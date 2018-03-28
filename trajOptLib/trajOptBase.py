#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
trajOptBase.py

Classes ready to be used for trajectory optimization.
"""
import sys, os, time
import numpy as np
from scipy.sparse import spmatrix, csr_matrix, csc_matrix, coo_matrix
import matplotlib.pyplot as plt
import logging
from .libsnopt import funBase, snoptConfig, result, probFun, solver


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class system(object):
    """Description of the dynamical system.

    To define a dynamical system, we need to specify dimension of state, control, and parameter.
    Optionally, integration approach can be selected.
    This function should be inherited and users are supposed to override dyn/Jdyn functions.

    """
    odes = ['RK4', 'Dis', 'Euler', 'BackEuler']
    def __init__(self, nx, nu, np=0, ode='RK4'):
        """Constructor for class.

        :param nx: int, dimension of state variable
        :param nu: int, dimension of control variable
        :param np: int, dimension of additional parameter
        :param ode: str, integration approach. Default RK4, options: Dis, Euler, BackEuler

        """
        self.nx = nx
        self.nu = nu
        self.np = np
        assert ode in system.odes
        self.ode = ode

    def setOde(self, method):
        """Set ode approach.

        :param method: str, name of ode approach.

        """
        assert method in system.odes
        self.ode = method

    def dyn(self, t, x, u, p=None, h=None):
        """Dynamics function without gradient information.

        It has to be overriden.

        :param t: float, time when evaluating system dynamics
        :param x: np.ndarray, (nx,) state variable
        :param u: np.ndarray, (nu,) control variable
        :param p: np.ndarray, (np,) additional optimizing variable
        :param h: float, used for discretized system. Integration step size
        :return: either dotx or x_k+1 depends on system type

        """
        raise NotImplementedError

    def Jdyn(self, t, x, u, p=None, h=None):
        """Dynamics function with Jacobian return.

        It has to be overriden.

        :param t, x, u, p, h: see dyn
        :returns: y: ndarray, either dotx or x_k+1
        :returns: J: ndarray/spmatrix, returning Jacobian of this function evaluation

        """
        raise NotImplementedError


class baseFun(funBase):
    """Base class for functions, including both objective and constraint.

    This function should be inherited to define your own functions.
    A function with user supplied gradient information, nx, nf, ng should be set.

    """
    grad = ['user', 'no']
    def __init__(self, nx, nf, gradmode, ng=None):
        """Constructor for base function

        :param nx: int, number of variables input
        :param nf: int, number of response output
        :param gradmode: string, mode of gradient
        :param ng: int, used only when gradmode == 'user', means number of nnz gradients

        """
        if gradmode == 'user':
            funBase.__init__(self, nx, nf, ng)
        elif gradmode == 'no':
            funBase.__init__(self, nx, nf)
        else:
            raise NotImplementedError
        self.gradmode = gradmode

    def __callf__(self, x, F):
        """Function call with no gradient information.

        :param x: ndarray, input to the function
        :param F: ndarray, output of the function which is written inplace

        """
        raise NotImplementedError

    def __callg__(self, x, F, G, row, col, rec, needg):
        """Function call with no gradient information.

        :param x: ndarray, input to the function
        :param F: ndarray, output of the function which is written inplace
        :param G: ndarray, gradient of the function in sparse form which stores the values
        :param row: ndarray of int, stores the rows of the sparse gradient matrix
        :param col: ndarray of int, stores the columns of the sparse gradient matrix
        :param rec: bool, determine if we write values to row and col
        :param needg: bool, determine if we need to calculate gradient and write to G

        """
        raise NotImplementedError


class linearObj(object):
    """Class for directly add linear objective function over the entire decision variable.

    It serves for objective of form :math:`y=Ax` where :math:`x` is the collected long vector.

    """
    def __init__(self, A):
        """Constructor for linear objective function using A

        :param A: np.ndarray or spmatrix, must of size equal to nsol

        """
        if isinstance(A, np.ndarray):
            assert A.ndim == 1
            A = coo_matrix(A)
        assert A.shape[0] == 1
        self.A = A


class linearPointObj(object):
    """Class for directly add linear objective function over the entire decision variable.

    It serves for objective function of the form :math:`y=Ax` where :math:`x` is the concatenated vector of state, 
    control and parameter at a selected index.

    """
    def __init__(self, index, A, nx, nu, np):
        """Constructor for linear objective function using A pointwise

        :param index: int, at which point is objective function evaluated
        :param A: np.ndarray or spmatrix, must of size equal to xdim
        :param nx: int, dimension of state
        :param nu: int, dimension of control
        :param np: int, dimension of parameter

        """
        xdim = 1 + nx + nu + np  # this x means length of variable at one point, (t, x, u, p)
        if isinstance(A, np.ndarray):
            assert A.ndim == 1 and xdim == len(A)
            A = csr_matrix(A)
        assert A.shape[0] == 1 and A.shape[1] == xdim
        self.A = A
        self.index = index


class nonLinObj(baseFun):
    """Class for general nonlinear objective function over the entire decision variables.

    The objective function is basically calculated by calling a nonlinear function.

    """
    def __init__(self, nsol, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.

        :param nsol: int, length of decision variable
        :param gradmode: str, how gradient is provided
        :param nG, int, number of nnz of Jacobian

        """
        baseFun.__init__(self, nsol, 1, gradmode, nG)


class nonPointObj(baseFun):
    """Class for defining point objective function.

    Similar to linear case. A function that takes the concatenated vector at a selected index is used.

    """
    def __init__(self, index, nx, nu, np=0, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.

        :param index: int, at which point is objective calculated
        :param nx, nu, np: int, dimensions
        :param gradmode: str, how gradient is provided
        :param nG, int, number of nnz of Jacobian

        """
        xdim = 1 + nx + nu + np
        baseFun.__init__(self, xdim, 1, gradmode, nG)
        self.index = index


class lqrObj(object):
    """Class for LQR objective since it is so common. It is treated independently with pathObj."""
    def __init__(self, F=None, Q=None, R=None, xfbase=None, xbase=None, ubase=None, tfweight=None, P=None, pbase=None):
        """Constructor for LQR objective function.

        :math:`c=\|x_f-x_{fbase}\|_F + \Sigma (\|x-x_{base}\|_Q + \|u-u_{base}\|_R + \|p-p_{base}\|_P) * h`

        :param F, Q, R: cost for terminal, path state, path ctrl.
        :param xfbase, xbase, ubase: the basis.
        :param P: might be None if we do not penalize p
        :param pbase: might be None if we do not penalize p

        """
        if F is not None:
            assert F.ndim == 1
            self.nx = len(F)
            self.F = coo_matrix(F)
        else:
            assert xfbase is None
            self.F = coo_matrix([])
        if Q is not None:
            assert Q.ndim == 1
            self.nx = len(Q)
            self.Q = coo_matrix(Q)
        else:
            assert xbase is None
            self.Q = coo_matrix([])
        if R is not None:
            self.nu = len(R)
            assert R.ndim == 1
            self.R = coo_matrix(R)
        else:
            assert ubase is None
            self.R = coo_matrix([])
        self.tfweight = tfweight
        if self.F.nnz > 0:
            self.xfbase = xfbase
            if self.xfbase is None:
                self.xfbase = np.zeros(self.nx)
        if self.Q.nnz > 0:
            self.xbase = xbase
            if self.xbase is None:
                self.xbase = np.zeros(self.nx)
        if self.R.nnz > 0:
            self.ubase = ubase
            if self.ubase is None:
                self.ubase = np.zeros(self.nu)
        if P is not None:
            assert P.ndim == 1
            self.P = coo_matrix(P)
            self.np = len(P)
            if pbase is None:
                self.pbase = np.zeros(self.np)
        else:
            self.P = None


class nonDiagLQRObj(object):
    """Class for LQR objective with non-diagonal entries"""
    # TODO: implement me
    pass


class pointConstr(baseFun):
    """Class for defining point constraint function."""
    def __init__(self, index, nc, nx, nu, np=0, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for nonlinear point constraint. Also serve as path constraint.

        :param index: int, at which point is objective calculated
        :param nc: int, dimension of constraint function
        :param nx, nu, np: int, dimensions
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param gradmode: str, how gradient is provided
        :param nG, int, number of nnz of Jacobian

        """
        xdim = 1 + nx + nu + np
        baseFun.__init__(self, xdim, nc, gradmode, nG)
        self.index = index
        self.lb = lb
        self.ub = ub


class nonLinConstr(baseFun):
    """Class for defining constraint function in a general form."""
    def __init__(self, nsol, nc, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for general nonlinear constraint.

        :param nsol: int, length of the solution vector, used to initialize baseFun
        :param nc: int, dimension of constraint function
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param gradmode: str, how gradient is provided
        :param nG, int, number of nnz of Jacobian

        """
        baseFun.__init__(self, nsol, nc, gradmode, nG)
        self.lb = lb
        self.ub = ub


def linearPointConstr(object):
    """Class for linear constraint at selected points."""
    def __init__(self, index, A, lb=None, ub=None):
        self.lb = lb
        self.ub = ub
        self.index = index
        self.A = A


def linearConstr(object):
    """Class for linear constraints based on the whole x length."""
    def __init__(self, A, lb=None, ub=None):
        self.A = A
        self.lb = lb
        self.ub = ub
