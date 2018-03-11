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
from scipy.sparse import spmatrix, csr_matrix, csc_matrix
import matplotlib.pyplot as plt
import logging
from libsnopt import funBase, snoptConfig, result, probFun, solver


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class system(object):
    """Description of the dynamical system."""
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
        """Set ode approach"""
        assert method in system.odes
        self.ode = method

    def dyn(self, t, x, u, p=None, h=None):
        """Dynamics function that has to be overriden.
        :param t: float, time when evaluating system dynamics
        :param x: np.ndarray, (nx,) state variable
        :param u: np.ndarray, (nu,) control variable
        :param p: np.ndarray, (np,) additional optimizing variable
        :param h: float, used for discretized system. Integration step size
        :rtype y: either dotx or x_k+1
        """
        raise NotImplementedError

    def Jdyn(self, t, x, u, p=None, h=None):
        """Dynamics function with Jacobian return.
        :param t, x, u, p, h: see dyn
        :rtype y: ndarray, either dotx or x_k+1
        :rtype J: ndarray/spmatrix, returning Jacobian of this function evaluation
        """
        raise NotImplementedError


class baseFun(funBase):
    """Base class for functions, including both objective and constraint."""
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
        raise NotImplementedError

    def __callg__(self, x, F, G, row, col, rec, needg):
        raise NotImplementedError


class linearObj(baseFun):
    """Class for directly add linear objective function over the entire decision variable."""
    def __init__(self, A):
        """Constructor for linear objective function using A
        :param A: np.ndarray or spmatrix, must of size equal to nsol
        """
        if isinstance(A, np.ndarray):
            assert A.ndim == 1
            A = csr_matrix(A)
        assert A.shape[0] == 1
        nx = A.shape[1]
        nG = A.nnz
        self.A = A
        baseFun.__init__(self, nx, 1, 'user', nG)


class linearPointObj(baseFun):
    """Class for directly add linear objective function over the entire decision variable."""
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
        nx = A.shape[1]
        nG = A.nnz
        self.A = A
        self.index = index
        baseFun.__init__(self, xdim, 1, 'user', nG)


class nonLinObj(baseFun):
    """Class for general nonlinear objective function over the entire decision variables."""
    def __init__(self, nsol, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.
        :param nsol: int, length of decision variable
        :param gradmode: str, how gradient is provided
        :param nG, int, number of nnz of Jacobian
        """
        baseFun.__init__(self, nsol, 1, gradmode, nG)


class nonPointObj(baseFun):
    """Class for defining point objective function"""
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


class pointConstr(baseFun):
    """Class for defining point constraint function."""
    def __init__(self, index, nc, nx, nu, np=0, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for nonlinear point constraint.
        :param index: int, at which point is objective calculated
        :param nc: int, dimension of constraint function
        :param nx, nu, np: int, dimensions
        :param gradmode: str, how gradient is provided
        :param nG, int, number of nnz of Jacobian
        """
        xdim = 1 + nx + nu + np
        baseFun.__init__(self, xdim, nc, gradmode, nG)
        self.index = index
        self.lb = lb
        self.ub = ub
