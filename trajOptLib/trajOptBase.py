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
import numpy as np
from scipy.sparse import spmatrix, csr_matrix, csc_matrix, coo_matrix
from .libsnopt import FunBase as funBase, SnoptConfig as snoptConfig, SnoptResult as result, probFun, solver


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


class daeSystem(object):
    """A DAE system."""
    def __init__(self, nx, nu, np, nf, nG):
        """Constructor for the problem.

        For a dae system described by :math `f(t, q, \dot{q}, \ddot{q}, u, p)=0`, nx=3dim(q), nf=dim(q).
        If it is described by :math `f(t, q, \dot{q}, u, p)=0`, nx=2dim(q), nf=dim(q).
        We define as order the highest time derivative of state q. Keeping this in mind, nf always equals dim(q),
        nx = (1+order)dim(q), nDefect = 2*order*nf
        Compared with system class, system dynamics can be given implicitly.
        For different order, we are gonna have different number of defect constraints and sizes.

        :param nx: int, dimension of states, it might also include acceleration
        :param nu: int, dimension of control
        :param np: int, dimension of parameter
        :param nf: int, dimension of dae system
        :param nG: int, nnz of Jacobian

        """
        self.nx = nx
        self.nu = nu
        self.np = np
        self.nf = nf
        self.nG = nG
        assert nx % nf == 0
        self.order = nx // nf - 1  # this is useful for detecting size
        self.autonomous = False
        self.timeindex = None

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        """Implementation of system dynamics expressed in a dae.

        It evaluates system dynamics like f(t, q, dq, ddq, p, u) = 0 and calculate gradients if necessary.
        :param t: float, time of evaluation
        :param x: ndarray, (nx,) state variable, it might contain q, dq, ddq or in more general case q and dq.
        :param u: ndarray, (nu,) control variable
        :param p: ndarray, (np,) parameter used such as reaction force from ground
        :param y: ndarray, (nq,) this constraint function. It is evaluated here.
        :param G, row, col: ndarray, (nG,) gradient of this constraints, row and col index
        :param rec: bool, if we need to write row and col
        :param needg: bool, if we have to evaluate gradients.
        """
        raise NotImplementedError

    def findTimeGradient(self, catx):
        """Detect if gradient is time related."""
        t = 0
        x = catx[:self.nx]
        u = catx[self.nx: self.nx + self.nu]
        p = catx[self.nx+self.nu: self.nx + self.nu + self.np]
        y = np.zeros(self.nf)
        G = np.zeros(self.nG)
        row = np.zeros(self.nG, dtype=int)
        col = np.zeros(self.nG, dtype=int)
        self.dyn(t, x, u, p, y, G, row, col, True, True)
        self.timeindex = np.where(col == 0)[0]
        if len(self.timeindex) == 0:
            self.autonomous = True


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
        self.timeindex = None
        self.autonomous = False

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

    def findTimeGradient(self, x):
        assert self.gradmode, 'Grad mode is off'
        tmpy = np.zeros(self.nf)
        tmpG = np.zeros(self.nG)
        tmprow = np.zeros(self.nG, dtype=int)
        tmpcol = np.zeros(self.nG, dtype=int)
        self.__callg__(x, tmpy, tmpG, tmprow, tmpcol, True, True)
        self.timeindex = np.where(tmpcol == 0)[0]  # the columns from time
        if len(self.timeindex) == 0:
            self.autonomous = True


class addX(object):
    """A description of additional optimizing parameter.

    It is intended to be used if the optimal control has like point constraint.
    In this class the user has to supply the size and bounds of those variables.
    """
    def __init__(self, n, lb=None, ub=None):
        """Constructor of this class.

        :param n: int, length of this variable.
        :param lb: ndarray, (n,) lower bounds for those variables. None means no bound
        :param ub: ndarray, (n,) uppere bounds for those variables. None means no bound

        """
        self.n = n
        if lb is None:
            self.lb = -1e20 * np.ones(n)
        else:
            self.lb = np.array(lb)
        if ub is None:
            self.ub = 1e20 * np.ones(n)
        else:
            self.ub = np.array(ub)


class _objectWithMatrix(object):
    """An abstract class that basically has a matrix in it."""
    def __init__(self):
        self.A = None
        self.timeindex = None
        self.autonomous = False

    def findTimeGradient(self):
        """For a matrix, find column 0 indice."""
        assert isinstance(self.A, coo_matrix)
        self.timeindex = np.where(self.A.col == 0)[0]  # the columns from time
        if len(self.timeindex) == 0:
            self.autonomous = True


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


class linearPointObj(_objectWithMatrix):
    """Class for directly add linear objective function over the entire decision variable.

    It serves for objective function of the form :math:`y=Ax` where :math:`x` is the concatenated vector of state, 
    control and parameter at a selected index.

    """
    def __init__(self, index, A, nx, nu, np_):
        """Constructor for linear objective function using A pointwise

        :param index: int, at which point is objective function evaluated
        :param A: np.ndarray or spmatrix, must of size equal to xdim
        :param nx: int, dimension of state
        :param nu: int, dimension of control
        :param np: int, dimension of parameter

        """
        xdim = 1 + nx + nu + np_  # this x means length of variable at one point, (t, x, u, p)
        if isinstance(A, np.ndarray):
            assert A.ndim == 1 and xdim == len(A)
            A = coo_matrix(A)
        assert A.shape[0] == 1 and A.shape[1] == xdim
        self.A = A
        self.index = index


class nonLinearObj(baseFun):
    """Class for general nonlinear objective function over the entire decision variables.

    The objective function is basically calculated by calling a nonlinear function.

    """
    def __init__(self, nsol, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.

        :param nsol: int, length of decision variable
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

        """
        baseFun.__init__(self, nsol, 1, gradmode, nG)


class nonLinearPointObj(baseFun):
    """Class for defining point objective function.

    Similar to linear case. A function that takes the concatenated vector at a selected index is used.

    """
    def __init__(self, index, nx, nu, np=0, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.

        :param index: int, at which point is objective calculated
        :param nx, nu, np: int, dimensions
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

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


class quadPenalty(nonLinearObj):
    """In many scenarios, we want to minimize the quadratic of some variables for some variables.

    This is generally different from LQR objective by that it is a point constraint and thus not integral one.
    To make it versatile, the user is allowed to pass indices so we can directly evaluate those variables.
    User friendly classes are also created so the indices are calculated internally.

    """
    def __init__(self, indices, weights):
        """Constructor for the class.

        :param indices: ndarray, indices of variables we aim to penalize.
        :param weights: float/ndarray, weights for terms

        """
        self.indices = indices
        self.weights = weights
        nonLinearObj.__init__(self, -1, 'user', nG=len(indices))

    def __callg__(self, x, y, G, row, col, rec, needg):
        y[0] = np.sum(self.weights * x[self.indices] ** 2)
        if needg:
            G[:] = 2 * self.weights * x[self.indices]
            if rec:
                row[:] = 0
                col[:] = self.indices


class nonDiagLQRObj(object):
    """Class for LQR objective with non-diagonal entries"""
    # TODO: implement me
    pass


class nonLinearPointConstr(baseFun):
    """Class for defining point constraint function."""
    def __init__(self, index, nc, nx, nu, np=0, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for nonlinear point constraint. Also serve as path constraint.

        :param index: int, at which point is objective calculated
        :param nc: int, dimension of constraint function
        :param nx, nu, np: int, dimensions
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

        """
        xdim = 1 + nx + nu + np
        baseFun.__init__(self, xdim, nc, gradmode, nG)
        self.index = index
        self.lb = lb
        self.ub = ub


class nonLinearConstr(baseFun):
    """Class for defining constraint function in a general form."""
    def __init__(self, nsol, nc, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for general nonlinear constraint.

        :param nsol: int, length of the solution vector, used to initialize baseFun
        :param nc: int, dimension of constraint function
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

        """
        baseFun.__init__(self, nsol, nc, gradmode, nG)
        self.lb = lb
        self.ub = ub


class linearPointConstr(_objectWithMatrix):
    """Class for linear constraint at selected points."""
    def __init__(self, index, A, lb=None, ub=None):
        self.lb = lb
        self.ub = ub
        self.index = index
        self.A = coo_matrix(A)


class linearConstr(object):
    """Class for linear constraints based on the whole x length."""
    def __init__(self, A, lb=None, ub=None):
        self.A = coo_matrix(A)
        self.lb = lb
        self.ub = ub
