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
import warnings
import numpy as np
import autograd
from scipy.sparse import spmatrix, csr_matrix, csc_matrix, coo_matrix


class System(object):
    """Description of the dynamical system.

    To define a dynamical system, we need to specify dimension of state, control, and parameter.
    Optionally, integration approach can be selected.
    This function should be inherited and users are supposed to override dyn/jac_dyn functions.

    """
    odes = ['RK4', 'Dis', 'Euler', 'BackEuler']
    fd_step = 1e-6
    def __init__(self, nx, nu, np=0, ode='Euler'):
        """Constructor for class.

        :param nx: int, dimension of state variable
        :param nu: int, dimension of control variable
        :param np: int, dimension of additional parameter
        :param ode: str, integration approach. Default RK4, options: Dis, Euler, BackEuler

        """
        self.nx = nx
        self.nu = nu
        self.np = np
        assert ode in System.odes
        self.ode = ode
        self._is_ad = self.__ad__.__func__ != System.__ad__
        if self._is_ad:
            self._jac_fun = autograd.jacobian(self.__ad__)

    def set_ode(self, method):
        """Set ode approach.

        :param method: str, name of ode approach.

        """
        assert method in System.odes
        self.ode = method

    def __ad__(self, txup):
        raise NotImplementedError

    def dyn(self, t, x, u, p=None):
        """Dynamics function without gradient information.

        It has to be overriden.

        :param t: float, time when evaluating system dynamics
        :param x: np.ndarray, (nx,) state variable
        :param u: np.ndarray, (nu,) control variable
        :param p: np.ndarray, (np,) additional optimizing variable
        :return: either dotx or x_k+1 depends on system type

        """
        raise NotImplementedError

    def jac_dyn(self, t, x, u, p=None):
        """Dynamics function with Jacobian return.

        It has to be overriden.

        :param t, x, u, p, h: see dyn
        :returns: y: ndarray, either dotx or x_k+1
        :returns: J: ndarray/spmatrix, returning Jacobian of this function evaluation

        """
        if self._is_ad:
            if self.ode == 'Dis':
                if p is not None:
                    xin = np.concatenate((x, u, p))
                else:
                    xin = np.concatenate((x, u))
            else:
                if p is not None:
                    xin = np.concatenate(([t], x, u, p))
                else:
                    xin = np.concatenate(([t], x, u))
            y0 = self.__ad__(xin)
            J = self._jac_fun(xin)
            return y0, J
        # implement a basic finite difference version
        if self.ode == 'Dis':
            J = np.zeros((self.nx, self.nx + self.nu + self.np))
            y0 = self.dyn(t, x, u, p)
            for i in range(self.nx):
                xp = x.copy()
                xp[i] += System.fd_step
                yi = self.dyn(t, xp, u, p)
                J[:, i] = (yi - y0) / System.fd_step
            # do it for u
            for i in range(self.nu):
                up = u.copy()
                up[i] += System.fd_step
                yi = self.dyn(t, x, up, p)
                J[:, self.nx + i] = (yi - y0) / System.fd_step
            # for p
            for i in range(self.np):
                pp = p.copy()
                pp[i] += System.fd_step
                yi = self.dyn(t, x, u, pp)
                J[:, self.nx + self.nu + i] = (yi - y0) / System.fd_step
            return y0, J
        raise NotImplementedError


class MultiStepDisSystem(System):
    """
    This class takes in a system and allow several steps of forward simulation
    
    :param system: The basic discrete system which is of type System and has ode == Dis
    :param k: int, the number of forward simulation steps
    """
    def __init__(self, system, k):
        assert isinstance(system, System) and system.ode == 'Dis'
        System.__init__(self, system.nx, nu, np=0, ode='Euler')
        self.system = system
        self.k = k
        self.ode = 'Dis'
        
    def dyn(self, t, x, u, p=None):
        """Dynamics function without gradient information.

        It has to be overriden.

        :param t: float, time when evaluating system dynamics
        :param x: np.ndarray, (nx,) state variable
        :param u: np.ndarray, (nu,) control variable
        :param p: np.ndarray, (np,) additional optimizing variable
        :return: either dotx or x_k+1 depends on system type

        """
        for i in range(self.k):
            x = system.dyn(i, x, u, p)
        return x

    def jac_dyn(self, t, x, u, p=None):
        """Dynamics function with Jacobian return.
        """
        J = np.c_[np.eye(self.nx), np.zeros((self.nx, self.nu + self.np))]  # this keeps the current jacobian estimation
        for i in range(self.k):
            x, Ji = system.jac_dyn(i, x, u, p)
            J[:, self.nx:] = Ji[:, self.nx:] + Ji[:, :self.nx].dot(J[:, self.nx:])
            J[:, :self.nx] = Ji[:, :self.nx].dot(J[:, :self.nx])
        return x, J


class MultiStepContSystem(System):
    """
    Given a continuous system and fixed step size, this class do forward integration for k steps.
    This greatly reduces the problem parameter numbers while keep integration accuracy higher.
    
    :param system: The basic continuous system which is of type System and has ode == 'Euler'
    :param h: float, the integration step size used in Euler integration
    :param k: int, the number of forward simulation steps, (hk) is the actual step size of optimized trajectory.
    """
    def __init__(self, system, h, k):
        self.system = system
        self.h = h
        self.k = k
        self.ode = 'Dis'
    
    def dyn(self, t, x, u, p=None):
        """
        Perform forward integration for k steps and return
        """
        for i in range(self.k):
            xdot = self.system.dyn(t + i * self.h, x, u, p)
            x = x + xdot * self.h
        return x
        
    def jac_dyn(self, t, x, u, p=None):
        """
        Perform forward integration for several steps, but with Jacobian information computed as well.
        """
        outj = np.zeros((self.nx, self.nx + self.nu + self.np))
        np.fill_diagonal(outj, 1)
        for i in range(self.k):
            xdot, jac = self.system.jac_dyn(-1, x, u, p)
            x = x + xdot * self.h
            tmp = np.eye(self.nx) + jac[:, 1: 1 + self.nx] * self.h
            outj[:, :self.nx] = tmp.dot(outj[:, :self.nx])
            outj[:, self.nx:] = tmp.dot(outj[:, self.nx:]) + self.h * jac[:, 1 + self.nx:]
        return x, outj
    


class DaeSystem(object):
    """A DAE system."""
    def __init__(self, nx, nu, np, nf, nG):
        r"""Constructor for the problem.

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
        if nx % nf != 0:
            warnings.warn(r"nx \% nf is not zero, make sure problem is defined correctly.")
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


class BaseFun(object):
    """Base class for functions, including both objective and constraint.

    This function should be inherited to define your own functions.
    A function with user supplied gradient information, nx, nf, ng should be set.

    """
    grad = ['user', 'no', 'ad', 'fd']
    def __init__(self, nx, nf, gradmode, ng=None):
        """Constructor for base function

        :param nx: int, number of variables input
        :param nf: int, number of response output
        :param gradmode: string, mode of gradient
        :param ng: int, used only when gradmode == 'user', means number of nnz gradients

        """
        self.nx = nx
        self.nf = nf
        if gradmode == 'user':
            self.nG = ng
            self.grad = True
        elif gradmode == 'no':
            self.grad = False
        elif gradmode == 'ad':
            self.grad = True
            self.nG = nx * nf
            self._jacobian = autograd.jacobian(self.__ad__)
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
        if self.gradmode == 'ad':  # the ad fun is called
            F[:] = self.__ad__(x)
            if needg:
                jac = self._jacobian(x)
                G[:] = jac.flat
                if rec:
                    row[:] = (np.arange(self.nf)[:, None] + np.zeros(self.nx)).flat
                    col[:] = (np.zeros(self.nf)[:, None] + np.arange(self.nx)).flat
        else:
            raise NotImplementedError

    def __ad__(self, x):
        """Function called to compute some stuff but with autodiff"""
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


class AddX(object):
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


class LinearObj(object):
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


class LinearPointObj(_objectWithMatrix):
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
            assert A.ndim == 1 and (xdim == len(A) or len(A) == xdim - 1)
            A = coo_matrix(A)
        assert A.shape[0] == 1 and (A.shape[1] == xdim or A.shape[1] == xdim - 1)
        self.A = A
        self.index = index


class LinearMultiPointObj(_objectWithMatrix):
    """Class for defining objective function that involves multiple points. It can be useful in cases where periodic motion is desired.
    Currently time is not supported.

    :param indexes: iterable, indexes of points we are considering
    :param As: iterable, np.ndarray or spmatrix, must of size equal to xdim
    :param nx: int, dimension of state
    :param nu: int, dimension of control
    :param np: int, dimension of parameter
    """
    def __init__(self, indexes, As, nx, nu, np_):
        self.nxup = (nx, nu, np_)
        self.idx_As = []
        xdim = 1 + nx + nu + np_  # length of the concatenated variable (t, x, u, p)
        for idx, A in zip(indexes, As):
            if isinstance(A, np.ndarray):
                assert A.ndim == 1 and xdim == len(A)
            A = coo_matrix(A)
            self.idx_As.append((idx, A))


class NonLinearObj(BaseFun):
    """Class for general nonlinear objective function over the entire decision variables.

    The objective function is basically calculated by calling a nonlinear function.

    """
    def __init__(self, nsol, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.

        :param nsol: int, length of decision variable, it can be set arbitrarily and won't be checked
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

        """
        BaseFun.__init__(self, nsol, 1, gradmode, nG)


class NonLinearMultiPointObj(BaseFun):
    """Class for defining point objective functions that require multiple points.

    The function is defined such that it can take a list of (t,x,u,p).

    """
    def __init__(self, indexes, nx, nu, np=0, gradmode='user', nG=None):
        """Constructor for nonlinear objective function.
        When defining __callf__, a list of points are passed as x
        When defining __callg__, a list of points are passed as x, the same applies to G, row, col

        :param index: int, at which point is objective calculated
        :param nx, nu, np: int, dimensions
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian, this has to be summed up for all points

        """
        self.indexes = indexes
        self.npoint = len(indexes)
        xdim = 1 + nx + nu + np
        BaseFun.__init__(self, xdim * self.npoint, 1, gradmode, nG)


class NonLinearPointObj(BaseFun):
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
        BaseFun.__init__(self, xdim, 1, gradmode, nG)
        self.index = index


class LqrObj(object):
    """Class for LQR objective since it is so common. It is treated independently with pathObj."""
    def __init__(self, F=None, Q=None, R=None, xfbase=None, xbase=None, ubase=None, tfweight=None, P=None, pbase=None):
        r"""Constructor for LQR objective function.

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


class QuadPenalty(NonLinearObj):
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
        NonLinearObj.__init__(self, -1, 'user', nG=len(indices))

    def __callg__(self, x, y, G, row, col, rec, needg):
        y[0] = np.sum(self.weights * x[self.indices] ** 2)
        if needg:
            G[:] = 2 * self.weights * x[self.indices]
            if rec:
                row[:] = 0
                col[:] = self.indices


class NonDiagLqrObj(object):
    """Class for LQR objective with non-diagonal entries"""
    # TODO: implement me
    pass


class NonLinearPointConstr(BaseFun):
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
        BaseFun.__init__(self, xdim, nc, gradmode, nG)
        self.index = index
        self.lb = lb
        self.ub = ub


class NonLinearMultiPointConstr(BaseFun):
    """Class for defining constraints where multiple points interacts.
    To define __call__, the inputs x, G, row, col now becomes list of appropriate sizes
    """
    def __init__(self, indexes, nc, nx, nu, np=0, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for nonlinear point constraint. Also serve as path constraint.

        :param indexes: iterable, at which point is objective calculated
        :param nc: int, dimension of constraint function
        :param nx, nu, np: int, dimensions
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

        """
        xdim = 1 + nx + nu + np
        BaseFun.__init__(self, xdim * len(indexes), nc, gradmode, nG)
        self.indexes = indexes
        self.npoint = len(indexes)
        self.lb = lb
        self.ub = ub
        self.create_buffer()

    def create_buffer(self):
        nG = self.nG
        if nG is not None and nG > 0:
            buffer = [(np.inf * np.ones(nG), np.zeros(nG, dtype=int), np.zeros(nG, dtype=int)) for _ in range(self.npoint)]
            self.buffer = list(zip(*buffer))

    def reset_buffer(self):
        for tmp in self.buffer[0]:
            tmp[:] = np.inf


class NonLinearConstr(BaseFun):
    """Class for defining constraint function in a general form."""
    def __init__(self, nsol, nc, lb=None, ub=None, gradmode='user', nG=None):
        """Constructor for general nonlinear constraint.

        :param nsol: int, length of the solution vector, used to initialize baseFun
        :param nc: int, dimension of constraint function
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param gradmode: str, how gradient is provided
        :param nG: int, number of nnz of Jacobian

        """
        BaseFun.__init__(self, nsol, nc, gradmode, nG)
        self.lb = lb
        self.ub = ub


class LinearPointConstr(_objectWithMatrix):
    """Class for linear constraint at selected points."""
    def __init__(self, index, A, lb=None, ub=None, offset=0):
        self.lb = lb
        self.ub = ub
        self.index = index
        self.A = coo_matrix(A)
        if offset:
            row, col = self.A.shape
            self.A = coo_matrix((self.A.data, (self.A.row, self.A.col + offset)),
                                shape=(row, col + offset))


class LinearMultiPointConstr(object):
    """Class for linear constraints that multiple points play a role"""
    def __init__(self, indexes, As, lb=None, ub=None):
        self.lb = lb
        self.ub = ub
        self.indexes = indexes
        self.As = [coo_matrix(A) for A in As]


class LinearConstr(object):
    """Class for linear constraints based on the whole x length."""
    def __init__(self, A, lb=None, ub=None, offset=0):
        self.A = coo_matrix(A)
        if offset:
            row, col = self.A.shape
            self.A = coo_matrix((self.A.data, (self.A.row, self.A.col + offset)), shape=(row, col + offset))
        self.lb = lb
        self.ub = ub
