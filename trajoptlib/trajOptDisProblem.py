#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Gao Tang <gaotang2@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
trajOptDisProblem.py

Class for describing the trajectory optimization problems.
But this is designed for discrete system and is handy when only a simulator is available.
"""
from __future__ import print_function, division
from collections.abc import Iterable
import numpy as np
from .trajOptBase import LinearObj, LinearPointObj, LinearMultiPointObj
from .trajOptBase import LinearPointConstr, LinearConstr, LinearMultiPointConstr
from .trajOptBase import NonLinearPointObj, NonLinearMultiPointObj, NonLinearObj
from .trajOptBase import NonLinearPointConstr, NonLinearConstr, NonLinearMultiPointConstr
from .trajOptBase import System, AddX as addX
from .trajOptBase import LqrObj
from pyoptsolver import OptProblem
from .utility import (random_gen_in_bound as randomGenInBound,
                      check_in_bounds as checkInBounds,
                      interp)
from scipy import sparse
from scipy.sparse import spmatrix, coo_matrix


class TrajOptDisProblem(OptProblem):
    """A class for definition of trajectory optimization problem.

    A general framework for using this class is to:

    1. Define a class inherited from system and write dyn/jac_dyn method.
    2. Optionally, write desired cost function by inheriting/creating from the list of available cost functions.
    3. Optionally, write desired constraint functions by inheriting from available constraints.
    4. Create this class with selected system, discretization, t0, tf range, gradient option
    5. Set bounds for state, control, parameters, x0 and xf
    6. Add objective functions and constraints to this class
    7. Call preProcess method explicitly
    8. Create snoptConfig instance and choose desired options
    9. Construct the solver
    10. Use the solver to solve with either automatic guess or user provided guess

    :currentmodule:
    .. exclude-members:: addLinearPointConstr
    .. exclude-members:: addLinearConstr
    .. exclude-members:: addNonLinearConstr
    .. exclude-members:: addNonLinearPointConstr
    .. exclude-members:: addNonLinearPointObj
    .. exclude-members:: addNonLinearObj
    .. exclude-members:: addLinearPointObj
    .. exclude-members:: addLinearObj

    """
    def __init__(self, sys, N, addx=None):
        """Initialize problem by system, discretization grid size, and allowable time

        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid, it has to be at least 2
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :param gradmode: bool, sets if we use gradient mode.
        :param addX: list of addX / one addX / None, additional optimization variables.

        """
        assert isinstance(sys, System) and sys.ode == 'Dis', "sys has to be System with Dis as ode type"
        self.sys = sys
        assert N >= 2 and isinstance(N, int), "N has to be integer and at least 2"
        self.N = N
        self.gradmode = True  # this controls if __callg__ will be called
        self.dimx = sys.nx  # I need N x, N - 1 u and N - 1 p
        self.dimu = sys.nu
        self.dimp = sys.np
        self.dimpoint = sys.nx + sys.nu + sys.np
        self.ubd = [None, None]
        self.xbd = [None, None]
        self.pbd = [None, None]
        self.x0bd = [None, None]
        self.xfbd = [None, None]
        # lqr cost function
        self.lqrObj = None
        # Linear cost function
        self.linearObj = []  # stores general linear cost
        self.linPointObj = []  # stores linear cost imposed at a point
        self.linPathObj = []  # stores Lagrange integral cost
        # nonlinear cost function
        self.nonLinObj = []  # stores general nonlinear cost
        self.nonPointObj = []  # stores nonlinear cost imposed at a point
        self.nonPathObj = []  # stores Lagrange integral cost. Includes LQR cost
        # nonlinear constraints. Linear constraints are treated as nonlinear
        self.pointConstr = []  # general constraint imposed at a certain point, such as initial and final point
        self.multiPointConstr = []
        self.pathConstr = []  # general constraint imposed everywhere such as collision avoidance
        self.nonLinConstr = []  # stores general nonlinear constraint
        self.linPointConstr = []
        self.linMultiPointConstr = []
        self.linPathConstr = []
        self.linearConstr = []
        # calculate number of variables to be optimized, time are always the last
        numX = self.N * self.dimx
        numU = self.N * self.dimu  # for convenience, I let U be there, even for the last one...
        numP = self.N * self.dimp  # the same applies for P
        if addx is None:
            self.lenAddX = 0
        else:
            if not isinstance(addx, list):
                addx = [addx]
            for tmp in addx:
                assert isinstance(tmp, addX)
            self.addX = addx
            self.lenAddX = sum([tmp.n for tmp in addx])
        numSol = numX + numU + numP
        self.numX = numX
        self.numU = numU
        self.numP = numP
        self.numTraj = numSol
        self.numSol = numSol + self.lenAddX

    def pre_process(self, *args):
        """Initialize the instances of probFun now we are ready.

        Call this function after the objectives and constraints have been set appropriately.
        It calculate the space required for SNOPT and allocates sparsity structure if necessary.

        """
        numDyn = self.dimx * (self.N - 1)  # constraints from system dynamics
        numC = 0
        for constr in self.pointConstr:
            numC += constr.nf
        for constr in self.multiPointConstr:
            numC += constr.nf
        for constr in self.pathConstr:
            numC += self.N * constr.nf
        for constr in self.nonLinConstr:
            numC += constr.nf
        nnonlincon = numC
        for constr in self.linPointConstr:
            numC += constr.A.shape[0]
        for constr in self.linMultiPointConstr:
            numC += constr.As[0].shape[0]
        for constr in self.linPathConstr:
            numC += constr.A.shape[0] * self.N  # this remains being argued
        for constr in self.linearConstr:
            numC += constr.A.shape[0]
        self.numF = 1 + numDyn + numC
        # time to analyze all objective functions in order to detect pattern for A, and additional variables for other nonlinear objective function
        spA, addn = self.__analyzeObj(self.numSol, self.numF)
        self.objaddn = addn  # this is important for multiple objective function support
        self.numSol += addn
        self.numF += addn
        OptProblem.__init__(self, self.numSol, self.numF)  # currently we do not know G yet
        self.__setAPattern(numDyn, nnonlincon, spA)
        self.__setXbound()
        self.__setFbound()
        # detect gradient information
        if self.gradmode:  # in this case, we randomly generate a guess and use it to initialize everything
            randX = self.randomGenX()
            self.__turnOnGrad(randX)

    def plot_jacobian(self, savefnm=None):
        """Plot the jacobian pattern of this problem.

        :param savefnm: str, the filename to save pattern into. If none, no figure is saved.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        if self.Acol.size > 0:
            ax.scatter(self.Acol, -self.Arow)
        randX = self.randomGenX()
        f, g, row, col = self.eval_g(randX)
        ax.scatter(col, -row)
        plt.show()
        if savefnm is not None:
            np.savez(savefnm, row=row, col=col, val=g, arow=self.Arow, acol=self.Acol, aval=self.Aval)

    def preProcess(self, *args):
        """Alias for pre_process"""
        self.pre_process(*args)

    def genGuessFromTraj(self, X=None, U=None, P=None, addx=None, obj=None, interp_kind='linear'):
        """Alias for gen_guess_from_traj"""
        return self.gen_guess_from_traj(X, U, P, addx, obj, interp_kind)

    def gen_guess_from_traj(self, X=None, U=None, P=None, addx=None, obj=None, interp_kind='linear'):
        """Generate an initial guess for the problem with user specified information.

        An amazing feature is the user does not have to give a solution of exactly the same time-stamped trajectory used internally.
        Interpolation approaches are used in such scenarios.
        The user is not required to provide them all, although we suggest they do so.

        :param X: ndarray, (x, x) each row corresponds to a state snapshot (if tstamp is None, assume equidistant grid). Column size can be dimx or dimx/sys.order
        :param U: ndarray, (x, dimu) each row corresponds to a control snapshot. Similar to X but with column size equals dimu
        :param P: ndarray, (x, dimp) each row corresponds to a parameter snapshot. Similar to X but with column size equals dimp
        :param t0/tf: float/array-like, initial/final time. If None, we randomly generate one
        :param addx: list of ndarray, guess of addx, if applicable
        :param tstamp: ndarray, (x,), None if the X/U/P are provided using equidistant grid.
        :param obj: ndarray, (x,) the objective part
        :param interp_kind: str, interpolation type for scipy.interpolate.interp1d, can be (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’)

        """
        randX = 2 * np.random.random(self.numSol) - 1
        Xtarget, Utarget, Ptarget = self.__parseX__(randX)
        if obj is not None:
            obj_ = self.__parseObj__(randX)
            obj_[:] = obj

        # interpolation for state variables
        nPoint = self.N
        dimx = self.dimx
        teval = np.arange(self.N)
        if X is not None:
            Xcol = X.shape[1]
            if not (Xcol == dimx):
                print('The column of X is not %d or %d, not use it' % (dimx))
                X = None
            else:  # use interpolation to do it
                if X.shape[0] == self.N:
                    Xtarget[:] = X
                else:
                    tstamp = np.linspace(0, self.N - 1, X.shape[0])
                    interp(tstamp, X, teval, Xtarget, interp_kind)
        else:
            # straight path go there
            for i in range(nPoint):
                Xtarget[i] = randomGenInBound(self.xbd, self.dimx)
            # randomly generate x0 and xf
            Xtarget[0, :dimx] = randomGenInBound(self.x0bd, self.dimx)
            Xtarget[-1, :dimx] = randomGenInBound(self.xfbd, self.dimx)

        # interpolation for control variable
        if U is not None:
            if U.shape[0] == self.N:
                Utarget[:] = U
            else:
                tstamp = np.linspace(0, self.N - 1, U.shape[0])
                interp(tstamp, U, teval, Utarget, interp_kind)
        else:
            for i in range(nPoint):
                Utarget[i] = randomGenInBound(self.ubd, self.dimu)
        if self.numP > 0:
            if P is not None:
                tstamp = np.linspace(0, self.N - 1, P.shape[0])
                interp(tstamp, P, teval, Ptarget, interp_kind)
            else:
                for i in range(nPoint):
                    Ptarget[i] = randomGenInBound(self.pbd, self.dimp)
        # generate for
        if self.lenAddX > 0:
            if addx is None:
                for field, addx_ in zip(self.__parseAddX__(randX), self.addX):
                    field[:] = randomGenInBound([addx_.lb, addx_.ub], addx_.n)
            else:
                for guess, field, addx_ in zip(addx, self.__parseAddX__(randX), self.addX):
                    field[:] = randomGenInBound(guess, addx_.n)
        # I do not have to worry about objaddn since they are linear
        return randX

    def __analyzeObj(self, numSol, numF):
        """Analyze the objective function.

        :param numSol: current estimation of free variables
        :param numF: current estimation of rows of constraints
        :return spA: coo sparse matrix, records first row of A, and last rows of A
        :return addn: int, additional nonlinear constraint. As long as nonlinear obj exists, this is non-zero

        """
        # detect how many nonlinear objective functions we have
        if self.lqrObj is None:
            nnlin = 0
        else:
            nnlin = 1
        nnlin += len(self.nonLinObj) + len(self.nonPathObj) + len(self.nonPointObj)
        addn = nnlin
        # analyze the linear objective functions in a good way
        A = np.zeros(numSol)
        for obj in self.linearObj:
            A[obj.A.col] += obj.A.data
        for obj in self.linPointObj:
            A[self.__patchCol__(obj.index, obj.A.col)] += obj.A.data
        for obj in self.linPathObj:  # this is not particularly useful, I have to say
            for i in range(self.numPoint):
                A[self.__patchCol__(i, obj.A.col)] += obj.A.data
        # get sparse representation of A
        nnzind = np.nonzero(A)[0]
        A_ = np.zeros(2 * addn)
        row_ = np.zeros(2 * addn, dtype=int)
        col_ = np.zeros(2 * addn, dtype=int)
        # for the addn
        for i in range(addn):
            A_[i] = 1
            A_[addn + i] = -1
            col_[i] = numSol + i
            col_[addn + i] = numSol + i
            row_[i] = 0
            row_[addn + i] = numF + i
        # concat them
        catA = np.concatenate((A[nnzind], A_))
        catArow = np.concatenate((np.zeros(len(nnzind), dtype=int), row_))
        catAcol = np.concatenate((nnzind, col_))
        spA = coo_matrix((catA, (catArow, catAcol)))
        return spA, addn

    def __setAPattern(self, ndyncon, nnonlincon, spA):
        """Set sparsity pattern for A. curRow is current row.

        It finds sparsity pattern from defect constraints.
        It also finds sparsity pattern from those linear constraints.

        :param ndyncon: int, describes how many dynamics constraints we have
        :param nnonlincon: int, describes how many nonlinear constraints we have
        :param spA: sparse matrix, how the objective function is described linearly.
        :return curRow: current row number after setting linear constraint
        :return curNA: current length of A

        """
        # find those linear constraints
        curRow = ndyncon + nnonlincon + 1  # since first row is objective
        lstCA = []
        lstCArow = []
        lstCAcol = []
        for constr in self.linPointConstr:
            lstCA.append(constr.A.data)
            lstCArow.append(constr.A.row + curRow)
            lstCAcol.append(self.__patchCol__(constr.index, constr.A.col))  # take care on here
            curRow += constr.A.shape[0]
        for constr in self.linMultiPointConstr:
            for idx, A in constr:
                lstCA.append(A.data)
                lstCArow.append(constr.A.row + curRow)
                lstCAcol.append(self.__patchCol__(idx, constr.A.col))
            curRow += constr.As[0].shape[0]
        for constr in self.linPathConstr:
            for index in range(self.N):
                lstCA.append(constr.A.data)
                lstCArow.append(constr.A.row + curRow)
                lstCAcol.append(self.__patchCol__(index, constr.A.col))
                curRow += constr.A.shape[0]
        for constr in self.linearConstr:
            lstCA.append(constr.A.data)
            lstCArow.append(constr.A.row + curRow)
            lstCAcol.append(constr.A.col)
            curRow += constr.A.shape[0]
        # for python3 compaticity, I extend the list
        lstCA.append(spA.data)
        lstCArow.append(spA.row)
        lstCAcol.append(spA.col)
        self.Aval = np.concatenate(lstCA)
        self.Arow = np.concatenate(lstCArow)
        self.Acol = np.concatenate(lstCAcol)
        self.set_a_by_triplet(self.Aval, self.Arow, self.Acol)
        curNA = len(self.Aval)  # this is just for bookkeeping
        return curRow, curNA

    def randomGenX(self):
        """Alias for :func:`~trajoptlib.TrajOptProblem.random_gen_guess`"""
        return self.random_gen_guess()

    def random_gen_guess(self):
        """A more reansonable approach to generate random guess for the problem.

        It considers bounds on initial and final states so this is satisfied.
        Then it linearly interpolate between states.
        Controls are randomly generated within control bound, if it presents. Otherwise [-1, 1]

        :return x: ndarray, (numSol, ) an initial guess of the solution
        """
        randX = np.zeros(self.numSol)
        X = np.reshape(randX[:self.numX], (self.N, self.dimx))
        U = np.reshape(randX[self.numX: self.numX + self.numU], (self.N, self.dimu))
        # randomly generate x0 and xf
        X[0] = randomGenInBound(self.x0bd, self.dimx)
        X[-1] = randomGenInBound(self.xfbd, self.dimx)
        # straight path go there
        for i in range(self.dimx):
            X[:, i] = np.linspace(X[0, i], X[-1, i], self.N)
        for i in range(self.N):
            U[i] = randomGenInBound(self.ubd, self.dimu)
        if self.numP > 0:
            P = np.reshape(randX[self.numX + self.numU: self.numX + self.numU + self.numP], (self.N, self.dimp))
            for i in range(self.N):
                P[i] = randomGenInBound(self.pbd, self.dimp)
        if self.lenAddX > 0:
            for field, addx in zip(self.__parseAddX__(randX), self.addX):
                field[:] = randomGenInBound([addx.lb, addx.ub], addx.n)
        return randX

    def __turnOnGrad(self, x0):
        """Turn on gradient, this is called after an initial x0 has been generated"""
        self.grad = True
        self.__getSparsity(x0)

    def __getSparsity(self, x0):
        """Detect sparsity of the problem with an initial guess."""
        numObjG = self.__getObjSparsity(x0)
        self.numObjG = numObjG
        # summarize number of pure linear constraints
        numDynG = self.__getDynSparsity(x0)
        numCG = 0  # G from C
        # I only care about those in numC
        for constr in self.pointConstr:
            numCG += constr.nG
        for constr in self.multiPointConstr:
            numCG += constr.nG
        for constr in self.pathConstr:
            numCG += self.N * constr.nG
        for constr in self.nonLinConstr:
            numCG += constr.nG
        numG = numObjG + numDynG + numCG
        self.numG = numG
        self.nG = numG

    def __getDynSparsity(self, x):
        """Set sparsity of the problem caused by system dynamics and other nonlinear constraints.

        :param x: ndarray, the guess/sol

        """
        # TODO: now sparsity is supported for Euler/ BackEuler/ Dis, RK4 needs to be done
        useX, useU, useP = self.__parseX__(x)
        # use one time/state/ctrl/para/h set to detect the sparsity pattern
        _, Jac = self.sys.jac_dyn(0, useX[0], useU[0], useP[0])
        nrow, ncol = Jac.shape
        assert nrow == self.dimx and ncol == self.dimpoint, "The shape of the dynamics Jacobian is incorrect"
        if isinstance(Jac, np.ndarray):
            nnz = nrow * ncol
        elif isinstance(Jac, spmatrix):
            nnz = Jac.nnz
        else:
            raise Exception("Jacobian from System.jac_dyn has to be np.ndarray or spmatrix")
        return (self.N - 1) * (nnz + self.dimx)

    def __getObjSparsity(self, x):
        """Set sparsity structure of the problem from objective function.

        The objective function pattern is composed of two parts:
        - linear parts. We sum all the coefficients and find sparsity pattern out of it
        - nonlinear parts. Each nonlinear objective is augmented with another row in jacobian
        and a another auxiliary optimization variable s.t. c(x) = y and J += y

        :param x: ndarray, the guess/sol
        :returns: nG: int, # Jacobian from nonlinear objective function

        """
        nG = 0
        # check sparseObj mode
        if self.lqrObj is not None:
            nG += self.LQRnG
        for obj in self.nonLinObj:
            nG += obj.nG
        for obj in self.nonPointObj:
            nG += obj.nG
        for obj in self.nonPathObj:
            nG += (self.N - 1) * obj.nG
        return nG

    def __setXbound(self):
        """Set bounds on decision variables."""
        # create bound on x
        dimpnt = self.dimpoint
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        xlb = np.zeros(self.numSol)
        xub = np.zeros(self.numSol)
        numDynVar = self.N * dimpnt
        Mxlb = np.reshape(xlb[:numDynVar], (self.N, dimpnt))
        Mxub = np.reshape(xub[:numDynVar], (self.N, dimpnt))
        Mulb = Mxlb[:, dimx:dimx + dimu]
        Muub = Mxub[:, dimx:dimx + dimu]
        Mplb = Mxlb[:, dimpnt - dimp:dimpnt]
        Mpub = Mxub[:, dimpnt - dimp:dimpnt]
        # set bounds for q and dq, agree with previous convention
        if self.xbd[0] is not None:
            Mxlb[:, :dimx] = self.xbd[0]
        else:
            Mxlb[:, :dimx] = -1e20
        # set lb for x0 and xf
        if self.x0bd[0] is not None:
            Mxlb[0, :dimx] = self.x0bd[0]
        if self.xfbd[0] is not None:
            Mxlb[-1, :dimx] = self.xfbd[0]

        if self.xbd[1] is not None:
            Mxub[:, :dimx] = self.xbd[1]
        else:
            Mxub[:, :dimx] = 1e20
        # set ub for x0 and xf
        if self.x0bd[1] is not None:
            Mxub[0, :dimx] = self.x0bd[1]
        if self.xfbd[1] is not None:
            Mxub[-1, :dimx] = self.xfbd[1]

        # set bounds for control variable
        if self.ubd[0] is not None:
            Mulb[:] = self.ubd[0]
        else:
            Mulb[:] = -1e20
        if self.ubd[1] is not None:
            Muub[:] = self.ubd[1]
        else:
            Muub[:] = 1e20
        if self.pbd[0] is not None and self.dimp > 0:
            Mplb[:] = self.pbd[0]
        else:
            Mplb[:] = -1e20
        if self.pbd[1] is not None and self.dimp > 0:
            Mpub[:] = self.pbd[1]
        else:
            Mpub[:] = 1e20

        # set bound on addX
        if self.lenAddX != 0:
            curN = self.numTraj
            for addx in self.addX:
                xlb[curN: curN + addx.n] = addx.lb
                xub[curN: curN + addx.n] = addx.ub

        # set bound on objaddn, this is obvious
        xlb[-self.objaddn:] = -1e20
        xub[-self.objaddn:] = 1e20

        # assign to where it should belong to
        self.set_xlb(xlb)
        self.set_xub(xub)

    def __setFbound(self):
        """Set bound on F"""
        # set bound on F
        numF = self.numF
        numDyn = self.dimx * (self.N - 1)  # constraints from system dynamics
        clb = np.zeros(numF)
        cub = np.zeros(numF)
        clb[0] = -1e20
        cub[0] = 1e20
        cind0 = 1 + numDyn
        for constr in self.pointConstr:
            if constr.lb is not None:
                clb[cind0: cind0 + constr.nf] = constr.lb
            if constr.ub is not None:
                cub[cind0: cind0 + constr.nf] = constr.ub
            cind0 += constr.nf
        for constr in self.multiPointConstr:
            if constr.lb is not None:
                clb[cind0: cind0 + constr.nf] = constr.lb
            if constr.ub is not None:
                cub[cind0: cind0 + constr.nf] = constr.ub
            cind0 += constr.nf
        for constr in self.pathConstr:
            tmplb = np.reshape(clb[cind0: cind0 + constr.nf * self.N], (self.N, constr.nf))
            tmpub = np.reshape(cub[cind0: cind0 + constr.nf * self.N], (self.N, constr.nf))
            cind0 += constr.nf * self.N
            if constr.lb is not None:
                tmplb[:] = constr.lb
            if constr.ub is not None:
                tmpub[:] = constr.ub
        for constr in self.nonLinConstr:
            if constr.lb is not None:
                clb[cind0: cind0 + constr.nf] = constr.lb
            if constr.ub is not None:
                cub[cind0: cind0 + constr.nf] = constr.ub
            cind0 += constr.nf
        # the rest are linear constraints and we should write those bounds, too
        for constr in self.linPointConstr:
            cindf = cind0 + constr.A.shape[0]
            clb[cind0: cindf] = constr.lb
            cub[cind0: cindf] = constr.ub
            cind0 = cindf
        for constr in self.linMultiPointConstr:
            cindf = cind0 + constr.As[0].shape[0]
            clb[cind0: cindf] = constr.lb
            cub[cind0: cindf] = constr.ub
            cind0 = cindf
        for constr in self.linPathConstr:
            cindf = cind0 + self.N * constr.A.shape[0]
            clb[cind0: cindf] = np.tile(constr.lb, (self.N, 1)).flatten()
            cub[cind0: cindf] = np.tile(constr.ub, (self.N, 1)).flatten()
            cind0 = cindf
        for constr in self.linearConstr:
            cindf = cind0 + constr.A.shape[0]
            clb[cind0: cindf] = constr.lb
            cub[cind0: cindf] = constr.ub
            cind0 = cindf
        # the bounds for objaddn is 0 so we are good so far
        # assign to where it should belong to
        self.lb = clb
        self.ub = cub

    def __parseX__(self, x):
        """Parse guess/sol into X, U, P"""
        X = np.reshape(x[:self.N * self.dimpoint], (self.N, self.dimpoint))
        useX = X[:, :self.dimx]
        useU = X[:, self.dimx:self.dimpoint - self.dimp]
        useP = X[:, self.dimpoint - self.dimp:]
        return useX, useU, useP

    def parse_f(self, guess):
        """Give an guess, evaluate it and parse into parts.

        :param guess: ndarray, (numSol, ) a guess or a solution to check
        :returns: dict, containing objective and parsed constraints

        """
        assert len(guess) == self.numSol
        N = self.N
        dimx = self.dimx
        y = np.ones(self.numF)
        if self.gradmode:
            self.__callg__(guess, y, np.zeros(1), np.zeros(1), np.zeros(1), False, False)
        else:
            self.__callf__(guess, y)
        obj = y[0]
        dynCon = np.reshape(y[1:(N - 1) * dimx + 1], (N - 1, dimx))
        curN = 1 + (N - 1) * dimx
        pointCon = []
        for constr in self.pointConstr:
            pointCon.append(y[curN: curN + constr.nf])
            curN += constr.nf
        multiPointCon = []
        for constr in self.multiPointConstr:
            multiPointCon.append(y[curN: curN + constr.nf])
            curN += constr.nf
        pathCon = []
        for constr in self.pathConstr:
            pathCon.append(np.reshape(y[curN: curN + N * constr.nf], (N, constr.nf)))
            curN += N * constr.nf
        nonLinCon = []
        for constr in self.nonLinConstr:
            nonLinCon.append(y[curN: curN + constr.nf])
            curN += constr.nf
        # check bounds, return a -1, 1 value for non-equality bounds, and 0 for equality bounds
        useX, useU, useP = self.__parseX__(guess)
        Xbound = checkInBounds(useX, self.xbd)
        x0bound = checkInBounds(useX[0], self.x0bd)
        xfbound = checkInBounds(useX[-1], self.xfbd)
        ubound = checkInBounds(useU, self.ubd)
        if self.dimp > 0:
            pbound = checkInBounds(useP, self.pbd)
        else:
            pbound = None
        if self.lenAddX > 0:
            addx = self.__parseAddX__(guess)
            addXbound = [checkInBounds(addx_, [addx__.lb, addx__.ub]) for addx_, addx__ in zip(addx, self.addX)]
        else:
            addXbound = None
        result = {'obj': obj, 'dyn': dynCon, 'Xbd': Xbound, 'Ubd': ubound, 'x0bd': x0bound, 'xfbd': xfbound, 'Pbd': pbound, 'addXbd': addXbound}
        if pointCon:
            result['point'] = pointCon
        if multiPointCon:
            result['mpoint'] = multiPointCon
        if pathCon:
            result['path'] = pathCon
        if nonLinCon:
            result['nonlin'] = nonLinCon
        return result

    def parse_sol(self, sol):
        """Call parseX function from utility and return a dict of solution.

        :param sol: ndarray, the solution.
        """
        X, U, P = self.__parseX__(sol)
        if self.dimp == 0:
            P = None
        obj = self.__parseObj__(sol)
        if self.lenAddX == 0:
            return {'t': np.arange(self.N), 'x': X, 'u': U, 'p': P, 'obj': obj}
        else:
            return {'t': np.arange(self.N), 'x': X, 'u': U, 'p': P, 'addx': self.__parseAddX__(sol), 'obj': obj}

    def __parseObj__(self, x):
        return x[self.numSol - self.objaddn:]

    def __parseAddX__(self, x):
        numTraj = self.numTraj
        addX = []
        for addx in self.addX:
            addX.append(x[numTraj: numTraj + addx.n])
            numTraj += addx.n
        return addX

    def __callf__(self, x, y):  # TODO: remove callf case
        """Evaluate those constraints and objective functions."""
        useX, useU, useP = self.__parseX__(x)
        # evaluate objective function
        self.__objModeF__(0, useX, useU, useP, x, y)
        # evaluate the system dynamics constraint
        curRow = 1
        curRow = self.__dynconstrModeF__(curRow, useX, useU, useP, y)
        # evaluate other constraints
        curRow = self.__constrModeF__(curRow, useX, useU, useP, x, y)
        return curRow

    def __objModeF__(self, curRow, useX, useU, useP, x, y):
        """Calculate objective function. F mode

        :param curRow: int, index from which we write on
        :param h, useT, useX, useU, useP: parsed solution
        :param x: ndarray, the sol, it is used for linear constraints
        :param y: ndarray, the F to be written. The first row stores the objective function

        """
        y[0] = 0.0
        tmpout = np.zeros(1)
        for obj in self.linPointObj:
            tmpx = np.concatenate((useX[obj.index], useU[obj.index], useP[obj.index]))
            y[0] += obj.A.dot(tmpx)
        for obj in self.linPathObj:
            for i in range(self.N - 1):
                tmpx = np.concatenate((useX[i], useU[i], useP[i]))
                obj.__callf__(tmpx, tmpout)
                y[0] += tmpout[0]
        for obj in self.linearObj:
            y[0] += obj.A.dot(x)
        for obj in self.nonPointObj:
            tmpx = np.concatenate((useX[obj.index], useU[obj.index], useP[obj.index]))
            obj.__callf__(tmpx, tmpout)
            y[0] += tmpout[0]
        for obj in self.nonPathObj:
            for i in range(self.N - 1):
                tmpx = np.concatenate((useX[i], useU[i], useP[i]))
                obj.__callf__(tmpx, tmpout)
                y[0] += tmpout[0]
        for obj in self.nonLinObj:
            if isinstance(obj, NonLinearObj):
                obj.__callf__(x, tmpout)
            else:  # NonLinearMultiPointObj
                xins = [np.concatenate((useX[idx], useU[idx], useP[idx])) for idx in obj.indexes]
                obj.__callf__(xins, tmpout)
            y[0] += tmpout[0]
        # add lqr cost, if applicable
        if self.lqrObj is not None:
            y[0] += self.lqrObj(useX, useU, useP)

    def __constrModeF__(self, curRow, useX, useU, useP, x, y):
        """Calculate constraint function. F mode

        :param curRow: int, index from which we write on
        :param h, useT, useX, useU, useP: parsed solution
        :param y: ndarray, the F to be written
        :returns: curRow: current row after we write on y

        """
        for constr in self.pointConstr:
            tmpx = np.concatenate((useX[constr.index], useU[constr.index], useP[constr.index]))
            constr.__evalf__(tmpx, y[curRow: curRow + constr.nf])
            curRow += constr.nf
        for constr in self.multiPointConstr:
            xin = [np.concatenate((useX[idx], useU[idx], useP[idx])) for idx in constr.indexes]
            constr.__evalf__(xin, y[curRow: curRow + constr.nf])
            curRow += constr.nf
        for constr in self.pathConstr:
            for i in range(self.N):
                tmpx = np.concatenate((useX[i], useU[i], useP[i]))
                constr.__evalf__(tmpx, y[curRow: curRow + constr.nf])
            self.numF
            curRow += constr.nf
        for constr in self.nonLinConstr:
            constr.__evalf__(x, y[curRow: curRow + constr.nf])
            curRow += constr.nf
        return curRow

    def __dynconstrModeF__(self, curRow, useX, useU, useP, y):
        """Calculate constraint from dynamics. F mode.

        :param curRow: int, index from which we write on
        :param h, useT, useX, useU, useP: the parsed sol
        :param y: ndarray, the F
        :returns: curRow: current row after we write on y

        """
        # loop over all system dynamics constraint
        cDyn = np.reshape(y[curRow:curRow + (self.N - 1) * self.dimx], (self.N - 1, self.dimx))
        for i in range(self.N - 1):
            # evaluate gradient of system dynamics
            ynext = self.sys.dyn(i, useX[i], useU[i], useP[i])
            cDyn[i] = ynext - useX[i + 1]
        return curRow

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix.

        :param x: ndarray, the solution to the problem
        :param y: ndarray, return F
        :param G, row, col: ndarray, information of gradient
        :param rec, needg: if we record/ if we need gradient

        """
        useX, useU, useP = self.__parseX__(x)
        # loop over all system dynamics constraint
        curRow = 1
        curNg = 0
        curRow, curNg = self.__dynconstrModeG(curRow, curNg, useX, useU, useP, y, G, row, col, rec, needg)
        curRow, curNg = self.__constrModeG(curRow, curNg, useX, useU, useP, x, y, G, row, col, rec, needg)
        curRow, curNg = self.__objModeG(curRow, curNg, useX, useU, useP, x, y, G, row, col, rec, needg)
        return curRow, curNg

    def __dynconstrModeG(self, curRow, curNg, useX, useU, useP, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        cDyn = np.reshape(y[curRow: curRow + (self.N - 1) * self.dimx], (self.N - 1, self.dimx))
        for i in range(self.N - 1):
            # evaluate gradient of system dynamics TODO: support other types of integration scheme
            ynext, Jac = self.sys.jac_dyn(i, useX[i], useU[i], useP[i])  # TODO: support in-place Jacobian
            cDyn[i] = ynext - useX[i + 1]
            if needg:
                if isinstance(Jac, np.ndarray):  # non-sparse mode
                    baseCol = i * self.dimpoint
                    # assign a block for x
                    G[curNg: curNg + dimx * self.dimpoint] = Jac.flatten()
                    if rec:
                        row[curNg: curNg + Jac.size] = curRow + (np.arange(dimx)[:, None] + np.zeros(self.dimpoint)).flatten()
                        col[curNg: curNg + Jac.size] = baseCol + (np.zeros(self.dimx)[:, None] + np.arange(self.dimpoint)).flatten()
                    curNg += Jac.size
                    # assign the diagonal block for x_{k+1}
                    G[curNg: curNg + dimx] = -1.0
                    if rec:
                        row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                        col[curNg: curNg + dimx] = baseCol + self.dimpoint + np.arange(dimx)
                    curNg += dimx
                else:
                    Jac = Jac.tocoo()  # what if some elements are eaten? Ignore it for now
                    rows = Jac.row
                    cols = Jac.col
                    data = Jac.data
                    G[curNg: curNg + Jac.nnz] = data
                    if rec:
                        row[curNg: curNg + Jac.nnz] = curRow + rows
                        col[curNg: curNg + Jac.nnz] = self.__patchCol__(i, cols)
                    curNg += Jac.nnz
                    # for x_{k+1}
                    G[curNg: curNg + dimx] = -1.0
                    if rec:
                        row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                        col[curNg: curNg + dimx] = np.arange(dimx) + (i + 1) * self.dimpoint
                    curNg += dimx
            curRow += dimx
        return curRow, curNg

    def __objModeG(self, curRow, curNg,  useX, useU, useP, x, y, G, row, col, rec, needg, first_row=False):
        """Calculate objective function. It just evaluates them and assign to correct position in y.

        See __constr_mode_g__ for arguments and output.
        """
        tmpout = np.zeros(1)
        y[0] = 0
        curRow = self.numF - self.objaddn
        if self.lqrObj is not None:  # the lqr obj
            Gpiece = G[curNg: curNg + self.LQRnG]
            rowpiece = row[curNg: curNg + self.LQRnG]
            colpiece = col[curNg: curNg + self.LQRnG]
            self.lqrObj(useX, useU, useP, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
            y[curRow] = tmpout[0]
            if rec:
                rowpiece[:] = curRow
            curNg += self.LQRnG
            curRow += 1

        # still in the point, path, obj order
        if self.nonPointObj:
            for obj in self.nonPointObj:
                tmpx = np.concatenate((useX[obj.index], useU[obj.index], useP[obj.index]))
                Gpiece = G[curNg: curNg + obj.nG]
                rowpiece = row[curNg: curNg + obj.nG]
                colpiece = col[curNg: curNg + obj.nG]
                obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                y[curRow] = tmpout[0]
                if rec:
                    rowpiece[:] = curRow
                    colpiece[:] = self.__patchCol__(obj.index, colpiece)
                curNg += obj.nG
                curRow += 1

        if self.nonPathObj:
            for obj in self.nonPathObj:
                y[curRow] = 0
                for i in range(self.N - 1):
                    tmpx = np.concatenate((useX[i], useU[i], useP[i]))
                    Gpiece = G[curNg: curNg + obj.nG]
                    rowpiece = row[curNg: curNg + obj.nG]
                    colpiece = col[curNg: curNg + obj.nG]
                    obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                    y[curRow] += tmpout[0]
                    if rec:
                        rowpiece[:] = curRow
                        colpiece[:] = self.__patchCol__(i, colpiece)
                    curNg += obj.nG
                curRow += 1

        if self.nonLinObj:
            for obj in self.nonLinObj:  # nonlinear cost function
                if isinstance(obj, NonLinearObj):
                    Gpiece = G[curNg: curNg + obj.nG]
                    rowpiece = row[curNg: curNg + obj.nG]
                    colpiece = col[curNg: curNg + obj.nG]
                    obj.__callg__(x, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                    y[curRow] = tmpout[0]
                    if rec:
                        rowpiece[:] = curRow
                    curNg += obj.nG
                    curRow += 1
                else:  # has to be NonLinearMultiPointObj
                    Gpiece = G[curNg: curNg + obj.nG]
                    rowpiece = row[curNg: curNg + obj.nG]
                    colpiece = col[curNg: curNg + obj.nG]
                    xins = [np.concatenate((useX[idx], useU[idx], useP[idx])) for idx in obj.indexes]
                    Gs = rows = cols = []
                    if needg:
                        Gs = [np.inf * np.ones(obj.nG) for _ in range(obj.npoint)]
                        if rec:
                            rows = [-1 * np.ones(obj.nG, dtype=int) for _ in range(obj.npoint)]
                            cols = [-1 * np.ones(obj.nG, dtype=int) for _ in range(obj.npoint)]
                    obj.__callg__(xins, tmpout, Gs, rows, cols, rec, needg)
                    if needg:
                        idx_ = 0
                        for i, idx in enumerate(obj.indexes):
                            good_len = np.sum(~np.isinf(Gs[i]))
                            Gpiece[idx_: idx_ + good_len] = Gs[i][:good_len]
                            if rec:
                                rowpiece[idx_: idx_ + good_len] = rows[i][:good_len]
                                colpiece[idx_: idx_ + good_len] = cols[i][:good_len] + idx * self.dimpoint
                            idx_ += good_len
        return curRow, curNg

    def __constrModeG(self, curRow, curNg, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate constraint function. G mode

        :param curRow: int, index from which we write on
        :param curNg: int, current index in G
        :param h, useT, useX, useU, useP: parsed solution
        :param y: ndarray, the F to be written
        :param G, row, col: ndarray, the G to be written and the locations
        :param rec: bool, if we record row and col
        :param needg: bool, if we need gradient information
        :return curRow: current row after we write on y
        :return curNg: current index in G after this

        """
        # loop over other constraints
        if self.pointConstr:
            for constr in self.pointConstr:
                tmpx = np.concatenate((useX[constr.index], useU[constr.index], useP[constr.index]))
                pieceG = G[curNg: curNg + constr.nG]
                pieceRow = row[curNg: curNg + constr.nG]
                pieceCol = col[curNg: curNg + constr.nG]
                constr.__callg__(tmpx, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
                if rec:
                    pieceRow += curRow
                    pieceCol[:] = self.__patchCol__(constr.index, pieceCol)
                curRow += constr.nf
                curNg += constr.nG
        if self.multiPointConstr:
            for constr in self.multiPointConstr:
                xin = [np.concatenate((useU[idx], useP[idx])) for idx in constr.indexes]
                constr.reset_buffer()
                constr.__callg__(xin, y[curRow: curRow + constr.nf], constr.buffer[0], constr.buffer[1], constr.buffer[2], rec, needg)
                if needg:
                    idx_ = 0
                    for i in range(constr.npoint):
                        leng = np.sum(~np.isinf(constr.buffer[0][i]))
                        pieceG[idx_: idx_ + leng] = constr.buffer[0][i][:leng]
                        if rec:
                            pieceRow[idx_: idx_ + leng] = curRow + constr.buffer[1][i][:leng]
                            pieceCol[idx_: idx_ + leng] = self.__patchCol__(self.indexes[i], constr.buffer[2][i][:leng])
                        idx_ += leng
                curRow += constr.nf
                curNg += constr.nG
        if self.pathConstr:
            for constr in self.pathConstr:
                for i in range(self.N):
                    tmpx = np.concatenate((useX[i], useU[i], useP[i]))
                    pieceG = G[curNg: curNg + constr.nG]
                    pieceRow = row[curNg: curNg + constr.nG]
                    pieceCol = col[curNg: curNg + constr.nG]
                    constr.__callg__(tmpx, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
                    if rec:
                        pieceRow += curRow
                        pieceCol[:] = self.__patchCol__(i, pieceCol)
                    curRow += constr.nf
                    curNg += constr.nG
        if self.nonLinConstr:
            for constr in self.nonLinConstr:
                pieceG = G[curNg: curNg + constr.nG]
                pieceRow = row[curNg: curNg + constr.nG]
                pieceCol = col[curNg: curNg + constr.nG]
                constr.__callg__(x, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
                if rec:
                    pieceRow += curRow
                curRow += constr.nf
                curNg += constr.nG
        return curRow, curNg

    # interface functions for ipopt
    def __cost__(self, x):
        """The eval_f function required by ipopt.

        :param x: a guess/solution of the problem
        :return f: float, objective function

        """
        y = np.zeros(1)
        useX, useU, useP = self.__parseX__(x)
        self.__objModeF__(0, useX, useU, useP, x, y)
        return y[0]

    def __gradient__(self, x, g):
        """Evaluation of the gradient of objective function.

        :param x: guess/solution to the problem
        :param g: the gradient of objective function w.r.t x to be written into

        """
        g[:] = 0
        mask = self.Arow == 0
        g[self.Acol[mask]] = self.Aval[mask]
        return True

    def __constraint__(self, x, f):
        """Evaluate constraint function.

        :param x: guess/solution to the problem
        :param f: constraints ready to be written upon
        """
        G = np.zeros(1)
        row = np.zeros(1, dtype=int)
        col = np.zeros(1, dtype=int)
        self.__callg__(x, f, G, row, col, False, False)
        return 0

    def __jacobian__(self, x, g, row, col, rec):
        """Evaluate the Jacobian of the problem.

        :param x: guess/solution to the problem
        :param g: the vector being written on for Jacobian entries
        """
        y = np.zeros(self.numF)
        self.__callg__(x, y, g, row, col, rec, True)
        return 1

    def __patchCol__(self, index, col, col_offset=0):
        """Find which indices it belongs to the original one for a local matrix at col"""
        if index < 0:
            index = index + self.N
        return col + col_offset + index * self.dimpoint

    def addLQRObj(self, lqrobj):
        """Alias for :func:`~trajOptLib.TrajOptProblem.add_lqr_obj`"""
        self.add_lqr_obj(lqrobj)

    def add_lqr_obj(self, lqrobj):
        """Add a lqr objective function to the problem. It changes lqrObj into a function being called in two modes...

        :param lqrobj: a lqrObj class.

        """
        Fcol = lqrobj.F.col
        Qcol = lqrobj.Q.col
        Rcol = lqrobj.R.col
        useF = len(Fcol)
        useQ = len(Qcol)
        useR = len(Rcol)
        if lqrobj.P is not None:
            Pcol = lqrobj.P.col
            useP = len(Pcol)
        else:
            useP = 0
            Pcol = []
        if lqrobj.tfweight is not None:
            tfweight = lqrobj.tfweight
        else:
            tfweight = 0.0
        self.LQRnG = (self.N - 1) * (useQ + useR + useP) + useF
        num1 = (self.N - 1) * (useQ + useR + useP)
        baseCol = self.dimpoint * np.arange(self.N - 1)[:, np.newaxis]  # a nPoint by 1 column matrix

        def __callf__(useX, useU, useP_):
            """Calculate the lqr cost.

            :param h: float, grid size
            :param useX: ndarray, parsed X
            :param useU: ndarray, parsed U
            :param useP: ndarray, parsed P, might be empty

            """
            y = 0
            if useF > 0:
                y += np.sum(lqrobj.F.data * ((useX[-1, Fcol] - lqrobj.xfbase[Fcol]) ** 2))
            if useQ > 0:
                y += np.sum(lqrobj.Q.data * np.sum((useX[:-1, Qcol] - lqrobj.xbase[Qcol]) ** 2, axis=0))
            if useR > 0:
                y += np.sum(lqrobj.R.data * np.sum((useU[:-1, Rcol] - lqrobj.ubase[Rcol]) ** 2, axis=0))
            if useP > 0:
                y += np.sum(lqrobj.P.data * np.sum((useP_[:-1, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0))
            return y

        def __callg__(useX, useU, useP_, y, G, row, col, rec, needg):
            """Calculate the lqr cost with gradient information.

            :param h: float, grid size
            :param useX, useU, useP: ndarray, parsed X, U, P from x, the same with __call__F
            :param y: ndarray, a location to write the objective function onto
            :param G, row, col: the gradient information
            :param rec: if we record row and col
            :param needg: if we need gradient information.

            """
            if rec:
                row[:] = 0
                col_ = np.reshape(col[:num1], (self.N - 1, -1))
            if not needg:
                y[0] = __callf__(useX, useU, useP_)
            else:
                yF = 0.0
                yQ = 0.0
                yR = 0.0
                yP = 0.0
                curG = 0
                h = 1
                G_ = G[:num1].reshape((self.N - 1, -1))
                if useQ > 0:
                    yQ = np.sum(lqrobj.Q.data * np.sum((useX[:-1, Qcol] - lqrobj.xbase[Qcol]) ** 2, axis=0)) * h
                    G_[:, :useQ] = 2 * h * ((useX[:-1, Qcol] - lqrobj.xbase[Qcol]) * lqrobj.Q.data)
                    if rec:
                        col_[:, :useQ] = Qcol + baseCol
                    curG += (self.N - 1) * useQ
                if useR > 0:
                    yR = np.sum(lqrobj.R.data * np.sum((useU[:-1, Rcol] - lqrobj.ubase[Rcol]) ** 2, axis=0)) * h
                    G_[:, useQ: useQ + useR] = 2 * h * ((useU[:-1, Rcol] - lqrobj.ubase[Rcol]) * lqrobj.R.data)
                    if rec:
                        col_[:, useQ: useQ + useR] = Rcol + baseCol + self.dimx
                    curG += (self.N - 1) * useR
                if useP > 0:
                    yP = np.sum(lqrobj.P.data * np.sum((useP_[:-1, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0)) * h
                    G_[:, useQ + useR: useQ + useR + useP] = 2 * h * (useP_[:-1, Pcol] - lqrobj.pbase[Pcol]) * lqrobj.P.data
                    if rec:
                        col_[:, useQ + useR: useQ + useR + useP] = Pcol + baseCol + self.dimx + self.dimu
                    curG += (self.N - 1) * useP
                if useF > 0:
                    yF = np.sum(lqrobj.F.data * ((useX[-1, Fcol] - lqrobj.xfbase[Fcol]) ** 2))
                    G[curG: curG + useF] = 2 * lqrobj.F.data * (useX[-1, Fcol] - lqrobj.xfbase[Fcol])
                    if rec:
                        col[curG: curG + useF] = self.dimpoint * (self.N - 1) + Fcol
                    curG += useF
                y[0] = yF + yQ + yR + yP

        if self.gradmode:
            self.lqrObj = __callg__
        else:
            self.lqrObj = __callf__

    def add_obj(self, obj, path=False):
        """A high level function that add objective function of any kind.

        :param obj: an objective instance or iterable of them, can be LinearObj/LinearPointObj/NonLinearObj/NonLinearPointObj, LqrObj
        :param path: bool, indicating if the constraint is path constraint.
        """
        if isinstance(obj, Iterable):
            for obj_ in obj:
                self.add_obj(obj_, path)
            return
        if isinstance(obj, LinearObj):
            self.addLinearObj(obj)
        elif isinstance(obj, LinearMultiPointObj):
            for idx, A in obj.idx_As:
                self.addLinearPointObj(LinearPointObj(idx, A, *obj.nxup), False)
        elif isinstance(obj, LinearPointObj):
            if path:
                self.linPathObj.append(obj)
            else:
                self.linPointObj.append(obj)
        elif isinstance(obj, (NonLinearObj, NonLinearMultiPointObj)):
            self.addNonLinearObj(obj)
        elif isinstance(obj, NonLinearPointObj):
            self.addNonLinearPointObj(obj, path)
        elif isinstance(obj, LqrObj):
            self.addLQRObj(obj)
        else:
            print("Skipping inappropriate type %s used as objective" % type(obj))

    def add_constr(self, constr, path=False):
        """Add a constraint to the problem.

        :param constr: a constraint object or iterable of them, can be LinearConstr/LinearPointConstr/NonLinearConstr/NonLinearPointConstr
        """
        if isinstance(constr, Iterable):
            for constr_ in constr:
                self.add_constr(constr_, path)
            return
        if isinstance(constr, LinearConstr):
            self.linearConstr.append(constr)
        elif isinstance(constr, LinearPointConstr):
            if path:
                self.linPathConstr.append(constr)
            else:
                self.linPointConstr.append(constr)
        elif isinstance(constr, LinearMultiPointConstr):
            self.linMultiPointConstr.append(constr)
        elif isinstance(constr, NonLinearConstr):
            self.nonLinConstr.append(constr)
        elif isinstance(constr, NonLinearPointConstr):
            if path:
                self.pathConstr.append(constr)
            else:
                self.pointConstr.append(constr)
        elif isinstance(constr, NonLinearMultiPointConstr):
            self.multiPointConstr.append(constr)
        else:
            print("Skipping inappropriate type %s used as constraint" % type(constr))

    def set_N(self, N):
        """Set N.

        :param N: the size of discretization.

        """
        self.N = N
