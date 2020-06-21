#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
trajOptProblem.py

Class for describing the trajectory optimization problems.
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
from scipy.sparse import spmatrix, coo_matrix, csr_matrix


class TrajOptProblem(OptProblem):
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
    def __init__(self, sys, N, t0, tf, gradmode=True, addx=None):
        """Initialize problem by system, discretization grid size, and allowable time

        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :param gradmode: bool, sets if we use gradient mode.
        :param addX: list of addX / one addX / None, additional optimization variables.

        """
        assert isinstance(sys, System)
        self.sys = sys
        self.N = N
        self.tf = tf
        self.t0 = t0
        if np.isscalar(tf):
            self.fixtf = True
        else:
            self.fixtf = False
            assert tf[0] <= tf[1]
        if np.isscalar(t0):
            self.fixt0 = True
        else:
            self.fixt0 = False
            assert t0[0] <= t0[1]
        self.gradmode = gradmode  # this controls if __callg__ will be called
        self.dimx = sys.nx
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
        numU = self.N * self.dimu
        numP = self.N * self.dimp
        if addx is None:
            self.lenAddX = 0
        else:
            if not isinstance(addx, list):
                addx = [addx]
            for tmp in addx:
                assert isinstance(tmp, addX)
            self.addX = addx
            self.lenAddX = sum([tmp.n for tmp in addx])
        numT = 2
        if self.fixt0:
            numT -= 1
        if self.fixtf:
            numT -= 1
        numSol = numX + numU + numP + numT
        self.numX = numX
        self.numU = numU
        self.numP = numP
        self.numT = numT
        self.numTraj = numSol
        self.numSol = numSol + self.lenAddX
        self.t0ind, self.tfind = self.__getTimeIndices()

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
        self._spA = coo_matrix((self.Aval, (self.Arow, self.Acol)), shape=(self.numF, self.numSol)).tocsr()
        self.__setXbound()
        self.__setFbound()
        # detect gradient information
        if self.gradmode:  # in this case, we randomly generate a guess and use it to initialize everything
            randX = self.randomGenX()
            self.__turnOnGrad(randX)
        else:
            raise Exception("Currently trajoptlib only supports gradient mode, it's more robust and efficient.\
                If you have trouble providing gradients, consider using finite difference or automatic differentiation.\
                These are naturally supported by trajoptlib and you can easily turn them on")

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

    def genGuessFromTraj(self, X=None, U=None, P=None, t0=None, tf=None, addx=None, tstamp=None, obj=None, interp_kind='linear'):
        """Alias for gen_guess_from_traj"""
        return self.gen_guess_from_traj(X, U, P, t0, tf, addx, tstamp, obj, interp_kind)

    def gen_guess_from_traj(self, X=None, U=None, P=None, t0=None, tf=None, addx=None, tstamp=None, obj=None, interp_kind='linear'):
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
        # generate t0 and tf, if applicable
        if self.t0ind > 0:
            if t0 is None:
                randX[self.t0ind] = randomGenInBound(self.t0)
            else:
                randX[self.t0ind] = randomGenInBound(t0)
            uset0 = randX[self.t0ind]
        else:
            uset0 = self.t0
        if self.tfind > 0:
            if tf is None:
                randX[self.tfind] = randomGenInBound(self.tf)
            else:
                randX[self.tfind] = randomGenInBound(tf)
            usetf = randX[self.tfind]
        else:
            usetf = self.tf
        teval = np.linspace(uset0, usetf, self.N)

        # interpolation for state variables
        nPoint = self.N
        dimx = self.dimx
        if X is not None:
            Xcol = X.shape[1]
            if not (Xcol == dimx):
                print('The column of X is not %d or %d, not use it' % (dimx))
                X = None
            else:  # use interpolation to do it
                interp(tstamp, X, teval, Xtarget, interp_kind)
        if X is None:
            # straight path go there
            for i in range(nPoint):
                Xtarget[i] = randomGenInBound(self.xbd, self.dimx)
            # randomly generate x0 and xf
            Xtarget[0, :dimx] = randomGenInBound(self.x0bd, self.dimx)
            Xtarget[-1, :dimx] = randomGenInBound(self.xfbd, self.dimx)

        # interpolation for control variable
        if U is not None:
            interp(tstamp, U, teval, Utarget, interp_kind)
        else:
            for i in range(nPoint):
                Utarget[i] = randomGenInBound(self.ubd, self.dimu)
        if self.numP > 0:
            if P is not None:
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
        dimx, dimu, dimp, dimpoint = self.dimx, self.dimu, self.dimp, self.dimpoint
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
        # compute first row mask
        mask = self.Arow == 0
        self.Aval_row0 = self.Aval[mask]
        self.Acol_row0 = self.Acol[mask]
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
        if self.t0ind > 0:
            randX[self.t0ind] = randomGenInBound(self.t0)
        if self.tfind > 0:
            randX[self.tfind] = randomGenInBound(self.tf)
        if self.lenAddX > 0:
            for field, addx in zip(self.__parseAddX__(randX), self.addX):
                field[:] = randomGenInBound([addx.lb, addx.ub], addx.n)
        # I do not have to worry about objaddn since they are linear
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
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        # use one time/state/ctrl/para/h set to detect the sparsity pattern
        ydot, Jac = self.sys.jac_dyn(useT[0], useX[0], useU[0], useP[0])
        nrow, ncol = Jac.shape
        if isinstance(Jac, np.ndarray):
            nnz = nrow * ncol
            self.dynSparse = False
        elif isinstance(Jac, spmatrix):
            nnz = Jac.nnz  # TODO: this is dangerous since nnz might be misleading
            self.dynSparse = True
        if not self.dynSparse:
            dynG = (self.N - 1) * self.dimx * (self.dimpoint + 1 + self.numT)  # G from dyn
            return dynG
        else:  # if sparse, then treat dis and Euler separately
            if self.sys.ode == 'Dis':
                timennz = Jac.getcol(0).nnz + Jac.getcol(-1).nnz
                dynnnz = nnz - timennz
                dynG = (self.N - 1) * (dynnnz + self.dimx + self.numT)
            else:
                timennz = Jac.getcol(0).nnz  # TODO: this is in fact wrong
                eyemat = sparse.diags(np.ones(self.dimx), offsets=1, shape=Jac.shape)
                Jac.data[Jac.data == 0] = 1e-10
                Jac = Jac.tocoo()
                if self.sys.ode == 'Euler':
                    Jac.data *= h
                    # add a block by identity
                    Jac = Jac + eyemat
                elif self.sys.ode == 'BackEuler':
                    Jac.data *= (-h)
                    Jac = Jac + eyemat
                dynnnz = Jac.nnz - timennz
                dynG = (self.N - 1) * (dynnnz + self.dimx + self.numT * self.dimx)
            return dynG

    def __getObjSparsity(self, x):
        """Set sparsity structure of the problem from objective function.

        The objective function pattern is composed of two parts:
        - linear parts. We sum all the coefficients and find sparsity pattern out of it
        - nonlinear parts. Each nonlinear objective is augmented with another row in jacobian
        and a another auxiliary optimization variable s.t. c(x) = y and J += y

        :param x: ndarray, the guess/sol
        :returns: nG: int, # Jacobian from nonlinear objective function

        """
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
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

    def __getTimeIndices(self):
        """Utility function for assigning sparsity structure."""
        t0ind = -1
        tfind = -1
        lenX = self.N * self.dimpoint
        if self.fixt0:
            if self.fixtf:
                pass
            else:
                tfind = lenX
        else:
            if self.fixtf:
                t0ind = lenX
            else:
                t0ind = lenX
                tfind = lenX + 1
        return t0ind, tfind

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

        # set bound on time
        if not self.fixt0:
            xlb[self.t0ind] = self.t0[0]
            xub[self.t0ind] = self.t0[1]
        if not self.fixtf:
            xlb[self.tfind] = self.tf[0]
            xub[self.tfind] = self.tf[1]

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

    def __getTimeGrid(self, x):
        """Based on initial guess x, get the time grid for discretization.

        :param x: ndarray, the guess/sol.
        :returns: h: float, grid size
        :returns: useT: the grid being used

        """
        if self.fixt0:
            uset0 = self.t0
        else:
            uset0 = x[self.t0ind]
        if self.fixtf:
            usetf = self.tf
        else:
            usetf = x[self.tfind]
        h = (1.0 * (usetf - uset0)) / (self.N - 1)
        if h <= 0 and not (self.fixt0 and self.fixtf):
            h = 1e-6
        useT = np.linspace(uset0, usetf, self.N)
        return h, useT

    def __parseX__(self, x):
        """Parse guess/sol into X, U, P"""
        n_var_pnt = self.dimx + self.dimu + self.dimp
        X = np.reshape(x[:self.N * n_var_pnt], (self.N, n_var_pnt))
        useX = X[:, :self.dimx]
        useU = X[:, self.dimx:self.dimpoint - self.dimp]
        useP = X[:, self.dimpoint - self.dimp:]
        return useX, useU, useP

    def parseF(self, guess):
        """Alias for :func:`~trajOptLib.TrajOptProblem.parse_f`"""
        return self.parse_f(guess)

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
        # add linear stuff
        for i, j, v in zip(self.Arow, self.Acol, self.Aval):
            y[i] += v * guess[j]
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
        if self.t0ind > 0:
            t0bound = checkInBounds(guess[self.t0ind], self.t0)
        else:
            t0bound = None
        if self.tfind > 0:
            tfbound = checkInBounds(guess[self.tfind], self.tf)
        else:
            tfbound = None
        if self.lenAddX > 0:
            addx = self.__parseAddX__(guess)
            addXbound = [checkInBounds(addx_, [addx__.lb, addx__.ub]) for addx_, addx__ in zip(addx, self.addX)]
        else:
            addXbound = None
        result = {'obj': obj, 'dyn': dynCon, 'Xbd': Xbound, 'Ubd': ubound, 'x0bd': x0bound, 'xfbd': xfbound, 'Pbd': pbound, 't0bd': t0bound, 'tfbd': tfbound, 'addXbd': addXbound}
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
        h, tgrid = self.__get_time_grid__(sol)
        obj = self.__parseObj__(sol)
        if self.lenAddX == 0:
            return {'t': tgrid, 'x': X, 'u': U, 'p': P, 'obj': obj}
        else:
            return {'t': tgrid, 'x': X, 'u': U, 'p': P, 'addx': self.__parseAddX__(sol), 'obj': obj}

    def parseSol(self, sol):
        """Alias for :func:`~trajOptLib.TrajOptProblem.parse_sol`"""
        return self.parse_sol(sol)

    def __get_time_grid__(self, x):
        """Based on initial guess x, get the time grid for discretization.

        :param x: ndarray, the guess/sol.
        :returns: h: float, grid size
        :returns: useT: the grid being used

        """
        if self.fixt0:
            uset0 = self.t0
        else:
            uset0 = x[self.t0ind]
        if self.fixtf:
            usetf = self.tf
        else:
            usetf = x[self.tfind]
        h = (usetf - uset0) / (self.N - 1)
        useT = np.linspace(uset0, usetf, self.N)
        return h, useT

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
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        # evaluate objective function
        self.__objModeF__(0, h, useT, useX, useU, useP, x, y)
        # evaluate the system dynamics constraint
        curRow = 1
        curRow = self.__dynconstrModeF__(curRow, h, useT, useX, useU, useP, y)
        # evaluate other constraints
        curRow = self.__constrModeF__(curRow, h, useT, useX, useU, useP, x, y)
        return curRow

    def __objModeF__(self, curRow, h, useT, useX, useU, useP, x, y):
        """Calculate objective function. F mode

        :param curRow: int, index from which we write on
        :param h, useT, useX, useU, useP: parsed solution
        :param x: ndarray, the sol, it is used for linear constraints
        :param y: ndarray, the F to be written. The first row stores the objective function

        """
        y[0] = 0.0
        tmpout = np.zeros(1)
        for obj in self.linPointObj:
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            y[0] += obj.A.dot(tmpx)
        for obj in self.linPathObj:
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                obj.__callf__(tmpx, tmpout)
                y[0] += tmpout[0] * h
        for obj in self.linearObj:
            y[0] += obj.A.dot(x)
        for obj in self.nonPointObj:
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            obj.__callf__(tmpx, tmpout)
            y[0] += tmpout[0]
        for obj in self.nonPathObj:
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                obj.__callf__(tmpx, tmpout)
                y[0] += tmpout[0] * h
        for obj in self.nonLinObj:
            if isinstance(obj, NonLinearObj):
                obj.__callf__(x, tmpout)
            else:  # NonLinearMultiPointObj
                xins = [np.concatenate(([useT[idx]], useX[idx], useU[idx], useP[idx])) for idx in obj.indexes]
                obj.__callf__(xins, tmpout)
            y[0] += tmpout[0]
        # add lqr cost, if applicable
        if self.lqrObj is not None:
            y[0] += self.lqrObjf(h, useX, useU, useP)

    def __constrModeF__(self, curRow, h, useT, useX, useU, useP, x, y):
        """Calculate constraint function. F mode

        :param curRow: int, index from which we write on
        :param h, useT, useX, useU, useP: parsed solution
        :param y: ndarray, the F to be written
        :returns: curRow: current row after we write on y

        """
        for constr in self.pointConstr:
            tmpx = np.concatenate(([useT[constr.index]], useX[constr.index], useU[constr.index], useP[constr.index]))
            constr.__evalf__(tmpx, y[curRow: curRow + constr.nf])
            curRow += constr.nf
        for constr in self.multiPointConstr:
            xin = [np.concatenate(([useT[idx]], useX[idx], useU[idx], useP[idx])) for idx in constr.indexes]
            constr.__evalf__(xin, y[curRow: curRow + constr.nf])
            curRow += constr.nf
        for constr in self.pathConstr:
            for i in range(self.N):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                constr.__evalf__(tmpx, y[curRow: curRow + constr.nf])
            self.numF
            curRow += constr.nf
        for constr in self.nonLinConstr:
            constr.__evalf__(x, y[curRow: curRow + constr.nf])
            curRow += constr.nf
        return curRow

    def __dynconstrModeF__(self, curRow, h, useT, useX, useU, useP, y):
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
            ydot = self.sys.dyn(useT[i], useX[i], useU[i], useP[i])
            cDyn[i] = useX[i] + h * ydot - useX[i + 1]
        return curRow

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix.

        :param x: ndarray, the solution to the problem
        :param y: ndarray, return F
        :param G, row, col: ndarray, information of gradient
        :param rec, needg: if we record/ if we need gradient

        """
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        # loop over all system dynamics constraint
        curRow = 1
        curNg = 0
        curRow, curNg = self.__dynconstrModeG(curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
        curRow, curNg = self.__constrModeG(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        curRow, curNg = self.__objModeG(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        return curRow, curNg

    def __dynconstrModeG(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        cDyn = np.reshape(y[curRow:curRow + (self.N - 1) * self.dimx], (self.N - 1, self.dimx))
        for i in range(self.N - 1):
            # evaluate gradient of system dynamics TODO: support other types of integration scheme
            ydot, Jac = self.sys.jac_dyn(useT[i], useX[i], useU[i], useP[i])  # TODO: support in-place Jacobian
            cDyn[i] = useX[i] + h * ydot - useX[i + 1]
            if needg:
                if not self.dynSparse:
                    Jac *= h  # always useful
                    baseCol = i * self.dimpoint
                    # assign a block for x
                    G[curNg: curNg + dimx * dimx] = (Jac[:, 1:1 + dimx] + np.eye(dimx)).flatten()
                    if rec:
                        tmpMat = np.tile(np.arange(dimx), (dimx, 1))
                        row[curNg: curNg + dimx * dimx] = curRow + tmpMat.T.flatten()
                        col[curNg: curNg + dimx * dimx] = baseCol + tmpMat.flatten()
                    curNg += dimx * dimx
                    # assign a block for u
                    G[curNg: curNg + dimx * dimu] = Jac[:, 1 + dimx:1 + dimx + dimu].flatten()
                    if rec:
                        row[curNg: curNg + dimx * dimu] = curRow + np.tile(np.arange(dimx), (dimu, 1)).T.flatten()
                        col[curNg: curNg + dimx * dimu] = baseCol + dimx + np.tile(np.arange(dimu), dimx).flatten()
                    curNg += dimx * dimu
                    # assign a block for p, if necessary
                    if dimp > 0:
                        G[curNg: curNg + dimx * dimp] = Jac[:, 1 + dimx + dimu:1 + dimx + dimu + dimp].flatten()
                        if rec:
                            row[curNg: curNg + dimx * dimp] = curRow + np.tile(np.arange(dimx), (dimp, 1)).T.flatten()
                            col[curNg: curNg + dimx * dimp] = baseCol + dimx + dimu + np.tile(np.arange(dimp), dimx).flatten()
                        curNg += dimx * dimp
                    # assign the diagonal block for x_{k+1}
                    G[curNg: curNg + dimx] = -1.0
                    if rec:
                        row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                        col[curNg: curNg + dimx] = baseCol + self.dimpoint + np.arange(dimx)
                    curNg += dimx
                    # assign a column for t0, if necessary
                    if self.t0ind > 0:
                        G[curNg: curNg + dimx] = -ydot / (self.N - 1) + Jac[:, 0] * (1 - float(i) / (self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.t0ind
                        curNg += dimx
                    # assign a column for tf, if necessary
                    if self.tfind > 0:
                        G[curNg: curNg + dimx] = ydot / (self.N - 1) + Jac[:, 0] * (float(i) / (self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.tfind
                        curNg += dimx
                else:
                    Jac.data[Jac.data == 0] = 1e-10
                    Jac.data *= h  # as always, no damage to this
                    eyemat = sparse.diags(np.ones(self.dimx), offsets=1, shape=Jac.shape)
                    Jac += eyemat  # TODO: this assume forward Euler, backward needs update
                    dynnnz = Jac.nnz - Jac.getcol(0).nnz  # TODO: this might cause serious issues since it changes structure
                    # convention is to assign all but first column, then assign first column
                    if not sparse.isspmatrix_coo(Jac):
                        Jac = Jac.tocoo()
                    rows = Jac.row
                    cols = Jac.col
                    data = Jac.data
                    timeone = cols == 0
                    nontime = np.logical_not(timeone)
                    G[curNg: curNg + dynnnz] = data[nontime]
                    if rec:
                        row[curNg: curNg + dynnnz] = curRow + rows[nontime]
                        col[curNg: curNg + dynnnz] = self.__patchCol__(i, cols[nontime])
                    curNg += dynnnz
                    # for x_{k+1}
                    G[curNg: curNg + dimx] = -1.0
                    if rec:
                        row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                        col[curNg: curNg + dimx] = np.arange(dimx) + (i + 1) * self.dimpoint
                    curNg += dimx
                    if self.t0ind > 0:
                        G[curNg: curNg + dimx] = -ydot / (self.N - 1) + Jac.getcol(0).toarray()[:, 0] * (1 - float(i) / (self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.t0ind
                        curNg += dimx
                    if self.tfind > 0:
                        G[curNg: curNg + dimx] = ydot / (self.N - 1) + Jac.getcol(0).toarray()[:, 0] * (float(i) / (self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.tfind
                        curNg += dimx
            curRow += dimx
        return curRow, curNg

    def __objModeG(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg, first_row=False):
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
            self.lqrObj(h, useX, useU, useP, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
            y[curRow] = tmpout[0]
            if rec:
                rowpiece[:] = curRow
            curNg += self.LQRnG
            curRow += 1

        # still in the point, path, obj order
        if self.nonPointObj:
            for obj in self.nonPointObj:
                tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
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
                    tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                    Gpiece = G[curNg: curNg + obj.nG]
                    rowpiece = row[curNg: curNg + obj.nG]
                    colpiece = col[curNg: curNg + obj.nG]
                    obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                    y[curRow] += tmpout[0] * h
                    Gpiece[:] *= h
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
                    xins = [np.concatenate(([useT[idx]], useX[idx], useU[idx], useP[idx])) for idx in obj.indexes]
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

    def __constrModeG(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
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
                tmpx = np.concatenate(([useT[constr.index]], useX[constr.index], useU[constr.index], useP[constr.index]))
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
                xin = [np.concatenate(([useT[idx]], useX[idx], useU[idx], useP[idx])) for idx in constr.indexes]
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
                    tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
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
        return np.dot(x[self.Acol_row0], self.Aval_row0)

    def __gradient__(self, x, g):
        """Evaluation of the gradient of objective function.

        :param x: guess/solution to the problem
        :param g: the gradient of objective function w.r.t x to be written into

        """
        g[:] = 0
        g[self.Acol_row0] = self.Aval_row0
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
        f += self._spA.dot(x)
        return 0

    def __jacobian__(self, x, g, row, col, rec):
        """Evaluate the Jacobian of the problem.

        :param x: guess/solution to the problem
        :param g: the vector being written on for Jacobian entries
        """
        y = np.zeros(self.numF)
        self.__callg__(x, y, g, row, col, rec, True)
        g[self.nG:] = self.Aval
        # append the linear parts here
        if rec:
            row[self.nG:] = self.Arow
            col[self.nG:] = self.Acol
        return 1

    def __patchCol__(self, index, col, col_offset=0):
        """Find which indices it belongs to the original one for a local matrix at col"""
        col = col[col > 0]
        if index < 0:
            index = index + self.N
        return col - 1 + col_offset + index * self.dimpoint

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
        self.LQRnG = (self.N - 1) * (useQ + useR + useP) + useF + self.numT
        num1 = (self.N - 1) * (useQ + useR + useP)
        baseCol = self.dimpoint * np.arange(self.N - 1)[:, np.newaxis]  # a nPoint by 1 column matrix

        def __callf__(h, useX, useU, useP_):
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
                y += np.sum(lqrobj.Q.data * np.sum((useX[:-1, Qcol] - lqrobj.xbase[Qcol]) ** 2, axis=0)) * h
            if useR > 0:
                y += np.sum(lqrobj.R.data * np.sum((useU[:-1, Rcol] - lqrobj.ubase[Rcol]) ** 2, axis=0)) * h
            if useP > 0:
                y += np.sum(lqrobj.P.data * np.sum((useP_[:-1, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0)) * h
            y += tfweight * (h * (self.N - 1))
            return y

        def __callg__(h, useX, useU, useP_, y, G, row, col, rec, needg):
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
                y[0] = __callf__(h, useX, useU, useP_)
            else:
                yF = 0.0
                yQ = 0.0
                yR = 0.0
                yP = 0.0
                yTf = tfweight * (h * (self.N - 1))
                curG = 0
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
                if self.t0ind > 0:
                    G[curG] = -(yQ + yR + yP) / h / (self.N - 1) - tfweight
                    if rec:
                        col[curG: curG + 1] = self.t0ind
                    curG += 1
                if self.tfind > 0:
                    G[curG] = (yQ + yR + yP) / h / (self.N - 1) + tfweight
                    if rec:
                        col[curG: curG + 1] = self.tfind
                    curG += 1
                y[0] = yF + yQ + yR + yP + yTf

        if self.gradmode:
            self.lqrObj = __callg__
            self.lqrObjf = __callf__
        else:
            self.lqrObj = __callf__

    def addLinearObj(self, linObj):
        """Add linear objective function. Obsolete, use add_obj instead.

        :param linObj: linearObj class

        """
        assert isinstance(linObj, LinearObj)
        self.linearObj.append(linObj)

    def addLinearPointObj(self, linPointObj, path=False):
        """Add linear point objective function. Obsolete, use add_obj instead.

        :param linPointObj: linearPointObj class
        :param path: bool, if this is path obj (at every point except for final one)

        """
        assert isinstance(linPointObj, LinearPointObj)
        if path:
            self.linPathObj.append(linPointObj)
        else:
            self.linPointObj.append(linPointObj)

    def addNonLinearObj(self, nonlinObj):
        """Add nonlinear objective function. Obsolete, use add_obj instead.

        :param nonLinObj: a nonLinObj class

        """
        assert isinstance(nonlinObj, (NonLinearObj, NonLinearPointObj))
        self.nonLinObj.append(nonlinObj)

    def addNonLinearPointObj(self, nonPntObj, path=False):
        """Add nonlinear point objective. Obsolete, use add_obj instead.

        :param nonPntObj: nonLinObj class
        :param path: bool, if this obj is pointwise

        """
        assert isinstance(nonPntObj, NonLinearPointObj)
        if path:
            self.nonPathObj.append(nonPntObj)
        else:
            self.nonPointObj.append(nonPntObj)

    def addNonLinearPointConstr(self, pntConstr, path=False):
        """Add point constraint. Obsolete, use add_constr instead.

        :param pntConstr: pointConstr class
        :param path: bool, if this obj

        """
        assert isinstance(pntConstr, NonLinearPointConstr)
        if path:
            self.pathConstr.append(pntConstr)
        else:
            self.pointConstr.append(pntConstr)

    def addNonLinearConstr(self, constr):
        """Add a general nonlinear constraint. Obsolete, use add_constr instead.

        :param constr: nonLinConstr class

        """
        assert isinstance(constr, NonLinearConstr)
        self.nonLinConstr.append(constr)

    def addLinearConstr(self, constr):
        assert isinstance(constr, LinearConstr)
        self.linearConstr.append(constr)

    def addLinearPointConstr(self, constr, path=False):
        assert isinstance(constr, LinearPointConstr)
        if path:
            self.linPathConstr.append(constr)
        else:
            self.linPointConstr.append(constr)

    def addObj(self, obj, path=False):
        """Alias for :func:`~trajOptLib.TrajOptProblem.add_obj`"""
        self.add_obj(obj, path)

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
            self.addLinearPointObj(obj, path)
        elif isinstance(obj, (NonLinearObj, NonLinearMultiPointObj)):
            self.addNonLinearObj(obj)
        elif isinstance(obj, NonLinearPointObj):
            self.addNonLinearPointObj(obj, path)
        elif isinstance(obj, LqrObj):
            self.addLQRObj(obj)
        else:
            print("Skipping inappropriate type %s used as objective" % type(obj))

    def addConstr(self, constr, path=False):
        """Alias for :func:`~trajOptLib.TrajOptProblem.add_constr`"""
        self.add_constr(constr, path)

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

    def set_t0_tf(self, t0, tf):
        """Set t0 and tf.

        :param t0: float/ndarray (2,) allowable t0
        :param tf: float/ndarray (2,) allowable tf

        """
        self.t0 = t0
        self.tf = tf
        if np.isscalar(tf):
            self.fixtf = True
        else:
            self.fixtf = False
            assert tf[0] <= tf[1]
        if np.isscalar(t0):
            self.fixt0 = True
        else:
            self.fixt0 = False
            assert t0[0] <= t0[1]
