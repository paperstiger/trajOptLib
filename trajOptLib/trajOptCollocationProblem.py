#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
trajOptCollocationProblem.py

This class implements the direct collocation approach for humanoid trajectory optimization
"""
import numpy as np
import logging
from .trajOptBase import system, linearObj, linearPointObj, nonLinObj, nonPointObj, pointConstr, nonLinConstr, lqrObj
from .libsnopt import snoptConfig, probFun, solver
from .utility import randomGenInBound, checkInBounds
from .trajOptProblem import trajOptProblem, addX
from scipy import sparse
from scipy.sparse import spmatrix, coo_matrix


class daeSystem(object):
    """A DAE system."""
    def __init__(self, nx, nu, np, nf, nG):
        """Constructor for the problem.

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


class trajOptCollocProblem(probFun):
    """A class for definition of trajectory optimization problem using collocation constraints.

    A general framework for using this class is to:

    1. Define a class which implements a DAE system, f(t, x, p, u)=0. This is a typical case for humanoid problem.
    2. Optionally, write desired cost function by inheriting/creating from the list of available cost functions.
    3. Optionally, write desired constraint functions by inheriting from available constraints.
    4. Create this class with selected system, discretization, t0, tf range, gradient option
    5. Set bounds for state, control, parameters, x0 and xf
    6. Add objective functions and constraints to this class
    7. Call preProcess method explicitly
    8. Create snoptConfig instance and choose desired options
    9. Construct the solver
    10. Use the solver to solve with either automatic guess or user provided guess

    The system dynamics constraints are imposed using direct collocation approach with some additional optimization
    variables as suggested by others.

    """
    def __init__(self, sys, N=20, t0=0.0, tf=1.0, addx=None):
        """Initialize problem by system, discretization grid size, and allowable time

        Change history: now I remove gradmode option since I require the gradient be provided analytically all the time.
        I remove unnecessary linear objective functions.
        I reorder variable so q, dq, ddq, u, p are at consecutive place.
        TODO: add linear constraints and enable it in SNOPT.

        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :param addX: list of addX / one addX / None, additional optimization variables.

        """
        assert isinstance(sys, daeSystem)
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
        self.dimx = sys.nx
        self.dimdyn = sys.nf  # dimension of dynamics constraint
        self.dimu = sys.nu
        self.dimp = sys.np
        self.dimpoint = sys.nx + sys.nu + sys.np  # each point have those variables
        self.ubd = [None, None]
        self.xbd = [None, None]
        self.pbd = [None, None]
        self.x0bd = [None, None]
        self.xfbd = [None, None]
        # lqr cost function
        self.lqrObj = None
        # nonlinear cost function
        self.nonLinObj = []  # stores general nonlinear cost
        self.nonPointObj = []  # stores nonlinear cost imposed at a point
        self.nonPathObj = []  # stores Lagrange integral cost. Includes LQR cost
        # nonlinear constraints. Linear constraints are treated as nonlinear
        self.pointConstr = []  # general constraint imposed at a certain point, such as initial and final point
        self.pathConstr = []  # general constraint imposed everywhere such as collision avoidance
        self.nonLinConstr = []  # stores general nonlinear constraint  TODO: implement linear constraints. They are cheap.
        # calculate number of variables to be optimized, time are always the last
        numX = (2 * self.N - 1) * self.dimx
        numU = (2 * self.N - 1) * self.dimu
        numP = (2 * self.N - 1) * self.dimp
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
        self.numDynVar = numX + numU + numP  # numDynVar includes q, dq, ddq, u, p
        self.nPoint = 2 * self.N - 1
        self.numTraj = numSol  # make it clear, numTraj contains numDynVar + time
        self.numSol = numSol + self.lenAddX
        self.t0ind, self.tfind = self.__getTimeIndices()

    def preProcess(self):
        """Initialize the instances of probFun now we are ready.

        Call this function after the objectives and constraints have been set appropriately.
        It calculate the space required for SNOPT and allocates sparsity structure if necessary.

        """
        numDyn = self.dimdyn * (2 * self.N - 1)  # constraints from system dynamics, they are imposed everywhere
        dynDefectSize = 4 * self.dimdyn  #TODO: this number 4 might change
        defectSize = dynDefectSize + self.dimu + self.dimp  # from x, dx, and u, p are average
        self.defectSize = defectSize
        defectDyn = (self.N - 1) * defectSize  # from enforcing those guys
        self.numDyn = numDyn + defectDyn
        numC = 0
        for constr in self.pointConstr:
            numC += constr.nf
        for constr in self.pathConstr:
            numC += self.N * constr.nf  # TODO: verify we do not have to impose those constraints on collocation points
        for constr in self.nonLinConstr:
            numC += constr.nf
        self.numF = 1 + numDyn + defectDyn + numC
        probFun.__init__(self, self.numSol, self.numF)  # not providing G means we use finite-difference
        # we are ready to write Aval, Arow, Acol for this problem. They are arranged right after dynamics
        curRow, curNA = self.__setAPattern(numDyn + 1)
        self.__setXbound()
        self.__setFbound()
        # detect gradient information
        randX = self.randomGenX()
        self.__turnOnGrad(randX)

    def __setAPattern(self, curRow):
        """Set sparsity pattern for A. curRow is current row."""
        # so rows are from numDyn + 1 to numDyn + 1 + defectDyn
        # the size of A is (4 * 5 * self.dimdyn + 3*(self.dimu + self.dimp)) * (self.N - 1)
        # sparse matrices for defect of states
        dimx, dimu, dimp, dimpoint, dimdyn = self.dimx, self.dimu, self.dimp, self.dimpoint, self.dimdyn
        dynDefectSize = 4 * self.dimdyn
        self.matL = np.zeros((dynDefectSize, self.dimx))
        self.matM = np.zeros((dynDefectSize, self.dimx))
        self.matR = np.zeros((dynDefectSize, self.dimx))
        two = 2*self.dimdyn
        one = self.dimdyn
        three = self.dimx
        h = (self.tf - self.t0) / (self.N - 1)
        np.fill_diagonal(self.matL[:two, :two], 0.5)
        np.fill_diagonal(self.matL[:two, one:], h/8)
        np.fill_diagonal(self.matL[two:, :two], -1.5/h)
        np.fill_diagonal(self.matL[two:, one:], -0.25)
        np.fill_diagonal(self.matM[:two, :two], -1)
        np.fill_diagonal(self.matM[two:, one:], -1)
        np.fill_diagonal(self.matR[:two, :two], 0.5)
        np.fill_diagonal(self.matR[:two, one:], -h/8)
        np.fill_diagonal(self.matR[two:, :two], 1.5/h)
        np.fill_diagonal(self.matR[two:, one:], -0.25)
        self.spL = coo_matrix(self.matL)
        self.spM = coo_matrix(self.matM)
        self.spR = coo_matrix(self.matR)
        lenA = (4*5*self.dimdyn + 3*(self.dimu + self.dimp)) * (self.N - 1)
        A = np.zeros(lenA)
        row = np.zeros(lenA, dtype=int)
        col = np.zeros(lenA, dtype=int)
        curNA = 0
        for i in range(self.N - 1):
            midi = 2*i + 1
            lefti = 2*i
            righti = 2*(i + 1)
            A[curNA: curNA + self.spL.nnz] = self.spL.data
            row[curNA: curNA + self.spL.nnz] = self.spL.row + curRow
            col[curNA: curNA + self.spL.nnz] = self.spL.col + lefti * dimpoint
            curNA += self.spL.nnz
            A[curNA: curNA + self.spM.nnz] = self.spM.data
            row[curNA: curNA + self.spM.nnz] = self.spM.row + curRow
            col[curNA: curNA + self.spM.nnz] = self.spM.col + midi * dimpoint
            curNA += self.spM.nnz
            A[curNA: curNA + self.spR.nnz] = self.spR.data
            row[curNA: curNA + self.spR.nnz] = self.spR.row + curRow
            col[curNA: curNA + self.spR.nnz] = self.spR.col + righti * dimpoint
            curNA += self.spR.nnz
            curRow += 4*dimdyn
        # then do the constraint of u and p on nodes and knots
        for i in range(self.N - 1):
            midi = 2*i + 1
            lefti = 2*i
            righti = 2*(i + 1)
            A[curNA: curNA + 2*dimu] = 0.5
            A[curNA + 2*dimu: curNA + 3*dimu] = -1
            if dimp > 0:
                A[curNA + 3*dimu: curNA + 3 * dimu + 2*dimp] = 0.5
                A[curNA + 3*dimu + 2 * dimp: curNA + 3*dimu + 3*dimp] = -1
            row[curNA: curNA + 3 * dimu] = curRow + np.tile(np.arange(dimu), (3, 1)).flatten()
            col[curNA: curNA + dimu] = lefti * dimpoint + dimx + np.arange(dimu)
            col[curNA + dimu: curNA + 2 * dimu] = righti * dimpoint + dimx + np.arange(dimu)
            col[curNA + 2 * dimu: curNA + 3 * dimu] = midi * dimpoint + dimx + np.arange(dimu)
            curNA_ = curNA + 3 * dimu
            if dimp > 0:
                row[curNA_: curNA_ + 3 * dimp] = curRow + dimu + np.tile(np.arange(dimp), (3, 1)).flatten()
                col[curNA_: curNA_ + dimp] = lefti * dimpoint + dimx + dimu + np.arange(dimp)
                col[curNA_ + dimp: curNA_ + 2 * dimp] = righti * dimpoint + dimx + dimu + np.arange(dimp)
                col[curNA_ + 2 * dimp: curNA_ + 3 * dimp] = midi * dimpoint + dimx + dimu + np.arange(dimp)
            curNA += 3 * (dimu + dimp)
            curRow += dimu + dimp
        self.Aval = A
        self.Arow = row
        self.Acol = col
        return curRow, curNA

    def randomGenX(self):
        """A more reansonable approach to generate random guess for the problem.

        It considers bounds on initial and final states so this is satisfied.
        Then it linearly interpolate between states.
        Controls are randomly generated within control bound, if it presents. Otherwise [-1, 1]

        :return x: ndarray, (numSol, ) an initial guess of the solution
        """
        nPoint = self.nPoint
        dimx = self.dimx
        randX = 2*np.random.random(self.numSol) - 1
        X, U, P = self.__parseX__(randX)
        # randomly generate x0 and xf
        X[0, :dimx] = randomGenInBound(self.x0bd, self.dimx)
        X[-1, :dimx] = randomGenInBound(self.xfbd, self.dimx)
        # straight path go there
        for i in range(self.dimx):
            X[:, i] = np.linspace(X[0, i], X[-1, i], nPoint)
        for i in range(nPoint):
            U[i] = randomGenInBound(self.ubd, self.dimu)
        if self.numP > 0:
            for i in range(nPoint):
                P[i] = randomGenInBound(self.pbd, self.dimp)
        if self.t0ind > 0:
            randX[self.t0ind] = randomGenInBound(self.t0)
        if self.tfind > 0:
            randX[self.tfind] = randomGenInBound(self.tf)
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
        arg1, objSparseMode = self.__getObjSparsity(x0)
        numObjG = arg1
        self.numObjG = numObjG
        self.objSparseMode = objSparseMode
        # summarize number of pure linear constraints
        numDynG = self.__getDynSparsity(x0)
        numCG = 0  # G from C
        # I only care about those in numC
        for constr in self.pointConstr:
            numCG += constr.nG
        for constr in self.pathConstr:
            numCG += self.N * constr.nG
        for constr in self.nonLinConstr:
            numCG += constr.nG
        numG = numObjG + numDynG + numCG
        self.numG = numG
        # we can change nG and grad
        self.nG = numG

    def __getDynSparsity(self, x):
        """Set sparsity of the problem caused by system dynamics and other nonlinear constraints.

        Sparsity pattern should be quite straight-forward to conclude since we introduce many auxiliary variables.
        The nG from dae system should be directly used and it is complete.
        The defect constraints introduce 5*dimq*4 gradients

        :param x: ndarray, the guess/sol

        """
        # this row is the case when defectDyn are imposed in G rather than A
        # dynG = self.sys.nG * self.nPoint + (20 * self.dimdyn + 3*(self.dimu + self.dimp)) * (self.N - 1)
        dynG = self.sys.nG * self.nPoint
        return dynG

    def __getObjSparsity(self, x):
        """Set sparsity structure of the problem from objective function.

        For objective function, we do calculation by assigning gradient to a vector of length numSol.
        The sparsity pattern is remembered which is a mask that is then used to assign G.
        This function also determines if it is safe to turn on aggressive mode where we directly assign G.
        We can safely turn it on in the following conditions:
        - only lqrObj is supplied
        - only one of nonLinObj, nonPointObj, nonPathObj, linearObj, linPointObj, linPathObj is non-empty
        In other cases, we still have to remember which maps to which and there is not much improvement compared with current implementation
        # TODO: add support such that user has option to safely insert an objective function. Such as path obj plus terminal cost. In long run, add linear obj support

        :param x: ndarray, the guess/sol
        :returns: indices: ndarray, the indices with non-zero objective gradient

        """
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        # check sparseObj mode
        sumLenObj = len(self.nonLinObj) + len(self.nonPointObj) + len(self.nonPathObj)
        if sumLenObj == 1 and self.lqrObj is None:
            self.sparseObjMode = True
        elif sumLenObj == 0 and self.lqrObj is not None:
            self.sparseObjMode = True
        else:
            self.sparseObjMode = False
        assert self.sparseObjMode  # TODO: support non-sparse obj mode
        # for sparseObjMode, we return #nnz as we will add to G, col will be added in __callg__
        if self.sparseObjMode:
            if sumLenObj == 1:
                nG = 0
                for obj in self.nonLinObj:
                    nG += obj.nG
                for obj in self.nonPointObj:
                    nG += obj.nG
                for obj in self.nonPathObj:
                    nG += (2 * self.N - 1) * obj.nG
                return nG, True
            else:
                return self.LQRnG, True

    def __assignTo__(self, index, nnz, value, target, selfadd):
        """Assign a block of value to a target vector. This helps us finding sparsity/assigning values.

        :param index: int, index of a point constraint or -1 indicates it is overall
        :param nnz: indices of nnz elements. For point constraint only, it is within xdim
        :param value: scalar / ndarray, the values we want to assign. If ndarray, its length has to equal nnz
        :param target: ndarray, the target array to assign value to
        :param selfadd: bool, indicates if we eadd value to target, otherwise we assign

        """
        if np.isscalar(value):
            value = np.ones_like(nnz, dtype=type(value)) * value
        if index == -1:  # simple case, we do not care too much
            if selfadd:
                target[nnz] += value
            else:
                target[nnz] = target[nnz] + value
        else:
            # for time
            mask = nnz < 1
            if np.sum(mask) > 0:
                if self.t0ind > -1:
                    if selfadd:
                        target[self.t0ind] -= value[mask]
                    else:
                        target[self.t0ind] = -value[mask]
                if self.tfind > -1:
                    if selfadd:
                        target[self.tfind] += value[mask]
                    else:
                        target[self.tfind] = value[mask]
            # for others
            nnzcol = self.__patchCol__(index, nnz)
            if selfadd:
                target[nnzcol] += value
            else:
                target[nnzcol] = value

    def __getTimeIndices(self):
        """Utility function for assigning sparsity structure."""
        t0ind = -1
        tfind = -1
        if self.fixt0:
            if self.fixtf:
                pass
            else:
                tfind = self.numTraj
        else:
            if self.fixxf:
                t0ind = self.numTraj
            else:
                t0ind = self.numTraj + 1
                tfind = self.numTraj + 2
        return t0ind, tfind

    def __setXbound(self):
        """Set bounds on decision variables."""
        # create bound on x
        dimpnt = self.dimpoint
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        xlb = np.zeros(self.numSol)
        xub = np.zeros(self.numSol)
        Mxlb = np.reshape(xlb[:self.numDynVar], (self.nPoint, dimpnt))
        Mxub = np.reshape(xub[:self.numDynVar], (self.nPoint, dimpnt))
        Mulb = Mxlb[:, dimx:dimx+dimu]
        Muub = Mxub[:, dimx:dimx+dimu]
        Mplb = Mxlb[:, dimpnt-dimp:dimpnt]
        Mpub = Mxub[:, dimpnt-dimp:dimpnt]
        # set bounds for q and dq, agree with previous convention
        if self.xbd[0] is not None:
            Mxlb[:, :dimx] = self.xbd[0]
            # set lb for x0 and xf
            if self.x0bd[0] is not None:
                Mxlb[0, :dimx] = self.x0bd[0]
            if self.x0bd[0] is not None:
                Mxlb[-1, :dimx] = self.xfbd[0]
        else:
            Mxlb[:, :dimx] = -1e20
        if self.xbd[1] is not None:
            Mxub[:, :dimx] = self.xbd[1]
            # set ub for x0 and xf
            if self.x0bd[1] is not None:
                Mxub[0, :dimx] = self.x0bd[1]
            if self.xfbd[1] is not None:
                Mxub[-1, :dimx] = self.xfbd[1]
        else:
            Mxub[:, :dimx] = 1e20
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
        if not self.fixt0:
            xlb[self.t0ind] = self.t0[0]
            xub[self.t0ind] = self.t0[1]
        # set bound on time
        if not self.fixtf:
            xlb[self.tfind] = self.tf[0]
            xub[self.tfind] = self.tf[1]
        # set bound on addX
        if self.lenAddX != 0:
            curN = self.numTraj
            for addx in self.addX:
                xlb[curN: curN + addx.n] = addx.lb
                xub[curN: curN + addx.n] = addx.ub
        # assign to where it should belong to
        self.xlb = xlb
        self.xub = xub

    def __setFbound(self):
        """Set bound on F"""
        # set bound on F
        numF = self.numF
        numDyn = self.numDyn
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
        h = (usetf - uset0) / (self.N - 1)
        useT = np.linspace(uset0, usetf, self.nPoint)
        return h, useT

    def __parseX__(self, x):
        """Parse guess/sol into X, U, P"""
        X = np.reshape(x[:self.numDynVar], (self.nPoint, self.dimpoint))
        useX = X[:, :self.dimx]
        useU = X[:, self.dimx:self.dimpoint - self.dimp]
        useP = X[:, self.dimpoint - self.dimp:]
        return useX, useU, useP

    def parseF(self, guess):
        """Give an guess, evaluate it and parse into parts.

        :param guess: ndarray, (numSol, ) a guess or a solution to check
        :returns: dict, containing objective and parsed constraints

        """
        assert len(guess) == self.numSol
        N = self.N
        dimx = self.dimx
        dimdyn = self.dimdyn
        y = np.ones(self.numF)
        if self.gradmode:
            self.__callg__(guess, y, np.zeros(1), np.zeros(1), np.zeros(1), False, False)
        else:
            self.__callf__(guess, y)
        obj = y[0]
        dynCon = np.reshape(y[1:(2*N-1)*dimdyn+1], (2*N - 1, dimdyn))
        curN = 1 + (2 * N - 1) * dimdyn
        curNf = curN + (N - 1) * 4 * dimdyn
        defectCon = np.reshape(y[curN: curNf], (N - 1, 4*dimdyn))
        curN = curNf
        pointCon = []
        for constr in self.pointConstr:
            pointCon.append(y[curN: curN + constr.nf])
            curN += constr.nf
        pathCon = []
        for constr in self.pathConstr:
            pathCon.append(np.reshape(y[curN: curN+N*constr.nf], (N, constr.nf)))
            curN += N*constr.nf
        nonLinCon = []
        for constr in self.nonLinConstr:
            nonLinCon.append(y[curN: curN+constr.nf])
            curN += constr.nf
        # check bounds, return a -1, 1 value for non-equality bounds, and 0 for equality bounds
        useX, useU, useP = self.__parseX__(guess)
        Xbound = checkInBounds(useX[:, :dimx], self.xbd)
        x0bound = checkInBounds(useX[0, :dimx], self.x0bd)
        xfbound = checkInBounds(useX[-1, :dimx], self.xfbd) # TODO: support acceleration bound check
        ubound = checkInBounds(useU, self.ubd)
        if self.dimp > 0:
            pbound = checkInBounds(useP, self.pbd)
        else:
            pbound = None
        if self.t0ind > 0:
            t0bound = checkInBounds([guess[self.t0ind]], self.t0)
        else:
            t0bound = None
        if self.tfind > 0:
            tfbound = checkInBounds([guess[self.tfind]], self.tf)
        else:
            tfbound = None
        if self.lenAddX > 0:
            addx = self.__parseAddX__(guess)
            addXbound = [checkInBounds(addx_, [addx__.lb, addx__.ub]) for addx_, addx__ in zip(addx, self.addX)]
        else:
            addXbound = None
        return {'obj': obj, 'dyn': dynCon, 'defect': defectCon, 'point': pointCon, 'path': pathCon, 'nonlin': nonLinCon,
                'Xbd': Xbound, 'Ubd': ubound, 'x0bd': x0bound, 'xfbd': xfbound, 'Pbd': pbound,
                't0bd': t0bound, 'tfbd': tfbound, 'addXbd': addXbound}

    def __parseAddX__(self, x):
        numTraj = self.numTraj
        addX = []
        for addx in self.addX:
            addX.append(x[numTraj: numTraj + addx.n])
            numTraj += addx.n
        return addX

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix.

        :param x: ndarray, the solution to the problem
        :param y: ndarray, return F
        :param G, row, col: ndarray, information of gradient
        :param rec, needg: if we record/ if we need gradient

        """
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        # loop over all the objective functions
        curRow, curNg = self.__objModeG__(0, 0, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        # loop over all system dynamics constraint
        curRow, curNg = self.__dynconstrModeG__(curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
        curRow, curNg = self.__constrModeG__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)

    def __dynconstrModeG__(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        dimpoint = self.dimpoint
        dimdyn = self.dimdyn  # this works for many cases
        nPoint = self.nPoint
        # first let's check the 2*N - 1 dimdyn constraints from dynamics
        cDyn = np.reshape(y[curRow:curRow + nPoint * dimdyn], (nPoint, dimdyn))
        for i in range(nPoint):
            Gpiece = G[curNg: curNg + self.sys.nG]
            rowpiece = row[curNg: curNg + self.sys.nG]
            colpiece = col[curNg: curNg + self.sys.nG]
            self.sys.dyn(useT[i], useX[i], useU[i], useP[i], cDyn[i], Gpiece, rowpiece, colpiece, rec, needg)
            if needg:
                curNg += self.sys.nG
                if rec:
                    rowpiece[:] += curRow
                    colpiece[:] = self.__patchCol__(i, colpiece[:])
            curRow += self.dimdyn
        # offset of row number due to defect dynamics constraint
        curRow += (4*dimdyn + dimu + dimp) * (self.N - 1)
        return curRow, curNg

    def __objModeG__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. G mode. See __constrModeG__ for arguments and output."""
        curRow, curNg = self.__objModeGKnown(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        return curRow, curNg

    def __objModeGKnown(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. G mode. Sparsity structure is known or easy to guess.
        See __constrModeG__ for arguments and output."""
        tmpout = np.zeros(1)
        y[0] = 0
        for obj in self.nonLinObj:  # nonlinear cost function
            Gpiece = G[curNg: curNg + obj.nG]
            rowpiece = row[curNg: curNg + obj.nG]
            colpiece = col[curNg: curNg + obj.nG]
            obj.__callg__(x, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
            y[0] += tmpout[0]
            if needg:
                if rec:
                    rowpiece[:] = curRow
                curNg += obj.nG
        # add lqr obj, which is similar to nonLinObj in this mode
        if self.lqrObj is not None:  # the lqr obj
            Gpiece = G[curNg: curNg + self.LQRnG]
            rowpiece = row[curNg: curNg + self.LQRnG]
            colpiece = col[curNg: curNg + self.LQRnG]
            self.lqrObj(h, useX, useU, useP, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
            y[0] += tmpout[0]
            if needg:
                if rec:
                    rowpiece[:] = curRow
                curNg += self.LQRnG
        return curRow + 1, curNg

    def __objDenseGradient__(self, h, useT, useX, useU, useP, x, y, needg):
        """Calculate a dense gradient even if every piece is given in sparse mode.

        return tmpG: ndarray, dense gradient
        """
        tmpG = np.zeros(self.numSol, dtype=float)  # those three arrays temporarily are used for gradient evaluation
        tmpout = np.zeros(1)
        tmpcallG = np.zeros(self.numSol, dtype=float)  # those three arrays temporarily are used for gradient evaluation
        tmpcallrow = np.zeros(self.numSol, dtype=int)
        tmpcallcol = np.zeros(self.numSol, dtype=int)
        for obj in self.linearObj:  # linear objective function
            y[0] += obj.A.dot(x)
            if needg:
                tmpG[obj.A.col] += obj.A.data
        for obj in self.linPointObj:  # linear point cost
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            y[0] += obj.A.dot(tmpx)
            if needg:
                self.__assignTo__(obj.index, obj.A.col, obj.A.data, tmpG, True)
        for obj in self.linPathObj:  # linear path cost
            # TODO: caveat, current implementation does not support obj on t_k and it is unnecessary
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                tmpval = obj.A.dot(tmpx)
                y[0] += tmpval * h
                if needg:
                    self.__assignTo__(i, obj.A.col, obj.A.data * h, tmpG, True)
                    if self.t0ind > 0:
                        tmpG[self.t0ind] -= tmpval / (self.N - 1)
                    if self.tfind > 0:
                        tmpG[self.t0ind] += tmpval / (self.N - 1)
        for obj in self.nonPointObj:  # nonlinear point cost
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            obj.__callg__(tmpx, tmpout, tmpcallG, tmpcallrow, tmpcallcol, True, needg)
            y[0] += tmpout[0]
            if needg:
                self.__assignTo__(obj.index, tmpcallcol[:obj.nG], tmpcallG[:obj.nG], tmpG, True)
        for obj in self.nonPathObj:  # nonlinear path cost
            # TODO: add support for gradient on tk
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                obj.__callg__(tmpx, tmpout, tmpcallG, tmpcallrow, tmpcallcol, True, needg)
                y[0] += tmpout[0] * h
                if needg:
                    self.__assignTo__(i, tmpcallcol[:obj.nG], tmpcallG[:obj.nG] * h, tmpG, True)
                    if self.t0ind > 0:
                        tmpG[self.t0ind] -= tmpout[0] / (self.N - 1)
                    if self.tfind > 0:
                        tmpG[self.t0ind] += tmpout[0] / (self.N - 1)
        for obj in self.nonLinObj:  # nonlinear cost function
            obj.__callg__(x, tmpout, tmpcallG, tmpcallrow, tmpcallcol, True, needg)
            y[0] += tmpout[0]
            if needg:
                self.__assignTo__(-1, tmpcallcol[:obj.nG], tmpcallG[:obj.nG], tmpG, True)
        # add lqr obj, which is similar to nonLinObj in this mode
        if self.lqrObj is not None:  # the lqr obj
            self.lqrObj(h, useX, useU, useP, tmpout, tmpcallG, tmpcallrow, tmpcallcol, True, needg)
            y[0] += tmpout[0]
            if needg:
                self.__assignTo__(-1, tmpcallcol[:self.LQRnG], tmpcallG[:self.LQRnG], tmpG, True)
        return tmpG

    def __constrModeG__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate constraint function. G mode

        :param curRow: int, index from which we write on
        :param curNg: int, current index in G
        :param h, useT, useX, useU, useP: parsed solution
        :param y: ndarray, the F to be written
        :param G, row, col: ndarray, the G to be written and the locations
        :param rec: bool, if we record row and col
        :param needg: bool, if we need gradient information
        :returns: curRow: current row after we write on y
        :returns: curNg: current index in G after this

        """
        # loop over other constraints
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
        for constr in self.pathConstr:
            for j in range(self.N):
                i = 2 * j
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
    def ipEvalF(self, x):
        """The eval_f function required by ipopt.

        :param x: a guess/solution of the problem
        :return f: float, objective function

        """
        y = np.zeros(1)
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        G = np.zeros(1)
        row = np.zeros(1, dtype=int)
        self.__objModeG__(0, 0, h, useT, useX, useU, useP, x, y, G, row, row, False, False)
        return y[0]

    def ipEvalGradF(self, x):
        """Evaluation of the gradient of objective function.

        :param x: guess/solution to the problem
        :return grad: gradient of objective function w.r.t x

        """
        h, useT = self.__getTimeGrid(x)
        useX, useU, useP = self.__parseX__(x)
        y = np.zeros(1)
        if not self.objSparseMode:
            tmpG = self.__objDenseGradient__(h, useT, useX, useU, useP, x, y, True)
        else:
            tmpG = np.zeros(self.numSol)
            spG = np.zeros(self.numObjG)
            spRow = np.zeros(self.numObjG, dtype=int)
            spCol = np.zeros(self.numObjG, dtype=int)
            self.__objModeGKnown__(0, 0, h, useT, useX, useU, useP, x, y, spG, spRow, spCol, True, True)
            tmpG[spCol] = spG
        return tmpG

    def ipEvalG(self, x):
        """Evaluation of the constraint function.

        :param x: ndarray, guess/solution to the problem.
        :return g: constraint function

        """
        y = np.zeros(self.numF)
        if self.gradmode:
            G = np.zeros(1)
            row = np.zeros(1, dtype=int)
            col = np.zeros(1, dtype=int)
            self.__callg__(x, y, G, row, col, False, False)
            return y
        else:
            self.__callf__(x, y)
        return y

    def ipEvalJacG(self, x, flag):
        """Evaluate jacobian of constraints. I simply call __callg__

        :param x: ndarray, guess / solution to the problem
        :param flag: bool, True return row/col, False return values

        """
        y = np.zeros(self.numF)
        G = np.zeros(self.nG)
        if flag:
            row = np.ones(self.nG, dtype=int)
            col = np.ones(self.nG, dtype=int)
            self.__callg__(x, y, G, row, col, True, True)
            return row, col
        else:
            row = np.ones(1, dtype=int)
            col = np.ones(1, dtype=int)
            self.__callg__(x, y, G, row, col, False, True)
            return G

    def __patchCol__(self, index, col):
        """Find which indices it belongs to the original one for a local matrix at col.

        Since we have changed how variables are arranged, now it should be quite straightforward to do so.
        """
        # TODO: with time, there is still issue
        return col - 1 + index * self.dimpoint

    def parseSol(self, sol):
        """Call parseX function from utility and return a dict of solution."""
        X, U, P = self.__parseX__(sol)
        if self.dimp == 0:
            P = None
        h, tgrid = self.__getTimeGrid(sol)
        if self.lenAddX == 0:
            return {'t': tgrid, 'x': X, 'u': U, 'p': P}
        else:
            return {'t': tgrid, 'x': X, 'u': U, 'p': P, 'addx': self.__parseAddX__(sol)}

    def addLQRObj(self, lqrobj):
        """Add a lqr objective function to the problem. It changes lqrObj into a function being called in two modes...

        :param lqrobj: a lqrObj class.

        """
        Fcol = lqrobj.F.col
        Qcol = lqrobj.Q.col
        Rcol = lqrobj.R.col
        useF = len(Fcol)
        useQ = len(Qcol)
        useR = len(Rcol)
        if useF > 0 and useQ > 0:
            assert np.allclose(Fcol, Qcol)
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
        nPoint = self.nPoint
        num1 = nPoint * (useQ + useR + useP)
        self.LQRnG = num1 + self.numT
        if useF > 0 and useQ == 0:
            self.LQRnG += useF
        weight = np.zeros((nPoint, 1))
        weight[1::2] = 2.0 / 3.0
        weight[0::2] = 1.0 / 3.0
        weight[0] = 1.0 / 6.0
        weight[-1] = 1.0 / 6.0
        baseCol = self.dimpoint * np.arange(self.nPoint)[:, np.newaxis]  # a nPoint by 1 column matrix

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
                col_ = np.reshape(col[:num1], (nPoint, -1))
            yF = 0.0
            yQ = 0.0
            yR = 0.0
            yP = 0.0
            yTf = tfweight * (h * (self.N - 1))
            curG = 0
            if useQ > 0:
                yQ = np.dot(weight[:, 0], np.sum(((useX[:, Qcol] - lqrobj.xbase[Qcol]) ** 2) * self.lqrobj.Q.data, axis=1)) * h
                if needg:
                    G[curG: curG + nPoint * useQ] = 2.0 * h * (weight * ((useX[:, Qcol] - lqrobj.xbase[Qcol]) * lqrobj.Q.data)).flatten()
                    if rec:
                        col_[:, :useQ] = Qcol + baseCol
                    curG += nPoint * useQ
            if useR > 0:
                yR = np.dot(weight[:, 0], np.sum(((useU[:, Rcol] - lqrobj.ubase[Rcol]) ** 2) * lqrobj.R.data, axis=1)) * h
                if needg:
                    G[curG: curG + nPoint * useR] = 2.0 * h * (weight * ((useU[:, Rcol] - lqrobj.ubase[Rcol]) * lqrobj.R.data)).flatten()
                    if rec:
                        col_[:, useQ: useQ+useR] = Rcol + baseCol + self.dimx
                    curG += nPoint * useR
            if useP > 0:
                yP = np.dot(weight[:, 0], lqrobj.P.data * np.sum((useP_[:, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0)) * h
                if needg:
                    G[curG: curG + nPoint * useP] = 2.0 * h * (weight * ((useP_[:, Pcol] - lqrobj.pbase[Pcol]) * lqrobj.P.data)).flatten()
                    if rec:
                        col_[:, useQ+useR: useQ+useR+useP] = Pcol + baseCol + self.dimx + self.dimu
                    curG += nPoint * useP
            if useF > 0:
                yF = np.sum(lqrobj.F.data * ((useX[-1, Fcol] - lqrobj.xfbase[Fcol]) ** 2))
                if useQ > 0:
                    if needg:
                        n0 = curG - numPoint * (useR + useP) - useF  # locate at the last row
                        G[n0: n0 + useF] += 2.0 * lqrobj.F.data * (useX[-1, Fcol] - lqrobj.xfbase[Fcol])
                else:
                    G[curG: curG + useF] = 2.0 * lqrobj.F.data * (useX[-1, Fcol] - lqrobj.xfbase[Fcol])
                    if rec:
                        col[curG: curG + useF] = Fcol + baseCol[-1]
                    curG += useF
            if needg:
                if self.t0ind > 0:
                    G[curG] = -(yQ + yR + yP) / h / (self.N - 1) - tfweight
                    if rec:
                        row[curG: curG+1] = 0
                        col[curG: curG+1] = self.t0ind
                    curG += 1
                if self.tfind > 0:
                    G[curG] = (yQ + yR + yP) / h / (self.N - 1) + tfweight
                    if rec:
                        row[curG: curG+1] = 0
                        col[curG: curG+1] = self.tfind
                    curG += 1
            y[0] = yF + yQ + yR + yP + yTf

        self.lqrObj = __callg__

    def addLinearObj(self, linObj):
        """Add linear objective function.

        :param linObj: linearObj class

        """
        assert isinstance(linObj, linearObj)
        self.linearObj.append(linObj)

    def addLinPointObj(self, linPointObj, path=False):
        """Add linear point objective function.

        :param linPointObj: linearPointObj class
        :param path: bool, if this is path obj (at every point except for final one)

        """
        assert isinstance(linPointObj, linearPointObj)
        if path:
            self.linPathObj.append(linPointObj)
        else:
            self.linPointObj.append(linPointObj)

    def addNonLinObj(self, nonlinObj):
        """Add nonlinear objective function.

        :param nonLinObj: a nonLinObj class

        """
        assert isinstance(nonlinObj, nonLinObj)
        self.nonLinObj.append(nonlinObj)

    def addNonPointObj(self, nonPntObj, path=False):
        """Add nonlinear point objective.

        :param nonPntObj: nonLinObj class
        :param path: bool, if this obj is pointwise

        """
        assert isinstance(nonPntObj, nonPointObj)
        if path:
            self.nonPathObj.append(nonPntObj)
        else:
            self.nonPointObj.append(nonPntObj)

    def addPointConstr(self, pntConstr, path=False):
        """Add point constraint.

        :param pntConstr: pointConstr class
        :param path: bool, if this obj

        """
        assert isinstance(pntConstr, pointConstr)
        if path:
            self.pathConstr.append(pntConstr)
        else:
            self.pointConstr.append(pntConstr)

    def addNonLinConstr(self, constr):
        """Add a general nonlinear constraint.

        :param constr: nonLinConstr class

        """
        assert isinstance(constr, nonLinConstr)
        self.nonLinConstr.append(constr)

    def setN(self, N):
        """Set N.

        :param N: the size of discretization.

        """
        self.N = N

    def sett0tf(self, t0, tf):
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
