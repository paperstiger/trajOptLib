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
import numpy as np
from .trajOptBase import linearObj, linearPointObj
from .trajOptBase import linearPointConstr, linearConstr
from .trajOptBase import nonLinearPointObj, nonLinearObj
from .trajOptBase import nonLinearPointConstr, nonLinearConstr
from .trajOptBase import system, addX
from .trajOptBase import lqrObj
from .libsnopt import snoptConfig, probFun, solver
from .utility import parseX, randomGenInBound, checkInBounds
from scipy import sparse
from scipy.sparse import spmatrix, coo_matrix


class trajOptProblem(probFun):
    """A class for definition of trajectory optimization problem.

    A general framework for using this class is to:

    1. Define a class inherited from system and write dyn/Jdyn method.
    2. Optionally, write desired cost function by inheriting/creating from the list of available cost functions.
    3. Optionally, write desired constraint functions by inheriting from available constraints.
    4. Create this class with selected system, discretization, t0, tf range, gradient option
    5. Set bounds for state, control, parameters, x0 and xf
    6. Add objective functions and constraints to this class
    7. Call preProcess method explicitly
    8. Create snoptConfig instance and choose desired options
    9. Construct the solver
    10. Use the solver to solve with either automatic guess or user provided guess

    """
    def __init__(self, sys, N=20, t0=0.0, tf=1.0, gradmode=True, addx=None):
        """Initialize problem by system, discretization grid size, and allowable time

        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :param gradmode: bool, sets if we use gradient mode.
        :param addX: list of addX / one addX / None, additional optimization variables.

        """
        assert isinstance(sys, system)
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
        self.pathConstr = []  # general constraint imposed everywhere such as collision avoidance
        self.nonLinConstr = []  # stores general nonlinear constraint
        self.linPointConstr = []
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

    def preProcess(self):
        """Initialize the instances of probFun now we are ready.

        Call this function after the objectives and constraints have been set appropriately.
        It calculate the space required for SNOPT and allocates sparsity structure if necessary.

        """
        numDyn = self.dimx * (self.N - 1)  # constraints from system dynamics
        numC = 0
        for constr in self.pointConstr:
            numC += constr.nf
        for constr in self.pathConstr:
            numC += self.N * constr.nf
        for constr in self.nonLinConstr:
            numC += constr.nf
        nnonlincon = numC
        for constr in self.linPointConstr:
            numC += constr.A.shape[0]
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
        probFun.__init__(self, self.numSol, self.numF)  # currently we do not know G yet
        self.__setAPattern(numDyn, nnonlincon, spA)
        self.__setXbound()
        self.__setFbound()
        # detect gradient information
        if self.gradmode:  # in this case, we randomly generate a guess and use it to initialize everything
            randX = self.randomGenX()
            self.__turnOnGrad(randX)

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
        catArow = np.concatenate((np.zeros(len(nnzind)), row_))
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
        curNA = len(self.Aval)  # this is just for bookkeeping
        return curRow, curNA


    def randomGenX(self):
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
        ydot, Jac = self.sys.Jdyn(useT[0], useX[0], useU[0], useP[0])
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
                Jac = Jac.tocoo()
                if self.sys.ode == 'Euler':
                    Jac.data *= h
                    # add a block by identity
                    Jac = Jac + eyemat
                elif self.sys.ode == 'BackEuler':
                    Jac.data *= (-h)
                    Jac = Jac + eyemat
                dynnnz = Jac.nnz - timennz
                dynG = (self.N - 1) * (dynnnz + self.dimx + self.numT)
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
            nG += self.nPoint * obj.nG
        return nG


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
        numX = self.numX
        numU = self.numU
        numP = self.numP
        xlb = np.zeros(self.numSol)
        xub = np.zeros(self.numSol)
        if self.xbd[0] is not None:
            stateLb = np.reshape(xlb[:numX], (self.N, self.dimx))
            stateLb[:] = self.xbd[0]
            # set lb for x0 and xf
            if self.x0bd[0] is not None:
                stateLb[0] = self.x0bd[0]
            if self.xfbd[0] is not None:
                stateLb[self.N - 1] = self.xfbd[0]
        if self.xbd[1] is not None:
            stateUb = np.reshape(xub[:numX], (self.N, self.dimx))
            stateUb[:] = self.xbd[1]
            # set ub for x0 and xf
            if self.x0bd[1] is not None:
                stateUb[0] = self.x0bd[1]
            if self.xfbd[1] is not None:
                stateUb[self.N - 1] = self.xfbd[1]
        if self.ubd[0] is not None:
            ctrlLb = np.reshape(xlb[numX:numX + numU], (self.N, self.dimu))
            ctrlLb[:] = self.ubd[0]
        if self.ubd[1] is not None:
            ctrlUb = np.reshape(xub[numX:numX + numU], (self.N, self.dimu))
            ctrlUb[:] = self.ubd[1]
        if self.pbd[0] is not None and self.dimp > 0:
            pLb = np.reshape(xlb[numX+numU:numX+numU+numP], (self.N, self.dimp))
            pLb[:] = self.pbd[0]
        if self.pbd[1] is not None and self.dimp > 0:
            pUb = np.reshape(xub[numX+numU:numX+numU+numP], (self.N, self.dimp))
            pUb[:] = self.pbd[1]
        if not self.fixt0:
            xlb[numX+numU+numP] = self.t0[0]
            xub[numX+numU+numP] = self.t0[1]
        # set bound on time
        if not self.fixtf:
            if not self.fixt0:
                xlb[numX+numU+numP+1] = self.tf[0]
                xub[numX+numU+numP+1] = self.tf[1]
            else:
                xlb[numX+numU+numP] = self.tf[0]
                xub[numX+numU+numP] = self.tf[1]
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
        self.xlb = xlb
        self.xub = xub

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
        h = (usetf - uset0) / (self.N - 1)
        useT = np.linspace(uset0, usetf, self.N)
        return h, useT

    def __parseX__(self, x):
        """Parse guess/sol into X, U, P"""
        numX = self.numX
        numU = self.numU
        numP = self.numP
        useX = np.reshape(x[:numX], (self.N, self.dimx))
        useU = np.reshape(x[numX: numX + numU], (self.N, self.dimu))
        useP = np.reshape(x[numX + numU: numX + numU + numP], (self.N, self.dimp))
        return useX, useU, useP

    def parseF(self, guess):
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
        dynCon = np.reshape(y[1:(N-1)*dimx+1], (N - 1, dimx))
        curN = 1 + (N - 1) * dimx
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
        Xbound = checkInBounds(useX, self.xbd)
        x0bound = checkInBounds(useX[0], self.x0bd)
        xfbound = checkInBounds(useX[-1], self.xfbd)
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
        return {'obj': obj, 'dyn': dynCon, 'point': pointCon, 'path': pathCon, 'nonlin': nonLinCon,
                'Xbd': Xbound, 'Ubd': ubound, 'x0bd': x0bound, 'xfbd': xfbound, 'Pbd': pbound,
                't0bd': t0bound, 'tfbd': tfbound, 'addXbd': addXbound}

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

    def __objModeF__(self, curRow, h, useT, useX, useU, useP, x, y):
        """Calculate objective function. F mode

        :param curRow: int, index from which we write on
        :param h, useT, useX, useU, useP: parsed solution
        :param x: ndarray, the sol, it is used for linear constraints
        :param y: ndarray, the F to be written. The first row stores the objective function

        """
        y[0] = 0.0
        tmpout = np.zeros(1)
        for obj in self.linearObj:
            y[0] += obj.A.dot(x)
        for obj in self.linPointObj:
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            y[0] += obj.A.dot(tmpx)
        for obj in self.nonLinObj:
            obj.__callf__(x, tmpout)
            y[0] += tmpout[0]
        for obj in self.linPathObj:
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                obj.__callf__(tmpx, tmpout)
                y[0] += tmpout[0] * h
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
            obj.__callf__(x, tmpout)
            y[0] += tmpout[0]
        # add lqr cost, if applicable
        if self.lqrObj is not None:
            y[0] += self.lqrObj(h, useX, useU, useP)

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
        for constr in self.pathConstr:
            for i in range(self.N):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                constr.__evalf__(tmpx, y[curRow: curRow + constr.nf])
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
        cDyn = np.reshape(y[curRow:curRow+(self.N - 1) * self.dimx], (self.N - 1, self.dimx))
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

    def __dynconstrModeG(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        numX, numU = self.numX, self.numU
        cDyn = np.reshape(y[curRow:curRow+(self.N - 1) * self.dimx], (self.N - 1, self.dimx))
        for i in range(self.N - 1):
            # evaluate gradient of system dynamics TODO: support other types of integration scheme
            ydot, Jac = self.sys.Jdyn(useT[i], useX[i], useU[i], useP[i])  # TODO: support in-place Jacobian
            cDyn[i] = useX[i] + h * ydot - useX[i + 1]
            if needg:
                if not self.dynSparse:
                    Jac *= h  # always useful
                    # assign a block for x
                    G[curNg: curNg + dimx*dimx] = (Jac[:, 1:1+dimx] + np.eye(dimx)).flatten()
                    if rec:
                        tmpMat = np.tile(np.arange(dimx), (dimx, 1))
                        row[curNg: curNg + dimx*dimx] = curRow + tmpMat.T.flatten()
                        col[curNg: curNg + dimx*dimx] = i * dimx + tmpMat.flatten()
                    curNg += dimx * dimx
                    # assign a block for u
                    G[curNg: curNg + dimx*dimu] = Jac[:, 1+dimx:1+dimx+dimu].flatten()
                    if rec:
                        row[curNg: curNg + dimx*dimu] = curRow + np.tile(np.arange(dimx), (dimu, 1)).T.flatten()
                        col[curNg: curNg + dimx*dimu] = numX + i * dimu + np.tile(np.arange(dimu), dimx).flatten()
                    curNg += dimx * dimu
                    # assign a block for p, if necessary
                    if dimp > 0:
                        G[curNg: curNg + dimx*dimp] = Jac[:, 1+dimx+dimu:1+dimx+dimu+dimp].flatten()
                        if rec:
                            row[curNg: curNg + dimx*dimp] = curRow + np.tile(np.arange(dimx), (dimp, 1)).T.flatten()
                            col[curNg: curNg + dimx*dimp] = numX + numU + i * dimp + np.tile(np.arange(dimp), dimx).flatten()
                        curNg += dimx * dimp
                    # assign the diagonal block for x_{k+1}
                    G[curNg: curNg + dimx] = -1.0
                    if rec:
                        row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                        col[curNg: curNg + dimx] = np.arange((i+1)*dimx, (i+2)*dimx)
                    curNg += dimx
                    # assign a column for t0, if necessary
                    if self.t0ind > 0:
                        G[curNg: curNg + dimx] = -ydot / (self.N - 1) + Jac[:, 0] * (1-float(i)/(self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.t0ind
                        curNg += dimx
                    # assign a column for tf, if necessary
                    if self.tfind > 0:
                        G[curNg: curNg + dimx] = ydot / (self.N - 1) + Jac[:, 0] * (float(i)/(self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.tfind
                        curNg += dimx
                else:
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
                        col[curNg: curNg + dimx] = np.arange((i+1)*dimx, (i+2)*dimx)
                    curNg += dimx
                    if self.t0ind > 0:
                        G[curNg: curNg + dimx] = -ydot/(self.N - 1) + Jac.getcol(0).toarray() * (1 - float(i)/(self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.t0ind
                        curNg += dimx
                    if self.tfind > 0:
                        G[curNg: curNg + dimx] = ydot/(self.N - 1) + Jac.getcol(0).toarray() * (float(i)/(self.N - 1.0))
                        if rec:
                            row[curNg: curNg + dimx] = curRow + np.arange(dimx)
                            col[curNg: curNg + dimx] = self.tfind
                        curNg += dimx
            curRow += dimx
        return curRow, curNg

    def __objModeG(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. It just evaluates them and assign to correct position in y.

        See __constrModeG for arguments and output.
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
        if len(self.nonPointObj) > 0:
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

        if len(self.nonPathObj) > 0:
            for obj in self.nonPathObj:
                y[curRow] = 0
                for i in range(self.nPoint):
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

        if(len(self.nonLinObj)) > 0:
            for obj in self.nonLinObj:  # nonlinear cost function
                Gpiece = G[curNg: curNg + obj.nG]
                rowpiece = row[curNg: curNg + obj.nG]
                colpiece = col[curNg: curNg + obj.nG]
                obj.__callg__(x, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                y[curRow] = tmpout[0]
                if rec:
                    rowpiece[:] = curRow
                curNg += obj.nG
                curRow += 1
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
        if len(self.pointConstr) > 0:
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
        if len(self.pathConstr) > 0:
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
        if len(self.nonLinConstr) > 0:
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
        """Find which indices it belongs to the original one for a local matrix at col"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        newcol = col[col > 0].copy()  # TODO: check if __patchCol__ causes error when time is involved.
        maskx = (col > 0) & (col < 1 + dimx)
        newcol[maskx] = index * self.dimx + col[maskx] - 1
        masku = (col >= 1 + dimx) & (col < 1 + dimx + dimu)
        newcol[masku] = self.numX + index * self.dimu + col[masku] - 1 - dimx
        if self.dimp > 0:
            maskp = (col >= 1 + dimx + dimu) & (col < 1 + dimx + dimu + dimp)
            newcol[maskp] = self.numX + self.numU + index * self.dimp + col[maskp] - 1 - dimx - dimu
        return newcol

    def parseSol(self, sol):
        """Call parseX function from utility and return a dict of solution."""
        return parseX(sol, self.N, self.dimx, self.dimu, self.dimp, self.t0, self.tf, True)

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
            if not needg:
                y[0] = __callf__(h, useX, useU, useP_)
            else:
                yF = 0.0
                yQ = 0.0
                yR = 0.0
                yP = 0.0
                yTf = tfweight * (h * (self.N - 1))
                curG = 0
                if useQ > 0:
                    yQ = np.sum(lqrobj.Q.data * np.sum((useX[:-1, Qcol] - lqrobj.xbase[Qcol]) ** 2, axis=0)) * h
                    G[curG: curG + (self.N - 1) * useQ] = 2.0 * h * ((useX[:-1, Qcol] - lqrobj.xbase[Qcol]) * lqrobj.Q.data).flatten()
                    if rec:
                        row[curG: curG + useQ * (self.N - 1)] = 0
                        col[curG: curG + useQ * (self.N - 1)] = (Qcol + self.dimx*(np.arange(0, self.N - 1)[:, np.newaxis])).flatten()
                    curG += (self.N - 1) * useQ
                if useR > 0:
                    yR = np.sum(lqrobj.R.data * np.sum((useU[:-1, Rcol] - lqrobj.ubase[Rcol]) ** 2, axis=0)) * h
                    G[curG: curG + (self.N - 1) * useR] = 2.0 * h * ((useU[:-1, Rcol] - lqrobj.ubase[Rcol]) * lqrobj.R.data).flatten()
                    if rec:
                        row[curG: curG + useR * (self.N - 1)] = 0
                        col[curG: curG + useR * (self.N - 1)] = (self.numX + Rcol + self.dimu*(np.arange(0, self.N - 1)[:, np.newaxis])).flatten()
                    curG += (self.N - 1) * useR
                if useP > 0:
                    yP = np.sum(lqrobj.P.data * np.sum((useP_[:-1, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0)) * h
                    G[curG: curG + (self.N - 1) * useP] = 2.0 * h * ((useP_[:-1, Pcol] - lqrobj.pbase[Pcol]) * lqrobj.P.data).flatten()
                    if rec:
                        row[curG: curG + useP * (self.N - 1)] = 0
                        col[curG: curG + useP * (self.N - 1)] = (self.numX + self.numU + Pcol + self.dimp*(np.arange(0, self.N - 1)[:, np.newaxis])).flatten()
                    curG += (self.N - 1) * useP
                if useF > 0:
                    yF = np.sum(lqrobj.F.data * ((useX[-1, Fcol] - lqrobj.xfbase[Fcol]) ** 2))
                    G[curG: curG + useF] = 2.0 * lqrobj.F.data * (useX[-1, Fcol] - lqrobj.xfbase[Fcol])
                    if rec:
                        row[curG: curG + useF] = 0
                        col[curG: curG + useF] = np.arange(self.numX - self.dimx, self.numX)
                    curG += useF
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

        if self.gradmode:
            self.lqrObj = __callg__
        else:
            self.lqrObj = __callf__

    def addLinearObj(self, linObj):
        """Add linear objective function.

        :param linObj: linearObj class

        """
        assert isinstance(linObj, linearObj)
        self.linearObj.append(linObj)

    def addLinearPointObj(self, linPointObj, path=False):
        """Add linear point objective function.

        :param linPointObj: linearPointObj class
        :param path: bool, if this is path obj (at every point except for final one)

        """
        assert isinstance(linPointObj, linearPointObj)
        if path:
            self.linPathObj.append(linPointObj)
        else:
            self.linPointObj.append(linPointObj)

    def addNonLinearObj(self, nonlinObj):
        """Add nonlinear objective function.

        :param nonLinObj: a nonLinObj class

        """
        assert isinstance(nonlinObj, nonLinearObj)
        self.nonLinObj.append(nonlinObj)

    def addNonLinearPointObj(self, nonPntObj, path=False):
        """Add nonlinear point objective.

        :param nonPntObj: nonLinObj class
        :param path: bool, if this obj is pointwise

        """
        assert isinstance(nonPntObj, nonLinearPointObj)
        if path:
            self.nonPathObj.append(nonPntObj)
        else:
            self.nonPointObj.append(nonPntObj)

    def addNonLinearPointConstr(self, pntConstr, path=False):
        """Add point constraint.

        :param pntConstr: pointConstr class
        :param path: bool, if this obj

        """
        assert isinstance(pntConstr, nonLinearPointConstr)
        if path:
            self.pathConstr.append(pntConstr)
        else:
            self.pointConstr.append(pntConstr)

    def addNonLinearConstr(self, constr):
        """Add a general nonlinear constraint.

        :param constr: nonLinConstr class

        """
        assert isinstance(constr, nonLinearConstr)
        self.nonLinConstr.append(constr)

    def addLinearConstr(self, constr):
        assert isinstance(constr, linearConstr)
        self.linearConstr.append(constr)

    def addLinearPointConstr(self, constr, path=False):
        assert isinstance(constr, linearPointConstr)
        if path:
            self.linPathConstr.append(constr)
        else:
            self.linPointConstr.append(constr)

    def addObj(self, obj, path=False):
        """A high level function that add objective function of any kind."""
        if isinstance(obj, linearObj):
            self.addLinearObj(obj)
        elif isinstance(obj, linearPointObj):
            self.addLinearPointObj(obj, path)
        elif isinstance(obj, nonLinearObj):
            self.addNonLinearObj(obj)
        elif isinstance(obj, nonLinearPointObj):
            self.addNonLinearPointObj(obj, path)
        elif isinstance(obj, lqrObj):
            self.addLQRObj(obj)

    def addConstr(self, constr, path=False):
        if isinstance(constr, linearConstr):
            self.addLinearConstr(constr)
        elif isinstance(constr, linearPointConstr):
            self.addLinearPointConstr(constr, path)
        elif isinstance(constr, nonLinearConstr):
            self.addNonLinearConstr(constr)
        elif isinstance(constr, nonLinearPointConstr):
            self.addNonLinearPointConstr(constr, path)

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
