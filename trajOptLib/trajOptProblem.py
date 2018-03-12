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
import logging
from .trajOptBase import system, linearObj, linearPointObj, nonLinObj, nonPointObj, pointConstr, nonLinConstr, lqrObj
from .libsnopt import snoptConfig, probFun, solver
from .utility import parseX
from scipy import sparse
from scipy.sparse import spmatrix, coo_matrix


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class trajOptProblem(probFun):
    """A class for definition of trajectory optimization problem.

    A general framework for using this class is to:

    1. Define a class inherited from system and write dyn/Jdyn method.
    2. Optionally, write desired cost function by inheriting/creating from the list of available cost functions.
    3. Optinally, write desired constraint functions by inheriting from available constraints.
    4. Create this class with selected system, discretization, t0, tf range, gradient option
    5. Set bounds for state, control, parameters, x0 and xf
    6. Add objective functions and constraints to this class
    7. Call preProcess method explicitly
    8. Create snoptConfig instance and choose desired options
    9. Construct the solver
    10. Use the solver to solve with either automatic guess or user provided guess

    :exclude-members: getSparsity

    """
    def __init__(self, sys, N=20, t0=0.0, tf=1.0, gradmode=True):
        """Initialize problem by system, discretization grid size, and allowable time

        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :param gradmode: bool, sets if we use gradient mode.

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
        # calculate number of variables to be optimized, time are always the last
        numX = self.N * self.dimx
        numU = self.N * self.dimu
        numP = self.N * self.dimp
        numT = 2
        if self.fixt0:
            numT -= 1
        if self.fixtf:
            numT -= 1
        self.t0ind, self.tfind = self.__getTimeIndices__()
        numSol = numX + numU + numP + numT
        self.numX = numX
        self.numU = numU
        self.numP = numP
        self.numT = numT
        self.numSol = numSol

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
        self.numF = 1 + numDyn + numC
        super(trajOptProblem, self).__init__(self.numSol, self.numF)  # not providing G means we use finite-difference
        self.setXbound()
        self.setFbound()
        if self.gradmode:  # in this case, we randomly generate a guess and use it to initialize everything
            randX = self.randomGenX()
            self.__turnOnGrad__(randX)

    def __turnOnGrad__(self, x0):
        """Turn on gradient, this is called after an initial x0 has been generated"""
        self.grad = True
        self.getSparsity(x0)

    def getSparsity(self, x0):
        arg1, objSparseMode = self.getObjSparsity(x0)
        if objSparseMode:
            numObjG = arg1
        else:
            self.indices = self.getObjSparsity(x0)
            numObjG = len(self.indices)  # G from objective, assume dense
        self.objSparseMode = objSparseMode
        # summarize number of pure linear constraints
        numDynG = self.getDynSparsity(x0)
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

    def getDynSparsity(self, x):
        """Set sparsity of the problem caused by system dynamics and other nonlinear constraints.

        :param x: ndarray, the guess/sol

        """
        # TODO: now sparsity is supported for Euler/ BackEuler/ Dis, RK4 needs to be done
        h, useT = self.__getTimeGrid__(x)
        useX, useU, useP = self.__parseX__(x)
        # use one time/state/ctrl/para/h set to detect the sparsity pattern
        ydot, Jac = self.sys.Jdyn(useT[0], useX[0], useU[0], useP[0])
        nrow, ncol = Jac.shape
        if isinstance(Jac, np.ndarray):
            nnz = nrow * ncol
            self.dynSparse = False
        elif isinstance(Jac, spmatrix):
            nnz = Jac.nnz
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
                timennz = Jac.getcol(0).nnz
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

    def getObjSparsity(self, x):
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
        h, useT = self.__getTimeGrid__(x)
        useX, useU, useP = self.__parseX__(x)
        # check sparseObj mode
        sumLenObj = len(self.nonLinObj) + len(self.nonPointObj) + len(self.nonPathObj) + len(self.linearObj) + len(self.linPointObj) + len(self.linPathObj)
        if sumLenObj == 1 and self.lqrObj is None:
            self.sparseObjMode = True
        elif sumLenObj == 0 and self.lqrObj is not None:
            self.sparseObjMode = True
        else:
            self.sparseObjMode = False
        # for sparseObjMode, we return #nnz as we will add to G, col will be added in __callg__
        if self.sparseObjMode:
            if sumLenObj == 1:
                nG = 0
                for obj in self.nonLinObj:
                    nG += obj.nG
                for obj in self.nonPointObj:
                    nG += obj.nG
                for obj in self.nonPathObj:
                    nG += (self.N - 1) * obj.nG
                for obj in self.linearObj:
                    nG += obj.nG
                for obj in self.linPointObj:
                    nG += obj.nG
                for obj in self.linPathObj:
                    nG += (self.N - 1) * obj.nG  # TODO: determine if it is better practice to merge lin and nonLin path/point together
                return nG, True
            else:
                return self.LQRnG, True
        # for non-sparseObjMode, we use a long vector to remember where the gradients of objective functions are
        # a long vector for temporary use
        G = np.zeros(self.numSol, dtype=float)
        row = np.zeros(self.numSol, dtype=int)
        col = np.zeros(self.numSol, dtype=int)
        tmpF = np.zeros(1)
        # for objective function
        mask = np.zeros(self.numSol, dtype=bool)
        t0ind, tfind = self.__getTimeIndices__()
        for obj in self.linearObj:  # for linear objective over all x
            if not isinstance(obj.A, coo_matrix):
                obj.A = obj.A.tocoo()
            mask[obj.A.nonzero()] = True
        for obj in self.linPointObj:  # for linear objective imposed at selected points
            if not isinstance(obj.A, coo_matrix):
                obj.A = obj.A.tocoo()
            ind = obj.index
            nnzind = obj.A.nonzero()
            self.__assignTo__(ind, nnzind, True, mask, False)
        for obj in self.linPathObj:  # for linear objective imposed at every point
            if not isinstance(obj.A, coo_matrix):
                obj.A = obj.A.tocoo()
            nnzind = obj.A.nonzero()
            for ind in range(self.N - 1):
                self.__assignTo__(ind, nnzind, True, mask, False)
            if self.t0ind > 0:
                mask[self.t0ind] = True  # for path obj, these two are always True
            if self.tfind > 0:
                mask[self.tfind] = True
        for obj in self.nonLinObj:  # for general nonlinear objective function
            obj.__callg__(x, tmpF, G, row, col, True, True)
            nG = obj.nG  # first nG elements are effective, col stores information
            self.__assignTo__(-1, col[:nG], True, mask, False)
        for obj in self.nonPointObj:  # for nonlinear constr imposed at selected points
            ind = obj.index
            nG = obj.nG
            piecex = np.concatenate((useT[ind:ind+1], useX[ind], useU[ind], useP[ind]))
            obj.__callg__(piecex, tmpF, G, row, col, True, True)
            self.__assignTo__(ind, col[:nG], True, mask, False)
        for obj in self.nonPathObj:  # for path, we only need to evaluate once
            nG = obj.nG
            ind = 0  # arbitrarily use the first one to get the sparsity structure
            piecex = np.concatenate((useT[ind:ind+1], useX[ind], useU[ind], useP[ind]))
            obj.__callg__(piecex, tmpF, G, row, col, True, True)
            for ind in range(self.N - 1):
                self.__assignTo__(ind, col[:nG], True, mask, False)
            if self.t0ind > 0:
                mask[self.t0ind] = True  # since it is path obj
            if self.tfind > 0:
                mask[self.tfind] = True
        # determine what lqr objective tells us
        if self.lqrObj is not None:
            nG = self.LQRnG
            self.lqrObj(h, useX, useU, useP, tmpF, G, row, col, True, True)
            self.__assignTo__(-1, col[:nG], True, mask, False)
            if self.t0ind > 0:
                mask[self.t0ind] = True  # since it is lqr obj
            if self.tfind > 0:
                mask[self.tfind] = True
        # summarize # nnz for obj
        nnzObj = np.sum(mask)
        logger.debug("#nnz from obj is %d" % nnzObj)
        indices = np.where(mask)[0]  # this means what indices we should choose from
        return indices, False

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

    def __getTimeIndices__(self):
        """Utility function for assigning sparsity structure."""
        t0ind = -1
        tfind = -1
        if self.fixt0:
            if self.fixtf:
                pass
            else:
                tfind = self.numSol - 1
        else:
            if self.fixxf:
                t0ind = self.numSol - 1
            else:
                t0ind = self.numSol - 2
                tfind = self.numSol - 1
        return t0ind, tfind

    def setXbound(self):
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
        # assign to where it should belong to
        self.xlb = xlb
        self.xub = xub

    def setFbound(self):
        """Set bound on F"""
        # set bound on F
        numF = self.numF
        numDyn = self.dimx * (self.N - 1)  # constraints from system dynamics
        clb = np.zeros(numF)
        cub = np.zeros(numF)
        cind0 = 1 + numDyn
        for constr in self.pointConstr:
            if constr.lb is not None:
                clb[cind0: cind0 + constr.nf] = constr.lb
            if constr.ub is not None:
                cub[cind0: cind0 + constr.nf] = constr.ub
            cind0 += constr.nf
        # assign to where it should belong to
        self.lb = clb
        self.ub = cub

    def __getTimeGrid__(self, x):
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

    def __callf__(self, x, y):
        """Evaluate those constraints and objective functions."""
        h, useT = self.__getTimeGrid__(x)
        useX, useU, useP = self.__parseX__(x)
        # evaluate objective function
        self.__objModeF__(0, h, useT, useX, useU, useP, x, y)
        # evaluate the system dynamics constraint
        curRow = 1
        curRow = self.__dynconstrModeF__(curRow, h, useT, useX, useU, useP, y)
        # evaluate other constraints
        curRow = self.__constrModeF__(curRow, h, useT, useX, useU, useP, x, y)
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("y is {}".format(y))

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
            for i in range(self.N - 1):
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
        h, useT = self.__getTimeGrid__(x)
        useX, useU, useP = self.__parseX__(x)
        # loop over all the objective functions
        curRow, curNg = self.__objModeG__(0, 0, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        # loop over all system dynamics constraint
        curRow, curNg = self.__dynconstrModeG__(curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
        curRow, curNg = self.__constrModeG__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)

    def __dynconstrModeG__(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics"""
        # TODO: clean up these code so sparse and dense are treated differently
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
                    dynnnz = Jac.nnz - Jac.getcol(0).nnz
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

    def __objModeG__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. G mode. See __constrModeG__ for arguments and output."""
        if not self.objSparseMode:
            curRow, curNg = self.__objModeGUnknown__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        else:
            curRow, curNg = self.__objModeGKnown__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        return curRow, curNg

    def __objModeGKnown__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. G mode. Sparsity structure is known or easy to guess. Use indices to record. 
        See __constrModeG__ for arguments and output."""
        tmpout = np.zeros(1)
        y[0] = 0
        for obj in self.linearObj:  # linear objective function
            y[0] += obj.A.dot(x)
            if needg:
                G[curNg: curNg + obj.nG] = obj.A.data
                if rec:
                    row[curNg: curNg + obj.nG] = curRow
                    col[curNg: curNg + obj.nG] = obj.A.col
                curNg += obj.nG
        for obj in self.linPointObj:  # linear point cost
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            y[0] += obj.A.dot(tmpx)
            if needg:
                G[curNg: curNg + obj.nG] = obj.A.col
                if rec:
                    row[curNg: curNg + obj.nG] = curRow
                    col[curNg: curNg + obj.nG] = self.__patchCol__(obj.index, obj.A.col)
                curNg += obj.nG
        for obj in self.linPathObj:  # linear path cost
            # TODO: caveat, current implementation does not support obj on t_k and it is unnecessary
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                tmpval = obj.A.dot(tmpx)
                y[0] += tmpval * h
                if needg:
                    G[curNg: curNg + obj.nG] = obj.A.data * h
                    if rec:
                        row[curNg: curNg + obj.nG] = curRow
                        col[curNg: curNg + obj.nG] = self.__patchCol__(i, obj.A.col)
                    if self.t0ind > 0:
                        print("Error! Current implementation does not support path objective dependent on time")
                        raise NotImplementedError  # in this case, we will have error
                    if self.tfind > 0:
                        print("Error! Current implementation does not support path objective dependent on time")
                        raise NotImplementedError  # in this case, we will have error
                    curNg += obj.nG
        for obj in self.nonPointObj:  # nonlinear point cost
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            Gpiece = G[curNg: curNg + obj.nG]
            rowpiece = row[curNg: curNg + obj.nG]
            colpiece = col[curNg: curNg + obj.nG]
            obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
            y[0] += tmpout[0]
            if needg:
                if rec:
                    rowpiece[:] = curRow
                    colpiece[:] = self.__patchCol__(obj.index, colpiece)
                curNg += obj.nG
        for obj in self.nonPathObj:  # nonlinear path cost
            # TODO: add support for gradient on tk
            for i in range(self.N - 1):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                Gpiece = G[curNg: curNg + obj.nG]
                rowpiece = row[curNg: curNg + obj.nG]
                colpiece = col[curNg: curNg + obj.nG]
                obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                y[0] += tmpout[0] * h
                if needg:
                    Gpiece *= h
                    if rec:
                        rowpiece[:] = curRow
                        colpiece[:] = self.__patchCol__(i, colpiece)
                    if self.t0ind > 0:
                        print("Error! Current implementation does not support path objective dependent on time")
                        raise NotImplementedError  # in this case, we will have error
                    if self.tfind > 0:
                        print("Error! Current implementation does not support path objective dependent on time")
                        raise NotImplementedError  # in this case, we will have error
                    curNg += obj.nG
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
        return 1, curNg

    def __objModeGUnknown__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. G mode. Sparsity structure not known a prior. Use indices to record.
        We directly accumulate those gradients and assume no violation.
        See __constrModeG__ for arguments and output."""
        tmpout = np.zeros(1)
        y[0] = 0
        tmpG = np.zeros(self.numSol, dtype=float)  # those three arrays temporarily are used for gradient evaluation
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
        # write those guys into G we need
        if needg:
            lenObjInd = len(self.indices)
            G[:lenObjInd] = tmpG[self.indices]
            if rec:
                row[:lenObjInd] = 0
                col[:lenObjInd] = self.indices
            curNg = lenObjInd  # it stores how many G we have collected. Wow, pretty heavy and tedious computation for cost function
        return 1, curNg

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
            constr.__evalg__(tmpx, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
            if rec:
                pieceRow += curRow
                pieceCol[:] = self.__patchCol__(constr.index, pieceCol)
            curRow += constr.nf
            curNg += constr.nG
        for constr in self.pathConstr:
            for i in range(self.N):
                tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                pieceG = G[curNg: curNg + constr.nG]
                pieceRow = row[curNg: curNg + constr.nG]
                pieceCol = col[curNg: curNg + constr.nG]
                constr.__evalg__(tmpx, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
                if rec:
                    pieceRow += curRow
                    pieceCol[:] = self.__patchCol__(i, pieceCol)
                curRow += constr.nf
                curNg += constr.nG
        for constr in self.nonLinConstr:
            pieceG = G[curNg: curNg + constr.nG]
            pieceRow = row[curNg: curNg + constr.nG]
            pieceCol = col[curNg: curNg + constr.nG]
            constr.__evalg__(x, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
            if rec:
                pieceRow += curRow
            curRow += curRow
            curNg += constr.nG
        return curRow, curNg

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
        if lqrobj.tfweight is not None:
            tfweight = lqrobj.tfweight
        else:
            tfweight = 0.0
        self.LQRnG = (self.N - 1) * (useQ + useR + useP) + useF + self.numT

        def __callf__(h, useX, useU, useP):
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
                y += np.sum(lqrobj.P.data * np.sum((useP[:-1, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0)) * h
            y += tfweight * (h * (self.N - 1))
            return y

        def __callg__(h, useX, useU, useP, y, G, row, col, rec, needg):
            """Calculate the lqr cost with gradient information.

            :param h: float, grid size
            :param useX, useU, useP: ndarray, parsed X, U, P from x, the same with __call__F
            :param y: ndarray, a location to write the objective function onto
            :param G, row, col: the gradient information
            :param rec: if we record row and col
            :param needg: if we need gradient information.

            """
            if not needg:
                y[0] = __callf__(h, useX, useU, useP)
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
                    yP = np.sum(lqrobj.P.data * np.sum((useP[:-1, Pcol] - lqrobj.pbase[Pcol]) ** 2, axis=0)) * h
                    G[curG: curG + (self.N - 1) * useP] = 2.0 * h * ((useU[:-1, Pcol] - lqrobj.ubase[Pcol]) * lqrobj.P.data).flatten()
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
