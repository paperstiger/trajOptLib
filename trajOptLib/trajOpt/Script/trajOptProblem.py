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
from trajOptBase import system, linearObj, linearPointObj, nonLinObj, nonPointObj, pointConstr
from libsnopt import snoptConfig, probFun, solver
from utility import parseX
from scipy import sparse
from scipy.sparse import spmatrix, coo_matrix


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class trajOptProblem(probFun):
    """A class for definition of trajectory optimization problem."""
    def __init__(self, sys, N=20, t0=0.0, tf=1.0, gradmode=True):
        """Initialize problem by system, discretization grid size, and allowable time
        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
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
        self.dimx = sys.nx
        self.dimu = sys.nu
        self.dimp = sys.np
        self.dimpoint = sys.nx + sys.nu + sys.np
        self.ubd = [None, None]
        self.xbd = [None, None]
        self.pbd = [None, None]
        self.x0bd = [None, None]
        self.xfbd = [None, None]
        self.linearObj = []
        self.linPointObj = []
        self.linPathObj = []
        self.nonLinObj = []
        self.nonPointObj = []
        self.nonPathObj = []
        self.pointConstr = []
        self.pathConstr = []
        self.gradmode = gradmode  # this controls if __callg__ will be called
        self.gradOn = False  # indicates if we have set up all the gradient information

    def preProcess(self):
        """Initialize the instances of probFun now we are ready."""
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
        numDyn = self.dimx * (self.N - 1)  # constraints from system dynamics
        numC = 0
        for constr in self.pointConstr:
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
        self.gradOn = True
        self.grad = True
        self.getSparsity(x0)

    def getSparsity(self, x0):
        self.indices = self.getObjSparsity(x0)
        numObjG = len(self.indices)  # G from objective, assume dense
        # summarize number of pure linear constraints
        # numDynG = (self.N - 1) * self.dimx * (self.dimpoint + numT)  # G from dyn
        numDynG = self.getDynSparsity(x0)
        numCG = 0  # G from C
        # I only care about those in numC
        for constr in self.pointConstr:
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
        :param x: ndarray, the guess/sol
        :rtype indices: ndarray, the indices with non-zero objective gradient
        """
        h, useT = self.__getTimeGrid__(x)
        useX, useU, useP = self.__parseX__(x)
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
        # summarize # nnz for obj
        nnzObj = np.sum(mask)
        logger.debug("#nnz from obj is %d" % nnzObj)
        indices = np.where(mask)[0]  # this means what indices we should choose from
        return indices

    def __assignTo__(self, index, nnz, value, target, selfadd):
        """Assign a block of value to a target vector. This helps us finding sparsity/assigning values.
        :param index: int, index of a point constraint or -1 indicates it is overall
        :param nnz: indices of nnz elements. For point constraint only, it is within xdim
        :param value: scalar / ndarray, the values we want to assign. If ndarray, its length has to equal nnz
        :param target: ndarray, the target array to assign value to
        :param selfadd: bool, indicates if we eadd value to target, otherwise we assign
        """
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
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
        :rtype: h: float, grid size
        :rtype: useT: the grid being used
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
        # loop over all the objective functions
        tmpout = np.zeros(1)
        y[0] = 0
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
        # loop over all system dynamics constraint
        cDyn = np.reshape(y[1:1+(self.N - 1) * self.dimx], (self.N - 1, self.dimx))
        for i in range(self.N - 1):
            # evaluate gradient of system dynamics
            ydot = self.sys.dyn(useT[i], useX[i], useU[i], useP[i])
            cDyn[i] = useX[i] + h * ydot - useX[i + 1]
        # loop over other constraints
        cInd = 1 + cDyn
        for constr in self.pointConstr:
            tmpx = np.concatenate(([useT[constr.index]], useX[constr.index], useU[constr.index], useP[constr.index]))
            constr.__evalf__(tmpx, y[cInd: cInd + constr.nf])

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix"""
        if not self.gradOn:
            self.__turnOnGrad__(x)
        h, useT = self.__getTimeGrid__(x)
        useX, useU, useP = self.__parseX__(x)
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        numX, numU, numP, numT = self.numX, self.numU, self.numP, self.numT
        # loop over all the objective functions
        tmpout = np.zeros(1)
        y[0] = 0
        tmpG = np.zeros(self.numSol, dtype=float)  # those three arrays temporarily are used for gradient evaluation
        tmpcallG = np.zeros(self.numSol, dtype=float)  # those three arrays temporarily are used for gradient evaluation
        tmpcallrow = np.zeros(self.numSol, dtype=int)
        tmpcallcol = np.zeros(self.numSol, dtype=int)
        for obj in self.linearObj:
            y[0] += obj.A.dot(x)
            if needg:
                tmpG[obj.A.col] += obj.A.data
        for obj in self.linPointObj:
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            y[0] += obj.A.dot(tmpx)
            if needg:
                self.__assignTo__(obj.index, obj.A.col, obj.A.data, tmpG, True)
        for obj in self.linPathObj:  # TODO: caveat, current implementation does not support obj on t_k and it is unnecessary
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
        for obj in self.nonLinObj:
            obj.__callg__(x, tmpout, tmpcallG, tmpcallrow, tmpcallcol, True, needg)
            y[0] += tmpout[0]
            if needg:
                self.__assignTo__(-1, tmpcallcol[:obj.nG], tmpcallG[:obj.nG], tmpG, True)
        for obj in self.nonPointObj:
            tmpx = np.concatenate(([useT[obj.index]], useX[obj.index], useU[obj.index], useP[obj.index]))
            obj.__callg__(tmpx, tmpout, tmpcallG, tmpcallrow, tmpcallcol, True, needg)
            y[0] += tmpout[0]
            if needg:
                self.__assignTo__(obj.index, tmpcallcol[:obj.nG], tmpcallG[:obj.nG], tmpG, True)
        for obj in self.nonPathObj:  # TODO: add support for gradient on tk
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
        # summarize the sparsity of objective function
        if needg:
            lenObjInd = len(self.indices)
            G[:lenObjInd] = tmpG[self.indices]
            if rec:
                row[:lenObjInd] = 0
                col[:lenObjInd] = self.indices
            curNg = lenObjInd  # it stores how many G we have collected. Wow, pretty heavy and tedious computation for cost function
        curRow = 1  # after objective
        # loop over all system dynamics constraint
        cDyn = np.reshape(y[1:1+(self.N - 1) * self.dimx], (self.N - 1, self.dimx))
        for i in range(self.N - 1):
            # evaluate gradient of system dynamics TODO: support other types of integration scheme
            ydot, Jac = self.sys.Jdyn(useT[i], useX[i], useU[i], useP[i])
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
        # loop over other constraints
        curRow = 1 + cDyn.shape[0]
        for constr in self.pointConstr:
            tmpx = np.concatenate(([useT[constr.index]], useX[constr.index], useU[constr.index], useP[constr.index]))
            pieceG = G[curNg: curNg + constr.nG]
            pieceRow = row[curNg: curNg + constr.nG]
            pieceCol = col[curNg: curNg + constr.nG]
            constr.__evalg__(tmpx, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
            if rec:
                pieceRow += curRow
            curRow += constr.nf
            curNg += constr.nG

    def __patchCol__(self, index, col):
        """Find which indices it belongs to the original one for a local matrix at col"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        newcol = col.copy()
        maskx = col < 1 + dimx
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
        :param: nonLinObj: a nonLinObj class
        """
        assert isinstance(nonlinObj, nonLinObj)
        self.nonLinObj.append(nonlinObj)

    def addNonPointObj(self, nonPntObj, path=False):
        """Add nonlinear point objective.
        :param nonPntObj: nonLinObj class
        :param path: bool, if this obj
        """
        assert isinstance(nonPntObj, nonPointObj)
        if path:
            self.nonPathObj.append(nonPntObj)
        else:
            self.nonPointObj.append(nonPntObj)

    def addPointConstr(self, pntConstr):
        """Add point constraint.
        :param: pntConstr: pointConstr class
        """
        assert isinstance(pntConstr, pointConstr)
        self.pointConstr.append(pntConstr)

    def setN(self, N):
        """Set N"""
        self.N = N

    def sett0tf(self, t0, tf):
        """Set t0 and tf"""
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
