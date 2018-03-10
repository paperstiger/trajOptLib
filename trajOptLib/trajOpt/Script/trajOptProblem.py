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
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
from trajOptBase import *
from libsnopt import snoptConfig, probFun, solver
from utility import parseX


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class trajOptProblem(probFun):
    """A class for definition of trajectory optimization problem."""
    def __init__(self, sys, N=20, t0=0.0, tf=1.0):
        """Initialize problem by system, discretization grid size, and allowable time
        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        """
        assert isinstance(sys, system)
        self.sys = sys
        self.N = 20
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
        numSol = numX + numU + numP + numT
        numObjG = numSol  # G from objective, assume dense
        # summarize number of pure linear constraints
        # TODO: implement linear objective and constraint
        numDyn = self.dimx * (self.N - 1)  # constraints from system dynamics
        # TODO: support sparsity in sys dyn constraints
        numDynG = (self.N - 1) * self.dimx * (self.dimpoint + numT)  # G from dyn
        numC = 0
        numCG = 0  # G from C
        # I only care about those in numC
        for constr in self.pointConstr:
            cind0 += constr.nf
            numC += constr.nf
            numG += constr.nG
        numF = 1 + numDyn + numC
        numG = numObjG + numDynG + numCG
        # TODO: implement gradient
        # currently try a case with only finite difference gradient
        super(trajOptProblem, self).__init__(numSol, numF)  # not providing G means we use finite-difference
        # create bound on x
        xlb = np.zeros(numSol)
        xub = np.zeros(numSol)
        if self.xbd[0] is not None:
            stateLb = np.reshape(xlb[:numX], (self.N, self.dimx))
            for i in range(self.dimx):
                stateLb[:, i] = self.xbd[0][i]
            # set lb for x0 and xf
            if self.x0bd[0] is not None:
                stateLb[0] = self.x0bd[0]
            if self.xfbd[0] is not None:
                stateLb[self.N - 1] = self.xfbd[0]
        if self.xbd[1] is not None:
            stateUb = np.reshape(xub[:numX], (self.N, self.dimx))
            for i in range(self.dimx):
                stateUb[:, i] = self.xbd[1][i]
            # set ub for x0 and xf
            if self.x0bd[1] is not None:
                stateUb[0] = self.x0bd[1]
            if self.xfbd[1] is not None:
                stateUb[self.N - 1] = self.xfbd[1]
        if self.ubd[0] is not None:
            ctrlLb = np.reshape(xlb[numX:numX + numU], (self.N, self.dimu))
            for i in range(self.dimu):
                ctrlLb[:, i] = self.ubd[0][i]
        if self.ubd[1] is not None:
            ctrlUb = np.reshape(xub[numX:numX + numU], (self.N, self.dimu))
            for i in range(self.dimu):
                ctrlUb[:, i] = self.ubd[1][i]
        if self.pbd[0] is not None and self.dimp > 0:
            pLb = np.reshape(xlb[numX+numU:numX+numU+numP], (self.N, self.dimp))
            for i in range(self.dimp):
                pLb[:, i] = self.pbd[0][i]
        if self.pbd[1] is not None and self.dimp > 0:
            pUb = np.reshape(xub[numX+numU:numX+numU+numP], (self.N, self.dimp))
            for i in range(self.dimp):
                pUb[:, i] = self.pbd[1][i]
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
        # set bound on F
        clb = np.zeros(numF)
        cub = np.zeros(numF)
        cind0 = 1 + numDyn
        for constr in self.pointConstr:
            if constr.lb is not None:
                clb[cind0: cind0 + constr.nf] = constr.lb
            if constr.ub is not None:
                cub[cind0: cind0 + constr.nf] = constr.ub
            cind0 += constr.nf
        # assign back to class
        self.lb = clb
        self.ub = cub
        self.xlb = xlb
        self.xub = xub
        self.numX = numX
        self.numU = numU
        self.numP = numP
        self.numT = numT

    def __callf__(self, x, y):
        """Evaluate those constraints and objective functions."""
        if self.fixt0:
            uset0 = self.t0
            if self.fixtf:
                usetf = self.tf
            else:
                usetf = x[self.numX + self.numU + self.numP]
        else:
            uset0 = x[self.numX + self.numU + self.numP]
            if self.fixtf:
                usetf = self.tf
            else:
                usetf = x[self.numX + self.numU + self.numP + 1]
        h = (usetf - uset0) / (self.N - 1)
        numX = self.numX
        numU = self.numU
        numP = self.numP
        useX = np.reshape(x[:numX], (self.N, self.dimx))
        useU = np.reshape(x[numX: numX + numU], (self.N, self.dimu))
        useP = np.reshape(x[numX + numU: numX + numU + numP], (self.N, self.dimp))
        useT = np.linspace(uset0, usetf, self.N)
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

    def addNonLinObj(self, nonLinObj):
        """Add nonlinear objective function.
        :param: nonLinObj: a nonLinObj class
        """
        assert isinstance(nonLinObj, nonLinObj)
        self.nonLinObj.append(nonLinObj)

    def addNonPointObj(self, nonPntObj, path=False):
        """Add nonlinear point objective.
        :param nonPntObj: nonLinObj class
        :param path: bool, if this obj
        """
        assert isinstance(nonPntObj, pointObj)
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


if __name__ == '__main__':
    main()
