#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
trajOptManifoldCollocationProblem.py

Direct collocation approach with manifold constraint.
We introduce additional variables to solve this issue.
Here the order of a system is really helpful.
For a constraint \phi(q)=0, by order we know up to which we should differentiate the constraints.
"""
from __future__ import division
import numpy as np
from .trajOptCollocationProblem import trajOptCollocProblem
from .libsnopt import probFun
from scipy.sparse import csr_matrix, coo_matrix
from .utility import randomGenInBound, checkInBounds, interp


class manifoldConstr(object):
    """An object which defines state manifold constraints. We might need other ugly hacks.

    Update 1, remove support for constr_order, only allow holonomic constraint, i.e. q
    However, user is allowed to control to any level, position / velocity / acceleration level.
    For different daeOrder, the implementation is slightly different

    """
    def __init__(self, nx, nc, order=2, nG=0, nnzJx=0, nnzJg=0):
        """Constructor for the class.

        Parameters
        ----------
        nx : the system state dimension, for daeSystem it is order * nf
        nc : the manifold constraint dimension, it is not multiplied by order yet
        order : the order of the daeSystem
        constr_order : which order of states (number of time derivative) does the constraint first appear
        nG : non-zero of the output Jacobian
        nnzJx : int, nnz for Jacobian of (J^T \gamma) w.r.t. x
        nnzJg : int, nnz for Jacobian of (J^T \gamma) w.r.t. gamma

        """
        self.nx = nx
        self.nc = nc
        self.nf = (order + 1) * nc
        self.nG = nG
        self.nnzJ_gamma = nnzJg
        self.nnzJ_x = nnzJx

    def __callg__(self, x, F, G, row, col, rec, needg):
        """Constraint evaluation up to order with no correction.

        This function serves as a basic path constraint imposed on every knot point
        It is different from nonLinearPointConstr in the following:
        1. It is autonomous and time is not considered
        2. Only state variables are passed in, no u or p is used

        Parameters
        ----------
        x : ndarray, the state variables including q, dq, ddq, etc
        gamma : ndarray, the correction term
        F : ndarray, it stores the calculation values
        G, row, col : ndarray, the triplet for Jacobian return
        rec, needg : control if G is calculated and row, col are modified

        Returns
        -------
        This subroutine has no return.

        """
        raise NotImplementedError

    def __calc_correction__(self, x, gamma, y, Gq, rowq, colq, Gg, rowg, colg, rec, needg):
        """Calculate the correction term

        It calculate J(qc)^T \gamma, and the sparse Jacobians associated with it

        Parameters
        ----------
        x : array-like, q, dq, ddq, etc
        gamma : array-like, the correction vector gamma
        y : array-like, the function calculation value, it is essentially J(qc)^T \gamma
        Gq, rowq, colq : array-like, store the Jacobian of correction term w.r.t. x
        Gg, rowg, colg : array-like, store the Jacobian of correction term w.r.t. gamma
        rec : if we record row and column information
        needg : if we calculate G

        Returns
        -------
        Calculations are done in-place so no return is required.

        """
        raise NotImplementedError


class trajOptManifoldCollocProblem(trajOptCollocProblem):
    """A class for solving state equality constrained problem.

    The difference between this class with trajOptCollocProblem is

    1. This class directly handles state constrained problem
    2. This class introduce additional variables at collocation points and change deflect constraints
    3. Additional parameters are linearly interpolated previously, now we let it unconstrained. This is in fact the same with adding
    4. How variables are arranged are different
    5. The state constraints should be defined using a special class which returns higher order time derivatives

    The optimization variables are arranged by
    1. Trajectory region, (2*N - 1) * [q, dq, ddq, ..., u, p]
    2. t0 and tf region, it might be empty
    3. addX region, of length lenAddX
    4. gamma region, of length (N-1) * (man_constr.nf)
    5. Auxiliary objective region

    """
    def __init__(self, sys, N, t0, tf, man_constr, addx=None):
        """Constructor for the class, the user is required to provide dynamical system, discretization size, 
        allowable times, state manifold constraint, and possibly auxiliary variables.

        Parameters
        ----------
        sys : the system instance which describes system dynamics
        N : int, discretization grid size
        t0 : float/array-like, allowable initial time
        tf : float/array-like, allowable final time
        man_constr : manifold constraint
        addx : list of addX / one addX / None, additional optimization variables.

        """
        trajOptCollocProblem.__init__(self, sys, N, t0, tf, addx)
        # modify a few parameters due to introduction of man_constr
        man_constr_dim = man_constr.nc
        ext_man_constr_dim = man_constr.nf
        self.numGamma = (N - 1) * man_constr_dim
        assert isinstance(man_constr, manifoldConstr)
        self.man_constr = man_constr
        self.man_constr_dim = man_constr_dim
        self.ext_man_constr_dim = ext_man_constr_dim
        self.numSol += self.numGamma  # add gamma to optimization variables

    def preProcess(self, defect_u=True, defect_p=False, gamma_bound=1):
        """Initialize the problem, allocating spaces.

        Compared with trajOptCollocProblem, it has new sets of variables, different constraints.
        It supports limitation of magnitude of gamma, as gamma_bound means.

        """
        self.colloc_constr_is_on = False  # we never need this, actually
        self.defectU = defect_u
        self.defectP = defect_p
        numDyn = self.dimdyn * self.nPoint  # constraints from system dynamics, they are imposed everywhere
        dynDefectSize = 2 * self.daeOrder * self.dimdyn
        self.dynDefectSize = dynDefectSize
        defectSize = dynDefectSize
        if defect_u:
            defectSize += self.dimu
        if defect_p:
            defectSize += self.dimp
        self.defectSize = defectSize
        numDefectDyn = (self.N - 1) * defectSize  # from enforcing those guys
        self.numDefectDyn = numDefectDyn
        self.numDyn = numDyn + numDefectDyn  # from both nonlinear dynamics and linear defect constraints
        # calculate number of constraint
        numC, nnonlincon, nlincon = self.__sumConstrNum__()
        nnonlincon += self.N * self.man_constr.nf
        # add constraint from manifold constraint
        numC += self.N * self.man_constr.nf

        self.numLinCon = nlincon
        self.numNonLinCon = nnonlincon
        trajOptCollocProblem.__findMaxNG__(self)
        self.numF = 1 + numDyn + numDefectDyn + numC

        # analyze all objective functions in order to detect pattern for A, and additional variables for other nonlinear objective function
        spA, addn = self.__analyzeObj__(self.numSol, self.numF)
        self.objaddn = addn  # this is important for multiple objective function support
        self.numSol += addn
        self.numF += addn
        probFun.__init__(self, self.numSol, self.numF)  # not providing G means we use finite-difference
        # we are ready to write Aval, Arow, Acol for this problem. They are arranged right after dynamics
        self.__setAPattern__(numDyn, nnonlincon, spA)
        self.__setXbound__(gamma_bound)
        self.__setFbound__()
        # detect gradient information
        randX = self.randomGenX()
        self.__turnOnGrad__(randX)


    def __getSparsity__(self, x0):
        """Detect the sparsity pattern with an initial guess."""
        trajOptCollocProblem.__getSparsity__(self, x0)
        # add nnz Jacobian from defect constraint, i.e. J(q)^T \gamma
        self.numG += (self.N - 1) * (self.man_constr.nnzJ_gamma + self.man_constr.nnzJ_x)
        # add nnz Jacobian from manifold constraint on knot points
        self.numG += self.N * self.man_constr.nG
        self.nG = self.numG

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix.

        :param x: ndarray, the solution to the problem
        :param y: ndarray, return F
        :param G, row, col: ndarray, information of gradient
        :param rec, needg: if we record/ if we need gradient

        """
        y[0] = 0  # since this row is purely linear
        h, useT = self.__get_time_grid__(x)
        useX, useU, useP = self.__parseX__(x)
        useGamma = self.__parseGamma__(x)
        # loop over all system dynamics constraint
        curRow = 1
        curNg = 0
        curRow, curNg = self.__dynconstr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, useGamma, y, G, row, col, rec, needg)
        curRow, curNg = self.__constr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        curRow, curNg = self.__manifold_constr_mode_g__(curRow, curNg, useX, y, G, row, col, rec, needg)
        curRow += self.numLinCon
        # loop over all the objective functions, I haven't checked if order is correct since some linear constraints are followed
        curRow, curNg = self.__obj_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        print(y[281:287])
        pass

    def __manifold_constr_mode_g__(self, curRow, curNg, useX, y, G, row, col, rec, needg):
        """Calculate the manifold constraints on knots with autonomous assumption."""
        for i in range(self.N):
            next_Ng = curNg + self.man_constr.nG
            next_Row = curRow + self.ext_man_constr_dim
            pieceG = G[curNg: next_Ng]
            pieceRow = row[curNg: next_Ng]
            pieceCol = col[curNg: next_Ng]
            self.man_constr.__callg__(useX[2*i], y[curRow: next_Row], pieceG, pieceRow, pieceCol, rec, needg)
            if rec:
                pieceRow += curRow
                pieceCol += self.getStateIndexByIndex(2*i)
            curRow = next_Row
            curNg = next_Ng
        return curRow, curNg

    def __dynconstr_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, useGamma, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics. Code are copied and modified.
        We implement this by call the trajOptCollocProblem version and append some variables
        """
        nPoint = self.nPoint
        dimdyn = self.dimdyn
        dimq = self.dimq
        if self.daeOrder == 2:
            cur_row = curRow + nPoint * dimdyn + dimq  # We add additional dimdyn since it where J^T\gamma comes into play
            row_step = self.dynDefectSize
        else:
            cur_row = curRow
            row_step = dimdyn
        # call ordinary dyn constr
        curRow, curNg = trajOptCollocProblem.__dynconstr_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
        # loop over defect constraints
        for i in range(self.N - 1):
            # this is done by calling the __calc_correction__ function implemented by the user
            # the only difference is it calculates two sets of sparse Jacobian, and it is autonomous
            # for second-order system, there is no problem, but for first order one, we implicitly assume the first dimq are
            Gx = G[curNg: curNg + self.man_constr.nnzJ_x]
            rowx = row[curNg: curNg + self.man_constr.nnzJ_x]
            colx = col[curNg: curNg + self.man_constr.nnzJ_x]
            curNg += self.man_constr.nnzJ_x
            Gg = G[curNg: curNg + self.man_constr.nnzJ_gamma]
            rowg = row[curNg: curNg + self.man_constr.nnzJ_gamma]
            colg = col[curNg: curNg + self.man_constr.nnzJ_gamma]
            curNg += self.man_constr.nnzJ_gamma
            self.man_constr.__calc_correction__(useX[2*i + 1, :dimq], useGamma[i], y[cur_row: cur_row+dimq],
                                                Gx, rowx, colx, Gg, rowg, colg, rec, needg)
            if rec:
                rowx += cur_row
                rowg += cur_row
                colg += self.getGammaIndexByIndex(i)
                colx += self.getStateIndexByIndex(2*i + 1)
            cur_row += row_step
        return curRow, curNg

    def parseF(self, guess):
        """Give an guess, evaluate it and parse into parts.

        :param guess: ndarray, (numSol, ) a guess or a solution to check
        :returns: dict, containing objective and parsed constraints

        """
        assert len(guess) == self.numSol
        N = self.N
        nPoint = self.nPoint
        dimx = self.dimx
        dimdyn = self.dimdyn
        nc = self.man_constr.nc
        nf = self.man_constr.nf
        y = np.zeros(self.numF)

        # call previous function
        rst = trajOptCollocProblem.parseF(self, guess, y)

        numC, nnonlincon, nlincon = self.__sumConstrNum__()
        curN = 1 + (2 * N - 1) * dimdyn + self.numDefectDyn + nnonlincon

        manCon = np.reshape(y[curN: curN + N*nf], (N, nf))
        rst['man'] = manCon
        # we ignore linear and auxVc constraints

        useGamma = self.__parseGamma__(guess)
        rst['gamma'] = useGamma

        return rst

    def __parseGamma__(self, x):
        """Return gamma as a N - 1 by correct size np array"""
        n1 = self.numTraj + self.lenAddX
        return np.reshape(x[n1: n1 + self.numGamma], (self.N - 1, -1))

    def getGammaIndexByIndex(self, i):
        """Return the index of gamma by its index."""
        return self.numTraj + self.lenAddX + i * self.man_constr_dim

    def __setXbound__(self, gamma_bound):
        """Set bounds on decision variables."""
        trajOptCollocProblem.__setXbound__(self)
        # we only need to set gamma which are unbounded
        curN = self.numTraj + self.lenAddX
        # gamma cannot be too large, I guess
        self.batchSetXlb(-gamma_bound*np.ones(self.numGamma), curN)
        self.batchSetXub(gamma_bound*np.ones(self.numGamma), curN)

    def __setFbound__(self):
        """Set bound on F"""
        # set bound on F
        numF = self.numF
        numDyn = self.numDyn
        clb = np.zeros(numF)
        cub = np.zeros(numF)
        clb[0] = -1e20
        cub[0] = 1e20
        cind0 = 1 + numDyn
        cind0 = self._setNonLinConstr(clb, cub, cind0)
        cind0 = self._setManifoldConstr(clb, cub, cind0)
        cind0 = self._setLinConstr(clb, cub, cind0)
        # assign to where it should belong to
        self.lb = clb
        self.ub = cub

    def _setManifoldConstr(self, clb, cub, cind0):
        """Set the manifold constraint."""
        cindf = cind0 + self.N * self.man_constr.nf
        clb[cind0: cindf] = 0
        cub[cind0: cindf] = 0
        return cindf
