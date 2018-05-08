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
    """An object which defines state manifold constraints. We might need other ugly hacks."""
    def __init__(self, nx, nc, order, constr_order=0, nG=0, nnzJx=0, nnzJg=0):
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
        self.order = order
        self.constr_order = constr_order
        self.nf = (order - constr_order + 1) * nc
        self.nG = nG
        self.nnzJ_gamma = nnzJg
        self.nnzJ_x = nnzJx

    def __callg__(self, x, F, G, row, col, rec, needg):
        """Constraint evaluation to order with no correction.

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
    5. auxiliary vc region, this is used for nasty conflict issue (N-1)*dimdyn # TODO: this needs verification and generalization
    4. Auxiliary objective region

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
        self.numAuxVc = (N - 1) * self.dimdyn
        assert isinstance(man_constr, manifoldConstr)
        self.man_constr = man_constr
        self.man_constr_dim = man_constr_dim
        self.ext_man_constr_dim = ext_man_constr_dim
        self.numSol += self.numGamma  # add those values to optimization variables
        self.numSol += self.numAuxVc  # add those values for inconvenient un-repeat issue between linear and nonlinear Jacobian

    def preProcess(self, defect_u=True, defect_p=False):
        """Initialize the problem, allocating spaces.

        Compared iwth trajOptCollocProblem, it has new sets of variables, different constraints.

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
        # add constraints from auxiliary velocities
        nauxvccon = self.numAuxVc
        numC += nauxvccon
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
        self.__setAPattern__(numDyn, nnonlincon, nlincon, spA)
        self.__setXbound__()
        self.__setFbound__()
        # detect gradient information
        randX = self.randomGenX()
        self.__turnOnGrad__(randX)

    def __setAPattern__(self, ndyncon, nnonlincon, nlincon, spA):
        """Set sparsity pattern from linear constraints and objective functions.

        see trajOptCollocProblem for details

        This overrides previous one by considering auxiliary Vc variables

        """
        curRow, A, row, col = self.__setDefectPattern__(ndyncon, self.defectU, self.defectP)
        curRow += nnonlincon
        # we are ready to parse linear constraints
        lstCA, lstCArow, lstCAcol = self.__parseLinearConstraints__(curRow)

        # append the auxiliary vc constraints right after linear constraints
        curRow += nlincon
        new_A = np.zeros((self.N - 1) * self.dimdyn * 2)
        new_A_row = np.zeros((self.N - 1) * self.dimdyn * 2, dtype=int)
        new_A_col = np.zeros((self.N - 1) * self.dimdyn * 2, dtype=int)
        accum_n = 0
        baseRange = np.arange(self.dimdyn)
        # set A for auxiliary velocity variables
        for i in range(self.N - 1):
            new_A[accum_n: accum_n + self.dimdyn] = 1.0
            new_A[accum_n + self.dimdyn: accum_n + 2*self.dimdyn] = -1.0
            new_A_row[accum_n: accum_n + self.dimdyn] = curRow + baseRange
            new_A_row[accum_n + self.dimdyn: accum_n + 2*self.dimdyn] = curRow + baseRange
            new_A_col[accum_n: accum_n + self.dimdyn] = self.getStateIndexByIndex(2*i + 1) + self.dimdyn + baseRange  # TODO: assume vc is dq term
            new_A_col[accum_n + self.dimdyn: accum_n + 2*self.dimdyn] = self.getAuxVcIndexByIndex(i) + baseRange  # TODO: assume vc has length dimdyn
            curRow += self.dimdyn
            accum_n += 2*self.dimdyn
        lstCA.append(new_A)
        lstCArow.append(new_A_row)
        lstCAcol.append(new_A_col)

        # concatenate all those things together
        lstCA.append(spA.data)
        lstCA.append(A)
        lstCArow.append(spA.row)
        lstCArow.append(row)
        lstCAcol.append(spA.col)
        lstCAcol.append(col)

        self.Aval = np.concatenate(lstCA)
        self.Arow = np.concatenate(lstCArow)
        self.Acol = np.concatenate(lstCAcol)
        self.spA = csr_matrix((self.Aval, (self.Arow, self.Acol)), shape=(self.nf, self.nx))
        self.spA_coo = self.spA.tocoo()


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
        useAuxVc = self.__parseAuxVc__(x)
        # loop over all system dynamics constraint
        curRow = 1
        curNg = 0
        curRow, curNg = self.__dynconstr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, useGamma, useAuxVc, y, G, row, col, rec, needg)
        curRow, curNg = self.__constr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        curRow, curNg = self.__manifold_constr_mode_g__(curRow, curNg, useX, y, G, row, col, rec, needg)
        curRow += self.numLinCon
        curRow += self.numAuxVc
        # loop over all the objective functions, I haven't checked if order is correct since some linear constraints are followed
        curRow, curNg = self.__obj_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
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

    def __dynconstr_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, useGamma, useAux, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics. Code are copied and modified.
        We implement this by call the trajOptCollocProblem version and append some variables
        """
        nPoint = self.nPoint
        dimdyn = self.dimdyn
        cur_row = curRow + nPoint * dimdyn + dimdyn*(2*self.man_constr.constr_order+1)  # from the 2*N - 1 dimdyn constraints from dynamics
        # explanation of why we use 2*constr_order + 1:
        # constr_order is for in case the constraint is based on control so we have to correct for acceleration
        # the 1 is because we are correcting for time derivative part, and in our (q, dq) pair q is leading
        curRow, curNg = trajOptCollocProblem.__dynconstr_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
        # loop over defect constraints
        for i in range(self.N - 1):
            # this is done by calling the __calc_correction__ function implemented by the user
            # the only difference is it calculates two sets of sparse Jacobian, and it is autonomous
            Gx = G[curNg: curNg + self.man_constr.nnzJ_x]
            rowx = row[curNg: curNg + self.man_constr.nnzJ_x]
            colx = col[curNg: curNg + self.man_constr.nnzJ_x]
            curNg += self.man_constr.nnzJ_x
            Gg = G[curNg: curNg + self.man_constr.nnzJ_gamma]
            rowg = row[curNg: curNg + self.man_constr.nnzJ_gamma]
            colg = col[curNg: curNg + self.man_constr.nnzJ_gamma]
            curNg += self.man_constr.nnzJ_gamma
            tmpX = useX[2*i + 1].copy()
            tmpX[dimdyn: 2*dimdyn] = useAux[i]
            self.man_constr.__calc_correction__(tmpX, useGamma[i], y[cur_row: cur_row+dimdyn],
                                                Gx, rowx, colx, Gg, rowg, colg, rec, needg)
            y[cur_row: cur_row + dimdyn] *= -1
            if needg:
                Gx[:] *= -1
                Gg[:] *= -1
            if rec:
                rowx += cur_row
                rowg += cur_row
                colg += self.getGammaIndexByIndex(i)
                # create a mask for colx filter out which are assoicated with vc, move it to aux
                mask = (colx >= dimdyn) & (colx < 2*dimdyn)
                colx[~mask] += self.getStateIndexByIndex(2*i + 1)
                colx[mask] += self.getAuxVcIndexByIndex(i) - dimdyn  # since this part only has vel, we have to subtract q part
            cur_row += self.dynDefectSize
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
        useAuxVc = self.__parseAuxVc__(guess)
        rst['gamma'] = useGamma
        rst['useVc'] = useAuxVc

        return rst

    def __parseGamma__(self, x):
        """Return gamma as a N - 1 by correct size np array"""
        n1 = self.numTraj + self.lenAddX
        return np.reshape(x[n1: n1 + self.numGamma], (self.N - 1, -1))

    def __parseAuxVc__(self, x):
        n1 = self.numTraj + self.lenAddX + self.numGamma
        return np.reshape(x[n1: n1+self.numAuxVc], (self.N - 1, -1))

    def getGammaIndexByIndex(self, i):
        """Return the index of gamma by its index."""
        return self.numTraj + self.lenAddX + i * self.man_constr_dim

    def getAuxVcIndexByIndex(self, i):
        """Return the index of gamma by its index."""
        return self.numTraj + self.lenAddX + self.numGamma + i * self.dimdyn

    def __setXbound__(self):
        """Set bounds on decision variables."""
        trajOptCollocProblem.__setXbound__(self)
        # we only need to set gamma which are unbounded
        curN = self.numTraj + self.lenAddX
        # gamma cannot be too large, I guess
        self.batchSetXlb(-1e20*np.ones(self.numGamma), curN)
        self.batchSetXub(1e20*np.ones(self.numGamma), curN)
        curN += self.numGamma
        finalN = self.numSol
        self.batchSetXlb(-1e20*np.ones(finalN - curN), curN)
        self.batchSetXub(1e20*np.ones(finalN - curN), curN)

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
        # the bounds for objaddn and auxVcCon is 0 so we are good so far
        # assign to where it should belong to
        self.lb = clb
        self.ub = cub

    def _setManifoldConstr(self, clb, cub, cind0):
        """Set the manifold constraint."""
        cindf = cind0 + self.N * self.man_constr.nf
        clb[cind0: cindf] = 0
        cub[cind0: cindf] = 0
        return cindf
