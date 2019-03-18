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
from __future__ import division
import numpy as np
from .trajOptBase import LinearObj as linearObj, LinearPointObj as linearPointObj
from .trajOptBase import LinearPointConstr as linearPointConstr, LinearConstr as linearConstr
from .trajOptBase import NonLinearPointObj as nonLinearPointObj, NonLinearObj as nonLinearObj
from .trajOptBase import NonLinearPointConstr as nonLinearPointConstr, NonLinearConstr as nonLinearConstr
from .trajOptBase import QuadPenalty
from .trajOptBase import LqrObj as lqrObj
from .trajOptBase import AddX as addX
from .trajOptBase import DaeSystem as daeSystem
from pyoptsolver import OptProblem
from .utility import random_gen_in_bound as randomGenInBound, check_in_bounds as checkInBounds, interp
from scipy.sparse import coo_matrix, csr_matrix


class TrajOptCollocProblem(OptProblem):
    """A class for definition of trajectory optimization problem using collocation constraints.

    A general framework for using this class is to:

    1. Define a class which implements a DAE system, f(t, x, p, u)=0. This is a typical case for humanoid problem, but also flexible enough for simpler first order system.
    2. Optionally, write desired cost function by subclass/creating from the list of available cost functions. This can be built incrementally by adding pieces.
    3. Optionally, write desired constraint functions by subclass from available constraints. This can be built incrementally by adding pieces. Our approach can correctly detect the Jacobian structure.
    4. Create this class with selected system, discretization, t0, tf range
    5. Set bounds for state, control, parameters, x0 and xf
    6. Add objective functions and constraints to this class
    7. Call preProcess method explicitly
    8. Create snoptConfig instance and choose desired options
    9. Construct the solver
    10. Use the solver to solve with either automatic guess or user provided guess

    The system dynamics constraints are imposed using direct collocation approach with some additional optimization
    variables as suggested by others.

    :currentmodule:
    .. exclude-members:: addLinearPointConstr
    .. exclude-members:: addLinearConstr
    .. exclude-members:: addNonLinearConstr
    .. exclude-members:: addNonLinearPointConstr
    .. exclude-members:: addNonLinearPointObj
    .. exclude-members:: addNonLinearObj
    .. exclude-members:: addLinearPointObj
    .. exclude-members:: addLinearObj
    .. exclude-members:: getAddXIndexByIndex
    .. exclude-members:: getStateIndexByIndex
    .. exclude-members:: getControlIndexByIndex
    .. exclude-members:: getParamIndexByIndex

    """
    def __init__(self, sys, N, t0, tf, addx=None):
        """Initialize problem by system, discretization grid size, and allowable time

        Change history: now I remove gradmode option since I require the gradient be provided analytically all the time.
        I remove unnecessary linear objective functions.
        I reorder variable so q, dq, ddq, u, p are at consecutive place.

        :param sys: system, describe system dynamics
        :param N: int, discretization grid size, a uniform grid
        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :param addX: list of addX / one addX / None, additional optimization variables.

        """
        assert isinstance(sys, daeSystem)
        self.sys = sys
        self.N = N
        self.nPoint = 2 * self.N - 1
        numT = self._handleTime(t0, tf)
        self.dimx = sys.nx
        self.dimdyn = sys.nf  # dimension of dynamics constraint
        self.dimu = sys.nu
        self.dimp = sys.np
        self.dimpoint = sys.nx + sys.nu + sys.np  # each point have those variables
        self.daeOrder = sys.order
        if self.daeOrder == 1:
            self.dimq = self.dimdyn // 2
        elif self.daeOrder == 2:
            self.dimq = self.dimdyn
        self.ubd = [None, None]
        self.xbd = [None, None]
        self.pbd = [None, None]
        self.x0bd = [None, None]
        self.xfbd = [None, None]
        # lqr cost function
        self.lqrObj = None
        self.LQRnG = 0
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
        self.pathConstrIndexPairs = []  # this records the indexes that path constraints are imposed
        self.nonLinConstr = []  # stores general nonlinear constraint
        self.linPointConstr = []
        self.linPathConstr = []
        self.linearConstr = []
        # calculate number of variables to be optimized, time are always the last
        numX = self.nPoint * self.dimx
        numU = self.nPoint * self.dimu
        numP = self.nPoint * self.dimp
        if addx is None:
            self.lenAddX = 0
        else:
            if not isinstance(addx, list):
                addx = [addx]
            for tmp in addx:
                assert isinstance(tmp, addX)
            self.addX = addx
            self.lenAddX = sum([tmp.n for tmp in addx])
        numSol = numX + numU + numP + numT
        self.numX = numX
        self.numU = numU
        self.numP = numP
        self.numT = numT
        self.numDynVar = numX + numU + numP  # numDynVar includes q, dq, ddq, u, p
        self.numTraj = numSol  # make it clear, numTraj contains numDynVar + time
        self.numSol = numSol + self.lenAddX
        self.t0ind, self.tfind = self.__getTimeIndices()
        self.colloc_constr_is_on = False

    def pre_process(self, colloc_constr_is_on=False, defect_u=True, defect_p=True):
        """Initialize the instances of probFun now we are ready.

        Call this function after the objectives and constraints have been set appropriately.
        It calculate the space required for SNOPT and allocates sparsity structure if necessary.

        :param colloc_constr_is_on: bool, if we also impose constraints on those collocation points.
        :param defect_u: bool, if we want to impose defect constraint on u, i.e. umid=(u0+uf)/2
        :param defect_p: bool, if we want to impose defect constraint on p, i.e. pmid=(p0+pf)/2

        **Caveat** it might make problem over-constrained, if the path constraints are equality constraints.

        """
        self.defectU = defect_u
        self.defectP = defect_p
        self.colloc_constr_is_on = colloc_constr_is_on
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
        # constrain t0 and tf
        self._set_t0_tf_constr()

        numC, nnonlincon, nlincon = self.__sumConstrNum__()

        self.numLinCon = nlincon
        self.numNonLinCon = nnonlincon
        self.__findMaxNG__()
        self.numF = 1 + numDyn + numDefectDyn + numC

        # analyze all objective functions in order to detect pattern for A, and additional variables for other nonlinear objective function
        spA, addn = self.__analyzeObj__(self.numSol, self.numF)
        self.objaddn = addn  # this is important for multiple objective function support
        self.numSol += addn
        self.numF += addn
        OptProblem.__init__(self, self.numSol, self.numF)  # not providing G means we use finite-difference
        # we are ready to write Aval, Arow, Acol for this problem. They are arranged right after dynamics
        self.__setAPattern__(numDyn, nnonlincon, spA)
        self.__setXbound__()
        self.__setFbound__()
        # detect gradient information
        randX = self.randomGenX()
        self.__turnOnGrad__(randX)

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

    def _handleTime(self, t0, tf):
        """Deal with time settings.

        If t0 and tf are both free, we want tf to be greater than t0

        :param t0: float/array like, allowable t0
        :param tf: float/array like, allowable tf
        :return: int, number of variable associated with time
        """
        self.tf = tf
        self.t0 = t0
        numT = 2
        if np.isscalar(tf):
            self.fixtf = True
            numT -= 1
        else:
            self.fixtf = False
            assert tf[0] <= tf[1]
        if np.isscalar(t0):
            self.fixt0 = True
            numT -= 1
        else:
            self.fixt0 = False
            assert t0[0] <= t0[1]
        if self.fixt0 and self.fixtf:
            self.fixTimeMode = True
        else:
            self.fixTimeMode = False
        return numT

    def _set_t0_tf_constr(self):
        """Based on occasions, we set constraints on time."""
        if self.fixt0:
            if not self.fixtf:
                self.tf[0] = max(self.tf[0], self.t0 + 1e-6)
        else:
            if self.fixtf:
                self.t0[1] = min(self.t0[1], self.tf - 1e-6)
            else:
                if self.t0[1] > self.tf[0]:  # we might have trouble
                    a = np.zeros((1, self.numSol))
                    a[0, self.t0ind] = 1
                    a[0, self.tfind] = -1
                    self.addConstr(linearConstr(a, -1e20*np.ones(1), 1e-6*np.ones(1)))

    def __sumConstrNum__(self):
        """It simply calculate constraint numbers."""
        numC = 0
        for constr in self.pointConstr:
            numC += constr.nf
        for constr, (start, end) in zip(self.pathConstr, self.pathConstrIndexPairs):
            tmpN = end - start
            if self.colloc_constr_is_on:
                numC += (2*tmpN - 1) * constr.nf
            else:
                numC += tmpN * constr.nf
        for constr in self.nonLinConstr:
            numC += constr.nf
        nnonlincon = numC
        for constr in self.linPointConstr:
            numC += constr.A.shape[0]
        for constr in self.linPathConstr:
            if self.colloc_constr_is_on:
                numC += constr.A.shape[0] * self.nPoint
            else:
                numC += constr.A.shape[0] * self.N
        for constr in self.linearConstr:
            numC += constr.A.shape[0]
        nlincon = numC - nnonlincon
        return numC, nnonlincon, nlincon

    def genGuessFromTraj(self, X=None, U=None, P=None, t0=None, tf=None, addx=None, tstamp=None,
                         obj=None, interp_kind='linear'):
        """Alias for :func:`~trajOptLib.TrajOptCollocProblem.gen_guess_from_traj`"""
        return self.gen_guess_from_traj(X, U, P, t0, tf, addx, tstamp, obj, interp_kind)

    def gen_guess_from_traj(self, X=None, U=None, P=None, t0=None, tf=None, addx=None, tstamp=None,
                         obj=None, interp_kind='linear'):
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
        teval = np.linspace(uset0, usetf, self.nPoint)

        # interpolation for state variables
        nPoint = self.nPoint
        dimx = self.dimx
        dimx_ = dimx // (self.sys.order + 1)
        if X is not None:
            Xcol = X.shape[1]
            if not (Xcol == dimx_ or Xcol == dimx):
                print('The column of X is not %d or %d, not use it' % (dimx, dimx_))
                X = None
            else:  # use interpolation to do it
                interp(tstamp, X, teval, Xtarget, interp_kind)
        if X is None:
            # straight path go there
            for i in range(self.nPoint):
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

    def genGuessFromSol(self, parsed_sol):
        """Alias for :func:`~trajOptLib.TrajOptCollocProblem.gen_guess_from_traj`"""
        return self.gen_guess_from_sol(parsed_sol)

    def gen_guess_from_sol(self, parsed_sol):
        """Generate an initial guess from a previous solution. Mainly change grid size or add perturbation. But determining structure is difficult

        :param parsed_sol: dictionary, output of calling parse_sol

        """
        t = parsed_sol['t']
        x = parsed_sol['x']
        u = parsed_sol['u']
        p = parsed_sol.get('p', None)
        addx = parsed_sol.get('addx', None)
        obj = parsed_sol.get('obj', None)
        return self.genGuessFromTraj(X=x, U=u, P=p, t0=t[0], tf=t[-1], addx=addx, tstamp=t, obj=obj, interp_kind='cubic')

    def __findMaxNG__(self):
        """Loop over all the constraints, find max NG. We then create temporary data for them."""
        maxnG = 0
        maxnG = max(maxnG, self.sys.nG)
        for constr in self.pointConstr:
            maxnG = max(maxnG, constr.nG)
        for constr in self.pathConstr:
            maxnG = max(maxnG, constr.nG)
        for constr in self.nonLinConstr:
            maxnG = max(maxnG, constr.nG)
        self.G = np.zeros(maxnG)
        self.row = np.zeros(maxnG, dtype=int)
        self.col = np.zeros(maxnG, dtype=int)

    def __analyzeObj__(self, numSol, numF):
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
            for i in range(self.nPoint):
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
        if catA.shape[0] == 0:
            spA = coo_matrix(([], ([], [])), shape=(addn + numF, numSol))
        else:
            spA = coo_matrix((catA, (catArow, catAcol)))
        return spA, addn

    def __setAPattern__(self, ndyncon, nnonlincon, spA):
        """Set sparsity pattern from linear constraints and objective functions.

        It finds sparsity pattern from defect constraints, linear constraints, and objective functions.
        It also finds sparsity pattern from those linear constraints.
        The A matrix from objective function is given in the sparse A and we just append it.
        The rows from linear constraints are straightforward.
        The constraints from defect constraints lie from rows 1 + ndyncon and occupies defectSize rows
        After nnonlincon rows (empty in A), we set linear constraints.

        The size of defect A: dynDefectSize = 2 * daeOrder * dimdyn and defectSize = dynDefectSize + dimu + dimp
        It sums up to 3*(dimu + dimp) + daeOrder*5*2*dimdyn nnz
        We have to do this for self.N - 1 mid-points.

        :param ndyncon: int, describes how many dynamics constraints we have
        :param nnonlincon: int, describes how many nonlinear constraints we have
        :param spA: sparse matrix, how the objective function is described linearly.

        """
        curRow, A, row, col = self.__setDefectPattern__(ndyncon, self.defectU, self.defectP)
        curRow += nnonlincon
        # we are ready to parse linear constraints
        lstCA, lstCArow, lstCAcol = self.__parseLinearConstraints__(curRow)
        # concatenate all those things together
        lstCA.append(spA.data)
        lstCA.append(A)
        lstCArow.append(spA.row)
        lstCArow.append(row)
        lstCAcol.append(spA.col)
        lstCAcol.append(col)
        Aval = np.concatenate(lstCA)
        Arow = np.concatenate(lstCArow)
        Acol = np.concatenate(lstCAcol)
        self.set_a_by_triplet(Aval, Arow, Acol)
        self.spA = csr_matrix((self.Aval, (self.Arow, self.Acol)), shape=(self.nf, self.nx))
        self.spA_coo = self.spA.tocoo()

    def __setDefectPattern__(self, ndyncon, withu=True, withp=True):
        """Set the sparse linear constraints from defect constraints.

        :param ndyncon: number of dynamical constraints. This sets starting row.
        :param withu: bool, determines if we set defect constraint on u
        :param withp: bool, determines if we set defect constraint on p

        """
        dimx, dimu, dimp, dimpoint, dimdyn = self.dimx, self.dimu, self.dimp, self.dimpoint, self.dimdyn
        if not withp:
            dimp = 0
        if not withu:
            dimu = 0
        if self.fixTimeMode:
            lenA = (10*self.daeOrder*self.dimdyn + 3*(dimu + dimp)) * (self.N - 1)
        else:
            lenA = (6*self.daeOrder*self.dimdyn + 3*(dimu + dimp)) * (self.N - 1)
        A = np.zeros(lenA)
        row = np.zeros(lenA, dtype=int)
        col = np.zeros(lenA, dtype=int)
        curNA = 0
        curRow = 1 + ndyncon
        # find those three matrix
        spL, spM, spR = self.__findMatLMRTemplate()
        for i in range(self.N - 1):
            midi = 2*i + 1
            lefti = 2*i
            righti = 2*(i + 1)
            for i in range(self.daeOrder):
                A[curNA: curNA + spL.nnz] = spL.data
                row[curNA: curNA + spL.nnz] = spL.row + curRow
                col[curNA: curNA + spL.nnz] = spL.col + lefti * dimpoint + i * dimdyn
                curNA += spL.nnz
                A[curNA: curNA + spM.nnz] = spM.data
                row[curNA: curNA + spM.nnz] = spM.row + curRow
                col[curNA: curNA + spM.nnz] = spM.col + midi * dimpoint + i * dimdyn
                curNA += spM.nnz
                A[curNA: curNA + spR.nnz] = spR.data
                row[curNA: curNA + spR.nnz] = spR.row + curRow
                col[curNA: curNA + spR.nnz] = spR.col + righti * dimpoint + i * dimdyn
                curNA += spR.nnz
                curRow += 2*dimdyn  # since size of spL, it is 2d by 2d
        # then do the constraint of u and p on nodes and knots, basically, midpoint is the average of two knots
        for i in range(self.N - 1):
            midi = 2*i + 1
            lefti = 2*i
            righti = 2*(i + 1)
            A[curNA: curNA + 2*dimu] = 0.5
            A[curNA + 2*dimu: curNA + 3*dimu] = -1
            A[curNA + 3*dimu: curNA + 3 * dimu + 2*dimp] = 0.5
            A[curNA + 3*dimu + 2 * dimp: curNA + 3*dimu + 3*dimp] = -1
            row[curNA: curNA + 3 * dimu] = curRow + np.tile(np.arange(dimu), (3, 1)).flatten()
            col[curNA: curNA + dimu] = lefti * dimpoint + dimx + np.arange(dimu)
            col[curNA + dimu: curNA + 2 * dimu] = righti * dimpoint + dimx + np.arange(dimu)
            col[curNA + 2 * dimu: curNA + 3 * dimu] = midi * dimpoint + dimx + np.arange(dimu)
            curNA_ = curNA + 3 * dimu
            row[curNA_: curNA_ + 3 * dimp] = curRow + dimu + np.tile(np.arange(dimp), (3, 1)).flatten()
            col[curNA_: curNA_ + dimp] = lefti * dimpoint + dimx + dimu + np.arange(dimp)
            col[curNA_ + dimp: curNA_ + 2 * dimp] = righti * dimpoint + dimx + dimu + np.arange(dimp)
            col[curNA_ + 2 * dimp: curNA_ + 3 * dimp] = midi * dimpoint + dimx + dimu + np.arange(dimp)
            curNA += 3 * (dimu + dimp)
            curRow += dimu + dimp
        return curRow, A, row, col

    def __parseLinearConstraints__(self, curRow):
        """Parse the linear constraints and form a sparse matrix.

        :param curRow: current row of accumulated constraints.
        """
        lstCA = []
        lstCArow = []
        lstCAcol = []
        for constr in self.linPointConstr:  # TODO: support for time to be done
            lstCA.append(constr.A.data)
            lstCArow.append(constr.A.row + curRow)
            lstCAcol.append(self.__patchCol__(constr.index, constr.A.col))  # take care on here
            curRow += constr.A.shape[0]
        for constr in self.linPathConstr:
            for j in range(self.nPoint):
                if not self.colloc_constr_is_on:
                    if j % 2 == 1:
                        continue
                index = j
                lstCA.append(constr.A.data)
                lstCArow.append(constr.A.row + curRow)
                lstCAcol.append(self.__patchCol__(index, constr.A.col))
                curRow += constr.A.shape[0]
        for constr in self.linearConstr:
            # the users have to be aware of the columns
            lstCA.append(constr.A.data)
            lstCArow.append(constr.A.row + curRow)
            lstCAcol.append(constr.A.col)
            curRow += constr.A.shape[0]
        return lstCA, lstCArow, lstCAcol

    def __findMatLMRTemplate(self):
        """Assume h is fixed, we find the L, M, R matrix for defect constraints.

        The goal is for a pair (q^(k), q^{(k+1)}) we want MatL*L+MatM*M+MatR*R=0
        where L, M, R are such pair at left point, mid-point, and right point.

        """
        d = self.dimdyn
        matL = np.zeros((2*d, 2*d))
        matM = np.zeros((2*d, 2*d))
        matR = np.zeros((2*d, 2*d))
        np.fill_diagonal(matL[:d, :d], 0.5)
        np.fill_diagonal(matL[d:, d:], -0.25)
        np.fill_diagonal(matM, -1)
        np.fill_diagonal(matR[:d, :d], 0.5)
        np.fill_diagonal(matR[d:, d:], -0.25)
        # time dependent parts
        if self.fixTimeMode:
            h = (self.tf - self.t0) / (self.N - 1)
            np.fill_diagonal(matL[:d, d:], h/8)
            np.fill_diagonal(matL[d:, :d], -1.5/h)
            np.fill_diagonal(matR[:d, d:], -h/8)
            np.fill_diagonal(matR[d:, :d], 1.5/h)
        spL = coo_matrix(matL)
        spM = coo_matrix(matM)
        spR = coo_matrix(matR)
        return spL, spM, spR

    def randomGenX(self):
        """Alias for :func:`trajOptLib.TrajOptCollocProblem.random_gen_guess`"""
        return self.random_gen_guess()

    def random_gen_guess(self):
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
        # I do not have to worry about objaddn since they are linear
        return randX

    def __turnOnGrad__(self, x0):
        """Turn on gradient, this is called after an initial x0 has been generated"""
        self.grad = True
        self.__getSparsity__(x0)

    def __getSparsity__(self, x0):
        """Detect sparsity of the problem with an initial guess."""
        numObjG = self.__getObjSparsity(x0)
        self.numObjG = numObjG
        # summarize number of pure linear constraints
        numDynG = self.__getDynSparsity(x0)
        numCG = 0  # G from C

        h, useT = self.__get_time_grid__(x0)
        useX, useU, useP = self.__parseX__(x0)
        i = np.random.randint(self.nPoint)
        tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))

        for constr in self.pointConstr:
            numCG += constr.nG
            constr.findTimeGradient(tmpx)
            if not constr.autonomous:
                n = len(constr.timeindex)
                numCG += (self.numT - 1) * n
        for constr, (start, end) in zip(self.pathConstr, self.pathConstrIndexPairs):
            tmpN = end - start
            if self.colloc_constr_is_on:
                numCG += (2*tmpN - 1) * constr.nG
            else:
                numCG += tmpN * constr.nG
            constr.findTimeGradient(tmpx)
            if not constr.autonomous:
                n = len(constr.timeindex)
                if self.colloc_constr_is_on:
                    numCG += (self.numT - 1) * n * self.nPoint
                else:
                    numCG += (self.numT - 1) * n * self.N
        for constr in self.nonLinConstr:
            numCG += constr.nG
        numG = numObjG + numDynG + numCG
        self.numG = numG
        self.nG = numG

    def __getDynSparsity(self, x):
        """Set sparsity of the problem caused by system dynamics and other nonlinear constraints.

        Sparsity pattern should be quite straight-forward to conclude since we introduce many auxiliary variables.
        The nG from dae system should be directly used and it is complete.
        The defect constraints introduce 5*dimq*4 gradients

        :param x: ndarray, the guess/sol

        """
        # this row is the case when numDefectDyn are imposed in G rather than A
        # dynG = self.sys.nG * self.nPoint + (20 * self.dimdyn + 3*(self.dimu + self.dimp)) * (self.N - 1)
        dynG = self.sys.nG * self.nPoint
        usex = x[:self.dimpoint]
        self.sys.findTimeGradient(usex)
        if not self.sys.autonomous:
            n = len(self.timeindex)
            dynG += self.nPoint * (self.numT - 1) * n
        # dynG arising from defect constraints
        if not self.fixTimeMode:
            dynG += (self.N - 1) * self.daeOrder * 4 * self.dimdyn  # those are purely from the defect matrix
            dynG += self.numT * (self.N - 1) * self.daeOrder * self.dimdyn * 2  # since time is free
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
        h, useT = self.__get_time_grid__(x)
        useX, useU, useP = self.__parseX__(x)
        i = np.random.randint(self.nPoint)
        tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
        nG = 0
        # check sparseObj mode
        if self.lqrObj is not None:
            nG += self.LQRnG  # this is always changing with time so do not worry
        for obj in self.nonLinObj:
            nG += obj.nG  # this assume user knows preciously what is going on so do not care
        for obj in self.nonPointObj:
            nG += obj.nG  # in this case, if not autonomous, depending on numT, we might have to
            obj.findTimeGradient(tmpx)  # detect pattern
            if not obj.autonomous:  # since it is objective, only one piece is time, we have counted it once
                nG += (self.numT - 1)  # so if numT=0, we remove it, 1 is fine, 2 we increase one more
        for obj in self.nonPathObj:
            nG += self.nPoint * obj.nG
            obj.findTimeGradient(tmpx)
            if not obj.autonomous:
                nG += self.nPoint * (self.numT - 1)  # this is okay, I guess
        return nG

    def __getTimeIndices(self):
        """Utility function for assigning sparsity structure."""
        t0ind = -1
        tfind = -1
        lenX = self.nPoint * self.dimpoint
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

    def __setXbound__(self):
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
        cind0 = self._setLinConstr(clb, cub, cind0)
        # the bounds for objaddn is 0 so we are good so far
        # assign to where it should belong to
        self.lb = clb
        self.ub = cub

    def _setNonLinConstr(self, clb, cub, cind0):
        for constr in self.pointConstr:
            if constr.lb is not None:
                clb[cind0: cind0 + constr.nf] = constr.lb
            if constr.ub is not None:
                cub[cind0: cind0 + constr.nf] = constr.ub
            cind0 += constr.nf
        for constr, (start, end) in zip(self.pathConstr, self.pathConstrIndexPairs):
            tmpN = end - start
            if self.colloc_constr_is_on:
                useN = 2 * tmpN - 1
            else:
                useN = tmpN
            tmplb = np.reshape(clb[cind0: cind0 + constr.nf * useN], (useN, constr.nf))
            tmpub = np.reshape(cub[cind0: cind0 + constr.nf * useN], (useN, constr.nf))
            cind0 += constr.nf * useN
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
        return cind0

    def _setLinConstr(self, clb, cub, cind0):
        # the rest are linear constraints and we should write those bounds, too
        for constr in self.linPointConstr:
            cindf = cind0 + constr.A.shape[0]
            clb[cind0: cindf] = constr.lb
            cub[cind0: cindf] = constr.ub
            cind0 = cindf
        for constr in self.linPathConstr:
            if self.colloc_constr_is_on:
                N = self.nPoint
            else:
                N = self.N
            cindf = cind0 + N * constr.A.shape[0]
            clb[cind0: cindf] = np.tile(constr.lb, (N, 1)).flatten()
            cub[cind0: cindf] = np.tile(constr.ub, (N, 1)).flatten()
            cind0 = cindf
        for constr in self.linearConstr:
            cindf = cind0 + constr.A.shape[0]
            clb[cind0: cindf] = constr.lb
            cub[cind0: cindf] = constr.ub
            cind0 = cindf
        return cind0

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
        if h <= 0 and not self.fixTimeMode:
            h = 1e-6
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
        """Alias for :func:`~trajOptLib.TrajOptCollocProblem.parse_f`"""
        return self.parse_f(guess)

    def parse_f(self, guess, y=None):
        """Give an guess, evaluate it and parse into parts.

        :param guess: ndarray, (numSol, ) a guess or a solution to check
        :param y: ndarray, (numF, ) if None, it stores the solution
        :returns: dict, containing objective and parsed constraints

        """
        assert len(guess) == self.numSol
        N = self.N
        nPoint = self.nPoint
        dimx = self.dimx
        dimdyn = self.dimdyn
        if y is None:
            y = np.zeros(self.numF)
        if self.grad:
            self.__callg__(guess, y, np.zeros(1), np.zeros(1), np.zeros(1), False, False)
        else:
            self.__callf__(guess, y)
        y += self.spA.dot(guess)
        obj = y[0]
        dynCon = np.reshape(y[1:(2*N-1)*dimdyn+1], (2*N - 1, dimdyn))
        curN = 1 + (2 * N - 1) * dimdyn
        curNf = curN + self.dynDefectSize * (N - 1)
        defectCon = dict()
        dynDefectCon = np.reshape(y[curN: curNf], (N - 1, -1))
        defectCon['dyn'] = dynDefectCon
        curN = curNf
        if self.defectU:
            dimu = self.dimu
        else:
            dimu = 0
        if self.defectP:
            dimp = self.dimp
        else:
            dimp = 0
        curNf += (dimu + dimp) * (N - 1)
        upDefectCon = np.reshape(y[curN: curNf], (N - 1, -1))
        if dimu > 0:
            defectCon['u'] = upDefectCon[:, :dimu]
        if dimp > 0:
            defectCon['p'] = upDefectCon[:, dimu:dimu+dimp]

        curN = curNf
        pointCon = []
        for constr in self.pointConstr:
            pointCon.append(y[curN: curN + constr.nf])
            curN += constr.nf
        pathCon = []
        for constr, (start, end) in zip(self.pathConstr, self.pathConstrIndexPairs):
            tmpN = end - start
            if self.colloc_constr_is_on:
                useN = 2*tmpN - 1
            else:
                useN = tmpN
            pathCon.append(np.reshape(y[curN: curN+useN*constr.nf], (useN, constr.nf)))
            curN += useN*constr.nf
        nonLinCon = []
        for constr in self.nonLinConstr:
            nonLinCon.append(y[curN: curN+constr.nf])
            curN += constr.nf
        # all linear constraints can be ignored, here we ignore them
        # check bounds, return a -1, 1 value for non-equality bounds, and 0 for equality bounds
        useX, useU, useP = self.__parseX__(guess)
        Xbound = checkInBounds(useX[:, :dimx], self.xbd)
        x0bound = checkInBounds(useX[0, :dimx], self.x0bd)
        xfbound = checkInBounds(useX[-1, :dimx], self.xfbd)
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
        useX, useU, useP = self.__parseX__(guess)
        objs = guess[self.numSol - self.objaddn:]
        rst = {'obj': obj, 'objs': objs, 'dyn': dynCon, 'defect': defectCon, 'point': pointCon, 'path': pathCon, 'nonlin': nonLinCon,
                'xbd': Xbound, 'ubd': ubound, 'x0bd': x0bound, 'xfbd': xfbound, 'pbd': pbound,
                't0bd': t0bound, 'tfbd': tfbound, 'addXbd': addXbound,
                'x': useX, 'u': useU, 'p': useP}
        if self.t0ind > 0:
            rst['t0'] = guess[self.t0ind]
        else:
            rst['t0'] = self.t0
        if self.tfind > 0:
            rst['tf'] = guess[self.tfind]
        else:
            rst['tf'] = self.tf
        rst['t'] = np.linspace(rst['t0'], rst['tf'], 2*N - 1)
        # parse addx
        if self.lenAddX > 0:
            addx = self.__parseAddX__(guess)
            for i, addx_ in enumerate(addx):
                rst['addx_%d' % i] = addx_
        return rst

    def __parseAddX__(self, x):
        numTraj = self.numTraj
        addX = []
        for addx in self.addX:
            addX.append(x[numTraj: numTraj + addx.n])
            numTraj += addx.n
        return addX

    def __parseObj__(self, x):
        return x[self.numSol - self.objaddn:]

    def getAddXIndexByIndex(self, i):
        """With i as index of addx, it returns the starting index in solution vector for this one.

        :param i: int, the index of addX we want to query.

        """
        index = self.numTraj
        for j in range(i):
            index += self.addX[j].n
        return index

    def getStateIndexByIndex(self, i):
        """With i as index for state variable, return the starting index in solution vector.

        :param i: int, the index of State we want to query

        """
        if i >= 0:
            return self.dimpoint * i
        else:
            return (self.nPoint + i) * self.dimpoint

    def getContrlIndexByIndex(self, i):
        """With i as index for control variable, return the starting index in solution vector.

        :param i: int, the index of control we want to query

        """
        if i >= 0:
            return self.dimpoint * i + self.dimx
        else:
            return (self.nPoint + i) * self.dimpoint + self.dimx

    def getParamIndexByIndex(self, i):
        """With i as index for parameter variable, return the starting index in solution vector.

        :param i: int, the index of parameter we want to query

        """
        if i >= 0:
            return self.dimpoint * i + self.dimx + self.dimu
        else:
            return (self.nPoint + i) * self.dimpoint + self.dimx + self.dimu

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix.

        :param x: ndarray, the solution to the problem
        :param y: ndarray, return F
        :param G/row/col: ndarray, information of gradient
        :param rec/needg: if we record/ if we need gradient

        """
        y[0] = 0  # since this row is purely linear
        h, useT = self.__get_time_grid__(x)
        useX, useU, useP = self.__parseX__(x)
        # loop over all system dynamics constraint
        curRow = 1
        curNg = 0
        curRow, curNg = self.__dynconstr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
        curRow, curNg = self.__constr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        curRow += self.numLinCon
        # loop over all the objective functions, I haven't checked if order is correct since some linear constraints are followed
        curRow, curNg = self.__obj_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
        return curRow, curNg

    def __dynconstr_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg):
        """Evaluate the constraints imposed by system dynamics"""
        dimx, dimu, dimp = self.dimx, self.dimu, self.dimp
        dimpoint = self.dimpoint
        dimdyn = self.dimdyn  # this works for many cases
        nPoint = self.nPoint
        # first let's check the 2*N - 1 dimdyn constraints from dynamics
        cDyn = np.reshape(y[curRow:curRow + nPoint * dimdyn], (nPoint, dimdyn))
        for i in range(nPoint):
            if self.sys.autonomous:
                Gpiece = G[curNg: curNg + self.sys.nG]
                rowpiece = row[curNg: curNg + self.sys.nG]
                colpiece = col[curNg: curNg + self.sys.nG]
                self.sys.dyn(useT[i], useX[i], useU[i], useP[i], cDyn[i], Gpiece, rowpiece, colpiece, rec, needg)
                if needg:
                    curNg += self.sys.nG
                    if rec:
                        rowpiece[:] += curRow
                        colpiece[:] = self.__patchCol__(i, colpiece[:])
            else:
                self.sys.dyn(useT[i], useX[i], useU[i], useP[i], cDyn[i], self.G, self.row, self.col, rec, needg)
                curNg = self.__copy_into_g__(i, G, row, col, curRow, curNg, self.sys.nG, self.sys.timeindex, False, rec,
                                             self.G, self.row, self.col)
            curRow += self.dimdyn
        # offset of row number due to defect dynamics constraint
        if self.fixTimeMode:
            curRow += self.numDefectDyn
        else:
            # manually set those defect dynamics and gradients, etc
            defectRow = (self.N - 1) * 2 * self.daeOrder * dimdyn
            cDefect = np.reshape(y[curRow: curRow + defectRow], (self.N - 1, self.daeOrder, 2, dimdyn))
            bscIndex = np.arange(dimdyn)
            for i in range(self.N - 1):
                lefti, righti = 2 * i, 2 * i + 2
                for j in range(self.daeOrder):
                    cDefect[i, j, 0, :] = h/8 * useX[lefti, (j+1)*dimdyn:(j+2)*dimdyn] - h/8 * useX[righti, (j+1)*dimdyn:(j+2)*dimdyn]
                    cDefect[i, j, 1, :] = -1.5/h * useX[lefti, j*dimdyn:(j+1)*dimdyn] + 1.5/h * useX[righti, j*dimdyn:(j+1)*dimdyn]
                    if needg:
                        G[curNg: curNg + dimdyn] = h / 8
                        G[curNg + dimdyn: curNg + 2*dimdyn] = -h / 8
                        G[curNg + 2*dimdyn: curNg + 3*dimdyn] = -1.5 / h
                        G[curNg + 3*dimdyn: curNg + 4*dimdyn] = 1.5 / h
                        if rec:
                            row[curNg:curNg + dimdyn] = curRow + bscIndex
                            row[curNg + dimdyn:curNg + 2*dimdyn] = curRow + bscIndex
                            col[curNg:curNg + dimdyn] = lefti * dimpoint + j * dimdyn + dimdyn + bscIndex
                            col[curNg + dimdyn:curNg + 2*dimdyn] = righti * dimpoint + j * dimdyn + dimdyn + bscIndex
                            row[curNg + 2*dimdyn:curNg + 3*dimdyn] = curRow + bscIndex + dimdyn
                            row[curNg + 3*dimdyn:curNg + 4*dimdyn] = curRow + bscIndex + dimdyn
                            col[curNg + 2*dimdyn:curNg + 3*dimdyn] = lefti * dimpoint + j * dimdyn + bscIndex
                            col[curNg + 3*dimdyn:curNg + 4*dimdyn] = righti * dimpoint + j * dimdyn + bscIndex
                        curNg += 4 * dimdyn
                        # if time related, we have to also consider them
                        if not self.fixTimeMode:
                            pc1pt = (useX[lefti, (j+1)*dimdyn:(j+2)*dimdyn] - useX[righti, (j+1)*dimdyn:(j+2)*dimdyn]) / 8
                            pc2pt = -(-useX[lefti, j*dimdyn:(j+1)*dimdyn] + useX[righti, j*dimdyn:(j+1)*dimdyn]) * 1.5 / h ** 2
                            if self.t0ind > 0:
                                G[curNg: curNg + dimdyn] = -pc1pt / (self.N - 1)
                                G[curNg + dimdyn: curNg + 2*dimdyn] = -pc2pt / (self.N - 1)
                                if rec:
                                    row[curNg: curNg + dimdyn] = curRow + bscIndex
                                    row[curNg + dimdyn: curNg + 2*dimdyn] = curRow + bscIndex + dimdyn
                                    col[curNg: curNg + 2*dimdyn] = self.t0ind
                                curNg += 2 * dimdyn
                            if self.tfind > 0:
                                G[curNg: curNg + dimdyn] = pc1pt / (self.N - 1)
                                G[curNg + dimdyn: curNg + 2*dimdyn] = pc2pt / (self.N - 1)
                                if rec:
                                    row[curNg: curNg + dimdyn] = curRow + bscIndex
                                    row[curNg + dimdyn: curNg + 2*dimdyn] = curRow + bscIndex + dimdyn
                                    col[curNg: curNg + 2*dimdyn] = self.tfind
                                curNg += 2*dimdyn
                    curRow += 2 * dimdyn
            if self.defectU:
                curRow += (self.N - 1) * self.dimu
            if self.defectP:
                curRow += (self.N - 1) * self.dimp
        return curRow, curNg

    def __copy_into_g__(self, index, G, row, col, curRow, curNg, nG, time_index, plus, rec,
                        G_src, row_src, col_src, col_offset=0):
        """With sparsity calculated in self.G, we assign to correct G.

        :param index: int, we are evaluating this at which point
        :param G, row, col: the G, row, col vector storing sparse Jacobian.
        :param curRow: accumulated row number.
        :param curNg: accumulated sparse Jacobian number
        :param nG: number of non-zero of Jacobian, this indicates how long of self.G we shall use
        :param timeindex: index indicating time related.
        :param plus: bool, if we plus value to time-related index (integral one)
        :param rec: bool, if we record index into row and col
        :param G_src/row_src/col_src: where we copy value from
        :param col_offset: int, offset of column, it is only used for multiple-phase problem
        :return curNg: updated occupied Ng
        """
        # use time index to build the mask for selecting data
        G_ = G_src[:nG]
        timemask = np.zeros(nG, dtype=bool)
        timemask[time_index] = True  # get a mask for time-related gradient
        statemask = np.logical_not(timemask)
        lenstate = nG - len(time_index)
        lenTime = nG - lenstate
        G[curNg: curNg + lenstate] = G_[statemask]
        if rec:
            col_ = col_src[:nG]
            row_ = row_src[:nG]
            col[curNg: curNg + lenstate] = col_[statemask] - 1 + index * self.dimpoint + col_offset
            row[curNg: curNg + lenstate] = row_[statemask] + curRow
        curNg += lenstate
        # for time related columns
        if self.t0ind > 0:
            ptpt0 = (self.N - 1 - index) / (self.N - 1)
            if plus:
                G[curNg: curNg + lenTime] += G_[time_index] * ptpt0
            else:
                G[curNg: curNg + lenTime] = G_[time_index] * ptpt0
            if rec:
                row[curNg: curNg + lenTime] = row_[time_index] + curRow
                col[curNg: curNg + lenTime] = self.t0ind + col_offset
            curNg += lenTime
        if self.tfind > 0:
            ptptf = index / (self.N - 1)
            if plus:
                G[curNg: curNg + lenTime] += G_[time_index] * ptptf
            else:
                G[curNg: curNg + lenTime] = G_[time_index] * ptptf
            if rec:
                row[curNg: curNg + lenTime] = row_[time_index] + curRow
                col[curNg: curNg + lenTime] = self.tfind + col_offset
            curNg += lenTime
        return curNg

    def __obj_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
        """Calculate objective function. Basically it moves all nonlinear objective function to the final rows
        See __constrModeG__ for arguments and output."""
        tmpout = np.zeros(1)
        # first lets do lqrobj
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
                if obj.autonomous:
                    Gpiece = G[curNg: curNg + obj.nG]
                    rowpiece = row[curNg: curNg + obj.nG]
                    colpiece = col[curNg: curNg + obj.nG]
                    obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                    if rec:
                        rowpiece[:] = curRow
                        colpiece[:] = self.__patchCol__(obj.index, colpiece)
                    curNg += obj.nG
                else:
                    obj.__callg__(tmpx, tmpout, self.G, self.row, self.col, rec, needg)
                    curNg = self.__copy_into_g__(obj.index, G, row, col, curRow, curNg, obj.nG, obj.timeindex, False, rec, self.G, self.row, self.col)
                y[curRow] = tmpout[0]
                curRow += 1

        if len(self.nonPathObj) > 0:
            weight = np.zeros((self.nPoint, 1))
            weight[1::2] = 2.0 / 3.0 * h
            weight[0::2] = 1.0 / 3.0 * h
            weight[0] = 1.0 / 6.0 * h
            weight[-1] = 1.0 / 6.0 * h
            for obj in self.nonPathObj:
                y[curRow] = 0
                if obj.autonomous:
                    for i in range(self.nPoint):
                        tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                        Gpiece = G[curNg: curNg + obj.nG]
                        rowpiece = row[curNg: curNg + obj.nG]
                        colpiece = col[curNg: curNg + obj.nG]
                        obj.__callg__(tmpx, tmpout, Gpiece, rowpiece, colpiece, rec, needg)
                        Gpiece[:] *= weight[i]
                        if rec:
                            rowpiece[:] = curRow
                            colpiece[:] = self.__patchCol__(i, colpiece)
                        curNg += obj.nG
                else:
                    for i in range(self.nPoint):
                        tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                        obj.__callg__(tmpx, tmpout, self.G, self.row, self.col, rec, needg)
                        self.G[:obj.nG] *= weight[i]
                        curNg = self.__copy_into_g__(i, G, row, col. curRow, curNg, obj.nG, obj.timeindex, True, rec, self.G, self.row, self.col)
                y[curRow] += weight[i] * tmpout[0]
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

    def __constr_mode_g__(self, curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg):
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
        if len(self.pointConstr) > 0:
            for constr in self.pointConstr:
                tmpx = np.concatenate(([useT[constr.index]], useX[constr.index], useU[constr.index], useP[constr.index]))
                if constr.autonomous:
                    pieceG = G[curNg: curNg + constr.nG]
                    pieceRow = row[curNg: curNg + constr.nG]
                    pieceCol = col[curNg: curNg + constr.nG]
                    constr.__callg__(tmpx, y[curRow: curRow + constr.nf], pieceG, pieceRow, pieceCol, rec, needg)
                    if rec:
                        pieceRow += curRow
                        pieceCol[:] = self.__patchCol__(constr.index, pieceCol)
                    curNg += constr.nG
                else:
                    constr.__callg__(tmpx, y[curRow: curRow + constr.nf], self.G, self.row, self.col, rec, needg)
                    curNg = self.__copy_into_g__(constr.index, G, row, col, curRow, curNg, constr.nG, constr.timeindex,
                                                 True, rec, self.G, self.row, self.col)
                curRow += constr.nf
        if len(self.pathConstr) > 0:
            for constr, (start, end) in zip(self.pathConstr, self.pathConstrIndexPairs):
                if constr.autonomous:
                    for j in range(2*start, 2*end - 1):
                        if not self.colloc_constr_is_on:
                            if j % 2 == 1:
                                continue
                        i = j
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
                else:
                    for j in range(2*start, 2*end - 1):
                        if not self.colloc_constr_is_on:
                            if j % 2 == 1:
                                continue
                        i = j
                        tmpx = np.concatenate(([useT[i]], useX[i], useU[i], useP[i]))
                        constr.__callg__(tmpx, y[curRow: curRow + constr.nf], self.G, self.row, self.col, rec, needg)
                        curNg = self.__copy_into_g__(i, G, row, col, curRow, curNg, constr.nG, constr.timeindex, True, rec,
                                                     self.G, self.row, self.col)
                        curRow += constr.nf
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
    def __cost__(self, x):
        """The eval_f function required by ipopt.

        :param x: a guess/solution of the problem
        :return f: float, objective function

        """
        row0 = self.spA.getrow(0)
        return np.dot(row0.data, x[row0.indices])

    def __gradient__(self, x, g):
        """Evaluation of the gradient of objective function.

        :param x: guess/solution to the problem
        :return grad: gradient of objective function w.r.t x

        """
        g[:] = self.spA.getrow(0).toarray().flatten()
        return True

    def __constraint__(self, x, y):
        """Evaluation of the constraint function.

        :param x: ndarray, guess/solution to the problem.
        :return g: constraint function

        """
        G = np.zeros(1)
        row = np.zeros(1, dtype=int)
        col = np.zeros(1, dtype=int)
        self.__callg__(x, y, G, row, col, False, False)
        # y should plus A times x
        y += self.spA.dot(x)
        return y

    def __jacobian__(self, x, g, row, col, rec):
        """Evaluate jacobian of constraints. I simply call __callg__

        :param x: ndarray, guess / solution to the problem
        :param flag: bool, True return row/col, False return values

        """
        y = np.zeros(self.nf)
        if rec:
            row = np.ones(self.nG + self.spA.nnz, dtype=int)
            col = np.ones(self.nG + self.spA.nnz, dtype=int)
            tmpx = self.randomGenX()
            self.__callg__(tmpx, y, g, row, col, True, True)
            # good news is there is no overlap of A and G
            row[self.nG:] = self.spA_coo.row
            col[self.nG:] = self.spA_coo.col
            return row, col
        else:
            row = np.ones(1, dtype=int)
            col = np.ones(1, dtype=int)
            self.__callg__(x, y, g, row, col, False, True)
            g[self.nG:] = self.spA_coo.data
            return g

    def __patchCol__(self, index, col, col_offset=0):
        """Find which indices it belongs to the original one for a local matrix at index with col.

        Since we have changed how variables are arranged, now it should be quite straightforward to do so.

        """
        col = col[col > 0]  # get rid of those with time
        return col - 1 + self.getStateIndexByIndex(index) + col_offset

    def parse_sol(self, sol):
        """Call parseX function from utility and return a dict of solution."""
        X, U, P = self.__parseX__(sol)
        if self.dimp == 0:
            P = None
        h, tgrid = self.__get_time_grid__(sol)
        obj = self.__parseObj__(sol)
        if self.lenAddX == 0:
            return {'t': tgrid, 'x': X, 'u': U, 'p': P, 'obj': obj}
        else:
            return {'t': tgrid, 'x': X, 'u': U, 'p': P, 'addx': self.__parseAddX__(sol), 'obj': obj}

    def addLQRObj(self, lqrobj):
        """Alias for :func:`trajOptLib.TrajOptCollocProblem.add_lqr_obj`"""
        return self.add_lqr_obj(lqrobj)

    def add_lqr_obj(self, lqrobj):
        """Add a lqr objective function to the problem. It changes lqrObj into a function being called.

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
            :param useX, useU, useP: ndarray, parsed X, U, P from x
            :param y: ndarray, a location to write the objective function onto
            :param G/row/col: the gradient information
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
            if needg:
                G_ = np.reshape(G[:num1], (nPoint, -1))  # use the same structure with col_
            if useQ > 0:
                yQ = np.dot(weight[:, 0], np.sum(((useX[:, Qcol] - lqrobj.xbase[Qcol]) ** 2) * lqrobj.Q.data, axis=1)) * h
                if needg:
                    G_[:, :useQ] = 2.0 * h * (weight * ((useX[:, Qcol] - lqrobj.xbase[Qcol]) * lqrobj.Q.data))
                    if rec:
                        col_[:, :useQ] = Qcol + baseCol
                    curG += nPoint * useQ
            if useR > 0:
                yR = np.dot(weight[:, 0], np.sum(((useU[:, Rcol] - lqrobj.ubase[Rcol]) ** 2) * lqrobj.R.data, axis=1)) * h
                if needg:
                    G_[:, useQ: useQ + useR] = (2.0 * h * (weight * ((useU[:, Rcol] - lqrobj.ubase[Rcol]) * lqrobj.R.data)))
                    if rec:
                        col_[:, useQ: useQ+useR] = Rcol + baseCol + self.dimx
                    curG += nPoint * useR
            if useP > 0:
                yP = np.dot(weight[:, 0], np.sum(((useP_[:, Pcol] - lqrobj.pbase[Pcol]) ** 2) * lqrobj.P.data, axis=1)) * h
                if needg:
                    G_[:, useQ+useR: useQ+useR+useP] = (2.0 * h * (weight * ((useP_[:, Pcol] - lqrobj.pbase[Pcol]) * lqrobj.P.data)))
                    if rec:
                        col_[:, useQ+useR: useQ+useR+useP] = Pcol + baseCol + self.dimx + self.dimu
                    curG += nPoint * useP
            if useF > 0:
                yF = np.sum(lqrobj.F.data * ((useX[-1, Fcol] - lqrobj.xfbase[Fcol]) ** 2))
                if useQ > 0:
                    if needg:
                        n0 = curG - (useR + useP) - useF  # locate at the last row
                        G[n0: n0 + useF] += 2.0 * lqrobj.F.data * (useX[-1, Fcol] - lqrobj.xfbase[Fcol])
                else:
                    G[curG: curG + useF] = 2.0 * lqrobj.F.data * (useX[-1, Fcol] - lqrobj.xfbase[Fcol])
                    if rec:
                        col[curG: curG + useF] = Fcol + baseCol[-1]
                    curG += useF
            if needg:
                if self.t0ind > 0:
                    if h <= 1e-8:
                        pass
                    G[curG] = -(yQ + yR + yP) / h / (self.N - 1) - tfweight
                    if rec:
                        col[curG: curG+1] = self.t0ind
                    curG += 1
                if self.tfind > 0:
                    G[curG] = (yQ + yR + yP) / h / (self.N - 1) + tfweight
                    if rec:
                        col[curG: curG+1] = self.tfind
                    curG += 1
            y[0] = yF + yQ + yR + yP + yTf

        self.lqrObj = __callg__

    def addVanillaQuadPenalty(self, indices, weights):
        """Add a quadratic penalty term based on indices in the solution vector and weights.

        This is dangerous and might cause unexpected trouble unless you are sure indices are correct.
        You can always use addStateQuadPenalty/addControlQuadPenalty/addParamQuadPenalty/addAddXQuadPenalty to finish these.

        :param indices: indices of variables to be penalized
        :param weights: float/ndarray weights associated with those variables.

        """
        self.addNonLinearObj(QuadPenalty(indices, weights))

    def addStateQuadPenalty(self, index, weights, mask=None):
        """Add a quadratic penalty on selected state variables.

        :param index: int/array-like, indices of state variables to be penalized.
        :param weights: float/array-like, weights of penalty
        :param mask: mask-like, filter for selecting subset of variables

        """
        stateindex = np.array(index)
        if np.isscalar(weights):
            useweight = weights
        else:
            useweight = []
        index = []
        for idx in stateindex:
            i0 = self.getStateIndexByIndex(idx)
            if mask is None:
                index.append(np.arange(i0, i0 + self.dimx))
            else:
                index.append(np.arange(i0, i0 + self.dimx)[mask])
            if not np.isscalar(weights):
                useweight.append(weights)
        indexes = np.concatenate(index)
        if isinstance(weights, list):
            useweight = np.concatenate(useweight)
        self.addVanillaQuadPenalty(indexes, useweight)

    def addControlQuadPenalty(self, index, weights, mask=None):
        """Add a quadratic penalty on selected control variables.

        :param index: int/array-like, indices of control variables to be penalized.
        :param weights: float/array-like, weights of penalty
        :param mask: filter to select subset

        """
        ctrlindex = np.array(index)
        if np.isscalar(weights):
            useweight = weights
        else:
            useweight = []
        index = []
        for idx in ctrlindex:
            i0 = self.getControlIndexByIndex(idx)
            if mask is None:
                index.append(np.arange(i0, i0 + self.dimu))
            else:
                index.append(np.arange(i0, i0 + self.dimu)[mask])
            if isinstance(weights, list):
                useweight.append(weights)
        indexes = np.concatenate(index)
        if isinstance(weights, list):
            useweight = np.concatenate(useweight)
        self.addVanillaQuadPenalty(indexes, useweight)

    def addParamQuadPenalty(self, index, weights, mask=None):
        """Add a quadratic penalty on selected parameter variables.

        :param index: int/array-like, indices of parameter variables to be penalized.
        :param weights: float/array-like, weights of penalty
        :param mask: filter of variables

        """
        paramindex = np.array(index)
        if np.isscalar(weights):
            useweight = weights
        else:
            useweight = []
        index = []
        for idx in index:
            i0 = self.getParamIndexByIndex(idx)
            if mask is None:
                index.append(np.arange(i0, i0 + self.dimp))
            else:
                index.append(np.arange(i0, i0 + self.dimp)[mask])
            if isinstance(weights, list):
                useweight.append(weights)
        indexes = np.concatenate(index)
        if isinstance(weights, list):
            useweight = np.concatenate(useweight)
        self.addVanillaQuadPenalty(indexes, useweight)

    def addAddXQuadPenalty(self, index, weights, mask=None):
        """Add quadratic penalty to addx variables.

        The API is slightly different from previous three since we cannot guarantee weights are of the same lenght.

        :param index: int, indices of parameter variables to be penalized.
        :param weights: float/array-like, weights of penalty
        :param mask: filter

        """
        assert isinstance(index, int)
        i0 = self.getAddXIndexByIndex(index)
        naddx = self.addX[index].n
        if mask is None:
            indexes = np.arange(i0, i0 + naddx)
        else:
            indexes = np.arange(i0, i0 + naddx)[mask]
        self.addVanillaQuadPenalty(indexes, weights)

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

    def addNonLinearPointConstr(self, pntConstr, path=False, **kwargs):
        """Add point constraint.

        :param pntConstr: pointConstr class
        :param path: bool, if this obj
        :kwargs: additional parameters, users can specify starting and ending indexes by specifying start and end

        """
        assert isinstance(pntConstr, nonLinearPointConstr)
        if path:
            start = kwargs.get('start', 0)
            end = kwargs.get('end', self.N)
            if end <= 0:
                end = self.N + end
            self.pathConstr.append(pntConstr)
            self.pathConstrIndexPairs.append((start, end))
        else:
            self.pointConstr.append(pntConstr)

    def addNonLinearConstr(self, constr):
        """Add a general nonlinear constraint.

        :param constr: nonLinConstr class

        """
        assert isinstance(constr, nonLinearConstr)
        self.nonLinConstr.append(constr)

    def addLinearConstr(self, constr):
        """Add a linear constraint to the problem.

        :param constr: a linearConstr object

        """
        assert isinstance(constr, linearConstr)
        self.linearConstr.append(constr)

    def addLinearPointConstr(self, constr, path=False):
        """Add a linear point constraint to the problem.

        :param constr: a linearPointConstr object
        :param path: if this constraint is path constraint

        """
        assert isinstance(constr, linearPointConstr)
        if path:
            self.linPathConstr.append(constr)
        else:
            self.linPointConstr.append(constr)

    def addObj(self, obj, path=False):
        """Alias for :func:`~trajOptLib.TrajOptCollocProblem.addObj`"""
        self.add_obj(obj, path)

    def add_obj(self, obj, path=False):
        """A high level function that add objective function of any kind.

        :param obj: an objective object.
        :param path: bool, if the point objective is an integral one.

        """
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
        else:
            print("Inappropriate type %s used as objective" % type(obj))

    def addConstr(self, constr, path=False, **kwargs):
        """Alias for :func:`~trajOptLib.TrajOptCollocProblem.add_constr`"""
        self.add_constr(constr, path, **kwargs)

    def add_constr(self, constr, path=False, **kwargs):
        """Add a constraint to the problem.

        :param constr: a constraint object.
        :param path: bool, if this constraint is a path constraint. Only applies for point constraint.

        """
        if isinstance(constr, linearConstr):
            self.addLinearConstr(constr)
        elif isinstance(constr, linearPointConstr):
            self.addLinearPointConstr(constr, path)
        elif isinstance(constr, nonLinearConstr):
            self.addNonLinearConstr(constr)
        elif isinstance(constr, nonLinearPointConstr):
            self.addNonLinearPointConstr(constr, path, **kwargs)
        else:
            print("Inappropriate type %s used as constraint" % type(constr))

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

    def set_x_bound(self, xlb, xub):
        """Set bounds on state variables.

        :param xlb: ndarray, (dimx,) lower bounds on state variables.
        :param xub: ndarray, (dimx,) upper bounds on state variables.

        """
        if len(xlb) != self.dimx:
            print('Incorrect length of xlb, it must be %d' % self.dimx)
        if len(xub) != self.dimx:
            print('Incorrect length of xub, it must be %d' % self.dimx)
        self.xbd = [np.array(xlb), np.array(xub)]

    def set_u_bound(self, ulb, uub):
        """Set bounds on control variables.

        :param ulb: ndarray, (dimu,) lower bounds on control variables.
        :param uub: ndarray, (dimu,) upper bounds on control variables.

        """
        if len(ulb) != self.dimu:
            print('Incorrect length of ulb, it must be %d' % self.dimu)
        if len(uub) != self.dimu:
            print('Incorrect length of uub, it must be %d' % self.dimu)
        self.ubd = [np.array(ulb), np.array(uub)]

    def set_p_bound(self, plb, pub):
        """Set bounds on parameter variables.

        :param plb: ndarray, (dimp,) lower bounds on parameter variables.
        :param pub: ndarray, (dimp,) upper bounds on parameter variables.

        """
        if len(plb) != self.dimp:
            print('Incorrect length of plb, it must be %d' % self.dimp)
        if len(pub) != self.dimp:
            print('Incorrect length of pub, it must be %d' % self.dimp)
        self.pbd = [np.array(plb), np.array(pub)]

    def set_x0_bound(self, x0lb, x0ub):
        """Set bounds on x0. This is optional but useful.

        :param x0lb: ndarray, (dimx,) lower bounds on x0 variables.
        :param x0ub: ndarray, (dimx,) upper bounds on x0 variables.

        """
        if len(x0lb) != self.dimx:
            print('Incorrect length of x0lb, it must be %d' % self.dimx)
        if len(x0ub) != self.dimx:
            print('Incorrect length of x0ub, it must be %d' % self.dimx)
        self.x0bd = [np.array(x0lb), np.array(x0ub)]

    def set_xf_bound(self, xflb, xfub):
        """Set bounds on xf. This is optional but useful.

        :param xflb: ndarray, (dimx,) lower bounds on xf variables.
        :param xfub: ndarray, (dimx,) upper bounds on xf variables.

        """
        if len(xflb) != self.dimx:
            print('Incorrect length of xflb, it must be %d' % self.dimx)
        if len(xfub) != self.dimx:
            print('Incorrect length of xfub, it must be %d' % self.dimx)
        self.xfbd = [np.array(xflb), np.array(xfub)]
