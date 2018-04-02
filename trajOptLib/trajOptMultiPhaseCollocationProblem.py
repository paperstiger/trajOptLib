#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
trajOptMultiPhaseCollocationProblem.py

This class implements methods for solving problem with multiple phases.
I will also change code style to suppress those warnings.
"""
from __future__ import division
import numpy as np
from .trajOptBase import linearPointObj
from .trajOptBase import linearPointConstr, linearConstr
from .trajOptBase import nonLinearPointObj, nonLinearObj
from .trajOptBase import nonLinearPointConstr, nonLinearConstr
from .trajOptBase import addX
from .libsnopt import snoptConfig, probFun, solver, result
from .utility import randomGenInBound
from scipy.sparse import coo_matrix
from .trajOptCollocationProblem import trajOptCollocProblem


class NonLinearConnectConstr(object):
    """Class for defining point constraint function."""
    def __init__(self, phase1, phase2, nc, lb=None, ub=None, nG=None, index1=-1, index2=0, addx_index=None):
        """Constructor for nonlinear point constraint. Also serve as path constraint.

        :param phase1/phase2: int, phase number that constraint is imposed. We read data from them
        :param nc: int, dimension of constraint function
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param nG, int, number of nnz of Jacobian
        :param index1/index2: int, by default we connect last point of phase 1 with first point of phase 2.
        However, we allow manual assignment of them.
        :param addx_index: int, which additional x we use to connect those two phases

        """
        self.phases = (phase1, phase2)
        self.indexes = (index1, index2)
        self.nf = nc
        self.nG = nG
        self.human = False  # if we call __callhumang__ instead of __callg__
        self.autonomous = (False, False)
        self.timeindex = (None, None)
        self.addx_index = addx_index
        if lb is None:
            self.lb = np.zeros(self.nf)
        else:
            self.lb = lb
        if ub is None:
            self.ub = np.zeros(self.nf)
        else:
            self.ub = ub

    def __callg__(self, x1, x2, F, G, row, col, rec, needg, addx=None):
        """Call and evaluate the constraint function. x1 and x2 are points in phase 1 and phase 2.
        
        I return them together in G, row, col since it might be derived from sympy so thay are together.
        But, human derived formulas love the previous one better
        """
        raise NotImplementedError

    def __callhumang__(self, x1, x2, needg, addx=None):
        """Call and evaluate the constraint function. x1 and x2 are points in phase 1 and phase 2.
        
        I return them together in G, row, col since it might be derived from sympy so thay are together.
        But, human derived formulas love the previous one better

        :param x1: ndarray, a (t, x, u, p) tuple from phase 1
        :param x2: ndarray, a (t, x, u, p) tuple from phase 2
        :param needg: if we return those G, too
        :param addx: additional parameters passed in
        :return F: where constraint functions are evaluated and stored
        :return G1, G2, G3: Jacobian of F w.r.t x1, x2, and optionally addx (G3 empty if addx is None)

        """
        raise NotImplementedError

    def find_time_gradient(self, x1, x2, addx=None):
        """Find where the gradients are w.r.t time."""
        if self.human:
            rvalue = self.__callhumang__(x1, x2, True, addx)
            if self.addxindex is None:
                F, G1, G2, = rvalue
                if self.nG is None:
                    self.nG = G1.nnz + G2.nnz
            else:
                F, G1, G2, Ga = rvalue
                if self.nG is None:
                    self.nG = G1.nnz + G2.nnz + Ga.nnz
            tindex1 = np.where(G1.col == 0)[0]
            tindex2 = np.where(G2.col == 0)[0]
        else:
            tmpy = np.zeros(self.nf)
            tmpG = np.zeros(self.nG)
            tmprow = np.zeros(self.nG, dtype=int)
            tmpcol = np.zeros(self.nG, dtype=int)
            self.__callg__(x1, x2, tmpy, tmpG, tmprow, tmpcol, True, True, addx)
            len1 = len(x1)
            tindex1 = np.where(tmpcol == 0)[0]
            tindex2 = np.where(tmpcol == len1)[0]
        auto1 = len(tindex1) == 0
        auto2 = len(tindex2) == 0
        self.timeindex = (tindex1, tindex2)
        self.autonomous = (auto1, auto2)
        if self.lb is None:
            self.lb = np.zeros(self.nf)
        if self.ub is None:
            self.ub = np.zeros(self.nf)


class LinearConnectConstr(object):
    """Class for defining linear constraint functions."""
    def __init__(self, phase1, phase2, a1, a2, lb=None, ub=None, index1=-1, index2=0, adda=None, addx_index=None):
        """Constructor for nonlinear point constraint. Also serve as path constraint.

        :param phase1/phase2: int, phase number that constraint is imposed. We read data from them
        :param a1/a2: ndarray, the constraint, :math: `l \le A_1 x_1 + A_2 x_2 \le u`
        :param lb, ub: lower and upper bound of the constraint function. None means equal to 0
        :param nG, int, number of nnz of Jacobian
        :param index1/index2: int, by default we connect last point of phase 1 with first point of phase 2.
        However, we allow manual assignment of them.
        :param adda: ndarray, if additional parameter comes into play, we use this
        :param addx_index: int, index of the additional parameters

        """
        self.phases = (phase1, phase2)
        self.indexes = (index1, index2)
        assert a1.shape[0] == a2.shape[0]
        self.a = (coo_matrix(a1), coo_matrix(a2))  # since everything is in pair
        self.nf = a1.shape[0]
        self.adda = coo_matrix(adda)
        self.addx_index = addx_index
        self.find_time_gradient()
        if lb is None:
            self.lb = np.zeros(self.nf)
        else:
            self.lb = lb
        if ub is None:
            self.ub = np.zeros(self.nf)
        else:
            self.ub = ub

    def find_time_gradient(self):
        """For a matrix, find column 0 indice."""
        timeindex1 = np.where(self.a[0].col == 0)[0]
        auto1 = len(timeindex1) == 0
        timeindex2 = np.where(self.a[1].col == 0)[0]
        auto2 = len(timeindex2) == 0
        self.timeindex = (timeindex1, timeindex2)
        self.autonomous = (auto1, auto2)


class TrajOptMultiPhaseCollocProblem(probFun):
    """A class for definition of trajectory optimization problem using collocation constraints with support for phase transition.

    The support of multiple phase is implemented by calling functions defined in trajOptCollocProblem since I do not want to reinvent the wheels.
    I can conveniently add support for other constraints that connects two phases.

    """
    def __init__(self, probs, addx=None):
        """Initialize a multi-phase problem by a list of trajOptColloProblem objects.

        :param probs: a list of trajOptCollocProblem
        :param addx: a list of addX

        """
        assert isinstance(probs, list)
        for prob in probs:
            assert isinstance(prob, trajOptCollocProblem)
            prob.preProcess()
        self.__parseAddX(addx)
        self.phases = probs
        self.num_phase = len(probs)
        # additional constraints, linear or nonlinear
        self.connect_nonlinear_constr = []  # general nonlinear constraint functions
        self.connect_linear_constr = []  # can be used for time
        self.nonlinear_constr = []  # a general constraint function, do not use it until necessary
        self.addx_nonlinear_constr = []  # constraints on addx, this might be useful
        self.addx_linear_constr = []  # linear constraints on addx, this is cheap and not a big deal
        # additional objective functions on addx or general
        self.addx_nonlinear_obj = []  # cost function associated addx
        self.addx_linear_obj = []  # always extends functions that probably no one uses
        self.nonlinear_obj = []  # a general objective function on the whole x
        # data for determining size of the problem
        vec_num_sol = np.array([prob.numSol for prob in self.phases])
        self.accum_num_sol = np.insert(np.cumsum(vec_num_sol), 0, 0)
        vec_num_f = np.array([prob.numF - 1 for prob in self.phases])  # -1 to remove objf row
        self.accum_num_f = np.insert(np.cumsum(vec_num_f), 0, 0)
        self.sum_num_f = self.accum_num_f[-1]
        self.__add_connect_time_constr()

    def pre_process(self):
        """Initialize the instances of probFun now we are ready.

        Call this function after the objectives and constraints have been set appropriately.
        It calculate the space required for SNOPT and allocates sparsity structure if necessary.

        """
        est_num_F = 1 + self.sum_num_f  # without additional param for objective, this is
        # get additional constraints number
        linear_F = 0
        for constr in self.connect_linear_constr:
            est_num_F += constr.nf
            linear_F += constr.nf
        for constr in self.connect_nonlinear_constr:
            est_num_F += constr.nf
        for constr in self.nonlinear_constr:
            est_num_F += constr.nf
        for constr in self.addx_nonlinear_constr:
            est_num_F += constr.nf
        for constr in self.addx_linear_constr:
            est_num_F += constr.nf
            linear_F += constr.nf
        est_num_sol = self.accum_num_sol[-1] + self.len_addX  # see previous line
        probFun.__init__(self, est_num_F, est_num_sol)  # we have to allocate A related variables
        self.__set_A_pattern(est_num_sol, est_num_F)
        # we can determine num_f and num_sol for this class
        self.num_f = est_num_F + self.num_aux_obj_var
        self.num_sol = est_num_sol + self.num_aux_obj_var
        self.num_linear_constr = linear_F
        self.nx = self.num_sol
        self.nf = self.num_f
        # find number of G
        rdx = self.random_gen_x()
        nG = self.__find_num_G(rdx)  # TODO: allow user not to give nG if given in human-mode?
        self.numG = nG
        self.nG = nG
        self.grad = True
        self.__set_x_bound()
        self.__set_f_bound()

    def random_gen_x(self):
        """Generate a random guess to the problem."""
        randX = np.zeros(self.num_sol)
        for i, phase in enumerate(self.phases):
            piecex = self.__get_phase_sol_piece(randX, i)
            piecex[:] = phase.randomGenX()
        # randomly generate for addX
        if self.len_addX > 0:
            for field, addx in zip(self.__parseAddX__(randX), self.addX):
                field[:] = randomGenInBound([addx.lb, addx.ub], addx.n)
        # the num_aux_obj_var auxiliary variables need no effort
        return randX

    def __callg__(self, x, y, G, row, col, rec, needg):
        """Evaluate those constraints, objective functions, and constraints. It simultaneously allocates sparsity matrix.

        We do this in several parts:
        - for each phase, we call the private functions explicitly. Ugly hack
        - call the connection nonlinear constraints
        - call the nonlinear objective functions associated with addx

        :param x: ndarray, the solution to the problem
        :param y: ndarray, return F
        :param G, row, col: ndarray, information of gradient
        :param rec, needg: if we record/ if we need gradient

        """
        curRow = 1
        curNg = 0
        for phase_num, phase in enumerate(self.phases):
            xi = self.__get_phase_sol_piece(x, phase_num)
            h, useT = phase.__get_time_grid__(xi)
            useX, useU, useP = phase.__parseX__(xi)
            # loop over all system dynamics constraint
            Ng0 = curNg
            curRow, curNg = phase.__dynconstr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, y, G, row, col, rec, needg)
            curRow, curNg = phase.__constr_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
            y[curRow: curRow + phase.numLinCon] = 0
            curRow += phase.numLinCon
            curRow, curNg = phase.__obj_mode_g__(curRow, curNg, h, useT, useX, useU, useP, x, y, G, row, col, rec, needg)
            col[Ng0: curNg] += self.accum_num_sol[phase_num]
        y[curRow: curRow + self.num_linear_constr] = 0
        curRow += self.num_linear_constr
        # evaluate the nonlinear constraints
        curRow, curNg = self.__calc_nonlinear_constr(curRow, curNg, x, y, G, row, col, rec, needg)
        # evaluate additional nonlinear objective functions
        if self.num_aux_obj_var > 0:
            curRow, curNg = self.__calc_nonlinear_obj(curRow, curNg, x, y, G, row, col, rec, needg)
        else:
            y[0] = 0  # just to make sure

    def parse_sol(self, sol):
        """Given a solution, we parse and return readable data structure.
        
        :param sol: ndarray or result object returned by SNOPT.
        :return traj: a dict of keys 't', 'x', 'u', 'p', 'phases'.
        - 't': is a concatenated time grid (might not be equidistant, but piecewise equidistant
        - 'x': ndarray, (*, dimx) is the concatenated state vector of all phases
        - 'u': ndarray, (*, dimu) is the concatenated control vector of all phases
        - 'p': ndarray, (*, dimp) is the concatenated parameter vector of all phases
        - 'phases': list of dictionaries composed of keys 't', 'x', 'u', 'p', 'addx' where 'addx' is
        additional optimization variables.

        """
        if isinstance(sol, np.ndarray):
            assert len(sol) == self.nx
            x = sol
        elif isinstance(sol, result):
            assert len(sol.sol) == self.nx
            x = sol.sol
        else:
            raise NotImplementedError
        traj = {'t': [], 'x': [], 'u': [], 'p': None}
        phases = []  # this is a list of parsed solution
        for phase_num, phase in enumerate(self.phases):
            soli = self.__get_phase_sol_piece(x, phase_num)
            rsti = phase.parseSol(soli)
            phases.append(rsti)
            traj['t'].append(rsti['t'])
            traj['x'].append(rsti['x'])
            traj['u'].append(rsti['u'])
            if rsti['p'] is not None:
                if traj['p'] is None:
                    traj['p'] = [rsti['p']]
                else:
                    traj['p'].append(rsti['p'])
        # concatenate them and get solution
        for key in traj.keys():
            if traj[key] is not None:
                traj[key] = np.concatenate(traj[key], axis=0)
        traj['phases'] = phases
        return traj

    def add_obj(self, obj):
        """Add a objective function to the problem.
        
        :param obj: a general objective function object.
        
        """
        if isinstance(obj, nonLinearObj):
            self.add_nonlinear_obj(obj)
        elif isinstance(obj, (nonLinearPointObj, linearPointObj)):
            self.add_addx_obj(obj)
        else:
            raise NotImplementedError

    def add_constr(self, constr):
        """Add a constraint to the problem.
        
        :param constr: a general constraint function object.

        """
        
        if isinstance(constr, nonLinearObj):
            self.add_nonlinear_constr(constr)
        elif isinstance(constr, (nonLinearPointConstr, linearPointConstr)):
            self.add_addx_constr(constr)
        elif isinstance(constr, (LinearConnectConstr, NonLinearConnectConstr)):
            self.add_connect_constr(constr)
        else:
            raise NotImplementedError

    def add_connect_constr(self, constr):
        """Add a linear constraint that connects two phases.

        :param constr: a connect constraint object.

        """
        if isinstance(constr, LinearConnectConstr):
            self.connect_linear_constr.append(constr)
        elif isinstance(constr, NonLinearConnectConstr):
            self.connect_nonlinear_constr.append(constr)
        else:
            raise NotImplementedError

    def add_connect_linear_constr(self, constr):
        """Add a linear connection constraint to the problem.
        
        :param constr: a LinearConnectConstr object.
        
        """
        assert isinstance(constr, LinearConnectConstr)
        self.connect_linear_constr.append(constr)

    def add_connect_nonlinear_constr(self, constr):
        """Add a nonlinear connect constraint to the problem.
        
        :param constr: a NonLinearConnectConstr object.
        
        """
        assert isinstance(constr, NonLinearConnectConstr)
        self.connect_nonlinear_constr.append(constr)

    def add_addx_constr(self, constr):
        """Add a linear or nonlinear constraint associated with addx.
        
        :param constr: a point constraint.

        """
        if isinstance(constr, linearPointConstr):
            self.addx_linear_constr.append(constr)
        elif isinstance(constr, nonLinearPointConstr):
            self.addx_nonlinear_constr.append(constr)
        else:
            raise NotImplementedError

    def add_nonlinear_constr(self, constr):
        """Add a general nonlinear constraint to the problem.
        
        :param constr: a nonLinearConstr object.

        """
        assert isinstance(constr, nonLinearConstr)
        self.nonlinear_constr.append(constr)

    def add_addx_obj(self, obj):
        """Add an objective evaluated at addx to the problem.
        
        :param obj: a point objective object

        """
        if isinstance(obj, linearPointObj):
            self.addx_linear_obj.append(obj)
        if isinstance(obj, nonLinearPointObj):
            self.addx_nonlinear_obj.append(obj)
        else:
            raise NotImplementedError

    def add_nonlinear_obj(self, obj):
        """Add a nonlinear objective to the problem.
        
        :param obj: a nonLinearObj object.
        
        """
        assert isinstance(obj, nonLinearObj)
        self.nonlinear_obj.append(obj)

    def change_connect_time_constr_bound(self, num, xl, xu):
        """Change the connect time constraint if other specification are required."""
        assert isinstance(num, int)
        assert isinstance(xl, (int, float)) and isinstance(xu, (int, float))
        if num >= len(self.phases):
            print('Possible range of num is 0 - %d' % (len(self.phases) - 1))
            return
        if xl > xu:
            print('xl should be smaller than xu')
            return
        self.connect_linear_constr[num].lb[:] = xl
        self.connect_linear_constr[num].ub[:] = xu

    def __parseAddX(self, addx):
        """Parse addx"""
        if addx is None:
            self.len_addX = 0
        else:
            if not isinstance(addx, list):
                addx = [addx]
            for tmp in addx:
                assert isinstance(tmp, addX)
            self.addX = addx
            vec_len_addX = [tmp.n for tmp in addx]
            self.len_addX = np.sum(vec_len_addX)  # total length of addx
            self.accum_len_addX = np.insert(np.cumsum(vec_len_addX), 0, 0)

    def __parseAddX__(self, x):
        """Return a list of addX values."""
        numTraj = self.__get_addx_leading_column(0)
        addX = []
        for addx in self.addX:
            addX.append(x[numTraj: numTraj + addx.n])
            numTraj += addx.n
        return addX

    def __get_leading_column(self, phase, index):
        """For the vector in phase at index, return column index of first element."""
        if index >= 0:
            return self.accum_num_sol[phase] + index * self.phases[phase].dimpoint
        else:
            return self.accum_num_sol[phase] + (self.phases[phase].nPoint + index) * self.phases[phase].dimpoint

    def __get_t0_ind(self, phase):
        """Return the index of t0 for selected phase."""
        return self.accum_num_sol[phase] + self.phases[phase].t0ind

    def __get_tf_ind(self, phase):
        """Return the index of t0 for selected phase."""
        return self.accum_num_sol[phase] + self.phases[phase].tfind

    def __get_addx_leading_column(self, addxindex):
        """Return the column index of selected addx."""
        return self.accum_len_addX[addxindex] + self.accum_num_sol[-1]

    def __get_phase_sol_piece(self, x, phase_num):
        """Return the long vector belonging to certain phase."""
        return x[self.accum_num_sol[phase_num]: self.accum_num_sol[phase_num + 1]]

    def __get_phase_f_piece(self, f, phase_num):
        """Return the piece that collects constraint function in certain phases"""
        n1 = self.accum_num_f[phase_num]
        n2 = self.accum_num_f[phase_num + 1]
        return f[1 + n1: 1 + n2]

    def __get_addx_by_index(self, x, index):
        """Return the piece of addx."""
        i0 = self.__get_addx_leading_column(index)
        return x[i0: i0 + self.addX[index].n]

    def __extract_func_arguments(self, x, constr):
        """Extract arguments for connect_nonlinear_constr"""
        phase1, phase2 = constr.phases
        index1, index2 = constr.indexes
        x1, x2 = self.__get_phase_sol_piece(x, phase1), self.__get_phase_sol_piece(x, phase2)
        h1, useT1 = self.phases[phase1].__get_time_grid__(x1)
        useX1, useU1, useP1 = self.phases[phase1].__parseX__(x1)
        h2, useT2 = self.phases[phase2].__get_time_grid__(x2)
        useX2, useU2, useP2 = self.phases[phase2].__parseX__(x2)
        if constr.addx_index is not None:
            n0 = self.__get_addx_leading_column(constr.addx_index)
            addx = x[n0: n0 + self.addX[constr.addx_index].n]
        else:
            addx = None
        x1 = np.concatenate(([[useT1[index1]], useX1[index1], useU1[index1], useP1[index1]]))
        x2 = np.concatenate(([[useT2[index2]], useX2[index2], useU2[index2], useP2[index2]]))
        return x1, x2, addx

    def __add_connect_time_constr(self):
        """Add linear constraints that limits t0 and tf of two connecting phases.

        This is done by automatically add connect linear constraints.
        It basically requires the continuity of time.

        """
        num_phase = len(self.phases)
        for i in range(num_phase - 1):
            a1 = np.ones(1)
            a2 = -a1
            constr = LinearConnectConstr(i, i + 1, a1, a2, lb=np.zeros(1), ub=np.zeros(1))
            self.add_connect_linear_constr(constr)

    def __set_x_bound(self):
        """Set xlb and xub for this new structure. Basically it copies and creates."""
        xlb = -1e20 * np.ones(self.num_sol)
        xub = -xlb
        for phase_num, phase in enumerate(self.phases):
            tmplb = self.__get_phase_sol_piece(xlb, phase_num)
            tmpub = self.__get_phase_sol_piece(xub, phase_num)
            tmplb[:] = phase.xlb
            tmpub[:] = phase.xub
        # set bounds for addx
        cur_col = self.accum_num_sol[-1]
        if self.len_addX > 0:
            for addx in self.addX:
                if addx.lb is None:
                    addx.lb = np.zeros(addx.n)
                if addx.ub is None:
                    addx.ub = np.zeros(addx.n)
                xlb[cur_col: cur_col + addx.n] = addx.lb
                xub[cur_col: cur_col + addx.n] = addx.ub
                cur_col += addx.n
        # I do not care about auxiliary variables since they are unbounded
        self.xlb = xlb
        self.xub = xub

    def __set_f_bound(self):
        """Set lb and ub for constraint functions."""
        lb = np.zeros(self.num_f)
        ub = np.zeros(self.num_f)
        lb[0] = -1e20
        ub[0] = 1e20
        # loop over all prob
        for phase_num, phase in enumerate(self.phases):
            tmplb = self.__get_phase_f_piece(lb, phase_num)
            tmpub = self.__get_phase_f_piece(ub, phase_num)
            tmplb[:] = phase.lb[1:]
            tmpub[:] = phase.ub[1:]
        # loop over linear constraints
        curRow = self.accum_num_f[-1] + 1
        for constr in self.connect_linear_constr:
            lb[curRow: curRow + constr.nf] = constr.lb
            ub[curRow: curRow + constr.nf] = constr.ub
            curRow += constr.nf
        for constr in self.addx_linear_constr:  # actually I can list + list for simplicity
            lb[curRow: curRow + constr.nf] = constr.lb
            ub[curRow: curRow + constr.nf] = constr.ub
            curRow += constr.nf
        # loop over nonlinear constraints
        for constr in self.connect_nonlinear_constr + self.addx_nonlinear_constr + self.nonlinear_constr:
            lb[curRow: curRow + constr.nf] = constr.lb
            ub[curRow: curRow + constr.nf] = constr.ub
            curRow += constr.nf
        # for auxiliary variables the bounds are zero so I am good
        self.lb = lb
        self.ub = ub

    def __set_A_pattern(self, est_num_F, est_num_sol):
        """Set the pattern of A. 
        
        We do three things:
        - analyze the linear constraints, generate a list of sparse-like matrix
        - analyze the additional cost functions such as penalty on addX
        - change row/col of the Aval/Arow/Acol for each phase

        :param est_num_F: estimated number of F (including objective) with all constraints considered
        :param est_num_sol: estimated length of x (including addx) but without auxiliary ones
        """
        # analyze the linear connect constraints
        lstA, lstArow, lstAcol = self.__analyze_linear_constr()
        lstA2, lstArow2, lstAcol2 = self.__analyze_cost_function(est_num_F, est_num_sol)
        lstA3, lstArow3, lstAcol3 = self.__analyze_existing_A()
        self.Aval = np.concatenate(lstA + lstA2 + lstA3)
        self.Arow = np.concatenate(lstArow + lstArow2 + lstArow3)
        self.Acol = np.concatenate(lstAcol + lstAcol2 + lstAcol3)

    def __analyze_linear_constr(self):
        """Detect the sparse A for linear constraints.
        
        returns lstA, lstArow, lstAcol: list of sparse matrix A from linear constraints.
        
        """
        lstA = []
        lstArow = []
        lstAcol = []
        curRow = 1 + self.sum_num_f

        def process_one_matrix(constr, auto1, phase1, index1, a1, tidx1, curRow):
            # for first constraint
            col1 = self.__get_leading_column(phase1, index1)
            if auto1:  # actually this part can be merged
                lstA.append(a1.data)
                lstArow.append(curRow + a1.row)
                lstAcol.append(col1 + a1.col)
            else:
                timemask = np.zeros(a1.nnz, dtype=bool)
                timemask[tidx1] = True
                statemask = np.logical_not(timemask)
                # whatever happens, state related are processed
                lstA.append(a1.data[statemask])
                lstArow.append(a1.row[statemask] + curRow)
                lstAcol.append(col1 + a1.col[statemask] - 1)
                phase = self.phases[phase1]
                if index1 >= 0:
                    t0factor = (phase.N - 1 - index1) / (phase.N - 1)
                    tffactor = (index1) / (phase.N - 1)
                else:
                    t0factor = (-index1 - 1) / (phase.N - 1)
                    tffactor = (phase.N + index1) / (phase.N - 1)
                if phase.t0ind > 0:
                    lstA.append(a1.data[timemask] * t0factor)
                    lstArow.append(a1.row[timemask] + curRow)
                    lstAcol.append(np.sum(timemask) * [self.__get_t0_ind(phase1)])
                else:  # we have constraints on it but it is fixed, we change lb and ub
                    value = -t0factor * a1.data[timemask] * self.phases[phase1].t0
                    constr.lb[a1.row[tidx1]] += value
                    constr.ub[a1.row[tidx1]] += value
                if phase.tfind > 0:
                    lstA.append(a1.data[timemask] * tffactor)
                    lstArow.append(a1.row[timemask] + curRow)
                    lstAcol.append(np.sum(timemask) * [self.__get_tf_ind(phase1)])
                else:
                    value = -tffactor * a1.data[timemask] * self.phases[phase1].tf
                    constr.lb[a1.row[tidx1]] += value
                    constr.ub[a1.row[tidx1]] += value

        def process_adda_matrix(constr):
            if constr.addx_index is not None:
                lstA.append(constr.adda.data)
                lstArow.append(constr.adda.row + curRow)
                lstAcol.append(constr.adda.col + self.__get_addx_leading_column(constr.addx_index))

        for constr in self.connect_linear_constr:
            # constr.find_time_gradient()
            phase1, phase2 = constr.phases
            index1, index2 = constr.indexes
            a1, a2 = constr.a
            tidx1, tidx2 = constr.timeindex
            auto1, auto2 = constr.autonomous
            process_one_matrix(constr, auto1, phase1, index1, a1, tidx1, curRow)
            process_one_matrix(constr, auto2, phase2, index2, a2, tidx2, curRow)
            process_adda_matrix(constr)
            curRow += constr.nf

        for constr in self.addx_linear_constr:
            lstA.append(constr.A.data)
            lstArow.append(constr.A.row + curRow)
            lstAcol.append(constr.A.col + self.__get_addx_leading_column(constr.index))
            curRow += constr.nf

        return lstA, lstArow, lstAcol

    def __analyze_cost_function(self, est_num_F, est_num_sol):
        """Analyze the cost functions. It has two things to do"""
        lstA = []
        lstArow = []
        lstAcol = []
        for phase_num, phase in enumerate(self.phases):
            col_index = phase.Acol[phase.Arow == 0]
            n_col_index = len(col_index)
            lstAcol.append(col_index + self.accum_num_sol[phase_num])
            lstA.append(np.ones(n_col_index))
            lstArow.append(np.zeros(n_col_index))
        # append linear objectives to first row
        if len(self.addx_linear_constr) > 0:
            tmpA = np.zeros(self.len_addX)
            basecol = self.accum_num_sol[-1]
            for obj in self.addx_linear_obj:
                leadcol = np.__get_addx_leading_column(obj.index)
                tmpA[leadcol - basecol + obj.A.col] += obj.A.data
            tmpA_ = coo_matrix(tmpA)
            lstA.append(tmpA_.data)
            lstArow.append(np.zeros(tmpA_.nnz))
            lstAcol.append(tmpA_.col + basecol)
        # check all nonlinear objectives on addx
        num_addx_obj = len(self.addx_nonlinear_obj)  # this is the number of auxiliary variables
        num_non_linear_obj = len(self.nonlinear_obj)
        self.num_aux_obj_var = num_addx_obj + num_non_linear_obj
        # the upper right corner
        lstAcol.append(np.arange(num_addx_obj) + est_num_sol)
        lstArow.append(np.zeros(num_addx_obj))
        lstA.append(np.ones(num_addx_obj))
        # the bottom right corner
        lstAcol.append(np.arange(num_addx_obj) + est_num_sol)
        lstArow.append(est_num_F + np.arange(num_addx_obj))
        lstA.append(-np.ones(num_addx_obj))
        return lstA, lstArow, lstAcol

    def __analyze_existing_A(self):
        """Loop over all phases, change row indexes of them."""
        lstA = []
        lstArow = []
        lstAcol = []
        sumf = 1  # record how many f we have used
        for phase_num, phase in enumerate(self.phases):
            # find non-zero row number
            mask = phase.Arow > 0
            lstA.append(phase.Aval[mask])
            lstArow.append(phase.Arow[mask] - 1 + sumf)
            sumf += phase.numF - 1
            lstAcol.append(phase.Acol[mask] + self.accum_num_sol[phase_num])
        return lstA, lstArow, lstAcol

    def __find_num_G(self, x):
        """Find number of G for this huge problem."""
        nG = 0
        maxnG = 0
        for phase in self.phases:
            nG += phase.nG
        for constr in self.connect_nonlinear_constr:
            nG += constr.nG
            maxnG = max(maxnG, constr.nG)
            x1, x2, addx = self.__extract_func_arguments(x, constr)
            constr.find_time_gradient(x1, x2, addx)
            n_time_related_0 = len(constr.timeindex[0])
            n_time_related_1 = len(constr.timeindex[1])
            nG += (self.phases[constr.phases[0]].numT - 1) * n_time_related_0
            nG += (self.phases[constr.phases[1]].numT - 1) * n_time_related_1
        for obj in self.addx_nonlinear_obj:
            nG += obj.nG
            maxnG = max(maxnG, obj.nG)
        self.maxnG = maxnG
        self.G = np.zeros(maxnG)
        self.row = np.zeros(maxnG, dtype=int)
        self.col = np.zeros(maxnG, dtype=int)
        return nG

    def __calc_nonlinear_obj(self, x, y, G, row, col, curRow, curNg, rec, needg):
        """Calculate the nonlinear functions. 
        
        See :func:`trajOptMultiPhaseCollocationProblem.trajOptMultiPhaseCollocProblem.__calc_nonlinear_constr`
        
        """
        tmpout = np.zeros(1)
        y[0] = 0  # since this row is purely linear
        curRow = self.numF - self.num_aux_obj_var
        # loop over addx_nonlinear_obj, it is a pointObj with index information inside it
        for obj in self.addx_nonlinear_obj:
            idx = obj.index
            i0 = self.__get_addx_leading_column(idx)
            addx = x[i0: i0 + self.addX[idx].n] 
            Gpiece = G[curNg: curNg + obj.nG]
            rowpiece = row[curNg: curNg + obj.nG]
            colpiece = col[curNg: curNg + obj.nG]
            obj.__callg__(self, addx, y[curRow: curRow + 1], Gpiece, rowpiece, colpiece, rec, needg)
            if rec:
                rowpiece[:] = curRow
                colpiece[:] += i0
            curNg += obj.nG
            curRow += 1
        # loop over nonlinear_obj
        for obj in self.nonlinear_obj:
            obj.__callg__(self, x, y[curRow: curRow + 1], G[curNg: curNg + obj.nG], row[curNg: curNg + obj.nG],
                        col[curNg: curNg + obj.nG], rec, needg)
            if rec:
                row[curNg: curNg + obj.nG] = curRow
            curRow += 1
            curNg += obj.nG
        return curRow, curNg  # TODO: add support for 

    def __calc_nonlinear_constr(self, curRow, curNg, x, y, G, row, col, rec, needg):
        """Calculation of nonlinear constraints.
        
        :param x: the long vector of solution
        :param y: where to store constraint function evaluation
        :param G, row, col: the sparse Jacobian is stored here
        :param curRow, curNg: accumulated number of y/G
        :return curRow, curNg: updated number of used y/G
        """
        # loop over connect_nonlinear_constr
        for constr in self.connect_nonlinear_constr:
            # get arguments for calling this function
            x1, x2, addx = self.__extract_func_arguments(x, constr)
            phase1 = self.phases[constr.phases[0]]
            phase2 = self.phases[constr.phases[1]]
            phase_num1, phase_num2 = constr.phases
            index1, index2 = constr.indexes
            offset1 = self.__get_leading_column(phase_num1, index1)
            offset2 = self.__get_leading_column(phase_num2, index2)
            if constr.human:
                rst = constr.__callhumang__(x1, x2, needg, addx)
                tmpy, G1, G2, G3 = rst
                y[curRow: curRow + constr.nf] = tmpy
                if needg:
                    # for phase 1
                    if constr.autonomous[0]:  # we are happy
                        G[curNg: curNg + G1.nnz] = G1.data  # np.concatenate((G1.data, G2.data, G3.data))
                        if rec:
                            row[curNg: curNg + G1.nnz] = np.concatenate((G1.row)) + curRow
                            col[curNg: curNg + G1.nnz] = phase1.__patchCol__(index1, G1.col, offset1)
                        curNg += G1.nnz
                    else:
                        curNg = phase1.__copy_into_g__(index1, G, row, col, curRow, curNg, G1.nnz, constr.timeindex[0], 
                                False, rec, G1.data, G1.row, G1.col, offset1)
                    # for phase 2
                    if constr.autonomous[1]:  # we are happy
                        G[curNg: curNg + G2.nnz] = G2.data
                        if rec:
                            row[curNg: curNg + G2.nnz] = np.concatenate((G2.row)) + curRow
                            col[curNg: curNg + G2.nnz] = phase1.__patchCol__(index2, G2.col, offset2)
                        curNg += G2.nnz
                    else:
                        curNg = phase2.__copy_into_g__(index2, G, row, col, curRow, curNg, G2.nnz, constr.timeindex[1],
                                False, rec, G2.data, G2.row, G2.col, offset2)
                    # for phase a, if applicable
                    if constr.addx_index is not None:
                        G[curNg: curNg + G3.nnz] = G3.data
                        if rec:
                            row[curNg: curNg + G3.nnz] = curRow + G3.row
                            col[curNg: curNg + G3.nnz] = self.__get_addx_leading_column(constr.addx_index) + G3.col
            else:
                if constr.autonomous[0] and constr.autonomous[1]:  # we are sure it is within bound and we are fine
                    Gpiece = G[curNg: curNg + constr.nG]
                    rowpiece = row[curNg: curNg + constr.nG]
                    colpiece = col[curNg: curNg + constr.nG]
                    tmpcol = self.col[:constr.nG]
                    constr.__callg__(x1, x2, y[curRow: curRow + constr.nf], Gpiece, rowpiece, tmpcol, rec, needg, addx)
                    if rec:
                        rowpiece[:] += curRow
                        # assign self.col to correct place
                        index1_ = np.where(tmpcol < len(x1))[0]
                        colpiece[index1_] = phase1.__patchCol__(0, tmpcol[index1_], offset1)  # why use 0? Because offset1 accounts for them
                        index2_ = np.where((tmpcol >= len(x1)) & (tmpcol < len(x1) + len(x2)))[0]
                        colpiece[index2_] = phase2.__patchCol__(0, tmpcol[index2_] - len(x1), offset2)
                        if constr.addx_index is not None:
                            index3_ = np.where(tmpcol >= len(x1) + len(x2))[0]
                            colpiece[index3_] = tmpcol[index3_] - len(x1) - len(x2) + self.__get_addx_leading_column(constr.addx_index)
                    curNg += constr.nG
                else:
                    constr.__callg__(x1, x2, y[curRow: curRow + constr.nf], self.G, self.row, self.col, rec, needg, addx)
                    mask1 = self.col < len(x1)
                    ng1 = np.sum(mask1)
                    curNg = phase1.__copy_into_g__(index1, G, row, col, curRow, curNg, ng1, constr.timeindex[0],
                            False, rec, self.G[mask1], self.row[mask1], self.col[mask1], offset1)
                    mask2 = (self.col >= len(x1)) & (self.col < len(x1) + len(x2))
                    ng2 = np.sum(mask2)
                    curNg = phase2.__copy_into_g__(index2, G, row, col, curRow, curNg, ng2, constr.timeindex[1],
                            False, rec, self.G[mask2], self.row[mask2], self.col[mask2] - len(x1), offset2)
                    if constr.timeindex is not None:
                        mask3 = self.col >= len(x1) + len(x2)
                        ng3 = np.sum(mask3)
                        G[curNg: curNg + ng3] = self.G[mask3]
                        if rec:
                            row[curNg: curNg + ng3] = self.row[mask3] + curRow
                            col[curNg: curNg + ng3] = self.col[mask3] - len(x1) - len(x2) + self.__get_addx_leading_column(constr.addx_index)
                        curNg += ng3
            curRow += constr.nf
        for constr in self.addx_nonlinear_constr:
            usex = self.__get_addx_by_index(x, constr.index)
            constr.__callg__(usex, y[curRow: curRow + constr.nf], G[curNg: curNg + constr.nG],
                            row[curNg: curNg + constr.nG], col[curNg: curNg + constr.nG], rec, needg)
            if rec:
                row[curNg: curNg + constr.nG] += curRow
                col[curNg: curNg + constr.nG] += self.__get_addx_leading_column(constr.index)
            curRow += constr.nf
            curNg += constr.nG
        for constr in self.nonlinear_constr:
            constr.__callg__(x, y[curRow: curRow + constr.nf], G[curNg: curNg + constr.nG], 
                    row[curNg: curNg + constr.nG], col[curNg: curNg + constr.nG], rec, needg)
            if rec:
                row[curNg: curNg + constr.nG] += curRow
            curRow += constr.nf
            curNg += constr.nG
        return curRow, curNg

    def __checkGError(self, G, row, col, curNg):
        """For debug purpose. Check if G and A overlap. Due to incorrect assignment of Jacobian row and column numbers.
        It uses bool operation.

        :param G, row, col: the Jacobian matrix to check
        :param curNg: accumulated G.
        """
        boolMat1 = np.zeros((self.nf, self.nx), dtype=bool)
        boolMat1[self.Arow, self.Acol] = True
        boolMat2 = np.zeros((self.nf, self.nx), dtype=bool)
        boolMat2[row[:curNg], col[:curNg]] = True
        # find both true
        badMat = boolMat1 & boolMat2
        nRepeat = np.sum(badMat)
        print('%d %d %d' % (np.sum(boolMat1), np.sum(boolMat2), nRepeat))
        if nRepeat > 0:
            print(np.sum(badMat))
        pass
