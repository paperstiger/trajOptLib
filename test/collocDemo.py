#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
collocDemo.py

Use the collocation version of problem.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
from trajoptlib.io import get_onoff_args
from trajoptlib import DaeSystem, TrajOptCollocProblem
from trajoptlib import NonLinearPointObj, LinearPointObj, LinearPointConstr
from trajoptlib import LqrObj
from trajoptlib import OptConfig, OptSolver
from trajoptlib.utility import show_sol
from scipy.sparse import coo_matrix


class OneDcase(DaeSystem):
    def __init__(self):
        DaeSystem.__init__(self, 3, 1, 0, 1, 2)  # ddx = u

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[2] - u[0]
        if needg:
            G[0] = 1
            G[1] = -1
            if rec:
                row[0] = 0
                row[1] = 0
                col[0] = 3
                col[1] = 4


class OrderOneOneD(DaeSystem):
    def __init__(self):
        DaeSystem.__init__(self, 4, 1, 0, 2, 4)  # ddx = u

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[2] - x[1]
        y[1] = x[3] - u[0]
        if needg:
            G[0] = 1
            G[1] = -1
            G[2] = 1
            G[3] = -1
            if rec:
                row[:2] = 0
                row[2:] = 1
                col[:] = [3, 2, 4, 5]


class Pendulum(DaeSystem):
    """Test pendulum nonlinearity."""
    def __init__(self):
        DaeSystem.__init__(self, 3, 1, 0, 1, 3)  # ddq = u/5 - .5*sin(q)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[2] - u[0] / 5. + 0.5 * np.sin(x[0])
        if needg:
            G[0] = 0.5 * np.cos(x[0])
            G[1] = 1
            G[2] = -0.2
            if rec:
                row[:3] = 0
                col[:3] = [1, 3, 4]


class OrderOnePendulum(DaeSystem):
    """Pendulum with order 1"""
    def __init__(self):
        DaeSystem.__init__(self, 4, 1, 0, 2, 5)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        theta, omega, dtheta, domega = x
        y[0] = dtheta - omega
        y[1] = domega - u[0] / 5.0 + 0.5 * np.sin(x[0])
        if needg:
            G[0] = 1
            G[1] = -1
            G[2] = 1
            G[3] = 0.5 * np.cos(x[0])
            G[4] = -0.2
            if rec:
                row[:2] = 0
                row[2:] = 1
                col[:] = [3, 2, 4, 1, 5]


def main():
    args = get_onoff_args('oned', 'pen', 'lqr', 'linear', 'orderone', 'backend ipopt')
    if args.oned:
        testOneD()
    if args.pen:
        testPen()
    if args.linear:
        testLinear()
    if args.orderone:
        testOrderOne(args)


def testOrderOne(args):
    """Test order one pendulum case, this is seen everywhere."""
    if args.pen:
        sys = OrderOnePendulum()
    else:
        sys = OrderOneOneD()
    N = 20
    t0 = 0.0
    tf = 20.0
    prob = TrajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20, -1e20]), np.array([0, 0, 1e20, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20, -1e20]), np.array([np.pi, 0, 1e20, 1e20])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(backend=args.backend)
    solver = OptSolver(prob, cfg)
    rst = solver.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        show_sol(sol)


def testLinear(args):
    """Test 1d problem with linear constraints and linear objective"""
    sys = OneDcase()
    N = 10
    t0 = 0.0
    tf = 2.0
    prob = TrajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([1, 0, -1e20]), np.array([1, 0, 1e20])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    A = np.zeros(5)
    A[1] = 1
    A[2] = 1  # so it basically does nothing
    linPntObj = LinearPointObj(0, A, 3, 1, 0)
    prob.add_obj(linPntObj)
    # add linear constraint that x is increasing
    A = np.zeros(5)
    A[1] = 1
    lb = np.zeros(1)
    ub = np.ones(1)
    linPntCon = LinearPointConstr(-1, A, lb, ub)
    prob.add_constr(linPntCon, True)
    # we want mid point to be close to 0.8
    wantState = np.array([0.8, 0])
    pntObj = PointObj(N, wantState)
    prob.addObj(pntObj)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(args.backend)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag, rst.sol)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


class PointObj(NonLinearPointObj):
    """A objective function to make mid point close to a selected point"""
    def __init__(self, N, state):
        NonLinearPointObj.__init__(self, 15, 3, 1, 0, 'user', 2)
        self.state = state
        self.weight = 100

    def __callg__(self, x, F, G, row, col, rec, needg):
        dx = x[1:3] - self.state
        F[0] = self.weight * np.sum(dx ** 2)
        if needg:
            G[:2] = self.weight * 2 * dx
            if rec:
                row[:2] = 0
                col[:2] = np.arange(1, 3)


def testOneD(args):
    """Test solving one-dim problem using collocation approach"""
    sys = OneDcase()
    N = 10
    t0 = [-1.0, 0]
    tf = [2.0, 3.0]
    prob = TrajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([1, 0, -1e20]), np.array([1, 0, 1e20])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(args.backend)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag, rst.sol)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


def testPen(args):
    """Test solving pendulum swing up problem using collocation approach"""
    sys = Pendulum()
    N = 20
    t0 = 0.0
    tf = 20.0
    prob = TrajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20]), np.array([np.pi, 0, 1e20])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(args.backend)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


if __name__ == '__main__':
    main()
