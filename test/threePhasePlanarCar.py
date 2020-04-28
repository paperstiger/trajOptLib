#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
The classical planar car problem composed of three segments.
And we solve a time optimal problem.
"""
from __future__ import print_function
import numpy as np
from trajoptlib import DaeSystem, TrajOptCollocProblem, LinearPointConstr, LinearObj
from trajoptlib import LqrObj
from trajoptlib import OptConfig, OptSolver
from trajoptlib import TrajOptMultiPhaseCollocProblem, LinearConnectConstr
from trajoptlib.utility import show_sol


class OrderTwoModel(DaeSystem):
    r"""A class with dynamics :math:`\ddot{x}=u; \ddot{y}=v`"""
    def __init__(self):
        DaeSystem.__init__(self, 6, 2, 0, 2, 4)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[4] - u[0]
        y[1] = x[5] - u[1]
        if needg:
            G[0:2] = 1
            G[2:4] = -1
            if rec:
                row[:] = [0, 1, 0, 1]
                col[:] = [5, 6, 7, 8]


class ConnectConstr(LinearConnectConstr):
    """Linear constraints on second order continuity"""
    def __init__(self, n, dimq):
        a1 = np.eye(1 + dimq)
        a2 = -np.eye(1 + dimq)
        LinearConnectConstr.__init__(self, n, n + 1, a1, a2)


class FirstPointConstr(LinearPointConstr):
    """Limit bounds on some variables"""
    def __init__(self, index, lb, ub):
        LinearPointConstr.__init__(self, index, np.eye(2), lb, ub)


def main():
    sys = OrderTwoModel()
    N = 10
    t0 = 0.0
    tf = 3.0
    prob1 = TrajOptCollocProblem(sys, N, t0, [0.1, tf])  # maybe I need to give tips on choosing times
    prob2 = TrajOptCollocProblem(sys, N, [0.1, tf], [0.11, tf])
    prob3 = TrajOptCollocProblem(sys, N, [0.2, tf], [0.21, tf])
    xlb = -1e20 * np.ones(6)
    xub = 1e20 * np.ones(6)
    ulb = -np.ones(2)
    uub = -ulb
    x0, xf = np.array([[0., 0.], [1.0, 0.5]])
    x0lb = np.concatenate((x0, -1e20 * np.ones(4)))
    x0ub = np.concatenate((x0, 1e20 * np.ones(4)))

    xflb = np.concatenate((xf, -1e20 * np.ones(4)))
    xfub = np.concatenate((xf, 1e20 * np.ones(4)))
    prob1.xbd = [xlb, xub]
    prob1.ubd = [ulb, uub]
    prob1.x0bd = [x0lb, x0ub]
    prob1.xfbd = [xlb, xub]
    prob2.xbd = [xlb, xub]
    prob2.ubd = [ulb, uub]
    prob2.x0bd = [xlb, xub]
    prob2.xfbd = [xlb, xub]
    prob3.xbd = [xlb, xub]
    prob3.ubd = [ulb, uub]
    prob3.x0bd = [xlb, xub]
    prob3.xfbd = [xflb, xfub]

    # set bounds constraints for prob1 and prob2 at final
    a = np.zeros((2, 9))  # 1 + 6 + 2
    np.fill_diagonal(a[:, 1:3], 1.0)
    prob1.add_constr(LinearPointConstr(-1, a, np.array([0.2, 0.2]), np.array([0.2, 1e20])))
    prob2.add_constr(LinearPointConstr(-1, a, np.array([0.8, -1e20]), np.array([0.8, 0.3])))

    # define objective function, #TODO: change to time optimal
    obj = LqrObj(R=0.01 * np.ones(2))
    prob1.add_obj(obj, path=True)
    prob2.add_obj(obj, path=True)
    prob3.add_obj(obj, path=True)
    # optimize time
    obj_a = np.zeros(prob3.tfind + 1)
    obj_a[-1] = 1.0
    prob3.add_obj(LinearObj(obj_a))
    # add constraints to some phases
    prob = TrajOptMultiPhaseCollocProblem([prob1, prob2, prob3], addx=None)
    constr1 = ConnectConstr(0, 6)
    constr2 = ConnectConstr(1, 6)
    prob.add_connect_constr(constr1)
    prob.add_connect_constr(constr2)

    # ready to solve
    prob.pre_process()
    # cfg = OptConfig(backend='snopt', deriv_check=1, print_file='tmp.out')
    cfg = OptConfig(print_level=5)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        sol = prob.parse_sol(rst)
        show_sol(sol)


if __name__ == '__main__':
    main()
