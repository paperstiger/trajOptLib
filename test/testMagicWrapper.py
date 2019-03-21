#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testMagicWrapper.py

We test if autograd can be used and use the toy problem to demonstrate so.

"""
import sys, os, time
import numpy as np
import autograd.numpy as numpy
import matplotlib.pyplot as plt
import logging
import autograd.numpy as numpy
from trajoptlib import DaeSystemWrapper
from trajoptlib import DaeSystem, TrajOptCollocProblem
from trajoptlib import LqrObj
from trajoptlib.utility import show_sol
from trajoptlib import OptConfig, OptSolver


def sysFunOrder2(t, X, u, p):
    x, dx, ddx = X
    return numpy.array([ddx - u[0]])


def sysFunOrder1(t, X, u, p):
    x, v, dx, dv = X
    y0 = dx - v
    y1 = dv - u[0]
    return numpy.array([y0, y1])


def main():
    # prob = constructOrderOne()
    prob = constructOrderTwo()
    # construct a solver for the problem
    cfg = OptConfig(print_level=5)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parse_sol(rst.sol)
        show_sol(sol)


def constructOrderOne():
    """Test the wrapper class for this naive problem"""
    sys = DaeSystemWrapper(sysFunOrder1, 4, 1, 0, 2)
    N = 20
    t0 = 0.0
    tf = 10.0
    prob = TrajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20, -1e20]), np.array([0, 0, 1e20, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20, -1e20]), np.array([np.pi, 0, 1e20, 1e20])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process()  # construct the problem
    return prob


def constructOrderTwo():
    """Test the wrapper class for yet another naive problem."""
    sys = DaeSystemWrapper(sysFunOrder2, 3, 1, 0, 1)
    N = 20
    t0 = 0.0
    tf = 10.0
    prob = TrajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20]), np.array([np.pi, 0, 1e20])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process()  # construct the problem
    return prob


if __name__ == '__main__':
    main()
