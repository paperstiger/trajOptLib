#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testManifoldColloc.py

Test the collocation approach with manifold support.
Basically I will define a second-order omni-directional car.
It is constrained to move on a circle. We move from (0, 0) to (1, 1) on a circle (x-0.5)^2+(y-0.5)^2=0.5
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
from trajoptlib import TrajOptManifoldCollocProblem, LqrObj, OptConfig, OptSolver, show_sol
from scipy.sparse import coo_matrix

from carCommon import CircleConstr, SecondOrderOmniCar


RX = 0
RY = 0
RADIUS = 1


def main():
    testAnalytic()


def testAnalytic():
    """Test car problem with analytic solution. See what's missing."""
    sys = SecondOrderOmniCar()  # simple car has no p
    N = 10
    t0 = 0
    tf = 1
    man_constr = CircleConstr()
    prob = TrajOptManifoldCollocProblem(sys, N, t0, tf, man_constr)
    setBound(prob)

    lqr = LqrObj(R=np.ones(2), P=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process(defect_u=True, defect_p=False)
    # construct solver
    cfg = OptConfig()
    slv = OptSolver(prob, cfg)
    # generate an guess
    x0 = np.zeros(prob.nx)
    # generate the random initial guess
    t = np.linspace(0, 1, 2*N - 1)
    theta = -np.pi * t ** 3 + 1.5 * np.pi * t ** 2
    dtheta = -3*np.pi*t**2 + 3*np.pi*t
    ddtheta = -6*np.pi*t + 3*np.pi
    # assign them to it
    useX, useU, useP = prob.__parseX__(x0)
    radius = 1
    stheta = np.sin(theta)
    ctheta = np.cos(theta)
    useX[:, 0] = radius * ctheta
    useX[:, 1] = radius * stheta
    useX[:, 2] = -radius * stheta * dtheta
    useX[:, 3] = radius * ctheta * dtheta
    useX[:, 4] = -radius * ctheta * dtheta ** 2 - radius * stheta * ddtheta
    useX[:, 5] = -radius * stheta * dtheta ** 2 + radius * ctheta * ddtheta
    useU[:, 0] = useX[:, 4]
    useU[:, 1] = useX[:, 5]
    # gamma is assumed zero
    psf0 = prob.parse_f(x0)
    rst = slv.solve_guess(x0)
    psf = prob.parse_f(rst.sol)
    # solve the problem
    rst = slv.solve_guess(x0)
    print(rst.flag)
    # parse the solution
    sol = prob.parse_sol(rst.sol)
    show_sol(sol)


def setBound(prob):
    x0lb = np.array([1, 0, 0, 0, -1e20, -1e20])
    x0ub = np.array([1, 0, 0, 0, 1e20, 1e20])
    prob.x0bd = [x0lb, x0ub]
    xflb = np.array([0, 1, 0, 0, -1e20, -1e20])
    xfub = np.array([0, 1, 0, 0, 1e20, 1e20])
    prob.xfbd = [xflb, xfub]
    prob.xbd = [-1e20*np.ones(6), 1e20*np.ones(6)]
    prob.ubd = [-1e20*np.ones(2), 1e20*np.ones(2)]
    if prob.numP > 0:
        prob.pbd = [-1e20*np.ones(1), 1e20*np.ones(1)]


if __name__ == '__main__':
    main()
