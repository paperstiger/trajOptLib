#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
 
"""
A cartpole problem. Contributed by my colleague.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
from trajoptlib.io import get_onoff_args
from trajoptlib import System, NonLinearPointObj, LqrObj
from trajoptlib import TrajOptProblem
from trajoptlib import OptConfig, OptSolver
from trajoptlib.utility import show_sol
from scipy.sparse import coo_matrix
import autograd.numpy as numpy


class CartPole(System):
    """Test pendulum nonlinearity."""
    def __init__(self):
        System.__init__(self, 4, 1, 0, 'Euler')
        self.gravity = 9.81
        self.mass_wheel = 0.42972
        self.mass_shin = 0.66162
        self.mass_thigh = 0.9061
        self.mass_body = 15.513/2
        self.length = 0.5
        self.dt = 0.002
        self.mass_cart = self.mass_wheel
        self.mass_pole = self.mass_shin + self.mass_thigh + self.mass_body
        self.mass_total = self.mass_cart + self.mass_pole

    def __ad__(self, txup):
        x = txup[1: 5]
        u = txup[5: 6]
        y1 = txup[3]
        y2 = txup[4]
        y3 = 1.0/(self.mass_cart+self.mass_pole*numpy.sin(x[1])**2)*(u[0]+self.mass_pole*numpy.sin(x[1])*(self.length*x[3]**2+self.gravity*numpy.cos(x[1])))
        y4 = 1.0/(self.length*(self.mass_cart+self.mass_pole*numpy.sin(x[1])**2))*(-u[0]*numpy.cos(x[1])-self.mass_pole*self.length*x[3]**2*numpy.cos(x[1])*numpy.sin(x[1])-(self.mass_cart+self.mass_pole)*self.gravity*numpy.sin(x[1]))
        return numpy.array([y1, y2, y3, y4])

    def dyn(self, t, x, u, p=None):
        y1 = x[2]
        y2 = x[3]
        y3 = 1.0/(self.mass_cart+self.mass_pole*np.sin(x[1])**2)*(u[0]+self.mass_pole*np.sin(x[1])*(self.length*x[3]**2+self.gravity*np.cos(x[1])))
        y4 = 1.0/(self.length*(self.mass_cart+self.mass_pole*np.sin(x[1])**2))*(-u[0]*np.cos(x[1])-self.mass_pole*self.length*x[3]**2*np.cos(x[1])*np.sin(x[1])-(self.mass_cart+self.mass_pole)*self.gravity*np.sin(x[1]))
        return np.array([y1, y2, y3, y4])

    # def jac_dyn(self, t, x, u, p=None):
    #     y1 = x[2]
    #     y2 = x[3]
    #     y3 = 1.0/(self.mass_cart+self.mass_pole*np.sin(x[1])**2)*(u[0]+self.mass_pole*np.sin(x[1])*(self.length*x[3]**2+self.gravity*np.cos(x[1])))
    #     y4 = 1.0/(self.length*(self.mass_cart+self.mass_pole*np.sin(x[1])**2))*(-u[0]*np.cos(x[1])-self.mass_pole*self.length*x[3]**2*np.cos(x[1])*np.sin(x[1])-(self.mass_cart+self.mass_pole)*self.gravity*np.sin(x[1]))
    #     y = np.array([y1, y2, y3, y4])

    #     J = np.zeros((4, 6))
    #     a1 = self.mass_pole * self.gravity / self.mass_cart
    #     a2 = self.mass_total * self.gravity / (self.length * self.mass_cart)
    #     J = np.array(
    #         [[0, 0, 0, 1, 0, 0], 
    #          [0, 0, 0, 0, 1, 0], 
    #          [0, 0, a1, 0, 0, 1 / self.mass_cart], 
    #          [0, 0, a2, 0, 0, 1 / (self.length * self.mass_cart)]]
    #     )
    #     return y, J


class QuadCost(NonLinearPointObj):
    """A quadratic cost on control."""
    def __init__(self):
        NonLinearPointObj.__init__(self, -1, 4, 1, 0, 'user', nG=1)
        self.R = 1.0

    def __callf__(self, x, y):
        u = x[5]
        y[0] = u * self.R * u
        return 0

    def __callg__(self, x, y, G, row, col, rec, needg):
        u = x[5]
        y[0] = u * self.R * u
        if needg:
            G[0] = 2 * self.R * u
            if rec:
                row[0] = 0
                col[0] = 5
        return (0, 0)


def main():
    args = get_onoff_args('cartpole', 'lqr', 'backend ipopt')
    if args.cartpole:
        gradmode(args)


def gradmode(args):
    """Run pendulum swing-up problem"""
    sys = CartPole()
    cost = QuadCost()
    N = 21
    t0 = 0.0
    tf = 4.0
    prob = TrajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [np.array([-1e20, -1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-100.0]), np.array([100.0])]
    prob.x0bd = [np.array([0, np.pi/8, 0, 0]), np.array([0, np.pi/8, 0, 0])]
    prob.xfbd = [np.array([-1e20, -1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20, 1e20])]
    prob.xfbd = [np.array([0, 0, 0, 0]), np.array([0, 0, 0, 0])]
    if not args.lqr:
        prob.add_obj(cost, True)  # add a path cost
    else:
        lqr = LqrObj(R=np.ones(1), F=10*np.ones(4), Q=np.ones(4))
        prob.add_lqr_obj(lqr)
    snopt_mode = args.backend == 'snopt'

    prob.preProcess(snopt_mode=snopt_mode)  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(args.backend, deriv_check = 0)  #, print_level=1)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        print(sol)
        show_sol(sol)

if __name__ == '__main__':
    main()

