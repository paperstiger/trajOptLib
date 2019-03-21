#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
quadDemo.py

A demo of quadrotor which is quite challenging.
Use the classical model used everywhere else.
"""
from math import sin, cos
import numpy as np
from trajoptlib import System, NonLinearObj, TrajOptProblem, LqrObj, OptConfig, OptSolver
from trajoptlib.utility import show_sol
from trajoptlib.io import get_onoff_args
from scipy.sparse import coo_matrix


class Rotor:
    def __init__(self):
        self.dimx = 12
        self.dimu = 4
        self.m = 0.5
        self.g = 9.81
        self.kF = 1
        self.kM = 0.0245
        self.L = 0.175
        self.In = [0.0023, 0.0023, 0.004]
        self.cfg0 = np.zeros(204)

    def dyn(self, t, x, u, f, df):
        self.drone_dyn(t, x, u)
        df[:] = 0
        f[:] = self.cfg0[:self.dimx]
        df[:, 1: 1 + self.dimx + self.dimu] = self.cfg0[self.dimx:].reshape((self.dimx, self.dimx + self.dimu), order='F')

    def drone_dyn(self, t, x, u):
        m, g, kF, kM, L, In, cg0 = self.m, self.g, self.kF, self.kM, self.L, self.In, self.cfg0
        phi = x[3]; theta = x[4]; psi = x[5]; xd = x[6]; yd = x[7]; zd = x[8]; p = x[9]; q = x[10]; r = x[11]
        t1 = cos(theta)
        t2 = sin(theta)
        t3 = p * t1 + r * t2
        t4 = sin(phi)
        t5 = cos(phi)
        t6 = 0.1e1 / t5
        t7 = t1 * r
        t8 = t2 * p
        t9 = t8 - t7
        t10 = t6 * t9
        t11 = cos(psi)
        t12 = sin(psi)
        t13 = t1 * t12
        t14 = t11 * t2
        t15 = (u[0] + u[1] + u[2] + u[3]) * kF
        t11 = t11 * t1
        t12 = t12 * t2
        t16 = -t11 * t4 + t12
        t17 = 0.1e1 / m
        t5 = t17 * t5
        t18 = t5 * t1
        t19 = -In[1] + In[2]
        t20 = q * t19
        t21 = In[0] - In[2]
        t22 = p * t21
        t23 = L * kF * (u[0] - u[2]) + t22 * r
        t24 = In[0] - In[1]
        t25 = p * t24
        t26 = (u[0] - u[1] + u[2] - u[3]) * kM + t25 * q
        t27 = t6 ** 0.2e1
        t28 = t27 * t4 ** 0.2e1 + 0.1e1
        t7 = -t7 * t28 + t8 * t28
        t8 = t6 * t3
        t28 = 0.1e1 / In[1]
        t29 = 0.1e1 / In[0]
        t30 = 0.1e1 / In[2]
        t31 = t18 * kF
        t32 = t17 * (t13 * t4 + t14)
        t33 = t32 * t15
        t32 = t32 * kF
        t34 = t17 * t16
        t35 = t34 * kF
        t36 = t28 * L * kF
        t37 = t29 * L * kF
        t38 = t30 * kM
        t39 = t1 * t6
        t6 = t2 * t6
        cg0[0] = xd
        cg0[1] = yd
        cg0[2] = zd
        cg0[3] = t3
        cg0[4] = t10 * t4 + q
        cg0[5] = -t10
        cg0[6] = t33
        cg0[7] = t34 * t15
        cg0[8] = t18 * t15 - g
        cg0[9] = -t29 * (-L * kF * (u[1] - u[3]) + t20 * r)
        cg0[10] = -t28 * t23
        cg0[11] = t30 * t26
        cg0[52] = t7
        cg0[53] = -t27 * t4 * t9
        cg0[54] = t5 * t13 * t15
        cg0[55] = -t5 * t11 * t15
        cg0[56] = -t17 * t4 * t1 * t15
        cg0[63] = -t9
        cg0[64] = t8 * t4
        cg0[65] = -t8
        cg0[66] = t17 * (-t12 * t4 + t11) * t15
        cg0[67] = t17 * (t14 * t4 + t13) * t15
        cg0[68] = -t5 * t2 * t15
        cg0[78] = -t17 * t16 * t15
        cg0[79] = t33
        cg0[84] = 1
        cg0[97] = 1
        cg0[110] = 1
        cg0[123] = t1
        cg0[124] = t6 * t4
        cg0[125] = -t6
        cg0[130] = -t28 * r * t21
        cg0[131] = t30 * q * t24
        cg0[136] = 1
        cg0[141] = -t29 * r * t19
        cg0[143] = t30 * t25
        cg0[147] = t2
        cg0[148] = -t39 * t4
        cg0[149] = t39
        cg0[153] = -t29 * t20
        cg0[154] = -t28 * t22
        cg0[162] = t32
        cg0[163] = t35
        cg0[164] = t31
        cg0[166] = -t36
        cg0[167] = t38
        cg0[174] = t32
        cg0[175] = t35
        cg0[176] = t31
        cg0[177] = t37
        cg0[179] = -t38
        cg0[186] = t32
        cg0[187] = t35
        cg0[188] = t31
        cg0[190] = t36
        cg0[191] = t38
        cg0[198] = t32
        cg0[199] = t35
        cg0[200] = t31
        cg0[201] = -t37
        cg0[203] = -t38


class QuadRotor(System, Rotor):
    """A class derived from system and Rotor"""
    def __init__(self):
        System.__init__(self, 12, 4, 0, 'Euler')
        Rotor.__init__(self)

    def jac_dyn(self, t, x, u, p=None):
        f = np.zeros(self.nx)
        J = np.zeros((self.nx, self.nx + self.nu + 1 + self.np), order='F')
        Rotor.dyn(self, t, x, u, f, J)
        J = np.ascontiguousarray(J)
        return f, J


class QuadCost(NonLinearObj):
    """A quadratic cost on control."""
    def __init__(self, N, dimx, dimu):
        lenSol = N * (dimx + dimu)
        NonLinearObj.__init__(self, lenSol, 'user', nG=N * dimu)
        self.R = 1.0
        self.N = N
        self.dimx = dimx
        self.dimu = dimu

    def __callf__(self, x, y):
        u = x[3]
        y[0] = u * self.R * u

    def __callg__(self, x, y, G, row, col, rec, needg):
        u = np.reshape(x[self.N * self.dimx:], (self.N, self.dimu))
        y[0] = np.sum(u ** 2)
        if needg:
            G[:self.N * self.dimu] = 2.0 * u.flatten()
            if rec:
                row[:self.N * self.dimu] = 0
                col[:self.N * self.dimu] = np.arange(self.N * self.dimx, self.N * (self.dimx + self.dimu))


def main():
    args = get_onoff_args('backend ipopt')
    sys = QuadRotor()
    N = 40
    dimx, dimu = sys.nx, sys.nu
    cost = QuadCost(N, sys.nx, sys.nu)
    t0 = 0.0
    tf = 5.0
    prob = TrajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [-1e20 * np.ones(sys.nx), 1e20 * np.ones(sys.nx)]
    prob.ubd = [0 * np.ones(sys.nu), 4 * np.ones(sys.nu)]
    prob.x0bd = [np.zeros(sys.nx), np.zeros(sys.nx)]
    prob.xfbd = [np.zeros(sys.nx), np.zeros(sys.nx)]
    prob.xfbd[0][:3] = 5
    prob.xfbd[1][:3] = 5
    if False:
        prob.add_obj(cost)
    else:
        lqr = LqrObj(R=np.ones(4))
        prob.add_lqr_obj(lqr)
    prob.pre_process()
    # construct a solver for the problem
    cfg = OptConfig(args.backend, print_level=5)
    slv = OptSolver(prob, cfg)
    guessx = np.zeros(prob.nx)
    straightx = np.reshape(guessx[:N * dimx], (N, dimx))
    for i in range(3):
        straightx[:, i] = np.linspace(0, prob.xfbd[0][i], N)
    guessx[N * dimx:-1] = np.random.random(N * dimu)
    rst = slv.solve_guess(guessx)
    print(rst.flag)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


if __name__ == '__main__':
    main()
