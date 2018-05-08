#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
carCommon.py

Common class files for the car problem.
"""
import numpy as np
from trajOptLib import daeSystem, trajOptCollocProblem
from trajOptLib import lqrObj, nonLinearObj
from trajOptLib import system
from trajOptLib import manifoldConstr, nonLinearPointConstr


class OmniCar(system):
    """A omni-car system in 2D space.

    It has simple dynamics \ddot{x}=u_x; \ddot{y}=u_y with support for constraint force

    """
    def __init__(self):
        system.__init__(self, 4, 2, 1)
        self.rx = 0
        self.ry = 0
        self.radius = 1

    def dyn(self, t, x, u, p, needg=False):
        y = np.zeros(self.nx)
        y[:2] = x[2:]
        y[2] = u[0] + 2*p[0] * (x[0] - self.rx)
        y[3] = u[1] + 2*p[0] * (x[1] - self.ry)
        return y
        if needg:
            G[:] = [1.0, 1.0, 2*p[0], 1.0, 2*(x[0] - self.rx), 2*p[0], 1.0, 2*(x[1] - self.ry)]
            if rec:
                row[:] = [0, 1, 2, 2, 2, 3, 3, 3]
                col[:] = [3, 4, 1, 5, 7, 2, 6, 7]


class FirstOrderOmniCar(daeSystem):
    """A car without constraint. you use control to guarantee it"""
    def __init__(self):
        daeSystem.__init__(self, 8, 2, 1, 4, 12)
        self.rx = 0
        self.ry = 0
        self.radius = 1

    def dyn(self, t, X, u, p, Y, G, row, col, rec, needg):
        x, y, vx, vy, dx, dy, dvx, dvy = X
        ux, uy = u
        Y[0] = dx - vx
        Y[1] = dy - vy
        Y[2] = ux + 2 * p[0] * (x - self.rx) - dvx
        Y[3] = uy + 2 * p[0] * (y - self.ry) - dvy
        if needg:
            G[:] = [1, -1, 1, -1, 1, 2*(x-self.rx), 2*p[0], -1, 1, 2*(y-self.ry), 2*p[0], -1]
            if rec:
                row[:] = [0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
                col[:] = [5, 3, 6, 4, 9, 11, 1, 7, 10, 11, 2, 8]


class CollocCircleConstr(nonLinearPointConstr):
    def __init__(self, with_vel=False, with_acce=False, order=2):
        nG = 2
        nf = 1
        self.order = order
        self.with_vel = with_vel
        self.with_acc = with_acce
        if with_vel:
            nG += 4
            nf += 1
            if with_acce:
                nG += 6
                nf += 1
        if order == 2:
            nx = 6
        else:
            nx = 8
        nu = 2
        nonLinearPointConstr.__init__(self, -1, nf, nx, nu, np=1, nG=nG)
        self.x = 0
        self.y = 0
        self.r = 1

    def fast_call(self, X):
        if self.order == 2:
            x, y, dx, dy, ddx, ddy = X
        else:
            x, y, vy, vy, dx, dy, ddx, ddy = X
        F = np.zeros(3)
        F[0] = (x-self.x)**2 + (y-self.y)**2 - self.r**2
        F[1] = (x-self.x)*dx + (y-self.y)*dy
        F[2] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy
        return F

    def __callg__(self, X, F, G, row, col, rec, needg):
        if self.order == 2:
            x, y, dx, dy, ddx, ddy = X[1:7]
        else:
            x, y, dx, dy, vx, vy, ddx, ddy = X[1:9]
        F[0] = (x-self.x)**2 + (y-self.y)**2 - self.r**2
        if self.with_vel:
            F[1] = (x-self.x)*dx + (y-self.y)*dy
            if self.with_acc:
                F[2] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy
        if needg:
            G[:2] = [2*(x-self.x), 2*(y-self.y)]
            if self.with_vel:
                if self.order == 2:
                    G[2:6] = [dx, dy, x - self.x, y - self.y]
                    if self.with_acc:
                        G[6:12] = [ddx, ddy, 2*dx, 2*dy, x - self.x, y - self.y]
                else:
                    G[2:6] = [dx, dy, x - self.x, y - self.y]
                    if self.with_acc:
                        G[6:12] = [ddx, ddy, 2*dx, 2*dy, x - self.x, y - self.y]
            if rec:
                row[:2] = [0, 0]
                col[:2] = [1, 2]
                if self.with_vel:
                    row[2:6] = [1, 1, 1, 1]
                    col[2:6] = [1, 2, 3, 4]
                    if self.with_acc:
                        row[6:12] = 2
                        if self.order == 2:
                            col[6:12] = np.arange(6) + 1
                        else:
                            col[6:12] = [1, 2, 3, 4, 7, 8]


class CircleConstr(manifoldConstr):
    def __init__(self):
        manifoldConstr.__init__(self, 8, 1, 2, constr_order=0, nG=12, nnzJx=2, nnzJg=2)
        self.x = 0
        self.y = 0
        self.r = 1

    def __callg__(self, X, F, G, row, col, rec, needg):
        """Calculate the constraints at knot points."""
        x, y, vx, vy, dx, dy, ddx, ddy = X
        F[0] = (x-self.x)**2 + (y-self.y)**2 - self.r**2
        F[1] = (x-self.x)*dx + (y-self.y)*dy
        F[2] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy
        if needg:
            G[:] = [2*(x-self.x), 2*(y-self.y), dx, dy, x - self.x, y - self.y, ddx, ddy, 2*dx, 2*dy, x - self.x, y-self.y]
            if rec:
                row[:] = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
                col[:] = [0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5]

    def __callf__(self, X, F, acce=False):
        x, y, vx, vy, dx, dy, ddx, ddy = X
        if not acce:
            F[0] = (x-self.x)**2 + (y-self.y)**2 - self.r**2
            F[1] = (x-self.x)*dx + (y-self.y)*dy
            F[2] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy
        else:
            F[0] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy

    def __calc_correction__(self, X, Gamma, Y, Gq, rowq, colq, Gg, rowg, colg, rec, needg):
        """Calculate the constraint J(q)^T \gamma and corresponding Jacobians"""
        x, y, vx, vy, dx, dy, ddx, ddy = X
        gamma = Gamma[0]
        Y[0] = 2 * (x - self.x) * gamma
        Y[1] = 2 * (y - self.x) * gamma
        if needg:
            Gg[:] = [2*(x-self.x), 2*(y-self.y)]
            Gq[:] = [2*gamma, 2*gamma]
            if rec:
                rowg[:] = [0, 1]
                colg[:] = [0, 0]
                rowq[:] = [0, 1]
                colq[:] = [0, 1]

    def __return_correction__(self, X, Gamma):
        x, y = X[:2]
        Y = np.array([2*(x - self.x) * Gamma[0],
                      2*(y - self.y) * Gamma[0]])
        return Y
