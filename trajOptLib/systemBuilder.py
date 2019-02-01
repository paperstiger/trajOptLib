#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 motion <motion@motion-MS-7971>
#
# Distributed under terms of the MIT license.

"""
Code to build dynamical systems
"""
from __future__ import print_function
import numpy as np
from scipy.sparse import coo_matrix

from .trajOptBase import daeSystem, system


class GeometricRobot(daeSystem):
    """A class that defines a geometric robot, no dynamics involved"""
    def __init__(self, dim, dimp=0):
        self.dim = dim
        self.dimp = dimp
        daeSystem.__init__(self, 3 * dim, dim, dimp, dim, 2 * dim)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        """Write dynamics function"""
        y[:] = x[2 * self.dim: 3 * self.dim] - u
        if needg:
            G[:self.dim] = 1
            G[self.dim:] = -1
            if rec:
                row[:self.dim] = np.arange(self.dim)
                col[:self.dim] = 1 + 2 * self.dim + np.arange(self.dim)
                row[self.dim:] = np.arange(self.dim)
                col[self.dim:] = 1 + 3 * self.dim + np.arange(self.dim)


class GeometricRobotSystem(system):
    """A class that implements a geometric robot, second order."""
    def __init__(self, dim, dimp=0, ode='Euler'):
        self.dim = dim
        self.dimp = dimp
        system.__init__(self, 2 * dim, dim, dimp, ode)

    def Jdyn(self, t, x, u, p, h=None):
        dx = np.concatenate((x[self.dim:], u))
        jac = coo_matrix((np.ones(2 * self.dim),
                          (np.arange(2 * self.dim), np.arange(2 * self.dim) + 1 + self.dim)))
        return dx, jac
