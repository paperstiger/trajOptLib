#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testBind.py

I want to use autograd to automatically get jacobians.
I test if I can write a decorator to return a wrapper.
"""
import sys, os, time
import matplotlib.pyplot as plt
from libsnopt import snoptConfig, probFun, solver
import logging
from autograd import jacobian, grad
import autograd.numpy as np
from pyLib.math import blockIndex


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class wrapper(probFun):
    """This class takes a function as input and returns a suitable probFun."""
    def __init__(self, fun, nx, nf, *args):
        probFun.__init__(self, nx, nf)
        self.grad = True
        self.nG = nx * nf
        self.fun = fun
        self.gfun = jacobian(fun)
        self.funargs = args

    def __callf__(self, x, y):
        y0 = self.fun(x, *self.funargs)
        y[:] = y0

    def __callg__(self, x, y, G, row, col, rec, needg):
        y0 = self.fun(x, *self.funargs)
        y[:] = y0
        if needg:
            g0 = self.gfun(x, *self.funargs)
            G[:] = g0.flat
            if rec:
                row[:], col[:] = blockIndex(0, 0, self.nf, self.nx)


def testFun(x, p):
    y0 = x[0] ** 2 + x[1] ** 2 + p
    y1 = x[0] + x[1]
    return np.array([y0, y1])


def main():
    prob = wrapper(testFun, 2, 2, 1)
    prob.lb = [-1e10, 1.0]
    prob.ub = [1e10, 1.0]
    prob.xlb = [0.6, -1]
    prob.xub = [1, 1]
    cfg = snoptConfig()
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag, rst.sol)


if __name__ == '__main__':
    main()
