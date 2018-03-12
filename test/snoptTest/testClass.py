#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testClass.py

Test if we can use the class provided by pybind11 to have more control of how to solve the problem.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
from libsnopt import snoptConfig, probFun, solver


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class problem(probFun):
    def __init__(self):
        probFun.__init__(self, 2, 2)
        self.grad = True
        self.nG = 4

    def __callf__(self, x, y):
        y[0] = x[0] ** 2 + x[1] ** 2
        y[1] = x[0] + x[1]

    def __callg__(self, x, y, G, row, col, rec, needg):
        y[0] = x[0] ** 2 + x[1] ** 2
        y[1] = x[0] + x[1]
        G[:2] = 2 * x
        G[2:] = 1.0
        if rec:
            row[:] = [0, 0, 1, 1]
            col[:] = [0, 1, 0, 1]


def main():
    prob = problem()
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
