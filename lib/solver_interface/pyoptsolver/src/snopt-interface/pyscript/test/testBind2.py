#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testBind2.py

This basically test if I do not use autograd.numpy module when defining function, what can I get.
"""
import sys, os, time
import matplotlib.pyplot as plt
from libsnopt import snoptConfig, probFun, solver
import logging
import numpy as np
import autograd.numpy as numpy
from pyLib.math import blockIndex
from testBind import wrapper


def testFun1(x, p):
    """This function actually uses no numpy function."""
    y = np.zeros(2)
    y0 = x[0] ** 2 + x[1] ** 2 + p
    y1 = x[0] + x[1]
    return np.array([y0, y1])


def testFun2(x, p):
    y = numpy.zeros(2)
    y0 = numpy.dot(x, x) + p
    y1 = x[0] + x[1]
    return numpy.array([y0, y1])


def main():
    prob = wrapper(testFun2, 2, 2, 1)
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
