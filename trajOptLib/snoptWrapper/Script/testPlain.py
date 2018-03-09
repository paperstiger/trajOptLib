#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testPlain.py

Test the wrapper for plain function
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
import libsnopt
from snoptWrapper import directSolve, inDirectSolve, gradSolve, inGradSolve, spGradSolve, inSpGradSolve
from scipy.sparse import csc_matrix
from pyLib.io import getOnOffArgs


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    args = getOnOffArgs('plain', 'grad', 'sparse', 'all')
    if args.plain or args.all:
        testPlain()
    if args.grad or args.all:
        testGrad()
    if args.sparse or args.all:
        testSparse()


def testPlain():
    x0 = np.random.random(2)
    xlb = np.array([0.6, 0.2])
    xub = np.array([1.6, 1.2])
    clb = np.array([0, 1.0])
    cub = np.array([0, 1.0])
    # test plain fun
    rst = directSolve(plainFun, x0, nf=None, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = directSolve(plainFun, x0, nf=None, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)
    # test inplain fun
    rst = inDirectSolve(inPlainFun, x0, nf=2, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = inDirectSolve(inPlainFun, x0, nf=2, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)


def testGrad():
    x0 = np.random.random(2)
    xlb = np.array([0.6, 0.2])
    xub = np.array([1.6, 1.2])
    clb = np.array([0, 1.0])
    cub = np.array([0, 1.0])
    # test grad fun
    rst = gradSolve(gradFun, x0, nf=None, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = gradSolve(gradFun, x0, nf=None, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)
    # test ingrad fun
    rst = inGradSolve(inGradFun, x0, nf=2, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = inGradSolve(inGradFun, x0, nf=2, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)


def testSparse():
    x0 = np.random.random(2)
    xlb = np.array([0.6, 0.2])
    xub = np.array([1.6, 1.2])
    clb = np.array([0, 1.0])
    cub = np.array([0, 1.0])
    # test grad fun
    rst = spGradSolve(spGradFun, x0, nf=None, nG=None, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = spGradSolve(spGradFun, x0, nf=None, nG=None, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)
    rst = inSpGradSolve(inSpGradFun, x0, nf=2, nG=4, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = inSpGradSolve(inSpGradFun, x0, nf=2, nG=4, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)


def plainFun(x):
    """The naive example"""
    y = x[0] **2 + x[1] ** 2
    c = x[0] + x[1]
    return np.array([y, c])


def inPlainFun(x, f):
    """A naive example with inplace fun"""
    f[0] = x[0] **2 + x[1] ** 2
    f[1] = x[0] + x[1]


def gradFun(x):
    """The naive example with grad"""
    y = x[0] **2 + x[1] ** 2
    c = x[0] + x[1]
    y = np.array([y, c])
    J = np.zeros((2, 2))
    J[0] = 2 * x
    J[1] = 1.0
    return y, J


def inGradFun(x, f, J):
    """A naive example with inplace fun"""
    f[0] = x[0] **2 + x[1] ** 2
    f[1] = x[0] + x[1]
    J[0] = 2 * x
    J[1] = 1.0


def spGradFun(x):
    """The naive example with sparse grad"""
    y = x[0] **2 + x[1] ** 2
    c = x[0] + x[1]
    y = np.array([y, c])
    J = np.zeros((2, 2))
    J[0] = 2 * x
    J[1] = 1.0
    spJ = csc_matrix(J)
    return y, spJ


def inSpGradFun(x, y, G, row, col, rec):
    """The naive example with sparse grad"""
    y[0] = x[0] **2 + x[1] ** 2
    y[1] = x[0] + x[1]
    G[:2] = 2 * x
    G[2:] = 1.0
    if rec:
        row[:] = [0, 0, 1, 1]
        col[:] = [0, 1, 0, 1]


if __name__ == '__main__':
    main()
