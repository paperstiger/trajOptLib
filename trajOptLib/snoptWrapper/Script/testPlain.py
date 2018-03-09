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
from snoptWrapper import directSolve


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def plainFun(x):
    """The naive example"""
    y = x[0] **2 + x[1] ** 2
    c = x[0] + x[1]
    return np.array([y, c])


def main():
    x0 = np.random.random(2)
    xlb = np.array([0.6, 0.2])
    xub = np.array([1.6, 1.2])
    clb = np.array([0, 1.0])
    cub = np.array([0, 1.0])
    rst = directSolve(plainFun, x0, nf=None, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)
    rst = directSolve(plainFun, x0, nf=None, xlb=None, xub=None, clb=clb, cub=cub)
    print(rst)


if __name__ == '__main__':
    main()
