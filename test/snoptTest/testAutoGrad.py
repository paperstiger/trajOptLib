#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testAutoGrad.py

Wrap a function and solve it using snopt.
Test if autograd can be used to provide gradient.
"""
import sys, os, time
import numpy
import matplotlib.pyplot as plt
import logging
import autograd.numpy as np
from autograd import jacobian, elementwise_grad, grad
from functools import partial
from trajOptLib import libsnopt
from trajOptLib.snoptWrapper import gradSolve


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def fun(x):
    d = np.sum(x)
    c = np.sum(x ** 2)
    return np.array([c, d])


def gradFun(x, f, J):
    """The naive example with grad"""
    y = f(x)
    g = J(x)
    return y, g


def main():
    J = jacobian(fun)

    def wrapper(x):
        return fun(x), J(x)

    xlb = np.array([0.6, 0.2])
    xub = np.array([1.6, 1.2])
    clb = np.array([0, 1.0])
    cub = np.array([0, 1.0])
    x0 = numpy.random.random(2)
    # test grad fun
    rst = gradSolve(wrapper, x0, nf=None, xlb=xlb, xub=xub, clb=clb, cub=cub)
    print(rst)


if __name__ == '__main__':
    main()
