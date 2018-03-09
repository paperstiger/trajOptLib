#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
snoptWrapper.py


Wrapper functions for calling snopt
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
import libsnopt


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parseResult(rst):
    return {'flag': rst.flag, 'obj': rst.obj, 'x': rst.sol, 'f': rst.fval}


def directSolve(fun, x0, nf=None, xlb=None, xub=None, clb=None, cub=None):
    """Directly solve the optimization problem described using fun with guess x0

    Usage::
    :param fun: A function like y = f(x) where x, y are np.ndarray
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :rtype rst: a dictionary containing the solution
    """
    nx = len(x0)
    if nf is None:
        if clb is not None and cub is not None:
            assert len(clb) == len(cub)
            nf = len(clb)
        else:
            y = fun(x0)
            nf = len(y)
    if xlb is None or xub is None:
        xlb = np.empty(0)
        xub = np.empty(0)
    if clb is None or cub is None:
        clb = np.empty(0)
        cub = np.empty(0)
    rst = libsnopt.directSolve(fun, x0, nx, nf, xlb, xub, clb, cub)
    return parseResult(rst)


if __name__ == '__main__':
    main()
