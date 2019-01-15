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
import pyoptsolver as libsnopt


def parseResult(rst):
    """Parse the results returned by snopt and convert to a dict."""
    return {'flag': rst.flag, 'obj': rst.obj, 'x': rst.sol, 'f': rst.fval}


def directSolve(fun, x0, nf=None, xlb=None, xub=None, clb=None, cub=None, cfg=None):
    """Directly solve the optimization problem described using fun with guess x0

    :param fun: A function like y = f(x) where x, y are np.ndarray
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param nf: int, length of y
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :param cfg: libsnopt.SnoptConfig, configuration of snopt solver
    :returns: a dictionary containing the solution

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
    if cfg is None:
        cfg = libsnopt.SnoptConfig()
    rst = libsnopt.directSolve(fun, x0, nx, nf, xlb, xub, clb, cub, cfg)
    return parseResult(rst)


def inDirectSolve(fun, x0, nf=None, xlb=None, xub=None, clb=None, cub=None, cfg=None):
    """Directly solve the optimization problem described using fun with guess x0

    :param fun: A function like f(x, y) where x, y are np.ndarray
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param nf: int, length of y
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :param cfg: libsnopt.SnoptConfig, configuration of snopt solver
    :returns: a dictionary containing the solution

    """
    nx = len(x0)
    if nf is None:
        if clb is not None and cub is not None:
            assert len(clb) == len(cub)
            nf = len(clb)
    assert nf is not None
    if xlb is None or xub is None:
        xlb = np.empty(0)
        xub = np.empty(0)
    if clb is None or cub is None:
        clb = np.empty(0)
        cub = np.empty(0)
    if cfg is None:
        cfg = libsnopt.SnoptConfig()
    rst = libsnopt.inDirectSolve(fun, x0, nx, nf, xlb, xub, clb, cub, cfg)
    return parseResult(rst)


def gradSolve(fun, x0, nf=None, xlb=None, xub=None, clb=None, cub=None, cfg=None):
    """Directly solve the optimization problem described using fun with guess x0

    :param fun: A function like y, J = f(x) where x, y, J are np.ndarray
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param nf: int, length of y
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :param cfg: libsnopt.SnoptConfig, configuration of snopt solver
    :returns: a dictionary containing the solution

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
    if cfg is None:
        cfg = libsnopt.SnoptConfig()
    rst = libsnopt.gradSolve(fun, x0, nx, nf, xlb, xub, clb, cub, cfg)
    return parseResult(rst)


def inGradSolve(fun, x0, nf=None, xlb=None, xub=None, clb=None, cub=None, cfg=None):
    """Directly solve the optimization problem described using fun with guess x0

    :param fun: A function like f(x, y, J) where x, y, J are np.ndarray
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param nf: int, length of y
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :param cfg: libsnopt.SnoptConfig, configuration of snopt solver
    :returns: a dictionary containing the solution

    """
    nx = len(x0)
    if nf is None:
        if clb is not None and cub is not None:
            assert len(clb) == len(cub)
            nf = len(clb)
    assert nf is not None
    if xlb is None or xub is None:
        xlb = np.empty(0)
        xub = np.empty(0)
    if clb is None or cub is None:
        clb = np.empty(0)
        cub = np.empty(0)
    if cfg is None:
        cfg = libsnopt.SnoptConfig()
    rst = libsnopt.inGradSolve(fun, x0, nx, nf, xlb, xub, clb, cub, cfg)
    return parseResult(rst)


def spGradSolve(fun, x0, nf=None, nG=None, xlb=None, xub=None, clb=None, cub=None, cfg=None):
    """Directly solve the optimization problem described using fun with guess x0

    :param fun: A function like y, spJ = f(x) where x, y are np.ndarray, J is scipy.sparse.csc_matrix
    :param nf: int, length of y
    :param nG: int, nnz of spJ
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :param cfg: libsnopt.SnoptConfig, configuration of snopt solver
    :returns: a dictionary containing the solution

    """
    nx = len(x0)
    if nf is None:
        if clb is not None and cub is not None:
            assert len(clb) == len(cub)
            nf = len(clb)
        else:
            y, spJ = fun(x0)
            nf = len(y)
            nG = spJ.nnz
    if nG is None:
        y, spJ = fun(x0)
        nG = spJ.nnz
    assert nf is not None
    assert nG is not None
    if xlb is None or xub is None:
        xlb = np.empty(0)
        xub = np.empty(0)
    if clb is None or cub is None:
        clb = np.empty(0)
        cub = np.empty(0)
    if cfg is None:
        cfg = libsnopt.SnoptConfig()
    rst = libsnopt.spGradSolve(fun, x0, nx, nf, nG, xlb, xub, clb, cub, cfg)
    return parseResult(rst)


def inSpGradSolve(fun, x0, nf=None, nG=None, xlb=None, xub=None, clb=None, cub=None, cfg=None):
    """Directly solve the optimization problem described using fun with guess x0

    :param fun: A function like f(x, y, G, row, col, rec) where x, y are np.ndarray, J is scipy.sparse.csc_matrix
    :param x0: np.ndarray (nx,) the initial guess to the solver
    :param nf: int, number of f
    :param nG: int number nonzero in Jacobian
    :param xlb: np.ndarray (nx,) lower bound on decision variable x
    :param xub: np.ndarray (nx,) upper bound on decision variable x
    :param clb: np.ndarray (nc,) lower bound on return function c
    :param cub: np.ndarray (nc,) upper bound on return function c
    :param cfg: libsnopt.SnoptConfig, configuration of snopt solver
    :returns: a dictionary containing the solution

    """
    nx = len(x0)
    if nf is None:
        if clb is not None and cub is not None:
            assert len(clb) == len(cub)
            nf = len(clb)
    assert nf is not None
    assert nG is not None
    if xlb is None or xub is None:
        xlb = np.empty(0)
        xub = np.empty(0)
    if clb is None or cub is None:
        clb = np.empty(0)
        cub = np.empty(0)
    if cfg is None:
        cfg = libsnopt.SnoptConfig()
    rst = libsnopt.inSpGradSolve(fun, x0, nx, nf, nG, xlb, xub, clb, cub, cfg)
    return parseResult(rst)
