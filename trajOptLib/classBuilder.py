#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
classBuilder.py

This module provides versatile functions that helps user build up classes quickly.
Specifically, it allows fast prototype of problems.
However, the users have to write functions that are autograd compatible.
Most basically, import autograd.numpy instead of numpy
"""
import autograd.numpy as np
from autograd import jacobian, grad
from .trajOptBase import system, daeSystem, baseFun, addX
from .trajOptBase import nonLinearPointConstr
from .libsnopt import snoptConfig, solver, probFun, result


def blockIndex(i, j, rows, cols, order='C'):
    """For a matrix block, we return the index of row and columns.

    For a matrix we choose a block using the upper left corner positioned 
    at (i, j) and size (row, col). Each element of the block has row and 
    col index, they are returned in two arrays. The order controls we use
    row or column major order

    For example, blockIndex(1, 3, 2, 3, 'C') returns 
    (array([1, 1, 1, 2, 2, 2]), array([3, 4, 5, 3, 4, 5]))

    Parameters
    ----------
    i: int
        the row of the upper left corner
    j: int
        the column of the upper left corner
    rows: int
        number of rows of the block
    cols: int
        number of columns of the block
    order: char
        ('C'/'F') if we return row or column major
    """
    if order == 'C':
        row = i + (np.arange(rows)[:, np.newaxis] + np.zeros(cols)).flatten()
        col = j + (np.zeros(rows)[:, np.newaxis] + np.arange(cols)).flatten()
    elif order == 'F':
        row = i + (np.zeros(cols)[:, np.newaxis] + np.arange(rows)).flatten()
        col = j + (np.arange(cols)[:, np.newaxis] + np.zeros(rows)).flatten()
    else:
        raise Exception("Unsupported order")
    return row, col


class systemWrapper(system):
    """This class takes a function and returns a system."""
    def __init__(self, fun, nx, nu, np, *args):
        """Constructor.

        Parameters
        ----------
        fun : callable
            a function that implements \dot{x}=f(t, x, u, p, \*args) but does not depend on t
        nx : int
            length of x
        nu : int
            length of u
        np : int
            length of p
        args : kwargs
            additional parameters to function
        """
        system.__init__(self, nx, nu, np)

        def wrapfun(X, *args):
            x = X[:nx]
            u = X[nx:nx+nu]
            p = X[nx+nu:nx+nu+np]
            return fun(0, x, u, p, *args)

        self.grad = True
        self.nG = nx * (nx + nu + np)
        self.fun = wrapfun
        self.gfun = jacobian(wrapfun)
        self.funargs = args

    def dyn(self, t, x, u, p):
        xin = np.concatenate((x, u, p))
        y = self.fun(xin, *self.funargs)
        return y

    def Jdyn(self, t, x, u, p):
        xin = np.concatenate((x, u, p))
        y = self.fun(xin, *self.funargs)
        G = self.gfun(xin, *self.funargs)
        return y, G


class daeSystemWrapper(daeSystem):
    """This class takes a function and returns a system."""
    def __init__(self, fun, nx, nu, np, nf, *args):
        """Constructor for the problem.

        Parameters
        ----------
        fun : callable
            a function that implements f(t, x, u, p, \*args) = 0 but does not depend on t
        nx : int
            length of x
        nu : int
            length of u
        np : int
            length of p
        nf : int
            length of output of f
        args : kwargs
            additional parameters
        """
        nG = nf * (nx + nu + np)
        daeSystem.__init__(self, nx, nu, np, nf, nG)

        def wrapfun(X, *args):
            x = X[:nx]
            u = X[nx: nx+nu]
            p = X[nx+nu:nx+nu+np]
            return fun(0, x, u, p, *args)

        self.fun = wrapfun
        self.gfun = jacobian(wrapfun)
        self.funargs = args

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        """Override class method"""
        X = np.concatenate((x, u, p))
        out = self.fun(X, *self.funargs)
        y[:] = out
        if needg:
            g = self.gfun(X, *self.funargs)
            G[:] = g.flat
            if rec:
                row[:], col[:] = blockIndex(0, 1, self.nf, self.nx + self.nu + self.np)


class nonLinearPointConstrWrapper(nonLinearPointConstr):
    """This class takes a function and wrap it as a point constraint."""
    def __init__(self, fun, nx, nu, np, nc, index, lb=None, ub=None, args=None):
        nG = nc * (nx + nu + np)
        nonLinearPointConstr.__init__(self, index, nc, nx, nu, np, lb, ub, nG=nG)

        def wrapfun(X, *args):
            x = X[:nx]
            u = X[nx:nx+nu]
            p = X[nx+nu:nx+nu+np]
            return fun(0, x, u, p, *args)

        self.fun = wrapfun
        self.gfun = jacobian(wrapfun)
        if args is not None:
            self.args = args
        else:
            self.args = ()

    def __callg__(self, X, F, G, row, col, rec, needg):
        """override a function"""
        # we will assume first value is useless
        y = self.fun(X[1:], *self.args)
        F[:] = y
        if needg:
            g = self.gfun(X[1:], *self.args)
            G[:] = g
            if rec:
                row[:], col[:] = blockIndex(0, 1, self.nf, self.nx - 1)
