#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
utility.py

Collection of utility functions such as solution parsing and showing.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline
from . import plot as pld


def parse_X(x, N, dimx, dimu, dimp, uset0, usetf, setfinal=True):
    """Parse the solution into a dict of time, state, control, parameter.

    :param x: ndarray, solution found by SNOPT
    :param N, dimx, dimu, dimp: discretization of the problem
    :param uset0, usetf: float/array like, specification of initial and final time
    :param setfinal: if we set the final element of u/p the same with second last
    :returns: a dictionary of keys 't', 'x', 'u', 'p'
    """
    if np.isscalar(uset0):
        fixt0 = True
    else:
        fixt0 = False
    if np.isscalar(usetf):
        fixtf = True
    else:
        fixtf = False
    lenX = N * dimx + N * dimu + N * dimp + 2 - (fixt0 + fixtf)
    assert len(x) >= lenX
    # parse state
    state = np.reshape(x[:N*dimx], (N, dimx))
    ctrl = np.reshape(x[N*dimx:N*(dimx+dimu)], (N, dimu))
    if setfinal:
        ctrl[-1] = ctrl[-2]
    if dimp > 0:
        param = np.reshape(x[N*(dimx+dimu):N*(dimx+dimp+dimu)], (N, dimp))
        if setfinal:
            param[-1] = param[-2]
    else:
        param = None
    if fixt0:
        t0 = uset0
        if fixtf:
            tf = usetf
        else:
            tf = x[-1]
    else:
        if fixtf:
            t0 = x[-1]
            tf = usetf
        else:
            t0 = x[-2]
            tf = x[-1]
    tgrid = np.linspace(t0, tf, N)
    return {'t': tgrid, 'x': state, 'u': ctrl, 'p': param}


def show_sol(solDict, xsplit=None, usplit=None, psplit=None, show=True):
    """Plot the parsed solution.

    :param solDict: dict, returned from solParse; or a list of such dict. Then they are compared
    :param xsplit: if state space is large, this splits it into several images it's a list of int like (3, 6). None means no split
    :param usplit: similar to xsplit, apply for control
    :param psplit: similar to xsplit, apply for parameter

    """
    assert isinstance(solDict, (dict, list))
    if 'phases' in solDict:
        mul_seg_mode = True
    else:
        mul_seg_mode = False
    if isinstance(solDict, dict):
        plott = [solDict['t']]
        plotx = solDict['x'][np.newaxis]
        plotu = solDict['u'][np.newaxis]
        if solDict['p'] is not None:
            plotp = solDict['p'][np.newaxis]
        else:
            plotp = None
    else:
        plott = [sol['t'] for sol in solDict]
        plotx = np.array([sol['x'] for sol in solDict])
        plotu = np.array([sol['u'] for sol in solDict])
        if solDict[0]['p'] is not None:
            plotp = np.array([sol['p'] for sol in solDict])
        else:
            plotp = None

    def cmpSplit(t, var, split):
        if split is None:
            axes = pld.compare(var, x=t)
        else:
            nsplit = len(split)
            ind0 = 0
            for i in range(nsplit):
                indf = ind0 + split[i]
                axes = pld.compare(var[:, :, ind0:indf], x=t)
                ind0 = indf

    # plot x
    cmpSplit(plott, plotx, xsplit)
    # plot u
    cmpSplit(plott, plotu, usplit)
    if plotp is not None:
        cmpSplit(plott, plotp, psplit)
    if show:
        plt.show()


def random_gen_in_bound(bds, n=None):
    """Randomly generate a vector within bounds / [-1, 1]

    :param bds: list/tuple of ndarray, the bounds of variable
    :param n: if bds are None, this is for determining size of vector
    :return x: ndarray, the random generated variable

    """
    if not isinstance(bds, list):
        return bds
    assert len(bds) == 2
    lb, ub = bds
    if lb is not None:
        lb = np.array(lb)
    if ub is not None:
        ub = np.array(ub)
    if (lb is not None) and (ub is not None):
        lb[lb < -1e19] = -1
        ub[ub > 1e19] = 1
        x = np.random.random(lb.shape) * (ub - lb) + lb
    elif lb is None:
        if ub is None:
            assert n is not None
            x = np.random.random(n) * 2 - 1
        else:
            x = ub - np.random.random(ub.shape)
    else:
        assert ub is None
        x = lb + np.random.random(lb.shape)
    return x


def check_in_bounds(x, bds):
    """Check the position of variables in a bound.

    Basically, it returns the position of actual value in the bounds [-1, 1].
    If both sides are bounded, it simply calculate this value.
    If one side is bounded, it check if value if within certain threshold with the bound
    Unbounded variables shall always return 0, the same applies to equality-bounded variables.

    """
    threshold = 1e-6
    lb, ub = bds
    if np.isscalar(x) and np.isscalar(lb) and np.isscalar(ub):
        x = np.atleast_1d(x)
        lb = np.atleast_1d(lb)
        ub = np.atleast_1d(ub)
    assert len(lb) == len(ub)
    lenbd = len(lb)
    if x.ndim == 1:
        assert len(x) % lenbd == 0
        X = np.reshape(x, (-1, lenbd))
    elif x.ndim == 2:
        assert x.shape[1] == lenbd
        X = x
    position = np.zeros_like(X)
    for i in range(lenbd):
        lb_ = np.atleast_1d(lb[i])
        ub_ = np.atleast_1d(ub[i])
        assert lb_ <= ub_
        if lb_ <= -1e19:
            if ub_ > 1e19:
                position[:, i] = 0
            else:
                position[:, i] = np.where(X[:, i] < ub_ - threshold, 0, 1)
        else:
            if ub_ > 1e19:
                position[:, i] = np.where(X[:, i] > lb_ + threshold, 0, -1)
            elif ub_ > lb_:
                position[:, i] = -1.0 + 2.0 * (X[:, i] - lb_) / (ub_ - lb_)
            else:
                position[:, i] = np.where(np.abs(X[:, i] - lb_) < threshold, 0, np.where(X[:, i] > lb_, 1, -1))
    if x.ndim == 1:
        return position.flatten()
    else:
        return position


def interp(t, X, teval, Xeval, kind):
    """Do interpolation on previous calculated solution.

    It handles case when t is None, in which case, teval is not used and we use a uniform grid in [0, 1]

    :param t: array-like, user-specified time stamps
    :param X: ndarray, (x, x), variable to be interpolated
    :param teval: array-like, where to evaluate for the interpolation
    :param Xeval: ndarray, (x, x), where to store the evaluation. It might has more columns than X, in which case we fill higher-order derivative
    :param kind: str, interpolation type for scipy.interpolate.interp1d, can be (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’)

    """
    targetN = Xeval.shape[0]
    giveN = X.shape[0]
    Xcol = X.shape[1]
    XevalCol = Xeval.shape[1]
    order = XevalCol // Xcol
    if t is None:
        teval = np.linspace(0, 1, targetN)
        t = np.linspace(0, 1, giveN)
    # construct object
    if order == 1:
        interp_ = interp1d(t, X, kind=kind, axis=0)
        Xeval[:] = interp_(teval)
    else:
        interp_ = CubicSpline(t, X)
        Xeval[:, :X.shape[1]] = interp_(teval)
        curind = 1
        while order > 1:
            interp_ = interp_.derivative()
            Xeval[:, curind*Xcol: (curind+1)*Xcol] = interp_(teval)
            order -= 1
            curind += 1
    return


def get_inf(n=None):
    """Return an inf array.

    :param n: int, size of the array to return, None means scalar
    :returns y: the scalar inf or array of inf

    """
    if n is None:
        return 1e20
    else:
        return 1e20 * np.ones(n)


class InfBuilder(object):
    """A class to help us return infinity"""
    def __init__(self, n=20):
        self.n = n
        self.inf = np.ones(n) * 1e20
        self.ninf = -self.inf
        self.inf.flags.writeable = False
        self.ninf.flags.writeable = False

    def __getitem__(self, n):
        if n > 0 and n < self.n:
            return self.inf[:n]
        if n < 0 and n > -self.n:
            return self.ninf[:-n]
        return np.ones(abs(n)) * 1e20 * np.sign(n)


class OneBuilder(object):
    """A class to return array of 1"""
    def __init__(self, n=20):
        self.n = n
        self.ones = np.ones(n)
        self.ones.flags.writeable = False

    def __getitem__(self, n):
        if n > 0 and n < self.n:
            return self.ones[:n]
        return np.ones(n)


class ZeroBuilder(object):
    """A class to return array of 1"""
    def __init__(self, n=20):
        self.n = n
        self.zeros = np.zeros(n)
        self.zeros.flags.writeable = False

    def __getitem__(self, n):
        if n > 0 and n < self.n:
            return self.zeros[:n]
        return np.zeros(n)
