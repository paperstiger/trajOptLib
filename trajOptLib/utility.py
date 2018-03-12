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
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
import pyLib.plot as pld


def parseX(x, N, dimx, dimu, dimp, uset0, usetf, setfinal=True):
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
    assert len(x) == lenX
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


def showSol(solDict, xsplit=None, usplit=None, psplit=None):
    """Plot the parsed solution.

    :param solDict: dict, returned from solParse; or a list of such dict. Then they are compared
    :param xsplit: if state space is large, this splits it into several images it's a list of int like (3, 6). None means no split
    :param usplit: similar to xsplit, apply for control
    :param psplit: similar to xsplit, apply for parameter

    """
    assert isinstance(solDict, (dict, list))
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
    plt.show()


if __name__ == '__main__':
    main()
