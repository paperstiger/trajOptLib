#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
compare.py

Extensively used when I want to compare several samples
"""
import numpy as np
import matplotlib.pyplot as plt
from .common import plotkwargs, getColorCycle, get3dAxis, getIndAlongAxis, scatterkwargs


def compare(arr, x=None, ax=None, transpose=False, show=False, **kwargs):
    """Given a matrix or 3d tensor, do some comparison

    arr is 3d tensor or list of 2d array
    x is the x-coordinate for each Cat, so it is a list of 1d nparray
    trnspose = True make the data be data by sample
    show = False means not plot immediately
    TODO:
    headtail = True means distinguish between 1st and last samples
    kwargs are allowed configurations, we have to pass in dict
    """
    colors = getColorCycle()
    if isinstance(arr, np.ndarray):
        nCat = arr.shape[0]
        if transpose:
            nFeature = arr.shape[1]
        else:
            nFeature = arr.shape[2]
    elif isinstance(arr, list):
        nCat = len(arr)
        for arr_ in arr:
            assert isinstance(arr_, np.ndarray)
        if transpose:
            nFeature = arr[0].shape[0]
        else:
            nFeature = arr[0].shape[1]
    # parse x for x-axis
    if x is not None:
        if isinstance(x, list):
            useList = True
        else:
            useList = False
    # get subplots
    nRow = int(np.floor(np.sqrt(nFeature)))
    if nFeature % nRow == 0:
        nCol = nFeature // nRow
    else:
        nCol = nFeature // nRow + 1
    # create figure
    if ax is None:
        fig, axes = plt.subplots(nRow, nCol)
        tight = True
    else:
        axes = ax  # we hope for the good
        tight = False
    for i in range(nFeature):
        row = i // nCol
        col = i % nCol
        try:
            ax = axes[row, col]
        except:
            try:
                ax = axes[col]
            except:
                ax = axes
        # plot for each one
        if len(kwargs) == 0 and x is None and isinstance(arr, np.ndarray):  # empty dict and no x information
            if transpose:
                ax.plot(arr[:, i, :].T)
            else:
                ax.plot(arr[:, :, i].T)
        else:
            for j in range(nCat):
                if transpose:
                    arr_ = arr[j][i, :]
                else:
                    arr_ = arr[j][:, i]
                dct = dict()
                try:
                    tmp = kwargs.iteritems()
                except:
                    tmp = kwargs.items()
                for key, item in tmp:
                        if isinstance(item, dict):
                            if j in item:
                                dct[key] = item[j]
                        elif key in plotkwargs:
                            dct[key] = item
                # dct = {key: item[j] for key, item in kwargs.iteritems() if j in item}
                if 'color' not in dct and 'c' not in dct:
                    dct['color'] = colors[j % len(colors)]  # avoid overflow
                if x is None:
                    ax.plot(arr_, **dct)
                else:
                    if useList:
                        ax.plot(x[j], arr_, **dct)
                    else:
                        ax.plot(x, arr_, **dct)
    if tight:
        fig.tight_layout()
    if show:
        plt.show()
    return axes


def compareXYZ(arr, ax=None, transpose=False, d3=False, scatter=False, show=False, **kwargs):
    """hybrid of compare, and plot. Assume we have a cat by N by dim dataset, we want to select a few col/row to plot in 2d/3d"""
    colors = getColorCycle()
    assert isinstance(arr, np.ndarray)
    assert arr.ndim == 3
    nCat = arr.shape[0]
    if transpose:
        alongaxis = 1
    else:
        alongaxis = 2
    # create figure
    if ax is None:
        if d3:
            fig, ax = get3dAxis()
        else:
            fig, ax = plt.subplots()
    # extract values
    xind = kwargs.get('xind', 0)
    yind = kwargs.get('xind', 1)
    # get values to plot
    tx = getIndAlongAxis(arr, alongaxis, xind)
    ty = getIndAlongAxis(arr, alongaxis, yind)
    if d3:
        zind = kwargs.get('zind', 2)  # get which column we should focus
        tz = getIndAlongAxis(arr, alongaxis, zind)
    # now we get a bunch of cat by N matrix, we plot cat by cat
    for j in range(nCat):
        # updated dct, properties can be set in bunch mode
        dct = dict()
        for key, item in kwargs.iteritems():
            if isinstance(item, dict):
                if j in item:
                    dct[key] = item[j]
            elif key in plotkwargs:
                dct[key] = item
        # dct = {key: item[j] for key, item in kwargs.iteritems() if j in item}
        if 'color' not in dct and 'c' not in dct:
            dct['color'] = colors[j % len(colors)]  # avoid overflow
        if d3:
            if scatter:
                ax.scatter(tx[j], ty[j], tz[j], **dct)
            else:
                ax.plot(tx[j], ty[j], tz[j], **dct)
        else:
            if scatter:
                ax.scatter(tx[j], ty[j], **dct)
            else:
                ax.plot(tx[j], ty[j], **dct)
    if show:
        plt.show()
    return ax
