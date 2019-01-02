#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
plot3d.py
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from .common import getColorCycle, getIndAlongAxis, get3dAxis


def plot(x, y=None, z=None, ax=None, axis=1, noz=False, show=False, scatter=False, **kwargs):
    """
    Easy-2-use function for plotting, either line or scatter is supported, 2d and 3d are supported
    Refer scatter for document
    """
    if y is not None:
        if y.ndim == 1:
            ty = y
        else:
            yaxis = kwargs.get('yaxis', 0)
            yind = kwargs.get('yind', 0)
            ty = getIndAlongAxis(y, yaxis, yind)
        # do the same for x
        xaxis = kwargs.get('xaxis', 0)
        xind = kwargs.get('xind', 0)
        tx = getIndAlongAxis(x, xaxis, xind)

    if z is not None:
        assert y is not None
        if z.ndim == 1:
            tz = z
        else:
            zaxis = kwargs.get('zaxis', 0)
            zind = kwargs.get('zind', 0)
            tz = getIndAlongAxis(z, zaxis, zind)
    else:
        if y is None:
            assert x.ndim > axis
            xind = kwargs.get('xind', 0)
            yind = kwargs.get('yind', 1)
            tx = getIndAlongAxis(x, axis, xind)
            ty = getIndAlongAxis(x, axis, yind)
            if x.shape[axis] >= 3 and not noz:
                zind = kwargs.get('zind', 2)
                tz = getIndAlongAxis(x, axis, zind)
            else:
                tz = None
        else:
            tz = None  # y is not None, z is None, then only 2d
    # construct a dict
    cfgDct = dict()
    allowKeys = ['color', 'c', 'label', 'linestyle', 'ls', 'linewidth', 'lw', 'marker']
    for key in allowKeys:
        if key in kwargs:
            cfgDct[key] = kwargs[key]
    if tz is not None:
        if ax is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
        if scatter:
            hdl = ax.scatter(tx, ty, tz, **cfgDct)
        else:
            hdl = ax.plot(tx, ty, tz, **cfgDct)
    else:
        if ax is None:
            fig, ax = plt.subplots()
        if scatter:
            hdl = ax.scatter(tx, ty, **cfgDct)
        else:
            hdl = ax.plot(tx, ty, **cfgDct)
    if show:
        plt.show()
    if scatter:
        return hdl
    else:
        return hdl[-1]  # return handle


def scatter(x, y=None, z=None, ax=None, axis=1, noz=False, show=False, **kwargs):
    """
    plot scatter figure. 
    if y and z is None, simply use M, col determines 2d (default) or 3d
    axis determines column-wise (default) or row-wise plot
    if y is not None, and M and y are both 1d, plot 2d, otherwise do it for 1st axis
    if z is not None, do the same for 3d case
    kwargs might be xaxis, xind, yaxis, yind, zaxis, zind
    """
    if y is not None:
        if y.ndim == 1:
            ty = y
        else:
            yaxis = kwargs.get('yaxis', 0)
            yind = kwargs.get('yind', 0)
            ty = getIndAlongAxis(y, yaxis, yind)
        # do the same for x
        xaxis = kwargs.get('xaxis', 0)
        xind = kwargs.get('xind', 0)
        tx = getIndAlongAxis(x, xaxis, xind)

    if z is not None:
        assert y is not None
        if z.ndim == 1:
            tz = z
        else:
            zaxis = kwargs.get('zaxis', 0)
            zind = kwargs.get('zind', 0)
            tz = getIndAlongAxis(z, zaxis, zind)
    else:
        if y is None:
            assert x.ndim > axis
            xind = kwargs.get('xind', 0)
            yind = kwargs.get('yind', 1)
            tx = getIndAlongAxis(x, axis, xind)
            ty = getIndAlongAxis(x, axis, yind)
            if x.shape[axis] >= 3 and not noz:
                zind = kwargs.get('zind', 2)
                tz = getIndAlongAxis(x, axis, zind)
            else:
                tz = None
        else:
            tz = None  # y is not None, z is None, then only 2d
    # construct a dict
    cfgDct = dict()
    allowKeys = ['color', 's', 'marker', 'c', 'cmap', 'norm', 'vmin', 'vmax', 'alpha']
    for key in allowKeys:
        if key in kwargs:
            cfgDct[key] = kwargs[key]
    if tz is not None:
        if ax is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
        hdl = ax.scatter(tx, ty, tz, **cfgDct)
    else:
        if ax is None:
            fig, ax = plt.subplots()
        hdl = ax.scatter(tx, ty, **cfgDct)
    if show:
        plt.show()
    return hdl


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
    if ax.name != "3d":
        print('Warning, axis is not 3d, exit')
        return

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def addSphere(ax, xc, yc, zc, r, elev=10, color='g', alpha=0.2, **kwargs):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = r * np.outer(np.cos(u), np.sin(v)) + xc
    y = r * np.outer(np.sin(u), np.sin(v)) + yc
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + zc
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color=color, alpha=alpha, linewidth=0, **kwargs)
    rot = 80.0 / 180.0 * np.pi
    a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
    b = np.array([0, 1, 0])
    b = b * np.cos(rot) + np.cross(a, b) * np.sin(rot) + a * np.dot(a, b) * (1 - np.cos(rot))
    ax.plot(xc + r * np.sin(u), yc + r*np.cos(u), zc,color='k', linestyle = 'dashed')
    horiz_front = np.linspace(0, np.pi, 100)
    ax.plot(xc + r * np.sin(horiz_front), yc + r * np.cos(horiz_front), zc, color='k')
    vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
    ax.plot(xc + r*(a[0] * np.sin(u) + b[0] * np.cos(u)), yc + r*( b[1] * np.cos(u)), zc + r*(a[2] * np.sin(u) + b[2] * np.cos(u)), color='k', linestyle='dashed')
    ax.plot(xc + r*(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front)), yc + r*(b[1] * np.cos(vert_front)), zc + r*(a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front)), color='k')
