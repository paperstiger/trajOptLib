#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""

"""
from .common import get3dAxis, getColorCycle, getIndAlongAxis, setFontSize
from .plot3d import  plot, scatter, set_axes_equal, addSphere
from .compare import compare, compareXYZ
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def show(*args, **kwargs):
    plt.show(*args, **kwargs)


def getSizeFont(size=16):
    font = mpl.font_manager.FontProperties(size=size)
    return font


def setTickSize(ax, size=16):
    ax.tick_params(axis='x', labelsize=size)
    ax.tick_params(axis='y', labelsize=size)
    if ax.name == '3d':
        ax.tick_params(axis='z', labelsize=size)


def alignRange(*args):
    """Given many axes, align the range"""
    xlims = np.array([arg.get_xlim() for arg in args])
    ylims = np.array([arg.get_ylim() for arg in args])
    left = np.amin(xlims[:, 0])
    right = np.amax(xlims[:, 1])
    up = np.amax(ylims[:, 1])
    down = np.amin(ylims[:, 0])
    for arg in args:
        arg.set_xlim((left, right))
        arg.set_ylim((down, up))
