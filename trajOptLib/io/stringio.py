#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
stringio.py

Provides subroutines to extract numbers from string
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
import re


#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


numRe = re.compile(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?")


def getNumber(string, mapfun=None):
    """Parse all numbers from a string"""
    if mapfun is None:
        return numRe.findall(string)
    else:
        return map(mapfun, numRe.findall(string))
