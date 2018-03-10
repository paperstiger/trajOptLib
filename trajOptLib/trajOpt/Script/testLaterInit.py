#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testLaterInit.py

Test if we can initialize a class later.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class A(object):
    def __init__(self, a):
        self.A = a


class B(A):
    def __init__(self, b):
        self.b = b

    def seta(self, a):
        self.a = a


def main():
    b = B(5)
    b.seta(10)
    print(b.a, b.b)


if __name__ == '__main__':
    main()
