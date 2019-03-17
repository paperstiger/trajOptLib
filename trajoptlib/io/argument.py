#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
argument.py

Easy arguments
"""
import argparse


def get_onoff_args(*args):
    """Get on-off arguments"""
    parser = argparse.ArgumentParser()
    for arg in args:
        assert isinstance(arg, str)
        parser.add_argument('-%s' % arg, action='store_true', default=False)
    return parser.parse_args()
