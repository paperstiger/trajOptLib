#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
printio.py

Control print
"""
import sys


class Mute(list):
    """Capture std.output."""
    def __init__(self, mute=True, record=False):
        self.mute = mute
        self.record = record
        if self.record:
            self.mute = True

    def __enter__(self):
        if self.mute:
            self._stdout = sys.stdout
            if not self.record:
                sys.stdout = open('/dev/null')
            else:
                sys.stdout = self._stringio = StringIO()

    def __exit__(self, *args):
        if self.mute:
            if self.record:
                del self._stringio    # free up some memory
            sys.stdout = self._stdout

    def getOutput(self):
        if not self.record:
            return None
        else:
            return self._stringio.getvalue()
