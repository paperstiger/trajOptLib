#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
fileio.py

Useful functions for file io, especially dict of dict type.
This type is heavily used for clustered training.
"""
import numpy as np
import os
try:
    import cPickle as pkl
except:
    import pickle as pickle
import json


def ddctParse(fnm, thekey=None):
    """Parse a file that is dict of dict.
    I wish this is faster than pickle
    """
    if '.pkl' in fnm:
        with open(fnm, 'rb') as f:
            data = pkl.load(f)
        if thekey is not None:
            assert thekey in data.keys()
            return data[thekey]
        else:
            return data
    tmp = np.load(fnm)
    keys = tmp.keys()
    if thekey is None:  # load all data
        try:
            rst = dict()
            for key in keys:
                rst[key] = tmp[key].item()
        except:
            rst = dict()
            for key in keys:
                rst[key] = tmp[key]
        return rst
    else:
        assert thekey in tmp.keys()
        try:
            return tmp[thekey].item()
        except:
            return tmp[thekey]


def ddctSave(fnm, arr, pklmode=False):
    """Save a dict of dict"""
    assert isinstance(arr, dict)
    if not pklmode:
        np.savez(fnm, **arr)
    else:
        with open(fnm, 'wb') as f:
            pkl.dump(arr, f)


def getLogPath(fnm, debug=False):
    dirname, filename = os.path.split(fnm)
    if dirname == '':
        dirname = '.'
    if debug:
        return '%s/logs/%s.log.debug' % (dirname, filename)
    else:
        return '%s/logs/%s.log' % (dirname, filename)


def getJson(fnm):
    """Parse a json file as a dict"""
    with open(fnm) as f:
        return json.load(f)


def dumpJson(obj, fnm):
    with open(fnm, 'w') as f:
        json.dump(obj, f, indent=4)


def assignModule(module, dct):
    """Assign values in dict to the module"""
    try:
        [setattr(module, key, v) for k, v in dct.iteritems()]
    except:
        [setattr(module, key, v) for k, v in dct.items()]


def assignModuleByJson(module, fnm):
    """Change contents of a module by json file"""
    dct = getJsonConfig(fnm)
    assignModule(module, dct)
