#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 motion <motion@motion-MS-7971>
#
# Distributed under terms of the MIT license.

"""
On a naive QP problem, test all solvers
"""
from __future__ import print_function
import numpy as np
from pyoptsolver import OptProblem, OptSolver, OptConfig
from test_ipopt import Demo


def main():
    problem = Demo(False)
    guess = np.random.random(problem.nx)
    # first test snopt
    snp_config = OptConfig('snopt', major_iter=100, print_level=1, history=True)
    snp_solver = OptSolver(problem, snp_config)
    snp_rst = snp_solver.solve_guess(guess)
    # then test ipopt
    print(snp_rst.obj, snp_rst.sol, snp_rst.history)


if __name__ == '__main__':
    main()
