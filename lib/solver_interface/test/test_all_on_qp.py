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
from test_ipopt import IpQP, Demo


def main():
    problem = Demo(False)
    guess = np.random.random(problem.nx)
    # first test snopt
    snp_config = OptConfig('snopt', major_iter=100, print_level=1)
    snp_solver = OptSolver(problem, snp_config)
    snp_rst = snp_solver.solve_guess(guess)
    # then test ipopt
    problem = Demo(True)
    ip_config = OptConfig('ipopt', major_iter=100, print_level=5, print_freq=1, linear_solver='ma27')
    ip_solver = OptSolver(problem, ip_config)
    ip_rst = ip_solver.solve_guess(guess)
    # finally scipy
    sp_config = OptConfig('scipy', major_iter=100, print_level=2, fea_tol=1e-4)
    sp_solver = OptSolver(problem, sp_config)
    sp_rst = sp_solver.solve_guess(guess)
    print(snp_rst.obj, ip_rst.obj, sp_rst.obj)
    print(snp_rst.sol, ip_rst.sol, sp_rst.sol)


if __name__ == '__main__':
    main()
