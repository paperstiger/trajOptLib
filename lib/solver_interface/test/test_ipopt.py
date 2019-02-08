#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 motion <motion@motion-MS-7971>
#
# Distributed under terms of the MIT license.

"""
Test ipopt solver here.
"""
import numpy as np
from pyoptsolver import OptProblem as IpoptProblem, solve_problem, IpoptConfig, IpoptSolver, set_verbosity


class BasicQP(IpoptProblem):
    def __init__(self):
        IpoptProblem.__init__(self, 2, 3, 4)
        self.batch_set_xlb(np.zeros(2), 0)
        self.batch_set_xub(2*np.ones(2), 0)
        lb = self.get_lb()
        ub = self.get_ub()
        lb[:] = [-1e20, 1.0, 1.0]
        ub[:] = [1e20, 3.0, 5.0]
        self.set_a_by_triplet(np.ones(2), [2, 2], [0, 1])

    def __callg__(self, x, F, G, row, col, rec, needg):
        x1, x2 = x
        F[0] = np.sum(x ** 2)
        F[1] = x1 * x2
        F[2] = 0
        if needg:
            G[:] = [2 * x1, 2 * x2, x2, x1]
            if rec:
                row[:] = [0, 0, 1, 1]
                col[:] = [0, 1, 0, 1]
        return 3, 4


class IpQP(IpoptProblem):
    def __init__(self):
        IpoptProblem.__init__(self, 2, 2, 4)
        self.set_xlb(np.zeros(2))
        self.set_xub(2*np.ones(2))
        self.set_lb([1.0, 1.0])
        self.set_ub([3.0, 5.0])
        self.ipopt_style()

    def __cost__(self, x):
        num = np.sum(x ** 2)
        return float(num)

    def __gradient__(self, x, g):
        g[:] = 2 * x
        return True

    def __constraint__(self, x, g):
        x1, x2 = x
        g[0] = x1 * x2
        g[1] = x1 + x2
        return 2

    def __jacobian__(self, x, g, row, col, rec):
        x1, x2 = x
        if rec:
            row[:] = [0, 0, 1, 1]
            col[:] = [0, 1, 0, 1]
            return 4
        else:
            g[:] = [x2, x1, 1.0, 1.0]
            return 4


class Demo(IpoptProblem):
    def __init__(self, ip=False):
        self.ip = ip
        if ip:
            IpoptProblem.__init__(self, 4, 2, 8)
            self.ipopt_style()
            self.set_lb([25, 40.])
            self.set_ub([2e19, 40.0])
        else:
            IpoptProblem.__init__(self, 4, 3, 12)
            self.set_lb([-2e19, 25, 40.])
            self.set_ub([2e19, 2e19, 40.0])
        self.set_xlb(np.ones(4))
        self.set_xub(5*np.ones(4))

    def __callg__(self, x, f, g, row, col, rec, needg):
        f[0] = self.eval_f(x)
        self.eval_constr(x, f[1:])
        if needg:
            self.eval_gradient(x, g[:4])
            self.eval_jacobian(x, g[4:], [0], [0], False)
            if rec:
                row[:], col[:] = np.unravel_index(np.arange(12), (3, 4))
        return 3, 12

    def eval_f(self, x):
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

    def eval_gradient(self, x, grad_f):
        grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
        grad_f[1] = x[0] * x[3];
        grad_f[2] = x[0] * x[3] + 1;
        grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
        return True

    def eval_constr(self, x, g):
        g[0] = x[0] * x[1] * x[2] * x[3];
        g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
        return 2

    def eval_jacobian(self, x, values, iRow, jCol, rec):
        if rec:
            iRow[:] = [0, 0, 0, 0, 1, 1, 1, 1]
            jCol[:] = [0, 1, 2, 3, 0, 1, 2, 3]
        else:
            values[0] = x[1]*x[2]*x[3];
            values[1] = x[0]*x[2]*x[3];
            values[2] = x[0]*x[1]*x[3];
            values[3] = x[0]*x[1]*x[2];
            values[4] = 2*x[0];
            values[5] = 2*x[1];
            values[6] = 2*x[2];
            values[7] = 2*x[3];
        return 8

    def eval_hessian(self, x, obj_factor, lmd, values, iRow, jCol, rec):
        if rec:
            idx = 0
            for row in range(4):
                for col in range(0, row + 1):
                    iRow[idx] = row
                    jCol[idx] = col
                    idx += 1
        else:
            values[0] = obj_factor * (2*x[3]);
            values[1] = obj_factor * (x[3]);
            values[2] = 0;
            values[3] = obj_factor * (x[3]);
            values[4] = 0; # 2,1
            values[5] = 0; # 2,2
            values[6] = obj_factor * (2*x[0] + x[1] + x[2]); # 3,0
            values[7] = obj_factor * (x[0]); # 3,1
            values[8] = obj_factor * (x[0]); # 3,2
            values[9] = 0; # 3,3

            values[1] += lmd[0] * (x[2] * x[3]); # 1,0
            values[3] += lmd[0] * (x[1] * x[3]); # 2,0
            values[4] += lmd[0] * (x[0] * x[3]); # 2,1
            values[6] += lmd[0] * (x[1] * x[2]); # 3,0
            values[7] += lmd[0] * (x[0] * x[2]); # 3,1
            values[8] += lmd[0] * (x[0] * x[1]); # 3,2

            values[0] += lmd[1] * 2; # 0,0
            values[2] += lmd[1] * 2; # 1,1
            values[5] += lmd[1] * 2; # 2,2
            values[9] += lmd[1] * 2; # 3,3
        return 10


def main():
    use_ip = False
    prob = Demo(use_ip)
    x0 = np.array([1.0, 5.0, 5.0, 1.0])
    config = IpoptConfig()
    # config.hessian_approximation = 'exact'
    solver = IpoptSolver(prob, config)
    result = solver.solve_guess(x0)
    print(result.sol)


def main2():
    prob = BasicQP()
    # prob = IpQP()
    x0 = np.ones(2)
    config = IpoptConfig()
    config.print_level = 5
    # config.linear_solver = 'ma57'
    # result = solve_problem(prob, config, x0)
    solver = IpoptSolver(prob, config)
    result = solver.solve_guess(x0)
    print(result.flag, result.obj, result.sol)


if __name__ == '__main__':
    main2()
