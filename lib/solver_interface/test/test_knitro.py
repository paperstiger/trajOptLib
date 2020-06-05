from pyoptsolver import OptProblem, OptSolver, OptConfig
import numpy as np


class TestFunction2(OptProblem):
    """The same problem but different style"""
    def __init__(self):
        OptProblem.__init__(self, 2, 2, 4)
        self.set_lb([1., 0.])
        self.set_ub([1e20, 1e20])
        self.set_xlb([-1e20, -1e20])
        self.set_xub([0.5, 1e20])
        
    def __cost__(self, x):
        """The eval_f function required by ipopt.

        :param x: a guess/solution of the problem
        :return f: float, objective function

        """
        return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2

    def __gradient__(self, x, g):
        """Evaluation of the gradient of objective function.

        :param x: guess/solution to the problem
        :param g: the gradient of objective function w.r.t x to be written into

        """
        v1 = 200 * (x[1] - x[0] ** 2)
        g[:] = [-2 * x[0] * v1 + 2 * x[0] - 1, v1]
        return True

    def __constraint__(self, x, f):
        """Evaluate constraint function.

        :param x: guess/solution to the problem
        :param f: constraints ready to be written upon
        """
        f[:] = np.array([x[0] * x[1], x[0] + x[1] * x[1]])
        return 1

    def __jacobian__(self, x, g, row, col, rec):
        """Evaluate the Jacobian of the problem.

        :param x: guess/solution to the problem
        :param g: the vector being written on for Jacobian entries
        """
        if rec:
            row[:] = [0, 0, 1, 1]
            col[:] = [0, 1, 0, 1]
        g[:] = [x[1], x[0], 1, 2 * x[1]]
        return 1
            
if __name__ == "__main__":

    
    print('\n\n Test another style\n\n')
    prob = TestFunction2()
    config = OptConfig(backend='knitro')
    print('solve with random guess')
    solver = OptSolver(prob, config)
    rst = solver.solve_rand()
    print(rst.flag, rst.obj, rst.sol)
    print('solve with provided guess')
    rst = solver.solve_guess([0.3, 0.4])
    print(rst.flag, rst.obj, rst.sol)
    print('solve with auto guess')
    rst = solver.solve_guess(None)
    print(rst.flag, rst.obj, rst.sol)

