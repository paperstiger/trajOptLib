"""
Test the discrete solver. I shall just try some simple system
"""
import numpy as np
from pyoptsolver import OptConfig, OptSolver
from trajoptlib import TrajOptDisProblem, LqrObj, System, inf_, zero_, one_
from trajoptlib.io import get_onoff_args
from trajoptlib.utility import show_sol


def main():
    args = get_onoff_args('linear', 'ad', 'pen', 'fd')
    if args.linear:
        test_linear(args.ad)
    if args.pen:
        test_pendulum(args.ad, args.fd)

class LinearSystem(System):
    """The second order problem..."""
    def __init__(self):
        System.__init__(self, 4, 2, 0, 'Dis')
        self.h = 0.2

    def jac_dyn(self, t, x, u, p=None):
        """Return value at next step"""
        dx = np.concatenate((x[2:], u))
        xkp = x + dx * self.h
        jac = np.zeros((4, 6))
        np.fill_diagonal(jac, 1)
        np.fill_diagonal(jac[:2, 2:4], self.h)
        np.fill_diagonal(jac[2:4, 4:6], self.h)
        return xkp, jac


class AdLinearSystem(System):
    """The second order problem..."""
    def __init__(self):
        System.__init__(self, 4, 2, 0, 'Dis')
        self.h = 0.2

    def __ad__(self, xin):
        """Return value at next step"""
        import autograd.numpy as np
        return np.array([xin[0] + xin[2] * self.h, xin[1] + xin[3] * self.h, xin[2] + self.h * xin[4], xin[3] + self.h * xin[5]])


class Pendulum(System):
    """Test pendulum nonlinearity."""
    def __init__(self):
        System.__init__(self, 2, 1, 0, 'Dis')
        self.h = 0.2

    def jac_dyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0] / 2 - 2 * np.sin(x[0])
        y = np.array([x[0] + y1 * self.h, x[1] + y2 * self.h])
        J = np.array([[1, self.h, 0], [self.h * (-2 * np.cos(x[0])), 1, self.h / 2]])
        return y, J


class AdPendulum(System):
    def __init__(self):
        System.__init__(self, 2, 1, 0, 'Dis')
        self.h = 0.2

    def __ad__(self, xin):
        import autograd.numpy as np
        return np.array([xin[0] + xin[1] * self.h, xin[1] + self.h * (xin[2] - np.sin(xin[0]))])


class FdPendulum(System):
    def __init__(self):
        System.__init__(self, 2, 1, 0, 'Dis')
        self.h = 0.2

    def dyn(self, t, x, u, p):
        return np.array([x[0] + self.h * x[1], x[1] + self.h * (u[0] / 2 - 2 * np.sin(x[0]))])


def test_pendulum(ad, fd):
    print("Test on pendulum problem")
    if ad:
        pen_sys = AdPendulum()
    elif fd:
        pen_sys = FdPendulum()
    else:
        pen_sys = Pendulum()

    N = 40
    prob = TrajOptDisProblem(pen_sys, N)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1.0]), np.array([1.0])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    delta = 0.2
    prob.xfbd = [np.array([np.pi - delta, -delta]), np.array([np.pi + delta, delta])]
    lqr = LqrObj(R=np.ones(1))
    prob.add_lqr_obj(lqr)
    prob.pre_process()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig('snopt', print_level=1, deriv_check=0, print_file='tmp.out', major_iter=200)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print('Solve flag is ', rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


def test_linear(ad):
    print("Test on linear system")
    if ad:
        lin_sys = AdLinearSystem()
    else:
        lin_sys = LinearSystem()

    N = 30
    prob = TrajOptDisProblem(lin_sys, N)
    lqr = LqrObj(R=np.ones(2))
    prob.add_lqr_obj(lqr)
    prob.xbd = [-inf_[4], inf_[4]]
    prob.ubd = [-one_[2], one_[2]]
    prob.x0bd = [0.5 * one_[4], 0.5 * one_[4]]
    prob.xfbd = [-0. * one_[4], 0. * one_[4]]
    prob.pre_process()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig('snopt', print_level=1, major_iter=100, deriv_check=0, print_file='tmp.out')
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


if __name__ == '__main__':
    main()
