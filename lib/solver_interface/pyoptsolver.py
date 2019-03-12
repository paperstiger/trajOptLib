from pyoptsolvercpp import IpoptConfig, IpoptSolver
from pyoptsolvercpp import SnoptConfig, SnoptSolver
from pyoptsolvercpp import OptProblem, OptResult
from pyoptsolvercpp import __with_snopt__, __with_ipopt__, __version__
import warnings
import six
import numpy as np
from scipy.optimize import NonlinearConstraint, minimize, Bounds
from scipy.sparse import coo_matrix


class OptConfig(object):
    """This class is used to keep track of various optimization configurations.

    It provides a unified interface for interacting with those solvers.
    Some custom options are also supported, but might not work depending on the backend.
    """
    shared_options = {'major_iter', 'opt_tol', 'fea_tol', 'print_level', 'deriv_check'}
    snopt_options = {'minor_iter', 'iter_limit'}
    ipopt_options = {'print_freq', 'linear_solver', 'exact_hessian'}
    scipy_kws = {'tol'}
    scipy_option_kws = {'grad', 'xtol', 'gtol', 'barrier_tol', 'initial_constr_penalty',
                        'initial_tr_radius', 'initial_barrier_parameter', 'initial_barrier_tolerance',
                        'factorization_method', 'disp'}
    def __init__(self, backend='ipopt', **kw):
        if backend not in ['ipopt', 'snopt', 'scipy']:
            warnings.warn('Backend should be in ipopt, snopt, scipy')
        if __with_ipopt__ and backend == 'ipopt':
            self.backend = 'ipopt'
            self.option = IpoptConfig()
        elif __with_snopt__ and backend == 'snopt':
            self.backend = 'snopt'
            self.option = SnoptConfig()
        else:
            self.backend = 'scipy'
            self.option = {'options': {'sparse_jacobian': True}, 'tol': 1e-6}
        self._process_kw(**kw)

    def _process_kw(self, **kws):
        is_ipopt = self.backend == 'ipopt'
        is_snopt = self.backend == 'snopt'
        is_scipy = self.backend == 'scipy'
        if is_scipy:
            options = self.option['options']
        for key, val in six.iteritems(kws):
            if key == 'major_iter':
                if not is_scipy:
                    self.option.set_major_iter(val)
                else:
                    options['maxiter'] = val
            elif key == 'opt_tol':
                if not is_scipy:
                    self.option.set_opt_tol(val)
                else:
                    options['gtol'] = val
            elif key == 'fea_tol':
                if not is_scipy:
                    self.option.set_fea_tol(val)
                else:
                    self.option['tol'] = val
            elif key == 'print_level':
                if not is_scipy:
                    self.option.set_print_level(val)
                else:
                    options['verbose'] = val
            elif key == 'deriv_check':
                if not is_scipy:
                    self.option.enable_deriv_check(val)
            elif key == 'minor_iter':
                if is_snopt:
                    self.option.set_minor_iter(val)
            elif key == 'iter_limit':
                if is_snopt:
                    self.option.set_iter_limit(val)
            elif key == 'print_freq':
                if is_ipopt:
                    self.option.set_print_freq(val)
            elif key == 'linear_solver':
                if is_ipopt:
                    self.option.set_linear_solver(val)
            elif key == 'exact_hessian':
                if is_ipopt and val:
                    self.option.enable_exact_hessian()
            elif key in self.scipy_kws:
                self.option[key] = val
            elif key in self.scipy_option_kws:
                options[key] = val
            else:  # considering we are adding several types
                if isinstance(val, int):
                    if not is_scipy:
                        self.option.add_int_option(key, val)
                elif isinstance(val, float):
                    if not is_scipy:
                        self.option.add_float_option(val)
                elif isinstance(val, str) or val is None:
                    if is_snopt:
                        self.option.add_string_option(key)
                    elif is_ipopt:
                        self.option.add_string_option(key, val)


class TrustConstrSolver(object):
    """A wrapper that builds on a already declared problem."""
    def __init__(self, problem, option):
        self.problem = problem
        self.option = option
        self.g = np.zeros(problem.nG)
        self.y = np.zeros(problem.nf)
        self.row, self.col = np.zeros((2, problem.nG), dtype=int)

    def _get_coo(self):
        return coo_matrix((self.g, (self.row, self.col)), shape=(self.problem.nf, self.problem.nx))

    def solve_rand(self):
        guess = self.problem.random_gen_x()
        return self.solve_guess(guess)

    def solve_guess(self, guess):
        if self.problem.ipstyle:
            def confun(x):
                self.problem.eval_constr(x, self.y)
                return self.y

            def jacfun(x):
                self.problem.__jacobian__(x, self.g, self.row, self.col, False)
                return self._get_coo()

            def costfun(x):
                g = np.zeros_like(x)
                self.problem.eval_gradient(x, g)
                return self.problem.__cost__(x), g
            # figure out the jacobian sparsity
            self.problem.__jacobian__(guess, self.g, self.row, self.col, True)
        else:
            def confun(x):
                if self.problem.grad:
                    self.problem.__callg__(x, self.y, self.g, self.row, self.col, False, False)
                else:
                    self.problem.__callf__(x, self.y)
                return self.y

            def jacfun(x):
                self.problem.__callg__(x, self.y, self.g, self.row, self.col, False, True)
                return self._get_coo()

            def costfun(x):
                self.problem.__callg__(x, self.y, self.g, self.row, self.col, False, True)
                if not hasattr(self, 'cost_row_mask'):
                    self.cost_row_mask = np.where(self.row == 0)[0]
                grad = np.zeros_like(x)
                grad[self.col[self.cost_row_mask]] = self.g[self.cost_row_mask]
                return self.y[0], grad
            # figure out the jacobian sparsity
            self.problem.__callg__(guess, self.y, self.g, self.row, self.col, True, True)
        # jac0 = self._get_coo()

        constr = NonlinearConstraint(confun, self.problem.get_lb(), self.problem.get_ub(), jac=jacfun)
        bounds = Bounds(self.problem.get_xlb(), self.problem.get_xub())
        res = minimize(costfun, guess, method='trust-constr', jac=True, bounds=bounds,
                       constraints=constr, **self.option)
        # use the attr trick to add contents to returned res
        setattr(res, 'obj', res.fun)
        setattr(res, 'flag', res.status != 0 and res.constr_violation < self.option['tol'])
        setattr(res, 'sol', res.x)
        setattr(res, 'fval', res.constr)
        setattr(res, 'lmd', res.v)
        setattr(res, 'xmul', res.v)
        return res


class OptSolver(object):
    def __init__(self, problem, config):
        if config.backend == 'snopt':
            self.solver = SnoptSolver(problem, config.option)
        elif config.backend == 'ipopt':
            self.solver = IpoptSolver(problem, config.option)
        else:
            self.solver = TrustConstrSolver(problem, config.option)

    def solve_rand(self):
        return self.solver.solve_rand()

    def solve_guess(self, guess):
        return self.solver.solve_guess(guess)
