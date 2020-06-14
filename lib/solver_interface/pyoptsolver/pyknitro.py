"""
Interface for knitro solver, it requires Python interface of knitro being installed.
With this interface, everything happens on Python side so it's quite handy.
"""
import numpy as np
from knitro import *
from pyoptsolver import OptResult, OptProblem


class KnitroSolver(object):
    def __init__(self, prob, config):
        """
        All we need is an instance of the problem in other formats 
        and a  specific configuration.
        """
        self.prob = prob
        self.config = config
        self._grad = np.zeros(self.prob.nx)
        self._y = np.zeros(self.prob.nf)
        nG = self.prob.nG + len(self.prob.Aval)
        self._jac = np.zeros(nG)
        self._row, self._col = np.zeros((2, nG), dtype=int)
        # prepare the solver and define some stuff
        kc = KN_new()
        xIndices = KN_add_vars(kc, prob.nx)
        KN_set_var_lobnds(kc, xLoBnds = prob.xlb)
        KN_set_var_upbnds(kc, xUpBnds = prob.xub)
        # Add the constraints and set the rhs and coefficients
        KN_add_cons(kc, prob.nf)
    
        # Add callback to evaluate nonlinear (non-quadratic) terms in the model:
        #    x0*x1*x2*x3  in the objective 
        #    x0^3         in first constraint c0
        #    x0^2*x3      in second constraint c1
        cb = KN_add_eval_callback (kc, evalObj = True, indexCons = range(prob.nf), funcCallback = self.callbackEvalFC)
        KN_set_con_lobnds(kc, cLoBnds = prob.lb)
        KN_set_con_upbnds(kc, cUpBnds = prob.ub)
        #
        rand_x = np.random.uniform(-np.ones(prob.nx), np.ones(prob.nx))
        rand_x = np.clip(rand_x, prob.xlb, prob.xub)
        self.prob.__jacobian__(rand_x, self._jac, self._row, self._col, True)
        # Set obj. gradient and nonlinear jac provided through callbacks.
        # Mark objective gradient as dense, and provide non-zero sparsity
        # structure for constraint Jacobian terms.
        cbjacIndexCons = self._row
        cbjacIndexVars = self._col
        KN_set_cb_grad (kc, cb, objGradIndexVars = KN_DENSE, jacIndexCons = cbjacIndexCons, 
                        jacIndexVars=cbjacIndexVars, gradCallback=self.callbackEvalGA)
        # Hessian is skipped
        # New point callback is skipped as well

        # Set option to print output after every iteration.
        # Just set some default values, may improve later
        # default options from Zherong
        use_direct = config.get('use_direct', False)
        KN_set_int_param (kc, KN_PARAM_ALGORITHM, KN_ALG_BAR_DIRECT if use_direct else KN_ALG_BAR_CG)
        KN_set_int_param (kc, KN_PARAM_BAR_FEASIBLE, KN_BAR_FEASIBLE_GET_STAY)
        KN_set_int_param (kc, KN_PARAM_CG_MAXIT, 50)
        KN_set_int_param (kc, KN_PARAM_OUTLEV, KN_OUTLEV_ITER)
        # set other options
        exclude_keys = {'use_direct', 'history', 'user_callback'}

        def get_key(key):
            if isinstance(key, int):
                return key
            elif isinstance(key, 'float'):
                if key.isnumeric():
                    return int(key)
                else:
                    # find the string at first
                    key = key.upper()
                    if not key.startswith('KN_PARAM'):
                        key = 'KN_PARAM_' + key
                    try:
                        return eval(key)
                    except:
                        raise Exception("Option %s is not defined in knitro" % key)
        
        def get_item(item):
            """Try to convert val to correct value"""
            if isinstance(item, (int, float)):
                return item
            elif isinstance(item, str):
                try:
                    uitem = item.upper()
                    if not uitem.startswith('KN_'):
                        uitem = 'KN_' + uitem
                    return eval(uitem)
                except:
                    return item

        for key, item in config.items():
            if key in exclude_keys:
                continue
            key = get_key(key)
            item = get_item(item)
            if isinstance(item, int):
                KN_set_int_param(kc, key, item)
            elif isinstance(item, float):
                KN_set_double_param(kc, key, item)
            elif isinstance(item, str):
                KN_set_char_param(kc, key, item)
        # default_options = {KN_PARAM_ALGORITHM: KN_
        # enable history
        if config.get('history', False):
            KN_set_newpt_callback(kc, self.callbackHistory)
            self.history = []
        if 'user_callback' in config:
            KN_set_newpt_callback(kc, config['user_callback'])

        # cache some results
        self.kc = kc
        self.xIndices = xIndices

    def callbackEvalFC(self, kc, cb, evalRequest, evalResult, userParams):
        if evalRequest.type != KN_RC_EVALFC:
            print ("*** callbackEvalFC incorrectly called with eval type %d" % evalRequest.type)
            return -1
        x = np.array(evalRequest.x)
        # Evaluate nonlinear term in objective
        evalResult.obj = self.prob.__cost__(x)
        # Evaluate nonlinear terms in constraints
        self._y[:] = 0
        self.prob.__constraint__(x, self._y)
        evalResult.c[:] = self._y
        return 0

    def callbackHistory(self, kc, x, lambda_, userParams):
        self.history.append(x.copy())
        return 0

    def callbackEvalGA(self, kc, cb, evalRequest, evalResult, userParams):
        if evalRequest.type != KN_RC_EVALGA:
            print ("*** callbackEvalGA incorrectly called with eval type %d" % evalRequest.type)
            return -1
        x = np.array(evalRequest.x)

        # Evaluate nonlinear term in objective gradient
        self._grad[:] = 0
        self.prob.__gradient__(x, self._grad)
        evalResult.objGrad[:] = self._grad

        # Evaluate nonlinear terms in constraint gradients (Jacobian)
        self.prob.__jacobian__(x, self._jac, self._row, self._col, False)
        evalResult.jac[:] = self._jac
        return 0

    def callbackEvalH(self, kc, cb, evalRequest, evalResult, userParams):
        raise NotImplementedError("Hessian is not implemented by default")

    # def callbackNewPoint(self, kc, x, lambda_, userParams):
    #     return 0
    def solve_rand(self):
        guess = np.random.random(self.prob.nx)
        guess = np.clip(guess, self.prob.xlb, self.prob.xub)
        return self.solve_guess(guess)

    def solve_guess(self, guess):
        if guess is not None:
            for x, v in zip(self.xIndices, guess):
                KN_set_var_primal_init_values(self.kc, x, v)
        KN_solve(self.kc)
        return self.get_solution()

    def get_solution(self):
        nStatus, objSol, x, lambda_ =  KN_get_solution(self.kc)
        print ("  feasibility violation    = %e" % KN_get_abs_feas_error(self.kc))
        print ("  KKT optimality violation = %e" % KN_get_abs_opt_error(self.kc))
        result = OptResult()
        result.flag = nStatus
        result.obj = objSol
        result.sol = x
        result.lmd = lambda_
        if hasattr(self, 'history'):
            result.history = np.array(self.history)
        return result

    def __del__(self):
        KN_free(self.kc)
