/*
 * solverwrapper.cpp
 * Copyright (C) 2017 Gao <gao.tang@duke.edu>
 *
 * Distributed under terms of the  license.
 */

/* Write a python wrapper for the problem */

#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#ifdef SNOPT
#include "snoptWrapper.h"
#include "funcStyle.h"
#include "classStyle.h"
#endif
#ifdef IPOPT
#include "ipoptWrapper.h"
#endif


#define VERSION_INFO "0.0.1"

double PYOPTSOLVER_FD_STEP = 1e-6;


class pyProbFun: public ProblemFun{
    public:
        using ProblemFun::ProblemFun;

        int operator()(cRefV x, RefV F) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "__callf__",
                    operator(),
                    x, F
                    );
        }

        pint operator()(cRefV x, RefV F, RefV G, RefVl row, RefVl col, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_NAME(
                    pint,
                    ProblemFun,
                    "__callg__",
                    operator(),
                    x, F, G, row, col, rec, needg
                    );
        }

#ifdef ENABLEIP
        double evalF(cRefV x) override{
            PYBIND11_OVERLOAD_NAME(
                    double,
                    ProblemFun,
                    "__cost__",
                    evalF,
                    x
                    );
        }

        bool evalGrad(cRefV x, RefV grad) override{
            PYBIND11_OVERLOAD_NAME(
                    bool,
                    ProblemFun,
                    "__gradient__",
                    evalGrad,
                    x, grad
                    );
        }

        int evalG(cRefV x, RefV g) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "__constraint__",
                    evalG,
                    x, g
                    );
        }

        int evalJac(cRefV x, RefV g, RefVl row, RefVl col, bool rec) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "__jacobian__",
                    evalJac,
                    x, g, row, col, rec
                    );
        }

        int evalHess(cRefV x, double sigma, cRefV lmd, RefV g, RefVl row, RefVl col, bool rec) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "__hessian__",
                    evalHess,
                    x, sigma, lmd, g, row, col, rec
                    );
        }
#endif
};


namespace py = pybind11;
PYBIND11_MODULE(pyoptsolvercpp, m){
    m.doc() = R"pbdoc(
        A Python interface for large scale optimization solver SNOPT and Ipopt
    )pbdoc";

    // this class was called result previously, consider changing __init__.py
    py::class_<optResult>(m, "OptResult", R"pbdoc(
        A class that contains the solution returned by either one.

        This class has the default constructor.

        Attributes:
            flag (int): the return flag by SNOPT.solve()
            obj (float): the cost function
            sol (ndarray): a copy of the solution
            fval (ndarray): a copy of F
            lmd (ndarray): a copy of Lagrangian multipliers
            xmul (ndarray): a copy of multipliers associated with bounds on x
    )pbdoc")
        .def(py::init<>(), R"pbdoc(
            Constructor for SnoptResult.

            It takes no arguments.
        )pbdoc")
        .def("get_sol", &optResult::get_sol, R"pbdoc(
            Return a reference to the solution without copying.
        )pbdoc")
        .def("get_fval", &optResult::get_fval, R"pbdoc(
            Return a reference to the function evaluation without copying.
        )pbdoc")
        .def("get_lambda", &optResult::get_lambda, R"pbdoc(
            Return a reference to the Lagrangian multiplier without copying.
        )pbdoc")
        .def("get_xmul", &optResult::get_xmul, R"pbdoc(
            Return a reference to the Lagrangian multiplier on bounds without copying.
        )pbdoc")
        .def_readonly("flag", &optResult::flag)
        .def_readonly("obj", &optResult::val)
        .def_readonly("sol", &optResult::sol)
        .def_readonly("fval", &optResult::c)
        .def_readonly("xmul", &optResult::xmul)
        .def_readonly("lmd", &optResult::lmd);


    py::class_<ProblemFun, pyProbFun>(m, "OptProblem", R"pbdoc(
        A class sub-classed by users to define a custom optimization problem.

        :currentmodule:
        .. automethod:: __callf__.
        .. automethod:: __callg__.
        .. automethod:: __cost__.
        .. automethod:: __gradient__.
        .. automethod:: __constraint__
        .. automethod:: __jacobian__.
        .. automethod:: __hessian__.

        This class defines an optimization problem to optimize :math:`y(0)` where :math:`y=f(x)+Ax`
        subject to constraints :math:`lb \le y \le ub` and :math:`xlb \le x \le xub`.

        Attributes:
            nx (int): number of variables to be optimized. The same with :class:`.FunBase`.
            nf (int): number of constraints. The same with :class:`.FunBase`.
            nG (int=0): nnz of the Jacobian matrix. The same with :class:`.FunBase`.
            grad (bool=False): if the function __callg__ has been implemented so we enable gradient.
            lb (ndarray): lower bound for f. It is readonly, use get_lb() to get a reference.
            ub (ndarray): upper bound for f. It is readonly, use get_ub() to get a reference.
            xlb (ndarray): lower bound for x. It is readonly, use get_xlb() to get a reference.
            xub (ndarray): upper bound for x. It is readonly, use get_xub() to get a reference.
            Aval (ndarray): the value ndarray for triplet A, use get_aval() to get a reference.
            Arow (ndarray): integer ndarray of rows for triplet A, use get_arow() to get a reference.
            Acol (ndarray): integer ndarray of columns for triplet A, use get_acol() to get a reference.

        Note:
            This class has 2 constructors. Refer to :class:`.FunBase` for details.
            The first entry of lb and ub should be -np.inf and np.inf, respectively.
            For linear function part, you have to set function return 0.
            Use ipopt_style() and snopt_style() to toggle function evaluation mode.
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<int, int, int>())
        .def("update_nf", &ProblemFun::updateNf, R"pbdoc(
            Update nf and reallocate space for lb and ub, if necessary.

            Args:
                nf (int): the desired nf
            )pbdoc")
        .def("update_ng", &ProblemFun::updateNg, R"pbdoc(
            Update nnz value, equivalently to directly writing to nG.

            Args:
                ng (int): the desired nG
            )pbdoc")
        .def("set_a_by_matrix", (void (ProblemFun::*)(crRefM)) &pyProbFun::setA, R"pbdoc(
            Set triplet matrix A by a dense row major matrix, i.e. 2d ndarray of size (m, n).

            Args:
                mat (ndarray): the 2d float A matrix in row major.
        )pbdoc")
        .def("set_a_by_triplet", (void (ProblemFun::*)(cRefV, cRefVi, cRefVi)) &pyProbFun::setA, R"pbdoc(
            Set triplet matrix A by the value, row, column triplet pairs.

            Args:
                val (ndarray): the value ndarray
                row (ndarray): the integer row ndarray
                col (ndarray): thet integer column ndarray
        )pbdoc")
        .def("set_a_by_triplet", (void (ProblemFun::*)(cRefV, cRefVl, cRefVl)) &pyProbFun::setA, R"pbdoc(
            Set triplet matrix A by the value, row, column triplet pairs.

            Args:
                val (ndarray): the value ndarray
                row (ndarray): the integer row ndarray
                col (ndarray): thet integer column ndarray
        )pbdoc")
        .def("check_overlap", &pyProbFun::overlap_check, R"pbdoc(
            Check if user provided linear and nonlinear Jacobian has overlap.
        )pbdoc")
        .def("detect_prob_size", (void (ProblemFun::*)(int, int)) &pyProbFun::detect_prob_size,
                R"pbdoc(
            Automatically detect nf and nG.

            This function requires __callf__ or __callg__ return non-trivial values.

            Args:
                nf (int): an estimated size for :math:`y`, it has to be larger.
                nG (int): an estimated size for Jacobian, it has to be larger. If set to 0, __callf__ is used and
                has to be defined. Otherwise __callg__ has to be defined.
                )pbdoc")
        .def("detect_prob_size", (void (ProblemFun::*)(cRefV, int, int)) &pyProbFun::detect_prob_size,
                R"pbdoc(
            Automatically detect nf and nG.

            This function requires __callf__ or __callg__ return non-trivial values.

            Args:
                x( ndarray): the guess being used. This will change nf too so please make sure its size is correct
                nf (int): an estimated size for :math:`y`, it has to be larger.
                nG (int): an estimated size for Jacobian, it has to be larger. If set to 0, __callf__ is used and
                has to be defined. Otherwise __callg__ has to be defined.
                )pbdoc")
        .def("enable_timer", &pyProbFun::enable_time_record, R"pbdoc(
            Turn on internal timer. By doing so, user is able to call get_timer function.
        )pbdoc")
        .def("set_debug_verbose", &pyProbFun::set_debug_verbose, R"pbdoc(
            Turn on or off verbosity.
        )pbdoc")
        .def("get_timer", &pyProbFun::get_time_obj_history, R"pbdoc(
            Return a tuple of time, cost, constr
        )pbdoc")
        .def("get_lb", &pyProbFun::get_lb)
        .def("get_ub", &pyProbFun::get_ub)
        .def("get_xlb", &pyProbFun::get_xlb)
        .def("get_xub", &pyProbFun::get_xub)
        .def("get_aval", &pyProbFun::get_aval)
        .def("get_arow", &pyProbFun::get_arow)
        .def("get_acol", &pyProbFun::get_acol)
        .def("set_lb", &pyProbFun::set_lb)
        .def("set_ub", &pyProbFun::set_ub)
        .def("set_xlb", &pyProbFun::set_xlb)
        .def("set_xub", &pyProbFun::set_xub)
        .def("batch_set_lb", &pyProbFun::batchSetLb, R"pydoc(
            Set lb in batch mode by specifying an ndarray and starting index.

            This value is not necessary since you can get a reference by get_lb and directly operate on that.
            We keep it here for backward compatibility.
            Other similar functions are not documented. Please refer to this.

            Args:
                lb (ndarray): a segment of lower bound to be set
                ind0 (int): the starting index for this segment.
        )pydoc")
        .def("batch_set_ub", &pyProbFun::batchSetUb)
        .def("batch_set_xlb", &pyProbFun::batchSetXlb)
        .def("batch_set_xub", &pyProbFun::batchSetXub)
        .def("random_gen_x", &pyProbFun::randomGenX, R"pydoc(
            Randomly generate an initial guess of x within lb and ub bounds.

            Returns:
                ndarray: the generated initial guess.
        )pydoc")
        .def("eval_f", &pyProbFun::eval_f, R"pydoc(
            Evaluate __callf__ once.

            This serves for debugging purpose.

            Args:
                x (ndarray): a candidate solution
            Returns:
                y (ndarray): the evaluated function.
        )pydoc")
        .def("eval_g", &pyProbFun::eval_g, R"pydoc(
            Evaluate __callf__ once.

            This serves for debugging purpose.

            Args:
                x (ndarray): a candidate solution
            Returns:
                y (ndarray): the evaluated function.
                G (ndarray): the value part of the returned Jacobian triplet.
                row (ndarray): the row part of the returned Jacobian triplet.
                col (ndarray): the column part of the returned Jacobian triplet.
        )pydoc")
        .def("__callf__", (int (ProblemFun::*)(cRefV, RefV)) &pyProbFun::operator())
        .def("__callg__", (std::pair<int, int> (ProblemFun::*)(cRefV, RefV, RefV, RefVl, RefVl, bool, bool)) &pyProbFun::operator())
#ifdef ENABLEIP
        .def("eval_cost", [](ProblemFun* self, cRefV x) {return self->evalF(x);})
        .def("eval_gradient", [](ProblemFun* self, cRefV x) {VX grad(self->nf); self->evalGrad(x, grad); return grad;})
        .def("eval_constr", [](ProblemFun* self, cRefV x) {VX grad(self->nf); self->evalGrad(x, grad); return grad;})
        .def("eval_jacobian", [](ProblemFun *self, cRefV x) {
              int nG = self->nG;
              VX g(nG);
              VXl row(nG), col(nG);
              self->evalJac(x, g, row, col, true);
              return std::make_tuple(g, row, col);
            })
        .def("__cost__", &ProblemFun::evalF)
        .def("__constraint__", &ProblemFun::evalG)
        .def("__gradient__", &ProblemFun::evalGrad)
        .def("__jacobian__", &ProblemFun::evalJac)
        .def("__hessian__", &ProblemFun::evalHess)
#endif
        .def("ipopt_style", [](ProblemFun *self) {self->ipStyle=true;})
        .def("snopt_style", [](ProblemFun *self) {self->ipStyle=false;})
        .def_readonly("ipstyle", &pyProbFun::ipStyle)
        .def_readwrite("nx", &pyProbFun::nx)
        .def_readwrite("nf", &pyProbFun::nf)
        .def_readwrite("nG", &pyProbFun::nG)
        .def_readwrite("nH", &pyProbFun::hess_nnz)
        .def_readwrite("grad", &pyProbFun::grad)
        .def_readwrite("Aval", &pyProbFun::Aval)
        .def_readwrite("Arow", &pyProbFun::Arow)
        .def_readwrite("Acol", &pyProbFun::Acol)
        .def_readwrite("lb", &pyProbFun::lb)
        .def_readwrite("ub", &pyProbFun::ub)
        .def_readwrite("xlb", &pyProbFun::xlb)
        .def_readwrite("xub", &pyProbFun::xub);


#ifdef SNOPT
    // this class was called snoptConfig previously, consider changing __init__.py
    py::class_<snoptConfig>(m, "SnoptConfig", R"pbdoc(
        A class that specifies configuration to snopt solver, keep for backward compatibility.

        This class has the default constructor.

        Attributes:
            name (str="Toy"): the name for the problem.
            print_file (str): the filename for printing results.
            print_level (int=0): the print level for major iteration; 0 disable iterations; 1 enables iterations.
            minor_print_level (int=0): print level for minor iteration.
            verify_level (int=0): level for gradient verification; 0 disable; 3 full verify.
            major_iter_limit (int=0): major iteration limit
            minor_iter_limit (int=0): minor iteration limit
            iter_limit (int=0): total iteration limit
            opt_tol (float=1e-6): tolerance for optimality
            fea_tol (float=1e-6): tolerance for feasibility
            int_options (list of (string, int) tuple): integer options
            float_options (list of (string, float) tuple): float options
            string_options (list of (string, string) tuple): string options
    )pbdoc")
        .def(py::init<>())
        .def("add_int_option", &snoptConfig::addIntOption, R"pbdoc(
            Add an integer option to configurations.

            Args:
                option (str): the option to configure
                value (int): the integer value for configuration
        )pbdoc")
        .def("add_float_option", &snoptConfig::addFloatOption, R"pbdoc(
            Add a float option to configurations.

            Args:
                option (str): the option to configure
                value (float): the float value for configuration
        )pbdoc")
        .def("add_string_option", &snoptConfig::addStringOption, R"pbdoc(
            Add a string option to configurations.

            Args:
                option (str): the option to configure
        )pbdoc")
        .def("set_major_iter", &snoptConfig::setMajorIter, R"pbdoc(
            Set up the major iteration limit.

            Args:
                iter (int): limit for major iteration 
        )pbdoc")
        .def("set_minor_iter", &snoptConfig::setMinorIter, R"pbdoc(
            Set up the minor iteration limit.

            Args:
                iter (int): limit for minor iteration 
        )pbdoc")
        .def("set_iter_limit", &snoptConfig::setIterLimit, R"pbdoc(
            Set up the total iteration limit.

            Args:
                iter (int): limit for total iteration 
        )pbdoc")
        .def("set_opt_tol", &snoptConfig::setOptTol, R"pbdoc(
            Set up the optimization tolerance.

            Args:
                tol (float): the tolerance
        )pbdoc")
        .def("set_fea_tol", &snoptConfig::setFeaTol, R"pbdoc(
            Set up the feasibility tolerance.

            Args:
                tol (float): the tolerance
        )pbdoc")
        .def("set_print_level", &snoptConfig::setPrintLevel, R"pbdoc(
            Set up the print level.

            Args:
                level (int): the print level
        )pbdoc")
        .def("enable_deriv_check", &snoptConfig::enableDerivCheck, py::arg("lvl")=3, R"pbdoc(
            Set up the derivative check param.

            Args:
                level (int): the derivative check level
        )pbdoc")
        .def("addIntOption", &snoptConfig::addIntOption)
        .def("addFloatOption", &snoptConfig::addFloatOption)
        .def("addStringOption", &snoptConfig::addStringOption)
        .def_readwrite("name", &snoptConfig::name)
        .def_readwrite("print_file", &snoptConfig::printFile)
        .def_readwrite("print_level", &snoptConfig::printlevel)
        .def_readwrite("minor_print_level", &snoptConfig::minorprintlevel)
        .def_readwrite("verify_level", &snoptConfig::verifylevel)
        .def_readwrite("major_iter_limit", &snoptConfig::majoriterlimit)
        .def_readwrite("minor_iter_limit", &snoptConfig::minoriterlimit)
        .def_readwrite("iter_limit", &snoptConfig::iterationslimit)
        .def_readwrite("opt_tol", &snoptConfig::optTol)
        .def_readwrite("fea_tol", &snoptConfig::feaTol)
        .def_readwrite("int_options", &snoptConfig::intOptions)
        .def_readwrite("float_options", &snoptConfig::floatOptions)
        .def_readwrite("string_options", &snoptConfig::stringOptions)
        .def_readwrite("printFile", &snoptConfig::printFile)
        .def_readwrite("printLevel", &snoptConfig::printlevel)
        .def_readwrite("minorPrintLevel", &snoptConfig::minorprintlevel)
        .def_readwrite("verifyLevel", &snoptConfig::verifylevel)
        .def_readwrite("majorIterLimit", &snoptConfig::majoriterlimit)
        .def_readwrite("minorIterLimit", &snoptConfig::minoriterlimit)
        .def_readwrite("iterLimit", &snoptConfig::iterationslimit)
        .def_readwrite("optTol", &snoptConfig::optTol)
        .def_readwrite("feaTol", &snoptConfig::feaTol)
        .def_readwrite("intOptions", &snoptConfig::intOptions)
        .def_readwrite("floatOptions", &snoptConfig::floatOptions)
        .def_readwrite("stringOptions", &snoptConfig::stringOptions);


        //.def("init", (void (pyServer::*)(const config &cfgin)) &pyServer::init, "init by config object")
    py::class_<pySnoptWrapper>(m, "SnoptSolver", R"pbdoc(
        The snopt solver class that is constructed by a SnoptProblem object and a SnoptConfig object.

        Args:
            prob (SnoptProblem): the problem object.
            cfg (SnoptConfig): the configuration.
        )pbdoc")
        .def(py::init<ProblemFun&, snoptConfig&>())
        .def("solve_rand", (optResult (pySnoptWrapper::*)()) &pySnoptWrapper::solve, R"pbdoc(
            Solve the problem using a random guess generated by the solver itself.

            Returns:
                SnoptResult: a SnoptResult object.
        )pbdoc")
        .def("solve_guess", (optResult (pySnoptWrapper::*)(cRefV x)) &pySnoptWrapper::solve, R"pbdoc(
            Solve the problem using a user-specified guess.

            Args:
                x (ndarray): an initial guess.

            Returns:
                SnoptResult: a SnoptResult object.
        )pbdoc")
        .def("solve_more", &pySnoptWrapper::solve_more, R"pbdoc(
            Use hot restart and continue a few major iterations.
            It can be used when SNOPT is solely used for backtracking.

            Args:
                iter (int): desired iteration number

            Returns:
                SnoptResult: a SnoptResult object
        )pbdoc")
        .def("obj_search", &pySnoptWrapper::obj_search, R"pbdoc(
            Given an initial, use search style to control how solver proceed.

            Args:
                guess (ndarray): an initial guess
                iter (int): initial iteration time
                step (int): iteration time for later search
                step_num (int): how many more step iterations do we play.
                abstol (float): absolute cost tolerance, exit if improvement is smaller than it
                reltol (float): relative cost tolerance, exit if improvement is smaller than reltol * last cost

            Returns:
                SnoptResult: a SnoptResult object
                costs: ndarray, recorded costs during search
                times: ndarray, recorded used times during search
        )pbdoc")
        .def("get_info", &pySnoptWrapper::getInfo, R"pbdoc(
            Get the inform variable from the solver. It indicates solving status.

            Returns:
                int: the info.
        )pbdoc")
        .def("set_workspace", &pySnoptWrapper::setWorkspace, R"pbdoc(
            Set the integer and real workspace size, this is necessary sometimes.

            Args:
                size (int): the desired integer workspace size.
                size (int): the desired real workspace size.
        )pbdoc")
        .def("fEval", &pySnoptWrapper::fEval)
        .def("setOptTol", &pySnoptWrapper::setOptTol)
        .def("setFeaTol", &pySnoptWrapper::setFeaTol)
        .def("setIntOption", &pySnoptWrapper::setIntOption)
        .def("setFloatOption", &pySnoptWrapper::setDoubleOption)
        .def("setWorkspace", &pySnoptWrapper::setWorkspace)
        .def("setMajorIter", &pySnoptWrapper::setMajorIter)
        .def("setPrintFile", &pySnoptWrapper::setPrintFile)
        .def("solveRand", (optResult (pySnoptWrapper::*)()) &pySnoptWrapper::solve)
        .def("solveGuess", (optResult (pySnoptWrapper::*)(cRefV x)) &pySnoptWrapper::solve)
        .def("getInfo", &pySnoptWrapper::getInfo)
        .def("fEval", &pySnoptWrapper::fEval);

    m.def("directSolve", &directSolve);
    m.def("inDirectSolve", &directInSolve);
    m.def("gradSolve", &gradFunSolve);
    m.def("inGradSolve", &inGradFunSolve);
    m.def("spGradSolve", &spGradFunSolve);
    m.def("inSpGradSolve", &inSpGradFunSolve);
    m.attr("__with_snopt__") = true;
#else
    m.attr("__with_snopt__") = false;
#endif

#ifdef IPOPT

    // this class was called snoptConfig previously, consider changing __init__.py
    py::class_<IpoptConfig>(m, "IpoptConfig", R"pbdoc(
        A class that specifies configuration to ipopt solver.

        This class has the default constructor.

        Attributes:
            print_level (int): level for printing = 4
            max_iter (int): maximum iteration = 1000
            linear_solver (string): linear solver = "mumps";
            hessian_approximation (string): Hessian approach = "limited-memory"
            jacobian_appximation (string): Jacobian approach = "exact";
            derivative_test (string): if test derivative = "no";
            tol (float): convergence tolerance = 1e-6;
            constr_vio_tol (float): tolerance on constraint = 1e-6;
            max_cpu_time (float): computation time = 1e6;
    )pbdoc")
        .def(py::init<>())
        .def("add_int_option", &IpoptConfig::addIntOption, R"pbdoc(
            Add an integer option to configurations.

            Args:
                option (str): the option to configure
                value (int): the integer value for configuration
        )pbdoc")
        .def("add_float_option", &IpoptConfig::addFloatOption, R"pbdoc(
            Add a float option to configurations.

            Args:
                option (str): the option to configure
                value (float): the float value for configuration
        )pbdoc")
        .def("add_string_option", &IpoptConfig::addPairStringOption, R"pbdoc(
            Add a string option to configurations.

            Args:
                option (str): the option to configure
                value (str): the string value for configuration
        )pbdoc")
        .def("set_major_iter", &IpoptConfig::setMajorIter, R"pbdoc(
            Set up the major iteration limit.

            Args:
                iter (int): limit for major iteration 
        )pbdoc")
        .def("set_opt_tol", &IpoptConfig::setOptTol, R"pbdoc(
            Set up the optimization tolerance.

            Args:
                tol (float): the tolerance
        )pbdoc")
        .def("set_fea_tol", &IpoptConfig::setFeaTol, R"pbdoc(
            Set up the feasibility tolerance.

            Args:
                tol (float): the tolerance
        )pbdoc")
        .def("set_print_level", &IpoptConfig::setPrintLevel, R"pbdoc(
            Set up the print level.

            Args:
                level (int): the print level
        )pbdoc")
        .def("enable_deriv_check", &IpoptConfig::enableDerivCheck, py::arg("lvl")=0, R"pbdoc(
            Set up the derivative check param.

            Args:
                level (int): the derivative check level
        )pbdoc")
        .def("set_print_freq", &IpoptConfig::setPrintFreq, R"pbdoc(
            Set up the print frequency.

            Args:
                freq (int): the print frequency 
        )pbdoc")
        .def("set_linear_solver", &IpoptConfig::setLinearSolver, R"pbdoc(
            Set the linear solver.

            Args:
                solver (str): the selected linear solver.
        )pbdoc")
        .def("enable_exact_hessian", &IpoptConfig::enableExactHessian, R"pbdoc(
            Enable usage of exact Hessian.
        )pbdoc")
        .def("enable_fd_jacobian", &IpoptConfig::enableFDJacobian, R"pbdoc(
            Enable usage of finite-difference approximation of Jacobian.
        )pbdoc")

        .def_readwrite("print_level", &IpoptConfig::print_level)
        .def_readwrite("print_frequency_iter", &IpoptConfig::print_frequency_iter)
        .def_readwrite("max_iter", &IpoptConfig::max_iter)
        .def_readwrite("linear_solver", &IpoptConfig::linear_solver)
        .def_readwrite("hessian_approximation", &IpoptConfig::hessian_approximation)
        .def_readwrite("jacobian_approximation", &IpoptConfig::jacobian_approximation)
        .def_readwrite("derivative_test", &IpoptConfig::derivative_test)
        .def_readwrite("tol", &IpoptConfig::tol)
        .def_readwrite("constr_vio_tol", &IpoptConfig::constr_viol_tol)
        .def_readwrite("max_cpu_time", &IpoptConfig::max_cpu_time)
         ;

    py::class_<IpoptSolver>(m, "IpoptSolver", R"pbdoc(
        The solver class.
        )pbdoc")
        .def(py::init<ProblemFun&, IpoptConfig &>())
        .def("solve_rand", &IpoptSolver::solve_rand)
        .def("solve_guess", &IpoptSolver::solve_guess)
        ;

    m.def("solve_problem", &solve_problem);

    m.def("set_verbosity", [](bool verbose) {VERBOSE = verbose;});
    m.attr("__with_ipopt__") = true;
#else
    m.attr("__with_ipopt__") = false;
#endif

    m.def("set_fd_step", [](double step) {PYOPTSOLVER_FD_STEP = step;});

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
