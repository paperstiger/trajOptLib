/*
 * pybindWrapper.cpp
 * Copyright (C) 2019 motion <motion@motion-MS-7971>
 *
 * Distributed under terms of the MIT license.
 */

#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "ipoptWrapper.h"


#define VERSION_INFO "0.0.1"

bool VERBOSE = false;


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

        pint operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_NAME(
                    pint,
                    ProblemFun,
                    "__callg__",
                    operator(),
                    x, F, G, row, col, rec, needg
                    );
        }

        double evalF(cRefV x) override{
            PYBIND11_OVERLOAD_NAME(
                    double,
                    ProblemFun,
                    "eval_cost",
                    evalF,
                    x
                    );
        }

        bool evalGrad(cRefV x, RefV grad) override{
            PYBIND11_OVERLOAD_NAME(
                    bool,
                    ProblemFun,
                    "eval_gradient",
                    evalGrad,
                    x, grad
                    );
        }

        int evalG(cRefV x, RefV g) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "eval_constr",
                    evalG,
                    x, g
                    );
        }

        int evalJac(cRefV x, RefV g, RefVi row, RefVi col, bool rec) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "eval_jacobian",
                    evalJac,
                    x, g, row, col, rec
                    );
        }

        int evalHess(cRefV x, double sigma, cRefV lmd, RefV g, RefVi row, RefVi col, bool rec) override{
            PYBIND11_OVERLOAD_NAME(
                    int,
                    ProblemFun,
                    "eval_hessian",
                    evalHess,
                    x, sigma, lmd, g, row, col, rec
                    );
        }
};


namespace py = pybind11;
using namespace pybind11::literals;
PYBIND11_MODULE(libpyipopt, m){
    m.doc() = R"pbdoc(
        A Python interface for IPOPT
    )pbdoc";

    // this class was called result previously, consider changing __init__.py
    py::class_<optResult>(m, "IpoptResult", R"pbdoc(
        A class that contains the solution returned by SNOPT.

        This class has the default constructor.

        Attributes:
            flag (int): the return flag by SNOPT.solve()
            obj (float): the cost function
            sol (ndarray): a copy of the solution
            sol (ndarray): a copy of F
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
        .def("add_string_option", &IpoptConfig::addStringOption, R"pbdoc(
            Add a string option to configurations.

            Args:
                option (str): the option to configure
                value (str): the string value for configuration
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

    py::class_<ProblemFun, pyProbFun>(m, "IpoptProblem", R"pbdoc(
        A class sub-classed by users to define a custom optimization problem.

        :currentmodule:
        .. automethod:: __callf__.
        .. automethod:: __callg__.

        This class defines an optimization problem to optimize :math:`y(0)` where :math:`y=f(x)+Ax`
        subject to constraints :math:`lb \le y \le ub` and :math:`xlb \le x \le xub`.

        Attributes:
            nx (int): number of variables to be optimized. The same with :class:`.FunBase`.
            nf (int): number of constraints. The same with :class:`.FunBase`.
            nG (int=0): nnz of the Jacobian matrix. The same with :class:`.FunBase`.
            grad (bool=False): if the function __callg__ has been implemented so we enable gradient.
            The same with :class:`.FunBase`.
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
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<int, int, int>())
        .def("update_nf", &ProblemFun::updateNf, R"pbdoc(
            Update nf and reallocate space for lb and ub, if necessary.

            Args:
                nf (int): the desired nf
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
        .def("__callg__", (std::pair<int, int> (ProblemFun::*)(cRefV, RefV, RefV, RefVi, RefVi, bool, bool)) &pyProbFun::operator())
        .def("eval_cost", &ProblemFun::evalF)
        .def("eval_constr", &ProblemFun::evalG)
        .def("eval_gradient", &ProblemFun::evalGrad)
        .def("eval_jacobian", &ProblemFun::evalJac)
        .def("eval_hessian", &ProblemFun::evalHess)
        .def("ipstyle", [](ProblemFun *self) {self->ipStyle=true;})
        .def_readwrite("nx", &pyProbFun::nx)
        .def_readwrite("nf", &pyProbFun::nf)
        .def_readwrite("nG", &pyProbFun::nG)
        .def_readwrite("nH", &pyProbFun::hess_nnz)
        .def_readwrite("grad", &pyProbFun::grad)
        .def_readonly("Aval", &pyProbFun::Aval)
        .def_readonly("Arow", &pyProbFun::Arow)
        .def_readonly("Acol", &pyProbFun::Acol)
        .def_readonly("lb", &pyProbFun::lb)
        .def_readonly("ub", &pyProbFun::ub)
        .def_readonly("xlb", &pyProbFun::xlb)
        .def_readonly("xub", &pyProbFun::xub);

    py::class_<IpoptSolver>(m, "IpoptSolver", R"pbdoc(
        The solver class.
        )pbdoc")
        .def(py::init<ProblemFun&, IpoptConfig &>())
        .def("solve_rand", &IpoptSolver::solve_rand)
        .def("solve_guess", &IpoptSolver::solve_guess)
        ;

    m.def("solve_problem", &solve_problem);

    m.def("set_verbose", [](bool verbose) {VERBOSE = verbose;});

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

