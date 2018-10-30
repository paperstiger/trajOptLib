/*
 * solverwrapper.cpp
 * Copyright (C) 2017 Gao <gao.tang@duke.edu>
 *
 * Distributed under terms of the  license.
 */

/* Write a python wrapper for the problem */

#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "snoptWrapper.h"
#include "funcStyle.h"
#include "classStyle.h"


#define VERSION_INFO "0.0.1"


class result : public optResult{
public:
    using optResult::optResult;
};

namespace py = pybind11;
PYBIND11_MODULE(libsnopt, m){
    m.doc() = R"pbdoc(
        A Python interface for SNOPT
    )pbdoc";

    // this class was called result previously, consider changing __init__.py
    py::class_<optResult>(m, "SnoptResult", R"pbdoc(
        A class that contains the solution returned by SNOPT.

        This class has the default constructor.

        Attributes:
            flag (int): the return flag by SNOPT.solve()
            obj (float): the cost function
            sol (ndarray): a copy of the solution
            sol (fval): a copy of F
            lmd (fval): a copy of Lagrangian multipliers
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
        .def_readonly("flag", &optResult::flag)
        .def_readonly("obj", &optResult::val)
        .def_readonly("sol", &optResult::sol)
        .def_readonly("fval", &optResult::c)
        .def_readonly("lmd", &optResult::lmd);

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
                value (str): the string value for configuration
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

    py::class_<funBase, pyFunBase>(m, "FunBase", R"pbdoc(
        A class sub-classed by users to define a custom optimization problem

        .. automethod:: __callf__.
        .. automethod:: __callg__.

        This function evaluates :math:`y=f(x)` where :math:`x` is optimization variables 
        and :math:`y` is the constraint where the first entry is cost function.

        Attributes:
            nx (int): number of variables to be optimized.
            nf (int): number of constraints, including the cost function itself at first entry.
            nG (int=0): nnz of the Jacobian matrix.
            grad (bool=False): if the function __callg__ has been implemented so we enable gradient.

        Note:
            This class has 2 constructors

        The first constructor accepts two integers and assume we do not know gradient information yet.

        Args:
            nx (int): the dimension of decision variables, i.e. nx.
            nf (int): the dimension of constraints, i.e. nf.

        The second constructor accepts an additional integers indicating nnz of Jacobian. This also sets
        grad to be true.

        Args:
            nx (int): the dimension of decision variables, i.e. nx.
            nf (int): the dimension of constraints, i.e. nf.
            nG (int): the nnz of Jacobian, i.e. nG.
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<int, int, int>())
        .def("__callf__", (int (funBase::*)(cRefV, RefV)) &pyFunBase::operator(), R"pbdoc(
            The function to be implemented by the user that evaluate constraints without Jacobian.

            The first argument is the candidate solution, it is 1d non-writable float ndarray of size (nx,)
            The second argument is to be written, it is 1d writable float ndarray of size (nf,). The first entry
            should be the cost function and others are constraints.

            Args:
                x (ndarray): a 1d ndarray of the candidate solution.
                f (ndarray): a 1d writable ndarray recording constraint function values.

            Returns:
                int: 0 or the length of f. This is necessary only if you want to let solver automatically detect length of f.
        )pbdoc")
        .def("__callg__", (std::pair<int, int> (funBase::*)(cRefV, RefV, RefV, RefVi, RefVi, bool, bool)) &pyFunBase::operator(), R"pbdoc(
            The function to be implemented by the user that evaluate constraints and Jacobian.

            This function has to be implemented if analytic gradient is used. This function is designed such that it does not 
            allocate any unnecessary memory so written to variables are done in place.
            The Jacobian is calculated using the triplet convention and the value, row, column are passed in as 1d ndarray with
            datatype of float, int, and int. The user can directly write to them.
            Furthermore, two flags controls written to the triplet. If needg is False, do not alter G.
            If rec is False, do not alter row and col. All indexes are 0 based.
            This function returns a tuple of two integers to indicate size of constraint and Jacobian.
            This is necessary is the user wants to automatically detect sizes of them, otherwise (0, 0) is acceptable.

            Args:
                x (ndarray): a 1d ndarray of the candidate solution.
                y (ndarray): a 1d writable ndarray recording constraint function values.
                G (ndarray): a 1d ndarray recording Jacobian matrix values.
                row (ndarray): an integer 1d ndarray recording rows of nnz of Jacobian.
                col (ndarray): an integer 1d ndarray recording columns of nnz of Jacobian.
                rec (bool): if True, row and col have to be modified; otherwise do not alter them.
                needg (bool): if True, modify G to record gradients; otherwise do not alter it.
            Retruns:
                tuple: a tuple of two integers, recording lengths of y and G. Can be (0, 0)
        )pbdoc")
        .def_readwrite("nx", &funBase::nx)
        .def_readwrite("nf", &funBase::nf)
        .def_readwrite("nG", &funBase::nG)
        .def_readwrite("grad", &funBase::grad);

    py::class_<ProblemFun, pyProbFun>(m, "SnoptProblem", R"pbdoc(
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
        .def("detect_prob_size", &pyProbFun::detect_prob_size, 
                R"pbdoc(
            Automatically detect nf and nG.

            This function requires __callf__ or __callg__ return non-trivial values.

            Args:
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
        .def("batchSetLb", &pyProbFun::batchSetLb)
        .def("batchSetUb", &pyProbFun::batchSetUb)
        .def("batchSetXlb", &pyProbFun::batchSetXlb)
        .def("batchSetXub", &pyProbFun::batchSetXub)
        .def("randomGenX", &pyProbFun::randomGenX)
        .def("detectNg", &pyProbFun::detectNg)
        .def("__callf__", (int (ProblemFun::*)(cRefV, RefV)) &pyProbFun::operator())
        .def("__callg__", (std::pair<int, int> (ProblemFun::*)(cRefV, RefV, RefV, RefVi, RefVi, bool, bool)) &pyProbFun::operator())
        .def_readwrite("nx", &pyProbFun::nx)
        .def_readwrite("nf", &pyProbFun::nf)
        .def_readwrite("nG", &pyProbFun::nG)
        .def_readwrite("grad", &pyProbFun::grad)
        .def_readwrite("Aval", &pyProbFun::Aval)
        .def_readwrite("Arow", &pyProbFun::Arow)
        .def_readwrite("Acol", &pyProbFun::Acol)
        .def_readwrite("lb", &pyProbFun::lb)
        .def_readwrite("ub", &pyProbFun::ub)
        .def_readwrite("xlb", &pyProbFun::xlb)
        .def_readwrite("xub", &pyProbFun::xub);

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
        .def("solve_guess", (optResult (pySnoptWrapper::*)(RefV x)) &pySnoptWrapper::solve, R"pbdoc(
            Solve the problem using a user-specified guess.

            Args:
                x (ndarray): an initial guess.

            Returns:
                SnoptResult: a SnoptResult object.
        )pbdoc")
        .def("get_info", &pySnoptWrapper::getInfo, R"pbdoc(
            Get the inform variable from the solver. It indicates solving status.

            Returns:
                int: the info.
        )pbdoc")
        .def("set_int_workspace", &pySnoptWrapper::setIntWorkspace, R"pbdoc(
            Set the integer workspace size, this is necessary sometimes.

            Args:
                size (int): the desired integer workspace size.
        )pbdoc")
        .def("set_real_workspace", &pySnoptWrapper::setRealWorkspace, R"pbdoc(
            Set the real workspace size, this is necessary sometimes.

            Args:
                size (int): the desired real workspace size.
        )pbdoc")
        .def("fEval", &pySnoptWrapper::fEval)
        .def("setOptTol", &pySnoptWrapper::setOptTol)
        .def("setFeaTol", &pySnoptWrapper::setFeaTol)
        .def("setIntOption", &pySnoptWrapper::setIntOption)
        .def("setFloatOption", &pySnoptWrapper::setDoubleOption)
        .def("setIntWorkspace", &pySnoptWrapper::setIntWorkspace)
        .def("setRealWorkspace", &pySnoptWrapper::setRealWorkspace)
        .def("setMajorIter", &pySnoptWrapper::setMajorIter)
        .def("setPrintFile", &pySnoptWrapper::setPrintFile)
        .def("solveRand", (optResult (pySnoptWrapper::*)()) &pySnoptWrapper::solve)
        .def("solveGuess", (optResult (pySnoptWrapper::*)(RefV x)) &pySnoptWrapper::solve)
        .def("getInfo", &pySnoptWrapper::getInfo)
        .def("fEval", &pySnoptWrapper::fEval);

    m.def("directSolve", &directSolve);
    m.def("inDirectSolve", &directInSolve);
    m.def("gradSolve", &gradFunSolve);
    m.def("inGradSolve", &inGradFunSolve);
    m.def("spGradSolve", &spGradFunSolve);
    m.def("inSpGradSolve", &inSpGradFunSolve);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
