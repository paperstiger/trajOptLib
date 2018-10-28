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
        A class sub-classed by users to define a custom optimization problem.

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
        .def("__callf__", (void (funBase::*)(cRefV, RefV)) &pyFunBase::operator(), R"pbdoc(
            The function to be implemented by the user that evaluate constraints without Jacobian.

            The first argument is the candidate solution, it is 1d non-writable float ndarray of size (nx,)
            The second argument is to be written, it is 1d writable float ndarray of size (nf,). The first entry
            should be the cost function and others are constraints.

            Args:
                x (ndarray): a 1d ndarray of the candidate solution.
                f (ndarray): a 1d writable ndarray recording function evaluation.

            Returns:
                This function has no returns.
        )pbdoc")
        .def("__callg__", (void (funBase::*)(cRefV, RefV, RefV, RefVi, RefVi, bool, bool)) &pyFunBase::operator())
        .def_readwrite("nx", &funBase::nx)
        .def_readwrite("nf", &funBase::nf)
        .def_readwrite("nG", &funBase::nG)
        .def_readwrite("grad", &funBase::grad);

    py::class_<ProblemFun, pyProbFun>(m, "probFun")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<int, int, int>())
        .def("batchSetLb", &pyProbFun::batchSetLb)
        .def("batchSetUb", &pyProbFun::batchSetUb)
        .def("batchSetXlb", &pyProbFun::batchSetXlb)
        .def("batchSetXub", &pyProbFun::batchSetXub)
        .def("randomGenX", &pyProbFun::randomGenX)
        .def("detectNg", &pyProbFun::detectNg)
        .def("__callf__", (void (ProblemFun::*)(cRefV, RefV)) &pyProbFun::operator())
        .def("__callg__", (void (ProblemFun::*)(cRefV, RefV, RefV, RefVi, RefVi, bool, bool)) &pyProbFun::operator())
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
    py::class_<pySnoptWrapper>(m, "solver")
        .def(py::init<ProblemFun&, snoptConfig&>())
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
