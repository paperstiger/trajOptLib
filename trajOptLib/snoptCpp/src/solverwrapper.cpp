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


#define NO_DEBUG


namespace py = pybind11;
PYBIND11_MODULE(libsnopt, m){
    py::class_<optResult>(m, "result")
        .def(py::init<>())
        .def_readwrite("flag", &optResult::flag)
        .def_readwrite("obj", &optResult::val)
        .def_readwrite("sol", &optResult::sol)
        .def_readwrite("fval", &optResult::c);

    py::class_<snoptConfig>(m, "snoptConfig")
        .def(py::init<>())
        .def("addIntOption", &snoptConfig::addIntOption)
        .def("addFloatOption", &snoptConfig::addFloatOption)
        .def("addStringOption", &snoptConfig::addStringOption)
        .def_readwrite("name", &snoptConfig::name)
        .def_readwrite("printFile", &snoptConfig::printFile)
        .def_readwrite("printLevel", &snoptConfig::printlevel)
        .def_readwrite("verifyLevel", &snoptConfig::verifylevel)
        .def_readwrite("majorIterLimit", &snoptConfig::majoriterlimit)
        .def_readwrite("minorIterLimit", &snoptConfig::minoriterlimit)
        .def_readwrite("iterLimit", &snoptConfig::iterationslimit)
        .def_readwrite("optTol", &snoptConfig::optTol)
        .def_readwrite("feaTol", &snoptConfig::feaTol)
        .def_readwrite("intOptions", &snoptConfig::intOptions)
        .def_readwrite("floatOptions", &snoptConfig::floatOptions)
        .def_readwrite("stringOptions", &snoptConfig::stringOptions);

    py::class_<funBase, pyFunBase>(m, "funBase")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<int, int, int>())
        .def("__callf__", (void (funBase::*)(cRefV, RefV)) &pyFunBase::operator())
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
}
