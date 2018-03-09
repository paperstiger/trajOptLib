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


#define NO_DEBUG


// define several general classes that might be passed in
// First three types support direct calling
// Last three need user specified problem size so I create temporary arrays for them
typedef std::function<VX(cRefV)> plainFun;  //function of type y=f(x)
typedef std::function<std::tuple<VX, rMX>(cRefV)> gradFun;  //function of type y, J = f(x)
typedef std::function<std::tuple<VX, SpMX>(cRefV)> spGradFun;  //function of type y, sJ = f(x)
typedef std::function<void(cRefV, RefV)> inPlainFun;  //function of type f(x, y). User have to specify size
typedef std::function<void(cRefV, RefV, rRefM)> inGradFun;  //function of type f(x, y, J). User have to specify size
typedef std::function<void(cRefV, RefV, RefSpMX)> inSpGradFun;  //function of type f(x, y, J). User have to specify size


enum FunType {PLAIN, GRAD, SPGRAD, INPLAIN, INGRAD, INSPGRAD};
class funcWrapper : public ProblemFun {
public:
    funcWrapper(const plainFun &f, int nx_, int nf_){
        // Initialize base class
        nx = nx_;
        nf = nf_;
        nG = 0;
        grad = false;
        // Initialize itself
        plainfun = f;
        mode = PLAIN;
    }

    void operator()(cRefV x, RefV F){
        switch(mode){
            case PLAIN:
                VX y = plainfun(x);
                F = y;
                break;
        }
    }
    void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, int rowadd, int nGadd, bool rec, bool needg){}
private:
    plainFun plainfun;
    gradFun gradfun;
    spGradFun spgradfun;
    inPlainFun inplainfun;
    inGradFun ingradfun;
    inSpGradFun inspgradfun;
    FunType mode;
};


class optResult{
    public:
        int flag;
        double val;
        VX sol;
        VX c;
        optResult(){}
};
// plain function to use default setting for solving a simple problem
// param: fun The function being called.
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult directSolve(plainFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap);
    int flag = snwrap.solve(x0.data());
    // copy the results back
    optResult rst;
    rst.flag = flag;
    // evaluate and final call
    VX c = snwrap.fEval(snwrap.getX());
    rst.c = c;
    MapV Mx(snwrap.getX(), snwrap.getVarNum());
    rst.sol = Mx;
    return rst;
}


namespace py = pybind11;
PYBIND11_MODULE(libsnopt, m){
    py::class_<optResult>(m, "result")
        .def(py::init<>())
        .def_readwrite("flag", &optResult::flag)
        .def_readwrite("obj", &optResult::val)
        .def_readwrite("sol", &optResult::sol)
        .def_readwrite("fval", &optResult::c);

    m.def("directSolve", &directSolve);
}
