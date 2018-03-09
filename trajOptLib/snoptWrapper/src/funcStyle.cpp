/*
 * funcStyle.cpp
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#include "funcStyle.h"

// plain function to use default setting for solving a simple problem
// param: fun The function being called.
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult directSolve(plainFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap, cfg);
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

// inPlain function to use default setting for solving a simple problem
// param: fun The function being called.
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult directInSolve(inPlainFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap, cfg);
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

// grad function to use default setting for solving a simple problem
// param: fun The function being called in form of f, J = fun(x)
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult gradFunSolve(gradFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap, cfg);
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


// grad function to use default setting for solving a simple problem
// param: fun The function being called in form of fun(x, f, J)
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult inGradFunSolve(inGradFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap, cfg);
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


// grad function to use default setting for solving a simple problem
// param: fun The function being called in form of f, J, spM = fun(x)
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult spGradFunSolve(spGradFun fun, RefV x0, int nx, int nf, int nG, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf, nG);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap, cfg);
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


// sp grad function to use default setting for solving a simple problem
// param: fun The function being called in form of fun(x, f, G, row, col, rec)
// param: x0 The initial guess provided
// param: xlb The lower bound on variable x
// param: xub The upper bound on variable x
// param: lb The lower bound on function f
// param: ub The upper bound on function f
optResult inSpGradFunSolve(inSpGradFun fun, RefV x0, int nx, int nf, int nG, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg){
    // construct the wrapper
    funcWrapper funWrap(fun, nx, nf, nG);
    // set the bounds
    funWrap.lb = lb;
    funWrap.ub = ub;
    funWrap.xlb = xlb;
    funWrap.xub = xub;
    // construct snopt wrapper
    snoptWrapper snwrap(&funWrap, cfg);
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
