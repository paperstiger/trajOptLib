/*
 * classStyle.h
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef CLASSSTYLE_H
#define CLASSSTYLE_H


#include "snoptWrapper.h"
#include "pybind11/functional.h"


class pyProbFun: public ProblemFun{
    public:
        using ProblemFun::ProblemFun;

        void operator()(cRefV x, RefV F) override{
            PYBIND11_OVERLOAD_PURE_NAME(
                    void,
                    ProblemFun,
                    "__callf__",
                    operator(),
                    x, F
                    );
        }

        void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, int rowadd, int nGadd, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_PURE_NAME(
                    void,
                    ProblemFun,
                    "__callg__",
                    operator(),
                    x, F, G, row, col, rowadd, nGadd, rec, needg
                    );
        }
};


class pySnoptWrapper: public snoptWrapper{
    public:
        pySnoptWrapper(ProblemFun &fun, snoptConfig &cfg): snoptWrapper(&fun, &cfg){};
        optResult solve(RefV x){
            int flag = snoptWrapper::solve(x.data());
            optResult rst;
            rst.flag = flag;
            // evaluate and final call
            VX c = snoptWrapper::fEval(getX());
            rst.c = c;
            MapV Mx(getX(), getVarNum());
            rst.sol = Mx;
            return rst;
        }
        optResult solve(){
            int flag = snoptWrapper::solve();
            optResult rst;
            rst.flag = flag;
            // evaluate and final call
            VX c = snoptWrapper::fEval(getX());
            rst.c = c;
            MapV Mx(getX(), getVarNum());
            rst.sol = Mx;
            return rst;
        }
        VX fEval(RefV x){
            return snoptWrapper::fEval(x.data());
        }
};

#endif /* !CLASSSTYLE_H */
