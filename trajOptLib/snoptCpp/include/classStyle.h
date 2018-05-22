/*
 * classStyle.h
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef CLASSSTYLE_H
#define CLASSSTYLE_H


#include "snoptWrapper.h"
#include "functionBase.h"
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

        void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_PURE_NAME(
                    void,
                    ProblemFun,
                    "__callg__",
                    operator(),
                    x, F, G, row, col, rec, needg
                    );
        }
};


class pySnoptWrapper: public snoptWrapper{
    public:
        pySnoptWrapper(ProblemFun &fun, snoptConfig &cfg): snoptWrapper(&fun, &cfg){
#ifdef DEBUG
            std::cout << "Entering construct of pySnoptWrapper\n";
#endif
        };

        optResult solve(RefV x){
            int flag = snoptWrapper::solve(x.data());
            optResult rst;
            rst.flag = flag;
            copyToResult(rst);
            return rst;
        }

        optResult solve(){
            srand(time(NULL));
            int flag = snoptWrapper::solve();
            optResult rst;
            rst.flag = flag;
            copyToResult(rst);
            return rst;
        }

        void copyToResult(optResult &rst){
            rst.c.resize(neF);
            rst.lmd.resize(neF);
            rst.sol.resize(n);
            // evaluate and final call
            snoptWrapper::copyX(rst.sol.data());
            snoptWrapper::copyF(rst.c.data());
            snoptWrapper::copyLmd(rst.lmd.data());
            rst.val = snoptWrapper::getObj();
        }

        VX fEval(RefV x){
            return snoptWrapper::fEval(x.data());
        }
};


class pyFunBase: public funBase{
    public:
        using funBase::funBase;

        void operator()(cRefV x, RefV F) override{
            PYBIND11_OVERLOAD_PURE_NAME(
                    void,
                    funBase,
                    "__callf__",
                    operator(),
                    x, F
                    );
        }

        void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_PURE_NAME(
                    void,
                    funBase,
                    "__callg__",
                    operator(),
                    x, F, G, row, col, rec, needg
                    );
        }
};

#endif /* !CLASSSTYLE_H */
