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
#include "pybind11/pybind11.h"


class pyProbFun: public ProblemFun{
    public:
        using ProblemFun::ProblemFun;

        int operator()(cRefV x, RefV F) override{
            PYBIND11_OVERLOAD_PURE_NAME(
                    int,
                    ProblemFun,
                    "__callf__",
                    operator(),
                    x, F
                    );
        }

        pint operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_PURE_NAME(
                    pint,
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
            rst.constr_vio.resize(neF - 1);
            rst.lmd.resize(neF);
            rst.sol.resize(n);
            // evaluate and final call
            snoptWrapper::copyX(rst.sol.data());
            snoptWrapper::copyF(rst.c.data());
            snoptWrapper::copyLmd(rst.lmd.data());
            rst.val = snoptWrapper::getObj();
            // write constraint violation
            for(int i = 1; i < neF; i++){
                double diff1 = rst.c(i) - prob->lb(i);
                double diff2 = rst.c(i) - prob->ub(i);
                if(diff1 < 0)
                    rst.constr_vio(i - 1) = diff1;
                else if(diff2 > 0)
                    rst.constr_vio(i - 1) = diff2;
                else
                    rst.constr_vio(i - 1) = 0;
            }
        }

        VX fEval(RefV x){
            return snoptWrapper::fEval(x.data());
        }
};


class pyFunBase: public funBase{
    public:
        using funBase::funBase;

        int operator()(cRefV x, RefV F) override{
            PYBIND11_OVERLOAD_PURE_NAME(
                    int,
                    funBase,
                    "__callf__",
                    operator(),
                    x, F
                    );
        }

        pint operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) override {
            PYBIND11_OVERLOAD_PURE_NAME(
                    pint,
                    funBase,
                    "__callg__",
                    operator(),
                    x, F, G, row, col, rec, needg
                    );
        }
};

#endif /* !CLASSSTYLE_H */
