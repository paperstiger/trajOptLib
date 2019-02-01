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
#include "pybind11/pybind11.h"


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

        optResult solve_more(int iter) {
            int flag = snoptWrapper::solve_more(iter);
            optResult rst;
            rst.flag = flag;
            copyToResult(rst);
            return rst;
        }

        std::tuple<optResult, VX, VX> obj_search(RefV x, int iter, int step, int step_num, double abstol, double reltol) {
            std::pair<VX, VX> rst = snoptWrapper::obj_search(x.data(), iter, step, step_num, abstol, reltol);
            optResult opt;
            opt.flag = getInfo();
            copyToResult(opt);
            return std::make_tuple(opt, rst.first, rst.second);
        }

        void copyToResult(optResult &rst){
            rst.c.resize(neF);
            rst.lmd.resize(neF);
            rst.sol.resize(n);
            rst.xmul.resize(n);
            // evaluate and final call
            snoptWrapper::copyX(rst.sol.data());
            snoptWrapper::copyF(rst.c.data());
            snoptWrapper::copyLmd(rst.lmd.data());
            snoptWrapper::copyXmul(rst.xmul.data());
            rst.val = snoptWrapper::getObj();
        }

        VX fEval(RefV x){
            return snoptWrapper::fEval(x.data());
        }
};


#endif /* !CLASSSTYLE_H */
