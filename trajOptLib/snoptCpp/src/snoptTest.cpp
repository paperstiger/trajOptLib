/*
 * snoptTest.cpp
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#include "snoptWrapper.h"

class userFun: public ProblemFun{
public:
    userFun(): ProblemFun(2, 2){
        setNg(4);
        setGrad(false);
        xlb = VX::Zero(2);
        xub = VX::Zero(2);
        xlb(0) = 0.6;
        xub(0) = 1.0;
        xlb(1) = 0.2;
        xub(1) = 1.0;
        lb = VX::Zero(2);
        ub = VX::Zero(2);
        lb(1) = 1.0;
        ub(1) = 1.0;
    }
    void operator()(cRefV x, RefV F){
        F(0) = x(0) * x(0) + x(1) * x(1);
        F(1) = x(0) + x(1);
    }
    void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg){
        F(0) = x(0) * x(0) + x(1) * x(1);
        F(1) = x(0) + x(1) - 1.0;
        if(rec){
            row(0) = 0;
            row(1) = 0;
            row(2) = 1;
            row(3) = 1;
            col(0) = 0;
            col(1) = 1;
            col(2) = 0;
            col(3) = 1;
        }
        if(needg){
            G(0) = 2 * x(0);
            G(1) = 2 * x(1);
            G(2) = 1;
            G(3) = 1;
        }
    }
};

int main(){
    srand(time(NULL));
    userFun fun;
    int nx = fun.getNx(), nf = fun.getNf();
    snoptWrapper snwrap(&fun);
    int flag = snwrap.solve();
    printf("flag is %d\n", flag);
    MapV mX(snwrap.getX(), snwrap.getVarNum());
    std::cout << mX << std::endl;
    return 0;
}
