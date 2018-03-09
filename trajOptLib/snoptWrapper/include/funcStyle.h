/*
 * funcStyle.h
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef FUNCSTYLE_H
#define FUNCSTYLE_H

#include "snoptWrapper.h"

// define several general classes that might be passed in
// First three types support direct calling
// Last three need user specified problem size so I create temporary arrays for them
typedef std::function<VX(cRefV)> plainFun;  //function of type y=f(x)
typedef std::function<std::tuple<VX, rMX>(cRefV)> gradFun;  //function of type y, J = f(x)
typedef std::function<std::tuple<VX, SpMX>(cRefV)> spGradFun;  //function of type y, sJ = f(x)
typedef std::function<void(cRefV, RefV)> inPlainFun;  //function of type f(x, y). User have to specify size
typedef std::function<void(cRefV, RefV, rRefM)> inGradFun;  //function of type f(x, y, J). User have to specify size
typedef std::function<void(cRefV, RefV, RefV, RefVi, RefVi, bool)> inSpGradFun;  //function of type f(x, y, J). User have to specify size


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

    funcWrapper(const inPlainFun &f, int nx_, int nf_){
        // Initialize base class
        nx = nx_;
        nf = nf_;
        nG = 0;
        grad = false;
        // Initialize itself
        inplainfun = f;
        mode = INPLAIN;
    }

    funcWrapper(const gradFun &f, int nx_, int nf_){
        // Initialize base class
        nx = nx_;
        nf = nf_;
        nG = nx * nf;
        grad = true;
        // Initialize itself
        gradfun = f;
        mode = GRAD;
    }

    funcWrapper(const inGradFun &f, int nx_, int nf_){
        // Initialize base class
        nx = nx_;
        nf = nf_;
        nG = nx * nf;
        grad = true;
        // Initialize itself
        ingradfun = f;
        mode = INGRAD;
    }

    funcWrapper(const spGradFun &f, int nx_, int nf_, int nG_){
        // Initialize base class
        nx = nx_;
        nf = nf_;
        nG = nG_;
        grad = true;
        // Initialize itself
        spgradfun = f;
        mode = SPGRAD;
    }

    funcWrapper(const inSpGradFun &f, int nx_, int nf_, int nG_){
        // Initialize base class
        nx = nx_;
        nf = nf_;
        nG = nG_;
        grad = true;
        // Initialize itself
        inspgradfun = f;
        mode = INSPGRAD;
    }

    void operator()(cRefV x, RefV F){
        if(mode == PLAIN){
            F = plainfun(x);
        }
        else if(mode == INPLAIN){
            inplainfun(x, F);
        }
    }

    void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, int rowadd, int nGadd, bool rec, bool needg){
        if(mode == GRAD){
            if(rec){
                assignG(row, col, rowadd, nGadd);
            }
            rMapM rMG(G.data(), nf, nx);
            std::tuple<VX, rMX> rst = gradfun(x);
            F = std::get<0>(rst);
            rMG = std::get<1>(rst);
        }
        else if(mode == INGRAD){
            if(rec){
                assignG(row, col, rowadd, nGadd);
            }
            rMapM rMG(G.data(), nf, nx);
            ingradfun(x, F, rMG);
        }
        else if(mode == SPGRAD){
            std::tuple<VX, SpMX> rst = spgradfun(x);
            F = std::get<0>(rst);
            SpMX mat = std::get<1>(rst);
            int nG = nGadd;
            for (int k=0; k<mat.outerSize(); ++k){
                for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
                    G(nG) = it.value();
                    if(rec){
                        row(nG) = rowadd + it.row();
                        col(nG) = it.col();
                    }
                    nG++;
                }
            }
        }
        else if(mode == INSPGRAD){
            inspgradfun(x, F, G, row, col, rec);
        }
    }

    void assignG(RefVi row, RefVi col, int rowadd, int nGadd){
        int nG = nGadd;
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < nf; j++){
                row(nG) = rowadd + j;
                col(nG) = i;
                nG++;
            }
        }
    }

private:
    plainFun plainfun;
    gradFun gradfun;
    spGradFun spgradfun;
    inPlainFun inplainfun;
    inGradFun ingradfun;
    inSpGradFun inspgradfun;
    FunType mode;
};


optResult directSolve(plainFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg);

optResult directInSolve(inPlainFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg);

optResult gradFunSolve(gradFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg);

optResult inGradFunSolve(inGradFun fun, RefV x0, int nx, int nf, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg);

optResult spGradFunSolve(spGradFun fun, RefV x0, int nx, int nf, int nG, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg);

optResult inSpGradFunSolve(inSpGradFun fun, RefV x0, int nx, int nf, int nG, cRefV xlb, cRefV xub, cRefV lb, cRefV ub, snoptConfig *cfg);
#endif /* !FUNCSTYLE_H */
