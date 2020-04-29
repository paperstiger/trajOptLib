/*
 * functionBase.h
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef FUNCTIONBASE_H
#define FUNCTIONBASE_H

#include "TigerEigen.h"

typedef std::pair<int, int> pint;
extern double PYOPTSOLVER_FD_STEP;

class funBase{
    public:
        int nx, nf;
        bool grad = false;
        int nG = 0;

        funBase(){
#ifdef DEBUG
            std::cout << "Entering construct of funBase\n";
#endif
        }
        funBase(int nx_, int nf_): nx(nx_), nf(nf_){
#ifdef DEBUG
            std::cout << "Entering construct of funBase\n";
#endif
        }
        funBase(int nx_, int nf_, int ng_): nx(nx_), nf(nf_), nG(ng_), grad(true){
#ifdef DEBUG
            std::cout << "Entering construct of funBase\n";
#endif
        }

#ifdef ENABLEIP
        virtual double evalF(cRefV x) = 0;
        virtual bool evalGrad(cRefV x, RefV grad) = 0;
        virtual int evalG(cRefV x, RefV g) = 0;
        virtual int evalJac(cRefV x, RefV G, RefVl row, RefVl col, bool rec) = 0;
#endif

        virtual int operator()(cRefV x, RefV F) = 0;  // A function to be overwritten by subclass, this is called to evaluate

        virtual std::pair<int, int> operator()(cRefV x, RefV F, RefV G, RefVl row, RefVl col, bool rec, bool needg)   // A function to be overwritten by subclass, this is called for both assigning structure.
        {
            int flag = operator()(x, F);
            if(needg) {
                VX tmpf(nf), tmpx(nx);
                for(int i = 0; i < nx; i++) {
                    tmpx = x;
                    tmpx(i) += PYOPTSOLVER_FD_STEP;
                    operator()(tmpx, tmpf);
                    // write to G, column by column
                    G.segment(i * nf, nf).array() = (tmpf - F).array() / PYOPTSOLVER_FD_STEP;
                    if(rec) {
                        row.segment(i * nf, nf).setLinSpaced(0, nf - 1);
                        col.segment(i * nf, nf).array() = i;
                    }
                }
            }
            return std::make_pair(flag, nx * nf);
        }

        int getNx() const {return nx;}

        int getNf() const {return nf;}

        bool getGrad() const {return grad;}

        int getNg() const {return nG;};

        void setNx(int nx_){
            nx = nx_;
        }
        void setNf(int nf_){
            nf = nf_;
        }
        void setGrad(bool grad_){
            grad = grad_;
        }
        void setNg(int ng_){
            nG = ng_;
        }
};

#endif /* !FUNCTIONBASE_H */
