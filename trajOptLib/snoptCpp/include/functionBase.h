/*
 * functionBase.h
 * Copyright (C) 2018 Gao Tang <gt70@duke.edu>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef FUNCTIONBASE_H
#define FUNCTIONBASE_H

#include "TigerTools/TigerEigen.h"

typedef std::pair<int, int> pint;

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

        virtual int operator()(cRefV x, RefV F) = 0;  // A function to be overwritten by subclass, this is called to evaluate

        virtual std::pair<int, int> operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) = 0;  // A function to be overwritten by subclass, this is called for both assigning structure.

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
