/* Construct the SNOPT problem
 */
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <tuple>
#include <string>
#include "TigerEigen.h"
#include "functionBase.h"


class optResult{
    public:
        int flag;
        double val;
        VX sol;
        VX c;
        VX lmd;
        VX xmul;
        optResult(){}

        RefV get_sol(){
            return RefV(sol);
        }

        RefV get_fval(){
            return RefV(c);
        }

        RefV get_lambda(){
            return RefV(lmd);
        }

        RefV get_xmul() {
            return RefV(xmul);
        }
};


// a functor which defines problems
class ProblemFun : public funBase{
public:
    VX lb, ub;
    VX xlb, xub;
    VX Aval;
    VXi Arow, Acol;
    bool ipStyle = false;
    int hess_nnz;

    ProblemFun():funBase(){}

    ProblemFun(int nx_):funBase(){nx = nx_;}

    ProblemFun(int nx_, int nf_): funBase(nx_, nf_){
        _allocate_space();
    }
    ProblemFun(int nx_, int nf_, int ng_): funBase(nx_, nf_, ng_){
        _allocate_space();
    }

    virtual int operator()(cRefV x, RefV F){return 0;};  // A function to be overwritten by subclass, this is called to evaluate

    virtual std::pair<int, int> operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg){return std::make_pair(0, 0);};  // A function to be overwritten by subclass, this is called for both assigning structure.

    virtual double evalF(cRefV x) {return 0;};
    virtual bool evalGrad(cRefV x, RefV grad) {return true;};
    virtual int evalG(cRefV x, RefV g) {return 0;};
    virtual int evalJac(cRefV x, RefV G, RefVi row, RefVi col, bool rec) {return 0;};
    virtual int evalHess(cRefV x, double sigma, cRefV lmd, RefV G, RefVi row, RefVi col, bool rec) {return 0;}

    void _allocate_space() {
        if(lb.size() < nf){
            lb.resize(nf);
            ub.resize(nf);
        }
        if(xlb.size() < nx) {
            xlb.resize(nx);
            xub.resize(nx);
        }
    }

    void updateNf(int nf_) {
        nf = nf_;
        _allocate_space();
    }

    void updateNg(int ng_) {
        nG = ng_;
    }

    void setA(crRefM A){
        int nrow = A.rows(), ncol = A.cols();
        int nnz = 0;
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                if(A(i, j) != 0)
                    nnz++;
            }
        }
        Aval.resize(nnz);
        Arow.resize(nnz);
        Acol.resize(nnz);
        nnz = 0;
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                if(A(i, j) != 0){
                    Aval(nnz) = A(i, j);
                    Arow(nnz) = i;
                    Acol(nnz) = j;
                    nnz++;
                }
            }
        }
    }

    void setA(cRefV val, cRefVi row, cRefVi col){
        Aval = val;
        Arow = row;
        Acol = col;
    }

    void set_lb(cRefV lb_) {
        if(lb_.size() != nf) {
            printf("Incorrect size when setting lb\n");
            return;
        }
        lb = lb_;
    }

    void set_ub(cRefV ub_) {
        if(ub_.size() != nf) {
            printf("Incorrect size when setting ub\n");
            return;
        }
        ub = ub_;
    }

    void set_xlb(cRefV xlb_) {
        if(xlb_.size() != nx) {
            printf("Incorrect size when setting xlb\n");
            return;
        }
        xlb = xlb_;
    }

    void set_xub(cRefV xub_) {
        if(xub_.size() != nx) {
            printf("Incorrect size when setting xub\n");
            return;
        }
        xub = xub_;
    }

    void batchSetLb(cRefV lb_, int ind0=0){
        if(lb_.size() + ind0 > lb.size())
            lb.segment(ind0, lb.size() - ind0) = lb_.head(lb.size() - ind0);
        else
            lb.segment(ind0, lb_.size()) = lb_;
    }

    void batchSetUb(cRefV ub_, int ind0=0){
        if(ub_.size() + ind0 > ub.size())
            ub.segment(ind0, ub.size() - ind0) = ub_.head(ub.size() - ind0);
        else
            ub.segment(ind0, ub_.size()) = ub_;
    }

    void batchSetXlb(cRefV lb_, int ind0=0){
        if(lb_.size() + ind0 > xlb.size())
            xlb.segment(ind0, xlb.size() - ind0) = lb_.head(xlb.size() - ind0);
        else
            xlb.segment(ind0, lb_.size()) = lb_;
    }

    void batchSetXub(cRefV ub_, int ind0=0){
        if(ub_.size() + ind0 > xub.size())
            xub.segment(ind0, xub.size() - ind0) = ub_.head(xub.size() - ind0);
        else
            xub.segment(ind0, ub_.size()) = ub_;
    }

    VX randomGenX(){
        VX mX = VX::Zero(nx);
        mX.setRandom();
        double *xlow = xlb.data(), *xupp = xub.data(), *x = mX.data();
        // move to bound
        if(xlb.size() > 0 && xub.size() > 0){
            for(int i = 0; i < nx; i++){
                if(x[i] < xlow[i] || x[i] > xupp[i]){
                    if(xlow[i] == -1e20)
                        x[i] = xupp[i] - 1e-3*static_cast<double>(rand())/RAND_MAX;
                    else if(xupp[i] == 1e20)
                        x[i] = xlow[i] + 1e-3*static_cast<double>(rand())/RAND_MAX;
                    else{
                        x[i] = (xlow[i] + xupp[i]) / 2.0;
                    }
                }
            }
        }
        return mX;
    }

    // use this function to automatically detect problem size, if they are not specified before.
    void detect_prob_size(int nF, int nG){
        VX rand_x = randomGenX();
        _detect_size(rand_x, nF, nG);
    }

    void detect_prob_size(cRefV x, int nF, int nG) {
        if(x.size() != nx)
            nx = x.size();
        _detect_size(x, nF, nG);
    }

    void _detect_size(cRefV x, int nF, int nG){
        VX F(nF);
        if(nG == 0){
            int nF = operator()(x, F);
            nf = nF;
            grad = false;
            nG = 0;
        }
        else{
            VX G(nG);
            VXi row(nG), col(nG);
            auto rst = operator()(x, F, G, row, col, false, true);  // we needg, but not recording.
            nf = rst.first;
            ProblemFun::nG = rst.second;
            grad = true;
        }
        _allocate_space();
    }

    // use this function to evaluate the problem function once, call __callf__
    VX eval_f(cRefV x){
        VX F(nf);
        operator()(x, F);
        _add_linear_part(F, x);
        return F;
    }

    // use this function to evaluate the problem function once, call __callf__
    std::tuple<VX, VX, VXi, VXi> eval_g(cRefV x){
        VX F(nf), G(nG);
        VXi row(nG), col(nG);
        operator()(x, F, G, row, col, true, true);
        _add_linear_part(F, x);
        return std::make_tuple(F, G, row, col);
    }

    void _add_linear_part(RefV F, cRefV x){
        for(int i = 0; i < Aval.size(); i++) {
            F[Arow[i]] += Aval[i] * x[Acol[i]];
        }
    }

    // use this function to detect the value of nG if the user has implemented __callg__ function
    // maxnG is an integer which should be larger than actual nG + 1 to store enough space
    void detectNg(int maxnG){
        VX x = randomGenX();
        VX F(nf);
        VX G(maxnG);
        VXi row = -1 * VXi::Ones(maxnG);
        VXi col = -1 * VXi::Ones(maxnG);
        operator()(x, F, G, row, col, true, true);
        int i = 0;
        bool found = false;
        for(; i < maxnG; i++){
            if(row(i) == -1 && col(i) == -1){
                found = true;
                break;
            }
        }
        if(found)
            nG = i;
        else
            nG = maxnG;  // this is the only possibility
        grad = true;
    }

    RefV get_lb(){
        return RefV(lb);
    }

    RefV get_ub(){
        return RefV(ub);
    }

    RefV get_xlb(){
        return RefV(xlb);
    }

    RefV get_xub(){
        return RefV(xub);
    }

    RefV get_aval(){
        return RefV(Aval);
    }

    RefVi get_arow(){
        return RefVi(Arow);
    }

    RefVi get_acol(){
        return RefVi(Acol);
    }
};


class SolverConfig{
    typedef std::pair<std::string, int> intOption;
    typedef std::pair<std::string, double> floatOption;
    typedef std::pair<std::string, std::string> pairStringOption;
    typedef std::string stringOption;
public:
    std::vector<intOption> intOptions;
    std::vector<floatOption> floatOptions;
    std::vector<stringOption> stringOptions;
    std::vector<pairStringOption> pairStringOptions;
    void addIntOption(const std::string &nm, int value){
        intOptions.push_back(std::make_pair(nm, value));
    }
    void addFloatOption(const std::string &nm, double value){
        floatOptions.push_back(std::make_pair(nm, value));
    }
    void addPairStringOption(const std::string &nm, const std::string &value){
        pairStringOptions.push_back(std::make_pair(nm, value));
    }
    void addStringOption(const std::string &nm){
        stringOptions.push_back(nm);
    }
};


#endif

