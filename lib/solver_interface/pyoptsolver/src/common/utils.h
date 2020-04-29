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
#include <set>
#include <tuple>
#include <string>
#include <chrono>
#include <system_error>
#include <exception>
#include "TigerEigen.h"
#include "functionBase.h"


class NotImplementedError : public std::logic_error{
public:
    NotImplementedError(std::string info) : std::logic_error("Not implemented error" + info) {};
};


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


typedef std::chrono::high_resolution_clock high_res_clock;
typedef high_res_clock::time_point time_point;
typedef std::chrono::milliseconds milliseconds;
typedef std::chrono::duration<double> dsec;


// a functor which defines problems
class ProblemFun : public funBase{
public:
    VX lb, ub;
    VX xlb, xub;
    VX Aval;
    VXi Arow, Acol;
    bool ipStyle = false;
    int hess_nnz;
    std::vector<std::tuple<double, double, double> > time_cost_constr_pair;
    VX tmp_y;
    std::vector<int> linear_obj_indices;
    time_point _start;
    bool is_record_time = false;
    bool first_call = true;
    bool DEBUG_VERBOSE = false;

    ProblemFun():funBase(){
        _start = high_res_clock::now();
    }

    ProblemFun(int nx_):funBase(){
        nx = nx_;
    }

    ProblemFun(int nx_, int nf_): funBase(nx_, nf_){
        _allocate_space();
    }

    ProblemFun(int nx_, int nf_, int ng_): funBase(nx_, nf_, ng_){
        _allocate_space();
    }

    void enable_time_record(){
        is_record_time = true;
        tmp_y.resize(nf);
        // scan linear cost and find if there is anything belonging there
        for(int i = 0; i < Aval.size(); i++) 
            if(Arow(i) == 0)
                linear_obj_indices.push_back(i);
    }

    // return the history of time-objective pairs.
    std::tuple<VX, VX, VX> get_time_obj_history() {
        int sz = time_cost_constr_pair.size();
        VX time(sz), cost(sz), constr(sz);
        for(int i = 0; i < sz; i++) {
            auto &pairs = time_cost_constr_pair[i];
            time(i) = std::get<0>(pairs);
            cost(i) = std::get<1>(pairs);
            constr(i) = std::get<2>(pairs);
        }
        return std::make_tuple(time, cost, constr);
    }

    double get_linear_cost(cRefV x) {
        double y = 0;
        for(int i = 0; i < linear_obj_indices.size(); i++)
            y += x(Acol(linear_obj_indices[i])) * Aval(linear_obj_indices[i]);
        return y;
    }

    void set_debug_verbose(bool verbose) {
        DEBUG_VERBOSE = verbose;
    }

    double get_constr_vio(cRefV x, cRefV F) {
        // evaluate violation of constraints
        tmp_y.setZero();
        tmp_y = F;  // first copy, then add linear term
        _add_linear_part(tmp_y, x);
        double vio_sqr = 0;
        for(int i = 0; i < nf; i++) {
            if(tmp_y(i) > ub(i)){
                vio_sqr += pow(tmp_y(i) - ub(i), 2.0);
            }
            else if(tmp_y(i) < lb(i)){
                vio_sqr += pow(lb(i) - tmp_y(i), 2.0);
            }
        }
        return sqrt(vio_sqr);
    }

    int callf(cRefV x, RefV F) {
        if(is_record_time) {
            if(first_call) {
                _start = high_res_clock::now();
                first_call = false;
            }
            int flag = operator()(x, F);
            dsec fs = (high_res_clock::now() - _start);
            double cost = F[0] + get_linear_cost(x);
            double constr_vio = get_constr_vio(x, F);
            if(flag < 0 && time_cost_constr_pair.size() > 0){
                cost = std::get<1>(time_cost_constr_pair.back());
                constr_vio = std::get<2>(time_cost_constr_pair.back());
            }
            time_cost_constr_pair.push_back(std::make_tuple(fs.count(), cost, constr_vio));
            return flag;
        }
        else
            return operator()(x, F);
    }

    std::pair<int, int> callg(cRefV x, RefV F, RefV G, RefVl row, RefVl col, bool rec, bool needg){
        if(is_record_time) {
            if(first_call) {
                _start = high_res_clock::now();
                first_call = false;
            }
            auto ret = operator()(x, F, G, row, col, rec, needg);
            dsec fs = (high_res_clock::now() - _start);
            double cost = F[0] + get_linear_cost(x);
            double constr_vio = get_constr_vio(x, F);
            if(ret.first < 0 && time_cost_constr_pair.size() > 0){
                cost = std::get<1>(time_cost_constr_pair.back());
                constr_vio = std::get<2>(time_cost_constr_pair.back());
            }
            time_cost_constr_pair.push_back(std::make_tuple(fs.count(), cost, constr_vio));
            return ret;
        }
        else
            return operator()(x, F, G, row, col, rec, needg);

    }

    virtual int operator()(cRefV x, RefV F){
#ifdef ENABLEIP
        double obj = evalF(x);
        F[0] = obj;
        MapV g(F.data() + 1, F.size() - 1);
        evalG(x, g);
        return 0;
#else
        throw NotImplementedError("callf");
        return 0;
#endif
    };  // A function to be overwritten by subclass, this is called to evaluate

    virtual std::pair<int, int> operator()(cRefV x, RefV F, RefV G, RefVl row, RefVl col, bool rec, bool needg){
        throw NotImplementedError("callg");
        return std::make_pair(0, 0);
    };  // A function to be overwritten by subclass, this is called for both assigning structure.

#ifdef ENABLEIP
    virtual double evalF(cRefV x) {throw NotImplementedError("evalF"); return 0;};
    virtual bool evalGrad(cRefV x, RefV grad) {throw NotImplementedError("evalGrad"); return true;};
    virtual int evalG(cRefV x, RefV g) {throw NotImplementedError("evalG"); return 0;};
    virtual int evalJac(cRefV x, RefV G, RefVl row, RefVl col, bool rec) {throw NotImplementedError("evalJac"); return 0;};
    virtual int evalHess(cRefV x, double sigma, cRefV lmd, RefV G, RefVl row, RefVl col, bool rec) {throw NotImplementedError("evalHess"); return 0;}
    // easy function calls for debugging purposes
    double ipEvalF(cRefV x) {
        return evalF(x);
    }
    VX ipEvalGrad(cRefV x) {
        VX grad(nx);
        evalGrad(x, grad);
        return grad;
    }
    VX ipEvalG(cRefV x) { // evaluate constraint function values only
        VX g(nf);
        evalG(x, g);
        return g;
    }
    std::tuple<VX, VXl, VXl> ipEvalJac(cRefV x) {  // evaluate Jacobian and return triplet
        VX g(nG);
        VXl row(nG), col(nG);
        evalJac(x, g, row, col, true);
        return std::make_tuple(g, row, col);
    }
#endif

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

    void setA(cRefV val, cRefVl row, cRefVl col){
        Aval = val;
        Arow = row.cast<int>();
        Acol = col.cast<int>();
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
#ifdef ENABLEIP
            if(ipStyle)
                nf = evalGrad(x, F);
            else
#endif
                nf = operator()(x, F);
            grad = false;
            this->nG = 0;
        }
        else{
            VX G(nG);
            VXl row(nG), col(nG);
#ifdef ENABLEIP
            if(ipStyle){
                nf = evalGrad(x, F);
                this->nG = evalJac(x, F, row, col, false);
            }
            else{
#endif
            auto rst = operator()(x, F, G, row, col, false, true);  // we needg, but not recording.
            nf = rst.first;
            ProblemFun::nG = rst.second;
#ifdef ENABLEIP
            }
#endif
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
    std::tuple<VX, VX, VXl, VXl> eval_g(cRefV x){
        VX F(nf), G(nG);
        VXl row(nG), col(nG);
        operator()(x, F, G, row, col, true, true);
        _add_linear_part(F, x);
        return std::make_tuple(F, G, row, col);
    }

    void _add_linear_part(RefV F, cRefV x){
        for(int i = 0; i < Aval.size(); i++) {
            F(Arow(i)) += Aval(i) * x(Acol(i));
        }
    }

    // use this function to detect the value of nG if the user has implemented __callg__ function
    // maxnG is an integer which should be larger than actual nG + 1 to store enough space
    void detectNg(int maxnG){
        VX x = randomGenX();
        VX F(nf);
        VX G(maxnG);
        VXl row = -1 * VXl::Ones(maxnG);
        VXl col = -1 * VXl::Ones(maxnG);
#ifdef ENABLEIP
        if(ipStyle) {
            evalJac(x, G, row, col, false);
        }
        else{
#endif
        operator()(x, F, G, row, col, true, true);
#ifdef ENABLEIP
        }
#endif
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

    // check if any A or G overlaps, this could be dangerous if it indeed occurs
    // I have to use vector and set for this
    // I shall still work in integer space, not long yet
    // The user has to guarantee ng is already set, like by detect_ng
    void overlap_check() {
        assert(nG >= 0);
        VX x = randomGenX();
        VX F(nf);
        VXl row = VXl::Zero(nG), col = row;
        VX val = VX::Zero(nG);
        operator()(x, F, val, row, col, true, true);
        // in most cases, we start from nonlinear and then check linear
        int A_nnz = Arow.size();
        if(true) {
            std::vector<std::set<int> > pool(nf);
            // loop over nonlinear entry
            for(int i = 0; i < nG; i++) {
                auto &seti = pool[col(i)];
                if(seti.count(row(i)) != 0)
                    printf("Nonlin row %d col %d overlaps\n", (int)row(i), (int)col(i));
                else
                    seti.insert(row(i));
            }
            // loop over linear entry
            for(int i = 0; i < A_nnz; i++) {
                auto &seti = pool[Acol(i)];
                if(seti.count(Arow(i)) != 0)
                    printf("Lin row %d col %d overlaps\n", (int)Arow(i), (int)Acol(i));
                else
                    seti.insert(Arow(i));
            }
        }
        /*
        else{
            std::vector<std::set<int> > pool(nx);
            // loop over nonlinear entry
            for(int i = 0; i < nG; i++) {
                auto &seti = pool[row(i)];
                if(seti.count(col(i)) != 0)
                    printf("Nonlin row %d col %d overlaps\n", (int)row(i), (int)col(i));
                else
                    seti.insert(col(i));
            }
            // loop over linear entry
            for(int i = 0; i < A_nnz; i++) {
                auto &seti = pool[Arow(i)];
                if(seti.count(Acol(i)) != 0)
                    printf("Lin row %d col %d overlaps\n", (int)Arow(i), (int)Acol(i));
                else
                    seti.insert(Acol(i));
            }
        }
        */
    }

    // return a list of variable reference so painless modifying is possible
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
    virtual void setMajorIter(int iter) {
        throw NotImplementedError("set major iter");
    }
    virtual void setOptTol(double tol) {
        throw NotImplementedError("set opt tol");
    }
    virtual void setFeaTol(double tol) {
        throw NotImplementedError("set fea tol");
    }
    virtual int setPrintLevel(int lvl) {
        throw NotImplementedError("set print level");
    }
    virtual void enableDerivCheck(int lvl=3) {
        throw NotImplementedError("set deriv check");
    }
};


#endif
