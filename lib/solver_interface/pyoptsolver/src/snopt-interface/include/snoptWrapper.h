/* Construct the SNOPT problem
 */
#ifndef SNPSOLVE_H
#define SNPSOLVE_H

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <tuple>
#include <string>
#include "time.h"
#include "TigerEigen.h"
#include "toyfunction.h"
#include "snoptProblem.hh"
#include "utils.h"


class snoptConfig : public SolverConfig{
public:
    std::string name = std::string("Toy");
    std::string printFile;
    int printlevel = 1;
    int minorprintlevel = 0;
    int verifylevel = 0;
    int majoriterlimit = 0;
    int minoriterlimit = 0;
    int iterationslimit = 0;
    double optTol = 1e-6;
    double feaTol = 1e-6;
    void setIterLimit(int iter) {
        iterationslimit = iter;
    }
    void setMinorIter(int iter) {
        minoriterlimit = iter;
    }
    virtual void setMajorIter(int iter) {
        majoriterlimit = iter;
    }
    virtual void setOptTol(double tol) {
        optTol = tol;
    }
    virtual void setFeaTol(double tol) {
        feaTol = tol;
    }
    virtual int setPrintLevel(int lvl) {
        printlevel = lvl;
    }
    virtual void enableDerivCheck(int lvl=3) {
        verifylevel = lvl;
    }
};


extern ProblemFun *PROB;

class snoptWrapper{
private:
    bool is_setup;
protected:
    integer n, neF, neA, neG, lenA, lenG;
    integer *iAfun, *jAvar, *iGfun, *jGvar, *xstate, *Fstate;
    doublereal *mA, *G, *x, *xlow, *xupp, *xmul, *F, *Flow, *Fupp, *Fmul;
    char *xnames, *Fnames;
    int cached_lw = -1, cached_rw = -1;
public:
    ProblemFun *prob;
    snoptConfig *snpcfg;
    snoptProblem2 ToyProb;
    double *getFmul() const {return Fmul;};
    double *getX() const {return x;};
    int getVarNum() const {return n;}
    int getFNum() const {return neF;}
    int getneF() const{return neF;}
    double getObj() const {return F[0];};
    double *getF() const{return F;}
    double *getFlow() const{return Flow;}
    double *getFupp() const{return Fupp;}
    //Construction function
    snoptWrapper(ProblemFun *pfun, snoptConfig *cfg=nullptr) : is_setup(false) {
        prob = pfun;
        snpcfg = cfg;
        PROB = pfun;
        n = pfun->getNx();
        neF = pfun->getNf();
        // linear constraint is basically given up here
        bool useA = false;
        if(pfun->getGrad()){
            if(pfun->Aval.size() == 0){
                lenA = 1;  //not sure if 0 works fine
            }
            else{
                lenA = pfun->Aval.size();
                useA = true;
            }
        }
        else
            lenA  = n * neF;
        iAfun = new integer[lenA];
        jAvar = new integer[lenA];
        mA  = new doublereal[lenA];
        if(pfun->getGrad() && pfun->Aval.size() > 0){
            MapV mapA(mA, lenA);
            MapVi mapiA(iAfun, lenA);
            MapVi mapjA(jAvar, lenA);
            mapA = pfun->Aval;
            mapiA = pfun->Arow;
            mapjA = pfun->Acol;
        }
        if(pfun->getGrad())
            lenG = pfun->getNg();
        else
            lenG = n * neF;
        iGfun = new integer[lenG];
        jGvar = new integer[lenG];
        G     = new doublereal[lenG];

        x      = new doublereal[n];
        xlow   = new doublereal[n];
        xupp   = new doublereal[n];
        xmul   = new doublereal[n];
        xstate = new    integer[n];

        F      = new doublereal[neF];
        Flow   = new doublereal[neF];
        Fupp   = new doublereal[neF];
        Fmul   = new doublereal[neF];
        Fstate = new integer[neF];

        integer nxnames = 1;
        integer nFnames = 1;
        xnames = new char[nxnames*8];
        Fnames = new char[nFnames*8];

        integer    ObjRow = 0;
        doublereal ObjAdd = 0;

        /***set bound on x and f***/
        setxbound();
        setfbound();

        if(pfun->getGrad()){
            neA = lenA;
            neG = lenG;
            if(!useA)
                neA = 0;
            ToyProb.setNeA         ( neA );
            ToyProb.setNeG         ( neG );
        }

        ToyProb.setProblemSize( n, neF );
        ToyProb.setObjective  ( ObjRow, ObjAdd );
        ToyProb.setA          ( lenA, iAfun, jAvar, mA );
        ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
        ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );
        ToyProb.setXNames     ( xnames, nxnames );
        ToyProb.setFNames     ( Fnames, nFnames );
        ToyProb.setProbName   ( "Gao" );
        ToyProb.setUserFun    ( toyusrf_ );

        if(snpcfg == nullptr){
            ToyProb.setIntParameter( "Verify level", 0 );
            ToyProb.setIntParameter( "Major print level", 0 );
            ToyProb.setIntParameter( "Major iterations limit", 3000);
            ToyProb.setIntParameter( "Minor print level", 0 );
            ToyProb.setRealParameter("Major optimality tolerance", 1e-6);
            ToyProb.setRealParameter("Major feasibility tolerance", 1e-6);
        }
        else{
            ToyProb.setProbName(snpcfg->name.c_str());
            ToyProb.setRealParameter("Major optimality tolerance", snpcfg->optTol);
            ToyProb.setRealParameter("Major feasibility tolerance", snpcfg->feaTol);
            ToyProb.setIntParameter( "Major print level", snpcfg->printlevel);
            ToyProb.setIntParameter( "Minor print level", snpcfg->minorprintlevel);
            ToyProb.setIntParameter( "Verify level", snpcfg->verifylevel);
            if(snpcfg->printlevel == 0)
                ToyProb.set_summary(0);  // turn off summary in this scenario
            if(snpcfg->printFile.size() > 0){
                setPrintFile(snpcfg->printFile);
            }
            if(snpcfg->majoriterlimit > 0){
                ToyProb.setIntParameter("Major iterations limit", snpcfg->majoriterlimit);
            }
            if(snpcfg->minoriterlimit > 0){
                ToyProb.setIntParameter("Minor iterations limit", snpcfg->minoriterlimit);
            }
            if(snpcfg->iterationslimit > 0){
                ToyProb.setIntParameter("Iterations limit", snpcfg->iterationslimit);
            }
            for(auto &icfg : snpcfg->intOptions){
                ToyProb.setIntParameter(std::get<0>(icfg).c_str(), std::get<1>(icfg));
            }
            for(auto &fcfg : snpcfg->floatOptions){
                ToyProb.setRealParameter(std::get<0>(fcfg).c_str(), std::get<1>(fcfg));
            }
            for(auto &option : snpcfg->stringOptions){
                ToyProb.setParameter(option.c_str());
            }
        }
    }

    void problem_setup() {
        /***Gradients, depends on how we are defining the problem***/
        if(prob->getGrad()){
            MapV mV(F, neF);
            MapVi row(iGfun, lenG), col(jGvar, lenG);
            MapV Gvalue(G, lenG);
            MapV mX(x, n);
            VXl lrow = row.cast<long>();
            VXl lcol = col.cast<long>();
            prob->operator()(mX, mV, Gvalue, lrow, lcol, true, true);
            row = lrow.cast<int>();
            col = lcol.cast<int>();
            ToyProb.setG          ( lenG, iGfun, jGvar );
        }
        if(prob->getGrad())
            ToyProb.setIntParameter( "Derivative option", 1 );
        else
            ToyProb.setIntParameter( "Derivative option", 0 );
        is_setup = true;
    }

    void setX(const double *xin){
        for (int i = 0; i < n; ++i)
            x[i] = xin[i];
    }

    void copyX(double *xout) const{
        for (int i = 0; i < n; ++i)
            xout[i] = x[i];
    }

    void copyF(double *fout) const{
        for(int i = 0; i < neF; i++)
            fout[i] = F[i];
    }

    void copyLmd(double *lmd) const{
        for(int i = 0; i < neF; i++)
            lmd[i] = Fmul[i];
    }

    void copyXmul(double *xmu) const{
        for(int i = 0; i < n; i++)
            xmu[i] = xmul[i];
    }

    void setxbound(){
        // Set the upper and lower bounds.
        if(prob->xlb.size() == 0){
            for (int i = 0; i < n; ++i){
                xlow[i] = -1e20; xupp[i] = 1e20; xstate[i] = 0;
            }
        }
        else{
            if(prob->xlb.size() != n){
                printf("Incorrect x bound size\n");
                exit(0);
            }
            for (int i = 0; i < n; ++i){
                xlow[i] = prob->xlb(i); xupp[i] = prob->xub(i); xstate[i] = 0;
            }
        }
    }

    // generate a random guess for certain problems, all variables are random
    void ranGenX(){
        MapV mX(x, n);
        mX.setRandom();
        // move to bound
        for(int i = 0; i < n; i++){
            if(x[i] < xlow[i] || x[i] > xupp[i]){
                if(xlow[i] == -1e20)
                    x[i] = xupp[i];
                else if(xupp[i] == 1e20)
                    x[i] = xlow[i];
                else{
                    x[i] = (xlow[i] + xupp[i]) / 2.0 + static_cast<double>(rand())/RAND_MAX * 1e-6;
                }
            }
        }
    }

    void setfbound(){
        //set upper and lower bound for constraint
        if(prob->lb.size() == 0){
            for(int i = 0; i < neF; i++){
                Flow[i] = 0;
                Fupp[i] = 0;
                Fmul[i] = 0;
                Fstate[i] = 0;
            }
            Flow[0] = -1e20;
            Fupp[0] = 1e20;
        }
        else{
            if(prob->lb.size() != neF){
                printf("Incorrect constraint size\n");
                exit(0);
            }
            for (int i = 0; i < neF; ++i){
                Flow[i] = prob->lb(i);
                Fupp[i] = prob->ub(i);
                Fmul[i] = 0;
                Fstate[i] = 0;
            }
        }
    }

    void setOptTol(double tol){
        ToyProb.setRealParameter("Major optimality tolerance", tol);
    }

    void setPrintFile(std::string &fnm){
        FILE *fp = fopen(fnm.c_str(), "wt");  // wipe out all contents
        if(!fp){
            std::cout << "Fail to reset file " << fnm << std::endl;
        }
        fclose(fp);
        ToyProb.setPrintFile  (fnm.c_str());
    }
    void setFeaTol(double tol){
        ToyProb.setRealParameter("Major feasibility tolerance", tol);
    }
    void setIntOption(std::string &fnm, int opt){
        ToyProb.setIntParameter(fnm.c_str(), opt);
    }
    void setDoubleOption(std::string &fnm, double opt){
        ToyProb.setRealParameter(fnm.c_str(), opt);
    }
    void setMajorIter(int val){
        ToyProb.setIntParameter( "Major iterations limit", val);
    }
    void setMinorIter(int val){
        ToyProb.setIntParameter( "Minor iterations limit", val);
    }

    void setWorkspace(int iws, int fws){
        ToyProb.setWorkspace(neF, n, neA, neG);
        ToyProb.reallocI(iws);
        ToyProb.reallocR(fws);
    }

    void update_problem() {
        setxbound();
        setfbound();
    }

    //Solve the problem. Modified Aug 13 2017, no need to warm start
    int solve(){
        int appendcol = 0;
        for (int i = 0; i < neF; ++i){
            F[i] = 0;
            Fstate[i] = 0;
            Fmul[i] = 0;
        }
        for (int i = 0; i < n; ++i){
            xstate[i] = 0;
            xmul[i] = 0;
        }
        //Generate initial guess using the method in TrajOpt.h
        ranGenX();
        if(!is_setup)
            problem_setup();
        ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );
        ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
        if(prob->getGrad())
            ToyProb.solve( 0 );
        else
            ToyProb.nograd_solve( 0 );
        return ToyProb.getInfo();
    };
    //Solve the problem. given the solution, lmdF
    int solve(const double *_x, double *_Fmul = NULL){
        for (int i = 0; i < neF; ++i){
            F[i] = 0;
            Fstate[i] = 0;
            if(_Fmul)
                Fmul[i] = _Fmul[i];
            else
                Fmul[i] = 0;
        }
        for (int i = 0; i < n; ++i){
            xstate[i] = 0;
            xmul[i] = 0;
            x[i] = _x[i];
        }
        if(!is_setup)
            problem_setup();
        ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );
        ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
        if(prob->getGrad())
            ToyProb.solve( 0 );
        else
            ToyProb.nograd_solve( 0 );
        return ToyProb.getInfo();
    };

    // In cases that do not alter problem structure, just solve for more iterations
    int solve_more(int iter) {
        setMajorIter(iter);
        ToyProb.solve( 2 );
        return ToyProb.getInfo();
    }

    // given an initial guess and a problem without nonlinear constraints, perform 
    // search using snopt. First an initial iteration number of search is performed,
    // then increamental steps are performed and cost are checked.
    // If either abstol or reltol criteria is satisfied, we are done.
    // Returns final cost, list of costs at search step, list of time spent
    std::pair<VX, VX> obj_search(double *_x, int iter, int step, int step_num, double abstol, double reltol) {
        clock_t tic, toc;
        std::vector<double> cost, timestamp;
        tic = clock();
        setMajorIter(iter);
        int flag = solve(_x);
        double cur_obj = getObj();  // this might be incorrect, I guess.
        toc = clock();
        cost.push_back(cur_obj);
        timestamp.push_back((double)(toc - tic) / CLOCKS_PER_SEC);
        if(flag == 32) {
            for(int i = 0; i < step_num; i++) {
                flag = solve_more(step);
                double new_obj = getObj();
                cost.push_back(new_obj);
                toc = clock();
                timestamp.push_back((double)(toc - tic) / CLOCKS_PER_SEC);
                if((cur_obj - new_obj < abstol) || (cur_obj - new_obj) < reltol * cur_obj)
                    break;
                if (flag != 32)
                    break;
                cur_obj = new_obj;
            }
        }
        VX costs(cost.size());
        for(int i = 0; i < cost.size(); i++)
            costs(i) = cost[i];
        VX times(timestamp.size());
        for(int i = 0; i < timestamp.size(); i++)
            times(i) = timestamp[i];
        return std::make_pair(costs, times);
    }

    int getInfo(){
        return ToyProb.getInfo();
    }

    /* Given a guess of solution, evaluate violation of constraints */
    VX fEval(double *_x){
        MapV c(F, neF);
        MapV value(G, lenG);
        MapV Mx(_x, n);
        if(prob->getGrad()){
            MapM g(G, 1, lenG);
            VXl row(1), col(1);
            prob->operator()(Mx, c, value, row, col, false, false);
        }
        else{
            prob->operator()(Mx, c);
        }
        // do not forget those linear terms, I do it myself
        for(int i = 0; i < neA; i++){
            F[iAfun[i]] += mA[i] * _x[jAvar[i]];
        }
        VX out = c;
        return out;
    }

    VX getMul(){
        MapV lmd(Fmul, neF);
    }

    ~snoptWrapper(){
        delete []iAfun;  delete []jAvar;  delete []mA;
        delete []iGfun;  delete []jGvar;  delete []G;

        delete []x;      delete []xlow;   delete []xupp;
        delete []xmul;   delete []xstate;

        delete []F;      delete []Flow;   delete []Fupp;
        delete []Fmul;   delete []Fstate;

        delete []xnames; delete []Fnames;
    };
};
#endif

