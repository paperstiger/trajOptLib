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
#include <string>
#include "TigerTools/TigerEigen.h"
#include "toyfunction.hh"
#include "snfilewrapper.hh"
#include "snopt.hh"
#include "functionBase.h"
#include "snoptProblem.hh"


class optResult{
    public:
        int flag;
        double val;
        VX sol;
        VX c;
        optResult(){}
};


// a functor which defines problems
class ProblemFun : public funBase{
    public:
        VX lb, ub;
        VX xlb, xub;
        VX Aval;
        VXi Arow, Acol;
        using funBase::funBase;

        virtual void operator()(cRefV x, RefV F) = 0;  // A function to be overwritten by subclass, this is called to evaluate

        virtual void operator()(cRefV x, RefV F, RefV G, RefVi row, RefVi col, bool rec, bool needg) = 0;  // A function to be overwritten by subclass, this is called for both assigning structure.

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
            /*
            int nnz = val.size();
            Aval.resize(nnz);
            Arow.resize(nnz);
            Acol.resize(nnz);
            */
            Aval = val;
            Arow = row;
            Acol = col;
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
                            x[i] = xupp[i] - 1e-3*rand()/RAND_MAX;
                        else if(xupp[i] == 1e20)
                            x[i] = xlow[i] + 1e-3*rand()/RAND_MAX;
                        else{
                            x[i] = (xlow[i] + xupp[i]) / 2.0;
                        }
                    }
                }
            }
            return mX;
        }
};


class snoptConfig{
    typedef std::pair<std::string, int> intOption;
    typedef std::pair<std::string, double> floatOption;
    typedef std::pair<std::string, std::string> stringOption;
public:
    std::string name = std::string("Toy");
    std::string printFile;
    std::vector<intOption> intOptions;
    std::vector<floatOption> floatOptions;
    std::vector<stringOption> stringOptions;
    int printlevel = 0;
    int verifylevel = 0;
    int majoriterlimit = 0;
    int minoriterlimit = 0;
    int iterationslimit = 0;
    double optTol = 1e-6;
    double feaTol = 1e-6;
    snoptConfig(){
#ifdef DEBUG
        std::cout << "Entering construction of snopt config\n";
#endif
    }
    void addIntOption(const std::string &nm, int value){
        intOptions.push_back(std::make_pair(nm, value));
    }
    void addFloatOption(const std::string &nm, double value){
        floatOptions.push_back(std::make_pair(nm, value));
    }
    void addStringOption(const std::string &nm, const std::string &value){
        stringOptions.push_back(std::make_pair(nm, value));
    }
};


extern ProblemFun *PROB;

class snoptWrapper{
private:
    integer n, neF, neA, neG, lenA, lenG;
    integer *iAfun, *jAvar, *iGfun, *jGvar, *xstate, *Fstate;
    doublereal *mA, *G, *x, *xlow, *xupp, *xmul, *F, *Flow, *Fupp, *Fmul;
    char *xnames, *Fnames;
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
    snoptWrapper(ProblemFun *pfun, snoptConfig *cfg=nullptr){
#ifdef DEBUG
        std::cout << "Entering construction\n";
#endif
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
#ifdef DEBUG
                std::cout << "A of size 0\n";
#endif
            }
            else{
                lenA = pfun->Aval.size();
                useA = true;
#ifdef DEBUG
                std::cout << "length of A is " << lenA << std::endl;
#endif
            }
        }
        else
            lenA  = n * neF;
#ifdef DEBUG
        std::cout << "Allocate space\n";
#endif
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
#ifdef DEBUG
        std::cout << "set bounds\n";
#endif
        setxbound();
        setfbound();

#ifdef DEBUG
        std::cout << "get gradient information\n";
#endif
        /***Gradients, depends on how we are defining the problem***/
        if(pfun->getGrad()){
            MapV mV(F, neF);
            MapVi row(iGfun, lenG), col(jGvar, lenG);
            MapV Gvalue(G, lenG);
            ranGenX();  // randomly generate an initial guess
            MapV mX(x, n);
            pfun->operator()(mX, mV, Gvalue, row, col, true, true);
        }
        // Set the problem, what's left is the initial guess
        if(pfun->getGrad()){
            neA = lenA;
            neG = lenG;
            if(!useA)
                neA = 0;
            ToyProb.setNeA         ( neA );
            ToyProb.setNeG         ( neG );
#ifdef DEBUG
            std::cout << "neA, neG = " << neA << " " <<  neG << std::endl;
#endif
        }
#ifdef DEBUG
        std::cout << "set properties\n";
#endif
        ToyProb.setProblemSize( n, neF );
        ToyProb.setObjective  ( ObjRow, ObjAdd );
        ToyProb.setA          ( lenA, iAfun, jAvar, mA );
        ToyProb.setG          ( lenG, iGfun, jGvar );
        ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
        ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );
        ToyProb.setXNames     ( xnames, nxnames );
        ToyProb.setFNames     ( Fnames, nFnames );
        ToyProb.setProbName   ( "Toy0" );
        ToyProb.setUserFun    ( toyusrf_ );
        if(!pfun->getGrad()){
            ranGenX();
            ToyProb.computeJac();
        }
        if(pfun->getGrad())
            ToyProb.setIntParameter( "Derivative option", 1 );
        else
            ToyProb.setIntParameter( "Derivative option", 0 );
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
            ToyProb.setIntParameter( "Minor print level", 0 );
            ToyProb.setIntParameter( "Verify level", snpcfg->verifylevel);
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
            for(auto &scfg : snpcfg->stringOptions){
                std::string option = std::get<0>(scfg) + std::string(" ") + std::get<1>(scfg);
                ToyProb.setParameter(option.c_str());
            }
        }
#ifdef DEBUG
        std::cout << "finish construct\n";
#endif
    }

    void setX(const double *xin){
        for (int i = 0; i < n; ++i)
            x[i] = xin[i];
    }
    void copyX(double *xout) const{
        for (int i = 0; i < n; ++i)
            xout[i] = x[i];
    }
    void setxbound(){
        // Set the upper and lower bounds.
        if(prob->xlb.size() == 0){
#ifdef DEBUG
            std::cout << "xlb of size 0 , n = " << n << std::endl;
#endif
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
#ifdef DEBUG
        MapV Mxlb(xlow, n);
        MapV Mxub(xupp, n);
        std::cout << "lowX, uppX " << Mxlb << "\n" << Mxub << std::endl;
#endif
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
                    x[i] = (xlow[i] + xupp[i]) / 2.0;
                }
            }
        }
#ifdef DEBUG
        MapV mXf(x, n);
        std::cout << "Ran x " << mXf << std::endl;
#endif
    }

    void setfbound(){
        //set upper and lower bound for constraint
        if(prob->lb.size() == 0){
#ifdef DEBUG
            std::cout << "lb of size 0 , neF = " << neF << std::endl;
#endif
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
#ifdef DEBUG
        MapV Mflb(Flow, neF);
        MapV Mfub(Fupp, neF);
        std::cout << "lowF, uppF " << Mflb << "\n" << Mfub << std::endl;
#endif
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
    void setIntWorkspace(int iws){
        ToyProb.setIntWs(iws);
    }
    void setRealWorkspace(int fws){
        ToyProb.setRealWs(fws);
    }
    //Solve the problem. Modified Aug 13 2017, no need to warm start
    int solve(){
        int appendcol = 0;
        //change constraints, here
        setxbound();
        for (int i = 0; i < neF; ++i){
            F[i] = 0;
            Fstate[i] = 0;
            Fmul[i] = 0;
        }
        for (int i = 0; i < n; ++i){
            xstate[i] = 0;
            xmul[i] = 0;
        }
#ifdef DEBUG
        std::cout << "solve from random guess\n";
#endif
        //Generate initial guess using the method in TrajOpt.h
        ranGenX();
        ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );
        ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
        ToyProb.solve( 0 );
        return ToyProb.getInfo();
    };
    //Solve the problem. given the solution, lmdF
    int solve(double *_x, double *_Fmul = NULL){
#ifdef DEBUG
        std::cout << "solve with guess\n";
#endif
        setxbound();//As always, do this
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
        ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );
        ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
        ToyProb.solve( 0 );
        return ToyProb.getInfo();
    };

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
            VXi row(1), col(1);
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

