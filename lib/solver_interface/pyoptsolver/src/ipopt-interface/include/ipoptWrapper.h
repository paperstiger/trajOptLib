/*
 * ipoptWrapper.h
 * Copyright (C) 2019 motion <motion@motion-MS-7971>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef IPOPTWRAPPER_H
#define IPOPTWRAPPER_H

#include <coin/IpTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include "utils.h"
#include <vector>
#include <string>


using std::string;


extern bool VERBOSE;


inline void print_array(const char msg[], const double *x, int n) {
    if(VERBOSE){
        printf("%s ", msg);
        for(int i = 0; i < n; i++)
            printf("%f ", x[i]);
        printf("\n");
    }
}

inline void print_array(const char msg[], const int *x, int n) {
    if(VERBOSE) {
        printf("%s ", msg);
        for(int i = 0; i < n; i++)
            printf("%d ", x[i]);
        printf("\n");
    }
}

namespace Ipopt{
/**
 * @brief An adaptor that takes ProblemFun type and wrap into Ipopt
 *
 * Given an optimization probleme defined accordinG to SNOPT style,
 * implement a class that is suitable to be solved by Ipopt
 *
 */
class IpoptAdapter : public TNLP {
public:
    IpoptAdapter(ProblemFun &prob_, const double *x0_) : prob(prob_), x0(x0_) {
        // at this stage, call to fill in rows, cols
        int n = prob.nx, m = prob.nf;
        F.resize(m);
        sol.resize(n);
        mu.resize(n);
        lmd.resize(m);
        c.resize(m);
        G.resize(prob.nG);
        rows.resize(prob.nG);
        cols.resize(prob.nG);

        cMapV Mx(x0_, n);

        if(prob.ipStyle) {
#ifdef ENABLEIP
            if(prob_.grad)
                prob.evalJac(Mx, G, rows, cols, true);
#endif
        }
        else{
            prob.operator()(Mx, F, G, rows, cols, true, true);
            std::cout << "nG = " << prob.nG << std::endl;
            for(int i = 0; i < prob.nG; i++)
                if(rows(i) == 0)
                    g_obj_index.push_back(i);
            std::cout << "nA = " << prob.Aval.size() << std::endl;
            for(int i = 0; i < prob.Aval.size(); i++)
                if(prob.Arow(i) == 0)
                    a_obj_index.push_back(i);
        }
    }

    VX sol, lmd, mu, c;
    double cost;


    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
            Index& nnz_h_lag, IndexStyleEnum& index_style) {
        n = prob.nx;
        m = prob.nf;
        nnz_jac_g = prob.nG + prob.Aval.size();
        if(prob.hess_nnz > 0)
            nnz_h_lag = prob.hess_nnz;
        else
            nnz_h_lag = n * (n + 1) / 2;
        index_style = C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, double* x_lower, double* x_upper,
            Index m, double* g_l, double* g_u)
    {
        for (uint c=0; c < prob.nx; ++c) {
            x_lower[c] = prob.xlb(c);
            x_upper[c] = prob.xub(c);
        }

        // specific bounds dependinG on equality and inequality constraints
        for (uint c=0; c<prob.nf; ++c) {
            g_l[c] = prob.lb(c);
            g_u[c] = prob.ub(c);
        }
        print_array("xlow", x_lower, prob.nx);
        print_array("xup", x_upper, prob.nx);
        print_array("cl", g_l, prob.nf);
        print_array("cu", g_u, prob.nf);
        return true;
    }

    virtual bool get_starting_point(Index n, bool init_x, double* x,
            bool init_z, double* z_L, double* z_U,
            Index m, bool init_lambda,
            double* lambda)
    {
        // Here, we assume we only have startinG values for x
        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        MapV(x, prob.nx) = cMapV(x0, prob.nx);

        print_array("start x0", x, n);

        return true;
    }

    virtual void receive_new_x(const double *x, bool needg) {
        if(VERBOSE)
            std::cout << "receive new x\n";
        cMapV Mx(x, prob.nx);
        VXl row(1), col(1);
        F.array() = 0;
        // prob(Mx, F, G, row, col, false, needg);
        prob.callg(Mx, F, G, row, col, false, needg);
        last_has_grad = needg;
        // plus the array
        for(int i = 0; i < prob.Aval.size(); i++) {
            F(prob.Arow(i)) += prob.Aval(i) * x[prob.Acol(i)];
        }
    }

    virtual bool eval_f(Index n, const double* x, bool new_x, double& obj_value)
    {
        if(VERBOSE)
            std::cout << "evalf\n";
        if(prob.ipStyle){
#ifdef ENABLEIP
            cMapV Mx(x, n);
            double value = prob.evalF(Mx);
            // printf("obj = %f\n", value);
            obj_value = value;
#endif
            return true;
        }
        else{
            if(new_x)
                receive_new_x(x, false);
            obj_value = F(0);
            print_array("x=", x, n);
            // printf("f = %f\n", obj_value);
            return true;
        }
    }

    virtual bool eval_grad_f(Index n, const double* x, bool new_x, double* grad_f)
    {
        if(VERBOSE)
            std::cout << "grad f\n";
        if(prob.ipStyle) {
#ifdef ENABLEIP
            cMapV Mx(x, n);
            MapV grad(grad_f, n);
            bool ret = prob.evalGrad(Mx, grad);
            print_array("grad=", grad_f, n);
            return ret;
#endif
            return true;
        }
        else{
            if(new_x)
                receive_new_x(x, false);
            for(int i = 0; i < n; i++)
                grad_f[i] = 0;
            for(auto &idx : g_obj_index)
                grad_f[cols[idx]] = G[idx];
            for(auto &idx : a_obj_index)
                grad_f[prob.Acol[idx]] = prob.Aval[idx];
            print_array("x=", x, n);
            print_array("grad=", grad_f, n);
            return true;
        }
    }

    virtual bool eval_g(Index n, const double* x, bool new_x, Index m, double* g)
    {
        if(VERBOSE)
            std::cout << "eval g\n";
        if(prob.ipStyle) {
#ifdef ENABLEIP
            cMapV Mx(x, n);
            MapV gfun(g, m);
            gfun.array() = 0;
            int ret = prob.evalG(Mx, gfun);
            for(int i = 0; i < prob.Aval.size(); i++)
                gfun(prob.Arow(i)) += prob.Aval(i) * x[prob.Acol(i)];
            print_array("g=", g, m);
#endif
            return true;
        }
        else{
            if(new_x || !last_has_grad)
                receive_new_x(x, true);
            MapV(g, m) = F;
            print_array("g=", g, m);
            return true;
        }
    }

    virtual bool eval_jac_g(Index n, const double* x, bool new_x,
            Index m, Index nele_jac, Index* iRow, Index *jCol,
            double* values)
    {
        if(VERBOSE)
            std::cout << "jac g\n";
        // defines the positions of the nonzero elements of the jacobian
        if (values == NULL) {
            // first copy nonlinear part
            for(int i = 0; i < prob.nG; i++) {
                iRow[i] = rows(i);
                jCol[i] = cols(i);
            }
            // then copy linear part, if any
            int nnz_a = prob.Aval.size();
            for(int i = 0; i < nnz_a; i++) {
                iRow[prob.nG + i] = prob.Arow(i);
                jCol[prob.nG + i] = prob.Acol(i);
            }
            print_array("row=", iRow, prob.nG + prob.Arow.size());
            print_array("col=", jCol, prob.nG + prob.Arow.size());
        }
        else {
            if(prob.ipStyle){
#ifdef ENABLEIP
                cMapV Mx(x, n);
                MapV G(values, prob.nG);
                VXl row(0), col(0);
                int ret = prob.evalJac(Mx, G, row, col, false);
                MapV(values + prob.nG, prob.Aval.size()) = prob.Aval;
                print_array("jac=", values, prob.nG + prob.Aval.size());
#endif
            }
            else{
                // only gets used if "jacobian_approximation finite-difference-values" is not set
                if(new_x || !last_has_grad)
                    receive_new_x(x, true);
                // first copy nonlinear part
                MapV(values, prob.nG) = G;
                MapV(values + prob.nG, prob.Aval.size()) = prob.Aval;
                print_array("jac=", values, prob.nG + prob.Aval.size());
            }
        }
        return true;
    }

    virtual bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number *lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values) {
        std::cout << "Entering hessian evaluation\n";
        cMapV Mx(x, n);
        cMapV lmd(lambda, m);
        MapV G(values, nele_hess);
        MapVi row(iRow, nele_hess);
        MapVi col(jCol, nele_hess);
#ifdef ENABLEIP
        if (values == NULL) {
            VXl hrows(nele_hess), hcols(nele_hess);
            prob.evalHess(Mx, obj_factor, lmd, G, hrows, hcols, true);
            row = hrows.cast<int>();
            col = hcols.cast<int>();
        }
        else {
            VXl hrows(0), hcols(0);
            prob.evalHess(Mx, obj_factor, lmd, G, hrows, hcols, false);
        }
#endif
        return true;
    }

    /** This method is called when the algorithm is complete so the TNLP can
    * store/write the solution */
    virtual void finalize_solution(SolverReturn status,
                                 Index n, const double* x, const double* z_L, const double* z_U,
                                 Index m, const double* g, const double* lambda,
                                 double obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq){
        sol = cMapV(x, n);
        lmd = cMapV(lambda, m);
        cost = obj_value;
        mu = cMapV(z_U, n) - cMapV(z_L, n);
        c = cMapV(g, m);
    }

private:
    ProblemFun &prob;  // this stores the vector
    const double *x0;  // store the pointer to initial guess
    VX F, G;
    std::vector<int> a_obj_index, g_obj_index;  // record which entries are in first row
    VXl rows, cols;
    bool last_has_grad = false;

};

};


class IpoptConfig : public SolverConfig{
public:
    int print_level = 5;
    int print_frequency_iter = 10;
    int max_iter = 1000;
    string linear_solver = "mumps";
    string hessian_approximation = "limited-memory";
    string jacobian_approximation = "exact";
    bool deriv_check_enabled = false;
    string derivative_test = "none";
    double tol = 1e-6;
    double constr_viol_tol = 1e-6;
    double max_cpu_time = 1e6;

    void setup_app(Ipopt::IpoptApplication *app) {
        auto opt = app->Options();
        opt->SetIntegerValue("print_level", print_level);
        opt->SetIntegerValue("print_frequency_iter", print_frequency_iter);
        opt->SetIntegerValue("max_iter", max_iter);
        opt->SetStringValue("linear_solver", linear_solver.c_str());
        opt->SetStringValue("hessian_approximation", hessian_approximation.c_str());
        opt->SetStringValue("jacobian_approximation", jacobian_approximation.c_str());
        if(deriv_check_enabled)
            opt->SetStringValue("derivative_test", "first-order");  //derivative_test.c_str());
        else
            opt->SetStringValue("derivative_test", "none");
        opt->SetStringValue("print_timing_statistics", "no");
        opt->SetStringValue("print_user_options", "no");
        opt->SetNumericValue("tol", tol);
        opt->SetNumericValue("constr_viol_tol", constr_viol_tol);
        opt->SetNumericValue("max_cpu_time", max_cpu_time);
        for(auto &icfg : intOptions){
            opt->SetIntegerValue(std::get<0>(icfg).c_str(), std::get<1>(icfg));
        }
        for(auto &fcfg : floatOptions){
            opt->SetNumericValue(std::get<0>(fcfg).c_str(), std::get<1>(fcfg));
        }
        for(auto &scfg : pairStringOptions){
            opt->SetStringValue(std::get<0>(scfg).c_str(), std::get<1>(scfg).c_str());
        }
    }

    virtual void setMajorIter(int iter) {
        max_iter = iter;
    }
    virtual void setOptTol(double tol) {
        this->tol = tol;
    }
    virtual void setFeaTol(double tol) {
        constr_viol_tol = tol;
    }
    virtual int setPrintLevel(int lvl) {
        print_level = lvl;
    }
    virtual void enableDerivCheck(int lvl=3) {
        if(lvl <= 0)
            deriv_check_enabled = false;
	else
            deriv_check_enabled = true;
    }

    void setPrintFreq(int freq) {
        print_frequency_iter = freq;
    }

    void setLinearSolver(std::string &stri) {
        linear_solver = stri;
    }

    void enableExactHessian() {
        hessian_approximation = "exact";
    }

    void enableFDJacobian() {
        jacobian_approximation = "finite-difference-values";
    }

};

optResult solve_problem(ProblemFun &prob, IpoptConfig &config, cRefV x0);


class IpoptSolver {
public:
    IpoptSolver(ProblemFun &prob_, IpoptConfig &config_) {
        prob = &prob_;
        config = &config_;
        if(!prob->grad) {
            config->enableFDJacobian();
        }
    }

    optResult solve_rand() {
        VX x0 = prob->randomGenX();
        return solve_problem(*prob, *config, x0);
    }

    optResult solve_guess(cRefV x0) {
        return solve_problem(*prob, *config, x0);
    }

private:
    ProblemFun *prob;
    IpoptConfig *config;
};

#endif /* !IPOPTWRAPPER_H */
