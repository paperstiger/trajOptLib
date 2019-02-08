/*
 * ipoptWrapper.cpp
 * Copyright (C) 2019 motion <motion@motion-MS-7971>
 *
 * Distributed under terms of the MIT license.
 */

#include "ipoptWrapper.h"
#include <vector>


using namespace Ipopt;
using namespace std;


optResult solve_problem(ProblemFun &prob, IpoptConfig &config, cRefV x0) {
    // create instance
    // IpoptAdapter *prob_in = new IpoptAdapter(prob, x0.data());
    // TNLP *mynlp = prob_in;
    SmartPtr<TNLP> mynlp = new IpoptAdapter(prob, x0.data());
    // create app
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);
    // some options waiting to be done
    app->Options()->SetNumericValue("tol", 1e-9);
    // Intialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    config.setup_app(&(*app));
    optResult result;
    if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        result.flag = (int) status + 1;
        return result;
    }
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);
    if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
        result.flag = (int) status + 1;
    }
    else {
        printf("\n\n*** The problem FAILED!\n");
        result.flag = 0;
    }
    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.
    IpoptAdapter *prob_in = (IpoptAdapter*)(Ipopt::GetRawPtr(mynlp));
    result.val = prob_in->cost;
    result.sol = prob_in->sol;
    result.c = prob_in->c;
    result.lmd = prob_in->lmd;
    result.xmul = prob_in->mu;
    return result;
}
