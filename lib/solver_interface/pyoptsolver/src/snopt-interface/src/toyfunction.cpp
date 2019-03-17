#include <stdio.h>
#include "stdlib.h"
#include <string.h>
#include "snoptWrapper.h"

ProblemFun *PROB;
//#define DEBUG

void toyusrf_(integer    *Status, integer *n,    doublereal x[],
         integer    *needF,  integer *neF,  doublereal F[],
         integer    *needG,  integer *neG,  doublereal G[],
         char       *cu,     integer *lencu,
         integer    iu[],    integer *leniu,
         doublereal ru[],    integer *lenru )
{
    //convert x to traj
    MapV c(F, *neF);
    MapV Mx(x, *n);
    if(PROB->getGrad()){
        MapV value(G, *neG);
        VXl row(1), col(1);
        bool rec = false;
        bool needg = true;
        if(*needG == 0)
            needg = false;
        std::pair<int, int> rst = PROB->callg(Mx, c, value, row, col, rec, needg);
        if(rst.first < 0)
            *Status = rst.first;
    }
    else{
        int nf = PROB->callf(Mx, c);
        if(nf < 0)
            *Status = nf;
    }
    return;
}
