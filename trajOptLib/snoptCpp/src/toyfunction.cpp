#include <stdio.h>
#include "stdlib.h"
#include <string.h>
#include "snoptWrapper.h"

ProblemFun *PROB;
//#define DEBUG

int toyusrf_(integer    *Status, integer *n,    doublereal x[],
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
        VXi row(1), col(1);
        PROB->operator()(Mx, c, value, row, col, false, true);
#ifdef DEBUG
        std::cout << "grad mode\n";
#endif
    }
    else{
        PROB->operator()(Mx, c);
    }
    return 0;
}
