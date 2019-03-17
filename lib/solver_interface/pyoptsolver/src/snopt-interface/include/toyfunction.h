#ifndef TOYFUNCTION
#define TOYFUNCTION

#include "snoptProblem.hh"
#ifdef __cplusplus
extern "C" {
#endif

  void toyusrf_ ( integer    *Status, integer *n,    doublereal x[],
		 integer    *needF,  integer *neF,  doublereal F[],
		 integer    *needG,  integer *neG,  doublereal G[],
		 char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru );

  void toyusrfg_( integer    *Status, integer *n,    doublereal x[],
		 integer    *needF,  integer *neF,  doublereal F[],
		 integer    *needG,  integer *neG,  doublereal G[],
		 char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru);

#ifdef __cplusplus
}
#endif

#endif
