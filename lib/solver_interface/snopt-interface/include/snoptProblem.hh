#ifndef SNOPTPROBLEM
#define SNOPTPROBLEM

#define SNOPT76
#ifdef SNOPT76
#include "snoptProblem.hpp"

typedef int integer;
typedef double doublereal;

typedef void (*My_fp)
     ( integer *Status, integer *n,
       doublereal x[],     integer *needF, integer *neF,  doublereal F[],
       integer    *needG,  integer *neG,  doublereal G[],
       char       *cu,     integer *lencu,
       integer    iu[],    integer *leniu,
       doublereal ru[],    integer *lenru );

// class snoptProblem performs problem set-up, initialization,
// and problem-specific variable storage.
namespace snopt{  // use namespace for protection
class snoptProblem : public snoptProblemA {
protected:
  //************************************************************
  // The variables below are set by the user with initialize
  char       Prob[200];

  integer     n, neF;
  integer     ObjRow;
  doublereal  ObjAdd;

  integer     lenA, neA;
  integer    *iAfun, *jAvar;
  doublereal *A;
  integer     lenG, neG;
  integer    *iGfun, *jGvar;

  doublereal *x, *xlow, *xupp, *xmul;
  doublereal *F, *Flow, *Fupp, *Fmul;

  integer    *xstate, *Fstate;

  char       *xnames, *Fnames;
  integer     nxnames, nFnames;

  My_fp usrfun;
  //***********************************************************
  //***********************************************************
  // void userDataSet();
  // void errMsgExit( const char *var );

  integer     inform, fortranStyleObj, fortranStyleAG;
  integer     initCalled;
  integer     minrw, miniw, mincw;
  integer     lenrw, leniw, lencw;
  doublereal *rw;
  integer    *iw;
  char       *cw;

  integer     iSpecs, iSumm, iPrint;
  char        specname[200], printname[200];
  integer     spec_len, prnt_len;

  /*
  void init2zero();
  void setMemory();
  void alloc    ( integer lencw, integer leniw, integer lenrw );
  void realloc  ( integer lencw, integer leniw, integer lenrw );
  void memcpyIn ( char *tcw, integer *tiw, doublereal *trw,
                  integer tlencw, integer tleniw, integer tlenrw);
  void memcpyOut( char *tcw, integer *tiw, doublereal *trw,
                  integer tlencw, integer tleniw, integer tlenrw);
  */
public:
  snoptProblem() {
      initialize("", 1);
      ObjRow = 0;
      ObjAdd = 0;
  };
  ~snoptProblem() {
  };
  /*
  void increment();
  void decrement();
  void computeJac() {
  };
  */

  // int snmema( integer &mincw, integer &miniw, integer &minrw);

  void init() {};

  void computeJac() {
      snoptProblemA::computeJac(neF, n, usrfun,
        x, xlow, xupp,
        iAfun, jAvar, A, neA,
        iGfun, jGvar, neG);
  }

  int solve           ( integer starttype ) {
      integer nS, nInf;
      doublereal sInf;
      inform = snoptProblemA::solve(starttype, neF, n, ObjAdd, ObjRow, usrfun, iAfun, jAvar, A, neA,
                iGfun, jGvar, neG, xlow, xupp, Flow, Fupp, 
                x, xstate, xmul,
                F, Fstate, Fmul,
                nS, nInf, sInf);
      return inform;
  }

  // Functions that set up the problem data:
  void setProblemSize( integer an, integer aneF ) {
      n = an;
      neF = aneF;
  }
  void setObjective  ( integer aObjRow, doublereal aObjAdd ){
      ObjRow = aObjRow;
      ObjAdd = aObjAdd;
  }
  void setA( integer alenA, integer *aiAfun,
                           integer *ajAvar, doublereal *aA )
  {
    //checkSet = checkSet+4;
    lenA  = alenA;
    iAfun = aiAfun;
    jAvar = ajAvar;
    A     = aA;
  }

  void setNeA        ( integer aneA ){
      neA = aneA;
  }
  void setG( integer alenG, integer *aiGfun, integer *ajGvar)
  {
    lenG  = alenG;
    iGfun = aiGfun;
    jGvar = ajGvar;
  }
  void setNeG        ( integer aneG ){
      neG = aneG;
  }
  void setX( doublereal *ax, doublereal *axlow, doublereal *axupp,
                           doublereal *axmul, integer *axstate )
  {
    //checkSet = checkSet+5;
    x      = ax;
    xlow   = axlow;
    xupp   = axupp;
    xmul   = axmul;
    xstate = axstate;
  }
  void setF( doublereal *aF, doublereal *aFlow, doublereal *aFupp,
                           doublereal *aFmul, integer *aFstate )
  {
    F      = aF;
    Flow   = aFlow;
    Fupp   = aFupp;
    Fmul   = aFmul;
    Fstate = aFstate;
  }
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  void setXNames( char *axnames, integer anxnames )
  {
    //checkSet = checkSet+2;
    xnames  = axnames;
    nxnames = anxnames;
  }
  
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  void setFNames( char *aFnames, integer anFnames )
  {
    //checkSet = checkSet+2;
    Fnames  = aFnames;
    nFnames = anFnames;
  }

  void setUserFun    ( My_fp usrfun_in ) {
      usrfun = usrfun_in;
  }
  int getInfo(){
    return inform;
  }
};

}  //namespace snopt
typedef snopt::snoptProblem snoptProblem2;
#else

#endif
#endif
