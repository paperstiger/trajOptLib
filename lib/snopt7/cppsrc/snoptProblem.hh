#ifndef SNOPTPROBLEM
#define SNOPTPROBLEM

#include "snopt.hh"


// class snoptProblem performs problem set-up, initialization,
// and problem-specific variable storage.
class snoptProblem {
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
  void userDataSet();
  void errMsgExit( const char *var );

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

  void init2zero();
  virtual void setMemory();
  void alloc    ( integer lencw, integer leniw, integer lenrw );
  void realloc  ( integer lencw, integer leniw, integer lenrw );
  void memcpyIn ( char *tcw, integer *tiw, doublereal *trw,
                  integer tlencw, integer tleniw, integer tlenrw);
  void memcpyOut( char *tcw, integer *tiw, doublereal *trw,
                  integer tlencw, integer tleniw, integer tlenrw);
public:
   snoptProblem();
  ~snoptProblem();
  void increment();
  void decrement();
  void computeJac();

  int snmema( integer &mincw, integer &miniw, integer &minrw);

  void init();
  void getParameter    ( const char *stroptin, char *stroptout );
  void getIntParameter ( const char *stropt,   integer    &opt );
  void getRealParameter( const char *stropt,   doublereal &opt );
  void setParameter    ( const char *stroptin );
  void setIntParameter ( const char *stropt,   integer     opt );
  void setRealParameter( const char *stropt,   doublereal  opt );
  void setPrintFile    ( const char printname[] );
  void setSpecsFile    ( const char specname[] );
  void solve           ( integer starttype );

  // Functions that set up the problem data:
  void setProblemSize( integer n, integer neF );
  void setObjective  ( integer ObjRow, doublereal ObjAdd );
  void setA          ( integer lenA, integer *iAfun, integer *jAvar, doublereal *A );
  void setNeA        ( integer neA );
  void setG          ( integer lenG, integer *iGfun, integer *jGvar );
  void setNeG        ( integer neG );
  void setX          ( doublereal *x, doublereal *xlow, doublereal *xupp,
                       doublereal *xmul, integer *xstate );
  void setF          ( doublereal *F, doublereal *Flow, doublereal *Fupp,
                       doublereal *Fmul, integer *Fstate );
  void setXNames     ( char *xnames, integer nxnames );
  void setFNames     ( char *Fnames, integer nFnames );
  void setProbName   ( const char *Prob );
  void setUserFun    ( My_fp usrfun );
};

class snoptProblem2: public snoptProblem{
private:
    int intws = 0, realws = 0;
public:
    int getInfo(){
        return inform;
    }

    void setIntWs(int iws){
        intws = iws;
    }

    void setRealWs(int fws){
        realws = fws;
    }

    void setMemory(){
      int memoryGuess;
      memoryGuess = this->snmema(mincw, miniw, minrw);
      if(miniw < intws)
          miniw = intws;
      if(minrw < realws)
          minrw = realws;
      if ( mincw > lencw | miniw > leniw | minrw > lenrw ) {
        // Reallocate memory while retaining the values set in sninit_
        this->realloc( mincw, miniw, minrw );
        // Save the lengths of the new work arrays.
        this->setIntParameter((char*)"Total real workspace   ", lenrw );
        this->setIntParameter((char*)"Total integer workspace", leniw );

        // Did we have to guess values of neA and neG for snmema()
        if ( memoryGuess == 1 ) {
          this->computeJac();
          memoryGuess = this->snmema(mincw, miniw, minrw);
          assert( memoryGuess == 0 );
          this->realloc( mincw, miniw, minrw );
          this->setIntParameter((char*)"Total real workspace   ", lenrw );
          this->setIntParameter((char*)"Total integer workspace", leniw );
        }
      }
    }
};
#endif
