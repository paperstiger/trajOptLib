#include <assert.h>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include "sqopt.h"
#include "sqoptEZ.h"

using namespace std;

// sqopt function part
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
sqoptProblem::sqoptProblem() {
  init2zero();
  sprintf(Prob, "%8s", "        ");

  allocI(500);
  allocR(500);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
sqoptProblem::sqoptProblem(const char *name) {
  init2zero();

  sprintf(Prob, "%8s", name);

  allocI(500);
  allocR(500);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
sqoptProblem::~sqoptProblem() {
  f_sqend(iw, leniw, rw, lenrw);

  delete []rw;  delete []iw;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::init2zero() {
  initCalled = 0; memCalled = 0;

  leniw = 0; lenrw = 0;
  iw    = 0; rw    = 0;

  leniu = 0; lenru = 0;
  iu    = 0; ru    = 0;

  // sqopt
  sqLog = 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::allocI(int aleniw) {
  // Reset work array lengths.
  // Allocate new memory for work arrays.
  leniw = aleniw;
  iw    = new int[leniw];
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::allocR(int alenrw) {
  // Reset work array lengths.
  // Allocate new memory for work arrays.
  lenrw = alenrw;
  rw    = new double[lenrw];
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::reallocI(int aleniw) {
  int  tleniw = leniw;
  int    *tiw = iw;

  // Allocate new memory
  allocI(aleniw);

  // Copy old workspace into new.
  int mleniw = leniw < tleniw ? leniw : tleniw;
  memcpy(iw, tiw, mleniw*sizeof(int));

  // Delete temporary work arrays
  delete []tiw;

  setIntParameter((char*)"Total int workspace", leniw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::reallocR(int alenrw) {
  int  tlenrw = lenrw;
  double *trw = rw;

  // Allocate new memory
  allocR(alenrw);

  // Copy old workspace into new.
  int mlenrw = lenrw < tlenrw ? lenrw : tlenrw;
  memcpy(rw, trw, mlenrw*sizeof(double));

  // Delete temporary work arrays
  delete []trw;

  setIntParameter((char*)"Total real workspace   ", lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setProbName(const char *name) {
  sprintf(Prob, "%8s", name);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setPrintFile(const char *prtname) {
  assert(initCalled == 1);

  int len = strlen(prtname);
  f_sqsetprint(prtname, len, iw, leniw, rw, lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::setParameter(const char *stropt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_sqset(stropt, stropt_len, &errors, iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::getParameter(const char *stroptin, char *stroptout) {
  assert(initCalled == 1);

  int errors;
  int stroptin_len  = strlen(stroptin);
  int stroptout_len = strlen(stroptout);

  f_sqgetc(stroptin, stroptin_len, stroptout, stroptout_len,
	    &errors, iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::setIntParameter(const char *stropt, int opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);

  f_sqseti(stropt, stropt_len, opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::getIntParameter(const char *stropt, int &opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_sqgeti(stropt, stropt_len, &opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::setRealParameter(const char *stropt, double opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_sqsetr(stropt, stropt_len, opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::getRealParameter(const char *stropt, double &opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_sqgetr(stropt, stropt_len, &opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setUserI(int *aiu, int aleniu) {
  leniu = aleniu;
  iu    = aiu;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setUserR(double *aru, int alenru) {
  lenru = alenru;
  ru    = aru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setUserspace  (int*aiu,     int aleniu,
				   double *aru, int alenru) {
  leniu = aleniu;
  iu    = aiu;

  lenru = alenru;
  ru    = aru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::initialize(const char*prtfile, int summOn) {
  int len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("  SQOPT  C++ interface  2.0.0  ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_sqinit(prtfile, len, summOn, iw, leniw, rw, lenrw);
  initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::setSpecsFile(const char *specname) {
  assert(initCalled == 1);

  int inform, len = strlen(specname);
  f_sqspec(specname, len, &inform, iw, leniw, rw, lenrw);
  if(inform != 101){
    printf("Warning: unable to find specs file %s \n", specname);
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setLog(isqLog asqLog) {
  sqLog  = asqLog;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void sqoptProblem::setWorkspace(int m, int n, int neA, int ncObj, int nnH) {
  assert(initCalled == 1);

  int inform, miniw, minrw;

  f_sqmem(&inform, m, n, neA, ncObj, nnH,
	   &miniw, &minrw, iw, leniw, rw, lenrw);

  if (miniw > leniw) { reallocI (miniw); }
  if (minrw > lenrw) { reallocR (minrw); }

  memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int sqoptProblem::solve(int starttype, sqFunHx qpHx,
			int m, int n, int neA, int ncObj, int nnH,
			int iObj, double ObjAdd,
			double *A, int *indA, int *locA,
			double *bl, double *bu, double *cObj,
			int *eType, int *hs, double *x, double *pi, double *rc,
			int &nS, int &nInf, double &sInf, double &objective) {
  assert(initCalled == 1);

  int inform, sniObj, miniw, minrw;

  if (memCalled == 0) { setWorkspace(m, n, neA, ncObj, nnH); }

  for (int i = 0; i < neA; i++) {
    indA[i]++;
  }
  for (int i = 0; i <= n; i++) {
    locA[i]++;
  }

  sniObj = iObj+1;

  f_snkerq(starttype, qpHx, sqLog,
       m, n, neA, ncObj, nnH,
       sniObj, ObjAdd, Prob,
       A, indA, locA, bl, bu, cObj,
       eType, hs, x, pi, rc,
       &inform, &nS, &nInf, &sInf, &objective,
       &miniw, &minrw, iu, leniu, ru, lenru,
       iw, leniw, rw, lenrw);

  for (int i = 0; i < neA; i++) {
    indA[i]--;
  }
  for (int i = 0; i <= n; i++) {
    locA[i]--;
  }

  return inform;
}


// initialize static member variable H
double* SqoptEZ::Hv = NULL; // values
int* SqoptEZ::Hr = NULL; // row indices
int* SqoptEZ::Hc = NULL; // column indices
double SqoptEZ::infBnd = 1e20; // infinity
