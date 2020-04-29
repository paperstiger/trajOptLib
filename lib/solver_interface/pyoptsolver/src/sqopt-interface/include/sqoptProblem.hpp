#ifndef SQOPTPROBLEM_H
#define SQOPTPROBLEM_H

#include "sqopt.h"

/* File sqoptProblem.hpp
 *   C++ interface for SQOPT
 *
 */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblem {
protected:
  snoptProblem();
  snoptProblem(const char*name);
  ~snoptProblem();

  void init2zero();

  char    Prob[30];

  int     initCalled, memCalled;

  int     leniw, lenrw;
  double *rw;
  int    *iw;

  int     lenru, leniu;
  double *ru;
  int    *iu;

  void allocI    (int leniw);
  void allocR    (int lenrw);
  void reallocI  (int leniw);
  void reallocR  (int lenrw);

public:
  void setProbName    (const char *Prob);
  void setPrintFile   (const char *prtname);

  int getParameter    (const char *stroptin, char *stroptout);
  int getIntParameter (const char *stropt,   int    &opt);
  int getRealParameter(const char *stropt,   double &opt);
  int setParameter    (const char *stroptin);
  int setIntParameter (const char *stropt,   int     opt);
  int setRealParameter(const char *stropt,   double  opt);

  void setUserI       (int    *iu, int leniu);
  void setUserR       (double *ru, int lenru);
  void setUserspace   (int    *iu, int leniu, double *ru, int lenru);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class sqoptProblem : public snoptProblem {
private:
  isqLog  sqLog;
  void init2zero();

public:
  sqoptProblem();
  sqoptProblem(const char*name);
  ~sqoptProblem();

  void initialize  (const char *prtfile, int summOn);
  int  setSpecsFile(const char *specname);
  void setLog      (isqLog sqLog);2
  void setWorkspace(int m, int n, int neA, int ncObj, int nnH);

  int solve(int starttype, sqFunHx qpHx,
	    int m, int n, int neA, int ncObj, int nnH,
	    int iObj, double ObjAdd,
	    double *A, int *indA, int *locA,
	    double *bl, double *bu, double *cObj,
	    int *eType, int *hs, double *x, double *pi, double *rc,
	    int &nS, int &nInf, double &sInf, double &objective);
};

#endif /* SQOPTPROBLEM_H */
