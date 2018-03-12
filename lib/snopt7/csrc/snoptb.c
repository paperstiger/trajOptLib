/* ./src/snoptb.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__130 = 130;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__3 = 3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  snoptb.f  --- the Basic interface for snOpt. */

/*     snOpt   snOptB   snKerB */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snopt_(char *start, integer *m, integer *n, integer *ne, 
	integer *nname, integer *nncon, integer *nnobj, integer *nnjac, 
	integer *iobj, doublereal *objadd, char *prob, U_fp funcon, U_fp 
	funobj, doublereal *jcol, integer *indj, integer *locj, doublereal *
	bl, doublereal *bu, char *names, integer *hs, doublereal *x, 
	doublereal *pi, doublereal *rc, integer *info, integer *mincw, 
	integer *miniw, integer *minrw, integer *ns, integer *ninf, 
	doublereal *sinf, doublereal *obj, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen start_len, ftnlen prob_len, ftnlen names_len, ftnlen cu_len, 
	ftnlen cw_len)
{
    extern /* Subroutine */ int snlog_(), sqlog_(), snlog2_();
    extern /* Subroutine */ int snkerb_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, char *, U_fp, U_fp, U_fp, U_fp, U_fp, U_fp, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    char *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    extern /* Subroutine */ int snstop_();

/*     ================================================================== */
/*     snOpt is a Fortran subroutine for constrained nonlinear */
/*     optimization.  The constraints take the form */

/*                            (   x  ) */
/*                      bl <= (      ) <= bu, */
/*                            ( F(x) ) */

/*     where bl and bu are constant lower and upper bounds. */

/*     o If all constraints are linear, F = J x for some sparse matrix J. */

/*     o If all constraints are nonlinear, F = fCon(x) for some vector */
/*       fCon of smooth functions. */

/*     o In general, there is a mixture of constraints of the form */
/*                      ( fCon(x1) +  J2 x2 ), */
/*                      (   J3 x1  +  J4 x2 ) */
/*       where the nonlinear variables x1 must appear first as shown. */

/*     o fCon(x1) and (optionally) its partial derivatives J1(x) are set */
/*       in subroutine funcon (see below). */

/*     o The matrices J2, J3, J4 and the sparsity pattern of J1(x) are */
/*       entered column-wise in the arrays Jcol, indJ, locJ (below). */

/*     o Internally, the constraints are converted into the form */

/*           fCon(x1) +  J2 x2  - s1      = 0,     bl <= ( x ) <= bu */
/*             J3 x1  +  J4 x2       - s2 = 0            ( s ) */

/*       where s = (s1,s2)  and the components of (x,s) are the */
/*       variables and slacks respectively. */

/*     ------------------------------------------------------------------ */
/*     NOTE: Before calling SNOPT, your calling program must call: */
/*     call snInit( iPrint, iSumm, */
/*    &             cw, lencw, iw, leniw, rw, lenrw ) */
/*     This sets the default values of the optional parameters. You can */
/*     also alter the default values of iPrint and iSumm before SNOPT */
/*     is used.  iPrint = 0, etc, is OK. */
/*     ------------------------------------------------------------------ */

/*     ON ENTRY: */

/*     Start   specifies how a starting basis (and certain other items) */
/*             are to be obtained. */
/*             Start = 'Cold' means that Crash should be used to choose */
/*                      an initial basis, unless a basis file is given */
/*                      via Old basis, Insert or Load in the Specs file. */
/*             Start = 'Basis file' means the same (but is more */
/*                      meaningful in the latter case). */
/*             Start = 'Warm' means that a basis is already defined in hs */
/*                      (probably from an earlier call). */
/*             Start = 'Hot' or 'Hot FHS' implies a hot start. */
/*                      hs defines a basis and an earlier call has */
/*                      defined certain other things that should also be */
/*                      kept.  The problem dimensions and the arrays */
/*                      cw(*), iw(*), and rw(*) must not have changed. */
/*                      F refers to the LU factors of the basis. */
/*                      H refers to the approximate reduced Hessian R. */
/*                      S refers to column and row scales. */
/*             Start = 'Hot H' (for example) means that only the Hessian */
/*                      is defined. */

/*     m       is the number of slacks (i.e., general constraints). */
/*             For LP, QP or LC  problems this means the number of rows */
/*             in the constraint matrix J. */
/*             m > 0. */

/*             For problems with no general constraints, set m = 1 and */
/*             impose the constraint that the sum of the variables */
/*             must lie between plus and minus infinity. This gives */
/*             J  one ``free row'' that will not alter the solution. */

/*     n       is the number of variables, excluding slacks. */
/*             For LP problems, this is the number of columns in J. */
/*             n > 0. */

/*     ne      is the number of nonzero entries in J (including the */
/*             Jacobian for any nonlinear constraints). */
/*             ne gt 0. */

/*     nName   is the number of column and row names provided in the */
/*             array  Names.  If nName = 1, there are NO names. */
/*             Generic names will be used in the printed solution. */
/*             Otherwise, nName = n+m and all names must be provided. */

/*     nnCon   is the number of nonlinear constraints. */
/*             nnCon ge 0. */

/*             a nonzero nnCon defines the row dimension of the */
/*             constraint Jacobian J1(x) defined in subroutine funcon. */

/*     nnObj   is the number of nonlinear Objective variables. */
/*             nnObj ge 0. */

/*     nnJac   is the number of nonlinear Jacobian variables. */
/*             If nnCon = 0, nnJac = 0. */
/*             if nnCon > 0, nnJac > 0. */

/*             a nonzero nnJac defines the column dimension of the */
/*             constraint Jacobian J1(x) defined in subroutine funcon. */

/*     iObj    says which row of J is a free row containing a linear */
/*             objective vector  c  (iObj = 0 if none). */
/*             iObj = 0  or  nnCon < iObj le m. */

/*     ObjAdd  is a constant that will be added to the objective. */
/*             Typically ObjAdd = 0.0d+0. */

/*     Prob    is an 8-character name for the problem, used in the */
/*             output.  A blank name can be assigned if necessary. */

/*     Jcol(ne) is the constraint Jacobian J, stored column-wise.  Every */
/*             element of Jcol(*) must be assigned a value.  Elements */
/*             in the nonlinear part (see NOTE 2 below) can be any dummy */
/*             value (e.g., zero) since they are initialized by SNOPT at */
/*             the first point that is feasible with respect to the */
/*             linear constraints.  The linear part of Jcol(*) must */
/*             contain the constant Jacobian elements. */

/*     indJ(ne)  is the list of row indices for each nonzero in Jcol(*). */

/*     locJ(n+1) is a set of pointers to the beginning of each column of */
/*             the constraint matrix within Jcol(*) and indJ(*). */
/*             Must have locJ(1) = 1 and locJ(n+1) = ne+1. */

/*  NOTES:  1. If the problem has a nonlinear objective, the first nnObj */
/*             columns of Jcol and indJ belong to the nonlinear objective */
/*             variables. Subroutine userfg deals with these variables. */

/*          2. If the problem has nonlinear constraints, the first nnJac */
/*             columns of Jcol and indJ belong to the nonlinear Jacobian */
/*             variables, and the first nnCon rows of Jcol and indJ */
/*             belong to the nonlinear constraints.  Subroutine funcon */
/*             deals with these variables and constraints. */

/*          3. If nnObj > 0 and nnJac > 0, the two sets of */
/*             nonlinear variables overlap.  The total number of */
/*             nonlinear variables is nnL = max( nnObj, nnJac ). */

/*          4. The Jacobian forms the top left corner of Jcol and indJ. */
/*             If a Jacobian column j (1 le j le nnJac) contains */
/*             any entries Jcol(k), indJ(k) associated with nonlinear */
/*             constraints (1 le indJ(k) le nnCon), those entries must */
/*             come before any other (linear) entries. */

/*          5. The row indices indJ(k) for a column may be in any order */
/*             (subject to Jacobian entries appearing first). */
/*             Subroutine funcon must define Jacobian entries in the */
/*             same order. */

/*     bl(n+m) is the lower bounds on each variable (x,s). */

/*     bu(n+m) is the upper bounds on each variable (x,s). */

/*     Names(nName) is an character*8 array. */
/*             If nName =  1, Names is not used.  The printed solution */
/*             will use generic names for the columns and row. */
/*             If nName = n+m, Names(j) should contain an 8 character */
/*             name of the jth variable (j = 1, n+m). */
/*             If j = n+i, the jth variable is the ith row. */

/*     hs(n+m) sometimes contains a set of initial states for each */
/*             variable (x, s).  See the following NOTES. */

/*     x(n+m)  is a set of initial values for each variable (x, s). */

/*  NOTES:  1. If Start = 'Cold' or 'Basis file' and a BASIS file */
/*             of some sort is to be input */
/*             (an OLD BASIS file, INSERT file or LOAD file), */
/*             hs and x need not be set at all. */

/*          2. Otherwise, hs(1:n) must be defined for a cold start. */
/*             If nothing special is known about the problem, or if */
/*             there is no wish to provide special information, */
/*             you may set hs(j) = 0, x(j) = 0.0d+0 for all j=1:n. */
/*             All variables will be eligible for the initial basis. */

/*             Less trivially, to say that variable j will probably */
/*             be equal to one of its bounds, */
/*             set hs(j) = 4 and x(j) = bl(j) */
/*             or  hs(j) = 5 and x(j) = bu(j) as appropriate. */

/*          3. For Cold starts with no basis file, a Crash procedure */
/*             is used to select an initial basis.  The initial basis */
/*             matrix will be triangular (ignoring certain small */
/*             entries in each column). */
/*             The values hs(j) = 0, 1, 2, 3, 4, 5 have the following */
/*             meaning: */

/*             hs(j)    State of variable j during Crash */

/*             0, 1, 3  Eligible for the basis.  3 is given preference. */
/*             2, 4, 5  Ignored. */

/*             After Crash, hs(j) = 2 entries are made superbasic. */
/*             Other entries not selected for the basis are made */
/*             nonbasic at the value x(j) if bl(j) <= x(j) <= bu(j), */
/*             or at the value bl(j) or bu(j) closest to x(j). */

/*          4. For Warm and Hot starts, all of hs(1:n+m) is assumed to be */
/*             set to the values 0, 1, 2 or 3 from some previous call. */

/*     pi(m)   contains an estimate of the vector of Lagrange multipliers */
/*             (shadow prices) for the NONLINEAR constraints.  The first */
/*             nnCon components must be defined.  If nothing is known */
/*             about lambda, set pi(i) = 0.0d+0, i = 1:nnCon. */

/*     nS      need not be specified for Cold starts, */
/*             but should retain its value from a previous call */
/*             when a Warm or Hot start is used. */


/*     ON EXIT: */

/*     hs(n+m) is the final state vector: */

/*                hs(j)    State of variable j    Normal value of x(j) */

/*                  0      nonbasic               bl(j) */
/*                  1      nonbasic               bu(j) */
/*                  2      superbasic             Between bl(j) and bu(j) */
/*                  3      basic                  ditto */

/*             Very occasionally there may be nonbasic variables for */
/*             which x(j) lies strictly between its bounds. */
/*             If nInf = 0, basic and superbasic variables may be outside */
/*             their bounds by as much as the Feasibility tolerance. */
/*             Note that if Scale is specified, the Feasibility tolerance */
/*             applies to the variables of the SCALED problem. */
/*             In this case, the variables of the original problem may be */
/*             as much as 0.1 outside their bounds, but this is unlikely */
/*             unless the problem is very badly scaled. */

/*     x(n+m)  contains the final variables and slacks (x, s). */

/*     pi(m)   is the vector of Lagrange multipliers (shadow prices) */
/*             for the general constraints. */

/*     rc(n+m) is a vector of reduced costs: rc = g - (J -I)'pi, where g */
/*             is the gradient of the objective if x is feasible */
/*             (or the gradient of the Phase-1 objective otherwise). */
/*             If nInf = 0, the last m entries are pi. */

/*     INFO    says what happened; see the User's Guide. */
/*             The possible values are as follows: */

/*             INFO    Meaning */

/*                0    finished successfully */
/*                1       optimality conditions satisfied */
/*                2       feasible point found */
/*                3       requested accuracy could not be achieved */

/*               10    the problem appears to be infeasible */
/*               11       infeasible linear constraints */
/*               12       infeasible linear equalities */
/*               13       nonlinear infeasibilities minimized */
/*               14       infeasibilities minimized */

/*               20    the problem appears to be unbounded */
/*               21       unbounded objective */
/*               22       constraint violation limit reached */

/*               30    resource limit error */
/*               31       iteration limit reached */
/*               32       major iteration limit reached */
/*               33       the superbasics limit is too small */

/*               40    terminated after numerical difficulties */
/*               41       current point cannot be improved */
/*               42       singular basis */
/*               43       cannot satisfy the general constraints */
/*               44       ill-conditioned null-space basis */

/*               50    error in the user-supplied functions */
/*               51       incorrect objective  derivatives */
/*               52       incorrect constraint derivatives */

/*               60    undefined user-supplied functions */
/*               61       undefined function at the first feasible point */
/*               62       undefined function at the initial point */
/*               63       unable to proceed into undefined region */

/*               70    user requested termination */
/*               71       terminated during function evaluation */
/*               72       terminated during constraint evaluation */
/*               73       terminated during objective evaluation */
/*               74       terminated from monitor routine */

/*               80    insufficient storage allocated */
/*               81       work arrays must have at least 500 elements */
/*               82       not enough character storage */
/*               83       not enough integer storage */
/*               84       not enough real storage */

/*               90    input arguments out of range */
/*               91       invalid input argument */
/*               92       basis file dimensions do not match this problem */

/*              140    system error */
/*              141       wrong no of basic variables */
/*              142       error in basis package */

/*     mincw   says how much character storage is needed to solve the */
/*             problem.  If INFO = 82, the work array cw(lencw) was */
/*             too small.  snOptA may be called again with lencw suitably */
/*             larger than mincw. */

/*     miniw   says how much integer storage is needed to solve the */
/*             problem.  If INFO = 83, the work array iw(leniw) was too */
/*             small.  snOptA may be called again with leniw suitably */
/*             larger than miniw.  (The bigger the better, since it is */
/*             not certain how much storage the basis factors need.) */

/*     minrw   says how much real storage is needed to solve the */
/*             problem.  If INFO = 84, the work array rw(lenrw) was too */
/*             small.  (See the comments above for miniw.) */

/*     nS      is the final number of superbasics. */

/*     nInf    is the number of infeasibilities. */

/*     sInf    is the sum    of infeasibilities. */

/*     Obj     is the value of the nonlinear part of the objective. */
/*             If nInf = 0, Obj includes the nonlinear objective if any. */
/*             If nInf > 0, Obj is just the linear objective if any. */

/*     cu(lencu), iu(leniu), ru(lenru)  are character, integer and real */
/*             arrays of USER workspace.  These arrays are available to */
/*             pass data to the user-defined routines funobj and funcon. */
/*             If no workspace is required, you can either use dummy */
/*             arrays for cu, iu and ru, or use cw, iw and rw */
/*             (see below). */

/*     cw(lencw), iw(leniw), rw(lenrw)  are character*8, integer and real */
/*             arrays of workspace used by SNOPT. */
/*             lencw  should be about at least 500. */
/*             leniw  should be about max( 500, 20(m+n) ) or larger. */
/*             lenrw  should be about max( 500, 40(m+n) ) or larger. */

/*     SNOPT is maintained by Philip E. Gill, */
/*     Dept of Mathematics, University of California, San Diego. */

/*     LUSOL is maintained by Michael A. Saunders, */
/*     Systems Optimization Laboratory, */
/*     Dept of Management Science & Engineering, Stanford University. */

/*     12 Nov 1994: Workspace separated into iw(*) and rw(*). */
/*     08 Aug 1996: First Min Sum version. */
/*     17 Jul 1997: first thread-safe version. */
/*     26 Jul 1997: User workspace added. */
/*     02 Oct 1997: Character workspace added. */
/*     11 Oct 1998: Added facility to combine funobj and funcon. */
/*     22 Dec 2002: Added input argument checking. */
/*     24 Dec 2002: Reincarnated as snOptB. */
/*     02 Jan 2003: Call sqOpt for LP's. */
/*     20 Jan 2003: QP-QN added. */
/*     31 Jul 2003: snPRNT and snEXIT adopted. */
/*     15 Oct 2004: snSTOP adopted. */
/*     01 Sep 2007: Sticky parameters added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --bu;
    --bl;
    --locj;
    --indj;
    --jcol;
    names -= 8;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    snkerb_(start, m, n, ne, nname, nncon, nnobj, nnjac, iobj, objadd, prob, (
	    U_fp)funcon, (U_fp)funobj, (U_fp)snlog_, (U_fp)snlog2_, (U_fp)
	    sqlog_, (U_fp)snstop_, &jcol[1], &indj[1], &locj[1], &bl[1], &bu[
	    1], names + 8, &hs[1], &x[1], &pi[1], &rc[1], info, mincw, miniw, 
	    minrw, ns, ninf, sinf, obj, cu + 8, lencu, &iu[1], leniu, &ru[1], 
	    lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, start_len, (
	    ftnlen)8, (ftnlen)8, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snopt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snOpt */
/* Subroutine */ int snoptb_(char *start, integer *m, integer *n, integer *ne,
	 integer *nname, integer *nncon, integer *nnobj, integer *nnjac, 
	integer *iobj, doublereal *objadd, char *prob, U_fp funcon, U_fp 
	funobj, doublereal *jcol, integer *indj, integer *locj, doublereal *
	bl, doublereal *bu, char *names, integer *hs, doublereal *x, 
	doublereal *pi, doublereal *rc, integer *info, integer *mincw, 
	integer *miniw, integer *minrw, integer *ns, integer *ninf, 
	doublereal *sinf, doublereal *obj, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen start_len, ftnlen prob_len, ftnlen names_len, ftnlen cu_len, 
	ftnlen cw_len)
{
    extern /* Subroutine */ int snlog_(), sqlog_(), snlog2_();
    extern /* Subroutine */ int snkerb_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, char *, U_fp, U_fp, U_fp, U_fp, U_fp, U_fp, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    char *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    extern /* Subroutine */ int snstop_();

/*     ================================================================== */
/*     snOptB is a Fortran subroutine for constrained nonlinear */
/*     optimization.  The constraints take the form */

/*                            (   x  ) */
/*                      bl <= (      ) <= bu, */
/*                            ( F(x) ) */


/*     24 Dec 2002: First version of snOptB. */
/*     01 Sep 2007: Sticky parameters added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --bu;
    --bl;
    --locj;
    --indj;
    --jcol;
    names -= 8;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    snkerb_(start, m, n, ne, nname, nncon, nnobj, nnjac, iobj, objadd, prob, (
	    U_fp)funcon, (U_fp)funobj, (U_fp)snlog_, (U_fp)snlog2_, (U_fp)
	    sqlog_, (U_fp)snstop_, &jcol[1], &indj[1], &locj[1], &bl[1], &bu[
	    1], names + 8, &hs[1], &x[1], &pi[1], &rc[1], info, mincw, miniw, 
	    minrw, ns, ninf, sinf, obj, cu + 8, lencu, &iu[1], leniu, &ru[1], 
	    lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, start_len, (
	    ftnlen)8, (ftnlen)8, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snoptb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snOptB */
/* Subroutine */ int snkerb_(char *start, integer *m, integer *n, integer *ne,
	 integer *nname, integer *nncon, integer *nnobj, integer *nnjac, 
	integer *iobj, doublereal *objadd, char *prob, U_fp funcon, U_fp 
	funobj, U_fp snlog, U_fp snlog2, U_fp sqlog, U_fp snstop, doublereal *
	jcol, integer *indj, integer *locj, doublereal *bl, doublereal *bu, 
	char *names, integer *hs, doublereal *x, doublereal *pi, doublereal *
	rc, integer *info, integer *mincw, integer *miniw, integer *minrw, 
	integer *ns, integer *ninf, doublereal *sinf, doublereal *obj, char *
	cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, 
	integer *lenru, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen start_len, ftnlen prob_len, 
	ftnlen names_len, ftnlen cu_len, ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal x0[1];
    static integer nb, lx0, nx0, nnh;
    static doublereal rhs[1];
    static integer lkx, nkx;
    static char str[80];
    static integer nnh0;
    static char str2[80];
    static doublereal fobj;
    static integer lenr, ngqp, maxr, maxs, nrhs;
    static logical gotr;
    extern /* Subroutine */ int sqhx_(), s0fgb_();
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer lenx0;
    extern /* Subroutine */ int s8map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer nrhs0;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer lgobj, ngobj, llocg, nlocg;
    extern /* Subroutine */ int s2mem0_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, ftnlen);
    static integer nlocj, maxcw;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxiw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxrw, ngobj0;
    extern /* Subroutine */ int s3argb_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    char *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen), s1file_(
	    integer *, integer *, integer *), s2bmap_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), s8gloc_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s1time_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *), s8dflt_(integer *, integer *,
	     integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen), s1perm_(integer *, 
	    integer *), s3prtb_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *), 
	    s8gsiz_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    extern /* Subroutine */ int s8qphx_();
    extern /* Subroutine */ int s5solv_(integer *, char *, integer *, U_fp, 
	    U_fp, U_fp, logical *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, doublereal *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), s8solv_(integer *, 
	    char *, integer *, U_fp, U_fp, U_fp, U_fp, U_fp, U_fp, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static integer negcon;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer nmajor, minmax, inform__, mqnmod, lvlhes;
    static logical prtmem;
    static integer lhetyp;
    static doublereal objtru;
    static logical fponly;
    static char usercw[8*130];
    static integer liwest;
    static char solver[6];
    static integer nextcw, errors, useriw[130], nextiw, lrwest;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal userrw[130];
    static integer nextrw, stkyop;

/*     ================================================================== */
/*     snKerB does the work for snOptB. (Kernel for snoptB) */

/*     Developers can call this version with customized versions of */
/*     snLog, snLog2  and  snSTOP. */

/*     17 Oct 2004: First version of snKerB. */
/*     01 Sep 2007: Sticky parameters added. */
/*     04 Jul 2010: mincw, miniw, minrw added to workspace. */
/*     21 Aug 2013: Ad hoc allocation of x0 made safe. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Problem name */
/* cold:warm:basis:hot start */
/*     ------------------------------------------------------------------ */
/* Current Hessian type */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --bu;
    --bl;
    --locj;
    --indj;
    --jcol;
    names -= 8;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNOPTB", (ftnlen)6, (ftnlen)6);
    *info = 0;
/*     ------------------------------------------------------------------ */
/*     Check memory limits and fetch the workspace starting positions. */
/*     ------------------------------------------------------------------ */
    s2mem0_(info, solver, lencw, leniw, lenrw, &iw[1], mincw, miniw, minrw, &
	    maxcw, &maxiw, &maxrw, &nextcw, &nextiw, &nextrw, (ftnlen)6);
    if (*info > 0) {
	goto L999;
    }
/*     Save the user's option choices  (weird choices get overwritten). */
/*     Initialize timers and the standard input file. */
/* Exit without printing */
    chcopy_(&c__130, cw + 408, &c__1, usercw, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__130, &rw[51], &c__1, userrw, &c__1);
    s1time_(&c__0, &c__0, &iw[1], leniw, &rw[1], lenrw);
    s1file_(&c__2, &iw[1], leniw);
/*     Check the input arguments. */
    s3argb_(&inform__, start, m, n, ne, nname, ns, nncon, nnobj, nnjac, iobj, 
	    &indj[1], &locj[1], &bl[1], &bu[1], names + 8, &hs[1], &pi[1], &
	    iw[69], &errors, &iw[1], leniw, &rw[1], lenrw, start_len, (ftnlen)
	    8);
    if (inform__ > 0) {
	*info = inform__;
	goto L800;
    }
/*     ------------------------------------------------------------------ */
/*     The obligatory call to snInit has already set the defaults. */
/*     Check that the optional parameters have sensible values. */
/*     Print the options. */
/*     ------------------------------------------------------------------ */
    s_copy(cw + 408, prob, (ftnlen)8, (ftnlen)8);
    s8dflt_(m, n, nncon, nnjac, nnobj, cw + 8, lencw, &iw[1], leniw, &rw[1], 
	    lenrw, (ftnlen)8);
    s3prtb_(m, n, nncon, nnjac, nnobj, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Compute the storage requirements for SNOPT  from the following */
/*     variables: */
/*         m,      n,     ne */
/*         lenR  , maxS */
/*         nnCon , nnJac, nnObj, */
/*         negCon */
/*     All have to be known before calling s2Mem. */
/*     The only one in doubt is negCon, the number of Jacobian elements. */
/*     Count them in s8Gsiz. */
/*     ------------------------------------------------------------------ */
    nb = *n + *m;
    nlocj = *n + 1;
    nkx = nb;
    s8gsiz_(m, nncon, nnjac, ne, &nlocj, &locj[1], &indj[1], &negcon);
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    lvlhes = iw[72];
/* 0, 1,  2  => LM, FM, Exact Hessian */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    nnh = max(*nnjac,*nnobj);
    ngobj = *nnobj;
/*     Load iw with various problem dimensions. */
/* Local nnObj is altered for FP */
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = *m;
/* copy of the number of rows */
    iw[17] = *ne;
/* copy of the number of nonzeros in Jcol */
    iw[20] = negcon;
/* # of nonzeros in gCon */
    iw[21] = *nnjac;
/* # of Jacobian  variables */
    iw[22] = *nnobj;
/* # of objective variables */
    iw[23] = *nncon;
/* # of nonlinear constraints */
    iw[24] = nnh;
/*   max( nnObj, nnJac ) */
    iw[28] = lenr;
/* R(lenR) is the reduced Hessian factor */
    iw[204] = *iobj;
/*     ------------------------------------------------------------------ */
/*     If only a feasible point is requested, save the base point for the */
/*     objective function:  1/2 || x - x0 ||^2 */

/*     s8Map does not allocate space for x0, so it is done here. */
/*     ------------------------------------------------------------------ */
/* position of the objective row in J */
    fponly = minmax == 0;
    if (fponly) {
	ngobj = nnh;
	lx0 = nextrw;
	nextrw = lx0 + nnh;
	iw[298] = lx0;
/* ad hoc allocated FP base point */
    }
/*     ------------------------------------------------------------------ */
/*     Allocate the local arrays for snOpt. */
/*     s8Map  maps snOpt integer and double arrays. */
/*     s2BMap maps the arrays for the LU routines. */
/*     s2Mem  checks what space is available and prints any messages. */
/*     ------------------------------------------------------------------ */
    s8map_(m, n, &negcon, &nkx, nncon, nnjac, &ngobj, &lenr, &maxr, &maxs, &
	    mqnmod, &lvlhes, &nextcw, &nextiw, &nextrw, &iw[1], leniw);
    s2bmap_(m, n, ne, &maxs, &nextiw, &nextrw, &maxiw, &maxrw, &liwest, &
	    lrwest, &iw[1], leniw);
    prtmem = TRUE_;
/* Print all messages in s2Mem */
    s2mem_(&inform__, &prtmem, &liwest, &lrwest, &nextcw, &nextiw, &nextrw, &
	    maxcw, &maxiw, &maxrw, lencw, leniw, lenrw, mincw, miniw, minrw, &
	    iw[1]);
    if (inform__ != 0) {
	*info = inform__;
	goto L800;
    }
/*     Its now safe to copy x to x0 */
    if (fponly) {
	dcopy_(&nnh, &x[1], &c__1, &rw[lx0], &c__1);
    }
/*     Define the row and column ordering for J. */
/*     SNOPT  uses natural order throughout, so kx = kxN. */
    iw[247] = nkx;
/* dimension of kx and its inverse, kxN */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    iw[252] = lkx;
/* jN = kxN(j ) => col j of Jcol is variable jN */
    s1perm_(n, &iw[lkx]);
    s1perm_(m, &iw[lkx + *n]);
/*     ------------------------------------------------------------------ */
/*     Construct column pointers for the nonlinear part of the  Jacobian. */
/*     ------------------------------------------------------------------ */
    if (*nncon > 0) {
	llocg = iw[260];
/* locG(nlocG) = column pointers for indG */
	nlocg = *nnjac + 1;
	s8gloc_(nncon, nnjac, ne, &nlocj, &locj[1], &indj[1], &negcon, &nlocg,
		 &iw[llocg]);
    }
/*     ------------------------------------------------------------------ */
/*     Solve the problem. */
/*     ------------------------------------------------------------------ */
    if (nnh == 0) {
/*        The problem is a linear program. */
	nrhs = 0;
/* No constraint rhs vector. */
	nx0 = 0;
/* No constant shift for x. */
	nrhs0 = max(nrhs,1);
	lenx0 = max(nx0,1);
	ngobj0 = max(ngobj,1);
	nnh0 = max(nnh,1);
	ngqp = max(ngobj,nnh);
	lhetyp = iw[283];
/* hEtype(nb) definition of elastic vars */
	lgobj = iw[297];
/* gObj(ngObj) = Objective gradient */
	iload_(&nb, &c__3, &iw[lhetyp], &c__1);
	s5solv_(info, solver, &iw[69], (U_fp)sqhx_, (U_fp)s8qphx_, (U_fp)
		sqlog, &gotr, m, n, &nb, &nnh0, &nnh, nname, &ngqp, &ngobj0, &
		ngobj, iobj, objadd, &fobj, &objtru, ninf, sinf, ne, &nlocj, &
		locj[1], &indj[1], &jcol[1], &bl[1], &bu[1], &rw[lgobj], 
		names + 8, &nrhs0, &nrhs, rhs, &lenx0, &nx0, x0, &iw[lhetyp], 
		&hs[1], &x[1], &pi[1], &rc[1], ns, cu + 8, lencu, &iu[1], 
		leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		lenrw, (ftnlen)6, (ftnlen)8, (ftnlen)8, (ftnlen)8);
    } else {
/*        The problem is nonlinear. */
/*        Define the type of initial Hessian. */
	if (iw[69] == 0) {
	    iw[202] = -1;
	} else if (iw[69] == 1) {
	    iw[202] = -1;
	} else if (iw[69] == 2) {
	    iw[202] = -1;
	} else if (iw[69] == 3) {
	    iw[202] = 0;
	}
	s8solv_(info, solver, &iw[69], (U_fp)s0fgb_, (U_fp)funcon, (U_fp)
		funobj, (U_fp)snlog, (U_fp)snlog2, (U_fp)snstop, &gotr, m, n, 
		&nb, nncon, nnjac, &ngobj, nname, iobj, objadd, &fobj, &
		objtru, ninf, sinf, ne, &nlocj, &locj[1], &indj[1], &jcol[1], 
		&bl[1], &bu[1], names + 8, &hs[1], &x[1], &pi[1], &rc[1], &
		nmajor, ns, cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
		8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)6, (ftnlen)8, 
		(ftnlen)8, (ftnlen)8);
    }
    *obj = fobj;
    *mincw = iw[47];
/* minimum length of cw */
    *miniw = iw[48];
/* minimum length of iw */
    *minrw = iw[49];
/*     If "sticky parameters no",  restore the user-defined options */
/* minimum length of rw */
    stkyop = iw[116];
    if (stkyop <= 0) {
	chcopy_(&c__130, usercw, &c__1, cw + 408, &c__1, (ftnlen)8, (ftnlen)8)
		;
	icopy_(&c__130, useriw, &c__1, &iw[51], &c__1);
	dcopy_(&c__130, userrw, &c__1, &rw[51], &c__1);
    }
/*     Print times for all clocks (if lvlTim > 0). */
    s1time_(&c__0, &c__2, &iw[1], leniw, &rw[1], lenrw);
    return 0;
/*     Local exit messages. */
L800:
    snwrap_(info, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* snkerb_ */

