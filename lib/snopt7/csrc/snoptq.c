/* ./src/snoptq.f -- translated by f2c (version 20090411).
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

static integer c__109 = 109;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__130 = 130;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  snoptq.f */

/*     snOptQ   snKerQ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snoptq_(char *start, U_fp qphx, integer *m, integer *n, 
	integer *ne, integer *nname, integer *lencobj, integer *ncolh, 
	integer *iobj, doublereal *objadd, char *prob, doublereal *acol, 
	integer *inda, integer *loca, doublereal *bl, doublereal *bu, 
	doublereal *cobj, char *names, integer *hetype, integer *hs, 
	doublereal *x, doublereal *pi, doublereal *rc, integer *info, integer 
	*mincw, integer *miniw, integer *minrw, integer *ns, integer *ninf, 
	doublereal *sinf, doublereal *obj, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen start_len, ftnlen prob_len, ftnlen names_len, ftnlen cu_len, 
	ftnlen cw_len)
{
    extern /* Subroutine */ int sqlog_();
    extern /* Subroutine */ int snkerq_(char *, U_fp, U_fp, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, char *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);

/*     ------------------------------------------------------------------ */
/*     snOptQ    solves the linear/quadratic programming problem */

/*     Min/max     ObjAdd + <a(iObj),x> + <cObj,x> + half <x,H*x> */
/*        x */
/*     subject to linear constraints, and upper and lower bounds on x */
/*     (see below). */

/*     ObjAdd  is a constant scalar. */
/*     a(iObj) is the iObj-th row of the constraint matrix A (see below). */
/*     cObj    is a constant vector. */
/*     H       is a constant symmetric matrix, defined implicitly in the */
/*             user-supplied subroutine qpHx. */
/*             o Subroutine qpHx must evaluate products of H with a given */
/*               vector x. This implies that H is defined as an operator, */
/*               and need never be defined explicitly. */
/*             o If H has zero rows and columns, it is most efficient to */
/*               order the variables so that */
/*                      Hx = ( H1  0 )( X1 ) = ( H1 x1 ), */
/*                           ( 0   0 )( X2)    (   0   ) */
/*               where the nonlinear variables x1 appear first as shown. */

/*     The constraints take the form */

/*                            (  x ) */
/*                      bl <= (    ) <= bu, */
/*                            ( Ax ) */

/*     where bl and bu are constant lower and upper bounds. */

/*     Internally, the constraints are rewritten as */

/*                                         ( x ) */
/*                   Ax - s = 0,     bl <= (   ) <= bu, */
/*                                         ( s ) */

/*     where s is an m-vector of slack variables. */
/*     components of (x,s) are the variables and slacks respectively. */
/*     The sparse matrix A is entered column-wise in the arrays */
/*     Acol, indA, locA (below). */

/*     ------------------------------------------------------------------ */
/*     NOTE: Before calling snOptQ, the calling program must call: */
/*     call sqInit( iPrint, iSumm, */
/*    &             cw, lencw, iw, leniw, rw, lenrw ) */
/*     This sets the default values of the optional parameters. You can */
/*     also alter the default values of iPrint and iSumm before snOptQ */
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

/*     m       is the number of rows in the constraint matrix A. */
/*             m > 0. */

/*     n       is the number of variables, excluding slacks. */
/*             For LP problems, this is the number of columns in A. */
/*             n > 0. */

/*     ne      is the number of nonzero entries in A. */
/*             ne > 0. */

/*     nName   is the number of column and row names provided in the */
/*             array  Names.  If nName = 1, there are NO names. */
/*             Generic names will be used in the printed solution. */
/*             Otherwise, nName = n+m and all names must be provided. */

/*     lencObj is the number of elements in the constant objective cObj. */
/*             lencObj ge 0. */

/*     ncolH   is the number of leading nonzero columns of the */
/*             QP Hessian. */
/*             If ncolH > 0, the subroutine qpHx must provide the */
/*             matrix-product Hx. */
/*             ncolH ge 0. */

/*     iObj    says which row of A is a free row containing a linear */
/*             objective vector  cObj  (iObj = 0 if none). */
/*             iObj = 0  or  iObj le m. */

/*     ObjAdd  is a constant that will be added to the objective before */
/*             printing.  Typically,  ObjAdd = 0.0d+0. */

/*     Prob    is an 8-character name for the problem, used in the */
/*             output.  A blank name can be assigned if necessary. */

/*     Acol(ne) is the constraint matrix, stored column-wise. */

/*     indA(ne) is the list of row indices for each nonzero in a(*). */

/*     locA(n+1) is a set of pointers to the beginning of each column of */
/*             the constraint matrix within Acol(*) and indA(*). */
/*             MUST HAVE locA(1) = 1 AND locA(n+1) = ne+1. */

/*  NOTES:  1. If lencObj > 0, the first lencObj columns of Acol and indA belong */
/*             to variables corresponding to the constant */
/*             objective term c. */

/*          2. If the problem has a quadratic objective, */
/*             The first ncolH columns of Acol and indA belong to variables */
/*             corresponding to the nonzero block of the QP Hessian. */
/*             Subroutine qpHx deals with these variables. */

/*          3. If lencObj > 0 and ncolH > 0, the two sets of */
/*             objective variables overlap.  The total number of */
/*             objective variables is nQP = max( lencObj, ncolH ). */

/*          4. The row indices indA(k) for a column may be in any order. */

/*     bl(n+m) is the lower bounds on the variables and slacks (x, s). */

/*     bu(n+m) is the upper bounds on (x, s). */

/*     Names(nName) is an character*8 array. */
/*             If nName =  1, Names is not used.  The printed solution */
/*             will use generic names for the columns and row. */
/*             If nName = n+m, Names(j) should contain an 8 character */
/*             name of the jth variable (j = 1, n+m). */
/*             If j = n+i, the jth variable is the ith row. */

/*     hEtype(n+m) indicate if the variable can violate its bound in */
/*             Elastic mode. */
/*             if hEtype(j) = 0, variable j is non-elastic and must not */
/*                               be violated. */
/*             if hEtype(j) = 1, variable j can violate its lower bound. */
/*             if hEtype(j) = 2, variable j can violate its upper bound. */
/*             if hEtype(j) = 3, variable j can violate either its */
/*                                          lower or upper bound. */

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

/*          4. For Warm starts, all of hs(1:n+m) is assumed to be */
/*             set to the values 0, 1, 2 or 3 from some previous call. */

/*     nS      need not be specified for Cold starts, */
/*             but should retain its value from a previous call */
/*             when a Warm start is used. */


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

/*     x(n+m)  is the final variables and slacks (x, s). */

/*     pi(m)   is the vector of Lagrange multipliers (shadow prices) */
/*             for the general constraints. */

/*     rc(n+m) is a vector of reduced costs: rc = g - (A -I)'*pi, where g */
/*             is the gradient of the objective if x is feasible */
/*             (or the gradient of the Phase 1 objective otherwise). */
/*             If nInf = 0, the last m entries are pi. */

/*     INFO    says what happened; see the User's Guide. */
/*             The possible values are as follows: */

/*             INFO    Meaning */

/*                0    finished successfully */
/*                1       optimality conditions satisfied */
/*                2       feasible point found */
/*                4       weak QP minimizer */

/*               10    the problem appears to be infeasible */
/*               11       infeasible linear constraints */
/*               12       infeasible linear equalities */
/*               14       infeasibilities minimized */

/*               20    the problem appears to be unbounded */
/*               21       unbounded objective */

/*               30    resource limit error */
/*               31       iteration limit reached */
/*               33       the superbasics limit is too small */

/*               40    terminated after numerical difficulties */
/*               42       singular basis */
/*               43       cannot satisfy the general constraints */
/*               44       ill-conditioned null-space basis */

/*               50    error in the user-supplied functions */
/*               53       the QP Hessian is indefinite */

/*               70    user requested termination */
/*               73       terminated during QP objective evaluation */
/*               74       terminated from monitor routine */

/*               80    insufficient storage allocated */
/*               81       work arrays must have at least 500 elements */
/*               82       not enough character storage */
/*               83       not enough integer storage */
/*               84       not enough real storage */

/*               90    input arguments out of range */
/*               91       invalid input argument */
/*               92       basis file dimensions do not match this problem */
/*               93       the QP Hessian is indefinite */

/*              140    system error */
/*              141       wrong no of basic variables */
/*              142       error in basis package */

/*     mincw   says how much character storage is needed to solve the */
/*             problem.  If INFO = 82, the work array cw(lencw) was */
/*             too small.  sqOpt may be called again with lencw suitably */
/*             larger than mincw. */

/*     miniw   says how much integer storage is needed to solve the */
/*             problem.  If INFO = 83, the work array iw(leniw) was too */
/*             small.  sqOpt may be called again with leniw suitably */
/*             larger than miniw.  (The bigger the better, since it is */
/*             not certain how much storage the basis factors need.) */

/*     minrw   says how much real storage is needed to solve the */
/*             problem.  If INFO = 84, the work array rw(lenrw) was too */
/*             small.  (See the comments above for miniw.) */

/*     nS      is the final number of superbasics. */

/*     nInf    is the number of infeasibilities. */

/*     sInf    is the sum    of infeasibilities. */

/*     Obj     is the value of the QP objective function. */
/*             Obj does NOT include ObjAdd or the objective row. */
/*             If nInf = 0, Obj includes the quadratic objective if any. */
/*             If nInf > 0, Obj is just the linear objective if any. */

/*     cu(lencu), iu(leniu), ru(lenru)  are character, integer and real */
/*             arrays of USER workspace.  These arrays are available to */
/*             pass data to the user-defined routine qpHx. */
/*             If no workspace is required, you can either use dummy */
/*             arrays for cu, iu and ru, or use cw, iw and rw */
/*             (see below). */

/*     cw(lencw), iw(leniw), rw(lenrw)  are character*8, integer and real */
/*             arrays of workspace used by SNOPTQ. */
/*             lencw  should be about at least 500. */
/*             leniw  should be about max( 500, 10(m+n) ) or larger. */
/*             lenrw  should be about max( 500, 20(m+n) ) or larger. */

/*     SNOPTQ is maintained by Philip E. Gill, */
/*     Dept of Mathematics, University of California, San Diego. */

/*     LUSOL is maintained by Michael A. Saunders, */
/*     Systems Optimization Laboratory, */
/*     Dept of Management Science & Engineering, Stanford University. */

/*     12 Nov 1994: Workspace separated into iw(*) and rw(*). */
/*     20 Jul 1996: Slacks changed to be the row value. */
/*     09 Aug 1996: First Min Sum version. */
/*     17 Jul 1997: Thread-safe version. */
/*     26 Jul 1997: User workspace added. */
/*     02 Oct 1997: Character workspace added. */
/*     29 Jul 2003: Variables marked correctly for final printing. */
/*     15 Oct 2003: snPRNT added. */
/*     01 Sep 2007: Sticky parameters added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --hetype;
    --bu;
    --bl;
    --loca;
    --inda;
    --acol;
    names -= 8;
    --cobj;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    snkerq_(start, (U_fp)qphx, (U_fp)sqlog_, m, n, ne, nname, lencobj, ncolh, 
	    iobj, objadd, prob, &acol[1], &inda[1], &loca[1], &bl[1], &bu[1], 
	    &cobj[1], names + 8, &hetype[1], &hs[1], &x[1], &pi[1], &rc[1], 
	    info, mincw, miniw, minrw, ns, ninf, sinf, obj, cu + 8, lencu, &
	    iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1],
	     lenrw, start_len, (ftnlen)8, (ftnlen)8, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snoptq_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snOptQ */
/* Subroutine */ int snkerq_(char *start, U_fp qphx, U_fp sqlog, integer *m, 
	integer *n, integer *ne, integer *nname, integer *lencobj, integer *
	ncolh, integer *iobj, doublereal *objadd, char *prob, doublereal *
	acol, integer *inda, integer *loca, doublereal *bl, doublereal *bu, 
	doublereal *cobj, char *names, integer *hetype, integer *hs, 
	doublereal *x, doublereal *pi, doublereal *rc, integer *info, integer 
	*mincw, integer *miniw, integer *minrw, integer *ns, integer *ninf, 
	doublereal *sinf, doublereal *obj, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen start_len, ftnlen prob_len, ftnlen names_len, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal x0[1];
    static integer nb, nx0, nnh;
    static doublereal rhs[1];
    static integer lkx, nkx;
    static char str[80];
    static integer nnh0;
    static char str2[80];
    static integer lenr, ngqp, maxr, maxs, nrhs;
    static logical gotr;
    extern /* Subroutine */ int sqhx_();
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s5map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer lenx0, nrhs0, nnjac, nloca, ngobj;
    extern /* Subroutine */ int s2mem0_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, ftnlen);
    static integer nncon, maxcw;
    static doublereal objqp;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxiw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxrw, ngobj0;
    extern /* Subroutine */ int s1file_(integer *, integer *, integer *), 
	    s2bmap_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), s3argq_(integer *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, char *, integer 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), s5dflt_(integer *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen), s1time_(integer *, 
	    integer *, integer *, integer *, doublereal *, integer *), 
	    s1perm_(integer *, integer *), s5solv_(integer *, char *, integer 
	    *, U_fp, U_fp, U_fp, logical *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), s3prtq_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *), chcopy_(integer *, char *, integer *, 
	    char *, integer *, ftnlen, ftnlen);
    static integer inform__;
    static logical prtmem;
    static doublereal objtru;
    static integer liwest;
    static char usercw[8*130];
    static integer nextcw;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer errors, nextiw, useriw[130], lrwest, nextrw;
    static doublereal userrw[130];
    static integer stkyop;

/*     ================================================================== */
/*     snKerQ does the work for snOptQ. (i.e., its the kernel for sqopt) */

/*     Developers can call snKerQ with customized versions of sqLog. */

/*     04 Jan 2005: First version of snKerQ. */
/*     01 Sep 2007: Sticky parameters added. */
/*     04 Jul 2010: mincw, miniw, minrw added to workspace. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Problem name */
/* cold:warm:basis:hot start */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --hetype;
    --bu;
    --bl;
    --loca;
    --inda;
    --acol;
    names -= 8;
    --cobj;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SQOPT ", (ftnlen)6, (ftnlen)6);
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
    chcopy_(&c__109, cw + 408, &c__1, usercw, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__109, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__109, &rw[51], &c__1, userrw, &c__1);
    s1time_(&c__0, &c__0, &iw[1], leniw, &rw[1], lenrw);
    s1file_(&c__2, &iw[1], leniw);
/*     Check the input arguments. */
    s3argq_(&inform__, start, m, n, ne, nname, ns, lencobj, iobj, ncolh, &
	    inda[1], &loca[1], &bl[1], &bu[1], names + 8, &hs[1], &pi[1], &iw[
	    69], &errors, &iw[1], leniw, &rw[1], lenrw, start_len, (ftnlen)8);
    if (inform__ > 0) {
	*info = inform__;
	goto L800;
    }
/*     ------------------------------------------------------------------ */
/*     The obligatory call to sqInit has already set the defaults. */
/*     Check that the optional parameters have sensible values. */
/*     Print the options if requested. */
/*     ------------------------------------------------------------------ */
    s_copy(cw + 408, prob, (ftnlen)8, (ftnlen)8);
    s5dflt_(m, n, lencobj, ncolh, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
	     (ftnlen)8);
    s3prtq_(m, n, lencobj, ncolh, &iw[1], leniw, &rw[1], lenrw);
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    nnh = *ncolh;
    nnh0 = max(nnh,1);
    ngobj = *lencobj;
    ngobj0 = max(ngobj,1);
    ngqp = max(ngobj,nnh);
    nb = *n + *m;
    nkx = nb;
    nncon = 0;
    nnjac = 0;
/*     Load iw with various problem dimensions. */
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = *m;
/* copy of the number of rows */
    iw[17] = *ne;
/* copy of the number of nonzeros in Jcol */
    iw[21] = nnjac;
/* # of Jacobian  variables */
    iw[22] = nnh;
/* # of objective variables */
    iw[23] = nncon;
/* # of nonlinear constraints */
    iw[24] = nnh;
/*   max( nnObj, nnJac ) */
    iw[26] = ngobj;
/* length of QP constant vector */
    iw[27] = nnh;
/* # QP Hessian columns */
    iw[28] = lenr;
/* R(lenR) is the reduced Hessian factor */
    iw[204] = *iobj;
/*     ------------------------------------------------------------------ */
/*     Allocate the local arrays for snOptQ. */
/*     s5Map  maps snOptQ integer and double arrays. */
/*     s2BMap maps the arrays for the LU routines. */
/*     s2Mem  checks what space is available and prints any messages. */
/*     ------------------------------------------------------------------ */
/* position of the objective row in J */
    s5map_(m, n, &nkx, &ngobj, &nnh, &lenr, &maxr, &maxs, &nextcw, &nextiw, &
	    nextrw, &iw[1], leniw);
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
/*     Define the row and column ordering for Acol. */
/*     SNOPTQ  uses natural order throughout, so kx = kxN. */
    iw[247] = nkx;
/* dimension of kx and its inverse, kxN */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    iw[252] = lkx;
/* jN = kxN(j ) => col j of Jcol is variable jN */
    s1perm_(n, &iw[lkx]);
    s1perm_(m, &iw[lkx + *n]);
/*     ------------------------------------------------------------------ */
/*     Set values for the unused features in snOptQ. */
/*     ------------------------------------------------------------------ */
    nrhs = 0;
/* No constraint rhs vector. */
    nrhs0 = 1;
    nx0 = 0;
/* No constant shift for x. */
    lenx0 = 1;
/*     ------------------------------------------------------------------ */
/*     Solve the problem. */
/*     ------------------------------------------------------------------ */
    nloca = *n + 1;
    s5solv_(info, solver, &iw[69], (U_fp)sqhx_, (U_fp)qphx, (U_fp)sqlog, &
	    gotr, m, n, &nb, &nnh0, &nnh, nname, &ngqp, &ngobj0, &ngobj, iobj,
	     objadd, &objqp, &objtru, ninf, sinf, ne, &nloca, &loca[1], &inda[
	    1], &acol[1], &bl[1], &bu[1], &cobj[1], names + 8, &nrhs0, &nrhs, 
	    rhs, &lenx0, &nx0, x0, &hetype[1], &hs[1], &x[1], &pi[1], &rc[1], 
	    ns, cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &
	    iw[1], leniw, &rw[1], lenrw, (ftnlen)6, (ftnlen)8, (ftnlen)8, (
	    ftnlen)8);
    *obj = objqp;
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
} /* snkerq_ */

