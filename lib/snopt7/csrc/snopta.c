/* ./src/snopta.f -- translated by f2c (version 20100827).
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
static doublereal c_b20 = 0.;
static integer c__3 = 3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  snopta.f --- the free format interface for snOpt */

/*     snOptA    snKerA */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snopta_(integer *start, integer *nf, integer *n, integer 
	*nxname, integer *nfname, doublereal *objadd, integer *objrow, char *
	prob, U_fp usrfun, integer *iafun, integer *javar, integer *lena, 
	integer *nea, doublereal *a, integer *igfun, integer *jgvar, integer *
	leng, integer *neg, doublereal *xlow, doublereal *xupp, char *xnames, 
	doublereal *flow, doublereal *fupp, char *fnames, doublereal *x, 
	integer *xstate, doublereal *xmul, doublereal *f, integer *fstate, 
	doublereal *fmul, integer *info, integer *mincw, integer *miniw, 
	integer *minrw, integer *ns, integer *ninf, doublereal *sinf, char *
	cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, 
	integer *lenru, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen prob_len, ftnlen xnames_len, 
	ftnlen fnames_len, ftnlen cu_len, ftnlen cw_len)
{
    extern /* Subroutine */ int snlog_(), sqlog_(), snlog2_();
    extern /* Subroutine */ int snkera_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, U_fp, U_fp,
	     U_fp, U_fp, U_fp, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, char *, doublereal *, doublereal *, 
	    char *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, char *, integer *, integer 
	    *, integer *, doublereal *, integer *, char *, integer *, integer 
	    *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen);
    extern /* Subroutine */ int snstop_();

/*     ================================================================== */
/*     snOptA  is a Fortran wrappper for the SNOPT solver. */
/*     snOptA   is a subroutine for constrained nonlinear */
/*     optimization.  The optimization problem involves m  functions */
/*     F(1), F(2), ... , F(nF), each a function of n variables */
/*     x(1), x(2), ... , x(n).  The  problem has the form: */

/*           minimize/maximize    ObjAdd + F(ObjRow) */

/*                            ( xlow <=  x  <= xupp, */
/*                 subject to ( */
/*                            ( Flow <=  F  <= Fupp, */

/*     where ObjAdd is a constant, ObjRow is a user-specified row of  F, */
/*     xlow, Flow, xupp and Fupp are constant lower and upper bounds. */

/*     ------------------------------------------------------------------ */
/*     NOTE: Before calling SNOPTA, your calling program MUST call the */
/*     initialization routine using the call: */
/*     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw ) */
/*     This sets the default values of the optional parameters. You can */
/*     also alter the default values of iPrint and iSumm before snOptA */
/*     is used.  iPrint = 0, etc, is OK. */
/*     ------------------------------------------------------------------ */

/*     o If ObjRow = 0, then snOptA will find a point satisfying the */
/*       constraints. */

/*     o If all functions are linear, F = A x for some sparse matrix A. */
/*       This defines a linear program (LP).  In this case,  the nonzero */
/*       elements of A can be input in coordinate form (i,j,A_ij) (see */
/*       below). */

/*     o If all functions are nonlinear, F = F(x) for some vector */
/*       F(x) of smooth functions.  In this case, the elements of  F  and */
/*       (optionally) their first and second partial derivatives must be */
/*       coded by the user in the subroutine usrfun  (see below). */

/*     o If some functions are linear and some are nonlinear, the user */
/*       can choose to set every component in usrfun.  It is usually more */
/*       efficient, however,  to supply the coefficients of the linear */
/*       functions via the sparse array  A (see below).   In this case, */
/*       the linear elements of  F  need not be assigned (SNOPTA will */
/*       figure out which components of  F  are needed). */

/*     o In the most general situation, the ith component of F(x) is the */
/*       sum of linear and nonlinear terms.  In this case, if F(x) can be */
/*       defined as a sum of "non-overlapping" linear and nonlinear */
/*       functions, then the nonlinear part of F can be defined in usrfun */
/*       and the linear part can be defined via the array A. */

/*       Suppose that the ith component of F(x) is of the form */
/*            F_i(x) = f_i(x) + sum (over j)  A_ij x_j, */
/*       where f_i(x) is a nonlinear function and the elements A_ij */
/*       are constant.   It is convenient to write  F_i(x)  in the */
/*       compact form  F_i(x) = f_i(x) + A_i' x, where A_i denotes a */
/*       column vector with components ( A_i1, A_i2, ..., A_in ), and */
/*       "A_i'" denotes the transpose of A_i. */

/*       Functions f_i and A_i are said to be "non-overlapping" if any */
/*       variable x_j  appearing explicitly in f_i(x) does not appear */
/*       explicitly in A_i'x, i.e., A_ij = 0.  (Equivalently, any */
/*       variable with a nonzero A_ij must not appear explicitly in */
/*       f_i(x).)  For example, the function */
/*         F_i(x) = 3x_1 + exp(x_2)x_4 + x_2^2 + 4x_4 - x_3 + x_5 */
/*       can be written as the sum of non-overlapping functions f_i and */
/*       A_i'x, such that */
/*           f_i(x) = exp(x_2)x_4 + x_2^2  + 4x_4  and */
/*           A_i'x  = 3x_1 - x_3 + x_5. */

/*       Given a non-overlapping sum for each component of F, we can */
/*       write  F(x) = f(x) + Ax, where f(x) is a vector-valued function */
/*       of x and A is a sparse matrix whose ith row is A_i'. */

/*       The nF by n  Jacobian of  F(x)  is the sum of two  nF by n */
/*       sparse matrices G and A,  i.e.,  J = G + A,  where G and A */
/*       contain the nonlinear and constant elements of J respectively. */
/*       The important property of non-overlapping functions is that */
/*       a nonzero entry of J is either an element of A, or an element */
/*       of G, but NOT BOTH (i.e., the nonzeros of  A  and  G  do not */
/*       overlap. */

/*       The nonzero elements of A and G must be provided in coordinate */
/*       form.  In coordinate form, a nonzero element G_ij of a matrix */
/*       G  is stored as the triple (i,j,G_ij).  The kth coordinate is */
/*       defined by iGfun(k) and jGvar(k)  (i.e., if i=iGfun(k) and */
/*       j=jGvar(k), then G(k) is the ijth element of G.)  Any known */
/*       values of G(k) must be assigned by the user in the routine */
/*       usrfun. */

/*       RESTRICTIONS: */
/*        1.  If the elements of G cannot be provided because they are */
/*            either too expensive or too complicated to evaluate,  it */
/*            is still necessary to specify the position of the nonzeros */
/*            as specified by the arrays iGfun and jGvar. */

/*        2.  If an element of G happens to be zero at a given point, */
/*            it must still be loaded in usrfun. (The order of the */
/*            list of coordinates (triples) is meaningful in snOptA.) */

/*       The elements of A and G can be stored in any order, (e.g., by */
/*       rows, by columns, or mixed). Duplicate entries are ignored. */

/*     ON ENTRY: */

/*     Start   specifies how a starting basis (and certain other items) */
/*             are to be obtained. */
/*             start =  0 (Cold) means that Crash should be used to */
/*                      choose an initial basis, unless a basis file is */
/*                      given by reference in the Specs file to an */
/*                      Old basis file. */
/*             start =  1 (Basis file) means the same (but is more */
/*                      meaningful in the latter case). */
/*             start =  2 (Warm) means that a basis is already defined */
/*                      in xstate and Fstate (probably from an earlier */
/*                      call). */

/*     nF      is the number  of problem functions in F, including the */
/*             objective function (if any) and the linear */
/*             and nonlinear constraints.  Simple upper and lower bound */
/*             constraints on the variables should not be included in  F. */
/*             nF > 0. */

/*     n       is the number of variables. */
/*             n > 0. */

/*     neA     is the number of nonzero entries in A. */
/*             neA >= 0. */

/*     nxname  is the number of 8-character column (i.e., variable) names */
/*             provided in the array xnames.  If nxname = 1,  then there */
/*             are NO column names (generic names will be used in the */
/*             printed solution).  Otherwise, nxname = n and every */
/*             column name must be provided. */

/*     nFname  is the number of 8-character row (i.e., constraint and */
/*             objective) names provided in the array Fnames. */
/*             If nFname = 1,  then there are NO row names (generic */
/*             names will be used in the printed solution).  Otherwise, */
/*             nFname = nF and every row name must be provided. */

/*     ObjAdd  is a constant that will be added to the objective. */
/*             Typically ObjAdd = 0.0d+0. */

/*     Prob    is an 8-character name for the problem, used in the */
/*             output.  A blank name can be assigned if necessary. */

/*     xlow(n) are the lower bounds on x. */

/*     xupp(n) are the upper bounds on x. */

/*     xnames(nxname) is an character*8 array of names for each x(j). */
/*             If nxname =  1, xnames is not used.  The printed solution */
/*             will use generic names for the variables. */
/*             If nxname = n, xnames(j) should contain an 8 character */
/*             name of the jth variable (j = 1:n). */

/*     Flow(n) are the lower bounds on F.  If component F(ObjRow) */
/*             is being optimized,  Flow(ObjRow) is ignored. */

/*     Fupp(n) are the upper bounds on F.  If component F(ObjRow) */
/*             is being optimized,  Fupp(ObjRow) is ignored. */

/*     Fnames(nFname) is an character*8 array of names for each F(i). */
/*             If nFname =  1, Fnames is not used.  The printed solution */
/*             will use generic names for the objective and constraints. */
/*             If nName = nF, Fnames(j) should contain an 8 character */
/*             name of the jth constraint (j=1:nF). */

/*     xstate(n) sometimes contains a set of initial states for each */
/*             variable x.  See the following NOTES. */

/*     x(n)    is a set of initial values for each variable  x. */

/*  NOTES:  1. If start = 0 (Cold) or 1 (Basis file) and an OLD BASIS */
/*             file is to be input, xstate and x need not be set at all. */

/*          2. Otherwise, xstate(1:n) must be defined for a cold start. */
/*             If nothing special is known about the problem, or if */
/*             there is no wish to provide special information, */
/*             you may set xstate(j) = 0, x(j) = 0.0d+0 for all j=1:n. */
/*             All variables will be eligible for the initial basis. */

/*             Less trivially, to say that variable j will probably */
/*             be equal to one of its bounds, */
/*             set xstate(j) = 4 and x(j) = bl(j) */
/*             or  xstate(j) = 5 and x(j) = bu(j) as appropriate. */

/*          3. For Cold starts with no basis file, a Crash procedure */
/*             is used to select an initial basis.  The initial basis */
/*             matrix will be triangular (ignoring certain small */
/*             entries in each column). */
/*             The values xstate(j) = 0, 1, 2, 3, 4, 5 have the following */
/*             meaning: */

/*             xstate(j)  State of variable j during Crash */
/*             ---------  -------------------------------- */
/*             0, 1, 3    Eligible for the basis.  3 is given preference. */
/*             2, 4, 5    Ignored. */

/*             After Crash, xstate(j) = 2 entries are made superbasic. */
/*             Other entries not selected for the basis are made */
/*             nonbasic at the value x(j) if bl(j) <= x(j) <= bu(j), */
/*             or at the value bl(j) or bu(j) closest to x(j). */

/*          4. For Warm starts, all of Fstate(1:nF) is assumed to be */
/*             set to the values 0, 1, 2 or 3 from some previous call. */

/*     Fmul(nF) contains an estimate of the Lagrange multipliers */
/*             (shadow prices) for the F- constraints.  They are used */
/*             to define the Lagrangian for the first major iteration. */
/*             If nothing is known about Fmul, set */
/*             Fmul(i) = 0.0d+0, i = 1:nF */

/*     ON EXIT: */

/*     xstate(n) is the final state vector for x: */

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

/*     x(n)    contains the final variables. */

/*     F(nF)   contains the final values of F. */

/*     xmul(nF) is the vector of Lagrange multipliers (shadow prices) */
/*             for the variables constraints. */

/*     Fmul(nF) is the vector of Lagrange multipliers (shadow prices) */
/*             for the general constraints. */

/*     INFO    says what happened; see the User's Guide. */
/*             The possible values are as follows: */

/*             INFO       Meaning */

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
/*               93       the QP Hessian is indefinite */

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

/*     cu(lencu), iu(leniu), ru(lenru)  are character, integer and real */
/*             arrays of USER workspace.  These arrays are available to */
/*             pass data to the user-defined routine usrfun. */
/*             If no workspace is required, you can either use dummy */
/*             arrays for cu, iu and ru, or use cw, iw and rw */
/*             (see below). */

/*     cw(lencw), iw(leniw), rw(lenrw)  are character*8, integer and real */
/*             arrays of workspace used by snOptA. */
/*             lencw should be at least 500, or nF+n if names are given. */
/*                              +. */
/*             leniw should be about max( 500, 20(nF+n) ) or larger. */
/*             lenrw should be about max( 500, 40(nF+n) ) or larger. */

/*     SNOPT package maintained by Philip E. Gill, */
/*     Dept of Mathematics, University of California, San Diego. */

/*     08 Nov 1998: First version based on the snopt of SNOPT 5.3-4. */
/*     25 Aug 1999: for SNOPT Version 6.0. */
/*     04 Nov 2001: LP's solved explicitly. */
/*     31 Jul 2003: snEXIT and snPRNT adopted. */
/*     02 May 2004: Call to base routine added. */
/*     01 Sep 2007: Sticky parameters added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fmul;
    --fstate;
    --f;
    --fupp;
    --flow;
    --xmul;
    --xstate;
    --x;
    --xupp;
    --xlow;
    xnames -= 8;
    fnames -= 8;
    --a;
    --javar;
    --iafun;
    --jgvar;
    --igfun;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    snkera_(start, nf, n, nxname, nfname, objadd, objrow, prob, (U_fp)usrfun, 
	    (U_fp)snlog_, (U_fp)snlog2_, (U_fp)sqlog_, (U_fp)snstop_, &iafun[
	    1], &javar[1], lena, nea, &a[1], &igfun[1], &jgvar[1], leng, neg, 
	    &xlow[1], &xupp[1], xnames + 8, &flow[1], &fupp[1], fnames + 8, &
	    x[1], &xstate[1], &xmul[1], &f[1], &fstate[1], &fmul[1], info, 
	    mincw, miniw, minrw, ns, ninf, sinf, cu + 8, lencu, &iu[1], leniu,
	     &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)8, (ftnlen)8, (ftnlen)8, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snopta_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snOptA */
/* Subroutine */ int snkera_(integer *start, integer *nf, integer *n, integer 
	*nxname, integer *nfname, doublereal *objadd, integer *objrow, char *
	prob, U_fp usrfun, U_fp snlog, U_fp snlog2, U_fp sqlog, U_fp snstop, 
	integer *iafun, integer *javar, integer *lena, integer *nea, 
	doublereal *a, integer *igfun, integer *jgvar, integer *leng, integer 
	*neg, doublereal *xlow, doublereal *xupp, char *xnames, doublereal *
	flow, doublereal *fupp, char *fnames, doublereal *x, integer *xstate, 
	doublereal *xmul, doublereal *f, integer *fstate, doublereal *fmul, 
	integer *info, integer *mincw, integer *miniw, integer *minrw, 
	integer *ns, integer *ninf, doublereal *sinf, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen prob_len, ftnlen xnames_len, ftnlen 
	fnames_len, ftnlen cu_len, ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer m;
    static doublereal x0[1];
    static integer mqnmodtmp, lvlhestmp, nb, ne, lx, lx0, nx0, lbl, lrc, lbu, 
	    nnh, lpi, lhs, lfx;
    static doublereal rhs[1];
    static integer liy, lkx, nkx;
    static char str[80];
    static integer nnh0;
    static char str2[80];
    static doublereal fobj;
    static integer iobj, lenr, ngqp, maxr, maxs;
    static logical gotr;
    static integer lkxn, nrhs;
    extern /* Subroutine */ int sqhx_(), s0fga_();
    extern /* Subroutine */ int s3ina_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, ftnlen, ftnlen, ftnlen), s2mem_(
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer lenx0;
    extern /* Subroutine */ int s8map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer nrhs0;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *);
    static integer nnjac, lgobj, lindg, nname, lindj, llocg, ljcol, llocj, 
	    nnobj, nncon, maxcw, nlocj, nlocg, maxiw;
    extern /* Subroutine */ int s2mem0_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, ftnlen);
    static integer maxrw;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s3blda_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *);
    static integer lgobj1, lgobj2;
    extern /* Subroutine */ int s3arga_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, char *, doublereal *, 
	    doublereal *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), s1file_(integer *, 
	    integer *, integer *), s3mapa_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), s2bmap_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static integer nnobj0;
    extern /* Subroutine */ int s1time_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *), s8dflt_(integer *, integer *,
	     integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen), s3prma_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), s3siza_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), s3prta_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *), s3outa_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
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
    static doublereal infbnd;
    static integer minbld, negcon, lnames, objspc, ligfun, ljgvar, lnglin, 
	    nmajor, inform__, mqnmod, lvlhes, objsav;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static logical prtmem;
    static integer lhetyp;
    static doublereal objtru;
    static integer liwest;
    static char usercw[8*130], solver[6];
    static integer nextcw, errors, useriw[130], nextiw, lrwest;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal userrw[130];
    static integer nextrw, stkyop;
    extern /* Subroutine */ int s3dflta_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer lenrtmp, maxrtmp, maxstmp;

/*     ================================================================== */
/*     snKerA does the work for snOptA. (Kernel for snoptA) */

/*     Developers can call this version with customized versions of */
/*     snLog, snLog2  and  snSTOP. */

/*     17 Oct 2004: First version of snKerA. */
/*     01 Sep 2007: Sticky parameters added. */
/*     24 Nov 2007: s3dfltA sets options dependent on problem dimensions. */
/*     04 Jul 2010: mincw, miniw, minrw added to workspace. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Problem name */
/* cold:warm:basis:hot start */
/* Current Hessian type */
/* QP user-routine call-status */
/* NP user-routine call-status */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fmul;
    --fstate;
    --f;
    --fupp;
    --flow;
    --xmul;
    --xstate;
    --x;
    --xupp;
    --xlow;
    xnames -= 8;
    fnames -= 8;
    --a;
    --javar;
    --iafun;
    --jgvar;
    --igfun;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNOPTA", (ftnlen)6, (ftnlen)6);
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
/*     Check the arguments of snOptA. */
    s3arga_(&inform__, start, nf, n, ns, nxname, nfname, objrow, nea, neg, &
	    xlow[1], &xupp[1], xnames + 8, &flow[1], &fupp[1], fnames + 8, &
	    xstate[1], &xmul[1], &fstate[1], &fmul[1], &iw[69], &errors, &iw[
	    1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
    if (inform__ > 0) {
	*info = inform__;
	goto L800;
    }
/*     ------------------------------------------------------------------ */
/*     Any objective row specified in the specs file overrides  ObjRow. */
/*     ------------------------------------------------------------------ */
/*     There is always an objective function, even if the user didn't */
/*     specify one. */
    objsav = *objrow;
    objspc = iw[103];
/* Optional parameter: objective row. */
    if (objspc != -11111) {
	*objrow = objspc;
    } else {
	iw[103] = *objrow;
    }
/*     Allocate temporary work arrays for s3sizA. */
    nkx = *n + *nf;
    nlocj = *n + 1;
/*     Permanent addresses first. */
    lkx = nextiw;
    minbld = lkx + nkx;
    if (minbld > maxiw) {
/*        --------------------------------------------------------------- */
/*        Not enough memory to build the problem. */
/*        Provide the user an (over) estimate of what is needed. */
/*        The problem dimensions have not been computed yet, so s3dfltA */
/*        assigns temporary estimates of  lvlHes, maxR, maxS  and  mQNmod */
/*        if they have not yet been specified by the user. */
/*        --------------------------------------------------------------- */
	ne = *nea + *neg;
	m = *nf;
	nlocg = *n + 1;
	if (*nxname == 1 && *nfname == 1) {
	    nname = 1;
	} else {
	    nname = *n + m;
	}
	nncon = m;
	nnjac = *n;
	nnobj = *n;
	negcon = ne;
	s3dflta_(n, &lvlhestmp, &maxrtmp, &maxstmp, &mqnmodtmp, &iw[1], leniw)
		;
	lenrtmp = maxrtmp * (maxrtmp + 1) / 2 + (maxstmp - maxrtmp);
	s8map_(&m, n, &negcon, &nkx, &nncon, &nnjac, &nnobj, &lenrtmp, &
		maxrtmp, &maxstmp, &mqnmodtmp, &lvlhestmp, &nextcw, &nextiw, &
		nextrw, &iw[1], leniw);
	s3mapa_(&m, n, &ne, nf, neg, &negcon, &nkx, &nnjac, &nname, &nextcw, &
		nextiw, &nextrw, &iw[1], leniw);
	s2bmap_(&m, n, &ne, &maxstmp, &nextiw, &nextrw, &maxiw, &maxrw, &
		liwest, &lrwest, &iw[1], leniw);
	prtmem = TRUE_;
/* Print all messages in s2Mem */
	s2mem_(&inform__, &prtmem, &liwest, &lrwest, &nextcw, &nextiw, &
		nextrw, &maxcw, &maxiw, &maxrw, lencw, leniw, lenrw, mincw, 
		miniw, minrw, &iw[1]);
	*info = inform__;
	goto L800;
    }
/*     Compute  m, negCon, ne, nnCon, nnJac, nnObj and iObj. */
/*     The integer array kx defines the order of the variables */
/*     and constraints given to SNOPT. */
    s3siza_(info, n, nf, &nkx, objrow, &iafun[1], &javar[1], lena, nea, &
	    igfun[1], &jgvar[1], leng, neg, &m, &negcon, &ne, &nncon, &nnjac, 
	    &nnobj, &iobj, &iw[lkx], leniw, &iw[1]);
    if (*info > 0) {
	goto L800;
    }
/*     The values of  ne,  nnCon,  nnJac  and  nnObj  are now known */
/*     Load the iw array with various problem dimensions. */
    nnh = max(nnjac,nnobj);
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = m;
/* copy of the number of rows */
    iw[17] = ne;
/* copy of the number of nonzeros in Jcol */
    iw[20] = negcon;
/* # of nonzero elems in J */
    iw[21] = nnjac;
/* # of Jacobian  variables */
    iw[22] = nnobj;
/* # of objective variables */
    iw[23] = nncon;
/* # of nonlinear constraints */
    iw[24] = nnh;
/*   max( nnObj, nnJac ) */
    iw[204] = iobj;
/* position of the objective row in J */
    iw[248] = *nf;
/* # of components of user-defined F */
    iw[249] = *neg;
/*     ------------------------------------------------------------------ */
/*     The obligatory call to snInit has already set the defaults. */
/*     All problem dimensions have been computed. */
/*     Check that the optional parameters have sensible values. */
/*     Print the options. */
/*     ------------------------------------------------------------------ */
/* # of components of user-defined G */
    s_copy(cw + 408, prob, (ftnlen)8, (ftnlen)8);
    s8dflt_(&m, n, &nncon, &nnjac, &nnobj, cw + 8, lencw, &iw[1], leniw, &rw[
	    1], lenrw, (ftnlen)8);
    s3prta_(&m, n, &nncon, &nnjac, &nnobj, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Compute the addresses of all work arrays. */
/*     ------------------------------------------------------------------ */
    nb = *n + m;
    nlocg = nnjac + 1;
    if (*nxname == 1 && *nfname == 1) {
	nname = 1;
    } else {
	nname = nb;
    }
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    iw[28] = lenr;
    if (*objrow == 0) {
	lx0 = nextrw;
	nextrw = lx0 + nnh;
	iw[298] = lx0;
/* ad hoc allocated FP base point */
    }
/*     ------------------------------------------------------------------ */
/*     Allocate the local arrays for snOptA. */
/*     s8Map  maps snOptA integer and double arrays. */
/*     s3mapA maps additional arrays for snOptA. */
/*     s2BMap maps the arrays for the LU routines. */
/*     s2Mem  checks what space is available and prints any messages. */
/*     ------------------------------------------------------------------ */
    s8map_(&m, n, &negcon, &nkx, &nncon, &nnjac, &nnobj, &lenr, &maxr, &maxs, 
	    &mqnmod, &lvlhes, &nextcw, &nextiw, &nextrw, &iw[1], leniw);
    s3mapa_(&m, n, &ne, nf, neg, &negcon, &nkx, &nnjac, &nname, &nextcw, &
	    nextiw, &nextrw, &iw[1], leniw);
    s2bmap_(&m, n, &ne, &maxs, &nextiw, &nextrw, &maxiw, &maxrw, &liwest, &
	    lrwest, &iw[1], leniw);
    prtmem = TRUE_;
/* OK to print messages in s2Mem */
    s2mem_(&inform__, &prtmem, &liwest, &lrwest, &nextcw, &nextiw, &nextrw, &
	    maxcw, &maxiw, &maxrw, lencw, leniw, lenrw, mincw, miniw, minrw, &
	    iw[1]);
    if (inform__ != 0) {
	*info = inform__;
	goto L800;
    }
/*     Allocate local work arrays. */
    lkxn = iw[252];
/* jN = kxN(j ) => col j of Jcol is variable jN */
    ljcol = iw[256];
/* Jcol(ne)    = Constraint Jacobian by columns */
    llocj = iw[257];
/* locJ(n+1)   = column pointers for indJ */
    lindj = iw[258];
/* indJ(ne) holds the row indices for Jij */
    llocg = iw[260];
/* locG(nlocG) = column pointers for indG */
    lindg = iw[261];
/* indG(neG) holds the row indices for gij */
    lnglin = iw[262];
/* nGlin(j) = # linear elems in col j of gCon */
    ligfun = iw[266];
/* iGfun(neG) row list of reordered G nonzeros */
    ljgvar = iw[267];
/* iGvar(neG) col list of reordered G nonzeros */
    lgobj = iw[297];
/* gObj(nnObj) = Objective gradient */
    lgobj1 = iw[324];
/* gObj1(nnObj) objective gradients at x1 */
    lgobj2 = iw[325];
/* gObj2(nnObj) work gObj */
    lbl = iw[271];
/* bl(nb)      = lower bounds */
    lbu = iw[272];
/* bu(nb)      = upper bounds */
    lpi = iw[279];
/* pi(m)       = the pi-vector */
    lrc = iw[280];
/* rc(nb)      = the reduced costs */
    lhs = iw[282];
/* the column state vector */
    lhetyp = iw[283];
/* hEtype(nb) definition of elastic vars */
    lx = iw[299];
/* x(nb)       = the solution (x,s) */
    liy = iw[308];
/* iy (nb)     =  integer work vector */
    lfx = iw[336];
/* Fx (nnCon)  = F(x) + A(linear)x */
    lnames = iw[359];
/*     ------------------------------------------------------------------ */
/*     Build the column-wise data structure for the Jacobian. */
/*     ------------------------------------------------------------------ */
/* Names(nName), row and column names */
    s3blda_(objrow, n, &nkx, &nncon, &nnjac, &iw[lkx], &iw[lnglin], &iafun[1],
	     &javar[1], lena, nea, &a[1], &igfun[1], &jgvar[1], leng, neg, &
	    ne, &nlocj, &iw[llocj], &iw[lindj], &rw[ljcol], &negcon, &nlocg, &
	    iw[llocg], &iw[lindg], &iw[liy]);
/*     ------------------------------------------------------------------ */
/*     Re-order the input data and invert the row and column orderings. */
/*     ------------------------------------------------------------------ */
    s3prma_(n, nf, &nkx, &igfun[1], &jgvar[1], &iw[ligfun], &iw[ljgvar], leng,
	     neg, &iw[lkx], &iw[lkxn]);
    iw[247] = nkx;
/*     ------------------------------------------------------------------ */
/*     Load the arrays used by SNOPTA. */
/*     These are for the data, */
/*              Jcol, indJ, locJ, bl, bu */
/*     and for the solution */
/*              hs, x, pi, rc, hs. */
/*     ------------------------------------------------------------------ */
/* dimension of kx and its inverse, kxN */
    infbnd = rw[70];
/* definition of an infinite bound. */
    s3ina_(&iw[69], &iobj, &m, n, &nb, &nncon, nf, &nkx, &iw[lkxn], &infbnd, 
	    xnames + 8, fnames + 8, cw + (lnames << 3), nxname, nfname, &
	    nname, &xlow[1], &xupp[1], &flow[1], &fupp[1], &rw[lbl], &rw[lbu],
	     &xstate[1], &fstate[1], &iw[lhs], &x[1], &f[1], &rw[lx], &rw[lfx]
	    , &fmul[1], &rw[lpi], (ftnlen)8, (ftnlen)8, (ftnlen)8);
    if (*objrow == 0) {
	dcopy_(&nb, &rw[lx], &c__1, &rw[lx0], &c__1);
    }
/*     ------------------------------------------------------------------ */
/*     Sparse obj. gradients are scattered into gObj, gObj1 and gObj2. */
/*     ------------------------------------------------------------------ */
    dload_(&nnobj, &c_b20, &rw[lgobj], &c__1);
    dload_(&nnobj, &c_b20, &rw[lgobj1], &c__1);
    dload_(&nnobj, &c_b20, &rw[lgobj2], &c__1);
/*     ------------------------------------------------------------------ */
/*     The problem has been built.  The call-status variables npstat and */
/*     qpstat are reset to ensure that the necessary housekeeping is done */
/*     on the first call to s0fgA. */

/*     In future versions, the function build will be separated from the */
/*     snoptA. This will allow the function build to be called from the */
/*     function wrapper s0fgA. */
/*     ------------------------------------------------------------------ */
    iw[235] = -11111;
    iw[236] = -11111;
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
	nnobj0 = max(nnobj,1);
	nnh0 = max(nnh,1);
	ngqp = max(nnobj,nnh);
	iload_(&nb, &c__3, &iw[lhetyp], &c__1);
	s5solv_(info, solver, &iw[69], (U_fp)sqhx_, (U_fp)s8qphx_, (U_fp)
		sqlog, &gotr, &m, n, &nb, &nnh0, &nnh, &nname, &ngqp, &nnobj0,
		 &nnobj, &iobj, objadd, &fobj, &objtru, ninf, sinf, &ne, &
		nlocj, &iw[llocj], &iw[lindj], &rw[ljcol], &rw[lbl], &rw[lbu],
		 &rw[lgobj], cw + (lnames << 3), &nrhs0, &nrhs, rhs, &lenx0, &
		nx0, x0, &iw[lhetyp], &iw[lhs], &rw[lx], &rw[lpi], &rw[lrc], 
		ns, cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)6, (ftnlen)8, (
		ftnlen)8, (ftnlen)8);
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
	s8solv_(info, solver, &iw[69], (U_fp)s0fga_, (U_fp)usrfun, (U_fp)
		usrfun, (U_fp)snlog, (U_fp)snlog2, (U_fp)snstop, &gotr, &m, n,
		 &nb, &nncon, &nnjac, &nnobj, &nname, &iobj, objadd, &fobj, &
		objtru, ninf, sinf, &ne, &nlocj, &iw[llocj], &iw[lindj], &rw[
		ljcol], &rw[lbl], &rw[lbu], cw + (lnames << 3), &iw[lhs], &rw[
		lx], &rw[lpi], &rw[lrc], &nmajor, ns, cu + 8, lencu, &iu[1], 
		leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		lenrw, (ftnlen)6, (ftnlen)8, (ftnlen)8, (ftnlen)8);
    }
/*     ------------------------------------------------------------------ */
/*     Restore and update the user data. */
/*     ------------------------------------------------------------------ */
    s3outa_(n, &nb, nf, &nncon, &nkx, &iw[lkxn], objrow, &iobj, &fobj, &
	    xstate[1], &fstate[1], &iw[lhs], &x[1], &f[1], &rw[lx], &rw[lfx], 
	    &xmul[1], &fmul[1], &rw[lrc]);
    *objrow = objsav;
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
} /* snkera_ */

