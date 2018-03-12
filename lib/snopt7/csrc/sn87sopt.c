/* ../snopt7/src/sn87sopt.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static doublereal c_b4 = 0.;
static integer c__4 = 4;
static integer c__0 = 0;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c_n2 = -2;
static integer c__13 = 13;
static integer c_n3 = -3;
static doublereal c_b129 = .33333;
static integer c__6 = 6;
static doublereal c_b168 = 1.;
static doublereal c_b169 = -1.;
static integer c__23 = 23;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn87sopt.f */

/*     s8solv */
/*     s8dflt   s8Map   s8SQP   s8Stat */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8solv_(integer *iexit, char *solver, integer *start, 
	S_fp fgwrap, U_fp funcon, U_fp funobj, U_fp mjrlog, U_fp mnrlog, U_fp 
	snstop, logical *gotr, integer *m, integer *n, integer *nb, integer *
	nncon, integer *nnjac, integer *nnobj, integer *nname, integer *iobj, 
	doublereal *objadd, doublereal *fobj, doublereal *objtru, integer *
	ninf, doublereal *sinf, integer *ne, integer *nlocj, integer *locj, 
	integer *indj, doublereal *jcol, doublereal *bl, doublereal *bu, char 
	*names, integer *hs, doublereal *x, doublereal *pi, doublereal *rc, 
	integer *nmajor, integer *ns, char *cu, integer *lencu, integer *iu, 
	integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen solver_len, ftnlen names_len, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1900[] = "(\002 No. of iterations\002,i20,2x,\002 Object"
	    "ive value\002,1p,e22.10)";
    static char fmt_1910[] = "(\002 No. of infeasibilities\002,i15,2x,\002 S"
	    "um of infeas\002,1p,e24.10)";
    static char fmt_1915[] = "(\002 Elastic weight            \002,1p,e11.1,"
	    "2x,\002 Scaled Merit \002,1p,e24.10)";
    static char fmt_1920[] = "(\002 No. of major iterations\002,i14,2x,\002 "
	    "Linear objective\002,1p,e21.10)";
    static char fmt_1930[] = "(\002 Penalty parameter\002,1p,e20.3,2x,\002 N"
	    "onlinear objective\002,1p,e18.10)";
    static char fmt_1950[] = "(\002 No. of calls to funobj\002,i15,2x,\002 N"
	    "o. of calls to funcon\002,i15)";
    static char fmt_1955[] = "(\002 Calls with modes 1,2 (known g)\002,i7,"
	    "2x,\002 Calls with modes 1,2 (known g)\002,i7)";
    static char fmt_1960[] = "(\002 Calls for forward differencing\002,i7,"
	    "2x,\002 Calls for forward differencing\002,i7)";
    static char fmt_1962[] = "(\002 Calls for central differencing\002,i7,"
	    "2x,\002 Calls for central differencing\002,i7)";
    static char fmt_1970[] = "(\002 No. of superbasics\002,i19,2x,\002 No. o"
	    "f basic nonlinears\002,i14)";
    static char fmt_1973[] = "(\002 No. of CG iterations\002,i17)";
    static char fmt_1975[] = "(\002 No. of degenerate steps\002,i14,2x,\002 "
	    "Percentage\002,f27.2)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    doublereal d__1;
    char ch__1[38];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, lr, ly, lx1, ly1, ly2, nx0, ldg, nnb, mbs, lrg, ldx, 
	    lfv, lfx, itn, nnl, liy, lkx, nkx;
    static char str[133];
    static integer lrg2, nnl0;
    static doublereal eps0;
    static integer liy1, liy2;
    extern /* Subroutine */ int s8fx_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static char str2[133];
    static integer lgbs, lkbs;
    static doublereal flin;
    static integer lenr, lhdx, lgqp, lpbs, lxbs, maxs, nrhs, ludx;
    static doublereal fmrt;
    static integer lxqp;
    static doublereal tolx;
    extern /* Subroutine */ int s6fdg_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    S_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static integer lenx0, nrhs0, lxqp0;
    static logical needb;
    extern /* Subroutine */ int s8sqp_(integer *, S_fp, U_fp, U_fp, U_fp, 
	    U_fp, U_fp, logical *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal degen;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *), 
	    dddiv_(integer *, doublereal *, integer *, doublereal *, integer *
	    ), ddscl_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer lgobj, lblbs, llocg, lfcon, lgcon, inewb, lbubs, nlocg;
    static char mprob[8];
    static integer numlc;
    static doublereal duinf, vilim;
    static integer lycon, lxpen, maxvi;
    static doublereal tolfp, vimax, virel, wtinf, tolqp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal xnorm, visup;
    static integer lgobj1, lgobj2, lfcon1, lfcon2, lgcon1, lgcon2, nnobj0, 
	    nnobj1, nncon0, nncon1;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    s2amat_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *), s8rand_(integer *, integer *, 
	    doublereal *), s5getb_(integer *, integer *, U_fp, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, char *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    extern doublereal dnrm1s_(integer *, doublereal *, integer *);
    static integer lycon1, lycon2;
    static doublereal wtinf0;
    extern /* Subroutine */ int s8feas_(integer *, U_fp, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen),
	     s2scla_(integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *), s7chkg_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    S_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen), s2scal_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *);
    static doublereal pnorm1, pnorm2;
    extern /* Subroutine */ int s8gcpy_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *), 
	    s8sclj_(integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *), s8sclg_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *), s2vmax_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s2crsh_(integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), s1time_(integer *, integer *
	    , integer *, integer *, doublereal *, integer *), s4stat_(integer 
	    *, char *, ftnlen), s5fixx_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    , s4newb_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, char *, char *, 
	    integer *, integer *, integer *, ftnlen, ftnlen), s4savb_(integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static integer lascal, ndegen, modefg, lhfeas;
    static doublereal infbnd;
    static integer icrash, negcon, lcrash;
    static doublereal sclobj;
    static logical lininf;
    static integer lhesta;
    static logical nlnobj, nlninf;
    static integer lblsav;
    static char istate[4*3];
    static logical nlncon, nonlin;
    static integer inform__, itnlim, lssave;
    static logical gotfun, fponly;
    static integer lvlsch, lbusav, ldycon, lgconu, lhetyp, lqprhs, minimz, 
	    mjrprt, mnrprt, nminor, numliq;
    static doublereal pennrm, pinorm, rgnorm, tcrash;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), snprnt_(integer *, 
	    char *, integer *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___117 = { 0, str, 0, fmt_1900, 133, 1 };
    static icilist io___118 = { 0, str, 0, fmt_1910, 133, 1 };
    static icilist io___119 = { 0, str, 0, fmt_1915, 133, 1 };
    static icilist io___120 = { 0, str, 0, fmt_1920, 133, 1 };
    static icilist io___121 = { 0, str, 0, fmt_1930, 133, 1 };
    static icilist io___122 = { 0, str, 0, fmt_1950, 133, 1 };
    static icilist io___123 = { 0, str, 0, fmt_1955, 133, 1 };
    static icilist io___124 = { 0, str, 0, fmt_1960, 133, 1 };
    static icilist io___125 = { 0, str, 0, fmt_1962, 133, 1 };
    static icilist io___126 = { 0, str, 0, fmt_1970, 133, 1 };
    static icilist io___127 = { 0, str, 0, fmt_1973, 133, 1 };
    static icilist io___128 = { 0, str, 0, fmt_1975, 133, 1 };


/*     ================================================================== */
/*     s8solv solves the current problem. */

/*     On entry, */
/*     the specs file has been read, */
/*     all data items have been loaded (including locJ, indJ, Jcol, ...), */
/*     and workspace has been allocated. */

/*     On exit, */
/*     iExit  =  0 if an optimal solution was found, */
/*            =  1 if the problem was infeasible, */
/*            =  2 if the problem was unbounded, */
/*            =  3 if the Iteration limit was exceeded, */
/*           ge  4 if iterations were terminated by some other */
/*                 error condition (see the SNOPT user's guide). */

/*     15 Nov 1991: First version based on Minos 5.4 routine misolv. */
/*     13 Feb 1994: Eliminated "Cycle" options. */
/*                  Simplified s4getb. */
/*     12 Nov 1994: Integer workspace added. */
/*     25 Jul 1996: Sign of the slacks changed. */
/*     28 Sep 1997: Character workspace added. */
/*     11 Nov 1997: Backtracking for undefined functions. */
/*     26 Dec 1997: Dummy Jacobian scaled in feasibility phase. */
/*     27 Aug 1998: Constant Jacobian elements handled correctly. */
/*     10 Oct 1998: Objective and constraint gradient checking merged. */
/*     11 Oct 1998: Facility to combine funobj and funcon added. */
/*     23 Dec 1999: Suboptimize option added. */
/*     30 Dec 2000: Housekeeping for first function call moved to snwrap. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     19 Mar 2006: Fx initialized for s8savB. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     18 Jun 2008: iy2, pBS and rg2 added to workspace. */
/*     23 Oct 2010: pinorm initialized at 1.0. */
/*     01 Dec 2012: heStat initialized in s8solv. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* = 0,1,2,3 or 4, deriv level */
/* 0,1,2  => LM, FM, Newton */
/* scale option */
/* 1, 0, -1  => MIN, FP, MAX */
/* =1(2) for forwd(cntrl) diffs */
/* > 0 => some exact derivs */
/* number of calls of fCon */
/* number of calls of fCon */
/* number of calls of fCon */
/* number of calls of fCon */
/* number of calls of fObj */
/* number of calls of fObj */
/* number of calls of fObj */
/* number of calls of fObj */
/* =1(0) for pd  QP Hessian */
/* # of LU factorizations */
/* # lines in log     file */
/* # lines in summary file */
/* NP user-routine call-status */
/* Number of symmlq iterations */
/* symmlq itns for last minor */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --bu;
    --bl;
    names -= 8;
    --jcol;
    --indj;
    --locj;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    inewb = iw[124];
/* new basis file */
    negcon = iw[20];
/* # of nonzero elems in J */
    lenr = iw[28];
/* R(lenR) is the reduced Hessian factor */
    maxs = iw[53];
/* max # of superbasics */
    lvlsch = iw[76];
/* >0     => use derivatives in the line search */
    icrash = iw[88];
/* Crash option */
    itnlim = iw[89];
/* limit on total iterations */
    mjrprt = iw[92];
/* Major print level */
    mnrprt = iw[93];
/* Minor print level */
    minimz = iw[199];
/* 1 (-1)    => minimize (maximize) */
    nkx = iw[247];
/* dimension of kx and its inverse, kxN */
    eps0 = rw[2];
    tolfp = rw[51];
/* Minor Phase 1 Opt tol */
    tolqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    tcrash = rw[62];
/* crash tolerance. */
    infbnd = rw[70];
/* definition of an infinite bound. */
    vilim = rw[81];
/* violation limit */
    wtinf0 = rw[88];
/* infeasibility weight */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/*     Addresses */
/* Problem name */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    llocg = iw[260];
/* locG(nnJac+1) = column pointers for indG */
    lblbs = iw[273];
/* blBS(mBS)   = lower bounds for xBS */
    lbubs = iw[274];
/* buBS(mBS)   = upper bounds for xBS */
    lblsav = iw[275];
/* blSav(nb)   = copy of bl */
    lbusav = iw[276];
/* buSav(nb)   = copy of bu */
    lpbs = iw[277];
/* pBS(nb)     = search direction */
    lqprhs = iw[278];
/* QPrhs(nnCon)=  QP constraint rhs */
    lhetyp = iw[283];
/* hEtype(nb) list of elastic vars */
    lhfeas = iw[284];
/* hfeas(mBS), feasibility types */
    lhesta = iw[285];
/* hEstat(nb), status of elastics */
    ldx = iw[287];
/* dx(nb)      = x1 - x */
    lhdx = iw[288];
/* Hdx(nnL)    = product of H with  x1 - x */
    ldg = iw[289];
/* dg(nnL)     = gradient difference */
    lgqp = iw[290];
/* gQP(ngQP)   = QP gradient */
    lgbs = iw[291];
/* gBS(mBS)    = BS components of g */
    lkbs = iw[292];
/* kBS(mBS), ( B  S ) list */
    lrg = iw[293];
/* rg (maxS)   = reduced gradient */
    lrg2 = iw[294];
/* rg2(maxS)   = reduced gradient */
    lr = iw[295];
/* R(lenR)     = factor of Z'HZ */
    lascal = iw[296];
/* Ascale(nb)  = row and column scales */
    lgobj = iw[297];
/* gObj(nnObj) = Objective gradient */
    lx1 = iw[300];
/* x1(nb)      = new x, used to store x0 */
    lxbs = iw[301];
/* xBS(mBS)    = basics, superbasics */
    lxpen = iw[304];
/* xPen(nnCon) = penalty params */
    lxqp = iw[305];
/* xQP(nb)     = QP solution */
    lxqp0 = iw[306];
/* xQP0(nb)    = QP feasible pt. */
    liy = iw[308];
/* iy(nb)      =  integer work vector */
    liy1 = iw[309];
/* iy1(nb)     =  integer work vector */
    liy2 = iw[310];
/* iy2(nb)     =  integer work vector */
    ly = iw[311];
/* y(nb)       =  real work vector */
    ly1 = iw[312];
/* y1(nb)      =  real work vector */
    ly2 = iw[313];
/* y2(nb)      =  real work vector */
    lfcon = iw[316];
/* fCon (nnCon) constraints at x */
    lfcon1 = iw[317];
/* fCon1(nnCon) constraints at x1 */
    lfcon2 = iw[318];
/* fCon2(nnCon) work vector */
    lgconu = iw[319];
/* record of unknown derivatives and constants */
    lgcon = iw[320];
/* gCon (negCon)   constraint gradients at x */
    lgcon1 = iw[321];
/* gCon1(negCon)   constraint gradients at x1 */
    lgcon2 = iw[322];
/* gCon2(negCon)   work vector */
    lgobj1 = iw[324];
/* gObj1(nnObj) objective gradients at x1 */
    lgobj2 = iw[325];
/* gObj2(nnObj) work gObj */
    lfx = iw[336];
/* Fx (nnCon)  = F(x) + A(linear)x */
    lfv = iw[337];
/* Fv          = F(x) + A(linear)x - sN */
    ludx = iw[345];
/* Udx(nnL)      = product of U with dx */
    lycon = iw[348];
/* yCon (nnCon)  = multipliers for F */
    lycon1 = iw[349];
/* yCon1(nnCon)  = yCon at x1 */
    lycon2 = iw[350];
/* yCon2(nnCon)  = work copy of yCon */
    ldycon = iw[351];
/* dyCon(nnCon)  = yCon1 - yCon */
    *iexit = 0;
    gotfun = FALSE_;
    fponly = iw[87] == 0;
    nnl = max(*nnjac,*nnobj);
    nlncon = *nncon > 0;
    nlnobj = *nnobj > 0;
    nonlin = nnl > 0;
    nnobj0 = max(*nnobj,1);
    nncon0 = max(*nncon,1);
    nnl0 = max(nnl,1);
    mbs = *m + maxs;
    nlocg = *nnjac + 1;
    numlc = *m - *nncon;
/*     Initialize yCon from pi. */
/*     Zap the pi(i) to prevent them being printed without being set. */
    if (nlncon) {
	dcopy_(nncon, &pi[1], &c__1, &rw[lycon], &c__1);
    }
    if (numlc > 0) {
	dload_(&numlc, &c_b4, &pi[*nncon + 1], &c__1);
    }
/*     Initialize a few things. */
/*     Define the Hessian type for the QP subproblem. */
    if (iw[72] == 0 || iw[72] == 1) {
	if (nnl < *n) {
	    iw[200] = 0;
	} else {
	    iw[200] = 1;
	}
    }
    iw[181] = 1;
    iw[210] = 0;
    iw[386] = 0;
    iw[387] = 0;
    *ninf = 0;
    wtinf = 1.;
    pinorm = 1.;
    duinf = 0.;
    fmrt = 0.;
    *fobj = 0.;
    *objtru = 0.;
    pennrm = 0.;
    vimax = 0.;
    virel = 0.;
    iw[220] = 0;
/* Line count for the print   file */
    iw[221] = 0;
/* Line count for the summary file */
    iload_(&c__4, &c__0, &iw[189], &c__1);
    iload_(&c__4, &c__0, &iw[194], &c__1);
    itn = 0;
    ndegen = 0;
    *nmajor = 0;
    nminor = 0;
    s1page_(&c__1, &iw[1], leniw);
/*     ------------------------------------------------------------------ */
/*     Print the matrix statistics before the nonlinear part of Jcol is */
/*     loaded with random elements. */
/*     Find the rowtypes for use in s5getB (they are held in iy2). */
/*     ------------------------------------------------------------------ */
    s2amat_(&c__1, &mjrprt, m, n, nb, nncon, nnjac, nnobj, iobj, &numlc, &
	    numliq, ne, nlocj, &locj[1], &indj[1], &jcol[1], &bl[1], &bu[1], &
	    iw[liy2], &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Make a permanent copy in gConu of the constant Jacobian elements */
/*     stored in J.  Load the nonlinear part of J with random elements. */
/*     ------------------------------------------------------------------ */
    if (nlncon) {
	s8gcpy_(nncon, nnjac, ne, nlocj, &locj[1], &indj[1], ne, nlocj, &locj[
		1], &jcol[1], &negcon, &nlocg, &iw[llocg], &rw[lgconu]);
	s8rand_(&negcon, &negcon, &rw[lgcon]);
	s8gcpy_(nncon, nnjac, ne, nlocj, &locj[1], &indj[1], &negcon, &nlocg, 
		&iw[llocg], &rw[lgcon], ne, nlocj, &locj[1], &jcol[1]);
    }
/* ================================================================= */
/* Find a basis kBS(1:m) for the linear constraints and bounds. */
/* ================================================================= */
/* s5getB does the following. */
/*  1. The linear constraints are (optionally) scaled. */
/*  2. The bounds bl and bu on the nonlinear rows are relaxed. */
/*  3. Elements x(n+1:n+m) of the initial x are assigned. */
/*  4. An LP is used to find a feasible x for the bounds and */
/*     linear equality constraints. */
/*  5. x(nb) is (optionally) scaled and saved in x1. */
    nrhs = 0;
/* No QP rhs when finding the first basis */
    nrhs0 = 1;
    nx0 = *nb;
/* elements in x1(nb), the base point */
    lenx0 = *nb;
    iload_(nb, &c__0, &iw[lhetyp], &c__1);
/* placeholder for calls in s5getB */
    iload_(nb, &c__0, &iw[lhesta], &c__1);
    s5getb_(&inform__, start, (U_fp)mnrlog, &needb, m, &maxs, &mbs, n, nb, 
	    nncon, nnjac, nnobj, nname, ns, &nminor, &itnlim, &itn, &ndegen, &
	    numlc, &numliq, &tolfp, &tolqp, &tolx, ninf, sinf, &wtinf, iobj, &
	    sclobj, &pinorm, &rgnorm, ne, nlocj, &locj[1], &indj[1], &jcol[1],
	     &iw[lhetyp], &iw[lhesta], &iw[liy2], &iw[lhfeas], &hs[1], &iw[
	    lkbs], names + 8, &rw[lascal], &bl[1], &bu[1], &rw[lblbs], &rw[
	    lbubs], &rw[lblsav], &rw[lbusav], &rw[lgbs], &pi[1], &rc[1], &
	    nrhs0, &nrhs, &rw[lqprhs], &lenx0, &nx0, &rw[lx1], &x[1], &rw[
	    lxbs], &iw[liy], &iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cw + 8, 
	    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
/* Possible inform values are = -3,-2,-1, 0, >0 */
    if (inform__ != 0) {
	if (inform__ > 0) {
	    *iexit = inform__;
/* fatal error */
	} else if (inform__ == -3) {
	    *iexit = 31;
/* too many iterations */
	} else {
	    *iexit = 12;
/* infeasible linear equalities */
	}
    }
/* ================================================================= */
/* An initial basis has been assigned to  kBS(1:m). */

/* Find a feasible point for all the linear constraints. */
/* The norm of x is minimized via a proximal-point QP. */
/* If there is no feasible point, the linear rows can be elastic. */
/* ================================================================= */
    if (numlc > 0) {
	if (*iexit == 0) {
/* the E rows are feasible, check LG rows */
	    s8feas_(iexit, (U_fp)mnrlog, &lenr, m, &maxs, &mbs, n, nb, &
		    nncon0, nncon, &nnl0, &nnl, &ndegen, ns, &numlc, &numliq, 
		    &itn, &itnlim, &nminor, &mnrprt, &sclobj, &tolqp, &tolx, 
		    ninf, sinf, &wtinf, &pinorm, &rgnorm, ne, nlocj, &locj[1],
		     &indj[1], &jcol[1], &iw[lhetyp], &iw[lhesta], &iw[lhfeas]
		    , &hs[1], &iw[lkbs], &rw[lascal], &bl[1], &bu[1], &rw[
		    lblsav], &rw[lbusav], &rw[lblbs], &rw[lbubs], &rw[lgbs], &
		    rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], &
		    rw[lrg], &rw[lqprhs], &rw[lx1], &x[1], &rw[lxbs], &iw[liy]
		    , &iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cw + 8, lencw, &
		    iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
	}
/* -------------------------------------------------------------- */
/* Do some housekeeping in case we have to exit */
/* -------------------------------------------------------------- */
/* Reinstate the scaled bounds on the nonlinear constraints. */
	if (nlncon) {
	    dcopy_(nncon, &rw[lblsav + *n], &c__1, &bl[*n + 1], &c__1);
	    dcopy_(nncon, &rw[lbusav + *n], &c__1, &bu[*n + 1], &c__1);
	}
/* Unscale the linear constraints. */
/* nlnCon */
	if (iw[75] > 0) {
	    s2scla_(&c__1, m, n, nb, iobj, &infbnd, &sclobj, ne, nlocj, &locj[
		    1], &indj[1], &jcol[1], &rw[lascal], &bl[1], &bu[1], &pi[
		    1], &x[1]);
	}
    }
/*     Exit if the linear constraints are infeasible. */
/* numLC > 0 */
    lininf = *iexit != 0;
    if (lininf) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Copy the constant Jacobian elements into gCon, gCon1 and gCon2. */
/*     Reset hEtype so that only nonlinear rows are elastic. */
/*     Make sure variables are not outside their bounds */
/*     (in particular, check the nonlinear slacks). */
/*     ------------------------------------------------------------------ */
    if (nlncon) {
	dcopy_(&negcon, &rw[lgconu], &c__1, &rw[lgcon], &c__1);
	dcopy_(&negcon, &rw[lgconu], &c__1, &rw[lgcon1], &c__1);
	dcopy_(&negcon, &rw[lgconu], &c__1, &rw[lgcon2], &c__1);
	iload_(nncon, &c__3, &iw[lhetyp + *n], &c__1);
    }
/* nlnCon */
    s5fixx_(&c__0, &c__1, nb, &tolx, &hs[1], &bl[1], &bu[1], &x[1]);
/*     ================================================================== */
/*     ================================================================== */
/*     The linear constraints have been satisfied! */
/*     Compute the problem functions at this all-important point. */
/*     No scaling yet. */
/*     ================================================================== */
/*     ================================================================== */
    if (nnl > 0) {
	lssave = iw[75];
	iw[75] = 0;
	modefg = 2;
	(*fgwrap)(&inform__, &modefg, &nlncon, &nlnobj, n, &negcon, &nncon0, 
		nncon, nnjac, &nnl, &nnobj0, nnobj, (U_fp)funcon, (U_fp)
		funobj, &x[1], ne, nlocj, &locj[1], &indj[1], &rw[lfcon], 
		fobj, &rw[lgcon], &rw[lgobj], cu + 8, lencu, &iu[1], leniu, &
		ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)8, (ftnlen)8);
	gotfun = inform__ == 0;
	if (! gotfun) {
	    if (inform__ < 0) {
		if (numlc > 0) {
		    *iexit = 61;
/* Undefined fun at first feasible point */
		} else {
		    *iexit = 62;
/* Undefined fun at initial point */
		}
	    } else {
		*iexit = inform__;
/* User wants to stop */
	    }
	    goto L900;
	}
	if (nlncon) {
/* Define Fx for s8savB */
	    dcopy_(nncon, &rw[lfcon], &c__1, &rw[lfx], &c__1);
	}
/*        --------------------------------------------------------------- */
/*        Check derivatives. */
/*        (One day, we will do this on the SCALED problem.) */
/*        --------------------------------------------------------------- */
	s7chkg_(iexit, n, &nncon0, nncon, nnjac, &nnl, &nnobj0, nnobj, (S_fp)
		fgwrap, (U_fp)funcon, (U_fp)funobj, &x[1], &rw[lx1], &bl[1], &
		bu[1], fobj, &rw[lgobj], ne, nlocj, &locj[1], &indj[1], &
		negcon, &nlocg, &iw[llocg], &rw[lfcon], &rw[lgcon], &rw[
		lgobj2], &rw[lfcon2], &rw[lgcon2], &rw[ly], &rw[ly1], cu + 8, 
		lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	if (*iexit != 0) {
	    goto L900;
	}
/*        --------------------------------------------------------------- */
/*        Compute any missing derivatives. */
/*        Load the Jacobian gCon in  J. */
/*        --------------------------------------------------------------- */
	s6fdg_(iexit, n, &negcon, &nncon0, nncon, nnjac, &nnl, &nnobj0, nnobj,
		 (S_fp)fgwrap, (U_fp)funcon, (U_fp)funobj, &bl[1], &bu[1], &x[
		1], ne, nlocj, &locj[1], &indj[1], &rw[lfcon], fobj, &rw[
		lgcon], &rw[lgobj], &rw[ly1], cu + 8, lencu, &iu[1], leniu, &
		ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)8, (ftnlen)8);
	if (*iexit != 0) {
	    goto L900;
	}
	if (nlncon) {
	    s8gcpy_(nncon, nnjac, ne, nlocj, &locj[1], &indj[1], &negcon, &
		    nlocg, &iw[llocg], &rw[lgcon], ne, nlocj, &locj[1], &jcol[
		    1]);
	    dcopy_(nncon, &rw[lfcon], &c__1, &rw[lfx], &c__1);
	}
	iw[75] = lssave;
    }
/*     ================================================================== */
/*     Scale the problem. */
/*     ================================================================== */
    if (iw[75] > 0) {
/*        --------------------------------------------------------------- */
/*        Recompute the vector of row types. */
/*        --------------------------------------------------------------- */
	s2amat_(&c__0, &mjrprt, m, n, nb, nncon, nnjac, nnobj, iobj, &numlc, &
		numliq, ne, nlocj, &locj[1], &indj[1], &jcol[1], &bl[1], &bu[
		1], &iw[liy2], &iw[1], leniw, &rw[1], lenrw);
	s2scal_(&mjrprt, m, n, nb, &nnl, nncon, nnjac, &iw[liy2], ne, nlocj, &
		locj[1], &indj[1], &jcol[1], &rw[lascal], &bl[1], &bu[1], &rw[
		ly], &rw[ly2], &iw[1], leniw, &rw[1], lenrw);
	s2scla_(&c__0, m, n, nb, iobj, &infbnd, &sclobj, ne, nlocj, &locj[1], 
		&indj[1], &jcol[1], &rw[lascal], &bl[1], &bu[1], &pi[1], &x[1]
		);
/*        --------------------------------------------------------------- */
/*        The objective and constraint functions haven't been scaled yet. */
/*        Scale the constant elements in gCon1 and gCon2. */
/*        Don't forget the initial pi. */
/*        --------------------------------------------------------------- */
	if (nlncon) {
	    dddiv_(nncon, &rw[lascal + *n], &c__1, &rw[lfcon], &c__1);
	    if (iw[184] > 0) {
		s8sclj_(nncon, nnjac, &negcon, n, &rw[lascal], ne, nlocj, &
			locj[1], &indj[1], &rw[lgcon], &rw[1], lenrw);
		dcopy_(&negcon, &rw[lgcon], &c__1, &rw[lgcon1], &c__1);
		dcopy_(&negcon, &rw[lgcon], &c__1, &rw[lgcon2], &c__1);
	    }
	    ddscl_(nncon, &rw[lascal + *n], &c__1, &rw[lycon], &c__1);
	}
	if (nlnobj && iw[184] > 0) {
	    s8sclg_(nnobj, &rw[lascal], &rw[lgobj], &rw[1], lenrw);
	}
    }
/*     ================================================================== */
/*     s8Fx computes the nonlinear constraint values Fx. */
/*     Copy these into the slacks x(n+i) and make sure they are feasible. */
/*     Crash uses them to decide which slacks to grab for the basis */
/*     If any nonbasic nonlinear slacks are close to a bound, */
/*     move them exactly onto the bound to avoid very small steps. */
/*     ================================================================== */
/* iw(lvlScl) > 0 */
    if (nlncon) {
	s8fx_(n, nncon, nnjac, &eps0, ne, nlocj, &locj[1], &indj[1], &jcol[1],
		 &rw[lfcon], &x[1], &rw[lfx]);
	s2vmax_(n, nncon, &maxvi, &vimax, &bl[1], &bu[1], &rw[lfx]);
/* Computing MAX */
	d__1 = vimax * 10.;
	visup = max(d__1,vilim);
	dcopy_(nncon, &rw[lfx], &c__1, &x[*n + 1], &c__1);
	i__1 = *n + 1;
	i__2 = *n + *nncon;
	s5fixx_(&c__1, &i__1, &i__2, &tolx, &hs[1], &bl[1], &bu[1], &x[1]);
/*        =============================================================== */
/*        Crash on the nonlinear rows. */
/*        hs(*) already defines a basis for the full problem,  but we */
/*        want to do better by not including all of the slacks. */
/*        =============================================================== */
	if (needb) {
/*           Load  iy2  with the row types. */
/*           s2crsh uses kBS as workspace.  It may alter x(n+i) for */
/*           nonlinear slacks. */
	    s2amat_(&c__0, &mjrprt, m, n, nb, nncon, nnjac, nnobj, iobj, &
		    numlc, &numliq, ne, nlocj, &locj[1], &indj[1], &jcol[1], &
		    bl[1], &bu[1], &iw[liy2], &iw[1], leniw, &rw[1], lenrw);
	    lcrash = 5;
	    s2crsh_(&lcrash, &mjrprt, m, n, nb, nncon, &icrash, &tcrash, ne, 
		    nlocj, &locj[1], &indj[1], &jcol[1], &iw[lkbs], &hs[1], &
		    iw[liy2], &bl[1], &bu[1], &x[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    needb = FALSE_;
	}
/* needB */
    }
/*     ------------------------------------------------------------------ */
/*     Solve the problem. */
/*     ------------------------------------------------------------------ */
/* nlnCon */
    s1page_(&c__1, &iw[1], leniw);
    s1time_(&c__2, &c__0, &iw[1], leniw, &rw[1], lenrw);
    s8sqp_(iexit, (S_fp)fgwrap, (U_fp)funcon, (U_fp)funobj, (U_fp)mjrlog, (
	    U_fp)mnrlog, (U_fp)snstop, gotr, &itn, &lenr, m, &maxs, &mbs, n, 
	    nb, ns, &nncon0, nncon, &nnobj0, nnobj, &nnl0, &nnl, nmajor, &
	    nminor, &ndegen, &duinf, &minimz, iobj, &sclobj, objadd, fobj, &
	    fmrt, &vimax, &virel, &visup, ninf, sinf, &wtinf0, &wtinf, &
	    pennrm, &pinorm, &xnorm, ne, nlocj, &locj[1], &indj[1], &jcol[1], 
	    &negcon, &nlocg, &iw[llocg], &iw[lhetyp], &iw[lhesta], &iw[lhfeas]
	    , &hs[1], &iw[lkbs], &rw[lascal], &bl[1], &bu[1], &rw[lblbs], &rw[
	    lbubs], &rw[lfv], &rw[lfx], &rw[lfcon], &rw[lgcon], &rw[lgobj], &
	    rw[lfcon1], &rw[lgcon1], &rw[lgobj1], &rw[lfcon2], &rw[lgcon2], &
	    rw[lgobj2], &rw[lgbs], &rw[lgqp], &rw[ldycon], &rw[ldx], &rw[ldg],
	     &rw[ludx], &rw[lhdx], &rw[lpbs], &rw[lycon], &rw[lycon1], &rw[
	    lycon2], &pi[1], &rw[lqprhs], &rw[lr], &rc[1], &rw[lrg], &rw[lrg2]
	    , &x[1], &rw[lx1], &rw[lxbs], &rw[lxqp0], &rw[lxqp], &rw[lxpen], &
	    iw[liy], &iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cu + 8, lencu, &
	    iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1],
	     lenrw, (ftnlen)8, (ftnlen)8);
    s1time_(&c_n2, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     ================================================================== */
/*     Exit. */
/*     Set output variables and print a summary of the final solution. */
/*     ================================================================== */
L900:
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)133, (
	    ftnlen)133);
    degen = ndegen * 100. / max(itn,1);
    if (*iobj == 0) {
	flin = *objadd;
    } else {
	flin = *objadd + x[*n + *iobj] * sclobj;
    }
    *objtru = flin + *fobj;
    nlninf = *ninf > 0;
    xnorm = dnrm1s_(n, &x[1], &c__1);
/*     Count basic nonlinear variables (used only for printing). */
    nnb = 0;
    i__1 = nnl;
    for (j = 1; j <= i__1; ++j) {
	if (hs[j] == 3) {
	    ++nnb;
	}
    }
    if (inewb > 0 && *iexit / 10 < 8) {
	k = *iexit / 10 + 1;
	s4stat_(&k, istate, (ftnlen)4);
	s4newb_(&c__1, &inewb, &minimz, m, n, nb, ns, &mbs, &itn, ninf, sinf, 
		fobj, &iw[lkbs], &hs[1], &rw[lascal], &bl[1], &bu[1], &x[1], &
		rw[lxbs], istate, cw + 8, lencw, &iw[1], leniw, (ftnlen)4, (
		ftnlen)8);
    }
/*     Print statistics. */
/* Writing concatenation */
    i__3[0] = 30, a__1[0] = " Problem name                 ";
    i__3[1] = 8, a__1[1] = mprob;
    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)38);
    snprnt_(&c__13, ch__1, &iw[1], leniw, (ftnlen)38);
    s_wsfi(&io___117);
    do_fio(&c__1, (char *)&itn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*objtru), (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
    if (nlninf) {
	s_wsfi(&io___118);
	do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	if (! lininf) {
	    s_wsfi(&io___119);
	    do_fio(&c__1, (char *)&wtinf, (ftnlen)sizeof(doublereal));
	    d__1 = fmrt / wtinf;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
    }
    if (gotfun) {
	if (nonlin) {
	    s_wsfi(&io___120);
	    do_fio(&c__1, (char *)&(*nmajor), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&flin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    s_wsfi(&io___121);
	    do_fio(&c__1, (char *)&pennrm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*fobj), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    s_wsfi(&io___122);
	    do_fio(&c__1, (char *)&iw[194], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iw[189], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	if (iw[70] < 3 || nonlin && lvlsch == 0) {
	    s_wsfi(&io___123);
	    do_fio(&c__1, (char *)&iw[195], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iw[190], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	if (iw[70] < 3) {
	    s_wsfi(&io___124);
	    do_fio(&c__1, (char *)&iw[196], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iw[191], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	    s_wsfi(&io___125);
	    do_fio(&c__1, (char *)&iw[197], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iw[192], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	if (*ns > 0) {
	    s_wsfi(&io___126);
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nnb, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	if (iw[386] > 0) {
	    s_wsfi(&io___127);
	    do_fio(&c__1, (char *)&iw[386], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
	}
	s_wsfi(&io___128);
	do_fio(&c__1, (char *)&ndegen, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&degen, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)133);
    }
/*     ------------------------------------------------------------------ */
/*     Unscale, compute nonlinear constraint violations, */
/*     save basis files and prepare to print the solution. */
/*     Clock 3 is "Output time". */
/*     ------------------------------------------------------------------ */
    s1time_(&c__3, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     Skip the functions if we don't have them. */
/*     Skip unscaling everything for infeasible linear constraints, they */
/*     have already been unscaled. */
    lssave = iw[75];
    if (! gotfun || fponly) {
	nncon1 = 0;
	nnobj1 = 0;
	if (lininf) {
	    iw[75] = 0;
	}
    } else {
	nncon1 = *nncon;
	nnobj1 = *nnobj;
    }
    s4savb_(iexit, &c__0, &minimz, m, n, nb, &nkx, &nncon0, &nncon1, &nnl0, &
	    nnobj1, nname, ns, &itn, ninf, sinf, &wtinf, &vimax, iobj, &
	    sclobj, objtru, &pnorm1, &pnorm2, &pinorm, &xnorm, ne, nlocj, &
	    locj[1], &indj[1], &jcol[1], &iw[lkx], &iw[lhesta], &hs[1], &rw[
	    lascal], &bl[1], &bu[1], &rw[lfx], &rw[lgobj], names + 8, &pi[1], 
	    &rc[1], &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)8, (ftnlen)8);
/*     If task = 'Print', s4savB prints the solution under the control */
/*     of lprSol (set by the  Solution  keyword in the SPECS file). */
/*     The printed solution may or may not be wanted, as follows: */

/*     lprSol = 0   means      No */
/*            = 2   means      Yes */
    s4savb_(iexit, &c__1, &minimz, m, n, nb, &nkx, &nncon0, &nncon1, &nnl0, 
	    nnobj, nname, ns, &itn, ninf, sinf, &wtinf, &vimax, iobj, &sclobj,
	     objtru, &pnorm1, &pnorm2, &pinorm, &xnorm, ne, nlocj, &locj[1], &
	    indj[1], &jcol[1], &iw[lkx], &iw[lhesta], &hs[1], &rw[lascal], &
	    bl[1], &bu[1], &rw[lfx], &rw[lgobj], names + 8, &pi[1], &rc[1], &
	    x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
	    ftnlen)8);
    iw[75] = lssave;
    s1time_(&c_n3, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     If the user hasn't already pulled the plug, */
/*     call the functions one last time with  Status .ge. 2. */
/*     Everything has been  unscaled, so we have to disable scaling. */
/*     modefg = 0  requests that no gradients are computed. */
/*     ------------------------------------------------------------------ */
    if (gotfun && *iexit / 10 != 7 && *iexit / 10 != 6) {
/* Computing MIN */
	i__1 = *iexit / 10;
	iw[236] = min(i__1,4) + 2;
	modefg = 0;
	lssave = iw[75];
	iw[75] = 0;
	(*fgwrap)(&inform__, &modefg, &nlncon, &nlnobj, n, &negcon, &nncon0, 
		nncon, nnjac, &nnl, &nnobj0, nnobj, (U_fp)funcon, (U_fp)
		funobj, &x[1], ne, nlocj, &locj[1], &indj[1], &rw[lfcon], 
		fobj, &rw[lgcon2], &rw[lgobj2], cu + 8, lencu, &iu[1], leniu, 
		&ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)8, (ftnlen)8);
	iw[75] = lssave;
    }
/* Save some things needed by solvers calling SNOPT */
/* L999: */
    rw[421] = *objtru;
/* The true objective */
    rw[422] = pinorm;
/* Lagrange multiplier norm */
    rw[423] = xnorm;
/* Norm of the variables (for GAMS) */
    rw[424] = wtinf;
/* Infeasibility weight */
    rw[433] = *sinf;
/* Sum of infeasibilities */
    rw[434] = flin;
/* Linear objective */
    rw[435] = *fobj;
/* Objective function */
    rw[436] = pennrm;
/* Norm of penalty parameters */
    iw[421] = itn;
/* Total iteration count */
    iw[422] = *nmajor;
/* Major iterations */
    iw[423] = maxs;
/* max # of superbasics */
    return 0;
} /* s8solv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8solv */
/* Subroutine */ int s8dflt_(integer *m, integer *n, integer *nncon, integer *
	nnjac, integer *nnobj, char *cw, integer *lencw, integer *iw, integer 
	*leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal c4, c6;
    static integer nnl;
    static doublereal eps, eps0, eps1, eps2, eps3, eps4;
    static integer npr1, npr2, kfac;
    static char mbnd[8];
    static integer kchk;
    static char mobj[8];
    static integer klog;
    static char mrng[8];
    static integer ksav, maxr, maxs;
    static char mrhs[8];
    static integer nout;
    static doublereal tolx, dens1, dens2, lmax1, lmax2, xpen0, utol1, utol2;
    static integer iback, ioldb;
    static doublereal bigdx, bigfx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal small, tolcg, xdlim;
    static char mprob[8];
    static integer idump, maxmn, never, mskip, isoln;
    static doublereal epsrf, tolfp, vilim;
    static integer ksumm;
    static doublereal tolqp, wtmax, fdint1, fdint2;
    static integer jverf1, jverf2, jverf3, jverf4, jverf5, jverf6;
    static doublereal toldj3, wtinf0, utol1m, hcndbd, utol2m;
    static integer iloadb, kdegen;
    static doublereal infbnd, zcndbd;
    static integer lemode;
    static doublereal chzbnd;
    static integer icrash;
    static logical linear;
    static integer lprdbg;
    static doublereal tolfac, uspace;
    static logical lincon;
    static integer maxcol;
    static doublereal tcrash;
    static integer mmajor;
    static doublereal toldcp;
    static logical nlncon;
    static integer lvlder, minmax, lvlinf, cgitmx, itnlim;
    static logical nonlin;
    static integer deropt, kreset, lprsch, lprscl, lvlhes, lvlpre, iprint, 
	    ireprt, lvlsch, iinsrt, lvlscl, lvlppm, lprsol, lprprm, lvlpiv, 
	    lvlver, mflush, minimz, minprc, mjrprt, mminor, mnrprt, mqnmod, 
	    mnewsb, nparpr, objrow, qpslvr, stkyop, tpivot, lvlsys;
    static doublereal scltol, tolcon, tolddp, toldpp, toldrp, toldup, tolnlp, 
	    tolpiv, tolrow, tolswp, tolupd, wolfeg;

/*     ================================================================== */
/*     s8dflt checks the optional parameter values and possibly changes */
/*     them to reasonable values. */

/*     Note that checking occurs before the amount of working storage has */
/*     been defined. */

/*     See  snworkspace.info  for full documentation of cw, iw and rw. */

/*     15 Nov 1991: first version. */
/*     27 Apr 2001: wtMax introduced. */
/*     10 Dec 2002: Added defaults for LU Rook and Diagonal Pivoting. */
/*     31 Dec 2002: Added default MPS character names. */
/*     30 Jul 2003: Added default CG tolerance. */
/*     22 Jun 2004: Added default LU mod singularity tol */
/*     20 Dec 2004: Default LU tols reduced. */
/*     21 Dec 2004: Default LU tols fixed up. */
/*     09 Jun 2005: Default tolpiv same as in s5dflt. */
/*     02 Jul 2005: Default Utol's set back to eps1. */
/*     02 May 2006: lvlTim removed. */
/*     01 Sep 2007: stkyOp added. */
/*     25 Nov 2007: Hessian options added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Set some local machine-dependent constants. */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    eps1 = rw[3];
/* eps**(2/3)          IEEE DP  3.67e-11 */
    eps2 = rw[4];
/* eps**(1/2)          IEEE DP  1.49e-08 */
    eps3 = rw[5];
/* eps**(1/3)          IEEE DP  6.05e-06 */
    eps4 = rw[6];
/*     ------------------------------------------------------------------ */
/*     rw(51)--rw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
/* eps**(1/4)          IEEE DP  1.22e-04 */
    tolfp = rw[51];
/* Minor Phase 1 Opt tol */
    tolqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolnlp = rw[53];
/* Major Optimality tolerance */
    tolcg = rw[54];
/* cg tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    tolcon = rw[57];
/* Major feasibility tolerance. */
    tolpiv = rw[60];
/* excludes small elements of y */
    tolrow = rw[61];
/* tolerance for the row error */
    tcrash = rw[62];
/* crash tolerance */
    utol1m = rw[63];
/* abs tol for small diag of U in LU mod */
    utol2m = rw[64];
/* rel tol for small diag of U in LU mod */
    tolswp = rw[65];
/* LU swap tolerance */
    tolfac = rw[66];
/* User-defined LU factor tolerance */
    tolupd = rw[67];
/* User-defined LU update tolerance */
    infbnd = rw[70];
/* definition of an infinite bound */
    bigfx = rw[71];
/* unbounded objective */
    bigdx = rw[72];
/* unbounded step */
    epsrf = rw[73];
/* relative function precision. */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    fdint2 = rw[77];
/* (2) cntrl  diff. interval */
    xdlim = rw[80];
/* Step limit */
    vilim = rw[81];
/* violation limit */
    etarg = rw[83];
/* Quasi-Newton QP rg tolerance */
    wolfeg = rw[84];
/* line search tolerance. */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    wtinf0 = rw[88];
/* infeasibility weight */
    xpen0 = rw[89];
/* initial penalty parameter. */
    wtmax = rw[90];
/* max     infeasibility weight */
    scltol = rw[92];
/*     ------------------------------------------------------------------ */
/*     rw(151)--rw(180) are parmLU parameters for LUSOL (some optional). */
/*     ------------------------------------------------------------------ */
/* scale tolerance. */
    small = rw[153];
/* defn of small real. */
    utol1 = rw[154];
/* abs tol for small diag of U. */
    utol2 = rw[155];
/* rel tol for small diag of U. */
    uspace = rw[156];
/* limit on waste space in U. */
    dens1 = rw[157];
/* switch to search maxcol columns and no rows. */
    dens2 = rw[158];
/*     ------------------------------------------------------------------ */
/*     rw(181)--rw(199) pass parameters into various routines. */
/*     ------------------------------------------------------------------ */
/*     toldj3    = rw(186) ! current optimality tol */
/*     ------------------------------------------------------------------ */
/*     iw(1)--iw(50): I/O file numbers and dimensions. */
/*     ------------------------------------------------------------------ */
/* switch to dense LU. */
    iprint = iw[12];
/*     ------------------------------------------------------------------ */
/*     iw(51)--iw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    qpslvr = iw[55];
/* 0(1) => QP(QN) QP solver */
    lemode = iw[56];
/* >0    => use elastic mode */
    kchk = iw[58];
/* check (row) frequency */
    kfac = iw[59];
/* factorization frequency */
    ksav = iw[60];
/* save basis map */
    klog = iw[61];
/* log/print frequency */
    ksumm = iw[62];
/* Summary print frequency */
    kdegen = iw[63];
/* max. expansions of featol */
    kreset = iw[64];
/* Hessian frequency */
    mflush = iw[66];
/* Hessian flush */
    mskip = iw[67];
/*     lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:basis:warm:hot start */
/* # largest value of nSkip */
    lvlder = iw[70];
/* = 0, 1 or 2, the derivative level */
    lvlsys = iw[71];
/* > 0   => print system info */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lvlinf = iw[73];
/* Elastic option */
    lvlscl = iw[75];
/* scale option */
    lvlsch = iw[76];
/* >0     => use derivatives in the line search */
    lvlpre = iw[77];
/* >0    => QN preconditioned CG */
    lvlver = iw[78];
/* Verify level */
    lvlppm = iw[79];
/* 1(2)-norm proximal point method for x0 */
    lvlpiv = iw[80];
/* 0(1 2 3) LU partial(rook complete diagonal) */
    lprprm = iw[81];
/* > 0    => parms are printed */
    lprsch = iw[82];
/* line search debug starting itn */
    lprscl = iw[83];
/* > 0    => print the scales */
    lprsol = iw[84];
/* > 0    => print the solution */
    lprdbg = iw[85];
/* > 0    => private debug print */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    icrash = iw[88];
/* Crash option */
    itnlim = iw[89];
/* limit on total iterations */
    mmajor = iw[90];
/* limit on major iterations */
    mminor = iw[91];
/* limit on minor iterations */
    mjrprt = iw[92];
/* Major print level */
    mnrprt = iw[93];
/* Minor print level */
    nparpr = iw[94];
/* # of partial pricing sections */
    mnewsb = iw[95];
/* maximum # of new superbasics per major */
    cgitmx = iw[97];
/* CG iteration limit */
    objrow = iw[103];
/* Objective row of user-defined F */
    deropt = iw[104];
/* 0, 1, 2 => derivative option */
    jverf1 = iw[110];
/* col # to start derivative checking */
    jverf2 = iw[111];
/* col # to stop  derivative checking */
    jverf3 = iw[112];
/* col # to start derivative checking */
    jverf4 = iw[113];
/* col # to stop  derivative checking */
    jverf5 = iw[114];
/* start col for Hessian checking */
    jverf6 = iw[115];
/* stop  col for Hessian checking */
    stkyop = iw[116];
/* > 0 => optional parameters are sticky */
    iback = iw[120];
/* backup file */
    idump = iw[121];
/* dump file */
    iloadb = iw[122];
/* load file */
    inewb = iw[124];
/* new basis file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    ipnch = iw[127];
/* punch file */
    ireprt = iw[130];
/* report file */
    isoln = iw[131];
/*     ------------------------------------------------------------------ */
/*     iw(151)--iw(180) are luparm parameters for LUSOL (some optional). */
/*     ------------------------------------------------------------------ */
/* solution file */
    maxcol = iw[153];
/*     ------------------------------------------------------------------ */
/*     Character  workspace. */
/*     cw(51)--cw(150): optional parameters */
/*     ------------------------------------------------------------------ */
/* lu1fac: max. # columns */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Problem name */
    s_copy(mobj, cw + 416, (ftnlen)8, (ftnlen)8);
/* Objective name */
    s_copy(mrhs, cw + 424, (ftnlen)8, (ftnlen)8);
/* rhs name */
    s_copy(mrng, cw + 432, (ftnlen)8, (ftnlen)8);
/* range name */
    s_copy(mbnd, cw + 440, (ftnlen)8, (ftnlen)8);
/*     ------------------------------------------------------------------ */
/* bounds name */
    c4 = max(1e-4,eps3);
    c6 = max(1e-6,eps2);
    never = 99999999;
/*     =============================================================== */
/*     Check the optional parameters. */
/*     =============================================================== */
    if (*nncon == 0) {
	*nnjac = 0;
    }
    if (*nnjac == 0) {
	*nncon = 0;
    }
    nnl = max(*nnjac,*nnobj);
    lincon = *nncon == 0;
    nlncon = *nncon > 0;
    linear = nnl == 0;
    nonlin = nnl > 0;
    if (iback == -11111) {
	iback = 0;
    }
    if (idump == -11111) {
	idump = 0;
    }
    if (iloadb == -11111) {
	iloadb = 0;
    }
    if (inewb == -11111) {
	inewb = 0;
    }
    if (iinsrt == -11111) {
	iinsrt = 0;
    }
    if (ioldb == -11111) {
	ioldb = 0;
    }
    if (ipnch == -11111) {
	ipnch = 0;
    }
    if (ireprt == -11111) {
	ireprt = 0;
    }
    if (isoln == -11111) {
	isoln = 0;
    }
/*     Set unspecified frequencies or silly values to defaults. */
    if (kchk == -11111) {
	kchk = 60;
    }
    if (kfac <= 0) {
	kfac = 100;
	if (nonlin) {
	    kfac = 50;
	}
    }
    if (klog == -11111) {
	klog = 100;
    }
    if (ksumm == -11111) {
	ksumm = 100;
    }
    if (ksav == -11111) {
	ksav = 100;
    }
    if (kdegen == -11111) {
	kdegen = 10000;
    }
    if (mflush == -11111) {
	mflush = 0;
    }
/*     Sometimes, frequency 0 means "almost never". */
    if (kchk <= 0) {
	kchk = never;
    }
    if (mflush <= 0) {
	mflush = never;
    }
    if (klog <= 0) {
	klog = never;
    }
    if (ksav <= 0) {
	ksav = never;
    }
    if (ksumm <= 0) {
	ksumm = never;
    }
    if (kdegen <= 0) {
	kdegen = never;
    }
    if (kreset <= 0) {
	kreset = never;
    }
    if (icrash < 0) {
	icrash = 3;
    }
    if (minmax == -11111) {
	minmax = 1;
    }
    if (minmax == -1) {
	minimz = -1;
    } else {
	minimz = 1;
    }
    if (objrow == 0) {
	minimz = 1;
	minmax = 0;
    } else if (objrow > 0 && minmax == 0) {
	minimz = 1;
	objrow = 0;
    }
    if (mjrprt == -11111) {
	mjrprt = 1;
    }
    if (mnrprt == -11111) {
	mnrprt = 1;
    }
/*     if (mMinor .lt. 0      ) mMinor = max( 1000,5*max( n,m ) ) */
    if (mminor < 0) {
	mminor = 500;
    }
    if (mmajor < 0) {
/* Computing MAX */
	i__1 = 1000, i__2 = max(*n,*m) * 3;
	mmajor = max(i__1,i__2);
    }
    if (mskip < 0 && lincon) {
	mskip = never;
    }
    if (mskip < 0 && nlncon) {
	mskip = 2;
    }
    if (mnewsb <= 0) {
	mnewsb = 99;
    }
    if (lprdbg < 0) {
	lprdbg = 0;
    }
    if (lprprm < 0) {
	lprprm = 1;
    }
    if (lprsch < 0) {
	lprsch = never;
    }
    if (lprscl < 0) {
	lprscl = 0;
    }
    if (lprsol < 0) {
	lprsol = 2;
    }
/*     lvlSrt is checked in s3argA or s3argB */
/*     if (lvlSrt .lt. 0      ) lvlSrt =  0 */
    if (deropt != -11111) {
	if (deropt < 0 || deropt > 1) {
	    deropt = 1;
	}
	if (deropt == 0) {
	    lvlder = 0;
	}
	if (deropt == 1) {
	    lvlder = 3;
	}
    } else if (lvlder != -11111) {
	if (lvlder < 0 || lvlder > 3) {
	    lvlder = 3;
	}
	deropt = 0;
	if (lvlder == 3) {
	    deropt = 1;
	}
    } else {
	deropt = 1;
	lvlder = 3;
    }
    if (lvlder == 2 && minmax == 0) {
	lvlder = 3;
    }
    if (lvlder == 0 && minmax == 0) {
	lvlder = 1;
    }
    if (lvlsys < 0) {
	lvlsys = 0;
    }
    if (lvlver == -11111) {
	lvlver = 0;
    }
    if (lvlver < 0) {
	lvlver = -1;
    }
    if (lvlver > 3) {
	lvlver = 0;
    }
    lvlinf = 2;
    if (lvlsch < 0) {
	if (deropt != 1) {
	    lvlsch = 0;
	}
	if (lvlder != 3) {
	    lvlsch = 0;
	}
	if (lvlsch < 0) {
	    lvlsch = 1;
	}
    }
    if (lvlppm < 0) {
	lvlppm = 1;
    }
    lemode = 1;
    if (stkyop < 0) {
	stkyop = 0;
    }
/*     Check superbasics limit maxS and Hessian dimension maxR. */
    if (nonlin) {
	if (maxr < 0) {
/* Computing MIN */
	    i__1 = 2000, i__2 = nnl + 1;
	    maxr = min(i__1,i__2);
	}
	if (maxs < 0) {
	    maxs = nnl + 1;
	}
/* Computing MAX */
	i__1 = min(maxr,*n);
	maxr = max(i__1,0);
/* Computing MAX */
	i__1 = min(maxs,*n);
	maxs = max(i__1,1);
    } else {
/* linear */
	if (maxs <= 0) {
	    maxs = 1;
	}
	if (maxr <= 0) {
	    maxr = 1;
	}
    }
    if (maxs < maxr) {
	maxr = maxs;
    }
    if (qpslvr < 0) {
	qpslvr = 0;
    }
    if (maxr == 0) {
	qpslvr = 1;
    }
    if (qpslvr == 2) {
	lvlpre = 0;
    }
    if (qpslvr == 0) {
	lvlpre = 0;
    }
    if (lvlpre > 1) {
	lvlpre = 1;
    }
    if (lvlpre < 0 && qpslvr == 1) {
	lvlpre = 0;
    }
    if (cgitmx < 0) {
	cgitmx = 100;
    }
    if (qpslvr == 1 || maxr < maxs) {
	if (lvlhes < 0) {
	    lvlhes = 0;
	}
	if (mqnmod < 0) {
	    mqnmod = 10;
	}
    } else {
	if (lvlhes < 0 && nnl > 75) {
	    lvlhes = 0;
	}
	if (lvlhes < 0 && nnl <= 75) {
	    lvlhes = 1;
	}
	if (lvlhes == 1) {
	    mqnmod = kreset;
	}
	if (mqnmod < 0) {
	    mqnmod = 10;
	}
    }
/*     --------------------------------- */
/*     CG QP optional parameters */
/*     --------------------------------- */
    if (etarg < 0. || etarg > 1.) {
	etarg = .1;
    }
/*     Check other options. */
    if (lvlscl < 0) {
	lvlscl = 2;
	if (nonlin) {
	    lvlscl = 1;
	}
    }
    lvlscl = min(lvlscl,2);
    if (lvlscl == 1 && nnl == *n) {
	lvlscl = 0;
    }
    if (nparpr <= 0) {
	nparpr = 10;
	if (nonlin) {
	    nparpr = 1;
	}
    }
    minprc = 10;
    npr1 = *n / nparpr;
    npr2 = *m / nparpr;
    if (max(npr1,npr2) < minprc) {
	maxmn = max(*m,*n);
	nparpr = maxmn / min(maxmn,minprc);
    }
/* Computing MAX */
    d__1 = 1. / (eps * 100.);
    chzbnd = max(d__1,1e6);
    if (infbnd < 0.) {
	infbnd = 1e20;
    }
    if (bigfx <= 0.) {
	bigfx = 1e15;
    }
    if (bigdx <= 0.) {
	bigdx = infbnd;
    }
    if (hcndbd <= 0.) {
	hcndbd = chzbnd;
    }
    if (tcrash < 0. || tcrash >= 1.) {
	tcrash = .1;
    }
    if (vilim <= 0.) {
	vilim = 1e6;
    }
    if (wolfeg < 0. || wolfeg > 1.) {
	wolfeg = .9;
    }
    if (wtmax < 0.) {
	wtmax = 1e10;
    }
    if (xdlim <= 0.) {
	xdlim = 2.;
    }
    if (xpen0 < 0.) {
	xpen0 = 0.;
    }
    if (zcndbd <= 0.) {
	if (qpslvr == 0) {
	    zcndbd = 1e4;
	} else {
	    zcndbd = 1e6;
	}
    }
/*     --------------------------------- */
/*     Set up the parameters for lu1fac. */
/*     --------------------------------- */
    if (maxcol < 0) {
	maxcol = 5;
    }
    nout = iprint;
    if (lvlsys == 0) {
	nout = 0;
    }
    if (lvlpiv <= 0) {
	lvlpiv = 0;
    }
    if (lvlpiv > 3) {
	lvlpiv = 0;
    }
    tpivot = lvlpiv;
    if (linear) {
	toldpp = 100.;
	toldrp = 10.;
	toldcp = 10.;
	tolddp = 10.;
	toldup = 10.;
    } else {
/* nonlinear */
	toldpp = 3.99;
	toldrp = 3.99;
	toldcp = 3.99;
	tolddp = 3.99;
	toldup = 3.99;
    }
    if (tolfac < 1.) {
	if (lvlpiv == 0) {
	    tolfac = toldpp;
	}
	if (lvlpiv == 1) {
	    tolfac = toldrp;
	}
	if (lvlpiv == 2) {
	    tolfac = toldcp;
	}
	if (lvlpiv == 3) {
	    tolfac = tolddp;
	}
    }
    if (tolupd < 1.) {
	tolupd = toldup;
    }
    lmax1 = tolfac;
    lmax2 = tolupd;
    if (utol1 <= 0.) {
	utol1 = eps1;
    }
    if (utol2 <= 0.) {
	utol2 = eps1;
    }
    if (utol1m <= 0.) {
	utol1m = eps1;
    }
    if (utol2m <= 0.) {
	utol2m = eps1;
    }
    if (dens2 < 0.) {
	dens2 = .6;
    }
    if (small <= 0.) {
	small = eps0;
    }
    if (uspace <= 0.) {
	uspace = 3.;
    }
    if (dens1 <= 0.) {
	dens1 = .3;
    }
/*     Set some SQP tolerances. */
/*     Set the minor and major optimality tolerances. */
/*     Solve the QP subproblems fairly accurately even if the */
/*     NLP Optimality Tolerance is big. */
    if (tolnlp <= 0.) {
	tolnlp = c6 * 2.;
	if (epsrf > 0.) {
/* Computing MAX */
	    d__1 = tolnlp, d__2 = sqrt(epsrf * 10.);
	    tolnlp = max(d__1,d__2);
	}
    }
    if (tolqp <= 0.) {
/* Computing MIN */
	d__1 = c6, d__2 = tolnlp / 2.;
	tolqp = min(d__1,d__2);
    }
    if (tolfp < 0.) {
	tolfp = c6;
    }
    if (tolcg <= 0.) {
	tolcg = .01;
    }
    if (tolrow <= 0.) {
	tolrow = c4;
    }
    if (tolswp <= 0.) {
	tolswp = eps4;
    }
    if (tolx <= 0.) {
	tolx = c6;
    }
    if (tolcon <= eps) {
	tolcon = c6;
    }
    toldj3 = tolqp;
    if (scltol <= 0.) {
	scltol = .9;
    }
    if (scltol >= 1.) {
	scltol = .99;
    }
    if (tolpiv <= 0.) {
	tolpiv = eps1;
    }
    if (lincon) {
	if (wtinf0 < 0.) {
	    wtinf0 = 1.;
	}
    } else {
	if (wtinf0 < 0.) {
	    wtinf0 = 1e4;
	}
    }
    if (epsrf <= 0.) {
	epsrf = eps0;
    }
    if (fdint1 <= 0.) {
	fdint1 = sqrt(epsrf);
    }
    if (fdint2 <= 0.) {
	fdint2 = pow_dd(&epsrf, &c_b129);
    }
/*     Check  START and STOP  column numbers for derivative checking. */
    if (jverf1 <= 0) {
	jverf1 = 1;
    }
    if (jverf2 < 0) {
	jverf2 = *n;
    }
    if (lvlver == 2 || lvlver == 0) {
	jverf2 = 0;
    }
    if (jverf3 <= 0) {
	jverf3 = 1;
    }
    if (jverf4 < 0) {
	jverf4 = *n;
    }
    if (lvlver == 1 || lvlver == 0) {
	jverf4 = 0;
    }
    if (jverf5 <= 0) {
	jverf5 = 1;
    }
    if (jverf6 < 0) {
	jverf6 = *n;
    }
    if (lvlver == 1 || lvlver == 0) {
	jverf6 = 0;
    }
    if (iback == inewb) {
	iback = 0;
    }
    if (itnlim < 0) {
/* Computing MAX */
	i__1 = 10000, i__2 = max(*n,*m) * 10;
	itnlim = max(i__1,i__2);
    }
/*     Set default names (they may be printed by the basis routines). */
    if (s_cmp(mprob, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mprob, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mobj, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mobj, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mrhs, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mrhs, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mrng, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mrng, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mbnd, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mbnd, "        ", (ftnlen)8, (ftnlen)8);
    }
/*     ------------------------------------------------------------------ */
/*     Done. */
/*     Re-assign the options to their respective work arrays. */
/*     ------------------------------------------------------------------ */
    rw[51] = tolfp;
    rw[52] = tolqp;
    rw[53] = tolnlp;
    rw[54] = tolcg;
    rw[56] = tolx;
    rw[57] = tolcon;
    rw[60] = tolpiv;
    rw[61] = tolrow;
    rw[62] = tcrash;
    rw[63] = utol1m;
    rw[64] = utol2m;
    rw[65] = tolswp;
    rw[66] = tolfac;
    rw[67] = tolupd;
    rw[70] = infbnd;
    rw[71] = bigfx;
    rw[72] = bigdx;
    rw[73] = epsrf;
    rw[76] = fdint1;
    rw[77] = fdint2;
    rw[80] = xdlim;
    rw[81] = vilim;
    rw[83] = etarg;
    rw[84] = wolfeg;
    rw[85] = hcndbd;
    rw[86] = zcndbd;
    rw[88] = wtinf0;
    rw[89] = xpen0;
    rw[90] = wtmax;
    rw[92] = scltol;
    rw[151] = lmax1;
    rw[152] = lmax2;
    rw[153] = small;
    rw[154] = utol1;
    rw[155] = utol2;
    rw[156] = uspace;
    rw[157] = dens1;
    rw[158] = dens2;
/*     Dependent parameters set in s8dflt. */
    rw[181] = toldpp;
    rw[182] = toldcp;
    rw[183] = toldup;
    rw[186] = toldj3;
    rw[187] = toldrp;
/*     Addresses for integer quantities. */
    iw[52] = maxr;
    iw[53] = maxs;
    iw[54] = mqnmod;
    iw[55] = qpslvr;
    iw[56] = lemode;
    iw[58] = kchk;
    iw[59] = kfac;
    iw[60] = ksav;
    iw[61] = klog;
    iw[62] = ksumm;
    iw[63] = kdegen;
    iw[64] = kreset;
    iw[66] = mflush;
    iw[67] = mskip;
/*     iw( 69)  =  lvlSrt */
    iw[70] = lvlder;
    iw[71] = lvlsys;
    iw[72] = lvlhes;
    iw[73] = lvlinf;
    iw[75] = lvlscl;
    iw[76] = lvlsch;
    iw[77] = lvlpre;
    iw[78] = lvlver;
    iw[79] = lvlppm;
    iw[80] = lvlpiv;
    iw[81] = lprprm;
    iw[82] = lprsch;
    iw[83] = lprscl;
    iw[84] = lprsol;
    iw[85] = lprdbg;
    iw[87] = minmax;
    iw[88] = icrash;
    iw[89] = itnlim;
    iw[90] = mmajor;
    iw[91] = mminor;
    iw[92] = mjrprt;
    iw[93] = mnrprt;
    iw[94] = nparpr;
    iw[95] = mnewsb;
    iw[97] = cgitmx;
    iw[103] = objrow;
    iw[104] = deropt;
    iw[110] = jverf1;
    iw[111] = jverf2;
    iw[112] = jverf3;
    iw[113] = jverf4;
    iw[114] = jverf5;
    iw[115] = jverf6;
    iw[116] = stkyop;
    iw[120] = iback;
    iw[121] = idump;
    iw[122] = iloadb;
    iw[124] = inewb;
    iw[125] = iinsrt;
    iw[126] = ioldb;
    iw[127] = ipnch;
    iw[130] = ireprt;
    iw[131] = isoln;
    iw[151] = nout;
    iw[153] = maxcol;
    iw[156] = tpivot;
/*     Dependent parameters set in s8dflt. */
    iw[199] = minimz;
/*     Addresses for character quantities. */
    s_copy(cw + 408, mprob, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 416, mobj, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 424, mrhs, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 432, mrng, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 440, mbnd, (ftnlen)8, (ftnlen)8);
    return 0;
} /* s8dflt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8dflt */
/* Subroutine */ int s8map_(integer *m, integer *n, integer *negcon, integer *
	nkx, integer *nncon, integer *nnjac, integer *nnobj, integer *lenr, 
	integer *maxr, integer *maxs, integer *mqnmod, integer *lvlhes, 
	integer *nextcw, integer *nextiw, integer *nextrw, integer *iw, 
	integer *leniw)
{
    static integer nb, lr, ls, lv, ly, lr1, lr2, ls1, ls2, ls3, lu0, lx1, ly1,
	     ly2, ly3, ldg, lhd, mbs, lfr, lrg, nnl, ldx, lfv, lfx, liy, lkx, 
	    lux, lrg2, liy1, liy2, lgbs, lkbs, lhdx, lpbs, lgqp, ngqp, lxbs, 
	    lxqp, lxqp0, lgobj, lblbs, llocg, lfcon, nlocg, lgcon, lbubs, 
	    lenfr, lycon, lxpen, lgobj1, lgobj2, lfcon1, lfcon2, lgcon1, 
	    lgcon2, lycon1, lycon2, lascal, lhfeas, lhesta, lgobju, lblsav, 
	    lxscal, lgconu, ldycon, lbusav, lhetyp, lqprhs;

/*     ================================================================== */
/*     s8Map   allocates all array storage for snopt, */
/*     using the values: */
/*        m    , n    , ne */
/*        maxS                    Set in s8dflt. */
/*        nnObj, nnCon, nnJac     Set by specs file or argument list. */
/*        lenR , negCon           Set in calling program */

/*     On exit, */
/*        nextcw, nextiw, nextrw are pointers to the next elements of */
/*                               free space in cw, iw, and rw. */

/*     29 Dec 2000: First version of s8Map. */
/*     24 Jan 2003: Added workspace for SYMMLQ */
/*     18 Jun 2008: Added space for iy2, pBS and rg2. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    nnl = max(*nnjac,*nnobj);
/*     All dimensions are computed from */
/*        m     , n    , ne */
/*        lenR  , maxS , mQMmod */
/*        nnObj , nnCon, nnJac */
/*        negCon */
    ngqp = nnl;
    mbs = *m + *maxs;
    nb = *n + *m;
/*     Nonlinear constraints. */
    nlocg = *nnjac + 1;
/*     Addresses for the integer arrays. */
    lkx = *nextiw;
    lhfeas = lkx + *nkx;
    lkbs = lhfeas + mbs;
    lhesta = lkbs + mbs;
    lhetyp = lhesta + nb;
    liy = lhetyp + nb;
    liy1 = liy + nb;
    liy2 = liy1 + nb;
    *nextiw = liy2 + nb;
/*     Addresses for the double precision arrays. */
    lascal = *nextrw;
    ly = lascal + nb;
    ly1 = ly + nb;
    ly2 = ly1 + nb;
    ly3 = ly2 + nb;
/* SYMMLQ workspace */
    ls1 = ly3 + nb;
/* SYMMLQ workspace */
    ls2 = ls1 + *maxs;
/* SYMMLQ workspace */
    ls3 = ls2 + *maxs;
/* SYMMLQ workspace */
    lr1 = ls3 + *maxs;
/* SYMMLQ workspace */
    lr2 = lr1 + *maxs;
/* SYMMLQ workspace */
    lblbs = lr2 + *maxs;
    lbubs = lblbs + mbs;
    lxbs = lbubs + mbs;
    lxscal = lxbs + mbs;
    lgbs = lxscal + nnl;
    lgqp = lgbs + mbs;
    lux = lgqp + ngqp;
    lhdx = lux + nnl;
    lpbs = lhdx + nnl;
    lu0 = lpbs + nb;
    ldg = lu0 + nnl;
    lr = ldg + nnl;
    lrg = lr + *lenr;
    lrg2 = lrg + *maxs;
    lblsav = lrg2 + *maxs;
    lbusav = lblsav + nb;
    *nextrw = lbusav + nb;
/*     Nonlinear Objective. */
    lgobj = *nextrw;
    lgobj1 = lgobj + *nnobj;
    lgobj2 = lgobj1 + *nnobj;
    lgobju = lgobj2 + *nnobj;
    *nextrw = lgobju + *nnobj;
/*     Nonlinear constraints. */
    llocg = *nextiw;
    *nextiw = llocg + nlocg;
    lfcon = *nextrw;
    lfcon1 = lfcon + *nncon;
    lfcon2 = lfcon1 + *nncon;
    lfx = lfcon2 + *nncon;
    lfv = lfx + *nncon;
    lycon = lfv + *nncon;
    lycon1 = lycon + *nncon;
    lycon2 = lycon1 + *nncon;
    ldycon = lycon2 + *nncon;
    lxpen = ldycon + *nncon;
    lgcon = lxpen + *nncon;
    lgcon1 = lgcon + *negcon;
    lgcon2 = lgcon1 + *negcon;
    lgconu = lgcon2 + *negcon;
    lqprhs = lgconu + *negcon;
    ldx = lqprhs + *m;
    lxqp = ldx + nb;
    lxqp0 = lxqp + nb;
    lx1 = lxqp0 + nb;
    *nextrw = lx1 + nb;
/*     Store the addresses in iw. */
    iw[251] = lkx;
    iw[260] = llocg;
    iw[273] = lblbs;
    iw[274] = lbubs;
    iw[275] = lblsav;
    iw[276] = lbusav;
    iw[277] = lpbs;
    iw[278] = lqprhs;
    iw[283] = lhetyp;
    iw[284] = lhfeas;
    iw[285] = lhesta;
    iw[287] = ldx;
    iw[288] = lhdx;
    iw[289] = ldg;
    iw[290] = lgqp;
    iw[291] = lgbs;
    iw[292] = lkbs;
    iw[293] = lrg;
    iw[294] = lrg2;
    iw[295] = lr;
    iw[296] = lascal;
    iw[297] = lgobj;
    iw[300] = lx1;
    iw[301] = lxbs;
    iw[302] = lxscal;
    iw[304] = lxpen;
    iw[305] = lxqp;
    iw[306] = lxqp0;
    iw[308] = liy;
    iw[309] = liy1;
    iw[310] = liy2;
    iw[311] = ly;
    iw[312] = ly1;
    iw[313] = ly2;
    iw[314] = ly3;
    iw[316] = lfcon;
    iw[317] = lfcon1;
    iw[318] = lfcon2;
    iw[319] = lgconu;
    iw[320] = lgcon;
    iw[321] = lgcon1;
    iw[322] = lgcon2;
    iw[323] = lgobju;
    iw[324] = lgobj1;
    iw[325] = lgobj2;
    iw[336] = lfx;
    iw[337] = lfv;
    iw[345] = lux;
    iw[346] = lu0;
    iw[348] = lycon;
    iw[349] = lycon1;
    iw[350] = lycon2;
    iw[351] = ldycon;
    iw[353] = lr1;
    iw[354] = lr2;
    iw[355] = ls1;
    iw[356] = ls2;
    iw[357] = ls3;
/*     Allocate space for an approximate Hessian. */
/*     The amount will depend on the method selected. */
    if (*lvlhes == 0) {
/*        --------------------------------------------------------------- */
/*        Compute the addresses of the limited-memory arrays. */
/*        These are saved and used for subsequent entries. */
/*        --------------------------------------------------------------- */
	lhd = *nextrw;
	ls = lhd + nnl;
	lv = ls + nnl * *mqnmod;
	*nextrw = lv + nnl * *mqnmod;
	iw[347] = lhd;
	iw[401] = ls;
	iw[402] = lv;
    } else if (*lvlhes == 1) {
	lenfr = nnl * (nnl + 1) / 2;
	lhd = *nextrw;
	lfr = lhd + nnl;
	*nextrw = lfr + lenfr;
	iw[347] = lhd;
	iw[391] = lfr;
	iw[392] = lenfr;
    }
    return 0;
} /* s8map_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Map */
/* Subroutine */ int s8sqp_(integer *iexit, S_fp fgwrap, U_fp funcon, U_fp 
	funobj, S_fp mjrlog, U_fp mnrlog, S_fp snstop, logical *gotr, integer 
	*itn, integer *lenr, integer *m, integer *maxs, integer *mbs, integer 
	*n, integer *nb, integer *ns, integer *nncon0, integer *nncon, 
	integer *nnobj0, integer *nnobj, integer *nnl0, integer *nnl, integer 
	*nmajor, integer *nminor, integer *ndegen, doublereal *duinf, integer 
	*minimz, integer *iobj, doublereal *sclobj, doublereal *objadd, 
	doublereal *fobj, doublereal *fmrt, doublereal *vimax, doublereal *
	virel, doublereal *visup, integer *ninf, doublereal *sinf, doublereal 
	*wtinf0, doublereal *wtinf, doublereal *pennrm, doublereal *pinorm, 
	doublereal *xnorm, integer *ne, integer *nlocj, integer *locj, 
	integer *indj, doublereal *jcol, integer *negcon, integer *nlocg, 
	integer *locg, integer *hetype, integer *hestat, integer *hfeas, 
	integer *hs, integer *kbs, doublereal *ascale, doublereal *bl, 
	doublereal *bu, doublereal *blbs, doublereal *bubs, doublereal *fv, 
	doublereal *fx, doublereal *fcon, doublereal *gcon, doublereal *gobj, 
	doublereal *fcon1, doublereal *gcon1, doublereal *gobj1, doublereal *
	fcon2, doublereal *gcon2, doublereal *gobj2, doublereal *gbs, 
	doublereal *gqp, doublereal *dycon, doublereal *dx, doublereal *dg, 
	doublereal *udx, doublereal *hdx, doublereal *pbs, doublereal *ycon, 
	doublereal *ycon1, doublereal *ycon2, doublereal *pi, doublereal *
	qprhs, doublereal *r__, doublereal *rc, doublereal *rg, doublereal *
	rg2, doublereal *x, doublereal *x1, doublereal *xbs, doublereal *xqp0,
	 doublereal *xqp, doublereal *xpen, integer *iy, integer *iy1, 
	doublereal *y, doublereal *y1, doublereal *y2, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Initialized data */

    static char line[4] = "----";
    static char msg[19*5] = "max step too small." "step too small.    " "no "
	    "minimizer.      " "too many functions." "uphill direction.  ";

    /* Format strings */
    static char fmt_1000[] = "(1x,29a4)";
    static char fmt_1010[] = "(\002 Start of major itn\002,i6)";
    static char fmt_3020[] = "(\002 Itn\002,i7,\002 -- Central differences i"
	    "nvoked.\002,\002 Small step length.\002)";
    static char fmt_1050[] = "(\002 Search exit\002,i3,\002 -- \002,a,\002  "
	    " Itn =\002,i7,\002  Dual Inf =\002,1p,e11.3)";

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j;
    static doublereal eps;
    extern /* Subroutine */ int s8h0_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static char str[132];
    extern /* Subroutine */ int s8fd_(integer *, integer *, integer *, 
	    integer *, integer *, logical *, logical *, logical *, logical *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    static doublereal u0ii, eps0, eps1, eps5;
    extern /* Subroutine */ int s8rc_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s8fx_(integer *, integer *, integer *, doublereal *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int s8hx_();
    static doublereal back;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical newb;
    static integer info[6];
    static logical newg;
    static integer imsg, jobj, klog, maxr, itqp;
    static doublereal gmrt, rnnl, step;
    static logical newx;
    static doublereal tolx, fobj1, fobj2;
    extern /* Subroutine */ int s6fdg_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    S_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static doublereal sinf1, sinf2;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static logical qppi0;
    static doublereal fmrt1, gmrt1, xpen0;
    extern /* Subroutine */ int s8hqn_(integer *, S_fp, U_fp, U_fp, logical *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, logical *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s8iqn_(integer *, integer *, U_fp, U_fp, U_fp, integer *, integer 
	    *, logical *, logical *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen), s8iqp_(integer *, integer *, U_fp, 
	    U_fp, U_fp, integer *, integer *, logical *, logical *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen), s8mrt_(integer *, doublereal *, doublereal *, doublereal 
	    *, logical *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , s8xhx_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *), 
	    dload_(integer *, doublereal *, doublereal *, integer *);
    static logical fdobj, badls, debug;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static logical fdcon;
    static integer nnjac;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical goodg;
    static doublereal wtfac;
    static logical usefd;
    static doublereal dxhdx, prinf, gnorm;
    static logical maxns;
    static integer lureq;
    static logical newlu;
    static integer isumm, ksumm, maxvi, nswap, nskip;
    static doublereal rviol, tolfp, tolqp, wtmax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2bfac_(integer *, integer *, logical *,
	     logical *, logical *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *), 
	    daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *), s1page_(integer *, integer *, integer *)
	    , s2aprd_(integer *, doublereal *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), s6rcnd_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal gnorm0;
    extern doublereal dnrm1s_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int s6srch_(integer *, S_fp, U_fp, U_fp, logical *
	    , logical *, logical *, logical *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s8sinf_(integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s8infs_(
	    logical *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s8gcpy_(integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *), 
	    s8winf_(integer *, logical *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), s2vmax_(integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    ), s2tols_(integer *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *);
    static doublereal utol1s, utol2s;
    extern /* Subroutine */ int s8hwrp_();
    static logical feasbl;
    extern /* Subroutine */ int s8step_(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    static integer modefg;
    extern /* Subroutine */ int s8sopt_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static logical feaslk, backtr, dufeas, elastc, needlu;
    static integer jduinf;
    static logical prfeas;
    static integer iabort;
    static logical ktcond[2], nlnobj, rowfea;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static logical boostd, centrl, incrun, maxits, nlncon, fponly, newtol, 
	    optiml, prtlog, usefls;
    static integer cditns, inform__, iprint;
    static logical frstqp;
    static integer jprinf, lprsch, lvlder;
    static logical restrt;
    static integer lvlhes, lvlpiv, lvlpre, lvlsch;
    static logical prtsum;
    static integer lvlsrt, minmax, mmajor, mstart, mjrprt, mnrprt;
    static logical tnystp;
    static integer ninfqp, nstart, qpslvr, rtrmod, typelu;
    static doublereal condhz, fobjqp, drzmax, drzmin, pendmp, penmax, phpmrt, 
	    sinfqp, sgnobj, steplm, stepmn, stepmx, tolqpk, tolnlp, tolcon, 
	    utolmn, weight, wolfeg, wtscal, xdnorm;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static logical nearopt;

    /* Fortran I/O blocks */
    static icilist io___397 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___399 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___400 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___401 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___458 = { 0, str, 0, fmt_3020, 132, 1 };
    static icilist io___461 = { 0, str, 0, fmt_1050, 132, 1 };
    static icilist io___463 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___464 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___465 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___466 = { 0, str, 0, fmt_1010, 132, 1 };


/*     ================================================================== */
/*     s8SQP  solves a nonlinear programming problem. */
/*     A basis is assumed to be specified by nS, hs, x and the */
/*     superbasic parts of kBS. */
/*     In particular, there must be nS values hs(j) = 2, and the */
/*     corresponding j's must be listed in kBS(m+1) thru kBS(m+ns). */
/*     The ordering in kBS(m+1:m+nS) matches the reduced Hessian R. */

/*     On entry, if there are nonlinear constraints, Fx contains */
/*     the true nonlinear slacks (i.e., constraint values) */
/*     Fx  =  fCon + (linear A)*x,   excluding slacks. */

/*     On exit, if  iExit .lt. 30  it is safe to save the final */
/*     basis files and print the solution.  Otherwise, a fatal error */
/*     condition exists and numerous items will be undefined. */
/*     The last basis map saved (if any) retains the only useful */
/*     information. */

/*     30 Dec 1991: First version based on npsol routine npcore. */
/*     23 Oct 1993: Proximal point FP added. */
/*     29 Oct 1993: Crash on LG rows moved outside s5QP. */
/*     24 Apr 1994: Nx columns no longer in Q. */
/*     26 May 1995: Column order of R defined by kBS. */
/*     04 Aug 1995: Limited memory update */
/*     11 Aug 1995: tolg changed from 0.1 to 1.0d-4. */
/*     09 Nov 1995: Updated multipliers used to define Lagrangian. */
/*     19 Dec 1995: Finite-differences added. */
/*     09 Oct 1996: First Min Sum version. */
/*     16 Jul 1997: First thread-safe version. */
/*     09 Jul 1998: Quasi-Newton updates implemented correctly. */
/*     24 Aug 1998: Fixed bug in s8x1 found by Alan Brown at Nag. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     16 Jan 1999: Name changed from s8core. */
/*     06 Apr 2001: For Hot Starts, don't mess with Utol. */
/*     27 Apr 2001: wtMax introduced as parameter to s8wInf. */
/*     15 Jan 2003: CG and QN  QP solvers added. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     04 Jul 2005: Switched to vanilla CG for QN with nS > maxR. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     18 Jun 2008: pBS added as argument. */
/*     01 Dec 2012: Moved initialization of heStat to s8solv. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* abs tol for small diag of U. */
/* rel tol for small diag of U. */
/* =1(2) forwd (cntrl) diffs */
/* number of Hx products */
/* Approximate Hessian type */
/* Current QP solver */
/* Current precon mode */
/* itns since last factorize */
/* >0 => Mnr heading for iPrint */
/* >0 => Mjr heading for iPrint */
/*     ------------------------------------------------------------------ */
/* >0 => Mjr heading for iSumm */
    /* Parameter adjustments */
    --r__;
    --qprhs;
    --pi;
    --rg2;
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --y2;
    --y1;
    --y;
    --iy1;
    --iy;
    --xqp;
    --xqp0;
    --x1;
    --x;
    --rc;
    --pbs;
    --dx;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --xpen;
    --ycon2;
    --ycon1;
    --ycon;
    --dycon;
    --fcon2;
    --fcon1;
    --fcon;
    --fx;
    --fv;
    --gobj2;
    --gobj1;
    --gobj;
    --hdx;
    --udx;
    --dg;
    --gqp;
    --jcol;
    --indj;
    --locj;
    --gcon2;
    --gcon1;
    --gcon;
    --locg;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    iprint = iw[12];
/* Print file */
    isumm = iw[13];
/* Summary file */
    nnjac = iw[21];
/* # nonlinear Jacobian variables */
    maxr = iw[52];
/* max columns of R. */
    qpslvr = iw[55];
/* = 0:1:2   => QPChol:CG:QN QP solver */
    klog = iw[61];
/* log/print frequency */
    ksumm = iw[62];
/* Summary print frequency */
    lvlsrt = iw[69];
/* = 0:1:2:3 => cold:warm:basis:hot start */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lvlsch = iw[76];
/* >0     => use derivatives in the line search */
    lvlpre = iw[77];
/* >0     => QN preconditioned CG */
    lvlpiv = iw[80];
/* 0/1 Threshold partial/complete pivoting: use */
    lprsch = iw[82];
/* line search debug starting itn */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    mmajor = iw[90];
/* limit on major iterations */
    mjrprt = iw[92];
/* Major print level */
    mnrprt = iw[93];
/* Minor print level */
    lvlder = iw[70];
/*     Constants */
/* = 0, 1, 2 or 3, the derivative level */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    eps1 = rw[3];
/* eps**(2/3)          IEEE DP  3.67e-11 */
    eps5 = rw[7];
/* eps**(1/5)          IEEE DP  7.40e-04 */
    tolfp = rw[51];
/* Minor Phase 1 Opt tol */
    tolqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolnlp = rw[53];
/* Major Optimality tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    tolcon = rw[57];
/* Major feasibility tolerance. */
    wolfeg = rw[84];
/* line search tolerance. */
    xpen0 = rw[89];
/* initial penalty parameter. */
    wtmax = rw[90];
/* max     infeasibility weight */
    nlncon = *nncon > 0;
    nlnobj = *nnobj > 0;
    fponly = minmax == 0;
    iw[188] = 0;
    iw[223] = 0;
    iw[224] = 0;
    iw[226] = 0;
    iw[208] = qpslvr;
/* Current QP solver */
    iw[209] = lvlpre;
/*     ------------------------------------------------------------------ */
/*     s8SQP  operates in either ``Normal'' or ``Elastic'' mode. */
/*     In elastic mode, the nonlinear slacks are allowed to be infeasible */
/*     while a weighted sum of the slack infeasibilities is minimized. */
/*     ------------------------------------------------------------------ */
/* Current precon mode */
    feaslk = TRUE_;
/*     Elastc =       FP  .and.  nlnCon */
    elastc = FALSE_;
    *ninf = 0;
    *sinf = 0.;
    sinf1 = 0.;
    sinf2 = 0.;
    *iexit = 0;
    lureq = 0;
    nskip = 0;
    nstart = 0;
    if (*nnl > 0) {
	mstart = 2;
    } else {
	mstart = 0;
    }
    rtrmod = 0;
    iload_(&c__6, &c__0, info, &c__1);
    info[0] = 1;
/* Suppresses first printing of 'n' */
    sgnobj = (doublereal) (*minimz);
    gnorm = 1.;
    if (*iobj > 0) {
	gnorm += *sclobj;
    }
    gnorm0 = 0.;
    if (nlnobj) {
/*        gNorm0 = dnormi( nnObj, gObj, 1 ) */
	gnorm0 = dnrm2_(nnobj, &gobj[1], &c__1);
    }
/*     Hwt    = 0.001d+0         ! Weight on the Hessian */
/*     if (gNorm0 .gt. zero) then */
/*     U0ii   = sgnObj*Hwt*gNorm0 */
/*      else */
/*     U0ii   = sgnObj */
/*     end if */
    if (*nnl == 0) {
	u0ii = 1.;
    } else {
	if (gnorm0 > 0.) {
	    rnnl = (doublereal) (*nnl);
	    u0ii = sqrt(gnorm0 / sqrt(rnnl));
	} else {
	    u0ii = 1.;
	}
    }
/* Computing MIN */
    d__1 = max(u0ii,.01);
    u0ii = min(d__1,10.);
    rviol = 0.;
    prinf = 0.;
    *duinf = 0.;
    *wtinf = *wtinf0;
    tolqpk = tolqp * 10.;
    gmrt = 0.;
    step = 0.;
    ktcond[0] = FALSE_;
    ktcond[1] = FALSE_;
    frstqp = TRUE_;
    qppi0 = FALSE_;
/*     QPpi0  = .true.       ! Use QP initial multipliers */
/* Use zero initial multipliers */
    condhz = 1.;
    fdobj = (lvlder == 0 || lvlder == 2) && *nnobj > 0;
    fdcon = (lvlder == 0 || lvlder == 1) && nnjac > 0;
    usefd = fdobj || fdcon;
    usefls = usefd || lvlsch == 0;
    if (mjrprt >= 10 || mnrprt >= 10) {
	prtlog = iprint > 0 && klog == 1;
	prtsum = isumm > 0 && ksumm == 1;
	if (prtlog) {
	    s_wsfi(&io___397);
	    for (j = 1; j <= 29; ++j) {
		do_fio(&c__1, line, (ftnlen)4);
	    }
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___399);
	    do_fio(&c__1, (char *)&(*nmajor), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
	if (prtsum && mnrprt >= 10) {
	    s_wsfi(&io___400);
	    for (j = 1; j <= 19; ++j) {
		do_fio(&c__1, line, (ftnlen)4);
	    }
	    e_wsfi();
	    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___401);
	    do_fio(&c__1, (char *)&(*nmajor), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    jobj = *n + *iobj;
    if (nlncon) {
/*        --------------------------------------------- */
/*        Initialize the penalty parameters. */
/*        Set an initial elastic weight. */
/*        --------------------------------------------- */
	incrun = TRUE_;
	pendmp = 1.;
	penmax = 1. / eps;
	*pennrm = xpen0;
	dload_(nncon, &xpen0, &xpen[1], &c__1);
    }
    if (*nnl > 0 && iw[202] == -1) {
/*        --------------------------------------------------------------- */
/*        The approximate Hessian needs to be initialized. */
/*        Use the identity matrix until something better comes along. */
/*        --------------------------------------------------------------- */
	s8h0_(&iw[202], nnl, &u0ii, &iw[1], leniw, &rw[1], lenrw);
    }
    if (*ns > maxr) {
	iw[208] = 1;
/* Use CG */
	iw[209] = 0;
/* with no preconditioning */
	*gotr = FALSE_;
    }
    dcopy_(nb, &x[1], &c__1, &xqp[1], &c__1);
    cditns = -1;
    newg = FALSE_;
/*     ======================Start of main loop========================== */
/*     Start of a Major Iteration. */
/*     ================================================================== */
/* +    do while (iExit .eq. 0) */
L100:
    if (*iexit == 0) {
	*nminor = 0;
/*        =============================================================== */
/*        Repeat                    (until an accurate gradient is found) */
L110:
	centrl = iw[181] == 2;
	if (newg) {
	    if (usefd) {
/*                 ------------------------------------------------------ */
/*                 Compute any missing derivatives. */
/*                 ------------------------------------------------------ */
		s6fdg_(iexit, n, negcon, nncon0, nncon, &nnjac, nnl, nnobj0, 
			nnobj, (S_fp)fgwrap, (U_fp)funcon, (U_fp)funobj, &bl[
			1], &bu[1], &x[1], ne, nlocj, &locj[1], &indj[1], &
			fcon[1], fobj, &gcon[1], &gobj[1], &y[1], cu + 8, 
			lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &
			iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		if (*iexit != 0) {
		    goto L100;
		}
/* Break */
	    }
/* useFD */
	    newg = FALSE_;
	}
	if (nlncon) {
/*              --------------------------------------------------------- */
/*              Load the scaled Jacobian in J. */
/*              Compute the QP right-hand side   QPrhs  =  Jx - fCon. */
/*              Find Fx the nonlinear constraint values. */
/*              --------------------------------------------------------- */
	    s8gcpy_(nncon, &nnjac, ne, nlocj, &locj[1], &indj[1], negcon, 
		    nlocg, &locg[1], &gcon[1], ne, nlocj, &locj[1], &jcol[1]);
	    dcopy_(nncon, &fcon[1], &c__1, &qprhs[1], &c__1);
	    s2aprd_(&c__0, &eps0, ne, nlocj, &locj[1], &indj[1], &jcol[1], &
		    c_b168, &x[1], &nnjac, &c_b169, &qprhs[1], nncon);
/*              --------------------------------------------------------- */
/*              s8sOpt  finds the nonlinear slacks  sN  that minimize the */
/*              merit function with  x(1:n)  and  yCon  held fixed. */
/*              The optimal slacks are loaded into  x(n+1:nb)  and the */
/*              violations are calculated: */
/*                    Fv = fCon  + A(linear)x - nonlinear slacks */
/*                       = Fx                 - sN */
/*              --------------------------------------------------------- */
	    gnorm = 1.;
	    if (*iobj > 0) {
		gnorm += *sclobj;
	    }
	    gnorm0 = 0.;
	    if (nlnobj) {
		gnorm0 = dnrm1s_(nnobj, &gobj[1], &c__1);
	    }
	    if (! elastc) {
		d__1 = gnorm + gnorm0;
		s8winf_(&c__0, &boostd, itn, &d__1, wtinf0, wtinf, &wtmax, &
			weight, &wtfac, &wtscal, &iw[1], leniw);
	    }
	    s8sopt_(n, nb, nncon, &hestat[1], pinorm, &eps0, wtinf, &bl[1], &
		    bu[1], &fv[1], &x[1], &ycon[1], &xpen[1], &fx[1]);
/*              call s8Fv */
/*    &            ( Elastc, n, nnCon, eps0, wtInf, */
/*    &              bl, bu, Fv, x, yCon, Fx ) */
	}
/*           ------------------------------------------------------------ */
/*           Prepare to (re-)solve the QP subproblem (possibly after the */
/*           elastic weight has been increased). */
/*           ------------------------------------------------------------ */
/*           Factorize the basis at x. */
/*           Compute xQP such that (J -I)*xQP = rhs. */
L300:
	if (frstqp) {
/*              --------------------------------------------------------- */
/*              First QP subproblem. */
/*              --------------------------------------------------------- */
	    needlu = TRUE_;
	    *gotr = FALSE_;
	    nswap = 0;
	    if (*ns == 0) {
		typelu = 0;
	    } else {
		typelu = 2;
	    }
	    utol1s = rw[154];
	    utol2s = rw[155];
/*              To avoid an unnecessarily ill-conditioned starting basis */
/*              for the first QP, use big singularity tols */
/*              (except if it's a Hot Start!). */
	    if (lvlsrt == 3) {
		utolmn = eps1;
	    } else {
		utolmn = eps5;
	    }
	    rw[154] = max(utol1s,utolmn);
	    rw[155] = max(utol2s,utolmn);
	} else {
/*              --------------------------------------------------------- */
/*              Subsequent factorizations. */
/*              --------------------------------------------------------- */
/*              For linearly constrained problems, the factors L, U and R */
/*              can be saved as long as a poor x does not force a */
/*              new factorization. (Even in this case, R can be saved if */
/*              there are no swaps.) */
	    needlu = nlncon;
	    typelu = 3;
/*              Reset the factor and update tols if they were changed */
/*              during the previous major iteration. */
	    s2tols_(&c__3, &newtol, itn, &iw[1], leniw, &rw[1], lenrw);
	}
	s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, &mjrprt, &
		lureq, m, mbs, n, nb, nnl, ns, &nswap, ne, nlocj, &locj[1], &
		indj[1], &jcol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &blbs[1], 
		&bubs[1], nncon0, nncon, &qprhs[1], &xqp[1], &xbs[1], &iy[1], 
		&iy1[1], &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    goto L100;
	}
/* Break */
	if (iw[208] == 0) {
	    *gotr = *gotr && ! newlu;
	}
	needlu = FALSE_;
	if (mjrprt >= 10) {
	    iw[224] = 1;
	}
	if (frstqp) {
	    iw[80] = lvlpiv;
/* Reset original TPP or TCP */
	    rw[154] = utol1s;
	    rw[155] = utol2s;
	}
/*           ------------------------------------------------------------ */
/*           Solve the QP subproblem to obtain kBS, xQP and pi. */
/*           The search direction will be dx = xQP - x. */
/*           Use x1 to store the first feasible point. */
/*           ------------------------------------------------------------ */
	if (*ns <= maxr) {
	    if (qpslvr == 0 || qpslvr == 2) {
		iw[208] = qpslvr;
	    }
	}
/* Unless CG was requested explicitly,  use maxS = maxR */
/* with Cholesky and QN. Then switch to CG if necessary. */
	inform__ = 0;
	if (iw[208] == 0) {
	    d__1 = *objadd + *fobj;
	    s8iqp_(&inform__, info, (U_fp)s8hwrp_, (U_fp)s8hx_, (U_fp)mnrlog, 
		    &iw[202], &iw[188], &elastc, gotr, itn, &itqp, lenr, m, &
		    maxr, mbs, n, nb, nncon0, nncon, nnobj0, nnobj, nnl0, nnl,
		     ns, ndegen, &mjrprt, &mnrprt, minimz, iobj, sclobj, &
		    d__1, &fobjqp, &tolfp, &tolqpk, &tolx, &ninfqp, &sinfqp, 
		    wtinf, &u0ii, pinorm, ne, nlocj, &locj[1], &indj[1], &
		    jcol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[1], &kbs[
		    1], &ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &gbs[
		    1], &gqp[1], &gobj[1], &hdx[1], &pbs[1], &pi[1], &r__[1], 
		    &rc[1], &rg[1], &qprhs[1], &x[1], &xbs[1], &xqp0[1], &xqp[
		    1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu,
		     &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	} else if (iw[208] == 2) {
	    d__1 = *objadd + *fobj;
	    s8iqn_(&inform__, info, (U_fp)s8hwrp_, (U_fp)s8hx_, (U_fp)mnrlog, 
		    &iw[202], &iw[188], &elastc, gotr, itn, &itqp, lenr, m, &
		    maxr, mbs, n, nb, nncon0, nncon, nnobj0, nnobj, nnl0, nnl,
		     ns, ndegen, &mjrprt, &mnrprt, minimz, iobj, &condhz, 
		    sclobj, &d__1, &fobjqp, &tolfp, &tolqpk, &tolx, &ninfqp, &
		    sinfqp, wtinf, &u0ii, pinorm, ne, nlocj, &locj[1], &indj[
		    1], &jcol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[1], &
		    kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &
		    gbs[1], &gqp[1], &gobj[1], &hdx[1], &pbs[1], &pi[1], &r__[
		    1], &rc[1], &rg[1], &rg2[1], &qprhs[1], &x[1], &xbs[1], &
		    xqp0[1], &xqp[1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], 
		    cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	}
	if (inform__ == -2 && maxr < *maxs) {
	    iw[208] = 1;
/* Switch to CG */
	    iw[209] = 0;
/* with no preconditioning */
	    *gotr = FALSE_;
	    info[4] = 0;
	}
	if (iw[208] == 1) {
	    d__1 = *objadd + *fobj;
	    s8iqn_(&inform__, info, (U_fp)s8hwrp_, (U_fp)s8hx_, (U_fp)mnrlog, 
		    &iw[202], &iw[188], &elastc, gotr, itn, &itqp, lenr, m, 
		    maxs, mbs, n, nb, nncon0, nncon, nnobj0, nnobj, nnl0, nnl,
		     ns, ndegen, &mjrprt, &mnrprt, minimz, iobj, &condhz, 
		    sclobj, &d__1, &fobjqp, &tolfp, &tolqpk, &tolx, &ninfqp, &
		    sinfqp, wtinf, &u0ii, pinorm, ne, nlocj, &locj[1], &indj[
		    1], &jcol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[1], &
		    kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &
		    gbs[1], &gqp[1], &gobj[1], &hdx[1], &pbs[1], &pi[1], &r__[
		    1], &rc[1], &rg[1], &rg2[1], &qprhs[1], &x[1], &xbs[1], &
		    xqp0[1], &xqp[1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], 
		    cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	}
/*           inform    Status */
/*           ------    ------ */
/*            >0       Fatal error */
/*             0       QP solution found */
/*            -1       Too many iterations */
/*            -2       Too many superbasics */
	*nminor += itqp;
	if (inform__ > 0) {
	    *iexit = inform__;
	    goto L100;
	}
/* Remember the QP inform until after the printing */
	maxns = inform__ == -2;
	maxits = inform__ == -1;
	if (frstqp) {
	    frstqp = FALSE_;
	}
	if (nlncon) {
	    if (qppi0) {
/* Initialize dyCon = 0, yCon = pi(qp) */
		dcopy_(nncon, &pi[1], &c__1, &ycon[1], &c__1);
		dload_(nncon, &c_b4, &dycon[1], &c__1);
	    } else {
/* Set        dyCon = yCon - pi(qp) */
		dcopy_(nncon, &pi[1], &c__1, &dycon[1], &c__1);
		daxpy_(nncon, &c_b169, &ycon[1], &c__1, &dycon[1], &c__1);
	    }
	    if (elastc && feaslk) {
/* Initialize the elastic slacks */
/* This only happens once. */
		dcopy_(nncon, &fx[1], &c__1, &x[*n + 1], &c__1);
		dcopy_(nncon, &pi[1], &c__1, &ycon[1], &c__1);
		dload_(nncon, &c_b4, &dycon[1], &c__1);
	    }
/* If yCon or x  changed, recompute Fv. */
	    if (qppi0 || elastc && feaslk) {
		s8sopt_(n, nb, nncon, &hestat[1], pinorm, &eps0, wtinf, &bl[1]
			, &bu[1], &fv[1], &x[1], &ycon[1], &xpen[1], &fx[1]);
/*                 call s8Fv  ( Elastc, n, nnCon, eps0, wtInf, */
/*    &                 bl, bu, Fv, x, yCon, Fx ) */
		if (elastc && feaslk) {
		    feaslk = FALSE_;
		}
		if (qppi0) {
		    qppi0 = FALSE_;
		}
	    }
/*              Find the sum of infeasibilities of the nonlinear slacks. */
	    s8sinf_(n, nb, nncon, &tolx, ninf, sinf, &bl[1], &bu[1], &x[1]);
	}
/*           Compute the search direction dx. */
	dcopy_(nb, &xqp[1], &c__1, &dx[1], &c__1);
	daxpy_(nb, &c_b169, &x[1], &c__1, &dx[1], &c__1);
	*xnorm = dnrm1s_(n, &x[1], &c__1);
	xdnorm = dnrm1s_(n, &dx[1], &c__1);
/*           Compute all the QP reduced costs. */
/*           (We could use yCon for the nonlinear pi's). */
/*           Compute the maximum dual infeasibility. */
	s8rc_(sclobj, minimz, iobj, m, n, nb, nnobj0, nnobj, nncon, &nnjac, 
		negcon, ne, nlocj, &locj[1], &indj[1], &jcol[1], &gobj[1], &
		gcon[1], &pi[1], &rc[1]);
	s8infs_(&elastc, n, nb, nncon0, nncon, &tolx, wtinf, &prinf, duinf, &
		jprinf, &jduinf, &bl[1], &bu[1], &fx[1], &rc[1], &x[1]);
/*           Compute the largest nonlinear row violation. */
	if (nlncon) {
	    rviol = dnormi_(nncon, &fv[1], &c__1);
	}
/*           ------------------------------------------------------------ */
/*           Test for convergence. */
/*           ------------------------------------------------------------ */
	if (*gotr) {
	    s6rcnd_(&maxr, ns, lenr, &r__[1], &drzmax, &drzmin, &condhz);
	}
	rviol /= *xnorm + 1.;
	prinf /= *xnorm + 1.;
	*duinf /= *pinorm + 1.;
	rowfea = rviol < tolcon && *ninf > 0;
	prfeas = prinf <= tolcon;
	dufeas = *duinf <= tolnlp;
	ktcond[0] = prfeas;
	ktcond[1] = dufeas;
	feasbl = prfeas || rowfea;
	optiml = dufeas && feasbl;
	nearopt = ! optiml && (ktcond[0] && *duinf < tolnlp * 10. || ktcond[1]
		 && prinf < tolcon * 10.);
	if (nlncon && optiml && *ninf > 0) {
	    d__1 = gnorm + gnorm0;
	    s8winf_(&c__1, &boostd, itn, &d__1, wtinf0, wtinf, &wtmax, &
		    weight, &wtfac, &wtscal, &iw[1], leniw);
	    if (boostd) {
		elastc = TRUE_;
		goto L300;
/* Solve the QP again */
	    }
	}
/*            if ( Elastc ) then */
/* !              tolQPk = tolQP */
/*               tolQpk = min(tolQPk, 0.1d+0/wtInf) */
/*            else */
	if (*nmajor == 0) {
	    tolqpk = min(*duinf,.001);
	}
	if (*nminor == 0 && tolqpk > tolqp) {
	    tolqpk *= .2;
	}
/* Computing MIN */
	d__1 = tolqpk * .5, d__2 = *duinf * .1;
	tolqpk = min(d__1,d__2);
	tolqpk = max(tolqpk,tolqp);
/*            end if */
/*           ------------------------------------------------------------ */
/*           Compute the current augmented Lagrangian merit function. */
/*           ObjAdd is added in the log routine. */
/*           ------------------------------------------------------------ */
	if (*iobj == 0) {
	    *fmrt = 0.;
	} else {
	    *fmrt = sgnobj * x[jobj] * *sclobj;
	}
	if (nlnobj) {
	    *fmrt += sgnobj * *fobj;
	}
	if (nlncon) {
	    dcopy_(nncon, &fv[1], &c__1, &y[1], &c__1);
	    ddscl_(nncon, &xpen[1], &c__1, &y[1], &c__1);
	    *fmrt = *fmrt - ddot_(nncon, &ycon[1], &c__1, &fv[1], &c__1) + 
		    ddot_(nncon, &y[1], &c__1, &fv[1], &c__1) * .5;
	    if (elastc) {
		*fmrt += *wtinf * *sinf;
	    }
	}
/*           ------------------------------------------------------------ */
/*           If the forward-difference estimate of the reduced gradient */
/*           of the Lagrangian is small,  prepare to: (i) switch to */
/*           central differences; (ii)  recompute the derivatives,  and */
/*           (iii) solve the QP again. */

/*           On the other hand, if central differences give a large */
/*           reduced-gradient norm, switch back to forward differences. */
/*           ------------------------------------------------------------ */
	s8fd_(nncon0, nncon, nnobj, itn, &cditns, &centrl, &goodg, &newg, &
		usefd, info, duinf, &fcon[1], fobj, &iw[1], leniw, &rw[1], 
		lenrw);
/*           ------------------------------------------------------------ */
/*           Print the details of this iteration. */
/*           ------------------------------------------------------------ */
	(*mjrlog)(&iabort, info, &iw[202], ktcond, &mjrprt, minimz, n, nb, 
		nncon0, ns, itn, nmajor, nminor, &nswap, &condhz, iobj, 
		sclobj, objadd, fmrt, pennrm, &step, &prinf, duinf, vimax, 
		virel, &hs[1], ne, nlocj, &locj[1], &indj[1], &jcol[1], &
		ascale[1], &bl[1], &bu[1], &fcon[1], &ycon[1], &x[1], cu + 8, 
		lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	(*snstop)(&iabort, info, &iw[202], ktcond, &mjrprt, minimz, m, maxs, 
		n, nb, nncon0, nncon, nnobj0, nnobj, ns, itn, nmajor, nminor, 
		&nswap, &condhz, iobj, sclobj, objadd, fmrt, pennrm, &step, &
		prinf, duinf, vimax, virel, &hs[1], ne, nlocj, &locj[1], &
		indj[1], &jcol[1], negcon, &ascale[1], &bl[1], &bu[1], &fcon[
		1], &gcon[1], &gobj[1], &ycon[1], &pi[1], &rc[1], &rg[1], &x[
		1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	if (iabort != 0) {
	    *iexit = 74;
/* User has aborted the run via snSTOP */
	    goto L100;
	}
	info[2] = 0;
/* +       until (.not. (useFD  .and.  .not.goodg)) */
	if (usefd && ! goodg) {
	    goto L110;
	}
/*        =============================================================== */
	if (prfeas && info[4] == 3) {
	    *iexit = 21;
	}
	if (maxits) {
	    *iexit = 31;
	}
	if (*nmajor >= mmajor) {
	    *iexit = 32;
	}
	if (maxns && *nmajor >= mmajor) {
	    *iexit = 33;
	}
	if (nearopt && *iexit != 0) {
	    *iexit = 3;
	}
	if (optiml) {
	    if (*ninf == 0) {
/* Optimal and feasible */
		if (fponly) {
		    *iexit = 2;
		} else {
		    *iexit = 1;
		}
	    } else {
/* Violations minimized */
		*iexit = 13;
	    }
	}
	if (*iexit != 0) {
	    goto L100;
	}
	step = 0.;
	nswap = 0;
/*        =============================================================== */
/*        Take a step in the right direction. */
/*        =============================================================== */
/*        Compute  dxHdx = s'Hs  and other directional derivatives. */
/*        Be prepared to fix up pHpMrt if there are linear variables. */
/*        --------------------------------------------------------------- */
	if (*nnl > 0) {
	    s8xhx_(nnl, &dx[1], &udx[1], &hdx[1], &dxhdx, &iw[1], leniw, &rw[
		    1], lenrw);
	} else {
	    dxhdx = 0.;
	}
	if (*nnl == *n) {
	    if (dxhdx == 0. && iw[202] != 2) {
		iw[202] = 2;
		s8h0_(&iw[202], nnl, &u0ii, &iw[1], leniw, &rw[1], lenrw);
		goto L100;
	    }
	    phpmrt = dxhdx;
	} else {
/* Computing MAX */
	    d__1 = eps1 * xdnorm * xdnorm;
	    phpmrt = max(d__1,dxhdx);
	}
/*        --------------------------------------------------------------- */
/*        Compute the contributions to the merit function and its */
/*        directional derivative from the nonlinear constraints. */
/*        The penalty parameters  xPen(j)  are increased if the */
/*        directional derivative is not sufficiently negative. */
/*        --------------------------------------------------------------- */
/*        First, compute the value and directional derivative of the */
/*        Lagrangian with respect to x and the multipliers. */
	if (*iobj == 0) {
	    *fmrt = 0.;
	    gmrt = 0.;
	} else {
	    *fmrt = sgnobj * x[jobj] * *sclobj;
	    gmrt = sgnobj * dx[jobj] * *sclobj;
	}
	if (nlnobj) {
	    *fmrt += sgnobj * *fobj;
	    gmrt += sgnobj * ddot_(nnobj, &gobj[1], &c__1, &dx[1], &c__1);
	}
	if (elastc) {
	    *fmrt += *sinf * *wtinf;
	    gmrt += (sinfqp - *sinf) * *wtinf;
	}
/*        --------------------------------------------------------------- */
/*        Compute the search direction for the multipliers and nonlinear */
/*        slacks, and the contributions to the merit function and its */
/*        directional derivative from the nonlinear constraints. */
/*        The penalty parameters  xPen(j)  are increased if the */
/*        directional derivative is not sufficiently negative. */
/*        --------------------------------------------------------------- */
	if (nlncon) {
	    *fmrt -= ddot_(nncon, &ycon[1], &c__1, &fv[1], &c__1);
	    gmrt += ddot_(nncon, &ycon[1], &c__1, &fv[1], &c__1);
	    gmrt -= ddot_(nncon, &dycon[1], &c__1, &fv[1], &c__1);
	    s8mrt_(nncon, fmrt, &gmrt, &phpmrt, &incrun, &pendmp, &penmax, 
		    pennrm, &fv[1], &xpen[1], &y[1], &rw[1], lenrw);
	}
/*        =============================================================== */
/*        Find  stepmn,  stepmx  and  step,  the maximum, minimum and */
/*        initial values for the line search step. */
/*        =============================================================== */
	s8step_(&centrl, &usefls, nb, nncon, nnobj, nmajor, &nskip, &step, &
		stepmn, &steplm, &stepmx, &eps0, &xdnorm, xnorm, &bl[1], &bu[
		1], &x[1], &dx[1], &iw[1], leniw, &rw[1], lenrw);
	debug = *nmajor >= lprsch;
	back = .1;
/*        =============================================================== */
/*        Prepare for the linesearch to find a better point */
/*           x1 = x + step*dx  and  yCon1 = yCon + step*dyCon. */
/*        where, on entry,  x1 = xQP and  yCon1 = pi. */

/*        fCon , gCon , gObj  and yCon  are defined at the current    x. */
/*        fCon1, gCon1, gObj1 and yCon1 are defined at the new point x1. */
/*        fCon2, gCon2, gObj2 and yCon2 are temporary work arrays. */

/*        s6srch returns the following values: */

/*        inform    Result */
/*        ------    ------ */
/*         >0      Fatal error */
/*          0      redo the search. */
/*         -1      The search is successful and step < stpmax. */
/*         -2      The search is successful and step = stpmax. */
/*         -3      A better point was found but no sufficient decrease. */
/*                 Most likely, the merit function is decreasing at the */
/*                 boundary, but there could be too many function calls. */
/*         -4      stpmax < tolabs (too small to do a search). */
/*         -5      step   < stepmn (lsrchq only -- maybe want to switch */
/*                 to central differences to get a better direction). */
/*         -6      No useful step. */
/*                 The interval of uncertainty is less than 2*tolabs. */
/*                 The minimizer is very close to step = zero */
/*                 or the gradients are not sufficiently accurate. */
/*         -7      there were too many function calls. */
/*         -8      the input parameters were bad */
/*                 (stpmax le toltny  or  oldg ge 0). */
/*        =============================================================== */
/*        x and sInf are saved in case we have to restart the search. */
/*        pBS  is used as x2 (the QP solution) in s6srch. */
/* backtracking factor */
L500:
	dcopy_(nb, &xqp[1], &c__1, &pbs[1], &c__1);
	if (nlncon) {
	    dcopy_(nncon, &pi[1], &c__1, &ycon2[1], &c__1);
	}
	if (elastc) {
	    sinf2 = sinfqp;
	}
	s6srch_(&inform__, (S_fp)fgwrap, (U_fp)funcon, (U_fp)funobj, &debug, &
		elastc, &usefls, &prfeas, iobj, sclobj, n, nb, nncon0, nncon, 
		&nnjac, nnl, nnobj0, nnobj, itn, &wolfeg, &sgnobj, &step, &
		stepmn, &stepmx, &xdnorm, xnorm, fmrt, &fmrt1, &gmrt, &gmrt1, 
		sinf, &sinf1, &sinf2, wtinf, ne, nlocj, &locj[1], &indj[1], &
		jcol[1], negcon, nlocg, &locg[1], &fobj1, &fcon1[1], &gcon1[1]
		, &gobj1[1], &fobj2, &fcon2[1], &gcon2[1], &gobj2[1], &dx[1], 
		&dycon[1], &x[1], &x1[1], &pbs[1], &ycon[1], &ycon1[1], &
		ycon2[1], &xpen[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[
		1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1]
		, lenrw, (ftnlen)8, (ftnlen)8);
	if (inform__ > 0) {
	    *iexit = inform__;
/* The user wants to stop */
	    goto L100;
	}
	restrt = nstart < mstart;
/* permission to restart */
	newx = step > 0.;
	backtr = inform__ == 0;
	badls = inform__ <= -4;
	tnystp = inform__ == -4 || inform__ == -6 || inform__ == -8;
	if (usefd && ! centrl) {
/*           ------------------------------------------------------------ */
/*           If the line search failed.  Switch to central differences */
/*           and solve the QP subproblem again. */
/*           ------------------------------------------------------------ */
	    if (badls) {
		cditns = 0;
		info[5] = 1;
		s_wsfi(&io___458);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
		iw[181] = 2;
		newg = TRUE_;
		goto L110;
/* Recompute the derivatives */
	    }
	}
	if (newx || backtr) {
/*           ------------------------------------------------------------ */
/*           See if the search needs to be redone with a smaller stepmx. */
/*           ------------------------------------------------------------ */
	    if (backtr) {
		info[2] = 2;
/* User rejected the search. */
	    } else if (nlncon) {
/* Check vimax ge viSup */
		s8fx_(n, nncon, &nnjac, &eps0, ne, nlocj, &locj[1], &indj[1], 
			&jcol[1], &fcon1[1], &x1[1], &fx[1]);
		s2vmax_(n, nncon, &maxvi, vimax, &bl[1], &bu[1], &fx[1]);
		*virel = *vimax / (*xnorm + 1.);
		if (*vimax > *visup) {
		    backtr = TRUE_;
		    info[2] = 1;
		}
	    }
	    if (backtr) {
		stepmx = back * step;
		step = stepmx;
		back *= back;
		goto L500;
/* Repeat the line search */
	    }
	}
	if (info[2] == 1) {
/*           ------------------------------------------------------------ */
/*           The line search is backing away from a violation limit. */
/*           If we are in elastc mode or the line search died, switch to */
/*           elastic mode with a bigger infeasibility weight. */
/*           ------------------------------------------------------------ */
	    if (nlncon && (*ninf > 0 || tnystp)) {
		d__1 = gnorm + gnorm0;
		s8winf_(&c__1, &boostd, itn, &d__1, wtinf0, wtinf, &wtmax, &
			weight, &wtfac, &wtscal, &iw[1], leniw);
		if (boostd) {
		    elastc = TRUE_;
		}
		if (tnystp) {
		    goto L300;
		}
/* Solve the QP again */
	    }
	}
	if (badls && ! newx) {
/*           ============================================================ */
/*           Deal with some obvious cases. */
/*           ============================================================ */
	    if (maxns) {
		*iexit = 33;
/* Superbasics limit */
	    } else if (nearopt) {
		*iexit = 3;
/* Requested accuracy could not be ... */
	    } else if (info[2] == 1) {
		*iexit = 22;
	    } else if (info[2] == 2) {
		*iexit = 63;
	    }
	    if (*iexit > 0) {
		goto L100;
	    }
	}
	if (badls) {
/*           ------------------------------------------------------------ */
/*           The line search failed to provide a sufficient decrease. */
/*           ------------------------------------------------------------ */
	    if (newx) {
/*              Relax. At least we got SOME decrease. */
	    } else {
/*              Desperate times. */
/*              If possible, reset everything and solve the QP again. */
		if (restrt) {
		    if (iw[202] != 2) {
			iw[202] = 2;
			s8h0_(&iw[202], nnl, &u0ii, &iw[1], leniw, &rw[1], 
				lenrw);
		    }
		    if (iw[215] > 0) {
/*                    --------------------------------------------------- */
/*                    Try and fix up the basis. */
/*                    --------------------------------------------------- */
			s2tols_(&c__2, &newtol, itn, &iw[1], leniw, &rw[1], 
				lenrw);
		    }
		    if (nlncon) {
			incrun = TRUE_;
			pendmp = 1.;
			penmax = 1. / eps;
			*pennrm = xpen0;
			dload_(nncon, &xpen0, &xpen[1], &c__1);
		    }
		    ++nstart;
		    frstqp = TRUE_;
		    goto L300;
/* Solve the QP again */
		} else {
/*                 ------------------------------------------------------ */
/*                 We have run out of things to try. Bummer. */
/*                 ------------------------------------------------------ */
		    imsg = -inform__;
		    s_wsfi(&io___461);
		    do_fio(&c__1, (char *)&imsg, (ftnlen)sizeof(integer));
		    do_fio(&c__1, msg + (imsg - 4) * 19, (ftnlen)19);
		    do_fio(&c__1, (char *)&(*nmajor), (ftnlen)sizeof(integer))
			    ;
		    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(
			    doublereal));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
		    *iexit = 41;
/* Current point cannot be improved... */
		    goto L100;
		}
	    }
	}
/*        =============================================================== */
/*        The new point  x1  has been computed. */
/*        =============================================================== */
	if (step >= steplm && *nnl > 0) {
	    info[2] = 3;
	}
	inform__ = 0;
	centrl = iw[181] == 2;
/*        --------------------------------------------------------------- */
/*        Some unknown derivatives may need to be calculated at x1. */
/*        --------------------------------------------------------------- */
	if (usefls && *nnl > 0) {
	    modefg = 1;
	    (*fgwrap)(iexit, &modefg, &nlncon, &nlnobj, n, negcon, nncon0, 
		    nncon, &nnjac, nnl, nnobj0, nnobj, (U_fp)funcon, (U_fp)
		    funobj, &x1[1], ne, nlocj, &locj[1], &indj[1], &fcon2[1], 
		    &fobj2, &gcon1[1], &gobj1[1], cu + 8, lencu, &iu[1], 
		    leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1]
		    , lenrw, (ftnlen)8, (ftnlen)8);
	    if (*iexit != 0) {
		goto L100;
	    }
/* Break */
	    if (usefd) {
		s6fdg_(iexit, n, negcon, nncon0, nncon, &nnjac, nnl, nnobj0, 
			nnobj, (S_fp)fgwrap, (U_fp)funcon, (U_fp)funobj, &bl[
			1], &bu[1], &x1[1], ne, nlocj, &locj[1], &indj[1], &
			fcon1[1], &fobj1, &gcon1[1], &gobj1[1], &y[1], cu + 8,
			 lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &
			iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		if (*iexit != 0) {
		    goto L100;
		}
/* Break */
	    }
	}
	inform__ = 0;
	++(*nmajor);
	if (centrl) {
	    ++cditns;
	}
	if (mjrprt >= 10 || mnrprt >= 10) {
	    prtlog = iprint > 0 && klog == 1;
	    prtsum = isumm > 0 && ksumm == 1;
	    if (prtlog) {
		s1page_(&c__0, &iw[1], leniw);
		s_wsfi(&io___463);
		for (j = 1; j <= 29; ++j) {
		    do_fio(&c__1, line, (ftnlen)4);
		}
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
		s_wsfi(&io___464);
		do_fio(&c__1, (char *)&(*nmajor), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	    }
	    if (prtsum && mnrprt >= 10) {
		s_wsfi(&io___465);
		for (j = 1; j <= 19; ++j) {
		    do_fio(&c__1, line, (ftnlen)4);
		}
		e_wsfi();
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
		s_wsfi(&io___466);
		do_fio(&c__1, (char *)&(*nmajor), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)132);
	    }
	}
/*        =============================================================== */
/*        The problem functions have been defined at the new x. */
/*        =============================================================== */
	if (*nnl > 0 && (lvlhes == 0 || lvlhes == 1)) {
/* Update a QN approximate Hessian. */
	    s8hqn_(&inform__, (S_fp)fgwrap, (U_fp)funcon, (U_fp)funobj, &
		    usefd, &iw[202], &iw[208], info, lenr, m, mbs, n, nb, 
		    nncon0, nncon, &nnjac, nnl, nnobj0, nnobj, ns, nmajor, &
		    nskip, &u0ii, &step, minimz, &dxhdx, &rtrmod, gotr, &
		    incrun, &pendmp, &penmax, fobj, &fcon[1], &gcon[1], &gobj[
		    1], &fcon1[1], &gcon1[1], &gobj1[1], ne, nlocj, &locj[1], 
		    &indj[1], &jcol[1], negcon, nlocg, &locg[1], &kbs[1], &bl[
		    1], &bu[1], &dx[1], &dg[1], &udx[1], &hdx[1], &ycon1[1], &
		    r__[1], &x[1], &x1[1], &xqp0[1], &xpen[1], &y[1], &y1[1], 
		    &y2[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 
		    8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		    ftnlen)8);
	    if (inform__ != 0) {
		*iexit = inform__;
		goto L100;
	    }
	}
/*        --------------------------------------------------------------- */
/*        Update the variables. */
/*        The QP solution, saved in xQP, is used to start the next QP. */
/*        (If a unit step was not taken last iteration, some more */
/*        nonbasics may be between their bounds. */
/*        Nov 10, 1994. Tried leaving the nonbasics between their */
/*        bounds after short step. In some cases, the number of minor */
/*        iterations increased dramatically with a very short step.) */
/*        --------------------------------------------------------------- */
	dcopy_(nb, &x1[1], &c__1, &x[1], &c__1);
	if (nlncon) {
	    dcopy_(negcon, &gcon1[1], &c__1, &gcon[1], &c__1);
	    dcopy_(nncon, &ycon1[1], &c__1, &ycon[1], &c__1);
	    dcopy_(nncon, &fcon1[1], &c__1, &fcon[1], &c__1);
	}
	if (nlnobj) {
	    *fobj = fobj1;
	    dcopy_(nnobj, &gobj1[1], &c__1, &gobj[1], &c__1);
	}
	*sinf = sinf1;
	*ninf = ninfqp;
/* Not updated by the line search */
	goto L100;
/* +    end while */
    }
/*     ======================end of main loop============================ */
    return 0;
} /* s8sqp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8SQP */
/* Subroutine */ int s8stat_(integer *status, integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 XXX  user-function call-status not recog"
	    "nized.\002,\002 Requested status =\002,i6)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___468 = { 0, str, 0, fmt_9999, 80, 1 };


/*     ================================================================== */
/*     s8Stat fetches the call-status for the snOpt user-defined */
/*     functions. */

/*     16 Jun 2008: First version of s8Stat. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* NP user-routine call-status */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (iw[236] == 0) {
/* Standard call */
	*status = 0;
    } else if (iw[236] < 0) {
/* First call */
	*status = 1;
	iw[236] = 0;
    } else if (iw[236] >= 2) {
/* Last orders please */
	*status = iw[236];
	iw[236] = -1;
    } else {
	*status = iw[236];
	s_wsfi(&io___468);
	do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s8stat_ */

