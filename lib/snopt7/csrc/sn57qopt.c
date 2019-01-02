/* ../snopt7/src/sn57qopt.f -- translated by f2c (version 20100827).
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

static doublereal c_b2 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n2 = -2;
static integer c__13 = 13;
static integer c__3 = 3;
static integer c_n3 = -3;
static integer c__22 = 22;
static integer c__26 = 26;
static integer c__23 = 23;
static integer c__21 = 21;
static integer c__25 = 25;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn57qopt.f */

/*     s5solv */
/*     s5dflt   s5Map   s5sLP    s5sQP   s5sQN   s5Stat */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5solv_(integer *iexit, char *solver, integer *start, 
	U_fp hprod, U_fp hprod1, U_fp qplog, logical *gotr, integer *m, 
	integer *n, integer *nb, integer *nnh0, integer *nnh, integer *nname, 
	integer *ngqp, integer *ngobj0, integer *ngobj, integer *iobj, 
	doublereal *objadd, doublereal *objqp, doublereal *objtru, integer *
	ninf, doublereal *sinf, integer *ne, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, doublereal *bl, doublereal *bu, 
	doublereal *gobj, char *names, integer *nrhs0, integer *nrhs, 
	doublereal *rhs, integer *lenx0, integer *nx0, doublereal *x0, 
	integer *hetype, integer *hs, doublereal *x, doublereal *pi, 
	doublereal *rc, integer *ns, char *cu, integer *lencu, integer *iu, 
	integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen solver_len, ftnlen names_len, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1900[] = "(\002 No. of iterations\002,i20,2x,\002 Object"
	    "ive value\002,1p,e22.10)";
    static char fmt_1910[] = "(\002 No. of infeasibilities\002,i15,2x,\002 S"
	    "um of infeas\002,1p,e24.10)";
    static char fmt_1920[] = "(\002 No. of Hessian products\002,i14,2x,\002 "
	    "Objective row\002,3x,1p,e21.10)";
    static char fmt_1930[] = "(40x,\002 Quadratic objective\002,1p,e18.10)";
    static char fmt_1970[] = "(\002 No. of superbasics\002,i19,2x,\002 No. o"
	    "f basic nonlinears\002,i14)";
    static char fmt_1975[] = "(\002 No. of degenerate steps\002,i14,2x,\002 "
	    "Percentage\002,f27.2)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    char ch__1[38];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, lr;
    static doublereal fx[1];
    static integer ly, ly1, ly2, nnb, mbs, lrg, itn, liy, lkx, nkx;
    static char str[132];
    static integer lrg2, liy1, liy2;
    static char str2[132];
    static integer lgbs, lkbs, lenr, lhdx, lpbs, lgqp, prob, maxr, lxbs, maxs,
	     itqp;
    static doublereal tolx;
    static integer ngqp0;
    extern /* Subroutine */ int s5slp_(integer *, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen), s5sqn_(integer *, U_fp, U_fp, U_fp, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen), s5sqp_(integer *, U_fp, U_fp, U_fp,
	     logical *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static logical needb;
    static doublereal degen;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *);
    static integer nnjac;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dddiv_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer lblbs, inewb, lbubs, nnobj;
    static doublereal objlp;
    static char mprob[8];
    static integer nncon, numlc;
    static doublereal tolfp, vimax;
    static logical useqp;
    static doublereal tolqp, xnorm;
    extern /* Subroutine */ int s2amat_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *), s5getb_(integer *, integer *,
	     U_fp, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    char *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static integer nncon0;
    extern /* Subroutine */ int s1time_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *), s4newb_(integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, char *, integer *, integer *, integer *, 
	    ftnlen, ftnlen), s4savb_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s5qpfg_(U_fp, U_fp, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal wtinf0, pnorm1, pnorm2;
    extern /* Subroutine */ int s5fixs_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), s4stat_(integer *, char *, ftnlen);
    static integer lascal, ndegen, lhfeas, lemode;
    static doublereal sclobj;
    static logical infsbl;
    static integer lhesta, lblsav;
    static char istate[4*3];
    static integer inform__, itnlim, lbusav, lssave, lvlinf, minimz, minmax, 
	    numliq;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static doublereal pinorm, rgnorm;
    static integer mjrprt, sqstat, mnrprt;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), snprnt_(integer *, 
	    char *, integer *, integer *, ftnlen);
    static integer qpslvr, itqpmax;

    /* Fortran I/O blocks */
    static icilist io___71 = { 0, str, 0, fmt_1900, 132, 1 };
    static icilist io___72 = { 0, str, 0, fmt_1910, 132, 1 };
    static icilist io___73 = { 0, str, 0, fmt_1920, 132, 1 };
    static icilist io___74 = { 0, str, 0, fmt_1930, 132, 1 };
    static icilist io___75 = { 0, str, 0, fmt_1970, 132, 1 };
    static icilist io___76 = { 0, str, 0, fmt_1975, 132, 1 };


/*     ================================================================== */
/*     s5solv solves the current problem. */

/*     On entry, */
/*     the SPECS file has been read, */
/*     all data items have been loaded (including Acol, indA, locA, ...), */
/*     and workspace has been allocated within cw, iw and rw. */
/*     Start = lvlSrt from s3argQ. */

/*     On exit, */
/*     iExit  =  0 if an optimal solution was found, */
/*            =  1 if the problem was infeasible, */
/*            =  2 if the problem was unbounded, */
/*            =  3 if the Iteration limit was exceeded, */
/*           ge  4 if iterations were terminated by some other */
/*                 error condition (see the SQOPT user's guide). */

/*     01 Oct 1994: First version of s5solv. */
/*     06 Aug 1996: Min Sum option added. */
/*     14 Jul 1997: Thread-safe version. */
/*     02 Aug 2003: snEXIT and snPRNT adopted. */
/*     08 Mar 2004: Hot starts implemented. */
/*     16 May 2006: Explicit target itQP added. */
/*     18 Jun 2008: Added space for iy2, pBS and rg2. */
/*     07 Mar 2013: MnrPrt changed to MjrPrt in call to s2Amat. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* scale option */
/* number of Hx products */
/* type of QP Hessian */
/* Current QP solver */
/* Current precon mode */
/* # of LU factorizations */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Mnr heading for iPrint */
/* >0 => Minor heading for iSumm */
/* Save the LU factors */
/* Save the reduced Hessian */
/* QP user-routine call-status */
/* Number of symmlq iterations */
/* symmlq itns for last minor */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --x;
    --hs;
    --hetype;
    --bu;
    --bl;
    names -= 8;
    --gobj;
    --acol;
    --inda;
    --loca;
    --rhs;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    tolfp = rw[51];
/* Minor Phase 1 Opt tol */
    tolqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    wtinf0 = rw[88];
/* infeasibility weight */
    nnobj = iw[22];
/* # of objective variables */
    lenr = iw[28];
/* R(lenR) is the reduced Hessian factor */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    qpslvr = iw[55];
/* = 0:1:2   => QPChol:CG:QN QP solver */
    lemode = iw[56];
/* >0    => use elastic mode */
    lvlinf = iw[73];
/* Elastic option */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    itnlim = iw[89];
/* limit on total iterations */
    mjrprt = iw[92];
/* Major print level */
    mnrprt = iw[93];
/* Minor print level */
    inewb = iw[124];
/* new basis file */
    minimz = iw[199];
/* 1 (-1)    => minimize (maximize) */
    nkx = iw[247];
/* dimension of kx and its inverse, kxN */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Addresses */
/* Problem name */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    lhfeas = iw[284];
/* hfeas(mBS)  = feasibility types */
    lhesta = iw[285];
/* hEstat(nb)  = status of elastics */
    lkbs = iw[292];
/* kBS(mBS)    = ( B  S ) list */
    lblbs = iw[273];
/* blBS(mBS)   = lower bounds for xBS */
    lbubs = iw[274];
/* buBS(mBS)   = upper bounds for xBS */
    lpbs = iw[277];
/* pBS(nb)     = search direction */
    lxbs = iw[301];
/* xBS(mBS)    = basics, superbasics */
    lgqp = iw[290];
/* gQP(ngQP)   = QP gradient */
    lgbs = iw[291];
/* gBS(mBS)    = BS components of g */
    lrg = iw[293];
/* rg (maxS)   = reduced gradient */
    lrg2 = iw[294];
/* rg2(maxS)   = reduced gradient copy */
    lr = iw[295];
/* R(lenR)     = factor of Z'HZ */
    lascal = iw[296];
/* Ascale(nb)  = row and column scales */
    liy = iw[308];
/* iy (nb)     = integer work vector */
    liy1 = iw[309];
/* iy1(nb)     = integer work vector */
    liy2 = iw[310];
/* iy2(nb)     = integer work vector */
    ly = iw[311];
/*  y (nb)     = real work vector */
    ly1 = iw[312];
/*  y1(nb)     = real work vector */
    ly2 = iw[313];
/*  y2(nb)     = real work vector */
    lhdx = iw[288];
/* Hdx(nnH)    = product of H with  x - x0 */
    lblsav = iw[275];
/* blSav(m)    = temp bounds */
    lbusav = iw[276];
/* buSav(m)    = temp bounds */
    *iexit = 0;
    mbs = *m + maxs;
/* Figure out what type of problem we have. */
    if (minmax == 0 || lemode == 2 && lvlinf == 2) {
	prob = 0;
    } else if (*ngqp == 0) {
/* No explicit objective. Must be an LP. */
	if (*iobj == 0) {
	    prob = 0;
	} else {
	    prob = 1;
	}
    } else {
/*  Explicit objective. Check for quadratic term. */
	if (*nnh > 0) {
	    prob = 2;
	} else {
	    prob = 1;
	}
    }
    iw[223] = 0;
/* Print the header for the Print   file */
    iw[225] = 0;
/* Print the header for the summary file */
    iw[220] = 0;
/* Line count for the print   file */
    iw[221] = 0;
/* Initialize counters based on gotHes and gotFac (set in s3prtQ) */
/* Line count for the summary file */
    if (iw[230] <= 0) {
	iw[210] = 0;
    }
    if (iw[231] <= 0) {
	iw[188] = 0;
    }
    iw[386] = 0;
    iw[387] = 0;
    itn = 0;
    itqp = 0;
    itqpmax = itnlim;
    ndegen = 0;
    nncon = 0;
    nncon0 = 1;
    nnjac = 0;
    numlc = *m;
    ngqp0 = max(*ngqp,1);
    iw[200] = 0;
/* QP Hessian may or may not be definite */
    iw[208] = qpslvr;
/* Local value of QPslvr */
    *objqp = 0.;
    sclobj = 1.;
/* Initialize quantities to avoid them being used before being set. */
    dload_(m, &c_b2, &pi[1], &c__1);
    dload_(&ngqp0, &c_b2, &rw[lgqp], &c__1);
    iload_(nb, &c__0, &iw[lhesta], &c__1);
/* ----------------------------------------------------------------- */
/* Print the matrix statistics. */
/* Find the rowtypes for use in s5getB (they are held in iy2). */
/* ----------------------------------------------------------------- */
    s2amat_(&c__1, &mjrprt, m, n, nb, &nncon, &nnjac, &nnobj, iobj, &numlc, &
	    numliq, ne, nloca, &loca[1], &inda[1], &acol[1], &bl[1], &bu[1], &
	    iw[liy2], &iw[1], leniw, &rw[1], lenrw);
/* ================================================================= */
/* Find a basis kBS(1:m) for the linear constraints and bounds. */
/* ================================================================= */
/* s5getB does the following. */
/*  1. The linear constraints are (optionally) scaled. */
/*  2. Elements x(n+1:n+m) of the initial x are assigned. */
/*  3. An LP is used to find a feasible x for the bounds and */
/*     linear equality constraints. */
/*  The base point x0 is not touched. */
    s5getb_(&inform__, start, (U_fp)qplog, &needb, m, &maxs, &mbs, n, nb, &
	    nncon, &nnjac, &nnobj, nname, ns, &itqp, &itqpmax, &itn, &ndegen, 
	    &numlc, &numliq, &tolfp, &tolqp, &tolx, ninf, sinf, &wtinf0, iobj,
	     &sclobj, &pinorm, &rgnorm, ne, nloca, &loca[1], &inda[1], &acol[
	    1], &hetype[1], &iw[lhesta], &iw[liy2], &iw[lhfeas], &hs[1], &iw[
	    lkbs], names + 8, &rw[lascal], &bl[1], &bu[1], &rw[lblbs], &rw[
	    lbubs], &rw[lblsav], &rw[lbusav], &rw[lgbs], &pi[1], &rc[1], 
	    nrhs0, nrhs, &rhs[1], &c__1, &c__0, &x0[1], &x[1], &rw[lxbs], &iw[
	    liy], &iw[liy1], &rw[ly], &rw[ly1], &rw[ly2], cw + 8, lencw, &iw[
	    1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
    if (*ngobj > 0 && iw[75] > 0) {
	ddscl_(ngobj, &rw[lascal], &c__1, &gobj[1], &c__1);
    }
/*     Possible inform values are = -3,-2,-1, 0, >0 */
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
    if (*iexit != 0) {
	goto L900;
    }
/*     ================================================================== */
/*     Solve the problem. */
/*     ================================================================== */
    s1time_(&c__2, &c__0, &iw[1], leniw, &rw[1], lenrw);
    useqp = *ngobj > 0 && prob == 1 || prob == 2;
    if (useqp) {
	if (iw[208] == 0) {
	    s5sqp_(&inform__, (U_fp)hprod, (U_fp)hprod1, (U_fp)qplog, gotr, &
		    prob, &lenr, m, &maxs, &mbs, n, nb, &ndegen, &iw[188], &
		    ngqp0, ngqp, ngobj0, ngobj, nnh0, nnh, ns, &itqp, &
		    itqpmax, &itn, &minimz, iobj, &sclobj, objadd, objqp, &
		    tolfp, &tolqp, &tolx, ninf, sinf, &wtinf0, &pinorm, ne, 
		    nloca, &loca[1], &inda[1], &acol[1], &hetype[1], &iw[
		    lhesta], &iw[lhfeas], &hs[1], &iw[lkbs], &rw[lascal], &bl[
		    1], &bu[1], &rw[lblbs], &rw[lbubs], &rw[lgbs], &gobj[1], &
		    rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], &
		    rw[lrg], nrhs0, nrhs, &rhs[1], lenx0, nx0, &x0[1], &x[1], 
		    &rw[lxbs], &x[1], &iw[liy], &iw[liy1], &rw[ly], &rw[ly1], 
		    &rw[ly2], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw 
		    + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		    ftnlen)8);
/* xFreez = x, not used */
	} else if (iw[208] == 2) {
	    *gotr = iw[231] > 0;
	    s5sqn_(&inform__, (U_fp)hprod, (U_fp)hprod1, (U_fp)qplog, gotr, &
		    prob, &lenr, m, &maxs, &mbs, n, nb, &ndegen, &iw[188], &
		    ngqp0, ngqp, ngobj0, ngobj, nnh0, nnh, ns, &itqp, &
		    itqpmax, &itn, &minimz, iobj, &sclobj, objadd, objqp, &
		    tolfp, &tolqp, &tolx, ninf, sinf, &wtinf0, &pinorm, ne, 
		    nloca, &loca[1], &inda[1], &acol[1], &hetype[1], &iw[
		    lhesta], &iw[lhfeas], &hs[1], &iw[lkbs], &rw[lascal], &bl[
		    1], &bu[1], &rw[lblbs], &rw[lbubs], &rw[lgbs], &gobj[1], &
		    rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], &
		    rw[lrg], &rw[lrg2], nrhs0, nrhs, &rhs[1], lenx0, nx0, &x0[
		    1], &x[1], &rw[lxbs], &x[1], &iw[liy], &iw[liy1], &rw[ly],
		     &rw[ly1], &rw[ly2], cu + 8, lencu, &iu[1], leniu, &ru[1],
		     lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		    ftnlen)8, (ftnlen)8);
/* xFreez = x, not used */
	}
	if (inform__ == 33 && maxr < maxs) {
	    if (iw[208] == 0) {
		iw[209] = 0;
/* with no preconditioning */
		*gotr = FALSE_;
	    } else if (iw[208] == 2) {
		iw[209] = 1;
/* with QN preconditioning */
	    }
	    iw[208] = 1;
/* Switch to CG */
	}
	if (iw[208] == 1) {
	    if (qpslvr == 1) {
		*gotr = FALSE_;
	    }
	    s5sqn_(&inform__, (U_fp)hprod, (U_fp)hprod1, (U_fp)qplog, gotr, &
		    prob, &lenr, m, &maxs, &mbs, n, nb, &ndegen, &iw[188], &
		    ngqp0, ngqp, ngobj0, ngobj, nnh0, nnh, ns, &itqp, &
		    itqpmax, &itn, &minimz, iobj, &sclobj, objadd, objqp, &
		    tolfp, &tolqp, &tolx, ninf, sinf, &wtinf0, &pinorm, ne, 
		    nloca, &loca[1], &inda[1], &acol[1], &hetype[1], &iw[
		    lhesta], &iw[lhfeas], &hs[1], &iw[lkbs], &rw[lascal], &bl[
		    1], &bu[1], &rw[lblbs], &rw[lbubs], &rw[lgbs], &gobj[1], &
		    rw[lgqp], &rw[lhdx], &rw[lpbs], &pi[1], &rw[lr], &rc[1], &
		    rw[lrg], &rw[lrg2], nrhs0, nrhs, &rhs[1], lenx0, nx0, &x0[
		    1], &x[1], &rw[lxbs], &x[1], &iw[liy], &iw[liy1], &rw[ly],
		     &rw[ly1], &rw[ly2], cu + 8, lencu, &iu[1], leniu, &ru[1],
		     lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		    ftnlen)8, (ftnlen)8);
/* xFreez = x, not used */
	}
	*iexit = inform__;
    } else {
	s5fixs_(&c__0, m, &maxs, &mbs, n, nb, ns, &hs[1], &iw[lkbs], &bl[1], &
		bu[1], &rw[lblbs], &rw[lbubs], &x[1], &rw[lxbs]);
	s5slp_(iexit, (U_fp)qplog, &prob, m, &mbs, n, nb, &ndegen, &itqp, &
		itqpmax, &itn, &minimz, iobj, &sclobj, objadd, &tolfp, &tolqp,
		 &tolx, ninf, sinf, &wtinf0, &pinorm, ne, nloca, &loca[1], &
		inda[1], &acol[1], &hetype[1], &iw[lhesta], &iw[lhfeas], &hs[
		1], &iw[lkbs], &rw[lascal], &bl[1], &bu[1], &rw[lblbs], &rw[
		lbubs], &rw[lgbs], &pi[1], &rc[1], nrhs0, nrhs, &rhs[1], &x[1]
		, &rw[lxbs], &x[1], &iw[liy], &iw[liy1], &rw[ly], &rw[ly1], 
		cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
/* xFreez = x, not used */
    }
    s1time_(&c_n2, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     ================================================================== */
/*     Exit. */
/*     Set output variables and print a summary of the final solution. */
/*     ObjTru is printed in s4newB */
/*     ================================================================== */
L900:
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)132, (
	    ftnlen)132);
    degen = ndegen * 100. / max(itn,1);
    objlp = 0.;
    if (prob == 0) {
	*objtru = 0.;
    } else if (prob == 1 || prob == 2) {
	*objtru = *objadd;
	if (*iobj > 0) {
	    objlp = x[*n + *iobj] * sclobj;
	    *objtru += objlp;
	}
	if (*ngqp > 0) {
	    *objtru += *objqp;
	}
    }
    infsbl = *ninf > 0;
    xnorm = dnormi_(n, &x[1], &c__1);
/* Count basic nonlinear variables (used only for printing). */
    nnb = 0;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	if (hs[j] == 3) {
	    ++nnb;
	}
    }
    if (inewb > 0 && *iexit / 10 < 8) {
	k = *iexit / 10 + 1;
	s4stat_(&k, istate, (ftnlen)4);
	s4newb_(&c__1, &inewb, &minimz, m, n, nb, ns, &mbs, &itn, ninf, sinf, 
		objtru, &iw[lkbs], &hs[1], &rw[lascal], &bl[1], &bu[1], &x[1],
		 &rw[lxbs], istate, cw + 8, lencw, &iw[1], leniw, (ftnlen)4, (
		ftnlen)8);
    }
/* Print statistics. */
/* Writing concatenation */
    i__2[0] = 30, a__1[0] = " Problem name                 ";
    i__2[1] = 8, a__1[1] = mprob;
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)38);
    snprnt_(&c__13, ch__1, &iw[1], leniw, (ftnlen)38);
    s_wsfi(&io___71);
    do_fio(&c__1, (char *)&itn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*objtru), (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    if (infsbl) {
	s_wsfi(&io___72);
	do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
    if (prob == 2) {
	s_wsfi(&io___73);
	do_fio(&c__1, (char *)&iw[188], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&objlp, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___74);
	do_fio(&c__1, (char *)&(*objqp), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
    if (*ns > 0) {
	s_wsfi(&io___75);
	do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nnb, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
    s_wsfi(&io___76);
    do_fio(&c__1, (char *)&ndegen, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&degen, (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
/*     ------------------------------------------------------------------ */
/*     Unscale, save basis files and prepare to print the solution. */
/*     Clock 3 is "Output time". */
/*     ------------------------------------------------------------------ */
    s1time_(&c__3, &c__0, &iw[1], leniw, &rw[1], lenrw);
    s4savb_(&inform__, &c__0, &minimz, m, n, nb, &nkx, &nncon0, &nncon, &
	    ngqp0, ngqp, nname, ns, &itn, ninf, sinf, &wtinf0, &vimax, iobj, &
	    sclobj, objtru, &pnorm1, &pnorm2, &pinorm, &xnorm, ne, nloca, &
	    loca[1], &inda[1], &acol[1], &iw[lkx], &iw[lhesta], &hs[1], &rw[
	    lascal], &bl[1], &bu[1], fx, &rw[lgqp], names + 8, &pi[1], &rc[1],
	     &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
	    ftnlen)8);
    if (*ngobj > 0 && iw[75] > 0) {
	dddiv_(ngobj, &rw[lascal], &c__1, &gobj[1], &c__1);
    }
/*     If task = PrintS, s4savB prints the solution under the control */
/*     of lprSol (set by the  Solution  keyword in the SPECS file). */
/*     The printed solution may or may not be wanted, as follows: */

/*     lprSol = 0   means      No */
/*            = 1   means      If optimal, infeasible or unbounded */
/*            = 2   means      Yes */
/*            = 3   means      If error condition */
    s4savb_(&inform__, &c__1, &minimz, m, n, nb, &nkx, &nncon0, &nncon, &
	    ngqp0, ngqp, nname, ns, &itn, ninf, sinf, &wtinf0, &vimax, iobj, &
	    sclobj, objtru, &pnorm1, &pnorm2, &pinorm, &xnorm, ne, nloca, &
	    loca[1], &inda[1], &acol[1], &iw[lkx], &iw[lhesta], &hs[1], &rw[
	    lascal], &bl[1], &bu[1], fx, &rw[lgqp], names + 8, &pi[1], &rc[1],
	     &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
	    ftnlen)8);
    s1time_(&c_n3, &c__0, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Set Obj for output. */
/*     Call  Hx  one last time with  nState .ge. 2. */
/*     Everything has been  unscaled, so we have to disable scaling. */
/*     ------------------------------------------------------------------ */
    lssave = iw[75];
    iw[75] = 0;
/* Computing MIN */
    i__1 = *iexit / 10;
    sqstat = min(i__1,4) + 2;
    iw[235] = sqstat;
    objlp = 0.;
    if (prob == 0) {
	*objtru = 0.;
    } else if (prob == 1 || prob == 2) {
	*objtru = *objadd;
	if (*iobj > 0) {
	    objlp = x[*n + *iobj] * sclobj;
	    *objtru += objlp;
	}
	if (*ngqp > 0) {
	    s5qpfg_((U_fp)hprod, (U_fp)hprod1, ngqp, ngobj0, ngobj, nnh, &
		    sqstat, &iw[188], objqp, &gobj[1], &rw[lgqp], lenx0, nx0, 
		    &x0[1], &x[1], &rw[ly], cu + 8, lencu, &iu[1], leniu, &ru[
		    1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		    ftnlen)8, (ftnlen)8);
	    *objtru += *objqp;
	}
    }
    iw[75] = lssave;
/*     Save some things needed by solvers calling SQOPT */
/* L999: */
    rw[421] = *objtru;
/* The true objective */
    rw[422] = pinorm;
/* Lagrange multiplier norm */
    rw[423] = xnorm;
    iw[421] = itn;
/* Total iteration count */
    iw[423] = maxs;
/* max # of superbasics */
    return 0;
} /* s5solv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5solv */
/* Subroutine */ int s5dflt_(integer *m, integer *n, integer *lencobj, 
	integer *ncolh, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal c4, c6;
    static logical qp;
    static doublereal eps, eps0, eps1, eps2, eps3, eps4;
    static integer npr1, npr2, kfac, kchk, klog, ksav, maxr, imps, maxs, nout;
    static doublereal tolx, dens1, dens2, lmax1, lmax2, utol1, utol2;
    static integer iback, ioldb;
    static doublereal bigdx, bigfx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal small, tolcg, xdlim;
    static integer idump, maxmn, never, mskip, isoln;
    static doublereal tolfp, rmaxs;
    static integer ksumm;
    static doublereal tolqp, toldj3, wtinf0, utol1m, hcndbd, utol2m;
    static integer iloadb, kdegen;
    static doublereal infbnd, zcndbd;
    static integer lemode;
    static doublereal chzbnd;
    static integer icrash;
    static logical linear;
    static integer lprdbg;
    static doublereal tolfac, uspace;
    static integer maxcol;
    static doublereal tcrash, toldcp;
    static integer precon;
    static doublereal tolddp;
    static integer minprc, minmax, lvlinf, cgitmx, itnlim, kreset, mflush, 
	    lprscl, lvlscl, mminor, mnewsb, minimz, lvlpre, iprint, ireprt, 
	    nparpr, iinsrt;
    static doublereal scltol, tolcon;
    static integer lprsol, lprprm, lvlpiv, mjrprt;
    static doublereal toldpp, toldrp, toldup;
    static integer mnrprt;
    static doublereal tolnlp;
    static integer luprnt, tpivot;
    static doublereal tolpiv;
    static integer qpslvr;
    static doublereal tolrow;
    static integer stkyop;
    static doublereal tolswp, tolupd;
    static integer lvlsys;

/*     ================================================================== */
/*     s5dflt checks and possibly prints the optional parameter values */
/*     for sqopt. */

/*     Optional parameters are checked and, if necessary,  changed to */
/*     reasonable values. */

/*     Note that parameters are checked before the amount of working */
/*     storage has been defined. */

/*     See  snworkspace.info  for full documentation of cw, iw and rw. */

/*     15 Nov 1991: first version. */
/*     02 Aug 2003: snPRNT adopted. */
/*     22 Jun 2004: Added default LU mod singularity tol */
/*     21 Dec 2004: Default LU tols fixed up. */
/*     02 May 2006: lvlTim removed. */
/*     01 Sep 2007: stkyOp added. */
/*     18 Jan 2010: PreCon initialized. */
/*     13 Jul 2013: Default lvldif set correctly. */
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
    tolcg = rw[54];
/* cg tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance */
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
/* LU factor tolerance */
    tolupd = rw[67];
/* LU update tolerance */
    infbnd = rw[70];
/* definition of plus infinity */
    bigfx = rw[71];
/* unbounded objective */
    bigdx = rw[72];
/* unbounded step */
    lvlpre = iw[77];
/* >0    => QN preconditioned CG */
    xdlim = rw[80];
/* Step limit */
    etarg = rw[83];
/* Quasi-Newton QP rg tolerance */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    wtinf0 = rw[88];
/* infeasibility weight */
    scltol = rw[92];
/*     ------------------------------------------------------------------ */
/*     rw(151)--rw(180) contain  parmLU  parameters for LUSOL. */
/*     ------------------------------------------------------------------ */
/* scale tolerance. */
    lmax1 = rw[151];
/* max L-multiplier in factor */
    lmax2 = rw[152];
/* max L-multiplier in update */
    small = rw[153];
/* defn of small real */
    utol1 = rw[154];
/* abs tol for small diag of U */
    utol2 = rw[155];
/* rel tol for small diag of U */
    uspace = rw[156];
/* limit on waste space in U */
    dens1 = rw[157];
/* switch to search maxcol columns and no rows */
    dens2 = rw[158];
/*     ------------------------------------------------------------------ */
/*     rw(181)--rw(199) pass parameters into various routines. */
/*     ------------------------------------------------------------------ */
/*     toldj3    = rw(186) ! current optimality tol */
/*     ------------------------------------------------------------------ */
/*     iw(1)--iw(50): I/O file numbers and dimensions. */
/*     ------------------------------------------------------------------ */
/* switch to dense LU */
    iprint = iw[12];
/*     ------------------------------------------------------------------ */
/*     iw(51)--iw(150): optional parameters set via the specs file. */
/*     ------------------------------------------------------------------ */
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
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
/*     lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start */
/* # largest value of nSkip */
    lvlsys = iw[71];
/* > 0   => print system info */
    lvlinf = iw[73];
/* Elastic option */
    lvlscl = iw[75];
/* scale option */
    lvlpiv = iw[80];
/* 0(1) LU threshold partial(complete) pivoting */
    lprprm = iw[81];
/* > 0    => parms are printed */
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
    mminor = iw[91];
/* limit on minor iterations */
    mnrprt = iw[93];
/* Minor print level */
    nparpr = iw[94];
/* # of partial pricing sections */
    mnewsb = iw[95];
/* # of working set changes */
    cgitmx = iw[97];
/* CG iteration limit */
    stkyop = iw[116];
/* > 0 => optional parameters are sticky */
    iback = iw[120];
/* backup file */
    idump = iw[121];
/* dump file */
    iloadb = iw[122];
/* load file */
    imps = iw[123];
/* MPS file */
    inewb = iw[124];
/* new basis file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    ipnch = iw[127];
/* punch file */
    ireprt = iw[130];
/* Report file */
    isoln = iw[131];
/*     ------------------------------------------------------------------ */
/*     iw(151)--iw(180) contain luparm parameters for LUSOL. */
/*     ------------------------------------------------------------------ */
/* Solution file */
    nout = iw[151];
/* unit # for printed messages */
    luprnt = iw[152];
/* print level in LU routines */
    maxcol = iw[153];
/*     ------------------------------------------------------------------ */
/* lu1fac: max. # columns */
    precon = iw[209];
/* Current precon mode (based on QPslvr) */
    c4 = max(1e-4,eps3);
    c6 = max(1e-6,eps2);
    never = 99999999;
    qp = *ncolh > 0;
    linear = ! qp;
/*     ================================================================== */
/*     Check the optional parameters. */
/*     ================================================================== */
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
	if (qp) {
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
/*     Sometimes, frequency 0 means "almost never". */
    if (kchk <= 0) {
	kchk = never;
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
    if (mminor < 0) {
/* Computing MAX */
	i__1 = 1000, i__2 = max(*n,*m) * 5;
	mminor = max(i__1,i__2);
    }
    if (mnewsb <= 0) {
	mnewsb = never;
    }
    if (lprdbg < 0) {
	lprdbg = 0;
    }
    if (lprprm < 0) {
	lprprm = 1;
    }
    if (lprscl < 0) {
	lprscl = 0;
    }
    if (lprsol < 0) {
	lprsol = 2;
    }
/*     lvlSrt is checked in s3argA or s3argB */
/*     if (lvlSrt .lt. 0      ) lvlSrt = 0 */
    if (mnrprt < 0) {
	mnrprt = 1;
    }
    mjrprt = mnrprt;
    if (lvlinf < 0 || lvlinf > 2) {
	lvlinf = -11111;
    }
    if (lvlinf == -11111) {
	lvlinf = 2;
    }
    if (lvlsys < 0) {
	lvlsys = 0;
    }
    if (lemode < 0 || lemode > 2) {
	lemode = -11111;
    }
    if (lemode == -11111) {
	lemode = 1;
    }
    if (stkyop < 0) {
	stkyop = 0;
    }
/*     Check superbasics limit and reduced Hessian size. */
    if (qp) {
	if (maxr < 0) {
/* Computing MIN */
	    i__1 = 2000, i__2 = *ncolh + 1;
	    maxr = min(i__1,i__2);
	}
	if (maxs < 0) {
	    maxs = *ncolh + 1;
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
	maxs = maxr;
    }
    if (qpslvr < 0) {
	qpslvr = 0;
    }
    if (maxr == 0) {
	qpslvr = 1;
    }
    if (lvlpre < 0 || lvlpre > 1) {
	lvlpre = 0;
	precon = 0;
    } else {
	precon = 1;
    }
    if (cgitmx < 0) {
	cgitmx = 100;
    }
    if (etarg < 0. || etarg > 1.) {
	etarg = .5;
    }
/*     Check other options. */
    if (lvlscl < 0) {
	lvlscl = 2;
    }
    lvlscl = min(lvlscl,2);
    if (nparpr <= 0) {
	nparpr = 10;
    }
    minprc = 10;
    npr1 = *n / nparpr;
    npr2 = *m / nparpr;
    if (max(npr1,npr2) < minprc) {
	maxmn = max(*m,*n);
	nparpr = maxmn / min(maxmn,minprc);
    }
    rmaxs = (doublereal) maxs;
/* Computing MAX */
    d__1 = 1. / (eps * 100. * rmaxs);
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
    if (xdlim <= 0.) {
	xdlim = 2.;
    }
    if (zcndbd <= 0.) {
	if (qpslvr == 0) {
	    zcndbd = 1e4;
	} else {
	    zcndbd = 1e6;
	}
    }
    if (tcrash < 0. || tcrash >= 1.) {
	tcrash = .1;
    }
/*     ------------------------------------ */
/*     Set up the parameters for lu1fac. */
/*     ------------------------------------ */
    if (maxcol < 0) {
	maxcol = 5;
    }
    if (luprnt == -11111) {
	luprnt = -1;
    }
    nout = iprint;
    if (lvlsys == 0) {
	nout = 0;
    }
    if (mnrprt > 10) {
	luprnt = 0;
    }
    if (lprdbg == 51) {
	luprnt = 1;
    }
    if (lprdbg == 52) {
	luprnt = 2;
    }
    if (iprint < 0) {
	luprnt = -1;
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
/* QP */
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
/*     Set some tolerances. */
/*     Set the optimality tolerance. */
/*     Solve the QP subproblems fairly accurately. */
    if (tolcg <= 0.) {
	tolcg = .01;
    }
    if (tolqp <= 0.) {
	tolqp = c6;
    }
    if (tolfp <= 0.) {
	tolfp = tolqp;
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
    if (wtinf0 < 0.) {
	wtinf0 = 1.;
    }
    if (iback == inewb) {
	iback = 0;
    }
    if (itnlim < 0) {
/* Computing MAX */
	i__1 = 10000, i__2 = max(*n,*m) * 10;
	itnlim = max(i__1,i__2);
    }
/*     Load tolerances used to mark variables during printing in s4SavB. */
    tolnlp = tolqp;
    tolcon = tolx;
/*     ------------------------------------------------------------------ */
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
    rw[65] = tolswp;
    rw[66] = tolfac;
    rw[67] = tolupd;
    rw[70] = infbnd;
    rw[71] = bigfx;
    rw[72] = bigdx;
    rw[80] = xdlim;
    rw[83] = etarg;
    rw[85] = hcndbd;
    rw[86] = zcndbd;
    rw[88] = wtinf0;
    rw[92] = scltol;
    rw[151] = lmax1;
/* max L-multiplier in factor */
    rw[152] = lmax2;
/* max L-multiplier in update */
    rw[153] = small;
/* defn of small real */
    rw[154] = utol1;
/* abs tol for small diag of U */
    rw[155] = utol2;
/* rel tol for small diag of U */
    rw[156] = uspace;
/* limit on waste space in U */
    rw[157] = dens1;
/* switch to search maxcol columns and no rows */
    rw[158] = dens2;
/* switch to dense LU */
    rw[181] = toldpp;
    rw[182] = toldcp;
    rw[183] = toldup;
    rw[186] = toldj3;
    rw[187] = toldrp;
    iw[52] = maxr;
    iw[53] = maxs;
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
/*     iw( 69) = lvlSrt */
    iw[71] = lvlsys;
    iw[73] = lvlinf;
    iw[75] = lvlscl;
    iw[77] = lvlpre;
    iw[80] = lvlpiv;
    iw[81] = lprprm;
    iw[83] = lprscl;
    iw[84] = lprsol;
    iw[85] = lprdbg;
    iw[87] = minmax;
    iw[88] = icrash;
    iw[89] = itnlim;
    iw[91] = mminor;
    iw[92] = mjrprt;
    iw[93] = mnrprt;
    iw[94] = nparpr;
    iw[95] = mnewsb;
    iw[97] = cgitmx;
    iw[116] = stkyop;
    iw[120] = iback;
    iw[121] = idump;
    iw[122] = iloadb;
    iw[123] = imps;
    iw[124] = inewb;
    iw[125] = iinsrt;
    iw[126] = ioldb;
    iw[127] = ipnch;
    iw[130] = ireprt;
    iw[131] = isoln;
    iw[151] = nout;
    iw[152] = luprnt;
    iw[153] = maxcol;
    iw[156] = tpivot;
    iw[199] = minimz;
    iw[209] = precon;
/* not optional parameters, but set here. */
    rw[186] = toldj3;
    return 0;
} /* s5dflt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5dflt */
/* Subroutine */ int s5map_(integer *m, integer *n, integer *nkx, integer *
	ngobj, integer *nnh, integer *lenr, integer *maxr, integer *maxs, 
	integer *nextcw, integer *nextiw, integer *nextrw, integer *iw, 
	integer *leniw)
{
    static integer nb, lr, ly, lr1, lr2, ls1, ls2, ls3, ly1, ly2, ly3, mbs, 
	    lrg, ldx, liy, lkx, lrg2, liy1, liy2, lgbs, lkbs, lhdx, lpbs, 
	    lgqp, ngqp, lxbs, lblbs, lbubs, lascal, lhfeas, lhesta, lblsav, 
	    lxscal, lbusav, lqprhs;

/*     ================================================================== */
/*     s5Map   allocates all array storage for sqopt, */
/*     using the values: */
/*        m    , n    , ne */
/*        maxS                    Set in s5dflt. */
/*        ngObj, nnH              Set from the argument list. */
/*        lenR                    Set in the calling program. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2core. */
/*     12 Nov 1994: Converted to integer and real storage. */
/*     06 Aug 1996: First min sum version. */
/*     14 Jul 1997: Thread-safe version. */
/*     01 May 1998: First version called by sqMem. This simplified */
/*                  version may slightly overestimate needed memory. */
/*     02 Aug 2003: snPRNT adopted. */
/*     13 May 2005: Bug fix: ly3 assigned to iw correctly */
/*     18 Jun 2008: Added space for iy2, pBS and rg2. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ngqp = max(*ngobj,*nnh);
    mbs = *m + *maxs;
    nb = *n + *m;
/*     sqopt can use all of cw, iw and rw */
/*     except the first user workspace partitions. */
    lkx = *nextiw;
    lhfeas = lkx + *nkx;
    lkbs = lhfeas + mbs;
    lhesta = lkbs + mbs;
    liy = lhesta + nb;
    liy1 = liy + nb;
    liy2 = liy1 + nb;
    *nextiw = liy2 + nb;
/*     Addresses for the double precision arrays. */
    lascal = *nextrw;
    ly = lascal + nb;
    ly1 = ly + nb;
    ly2 = ly1 + nb;
    if (*maxr < *maxs) {
/* Define SYMMLQ workspace */
	ly3 = ly2 + nb;
	ls1 = ly3 + nb;
	ls2 = ls1 + *maxs;
	ls3 = ls2 + *maxs;
	lr1 = ls3 + *maxs;
	lr2 = lr1 + *maxs;
	lblbs = lr2 + *maxs;
    } else {
	ly3 = ly2 + nb;
	ls1 = ly3;
	ls2 = ls1;
	ls3 = ls2;
	lr1 = ls3;
	lr2 = lr1;
	lblbs = lr2;
    }
    lbubs = lblbs + mbs;
    lxbs = lbubs + mbs;
    lxscal = lxbs + mbs;
    lhdx = lxscal + *nnh;
    lpbs = lhdx + *nnh;
    lgqp = lpbs + nb;
    lgbs = lgqp + ngqp;
    lr = lgbs + mbs;
    lrg = lr + *lenr;
    lrg2 = lrg + *maxs;
    lblsav = lrg2 + *maxs;
    lbusav = lblsav + nb;
    lqprhs = lbusav + nb;
    ldx = lqprhs + *m;
    *nextrw = ldx + ngqp;
/*     --------------------------- */
/*     Store the addresses in iw. */
/*     --------------------------- */
    iw[251] = lkx;
    iw[273] = lblbs;
    iw[274] = lbubs;
    iw[275] = lblsav;
    iw[276] = lbusav;
    iw[277] = lpbs;
    iw[278] = lqprhs;
    iw[284] = lhfeas;
    iw[285] = lhesta;
    iw[287] = ldx;
    iw[288] = lhdx;
    iw[290] = lgqp;
    iw[291] = lgbs;
    iw[292] = lkbs;
    iw[293] = lrg;
    iw[294] = lrg2;
    iw[295] = lr;
    iw[296] = lascal;
    iw[301] = lxbs;
    iw[302] = lxscal;
    iw[308] = liy;
    iw[309] = liy1;
    iw[310] = liy2;
    iw[311] = ly;
    iw[312] = ly1;
    iw[313] = ly2;
    iw[314] = ly3;
    iw[353] = lr1;
    iw[354] = lr2;
    iw[355] = ls1;
    iw[356] = ls2;
    iw[357] = ls3;
    return 0;
} /* s5map_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Map */
/* Subroutine */ int s5slp_(integer *iexit, U_fp qplog, integer *prob, 
	integer *m, integer *mbs, integer *n, integer *nb, integer *ndegen, 
	integer *itqp, integer *itqpmax, integer *itn, integer *minimz, 
	integer *iobj, doublereal *sclobj, doublereal *objadd, doublereal *
	tolfp, doublereal *tolqp, doublereal *tolx, integer *ninf, doublereal 
	*sinf, doublereal *wtinf, doublereal *pinorm, integer *ne, integer *
	nloca, integer *loca, integer *inda, doublereal *acol, integer *
	hetype, integer *hestat, integer *hfeas, integer *hs, integer *kbs, 
	doublereal *ascale, doublereal *bl, doublereal *bu, doublereal *blbs, 
	doublereal *bubs, doublereal *gbs, doublereal *pi, doublereal *rc, 
	integer *nrhs0, integer *nrhs, doublereal *rhs, doublereal *x, 
	doublereal *xbs, doublereal *xfreez, integer *iy, integer *iy1, 
	doublereal *y, doublereal *y1, char *cw, integer *lencw, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in feasibility\002,\002 phase\002)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Feasible non-elastic"
	    "s\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Feasible constraints"
	    "\002)";
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in LP optimality\002,\002 phase\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer ns, nnh;
    static char str[120];
    extern /* Subroutine */ int s5lp_(integer *, integer *, char *, logical *,
	     integer *, U_fp, logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static logical done, newb, luok, badlu;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static logical needx, found;
    static integer lureq;
    static logical newlu;
    static integer nswap;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static integer lemod0, lvlin0, gotfac;
    static logical infeas;
    static integer lemode;
    static logical elastc, itnbig, needlu;
    static integer inform__, lvlinf;
    static logical getopt;
    static doublereal rgnorm;
    static logical fponly;
    static integer subopt, mnrprt, prtlvl, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static logical unbnded;
    extern /* Subroutine */ int s2trylu_(integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *);
    static char probtag[20];

    /* Fortran I/O blocks */
    static icilist io___252 = { 0, str, 0, fmt_1000, 120, 1 };
    static icilist io___253 = { 0, str, 0, fmt_1200, 120, 1 };
    static icilist io___254 = { 0, str, 0, fmt_1300, 120, 1 };
    static icilist io___255 = { 0, str, 0, fmt_1100, 120, 1 };


/* ================================================================= */
/* s5sLP   solves the current problem.  An initial basis is assumed */
/* to be specified by nS, hs, x and the superbasic parts of kBS. */

/* 25 Oct 2003: First version of s5sLP based on s5SQP. */
/* 08 Mar 2004: gotFac implemented. */
/* 18 Jun 2008: xFreez added as argument. */
/* 29 Apr 2011: Some minor cosmetic changes. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* >0 => Mnr heading for iPrint */
/* ----------------------------------------------------------------- */
/* >0 => Mnr heading for iSumm */
    /* Parameter adjustments */
    --pi;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --y1;
    --y;
    --iy1;
    --iy;
    --xfreez;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --acol;
    --inda;
    --loca;
    --rhs;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lemode = iw[56];
/* >0    => use elastic mode */
    lvlinf = iw[73];
/* Elastic option */
    mnrprt = iw[93];
/* Minor print level */
    gotfac = iw[230];
/* >0 => Save the LU factors */
    prtlvl = mnrprt;
    s_copy(probtag, "linear constraints", (ftnlen)20, (ftnlen)18);
    iload_(nb, &c__0, &hestat[1], &c__1);
    elastc = lemode == 2;
    getopt = *prob == 1 || *prob == 2;
    fponly = *prob == 0;
    ns = 0;
/* Local value */
    nnh = 0;
/* ----------------------------------------------------------------- */
/* Call s5LP with argument "FP" to find a feasible point. */
/* ----------------------------------------------------------------- */
/* Local value */
    if (elastc) {
/* Phase 2 will start in elastic mode with elastic objective */
/* determined by lvlInf. For FP, make the nonelastics feasible. */
	lemod0 = lemode;
	lvlin0 = lvlinf;
    } else {
/* FP will make all constraints feasible. */
/* Enter elastic mode if infeasible and minimize the elastic */
/* infeasibilities. */
	lemod0 = 1;
	lvlin0 = 2;
    }
    subopt = -1;
    if (gotfac == 0) {
	lureq = 1;
    } else {
	lureq = 0;
    }
    typelu = 3;
    luok = TRUE_;
    done = FALSE_;
/*     ================================================================== */
/* +    while (.not. done  .and.  LUok) do */
L500:
    if (! done && luok) {
	needlu = lureq > 0;
	s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, &prtlvl, &
		lureq, m, mbs, n, nb, &nnh, &ns, &nswap, ne, nloca, &loca[1], 
		&inda[1], &acol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &blbs[1],
		 &bubs[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], &iy[1], &iy1[
		1], &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    goto L900;
	}
	needx = needlu;
	s5lp_(&inform__, &c__0, probtag, &elastc, &subopt, (U_fp)qplog, &
		needlu, &needx, m, n, nb, ndegen, itqp, itqpmax, itn, &lemod0,
		 &lvlin0, &prtlvl, minimz, iobj, sclobj, objadd, tolfp, tolqp,
		 tolx, ninf, sinf, wtinf, pinorm, &rgnorm, ne, nloca, &loca[1]
		, &inda[1], &acol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[
		1], &kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &
		gbs[1], &pi[1], &rc[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], 
		&xfreez[1], &iy[1], &iy1[1], &y[1], &y1[1], cw + 8, lencw, &
		iw[1], leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8);
/* Check for trouble in s5LP.  Here are the possibilities: */
	badlu = inform__ > 0;
/* Fatal LU error */
	found = inform__ == 0;
/* LP solution found */
	infeas = inform__ == -1;
/* LP is infeasible */
	unbnded = inform__ == -2;
/* LP is unbounded */
	itnbig = inform__ == -3;
/* Too many iterations */
	if (badlu) {
	    *iexit = inform__;
/* Fatal LU error */
	    goto L900;
	}
	done = found || infeas || itnbig;
	if (! done) {
/* =========================================================== */
/* Trouble. */
/* The phase 1 was unbounded, which can only occur if a bad */
/* basis has given a large search direction. */
/* =========================================================== */
/* Assume infeasible constraints.  Factor with tighter tols. */
	    infeas = TRUE_;
	    s_wsfi(&io___252);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    s2trylu_(itn, &c__22, &ns, &lureq, &luok, &typelu, &iw[1], leniw, 
		    &rw[1], lenrw);
	}
	goto L500;
    }
/* +    end while */
/*     ================================================================== */
    if (! found) {
	goto L800;
    }
/* Itns or inf */
    if (*ninf > 0 && lemode != 1) {
	goto L800;
    }
/* Infeas */
    if (getopt) {
/*        call dcopy ( nb, x, 1, xFreez, 1 ) ! Save feasible nonelastics */
/*        if (PrtLvl .ge. 1  .and.  PrtLvl .lt. 10) then */
	if (prtlvl >= 10) {
	    if (elastc) {
		s_wsfi(&io___253);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    } else {
		s_wsfi(&io___254);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    }
	}
    }
/* ================================================================= */
/* Get optimal. */
/* ================================================================= */
    if (elastc) {
/* Relax */
    } else if (lemode == 1) {
/* Elastc mode never needed */
	lvlinf = 1;
    }
    lureq = 0;
    typelu = 3;
    luok = TRUE_;
    done = fponly;
/*     ================================================================== */
/* +    while (.not. done  .and.  LUok) do */
L600:
    if (! done && luok) {
	newlu = lureq > 0;
	needx = needlu;
/* -------------------------------------------------------------- */
/* LP with objective row in A. */
/* -------------------------------------------------------------- */
	iw[223] = 0;
	iw[225] = 0;
	s5lp_(&inform__, prob, probtag, &elastc, &subopt, (U_fp)qplog, &
		needlu, &needx, m, n, nb, ndegen, itqp, itqpmax, itn, &lemode,
		 &lvlinf, &prtlvl, minimz, iobj, sclobj, objadd, tolfp, tolqp,
		 tolx, ninf, sinf, wtinf, pinorm, &rgnorm, ne, nloca, &loca[1]
		, &inda[1], &acol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[
		1], &kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &
		gbs[1], &pi[1], &rc[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], 
		&xfreez[1], &iy[1], &iy1[1], &y[1], &y1[1], cw + 8, lencw, &
		iw[1], leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8);
/* Check for trouble in s5LP.  Here are the possibilities: */
	badlu = inform__ > 0;
/* Fatal LU error */
	found = inform__ == 0;
/* LP solution found */
	infeas = inform__ == -1;
/* LP is infeasible */
	unbnded = inform__ == -2;
/* LP is unbounded */
	itnbig = inform__ == -3;
/* Too many iterations */
	if (badlu) {
	    *iexit = inform__;
	    goto L900;
	}
	done = found || unbnded || itnbig;
	if (done) {
/* =========================================================== */
/* Relax, we are finished. */
/* =========================================================== */
	} else {
/* ----------------------------------------------------------- */
/* The non-elastics are infeasible. This should not happen. */
/* Phase 1 has already found a feasible point for the */
/* nonelastics, so the basis must be ill-conditioned. */
/* Refactorize with tighter tols and restart at the known */
/* feasible point.  Reduce the feasibility tol to try and */
/* prevent repeats. */
/* ----------------------------------------------------------- */
	    s_wsfi(&io___255);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    s2trylu_(itn, &c__22, &ns, &lureq, &luok, &typelu, &iw[1], leniw, 
		    &rw[1], lenrw);
	}
	goto L600;
    }
/* +    end while */
/*     ================================================================== */
L800:
    if (found) {
	if (*ninf == 0) {
	    if (fponly) {
		*iexit = 2;
/* Feasible */
	    } else {
		*iexit = 1;
/* optimal */
	    }
	} else if (lemode == 0) {
	    *iexit = 11;
/* infeasible linear constraints */
	} else {
	    *iexit = 14;
/* infeasibilites minimized */
	}
    } else if (infeas) {
	*iexit = 11;
/* infeasible nonelastics */
    } else if (unbnded) {
	*iexit = 21;
/* unbounded */
    } else if (itnbig) {
	*iexit = 31;
/* too many iterations */
    }
L900:
    return 0;
} /* s5slp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5sLP */
/* Subroutine */ int s5sqp_(integer *iexit, U_fp hprod, U_fp hprod1, U_fp 
	qplog, logical *gotr, integer *prob, integer *lenr, integer *m, 
	integer *maxs, integer *mbs, integer *n, integer *nb, integer *ndegen,
	 integer *hvcalls, integer *ngqp0, integer *ngqp, integer *ngobj0, 
	integer *ngobj, integer *nnh0, integer *nnh, integer *ns, integer *
	itqp, integer *itqpmax, integer *itn, integer *minimz, integer *iobj, 
	doublereal *sclobj, doublereal *objadd, doublereal *objqp, doublereal 
	*tolfp, doublereal *tolqp, doublereal *tolx, integer *ninf, 
	doublereal *sinf, doublereal *wtinf, doublereal *pinorm, integer *ne, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *hetype, integer *hestat, integer *hfeas, integer *hs, 
	integer *kbs, doublereal *ascale, doublereal *bl, doublereal *bu, 
	doublereal *blbs, doublereal *bubs, doublereal *gbs, doublereal *gobj,
	 doublereal *gqp, doublereal *hdx, doublereal *pbs, doublereal *pi, 
	doublereal *r__, doublereal *rc, doublereal *rg, integer *nrhs0, 
	integer *nrhs, doublereal *rhs, integer *lenx0, integer *nx0, 
	doublereal *x0, doublereal *x, doublereal *xbs, doublereal *xfreez, 
	integer *iy, integer *iy1, doublereal *y, doublereal *y1, doublereal *
	y2, char *cu, integer *lencu, integer *iu, integer *leniu, doublereal 
	*ru, integer *lenru, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in feasibility\002,\002 phase\002)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Feasible non-elastic"
	    "s\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Feasible constraints"
	    "\002)";
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP optimality\002,\002 phase\002)";
    static char fmt_1600[] = "(\002 Itn\002,i7,\002: Indefinite reduced QP H"
	    "essian --\002,\002 refactor the basis\002)";
    static char fmt_1700[] = "(\002 Itn\002,i7,\002: Large reduced gradien"
	    "t\002)";
    static char fmt_1800[] = "(\002 Itn\002,i7,\002: Ill-conditioned null sp"
	    "ace\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer itqptargt;
    static doublereal eps;
    static char str[120];
    extern /* Subroutine */ int s5qp_(integer *, integer *, char *, logical *,
	     integer *, U_fp, U_fp, U_fp, logical *, logical *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static logical badz, done, newb, luok;
    static integer znnh;
    static logical badlu, badzg;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static doublereal objfp;
    static logical nsbig, needx;
    static doublereal flmax;
    static logical found;
    static integer lureq;
    static logical newlu;
    static integer nswap, zngqp;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static integer lemod0, lvlin0;
    extern /* Subroutine */ int s5getr_(integer *, U_fp, U_fp, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static doublereal hcndbd;
    static logical indefh;
    static integer gotfac;
    static doublereal zcndbd;
    static logical infeas;
    static integer lemode;
    static logical elastc, itnbig, needlu;
    static integer inform__, lvlinf;
    static doublereal targth, plinfy;
    static logical getopt;
    static doublereal rgnorm;
    static logical fponly;
    static doublereal targtz;
    static integer subopt, mnrprt, prtlvl, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static logical unbnded;
    extern /* Subroutine */ int s2trylu_(integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *);
    static logical indefqp, weakmin;
    static char probtag[20];

    /* Fortran I/O blocks */
    static icilist io___296 = { 0, str, 0, fmt_1000, 120, 1 };
    static icilist io___297 = { 0, str, 0, fmt_1200, 120, 1 };
    static icilist io___298 = { 0, str, 0, fmt_1300, 120, 1 };
    static icilist io___305 = { 0, str, 0, fmt_1100, 120, 1 };
    static icilist io___306 = { 0, str, 0, fmt_1600, 120, 1 };
    static icilist io___307 = { 0, str, 0, fmt_1700, 120, 1 };
    static icilist io___308 = { 0, str, 0, fmt_1800, 120, 1 };


/* ================================================================= */
/* s5sQP   solves the current problem.  An initial basis is assumed */
/* to be specified by nS, hs, x and the superbasic parts of kBS. */
/* In particular, there must be nS values hs(j) = 2, and the */
/* corresponding j's must be listed in kBS(m+1) thru kBS(m+nS). */
/* The ordering in kBS matches the reduced Hessian R (if any). */

/* 05 Oct 1994: First version of s5sQP. */
/* 06 Aug 1996: Min Sum option added. */
/* 14 Jul 1997: Thread-safe version. */
/* 28 Apr 2002: Updated to reflect LU changes. */
/* 02 Aug 2003: snEXIT and snPRNT adopted. */
/* 08 Mar 2004: gotFac implemented. */
/* 16 May 2006: Explicit target itQP added. */
/* 27 Apr 2011: Finds indefinite Hessian with one BS factorize. */
/* 28 Apr 2011: Do one refactorize after detecting indefiniteness. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* >0 => Mnr heading for iPrint */
/* >0 => Mnr heading for iSumm */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;
    --pi;
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
    --xfreez;
    --x;
    --rc;
    --pbs;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --gqp;
    --gobj;
    --hdx;
    --acol;
    --inda;
    --loca;
    --rhs;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lemode = iw[56];
/* >0    => use elastic mode */
    lvlinf = iw[73];
/* Elastic option */
    mnrprt = iw[93];
/* Minor print level */
    gotfac = iw[230];
/* >0 => Save the LU factors */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    flmax = rw[8];
/* est. of the largest pos. real */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    prtlvl = mnrprt;
    subopt = -1;
    plinfy = flmax;
    *objqp = 0.;
    targth = hcndbd;
    targtz = zcndbd;
    itqptargt = *itqpmax;
    s_copy(probtag, "linear constraints", (ftnlen)20, (ftnlen)18);
    getopt = *prob == 1 || *prob == 2;
    fponly = *prob == 0;
    elastc = lemode == 2;
    iload_(nb, &c__0, &hestat[1], &c__1);
/* ----------------------------------------------------------------- */
/* Find a feasible point. */
/* If the constraints are infeasible, minimize the sum of the */
/* elastic variables, subject to keeping the non-elastic variables */
/* feasible.  Elastic variables can move outside their bounds. */
/* ----------------------------------------------------------------- */
    if (elastc) {
/* make the nonelastics feasible */
	lemod0 = lemode;
	lvlin0 = lvlinf;
    } else {
/* make everything feas or min sum */
	lemod0 = 1;
	lvlin0 = 2;
    }
    *gotr = FALSE_;
    zngqp = 0;
/* No objective in phase 1 */
    znnh = 0;
/* No Hessian either */
    iw[223] = 1;
/* Refresh print   heading. */
    iw[224] = 1;
/* Refresh summary heading */
    if (gotfac == 0) {
	lureq = 1;
    } else {
	lureq = 0;
    }
    typelu = 2;
    luok = TRUE_;
    done = FALSE_;
/*     ================================================================== */
/* +    while (.not. done  .and.  LUok) do */
L500:
    if (! done && luok) {
	needlu = lureq > 0;
	s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, &prtlvl, &
		lureq, m, mbs, n, nb, nnh, ns, &nswap, ne, nloca, &loca[1], &
		inda[1], &acol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &blbs[1], 
		&bubs[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], &iy[1], &iy1[
		1], &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    goto L900;
	}
	needx = needlu;
	s5qp_(&inform__, &c__0, probtag, &elastc, &subopt, (U_fp)hprod, (U_fp)
		hprod1, (U_fp)qplog, gotr, &needlu, &typelu, &needx, lenr, m, 
		maxs, mbs, n, nb, ndegen, hvcalls, ngqp0, &zngqp, ngobj0, 
		ngobj, nnh0, &znnh, ns, itqp, itqpmax, &itqptargt, itn, &
		lemod0, &lvlin0, &prtlvl, minimz, iobj, sclobj, objadd, &
		objfp, &targth, &targtz, tolfp, tolqp, tolx, ninf, sinf, 
		wtinf, pinorm, &rgnorm, ne, nloca, &loca[1], &inda[1], &acol[
		1], &hetype[1], &hestat[1], &hfeas[1], &hs[1], &kbs[1], &
		ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &gbs[1], &gobj[
		1], &gqp[1], &hdx[1], &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1]
		, nrhs0, nrhs, &rhs[1], lenx0, nx0, &x0[1], &x[1], &xbs[1], &
		xfreez[1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, 
		lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* Check for trouble in s5QP.  Here are the possibilities: */
	badlu = inform__ > 0;
/* Fatal LU error */
	found = inform__ == 0;
/* QP solution found */
	infeas = inform__ == -1;
/* QP is infeasible */
	unbnded = inform__ == -2;
/* QP is unbounded */
	itnbig = inform__ == -3;
/* Too many iterations */
	if (badlu) {
	    *iexit = inform__;
/* Fatal LU error */
	    goto L900;
	}
	done = found || infeas || itnbig;
	if (! done) {
/* =========================================================== */
/* Trouble. */
/* Looks like the phase 1 problem is unbounded, which can */
/* occur only if a bad basis has given a large direction. */
/* =========================================================== */
/* Assume infeasible and factor with tighter tols. */
	    infeas = TRUE_;
	    s_wsfi(&io___296);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], leniw, &
		    rw[1], lenrw);
	}
	goto L500;
    }
/* +    end while */
/*     ================================================================== */
    if (! found) {
	goto L800;
    }
/* Itns or infs */
    if (*ninf > 0 && lemode != 1) {
	goto L800;
    }
/* Infeas */
    if (getopt) {
/*        call dcopy ( nb, x, 1, xFreez, 1 ) ! Save feasible nonelastics */
/*        if (PrtLvl .ge. 1  .and.  PrtLvl .lt. 10) then */
	if (prtlvl >= 10) {
	    if (elastc) {
		s_wsfi(&io___297);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    } else {
		s_wsfi(&io___298);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    }
	}
    }
/* ================================================================= */
/* Get optimal. */
/* ================================================================= */
    if (elastc) {
/* Relax */
    } else if (lemode == 1) {
/* Elastc mode never needed */
	lvlinf = 1;
    }
    iw[223] = 1;
/* Refresh print heading. */
    iw[224] = 1;
    indefqp = FALSE_;
    lureq = 0;
    typelu = 3;
    luok = TRUE_;
    done = fponly;
/*     ================================================================== */
/* +    while (.not. done  .and.  LUok) do */
L600:
    if (! done && luok) {
/* -------------------------------------------------------------- */
/* Compute and factorize the initial Z'HZ. */
/* The basis is refactorized if necessary. */
/* -------------------------------------------------------------- */
	if (! (*gotr) || lureq > 0) {
	    s5getr_(iexit, (U_fp)hprod, (U_fp)hprod1, hvcalls, gotr, &typelu, 
		    &lureq, itn, lenr, m, mbs, n, nb, nnh, ns, &prtlvl, 
		    minimz, iobj, &targth, &targtz, ne, nloca, &loca[1], &
		    inda[1], &acol[1], &hs[1], &kbs[1], &bl[1], &bu[1], &blbs[
		    1], &bubs[1], &r__[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[
		    1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu,
		     &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	    if (*iexit != 0) {
		goto L900;
	    }
	}
/* -------------------------------------------------------------- */
/* Solve the QP. */
/* -------------------------------------------------------------- */
	lureq = 0;
	s5qp_(&inform__, prob, probtag, &elastc, &subopt, (U_fp)hprod, (U_fp)
		hprod1, (U_fp)qplog, gotr, &needlu, &typelu, &needx, lenr, m, 
		maxs, mbs, n, nb, ndegen, hvcalls, ngqp0, ngqp, ngobj0, ngobj,
		 nnh0, nnh, ns, itqp, itqpmax, &itqptargt, itn, &lemode, &
		lvlinf, &prtlvl, minimz, iobj, sclobj, objadd, objqp, &targth,
		 &targtz, tolfp, tolqp, tolx, ninf, sinf, wtinf, pinorm, &
		rgnorm, ne, nloca, &loca[1], &inda[1], &acol[1], &hetype[1], &
		hestat[1], &hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1], &
		bu[1], &blbs[1], &bubs[1], &gbs[1], &gobj[1], &gqp[1], &hdx[1]
		, &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], nrhs0, nrhs, &rhs[
		1], lenx0, nx0, &x0[1], &x[1], &xbs[1], &xfreez[1], &iy[1], &
		iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], leniu, &
		ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)20, (ftnlen)8, (ftnlen)8);
/* Check for trouble in s5QP.  Here are the possibilities: */
	badlu = inform__ > 0;
/* Fatal LU error */
	found = inform__ == 0;
/* QP solution found */
	infeas = inform__ == -1;
/* QP is infeasible */
	unbnded = inform__ == -2;
/* QP is unbounded */
	itnbig = inform__ == -3;
/* Too many iterations */
	weakmin = inform__ == -4;
/* Weak QP minimizer */
	nsbig = inform__ == -5;
/* Too many superbasics */
	indefh = inform__ == -6;
/* QP Hessian not positive semidefinite */
	badzg = inform__ == -7;
/* Z'g could not be made small enough */
	badz = inform__ == -8;
/* Ill-conditioned Z */
	if (badlu) {
	    *iexit = inform__;
	    goto L900;
	}
	indefqp = indefh && indefqp;
	done = found || unbnded || infeas || nsbig || indefqp || itnbig || 
		weakmin;
	if (done) {
/* ----------------------------------------------------------- */
/* Relax, we are finished. */
/* ----------------------------------------------------------- */
	} else {
/* ----------------------------------------------------------- */
/* Numerical trouble in s5QP. */
/* ----------------------------------------------------------- */
	    *gotr = FALSE_;
	    if (infeas) {
/* -------------------------------------------------------- */
/* The nonelastics are infeasible. This should not happen. */
/* Phase 1 has already found a feasible point for the */
/* nonelastics, so the basis must be ill-conditioned. */
/* Refactorize with tighter tols and restart at the known */
/* feasible point.  Reduce the feasibility tol to try and */
/* prevent repeats. */
/* ------------------------------------------------------- */
		s_wsfi(&io___305);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
		s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		elastc = FALSE_;
/*              call dcopy  ( nb, xFreez, 1, x, 1 ) */
	    } else if (indefh) {
/* -------------------------------------------------------- */
/* Indefinite Z'HZ.  Could be an ill-conditioned Z'HZ. */
/* -------------------------------------------------------- */
		s_wsfi(&io___306);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
		indefqp = indefh;
		targth = 1. / (eps * eps);
		s2trylu_(itn, &c__26, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
	    } else if (badzg) {
/* -------------------------------------------------------- */
/* Large Z'g.        Probably an ill-conditioned Z. */
/* -------------------------------------------------------- */
		s_wsfi(&io___307);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		s2trylu_(itn, &c__21, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
	    } else if (badz) {
/* -------------------------------------------------------- */
/* condZ > targtZ  while forming Z'HZ for a freq. check. */
/* Refactorize B, possibly with a reduced factor tol. If */
/* the factor tol is already tight, accept Z, however bad. */
/* -------------------------------------------------------- */
		s_wsfi(&io___308);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
		s2trylu_(itn, &c__25, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    targtz = plinfy;
		    luok = TRUE_;
		}
	    }
	}
	goto L600;
    }
/* +    end while */
/*     ================================================================== */
L800:
    if (found) {
	if (*ninf == 0) {
	    if (fponly) {
		*iexit = 2;
/* feasible */
	    } else {
		*iexit = 1;
/* optimal */
	    }
	} else {
	    *iexit = 14;
/* optimal in elastic mode */
	}
    } else if (infeas) {
	*iexit = 11;
/* infeasible nonelastics */
    } else if (unbnded) {
	*iexit = 21;
/* unbounded */
    } else if (itnbig) {
	*iexit = 31;
/* too many iterations */
    } else if (weakmin) {
	*iexit = 4;
/* weak minimizer */
    } else if (nsbig) {
	*iexit = 33;
/* too many superbasics */
    } else if (indefqp) {
	*iexit = 53;
/* Hessian indefinite */
    } else {
	*iexit = 41;
/* Current point cannot be improved */
    }
L900:
    return 0;
} /* s5sqp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5sQP */
/* Subroutine */ int s5sqn_(integer *iexit, U_fp hprod, U_fp hprod1, U_fp 
	qplog, logical *gotr, integer *prob, integer *lenr, integer *m, 
	integer *maxs, integer *mbs, integer *n, integer *nb, integer *ndegen,
	 integer *hvcalls, integer *ngqp0, integer *ngqp, integer *ngobj0, 
	integer *ngobj, integer *nnh0, integer *nnh, integer *ns, integer *
	itqp, integer *itqpmax, integer *itn, integer *minimz, integer *iobj, 
	doublereal *sclobj, doublereal *objadd, doublereal *objqp, doublereal 
	*tolfp, doublereal *tolqp, doublereal *tolx, integer *ninf, 
	doublereal *sinf, doublereal *wtinf, doublereal *pinorm, integer *ne, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *hetype, integer *hestat, integer *hfeas, integer *hs, 
	integer *kbs, doublereal *ascale, doublereal *bl, doublereal *bu, 
	doublereal *blbs, doublereal *bubs, doublereal *gbs, doublereal *gobj,
	 doublereal *gqp, doublereal *hdx, doublereal *pbs, doublereal *pi, 
	doublereal *r__, doublereal *rc, doublereal *rg, doublereal *rg2, 
	integer *nrhs0, integer *nrhs, doublereal *rhs, integer *lenx0, 
	integer *nx0, doublereal *x0, doublereal *x, doublereal *xbs, 
	doublereal *xfreez, integer *iy, integer *iy1, doublereal *y, 
	doublereal *y1, doublereal *y2, char *cu, integer *lencu, integer *iu,
	 integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in feasibility\002,\002 phase\002)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Feasible non-elastic"
	    "s\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Feasible constraints"
	    "\002)";
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QN optimality\002,\002 phase\002)";
    static char fmt_1800[] = "(\002 Itn\002,i7,\002: Ill-conditioned null sp"
	    "ace\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer itqptargt;
    static char str[120];
    extern /* Subroutine */ int s5qn_(integer *, integer *, char *, logical *,
	     integer *, U_fp, U_fp, U_fp, logical *, logical *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical badz, done, newb, luok, badlu, badzg;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static doublereal objfp;
    static logical nsbig, needx;
    static doublereal flmax;
    static logical found;
    static integer lureq;
    static logical newlu;
    static integer nswap, zngqp;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static integer lemod0;
    static logical elast0;
    static integer lvlin0;
    static logical indefh;
    static integer gotfac;
    static doublereal zcndbd;
    static logical infeas;
    static integer lemode;
    static logical elastc, itnbig, needlu;
    static doublereal condhz;
    static integer inform__, lvlinf;
    static doublereal plinfy;
    static logical getopt, fponly;
    static doublereal targtz;
    static integer subopt, mnrprt, prtlvl, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static logical unbnded;
    extern /* Subroutine */ int s2trylu_(integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *);
    static logical weakmin;
    static char probtag[20];
    static logical zitnbig;

    /* Fortran I/O blocks */
    static icilist io___346 = { 0, str, 0, fmt_1000, 120, 1 };
    static icilist io___347 = { 0, str, 0, fmt_1200, 120, 1 };
    static icilist io___348 = { 0, str, 0, fmt_1300, 120, 1 };
    static icilist io___355 = { 0, str, 0, fmt_1100, 120, 1 };
    static icilist io___356 = { 0, str, 0, fmt_1800, 120, 1 };


/* ================================================================= */
/* s5sQN   solves the current problem.  An initial basis is assumed */
/* to be specified by nS, hs, x and the superbasic parts of kBS. */
/* In particular, there must be nS values hs(j) = 2, and the */
/* corresponding j's must be listed in kBS(m+1) thru kBS(m+nS). */
/* The ordering in kBS matches the reduced Hessian R (if any). */

/* 01 Nov 2003: First version of s5sQN based on s5sQP. */
/* 08 Mar 2004: gotFac, gotHes implemented. */
/* 16 May 2006: Explicit target itQP added. */
/* 28 Apr 2011: Some cosmetic restructuring. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* >0 => Mnr heading for iPrint */
/* >0 => Mnr heading for iSumm */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;
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
    --xfreez;
    --x;
    --rc;
    --pbs;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --gqp;
    --gobj;
    --hdx;
    --acol;
    --inda;
    --loca;
    --rhs;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lemode = iw[56];
/* >0    => use elastic mode */
    lvlinf = iw[73];
/* Elastic option */
    mnrprt = iw[93];
/* Minor print level */
    gotfac = iw[230];
/* >0 => Save the LU factors */
    flmax = rw[8];
/* est. of the largest pos. real */
    zcndbd = rw[86];
/* bound on the condition of Z */
    prtlvl = mnrprt;
    subopt = -1;
    plinfy = flmax;
    *objqp = 0.;
    targtz = zcndbd;
    condhz = 0.;
    itqptargt = *itqpmax;
    s_copy(probtag, "linear constraints", (ftnlen)20, (ftnlen)18);
    getopt = *prob == 1 || *prob == 2;
    fponly = *prob == 0;
    elastc = lemode == 2;
    iload_(nb, &c__0, &hestat[1], &c__1);
/* ----------------------------------------------------------------- */
/* Find a feasible point. */
/* If the constraints are infeasible, minimize the sum of the */
/* elastic variables, subject to keeping the non-elastic variables */
/* feasible.  Elastic variables can move outside their bounds. */
/* ----------------------------------------------------------------- */
    lvlin0 = 2;
/* local value of lvlInf */
    lemod0 = 1;
/* local value of lEmode */
    elast0 = FALSE_;
/* local value of Elastc */
    zngqp = 0;
/* No objective in phase 1 */
    needx = FALSE_;
/* needLU will ask for new factors. */
    if (gotfac == 0) {
	lureq = 1;
    } else {
	lureq = 0;
    }
    typelu = 2;
    luok = TRUE_;
    done = FALSE_;
/*     ================================================================== */
/* +    while (.not. done  .and.  LUok) do */
L500:
    if (! done && luok) {
	needlu = lureq > 0;
	s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, &prtlvl, &
		lureq, m, mbs, n, nb, nnh, ns, &nswap, ne, nloca, &loca[1], &
		inda[1], &acol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &blbs[1], 
		&bubs[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], &iy[1], &iy1[
		1], &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    goto L900;
	}
	needx = needlu;
	s5qn_(&inform__, &c__0, probtag, &elast0, &subopt, (U_fp)hprod, (U_fp)
		hprod1, (U_fp)qplog, gotr, &needlu, &typelu, &needx, lenr, m, 
		maxs, mbs, n, nb, ndegen, hvcalls, ngqp0, &zngqp, ngobj0, 
		ngobj, nnh0, nnh, ns, itqp, itqpmax, &itqptargt, itn, &lemod0,
		 &lvlin0, &prtlvl, minimz, iobj, sclobj, objadd, &objfp, &
		condhz, &targtz, tolfp, tolqp, tolx, ninf, sinf, wtinf, 
		pinorm, ne, nloca, &loca[1], &inda[1], &acol[1], &hetype[1], &
		hestat[1], &hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1], &
		bu[1], &blbs[1], &bubs[1], &gbs[1], &gobj[1], &gqp[1], &hdx[1]
		, &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], &rg2[1], nrhs0, 
		nrhs, &rhs[1], lenx0, nx0, &x0[1], &x[1], &xbs[1], &xfreez[1],
		 &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1]
		, leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* Check for trouble in s5QN.  Here are the possibilities: */
	badlu = inform__ > 0;
/* Fatal LU error */
	found = inform__ == 0;
/* QP solution found */
	infeas = inform__ == -1;
/* QP is infeasible */
	unbnded = inform__ == -2;
/* QP is unbounded */
	itnbig = inform__ == -3;
/* Too many iterations */
	if (badlu) {
	    *iexit = inform__;
/* Fatal LU error */
	    goto L900;
	}
	done = found || infeas || itnbig;
	if (! done) {
/* =========================================================== */
/* Trouble. */
/* The phase 1 was unbounded, which can only occur if a bad */
/* basis has given a large search direction. */
/* =========================================================== */
/* Assume infeasible constraints and factor with tighter tols */
	    infeas = TRUE_;
	    s_wsfi(&io___346);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], leniw, &
		    rw[1], lenrw);
	}
	goto L500;
    }
/* +    end while */
/*     ================================================================== */
    if (! found) {
	goto L800;
    }
/* Itns or infeasible */
    if (getopt) {
/*        call dcopy ( nb, x, 1, xFreez, 1 ) ! Save feasible nonelastics */
/*        if (PrtLvl .ge. 1  .and.  PrtLvl .lt. 10) then */
	if (prtlvl >= 10) {
	    if (elastc) {
		s_wsfi(&io___347);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    } else {
		s_wsfi(&io___348);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
	    }
	}
    }
/* ================================================================= */
/* Get optimal. */
/* ================================================================= */
    lureq = 0;
    typelu = 3;
    luok = TRUE_;
    done = fponly;
/*     ================================================================== */
/* +    while (.not. done  .and.  LUok) do */
L600:
    if (! done && luok) {
/* -------------------------------------------------------------- */
/* Refactorize the basis if necessary. */
/* -------------------------------------------------------------- */
	needlu = lureq > 0;
	if (needlu) {
	    s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, &
		    prtlvl, &lureq, m, mbs, n, nb, nnh, ns, &nswap, ne, nloca,
		     &loca[1], &inda[1], &acol[1], &kbs[1], &hs[1], &bl[1], &
		    bu[1], &blbs[1], &bubs[1], nrhs0, nrhs, &rhs[1], &x[1], &
		    xbs[1], &iy[1], &iy1[1], &y[1], &y1[1], &iw[1], leniw, &
		    rw[1], lenrw);
	    if (*iexit != 0) {
		goto L900;
	    }
	}
/* How do we update R if the superbasics change? */
/* -------------------------------------------------------------- */
/* Solve the QP using a quasi-Newton method. */
/* -------------------------------------------------------------- */
	iw[223] = 1;
/* Refresh print heading. */
	iw[224] = 1;
	lureq = 0;
	if (elastc) {
	    lvlinf = 1;
/* W1 = 1.0, W2 = wtInf  Elastic Phase 2 composite */
	}
	s5qn_(&inform__, prob, probtag, &elastc, &subopt, (U_fp)hprod, (U_fp)
		hprod1, (U_fp)qplog, gotr, &needlu, &typelu, &needx, lenr, m, 
		maxs, mbs, n, nb, ndegen, hvcalls, ngqp0, ngqp, ngobj0, ngobj,
		 nnh0, nnh, ns, itqp, itqpmax, &itqptargt, itn, &lemode, &
		lvlinf, &prtlvl, minimz, iobj, sclobj, objadd, objqp, &condhz,
		 &targtz, tolfp, tolqp, tolx, ninf, sinf, wtinf, pinorm, ne, 
		nloca, &loca[1], &inda[1], &acol[1], &hetype[1], &hestat[1], &
		hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[
		1], &bubs[1], &gbs[1], &gobj[1], &gqp[1], &hdx[1], &pbs[1], &
		pi[1], &r__[1], &rc[1], &rg[1], &rg2[1], nrhs0, nrhs, &rhs[1],
		 lenx0, nx0, &x0[1], &x[1], &xbs[1], &xfreez[1], &iy[1], &iy1[
		1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], leniu, &ru[
		1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)20, (ftnlen)8, (ftnlen)8);
/* Check for trouble in s5QN.  Here are the possibilities: */
	badlu = inform__ > 0;
/* Fatal LU error */
	found = inform__ == 0;
/* QP solution found */
	infeas = inform__ == -1;
/* QP is infeasible */
	unbnded = inform__ == -2;
/* QP is unbounded */
	itnbig = inform__ == -3;
/* Too many iterations */
	weakmin = inform__ == -4;
/* Void, should not happen */
	nsbig = inform__ == -5;
/* Too many superbasics */
	indefh = inform__ == -6;
/* Void, should not happen */
	badzg = inform__ == -7;
/* Void, should not happen */
	badz = inform__ == -8;
/* Ill-conditioned Z */
	zitnbig = inform__ == -9;
/* Too many subspace iterations */
	if (badlu) {
	    *iexit = inform__;
	    goto L900;
	}
	done = found || unbnded || infeas || nsbig || itnbig || weakmin || 
		indefh || badzg || zitnbig;
	if (done) {
/* ----------------------------------------------------------- */
/* Relax, we are finished. */
/* ----------------------------------------------------------- */
	} else {
/* ----------------------------------------------------------- */
/* Numerical trouble in s5QN. */
/* ----------------------------------------------------------- */
	    *gotr = FALSE_;
	    if (infeas) {
/* -------------------------------------------------------- */
/* The nonelastics are infeasible. This should not happen. */
/* Phase 1 has already found a feasible point for the */
/* nonelastics, so the basis must be ill-conditioned. */
/* Refactorize with tighter tols and restart at the known */
/* feasible point.  Reduce the feasibility tol to try and */
/* prevent repeats. */
/* -------------------------------------------------------- */
		s_wsfi(&io___355);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
		s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		elastc = FALSE_;
/*              call dcopy  ( nb, xFreez, 1, x, 1 ) */
	    } else if (badz) {
/* -------------------------------------------------------- */
/* Numerical error.  Probably a bad basis. */
/* Refactorize B, possibly with a reduced factor tol. If */
/* the factor tol is already tight, accept Z, however bad. */
/* -------------------------------------------------------- */
		s_wsfi(&io___356);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)120);
		s2trylu_(itn, &c__25, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    targtz = plinfy;
		    luok = TRUE_;
		}
	    }
	}
	goto L600;
    }
/* +    end while */
/*     ================================================================== */
L800:
    if (found) {
	if (*ninf == 0) {
	    *iexit = 1;
/* optimal */
	} else {
	    *iexit = 14;
/* optimal */
	}
    } else if (infeas) {
	*iexit = 11;
/* infeasible nonelastics */
    } else if (unbnded) {
	*iexit = 21;
/* unbounded */
    } else if (itnbig) {
	*iexit = 31;
/* too many iterations */
    } else if (nsbig) {
	*iexit = 33;
/* too many superbasics */
    } else {
	*iexit = 41;
/* Current point cannot be improved */
    }
L900:
    return 0;
/* L1600: */
/* L1700: */
} /* s5sqn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5sQN */
/* Subroutine */ int s5stat_(integer *status, integer *iw, integer *leniw)
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
    static icilist io___358 = { 0, str, 0, fmt_9999, 80, 1 };


/* ================================================================= */
/* s5Stat fetches the call-status for the sqOpt user-defined */
/* matrix-vector product. */

/* 16 Jun 2008: First version of s5Stat. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* QP user-routine call-status */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (iw[235] == 0) {
/* Standard call */
	*status = 0;
    } else if (iw[235] < 0) {
/* First call */
	*status = 1;
	iw[235] = 0;
    } else if (iw[235] >= 2) {
/* Last orders please */
	*status = iw[235];
	iw[235] = -1;
    } else {
	*status = iw[235];
	s_wsfi(&io___358);
	do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s5stat_ */

