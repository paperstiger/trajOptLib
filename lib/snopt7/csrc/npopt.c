/* ./src/npopt.f -- translated by f2c (version 20100827).
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

static integer c__500 = 500;
static integer c__130 = 130;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  npOpt.f --- the NPSOL wrapper for SNOPT. */

/*     npOpt   npKerN */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int npopt_(integer *n, integer *nclin, integer *ncnln, 
	integer *lda, integer *ldcj, integer *ldh, doublereal *a, doublereal *
	bl, doublereal *bu, U_fp funcon, U_fp funobj, integer *info, integer *
	majits, integer *istate, doublereal *ccon, doublereal *cjac, 
	doublereal *cmul, doublereal *objf, doublereal *grad, doublereal *
	hess, doublereal *x, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* System generated locals */
    integer a_dim1, a_offset, cjac_dim1, cjac_offset, hess_dim1, hess_offset;

    /* Local variables */
    extern /* Subroutine */ int snlog_(), snlog2_();
    extern /* Subroutine */ int npkern_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, U_fp, U_fp, U_fp, U_fp, U_fp, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    extern /* Subroutine */ int snstop_();

/*     ================================================================== */
/*     npOpt  solves the nonlinear programming problem */

/*            minimize                   f(x) */

/*                                    (      x  ) */
/*            subject to    bl  .le.  (    A*x  )  .le.  bu */
/*                                    ( cCon(x) ) */

/*     where  f(x)  is a smooth scalar function,  A  is a constant matrix */
/*     and  cCon(x)  is a vector of smooth nonlinear functions. */
/*     The feasible region is defined by a mixture of linear and */
/*     nonlinear equality or inequality constraints on  x. */

/*     The calling sequence of NPOPT and the user-defined functions */
/*     funcon and funobj are identical to those of the dense code NPSOL */
/*     (see the User's Guide for NPSOL (Version 4.0): a Fortran Package */
/*     for Nonlinear Programming, Systems Optimization Laboratory Report */
/*     SOL 86-2, Department of Operations Research, Stanford University, */
/*     1986.) */

/*     The dimensions of the problem are... */

/*     n        the number of variables (dimension of  x), */

/*     nclin    the number of linear constraints (rows of the matrix  A), */

/*     ncnln    the number of nonlinear constraints (dimension of  c(x)), */

/*     NPOPT  is maintained by Philip E. Gill, */
/*     Dept of Mathematics, University of California, San Diego. */

/*     LUSOL is maintained by Michael A. Saunders, */
/*     Systems Optimization Laboratory, */
/*     Dept of Management Science & Engineering, Stanford University. */

/*     22 Mar 1997: First   version of npOpt. */
/*     31 Jul 2003: snEXIT and snPRNT adopted. */
/*     15 Oct 2004: snSTOP adopted. */
/*     01 Sep 2007: Sticky parameters added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --grad;
    --cmul;
    --istate;
    --bu;
    --bl;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    cjac_dim1 = *ldcj;
    cjac_offset = 1 + cjac_dim1;
    cjac -= cjac_offset;
    hess_dim1 = *ldh;
    hess_offset = 1 + hess_dim1;
    hess -= hess_offset;
    --ccon;
    --iw;
    --rw;

    /* Function Body */
    npkern_(n, nclin, ncnln, lda, ldcj, ldh, &a[a_offset], &bl[1], &bu[1], (
	    U_fp)funcon, (U_fp)funobj, (U_fp)snlog_, (U_fp)snlog2_, (U_fp)
	    snstop_, info, majits, &istate[1], &ccon[1], &cjac[cjac_offset], &
	    cmul[1], objf, &grad[1], &hess[hess_offset], &x[1], &iw[1], leniw,
	     &rw[1], lenrw);
    return 0;
} /* npopt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npOpt */
/* Subroutine */ int npkern_(integer *n, integer *nclin, integer *ncnln, 
	integer *lda, integer *ldcj, integer *ldh, doublereal *a, doublereal *
	bl, doublereal *bu, U_fp funcon, U_fp funobj, U_fp snlog, U_fp snlog2,
	 U_fp snstop, integer *info, integer *majits, integer *istate, 
	doublereal *ccon, doublereal *cjac, doublereal *cmul, doublereal *
	objf, doublereal *grad, doublereal *hess, doublereal *x, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer a_dim1, a_offset, cjac_dim1, cjac_offset, hess_dim1, hess_offset, 
	    i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, m, nb, ne;
    static char cw[8*500];
    static integer ns, lx, lx0, lbl, lrc, lbu, nnh, lpi, lhs, lkx, nkx;
    static char str[80], str2[80];
    static doublereal fobj;
    static integer iobj, ninf, ncon;
    static doublereal sinf;
    static integer lenr, maxr, maxs;
    static logical gotr;
    extern /* Subroutine */ int s0fgn_();
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s8map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s3inn_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *);
    static integer nnjac, lgobj, ngobj, lindj, lfcon, lgcon, ljcol, llocg, 
	    llocj, nlocj, nlocg, nncol, nncon, mincw, maxcw, nnobj, miniw, 
	    maxiw;
    extern /* Subroutine */ int s2mem0_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, ftnlen),
	     icopy_(integer *, integer *, integer *, integer *, integer *), 
	    dcopy_(integer *, doublereal *, integer *, doublereal *, integer *
	    );
    static integer minrw, start, maxrw;
    static doublereal xnorm;
    extern /* Subroutine */ int s1file_(integer *, integer *, integer *), 
	    s2bmap_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    static integer nncon0;
    extern /* Subroutine */ int s3argn_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, char *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, ftnlen, ftnlen), s8gloc_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     s1time_(integer *, integer *, integer *, integer *, doublereal *,
	     integer *), s3inin_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    extern doublereal dnrm1s_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int s8dflt_(integer *, integer *, integer *, 
	    integer *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), s1perm_(integer *, integer *), 
	    s3prtb_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    static doublereal objadd;
    extern /* Subroutine */ int s3prtn_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s3outn_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s8solv_(integer *, char *, integer *, 
	    U_fp, U_fp, U_fp, U_fp, U_fp, U_fp, logical *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static integer negcon, lnames, lgsave, nmajor, inform__, minmax, mqnmod;
    static char cstart[8];
    static logical prtmem;
    static doublereal objtru;
    static logical fponly;
    static integer liwest, lprsav;
    static char solver[6];
    static integer errors, nextcw, useriw[130], nextiw, lrwest;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal userrw[130];
    static integer nextrw, stkyop;

/*     ================================================================== */
/*     npKerN does the work for npOpt. (Kernel for npOpt) */

/*     Developers can call this version with customized versions of */
/*     snLog, snLog2  and  snSTOP. */

/*     17 Oct 2004: First version of npKerN. */
/*     01 Sep 2007: Sticky parameters added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* LU factor tolerance. */
/* LU update tolerance. */
/* = 0(1) => cold(warm) start */
/* 0,1,2  => LM, FM, SD Hessian */
/* > 0    => print the solution */
/* 0(1) LU part(complete) piv */
/* Current Hessian type */
/* # of row and col. names */
/* Problem name */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --grad;
    --cmul;
    --istate;
    --bu;
    --bl;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    cjac_dim1 = *ldcj;
    cjac_offset = 1 + cjac_dim1;
    cjac -= cjac_offset;
    hess_dim1 = *ldh;
    hess_offset = 1 + hess_dim1;
    hess -= hess_offset;
    --ccon;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "NPOPT ", (ftnlen)6, (ftnlen)6);
    *info = 0;
/*     ------------------------------------------------------------------ */
/*     Check memory limits and fetch the workspace starting positions. */
/*     ------------------------------------------------------------------ */
    s2mem0_(info, solver, &c__500, leniw, lenrw, &iw[1], &mincw, &miniw, &
	    minrw, &maxcw, &maxiw, &maxrw, &nextcw, &nextiw, &nextrw, (ftnlen)
	    6);
    if (*info > 0) {
	goto L999;
    }
/*     Save the user's option choices  (weird choices get overwritten). */
/*     Initialize timers and the standard input file. */
/* Exit without printing */
    icopy_(&c__130, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__130, &rw[51], &c__1, userrw, &c__1);
    s1time_(&c__0, &c__0, &iw[1], leniw, &rw[1], lenrw);
    s1file_(&c__2, &iw[1], leniw);
    lnames = nextcw - 1;
/*     Check the arguments of npOpt. */
/* No names */
    s_copy(cstart, "Cold", (ftnlen)8, (ftnlen)4);
/* Preempted by lvlSrt */
    s3argn_(&inform__, cstart, lda, ldcj, ldh, n, nclin, ncnln, &iw[233], &bl[
	    1], &bu[1], cw + (lnames - 1 << 3), &istate[1], &cmul[1], &iw[69],
	     &errors, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
    if (inform__ > 0) {
	*info = inform__;
	goto L800;
    }
/*     Load the local problem dimensions. */
    ncon = *nclin + *ncnln;
    if (ncon == 0) {
/*        The problem is unconstrained. */
/*        Include a dummy row of zeros. */
	nncol = 0;
	m = 1;
	ne = 1;
    } else {
	nncol = *n;
	m = ncon;
	ne = m * *n;
    }
/*     Load the iw array with various problem dimensions. */
/*     First record problem dimensions for smart users to access in iw. */
    iobj = 0;
    nncon = *ncnln;
    nnjac = nncol;
    nnobj = *n;
/*     ------------------------------------------------------------------ */
/*     The obligatory call to npInit has already ``unset'' */
/*     the optional parameters.  However, it could not undefine */
/*     the char*8 options.  Do it now. */
/*     ------------------------------------------------------------------ */
    for (i__ = 51; i__ <= 180; ++i__) {
	s_copy(cw + (i__ - 1 << 3), "-1111111", (ftnlen)8, (ftnlen)8);
    }
/*     Set default options that relate specially to npOpt. */
/*     (Mainly, LU complete pivoting with small threshold pivot tol). */
    if (rw[66] < 0.) {
	rw[66] = 1.1;
    }
    if (rw[67] < 0.) {
	rw[67] = 1.1;
    }
    if (iw[156] < 0) {
	iw[156] = 2;
    }
    if (iw[72] < 0) {
	iw[72] = 1;
    }
    objadd = 0.;
/*     ------------------------------------------------------------------ */
/*     Load a generic problem name. */
/*     Check that the optional parameters have sensible values. */
/*     Delay printing the options until the arguments have been checked. */
/*     ------------------------------------------------------------------ */
    s_copy(cw + 400, "     NLP", (ftnlen)8, (ftnlen)8);
    s8dflt_(&m, n, &nncon, &nnjac, &nnobj, cw, &c__500, &iw[1], leniw, &rw[1],
	     lenrw, (ftnlen)8);
    s3prtb_(&m, n, &nncon, &nnjac, &nnobj, &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Determine storage requirements using the */
/*     following variables: */
/*         m,      n,     ne */
/*         lenR  , maxS , nnL */
/*         nnObj , nnCon, nnJac */
/*         negCon */
/*     All have to be known before calling s8Map. */
/*     ------------------------------------------------------------------ */
    nb = *n + m;
    nlocj = *n + 1;
    nkx = nb;
/* Computing MAX */
    i__1 = *ncnln * *n;
    negcon = max(i__1,1);
/*     Allocate arrays that are arguments of s8solv. */
/*     These are for the data, */
/*              locJ, indJ, Jcol, bl, bu, Names, */
/*     and for the solution */
/*              hs, x, pi, rc, hs. */
    lindj = nextiw;
    llocj = lindj + ne;
    lhs = llocj + nlocj;
    nextiw = lhs + nb;
    ljcol = nextrw;
    lbl = ljcol + ne;
    lbu = lbl + nb;
    lx = lbu + nb;
    lpi = lx + nb;
    lrc = lpi + m;
    nextrw = lrc + nb;
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    iw[20] = negcon;
    iw[28] = lenr;
/*     Load the iw array with various problem dimensions. */
    nnh = max(nnjac,nnobj);
    ngobj = nnobj;
/* Local nnObj is altered for FP */
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = m;
/* copy of the number of rows */
    iw[17] = ne;
/* copy of the number of nonzeros in Jcol */
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
    iw[233] = 1;
/*     ------------------------------------------------------------------ */
/*     If only a feasible point is requested, save the base point for the */
/*     objective function:  1/2 || x - x0 ||^2 */
/*     ------------------------------------------------------------------ */
    fponly = minmax == 0;
    if (fponly) {
	ngobj = nnh;
	lx0 = nextrw;
	lgsave = lx0 + nnh;
	nextrw = lgsave + nnobj;
	iw[298] = lx0;
/* ad hoc allocated FP base point */
	iw[339] = lgsave;
/* ad hoc copy of true obj grad for FP */
    }
/*     ------------------------------------------------------------------ */
/*     Allocate the local arrays for npOpt. */
/*     s8Map  maps snopt integer and double arrays. */
/*     s2BMap maps the arrays for the LU routines. */
/*     s2Mem  checks what space is available and prints any messages. */
/*     ------------------------------------------------------------------ */
    s8map_(&m, n, &negcon, &nkx, &nncon, &nnjac, &ngobj, &lenr, &maxr, &maxs, 
	    &mqnmod, &iw[72], &nextcw, &nextiw, &nextrw, &iw[1], leniw);
    s2bmap_(&m, n, &ne, &maxs, &nextiw, &nextrw, &maxiw, &maxrw, &liwest, &
	    lrwest, &iw[1], leniw);
    prtmem = TRUE_;
/* Print all messages in s2Mem */
    s2mem_(&inform__, &prtmem, &liwest, &lrwest, &nextcw, &nextiw, &nextrw, &
	    maxcw, &maxiw, &maxrw, &c__500, leniw, lenrw, &mincw, &miniw, &
	    minrw, &iw[1]);
    if (inform__ != 0) {
	*info = inform__;
	goto L800;
    }
    iw[256] = ljcol;
/* Jcol(neJ)   = Constraint Jacobian by columns */
    iw[257] = llocj;
/* locJ(n+1)   = column pointers for indJ */
    iw[258] = lindj;
/* indJ(ne) holds the row indices for Jij */
    iw[271] = lbl;
/* bl(nb)      = lower bounds */
    iw[272] = lbu;
/* bu(nb)      = upper bounds */
    iw[299] = lx;
/* x(nb)       = the solution (x,s) */
    iw[279] = lpi;
/* pi(m)       = the pi-vector */
    iw[280] = lrc;
/* rc(nb)      = the reduced costs */
    iw[282] = lhs;
/* the column state vector */
    iw[359] = lnames;
/* Names(nName) */
    lgobj = iw[297];
/* gObj(nnObj) = Objective gradient */
    lfcon = iw[316];
/* fCon (nnCon)  constraints at x */
    lgcon = iw[320];
/*     Define the row and column ordering for J. */
/*     NPOPT  uses natural order throughout, so kx = kxN. */
/* gCon (negCon) constraint gradients at x */
    iw[247] = nkx;
/* dimension of kx and its inverse, kxN */
    lkx = iw[251];
/* j  = kx (jN) => col j of Jcol is variable jN */
    iw[252] = lkx;
/* jN = kxN(j ) => col j of Jcol is variable jN */
    s1perm_(n, &iw[lkx]);
    s1perm_(&m, &iw[lkx + *n]);
    if (iw[69] == 0) {
	start = 0;
    } else {
	start = 2;
    }
/*     ------------------------------------------------------------------ */
/*     Initialize some SNOPT arrays that are copied to NPOPT arrays. */
/*     Build the Jacobian and load the SNOPT arrays. */
/*     ------------------------------------------------------------------ */
    nncon0 = max(nncon,1);
    s3inin_(&start, n, &nb, &nncon0, &nncon, &negcon, &iw[lhs], &rw[lfcon], &
	    rw[lgcon], &rw[lgobj], &rw[lrc], &rw[lx], &iw[1], leniw, &rw[1], 
	    lenrw);
    s3inn_(&start, lda, ldh, &m, n, nclin, &ncon, &nncol, &nb, &nncon0, &
	    nncon, &iw[lhs], &istate[1], &a[a_offset], &ne, &nlocj, &iw[llocj]
	    , &iw[lindj], &rw[ljcol], &bl[1], &bu[1], &rw[lbl], &rw[lbu], &
	    ccon[1], &cmul[1], &hess[hess_offset], &rw[lpi], &x[1], &rw[lx], &
	    iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Construct column pointers for the nonlinear part of the  Jacobian. */
/*     ------------------------------------------------------------------ */
    if (nncon > 0) {
	llocg = iw[260];
/* locG(nlocG) = column pointers for indG */
	nlocg = nnjac + 1;
	s8gloc_(&nncon, &nnjac, &ne, &nlocj, &iw[llocj], &iw[lindj], &negcon, 
		&nlocg, &iw[llocg]);
    }
    if (fponly) {
	dcopy_(&nnh, &rw[lx], &c__1, &rw[lx0], &c__1);
    }
/*     ------------------------------------------------------------------ */
/*     Solve the problem. */
/*     Tell s8solv that we don't have an initial Hessian. */
/*     ------------------------------------------------------------------ */
    iw[202] = -1;
    lprsav = iw[84];
    iw[84] = 0;
    s8solv_(info, solver, &c__0, (U_fp)s0fgn_, (U_fp)funcon, (U_fp)funobj, (
	    U_fp)snlog, (U_fp)snlog2, (U_fp)snstop, &gotr, &m, n, &nb, &nncon,
	     &nnjac, &ngobj, &iw[233], &iobj, &objadd, &fobj, &objtru, &ninf, 
	    &sinf, &ne, &nlocj, &iw[llocj], &iw[lindj], &rw[ljcol], &rw[lbl], 
	    &rw[lbu], cw + (lnames - 1 << 3), &iw[lhs], &rw[lx], &rw[lpi], &
	    rw[lrc], &nmajor, &ns, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, 
	    cw, &c__500, &iw[1], leniw, &rw[1], lenrw, (ftnlen)6, (ftnlen)8, (
	    ftnlen)8, (ftnlen)8);
    iw[84] = lprsav;
    *objf = fobj;
    *majits = nmajor;
/*     ------------------------------------------------------------------ */
/*     Unload the SNOPT arrays. */
/*     ------------------------------------------------------------------ */
    if (fponly && nnobj > 0) {
	dcopy_(&nnobj, &rw[lgsave], &c__1, &rw[lgobj], &c__1);
    }
    s3outn_(ldcj, ldh, n, nclin, &ncon, &nb, &nncon0, &nncon, &iw[lhs], &
	    istate[1], &ccon[1], &cjac[cjac_offset], &cmul[1], &rw[lfcon], &
	    rw[lgcon], &rw[lgobj], &grad[1], &hess[hess_offset], &rw[lrc], &x[
	    1], &rw[lx], &iw[1], leniw, &rw[1], lenrw);
    xnorm = dnrm1s_(n, &rw[lx], &c__1);
    i__1 = *n + ncon;
    s3prtn_(n, &i__1, nclin, &nncon0, lda, &iw[84], &xnorm, &istate[1], &a[
	    a_offset], &bl[1], &bu[1], &ccon[1], &cmul[1], &x[1], &rw[lx], &
	    iw[1], leniw, &rw[1], lenrw);
/*     If "sticky parameters no",  restore the user-defined options */
    stkyop = iw[116];
    if (stkyop <= 0) {
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
} /* npkern_ */

