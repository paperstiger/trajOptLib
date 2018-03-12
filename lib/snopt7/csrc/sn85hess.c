/* ./src/sn85hess.f -- translated by f2c (version 20100827).
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
static doublereal c_b11 = 1.;
static doublereal c_b16 = -1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn85Hess.f */

/*     s8getH   s8H0     s8Hfix   s8HQN    s8Hwrp   s8Hx   s8xHx */
/*     s8Hupd   s8x1 */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8geth_(integer *nnh, integer *lenh, doublereal *u, 
	doublereal *h__, doublereal *y, doublereal *y1)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), s6rprd_(integer *, integer *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *);
    static integer jthcol;

/*     ================================================================== */
/*     s8getH  computes the product H = U'U, where  U is the Cholesky */
/*     factor of the approximate Hessian of the Lagrangian.  The matrix */
/*     U is stored by rows in the one-dimensional array  U. */
/*     lenH defines the length of U.  lenH must be at least */
/*     nnH*(nnH + 1)/2.  The result is stored by columns in the upper */
/*     triangular array H. */

/*     03 Sep 2006: First version of s8getH. */
/*     03 Sep 2006: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Compute the product y1 = U'Uy, where  U is an upper- */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y1;
    --y;
    --h__;
    --u;

    /* Function Body */
    jthcol = 1;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	jthcol = jthcol + j - 1;
	dload_(nnh, &c_b2, &y[1], &c__1);
	y[j] = 1.;
	s6rprd_(&c__0, nnh, nnh, lenh, &u[1], &y[1], &y1[1]);
	s6rprd_(&c__1, nnh, nnh, lenh, &u[1], &y1[1], &y[1]);
	dcopy_(&j, &y[1], &c__1, &h__[jthcol], &c__1);
    }
    return 0;
} /* s8geth_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8getH */
/* Subroutine */ int s8h0_(integer *hqntype, integer *nnh, doublereal *u0pre, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    static integer lu, lu0, lhd, lenu;
    extern /* Subroutine */ int s9fmh0_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    s9lmh0_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal hcndbd;
    static integer lvlhes;

/*     ================================================================== */
/*     s8H0    initializes the BFGS approximate Hessian. */

/*     s8H0   calls one of the Hessian routines s9LMH0, s9FMH0, s9SDH0 */
/*     according to the value of the option lvlHes. */
/*     Each of these routines defines a particular form of the Hessian. */
/*     At the moment the options are: */
/*        lvlHes = 0   Limited-Memory (LM) BFGS Hessian  (the default). */
/*        lvlHes = 1   Full-Memory    (FM) BFGS Hessian. */

/*     On entry, the value of HQNType is as follows: */

/*      HQNType */
/*      ------- */
/*       -1  (= HUnset)  H is undefined. */
/*        0  (= HNorml)  H is of the form defined by lvlHes. */
/*        1  (= HDiag )  H is a diagonal matrix. */
/*        2  (= HUnit )  H is an identity matrix. */

/*     19 Jul 1995: First version of s8H0. */
/*     12 Jan 1996: Full memory Hessian option added. */
/*     18 Feb 2001: H stored in product form. */
/*     23 Jun 2008: Exact optiona added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* # of updates since last reset */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    lvlhes = iw[72];
/* LM, FM or Exact Hessian */
    if (lvlhes == 0) {
/*        ----------------------- */
/*        Limited-memory Hessian. */
/*        ----------------------- */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	lhd = iw[347];
/* Diagonal of BFGS Hessian */
	s9lmh0_(hqntype, nnh, &hcndbd, u0pre, &rw[lhd], &rw[lu0], &iw[381]);
    } else if (lvlhes == 1) {
/*        ----------------------- */
/*        Full-memory Hessian. */
/*        ----------------------- */
	lhd = iw[347];
/* Diagonal of BFGS Hessian */
	lu = iw[391];
/* U(lenU), full-memory BFGS H = U'U */
	lenu = iw[392];

	s9fmh0_(hqntype, nnh, &hcndbd, u0pre, &rw[lhd], &lenu, &rw[lu], &iw[
		381]);
    }
    return 0;
} /* s8h0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8H0 */
/* Subroutine */ int s8hfix_(integer *nncon, integer *nnjac, doublereal *tolz,
	 integer *ne, integer *nlocj, integer *locj, integer *indj, integer *
	negcon, integer *nlocg, integer *locg, doublereal *ydx, doublereal *
	ydxmin, doublereal *penunm, doublereal *fcon, doublereal *fcon1, 
	doublereal *gcon, doublereal *gcon1, doublereal *dx, doublereal *gd, 
	doublereal *penu, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal wi, diff, beta;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static doublereal peni, wmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal wnorm;
    extern /* Subroutine */ int s8gprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static logical gotpen, overfl;

/*     ================================================================== */
/*     s8Hfix  attempts to find the a vector xPen  of minimum two-norm */
/*     such that there exists a BFGS update for the modified Lagrangian */
/*       La   = f(x) - lambda'(fCon1 - LfCon) */
/*                   + 1/2 (fCon1 - LfCon)'*diag(PenU)*(fCon1 - LfCon), */

/*     where  LfCon = fCon + J(x1)*dx. */

/*     On entry, */
/*     dx     is the nonlinear part of the search direction x2 - x1. */
/*     gd     is the Lagrangian gradient difference. */
/*     gCon    is the Jacobian at the old x. */
/*     gCon1    is the Jacobian at the new x. */
/*     ydx    is the approximate curvature of the Lagrangian. */
/*     ydxmin   (ydx < ydxmin) is the smallest acceptable approximate */
/*              curvature. */

/*     On exit, */
/*     gd     is the augmented Lagrangian gradient difference. */
/*     PenU     are the penalty parameters. */
/*     ydx    is unchanged unless gotPen is true, in which case */
/*              ydx = ydxmin. */

/*     08 Dec 1991: First version based on  Npsol  routine npupdt. */
/*     26 Oct 2000: Current version of s8Hfix. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* -->  parameter         (PenMax = 1.0d+5) */
/* -->  parameter         (PenMax = 1.0d+16) */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --w;
    --v;
    --penu;
    --fcon1;
    --fcon;
    --gd;
    --dx;
    --indj;
    --locj;
    --gcon1;
    --gcon;
    --locg;

    /* Function Body */
    overfl = FALSE_;
/*     Try an augmented Lagrangian term to increase ydx. */
    *penunm = 0.;
/*     Compute  v = J1*dx and w = (J2 - J1)*dx = J2*dx - v. */
    s8gprd_(&c__0, tolz, ne, nlocj, &locj[1], &indj[1], negcon, nlocg, &locg[
	    1], &gcon[1], &c_b11, &dx[1], nnjac, &c_b2, &v[1], nncon);
    s8gprd_(&c__0, tolz, ne, nlocj, &locj[1], &indj[1], negcon, nlocg, &locg[
	    1], &gcon1[1], &c_b11, &dx[1], nnjac, &c_b2, &w[1], nncon);
    daxpy_(nncon, &c_b16, &v[1], &c__1, &w[1], &c__1);
/*     Compute the difference between c and its linearization. */
/*     v  =  c - cL = fCon1 - (fCon + J1*s) = fCon1 - fCon - J1*s. */
    daxpy_(nncon, &c_b16, &fcon1[1], &c__1, &v[1], &c__1);
    daxpy_(nncon, &c_b11, &fcon[1], &c__1, &v[1], &c__1);
    dscal_(nncon, &c_b16, &v[1], &c__1);
/*     --------------------------------------------------------- */
/*     Compute the minimum-length vector of penalty parameters */
/*     that makes the approximate curvature equal to  ydxmin. */
/*     --------------------------------------------------------- */
/*     Use w to hold the constraint on PenU. */
/*     Minimize            norm(PenU) */
/*     subject to   ( Sum( w(i)*PenU(i) )  =   const, */
/*                  (           PenU(i)   .ge. 0. */
    wmax = 0.;
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__] * v[i__];
	wmax = max(wmax,wi);
	w[i__] = max(0.,wi);
    }
    wnorm = dnrm2_(nncon, &w[1], &c__1);
    diff = *ydxmin - *ydx;
    d__1 = wmax * diff;
/* Computing 2nd power */
    d__3 = wnorm;
    d__2 = d__3 * d__3;
    beta = ddiv_(&d__1, &d__2, &overfl);
    gotpen = ! overfl && wmax > 0. && beta < 1e5;
    if (gotpen) {
/* Computing 2nd power */
	d__1 = wnorm;
	beta = diff / (d__1 * d__1);
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wi = w[i__];
	    peni = beta * wi;
	    v[i__] = peni * v[i__];
	    *ydx += peni * wi;
	    penu[i__] = peni;
	}
	*ydx = max(*ydx,*ydxmin);
	*penunm = dnrm2_(nncon, &penu[1], &c__1);
/*        Update  gd  by the term  (J2' - J1')*v, */
/*        with v = diag(PenU)*(fCon1 - fCon - J1*s) from above. */
	s8gprd_(&c__1, tolz, ne, nlocj, &locj[1], &indj[1], negcon, nlocg, &
		locg[1], &gcon1[1], &c_b11, &v[1], nncon, &c_b11, &gd[1], 
		nnjac);
	s8gprd_(&c__1, tolz, ne, nlocj, &locj[1], &indj[1], negcon, nlocg, &
		locg[1], &gcon[1], &c_b16, &v[1], nncon, &c_b11, &gd[1], 
		nnjac);
    }
/* gotPen */
    return 0;
} /* s8hfix_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Hfix */
/* Subroutine */ int s8hqn_(integer *iexit, U_fp fgwrap, U_fp fgcon, U_fp 
	fgobj, logical *usefd, integer *hqntype, integer *qptype, integer *
	info, integer *lenr, integer *m, integer *mbs, integer *n, integer *
	nb, integer *nncon0, integer *nncon, integer *nnjac, integer *nnl, 
	integer *nnobj0, integer *nnobj, integer *ns, integer *nmajor, 
	integer *nskip, doublereal *u0ii, doublereal *step, integer *minimz, 
	doublereal *dxhdx, integer *rtrmod, logical *gotr, logical *incrun, 
	doublereal *pendmp, doublereal *penmax, doublereal *fobj, doublereal *
	fcon, doublereal *gcon, doublereal *gobj, doublereal *fcon1, 
	doublereal *gcon1, doublereal *gobj1, integer *ne, integer *nlocj, 
	integer *locj, integer *indj, doublereal *jcol, integer *negcon, 
	integer *nlocg, integer *locg, integer *kbs, doublereal *bl, 
	doublereal *bu, doublereal *dx, doublereal *dg, doublereal *udx, 
	doublereal *hdx, doublereal *ycon1, doublereal *r__, doublereal *x, 
	doublereal *x1, doublereal *xqp0, doublereal *xpen, doublereal *y, 
	doublereal *y1, doublereal *y2, char *cu, integer *lencu, integer *iu,
	 integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer nbs;
    static doublereal eps;
    extern /* Subroutine */ int s8h0_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), s8x1_(integer *, 
	    U_fp, U_fp, U_fp, logical *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static doublereal ydx, eps0, eps1;
    static integer kfac;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), ddiv_(doublereal *, doublereal *, logical *);
    static integer maxr;
    static doublereal rnnl, rydx;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal xpen0;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer mskip;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer nzero;
    static doublereal snorm, u0scal;
    extern /* Subroutine */ int s8gprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), s8hfix_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s8hupd_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s6rupd_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static logical updatd;
    static doublereal sgnobj;
    static logical nlncon;
    static integer inform__;
    static logical overfl;
    static doublereal glnorm, rdxhdx, penunm, ydxmin;
    static integer lvlsrt;

/*     ================================================================== */
/*     s8HQN  does the quasi-Newton update with vectors */
/*        dx = x1 - x   and   dg = gL(x1) - gL(x). */

/*     On entry: */
/*      xQP is the QP solution. */

/*     23 Apr 1999: First version of s8HQN, */
/*     18 Feb 2001: LM H stored in product form. */
/*     12 Oct 2003: snEXIT and SNPRNT adopted */
/*     10 Jan 2005: FM H stored in product form. */
/*     23 Jun 2008: y3 no longer an argument. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --info;
    --r__;
    --kbs;
    --x;
    --y2;
    --y1;
    --y;
    --xqp0;
    --bu;
    --bl;
    --xpen;
    --ycon1;
    --fcon1;
    --fcon;
    --x1;
    --hdx;
    --udx;
    --dg;
    --dx;
    --gobj1;
    --gobj;
    --jcol;
    --indj;
    --locj;
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
    maxr = iw[52];
/* max columns of R */
    kfac = iw[59];
/* factorization frequency */
    mskip = iw[67];
/* # largest allowable  nSkip */
    lvlsrt = iw[69];
/* = 0:1:2:3 => cold:warm:basis:hot start */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5) */
    eps1 = rw[3];
/* eps**(2/3) */
    xpen0 = rw[89];
/* initial penalty parameter. */
    *iexit = 0;
    overfl = FALSE_;
    nlncon = *nncon > 0;
    nbs = *m + *ns;
    sgnobj = (doublereal) (*minimz);
    info[1] = 0;
    info[2] = 0;
    ydx = 0.;
/*     --------------------------------------------------------------- */
/*     Compute  dx = x1 - x  and  dg = gL1 - gL. */
/*     Compute the approx. curvature ydx and new scale factor U0. */
/*     --------------------------------------------------------------- */
    dcopy_(nnl, &x1[1], &c__1, &dx[1], &c__1);
    daxpy_(nnl, &c_b16, &x[1], &c__1, &dx[1], &c__1);
    dscal_(nnl, step, &hdx[1], &c__1);
    dscal_(nnl, step, &udx[1], &c__1);
    *dxhdx = *dxhdx * *step * *step;
    if (*nnobj > 0) {
	dcopy_(nnobj, &gobj1[1], &c__1, &dg[1], &c__1);
	if (*minimz < 0) {
	    dscal_(nnobj, &sgnobj, &dg[1], &c__1);
	}
    }
    nzero = *nnl - *nnobj;
    if (nzero > 0) {
	dload_(&nzero, &c_b2, &dg[*nnobj + 1], &c__1);
    }
    if (*nncon > 0) {
	s8gprd_(&c__1, &eps0, ne, nlocj, &locj[1], &indj[1], negcon, nlocg, &
		locg[1], &gcon1[1], &c_b16, &ycon1[1], nncon, &c_b11, &dg[1], 
		nnjac);
    }
/*     gLnorm = dnormi( nnL, dg, 1 ) */
    glnorm = dnrm2_(nnl, &dg[1], &c__1);
    if (*nnobj > 0) {
	d__1 = -sgnobj;
	daxpy_(nnobj, &d__1, &gobj[1], &c__1, &dg[1], &c__1);
    }
    if (*nncon > 0) {
	s8gprd_(&c__1, &eps0, ne, nlocj, &locj[1], &indj[1], negcon, nlocg, &
		locg[1], &gcon[1], &c_b11, &ycon1[1], nncon, &c_b11, &dg[1], 
		nnjac);
    }
    ydx = ddot_(nnl, &dg[1], &c__1, &dx[1], &c__1);
    if (*nmajor == 1 && lvlsrt != 3) {
/*        =============================================================== */
/*        First iteration. Do not attempt a BFGS update. */
/*        Use latest curvature information to get a better scaled H. */
/*        =============================================================== */
	if (glnorm > 0.) {
	    rnnl = (doublereal) (*nnl);
	    *u0ii = sqrt(glnorm / sqrt(rnnl));
	} else {
	    *u0ii = 1.;
	}
	*gotr = FALSE_;
/* Computing MIN */
	d__1 = max(*u0ii,.01);
	*u0ii = min(d__1,10.);
	s8h0_(hqntype, nnl, u0ii, &iw[1], leniw, &rw[1], lenrw);
    } else {
/*        =============================================================== */
/*        Except on the first iteration, attempt a BFGS update. */
/*        Compute the smallest allowable curvature. */
/*        If the update cannot be done, s8x1 attempts to find a */
/*        modified update using  dx = x1 - x defined with a new x. */
/*        Arrays fCon, gCon and gObj must be redefined at the new x. */
/*        =============================================================== */
	snorm = dnrm2_(nnl, &dx[1], &c__1);
	d__1 = sqrt((abs(ydx)));
	*u0ii = ddiv_(&d__1, &snorm, &overfl);
/* Computing MIN */
	d__1 = max(*u0ii,.01);
	*u0ii = min(d__1,10.);
	penunm = 0.;
	ydxmin = *dxhdx * .001;
	updatd = *dxhdx > 0. && (ydx >= ydxmin || ydx >= eps1);
	if (nlncon && ! updatd) {
/*           ------------------------------------------------------------ */
/*           Redefine  x, dx, Hdx and dg. */
/*           The problem functions are recomputed at x. */
/*           ------------------------------------------------------------ */
	    s8x1_(&inform__, (U_fp)fgwrap, (U_fp)fgcon, (U_fp)fgobj, usefd, n,
		     nb, nncon0, nncon, nnjac, nnobj0, nnobj, nnl, minimz, 
		    step, dxhdx, &ydx, fobj, &fcon[1], &gcon[1], &gobj[1], &
		    gcon1[1], &gobj1[1], ne, nlocj, &locj[1], &indj[1], 
		    negcon, nlocg, &locg[1], &bl[1], &bu[1], &dx[1], &dg[1], &
		    udx[1], &hdx[1], &ycon1[1], &y1[1], &x[1], &xqp0[1], &y[1]
		    , cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	    if (inform__ > 0) {
		*iexit = inform__;
/* User wants to stop */
		goto L999;
	    }
	    ydxmin = *dxhdx * .001;
	    updatd = *dxhdx > 0. && (ydx >= ydxmin || ydx >= eps1);
	    if (updatd) {
		info[2] = 1;
	    }
	    if (! updatd && *dxhdx > 0.) {
/*              --------------------------------------------------------- */
/*              If all else fails, attempt to update the Hessian of */
/*              the augmented Lagrangian. */
/*              The target ydx is defined via tolg2. */
/*              --------------------------------------------------------- */
		ydxmin = *dxhdx * .1;
		s8hfix_(nncon, nnjac, &eps0, ne, nlocj, &locj[1], &indj[1], 
			negcon, nlocg, &locg[1], &ydx, &ydxmin, &penunm, &
			fcon[1], &fcon1[1], &gcon[1], &gcon1[1], &dx[1], &dg[
			1], &y[1], &y1[1], &y2[1]);
		updatd = ydx >= ydxmin;
		if (updatd) {
		    info[2] = 2;
		}
	    }
	}
/* nlnCon */
	if (updatd) {
/*           ------------------------------------------------------------ */
/*           Update the approximate Hessian using (dg,Hdx). */
/*           If there are no nonlinear constraints,  apply the update */
/*           to the reduced Hessian. */
/*           ------------------------------------------------------------ */
	    *nskip = 0;
	    if (ydx >= ydxmin && *hqntype == 0) {
		info[1] = 1;
/* conventional BFGS */
	    } else {
		info[1] = 2;
/* self-scaled  BFGS */
	    }
	    rydx = sqrt(ydx);
	    rdxhdx = sqrt(*dxhdx);
	    s8hupd_(&info[1], hqntype, nnl, u0ii, &u0scal, &rydx, &rdxhdx, &
		    dx[1], &hdx[1], &dg[1], &iw[1], leniw, &rw[1], lenrw);
	    if (*qptype == 0) {
		*gotr = *gotr && *nncon == 0 && *ns > 0 && *hqntype == 0 && *
			rtrmod < kfac;
	    } else {
		*gotr = *gotr && *ns > 0 && *hqntype == 0;
	    }
	    if (*gotr) {
		s6rupd_(&info[1], &maxr, lenr, m, n, &nbs, nnl, ns, &u0scal, &
			rdxhdx, ne, nlocj, &locj[1], &indj[1], &jcol[1], &kbs[
			1], &dg[1], &hdx[1], &r__[1], &y[1], &y1[1], &y2[1], &
			iw[1], leniw, &rw[1], lenrw);
		++(*rtrmod);
	    } else {
		*rtrmod = 0;
	    }
	} else {
/*           ------------------------------------------------------------ */
/*           No suitable update pair (dg,Hdx) could be found. */
/*           Skip the update.  Too many skips and we reset. */
/*           ------------------------------------------------------------ */
	    ++(*nskip);
/*           Apply all updates to H and discard the off-diagonals. */
	    if (*nskip % mskip == 0) {
		s8h0_(hqntype, nnl, u0ii, &iw[1], leniw, &rw[1], lenrw);
		if (*nskip % (mskip << 1) == 0) {
/*                 ------------------------------------------------------ */
/*                 Reset the multipliers and penalty parameters */
/*                 ------------------------------------------------------ */
		    *incrun = TRUE_;
		    *pendmp = 1.;
		    *penmax = 1. / eps;
		    dload_(nncon, &xpen0, &xpen[1], &c__1);
		    dload_(nncon, &c_b2, &ycon1[1], &c__1);
		}
	    }
	}
    }
/* nMajor > 1 */
L999:
    return 0;
} /* s8hqn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8HQN */
/* Subroutine */ int s8hwrp_(S_fp hprod, integer *nnh, doublereal *x, 
	doublereal *hx, integer *status, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
/*     ================================================================== */
/*     s8Hwrp wraps Hprod, which multiplies the QP Hessian H by the */
/*     vector  x.   It is called by the QP solver. */

/*     On entry: */
/*        Status  = 0  => a normal call for H*x. */
/*        Status  = 1  => the first entry for a given QP. */
/*        Status ge 2  => last call for a given QP. Status = 2+iExit. */

/*     On exit: */
/*        Status lt 0   the user wants to stop. */

/*     03 Nov 2000: First version of s8Hwrp. */
/*     03 Nov 2000: Current version. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hx;
    --x;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    (*hprod)(nnh, &x[1], &hx[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, 
	    (ftnlen)8);
    return 0;
} /* s8hwrp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Hwrp */
/* Subroutine */ int s8hx_(integer *nnh, doublereal *x, doublereal *hx, char *
	cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw, ftnlen cw_len)
{
    static integer ls, lu, lv, lu0, lux, lenu;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s9fmhx_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *), s9lmhx_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal sgnobj;
    static integer mqnmod, lvlhes, minimz;

/*     ================================================================== */
/*     s8Hx  multiplies the QP Hessian  H by the vector  x. */
/*     It is used to define Hx for the QP subproblem. */

/*     This routine is called by a general QP solver, which will rescale */
/*     Hx by  sgnObj  when maximizing. */

/*     s8Hx calls one of the Hessian routines s9LMH, s9FMH, s9SDH, ... */
/*     according to the value of the options lvlDer and lvlHes. */
/*     Each of these routines defines a particular form of the Hessian. */
/*     At the moment the options are: */

/*        lvlHes = LM      Limited-Memory (LM) BFGS  (the default). */
/*        lvlHes = FM      Full-Memory    (FM) BFGS */
/*        lvlHes = Exact   FD or exact Hessian */

/*     30 Dec 1991: First version of s8Hx. */
/*     12 Jan 1996: Full memory Hessian option added. */
/*     04 Apr 1999: Exact and FD Hessian option added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     10 Jan 2005: FM H stored in product form. */
/*     23 Jun 2008: Exact option added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --x;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lvlhes = iw[72];
/* LM, FM or Exact Hessian */
    minimz = iw[199];
/* (-1)(+1)    => (max)(min) */
    lux = iw[345];
/* Ux(nnL)     = product of U with x */
    sgnobj = (doublereal) minimz;
    if (lvlhes == 0) {
/*        ----------------------- */
/*        Limited-memory Hessian. */
/*        ----------------------- */
	mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	ls = iw[401];
/* sk's for BFGS products: (I + sk*vk') */
	lv = iw[402];
/* vk's for BFGS products: (I + sk*vk') */
	s9lmhx_(nnh, &x[1], &rw[lux], &hx[1], &mqnmod, &iw[381], &rw[lu0], &
		rw[ls], &rw[lv]);
    } else if (lvlhes == 1) {
/*        ----------------------- */
/*        Full-memory Hessian. */
/*        ----------------------- */
	lu = iw[391];

	lenu = iw[392];

	s9fmhx_(nnh, &x[1], &rw[lux], &hx[1], &lenu, &rw[lu]);
    }
    if (minimz < 0) {
	dscal_(nnh, &sgnobj, &hx[1], &c__1);
    }
    return 0;
} /* s8hx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Hx */
/* Subroutine */ int s8xhx_(integer *nnh, doublereal *x, doublereal *ux, 
	doublereal *hx, doublereal *xhx, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    static integer ls, lu, lv, lu0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer lenu;
    extern /* Subroutine */ int s9fmhx_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *), s9lmhx_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer mqnmod, lvlhes;

/*     ================================================================== */
/*     s8xHx  computes x'Hx and Hx, where H = U'U. */

/*     s8xHx calls one of the Hessian routines s9LMH, s9FMH, s9SDH, ... */
/*     according to the value of the options lvlDer and lvlHes. */
/*     Each of these routines defines a particular form of the Hessian. */
/*     At the moment the options are: */

/*        lvlHes = LM      Limited-Memory (LM) BFGS  (the default). */
/*        lvlHes = FM      Full-Memory    (FM) BFGS */
/*        lvlHes = Exact   FD or exact Hessian */

/*     10 Jan 2005: First version of s8Hx based on s8Hx */
/*     18 Feb 2001: LM H stored in product form. */
/*     10 Jan 2005: FM H stored in product form. */
/*     23 Jun 2008: Exact option added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --ux;
    --x;
    --iw;
    --rw;

    /* Function Body */
    lvlhes = iw[72];
/* LM, FM or Exact Hessian */
    if (lvlhes == 0) {
/*        ----------------------- */
/*        Limited memory Hessian. */
/*        ----------------------- */
	mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	ls = iw[401];
/* sk's for BFGS products: (I + sk*vk') */
	lv = iw[402];
/* vk's for BFGS products: (I + sk*vk') */
	s9lmhx_(nnh, &x[1], &ux[1], &hx[1], &mqnmod, &iw[381], &rw[lu0], &rw[
		ls], &rw[lv]);
    } else if (lvlhes == 1) {
/*        ----------------------- */
/*        Full memory Hessian. */
/*        ----------------------- */
	lu = iw[391];
/* U(lenU), dense Hessian factor */
	lenu = iw[392];

	s9fmhx_(nnh, &x[1], &ux[1], &hx[1], &lenu, &rw[lu]);
    }
    *xhx = ddot_(nnh, &ux[1], &c__1, &ux[1], &c__1);
    return 0;
} /* s8xhx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8xHx */
/* Subroutine */ int s8hupd_(integer *update, integer *hqntype, integer *nnh, 
	doublereal *u0pre, doublereal *u0scal, doublereal *rydx, doublereal *
	rdxhdx, doublereal *dx, doublereal *hdx, doublereal *y, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    static integer ls, lu, lv, lu0, lhd, lenu, ludx;
    static doublereal hcndbd;
    extern /* Subroutine */ int s9fmup_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *), s9lmup_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer mqnmod, lvlhes;

/*     ================================================================== */
/*     s8Hupd  applies the pair of vectors that define the BFGS update */
/*     or self-scaled BFGS update. */

/*     On entry: */

/*     Hdx   contains  H  times the difference x1 - x. */
/*     y     contains the gradient  difference g1 - g. */

/*     On exit: */
/*     Hdx   contains  H  times the difference x1 - x(new). */
/*     y     contains the gradient  difference g1 - g(new). */

/*     s8Hupd calls one of the Hessian routines s9LMH, s9FMH, s9SDH, ... */
/*     according to the value of the option lvlHes. */
/*     At the moment the options are: */

/*        lvlHes = LM      Limited-Memory (LM) BFGS  (the default). */
/*        lvlHes = FM      Full-Memory    (FM) BFGS */
/*        lvlHes = Exact   FD or exact Hessian */

/*     19 Jul 1995: First version of s8Hupd. */
/*     12 Jan 1996: Full-memory Hessian option added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     12 Jan 2005: FM H stored in product form. */
/*     16 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* # of updates since last reset */
    /* Parameter adjustments */
    --y;
    --hdx;
    --dx;
    --iw;
    --rw;

    /* Function Body */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    lvlhes = iw[72];
/* LM, FM or Exact Hessian */
    lhd = iw[347];
/* Diagonal of BFGS Hessian */
    if (lvlhes == 0) {
	mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
	lu0 = iw[346];
/* Square root of initial BFGS diagonal */
	ls = iw[401];
/* sk's for BFGS products: (I + sk*vk') */
	lv = iw[402];
/* vk's for BFGS products: (I + sk*vk') */
	s9lmup_(update, hqntype, nnh, &mqnmod, &iw[381], &hcndbd, u0pre, 
		u0scal, rydx, rdxhdx, &rw[lhd], &hdx[1], &y[1], &dx[1], &rw[
		lu0], &rw[ls], &rw[lv]);
    } else if (lvlhes == 1) {
	mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
	ludx = iw[345];
/* Ux(nnL)      = product of U with x */
	lu = iw[391];
/* U(lenU), full-memory BFGS Hessian H = U'U */
	lenu = iw[392];

	s9fmup_(update, hqntype, nnh, &mqnmod, &iw[381], &hcndbd, u0pre, 
		u0scal, rydx, rdxhdx, &rw[lhd], &hdx[1], &y[1], &rw[ludx], &
		lenu, &rw[lu]);
    }
    return 0;
} /* s8hupd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Hupd */
/* Subroutine */ int s8x1_(integer *iexit, S_fp fgwrap, U_fp fgcon, U_fp 
	fgobj, logical *usefd, integer *n, integer *nb, integer *nncon0, 
	integer *nncon, integer *nnjac, integer *nnobj0, integer *nnobj, 
	integer *nnl, integer *minimz, doublereal *step, doublereal *dxhdx, 
	doublereal *ydx, doublereal *fobj, doublereal *fcon, doublereal *gcon,
	 doublereal *gobj, doublereal *gcon1, doublereal *gobj1, integer *ne, 
	integer *nlocj, integer *locj, integer *indj, integer *negcon, 
	integer *nlocg, integer *locg, doublereal *bl, doublereal *bu, 
	doublereal *dx, doublereal *dg, doublereal *udx, doublereal *hdx, 
	doublereal *ycon1, doublereal *tdx, doublereal *x, doublereal *xqp0, 
	doublereal *y, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal eps, eps0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int s6fdg_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    S_fp, U_fp, U_fp, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen), s8xhx_(integer *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *), 
	    dload_(integer *, doublereal *, doublereal *, integer *), dscal_(
	    integer *, doublereal *, doublereal *, integer *), dcopy_(integer 
	    *, doublereal *, integer *, doublereal *, integer *), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nzero;
    extern /* Subroutine */ int s8gprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static integer modefg;
    static logical nlnobj;
    static doublereal sgnobj;
    static logical nlncon;

/*     ================================================================== */
/*     s8x1   redefines the quantities x and dx, Hdx and dg  used for the */
/*     quasi-Newton update.  The problem functions are recomputed at x. */

/*     The new  x1  is  x1 + step*(xQP0 - x1),  where xQP0 is a */
/*     (nonelastic) feasible point from the QP subproblem. */

/*     s8x1 is always called with nnL > 0. */

/*     02 Dec 1994: First version of s8x1. */
/*     20 Jul 1998: s8x1 made self-contained */
/*     24 Aug 1998: Fixed bug found by Alan Brown at Nag. */
/*                  FD derivatives now computed correctly. */
/*                  Parameter useFD added. */
/*     11 Oct 1998: Facility to combine funobj and funcon added. */
/*     12 Oct 2003: snEXIT and snPRNT adopted. */
/*     14 Jan 2005: Argument Udx added for call to s8xHx. */
/*     16 Jun 2008: Call-status implemented correctly */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --y;
    --ycon1;
    --fcon;
    --gobj1;
    --gobj;
    --xqp0;
    --tdx;
    --hdx;
    --udx;
    --dg;
    --dx;
    --indj;
    --locj;
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
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5) */
    *iexit = 0;
    modefg = 2;
    sgnobj = (doublereal) (*minimz);
    nlncon = *nncon > 0;
    nlnobj = *nnobj > 0;
/*     Save dx in case a better update pair cannot be found. */
    dcopy_(nnl, &dx[1], &c__1, &tdx[1], &c__1);
/*     ------------------------------------------------- */
/*     dx =  dx - step*y,  with  y = xQP0 - x */
/*     ------------------------------------------------- */
    daxpy_(nnl, &c_b16, &x[1], &c__1, &xqp0[1], &c__1);
    d__1 = -(*step);
    daxpy_(nnl, &d__1, &xqp0[1], &c__1, &dx[1], &c__1);
/*     ------------------------------------------------- */
/*     Compute the minimum curvature. */
/*     If nnL < n, dxHdx may be zero (or negative fuzz). */
/*     ------------------------------------------------- */
    s8xhx_(nnl, &dx[1], &udx[1], &hdx[1], dxhdx, &iw[1], leniw, &rw[1], lenrw)
	    ;
    if (*dxhdx >= eps) {
/*        ----------------------------------------------- */
/*        Redefine  x  as   x + step*y  (y held in xQP0.) */
/*        Evaluate the functions at the new x. */
/*        ----------------------------------------------- */
	daxpy_(nnl, step, &xqp0[1], &c__1, &x[1], &c__1);
	(*fgwrap)(iexit, &modefg, &nlncon, &nlnobj, n, negcon, nncon0, nncon, 
		nnjac, nnl, nnobj0, nnobj, (U_fp)fgcon, (U_fp)fgobj, &x[1], 
		ne, nlocj, &locj[1], &indj[1], &fcon[1], fobj, &gcon[1], &
		gobj[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	if (*iexit > 0) {
	    return 0;
	}
	if (*iexit == 0 && *usefd) {
	    s6fdg_(iexit, n, negcon, nncon0, nncon, nnjac, nnl, nnobj0, nnobj,
		     (S_fp)fgwrap, (U_fp)fgcon, (U_fp)fgobj, &bl[1], &bu[1], &
		    x[1], ne, nlocj, &locj[1], &indj[1], &fcon[1], fobj, &
		    gcon[1], &gobj[1], &y[1], cu + 8, lencu, &iu[1], leniu, &
		    ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
		     (ftnlen)8, (ftnlen)8);
	    if (*iexit > 0) {
		return 0;
	    }
	}
	if (*iexit == 0) {
/*           ------------------------------------------------------------ */
/*           The functions have been computed at x. */
/*           ------------------------------------------------------------ */
	    if (*nnobj > 0) {
		dcopy_(nnobj, &gobj1[1], &c__1, &dg[1], &c__1);
		daxpy_(nnobj, &c_b16, &gobj[1], &c__1, &dg[1], &c__1);
		if (*minimz < 0) {
		    dscal_(nnobj, &sgnobj, &dg[1], &c__1);
		}
	    }
	    nzero = *nnl - *nnobj;
	    if (nzero > 0) {
		dload_(&nzero, &c_b2, &dg[*nnobj + 1], &c__1);
	    }
	    if (*nncon > 0) {
		s8gprd_(&c__1, &eps0, ne, nlocj, &locj[1], &indj[1], negcon, 
			nlocg, &locg[1], &gcon1[1], &c_b16, &ycon1[1], nncon, 
			&c_b11, &dg[1], nnjac);
		s8gprd_(&c__1, &eps0, ne, nlocj, &locj[1], &indj[1], negcon, 
			nlocg, &locg[1], &gcon[1], &c_b11, &ycon1[1], nncon, &
			c_b11, &dg[1], nnjac);
	    }
	    *ydx = ddot_(nnl, &dg[1], &c__1, &dx[1], &c__1);
	}
    }
    if (*dxhdx < eps || *iexit != 0) {
	dcopy_(nnl, &tdx[1], &c__1, &dx[1], &c__1);
	s8xhx_(nnl, &dx[1], &udx[1], &hdx[1], dxhdx, &iw[1], leniw, &rw[1], 
		lenrw);
    }
    return 0;
} /* s8x1_ */

