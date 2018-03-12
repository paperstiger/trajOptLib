/* ./src/sn95fmqn.f -- translated by f2c (version 20090411).
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

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn95fmqn.f.   Full memory BFGS routines. */

/*     s9FMH0   s9FMup   s9FMHx */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s9fmh0_(integer *hqntype, integer *nnh, doublereal *
	hcndbd, doublereal *u0pre, doublereal *hd, integer *lenu, doublereal *
	u, integer *nqnmod)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static integer incr;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal hdmin, hdmax;
    static integer nzero;
    static doublereal condhd;
    static logical overfl;

/*     ================================================================== */
/*     s9FMH0 resets  U  such that  H = U'U  to the square root of Hd. */
/*     If U is already diagonal it is set to the identity matrix. */
/*     On entry, the value of HQNType is as follows: */

/*       HQNType */
/*       ------- */
/*       HUnset (-1)      H not set. */
/*       HNorml ( 0)      H is a Hessian of the form defined by  lvlHes. */
/*       HDiag  ( 1)      H is a diagonal matrix. */
/*       HUnit  ( 2)      H is an identity matrix. */

/*     19 Jul 1995: First version of s9FMH0 written by PEG. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     13 Jan 2005: Hd always positive semidefinite. */
/*     14 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hd;
    --u;

    /* Function Body */
    if (*hqntype == 0) {
/*        --------------------------------------------------------------- */
/*        Try and set U so that U'U = Hd, where Hd is the diagonal of the */
/*        Hessian.  First, check the condition of Hd and reset it to */
/*        U0pre*U0pre*I  if its ill-conditioned or not positive definite. */
/*        --------------------------------------------------------------- */
	overfl = FALSE_;
	hdmin = hd[1];
/* strictly positive in exact arithmetic */
	hdmax = hdmin;
	i__1 = *nnh;
	for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MIN */
	    d__1 = hd[i__];
	    hdmin = min(d__1,hdmin);
/* Computing MAX */
	    d__1 = hd[i__];
	    hdmax = max(d__1,hdmax);
	}
	condhd = ddiv_(&hdmax, &hdmin, &overfl);
	if (hdmin <= 0. || condhd >= *hcndbd) {
	    *hqntype = 2;
	}
    }
    if (*hqntype == 0) {
	*hqntype = 1;
/* Set U0 to sqrt(Hd) */
    } else {
	*hqntype = 2;
/* Set U0 to U0pre */
    }
/*     ------------------------------------------------------------ */
/*     Zero the off-diagonal elements of U. */
/*     ------------------------------------------------------------ */
    l = 1;
    incr = *nnh;
    nzero = *nnh - 1;
    i__1 = *nnh;
    for (j = 1; j <= i__1; ++j) {
	if (*hqntype == 2) {
	    u[l] = *u0pre;
	    hd[j] = *u0pre * *u0pre;
	} else {
/* HQNType .eq. HDiag */
	    u[l] = sqrt(hd[j]);
	}
	if (j < *nnh) {
	    dload_(&nzero, &c_b2, &u[l + 1], &c__1);
	    l += incr;
	    --incr;
	    --nzero;
	}
    }
    *nqnmod = 0;
    return 0;
} /* s9fmh0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s9FMH0 */
/* Subroutine */ int s9fmup_(integer *update, integer *hqntype, integer *nnh, 
	integer *mqnmod, integer *nqnmod, doublereal *hcndbd, doublereal *
	u0pre, doublereal *u0scal, doublereal *rydx, doublereal *rdxhdx, 
	doublereal *hd, doublereal *hdx, doublereal *y, doublereal *udx, 
	integer *lenu, doublereal *u)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal t, yi;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal hdxi, told;
    static integer numu;
    static doublereal tolz;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s9fmh0_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer iexit;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal ulast, h0scal;
    extern /* Subroutine */ int s6rmod_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer lastnz;

/*     ================================================================== */
/*     s9FMup applies the full-memory BFGS update to H = U'U. */
/*     If defined, the self-scaling BFGS update parameter is saved. */
/*     It is needed to update the reduced Hessian when there are only */
/*     linear constraints. */

/*     19 Jul 1995: First version of s9FMup. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     13 Jan 2005: FM H stored in product form. */
/*     15 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --udx;
    --y;
    --hdx;
    --hd;
    --u;

    /* Function Body */
    told = 0.;
    tolz = 0.;
    if (*update > 1) {
	*u0scal = *rydx / *rdxhdx;
	h0scal = *u0scal * *u0scal;
/* Self-scaling parameter */
	*rdxhdx *= *u0scal;
/* Used again for the LC update. */
	numu = *nnh * (*nnh + 1) / 2;
	dscal_(&numu, u0scal, &u[1], &c__1);
/* multiplies U by U0scal. */
	dscal_(nnh, u0scal, &udx[1], &c__1);
	dscal_(nnh, &h0scal, &hd[1], &c__1);
	dscal_(nnh, &h0scal, &hdx[1], &c__1);
    }
/*     Include the latest update in the Hessian diagonal. */
    i__1 = *nnh;
    for (i__ = 1; i__ <= i__1; ++i__) {
	hdxi = hdx[i__] / *rdxhdx;
	yi = y[i__] / *rydx;
/* Computing 2nd power */
	d__1 = hdxi;
/* Computing 2nd power */
	d__2 = yi;
	hd[i__] = hd[i__] - d__1 * d__1 + d__2 * d__2;
    }
    if (*nqnmod >= *mqnmod) {
/*        --------------------------------------------------------------- */
/*        Reset H = U'U to be the diagonal of the current H. */
/*        Discard any updates accumulated so far. */
/*        --------------------------------------------------------------- */
	s9fmh0_(hqntype, nnh, hcndbd, u0pre, &hd[1], lenu, &u[1], nqnmod);
    } else {
/*        --------------------------------------------------------------- */
/*        Overwrite (Udx,y) with the vectors (Us,v)  such that */
/*        Us  = Udx / rdxHdx,    v = (1/rydx) gdif - (1/rdxHdx) Hdx. */

/*        Then, U(new) = U + Us v',  with H = U'U. */

/*        Hdx and v  are saved to update R for LC problems. */
/*        --------------------------------------------------------------- */
	++(*nqnmod);
	t = ddot_(nnh, &y[1], &c__1, &hdx[1], &c__1);
	if (t >= 0.) {
	    d__1 = 1. / *rydx;
	    dscal_(nnh, &d__1, &y[1], &c__1);
	} else {
	    d__1 = -1. / *rydx;
	    dscal_(nnh, &d__1, &y[1], &c__1);
	}
	d__1 = -1. / *rdxhdx;
	daxpy_(nnh, &d__1, &hdx[1], &c__1, &y[1], &c__1);
	d__1 = 1. / *rdxhdx;
	dscal_(nnh, &d__1, &udx[1], &c__1);
/*        --------------------------------------------------------------- */
/*        Restore  U + Us y' to triangular form  (overwriting Udx). */
/*        --------------------------------------------------------------- */
	ulast = 0.;
	lastnz = *nnh;
	s6rmod_(&iexit, nnh, nnh, lenu, &u[1], &udx[1], &y[1], &lastnz, &
		ulast, &told, &tolz);
	*hqntype = 0;
    }
    return 0;
} /* s9fmup_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s9FMup */
/* Subroutine */ int s9fmhx_(integer *nnh, doublereal *x, doublereal *ux, 
	doublereal *hx, integer *lenu, doublereal *u)
{
    extern /* Subroutine */ int s6rprd_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s9FMHx  computes the product Hx = U'Ux, where  U is an upper- */
/*     triangular matrix stored by rows in the one-dimensional array  U. */
/*     lenU defines the length of U.  lenU must be at least */
/*     nnH*(nnH + 1)/2. */

/*     12 Jan 1996: First version of s9FMHx */
/*     12 Jan 2005: H held as U'U. */
/*     15 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hx;
    --ux;
    --x;
    --u;

    /* Function Body */
    s6rprd_(&c__0, nnh, nnh, lenu, &u[1], &x[1], &ux[1]);
    s6rprd_(&c__1, nnh, nnh, lenu, &u[1], &ux[1], &hx[1]);
    return 0;
} /* s9fmhx_ */

