/* ./src/sn90lmqn.f -- translated by f2c (version 20090411).
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

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn90lmqn.f.   Limited-memory BFGS routines. */

/*     s9LMH0   s9LMup   s9LMHx */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s9lmh0_(integer *hqntype, integer *nnh, doublereal *
	hcndbd, doublereal *u0pre, doublereal *hd, doublereal *u0, integer *
	nqnmod)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal hdmin, hdmax, condhd;
    static logical overfl;

/*     ================================================================== */
/*     s9LMH0 sets the diagonal U0  such that  H = U0'U0. */
/*     If the BFGS diagonal Hd is not ill-conditioned, U0 = sqrt(Hd). */
/*     Otherwise, U0 =  U0pre*I and Hd = U0'U0. */

/*     On entry, the value of HQNType is as follows: */

/*       HQNType */
/*       ------- */
/*       HUnset (-1)      H not set. */
/*       HNorml ( 0)      H is a Hessian of the form defined by  lvlHes. */
/*       HDiag  ( 1)      H is a diagonal matrix. */
/*       HUnit  ( 2)      H is a scaled identity matrix. */

/*     19 Jul 1995: First version of s9LMH0. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     18 Feb 2001: H stored in product form. */
/*     13 Jan 2005: Hd always positive semidefinite. */
/*     15 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --u0;
    --hd;

    /* Function Body */
    if (*hqntype == 0) {
/*        --------------------------------------------------------------- */
/*        Check  Hd, the diagonal of the current BFGS  H. */
/*        Reset  Hd = U0pre*I  if Hd is not positive definite or is ill */
/*        conditioned */
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
/*        --------------------------------------------------------------- */
/*        Set U0 to sqrt(Hd) */
/*        --------------------------------------------------------------- */
	*hqntype = 1;
	i__1 = *nnh;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u0[i__] = sqrt(hd[i__]);
	}
    } else {
/*        --------------------------------------------------------------- */
/*        Set U0 to a multiple of the identity. */
/*        We come straight here for  HUnset. */
/*        --------------------------------------------------------------- */
	*hqntype = 2;
	dload_(nnh, u0pre, &u0[1], &c__1);
	d__1 = *u0pre * *u0pre;
	dload_(nnh, &d__1, &hd[1], &c__1);
    }
    *nqnmod = 0;
    return 0;
} /* s9lmh0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s9LMH0 */
/* Subroutine */ int s9lmup_(integer *update, integer *hqntype, integer *nnh, 
	integer *mqnmod, integer *nqnmod, doublereal *hcndbd, doublereal *
	u0pre, doublereal *u0scal, doublereal *rydx, doublereal *rdxhdx, 
	doublereal *hd, doublereal *hdx, doublereal *y, doublereal *dx, 
	doublereal *u0, doublereal *s, doublereal *v)
{
    /* System generated locals */
    integer s_dim1, s_offset, v_dim1, v_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal t, yi;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal hdxi;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s9lmh0_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), dcopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal h0scal;

/*     ================================================================== */
/*     s9LMup does almost everything associated with the limited-memory */
/*     quasi-Newton update. */

/*     If defined, the self-scaling BFGS parameter U0scl is saved. */
/*     It is needed to update the reduced Hessian when there are only */
/*     linear constraints. */

/*     19 Jul 1995: First version of s9LMup. */
/*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added. */
/*     18 Feb 2001: LM H stored in product form. */
/*     13 Jan 2005: FM H stored in product form. */
/*     16 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --u0;
    --dx;
    --y;
    --hdx;
    --hd;
    v_dim1 = *nnh;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    s_dim1 = *nnh;
    s_offset = 1 + s_dim1;
    s -= s_offset;

    /* Function Body */
    if (*update > 1) {
	*u0scal = *rydx / *rdxhdx;
	h0scal = *u0scal * *u0scal;
	*rdxhdx *= *u0scal;
/* Used later for the LC update. */
	dscal_(nnh, u0scal, &u0[1], &c__1);
/* multiplies U by U0scal. */
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
/*        Insufficient space for storing the new update pair (s,v). */
/*        Reset U0 to be the root of the diagonal of the current H. */
/*        Change nQNmod to discard any updates accumulated so far. */
/*        --------------------------------------------------------------- */
	s9lmh0_(hqntype, nnh, hcndbd, u0pre, &hd[1], &u0[1], nqnmod);
    } else {
/*        --------------------------------------------------------------- */
/*        Space remains. Store s and v, where */
/*        U(new) = U(I + sv'), with H = U'U. */
/*        S  =  dx / rdxHdx,    V =  (1/rydx) gdif - (1/rdxHdx) Hdx. */

/*        Hdx and the modified y (= v) are used to update the reduced */
/*        Hessian for LC problems. */
/*        --------------------------------------------------------------- */
	++(*nqnmod);
	dcopy_(nnh, &dx[1], &c__1, &s[*nqnmod * s_dim1 + 1], &c__1);
	d__1 = 1. / *rdxhdx;
	dscal_(nnh, &d__1, &s[*nqnmod * s_dim1 + 1], &c__1);
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
	dcopy_(nnh, &y[1], &c__1, &v[*nqnmod * v_dim1 + 1], &c__1);
	*hqntype = 0;
    }
    return 0;
} /* s9lmup_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s9LMup */
/* Subroutine */ int s9lmhx_(integer *nnh, doublereal *x, doublereal *ux, 
	doublereal *hx, integer *mqnmod, integer *nqnmod, doublereal *u0, 
	doublereal *s, doublereal *v)
{
    /* System generated locals */
    integer s_dim1, s_offset, v_dim1, v_offset, i__1;

    /* Local variables */
    static doublereal c__;
    static integer k;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s9LMHx forms the product  Hx  for the limited-memory s8Hx. */
/*     H = U'U, where U = U0*(I + s1*v1')*(I + s2*v2')...(I + sk*vk'). */
/*     with  k = nQNmod */

/*     19 Jul 1995: First version of s9LMHx */
/*     18 Feb 2001: H stored in product form. */
/*     12 Jan 2005: Ux added as argument. */
/*     16 Jan 2005: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --u0;
    --hx;
    --ux;
    --x;
    v_dim1 = *nnh;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    s_dim1 = *nnh;
    s_offset = 1 + s_dim1;
    s -= s_offset;

    /* Function Body */
    dcopy_(nnh, &x[1], &c__1, &ux[1], &c__1);
/*     Multiply by U. */
    for (k = *nqnmod; k >= 1; --k) {
	c__ = ddot_(nnh, &v[k * v_dim1 + 1], &c__1, &ux[1], &c__1);
	daxpy_(nnh, &c__, &s[k * s_dim1 + 1], &c__1, &ux[1], &c__1);
    }
    ddscl_(nnh, &u0[1], &c__1, &ux[1], &c__1);
/*     Multiply by U'. */
    dcopy_(nnh, &ux[1], &c__1, &hx[1], &c__1);
    ddscl_(nnh, &u0[1], &c__1, &hx[1], &c__1);
    i__1 = *nqnmod;
    for (k = 1; k <= i__1; ++k) {
	c__ = ddot_(nnh, &s[k * s_dim1 + 1], &c__1, &hx[1], &c__1);
	daxpy_(nnh, &c__, &v[k * v_dim1 + 1], &c__1, &hx[1], &c__1);
    }
    return 0;
} /* s9lmhx_ */

