/* ./src/sn17util.f -- translated by f2c (version 20100827).
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

static real c_b14 = 1.f;
static doublereal c_b16 = .8;
static doublereal c_b17 = .67;
static doublereal c_b18 = .5;
static doublereal c_b19 = .33;
static doublereal c_b20 = .25;
static doublereal c_b21 = .2;
static integer c__0 = 0;
static integer c__5 = 5;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn17util.f */

/*     Various utility routines for SNOPT. */
/*     The Level 1 and Level 2  BLAS routines are in sn15blas.f */

/*     ddiv     dddiv    ddscl    dnormi   dnormj   dnrm1s   dload */
/*     chcopy   icopy    iload    jdamax */

/*     These could be tuned to the machine being used. */
/*     dload  is used the most. */

/*     ddrand */

/*     Bibs and bobs */

/*     s1init   s1time   s1timp */

/*     s1perm   s1trim */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
doublereal ddiv_(doublereal *a, doublereal *b, logical *fail)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal div, absb, flmin, flmax;
    extern doublereal s1flmn_(void), s1flmx_(void);

/*     ================================================================== */
/*     ddiv  returns the value div given by */

/*     div = ( a/b                 if a/b does not overflow, */
/*           ( */
/*           ( 0.0                 if a .eq. 0.0, */
/*           ( */
/*           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow, */

/*     where  flmax  is a large value, via the function name. */
/*     In addition, if  a/b  would overflow then  fail  is returned as */
/*     true, otherwise  fail  is returned as false. */

/*     Note that when  a and b  are both zero, fail is returned as true, */
/*     but  div  is returned as  0.0. In all other cases of overflow */
/*     div is such that  abs( div ) = flmax. */

/*     When  b = 0  then  sign( a/b )  is taken as  sign( a ). */

/*     15 Nov 1991: First version based on Nag routine f06. */
/*     01 Jun 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    flmax = s1flmx_();
    flmin = s1flmn_();
    if (*a == 0.) {
	div = 0.;
	if (*b == 0.) {
	    *fail = TRUE_;
	} else {
	    *fail = FALSE_;
	}
    } else {
	if (*b == 0.) {
	    div = d_sign(&flmax, a);
	    *fail = TRUE_;
	} else {
	    absb = abs(*b);
	    if (absb >= 1.) {
		*fail = FALSE_;
		if (abs(*a) >= absb * flmin) {
		    div = *a / *b;
		} else {
		    div = 0.;
		}
	    } else {
		if (abs(*a) <= absb * flmax) {
		    *fail = FALSE_;
		    div = *a / *b;
		} else {
		    *fail = TRUE_;
		    div = flmax;
		    if (*a < 0. && *b > 0. || *a > 0. && *b < 0.) {
			div = -div;
		    }
		}
	    }
	}
    }
    ret_val = div;
    return ret_val;
} /* ddiv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* double precision function ddiv */
/* Subroutine */ int dddiv_(integer *n, doublereal *d__, integer *incd, 
	doublereal *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, id, ix;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);

/*     dddiv  performs the diagonal scaling  x  =  x / d. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --d__;

    /* Function Body */
    if (*n > 0) {
	if (*incd == 0 && *incx != 0) {
	    d__1 = 1. / d__[1];
	    i__1 = abs(*incx);
	    dscal_(n, &d__1, &x[1], &i__1);
	} else if (*incd == *incx && *incd > 0) {
	    i__1 = (*n - 1) * *incd + 1;
	    i__2 = *incd;
	    for (id = 1; i__2 < 0 ? id >= i__1 : id <= i__1; id += i__2) {
		x[id] /= d__[id];
	    }
	} else {
	    if (*incx >= 0) {
		ix = 1;
	    } else {
		ix = 1 - (*n - 1) * *incx;
	    }
	    if (*incd > 0) {
		i__2 = (*n - 1) * *incd + 1;
		i__1 = *incd;
		for (id = 1; i__1 < 0 ? id >= i__2 : id <= i__2; id += i__1) {
		    x[ix] /= d__[id];
		    ix += *incx;
		}
	    } else {
		id = 1 - (*n - 1) * *incd;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    x[ix] /= d__[id];
		    id += *incd;
		    ix += *incx;
		}
	    }
	}
    }
    return 0;
} /* dddiv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine dddiv */
/* Subroutine */ int ddscl_(integer *n, doublereal *d__, integer *incd, 
	doublereal *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, id, ix;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);

/*     ------------------------------------------------------------------ */
/*     ddscl  performs the diagonal scaling  x  =  d * x. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --d__;

    /* Function Body */
    if (*n > 0) {
	if (*incd == 0 && *incx != 0) {
	    i__1 = abs(*incx);
	    dscal_(n, &d__[1], &x[1], &i__1);
	} else if (*incd == *incx && *incd > 0) {
	    i__1 = (*n - 1) * *incd + 1;
	    i__2 = *incd;
	    for (id = 1; i__2 < 0 ? id >= i__1 : id <= i__1; id += i__2) {
		x[id] = d__[id] * x[id];
	    }
	} else {
	    if (*incx >= 0) {
		ix = 1;
	    } else {
		ix = 1 - (*n - 1) * *incx;
	    }
	    if (*incd > 0) {
		i__2 = (*n - 1) * *incd + 1;
		i__1 = *incd;
		for (id = 1; i__1 < 0 ? id >= i__2 : id <= i__2; id += i__1) {
		    x[ix] = d__[id] * x[ix];
		    ix += *incx;
		}
	    } else {
		id = 1 - (*n - 1) * *incd;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    x[ix] = d__[id] * x[ix];
		    id += *incd;
		    ix += *incx;
		}
	    }
	}
    }
    return 0;
} /* ddscl_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine ddscl */
doublereal dnormi_(integer *n, doublereal *x, integer *incx)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static integer kmax;
    extern integer idamax_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     dnormi  returns the infinity-norm of the vector  x. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1) {
	ret_val = 0.;
    } else {
	kmax = idamax_(n, &x[1], incx);
	ret_val = (d__1 = x[kmax], abs(d__1));
    }
    return ret_val;
} /* dnormi_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* double precision function dnormi */
doublereal dnormj_(integer *n, doublereal *x, integer *incx)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static integer kmax;
    extern doublereal s1flmx_(void);
    extern integer jdamax_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     dnormj returns the infinity-norm of the vector  x  in most cases. */
/*     flmax is returned if x(*) contains any NaNs or Infs. */

/*     29 Jul 2003: First version of dnormj. */
/*                  Realized that if x(*) contains NaN but not Inf, */
/*                  idamax may return a normal-sized entry, not the NaN. */
/*                  Implemented jdamax and dnormj to test more carefully. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    kmax = jdamax_(n, &x[1], incx);
    if (kmax == 0) {
	ret_val = 0.;
    } else if (kmax > 0) {
	ret_val = (d__1 = x[kmax], abs(d__1));
    } else {
	ret_val = s1flmx_();
    }
    return ret_val;
} /* dnormj_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* double precision function dnormj */
doublereal dnrm1s_(integer *n, doublereal *x, integer *incx)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    extern doublereal dasum_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     dnrm1s  returns the 1-norm of the vector  x,  scaled by root(n). */
/*     This approximates the two-norm of x without the expense. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1) {
	ret_val = 0.;
    } else {
	d__ = (doublereal) (*n);
	d__ = dasum_(n, &x[1], incx) / sqrt(d__);
	ret_val = d__;
    }
    return ret_val;
} /* dnrm1s_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* double precision function dnrm1s */
/* Subroutine */ int dload_(integer *n, doublereal *const__, doublereal *x, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, ix;

/*     ================================================================== */
/*     dload  loads every component of a vector x with the constant. */
/*     Special attention is given to the case incx = 1, const = zero. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n > 0) {
	if (*const__ == 0. && *incx == 1) {
	    m = *n % 7;
	    i__1 = m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = 0.;
	    }
	    i__1 = *n;
	    for (i__ = m + 1; i__ <= i__1; i__ += 7) {
		x[i__] = 0.;
		x[i__ + 1] = 0.;
		x[i__ + 2] = 0.;
		x[i__ + 3] = 0.;
		x[i__ + 4] = 0.;
		x[i__ + 5] = 0.;
		x[i__ + 6] = 0.;
	    }
	} else if (*const__ == 0.) {
	    i__1 = (*n - 1) * *incx + 1;
	    i__2 = *incx;
	    for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
		x[ix] = 0.;
	    }
	} else {
	    i__2 = (*n - 1) * *incx + 1;
	    i__1 = *incx;
	    for (ix = 1; i__1 < 0 ? ix >= i__2 : ix <= i__2; ix += i__1) {
		x[ix] = *const__;
	    }
	}
    }
    return 0;
} /* dload_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine dload */
/* Subroutine */ int chcopy_(integer *n, char *x, integer *incx, char *y, 
	integer *incy, ftnlen x_len, ftnlen y_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, ix, iy;

/*     ================================================================== */
/*     chcopy  is the character*8 version of dcopy. */
/*     ================================================================== */
    /* Parameter adjustments */
    y -= 8;
    x -= 8;

    /* Function Body */
    if (*n > 0) {
	if (*incx == *incy && *incy > 0) {
	    i__1 = (*n - 1) * *incy + 1;
	    i__2 = *incy;
	    for (iy = 1; i__2 < 0 ? iy >= i__1 : iy <= i__1; iy += i__2) {
		s_copy(y + (iy << 3), x + (iy << 3), (ftnlen)8, (ftnlen)8);
	    }
	} else {
	    if (*incx >= 0) {
		ix = 1;
	    } else {
		ix = 1 - (*n - 1) * *incx;
	    }
	    if (*incy > 0) {
		i__2 = (*n - 1) * *incy + 1;
		i__1 = *incy;
		for (iy = 1; i__1 < 0 ? iy >= i__2 : iy <= i__2; iy += i__1) {
		    s_copy(y + (iy << 3), x + (ix << 3), (ftnlen)8, (ftnlen)8)
			    ;
		    ix += *incx;
		}
	    } else {
		iy = 1 - (*n - 1) * *incy;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    s_copy(y + (iy << 3), x + (ix << 3), (ftnlen)8, (ftnlen)8)
			    ;
		    iy += *incy;
		    ix += *incx;
		}
	    }
	}
    }
    return 0;
} /* chcopy_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine chcopy */
/* Subroutine */ int icopy_(integer *n, integer *x, integer *incx, integer *y,
	 integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, ix, iy;

/*     ================================================================== */
/*     icopy  is the integer version of dcopy. */
/*     ================================================================== */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    if (*n > 0) {
	if (*incx == *incy && *incy > 0) {
	    i__1 = (*n - 1) * *incy + 1;
	    i__2 = *incy;
	    for (iy = 1; i__2 < 0 ? iy >= i__1 : iy <= i__1; iy += i__2) {
		y[iy] = x[iy];
	    }
	} else {
	    if (*incx >= 0) {
		ix = 1;
	    } else {
		ix = 1 - (*n - 1) * *incx;
	    }
	    if (*incy > 0) {
		i__2 = (*n - 1) * *incy + 1;
		i__1 = *incy;
		for (iy = 1; i__1 < 0 ? iy >= i__2 : iy <= i__2; iy += i__1) {
		    y[iy] = x[ix];
		    ix += *incx;
		}
	    } else {
		iy = 1 - (*n - 1) * *incy;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = x[ix];
		    iy += *incy;
		    ix += *incx;
		}
	    }
	}
    }
    return 0;
} /* icopy_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine icopy */
/* Subroutine */ int iload_(integer *n, integer *const__, integer *x, integer 
	*incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ix;

/*     iload  loads elements of x with const. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n > 0) {
	if (*incx == 1 && *const__ == 0) {
	    i__1 = *n;
	    for (ix = 1; ix <= i__1; ++ix) {
		x[ix] = 0;
	    }
	} else {
	    i__1 = (*n - 1) * *incx + 1;
	    i__2 = *incx;
	    for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
		x[ix] = *const__;
	    }
	}
    }
    return 0;
} /* iload_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine iload */
integer jdamax_(integer *n, doublereal *x, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, ix;
    static doublereal xi, dmax__;
    static integer kmax;
    static doublereal flmax;
    extern doublereal s1flmx_(void);

/*     ================================================================== */
/*     jdamax does the same as idamax in most cases. */
/*     jdamax > 0 if x contains normal values. */
/*     jdamax = 0 if n = 0. */
/*     jdamax < 0 means x(-jdamax) contains the first NaN or Inf. */

/*     29 Jul 2003: First version of jdamax implemented for s5setx. */
/*     29 Jul 2003: Current version of jdamax */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1) {
	ret_val = 0;
	return ret_val;
    }
    flmax = s1flmx_();
    dmax__ = 0.;
    ix = 1;
    kmax = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = (d__1 = x[ix], abs(d__1));
	if (xi <= flmax) {
/* false if xi = Nan or Inf */
	    if (dmax__ < xi) {
		dmax__ = xi;
		kmax = ix;
	    }
	} else {
	    goto L800;
	}
	ix += *incx;
    }
    ret_val = kmax;
    return ret_val;
L800:
    ret_val = -ix;
    return ret_val;
} /* jdamax_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function jdamax */
/* Subroutine */ int ddrand_(integer *n, doublereal *x, integer *incx, 
	integer *seeds)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double r_mod(real *, real *);

    /* Local variables */
    static integer ix;

/*     ------------------------------------------------------------------ */
/*     ddrand fills a vector x with uniformly distributed random numbers */
/*     in the interval (0, 1) using a method due to  Wichman and Hill. */

/*     seeds(1:3) should be set to integer values */
/*     between 1 and 30000 before the first entry. */

/*     Integer arithmetic up to 30323 is required. */

/*     Blatantly copied from Wichman and Hill 19-January-1987. */
/*     14-Feb-94. Original version. */
/*     30 Jun 1999. seeds stored in an array. */
/*     30 Jun 1999. This version of ddrand. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --seeds;
    --x;

    /* Function Body */
    if (*n < 1) {
	return 0;
    }
    i__1 = (*n - 1) * *incx + 1;
    i__2 = *incx;
    for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	seeds[1] = seeds[1] % 177 * 171 - (seeds[1] / 177 << 1);
	seeds[2] = seeds[2] % 176 * 172 - seeds[2] / 176 * 35;
	seeds[3] = seeds[3] % 178 * 170 - seeds[3] / 178 * 63;
	if (seeds[1] < 0) {
	    seeds[1] += 30269;
	}
	if (seeds[2] < 0) {
	    seeds[2] += 30307;
	}
	if (seeds[3] < 0) {
	    seeds[3] += 30323;
	}
	r__1 = (real) seeds[1] / 30269.f + (real) seeds[2] / 30307.f + (real) 
		seeds[3] / 30323.f;
	x[ix] = r_mod(&r__1, &c_b14);
    }
    return 0;
} /* ddrand_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine ddrand */
/* Subroutine */ int s1init_(char *title, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen title_len)
{
    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal eps, eps0, eps1, eps2, eps3, eps4, eps5;
    extern doublereal s1eps_(void);
    static doublereal flmin, flmax;
    extern doublereal s1flmn_(void), s1flmx_(void);
    extern /* Subroutine */ int s1envt_(integer *, integer *, integer *);
    static doublereal rtundf;

/*     ================================================================== */
/*     s1init saves some machine-dependent constants in iw and rw. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    eps = s1eps_();
    flmax = s1flmx_();
    flmin = s1flmn_();
    rtundf = sqrt(flmin);
/*     Use eps to set other machine precision constants. */
    eps0 = pow_dd(&eps, &c_b16);
    eps1 = pow_dd(&eps, &c_b17);
    eps2 = pow_dd(&eps, &c_b18);
    eps3 = pow_dd(&eps, &c_b19);
    eps4 = pow_dd(&eps, &c_b20);
    eps5 = pow_dd(&eps, &c_b21);
    rw[1] = eps;
/* machine precision */
    rw[2] = eps0;
/* eps**(4/5) */
    rw[3] = eps1;
/* eps**(2/3) */
    rw[4] = eps2;
/* eps**(1/2) */
    rw[5] = eps3;
/* eps**(1/3) */
    rw[6] = eps4;
/* eps**(1/4) */
    rw[7] = eps5;
/* eps**(1/5) */
    rw[8] = flmax;
/* est. of the largest pos. real. */
    rw[9] = flmin;
/* smallest positive real. */
    rw[10] = rtundf;
/*     Set the environment (for later use). */
/* sqrt of flmin. */
    s1envt_(&c__0, &iw[1], leniw);
    return 0;
} /* s1init_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1init */
/* Subroutine */ int s1time_(integer *clock, integer *prtopt, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    extern /* Subroutine */ int s1body_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *);

/*     ================================================================== */
/*     s1time, s1timp and s1cpu are derived from timer, timout and nowcpu */
/*     written for DEC VAX/VMS systems by Irvin Lustig, */
/*     Department of Operations Research, Stanford University, 1987. */

/*     SNOPT  calls s1time only.  s1time calls s1cpu  and  s1timp. */
/*     Only s1cpu is intrinsically machine dependent. */

/*     If a timer is available, call it in s1cpu  and arrange that */
/*     s1cpu  returns the current CPU time in seconds. */

/*     If a timer is not available or not wanted, set time = 0.0 in s1cpu. */
/*     Timing will be turned off and s1timp will not be called. */
/*     ------------------------------------------------------------------ */

/*     s1time turns on or off a selected clock and optionally prints */
/*     statistics regarding all clocks or just the clock chosen. */

/*     The value of abs(clock) is which clock to use. */
/*     If clock  = 0 and prtopt = 0, all clocks and statistics are reset. */
/*     If clock  > 0, the clock is reset to start timing at the */
/*                    current time (determined by calling the */
/*                    machine-dependent subroutine s1cpu). */
/*     If clock  < 0, the clock is turned off and the statistic is */
/*                    recorded for the amount of time since the clock */
/*                    was turned on. */

/*     prtopt is the print option. */
/*     If lvlTim < 0, nothing is printed.  Otherwise, */
/*     prtopt =  0 indicates print nothing, */
/*            =  1 indicates print last time for this clock, */
/*                 only if clock < 0 (it has just been turned off), */
/*            =  2 indicates print total time for all clocks, */
/*            =  3 indicates print mean  time for all clocks. */

/*     The procedure for adding a new timer n is as follows: */
/*     1)  Change ntime to n in the parameter statement in s1time. */
/*     2)  Expand the array "label" to length n in subroutine s1timp. */

/*     04 Jun 1989: Irv's VMS/VAXC version of s1cpu installed, */
/*                  with changes to return time in seconds. */
/*     10 Jul 1992: More clocks added for use in AMPL (and elsewhere). */
/*     31 Jul 2003: snPRNT adopted. */
/*     31 Jul 2003: Current version of s1time. */
/*     ------------------------------------------------------------------ */

/*        Clock 1 is for input time. */
/*        Clock 2 is for solve time. */
/*        Clock 3 is for output time. */
/*        Clock 4 is for the nonlinear functions. */

/*        numt(i)  is the number of times clock i has been turned on. */
/*        tlast(i) is the time at which clock i was last turned on. */
/*        tsum(i)  is the total time elapsed while clock i was on. */
/*        lvlTim   is the Timing level set in the Specs file. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* 1st element of numt(10) */
/* 1st element of tlast(10) */
/* 1st element of tsum(10) */
/*     ------------------------------------------------------------------ */
/* Timing level */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s1body_(clock, prtopt, &c__5, &iw[182], &rw[451], &rw[461], &iw[451], &iw[
	    1], leniw);
    return 0;
} /* s1time_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1time */
/* Subroutine */ int s1body_(integer *clock, integer *prtopt, integer *ntime, 
	integer *lvltim, doublereal *tlast, doublereal *tsum, integer *numt, 
	integer *iw, integer *leniw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ihi, ilo;
    static real time;
    static doublereal stat;
    extern /* Subroutine */ int s1cpu_(integer *, real *);
    static doublereal dtime;
    static integer istat;
    extern /* Subroutine */ int s1timp_(integer *, char *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer iclock;

/*     ================================================================== */
/*     s1body does the work for s1time. */

/*     31 Jul 2003: snPRNT adopted. */
/*     31 Jul 2003: Current version of s1body. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --numt;
    --tsum;
    --tlast;
    --iw;

    /* Function Body */
    if (*lvltim == 0) {
	return 0;
    }
    iclock = abs(*clock);
    if (*clock == 0) {
	if (*prtopt == 0) {
/*           clock = 0, prtopt = 0.  Reset everything. */
	    s1cpu_(&c__1, &time);
	    s1cpu_(&c__0, &time);
	    i__1 = *ntime;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dtime = (doublereal) time;
		tlast[i__] = dtime;
		tsum[i__] = 0.;
		numt[i__] = 0;
	    }
/*           If the s1cpu( 0, time ) gave time < 0.0, we assume that */
/*           the clock is a dummy.  Turn off future timing. */
	    if (time < 0.f) {
		*lvltim = 0;
	    }
	}
    } else {
	s1cpu_(&c__0, &time);
	dtime = (doublereal) time;
	if (*clock > 0) {
	    tlast[iclock] = dtime;
	} else {
	    stat = dtime - tlast[iclock];
	    tsum[iclock] += stat;
	    ++numt[iclock];
	}
    }
/*     Now deal with print options. */
    if (*prtopt == 0 || *lvltim < 0) {
/*        Do nothing. */
    } else if (*prtopt == 1) {
/*        Print statistic for last clock if just turned off. */
	if (*clock < 0) {
	    s1timp_(&iclock, "Last time", &stat, &iw[1], leniw, (ftnlen)9);
	}
    } else {
/*        prtopt >= 2.  Print all statistics if clock = 0, */
/*        or print statistic for individual clock. */
	if (*clock == 0) {
	    s1cpu_(&c_n1, &time);
	    ilo = 1;
	    ihi = *ntime;
	} else {
	    ilo = iclock;
	    ihi = iclock;
	}
	i__1 = ihi;
	for (i__ = ilo; i__ <= i__1; ++i__) {
	    stat = tsum[i__];
	    if (*prtopt == 2) {
		s1timp_(&i__, "Time", &stat, &iw[1], leniw, (ftnlen)4);
	    } else if (*prtopt == 3) {
		istat = numt[i__];
		if (istat > 0) {
		    stat /= istat;
		}
		s1timp_(&i__, "Mean time", &stat, &iw[1], leniw, (ftnlen)9);
	    }
	}
    }
    return 0;
} /* s1body_ */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1body */
/* Subroutine */ int s1timp_(integer *iclock, char *lstat, doublereal *stat, 
	integer *iw, integer *leniw, ftnlen lstat_len)
{
    /* Initialized data */

    static char label[24*5] = "for MPS input           " "for solving proble"
	    "m     " "for solution output     " "for constraint functions" 
	    "for objective function  ";

    /* Format strings */
    static char fmt_1000[] = "(1x,a,f13.2,\002 seconds\002)";

    /* System generated locals */
    address a__1[3];
    integer i__1[3];

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[60], tabby[38];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___52 = { 0, str, 0, fmt_1000, 60, 1 };


/*     ================================================================== */
/*     s1timp  prints CPU time for s1time on file iPrint and/or iSumm. */
/*     It is not intrinsically machine dependent. */

/*     iclock  selects the correct label. */
/*     lstat   is a string to print to tell which type of statistic. */
/*     stat    is the statistic to print out. */
/*             If it is zero, we figure it was never timed, so no print. */

/*     31 Jul 2003: snPRNT adopted. */
/*     05 Dec 2004: Replaced tab edit descriptor. */
/*     05 Dec 2004: Current version of s1timp. */
/*     ================================================================== */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    if (*iclock == 1) {
	snprnt_(&c__3, " ", &iw[1], leniw, (ftnlen)1);
    }
/* Writing concatenation */
    i__1[0] = i_len(lstat, lstat_len), a__1[0] = lstat;
    i__1[1] = 1, a__1[1] = " ";
    i__1[2] = 24, a__1[2] = label + (*iclock - 1) * 24;
    s_cat(tabby, a__1, i__1, &c__3, (ftnlen)38);
    s_wsfi(&io___52);
    do_fio(&c__1, tabby, (ftnlen)38);
    do_fio(&c__1, (char *)&(*stat), (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)60);
    return 0;
} /* s1timp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1timp */
/* Subroutine */ int s1perm_(integer *nkx, integer *kx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;

/*     ================================================================== */
/*     s1Perm  sets the default row and column order for the Jacobian. */

/*     28 Dec 2000: First version written for snopt 5.4. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --kx;

    /* Function Body */
    i__1 = *nkx;
    for (j = 1; j <= i__1; ++j) {
	kx[j] = j;
    }
    return 0;
} /* s1perm_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1Perm */
/* Subroutine */ int s1trim_(char *buffer, integer *lenbuf, ftnlen buffer_len)
{
    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer k;

/*     ================================================================== */
/*     s1trim  returns the length of buffer with trailing blanks omitted. */

/*     02 Dec 2000: First version written for snLog and snLog2. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    *lenbuf = i_len(buffer, buffer_len);
    for (k = *lenbuf; k >= 2; --k) {
	if (*(unsigned char *)&buffer[k - 1] != ' ') {
	    goto L100;
	}
	--(*lenbuf);
    }
L100:
    return 0;
} /* s1trim_ */

