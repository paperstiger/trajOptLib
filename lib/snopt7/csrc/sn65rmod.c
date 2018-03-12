/* ./src/sn65rmod.f -- translated by f2c (version 20100827).
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
static doublereal c_b13 = 1.;
static doublereal c_b29 = 0.;
static integer c__2 = 2;
static doublereal c_b50 = -1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  sn65rmod.f */

/*     s6Radd   s6RBFS   s6chol   s6mchl   s6pchl   s6Rcnd   s6Rcol */
/*     s6Rdel   s6Rfix   s6Rmod   s6Rprd   s6Rrow   s6Rset   s6Rsol */
/*     s6Rswp   s6Rupd   s6Rqn */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s6radd_(integer *maxr, integer *lenr, integer *ns, 
	doublereal *r__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, lr, incr;

/* ================================================================= */
/* s6Radd  adds a unit column nS to the upper triangular matrix R. */

/* 12 Jun 2001: First version based on MINOS routine m6radd. */
/* 12 Jun 2001: Current version of s6Radd. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;

    /* Function Body */
    if (*ns <= *maxr) {
	lr = *ns;
	incr = *maxr;
	i__1 = *ns - 1;
	for (k = 1; k <= i__1; ++k) {
	    r__[lr] = 0.;
	    --incr;
	    lr += incr;
	}
    } else {
	lr = *maxr * (*maxr + 1) / 2 + (*ns - *maxr);
    }
    r__[lr] = 1.;
    return 0;
} /* s6radd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Radd */
/* Subroutine */ int s6rbfs_(integer *iexit, integer *maxr, integer *ns, 
	integer *n, integer *lenr, doublereal *told, doublereal *tolz, 
	doublereal *r__, doublereal *u, doublereal *v, doublereal *delta, 
	doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, l, j1, nr;
    static doublereal unz;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal ulast;
    extern /* Subroutine */ int s6rmod_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer lastnz;

/* ================================================================= */
/* s6RBFS applies the quasi-Newton update to the first nS rows and */
/* columns of the n x n upper-triangular  R  such that H = R'R. */

/* On entry: */
/* R   contains the first nS rows and columns of an n x n */
/*     upper-triangular matrix R such that R'R = Q'HQ. */

/* u   contains the first  nS  components of  Rdx/norm(Rdx), where */
/*     R'*(Rdx) = Q'Hdx.  Note that norm(u)=1  if  nS = n. */
/*     It is overwritten. */

/* v   contains the first nS components of the BFGS update vector */
/*     v such that U(new) = U(I + sv'), with H = U'U. */
/*     s  = dx / rdxHdx, v = (1/rydx) gdif - (1/rdxHdx) Hdx. */

/* (delta, w)  is a scalar-vector pair such that */
/*       delta*w = (1/rdxHdx) Hdx */
/*               = U'U. */

/* 30 Dec 1991: First version based on NPSOL routine npbfgs. */
/* 05 May 1999: Last column-wise version. */
/* 03 Dec 2000: Converted to row-wise storage. */
/* 18 Feb 2001: H stored in product form. */
/* 23 Jul 2001: Excess elements of R implemented. */
/* 17 Nov 2001: Current version of s6RBFS. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --w;
    --v;
    --u;
    --r__;

    /* Function Body */
    nr = min(*ns,*maxr);
/* Apply the update in the form  R + u v',  where  u  is held */
/* in u and ulast. */
    ulast = 0.;
    if (*ns < *n) {
/* Computing MAX */
/* Computing 2nd power */
	d__3 = dnrm2_(ns, &u[1], &c__1);
	d__1 = 0., d__2 = 1. - d__3 * d__3;
	ulast = sqrt((max(d__1,d__2)));
    }
    if (*ns > *maxr) {
	i__1 = *ns - *maxr;
/* Computing 2nd power */
	d__1 = dnrm2_(&i__1, &u[*maxr + 1], &c__1);
/* Computing 2nd power */
	d__2 = ulast;
	ulast = sqrt(d__1 * d__1 + d__2 * d__2);
    }
/* -------------------------- */
/* Find the last nonzero in u. */
/* -------------------------- */
    unz = ulast;
    lastnz = nr + 1;
/* +    while (lastnz .gt. 1  .and.  unz .le. tolz) do */
L100:
    if (lastnz > 1 && unz <= *tolz) {
	--lastnz;
	unz = (d__1 = u[lastnz], abs(d__1));
	goto L100;
/* +    end while */
    }
/* ------------------------------------------------------ */
/* Restore  R + u v'  to triangular form  (overwriting u). */
/* ------------------------------------------------------ */
    s6rmod_(iexit, maxr, &nr, lenr, &r__[1], &u[1], &v[1], &lastnz, &ulast, 
	    told, tolz);
/* Deal with surplus diagonals of  R. */
/* They are defined to be the diagonals of the rank-two */
/* modification:  v w' + w v' + v v' (where w is already */
/* scaled by delta). */
    if (*ns > *maxr) {
	j1 = *maxr + 1;
	l = *maxr * j1 / 2;
	i__1 = *ns;
	for (j = j1; j <= i__1; ++j) {
	    ++l;
/* Rdsq = max(told,(R(l)**2 + v(j)**2 + two*delta*w(j)*v(j))) */
/* R(l) = sqrt(Rdsq) */
	    r__[l] = 1.;
	}
    }
    return 0;
} /* s6rbfs_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6RBFS */
/* Subroutine */ int s6chol_(integer *iexit, integer *swap, integer *maxh, 
	integer *nh, integer *lenh, doublereal *h__, doublereal *hdmin, 
	doublereal *dmax__, integer *irank, integer *perm)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublereal s;
    static integer is, js, ls, jx, lcj, lck, lrj, lrk;
    static doublereal hjx;
    static integer incj, inck, incr, incs, kmax, lmax, ldiag;
    static doublereal dsmall;
    static logical permut;

/* ================================================================= */
/* s6chol  forms the upper-triangular Cholesky factor R such that */
/* H = R'R  or  P H P' = R'R for some permutation P. */

/* On entry, */
/*  Swap   specifies the pivoting strategy. */
/*  Swap=0 means let P = I (no interchanges). */
/*         Otherwise P chooses the maximum diagonal at each stage. */

/*   Hdmin is the square of the smallest allowable diagonal of R. */

/* On exit, */
/*  perm   contains details of the permutation matrix P, such that */
/*         perm(k) = k  if no interchange occurred at the kth step, */
/*         perm(k) = j  (k < j <= nH)  if rows k and j were */
/*                        interchanged at the kth step. */

/* Only the diagonal and super-diagonal elements of H are used, */
/* and they are overwritten by R. */

/* 28 May 1994: First version of s6chol based on routine chlfac. */
/* 15 Oct 2000: Last version with column-wise storage for R. */
/* 05 Dec 2000: Converted to row-wise storage. */
/* 06 Dec 2000: Mild threshold diagonal pivoting implemented */
/*              to reduce frequency of symmetric interchanges. */
/*              TDPtol < 1.0 disables TDP. */
/* 25 Mar 2001: Current version of s6chol. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --perm;
    --h__;

    /* Function Body */
    *iexit = 0;
    *irank = 0;
    if (*nh == 0) {
	return 0;
    }
    permut = *swap != 0;
    dsmall = *hdmin;
/* ----------------------------------------------------------------- */
/* Form the Cholesky factorization R'R = H or P H P'. */
/* Process the first nH rows of H. */

/* ldiag is the location of H(j,j). */
/* lmax  is the location of H(k,k) if j and k will be interchanged. */
/* ----------------------------------------------------------------- */
    ldiag = 1;
    incr = *maxh;
    i__1 = *nh;
    for (j = 1; j <= i__1; ++j) {
	*dmax__ = h__[ldiag];
	kmax = j;
	lmax = ldiag;
	if (permut) {
/* Find the largest diagonal of the Schur complement. */
	    ls = ldiag + incr;
	    incs = incr - 1;
	    i__2 = *nh;
	    for (k = j + 1; k <= i__2; ++k) {
		if (*dmax__ < h__[ls]) {
		    *dmax__ = h__[ls];
		    kmax = k;
		    lmax = ls;
		}
		ls += incs;
		--incs;
	    }
	}
	if (*dmax__ <= dsmall) {
	    goto L800;
	}
/* The diagonal is too small. */
	if (h__[ldiag] * 1.1 >= *dmax__) {
	    kmax = j;
/* Don't interchange after all. */
	    lmax = ldiag;
	    *dmax__ = h__[ldiag];
	}
	perm[j] = kmax;
	if (kmax != j) {
/* Perform a symmetric interchange. */
/* First swap row H(j,j+1:kmax) with col H(j+1:kmax,kmax). */
	    lrj = ldiag + 1;
	    lck = ldiag + kmax - j + incr - 1;
	    inck = *maxh - j - 1;
	    i__2 = kmax;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		s = h__[lrj];
		h__[lrj] = h__[lck];
		h__[lck] = s;
		++lrj;
		lck += inck;
		--inck;
	    }
/* Now swap col H(1:j,j) with col H(1:j,kmax) */
	    lcj = j;
	    lck = kmax;
	    incj = *maxh - 1;
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s = h__[lcj];
		h__[lcj] = h__[lck];
		h__[lck] = s;
		lcj += incj;
		lck += incj;
		--incj;
	    }
/* Finally swap row H(j,kmax:n) with row H(kmax,kmax:n). */
	    lrj = ldiag + kmax - j;
	    lrk = lmax;
	    i__2 = *nh;
	    for (i__ = kmax; i__ <= i__2; ++i__) {
		s = h__[lrj];
		h__[lrj] = h__[lrk];
		h__[lrk] = s;
		++lrj;
		++lrk;
	    }
	}
/* Set the diagonal of  R. */
/* kmax ne j */
	d__ = sqrt(*dmax__);
	h__[ldiag] = d__;
	++(*irank);
	if (j < *nh) {
/* Set the super-diagonal elements of the jth row of R. */
	    jx = ldiag + 1;
	    i__2 = *nh;
	    for (k = j + 1; k <= i__2; ++k) {
		h__[jx] /= d__;
		++jx;
	    }
/* Do a rank-one update to the Schur complement. */
/* Form the upper-triangular part of H = H - x x', */
/* where x is the row H(j,j+1:nH). */
/* H(js,:) = H(js,:) - H(j,js)*H(j,:). */
	    jx = ldiag + 1;
	    ls = ldiag + incr;
	    incs = incr - 1;
	    i__2 = *nh;
	    for (js = j + 1; js <= i__2; ++js) {
		hjx = h__[jx];
		if (hjx != 0.) {
		    i__ = jx;
		    k = ls;
		    i__3 = *nh;
		    for (is = js; is <= i__3; ++is) {
			h__[k] -= hjx * h__[i__];
			++i__;
			++k;
		    }
		}
		++jx;
		ls += incs;
		--incs;
	    }
	}
/* rank-one update */
	ldiag += incr;
	--incr;
    }
/* ================================================================ */
/* Test if  H  is not positive definite. */
/* ================================================================ */
L800:
    if (*irank < *nh) {
	*iexit = 1;
    }
    return 0;
} /* s6chol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6chol */
/* Subroutine */ int s6mchl_(integer *iexit, integer *swap, integer *maxh, 
	integer *nh, integer *lenh, doublereal *h__, doublereal *hdmin, 
	doublereal *eps, doublereal *dmax__, integer *irank, integer *nmodh, 
	integer *perm, doublereal *e)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublereal s;
    static integer is, jr, js, lr, ls, jx, lcj, lck, lrj, lrk;
    static doublereal hjx;
    static integer incj, inck, incr, incs, kmax, lmax;
    static doublereal supj;
    static integer ldiag;
    static doublereal dsmall, betasq;
    static logical permut;
    static doublereal supmax;

/* ================================================================= */
/* s6mchl  forms a modified Cholesky factorization of the matrix H, */
/* such that */
/*            H + E = R'R  or  P (H + E) P' = R'R, */
/* where */
/*            E is diagonal, */
/*            R is upper triangular, */
/*            P is a permutation. */
/* If H is sufficiently positive definite, E will be zero. */

/* On entry, */
/*  Swap   specifies the pivoting strategy. */
/*  Swap   = 0 means let P = I (no interchanges). */
/*         Otherwise P chooses the maximum diagonal at each stage. */

/*  Hdmin  is the square of the smallest allowable diagonal of R. */

/*  eps    causes termination if the remaining diagonals are */
/*         that small.  Intended to leave singularities unaltered. */

/* On exit, */
/*  iExit = 0 if iRank = nH, */
/*        = 1 if iRank < nH (so the returned factor is singular). */

/*  dmax   returns the largest diagonal AT THE TIME OF TERMINATION. */
/*         (More explanation needed.) */

/*  iRank  is the number of rows in R. */

/*  nmodH  is the number of modified diagonals out of 1:iRank. */
/*         Note that diagonals iRank+1:nH probably need to be */
/*         taken care of by the calling routine. */

/*  perm   contains details of the permutation matrix P, such that */
/*         perm(j) = j  if no interchange occurred at the kth step, */
/*         perm(j) = k  (j < k <= nH)  if rows j and k were */
/*                       interchanged at the j-th step. */

/*  e      contains the modifications to each diagonal: */
/*         e(1:nH) >= 0.0,  E = diag(e). */

/* Only the diagonal and super-diagonal elements of H are used, */
/* and they are overwritten by R. */

/* 21 Oct 1993: s6mchl: First version based on routine chlfac. */
/* 15 Oct 1994: s6mchl: Symmetric interchanges added. */
/* 28 May 1999: s6mchl: Last version with column-wise R. */
/* 05 Dec 2000: s6chol: "Unmodified Cholesky" routine */
/*                      converted to row-wise storage. */
/* 06 Dec 2000: s6chol: Theshold Diagonal Pivoting implemented. */
/* 20 Dec 2000: s6mchl: Derived from s6chol and previous s6mchl. */
/* 25 Jan 2005: m6mchl = s6mchl from SNOPT. */
/* 25 Jan 2005: perm not touched if Permut = .false. */
/* 12 Mar 2005: nmodH is an output parameter (no longer input). */
/* 15 Mar 2005: Bug found by MPF in first loop.  incR is needed. */
/*              dmax and supmax were not being set properly. */
/* 23 Mar 2005: s6mchl = m6mchl from MINOS6. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --e;
    --perm;
    --h__;

    /* Function Body */
    *iexit = 0;
    *irank = 0;
    *nmodh = 0;
    if (*nh == 0) {
	return 0;
    }
    permut = *swap != 0;
/* ----------------------------------------------------------------- */
/* Find the maximum diagonal and super-diagonal elements of H. */
/* TDPtol = 1.1 (say) implements mild */
/* threshold diagonal pivoting (TDP) */
/* to reduce frequency of symmetric interchanges. */
/* TDPtol < 1.0 disables TDP. */
/* ----------------------------------------------------------------- */
    dsmall = *hdmin;
    *dmax__ = dsmall;
    supmax = 0.;
    ldiag = 1;
    incr = *maxh;
    i__1 = *nh;
    for (j = 1; j <= i__1; ++j) {
	e[j] = 0.;
	lr = ldiag;
/* Computing MAX */
	d__2 = *dmax__, d__3 = (d__1 = h__[lr], abs(d__1));
	*dmax__ = max(d__2,d__3);
	i__2 = *nh;
	for (k = j + 1; k <= i__2; ++k) {
	    ++lr;
/* Computing MAX */
	    d__2 = supmax, d__3 = (d__1 = h__[lr], abs(d__1));
	    supmax = max(d__2,d__3);
	}
	ldiag += incr;
	--incr;
    }
/* Computing MAX */
    d__1 = *dmax__, d__2 = supmax / *nh;
    betasq = max(d__1,d__2);
/* ----------------------------------------------------------------- */
/* Form the Cholesky factorization R'R = H or P H P'. */
/* Process the first nH rows of H. */

/* ldiag is the location of H(j,j). */
/* lmax  is the location of H(k,k) if j and k will be interchanged. */
/* ----------------------------------------------------------------- */
/* Bound on the off-diagonals */
    ldiag = 1;
    incr = *maxh;
    i__1 = *nh;
    for (j = 1; j <= i__1; ++j) {
	*dmax__ = (d__1 = h__[ldiag], abs(d__1));
	kmax = j;
	lmax = ldiag;
	if (permut) {
/* Find the diagonal of the Schur complement with */
/* maximum absolute value. */
	    ls = ldiag + incr;
	    incs = incr - 1;
	    i__2 = *nh;
	    for (k = j + 1; k <= i__2; ++k) {
		if (*dmax__ < (d__1 = h__[ls], abs(d__1))) {
		    *dmax__ = (d__1 = h__[ls], abs(d__1));
		    kmax = k;
		    lmax = ls;
		}
		ls += incs;
		--incs;
	    }
	}
	if (*dmax__ <= *eps) {
	    goto L800;
	}
/* The diagonal is too small. */
	if ((d__1 = h__[ldiag], abs(d__1)) * 1.1 >= *dmax__) {
	    kmax = j;
/* Don't interchange after all. */
	    lmax = ldiag;
	}
	*dmax__ = h__[lmax];
/* NOTE: not abs(.) anymore. */
	if (permut) {
	    perm[j] = kmax;
	}
/* --------------------------------------------------------- */
/* Perform a symmetric interchange. */
/* --------------------------------------------------------- */
	if (kmax != j) {
/* First swap row H(j,j+1:kmax) with col H(j+1:kmax,kmax). */
	    lrj = ldiag + 1;
	    lck = ldiag + kmax - j + incr - 1;
	    inck = *maxh - j - 1;
	    i__2 = kmax;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		s = h__[lrj];
		h__[lrj] = h__[lck];
		h__[lck] = s;
		++lrj;
		lck += inck;
		--inck;
	    }
/* Now swap col H(1:j,j) with col H(1:j,kmax) */
	    lcj = j;
	    lck = kmax;
	    incj = *maxh - 1;
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s = h__[lcj];
		h__[lcj] = h__[lck];
		h__[lck] = s;
		lcj += incj;
		lck += incj;
		--incj;
	    }
/* Finally swap row H(j,kmax:n) with row H(kmax,kmax:n). */
	    lrj = ldiag + kmax - j;
	    lrk = lmax;
	    i__2 = *nh;
	    for (i__ = kmax; i__ <= i__2; ++i__) {
		s = h__[lrj];
		h__[lrj] = h__[lrk];
		h__[lrk] = s;
		++lrj;
		++lrk;
	    }
	}
/* --------------------------------------------------------- */
/* Find the largest super-diagonal in the jth row. */
/* --------------------------------------------------------- */
/* kmax ne j */
	supj = 0.;
	jr = ldiag + 1;
	i__2 = *nh;
	for (k = j + 1; k <= i__2; ++k) {
/* Computing MAX */
	    d__2 = supj, d__3 = (d__1 = h__[jr], abs(d__1));
	    supj = max(d__2,d__3);
	    ++jr;
	}
	if (supj < dsmall) {
	    supj = 0.;
	}
/* --------------------------------------------------------- */
/* Set the diagonal of  R. */
/* --------------------------------------------------------- */
/* Computing MAX */
/* Computing 2nd power */
	d__3 = supj;
	d__1 = dsmall, d__2 = abs(*dmax__), d__1 = max(d__1,d__2), d__2 = 
		d__3 * d__3 / betasq;
	d__ = max(d__1,d__2);
	e[j] = d__ - *dmax__;
	if (e[j] > 0.) {
	    ++(*nmodh);
	}
	d__ = sqrt(d__);
	h__[ldiag] = d__;
	++(*irank);
	if (j < *nh) {
/* Set the super-diagonal elements of the jth row of R. */
	    jr = ldiag + 1;
	    i__2 = *nh;
	    for (k = j + 1; k <= i__2; ++k) {
		h__[jr] /= d__;
		++jr;
	    }
/* Do a rank-one update to the Schur complement. */
/* Form the upper-triangular part of H = H - x x', */
/* where x is the row H(j,j+1:nH). */
/* H(js,:) = H(js,:) - H(j,js)*H(j,:). */
	    jx = ldiag + 1;
	    ls = ldiag + incr;
	    incs = incr - 1;
	    i__2 = *nh;
	    for (js = j + 1; js <= i__2; ++js) {
		hjx = h__[jx];
		if (hjx != 0.) {
		    i__ = jx;
		    k = ls;
		    i__3 = *nh;
		    for (is = js; is <= i__3; ++is) {
			h__[k] -= hjx * h__[i__];
			++i__;
			++k;
		    }
		}
		++jx;
		ls += incs;
		--incs;
	    }
	}
/* rank-one update */
	ldiag += incr;
	--incr;
    }
/* ================================================================ */
/* Test if  H  is not positive definite. */
/* ================================================================ */
L800:
    if (*irank < *nh) {
	*iexit = 1;
    }
    return 0;
} /* s6mchl_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6mchl */
/* Subroutine */ int s6pchl_(integer *iexit, doublereal *tolpiv, integer *
	maxh, integer *nh, integer *lenh, doublereal *h__, integer *npos, 
	integer *perm)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublereal s;
    static integer is, jr, js, ls, jx, lcj, lck, lrj, lrk;
    static doublereal hjx;
    static integer incj, inck;
    static doublereal dmax__;
    static integer incr, incs, kmax, lmax, ldiag;
    static doublereal supmax;

/* ================================================================= */
/* s6pchl  forms the partial Cholesky factorization: */
/*            P H P' = R'DR, */
/* where */
/*            H is symmetric, */
/*            D is the block diagonal, with  D = (  I   0   ) */
/*                                               (  0   B   ), */
/*            R is upper triangular, with    R = ( R11 R12 ) */
/*                                               (  0   I  ), */
/*            P is a permutation. */
/* If H is sufficiently positive definite, D  will be the identity. */

/* On entry, */

/*  H      is the upper-triangular part of the maxH by maxH */
/*         matrix  H, stored by rows in H(1:lenH). */

/*  tolPiv is the relative pivot tolerance. */

/* On exit, */

/*  iExit = 0 if H is sufficiently pos. definite with nPos = nH, */
/*        = 1 if H has at least one non-positive eigenvalue. */

/*  nPos   is the dimension of the non-unit part of R. */
/*         If H is sufficiently positive definite,  nPos = nH. */
/*         Otherwise, B has at least one negative eigenvalue. */

/*  perm   contains details of the permutation matrix P, such that */
/*         perm(j) = j  if no interchange occurred at the kth step, */
/*         perm(j) = k  (j < k <= nH)  if rows j and k were */
/*                        interchanged at the jth step. */

/* The elements of  H are overwritten by  R  and  D. */

/* 13 Dec 2005: s6pchl: First version based on routine s6mdcl. */
/* 14 Dec 2005: s6pchl: this verson. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/*                  TDPtol = 1.1 (say) implements mild */
/*                  threshold diagonal pivoting (TDP) */
/*                  to reduce frequency of symmetric interchanges. */
/*                  TDPtol < 1.0 disables TDP. */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --perm;
    --h__;

    /* Function Body */
    *iexit = 0;
    *npos = 0;
    if (*nh == 0) {
	return 0;
    }
/* ----------------------------------------------------------------- */
/* Start the Cholesky factorization R'R = P H P'. */
/* Process the first nH rows of H. */
/* ----------------------------------------------------------------- */
    ldiag = 1;
/* Location of H(j,j) */
    incr = *maxh;
    i__1 = *nh;
    for (j = 1; j <= i__1; ++j) {
	dmax__ = h__[ldiag];
/* H(k,k), the largest diagonal */
	kmax = j;
/* kmax is k */
	lmax = ldiag;
/* Find the largest diagonal */
/* Location of H(k,k) */
	ls = ldiag + incr;
	incs = incr - 1;
	i__2 = *nh;
	for (k = j + 1; k <= i__2; ++k) {
	    if (dmax__ < h__[ls]) {
		dmax__ = h__[ls];
		kmax = k;
		lmax = ls;
	    }
	    ls += incs;
	    --incs;
	}
	supmax = 0.;
/* Search column  H(j:kmax-1,kmax) */
/* Element of largest abs in the kth row */
	lck = ldiag + kmax - j;
	inck = *maxh - j;
	i__2 = kmax - 1;
	for (i__ = j; i__ <= i__2; ++i__) {
	    if (supmax < (d__1 = h__[lck], abs(d__1))) {
		supmax = (d__1 = h__[lck], abs(d__1));
	    }
	    lck += inck;
	    --inck;
	}
/* Search row  H(kmax,kmax+1:n) */
	lrk = lmax + 1;
	i__2 = *nh;
	for (i__ = kmax + 1; i__ <= i__2; ++i__) {
	    if (supmax < (d__1 = h__[lrk], abs(d__1))) {
		supmax = (d__1 = h__[lrk], abs(d__1));
	    }
	    ++lrk;
	}
	if (dmax__ <= 0. || dmax__ < *tolpiv * supmax) {
	    goto L800;
/* Break */
	}
	dmax__ = h__[lmax];
	perm[j] = kmax;
/* --------------------------------------------------------- */
/* Perform a symmetric interchange. */
/* --------------------------------------------------------- */
	if (kmax != j) {
/* First swap row H(j,j+1:kmax) with col H(j+1:kmax,kmax). */
	    lrj = ldiag + 1;
	    lck = ldiag + kmax - j + incr - 1;
	    inck = *maxh - j - 1;
	    i__2 = kmax;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		s = h__[lrj];
		h__[lrj] = h__[lck];
		h__[lck] = s;
		++lrj;
		lck += inck;
		--inck;
	    }
/* Now swap col H(1:j,j) with col H(1:j,kmax) */
	    lcj = j;
	    lck = kmax;
	    incj = *maxh - 1;
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s = h__[lcj];
		h__[lcj] = h__[lck];
		h__[lck] = s;
		lcj += incj;
		lck += incj;
		--incj;
	    }
/* Finally swap row H(j,kmax:n) with row H(kmax,kmax:n). */
	    lrj = ldiag + kmax - j;
	    lrk = lmax;
	    i__2 = *nh;
	    for (i__ = kmax; i__ <= i__2; ++i__) {
		s = h__[lrj];
		h__[lrj] = h__[lrk];
		h__[lrk] = s;
		++lrj;
		++lrk;
	    }
	}
/* --------------------------------------------------------- */
/* Set the diagonal of  R. */
/* --------------------------------------------------------- */
/* kmax ne j */
	d__ = sqrt(dmax__);
	h__[ldiag] = d__;
	++(*npos);
	if (j < *nh) {
/* Set the super-diagonal elements of the jth row of R. */
	    jr = ldiag + 1;
	    i__2 = *nh;
	    for (k = j + 1; k <= i__2; ++k) {
		h__[jr] /= d__;
		++jr;
	    }
/* Do a rank-one update to the Schur complement. */
/* Form the upper-triangular part of H = H - x x', */
/* where x is the row H(j,j+1:nH). */
/* H(jS,:) = H(jS,:) - H(j,jS)*H(j,:). */
	    jx = ldiag + 1;
	    ls = ldiag + incr;
	    incs = incr - 1;
	    i__2 = *nh;
	    for (js = j + 1; js <= i__2; ++js) {
		hjx = h__[jx];
		if (hjx != 0.) {
		    i__ = jx;
		    k = ls;
		    i__3 = *nh;
		    for (is = js; is <= i__3; ++is) {
			h__[k] -= hjx * h__[i__];
			++i__;
			++k;
		    }
		}
		++jx;
		ls += incs;
		--incs;
	    }
	}
/* rank-one update */
	ldiag += incr;
	--incr;
    }
/* ================================================================ */
/* Test if  H  is not positive definite. */
/* ================================================================ */
L800:
    if (*npos < *nh) {
	*iexit = 1;
    }
    return 0;
} /* s6pchl_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6pchl */
/* Subroutine */ int s6rcnd_(integer *maxr, integer *ns, integer *lenr, 
	doublereal *r__, doublereal *drmax, doublereal *drmin, doublereal *
	condh)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal d__;
    static integer j, l, j1, nr;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static integer incr;
    static logical overfl;

/* ================================================================= */
/* s6Rcnd  finds the largest and smallest diagonals of the */
/* upper-triangular R, and returns the square of their ratio. */
/* This is a lower bound on the condition number of R' R. */

/* 03 Dec 2000: Converted to row-wise storage. */
/*              Rmax is NO LONGER OUTPUT. */
/* 21 Feb 2004: Current version of s6Rcnd. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;

    /* Function Body */
    overfl = FALSE_;
    nr = min(*ns,*maxr);
    if (*ns == 0) {
	*condh = 0.;
	*drmin = 1.;
	*drmax = 1.;
    } else {
	*drmax = abs(r__[1]);
	*drmin = *drmax;
	if (*ns == 1) {
/* Computing 2nd power */
	    d__2 = *drmin;
	    d__1 = d__2 * d__2;
	    *condh = ddiv_(&c_b13, &d__1, &overfl);
	    if (*condh < 1.) {
		*condh = 1. / *condh;
	    }
	} else {
	    l = 1;
	    incr = *maxr;
	    i__1 = nr;
	    for (j = 2; j <= i__1; ++j) {
		l += incr;
		--incr;
		d__ = (d__1 = r__[l], abs(d__1));
		*drmin = min(*drmin,d__);
		*drmax = max(*drmax,d__);
	    }
/* Deal with surplus diagonals of  R. */
	    if (*ns > *maxr) {
		j1 = *maxr + 1;
		l = *maxr * j1 / 2;
		i__1 = *ns;
		for (j = j1; j <= i__1; ++j) {
		    ++l;
		    d__ = (d__1 = r__[l], abs(d__1));
		    *drmin = min(*drmin,d__);
		    *drmax = max(*drmax,d__);
		}
	    }
/* Computing 2nd power */
	    d__2 = *drmax;
	    d__1 = d__2 * d__2;
/* Computing 2nd power */
	    d__4 = *drmin;
	    d__3 = d__4 * d__4;
	    *condh = ddiv_(&d__1, &d__3, &overfl);
	}
    }
    return 0;
} /* s6rcnd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rcnd */
/* Subroutine */ int s6rcol_(integer *jr, integer *maxr, integer *nr, integer 
	*lenr, doublereal *r__, doublereal *v, integer *ldiag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, l, incr;

/*     ================================================================== */
/*     s6Rcol  inserts a column of the upper-triangular matrix R. */
/*     v(1:jR) becomes column jR. */
/*     ldiag   returns the location of the new diagonal element. */

/*     03 Dec 2000: s6Radd derived from MINOS routine m6radd. */
/*     06 Dec 2000: s6Radd changed to s6Rcol because we also need s6Rrow. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --v;
    --r__;

    /* Function Body */
    l = *jr;
    incr = *maxr;
    i__1 = *jr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[l] = v[i__];
	--incr;
	l += incr;
    }
    *ldiag = l - incr;
    return 0;
} /* s6rcol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rcol */
/* Subroutine */ int s6rdel_(integer *jq, integer *maxr, integer *ns, integer 
	*lenr, doublereal *r__, doublereal *tolz)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b;
    static integer i__, j, k;
    static doublereal cs;
    static integer jr, js, lr, ls, nr;
    static doublereal rk, sk, sn, diag;
    static integer incr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer insav, nmove, lrsav;

/* ================================================================ */
/* s6Rdel  deletes the jq-th column from R. */
/* The dimension of R decreases from nS to nS - 1. */

/* 03 Dec 2000: s6Rdel derived from MINOS routine m6rdel. */
/* 17 Apr 1994: m6rdel converted to row-wise storage. */
/* 30 Jul 1994: Bottom part of R moved north-west AFTER the sweep */
/*              of rotations has eliminated the jq-th row. */
/*              This often means double handling of that part of R, */
/*              but it's the only way to skip identity rotations. */
/*              (Who knows if there are any.) */
/* 13 Jun 2001: Implemented surplus elements of R. */
/* 14 Jun 2001: Current version of s6Rdel. */
/* ================================================================ */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;

    /* Function Body */
    if (*jq == *ns) {
	return 0;
    }
/* ----------------------------------------------------------------- */
/* Delete the jq-th column of R from the top rows. */
/* For the first jq-1 rows, elements R(*,jq+1:nR) of each row */
/* are shifted 1 place to the left. */
/* ----------------------------------------------------------------- */
    nr = min(*maxr,*ns);
    lr = *jq;
    incr = *maxr;
    nmove = nr - *jq;
    i__1 = *jq - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = lr + nmove - 1;
	for (k = lr; k <= i__2; ++k) {
	    r__[k] = r__[k + 1];
	}
	--incr;
	lr += incr;
    }
/* ----------------------------------------------------------------- */
/* Triangularize the remaining rows of R, */
/* using a partial forward sweep of rotations. */

/* x x x x x x     becomes   x x x x x x */
/*   x x x x x                 x x x x x */
/*     . - - - -                 . 0 0 0 0 */
/*       x x x x                   + + + + */
/*         x x x                     + + + */
/*           x x                       + + */
/*             x                         + */
/*     |                         | */
/*    jq                        jq */

/* The . is not touched because it is later overwritten. */
/* ls marks the - being eliminated. */
/* lr marks the start of the next + + + row. */
/* ----------------------------------------------------------------- */
    lrsav = lr;
    insav = incr;
    ls = lr;
    i__1 = nr;
    for (j = *jq + 1; j <= i__1; ++j) {
	++ls;
	lr += incr;
	--incr;
	b = r__[ls];
	if (abs(b) > *tolz) {
	    a = r__[lr];
/* Computing 2nd power */
	    d__1 = a;
/* Computing 2nd power */
	    d__2 = b;
	    diag = sqrt(d__1 * d__1 + d__2 * d__2);
	    r__[lr] = diag;
	    if (j < nr) {
		cs = a / diag;
		sn = b / diag;
		jr = lr;
		js = ls;
		i__2 = nr;
		for (k = j + 1; k <= i__2; ++k) {
		    ++jr;
		    ++js;
		    rk = r__[jr];
		    sk = r__[js];
		    r__[jr] = cs * rk + sn * sk;
		    r__[js] = sn * rk - cs * sk;
		}
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Shift the + + + triangle up and left. */
/* lr marks the start of each + + + row being moved. */
/* ls marks the start of its final position. */
/* ----------------------------------------------------------------- */
    lr = lrsav;
    incr = insav;
    nmove = nr - *jq;
    i__1 = nr;
    for (j = *jq + 1; j <= i__1; ++j) {
	ls = lr;
	lr += incr;
	dcopy_(&nmove, &r__[lr], &c__1, &r__[ls], &c__1);
	--incr;
	--nmove;
    }
/* ----------------------------------------------------------------- */
/* Deal with surplus diagonals of R. */
/* ----------------------------------------------------------------- */
    if (*ns > *maxr) {
	if (*jq <= *maxr) {
/* Clear out the last column of R. */
	    lr = *maxr;
	    incr = *maxr;
	    i__1 = *maxr;
	    for (k = 1; k <= i__1; ++k) {
		r__[lr] = 0.;
		--incr;
		lr += incr;
	    }
	}
/* Shift surplus diagonals of R to the left. */
	k = max(*maxr,*jq);
	lr = *maxr * (*maxr + 1) / 2 + (k - *maxr);
	nmove = *ns - k;
	i__1 = lr + nmove - 1;
	for (k = lr; k <= i__1; ++k) {
	    r__[k] = r__[k + 1];
	}
    }
    return 0;
} /* s6rdel_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rdel */
/* Subroutine */ int s6rfix_(integer *maxr, integer *ns, integer *lenr, 
	doublereal *r__, doublereal *drmax, doublereal *drmin, doublereal *
	condh, doublereal *hcndbd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer j, l, j1, nr, incr;

/* ================================================================ */
/* s6Rfix  is called after s6Rcnd if R is ill-conditioned. */
/* It increases the magnitude of small diagonals. */
/* dRmax, dRmin, Hcndbd are input. */
/*        dRmin, condH  are output. */

/* 29 Jul 2003: First version of s6Rfix. */
/* 29 Jul 2003: Current version of s6Rfix. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;

    /* Function Body */
    if (*ns == 0) {
	return 0;
    }
    nr = min(*ns,*maxr);
/* Computing MAX */
    d__1 = *drmin, d__2 = *drmax / *hcndbd;
    *drmin = max(d__1,d__2);
    if (*ns == 1) {
	r__[1] = 1.;
/* Reset! */
	*condh = 1.;
    } else {
	l = 1;
	incr = *maxr;
	i__1 = nr;
	for (j = 2; j <= i__1; ++j) {
	    l += incr;
	    --incr;
	    d__ = r__[l];
	    if (abs(d__) < *drmin) {
		if (d__ >= 0.) {
		    d__ = *drmin;
		} else {
		    d__ = -(*drmin);
		}
		r__[l] = d__;
	    }
	}
/* Deal with surplus diagonals of  R. */
	if (*ns > *maxr) {
	    j1 = *maxr + 1;
	    l = *maxr * j1 / 2;
	    i__1 = *ns;
	    for (j = j1; j <= i__1; ++j) {
		++l;
		d__ = r__[l];
		if (abs(d__) < *drmin) {
		    if (d__ >= 0.) {
			d__ = *drmin;
		    } else {
			d__ = -(*drmin);
		    }
		    r__[l] = d__;
		}
	    }
	}
	*condh = *hcndbd;
    }
    return 0;
} /* s6rfix_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rfix */
/* Subroutine */ int s6rmod_(integer *iexit, integer *maxr, integer *ns, 
	integer *lenr, doublereal *r__, doublereal *u, doublereal *v, integer 
	*lastnz, doublereal *ulast, doublereal *told, doublereal *tolz)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l;
    static doublereal s, t, u2, cs, sn;
    static integer lm1, incr;
    static doublereal root;
    static integer ldiag;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nmove, lastr;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ================================================================= */
/* s6Rmod  modifies the (nS+1) x nS upper-triangular matrix R so */
/* that Q*(R + u v') is upper triangular,  where  Q  is orthogonal. */
/* The arrays v and u hold v and the first nS elements of u */
/* respectively.  The (nS+1)th component of u is held in ulast. */
/* The new R overwrites the old. */

/* Q is the product of two sweeps of plane rotations (not stored). */
/* These affect the (lastnz)th row of R, which is temporarily held */
/* in the array u.  Thus, u is overwritten.  v is not altered. */

/* ulast  holds  u(nS+1).   It is overwritten. */

/* lastnz points to the last nonzero of the vector ( u, ulast ). */
/*        The value lastnz = nS+1 is always safe, but sometimes */
/*        a value of lastnz less than nS+1 is known, in which case */
/*        Q reduces to two partial sweeps. */

/* told   is a tolerance on the lastv-th diagonal of R. */
/* tolz   is a tolerance for negligible elements in  u. */

/* On exit, */
/* iExit = 1  if the diagonal of R is larger than told, */
/*       = 2  if not (the diagonal is not modified). */

/* 06 Sep 1991: First version based on Minos routine m6rmod. */
/* 11 Sep 1994: Modified to update a principal submatrix of R. */
/* 03 Dec 2000: Converted to row-wise storage. */
/* 17 Jun 2001: Current version of s6Rmod. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --v;
    --u;
    --r__;

    /* Function Body */
    if (*lastnz <= *ns) {
	*ulast = u[*lastnz];
    }
/* Copy the (lastnz)th row of R into the end of u. */
    lm1 = *lastnz - 1;
    lastr = lm1 * *maxr + (3 - *lastnz) * *lastnz / 2;
    nmove = *ns - lm1;
    if (nmove > 0) {
	dcopy_(&nmove, &r__[lastr], &c__1, &u[*lastnz], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Reduce u to a multiple of e(lastnz) using a partial backward swe */
/* of rotations.  This fills in the (lastnz)th row of R (held in u) */
/* ----------------------------------------------------------------- */
    if (*lastnz > 1) {
/* Computing 2nd power */
	d__1 = *ulast;
	u2 = d__1 * d__1;
	ldiag = lastr;
	incr = *maxr - lm1;
	for (i__ = lm1; i__ >= 1; --i__) {
	    ++incr;
	    ldiag -= incr;
	    s = u[i__];
	    u[i__] = 0.;
	    if (abs(s) > *tolz) {
/* Computing 2nd power */
		d__1 = s;
		u2 = d__1 * d__1 + u2;
		root = sqrt(u2);
		cs = *ulast / root;
		sn = s / root;
		*ulast = root;
		l = ldiag;
		i__1 = *ns;
		for (j = i__; j <= i__1; ++j) {
		    s = u[j];
		    t = r__[l];
		    u[j] = cs * s + sn * t;
		    r__[l] = sn * s - cs * t;
		    ++l;
		}
	    }
	}
    }
    daxpy_(ns, ulast, &v[1], &c__1, &u[1], &c__1);
/* ----------------------------------------------------------------- */
/* Eliminate the front of the (lastnz)th row of R (held in u) using */
/* partial forward sweep of rotations. */
/* ----------------------------------------------------------------- */
/* Set u = u  +  ulast*v */
    if (*lastnz > 1) {
	ldiag = 1;
	incr = *maxr;
	i__1 = lm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = u[i__];
	    if (abs(t) > *tolz) {
		s = r__[ldiag];
/* Computing 2nd power */
		d__1 = s;
/* Computing 2nd power */
		d__2 = t;
		root = sqrt(d__1 * d__1 + d__2 * d__2);
		cs = s / root;
		sn = t / root;
		r__[ldiag] = root;
		l = ldiag;
		i__2 = *ns;
		for (j = i__ + 1; j <= i__2; ++j) {
		    ++l;
		    s = r__[l];
		    t = u[j];
		    r__[l] = cs * s + sn * t;
		    u[j] = sn * s - cs * t;
		}
	    }
	    ldiag += incr;
	    --incr;
	}
    }
/* Insert the new (lastnz)th row of  R. */
    if (nmove > 0) {
	dcopy_(&nmove, &u[*lastnz], &c__1, &r__[lastr], &c__1);
/* -------------------------------------------------------------- */
/* Test for (unlikely) singularity. */
/* -------------------------------------------------------------- */
	*iexit = 1;
	if ((d__1 = r__[lastr], abs(d__1)) <= *told) {
	    *iexit = 2;
/* r(lastr) = told */
	}
    }
    return 0;
} /* s6rmod_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rmod */
/* Subroutine */ int s6rprd_(integer *task, integer *maxr, integer *ns, 
	integer *lenr, doublereal *r__, doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, l, nr;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer incr, numr;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);

/* ================================================================= */
/* s6Rprd  forms y = Rx or y = R'x, where R is an upper-triangular */
/* matrix of dimension nS stored by rows in R(lenR). */

/* 03 Dec 2000: Converted to row-wise storage. */
/* 12 Jun 2001: Surplus diagonals allowed in R. */
/* 31 Jul 2003: Allow maxR = 0. */
/* 20 Jun 2004: Current version of s6Rprd. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --y;
    --x;
    --r__;

    /* Function Body */
    if (*maxr > 0) {
	nr = min(*maxr,*ns);
	if (*task == 0) {
/* Form y = R x. */
	    numr = nr;
	    l = 1;
	    incr = *maxr;
	    i__1 = nr;
	    for (j = 1; j <= i__1; ++j) {
		y[j] = ddot_(&numr, &r__[l], &c__1, &x[j], &c__1);
		l += incr;
		--numr;
		--incr;
	    }
	} else if (*task == 1) {
/* Form y = Rtranspose x. */
	    dload_(&nr, &c_b29, &y[1], &c__1);
	    l = 1;
	    incr = *maxr;
	    numr = nr;
	    i__1 = nr;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		daxpy_(&numr, &x[i__], &r__[l], &c__1, &y[i__], &c__1);
		l += incr;
		--incr;
		--numr;
	    }
	}
    }
/* Deal with surplus diagonals of R. */
    if (*ns > *maxr) {
	l = *maxr * (*maxr + 1) / 2;
	i__1 = *ns;
	for (j = *maxr + 1; j <= i__1; ++j) {
	    ++l;
	    y[j] = r__[l] * x[j];
	}
    }
    return 0;
} /* s6rprd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rprd */
/* Subroutine */ int s6rrow_(integer *jr, integer *maxr, integer *nr, integer 
	*lenr, doublereal *r__, doublereal *v, integer *ldiag)
{
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nmove;

/* ================================================================ */
/* s6Rrow     inserts a row of the upper-triangular matrix R. */
/* v(jR+1:nR) becomes row jR. */
/* ldiag      returns the location of the new diagonal element. */

/* 06 Dec 2000: First version of s6Rrow, called from s5HZ. */
/* ================================================================ */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --v;
    --r__;

    /* Function Body */
    *ldiag = (*jr - 1) * *maxr + (3 - *jr) * *jr / 2;
/* Magic formula! */
    nmove = *nr - *jr + 1;
    dcopy_(&nmove, &v[*jr], &c__1, &r__[*ldiag], &c__1);
    return 0;
} /* s6rrow_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rrow */
/* Subroutine */ int s6rset_(logical *gotr, integer *maxr, integer *ns, 
	integer *lenr, doublereal *r__, doublereal *w, doublereal *condr)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, l, nr;
    static doublereal diag, dmin__, dmax__;
    static integer incr, ldiag;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer nzero;

/* ================================================================= */
/* s6Rset  alters R, the upper-triangular factor of the approximate */
/* reduced Hessian. */

/* If  gotR = .false., R does not exist. */
/* In this case, R is initialized to the identity matrix. */

/* Otherwise, R already exists and we attempt to make it better */
/* conditioned by scaling its columns by the square roots of its */
/* current diagonals. */

/* 12 Jun 2001: First version based on MINOS routine m6rset. */
/* 21 Feb 2004: Current version of s6Rset. */
/* ================================================================ */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --w;
    --r__;

    /* Function Body */
    nr = min(*maxr,*ns);
    if (nr == 0) {
	*condr = 0.;
	*gotr = TRUE_;
    } else if (*gotr) {
/* -------------------------------------------------------------- */
/* Scale the columns of R. */
/* -------------------------------------------------------------- */
/* Find dmin and dmax, and set w = set of scale factors. */
	dmax__ = abs(r__[1]);
	dmin__ = dmax__;
	ldiag = 1;
	incr = *maxr;
	i__1 = nr;
	for (k = 1; k <= i__1; ++k) {
	    diag = (d__1 = r__[ldiag], abs(d__1));
	    dmax__ = max(dmax__,diag);
	    dmin__ = min(dmin__,diag);
	    w[k] = 1. / sqrt(diag);
	    ldiag += incr;
	    --incr;
	}
/* Apply column scales to each row of R. */
	ldiag = 1;
	incr = *maxr;
	i__1 = nr;
	for (k = 1; k <= i__1; ++k) {
	    l = ldiag;
	    i__2 = nr;
	    for (j = k; j <= i__2; ++j) {
		r__[l] *= w[j];
		++l;
	    }
	    ldiag += incr;
	    --incr;
	}
	*condr = dmax__ / dmin__;
    } else {
/* Set R = the identity. */
	ldiag = 1;
	incr = *maxr;
	nzero = nr - 1;
	i__1 = nr - 1;
	for (k = 1; k <= i__1; ++k) {
	    r__[ldiag] = 1.;
	    dload_(&nzero, &c_b29, &r__[ldiag + 1], &c__1);
	    ldiag += incr;
	    --incr;
	    --nzero;
	}
	r__[ldiag] = 1.;
	i__1 = *ns;
	for (k = *maxr + 1; k <= i__1; ++k) {
	    ++ldiag;
	    r__[ldiag] = 1.;
	}
	*condr = 1.;
	*gotr = TRUE_;
    }
    return 0;
} /* s6rset_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rset */
/* Subroutine */ int s6rsol_(integer *task, integer *maxr, integer *ns, 
	integer *lenr, doublereal *r__, doublereal *y)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l;
    static doublereal t;
    static integer nr;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer incr, numr;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ================================================================ */
/* s6Rsol  solves Rx = y or R'x = y, where R is an */
/* upper-triangular matrix of dimension nS stored by rows in R(lenR */
/* The solution x overwrites y. */

/* 03 Dec 2000: Converted to row-wise storage. */
/* 12 Jun 2001: Surplus diagonals allowed in R. */
/* 31 Jul 2003: Allow maxR = 0. */
/* 31 Jul 2003: Current version of s6Rsol. */
/* ================================================================ */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --y;
    --r__;

    /* Function Body */
    if (*maxr > 0) {
	nr = min(*maxr,*ns);
	if (*task == 0) {
/* Solve R y = y. */
	    l = (nr - 1) * *maxr + (3 - nr) * nr / 2;
	    y[nr] /= r__[l];
	    incr = *maxr + 1 - nr;
	    numr = 0;
	    for (j = nr - 1; j >= 1; --j) {
		++numr;
		++incr;
		l -= incr;
		t = ddot_(&numr, &r__[l + 1], &c__1, &y[j + 1], &c__1);
		y[j] = (y[j] - t) / r__[l];
	    }
	} else if (*task == 1) {
/* Solve Rtranspose y = y. */
	    l = 1;
	    incr = *maxr;
	    numr = nr - 1;
	    i__1 = nr - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] /= r__[l];
		d__1 = -y[i__];
		daxpy_(&numr, &d__1, &r__[l + 1], &c__1, &y[i__ + 1], &c__1);
		l += incr;
		--incr;
		--numr;
	    }
	    y[nr] /= r__[l];
	}
    }
/* Deal with surplus diagonals of R. */
    if (*ns > *maxr) {
	l = *maxr * (*maxr + 1) / 2;
	i__1 = *ns;
	for (j = *maxr + 1; j <= i__1; ++j) {
	    ++l;
	    y[j] /= r__[l];
	}
    }
    return 0;
} /* s6rsol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rsol */
/* Subroutine */ int s6rswp_(integer *maxr, integer *nr, integer *lenr, 
	doublereal *r__, doublereal *v, doublereal *w, integer *lastnz, 
	doublereal *tolz)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, l, incr;
    static doublereal told, vlast, vnorm;
    extern /* Subroutine */ int s6rmod_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__;

/* ================================================================ */
/* s6Rswp  modifies the upper-triangular matrix R to account for a */
/* basis exchange in which the (lastnz)th superbasic becomes basic. */
/* R is changed to R + vw', which is triangularized by  s6Rmod. */
/* v is the (lastnz)th column of R, and w is input. */

/* 01 Dec 1991: First version based on Minos routine m6bswp. */
/* 24 Jan 1996: v left unscaled. */
/* 03 Dec 2000: Converted to row-wise storage. */
/* 17 Jun 2001: Current version of s6Rswp. */
/* ================================================================ */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --w;
    --v;
    --r__;

    /* Function Body */
    told = 0.;
/* Set  v =  (lastnz)-th column of  R  and find its norm. */
/* Singularity checked elsewhere */
    l = *lastnz;
    incr = *maxr - 1;
    i__1 = *lastnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__] = r__[l];
	l += incr;
	--incr;
    }
    vnorm = dnormi_(lastnz, &v[1], &c__1);
    vlast = 0.;
/* v(nR+1) = 0 */
    d__1 = vnorm * *tolz;
    s6rmod_(&inform__, maxr, nr, lenr, &r__[1], &v[1], &w[1], lastnz, &vlast, 
	    &told, &d__1);
    return 0;
} /* s6rswp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rswp */
/* Subroutine */ int s6rupd_(integer *update, integer *maxr, integer *lenr, 
	integer *m, integer *n, integer *nbs, integer *nnl, integer *ns, 
	doublereal *u0scal, doublereal *rdxhdx, integer *ne, integer *nloca, 
	integer *loca, integer *inda, doublereal *acol, integer *kbs, 
	doublereal *v, doublereal *hdx, doublereal *r__, doublereal *w, 
	doublereal *y, doublereal *y1, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l, nr;
    static doublereal eps0;
    static integer incr;
    static doublereal told;
    static integer numr;
    static doublereal tolz;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), s2bprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), s2bsol_(integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s6rbfs_(integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *), s6rsol_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static integer inform__;
    extern /* Subroutine */ int s2gathr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);

/* ================================================================= */
/* s6Rupd computes the effect of the full BFGS update on the */
/* Cholesky factor of the reduced Hessian. */


/* On entry: */
/* R   contains the first nS rows and columns of an nnL x nnL */
/*     upper-triangular matrix R such that R'R = Q'HQ, Q = ( Z Y ). */

/* v   contains the nnL components of the BFGS update */
/*     U(new) = U(I + sv'), with H = U'U. */
/*     s  = dx / rdxHdx, v = (1/rydx) gdif - (1/rdxHdx) Hdx. */

/* Hdx contains the product H dx. */

/* w, y, y1 are work vectors. */

/* 04 Aug 1995: First version based on NPSOL routine npupdt. */
/* 03 Dec 2000: Converted to row-wise storage. */
/* 18 Feb 2001: H stored in product form. */
/* 17 Nov 2001: Current version of s6Rupd. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --r__;
    --y1;
    --y;
    --w;
    --kbs;
    --hdx;
    --v;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
    told = 0.;
    tolz = 0.;
    nr = min(*maxr,*ns);
/* ----------------------------------------------------------------- */
/* Find  y (B) =   v(B),  y (S) =   v(S) - S'Binv v(B)   = Z'v. */
/* Find  y1(B) = Hdx(B),  y1(S) = Hdx(S) - S'Binv Hdx(B) = Z'Hdx. */
/* ----------------------------------------------------------------- */
    s2gathr_(nnl, nbs, &kbs[1], &c_b13, &v[1], &y[1]);
    s2bsol_(&inform__, &c__2, m, &y[1], &w[1], &iw[1], leniw, &rw[1], lenrw);
    if (*ns > 0) {
	s2bprd_(&c__1, &eps0, n, ns, &kbs[*m + 1], ne, nloca, &loca[1], &inda[
		1], &acol[1], &c_b50, &w[1], m, &c_b13, &y[*m + 1], ns);
    }
    s2gathr_(nnl, nbs, &kbs[1], &c_b13, &hdx[1], &y1[1]);
    s2bsol_(&inform__, &c__2, m, &y1[1], &w[1], &iw[1], leniw, &rw[1], lenrw);
    if (*ns > 0) {
	s2bprd_(&c__1, &eps0, n, ns, &kbs[*m + 1], ne, nloca, &loca[1], &inda[
		1], &acol[1], &c_b50, &w[1], m, &c_b13, &y1[*m + 1], ns);
    }
    if (*update > 1) {
/* -------------------- */
/* Self-scaled BFGS. */
/* -------------------- */
	l = 1;
	incr = *maxr;
	numr = nr;
	i__1 = nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dscal_(&numr, u0scal, &r__[l], &c__1);
	    l += incr;
	    --incr;
	    --numr;
	}
/* Deal with surplus diagonals of R. */
	if (*ns > *maxr) {
	    l = *maxr * (*maxr + 1) / 2;
	    i__1 = *ns;
	    for (j = *maxr + 1; j <= i__1; ++j) {
		++l;
		r__[l] = *u0scal * r__[l];
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Apply the BFGS update to the first nS columns of R. */
/* y(S) = Z'v  and  y1(S) = u/norm(Rdx). */
/* ----------------------------------------------------------------- */
    inform__ = 0;
    dcopy_(ns, &y1[*m + 1], &c__1, &w[*m + 1], &c__1);
    s6rsol_(&c__1, maxr, ns, lenr, &r__[1], &y1[*m + 1]);
    d__1 = 1. / *rdxhdx;
    dscal_(ns, &d__1, &y1[*m + 1], &c__1);
    d__1 = 1. / *rdxhdx;
    s6rbfs_(&inform__, maxr, ns, nnl, lenr, &told, &tolz, &r__[1], &y1[*m + 1]
	    , &y[*m + 1], &d__1, &w[*m + 1]);
    return 0;
} /* s6rupd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6Rupd */
/* Subroutine */ int s6rqn_(integer *iexit, logical *exactp, integer *maxr, 
	integer *ns, integer *lenr, doublereal *r__, doublereal *gtp, 
	doublereal *g, doublereal *g2, doublereal *p, doublereal *w, 
	doublereal *step, doublereal *told, doublereal *tolz)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer j;
    static doublereal pj, gtp2;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal delta1, delta2;
    extern /* Subroutine */ int s6rbfs_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s6rprd_(
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);

/* ================================================================= */
/* s6Rqn  applies the BFGS update to the upper-triangular matrix R, */
/* which holds the factor of the quasi-Newton approximation of the */
/* reduced Hessian. */

/* R  contains a triangle of size  nR = min( nS, maxR ). */
/* If  nS .gt. maxR,  R also contains a diagonal of size nS - maxR. */

/* g, g2   hold the gradients.                  g is overwritten. */
/* p       holds the search direction.         It is overwritten. */
/* w       must satisfy  w = -Rp.              It is overwritten. */
/* Exactp  = true, means that  w  also satisfies   R' w = g. */
/*         This holds if p satisfies R'Rp = -g. */

/* On exit, */
/* iExit = 0  if no update was performed, */
/*       = 1  if the update was successful, */
/*       = 2  if it was nearly singular. */

/* 12 Jun 2001: First version based on MINOS routine m6bfgs. */
/* 14 Apr 2004: Added option allowing arbitrary directions. */
/* 15 Apr 2004: Current version of s6Rqn. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --w;
    --p;
    --g2;
    --g;
    --r__;

    /* Function Body */
    *iexit = 0;
    gtp2 = ddot_(ns, &g2[1], &c__1, &p[1], &c__1);
    if (gtp2 <= *gtp * .91) {
	return 0;
    }
    delta2 = 1. / sqrt(*step * (gtp2 - *gtp));
/* 1/sqrt(y's) */
    if (*exactp) {
/* -------------------------------------------------------------- */
/* The vector w satisfies R'w = g. */
/* -------------------------------------------------------------- */
	delta1 = 1. / sqrt((abs(*gtp)));
/* Normalize  w  and change its sign. Now w = Rp/norm(Rp) */
/* 1/sqrt(p'Bp) */
	d__1 = -delta1;
	dscal_(ns, &d__1, &w[1], &c__1);
/* Avoid cancellation error in forming the new vector  p. */
	if ((d__1 = delta1 / delta2 - 1., abs(d__1)) >= .5) {
	    i__1 = *ns;
	    for (j = 1; j <= i__1; ++j) {
		p[j] = delta2 * (g2[j] - g[j]) + delta1 * g[j];
	    }
	} else {
	    d__ = delta1 - delta2;
	    i__1 = *ns;
	    for (j = 1; j <= i__1; ++j) {
		p[j] = delta2 * g2[j] + d__ * g[j];
	    }
	}
    } else {
/* -------------------------------------------------------------- */
/* The product  R'w = -R'Rp = -Bp is stored in g. */
/* -------------------------------------------------------------- */
	delta1 = 1. / dnrm2_(ns, &w[1], &c__1);
/* Form   R'*w = -R'Rp in p */
/* 1/sqrt(p'Bp) */
	s6rprd_(&c__1, maxr, ns, lenr, &r__[1], &w[1], &p[1]);
/* Normalize  w  and change its sign. Now w = Rp/norm(Rp) */
	d__1 = -delta1;
	dscal_(ns, &d__1, &w[1], &c__1);
	i__1 = *ns;
	for (j = 1; j <= i__1; ++j) {
	    pj = p[j];
	    p[j] = delta2 * (g2[j] - g[j]) + delta1 * pj;
	    g[j] = pj;
	}
    }
/* Triangularize   R  +  w p'. */
/* Deal with surplus diagonals of  R. */
    d__1 = -delta1;
    s6rbfs_(iexit, maxr, ns, ns, lenr, told, tolz, &r__[1], &w[1], &p[1], &
	    d__1, &g[1]);
    return 0;
} /* s6rqn_ */

