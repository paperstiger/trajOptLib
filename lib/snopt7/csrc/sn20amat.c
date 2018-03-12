/* ../snopt7/src/sn20amat.f -- translated by f2c (version 20100827).
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

static integer c__11 = 11;
static integer c__1 = 1;
static integer c__12 = 12;
static integer c__2 = 2;
static doublereal c_b94 = 0.;
static integer c__31 = 31;
static integer c__0 = 0;
static doublereal c_b109 = 1.;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__21 = 21;
static integer c__14 = 14;
static integer c__13 = 13;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn20amat.f */

/*     s2Amat   s2Aprd   s2Bprd   s2bInf   s2ColN   s2RowN   s2VarN */
/*     s2dInf   s2gathr  s2scatr  s2crsh   s2Mem0   s2Mem    s2rcA */
/*     s2Scal   s2SclA   s2unpk   s2vmax   s2dmat   s2xmat */

/*     27 Dec 2013: For Scale option 1, s2Scal sets jStart = nnJac + 1, not nnL + 1. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s2amat_(integer *task, integer *lprint, integer *m, 
	integer *n, integer *nb, integer *nncon, integer *nnjac, integer *
	nnobj, integer *iobj, integer *numlc, integer *numliq, integer *ne, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	doublereal *bl, doublereal *bu, integer *hrtype, integer *iw, integer 
	*leniw, doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_2300[] = "(a,5i12)";
    static char fmt_2400[] = "(\002 No. of matrix elements\002,i21,5x,\002De"
	    "nsity\002,f12.3)";
    static char fmt_2410[] = "(\002 Biggest \002,1p,e35.4,\002  (excluding f"
	    "ixed columns,\002)";
    static char fmt_2420[] = "(\002 Smallest\002,1p,e35.4,\002   free rows, "
	    "and RHS)\002)";
    static char fmt_2411[] = "(\002 Biggest  constant element\002,1p,e18.4"
	    ",\002  (excluding fixed columns,\002)";
    static char fmt_2421[] = "(\002 Smallest constant element\002,1p,e18.4"
	    ",\002   free rows, and RHS)\002)";
    static char fmt_2430[] = "(\002 No. of objective coefficients\002,i14)";
    static char fmt_2450[] = "(\002 Biggest \002,1p,e35.4,\002  (excluding f"
	    "ixed columns)\002)";
    static char fmt_2460[] = "(\002 Smallest\002,1p,e35.4)";
    static char fmt_2500[] = "(\002 Nonlinear constraints\002,i8,5x,\002Line"
	    "ar constraints\002,i8)";
    static char fmt_2510[] = "(\002 Nonlinear variables  \002,i8,5x,\002Line"
	    "ar variables  \002,i8)";
    static char fmt_2520[] = "(\002 Jacobian  variables  \002,i8,5x,\002Obje"
	    "ctive variables\002,i7)";
    static char fmt_2530[] = "(\002 Total constraints    \002,i8,5x,\002Tota"
	    "l variables   \002,i8)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, l;
    static doublereal b1, b2, aij;
    static integer nnl;
    static char str[80];
    static integer free;
    static doublereal cmin, cmax;
    static integer lcon, lvar, cost;
    static logical summ1, prnt1;
    static integer bnded, fixed;
    static doublereal bplus;
    static integer norml;
    static doublereal infbnd, aijmin, aijmax, bminus;
    static integer numleq;
    static doublereal adnsty;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___17 = { 0, str, 0, fmt_2300, 80, 1 };
    static icilist io___18 = { 0, str, 0, fmt_2300, 80, 1 };
    static icilist io___27 = { 0, str, 0, fmt_2400, 80, 1 };
    static icilist io___28 = { 0, str, 0, fmt_2410, 80, 1 };
    static icilist io___29 = { 0, str, 0, fmt_2420, 80, 1 };
    static icilist io___30 = { 0, str, 0, fmt_2411, 80, 1 };
    static icilist io___31 = { 0, str, 0, fmt_2421, 80, 1 };
    static icilist io___32 = { 0, str, 0, fmt_2430, 80, 1 };
    static icilist io___33 = { 0, str, 0, fmt_2450, 80, 1 };
    static icilist io___34 = { 0, str, 0, fmt_2460, 80, 1 };
    static icilist io___37 = { 0, str, 0, fmt_2500, 80, 1 };
    static icilist io___38 = { 0, str, 0, fmt_2510, 80, 1 };
    static icilist io___39 = { 0, str, 0, fmt_2520, 80, 1 };
    static icilist io___40 = { 0, str, 0, fmt_2530, 80, 1 };


/*     ================================================================== */
/*     s2Amat defines hrtype, the set of row types. */

/*     If Task = Rowtyp (= 0), only the row types are computed. */
/*     If Task = Stats  (= 1), matrix statistics are also computed. */

/*     The vector of row-types is as follows: */
/*        hrtype(i) = 0  for E rows      (equalities) */
/*        hrtype(i) = 1  for L or G rows (inequalities) */
/*        hrtype(i) = 2  for N rows      (free rows) */
/*     They are used in s2Scal and s2crsh. */

/*     15 Feb 1991: First version based on Minos 5.4 routine m2Amat. */
/*     04 Apr 1999: Objective stored in A. */
/*     31 Jul 2003: snPRNT adopted. */
/*     14 Jun 2008: numLC and numLIQ  computed here instead of s5getB. */
/*                  Constant elements used for the matrix statistics. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hrtype;
    --bu;
    --bl;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    nnl = max(*nnobj,*nnjac);
    bplus = infbnd * .9;
    bminus = -bplus;
    prnt1 = *lprint >= 1;
    summ1 = *lprint >= 1;
    if (*task == 0 || *task == 1) {
/* Construct the vector of row-types. */
	fixed = 0;
	free = 0;
	norml = 0;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = *n + i__;
	    b1 = bl[j];
	    b2 = bu[j];
	    if (b1 == b2) {
		hrtype[i__] = 0;
		++fixed;
	    } else if (b1 <= bminus && b2 >= bplus) {
		hrtype[i__] = 2;
		++free;
	    } else {
		hrtype[i__] = 1;
		if (b1 <= bminus || b2 >= bplus) {
		    ++norml;
		}
	    }
	}
    }
    if (*task == 1) {
/* Count the linear equality and inequality constraints */
/* (E and LG rows). */
	numleq = 0;
	i__1 = *nb;
	for (j = *n + *nncon + 1; j <= i__1; ++j) {
	    if (bl[j] == bu[j]) {
		++numleq;
	    }
	}
	*numliq = *numlc - numleq;
	if (prnt1) {
	    bnded = *m - fixed - free - norml;
	    snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
	    snprnt_(&c__11, " Matrix statistics", &iw[1], leniw, (ftnlen)18);
	    snprnt_(&c__1, " -----------------", &iw[1], leniw, (ftnlen)18);
	    snprnt_(&c__1, "               Total      Normal        Free    "
		    "   Fixed     Bounded", &iw[1], leniw, (ftnlen)68);
	    s_wsfi(&io___17);
	    do_fio(&c__1, " Rows   ", (ftnlen)8);
	    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&norml, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&free, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fixed, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bnded, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    fixed = 0;
	    free = 0;
	    norml = 0;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		b1 = bl[j];
		b2 = bu[j];
		if (b1 == b2) {
		    ++fixed;
		} else {
		    if (b1 == 0.) {
			if (b2 >= bplus) {
			    ++norml;
			}
		    } else if (b1 <= bminus) {
			if (b2 == 0.) {
			    ++norml;
			} else if (b2 >= bplus) {
			    ++free;
			}
		    }
		}
	    }
	    bnded = *n - fixed - free - norml;
	    s_wsfi(&io___18);
	    do_fio(&c__1, " Columns", (ftnlen)8);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&norml, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&free, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fixed, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bnded, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
/* Find the biggest and smallest constant elements in Acol, */
/* excluding free rows and fixed columns.  Also find the */
/* largest objective coefficient in Acol. */
	    aijmax = 0.;
	    aijmin = bplus;
	    cmax = 0.;
	    cmin = bplus;
	    cost = 0;
	    i__1 = nnl;
	    for (j = 1; j <= i__1; ++j) {
		if (bl[j] < bu[j]) {
		    i__2 = loca[j + 1] - 1;
		    for (l = loca[j]; l <= i__2; ++l) {
			i__ = inda[l];
			if (i__ > *nncon) {
			    if (hrtype[i__] == 2) {
				if (i__ == *iobj) {
				    aij = (d__1 = acol[l], abs(d__1));
				    if (aij > 0.) {
					++cost;
					cmax = max(cmax,aij);
					cmin = min(cmin,aij);
				    }
				}
			    } else {
				aij = (d__1 = acol[l], abs(d__1));
				aijmax = max(aijmax,aij);
				aijmin = min(aijmin,aij);
			    }
			}
		    }
		}
	    }
	    if (aijmin == bplus) {
		aijmin = 0.;
	    }
	    adnsty = *ne * 100. / (*m * *n);
	    s_wsfi(&io___27);
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&adnsty, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    if (*nncon == 0) {
		s_wsfi(&io___28);
		do_fio(&c__1, (char *)&aijmax, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
		s_wsfi(&io___29);
		do_fio(&c__1, (char *)&aijmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    } else {
		s_wsfi(&io___30);
		do_fio(&c__1, (char *)&aijmax, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
		s_wsfi(&io___31);
		do_fio(&c__1, (char *)&aijmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	    s_wsfi(&io___32);
	    do_fio(&c__1, (char *)&cost, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    if (cost > 0) {
		s_wsfi(&io___33);
		do_fio(&c__1, (char *)&cmax, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
		s_wsfi(&io___34);
		do_fio(&c__1, (char *)&cmin, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
/* Print a few things that may be gathered as statistics */
/* from a bunch of test runs. */
	if (prnt1 || summ1) {
	    lvar = *n - nnl;
	    lcon = *m - *nncon;
	    s_wsfi(&io___37);
	    do_fio(&c__1, (char *)&(*nncon), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lcon, (ftnlen)sizeof(integer));
	    e_wsfi();
	    if (prnt1) {
		snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    }
	    if (summ1) {
		snprnt_(&c__12, str, &iw[1], leniw, (ftnlen)80);
	    }
	    s_wsfi(&io___38);
	    do_fio(&c__1, (char *)&nnl, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lvar, (ftnlen)sizeof(integer));
	    e_wsfi();
	    if (prnt1) {
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	    if (summ1) {
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)80);
	    }
	    s_wsfi(&io___39);
	    do_fio(&c__1, (char *)&(*nnjac), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nnobj), (ftnlen)sizeof(integer));
	    e_wsfi();
	    if (prnt1) {
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	    if (summ1) {
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)80);
	    }
	    s_wsfi(&io___40);
	    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    e_wsfi();
	    if (prnt1) {
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	    if (summ1) {
		snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
    return 0;
} /* s2amat_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Amat */
/* Subroutine */ int s2aprd_(integer *task, doublereal *tolz, integer *ne, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	doublereal *alpha, doublereal *x, integer *lenx, doublereal *beta, 
	doublereal *y, integer *leny)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal sum, alphxj;

/*     ================================================================== */
/*     s2Aprd computes matrix-vector products involving  x  and a sparse */
/*     matrix A  stored by columns. */
/*     The parameter Task specifies the operation to be done as follows: */
/*       Task = Normal (=0)    y := alpha*A *x + beta*y, */
/*       Task = Transp (=1)    y := alpha*A'*x + beta*y, */
/*     where alpha and beta are scalars, x and y are vectors and A is a */
/*     sparse matrix stored by columns. */

/*     28 Jul 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --acol;
    --inda;
    --loca;
    --x;
    --y;

    /* Function Body */
    if (*alpha == 0. && *beta == 1.) {
	return 0;
    }
/*     First form  y := beta*y. */
    if (*beta != 1.) {
	if (*beta == 0.) {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = 0.;
	    }
	} else {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = *beta * y[i__];
	    }
	}
    }
    if (*alpha == 0.) {
/*        Relax */
    } else if (*task == 0) {
	i__1 = *lenx;
	for (j = 1; j <= i__1; ++j) {
	    alphxj = *alpha * x[j];
	    if (abs(alphxj) > *tolz) {
		i__2 = loca[j + 1] - 1;
		for (l = loca[j]; l <= i__2; ++l) {
		    i__ = inda[l];
		    if (i__ <= *leny) {
			y[i__] += acol[l] * alphxj;
		    }
		}
	    }
	}
    } else if (*task == 1) {
	i__1 = *leny;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
	    i__2 = loca[j + 1] - 1;
	    for (l = loca[j]; l <= i__2; ++l) {
		i__ = inda[l];
		if (i__ <= *lenx) {
		    sum += acol[l] * x[i__];
		}
	    }
	    y[j] += *alpha * sum;
	}
    }
    return 0;
} /* s2aprd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Aprd */
/* Subroutine */ int s2bprd_(integer *task, doublereal *tolz, integer *n, 
	integer *lenkbs, integer *kbs, integer *ne, integer *nloca, integer *
	loca, integer *inda, doublereal *acol, doublereal *alpha, doublereal *
	x, integer *lenx, doublereal *beta, doublereal *y, integer *leny)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t, alphxj;

/*     ================================================================== */
/*     s2Bprd computes various matrix-vector products involving */
/*     B  and  S,  the basic and superbasic columns of  A. The variable */
/*     Task specifies the operation to be performed as follows: */
/*         Task = Normal (= 0)         y := alpha*A *x + beta*y, */
/*         Task = Transp (= 1)         y := alpha*A'*x + beta*y, */
/*     where alpha and beta are scalars, x and y are vectors, and A is a */
/*     sparse matrix whose columns are in the basic-superbasic order */
/*     A = ( B  S ). */

/*     23 Nov 1991: First version of s2Bprd. */
/*     17 Jul 1996: Standard implementation for slacks. */
/*     22 Mar 1999: Reverted to integer Task control. */
/*     01 Apr 1999: Acol includes the objective as row 0. */
/*     28 Jul 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --kbs;
    --acol;
    --inda;
    --loca;
    --x;
    --y;

    /* Function Body */
    if (*alpha == 0. && *beta == 1.) {
	return 0;
    }
/*     First form  y := beta*y. */
    if (*beta != 1.) {
	if (*beta == 0.) {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = 0.;
	    }
	} else {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = *beta * y[i__];
	    }
	}
    }
    if (*alpha == 0.) {
/*         Relax */
    } else if (*task == 0) {
	i__1 = *lenx;
	for (k = 1; k <= i__1; ++k) {
	    alphxj = *alpha * x[k];
	    if (abs(alphxj) > *tolz) {
		j = kbs[k];
		if (j <= *n) {
/*                 ------------------- */
/*                 Column of A. */
/*                 ------------------- */
		    i__2 = loca[j + 1] - 1;
		    for (l = loca[j]; l <= i__2; ++l) {
			i__ = inda[l];
			y[i__] += acol[l] * alphxj;
		    }
		} else {
/*                 -------------------- */
/*                 Slack column. */
/*                 -------------------- */
		    i__ = j - *n;
		    y[i__] -= alphxj;
		}
	    }
	}
    } else if (*task == 1) {
	i__1 = *leny;
	for (k = 1; k <= i__1; ++k) {
	    t = 0.;
	    j = kbs[k];
	    if (j <= *n) {
/*              ------------------- */
/*              Column of A. */
/*              ------------------- */
		i__2 = loca[j + 1] - 1;
		for (l = loca[j]; l <= i__2; ++l) {
		    i__ = inda[l];
		    t += acol[l] * x[i__];
		}
	    } else {
/*              ------------------- */
/*              Slack column. */
/*              ------------------- */
		t = -x[j - *n];
	    }
	    y[k] += *alpha * t;
	}
    }
    return 0;
} /* s2bprd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bprd */
/* Subroutine */ int s2binf_(integer *nb, doublereal *bl, doublereal *bu, 
	doublereal *x, doublereal *binf, integer *jbinf)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal d1, d2;

/*     ================================================================== */
/*     s2bInf  computes the maximum infeasibility with respect to */
/*     the bounds on x. */
/*     s2bInf  is called by s5savB and s8savB before and after unscaling. */

/*     On exit, */
/*      bInf  is the maximum bound infeasibility. */
/*     jbInf  is the corresponding variable. */

/*     28 Jul 1999: First version based on Minos routine m2bInf. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;

    /* Function Body */
    *jbinf = 0;
    *binf = 0.;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	d1 = bl[j] - x[j];
	d2 = x[j] - bu[j];
	if (*binf < d1) {
	    *binf = d1;
	    *jbinf = j;
	}
	if (*binf < d2) {
	    *binf = d2;
	    *jbinf = j;
	}
    }
    return 0;
} /* s2binf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2bInf */
integer s2coln_(integer *j, integer *leniw, integer *iw)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer jn, lkxn;

/*     ================================================================== */
/*     s2ColN  gives the natural column number of Jacobian column j. */

/*     28 Dec 2000: First version written for snoptf */
/*     07 Jun 2001: Current version */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    lkxn = iw[252];
/* jN = kxN(j ) => col j of Jcol is variable iN */
    if (*j > 0) {
	jn = iw[lkxn + *j - 1];
    } else {
	jn = 0;
    }
    ret_val = jn;
    return ret_val;
} /* s2coln_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function s2ColN */
integer s2rown_(integer *i__, integer *leniw, integer *iw)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer j, n, in, lkxn;

/*     ================================================================== */
/*     s2RowN  gives the natural row number of Jacobian row i. */

/*     28 Dec 2000: First version written for snoptf */
/*     07 Jun 2001: Current version */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    n = iw[15];
/* copy of the number of columns */
    lkxn = iw[252];
/* jN = kxN(j ) => var j is natural variable iN */
    if (*i__ > 0) {
	j = n + *i__;
	in = iw[lkxn + j - 1];
    } else {
	in = 0;
    }
    ret_val = in;
    return ret_val;
} /* s2rown_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function s2RowN */
integer s2varn_(integer *j, integer *leniw, integer *iw)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer n, jj, jn, lkxn;

/*     ================================================================== */
/*     s2VarN  gives the natural variable number corresponding to var i. */

/*     28 Dec 2000: First version written for snoptf */
/*     17 Nov 2001: Current version */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    n = iw[15];
/* copy of the number of columns */
    lkxn = iw[252];
/* jN = kxN(j ) => var j is natural variable iN */
    if (*j >= 0) {
	jj = *j;
    } else {
	jj = -(*j);
    }
    if (jj > n) {
	jn = n + iw[lkxn + jj - 1];
    } else if (jj > 0) {
	jn = iw[lkxn + jj - 1];
    } else {
	jn = 0;
    }
    if (*j >= 0) {
	ret_val = jn;
    } else {
	ret_val = -jn;
    }
    return ret_val;
} /* s2varn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function s2VarN */
/* Subroutine */ int s2dinf_(integer *n, integer *nb, integer *iobj, 
	doublereal *tol, doublereal *bl, doublereal *bu, doublereal *rc, 
	doublereal *x, doublereal *dinf, integer *jdinf)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal dj;
    static integer jobj;
    static doublereal blobj, tolrel;

/*     ================================================================== */
/*     s2dInf  computes the maximum complementarity. */
/*     s2dInf  is called by s4savB before and after unscaling. */

/*     On exit, */
/*     dInf  is the maximum dual infeasibility. */
/*     jdInf  is the corresponding variable. */

/*     05 Apr 1996: First version based on Minos routine m2dInf. */
/*     29 Nov 2002: Switched to max complementarity. */
/*     26 Dec 2002: This version */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --rc;
    --bu;
    --bl;

    /* Function Body */
    jobj = *n + *iobj;
    if (*iobj > 0) {
	blobj = bl[jobj];
	bl[jobj] = bu[jobj];
    }
    *jdinf = 0;
    *dinf = 0.;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	if (bl[j] < bu[j]) {
	    tolrel = *tol * ((d__1 = x[j], abs(d__1)) + 1.);
	    dj = rc[j];
	    if (x[j] <= bl[j] + tolrel) {
		dj = -dj;
	    } else if (x[j] >= bu[j] - tolrel) {
/*              dj  = + dj */
	    } else {
		dj = abs(dj);
	    }
	    if (*dinf < dj) {
		*dinf = dj;
		*jdinf = j;
	    }
	}
    }
    if (*iobj > 0) {
	bl[jobj] = blobj;
    }
    return 0;
} /* s2dinf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2dInf */
/* Subroutine */ int s2gathr_(integer *neg, integer *nbs, integer *kbs, 
	doublereal *alpha, doublereal *g, doublereal *gbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;

/*     ================================================================== */
/*     s2gathr   performs the gather operation alpha*g  --> gBS */
/*     between vectors  g(neg)  and  gBS(nBS) according to */
/*     the index array  kBS. */

/*     The case alpha = 1 is treated specially. */

/*     On entry: 0 < nBS <= neg. */

/*     20 Jul 2000: First version of s2gathr */
/*     19 Dec 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --g;
    --gbs;
    --kbs;

    /* Function Body */
    if (*alpha == 1.) {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *neg) {
		gbs[k] = g[j];
	    } else {
		gbs[k] = 0.;
	    }
	}
    } else if (*alpha == -1.) {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *neg) {
		gbs[k] = -g[j];
	    } else {
		gbs[k] = 0.;
	    }
	}
    } else {
/* general alpha */
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *neg) {
		gbs[k] = *alpha * g[j];
	    } else {
		gbs[k] = 0.;
	    }
	}
    }
    return 0;
} /* s2gathr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2gathr */
/* Subroutine */ int s2scatr_(integer *neg, integer *nbs, integer *kbs, 
	doublereal *alpha, doublereal *gbs, doublereal *g)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);

/*     ================================================================== */
/*     s2scatr   performs the scatter operation alpha*gBS  --> g */
/*     between vectors  gBS(nBS)  and  g(neg)  according to */
/*     the index array  kBS. */

/*     The case alpha = 1 is treated specially. */

/*     On entry: 0 < nBS <= neg. */

/*     20 Jul 2000: First version of s2scatr */
/*     19 Dec 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --g;
    --gbs;
    --kbs;

    /* Function Body */
    dload_(neg, &c_b94, &g[1], &c__1);
    if (*alpha == 1.) {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *neg) {
		g[j] = gbs[k];
	    }
	}
    } else if (*alpha == -1.) {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *neg) {
		g[j] = -gbs[k];
	    }
	}
    } else {
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *neg) {
		g[j] = gbs[k] * *alpha;
	    }
	}
    }
    return 0;
} /* s2scatr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2scatr */
/* Subroutine */ int s2crsh_(integer *lcrash, integer *lprint, integer *m, 
	integer *n, integer *nb, integer *nncon, integer *icrash, doublereal *
	tcrash, integer *ne, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *hpiv, integer *hs, integer *hrtype, 
	doublereal *bl, doublereal *bu, doublereal *x, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_1200[] = "(\002 Slacks\002,i6,\002  Free cols\002,i6,"
	    "\002  Preferred\002,i6)";
    static char fmt_1210[] = "(\002 Unit  \002,i6,\002  Double   \002,i6,"
	    "\002  Triangle \002,i6,\002  Pad\002,i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k;
    static doublereal d1, d2;
    static integer i1, i2, k1, k2;
    static doublereal ai;
    static integer ip, js, nz, num[6];
    static char str[80];
    static doublereal eps0;
    static logical free;
    static integer npad;
    static doublereal apiv[2];
    static integer ipiv[2], npiv;
    static logical prnt1;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static doublereal aimax;
    static integer stage;
    static doublereal aitol;
    static integer nrows;
    static logical stage2, stage3, stage4, stage5;
    extern /* Subroutine */ int s2aprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer nbasic;
    static logical addslk, prefer, gotslk;
    static doublereal tolslk;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___77 = { 0, str, 0, "(a, i3)", 80, 1 };
    static icilist io___110 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___111 = { 0, str, 0, fmt_1210, 80, 1 };


/*     ================================================================== */
/*     s2crsh  looks for a basis in the columns of ( A  -I ). */

/*     ON ENTRY */

/*     iCrash    = the Crash option used by s5getB to set lCrash. */
/*     tCrash    = the Crash tolerance.  Default = 0.1 */

/*     lCrash      specifies the action to be taken by Crash. */
/*        0,1,2,3  The call is from s5getB. */
/*        4,5      The call is from s8solv. */

/*        0        The all-slack basis is set up. */

/*        1        A triangular Crash is applied to the columns of A. */
/*                 hs(1:n) is used to help select columns. */
/*                 tCrash is used to ignore small entries in each column. */
/*                 Depending on the size of tCrash, the resulting basis */
/*                 will be nearly (but not strictly) lower triangular. */

/*        2        As for 1, but nonlinear rows are ignored. */

/*        3        As for 2, but linear LG rows are also ignored. */

/*        4        Linear LG rows are now included. */
/*                 All hs(1:nb) and x(n+i) are defined. */
/*                 Slack values of x(n+i) are used to select LG rows. */

/*        5        Nonlinear rows are now included. */

/*     hrtype(*)   should be defined as described in s2Amat: */
/*     hrtype(i) = 0  for E rows      (equalities) */
/*     hrtype(i) = 1  for L or G rows (inequalities) */
/*     hrtype(i) = 2  for N rows      (objective or free rows) */

/*     x           If lCrash <= 4, x(1:n) is used to initialize */
/*                 slacks as x(n+1:nb) = A*x. */
/*                 Used to select slacks from LG rows to be in B (basis). */
/*                 If lCrash  = 5, x(n+1:n+nnCon) contains slack values */
/*                 evaluated from x(1:n) and Fx(*). */
/*                 Used to select slacks from nonlinear rows to be in B. */

/*     hs          If lCrash = 1, 2 or 3, hs(1:n)  is used. */
/*                 If lCrash =    4 or 5, hs(1:nb) is used. */
/*                 If hs(j) =  0, 1 or 3, column j is eligible for B, */
/*                                        with 3 being "preferred". */
/*                 If hs(j) =  2, 4 or 5, column j is ignored. */


/*     Crash has several stages. */

/*     Stage 1: Insert any slacks (N, L or G rows, hrtype = 1 or 2). */

/*     Stage 2: Do triangular Crash on any free columns (wide bounds) */

/*     Stage 3: Do triangular Crash on "preferred" columns (hs(j) < 0). */
/*              For the linear Crash, this includes variables set */
/*              between their bounds in the MPS file via FR INITIAL. */
/*              For the nonlinear Crash, it includes nonbasics */
/*              between their bounds. */
/*              (That is, "pegged" variables in both cases.) */

/*     Stage 4: Grab unit columns. */

/*     Stage 5: Grab double columns. */

/*     Stage 6: Do triangular Crash on all columns. */

/*     Slacks are then used to pad the basis. */


/*     ON EXIT */

/*     hs          is set to denote an initial (B S N) partition. */
/*                 hs(j) = 3 denotes variables for the initial basis. */
/*                 If hs(j) = 2 still, variable j will be superbasic. */
/*                 If hs(j) = 4 or 5 still, it will be changed to 0 or 1 */
/*                 by s4chek and variable j will be nonbasic. */

/*     x          If lCrash <= 4, slacks x(n+1:nb) are initialized. */

/*     ------------------------------------------------------------------ */
/*        Nov 1986: Essentially the same as in 1976. */
/*                  Crash tolerance added. */
/*                  Attention paid to various hs values. */

/*     12 Nov 1988: After free rows and columns have been processed */
/*                  (stage 1 and 2), slacks on L or G rows are inserted */
/*                  if their rows have not yet been assigned. */

/*     28 Dec 1988: Refined as follows. */
/*                  Stage 1 inserts free and preferred rows (slacks). */
/*                  Stage 2 performs a triangular Crash on free or */
/*                          preferred columns, ignoring free rows. */
/*                          Unpivoted L or G slacks are then inserted. */
/*                  Stage 3 performs a triangular Crash on the other */
/*                          columns, ignoring rows whose slack is basic. */
/*                          (Such rows form the top part of U.  The */
/*                          remaining rows form the beginning of L.) */

/*     30 Apr 1989: Stage 1 now also looks for singleton columns */
/*                  (ignoring free and preferred rows). */
/*     05 May 1989: Stage 2 doesn't insert slacks if Crash option < 0. */

/*     06 Dec 1989: Stage 2, 3, 4 modified.  Columns of length 2 are */
/*                  now treated specially. */

/*     20 Dec 1989: Stage 2 thru 5 modified.  Free columns done before */
/*                  unit and double columns. */

/*     19 May 1992: x now used to help initialize slacks. */
/*                  Stage 1 thru 7 redefined as above. */

/*     01 Jun 1992: abs used to define closeness of slacks to bounds. */
/*                  Unfortunately, x(1:n) seldom has meaningful values. */

/*     02 Jun 1992: Poor performance on the larger problems. */
/*                  Reverted to simple approach: all slacks grabbed. */

/*     04 Jun 1992: Compromise -- Crash 3 now has 3 phases: */
/*                  (a) E rows. */
/*                  (b) LG rows. */
/*                  (c) Nonlinear rows. */
/*                  x(1:n) should then define the slack values better */
/*                  for (b) and (c). */
/*     17 Jul 1996: Standard implementation for slacks. */
/*     28 Jul 2003: Print level reduced to 1. */
/*     31 Jul 2003: snPRNT adopted. */
/*     31 Jul 2003: Current version of s2crsh. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hrtype;
    --hpiv;
    --x;
    --bu;
    --bl;
    --hs;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5) */
    prnt1 = *lprint >= 1;
    if (prnt1) {
	if (*lcrash <= 3) {
	    s_wsfi(&io___77);
	    do_fio(&c__1, " Crash option", (ftnlen)13);
	    do_fio(&c__1, (char *)&(*icrash), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
	}
	if (*lcrash == 3 && *nncon < *m) {
	    snprnt_(&c__31, " Crash on linear E  rows:", &iw[1], leniw, (
		    ftnlen)25);
	}
	if (*lcrash == 4) {
	    snprnt_(&c__31, " Crash on linear LG rows:", &iw[1], leniw, (
		    ftnlen)25);
	}
	if (*lcrash == 5) {
	    snprnt_(&c__31, " Crash on nonlinear rows:", &iw[1], leniw, (
		    ftnlen)25);
	}
    }
    if (*lcrash <= 4) {
/*        Sets slacks x(n+1:nb) = A*x. */
/*        This is where the slacks are initialized. */
/*        They may be altered later (see the end of Crash). */
	s2aprd_(&c__0, &eps0, ne, nloca, &loca[1], &inda[1], &acol[1], &
		c_b109, &x[1], n, &c_b94, &x[*n + 1], m);
    }
/*     ------------------------------------------------------------------ */
/*     For Crash option 0, set hs(j) = 3 for all slacks and quit. */
/*     ------------------------------------------------------------------ */
    if (*lcrash == 0) {
	iload_(n, &c__0, &hs[1], &c__1);
	iload_(m, &c__3, &hs[*n + 1], &c__1);
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Crash option 1, 2 or 3.   lCrash = 1, 2, 3, 4, or 5. */
/*     tolslk measures closeness of slacks to bounds. */
/*     i1,i2  are the first and last rows of A involved in Stage 1. */
/*     ------------------------------------------------------------------ */
/* -->  tolslk = 0.25 */
    tolslk = .01;
    iload_(&c__6, &c__0, num, &c__1);
    if (*lcrash <= 3) {
/*        --------------------------------------------------------------- */
/*        First call.   lCrash = 1, 2 or 3. */
/*        Initialize hpiv(*) for all rows and hs(*) for all slacks. */
/*        --------------------------------------------------------------- */
	i1 = 1;
	if (*lcrash >= 2) {
	    i1 = *nncon + 1;
	}
	i2 = *m;
	nrows = i2 - i1 + 1;
/*        Make sure there are no basic columns already (hs(j) = 3). */
/*        If there are, make them "preferred". */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] == 3) {
		hs[j] = -1;
	    }
	}
/*        Make relevant rows available:  hpiv(i) = 1, hs(n+i) = 0. */
	if (nrows > 0) {
	    iload_(&nrows, &c__1, &hpiv[i1], &c__1);
	    iload_(&nrows, &c__0, &hs[*n + i1], &c__1);
	}
	if (*lcrash == 1) {
	    nbasic = 0;
	} else {
/*           lCrash = 2 or 3:  Insert nonlinear slacks. */
	    nbasic = *nncon;
	    if (*nncon > 0) {
		iload_(nncon, &c__3, &hpiv[1], &c__1);
		iload_(nncon, &c__3, &hs[*n + 1], &c__1);
	    }
	}
	if (*lcrash == 3) {
/*           Insert linear inequality slacks (including free rows). */
	    i__1 = *m;
	    for (i__ = i1; i__ <= i__1; ++i__) {
		if (hrtype[i__] >= 1) {
		    ++nbasic;
		    --nrows;
		    hpiv[i__] = 3;
		    hs[*n + i__] = 3;
		}
	    }
	}
/*        We're done if there are no relevant rows. */
	if (nrows == 0) {
	    goto L800;
	}
    } else {
/*        --------------------------------------------------------------- */
/*        Second or third call.  lCrash = 4 or 5. */
/*        Initialize hpiv(*) for all rows. */
/*        hs(*) already defines a basis for the full problem, */
/*        but we want to do better by including only some of the slacks. */
/*        --------------------------------------------------------------- */
	if (*lcrash == 4) {
/*           ------------------------------------------------------------ */
/*           Crash on linear LG rows. */
/*           ------------------------------------------------------------ */
	    if (*nncon == *m) {
		goto L900;
	    }
	    i1 = *nncon + 1;
	    i2 = *m;
/*           Mark nonlinear rows as pivoted: hpiv(i) = 3. */
	    nbasic = *nncon;
	    if (nbasic > 0) {
		iload_(&nbasic, &c__3, &hpiv[1], &c__1);
	    }
/*           Mark linear E  rows as pivoted: hpiv(i) = 3 */
/*           Make linear LG rows available:  hpiv(i) = 1, hs(n+i) = 0. */
	    i__1 = *m;
	    for (i__ = i1; i__ <= i__1; ++i__) {
		if (hrtype[i__] == 0) {
		    ++nbasic;
		    hpiv[i__] = 3;
		} else {
		    hpiv[i__] = 1;
		    hs[*n + i__] = 0;
		}
	    }
/*           Mark linear LG rows with hpiv(i) = 2 */
/*           if any basic columns contain a nonzero in row i. */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (hs[j] == 3) {
		    i__2 = loca[j + 1] - 1;
		    for (k = loca[j]; k <= i__2; ++k) {
			i__ = inda[k];
			if (hrtype[i__] == 1) {
			    if (i__ > *nncon) {
				if (acol[k] != 0.) {
				    hpiv[i__] = 2;
				}
			    }
			}
		    }
		}
	    }
	} else {
/*           ------------------------------------------------------------ */
/*           lCrash = 5.  Crash on nonlinear rows. */
/*           ------------------------------------------------------------ */
	    i1 = 1;
	    i2 = *nncon;
/*           Mark all linear rows as pivoted: hpiv(i) = 3 */
	    nbasic = *m - *nncon;
	    if (nbasic > 0) {
		iload_(&nbasic, &c__3, &hpiv[*nncon + 1], &c__1);
	    }
/*           Make nonlinear rows available:  hpiv(i) = 1, hs(n+i) = 0. */
	    iload_(nncon, &c__1, &hpiv[1], &c__1);
	    iload_(nncon, &c__0, &hs[*n + 1], &c__1);
/*           Mark nonlinear rows with hpiv(i) = 2 */
/*           if any basic columns contain a nonzero in row i. */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (hs[j] == 3) {
		    i__2 = loca[j + 1] - 1;
		    for (k = loca[j]; k <= i__2; ++k) {
			i__ = inda[k];
			if (i__ <= *nncon) {
			    if (acol[k] != 0.) {
				hpiv[i__] = 2;
			    }
			}
		    }
		}
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     Stage 1: Insert relevant slacks (N, L or G rows, hrtype = 1 or 2). */
/*              If lCrash = 4 or 5, grab them only if they are more than */
/*              tolslk from their bound. */
/*     ------------------------------------------------------------------ */
    stage = 1;
    gotslk = *lcrash == 4 || *lcrash == 5;
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
	j = *n + i__;
	if (hs[j] <= 1 && hrtype[i__] > 0) {
	    addslk = TRUE_;
	    if (gotslk) {
		d1 = x[j] - bl[j];
		d2 = bu[j] - x[j];
		if (min(d1,d2) <= tolslk) {
/*                 The slack is close to a bound or infeasible. */
/*                 Move it exactly onto the bound. */
		    addslk = FALSE_;
		    if (d1 <= d2) {
			x[j] = bl[j];
			hs[j] = 0;
		    } else {
			x[j] = bu[j];
			hs[j] = 1;
		    }
		}
	    }
	    if (addslk) {
		++nbasic;
		++num[stage - 1];
		hpiv[i__] = 3;
		hs[j] = 3;
	    }
	}
    }
    if (nbasic == *m) {
	goto L700;
    }
/*     ------------------------------------------------------------------ */
/*     Apply a triangular Crash to various subsets of the columns of A. */

/*        hpiv(i) = 1  if row i is unmarked (initial state). */
/*        hpiv(i) = 3  if row i has been given a pivot */
/*                     in one of a set of triangular columns. */
/*        hpiv(i) = 2  if one of the triangular columns contains */
/*                     a nonzero in row i below the triangle. */
/*     ------------------------------------------------------------------ */
    for (stage = 2; stage <= 6; ++stage) {
	stage2 = stage == 2;
	stage3 = stage == 3;
	stage4 = stage == 4;
	stage5 = stage == 5;
/*        --------------------------------------------------------------- */
/*        Main loop for triangular Crash. */
/*        --------------------------------------------------------------- */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    js = hs[j];
	    if (js > 1) {
		goto L200;
	    }
	    if (bl[j] == bu[j]) {
		goto L200;
	    }
	    if (stage2) {
		free = bl[j] <= -1e4 && bu[j] >= 1e4;
		if (! free) {
		    goto L200;
		}
	    } else if (stage3) {
		prefer = js < 0;
		if (! prefer) {
		    goto L200;
		}
	    }
/*           Find the biggest aij, ignoring free rows. */
	    k1 = loca[j];
	    k2 = loca[j + 1] - 1;
	    aimax = 0.;
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
		i__ = inda[k];
		if (hrtype[i__] != 2) {
		    ai = acol[k];
/* Computing MAX */
		    d__1 = aimax, d__2 = abs(ai);
		    aimax = max(d__1,d__2);
		}
	    }
/*           Prevent small pivots if Crash tol is too small. */
	    if (aimax <= .001) {
		goto L200;
	    }
/*           Find the biggest pivots in rows that are still */
/*           unpivoted and unmarked.  Ignore smallish elements. */
/*           nz counts the number of relevant nonzeros. */
	    aitol = aimax * *tcrash;
	    nz = 0;
	    npiv = 0;
	    ipiv[0] = 0;
	    ipiv[1] = 0;
	    apiv[0] = 0.;
	    apiv[1] = 0.;
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
		i__ = inda[k];
		if (hs[*n + i__] != 3) {
		    ai = (d__1 = acol[k], abs(d__1));
		    if (ai > aitol) {
			++nz;
			ip = hpiv[i__];
			if (ip <= 2) {
			    if (apiv[ip - 1] < ai) {
				apiv[ip - 1] = ai;
				ipiv[ip - 1] = i__;
			    }
			} else {
			    ++npiv;
			}
		    }
		}
	    }
/*           Grab unit or double columns. */
	    if (stage4) {
		if (nz != 1) {
		    goto L200;
		}
	    } else if (stage5) {
		if (nz != 2) {
		    goto L200;
		}
	    }
/*           See if the column contained a potential pivot. */
/*           An unmarked row is favored over a marked row. */
	    ip = 1;
	    if (ipiv[0] == 0 && npiv == 0) {
		ip = 2;
	    }
	    i__ = ipiv[ip - 1];
	    if (i__ > 0) {
		++nbasic;
		++num[stage - 1];
		hpiv[i__] = 3;
		hs[j] = 3;
		if (nbasic >= *m) {
		    goto L700;
		}
/*              Mark off all relevant unmarked rows. */
		i__2 = k2;
		for (k = k1; k <= i__2; ++k) {
		    i__ = inda[k];
		    if (hs[*n + i__] != 3) {
			ai = (d__1 = acol[k], abs(d__1));
			if (ai > aitol) {
			    if (hpiv[i__] == 1) {
				hpiv[i__] = 2;
			    }
			}
		    }
		}
	    }
L200:
	    ;
	}
    }
/*     ------------------------------------------------------------------ */
/*     All stages finished. */
/*     Fill remaining gaps with slacks. */
/*     ------------------------------------------------------------------ */
L700:
    npad = *m - nbasic;
    if (prnt1) {
	s_wsfi(&io___110);
	do_fio(&c__1, (char *)&num[0], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&num[1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&num[2], (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	s_wsfi(&io___111);
	do_fio(&c__1, (char *)&num[3], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&num[4], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&num[5], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&npad, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
    }
    if (npad > 0) {
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (hpiv[i__] < 3) {
		++nbasic;
		hs[*n + i__] = 3;
		if (nbasic >= *m) {
		    goto L800;
		}
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     Make sure there aren't lots of nonbasic slacks floating off */
/*     their bounds.  They could take lots of iterations to move. */
/*     ------------------------------------------------------------------ */
L800:
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
	j = *n + i__;
	if (hs[j] <= 1 && hrtype[i__] > 0) {
	    d1 = x[j] - bl[j];
	    d2 = bu[j] - x[j];
	    if (min(d1,d2) <= tolslk) {
/*              The slack is close to a bound or infeasible. */
/*              Move it exactly onto the bound. */
		if (d1 <= d2) {
		    x[j] = bl[j];
		    hs[j] = 0;
		} else {
		    x[j] = bu[j];
		    hs[j] = 1;
		}
	    }
	}
    }
L900:
    return 0;
} /* s2crsh_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2crsh */
/* Subroutine */ int s2mem0_(integer *iexit, char *solver, integer *lencw, 
	integer *leniw, integer *lenrw, integer *iw, integer *mincw, integer *
	miniw, integer *minrw, integer *maxcw, integer *maxiw, integer *maxrw,
	 integer *nextcw, integer *nextiw, integer *nextrw, ftnlen solver_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static char str[80], str2[80];
    static integer mxcu, mxcw, mxiu, mxiw, mxru, mxrw, maxcu, maxiu, maxru, 
	    mincu1, mincu2, maxcu1, maxcu2, miniu1, miniu2, maxiu1, maxiu2, 
	    minru1, minru2, maxru1, maxru2;
    static logical fixdup;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), snprnt_(integer *, 
	    char *, integer *, integer *, ftnlen);

/*     ================================================================== */
/*     s2Mem0   checks the memory requirements for sqopt/snopt and sets */
/*     the pointers to the beginning and end of sqopt/snopt workspace. */

/*     Note: cw, iw and rw hold constants and work-space addresses. */
/*        They must have dimension at least 500. */

/*     The SPECS file has been read, and values are known for */
/*     maxcu, maxiu, maxru  (upper limit of user  partition 1) */
/*     maxcw, maxiw, maxrw  (upper limit of SNOPT partition) */

/*     The default values for these values are */
/*     maxcu = 500  ,   maxiu = 500  ,   maxru = 500, */
/*     maxcw = lencw,   maxiw = leniw,   maxrw = lenrw, */
/*     which are set in snInit: */

/*     The user can alter these in the SPECS file via */
/*     lines of the form */

/*        User  character workspace      10000    (Sets maxcu) */
/*        User  integer   workspace      10000    (Sets maxiu) */
/*        User  real      workspace      10000    (Sets maxru) */
/*        Total character workspace      90000    (Sets maxcw) */
/*        Total integer   workspace      90000    (Sets maxiw) */
/*        Total real      workspace      90000    (Sets maxrw) */

/*     SNOPT will use only rw(1:500) and rw(maxru+1:maxrw). */
/*     Hence, rw(501:maxru) and possibly rw(maxrw+1:lenrw) may be used as */
/*     workspace by the user during solution of the problem (e.g., within */
/*     the user-supplied function routines).  Similarly for the integer */
/*     and character work arrays iw and cw. */

/*     Setting maxiw and maxrw less than leniw and lenrw may serve to */
/*     reduce paging activity on a machine with virtual memory, by */
/*     confining SNOPT (in particular the LU-factorization routines) to */
/*     an area of memory that is sensible for the current problem.  This */
/*     often allows cw(*), iw(*) and rw(*) to be declared arbitrarily */
/*     large at compile time. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2core. */
/*     29 Mar 1998: First version called by snMem. This simplified */
/*                  version may slightly overestimate needed memory. */
/*     27 Apr 1999: Renamed s8Mem. */
/*     16 Jul 2003: Standard output added */
/*     31 Jul 2003: snEXIT and snPRNT adopted. */
/*     27 Oct 2003: Renamed s2Mem0. */
/*     15 Oeb 2004: Current version of s2Mem0. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    *iexit = 0;
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	*mincw = 500;
	*miniw = 500;
	*minrw = 500;
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
    mxcu = iw[6];
/* maxcu+1 is the start of SNOPT part of cw */
    mxiu = iw[4];
/* maxiu+1 is the start of SNOPT part of iw */
    mxru = iw[2];
/* maxru+1 is the start of SNOPT part of rw */
    mxcw = iw[7];
/* end of SNOPT part of cw */
    mxiw = iw[5];
/* end of SNOPT part of iw */
    mxrw = iw[3];
/*     Check for silly values. */
/* end of SNOPT part of rw */
/* Computing MIN */
    i__1 = max(mxcu,500);
    maxcu = min(i__1,*lencw);
/* Computing MIN */
    i__1 = max(mxiu,500);
    maxiu = min(i__1,*leniw);
/* Computing MIN */
    i__1 = max(mxru,500);
    maxru = min(i__1,*lenrw);
/* Computing MIN */
    i__1 = max(mxcw,500);
    *maxcw = min(i__1,*lencw);
/* Computing MIN */
    i__1 = max(mxiw,500);
    *maxiw = min(i__1,*leniw);
/* Computing MIN */
    i__1 = max(mxrw,500);
    *maxrw = min(i__1,*lenrw);
    maxcu = min(maxcu,*maxcw);
    maxiu = min(maxiu,*maxiw);
    maxru = min(maxru,*maxrw);
    fixdup = mxcu != maxcu || mxcw != *maxcw || mxiu != maxiu || mxiw != *
	    maxiw || mxru != maxru || mxrw != *maxrw;
    if (fixdup) {
	snprnt_(&c__14, " XXX  User workspace parameters had to be modified", 
		&iw[1], leniw, (ftnlen)50);
    }
/*     Save the checked values */
    iw[6] = maxcu;
/* maxcu+1 is the start of SNOPT part of cw */
    iw[4] = maxiu;
/* maxiu+1 is the start of SNOPT part of iw */
    iw[2] = maxru;
/* maxru+1 is the start of SNOPT part of rw */
    iw[7] = *maxcw;
/* end of SNOPT part of cw */
    iw[5] = *maxiw;
/* end of SNOPT part of iw */
    iw[3] = *maxrw;
/*     ------------------------------------------------------------------ */
/*     Save the limits of the two user-accessible workspace partitions */
/*     to allow the user can grab them from iw. */

/*     Upper limits of partition 1 may be set in the specs file. */
/*     lower limits of partition 2 are set here. */
/*     ------------------------------------------------------------------ */
/* end of SNOPT part of rw */
    mincu1 = 501;
/* Lower limits on partition 1 */
    miniu1 = 501;
    minru1 = 501;
    maxcu1 = maxcu;
/* User-defined upper limits on */
    maxiu1 = maxiu;
/* partition 1 (default 500). */
    maxru1 = maxru;
    maxcu2 = *lencw;
/* Upper limits on partition 2 */
    maxiu2 = *leniw;
    maxru2 = *lenrw;
    iw[31] = mincu1;
/* Start of first  user partition of cw */
    iw[36] = miniu1;
/* Start of first  user partition of iw */
    iw[41] = minru1;
/* Start of first  user partition of rw */
    iw[32] = maxcu1;
/* End   of first  user partition of cw */
    iw[37] = maxiu1;
/* End   of first  user partition of iw */
    iw[42] = maxru1;
/* End   of first  user partition of rw */
    iw[34] = maxcu2;
/* End   of second user partition of cw */
    iw[39] = maxiu2;
/* End   of second user partition of iw */
    iw[44] = maxru2;
/* End   of second user partition of rw */
    mincu2 = *maxcw + 1;
/* Lower limits on partition 2 */
    miniu2 = *maxiw + 1;
    minru2 = *maxrw + 1;
    iw[33] = mincu2;
/* Start of second user partition of cw */
    iw[38] = miniu2;
/* Start of second user partition of iw */
    iw[43] = minru2;
/*     ------------------------------------------------------------------ */
/*     Define the start of the character, integer and real workspace. */
/*     ------------------------------------------------------------------ */
/* Start of second user partition of rw */
    *nextcw = maxcu + 1;
    *nextiw = maxiu + 1;
    *nextrw = maxru + 1;
L999:
    return 0;
} /* s2mem0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Mem0 */
/* Subroutine */ int s2mem_(integer *iexit, logical *prtmem, integer *liwest, 
	integer *lrwest, integer *nextcw, integer *nextiw, integer *nextrw, 
	integer *maxcw, integer *maxiw, integer *maxrw, integer *lencw, 
	integer *leniw, integer *lenrw, integer *mincw, integer *miniw, 
	integer *minrw, integer *iw)
{
    /* Format strings */
    static char fmt_1110[] = "(\002 Total char*8  workspace\002,i10,6x,\002 "
	    "Total integer workspace\002,i10,6x,\002 Total real    workspac"
	    "e\002,i10)";
    static char fmt_1120[] = "(\002 Total char*8  (minimum)\002,i10,6x,\002 "
	    "Total integer (minimum)\002,i10,6x,\002 Total real    (minimum"
	    ")\002,i10)";
    static char fmt_1201[] = "(\002 Elements cw(\002,i10,\002:\002,i10,\002"
	    ")\002,6x,\002are free\002,\002 for USER CHAR*8  WORKSPACE\002)";
    static char fmt_1202[] = "(\002 Elements iw(\002,i10,\002:\002,i10,\002"
	    ")\002,6x,\002are free\002,\002 for USER INTEGER WORKSPACE\002)";
    static char fmt_1203[] = "(\002 Elements rw(\002,i10,\002:\002,i10,\002"
	    ")\002,6x,\002are free\002,\002 for USER REAL    WORKSPACE\002)";
    static char fmt_9420[] = "(\002 Total character workspace should be sign"
	    "ificantly\002,\002 more than\002,i8)";
    static char fmt_9430[] = "(\002 Total integer   workspace  should be sig"
	    "nificantly\002,\002 more than\002,i8)";
    static char fmt_9440[] = "(\002 Total real      workspace  should be sig"
	    "nificantly\002,\002 more than\002,i8)";
    static char fmt_9500[] = "(24x,\002        Current    Recommended\002)";
    static char fmt_9510[] = "(\002 Total integer workspace\002,2i15)";
    static char fmt_9520[] = "(\002 Total real    workspace\002,2i15)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[132];
    static integer mincu1, maxcu1, miniu1, maxiu1, minru1, maxru1, lenalu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___143 = { 0, str, 0, fmt_1110, 132, 1 };
    static icilist io___144 = { 0, str, 0, fmt_1120, 132, 1 };
    static icilist io___145 = { 0, str, 0, fmt_1201, 132, 1 };
    static icilist io___146 = { 0, str, 0, fmt_1202, 132, 1 };
    static icilist io___147 = { 0, str, 0, fmt_1203, 132, 1 };
    static icilist io___148 = { 0, str, 0, fmt_9420, 132, 1 };
    static icilist io___149 = { 0, str, 0, fmt_9430, 132, 1 };
    static icilist io___150 = { 0, str, 0, fmt_9440, 132, 1 };
    static icilist io___152 = { 0, str, 0, fmt_9500, 132, 1 };
    static icilist io___153 = { 0, str, 0, fmt_9510, 132, 1 };
    static icilist io___154 = { 0, str, 0, fmt_9520, 132, 1 };


/*     ================================================================== */
/*     s2Mem  prints details of the workspace. */

/*     On exit. */
/*        If iExit = 0,  mincw, miniw, minrw give the amounts of */
/*        character, integer and real workspace needed to hold */
/*        the problem. (The LU factorization routines may */
/*        subsequently ask for more.) */

/*        If iExit > 0,  insufficient storage is provided to hold the */
/*        problem.  In this case, mincw, miniw and minrw give estimates */
/*        of reasonable lengths for cw(*), iw(*) and rw(*). */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2core. */
/*     29 Mar 1998: First version called by snMem. This simplified */
/*                  version may slightly overestimate needed memory. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     04 Jul 2010: mincw, miniw, minrw added to workspace. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    *iexit = 0;
/*     Compute the minimum storage required */
    *mincw = *nextcw - 1;
    *miniw = *nextiw - 1;
    *minrw = *nextrw - 1;
/*     ------------------------------------------------------------------ */
/*     Print details of the workspace. */
/*     ------------------------------------------------------------------ */
    mincu1 = iw[31];
/* Start of first  user partition of cw */
    miniu1 = iw[36];
/* Start of first  user partition of iw */
    minru1 = iw[41];
/* Start of first  user partition of rw */
    maxcu1 = iw[32];
/* End   of first  user partition of cw */
    maxiu1 = iw[37];
/* End   of first  user partition of iw */
    maxru1 = iw[42];
/* End   of first  user partition of rw */
    if (*prtmem) {
	s_wsfi(&io___143);
	do_fio(&c__1, (char *)&(*maxcw), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*maxiw), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*maxrw), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___144);
	do_fio(&c__1, (char *)&(*mincw), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*miniw), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*minrw), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__21, " ", &iw[1], leniw, (ftnlen)1);
	if (maxcu1 >= mincu1) {
	    s_wsfi(&io___145);
	    do_fio(&c__1, (char *)&mincu1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&maxcu1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
	}
	if (maxiu1 >= miniu1) {
	    s_wsfi(&io___146);
	    do_fio(&c__1, (char *)&miniu1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&maxiu1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
	}
	if (maxru1 >= minru1) {
	    s_wsfi(&io___147);
	    do_fio(&c__1, (char *)&minru1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&maxru1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    iw[47] = *mincw;
/* minimum length of cw */
    iw[48] = *miniw;
/* minimum length of iw */
    iw[49] = *minrw;
/* minimum length of rw */
    if (*mincw > *maxcw || *miniw > *maxiw || *minrw > *maxrw) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to solve the problem. */
/*        --------------------------------------------------------------- */
	if (*prtmem) {
	    snprnt_(&c__3, " XXX  Not enough storage to start the problem...",
		     &iw[1], leniw, (ftnlen)48);
	}
    }
    if (*mincw > *lencw) {
	*iexit = 82;
/* Not enough character workspace. */
	if (*prtmem) {
	    s_wsfi(&io___148);
	    do_fio(&c__1, (char *)&(*mincw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*miniw > *leniw) {
	*iexit = 83;
/* Not enough integer workspace. */
	*miniw = *liwest;
	if (*prtmem) {
	    s_wsfi(&io___149);
	    do_fio(&c__1, (char *)&(*miniw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*minrw > *lenrw) {
	*iexit = 84;
/* Not enough real workspace. */
	*minrw = *lrwest;
	if (*prtmem) {
	    s_wsfi(&io___150);
	    do_fio(&c__1, (char *)&(*minrw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/*     LUSOL stores  indc(*), indr(*) and A(*) in  iw(miniw:maxiw) and */
/*     rw(minrw:maxrw), with miniw pointing to the start of */
/*     indc(*), indr(*), and minrw pointing to the start of A(*). */
    lenalu = iw[213];
    if (*iexit == 0 && lenalu == 0) {
/*        ------------------------------------------- */
/*        Insufficient storage to factorize B. */
/*        ------------------------------------------- */
	*iexit = 82;
	*miniw = *liwest;
	*minrw = *lrwest;
	if (*prtmem) {
	    snprnt_(&c__3, " XXX  Not enough storage for the basis factors", &
		    iw[1], leniw, (ftnlen)46);
	    s_wsfi(&io___152);
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___153);
	    do_fio(&c__1, (char *)&(*maxiw), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*liwest), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	    s_wsfi(&io___154);
	    do_fio(&c__1, (char *)&(*maxrw), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*lrwest), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    return 0;
} /* s2mem_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Mem */
/* Subroutine */ int s2rca_(logical *feasbl, doublereal *featol, integer *
	iobj, integer *minimz, doublereal *wtinf, integer *m, integer *n, 
	integer *nb, integer *leng, integer *neg, integer *ne, integer *nloca,
	 integer *loca, integer *inda, doublereal *acol, integer *hestat, 
	integer *hs, doublereal *bl, doublereal *bu, doublereal *g, 
	doublereal *pi, doublereal *rc, doublereal *x)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal d1, d2, dj;
    static integer jobj;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal sgnobj;

/*     ================================================================== */
/*     s2rcA  computes reduced costs rc(*) for all columns of ( A  -I ). */
/*     If x is feasible, the gradient includes the explicit objective */
/*     vector  g.  Otherwise, the phase 1 objective is included. */

/*     s2rcA  is called by s4savb and s8SQP for scaled data. */
/*     External values of hs(*) are used (0, 1, 2, 3), */
/*     but internal would be ok too since we only test for > 1. */

/*     19 Feb 1994: First version based on Minos routine m4rc. */
/*     17 Jul 1996: Standard implementation for slacks. */
/*     31 Jul 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --x;
    --rc;
    --bu;
    --bl;
    --hs;
    --hestat;
    --g;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    sgnobj = (doublereal) (*minimz);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dj = 0.;
	i__2 = loca[j + 1] - 1;
	for (l = loca[j]; l <= i__2; ++l) {
	    i__ = inda[l];
	    dj += pi[i__] * acol[l];
	}
	rc[j] = -dj;
    }
    dcopy_(m, &pi[1], &c__1, &rc[*n + 1], &c__1);
    if (*feasbl) {
/*        Include the nonlinear objective gradient. */
/*        Include the gradient of the linear term. */
	sgnobj = (doublereal) (*minimz);
	if (*neg > 0) {
	    daxpy_(neg, &sgnobj, &g[1], &c__1, &rc[1], &c__1);
	}
	if (*iobj > 0) {
	    jobj = *n + *iobj;
	    rc[jobj] += sgnobj;
	}
    } else {
/*        Include the Phase 1 objective. */
/*        Only basics and superbasics can be infeasible. */
/*        Check that this works for scaling. */
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] > 1) {
		d1 = bl[j] - x[j];
		d2 = x[j] - bu[j];
		if (hestat[j] == 0) {
		    if (d1 > *featol) {
			rc[j] += -1.;
		    }
		    if (d2 > *featol) {
			rc[j] += 1.;
		    }
		} else {
		    if (d1 > *featol) {
			rc[j] -= *wtinf;
		    }
		    if (d2 > *featol) {
			rc[j] += *wtinf;
		    }
		}
	    }
	}
    }
    return 0;
} /* s2rca_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2rcA */
/* Subroutine */ int s2scal_(integer *lprint, integer *m, integer *n, integer 
	*nb, integer *nnl, integer *nncon, integer *nnjac, integer *hrtype, 
	integer *ne, integer *nloca, integer *loca, integer *inda, doublereal 
	*acol, doublereal *ascale, doublereal *bl, doublereal *bu, doublereal 
	*rmin, doublereal *rmax, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Format strings */
    static char fmt_1200[] = "(\002 After\002,i3,1p,e12.2,e12.2,0p,f20.2)";
    static char fmt_1310[] = "(1x,a,i7,1p,e10.1,12x,a,i7,e10.1,i17,0p,f8.1)";
    static char fmt_1400[] = "(\002 Norm of fixed columns and slacks\002,1p,"
	    "e20.1)";
    static char fmt_1410[] = "(\002 (before and after row scaling)  \002,1p,"
	    "e20.1)";
    static char fmt_1450[] = "(\002 Scales are large --- reduced by \002,1p,"
	    "e20.1)";
    static char fmt_1500[] = "(i6,g16.5)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, l;
    static doublereal b1, b2, ac, ar, big, bnd;
    static char str[80];
    static doublereal amin, amax, cmin, cmax;
    static integer imin, jmin, imax, jmax;
    static doublereal tolx;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), dddiv_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal sigma, small;
    static logical prnt10;
    static integer npass;
    static doublereal bplus;
    static logical lonly;
    static doublereal close1, close2, infbnd, aratio, cratio;
    static integer nclose;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer lprscl, lvlscl;
    static doublereal sratio, scltol;
    static integer istart, jstart, mxpass;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___189 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___195 = { 0, str, 0, fmt_1310, 80, 1 };
    static icilist io___198 = { 0, str, 0, fmt_1310, 80, 1 };
    static icilist io___203 = { 0, str, 0, fmt_1400, 80, 1 };
    static icilist io___204 = { 0, str, 0, fmt_1410, 80, 1 };
    static icilist io___206 = { 0, str, 0, fmt_1450, 80, 1 };
    static icilist io___207 = { 0, str, 0, fmt_1500, 80, 1 };
    static icilist io___208 = { 0, str, 0, fmt_1500, 80, 1 };


/*     ================================================================== */
/*     s2Scal computes scale factors  Ascale  from  A, bl, bu. */

/*     In phase 1, an iterative procedure based on geometric means is */
/*     used to compute scales from a alone.  This procedure is derived */
/*     from a routine written by Robert Fourer, 1979.  The main steps */
/*     are: */

/*        (1) Compute aRatio = max(i1,i2,j)  |A(i1,j)| / |A(i2,j)|. */
/*            i.e., the max ratio of the magnitudes of any two nonzeros */
/*            elements in the same column of A. */

/*        (2) Divide each row i by */
/*               [min_j A(i,j)*max_j A(i,j)]**1/2. */
/*            taking the min over all nonzero elements in row i. */

/*        (3) Divide each column j by */
/*               [min_i A(i,j)*max_i A(i,j)]**1/2. */
/*            taking the min over all nonzero elements in column i. */

/*        (4) Compute sRatio as the max ratio of the magnitudes of any */
/*            two nonzeros elements in the same column of A, as scaled. */

/*        (5) If sRatio .lt. tolScale*aRatio, repeat from step (1). */

/*        Free rows (hrtype=2) and fixed columns (bl=bu) are not used */
/*        at this stage. */

/*     In phase 2, the scales for free rows are set to be their largest */
/*     element. */

/*     In phase 3, fixed columns are summed in order to compute */
/*     a scale factor sigma that allows for the effective rhs of the */
/*     constraints.  All scales are then multiplied by sigma. */


/*     If lvlScl = 1, the first nnCon rows and the first nnJac columns will */
/*     retain scales of 1.0 during phases 1-2, and phase 3 will not be */
/*     performed.  (This takes effect if the problem is nonlinear but */
/*     the user has specified 'scale linear variables' only.) */
/*     However, all rows    contribute to the linear column scales, */
/*     and      all columns contribute to the linear row    scales. */

/*     If lvlScl = 2, all rows and columns are scaled.  To guard against */
/*     misleadingly small Jacobians, if the maximum element in any of */
/*     the first nnCon rows or the first nnL columns is less than */
/*     smallj, the corresponding row or column retains a scale of 1.0. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2scal. */
/*     24 Jan 2001: Don't scale BIG free rows anymore! */
/*     31 Jul 2003: snPRNT adopted. */
/*     26 Nov 2013: Current version of s2Scal. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --rmax;
    --rmin;
    --hrtype;
    --bu;
    --bl;
    --ascale;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    infbnd = rw[70];
/* definition of an infinite bound */
    lvlscl = iw[75];
/* scale option */
    lprscl = iw[83];
/* > 0    => print the scales */
    scltol = rw[92];
/* scale tolerance. */
    prnt10 = *lprint >= 10;
    if (prnt10) {
	snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
	snprnt_(&c__1, " Scaling", &iw[1], leniw, (ftnlen)8);
	snprnt_(&c__1, " -------", &iw[1], leniw, (ftnlen)8);
	snprnt_(&c__1, "             Min elem    Max elem       Max col ratio"
		, &iw[1], leniw, (ftnlen)53);
    }
    bplus = infbnd * .1;
    aratio = bplus;
    mxpass = 10;
    lonly = lvlscl == 1;
    dload_(nb, &c_b109, &ascale[1], &c__1);
    if (lonly) {
/* if (nnL .ge. n) return */
	istart = *nncon + 1;
/* jStart = nnL   + 1 */
	jstart = *nnjac + 1;
    } else {
	istart = 1;
	jstart = 1;
    }
/*     ------------------------------------------------------------------ */
/*     Main loop for phase 1. */
/*     Only the following row-types are used: */
/*        hrtype(i) = 2       for type N rows (objective or free rows), */
/*        hrtype(i) = 0 or 1  otherwise. */
/*     ------------------------------------------------------------------ */
    i__1 = mxpass;
    for (npass = 0; npass <= i__1; ++npass) {
/*        Find the largest column ratio. */
/*        Also set new column scales (except on pass 0). */
	amin = bplus;
	amax = 0.;
	small = .01;
	sratio = 1.;
	i__2 = *n;
	for (j = jstart; j <= i__2; ++j) {
	    if (bl[j] < bu[j]) {
		cmin = bplus;
		cmax = 0.;
		i__3 = loca[j + 1] - 1;
		for (l = loca[j]; l <= i__3; ++l) {
		    i__ = inda[l];
		    if (hrtype[i__] != 2) {
			ar = (d__1 = acol[l], abs(d__1));
			if (ar > 0.) {
			    ar /= ascale[*n + i__];
			    cmin = min(cmin,ar);
			    cmax = max(cmax,ar);
			}
		    }
		}
/* Computing MAX */
		d__1 = cmin, d__2 = cmax * 1e-4;
		ac = max(d__1,d__2);
		ac = sqrt(ac) * sqrt(cmax);
		if (j > *nnjac) {
		    small = 0.;
		}
		if (cmax <= small) {
		    ac = 1.;
		}
		if (npass > 0) {
		    ascale[j] = ac;
		}
/* Computing MIN */
		d__1 = amin, d__2 = cmin / ascale[j];
		amin = min(d__1,d__2);
/* Computing MAX */
		d__1 = amax, d__2 = cmax / ascale[j];
		amax = max(d__1,d__2);
		cratio = cmax / cmin;
		sratio = max(sratio,cratio);
	    }
	}
	if (prnt10) {
	    s_wsfi(&io___189);
	    do_fio(&c__1, (char *)&npass, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&amin, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&amax, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&sratio, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
/*        --------------------------------------------------------------- */
/*        Test for convergence. */
/*        --------------------------------------------------------------- */
	if (npass >= 3 && sratio >= aratio * scltol) {
	    goto L420;
	}
/* Break */
	if (npass < mxpass) {
	    aratio = sratio;
/*           Set new row scales for the next pass. */
	    if (istart <= *m) {
		i__2 = *m - istart + 1;
		dload_(&i__2, &bplus, &rmin[istart], &c__1);
		i__2 = *m - istart + 1;
		dload_(&i__2, &c_b94, &rmax[istart], &c__1);
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    if (bl[j] < bu[j]) {
			ac = ascale[j];
			i__3 = loca[j + 1] - 1;
			for (l = loca[j]; l <= i__3; ++l) {
			    i__ = inda[l];
			    if (i__ >= istart) {
				ar = (d__1 = acol[l], abs(d__1));
				if (ar > 0.) {
				    ar /= ac;
/* Computing MIN */
				    d__1 = rmin[i__];
				    rmin[i__] = min(d__1,ar);
/* Computing MAX */
				    d__1 = rmax[i__];
				    rmax[i__] = max(d__1,ar);
				}
			    }
			}
		    }
		}
		i__2 = *m;
		for (i__ = istart; i__ <= i__2; ++i__) {
		    j = *n + i__;
		    ar = rmax[i__];
		    if (i__ <= *nncon && ar <= .01) {
			ascale[j] = 1.;
		    } else {
/* Computing MAX */
			d__1 = rmin[i__], d__2 = ar * 1e-4;
			ac = max(d__1,d__2);
			ascale[j] = sqrt(ac) * sqrt(ar);
		    }
		}
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     End of main loop. */
/*     ------------------------------------------------------------------ */
/*     Invert the column scales, so that structurals and logicals */
/*     can be treated the same way during subsequent unscaling. */
/*     Find the min and max column scales while we're at it. */
/*     Nov 1989: nclose counts how many are "close" to 1. */
/*     For problems that are already well-scaled, it seemed sensible to */
/*     set the "close" ones exactly equal to 1. */
/*     Tried "close" = (0.5,2.0) and (0.9,1.1), but they helped only */
/*     occasionally.  Have decided not to interfere. */
L420:
    amin = bplus;
    amax = 0.;
    close1 = .5;
    close2 = 2.;
    nclose = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ac = 1. / ascale[j];
	if (amin > ac) {
	    amin = ac;
	    jmin = j;
	}
	if (amax < ac) {
	    amax = ac;
	    jmax = j;
	}
	if (ac > close1 && ac < close2) {
	    ++nclose;
/* ----        Ac     =  one */
	}
	ascale[j] = ac;
    }
/*     Remember, column scales are upside down. */
    amax = 1. / amax;
    amin = 1. / amin;
    if (prnt10) {
	snprnt_(&c__11, "            Min scale                       Max sca"
		"le      Between 0.5 and 2.0", &iw[1], leniw, (ftnlen)78);
	s_wsfi(&io___195);
	do_fio(&c__1, "Col", (ftnlen)3);
	do_fio(&c__1, (char *)&jmax, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&amax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "Col", (ftnlen)3);
	do_fio(&c__1, (char *)&jmin, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&amin, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&nclose, (ftnlen)sizeof(integer));
	d__1 = nclose * 100. / *n;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    }
/*     ------------------------------------------------------------------ */
/*     Phase 2.  Deal with empty rows and free rows. */
/*     Find the min and max row scales while we're at it. */
/*     ------------------------------------------------------------------ */
    amin = bplus;
    amax = 0.;
    imin = 0;
    imax = 0;
    nclose = 0;
    i__1 = *m;
    for (i__ = istart; i__ <= i__1; ++i__) {
	j = *n + i__;
	if (hrtype[i__] == 2) {
/* Computing MIN */
	    d__1 = rmax[i__];
	    ar = min(d__1,1.);
	    if (ar == 0.) {
		ar = 1.;
	    }
	} else {
	    ar = ascale[j];
	    if (ar == 0.) {
		ar = 1.;
	    }
	    if (amin > ar) {
		amin = ar;
		imin = i__;
	    }
	    if (amax < ar) {
		amax = ar;
		imax = i__;
	    }
	    if (ar > close1 && ar < close2) {
		++nclose;
/* ----           Ar     =  one */
	    }
	}
	ascale[j] = ar;
    }
    if (imin == 0) {
	amin = 0.;
	amax = 0.;
    }
    if (prnt10) {
	s_wsfi(&io___198);
	do_fio(&c__1, "Row", (ftnlen)3);
	do_fio(&c__1, (char *)&imin, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&amin, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "Row", (ftnlen)3);
	do_fio(&c__1, (char *)&imax, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&amax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&nclose, (ftnlen)sizeof(integer));
	d__1 = nclose * 100. / *m;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    }
/*     ------------------------------------------------------------------ */
/*     Phase 3. */
/*     Compute what is effectively the rhs for the constraints. */
/*     We set  rmax  =  ( A  -I )*x  for fixed columns and slacks, */
/*     including positive lower bounds and negative upper bounds. */
/*     ------------------------------------------------------------------ */
    if (! lonly) {
	dload_(m, &c_b94, &rmax[1], &c__1);
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    bnd = 0.;
	    b1 = bl[j];
	    b2 = bu[j];
	    if (b1 == b2) {
		bnd = b1;
	    }
	    if (b1 > 0.) {
		bnd = b1;
	    }
	    if (b2 < 0.) {
		bnd = b2;
	    }
	    if (bnd != 0.) {
		if (j <= *n) {
		    i__2 = loca[j + 1] - 1;
		    for (l = loca[j]; l <= i__2; ++l) {
			i__ = inda[l];
			rmax[i__] += acol[l] * bnd;
		    }
		} else {
/*                 Free slacks never get here, no need to skip them. */
		    i__ = j - *n;
		    rmax[i__] -= bnd;
		}
	    }
	}
/*        We don't want nonzeros in free rows to interfere. */
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (hrtype[i__] == 2) {
		rmax[i__] = 0.;
	    }
	}
/*        Scale rmax = rmax / (row scales),  and use its norm sigma */
/*        to adjust both row and column scales. */
	ac = dnormi_(m, &rmax[1], &c__1);
	dddiv_(m, &ascale[*n + 1], &c__1, &rmax[1], &c__1);
	sigma = dnormi_(m, &rmax[1], &c__1);
	if (prnt10) {
	    s_wsfi(&io___203);
	    do_fio(&c__1, (char *)&ac, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    s_wsfi(&io___204);
	    do_fio(&c__1, (char *)&sigma, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	sigma = max(sigma,1.);
	dscal_(nb, &sigma, &ascale[1], &c__1);
/*        Big scales might lead to excessive infeasibility when the */
/*        problem is unscaled.  If any are too big, scale them down. */
	amax = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = amax, d__2 = ascale[j];
	    amax = max(d__1,d__2);
	}
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (hrtype[i__] != 2) {
/* Computing MAX */
		d__1 = amax, d__2 = ascale[*n + i__];
		amax = max(d__1,d__2);
	    }
	}
	big = .1 / tolx;
	sigma = big / amax;
	if (sigma < 1.) {
	    dscal_(nb, &sigma, &ascale[1], &c__1);
	    if (prnt10) {
		s_wsfi(&io___206);
		do_fio(&c__1, (char *)&sigma, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
    if (lprscl > 0) {
	snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
	snprnt_(&c__1, " Row scales  r(i)      a(i,j)  =  r(i) * scaled a(i,"
		"j) / c(j)", &iw[1], leniw, (ftnlen)61);
	snprnt_(&c__1, " ----------------", &iw[1], leniw, (ftnlen)17);
	i__1 = *m;
	for (i__ = istart; i__ <= i__1; ++i__) {
	    s_wsfi(&io___207);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ascale[*n + i__], (ftnlen)sizeof(
		    doublereal));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
	snprnt_(&c__1, " Column scales  c(j)      x(j)    =  c(j) * scaled x"
		"(j)", &iw[1], leniw, (ftnlen)55);
	snprnt_(&c__1, " -------------------", &iw[1], leniw, (ftnlen)20);
	i__1 = *n;
	for (j = jstart; j <= i__1; ++j) {
	    s_wsfi(&io___208);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ascale[j], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    return 0;
} /* s2scal_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Scal */
/* Subroutine */ int s2scla_(integer *task, integer *m, integer *n, integer *
	nb, integer *iobj, doublereal *infbnd, doublereal *sclobj, integer *
	ne, integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	doublereal *ascale, doublereal *bl, doublereal *bu, doublereal *pi, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dddiv_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal bplus, cscale;

/*     ================================================================== */
/*     s2SclA scales or unscales A, bl, bu, pi, x using row and column */
/*     scales Sr and Sc stored as  Ascale = ( Sc  Sr ). */

/*       A(scaled)  = inv(Sr)*A*Sc */
/*       x(scaled)  = inv(Sc)*x,    s(scaled) = inv(Sr)*s, */
/*       pi(scaled) =     Sr *pi. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2scla. */
/*     24 Mar 2000: Current version of s2SclA */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --x;
    --bu;
    --bl;
    --ascale;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    bplus = *infbnd * .1;
    if (*task == 0) {
/*        --------------------------------------------------------------- */
/*        Scale A, bl, bu, x and pi. */
/*        --------------------------------------------------------------- */
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    cscale = ascale[j];
	    if (j <= *n) {
		i__2 = loca[j + 1] - 1;
		for (l = loca[j]; l <= i__2; ++l) {
		    i__ = inda[l];
		    acol[l] *= cscale / ascale[*n + i__];
		}
	    }
	    x[j] /= cscale;
	    if (bl[j] > -bplus) {
		bl[j] /= cscale;
	    }
	    if (bu[j] < bplus) {
		bu[j] /= cscale;
	    }
	}
	ddscl_(m, &ascale[*n + 1], &c__1, &pi[1], &c__1);
	if (*iobj > 0) {
	    *sclobj = ascale[*n + *iobj];
	}
    } else if (*task == 1) {
/*        --------------------------------------------------------------- */
/*        Unscale everything. */
/*        --------------------------------------------------------------- */
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    cscale = ascale[j];
	    if (j <= *n) {
		i__2 = loca[j + 1] - 1;
		for (l = loca[j]; l <= i__2; ++l) {
		    i__ = inda[l];
		    acol[l] *= ascale[*n + i__] / cscale;
		}
	    }
	    x[j] *= cscale;
	    if (bl[j] > -bplus) {
		bl[j] *= cscale;
	    }
	    if (bu[j] < bplus) {
		bu[j] *= cscale;
	    }
	}
	dddiv_(m, &ascale[*n + 1], &c__1, &pi[1], &c__1);
	*sclobj = 1.;
    }
    return 0;
} /* s2scla_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2SclA */
/* Subroutine */ int s2unpk_(integer *jq, integer *m, integer *n, integer *ne,
	 doublereal *colnrm, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, doublereal *y)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, l;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer islack;

/*     ================================================================== */
/*     s2unpk  expands the jq-th column of  ( A  -I )  into  y. */

/*     21 Jun 2004: Norm of unpacked column computed. */
/*     21 Jun 2004: Current version of s2unpk */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    dload_(m, &c_b94, &y[1], &c__1);
    *colnrm = 1.;
    if (*jq <= *n) {
	i__1 = loca[*jq + 1] - 1;
	for (l = loca[*jq]; l <= i__1; ++l) {
	    i__ = inda[l];
	    y[i__] = acol[l];
/* Computing MAX */
	    d__2 = *colnrm, d__3 = (d__1 = y[i__], abs(d__1));
	    *colnrm = max(d__2,d__3);
	}
    } else {
	islack = *jq - *n;
	y[islack] = -1.;
    }
    return 0;
} /* s2unpk_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2unpk */
/* Subroutine */ int s2vmax_(integer *n, integer *nncon, integer *maxvi, 
	doublereal *vimax, doublereal *bl, doublereal *bu, doublereal *fx)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal viol, slacki;

/*     ================================================================== */
/*     s2vmax  finds the largest nonlinear constraint violation. */
/*     maxvi   points to the biggest violation in Fx. */
/*     vimax   is the biggest violation in Fx. */

/*     13 Apr 1999: First version of s2viol. */
/*     21 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     See how much  Fx  violates the bounds on the nonlinear slacks. */
    /* Parameter adjustments */
    --fx;
    --bu;
    --bl;

    /* Function Body */
    *vimax = 0.;
    *maxvi = 1;
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	slacki = fx[i__];
/* Computing MAX */
	d__1 = 0., d__2 = bl[j] - slacki, d__1 = max(d__1,d__2), d__2 = 
		slacki - bu[j];
	viol = max(d__1,d__2);
	if (*vimax < viol) {
	    *vimax = viol;
	    *maxvi = i__;
	}
    }
    return 0;
} /* s2vmax_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2vmax */
/* Subroutine */ int s2dmat_(integer *matfil, integer *n, integer *nb, 
	integer *ne, integer *nloca, doublereal *acol, integer *inda, integer 
	*loca, integer *hs)
{
    /* Format strings */
    static char fmt_1000[] = "(1p,i10,i10,e24.14)";

    /* System generated locals */
    integer i__1, i__2;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer f_rew(alist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, k, l, ls;
    static doublereal aik;
    static integer ncol;

    /* Fortran I/O blocks */
    static cilist io___228 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_1000, 0 };


/*     ================================================================== */
/*     s2dmat  outputs A or B or (B S) to file matfil. */
/*     A triple (i, j, aij) is output for each nonzero entry, */
/*     intended for input to Matlab. */

/*     matfil    Matrix    hs(j)    j */
/*       91        A       >= 0    1:n */
/*       92        B       >= 3    1:nb */
/*       93       (B S)    >= 2    1:nb */

/*     10 Oct 2000: First version of m2xmat in MINOS, */
/*                  intended for output at the end of a successful run. */
/*                  Option "Report file 91" etc sets ireprt. */
/*                  m4savb calls m2xmat( ireprt, ... ) */
/*                  if ireprt = {91,92,93}. */

/*     09 Jun 2013: First version of s2dmat, intended for output before */
/*                  an unsuccessful B = LU factorization (Report file 92). */
/*                  If multiple LU factorizations occur, each B will be */
/*                  written to matfil, overwriting any previous B. */
/*                  Only the last B will be saved before the last LU. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hs;
    --inda;
    --acol;
    --loca;

    /* Function Body */
    if (*matfil == 91) {
	ls = 0;
	ncol = *n;
    } else if (*matfil == 92) {
	ls = 3;
	ncol = *nb;
    } else if (*matfil == 93) {
	ls = 2;
	ncol = *nb;
    } else {
	return 0;
    }
/* Output certain columns of A. */
    al__1.aerr = 0;
    al__1.aunit = *matfil;
    f_rew(&al__1);
    k = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (hs[j] >= ls) {
	    ++k;
	    i__2 = loca[j + 1] - 1;
	    for (l = loca[j]; l <= i__2; ++l) {
		i__ = inda[l];
		aik = acol[l];
		if (aik != 0.) {
		    io___228.ciunit = *matfil;
		    s_wsfe(&io___228);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&aik, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }
/* Output certain columns of -I. */
    aik = -1.;
    i__1 = ncol;
    for (j = *n + 1; j <= i__1; ++j) {
	if (hs[j] >= ls) {
	    ++k;
	    i__ = j - *n;
	    io___229.ciunit = *matfil;
	    s_wsfe(&io___229);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&aik, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = *matfil;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* s2dmat_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2dmat */
/* Subroutine */ int s2xmat_(integer *matfil, integer *n, integer *nb, 
	integer *ne, integer *nloca, integer *loca, integer *inda, doublereal 
	*acol, integer *hs)
{
    /* Format strings */
    static char fmt_1000[] = "(1p,i10,i10,e24.14)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, l, js, ls;
    static doublereal aik;
    static integer ncol;

    /* Fortran I/O blocks */
    static cilist io___238 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___239 = { 0, 0, 0, fmt_1000, 0 };


/*     ================================================================== */
/*     s2xmat  outputs A or B or (B S) to file matfil. */
/*     A triple (i, j, aij) is output for each nonzero entry, */
/*     intended for input to Matlab. */

/*     matfil    Matrix    hs(j)    j */
/*       91        A       >= 0    1:n */
/*       92        B       >= 3    1:nb */
/*       93       (B S)    >= 2    1:nb */

/*     10 Oct 2000: First version of s2xmat. */
/*     26 Jul 2003: Allow for hs(*) having "internal" values. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --hs;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    if (*matfil == 91) {
	ls = 0;
	ncol = *n;
    } else if (*matfil == 92) {
	ls = 3;
	ncol = *nb;
    } else if (*matfil == 93) {
	ls = 2;
	ncol = *nb;
    } else {
	return 0;
    }
/* Output certain columns of A. */
/* Treat internal values hs(j) = -1 or 4 as 0 (nonbasic). */
    k = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	js = hs[j];
	if (js < 0 || js == 4) {
	    js = 0;
	}
	if (js >= ls) {
	    ++k;
	    i__2 = loca[j + 1] - 1;
	    for (l = loca[j]; l <= i__2; ++l) {
		i__ = inda[l];
		aik = acol[l];
		if (aik != 0.) {
		    io___238.ciunit = *matfil;
		    s_wsfe(&io___238);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&aik, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }
/* Output certain columns of -I. */
    aik = -1.;
    i__1 = ncol;
    for (j = *n + 1; j <= i__1; ++j) {
	js = hs[j];
	if (js < 0 || js == 4) {
	    js = 0;
	}
	if (js >= ls) {
	    ++k;
	    i__ = j - *n;
	    io___239.ciunit = *matfil;
	    s_wsfe(&io___239);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&aik, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* s2xmat_ */

