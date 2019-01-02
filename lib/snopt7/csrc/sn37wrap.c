/* ./src/sn37wrap.f -- translated by f2c (version 20090411).
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
static integer c__13 = 13;
static integer c__3 = 3;
static integer c__11 = 11;
static integer c__22 = 22;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b367 = 0.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn37wrap.f --- auxiliaries for SNOPT wrappers. */

/*     SNOPTA  Arbitrary order for variables and constraints */
/*     s3argA   s3dfltA  s3prtA   s3sizA   s3bldA   s3inA    s3outA */
/*     s3prmA   s3mapA */

/*     SNOPTB  Basic format */
/*     s3argB   s3prtB */

/*     SQOPT   QP wrapper */
/*     s3argQ   s3prtQ */

/*     NPOPT   auxiliaries for NPSOL-like calls. */
/*     s3argN   s3HesN   s3inN   s3outN   s3prtN */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s3arga_(integer *iexit, integer *start, integer *nf, 
	integer *n, integer *ns, integer *nxname, integer *nfname, integer *
	objrow, integer *nea, integer *neg, doublereal *xlow, doublereal *
	xupp, char *xnames, doublereal *flow, doublereal *fupp, char *fnames, 
	integer *xstate, doublereal *xmul, integer *fstate, doublereal *fmul, 
	integer *lvlsrt, integer *errors, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen xnames_len, ftnlen fnames_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 XXX Start parameter not recognized:  "
	    "\002,i5)";
    static char fmt_1100[] = "(\002 XXX  Argument out of range:  \002,a6,"
	    "\002 = \002,i6)";
    static char fmt_1200[] = "(\002 XXX  Invalid argument xstate: \002,i6"
	    ",\002 elements modified to be in range.\002)";
    static char fmt_1210[] = "(\002 XXX  Invalid argument Fstate: \002,i6"
	    ",\002 elements modified to be in range.\002)";
    static char fmt_1310[] = "(\002 XXX  The equal bounds  \002,a8,\002  are"
	    " infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16.7)";
    static char fmt_1315[] = "(\002 XXX  The bounds on \002,a8,\002  are inc"
	    "onsistent.  xlow =\002,g16.7,\002  xupp =\002,g16.7)";
    static char fmt_1300[] = "(\002 XXX  The equal bounds on variable \002,i"
	    "6,\002  are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g1"
	    "6.7)";
    static char fmt_1305[] = "(\002 XXX  The bounds on variable \002,i6,\002"
	    "  are inconsistent.  xlow =\002,g16.7,\002  xupp =\002,g16.7)";
    static char fmt_1600[] = "(\002 XXX  Invalid arguments xlow, xupp: \002,"
	    "i6,\002 inconsistent bounds or infinite equal bounds.\002)";
    static char fmt_1510[] = "(\002 XXX  The equal bounds on  \002,a8,\002  "
	    "are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16.7)";
    static char fmt_1515[] = "(\002 XXX  The bounds on  \002,a8,\002  are in"
	    "consistent.  Flow =\002,g16.7,\002  Fupp =\002,g16.7)";
    static char fmt_1610[] = "(\002 XXX  Invalid arguments Flow, Fupp: \002,"
	    "i6,\002 inconsistent bounds or infinite equal bounds.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k;
    static doublereal b1, b2;
    static integer n1;
    static logical ok;
    static integer is, js;
    static char str[132];
    static integer errs;
    static logical named;
    static integer fmods, xmods;
    static doublereal infbnd;
    static integer argerr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___7 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___8 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___9 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___12 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___13 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___14 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___15 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___16 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___17 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___22 = { 0, str, 0, fmt_1200, 132, 1 };
    static icilist io___23 = { 0, str, 0, fmt_1210, 132, 1 };
    static icilist io___28 = { 0, str, 0, fmt_1310, 132, 1 };
    static icilist io___29 = { 0, str, 0, fmt_1315, 132, 1 };
    static icilist io___30 = { 0, str, 0, fmt_1300, 132, 1 };
    static icilist io___31 = { 0, str, 0, fmt_1305, 132, 1 };
    static icilist io___32 = { 0, str, 0, fmt_1600, 132, 1 };
    static icilist io___33 = { 0, str, 0, fmt_1510, 132, 1 };
    static icilist io___34 = { 0, str, 0, fmt_1515, 132, 1 };
    static icilist io___35 = { 0, str, 0, fmt_1510, 132, 1 };
    static icilist io___36 = { 0, str, 0, fmt_1515, 132, 1 };
    static icilist io___37 = { 0, str, 0, fmt_1610, 132, 1 };


/* ================================================================= */
/* s3argA   checks the arguments for snOptA. */

/* On exit, Errors says how many errors were encountered. */
/* lvlSrt is an integer version of Start: */
/*    Start   lvlSrt */
/*    'Cold'       0 */
/*    'Basis'      1 */
/*    'Warm'       2 */
/*    'Hot'        3 */
/* iw(gotFac,gotHes,gotScl) are set to 0 or 1. */
/*      33     Keep Factors of basis (LU) */
/*     303     Keep Hessian */
/*    3003     Keep Scales */
/*     333     Keep Factors and Hessian */
/*    3333     etc. */
/*       .     . */
/*       3     is treated as 3333 */

/* 21 Dec 2002: First version of s3argA. */
/* 03 Aug 2003: snPRNT adopted. */
/* 13 Mar 2004: Start = 3111 options decoded. */
/* 03 Apr 2004: Prevent n*nF overflow. */
/* 02 Sep 2004: Length of str increased to 132. */
/* 04 Dec 2004: Current version of s3argA. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --fmul;
    --fstate;
    --fupp;
    --flow;
    --xmul;
    --xstate;
    --xupp;
    --xlow;
    xnames -= 8;
    fnames -= 8;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    argerr = iw[106];
/* maximum # errors in MPS data */
    *iexit = 0;
    *errors = 0;
/* The options haven't been checked yet. */
    if (infbnd < 0.) {
	infbnd = 1e20;
    }
    if (argerr < 0) {
	argerr = 20;
    }
/* Print 20 lines max */
    named = *nxname == *n;
/* ================================================================= */
/* Decode 'Start'. */
/* ================================================================= */
    iw[230] = 0;
    iw[231] = 0;
    iw[232] = 0;
/* Determine the type of start. */
/* Optional parameters take precedence. */
    if (*lvlsrt != -11111) {
/* lvlSrt was set as an option */
	if (*lvlsrt == 0 || *lvlsrt == 1 || *lvlsrt == 2) {
/* Relax */
	} else if (*lvlsrt == 3) {
	    iw[230] = 1;
	    iw[231] = 1;
	    iw[232] = 1;
	} else {
/* lvlSrt is an unrecognized option */
	    *lvlsrt = -11111;
	}
    }
    if (*lvlsrt == -11111) {
/* lvlSrt is unset */
	if (*start == 0 || *start == 1 || *start == 2) {
	    *lvlsrt = *start;
	} else {
	    js = *start;
	    j = js - js / 10 * 10;
	    js -= j;
	    if (j == 3) {
		*lvlsrt = j;
		iw[230] = 1;
		iw[231] = 1;
		iw[232] = 1;
		j = js - js / 100 * 100;
		js -= j;
		if (j == 0) {
		    iw[232] = 0;
		}
		j = js - js / 10 * 10;
		js -= j;
		if (j == 0) {
		    iw[231] = 0;
		}
		if (js == 0) {
		    iw[230] = 0;
		}
	    } else {
		*lvlsrt = 0;
		s_wsfi(&io___7);
		do_fio(&c__1, (char *)&(*start), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	    }
	}
    }
/* ================================================================= */
/* Check the other arguments. */
/* ================================================================= */
    if (*nf < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___8);
	    do_fio(&c__1, "nF    ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*n < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___9);
	    do_fio(&c__1, "n     ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/* if (nS .lt. 0) then */
/*   Errors = Errors + 1 */
/*   if (Errors .le. argerr) then */
/*      write(str, 1100) 'nS    ', nS */
/*      call snPRNT( 13, str, iw, leniw ) */
/*   end if */
/* end if */
/* Test if neA or neG are bigger than n*nF. Avoid integer overflow. */
    n1 = max(*n,1);
    k = *nea % n1;
    if (*nea < 0 || k == 0 && *nea / n1 > *nf || k > 0 && *nea / n1 >= *nf) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___12);
	    do_fio(&c__1, "neA   ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nea), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    k = *neg % n1;
    if (*neg < 0 || k == 0 && *neg / n1 > *nf || k > 0 && *neg / n1 >= *nf) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___13);
	    do_fio(&c__1, "neG   ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*neg), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nxname != 1 && *nxname != *n) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___14);
	    do_fio(&c__1, "nxname ", (ftnlen)7);
	    do_fio(&c__1, (char *)&(*nxname), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nfname != 1 && *nfname != *nf) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___15);
	    do_fio(&c__1, "nFname ", (ftnlen)7);
	    do_fio(&c__1, (char *)&(*nfname), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nxname == 1 && *nfname != 1 || *nxname != 1 && *nfname == 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___16);
	    do_fio(&c__1, "nxname ", (ftnlen)7);
	    do_fio(&c__1, (char *)&(*nxname), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*objrow < 0 || *objrow > *nf) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___17);
	    do_fio(&c__1, "ObjRow", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*objrow), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    xmods = 0;
    fmods = 0;
    if (*lvlsrt == 0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    js = xstate[j];
	    if (js < -1 || js > 5) {
		++xmods;
		xstate[j] = 0;
	    }
	}
    } else if (*lvlsrt == 2 || *lvlsrt == 3) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    js = xstate[j];
	    if (js < 0 || js > 3) {
		++xmods;
		xstate[j] = 0;
	    }
	}
	i__1 = *nf;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is = fstate[i__];
	    if (is < 0 || is > 3) {
		++fmods;
		fstate[i__] = 0;
	    }
	}
    }
    if (xmods > 0) {
	s_wsfi(&io___22);
	do_fio(&c__1, (char *)&xmods, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    if (fmods > 0) {
	s_wsfi(&io___23);
	do_fio(&c__1, (char *)&fmods, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------------------------------------------------- */
/* Check the bounds on all variables and constraints. */
/* ----------------------------------------------------------------- */
    errs = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	b1 = xlow[j];
	b2 = xupp[j];
	ok = b1 <= b2 && b1 < infbnd && b2 > -infbnd;
	if (! ok) {
	    ++errs;
	    ++(*errors);
	    if (*errors <= argerr) {
		if (named) {
		    if (b1 == b2) {
			s_wsfi(&io___28);
			do_fio(&c__1, xnames + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___29);
			do_fio(&c__1, xnames + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (b1 == b2) {
			s_wsfi(&io___30);
			do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___31);
			do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	}
    }
    if (errs > 0) {
	s_wsfi(&io___32);
	do_fio(&c__1, (char *)&errs, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    errs = 0;
    i__1 = *nf;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b1 = flow[i__];
	b2 = fupp[i__];
	ok = b1 <= b2 && b1 < infbnd && b2 > -infbnd;
	if (! ok) {
	    ++errs;
	    ++(*errors);
	    if (*errors <= argerr) {
		if (named) {
		    if (b1 == b2) {
			s_wsfi(&io___33);
			do_fio(&c__1, fnames + (i__ << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___34);
			do_fio(&c__1, fnames + (i__ << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (b1 == b2) {
			s_wsfi(&io___35);
			do_fio(&c__1, "F(i) ", (ftnlen)5);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___36);
			do_fio(&c__1, "F(i) ", (ftnlen)5);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	}
    }
    if (errs > 0) {
	s_wsfi(&io___37);
	do_fio(&c__1, (char *)&errs, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    if (*errors >= argerr) {
	snprnt_(&c__3, " and so on ...", &iw[1], leniw, (ftnlen)14);
    }
    if (*errors > 0) {
/* ---------------------------- */
/* Invalid arguments. */
/* ---------------------------- */
	*iexit = 91;
/* invalid input argument */
    }
    return 0;
/* L1500: */
/* L1505: */
} /* s3arga_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3argA */
/* Subroutine */ int s3dflta_(integer *n, integer *lvlhes, integer *maxr, 
	integer *maxs, integer *mqnmod, integer *iw, integer *leniw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer never, kreset, qpslvr;

/* ================================================================= */
/* s3dfltA checks the optional parameters that depend on dimensions */
/* that are not yet known. */

/* The parameters will be re-checked later. */

/* See  snworkspace.info  for full da description of iw. */

/* 24 Nov 2007: first version of s3dfltA. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    *maxr = iw[52];
/* max columns of R. */
    *maxs = iw[53];
/* max # of superbasics */
    *mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    qpslvr = iw[55];
/* 0(1) => QP(QN) QP solver */
    kreset = iw[64];
/* Hessian frequency */
    *lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    never = 99999999;
    if (kreset <= 0) {
	kreset = never;
    }
/* Check superbasics limit maxS and Hessian dimension maxR. */
    if (*maxr < 0) {
	*maxr = min(2000,*n);
    }
    if (*maxs < 0) {
	*maxs = *n;
    }
/* Computing MAX */
    i__1 = min(*maxr,*n);
    *maxr = max(i__1,0);
/* Computing MAX */
    i__1 = min(*maxs,*n);
    *maxs = max(i__1,1);
    if (qpslvr < 0) {
	qpslvr = 0;
    }
    if (*maxr == 0) {
	qpslvr = 1;
    }
    if (qpslvr == 1 || *maxr < *maxs) {
	if (*lvlhes < 0) {
	    *lvlhes = 0;
	}
	if (*mqnmod < 0) {
	    *mqnmod = 10;
	}
    } else {
	if (*lvlhes < 0 && *n > 75) {
	    *lvlhes = 0;
	}
	if (*lvlhes < 0 && *n <= 75) {
	    *lvlhes = 1;
	}
	if (*lvlhes == 1) {
	    *mqnmod = kreset;
	}
	if (*mqnmod < 0) {
	    *mqnmod = 10;
	}
    }
    return 0;
} /* s3dflta_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3dfltA */
/* Subroutine */ int s3prta_(integer *m, integer *n, integer *nncon, integer *
	nnjac, integer *nnobj, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Initialized data */

    static char hestyp[24*3] = " Limited-Memory Hessian." " Full-Memory Hess"
	    "ian...." " Exact Hessian..........";
    static char lsrch[24*2] = " Nonderiv.  linesearch.." " Derivative linese"
	    "arch..";
    static char pivtyp[24*4] = " LU partial  pivoting..." " LU rook     pivo"
	    "ting..." " LU complete pivoting..." " LU diagonal pivoting...";
    static char prbtyp[24*3] = " Maximize..............." " Feasible point o"
	    "nly...." " Minimize...............";
    static char qptype[24*3] = " QPsolver Cholesky......" " QPsolver CG....."
	    "......." " QPsolver QN............";
    static char srttyp[24*4] = " Cold start............." " Basis file......"
	    "......." " Warm start............." " Hot start..............";
    static char noyes[3*2] = " No" "Yes";

    /* Format strings */
    static char fmt_2110[] = "(\002 Solution file..........\002,i10,6x,\002 "
	    "Old basis file ........\002,i10,6x,\002 Standard input........"
	    ".\002,i10)";
    static char fmt_2120[] = "(\002 Insert file............\002,i10,6x,\002 "
	    "New basis file ........\002,i10,6x,\002 (Printer)............."
	    ".\002,i10)";
    static char fmt_2130[] = "(\002 Punch file.............\002,i10,6x,\002 "
	    "Backup basis file......\002,i10,6x,\002 (Specs file).........."
	    ".\002,i10)";
    static char fmt_2140[] = "(\002 Load file..............\002,i10,6x,\002 "
	    "Dump file..............\002,i10,6x,\002 Standard output......."
	    ".\002,i10)";
    static char fmt_2210[] = "(\002 Print frequency........\002,i10,6x,\002 "
	    "Check frequency........\002,i10,6x,\002 Save new basis map...."
	    ".\002,i10)";
    static char fmt_2220[] = "(\002 Summary frequency......\002,i10,6x,\002 "
	    "Factorization frequency\002,i10,6x,\002 Expand frequency......"
	    ".\002,i10)";
    static char fmt_2310[] = "(a24)";
    static char fmt_2320[] = "(\002 Scale tolerance........\002,0p,f10.3,6x"
	    ",\002 Minor feasibility tol..\002,1p,e10.2,6x,\002 Iteration lim"
	    "it........\002,i10)";
    static char fmt_2330[] = "(\002 Scale option...........\002,i10,6x,\002 "
	    "Minor optimality  tol..\002,1p,e10.2,6x,\002 Minor print level.."
	    "....\002,i10)";
    static char fmt_2340[] = "(\002 Crash tolerance........\002,0p,f10.3,6x"
	    ",\002 Pivot tolerance........\002,1p,e10.2,6x,\002 Partial price"
	    "..........\002,i10)";
    static char fmt_2350[] = "(\002 Crash option...........\002,i10,6x,\002 "
	    "Elastic weight.........\002,1p,e10.2,6x,\002 Prtl price section "
	    "( A)\002,i10)";
    static char fmt_2360[] = "(40x,\002 New superbasics........\002,i10,6x"
	    ",\002 Prtl price section (-I)\002,i10)";
    static char fmt_2380[] = "(\002 Subspace tolerance.....\002,0p,f10.5,6x"
	    ",\002 CG tolerance...........\002,1p,e10.2,6x,\002 CG Iterations"
	    "..........\002,i10)";
    static char fmt_2390[] = "(80x,\002 CG preconditioning.....\002,i10)";
    static char fmt_2410[] = "(a24,16x,a24,16x,\002 Proximal Point method."
	    ".\002,i10)";
    static char fmt_2420[] = "(\002 Nonlinear objectiv vars\002,i10,6x,\002 "
	    "Objective Row..........\002,i10,6x,\002 Function precision...."
	    ".\002,1p,e10.2)";
    static char fmt_2430[] = "(\002 Unbounded step size....\002,1p,e10.2,6x"
	    ",\002 Superbasics limit......\002,i10,6x,\002 Difference interva"
	    "l....\002,1p,e10.2)";
    static char fmt_2440[] = "(\002 Unbounded objective....\002,1p,e10.2,6x"
	    ",\002 Reduced Hessian dim....\002,i10,6x,\002 Central difference"
	    " int.\002,1p,e10.2)";
    static char fmt_2450[] = "(\002 Major step limit.......\002,1p,e10.2,6x,"
	    "a11,\002 linesearch..\002,16x,\002 Derivative option......\002,i"
	    "10)";
    static char fmt_2460[] = "(\002 Major iterations limit.\002,i10,6x,\002 "
	    "Linesearch tolerance...\002,0p,f10.5,6x,\002 Verify level......."
	    "....\002,i10)";
    static char fmt_2470[] = "(\002 Minor iterations limit.\002,i10,6x,\002 "
	    "Penalty parameter......\002,1p,e10.2,6x,\002 Major Print Level.."
	    "....\002,i10)";
    static char fmt_2480[] = "(40x,\002 Major optimality tol...\002,1p,e10.2)"
	    ;
    static char fmt_2510[] = "(a24,16x,\002 Hessian updates........\002,i10,"
	    "6x,\002 Hessian frequency......\002,i10)";
    static char fmt_2520[] = "(80x,\002 Hessian flush..........\002,i10)";
    static char fmt_2610[] = "(\002 Nonlinear constraints..\002,i10,6x,\002 "
	    "Major feasibility tol..\002,1p,e10.2,6x,\002 Violation limit...."
	    "....\002,e10.2)";
    static char fmt_2620[] = "(\002 Nonlinear Jacobian vars\002,i10)";
    static char fmt_2710[] = "(\002 LU factor tolerance....\002,0p,f10.2,6x"
	    ",\002 LU singularity tol.....\002,1p,e10.2,6x,\002 Timing level."
	    "..........\002,i10)";
    static char fmt_2720[] = "(\002 LU update tolerance....\002,0p,f10.2,6x"
	    ",\002 LU swap tolerance......\002,1p,e10.2,6x,\002 Debug level.."
	    "..........\002,i10)";
    static char fmt_2730[] = "(a24,16x,\002 eps (machine precision)\002,1p,e"
	    "10.2,6x,\002 System information.....\002,7x,a3)";
    static char fmt_2740[] = "(80x,\002 Sticky parameters......\002,7x,a3)";
    static char fmt_3000[] = "(\002 Scale option\002,i3,\002,    Partial pri"
	    "ce\002,i4)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer nnl;
    static doublereal eps;
    static integer npr1, npr2;
    static char str1[132], str2[132], str3[132], str4[132], str5[132], str6[
	    132];
    static integer kfac, kchk, klog, ksav, maxr, maxs;
    static doublereal tolx, xpen0, utol1;
    static integer iback, ioldb;
    static doublereal bigdx, bigfx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal tolcg;
    static integer istdi;
    static doublereal xdlim;
    static integer idump, maxmn;
    static doublereal epsrf;
    static integer istdo;
    static doublereal vilim;
    static integer isoln, ksumm;
    static doublereal tolqp;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *);
    static doublereal fdint1, fdint2, wtinf0;
    static integer iloadb, kdegen;
    static doublereal tolfac;
    static integer icrash, lprdbg;
    static doublereal wolfeg, tcrash;
    static integer mmajor, ispecs, minprc, minmax, cgitmx, itnlim, deropt, 
	    kreset, lvlhes, lvlsch, lvlscl, mflush, mminor, lvlpre, iprint, 
	    mqnmod, lvltim, iinsrt, mnewsb, lvlppm, lvlver, lprprm, lvlpiv, 
	    mjrprt, nparpr, objrow;
    static doublereal scltol;
    static integer mnrprt;
    static doublereal tolcon, tolnlp, tolpiv;
    static integer lvlsrt, qpslvr;
    static doublereal tolswp;
    static integer stkyop;
    static doublereal tolupd;
    static integer lvlsys;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___127 = { 0, str1, 0, fmt_2110, 132, 1 };
    static icilist io___129 = { 0, str2, 0, fmt_2120, 132, 1 };
    static icilist io___131 = { 0, str3, 0, fmt_2130, 132, 1 };
    static icilist io___133 = { 0, str4, 0, fmt_2140, 132, 1 };
    static icilist io___134 = { 0, str1, 0, fmt_2210, 132, 1 };
    static icilist io___135 = { 0, str2, 0, fmt_2220, 132, 1 };
    static icilist io___136 = { 0, str1, 0, fmt_2310, 132, 1 };
    static icilist io___137 = { 0, str2, 0, fmt_2320, 132, 1 };
    static icilist io___138 = { 0, str3, 0, fmt_2330, 132, 1 };
    static icilist io___139 = { 0, str4, 0, fmt_2340, 132, 1 };
    static icilist io___141 = { 0, str5, 0, fmt_2350, 132, 1 };
    static icilist io___143 = { 0, str6, 0, fmt_2360, 132, 1 };
    static icilist io___144 = { 0, str1, 0, fmt_2380, 132, 1 };
    static icilist io___145 = { 0, str2, 0, fmt_2390, 132, 1 };
    static icilist io___146 = { 0, str1, 0, fmt_2410, 132, 1 };
    static icilist io___147 = { 0, str2, 0, fmt_2420, 132, 1 };
    static icilist io___148 = { 0, str3, 0, fmt_2430, 132, 1 };
    static icilist io___149 = { 0, str4, 0, fmt_2440, 132, 1 };
    static icilist io___150 = { 0, str1, 0, fmt_2450, 132, 1 };
    static icilist io___151 = { 0, str2, 0, fmt_2460, 132, 1 };
    static icilist io___152 = { 0, str3, 0, fmt_2470, 132, 1 };
    static icilist io___153 = { 0, str4, 0, fmt_2480, 132, 1 };
    static icilist io___154 = { 0, str1, 0, fmt_2510, 132, 1 };
    static icilist io___155 = { 0, str2, 0, fmt_2520, 132, 1 };
    static icilist io___156 = { 0, str1, 0, fmt_2610, 132, 1 };
    static icilist io___157 = { 0, str2, 0, fmt_2620, 132, 1 };
    static icilist io___158 = { 0, str1, 0, fmt_2710, 132, 1 };
    static icilist io___159 = { 0, str2, 0, fmt_2720, 132, 1 };
    static icilist io___160 = { 0, str3, 0, fmt_2730, 132, 1 };
    static icilist io___161 = { 0, str4, 0, fmt_2740, 132, 1 };
    static icilist io___162 = { 0, str1, 0, fmt_3000, 132, 1 };


/* ================================================================= */
/* s3prtA prints the settings of the optional parameters for the */
/* snoptA wrapper. */

/* See  snworkspace.info  for full documentation of cw, iw and rw. */

/* 15 Nov 1991: First version (s8prnt). */
/* 10 Dec 2002: More LU pivoting options. */
/* 03 Aug 2003: snPRNT adopted. */
/* 01 Sep 2007: Sticky parameters added. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
/* Set some local machine-dependent constants. */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
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
/* excludes small elements of y. */
    tcrash = rw[62];
/* crash tolerance. */
    tolswp = rw[65];
/* LU swap tolerance. */
    tolfac = rw[66];
/* LU factor tolerance. */
    tolupd = rw[67];
/* LU update tolerance. */
    bigfx = rw[71];
/* unbounded objective. */
    bigdx = rw[72];
/* unbounded step. */
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
    wtinf0 = rw[88];
/* infeasibility weight */
    xpen0 = rw[89];
/* initial penalty parameter. */
    scltol = rw[92];
/* scale tolerance. */
    utol1 = rw[154];
/* abs tol for small diag of U. */
    istdi = iw[9];
/* Standard Input */
    istdo = iw[10];
/* Standard Output */
    ispecs = iw[11];
/* Specs (options) file */
    iprint = iw[12];
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    qpslvr = iw[55];
/* 0(1/2) => QPsolver */
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
    lvlsrt = iw[69];
/* = 0:1:2:3 => cold:warm:basis:hot start */
    lvlsys = iw[71];
/* > 0   => print system info */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lvlscl = iw[75];
/* scale option */
    lvlsch = iw[76];
/* >0     => use derivatives in the line search */
    lvlpre = iw[77];
/* >0    => QN preconditioned CG */
    lvlver = iw[78];
/* Verify level */
    lvlppm = iw[79];
/* Proximal Point method for x0 */
    lvlpiv = iw[80];
/* 0(1) LU threshold partial(complete) pivoting */
    lprprm = iw[81];
/* > 0    => parms are printed */
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
    isoln = iw[131];
/* solution file */
    lvltim = iw[182];
/* ----------------------------------------------------------------- */
/* Timing level */
    if (iprint <= 0 || mjrprt == 0 || lprprm == 0) {
	return 0;
    }
    nnl = max(*nnjac,*nnobj);
    minprc = 10;
    npr1 = *n / nparpr;
    npr2 = *m / nparpr;
    if (max(npr1,npr2) < minprc) {
	maxmn = max(*m,*n);
	nparpr = maxmn / min(maxmn,minprc);
	npr1 = *n / nparpr;
	npr2 = *m / nparpr;
    }
/* ================================================================= */
/* Print parameters except if PRINT LEVEL = 0 */
/* or SUPPRESS PARAMETERS was specified. */
/* ================================================================= */
    s1page_(&c__1, &iw[1], leniw);
    snprnt_(&c__1, " Parameters", &iw[1], leniw, (ftnlen)11);
    snprnt_(&c__1, " ==========", &iw[1], leniw, (ftnlen)11);
/* -------------------- */
/* Files. */
/* -------------------- */
    snprnt_(&c__11, " Files", &iw[1], leniw, (ftnlen)6);
    snprnt_(&c__1, " -----", &iw[1], leniw, (ftnlen)6);
    s_wsfi(&io___127);
    do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ioldb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&istdi, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___129);
    do_fio(&c__1, (char *)&iinsrt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inewb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iprint, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___131);
    do_fio(&c__1, (char *)&ipnch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iback, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ispecs, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___133);
    do_fio(&c__1, (char *)&iloadb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&idump, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&istdo, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
/* -------------------- */
/* Frequencies. */
/* -------------------- */
    snprnt_(&c__11, " Frequencies", &iw[1], leniw, (ftnlen)12);
    snprnt_(&c__1, " -----------", &iw[1], leniw, (ftnlen)12);
    s_wsfi(&io___134);
    do_fio(&c__1, (char *)&klog, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kchk, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ksav, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___135);
    do_fio(&c__1, (char *)&ksumm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kfac, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kdegen, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
/* -------------------- */
/* QP subproblems. */
/* -------------------- */
    snprnt_(&c__11, " QP subproblems", &iw[1], leniw, (ftnlen)15);
    snprnt_(&c__1, " --------------", &iw[1], leniw, (ftnlen)15);
    s_wsfi(&io___136);
    do_fio(&c__1, qptype + qpslvr * 24, (ftnlen)24);
    e_wsfi();
    s_wsfi(&io___137);
    do_fio(&c__1, (char *)&scltol, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&itnlim, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___138);
    do_fio(&c__1, (char *)&lvlscl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&tolqp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&mnrprt, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___139);
    do_fio(&c__1, (char *)&tcrash, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolpiv, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nparpr, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___141);
    do_fio(&c__1, (char *)&icrash, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&wtinf0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&npr1, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___143);
    do_fio(&c__1, (char *)&mnewsb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npr2, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str5, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str6, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* CG QP solver. */
/* ----------------------- */
    if (qpslvr == 1 || maxr < maxs) {
	snprnt_(&c__11, " Conjugate-gradient QP solver", &iw[1], leniw, (
		ftnlen)29);
	snprnt_(&c__1, " ----------------------------", &iw[1], leniw, (
		ftnlen)29);
	s_wsfi(&io___144);
	do_fio(&c__1, (char *)&etarg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&tolcg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cgitmx, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___145);
	do_fio(&c__1, (char *)&lvlpre, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    }
/* -------------------- */
/* SQP method. */
/* -------------------- */
    snprnt_(&c__11, " The SQP Method", &iw[1], leniw, (ftnlen)15);
    snprnt_(&c__1, " --------------", &iw[1], leniw, (ftnlen)15);
    s_wsfi(&io___146);
    do_fio(&c__1, prbtyp + (minmax + 1) * 24, (ftnlen)24);
    do_fio(&c__1, srttyp + lvlsrt * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___147);
    do_fio(&c__1, (char *)&(*nnobj), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&objrow, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&epsrf, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___148);
    do_fio(&c__1, (char *)&bigdx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&maxs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&fdint1, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___149);
    do_fio(&c__1, (char *)&bigfx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&maxr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&fdint2, (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    s_wsfi(&io___150);
    do_fio(&c__1, (char *)&xdlim, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, lsrch + lvlsch * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&deropt, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___151);
    do_fio(&c__1, (char *)&mmajor, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&wolfeg, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lvlver, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___152);
    do_fio(&c__1, (char *)&mminor, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&xpen0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&mjrprt, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___153);
    do_fio(&c__1, (char *)&tolnlp, (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
/* -------------------- */
/* Hessian approximation. */
/* -------------------- */
    if (nnl > 0) {
	snprnt_(&c__11, " Hessian Approximation", &iw[1], leniw, (ftnlen)22);
	snprnt_(&c__1, " ---------------------", &iw[1], leniw, (ftnlen)22);
	s_wsfi(&io___154);
	do_fio(&c__1, hestyp + lvlhes * 24, (ftnlen)24);
	do_fio(&c__1, (char *)&mqnmod, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kreset, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___155);
	do_fio(&c__1, (char *)&mflush, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    }
/* -------------------- */
/* Nonlinear constraints. */
/* -------------------- */
    if (*nncon > 0) {
	snprnt_(&c__11, " Nonlinear constraints", &iw[1], leniw, (ftnlen)22);
	snprnt_(&c__1, " ---------------------", &iw[1], leniw, (ftnlen)22);
	s_wsfi(&io___156);
	do_fio(&c__1, (char *)&(*nncon), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&tolcon, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vilim, (ftnlen)sizeof(doublereal));
	e_wsfi();
	s_wsfi(&io___157);
	do_fio(&c__1, (char *)&(*nnjac), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    }
/* -------------------- */
/* Miscellaneous */
/* -------------------- */
    snprnt_(&c__11, " Miscellaneous", &iw[1], leniw, (ftnlen)14);
    snprnt_(&c__1, " -------------", &iw[1], leniw, (ftnlen)14);
    s_wsfi(&io___158);
    do_fio(&c__1, (char *)&tolfac, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&utol1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lvltim, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___159);
    do_fio(&c__1, (char *)&tolupd, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolswp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lprdbg, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___160);
    do_fio(&c__1, pivtyp + lvlpiv * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, noyes + lvlsys * 3, (ftnlen)3);
    e_wsfi();
    s_wsfi(&io___161);
    do_fio(&c__1, noyes + stkyop * 3, (ftnlen)3);
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    s_wsfi(&io___162);
    do_fio(&c__1, (char *)&lvlscl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nparpr, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__22, str1, &iw[1], leniw, (ftnlen)132);
    return 0;
} /* s3prta_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3prtA */
/* Subroutine */ int s3siza_(integer *iexit, integer *n, integer *nf, integer 
	*nkx, integer *objrow, integer *iafun, integer *javar, integer *lena, 
	integer *nea, integer *igfun, integer *jgvar, integer *leng, integer *
	neg, integer *m, integer *negcon, integer *ne, integer *nncon, 
	integer *nnjac, integer *nnobj, integer *iobj, integer *kx, integer *
	leniw, integer *iw)
{
    /* Format strings */
    static char fmt_2000[] = "(\002 XXX  Nonlinear derivative element    k"
	    " = \002,i6,\002,  row \002,i6,\002 column \002,i6,\002 is out of"
	    " range.\002)";
    static char fmt_2010[] = "(\002 XXX  Linear    derivative element    k"
	    " = \002,i6,\002,  row \002,i6,\002 column \002,i6,\002 is out of"
	    " range.\002)";
    static char fmt_1200[] = "(\002 ===>  WARNING - Column \002,i6,\002 of t"
	    "he Jacobian is empty.\002)";
    static char fmt_1300[] = "(\002 ===>  WARNING - Row    \002,i6,\002 of t"
	    "he Jacobian is empty.\002)";
    static char fmt_1000[] = "(\002 ===>  WARNING - too many rows.\002)";
    static char fmt_1500[] = "(\002 ===>  No objective row specified --"
	    "-\002,\002 finding a feasible point.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k, in, jn, js, col, nnl, row;
    static char str[132];
    static integer njac, nobj, ncol, nrow, jsec1, jsec2, jsec3, jsec4;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer nboth;
    static logical first;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___172 = { 0, str, 0, fmt_2000, 132, 1 };
    static icilist io___176 = { 0, str, 0, fmt_2010, 132, 1 };
    static icilist io___178 = { 0, str, 0, fmt_1200, 132, 1 };
    static icilist io___179 = { 0, str, 0, fmt_1300, 132, 1 };
    static icilist io___180 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___188 = { 0, str, 0, fmt_1500, 132, 1 };


/* ================================================================= */
/* s3sizA runs through the data once in order to compute the */
/* dimensions of the quantities: */
/* Jcol    is the      m x n      constraint Jacobian. */
/* gCon    is the  nnCon x nnJac  nonlinear part of the Jacobian. */
/* gObj    is the  nnObj x 1      objective gradient. */

/*         negCon,  ne,  nnCon, nnJac,  nnObj  and  iObj. */

/* On entry: */
/* --------- */
/* ObjRow  is the (natural) row number of the objective. */

/* nkx     (= n + nF) is the length of kx. */

/* On exit: */
/* -------- */
/* m       is the row dimension of Jcol.  m = nF. */

/* iObj    is the position of the linear objective row in J. */
/*         iObj = 0 means that there is no linear objective row. */

/* kx(nkx) is a list of the indices that defines the position of */
/*         each user-defined F(i) and x(j) in  Jcol, i.e., */
/*         k = kx(j)   => variable   j is the k-th column of Jcol. */
/*         l = kx(n+i) => constraint i is the l-th row    of Jcol. */

/* 08 Nov 1998: First version of s3sizA. */
/* 18 Jul 2000: First version for SNOPT 6.1. */
/* 03 Aug 2003: snPRNT adopted. */
/* 06 Apr 2008: ObjRow > 0 implies iObj > 0. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --kx;
    --javar;
    --iafun;
    --jgvar;
    --igfun;
    --iw;

    /* Function Body */
    *iexit = 0;
/* ================================================================= */
/* Count the number of nonlinear objective and Jacobian variables. */
/* put the nonlinear rows first, recording their positions in kx. */
/* Check if variable is nonlinear in the objective, the Jacobian, */
/* or both. */

/*  kx             type of variable */
/*  --             ---------------- */
/*  Linear         linear    Jacobian,   linear    Objective */
/*  Jac            nonlinear Jacobian,   linear    Objective */
/*  Obj            linear    Jacobian,   nonlinear Objective */
/*  Both           nonlinear Jacobian,   nonlinear Objective */
/* ================================================================= */
    nobj = 0;
    njac = 0;
    nboth = 0;
    iload_(nkx, &c_n1, &kx[1], &c__1);
    *iobj = 0;
    *nncon = 0;
    *negcon = 0;
    nrow = 0;
    ncol = 0;
/* ----------------------------------------------------------------- */
/* Process the nonlinear Jacobian elements.  The array G holds the */
/* nonzero nonlinear elements of the derivative matrix G(x) that */
/* includes the objective gradient (if there is one). */
/* kx(1:n) and kx(n+1:nkx) give the coordinates of J(i,j). */
/* ----------------------------------------------------------------- */
    i__1 = *neg;
    for (k = 1; k <= i__1; ++k) {
	in = igfun[k];
	jn = jgvar[k];
	if (jn < 0 || jn > *n || in < 0 || in > *nf) {
	    s_wsfi(&io___172);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&in, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	    *iexit = 91;
/* Invalid input argument */
	    goto L900;
	}
	col = jn;
	if (in == *objrow) {
/* ---------------------------- */
/* Element of gObj. */
/* ---------------------------- */
	    if (kx[col] == -1) {
		kx[col] = 2;
		++nobj;
		++ncol;
	    } else if (kx[col] == 1) {
		kx[col] = 3;
		++nobj;
		++nboth;
	    }
	} else {
/* -------------------------- */
/* Element of Jcol and gCon. */
/* -------------------------- */
	    ++(*negcon);
/* See if this element starts a new nonlinear column. */
	    if (kx[col] == -1) {
		kx[col] = 1;
		++njac;
		++ncol;
	    } else if (kx[col] == 2) {
		kx[col] = 3;
		++njac;
		++nboth;
	    }
/* See if this is the first element in a new nonlinear row. */
	    row = *n + in;
	    if (kx[row] == -1) {
		++(*nncon);
		++nrow;
		kx[row] = nrow;
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Now we know all the dimensions of the nonlinear part. */
/* Figure out the variable order in the first nnL columns. */
/* Set  nnObj,  nnCon,  nnJac  and  nnL. */
/* ----------------------------------------------------------------- */
    nnl = ncol;
    *nnjac = njac;
    if (nobj > nboth) {
	*nnobj = njac + (nobj - nboth);
    } else {
	*nnobj = nboth;
    }
/* ----------------------------------------------------------------- */
/* Process the constant part of the Jacobian. */
/* The linear rows follow the nonlinear rows. */
/* ----------------------------------------------------------------- */
    *ne = *negcon;
    i__1 = *nea;
    for (k = 1; k <= i__1; ++k) {
	in = iafun[k];
	jn = javar[k];
	if (jn < 0 || jn > *n || in < 0 || in > *nf) {
	    s_wsfi(&io___176);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&in, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	    *iexit = 91;
/* Invalid input argument */
	    goto L900;
	}
/* See if this element starts a new linear column. */
	col = jn;
	if (kx[col] == -1) {
	    ++ncol;
	    kx[col] = 0;
	}
	++(*ne);
/* See if this element starts a new linear row. */
/* All linear terms are held in J. */
	row = *n + in;
	if (kx[row] == -1) {
	    ++nrow;
	    kx[row] = nrow;
	    if (in == *objrow) {
		*iobj = nrow;
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* If necessary, add a dummy free linear row. */
/* ----------------------------------------------------------------- */
    if (*objrow > 0 && *iobj == 0) {
	++nrow;
	*iobj = nrow;
	row = *n + *objrow;
	kx[row] = *iobj;
    }
/* ------------------------ */
/* Check for empty columns. */
/* ------------------------ */
    first = TRUE_;
    i__1 = *n;
    for (jn = 1; jn <= i__1; ++jn) {
	if (kx[jn] == -1) {
	    if (first) {
		snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
		first = FALSE_;
	    }
	    s_wsfi(&io___178);
	    do_fio(&c__1, (char *)&jn, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	    ++ncol;
	    kx[jn] = 0;
	}
    }
/* --------------------- */
/* Check for empty rows. */
/* --------------------- */
    i__1 = *nf;
    for (in = 1; in <= i__1; ++in) {
	row = *n + in;
	if (kx[row] == -1) {
	    if (first) {
		snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
		first = FALSE_;
	    }
	    s_wsfi(&io___179);
	    do_fio(&c__1, (char *)&in, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
	    ++nrow;
	    kx[row] = nrow;
	}
    }
/* ----------------------------------------------------------------- */
/* Set the row dimension of Jcol. */
/* ----------------------------------------------------------------- */
    *m = *nf;
    if (nrow != *m) {
	if (first) {
	    snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
	    first = FALSE_;
	}
	s_wsfi(&io___180);
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
/* ================================================================= */
/*  Re-order the variables into four sections: */
/*     Section 1 - type  Both */
/*     Section 2 - type  Jac */
/*     Section 3 - type  Obj */
/*     Section 4 - type  Linear */

/*  jSecn+1  points to the next free slot in section n. */
/* ================================================================= */
    jsec1 = 0;
    jsec2 = nboth;
    jsec3 = *nnjac;
    jsec4 = nnl;
/* Loop through the columns and initialize kx, which holds */
/* the final column order. */
    i__1 = *n;
    for (jn = 1; jn <= i__1; ++jn) {
	if (kx[jn] == 3) {
	    ++jsec1;
	    js = jsec1;
	} else if (kx[jn] == 1) {
	    ++jsec2;
	    js = jsec2;
	} else if (kx[jn] == 2) {
	    ++jsec3;
	    js = jsec3;
	} else {
	    ++jsec4;
	    js = jsec4;
	}
	kx[jn] = js;
    }
/* ----------------------------------------------------------------- */
/* Check for constant elements in the nonlinear part of Jcol. */
/* ----------------------------------------------------------------- */
    i__1 = *nea;
    for (k = 1; k <= i__1; ++k) {
	in = iafun[k];
	jn = javar[k];
	row = *n + in;
	col = jn;
	i__ = kx[row];
	j = kx[col];
	if (i__ <= *nncon && j <= *nnjac) {
	    ++(*negcon);
	}
    }
    if (*ne == 0) {
/* no Jacobian elements */
/* Add a dummy row */
	++(*ne);
    }
    if (*objrow == 0) {
/* ------------------------------------------- */
/* No objective row.  SNOPT  will provide one. */
/* ------------------------------------------- */
	if (first) {
	    snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
	    first = FALSE_;
	}
	s_wsfi(&io___188);
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)132);
    }
L900:
    return 0;
} /* s3siza_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3sizA */
/* Subroutine */ int s3blda_(integer *objrow, integer *n, integer *nkx, 
	integer *nncon, integer *nnjac, integer *kx, integer *nglin, integer *
	iafun, integer *javar, integer *lena, integer *nea, doublereal *a, 
	integer *igfun, integer *jgvar, integer *leng, integer *neg, integer *
	ne, integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	integer *negcon, integer *nlocg, integer *locg, integer *indg, 
	integer *iy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, in, jn, nej, col, row;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer ncolj;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer nextg, nextj;

/* ================================================================= */
/* s3bldA builds the data structure (locJ, indJ, Jcol) for the full */
/* sparse Jacobian. */

/* On entry, */
/*  neG             holds the number of nozero nonlinear derivative */
/*                  elements. */
/*  G, iGfun, jGvar hold neG user-defined coordinates (Gij,i,j). */
/*                  neG  must be nonnegative. */
/*  A, iAfun, jAvar hold neA user-defined coordinates (Aij,i,j). */
/*                  neA  must be nonnegative. */

/* 25 Sep 1998: First version of s3bldA. */
/* 18 Jul 2000: First version for SNOPT 6.1. */
/* 06 Apr 2008: ObjRow > 0 implies iObj > 0. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* Count the number of entries in each column (again). */
/* Objective elements are excluded from  gCon and Jcol. */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --iy;
    --kx;
    --a;
    --javar;
    --iafun;
    --jgvar;
    --igfun;
    --jcol;
    --indj;
    --locj;
    --indg;
    --locg;
    --nglin;

    /* Function Body */
    iload_(n, &c__0, &locj[1], &c__1);
    if (*nnjac > 0) {
	i__1 = *neg;
	for (k = 1; k <= i__1; ++k) {
	    in = igfun[k];
	    jn = jgvar[k];
	    col = jn;
	    j = kx[col];
	    if (in != *objrow) {
		++locj[j];
	    }
	}
	icopy_(nnjac, &locj[1], &c__1, &locg[1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Now count the constant derivatives. */
/* nGlin(j) = the number of constant elements in column j of gCon. */
/* ----------------------------------------------------------------- */
    if (*nnjac > 0) {
	iload_(nnjac, &c__0, &nglin[1], &c__1);
	iload_(nnjac, &c__0, &iy[1], &c__1);
    }
    i__1 = *nea;
    for (k = 1; k <= i__1; ++k) {
	in = iafun[k];
	jn = javar[k];
	row = *n + in;
	col = jn;
	i__ = kx[row];
	j = kx[col];
	++locj[j];
	if (i__ <= *nncon && j <= *nnjac) {
	    ++locg[j];
	    ++nglin[j];
	}
    }
/* ----------------------------------------------------------------- */
/* Initialize the column pointers. */
/* ----------------------------------------------------------------- */
    locj[*nlocj] = *ne + 1;
    for (j = *n; j >= 1; --j) {
	locj[j] = locj[j + 1] - locj[j];
    }
    nej = 0;
    if (*negcon > 0) {
	locg[*nlocg] = *negcon + 1;
	for (i__ = *nnjac; i__ >= 1; --i__) {
	    locg[i__] = locg[i__ + 1] - locg[i__];
	}
/* -------------------------------------------------------------- */
/* Assign the nonlinear constraint derivatives to indJ and Jcol. */
/* locJ(j) points to the next nonzero in column i. */
/* -------------------------------------------------------------- */
	i__1 = *neg;
	for (k = 1; k <= i__1; ++k) {
	    in = igfun[k];
	    jn = jgvar[k];
	    col = jn;
	    j = kx[col];
	    if (in != *objrow) {
		row = *n + in;
		i__ = kx[row];
		nextj = locj[j];
		indj[nextj] = i__;
		jcol[nextj] = 0.;
		locj[j] = nextj + 1;
		++nej;
		nextg = locg[j];
		indg[nextg] = i__;
		locg[j] = nextg + 1;
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Place any linear elements after the nonlinear elements. */
/* ----------------------------------------------------------------- */
/* locJ (j) points to the next nonzero in the    linear part. */
/* iy(j)    points to the next nonzero in the nonlinear part. */
    if (*nnjac > 0) {
	i__1 = *nnjac;
	for (j = 1; j <= i__1; ++j) {
	    ncolj = nglin[j];
	    if (ncolj > 0) {
		nextj = locj[j];
		locj[j] += ncolj;
		iy[j] = nextj;
	    }
	}
    }
    i__1 = *nea;
    for (k = 1; k <= i__1; ++k) {
	in = iafun[k];
	jn = javar[k];
	row = *n + in;
	col = jn;
	i__ = kx[row];
	j = kx[col];
/* Assign the elements of  indJ  and  Jcol. */
	if (i__ <= *nncon && j <= *nnjac) {
	    nextj = iy[j];
	    indj[nextj] = i__;
	    jcol[nextj] = a[k];
	    iy[j] = nextj + 1;
	    ++nej;
	    nextg = locg[j];
	    indg[nextg] = i__;
	    locg[j] = nextg + 1;
	} else {
	    nextj = locj[j];
	    indj[nextj] = i__;
	    jcol[nextj] = a[k];
	    locj[j] = nextj + 1;
	    ++nej;
	}
    }
/* Check for an empty Jacobian. */
/* In this case, ne = 1 was set in s3SizA. */
    if (nej == 0) {
	jcol[*ne] = 0.;
	indj[*ne] = 1;
    }
/* Reset the column pointers. */
    for (j = *n; j >= 2; --j) {
	locj[j] = locj[j - 1];
    }
    locj[1] = 1;
    for (j = *nnjac; j >= 2; --j) {
	locg[j] = locg[j - 1];
    }
    locg[1] = 1;
    return 0;
} /* s3blda_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3bldA */
/* Subroutine */ int s3ina_(integer *start, integer *iobj, integer *m, 
	integer *n, integer *nb, integer *nncon, integer *nf, integer *nkx, 
	integer *kxn, doublereal *infbnd, char *xnames, char *fnames, char *
	names, integer *nxname, integer *nfname, integer *nname, doublereal *
	xlow, doublereal *xupp, doublereal *flow, doublereal *fupp, 
	doublereal *bl, doublereal *bu, integer *xstate, integer *fstate, 
	integer *hs, doublereal *xn, doublereal *f, doublereal *x, doublereal 
	*fx, doublereal *fmul, doublereal *pi, ftnlen xnames_len, ftnlen 
	fnames_len, ftnlen names_len)
{
    /* Initialized data */

    static char colnm[1] = "x";
    static char rownm[1] = "r";

    /* System generated locals */
    integer i__1;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, in, jn;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), chcopy_(integer *, char *, integer *, char *, integer 
	    *, ftnlen, ftnlen);

/* ================================================================= */
/* s3inA  loads re-ordered copies of the user defined values of */
/* bl, bu, x, pi  and  hs into the SNOPT arrays. */

/* 25 Sep 1998: First version of s3inA. */
/* 28 Jun 2000: First version for snoptA. */
/* 06 Apr 2008: ObjRow > 0 implies iObj > 0. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --pi;
    --xn;
    --xstate;
    --xupp;
    --xlow;
    --x;
    --hs;
    --bu;
    --bl;
    --fmul;
    --fx;
    --f;
    --fstate;
    --fupp;
    --flow;
    --kxn;
    xnames -= 8;
    fnames -= 8;
    names -= 8;

    /* Function Body */
/* ----------------------------------------------------------------- */
/* Load the bounds on the variables first. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jn = kxn[j];
	bl[j] = xlow[jn];
	bu[j] = xupp[jn];
	x[j] = xn[jn];
	hs[j] = xstate[jn];
    }
/* Bounds on the linear and nonlinear rows. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	in = kxn[j];
	if (i__ <= *nncon) {
	    pi[i__] = fmul[in];
	}
/* The linear objective row is always free. */
	if (i__ == *iobj) {
	    bl[j] = -(*infbnd);
	    bu[j] = *infbnd;
	    if (*start == 2) {
		hs[j] = 0;
		x[j] = 0.;
	    }
	} else {
	    bl[j] = flow[in];
	    bu[j] = fupp[in];
	    if (*start == 2) {
		hs[j] = fstate[in];
		x[j] = f[in];
	    } else {
		hs[j] = 0;
	    }
	}
    }
    dload_(m, &c_b367, &pi[1], &c__1);
    if (*nncon > 0) {
	dload_(nncon, &c_b367, &fx[1], &c__1);
    }
    if (*nxname > 1 || *nfname > 1) {
/* User has supplied some row or column names. */
	if (*nxname > 1) {
	    chcopy_(nxname, xnames + 8, &c__1, names + 8, &c__1, (ftnlen)8, (
		    ftnlen)8);
	} else {
/* Provide generic variable names. */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 8;
		ici__1.iciunit = names + (j << 3);
		ici__1.icifmt = "(a1,i7)";
		s_wsfi(&ici__1);
		do_fio(&c__1, colnm, (ftnlen)1);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		e_wsfi();
	    }
	}
	if (*nfname > 1) {
	    chcopy_(nfname, fnames + 8, &c__1, names + (*n + 1 << 3), &c__1, (
		    ftnlen)8, (ftnlen)8);
	} else {
/* Provide generic constraint names. */
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 8;
		ici__1.iciunit = names + (i__ << 3);
		ici__1.icifmt = "(a1,i7)";
		s_wsfi(&ici__1);
		do_fio(&c__1, rownm, (ftnlen)1);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfi();
	    }
	}
    }
    return 0;
} /* s3ina_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3inA */
/* Subroutine */ int s3outa_(integer *n, integer *nb, integer *nf, integer *
	nncon, integer *nkx, integer *kxn, integer *objrow, integer *iobj, 
	doublereal *fobj, integer *xstate, integer *fstate, integer *hs, 
	doublereal *xn, doublereal *f, doublereal *x, doublereal *fx, 
	doublereal *xmul, doublereal *fmul, doublereal *rc)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, in, jn;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *);

/* ================================================================= */
/* s3outA prepares the data for final exit. Various arrays are */
/* un-permuted and the working versions of x, rc  and  hs are */
/* copied into their originals. */

/* 25 Sep 1998: First version of s3outA. */
/* 18 Jul 2000: First version for snoptA. */
/* 17 Jul 2005: fObj assigned to F + x(n+iObj) instead of ObjTru. */
/* 06 Apr 2008: ObjRow > 0 implies iObj > 0. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --xmul;
    --xn;
    --xstate;
    --rc;
    --x;
    --hs;
    --fmul;
    --fx;
    --f;
    --fstate;
    --kxn;

    /* Function Body */
    iload_(n, &c__0, &xstate[1], &c__1);
    dload_(n, &c_b367, &xmul[1], &c__1);
    iload_(nf, &c__0, &fstate[1], &c__1);
    dload_(nf, &c_b367, &fmul[1], &c__1);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jn = kxn[j];
	xn[jn] = x[j];
	xstate[jn] = hs[j];
	xmul[jn] = rc[j];
    }
    i__1 = *nkx;
    for (j = *n + 1; j <= i__1; ++j) {
	i__ = j - *n;
	in = kxn[j];
	if (in == *objrow) {
	    if (*iobj > 0) {
		f[*objrow] = *fobj + x[*n + *iobj];
	    } else {
		f[*objrow] = *fobj;
	    }
	} else {
	    if (i__ <= *nncon) {
		f[in] = fx[i__];
	    } else {
		f[in] = x[j];
	    }
	    fstate[in] = hs[j];
	    fmul[in] = rc[j];
	}
    }
    return 0;
} /* s3outa_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3outA */
/* Subroutine */ int s3prma_(integer *n, integer *nf, integer *nkx, integer *
	igfunn, integer *jgvarn, integer *igfun, integer *jgvar, integer *
	leng, integer *neg, integer *kx, integer *kxn)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, jn, col, row, rown;

/* ================================================================= */
/* s3prmA makes a reordered copy of the G coordinate array, */
/* and inverts the row and column permutations. */

/* 29 Dec 2000: First version of s3prmA. */
/* 06 Jun 2001: Current version. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* Assign the reordered constant Jacobian elements. */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --kxn;
    --kx;
    --jgvar;
    --igfun;
    --jgvarn;
    --igfunn;

    /* Function Body */
    i__1 = *neg;
    for (k = 1; k <= i__1; ++k) {
	row = *n + igfunn[k];
	col = jgvarn[k];
	igfun[k] = kx[row];
	jgvar[k] = kx[col];
    }
/* ----------------------------------------------------------------- */
/* Invert the row and column lists. */
/* Note: if iObj = 0, then kx(ObjRow) = 0, kxN(nF) = ObjRow. */
/* ----------------------------------------------------------------- */
    i__1 = *n;
    for (jn = 1; jn <= i__1; ++jn) {
	j = kx[jn];
	kxn[j] = jn;
    }
    i__1 = *nkx;
    for (rown = *n + 1; rown <= i__1; ++rown) {
	i__ = kx[rown];
	kxn[*n + i__] = rown - *n;
    }
    return 0;
} /* s3prma_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3prmA */
/* Subroutine */ int s3mapa_(integer *m, integer *n, integer *ne, integer *nf,
	 integer *neg, integer *negcon, integer *nkx, integer *nnjac, integer 
	*nname, integer *nextcw, integer *nextiw, integer *nextrw, integer *
	iw, integer *leniw)
{
    static integer nb, lf, lg, lx, lbl, lrc, lbu, lpi, lhs, lxn, lkxn, lindg, 
	    lindj, nlocg, ljcol, llocj, nlocj, lfmul, lennam, lnames, lnglin, 
	    ligfun, ljgvar;

/* ================================================================= */
/* s3mapA   allocates additional array storage for snoptA. */

/* 29 Dec 2000: First version of s3mapA. */
/* 27 Dec 2002: Current version. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    nb = *n + *m;
/* Nonlinear constraints. */
    nlocg = *nnjac + 1;
    nlocj = *n + 1;
    lennam = 0;
    if (*nname == nb) {
	lennam = *nname;
    }
/* Define the addresses. */
    if (lennam == 0) {
	lnames = *nextcw - 1;
/* cw(lNames) unused */
    } else {
	lnames = *nextcw;
	*nextcw = lnames + lennam;
    }
    lhs = *nextiw;
    lkxn = lhs + nb;
    llocj = lkxn + *nkx;
    lindj = llocj + nlocj;
    lindg = lindj + *ne;
    lnglin = lindg + *negcon;
    ligfun = lnglin + nlocg;
    ljgvar = ligfun + *neg;
    *nextiw = ljgvar + *neg;
    ljcol = *nextrw;
    lbl = ljcol + *ne;
    lbu = lbl + nb;
    lpi = lbu + nb;
    lrc = lpi + *m;
    lx = lrc + nb;
    lxn = lx + nb;
    lf = lxn + *n;
    lfmul = lf + *nf;
    lg = lfmul + *nf;
    *nextrw = lg + *neg;
/* Store the addresses in iw. */
    iw[252] = lkxn;
/* jN = kxN(j ) => col j of Jcol is variable jN */
    iw[256] = ljcol;
/* Jcol(ne)    = Constraint Jacobian by columns */
    iw[257] = llocj;
/* locJ(n+1)   = column pointers for indJ */
    iw[258] = lindj;
/* indJ(ne) holds the row indices for Jij */
    iw[261] = lindg;
/* indG(negCon) holds the row indices for gij */
    iw[262] = lnglin;
/* nGlin(j) = # linear elems in col j of gCon */
    iw[266] = ligfun;
/* iGfun(neG) row list of reordered G nonzeros */
    iw[267] = ljgvar;
/* iGvar(neG) col list of reordered G nonzeros */
    iw[271] = lbl;
/* bl(nb)      = lower bounds */
    iw[272] = lbu;
/* bu(nb)      = upper bounds */
    iw[279] = lpi;
/* pi(m)       = the pi-vector */
    iw[280] = lrc;
/* rc(nb)      = the reduced costs */
    iw[282] = lhs;
/* the column state vector */
    iw[299] = lx;
/* x(nb)       = the solution (x,s) */
    iw[327] = lxn;
/* xN(n)       = copy of user-defined x */
    iw[328] = lf;
/* F(nF)       = user-defined F */
    iw[329] = lfmul;
/* Fmul(nF)    = user-defined multipliers */
    iw[330] = lg;
/* G (lenG)    = problem derivatives */
    iw[359] = lnames;
/* Names(nName), row and column names */
    return 0;
} /* s3mapa_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3mapA */
/* Subroutine */ int s3argb_(integer *iexit, char *start, integer *m, integer 
	*n, integer *ne, integer *nname, integer *ns, integer *nncon, integer 
	*nnobj, integer *nnjac, integer *iobj, integer *indj, integer *locj, 
	doublereal *bl, doublereal *bu, char *names, integer *hs, doublereal *
	pi, integer *lvlsrt, integer *errors, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen start_len, ftnlen names_len)
{
    /* Initialized data */

    static char id[5*3] = "varbl" "lncon" "nlcon";

    /* Format strings */
    static char fmt_1000[] = "(\002 XXX Start parameter not recognized:  "
	    "\002,a)";
    static char fmt_1100[] = "(\002 XXX  Argument out of range:  \002,a6,"
	    "\002 = \002,i6)";
    static char fmt_1200[] = "(\002 XXX  Invalid locA(1), locA(n+1) =\002,2i"
	    "8)";
    static char fmt_1300[] = "(\002 XXX  Invalid argument indA: linear row"
	    " \002,i6,\002 in column \002,i6,\002 appears before nonlinear ro"
	    "w \002,i6)";
    static char fmt_1400[] = "(\002 XXX  Invalid argument hs: \002,i6,\002 e"
	    "lements modified to be in range.\002)";
    static char fmt_1510[] = "(\002 XXX  The equal bounds on  \002,a8,\002  "
	    "are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16.7)";
    static char fmt_1515[] = "(\002 XXX  The bounds on  \002,a8,\002  are in"
	    "consistent.   bl =\002,g16.7,\002   bu =\002,g16.7)";
    static char fmt_1500[] = "(\002 XXX  The equal bounds on  \002,a5,i6,"
	    "\002  are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16."
	    "7)";
    static char fmt_1505[] = "(\002 XXX  The bounds on  \002,a5,i6,\002  are"
	    " inconsistent.   bl =\002,g16.7,\002   bu =\002,g16.7)";
    static char fmt_1600[] = "(\002 XXX  Invalid arguments bl, bu: \002,i6"
	    ",\002 inconsistent bounds or infinite equal bounds.\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, l;
    static doublereal b1, b2;
    static integer n1, nb;
    static logical ok;
    static integer ir, js;
    static char ch1[1], str[132];
    static integer mods, errs;
    static logical named;
    static integer nchar;
    extern /* Subroutine */ int s1trim_(char *, integer *, ftnlen);
    static doublereal infbnd;
    static integer argerr, linrow;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___249 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___250 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___251 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___254 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___255 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___256 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___257 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___258 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___259 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___260 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___261 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___262 = { 0, str, 0, fmt_1200, 132, 1 };
    static icilist io___265 = { 0, str, 0, fmt_1300, 132, 1 };
    static icilist io___268 = { 0, str, 0, fmt_1400, 132, 1 };
    static icilist io___274 = { 0, str, 0, fmt_1510, 132, 1 };
    static icilist io___275 = { 0, str, 0, fmt_1515, 132, 1 };
    static icilist io___276 = { 0, str, 0, fmt_1500, 132, 1 };
    static icilist io___277 = { 0, str, 0, fmt_1505, 132, 1 };
    static icilist io___278 = { 0, str, 0, fmt_1600, 132, 1 };


/* ================================================================= */
/* s3argB   checks the arguments for snOpt. */

/* On exit, Errors says how many errors were encountered. */

/* 21 Dec 2002: First version of s3argB. */
/* 03 Aug 2003: snPRNT adopted. */
/* 14 Sep 2003: hs = -1 accepted. */
/* 13 Mar 2004: Start = 'Hot FHS' options decoded. */
/* 26 Mar 2004: Prevent n*m overflow (noticed by Anders Goran). */
/* 04 Sep 2004: Length of str increased to 132. */
/* 04 Dec 2004: Current version of s3argB. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --pi;
    --hs;
    --bu;
    --bl;
    --locj;
    --indj;
    names -= 8;
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
    infbnd = rw[70];
/* definition of an infinite bound */
    argerr = iw[106];
/* maximum # errors in MPS data */
    *iexit = 0;
    *errors = 0;
/* The options haven't been checked yet. */
    if (infbnd < 0.) {
	infbnd = 1e20;
    }
    if (argerr < 0) {
	argerr = 20;
    }
/* Print 20 lines max */
    nb = *n + *m;
    named = *nname == nb;
/* ================================================================= */
/* Decode 'Start'. */
/* ================================================================= */
    iw[230] = 0;
    iw[231] = 0;
    iw[232] = 0;
/* Determine the type of start. */
/* Optional parameters take precedence. */
    if (*lvlsrt != -11111) {
/* lvlSrt was set as an option */
	if (*lvlsrt == 0 || *lvlsrt == 1 || *lvlsrt == 2) {
/* Relax */
	} else if (*lvlsrt == 3) {
	    iw[230] = 1;
	    iw[231] = 1;
	    iw[232] = 1;
	} else {
/* lvlSrt is an unrecognized option */
	    *lvlsrt = -11111;
	}
    }
    if (*lvlsrt == -11111) {
/* lvlSrt is unset */
	*(unsigned char *)ch1 = *(unsigned char *)start;
	if (*(unsigned char *)ch1 == 'C' || *(unsigned char *)ch1 == 'c') {
	    *lvlsrt = 0;
	} else if (*(unsigned char *)ch1 == 'B' || *(unsigned char *)ch1 == 
		'b') {
	    *lvlsrt = 1;
	} else if (*(unsigned char *)ch1 == 'W' || *(unsigned char *)ch1 == 
		'w') {
	    *lvlsrt = 2;
	} else if (*(unsigned char *)ch1 == 'H' || *(unsigned char *)ch1 == 
		'h') {
	    *lvlsrt = 3;
/* nchar  = len_trim(Start)         ! An F90 intrinsic */
	    s1trim_(start, &nchar, start_len);
/* Decode    Start = 'HOT .. */
/* The F77 equivalent */
	    if (nchar <= 4) {
/* 'Hot' or 'Hot ' = 'Hot FH */
		iw[230] = 1;
		iw[231] = 1;
		iw[232] = 1;
	    } else {
		i__1 = nchar;
		for (j = 5; j <= i__1; ++j) {
/* Decode 1 or more of FHS */
		    *(unsigned char *)ch1 = *(unsigned char *)&start[j - 1];
		    if (*(unsigned char *)ch1 == 'F' || *(unsigned char *)ch1 
			    == 'f') {
			iw[230] = 1;
		    }
		    if (*(unsigned char *)ch1 == 'H' || *(unsigned char *)ch1 
			    == 'h') {
			iw[231] = 1;
		    }
		    if (*(unsigned char *)ch1 == 'S' || *(unsigned char *)ch1 
			    == 's') {
			iw[232] = 1;
		    }
		}
	    }
	} else {
	    *lvlsrt = 0;
	    s_wsfi(&io___249);
	    do_fio(&c__1, start, start_len);
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/* ================================================================= */
/* Check the other arguments. */
/* ================================================================= */
    if (*m < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___250);
	    do_fio(&c__1, "m     ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*n < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___251);
	    do_fio(&c__1, "n     ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/* 26 Mar 2004: Test if ne > n*m without causing overflow. */
    n1 = max(*n,1);
    k = *ne % n1;
    if (*ne < 1 || k == 0 && *ne / n1 > *m || k > 0 && *ne / n1 >= *m) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___254);
	    do_fio(&c__1, "ne    ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nnobj < 0 || *nnobj > *n) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___255);
	    do_fio(&c__1, "nnObj ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nnobj), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nncon < 0 || *nncon > *m) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___256);
	    do_fio(&c__1, "nnCon ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nncon), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nnjac < 0 || *nnjac > *n) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___257);
	    do_fio(&c__1, "nnJac ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nnjac), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nncon == 0 && *nnjac > 0) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___258);
	    do_fio(&c__1, "nnJac ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nnjac), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nncon > 0 && *nnjac == 0) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___259);
	    do_fio(&c__1, "nnJac ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nnjac), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nname != 1 && *nname != nb) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___260);
	    do_fio(&c__1, "nName ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nname), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*iobj > 0 && *iobj < *nncon || *iobj < 0 || *iobj > *m) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___261);
	    do_fio(&c__1, "iObj  ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*iobj), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (locj[1] != 1 || locj[*n + 1] != *ne + 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___262);
	    do_fio(&c__1, (char *)&locj[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&locj[*n + 1], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	linrow = 0;
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
/* Linear row */
		linrow = ir;
	    } else if (linrow > 0) {
/* Nonlinear row occurs after a linear row. */
		++(*errors);
		if (*errors <= argerr) {
		    s_wsfi(&io___265);
		    do_fio(&c__1, (char *)&linrow, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		}
	    }
	}
    }
    mods = 0;
    if (*lvlsrt == 0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    js = hs[j];
	    if (js < -1 || js > 5) {
		++mods;
		hs[j] = 0;
	    }
	}
    } else if (*lvlsrt == 2 || *lvlsrt == 3) {
	i__1 = nb;
	for (j = 1; j <= i__1; ++j) {
	    js = hs[j];
	    if (js < 0 || js > 3) {
		++mods;
		hs[j] = 0;
	    }
	}
    }
    if (mods > 0) {
	s_wsfi(&io___268);
	do_fio(&c__1, (char *)&mods, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------------------------------------------------- */
/* Check the bounds on all variables and constraints. */
/* ----------------------------------------------------------------- */
    errs = 0;
    i__1 = nb;
    for (j = 1; j <= i__1; ++j) {
	b1 = bl[j];
	b2 = bu[j];
	ok = b1 <= b2 && b1 < infbnd && b2 > -infbnd;
	if (! ok) {
	    ++errs;
	    ++(*errors);
	    if (j > *n + *nncon) {
/* Linear    constraint */
		k = j - (*n + *nncon);
		l = 2;
	    } else if (j > *n) {
/* Nonlinear constraint */
		k = j - *n;
		l = 3;
	    } else {
/* Variable */
		k = j;
		l = 1;
	    }
	    if (*errors <= argerr) {
		if (named) {
		    if (b1 == b2) {
			s_wsfi(&io___274);
			do_fio(&c__1, names + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___275);
			do_fio(&c__1, names + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (b1 == b2) {
			s_wsfi(&io___276);
			do_fio(&c__1, id + (l - 1) * 5, (ftnlen)5);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___277);
			do_fio(&c__1, id + (l - 1) * 5, (ftnlen)5);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	}
    }
    if (errs > 0) {
	s_wsfi(&io___278);
	do_fio(&c__1, (char *)&errs, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    if (*errors >= argerr) {
	snprnt_(&c__3, " and so on ...", &iw[1], leniw, (ftnlen)14);
    }
    if (*errors > 0) {
	*iexit = 91;
/* Invalid arguments for snOptB. */
    }
    return 0;
} /* s3argb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3argB */
/* Subroutine */ int s3prtb_(integer *m, integer *n, integer *nncon, integer *
	nnjac, integer *nnobj, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Initialized data */

    static char hestyp[24*3] = " Limited-Memory Hessian." " Full-Memory Hess"
	    "ian...." " Exact Hessian..........";
    static char lsrch[24*2] = " Nonderiv.  linesearch.." " Derivative linese"
	    "arch..";
    static char pivtyp[24*4] = " LU partial  pivoting..." " LU rook     pivo"
	    "ting..." " LU complete pivoting..." " LU diagonal pivoting...";
    static char prbtyp[24*3] = " Maximize..............." " Feasible point o"
	    "nly...." " Minimize...............";
    static char qptype[24*3] = " QPsolver Cholesky......" " QPsolver CG....."
	    "......." " QPsolver QN............";
    static char srttyp[24*4] = " Cold start............." " Basis file......"
	    "......." " Warm start............." " Hot start..............";
    static char noyes[3*2] = " No" "Yes";

    /* Format strings */
    static char fmt_2110[] = "(\002 Solution file..........\002,i10,6x,\002 "
	    "Old basis file ........\002,i10,6x,\002 Standard input........"
	    ".\002,i10)";
    static char fmt_2120[] = "(\002 Insert file............\002,i10,6x,\002 "
	    "New basis file ........\002,i10,6x,\002 (Printer)............."
	    ".\002,i10)";
    static char fmt_2130[] = "(\002 Punch file.............\002,i10,6x,\002 "
	    "Backup basis file......\002,i10,6x,\002 (Specs file).........."
	    ".\002,i10)";
    static char fmt_2140[] = "(\002 Load file..............\002,i10,6x,\002 "
	    "Dump file..............\002,i10,6x,\002 Standard output......."
	    ".\002,i10)";
    static char fmt_2210[] = "(\002 Print frequency........\002,i10,6x,\002 "
	    "Check frequency........\002,i10,6x,\002 Save new basis map...."
	    ".\002,i10)";
    static char fmt_2220[] = "(\002 Summary frequency......\002,i10,6x,\002 "
	    "Factorization frequency\002,i10,6x,\002 Expand frequency......"
	    ".\002,i10)";
    static char fmt_2310[] = "(a24)";
    static char fmt_2320[] = "(\002 Scale tolerance........\002,0p,f10.3,6x"
	    ",\002 Minor feasibility tol..\002,1p,e10.2,6x,\002 Iteration lim"
	    "it........\002,i10)";
    static char fmt_2330[] = "(\002 Scale option...........\002,i10,6x,\002 "
	    "Minor optimality  tol..\002,1p,e10.2,6x,\002 Minor print level.."
	    "....\002,i10)";
    static char fmt_2340[] = "(\002 Crash tolerance........\002,0p,f10.3,6x"
	    ",\002 Pivot tolerance........\002,1p,e10.2,6x,\002 Partial price"
	    "..........\002,i10)";
    static char fmt_2350[] = "(\002 Crash option...........\002,i10,6x,\002 "
	    "Elastic weight.........\002,1p,e10.2,6x,\002 Prtl price section "
	    "( A)\002,i10)";
    static char fmt_2360[] = "(40x,\002 New superbasics........\002,i10,6x"
	    ",\002 Prtl price section (-I)\002,i10)";
    static char fmt_2380[] = "(\002 Subspace tolerance.....\002,0p,f10.5,6x"
	    ",\002 CG tolerance...........\002,1p,e10.2,6x,\002 CG Iterations"
	    "..........\002,i10)";
    static char fmt_2390[] = "(80x,\002 CG preconditioning.....\002,i10)";
    static char fmt_2410[] = "(a24,16x,a24,16x,\002 Proximal Point method."
	    ".\002,i10)";
    static char fmt_2420[] = "(\002 Nonlinear objectiv vars\002,i10,6x,\002 "
	    "Major optimality tol...\002,1p,e10.2,6x,\002 Function precision."
	    "....\002,1p,e10.2)";
    static char fmt_2430[] = "(\002 Unbounded step size....\002,1p,e10.2,6x"
	    ",\002 Superbasics limit......\002,i10,6x,\002 Difference interva"
	    "l....\002,1p,e10.2)";
    static char fmt_2440[] = "(\002 Unbounded objective....\002,1p,e10.2,6x"
	    ",\002 Reduced Hessian dim....\002,i10,6x,\002 Central difference"
	    " int.\002,1p,e10.2)";
    static char fmt_2450[] = "(\002 Major step limit.......\002,1p,e10.2,6x,"
	    "a11,\002 linesearch..\002,16x,\002 Derivative level.......\002,i"
	    "10)";
    static char fmt_2460[] = "(\002 Major iterations limit.\002,i10,6x,\002 "
	    "Linesearch tolerance...\002,0p,f10.5,6x,\002 Verify level......."
	    "....\002,i10)";
    static char fmt_2470[] = "(\002 Minor iterations limit.\002,i10,6x,\002 "
	    "Penalty parameter......\002,1p,e10.2,6x,\002 Major Print Level.."
	    "....\002,i10)";
    static char fmt_2510[] = "(a24,16x,\002 Hessian updates........\002,i10,"
	    "6x,\002 Hessian frequency......\002,i10)";
    static char fmt_2520[] = "(80x,\002 Hessian flush..........\002,i10)";
    static char fmt_2610[] = "(\002 Nonlinear constraints..\002,i10,6x,\002 "
	    "Major feasibility tol..\002,1p,e10.2,6x,\002 Violation limit...."
	    "....\002,e10.2)";
    static char fmt_2620[] = "(\002 Nonlinear Jacobian vars\002,i10)";
    static char fmt_2710[] = "(\002 LU factor tolerance....\002,0p,f10.2,6x"
	    ",\002 LU singularity tol.....\002,1p,e10.2,6x,\002 Timing level."
	    "..........\002,i10)";
    static char fmt_2720[] = "(\002 LU update tolerance....\002,0p,f10.2,6x"
	    ",\002 LU swap tolerance......\002,1p,e10.2,6x,\002 Debug level.."
	    "..........\002,i10)";
    static char fmt_2730[] = "(a24,16x,\002 eps (machine precision)\002,1p,e"
	    "10.2,6x,\002 System information.....\002,7x,a3)";
    static char fmt_2740[] = "(80x,\002 Sticky parameters......\002,7x,a3)";
    static char fmt_3000[] = "(\002 Scale option\002,i3,\002,    Partial pri"
	    "ce\002,i4)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer nnl;
    static doublereal eps;
    static integer npr1, npr2;
    static char str1[132], str2[132], str3[132], str4[132], str5[132], str6[
	    132];
    static integer kfac, kchk, klog, ksav, maxr, maxs;
    static doublereal tolx, xpen0, utol1;
    static integer iback, ioldb;
    static doublereal bigdx, bigfx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal tolcg;
    static integer istdi;
    static doublereal xdlim;
    static integer idump, maxmn;
    static doublereal epsrf;
    static integer istdo;
    static doublereal vilim;
    static integer isoln, ksumm;
    static doublereal tolqp;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *);
    static doublereal fdint1, fdint2, wtinf0;
    static integer iloadb, kdegen;
    static doublereal tolfac;
    static integer icrash, lprdbg;
    static doublereal wolfeg, tcrash;
    static integer mmajor, ispecs, lvlder, minmax, minprc, cgitmx, itnlim, 
	    kreset, lvlhes, lvlsch, lvlscl, mflush, mminor, mnewsb, lvlpre, 
	    iprint, mqnmod, lvltim, iinsrt, nparpr, lvlppm, lvlver, lprprm, 
	    lvlpiv, mjrprt;
    static doublereal scltol, tolcon, tolnlp;
    static integer mnrprt;
    static doublereal tolpiv, tolupd;
    static integer lvlsrt, qpslvr;
    static doublereal tolswp;
    static integer stkyop, lvlsys;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___364 = { 0, str1, 0, fmt_2110, 132, 1 };
    static icilist io___366 = { 0, str2, 0, fmt_2120, 132, 1 };
    static icilist io___368 = { 0, str3, 0, fmt_2130, 132, 1 };
    static icilist io___370 = { 0, str4, 0, fmt_2140, 132, 1 };
    static icilist io___371 = { 0, str1, 0, fmt_2210, 132, 1 };
    static icilist io___372 = { 0, str2, 0, fmt_2220, 132, 1 };
    static icilist io___373 = { 0, str1, 0, fmt_2310, 132, 1 };
    static icilist io___374 = { 0, str2, 0, fmt_2320, 132, 1 };
    static icilist io___375 = { 0, str3, 0, fmt_2330, 132, 1 };
    static icilist io___376 = { 0, str4, 0, fmt_2340, 132, 1 };
    static icilist io___378 = { 0, str5, 0, fmt_2350, 132, 1 };
    static icilist io___380 = { 0, str6, 0, fmt_2360, 132, 1 };
    static icilist io___381 = { 0, str1, 0, fmt_2380, 132, 1 };
    static icilist io___382 = { 0, str2, 0, fmt_2390, 132, 1 };
    static icilist io___383 = { 0, str1, 0, fmt_2410, 132, 1 };
    static icilist io___384 = { 0, str2, 0, fmt_2420, 132, 1 };
    static icilist io___385 = { 0, str3, 0, fmt_2430, 132, 1 };
    static icilist io___386 = { 0, str4, 0, fmt_2440, 132, 1 };
    static icilist io___387 = { 0, str1, 0, fmt_2450, 132, 1 };
    static icilist io___388 = { 0, str2, 0, fmt_2460, 132, 1 };
    static icilist io___389 = { 0, str3, 0, fmt_2470, 132, 1 };
    static icilist io___390 = { 0, str1, 0, fmt_2510, 132, 1 };
    static icilist io___391 = { 0, str2, 0, fmt_2520, 132, 1 };
    static icilist io___392 = { 0, str1, 0, fmt_2610, 132, 1 };
    static icilist io___393 = { 0, str2, 0, fmt_2620, 132, 1 };
    static icilist io___394 = { 0, str1, 0, fmt_2710, 132, 1 };
    static icilist io___395 = { 0, str2, 0, fmt_2720, 132, 1 };
    static icilist io___396 = { 0, str3, 0, fmt_2730, 132, 1 };
    static icilist io___397 = { 0, str4, 0, fmt_2740, 132, 1 };
    static icilist io___398 = { 0, str1, 0, fmt_3000, 132, 1 };


/* ================================================================= */
/* s3prtB prints the settings of the optional parameters for the */
/* standard SNOPT wrapper. */

/* See  snworkspace.info  for full documentation of cw, iw and rw. */

/* 15 Nov 1991: First version of s3prtB (s8prnt). */
/* 10 Dec 2002: More LU pivoting options. */
/* 01 Jul 2003: QN QP options added. */
/* 30 Jul 2003: QN CG options added. */
/* 03 Aug 2003: snPRNT adopted. */
/* 01 Sep 2007: Sticky parameters added. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
/* Set some local machine-dependent constants. */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
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
/* excludes small elements of y. */
    tcrash = rw[62];
/* crash tolerance. */
    tolswp = rw[65];
/* LU swap tolerance. */
    tolfac = rw[66];
/* LU factor tolerance. */
    tolupd = rw[67];
/* LU update tolerance. */
    bigfx = rw[71];
/* unbounded objective. */
    bigdx = rw[72];
/* unbounded step. */
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
    wtinf0 = rw[88];
/* infeasibility weight */
    xpen0 = rw[89];
/* initial penalty parameter. */
    scltol = rw[92];
/* scale tolerance. */
    utol1 = rw[154];
/* abs tol for small diag of U. */
    istdi = iw[9];
/* Standard Input */
    istdo = iw[10];
/* Standard Output */
    ispecs = iw[11];
/* Specs (options) file */
    iprint = iw[12];
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    qpslvr = iw[55];
/* 0(1) => QP(QN) QP solver */
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
    lvlsrt = iw[69];
/* = 0:1:2:3 => cold:warm:basis:hot start */
    lvlder = iw[70];
/* = 0, 1 or 2, the derivative level */
    lvlsys = iw[71];
/* > 0   => print system info */
    lvlhes = iw[72];
/* 0,1,2 => LM, FM, Exact Hessian */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lvlscl = iw[75];
/* scale option */
    lvlsch = iw[76];
/* >0    => use derivatives in the line search */
    lvlpre = iw[77];
/* >0    => QN preconditioned CG */
    lvlver = iw[78];
/* Verify level */
    lvlppm = iw[79];
/* Proximal Point method for x0 */
    lvlpiv = iw[80];
/* 0(1) LU threshold partial(complete) pivoting */
    lprprm = iw[81];
/* > 0    => parms are printed */
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
    isoln = iw[131];
/* solution file */
    lvltim = iw[182];
/* ----------------------------------------------------------------- */
/* Timing level */
    if (iprint <= 0 || mjrprt == 0 || lprprm == 0) {
	return 0;
    }
    nnl = max(*nnjac,*nnobj);
    minprc = 10;
    npr1 = *n / nparpr;
    npr2 = *m / nparpr;
    if (max(npr1,npr2) < minprc) {
	maxmn = max(*m,*n);
	nparpr = maxmn / min(maxmn,minprc);
	npr1 = *n / nparpr;
	npr2 = *m / nparpr;
    }
/* ================================================================= */
/* Print parameters except if PRINT LEVEL = 0 */
/* or SUPPRESS PARAMETERS was specified. */
/* ================================================================= */
    s1page_(&c__1, &iw[1], leniw);
    snprnt_(&c__1, " Parameters", &iw[1], leniw, (ftnlen)11);
    snprnt_(&c__1, " ==========", &iw[1], leniw, (ftnlen)11);
/* ----------------------- */
/* Files. */
/* ----------------------- */
    snprnt_(&c__11, " Files", &iw[1], leniw, (ftnlen)6);
    snprnt_(&c__1, " -----", &iw[1], leniw, (ftnlen)6);
    s_wsfi(&io___364);
    do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ioldb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&istdi, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___366);
    do_fio(&c__1, (char *)&iinsrt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inewb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iprint, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___368);
    do_fio(&c__1, (char *)&ipnch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iback, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ispecs, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___370);
    do_fio(&c__1, (char *)&iloadb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&idump, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&istdo, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* Frequencies. */
/* ----------------------- */
    snprnt_(&c__11, " Frequencies", &iw[1], leniw, (ftnlen)12);
    snprnt_(&c__1, " -----------", &iw[1], leniw, (ftnlen)12);
    s_wsfi(&io___371);
    do_fio(&c__1, (char *)&klog, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kchk, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ksav, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___372);
    do_fio(&c__1, (char *)&ksumm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kfac, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kdegen, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* QP subproblems. */
/* ----------------------- */
    snprnt_(&c__11, " QP subproblems", &iw[1], leniw, (ftnlen)15);
    snprnt_(&c__1, " --------------", &iw[1], leniw, (ftnlen)15);
    s_wsfi(&io___373);
    do_fio(&c__1, qptype + qpslvr * 24, (ftnlen)24);
    e_wsfi();
    s_wsfi(&io___374);
    do_fio(&c__1, (char *)&scltol, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&itnlim, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___375);
    do_fio(&c__1, (char *)&lvlscl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&tolqp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&mnrprt, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___376);
    do_fio(&c__1, (char *)&tcrash, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolpiv, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nparpr, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___378);
    do_fio(&c__1, (char *)&icrash, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&wtinf0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&npr1, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___380);
    do_fio(&c__1, (char *)&mnewsb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npr2, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str5, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str6, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* CG QP solver. */
/* ----------------------- */
    if (qpslvr == 1 || maxr < maxs) {
	snprnt_(&c__11, " Conjugate-gradient QP solver", &iw[1], leniw, (
		ftnlen)29);
	snprnt_(&c__1, " ----------------------------", &iw[1], leniw, (
		ftnlen)29);
	s_wsfi(&io___381);
	do_fio(&c__1, (char *)&etarg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&tolcg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cgitmx, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___382);
	do_fio(&c__1, (char *)&lvlpre, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------- */
/* SQP method. */
/* ----------------------- */
    snprnt_(&c__11, " The SQP Method", &iw[1], leniw, (ftnlen)15);
    snprnt_(&c__1, " --------------", &iw[1], leniw, (ftnlen)15);
    s_wsfi(&io___383);
    do_fio(&c__1, prbtyp + (minmax + 1) * 24, (ftnlen)24);
    do_fio(&c__1, srttyp + lvlsrt * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___384);
    do_fio(&c__1, (char *)&(*nnobj), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&tolnlp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&epsrf, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___385);
    do_fio(&c__1, (char *)&bigdx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&maxs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&fdint1, (ftnlen)sizeof(doublereal));
    e_wsfi();
    s_wsfi(&io___386);
    do_fio(&c__1, (char *)&bigfx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&maxr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&fdint2, (ftnlen)sizeof(doublereal));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    s_wsfi(&io___387);
    do_fio(&c__1, (char *)&xdlim, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, lsrch + lvlsch * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&lvlder, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___388);
    do_fio(&c__1, (char *)&mmajor, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&wolfeg, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lvlver, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___389);
    do_fio(&c__1, (char *)&mminor, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&xpen0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&mjrprt, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* Hessian approximation. */
/* ----------------------- */
    if (nnl > 0) {
	snprnt_(&c__11, " Hessian Approximation", &iw[1], leniw, (ftnlen)22);
	snprnt_(&c__1, " ---------------------", &iw[1], leniw, (ftnlen)22);
	s_wsfi(&io___390);
	do_fio(&c__1, hestyp + lvlhes * 24, (ftnlen)24);
	do_fio(&c__1, (char *)&mqnmod, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kreset, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___391);
	do_fio(&c__1, (char *)&mflush, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------- */
/* Nonlinear constraints. */
/* ----------------------- */
    if (*nncon > 0) {
	snprnt_(&c__11, " Nonlinear constraints", &iw[1], leniw, (ftnlen)22);
	snprnt_(&c__1, " ---------------------", &iw[1], leniw, (ftnlen)22);
	s_wsfi(&io___392);
	do_fio(&c__1, (char *)&(*nncon), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&tolcon, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vilim, (ftnlen)sizeof(doublereal));
	e_wsfi();
	s_wsfi(&io___393);
	do_fio(&c__1, (char *)&(*nnjac), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------- */
/* Miscellaneous */
/* ----------------------- */
    snprnt_(&c__11, " Miscellaneous", &iw[1], leniw, (ftnlen)14);
    snprnt_(&c__1, " -------------", &iw[1], leniw, (ftnlen)14);
    s_wsfi(&io___394);
    do_fio(&c__1, (char *)&tolfac, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&utol1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lvltim, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___395);
    do_fio(&c__1, (char *)&tolupd, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolswp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lprdbg, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___396);
    do_fio(&c__1, pivtyp + lvlpiv * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, noyes + lvlsys * 3, (ftnlen)3);
    e_wsfi();
    s_wsfi(&io___397);
    do_fio(&c__1, noyes + stkyop * 3, (ftnlen)3);
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    s_wsfi(&io___398);
    do_fio(&c__1, (char *)&lvlscl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nparpr, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__22, str1, &iw[1], leniw, (ftnlen)132);
    return 0;
} /* s3prtb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3prtB */
/* Subroutine */ int s3argq_(integer *iexit, char *start, integer *m, integer 
	*n, integer *ne, integer *nname, integer *ns, integer *lencobj, 
	integer *iobj, integer *ncolh, integer *inda, integer *loca, 
	doublereal *bl, doublereal *bu, char *names, integer *hs, doublereal *
	pi, integer *lvlsrt, integer *errors, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen start_len, ftnlen names_len)
{
    /* Initialized data */

    static char id[5*2] = "varbl" "lncon";

    /* Format strings */
    static char fmt_1000[] = "(\002 XXX Start parameter not recognized:  "
	    "\002,a)";
    static char fmt_1100[] = "(\002 XXX  Argument out of range:  \002,a,\002"
	    " = \002,i6)";
    static char fmt_1200[] = "(\002 XXX  Invalid locA(1), locA(n+1) =\002,2i"
	    "8)";
    static char fmt_1300[] = "(\002 XXX  Invalid argument indA: row index"
	    " \002,i6,\002 in column \002,i6,\002 is out of range \002,i6)";
    static char fmt_1400[] = "(\002 XXX  Invalid argument hs: \002,i6,\002 e"
	    "lements modified to be in range.\002)";
    static char fmt_1510[] = "(\002 XXX  The equal bounds on  \002,a8,\002  "
	    "are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16.7)";
    static char fmt_1515[] = "(\002 XXX  The bounds on  \002,a8,\002  are in"
	    "consistent.   bl =\002,g16.7,\002   bu =\002,g16.7)";
    static char fmt_1500[] = "(\002 XXX  The equal bounds on  \002,a5,i6,"
	    "\002  are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16."
	    "7)";
    static char fmt_1505[] = "(\002 XXX  The bounds on  \002,a5,i6,\002  are"
	    " inconsistent.   bl =\002,g16.7,\002   bu =\002,g16.7)";
    static char fmt_1600[] = "(\002 XXX  Invalid arguments bl, bu: \002,i6"
	    ",\002 inconsistent bounds or infinite equal bounds.\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, l;
    static doublereal b1, b2;
    static integer n1, nb;
    static logical ok;
    static integer ir, js, nz;
    static char ch1[1], str[132];
    static integer mods, errs;
    static logical named;
    static integer nchar;
    extern /* Subroutine */ int s1trim_(char *, integer *, ftnlen);
    static doublereal infbnd;
    static integer argerr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___408 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___409 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___410 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___413 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___414 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___415 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___416 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___417 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___418 = { 0, str, 0, fmt_1200, 132, 1 };
    static icilist io___421 = { 0, str, 0, fmt_1300, 132, 1 };
    static icilist io___422 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___425 = { 0, str, 0, fmt_1400, 132, 1 };
    static icilist io___431 = { 0, str, 0, fmt_1510, 132, 1 };
    static icilist io___432 = { 0, str, 0, fmt_1515, 132, 1 };
    static icilist io___433 = { 0, str, 0, fmt_1500, 132, 1 };
    static icilist io___434 = { 0, str, 0, fmt_1505, 132, 1 };
    static icilist io___435 = { 0, str, 0, fmt_1600, 132, 1 };


/* ================================================================= */
/* s3argQ   checks the arguments for sqOpt. */

/* On exit, */
/* Errors says how many errors were encountered. */
/* lvlSrt is an integer version of Start: */
/*    Start   lvlSrt */
/*    'Cold'       0 */
/*    'Basis'      1 */
/*    'Warm'       2 */
/*    'Hot'        3 */
/* iw(gotFac,gotHes,gotScl) are set to 0 or 1. */
/*    'Hot F'     Keep Factors of basis (LU) */
/*    'Hot H'     Keep Hessian */
/*    'Hot S'     Keep Scales */
/*    'Hot FH'    Keep Factors and Hessian */
/*    'Hot FHS'   etc. */

/* 18 Oct 2003: First version of s3argQ. */
/* 08 Mar 2004: Start = 'Hot FHS' options decoded. */
/* 26 Mar 2004: Prevent n*m overflow (noticed by Anders Goran). */
/* 02 Sep 2004: Length of str increased to 132. */
/* 04 Dec 2004: Current version of s3argQ. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --pi;
    --hs;
    --bu;
    --bl;
    --loca;
    --inda;
    names -= 8;
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
    infbnd = rw[70];
/* definition of an infinite bound */
    argerr = iw[106];
/* maximum # errors in MPS data */
    *iexit = 0;
    *errors = 0;
/* The options haven't been checked yet. */
    if (infbnd <= 0.) {
	infbnd = 1e20;
    }
    if (argerr < 0) {
	argerr = 20;
    }
/* Print 20 lines max */
    nb = *n + *m;
    named = *nname == nb;
/* ================================================================= */
/* Decode 'Start'. */
/* ================================================================= */
    iw[230] = 0;
    iw[231] = 0;
    iw[232] = 0;
/* Determine the type of start. */
/* Preset optional parameters take precedence. */
    if (*lvlsrt != -11111) {
/* lvlSrt was set as an option */
	if (*lvlsrt == 0 || *lvlsrt == 1 || *lvlsrt == 2) {
/* Relax */
	} else if (*lvlsrt == 3) {
	    iw[230] = 1;
	    iw[231] = 1;
	    iw[232] = 1;
	} else {
/* lvlSrt is an unrecognized option */
	    *lvlsrt = -11111;
	}
    }
    if (*lvlsrt == -11111) {
/* lvlSrt is unset */
	*(unsigned char *)ch1 = *(unsigned char *)start;
	if (*(unsigned char *)ch1 == 'C' || *(unsigned char *)ch1 == 'c') {
	    *lvlsrt = 0;
	} else if (*(unsigned char *)ch1 == 'B' || *(unsigned char *)ch1 == 
		'b') {
	    *lvlsrt = 1;
	} else if (*(unsigned char *)ch1 == 'W' || *(unsigned char *)ch1 == 
		'w') {
	    *lvlsrt = 2;
	} else if (*(unsigned char *)ch1 == 'H' || *(unsigned char *)ch1 == 
		'h') {
	    *lvlsrt = 3;
/* nchar  = len_trim(start)       ! An F90 intrinsic */
	    s1trim_(start, &nchar, start_len);
/* Decode    Start = 'HOT ...' */
/* The F77 equivalent */
	    if (nchar <= 4) {
/* 'Hot' or 'Hot ' = 'Hot FHS' */
		iw[230] = 1;
		iw[231] = 1;
		iw[232] = 1;
	    } else {
		i__1 = nchar;
		for (j = 5; j <= i__1; ++j) {
/* Decode 1 or more of FHS */
		    *(unsigned char *)ch1 = *(unsigned char *)&start[j - 1];
		    if (*(unsigned char *)ch1 == 'F' || *(unsigned char *)ch1 
			    == 'f') {
			iw[230] = 1;
		    }
		    if (*(unsigned char *)ch1 == 'H' || *(unsigned char *)ch1 
			    == 'h') {
			iw[231] = 1;
		    }
		    if (*(unsigned char *)ch1 == 'S' || *(unsigned char *)ch1 
			    == 's') {
			iw[232] = 1;
		    }
		}
	    }
	} else {
	    *lvlsrt = 0;
	    s_wsfi(&io___408);
	    do_fio(&c__1, start, start_len);
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/* ================================================================= */
/* Check the other arguments. */
/* ================================================================= */
    if (*m < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___409);
	    do_fio(&c__1, "m     ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*n < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___410);
	    do_fio(&c__1, "n     ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/* 26 Mar 2004: Test if ne > n*m without causing overflow. */
    n1 = max(*n,1);
    k = *ne % n1;
    if (*ne < 1 || k == 0 && *ne / n1 > *m || k > 0 && *ne / n1 >= *m) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___413);
	    do_fio(&c__1, "ne    ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*lencobj < 0 || *lencobj > *n) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___414);
	    do_fio(&c__1, "lencObj", (ftnlen)7);
	    do_fio(&c__1, (char *)&(*lencobj), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*ncolh < 0 || *ncolh > *n) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___415);
	    do_fio(&c__1, "ncolH ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*ncolh), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nname != 1 && *nname != nb) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___416);
	    do_fio(&c__1, "nName ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nname), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*iobj < 0 || *iobj > *m) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___417);
	    do_fio(&c__1, "iObj  ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*iobj), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (loca[1] != 1 || loca[*n + 1] != *ne + 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___418);
	    do_fio(&c__1, (char *)&loca[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&loca[*n + 1], (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    } else {
	nz = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = loca[j + 1] - 1;
	    for (k = loca[j]; k <= i__2; ++k) {
		ir = inda[k];
		++nz;
		if (ir > *m || ir <= 0) {
/* Row index out of range. */
		    ++(*errors);
		    if (*errors <= argerr) {
			s_wsfi(&io___421);
			do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	}
	if (nz != *ne) {
	    ++(*errors);
	    if (*errors <= argerr) {
		s_wsfi(&io___422);
		do_fio(&c__1, "ne    ", (ftnlen)6);
		do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	    }
	}
    }
    mods = 0;
    if (*lvlsrt == 0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    js = hs[j];
	    if (js < -1 || js > 5) {
		++mods;
		hs[j] = 0;
	    }
	}
    } else if (*lvlsrt == 2 || *lvlsrt == 3) {
	i__1 = nb;
	for (j = 1; j <= i__1; ++j) {
	    js = hs[j];
	    if (js < 0 || js > 3) {
		++mods;
		hs[j] = 0;
	    }
	}
    }
    if (mods > 0) {
	s_wsfi(&io___425);
	do_fio(&c__1, (char *)&mods, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------------------------------------------------- */
/* Check the bounds on all variables and constraints. */
/* ----------------------------------------------------------------- */
    errs = 0;
    i__1 = nb;
    for (j = 1; j <= i__1; ++j) {
	b1 = bl[j];
	b2 = bu[j];
	ok = b1 <= b2 && b1 < infbnd && b2 > -infbnd;
	if (! ok) {
	    ++errs;
	    ++(*errors);
	    if (j > *n) {
/* Linear    constraint */
		k = j - *n;
		l = 2;
	    } else {
/* Variable */
		k = j;
		l = 1;
	    }
	    if (*errors <= argerr) {
		if (named) {
		    if (b1 == b2) {
			s_wsfi(&io___431);
			do_fio(&c__1, names + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___432);
			do_fio(&c__1, names + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (b1 == b2) {
			s_wsfi(&io___433);
			do_fio(&c__1, id + (l - 1) * 5, (ftnlen)5);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___434);
			do_fio(&c__1, id + (l - 1) * 5, (ftnlen)5);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	}
    }
    if (errs > 0) {
	s_wsfi(&io___435);
	do_fio(&c__1, (char *)&errs, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    if (*errors >= argerr) {
	snprnt_(&c__3, " and so on ...", &iw[1], leniw, (ftnlen)14);
    }
    if (*errors > 0) {
	*iexit = 91;
/* Invalid arguments for sQOpt. */
    }
    return 0;
} /* s3argq_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3argQ */
/* Subroutine */ int s3prtq_(integer *m, integer *n, integer *ngobj, integer *
	nnh, integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* Initialized data */

    static char pivtyp[24*4] = " LU partial  pivoting..." " LU rook     pivo"
	    "ting..." " LU complete pivoting..." " LU diagonal pivoting...";
    static char prbtyp[24*3] = " Maximize..............." " Feasible point o"
	    "nly...." " Minimize...............";
    static char qptype[24*3] = " QPsolver Cholesky......" " QPsolver CG....."
	    "......." " QPsolver QN............";
    static char srttyp[24*4] = " Cold start............." " Basis file......"
	    "......." " Warm start............." " Hot start..............";
    static char noyes[3*2] = " No" "Yes";

    /* Format strings */
    static char fmt_2110[] = "(\002 Solution file..........\002,i10,6x,\002 "
	    "Old basis file ........\002,i10,6x,\002 Standard input........"
	    ".\002,i10)";
    static char fmt_2120[] = "(\002 Insert file............\002,i10,6x,\002 "
	    "New basis file ........\002,i10,6x,\002 (Printer)............."
	    ".\002,i10)";
    static char fmt_2130[] = "(\002 Punch file.............\002,i10,6x,\002 "
	    "Backup basis file......\002,i10,6x,\002 (Specs file).........."
	    ".\002,i10)";
    static char fmt_2140[] = "(\002 Load file..............\002,i10,6x,\002 "
	    "Dump file..............\002,i10,6x,\002 Standard output......."
	    ".\002,i10)";
    static char fmt_2210[] = "(\002 Print frequency........\002,i10,6x,\002 "
	    "Check frequency........\002,i10,6x,\002 Save new basis map...."
	    ".\002,i10)";
    static char fmt_2220[] = "(\002 Summary frequency......\002,i10,6x,\002 "
	    "Factorization frequency\002,i10,6x,\002 Expand frequency......"
	    ".\002,i10)";
    static char fmt_2310[] = "(a24,16x,a24,16x,a24)";
    static char fmt_2320[] = "(\002 Scale tolerance........\002,0p,f10.3,6x"
	    ",\002 Feasibility tolerance..\002,1p,e10.2,6x,\002 Iteration lim"
	    "it........\002,i10)";
    static char fmt_2330[] = "(\002 Scale option...........\002,i10,6x,\002 "
	    "Optimality tolerance...\002,1p,e10.2,6x,\002 Print level........"
	    "....\002,i10)";
    static char fmt_2340[] = "(\002 Crash tolerance........\002,0p,f10.3,6x"
	    ",\002 Pivot tolerance........\002,1p,e10.2,6x,\002 Partial price"
	    "..........\002,i10)";
    static char fmt_2350[] = "(\002 Crash option...........\002,i10,6x,\002 "
	    "Elastic weight.........\002,1p,e10.2,6x,\002 Prtl price section "
	    "( A)\002,i10)";
    static char fmt_2360[] = "(\002 Elastic mode...........\002,i10,6x,\002 "
	    "Elastic objective......\002,i10,6x,\002 Prtl price section (-I"
	    ")\002,i10)";
    static char fmt_2410[] = "(\002 Objective variables....\002,i10,6x,\002 "
	    "Hessian columns........\002,i10,6x,\002 Superbasics limit....."
	    ".\002,i10)";
    static char fmt_2420[] = "(\002 Nonlin Objective vars..\002,i10,6x,\002 "
	    "Unbounded step size....\002,1p,e10.2)";
    static char fmt_2430[] = "(\002 Linear Objective vars..\002,i10)";
    static char fmt_2380[] = "(\002 Subspace tolerance.....\002,0p,f10.5,6x"
	    ",\002 CG tolerance...........\002,1p,e10.2)";
    static char fmt_2710[] = "(\002 LU factor tolerance....\002,0p,f10.2,6x"
	    ",\002 LU singularity tol.....\002,1p,e10.2,6x,\002 Timing level."
	    "..........\002,i10)";
    static char fmt_2720[] = "(\002 LU update tolerance....\002,0p,f10.2,6x"
	    ",\002 LU swap tolerance......\002,1p,e10.2,6x,\002 Debug level.."
	    "..........\002,i10)";
    static char fmt_2730[] = "(a24,16x,\002 eps (machine precision)\002,1p,e"
	    "10.2,6x,\002 System information.....\002,7x,a3)";
    static char fmt_2740[] = "(80x,\002 Sticky parameters......\002,7x,a3)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static logical qp;
    static doublereal eps;
    static integer nqp, npr1, npr2;
    static char str1[132], str2[132], str3[132], str4[132], str5[132], str6[
	    132];
    static integer kfac, kchk, klog, ksav, maxr, maxs;
    static doublereal tolx, utol1;
    static integer iback, ioldb;
    static doublereal bigdx;
    static integer ipnch;
    static doublereal etarg;
    static integer inewb;
    static doublereal tolcg;
    static integer istdi, idump, maxmn, istdo, isoln, ksumm;
    static doublereal tolqp;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *);
    static doublereal wtinf0;
    static integer iloadb, kdegen, lemode;
    static doublereal tolfac;
    static integer icrash, lprdbg;
    static doublereal tcrash;
    static integer ispecs, minprc, minmax, lvlinf, itnlim, lvlscl;
    static doublereal scltol;
    static integer nparpr, iprint, lvltim, iinsrt;
    static doublereal tolupd;
    static integer lprprm, lvlpiv;
    static doublereal tolpiv;
    static integer mnrprt;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer lvlsrt, qpslvr;
    static doublereal tolswp;
    static integer stkyop, lvlsys;

    /* Fortran I/O blocks */
    static icilist io___498 = { 0, str1, 0, fmt_2110, 132, 1 };
    static icilist io___500 = { 0, str2, 0, fmt_2120, 132, 1 };
    static icilist io___502 = { 0, str3, 0, fmt_2130, 132, 1 };
    static icilist io___504 = { 0, str4, 0, fmt_2140, 132, 1 };
    static icilist io___505 = { 0, str1, 0, fmt_2210, 132, 1 };
    static icilist io___506 = { 0, str2, 0, fmt_2220, 132, 1 };
    static icilist io___507 = { 0, str1, 0, fmt_2310, 132, 1 };
    static icilist io___508 = { 0, str2, 0, fmt_2320, 132, 1 };
    static icilist io___509 = { 0, str3, 0, fmt_2330, 132, 1 };
    static icilist io___510 = { 0, str4, 0, fmt_2340, 132, 1 };
    static icilist io___512 = { 0, str5, 0, fmt_2350, 132, 1 };
    static icilist io___514 = { 0, str6, 0, fmt_2360, 132, 1 };
    static icilist io___515 = { 0, str1, 0, fmt_2410, 132, 1 };
    static icilist io___516 = { 0, str2, 0, fmt_2420, 132, 1 };
    static icilist io___517 = { 0, str3, 0, fmt_2430, 132, 1 };
    static icilist io___518 = { 0, str1, 0, fmt_2380, 132, 1 };
    static icilist io___519 = { 0, str1, 0, fmt_2710, 132, 1 };
    static icilist io___520 = { 0, str2, 0, fmt_2720, 132, 1 };
    static icilist io___521 = { 0, str3, 0, fmt_2730, 132, 1 };
    static icilist io___522 = { 0, str4, 0, fmt_2740, 132, 1 };


/* ================================================================= */
/* s3prtQ prints the settings of the optional parameters for the */
/* standard SNOPT wrapper. */

/* See  snworkspace.info  for full documentation of cw, iw and rw. */

/* 15 Nov 1991: First version of s3prtQ (s8prnt). */
/* 10 Dec 2002: More LU pivoting options. */
/* 01 Jul 2003: QN QP options added. */
/* 30 Jul 2003: QN CG options added. */
/* 03 Aug 2003: snPRNT adopted. */
/* 01 Sep 2007: Sticky parameters added. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
/* Set some local machine-dependent constants. */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    tolqp = rw[52];
/* Minor Phase 2 Opt tol */
    tolcg = rw[54];
/* cg tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance */
    tolpiv = rw[60];
/* excludes small elements of y */
    tcrash = rw[62];
/* crash tolerance */
    tolswp = rw[65];
/* LU swap tolerance */
    tolfac = rw[66];
/* LU factor tolerance */
    tolupd = rw[67];
/* LU update tolerance */
    bigdx = rw[72];
/* unbounded step */
    etarg = rw[83];
/* Quasi-Newton QP rg tolerance */
    wtinf0 = rw[88];
/* infeasibility weight */
    scltol = rw[92];
/* scale tolerance. */
    utol1 = rw[154];
/* abs tol for small diag of U */
    istdi = iw[9];
/* Standard Input */
    istdo = iw[10];
/* Standard Output */
    ispecs = iw[11];
/* Specs (options) file */
    iprint = iw[12];
/* Print file */
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    qpslvr = iw[55];
/* = 0:1:2   => QPChol:CG:QN QP solver */
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
    lvlsrt = iw[69];
/* = 0:1:2:3 => cold:warm:basis:hot start */
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
    lprdbg = iw[85];
/* > 0    => private debug print */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    icrash = iw[88];
/* Crash option */
    itnlim = iw[89];
/* limit on total iterations */
    mnrprt = iw[93];
/* Minor print level */
    nparpr = iw[94];
/* # of partial pricing sections */
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
    isoln = iw[131];
/* Solution file */
    lvltim = iw[182];
/* ----------------------------------------------------------------- */
/* Timing level */
    if (lprprm == 0 || mnrprt == 0) {
	return 0;
    }
    qp = *nnh > 0;
    nqp = max(*ngobj,*nnh);
    minprc = 10;
    npr1 = *n / nparpr;
    npr2 = *m / nparpr;
    if (max(npr1,npr2) < minprc) {
	maxmn = max(*m,*n);
	nparpr = maxmn / min(maxmn,minprc);
	npr1 = *n / nparpr;
	npr2 = *m / nparpr;
    }
/* ================================================================= */
/* Print parameters except if PRINT LEVEL = 0 */
/* or SUPPRESS PARAMETERS was specified. */
/* ================================================================= */
    s1page_(&c__1, &iw[1], leniw);
    snprnt_(&c__1, " Parameters", &iw[1], leniw, (ftnlen)11);
    snprnt_(&c__1, " ==========", &iw[1], leniw, (ftnlen)11);
/* ----------------------- */
/* Files. */
/* ----------------------- */
    snprnt_(&c__11, " Files", &iw[1], leniw, (ftnlen)6);
    snprnt_(&c__1, " -----", &iw[1], leniw, (ftnlen)6);
    s_wsfi(&io___498);
    do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ioldb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&istdi, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___500);
    do_fio(&c__1, (char *)&iinsrt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inewb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iprint, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___502);
    do_fio(&c__1, (char *)&ipnch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&iback, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ispecs, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___504);
    do_fio(&c__1, (char *)&iloadb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&idump, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&istdo, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* Frequencies. */
/* ----------------------- */
    snprnt_(&c__11, " Frequencies", &iw[1], leniw, (ftnlen)12);
    snprnt_(&c__1, " -----------", &iw[1], leniw, (ftnlen)12);
    s_wsfi(&io___505);
    do_fio(&c__1, (char *)&klog, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kchk, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ksav, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___506);
    do_fio(&c__1, (char *)&ksumm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kfac, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&kdegen, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* LP/QP parameters. */
/* ----------------------- */
    snprnt_(&c__11, " LP/QP Parameters", &iw[1], leniw, (ftnlen)17);
    snprnt_(&c__1, " ----------------", &iw[1], leniw, (ftnlen)17);
    s_wsfi(&io___507);
    do_fio(&c__1, prbtyp + (minmax + 1) * 24, (ftnlen)24);
    do_fio(&c__1, qptype + qpslvr * 24, (ftnlen)24);
    do_fio(&c__1, srttyp + lvlsrt * 24, (ftnlen)24);
    e_wsfi();
    s_wsfi(&io___508);
    do_fio(&c__1, (char *)&scltol, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&itnlim, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___509);
    do_fio(&c__1, (char *)&lvlscl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&tolqp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&mnrprt, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___510);
    do_fio(&c__1, (char *)&tcrash, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolpiv, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nparpr, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___512);
    do_fio(&c__1, (char *)&icrash, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&wtinf0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&npr1, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___514);
    do_fio(&c__1, (char *)&lemode, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&lvlinf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npr2, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str5, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str6, &iw[1], leniw, (ftnlen)132);
/* ----------------------- */
/* QP objective */
/* ----------------------- */
    if (qp) {
	snprnt_(&c__11, " QP objective", &iw[1], leniw, (ftnlen)13);
	snprnt_(&c__1, " ------------", &iw[1], leniw, (ftnlen)13);
	s_wsfi(&io___515);
	do_fio(&c__1, (char *)&nqp, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nnh), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&maxs, (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___516);
	do_fio(&c__1, (char *)&(*nnh), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&bigdx, (ftnlen)sizeof(doublereal));
	e_wsfi();
	s_wsfi(&io___517);
	do_fio(&c__1, (char *)&(*ngobj), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------- */
/* Quasi-Newton QP solver. */
/* ----------------------- */
    if (qpslvr == 1 || maxr < maxs) {
	snprnt_(&c__11, " Conjugate-gradient QP solver", &iw[1], leniw, (
		ftnlen)29);
	snprnt_(&c__1, " ----------------------------", &iw[1], leniw, (
		ftnlen)29);
	s_wsfi(&io___518);
	do_fio(&c__1, (char *)&etarg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&tolcg, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------- */
/* Miscellaneous */
/* ----------------------- */
    snprnt_(&c__11, " Miscellaneous", &iw[1], leniw, (ftnlen)14);
    snprnt_(&c__1, " -------------", &iw[1], leniw, (ftnlen)14);
    s_wsfi(&io___519);
    do_fio(&c__1, (char *)&tolfac, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&utol1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lvltim, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___520);
    do_fio(&c__1, (char *)&tolupd, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tolswp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&lprdbg, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___521);
    do_fio(&c__1, pivtyp + lvlpiv * 24, (ftnlen)24);
    do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, noyes + lvlsys * 3, (ftnlen)3);
    e_wsfi();
    s_wsfi(&io___522);
    do_fio(&c__1, noyes + stkyop * 3, (ftnlen)3);
    e_wsfi();
    snprnt_(&c__1, str1, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str2, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str3, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, str4, &iw[1], leniw, (ftnlen)132);
    return 0;
} /* s3prtq_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3prtQ */
/* Subroutine */ int s3argn_(integer *iexit, char *start, integer *lda, 
	integer *ldcj, integer *ldh, integer *n, integer *nclin, integer *
	ncnln, integer *nname, doublereal *bl, doublereal *bu, char *names, 
	integer *istate, doublereal *cmul, integer *lvlsrt, integer *errors, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	start_len, ftnlen names_len)
{
    /* Initialized data */

    static char id[5*3] = "varbl" "lncon" "nlcon";

    /* Format strings */
    static char fmt_1000[] = "(\002 XXX Start parameter not recognized:  "
	    "\002,a)";
    static char fmt_1100[] = "(\002 XXX  Argument out of range:  \002,a6,"
	    "\002 = \002,i6)";
    static char fmt_1400[] = "(\002 XXX  Invalid argument iState: \002,i6"
	    ",\002 elements modified to be in range.\002)";
    static char fmt_1510[] = "(\002 XXX  The equal bounds on  \002,a8,\002  "
	    "are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16.7)";
    static char fmt_1515[] = "(\002 XXX  The bounds on  \002,a8,\002  are in"
	    "consistent.   bl =\002,g16.7,\002   bu =\002,g16.7)";
    static char fmt_1500[] = "(\002 XXX  The equal bounds on  \002,a5,i6,"
	    "\002  are infinite.   Bounds =\002,g16.7,\002  infBnd =\002,g16."
	    "7)";
    static char fmt_1505[] = "(\002 XXX  The bounds on  \002,a5,i6,\002  are"
	    " inconsistent.   bl =\002,g16.7,\002   bu =\002,g16.7)";
    static char fmt_1600[] = "(\002 XXX  Invalid arguments bl, bu: \002,i6"
	    ",\002 inconsistent bounds or infinite equal bounds.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal b1, b2;
    static integer nb;
    static logical ok;
    static integer is;
    static char ch1[1];
    static doublereal mul;
    static char str[132];
    static integer mods, errs;
    static logical named;
    static integer nchar;
    extern /* Subroutine */ int s1trim_(char *, integer *, ftnlen);
    static doublereal infbnd;
    static integer argerr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___532 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___533 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___534 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___535 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___536 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___537 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___538 = { 0, str, 0, fmt_1100, 132, 1 };
    static icilist io___543 = { 0, str, 0, fmt_1400, 132, 1 };
    static icilist io___550 = { 0, str, 0, fmt_1510, 132, 1 };
    static icilist io___551 = { 0, str, 0, fmt_1515, 132, 1 };
    static icilist io___552 = { 0, str, 0, fmt_1500, 132, 1 };
    static icilist io___553 = { 0, str, 0, fmt_1505, 132, 1 };
    static icilist io___554 = { 0, str, 0, fmt_1600, 132, 1 };


/* ================================================================= */
/* s3argN   checks the arguments for npOpt. */

/* On exit, Errors says how many errors were encountered. */

/* 17 Mar 2002: First version of s3argN. */
/* 17 Jun 2004: Current version of s3argN. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --cmul;
    --istate;
    --bu;
    --bl;
    names -= 8;
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
    infbnd = rw[70];
/* definition of an infinite bound */
    argerr = iw[106];
/* maximum # errors in MPS data */
    *iexit = 0;
    *errors = 0;
/* The options haven't been checked yet. */
    if (infbnd < 0.) {
	infbnd = 1e20;
    }
    if (argerr < 0) {
	argerr = 20;
    }
/* Print 20 lines max */
    nb = *n + *nclin + *ncnln;
    named = *nname == nb;
/* ================================================================= */
/* Decode 'Start'. */
/* ================================================================= */
    iw[230] = 0;
    iw[231] = 0;
    iw[232] = 0;
/* Determine the type of start. */
/* Optional parameters take precedence. */
    if (*lvlsrt != -11111) {
/* lvlSrt was set as an option */
	if (*lvlsrt == 0 || *lvlsrt == 1 || *lvlsrt == 2) {
/*           Relax */
	} else if (*lvlsrt == 3) {
	    iw[230] = 1;
	    iw[231] = 1;
	    iw[232] = 1;
	} else {
/* lvlSrt is an unrecognized option */
	    *lvlsrt = -11111;
	}
    }
    if (*lvlsrt == -11111) {
/* lvlSrt is unset */
	*(unsigned char *)ch1 = *(unsigned char *)start;
	if (*(unsigned char *)ch1 == 'C' || *(unsigned char *)ch1 == 'c') {
	    *lvlsrt = 0;
	} else if (*(unsigned char *)ch1 == 'B' || *(unsigned char *)ch1 == 
		'b') {
	    *lvlsrt = 1;
	} else if (*(unsigned char *)ch1 == 'W' || *(unsigned char *)ch1 == 
		'w') {
	    *lvlsrt = 2;
	} else if (*(unsigned char *)ch1 == 'H' || *(unsigned char *)ch1 == 
		'h') {
	    *lvlsrt = 3;
/* nchar  = len_trim(Start)       ! An F90 intrinsic */
	    s1trim_(start, &nchar, start_len);
/* Decode    Start = 'HOT ...' */
/* The F77 equivalent */
	    if (nchar <= 4) {
/* 'Hot' or 'Hot ' = 'Hot FHS' */
		iw[230] = 1;
		iw[231] = 1;
		iw[232] = 1;
	    } else {
		i__1 = nchar;
		for (j = 5; j <= i__1; ++j) {
/* Decode 1 or more of FHS */
		    *(unsigned char *)ch1 = *(unsigned char *)&start[j - 1];
		    if (*(unsigned char *)ch1 == 'F' || *(unsigned char *)ch1 
			    == 'f') {
			iw[230] = 1;
		    }
		    if (*(unsigned char *)ch1 == 'H' || *(unsigned char *)ch1 
			    == 'h') {
			iw[231] = 1;
		    }
		    if (*(unsigned char *)ch1 == 'S' || *(unsigned char *)ch1 
			    == 's') {
			iw[232] = 1;
		    }
		}
	    }
	} else {
	    *lvlsrt = 0;
	    s_wsfi(&io___532);
	    do_fio(&c__1, start, start_len);
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
/* ================================================================= */
/* Check the other arguments. */
/* ================================================================= */
    if (*n < 1) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___533);
	    do_fio(&c__1, "n     ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*nclin < 0) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___534);
	    do_fio(&c__1, "nclin ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*nclin), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*ncnln < 0) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___535);
	    do_fio(&c__1, "ncnln ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*ncnln), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*lda < 0 || *nclin > *lda) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___536);
	    do_fio(&c__1, "ldA   ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*lda), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*ldcj < 0 || *ncnln > *ldcj) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___537);
	    do_fio(&c__1, "ldcJ  ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*ldcj), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*ldh < 0 || *n > *ldh) {
	++(*errors);
	if (*errors <= argerr) {
	    s_wsfi(&io___538);
	    do_fio(&c__1, "ldH   ", (ftnlen)6);
	    do_fio(&c__1, (char *)&(*ldh), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    mods = 0;
    if (*lvlsrt == 2 || *lvlsrt == 3) {
	i__1 = nb;
	for (j = 1; j <= i__1; ++j) {
	    is = istate[j];
	    if (is < -2 || is > 4) {
		++mods;
		istate[j] = 0;
	    }
	}
	i__1 = *ncnln;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = *n + *nclin + i__;
	    is = istate[j];
	    mul = cmul[j];
	    if (is == 0) {
		mul = 0.;
	    } else if (is == 1) {
		if (bl[j] <= -infbnd) {
		    is = 0;
		}
		if (mul < 0. || is == 0) {
		    mul = 0.;
		}
	    } else if (is == 2) {
		if (bu[j] >= infbnd) {
		    is = 0;
		}
		if (mul > 0. || is == 0) {
		    mul = 0.;
		}
	    } else if (is == 3) {
		if (bl[j] < bu[j]) {
		    is = 0;
		}
	    }
	    istate[j] = is;
	    cmul[j] = mul;
	}
    }
    if (mods > 0) {
	s_wsfi(&io___543);
	do_fio(&c__1, (char *)&mods, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
/* ----------------------------------------------------------------- */
/* Check the bounds on all variables and constraints. */
/* ----------------------------------------------------------------- */
    errs = 0;
    i__1 = nb;
    for (j = 1; j <= i__1; ++j) {
	b1 = bl[j];
	b2 = bu[j];
	ok = b1 <= b2 && b1 < infbnd && b2 > -infbnd;
	if (! ok) {
	    ++errs;
	    ++(*errors);
	    if (j > *n + *nclin) {
		k = j - *n - *nclin;
		l = 3;
	    } else if (j > *n) {
		k = j - *n;
		l = 2;
	    } else {
		k = j;
		l = 1;
	    }
	    if (*errors <= argerr) {
		if (named) {
		    if (b1 == b2) {
			s_wsfi(&io___550);
			do_fio(&c__1, names + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___551);
			do_fio(&c__1, names + (j << 3), (ftnlen)8);
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (b1 == b2) {
			s_wsfi(&io___552);
			do_fio(&c__1, id + (l - 1) * 5, (ftnlen)5);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&infbnd, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    } else {
			s_wsfi(&io___553);
			do_fio(&c__1, id + (l - 1) * 5, (ftnlen)5);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal))
				;
			e_wsfi();
			snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	}
    }
    if (errs > 0) {
	s_wsfi(&io___554);
	do_fio(&c__1, (char *)&errs, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    if (*errors >= argerr) {
	snprnt_(&c__3, " and so on ...", &iw[1], leniw, (ftnlen)14);
    }
    if (*errors > 0) {
	*iexit = 91;
/* Invalid arguments for npOptN. */
    }
    return 0;
} /* s3argn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3argN */
/* Subroutine */ int s3hesn_(integer *task, integer *ldh, integer *lenh, 
	integer *n, doublereal *h__, doublereal *hess)
{
    /* System generated locals */
    integer hess_dim1, hess_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l;

/* ================================================================= */
/* s3HesN loads the problem into SNOPT format. */

/* 07 Jul 1998: First version of s3HesN. */
/* 04 Nov 2000: Current version. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    hess_dim1 = *ldh;
    hess_offset = 1 + hess_dim1;
    hess -= hess_offset;
    --h__;

    /* Function Body */
    if (*task == 0) {
/* -------------------------------------------------------------- */
/* Load the user-supplied Hessian Hess into H. */
/* -------------------------------------------------------------- */
	l = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = i__; j <= i__2; ++j) {
		++l;
		h__[l] = hess[i__ + j * hess_dim1];
	    }
	}
    } else if (*task == 1) {
/* -------------------------------------------------------------- */
/* Down load the SNOPT approximate Hessian into Hess. */
/* -------------------------------------------------------------- */
	l = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = i__; j <= i__2; ++j) {
		++l;
		hess[i__ + j * hess_dim1] = h__[l];
		hess[j + i__ * hess_dim1] = h__[l];
	    }
	}
    }
    return 0;
} /* s3hesn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3HesN */
/* Subroutine */ int s3inin_(integer *start, integer *n, integer *nb, integer 
	*nncon0, integer *nncon, integer *negcon, integer *hs, doublereal *
	fcon, doublereal *gcon, doublereal *gobj, doublereal *rc, doublereal *
	x, integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, lfh;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lenfh, lvlhes;

/* ================================================================= */
/* s3iniN initializes SNOPT variables that are ultimately copied to */
/* NPSOL format in s3outN. */

/* 04 Dec 2004: First version of s3iniN. */
/* 04 Dec 2004: Current version of s3iniN. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --gobj;
    --x;
    --rc;
    --hs;
    --fcon;
    --gcon;
    --iw;
    --rw;

    /* Function Body */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lfh = iw[391];
/* H(lenfH), full-memory BFGS Hessian */
    lenfh = iw[392];

    if (*start == 0) {
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    hs[j] = 0;
	    x[j] = 0.;
	    rc[j] = 0.;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    gobj[j] = 0.;
	}
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fcon[i__] = 0.;
	}
	i__1 = *nncon * *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    gcon[i__] = 0.;
	}
	if (lvlhes == 1) {
	    dload_(&lenfh, &c_b367, &rw[lfh], &c__1);
	}
    }
    return 0;
} /* s3inin_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3iniN */
/* Subroutine */ int s3inn_(integer *start, integer *lda, integer *ldh, 
	integer *m, integer *n, integer *nclin, integer *ncon, integer *nncol,
	 integer *nb, integer *nncon0, integer *nncon, integer *hs, integer *
	istate, doublereal *alin, integer *ne, integer *nlocj, integer *locj, 
	integer *indj, doublereal *jcol, doublereal *bl, doublereal *bu, 
	doublereal *bbl, doublereal *bbu, doublereal *c__, doublereal *cmul, 
	doublereal *hess, doublereal *pi, doublereal *x, doublereal *xs, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer alin_dim1, alin_offset, hess_dim1, hess_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, ij, is, js;
    static doublereal xj;
    static integer lh0, lfh;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lenfh;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), s3hesn_(
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static doublereal infbnd;
    static integer lvlhes;

/* ================================================================= */
/* s3inN loads a problem in NPSOL format into SNOPT format. */

/* 22 Mar 1997: First version of s3inN. */
/* 04 Jan 2001: Current version. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    alin_dim1 = *lda;
    alin_offset = 1 + alin_dim1;
    alin -= alin_offset;
    hess_dim1 = *ldh;
    hess_offset = 1 + hess_dim1;
    hess -= hess_offset;
    --x;
    --cmul;
    --bu;
    --bl;
    --istate;
    --xs;
    --bbu;
    --bbl;
    --hs;
    --pi;
    --c__;
    --jcol;
    --indj;
    --locj;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lh0 = iw[346];
/* ================================================================= */
/* Load the snopt arrays. */
/* Copy the bounds, x's first, linears next, then nonlinears. */
/* ================================================================= */
/* Initial diagonal Hessian */
    dcopy_(n, &bl[1], &c__1, &bbl[1], &c__1);
    dcopy_(n, &bu[1], &c__1, &bbu[1], &c__1);
    if (*nclin > 0) {
	dcopy_(nclin, &bl[*n + 1], &c__1, &bbl[*n + *nncon + 1], &c__1);
	dcopy_(nclin, &bu[*n + 1], &c__1, &bbu[*n + *nncon + 1], &c__1);
    }
    if (*nncon > 0) {
	dcopy_(nncon, &bl[*n + *nclin + 1], &c__1, &bbl[*n + 1], &c__1);
	dcopy_(nncon, &bu[*n + *nclin + 1], &c__1, &bbu[*n + 1], &c__1);
    }
    if (*ncon == 0) {
	bbl[*nb] = -infbnd;
	bbu[*nb] = infbnd;
    }
    if (*start == 0) {
/* -------------------------------------------- */
/* Cold Start. */
/* -------------------------------------------- */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    xj = x[j];
	    if (xj <= bl[j]) {
		hs[j] = 4;
	    } else if (xj >= bu[j]) {
		hs[j] = 5;
	    } else {
		hs[j] = 0;
	    }
	    xs[j] = xj;
	}
	if (*nncon > 0) {
	    dload_(nncon, &c_b367, &pi[1], &c__1);
	}
    } else if (*start == 2) {
/* ----------------------------------------------- */
/* Warm Start. */
/* Input values of x, cMul, hs and Hess are used. */
/* The use of Hess is unique to snoptn. */
/* ----------------------------------------------- */
	dcopy_(n, &x[1], &c__1, &xs[1], &c__1);
	if (*nncon > 0) {
	    dcopy_(nncon, &c__[1], &c__1, &xs[*n + 1], &c__1);
	    dcopy_(nncon, &cmul[*n + *nclin + 1], &c__1, &pi[1], &c__1);
	}
	if (*nclin > 0) {
	    dload_(nclin, &c_b367, &xs[*n + *nncon + 1], &c__1);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(nclin, &xs[j], &alin[j * alin_dim1 + 1], &c__1, &xs[*n 
			+ *nncon + 1], &c__1);
	    }
	}
	l = 1;
	i__1 = *n + *ncon;
	for (j = 1; j <= i__1; ++j) {
	    is = istate[j];
	    js = 0;
	    if (is == 1) {
		js = 0;
	    } else if (is == 2) {
		js = 1;
	    }
	    hs[l] = js;
	    if (j == *n && *nclin > 0) {
		l = *n + *nncon;
	    } else if (j == *n + *nclin) {
		l = *n + 1;
	    } else {
		++l;
	    }
	}
	if (lvlhes == 0) {
	    i__1 = *ldh + 1;
	    dcopy_(n, &hess[hess_dim1 + 1], &i__1, &rw[lh0], &c__1);
	} else if (lvlhes == 1) {
	    lfh = iw[391];
/* H(lenfH), full-memory BFGS Hessian */
	    lenfh = iw[392];

	    s3hesn_(&c__0, ldh, &lenfh, n, &rw[lfh], &hess[hess_offset]);
	}
    }
/* -------------------------------------------------------------- */
/* Load the linear part of A with the linear constraints. */
/* -------------------------------------------------------------- */
/* cold start */
    if (*nncol == 0) {
/* Sparse dummy row */
	jcol[1] = 0.;
	indj[1] = 1;
	locj[1] = 1;
	i__1 = *n + 1;
	for (j = 2; j <= i__1; ++j) {
	    locj[j] = 2;
	}
    } else {
	ij = 1;
	locj[1] = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		indj[ij] = i__;
		if (i__ <= *nncon) {
		    jcol[ij] = 0.;
		} else if (i__ <= *ncon) {
		    jcol[ij] = alin[i__ - *nncon + j * alin_dim1];
		}
		++ij;
	    }
	    locj[j + 1] = ij;
	}
    }
    return 0;
} /* s3inn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3inN */
/* Subroutine */ int s3outn_(integer *ldcj, integer *ldh, integer *n, integer 
	*nclin, integer *ncon, integer *nb, integer *nncon0, integer *nncon, 
	integer *hs, integer *istate, doublereal *c__, doublereal *cjac, 
	doublereal *cmul, doublereal *fcon, doublereal *gcon, doublereal *
	gobj, doublereal *grad, doublereal *hess, doublereal *rc, doublereal *
	x, doublereal *xs, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* System generated locals */
    integer cjac_dim1, cjac_offset, gcon_dim1, gcon_offset, hess_dim1, 
	    hess_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, is, js, lfh, lenfh;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s3hesn_(integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *);
    static integer lvlhes;

/* ================================================================= */
/* s3outN changes the problem from SNOPT to NPSOL format. */

/* 22 Mar 1997: First version of s3outN. */
/* 04 Jan 2001: Current version. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    cjac_dim1 = *ldcj;
    cjac_offset = 1 + cjac_dim1;
    cjac -= cjac_offset;
    hess_dim1 = *ldh;
    hess_offset = 1 + hess_dim1;
    hess -= hess_offset;
    --x;
    --grad;
    --gobj;
    --cmul;
    --istate;
    --xs;
    --rc;
    --hs;
    gcon_dim1 = *nncon0;
    gcon_offset = 1 + gcon_dim1;
    gcon -= gcon_offset;
    --fcon;
    --c__;
    --iw;
    --rw;

    /* Function Body */
    lvlhes = iw[72];
/* ================================================================= */
/* Unload the SNOPT solution into the snoptn arrays */
/* Copy gCon, gObj into cJac and grad, */
/* ================================================================= */
/* 0,1,2  => LM, FM, Exact Hessian */
    l = 1;
    i__1 = *n + *ncon;
    for (j = 1; j <= i__1; ++j) {
	js = hs[l];
	is = 0;
	if (js == 0) {
	    is = 1;
	} else if (js == 1) {
	    is = 2;
	}
	istate[j] = is;
	if (j == *n && *nclin > 0) {
	    l = *n + *nncon;
	} else if (j == *n + *nclin) {
	    l = *n + 1;
	} else {
	    ++l;
	}
    }
/* ----------------------------------------------------------------- */
/* Copy gCon, gObj into cJac and grad, */
/* ----------------------------------------------------------------- */
    dcopy_(n, &xs[1], &c__1, &x[1], &c__1);
    dcopy_(n, &gobj[1], &c__1, &grad[1], &c__1);
    dcopy_(n, &rc[1], &c__1, &cmul[1], &c__1);
    if (*nclin > 0) {
	dcopy_(nclin, &rc[*n + *nncon + 1], &c__1, &cmul[*n + 1], &c__1);
    }
    if (*nncon > 0) {
	dcopy_(nncon, &fcon[1], &c__1, &c__[1], &c__1);
	dcopy_(nncon, &rc[*n + 1], &c__1, &cmul[*n + *nclin + 1], &c__1);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nncon;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		cjac[i__ + j * cjac_dim1] = gcon[i__ + j * gcon_dim1];
	    }
	}
    }
    if (lvlhes == 1) {
	lfh = iw[391];
/* H(lenfH), full-memory BFGS Hessian */
	lenfh = iw[392];

	s3hesn_(&c__1, ldh, &lenfh, n, &rw[lfh], &hess[hess_offset]);
    }
    return 0;
} /* s3outn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3outN */
/* Subroutine */ int s3prtn_(integer *n, integer *nb, integer *nclin, integer 
	*nncon0, integer *lda, integer *lprsol, doublereal *xnorm, integer *
	istate, doublereal *a, doublereal *bl, doublereal *bu, doublereal *
	c__, doublereal *cmul, doublereal *x, doublereal *r__, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* Initialized data */

    static char lstate[2*7] = "--" "++" "FR" "LL" "UL" "EQ" "TF";

    /* Format strings */
    static char fmt_1000[] = "(1x,a15,2x,\002State\002,6x,\002Value\002,7x"
	    ",\002Lower bound\002,5x,\002Upper bound\002,3x,\002Lagr multipli"
	    "er\002,4x,\002   Slack\002)";
    static char fmt_2000[] = "(1x,a8,i6,3x,a1,1x,a2,4g16.7,g16.4)";

    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal b1, b2;
    static integer is;
    static doublereal rj;
    static char key[1];
    static doublereal slk, tol;
    static char str[132];
    static doublereal slk1, slk2;
    static char name__[8], line[132];
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal wlam, tolx;
    static integer nplin;
    static char state[2];
    static doublereal infbnd;
    static integer number, iprint;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___590 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___603 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___604 = { 0, str, 0, fmt_1000, 132, 1 };
    static icilist io___608 = { 0, line, 0, fmt_2000, 132, 1 };


/* ================================================================= */
/* s3prtN  prints  x,  A*x, c(x), the bounds, the */
/* multipliers, and the slacks (distance to the nearer bound). */

/* 22 Mar 1997: First version of s3prtN. */
/* 03 Aug 2003: snPRNT adopted. */
/* 22 Nov 2003: Current version of s3prtN. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --x;
    --r__;
    --cmul;
    --bu;
    --bl;
    --istate;
    --c__;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --iw;
    --rw;

    /* Function Body */
/* ----------------------------------------------------------------- */
    iprint = iw[12];
/* Print file */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    infbnd = rw[70];
/* definition of an infinite bound */
    if (iprint == 0 || *lprsol == 0) {
	return 0;
    }
    nplin = *n + *nclin;
    tol = tolx;
    s_wsfi(&io___590);
    do_fio(&c__1, "Variable       ", (ftnlen)15);
    e_wsfi();
    snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    s_copy(name__, "variable", (ftnlen)8, (ftnlen)8);
    nplin = *n + *nclin;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	b1 = bl[j];
	b2 = bu[j];
	wlam = cmul[j];
	if (j <= *n) {
	    rj = x[j];
	    tol = tolx;
	} else {
	    tol = tolx * *xnorm;
	    if (j <= nplin) {
		i__ = j - *n;
		rj = ddot_(n, &a[i__ + a_dim1], lda, &x[1], &c__1);
	    } else {
		i__ = j - nplin;
		rj = c__[i__];
	    }
	}
	slk1 = rj - b1;
	slk2 = b2 - rj;
	r__[j] = rj;
/* Reset istate if necessary. */
	is = istate[j];
	if (slk1 < -tol) {
	    is = -2;
	}
	if (slk2 < -tol) {
	    is = -1;
	}
	if (is == 1 && slk1 > tol) {
	    is = 0;
	}
	if (is == 2 && slk2 > tol) {
	    is = 0;
	}
	istate[j] = is;
	s_copy(state, lstate + (is + 2 << 1), (ftnlen)2, (ftnlen)2);
	if (j <= *n) {
	    number = j;
	} else if (j <= nplin) {
	    number = j - *n;
	    if (number == 1) {
		s_wsfi(&io___603);
		do_fio(&c__1, "Linear constrnt", (ftnlen)15);
		e_wsfi();
		snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
		snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
		s_copy(name__, "lincon  ", (ftnlen)8, (ftnlen)8);
	    }
	} else {
	    number = j - nplin;
	    if (number == 1) {
		s_wsfi(&io___604);
		do_fio(&c__1, "Nonlin constrnt", (ftnlen)15);
		e_wsfi();
		snprnt_(&c__11, " ", &iw[1], leniw, (ftnlen)1);
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
		snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
		s_copy(name__, "nlncon  ", (ftnlen)8, (ftnlen)8);
	    }
	}
/* ------------------------------------------------- */
/* Print a line for the jth variable or constraint. */
/* ------------------------------------------------- */
	if (abs(slk1) < abs(slk2)) {
	    slk = slk1;
	    if (b1 <= -infbnd) {
		slk = slk2;
	    }
	} else {
	    slk = slk2;
	    if (b2 >= infbnd) {
		slk = slk1;
	    }
	}
/* Flag infeasibilities, primal and dual degeneracies, */
/* and active QP constraints that are loose in NP. */
	*(unsigned char *)key = ' ';
	if (slk1 < -tol || slk2 < -tol) {
	    *(unsigned char *)key = 'I';
	}
	if (is == 0 && abs(slk) <= tol) {
	    *(unsigned char *)key = 'D';
	}
	if (is >= 1 && abs(wlam) <= tol) {
	    *(unsigned char *)key = 'A';
	}
	s_wsfi(&io___608);
	do_fio(&c__1, name__, (ftnlen)8);
	do_fio(&c__1, (char *)&number, (ftnlen)sizeof(integer));
	do_fio(&c__1, key, (ftnlen)1);
	do_fio(&c__1, state, (ftnlen)2);
	do_fio(&c__1, (char *)&rj, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&wlam, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&slk, (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Reset special cases: */
/*  Infinite bounds */
/*  Zero bounds */
/*  Lagrange multipliers for inactive constraints */
/*  Lagrange multipliers for infinite bounds */
/*  Infinite slacks */
/*  Zero slacks */
	if (b1 <= -infbnd) {
	    s_copy(line + 38, "      None      ", (ftnlen)16, (ftnlen)16);
	}
	if (b2 >= infbnd) {
	    s_copy(line + 54, "      None      ", (ftnlen)16, (ftnlen)16);
	}
	if (b1 == 0.) {
	    s_copy(line + 38, "        .       ", (ftnlen)16, (ftnlen)16);
	}
	if (b2 == 0.) {
	    s_copy(line + 54, "        .       ", (ftnlen)16, (ftnlen)16);
	}
	if (is == 0 || wlam == 0.) {
	    s_copy(line + 70, "        .       ", (ftnlen)16, (ftnlen)16);
	}
	if (b1 <= -infbnd && b2 >= infbnd) {
	    s_copy(line + 70, "                ", (ftnlen)16, (ftnlen)16);
	    s_copy(line + 86, "                ", (ftnlen)16, (ftnlen)16);
	}
	if (slk == 0.) {
	    s_copy(line + 86, "        .       ", (ftnlen)16, (ftnlen)16);
	}
	snprnt_(&c__1, line, &iw[1], leniw, (ftnlen)132);
    }
    return 0;
} /* s3prtn_ */

