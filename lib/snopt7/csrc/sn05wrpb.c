/* ./src/sn05wrpb.f -- translated by f2c (version 20100827).
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

static integer c__13 = 13;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__0 = 0;
static integer c_n4 = -4;
static doublereal c_b18 = -1.;
static integer c__5 = 5;
static integer c_n5 = -5;
static integer c__3 = 3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn05wrpb.f  --- user function wrapper for snOptB (aka snOpt) */

/*     s0fgB */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s0fgb_(integer *iexit, integer *modefg, logical *getcon, 
	logical *getobj, integer *n, integer *negcon, integer *nncon0, 
	integer *nncon, integer *nnjac, integer *nnl, integer *ngobj0, 
	integer *ngobj, S_fp fgcon, S_fp fgobj, doublereal *x, integer *ne, 
	integer *nlocj, integer *locj, integer *indj, doublereal *fcon, 
	doublereal *fobj, doublereal *gcon, doublereal *gobj, char *cu, 
	integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer *
	lenru, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1100[] = "(\002 The user has defined\002,i8,\002   out"
	    " of\002,i8,\002   constraint gradients.\002)";
    static char fmt_2010[] = "(\002 SnOptB  will define \002,i8,\002   gradi"
	    "ents for the \002,\002 FP objective.\002)";
    static char fmt_2000[] = "(\002 The user has defined\002,i8,\002   out"
	    " of\002,i8,\002   objective  gradients.\002)";
    static char fmt_2100[] = "(\002 XXX  Some objective  derivatives are mis"
	    "sing ---\002,\002 derivative level reduced to\002,i3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer l, lg, lx0, statususer;
    static char str[80];
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer modec, modef, ngrad;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dddiv_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer nnobj;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), s1time_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *), s8sclg_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *), s8sclj_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *), s8stat_(integer *, 
	    integer *, integer *);
    static logical scaled;
    static integer lascal, lgobju, nnglin, lxscal, lgconu, minmax, lvlscl;
    static doublereal gdummy;
    static logical fponly;
    static integer lvltim, status;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___21 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___23 = { 0, str, 0, fmt_2010, 80, 1 };
    static icilist io___24 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___25 = { 0, str, 0, fmt_2100, 80, 1 };


/*     ================================================================== */
/*     s0fgB  is an instance of fgwrap that calls the user-written */
/*     routines  fgcon  and  fgobj  to evaluate the problem functions */
/*     and possibly their gradients. */

/*     Arguments  fgcon  and  fgobj  are called using modefg to control */
/*     the gradients as follows: */

/*     modefg        Task */
/*     ------        ---- */
/*       2     Assign fCon, fObj and all known elements of gCon and gObj. */
/*       1     Assign all known elements of gCon and gObj. */
/*             (fObj and fCon are ignored). */
/*       0     Assign fObj, fCon.  (gCon and gObj are ignored). */

/*     If s0fgB is called with minmax = 0 (feasible point only) then */
/*     ngObj = max(nnJac,nnObj)  and the user objective is not used. */

/*     09-Jan 1992: First version based on Minos 5.4 routine m6fobj. */
/*     28-Sep 1993: Test for mode = -1 exit from funobj. */
/*     22-Dec 1994: Simplified to conform to npfobj in npopt. */
/*     16-Jan 1997: fgobj added as argument. */
/*     19-Mar 1997: First instantiated version for SNOPT 5.1. */
/*     17-Jul 1997: First thread-safe version. */
/*     02-Oct 1997: Character workspace added. */
/*     11 Oct 1998: Logicals getCon and getObj added. */
/*     30 Dec 2000: Housekeeping for Status = 1 included. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     15 Nov 2010: Call-status removed from the argument list. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* = 0,1,2 or 3, deriv. level */
/* > 0 => some differences needed */
/* > 0 => some exact derivs */
/* > 0 => constant Jacob elements */
/* calls to fCon: mode = 0 */
/* calls to fCon  mode > 0 */
/* calls to fObj: mode = 0 */
/* calls to fObj: mode > 0 */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --gcon;
    --fcon;
    --gobj;
    --indj;
    --locj;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    nnobj = iw[22];
/* # of objective variables */
    lvlscl = iw[75];
/* scale option */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    lvltim = iw[182];
/* Timing level */
    lascal = iw[296];
/* Ascale(nb) = row and column scales */
    lx0 = iw[298];
/* x0(nnL)    = proximal point starting point */
    lxscal = iw[302];
/* xScal(n)   = copy of scaled  x */
    lgconu = iw[319];
/* record of unknown derivatives and constants */
    lgobju = iw[323];
/* record of unknown derivatives */
    gdummy = rw[69];
/* definition of an 'unset' value */
    *iexit = 0;
    fponly = minmax == 0;
    scaled = lvlscl == 2;
    modec = *modefg;
    modef = *modefg;
/* Determine the user-function call-status. */
    s8stat_(&status, &iw[1], leniw);
    if (status == 1) {
/* -------------------------------------------------------------- */
/* First evaluation of the problem functions in snOptB */
/* On entry, lvlScl = 0. */
/* -------------------------------------------------------------- */
	iw[183] = 0;
	iw[184] = 0;
	iw[185] = 0;
	snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
	dload_(negcon, &gdummy, &gcon[1], &c__1);
	dload_(ngobj, &gdummy, &gobj[1], &c__1);
    }
/* ----------------------------------------------------------------- */
/* Unscale x (never required for Status = 1) */
/* ----------------------------------------------------------------- */
    if (scaled) {
	dcopy_(nnl, &x[1], &c__1, &rw[lxscal], &c__1);
	ddscl_(nnl, &rw[lascal], &c__1, &x[1], &c__1);
/* If the Jacobian has some constant elements, they are wrecked */
/* by the scaling.  Restore them from gConu. */
	if (*getcon) {
	    if (*modefg > 0 && iw[185] > 0) {
		dcopy_(negcon, &rw[lgconu], &c__1, &gcon[1], &c__1);
	    }
	}
    }
/* ================================================================= */
/* Compute the constraint functions. */
/* ================================================================= */
    if (*getcon) {
	statususer = status;
/* In case fgCon alters Status */
	if (lvltim >= 2) {
	    s1time_(&c__4, &c__0, &iw[1], leniw, &rw[1], lenrw);
	}
	(*fgcon)(&modec, nncon, nnjac, negcon, &x[1], &fcon[1], &gcon[1], &
		statususer, cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (
		ftnlen)8);
	if (lvltim >= 2) {
	    s1time_(&c_n4, &c__0, &iw[1], leniw, &rw[1], lenrw);
	}
	++iw[189];
	if (*modefg > 0) {
	    ++iw[190];
	}
    }
/* ================================================================= */
/* Compute the objective function. */
/* ================================================================= */
    if (*getobj && modec >= 0) {
	if (fponly) {
	    dcopy_(ngobj, &x[1], &c__1, &gobj[1], &c__1);
	    daxpy_(ngobj, &c_b18, &rw[lx0], &c__1, &gobj[1], &c__1);
	    *fobj = ddot_(ngobj, &gobj[1], &c__1, &gobj[1], &c__1) * .5;
	} else {
/* ngObj = nnObj */
	    statususer = status;
/* In case fgobj alters Status */
	    if (lvltim >= 2) {
		s1time_(&c__5, &c__0, &iw[1], leniw, &rw[1], lenrw);
	    }
	    (*fgobj)(&modef, ngobj, &x[1], fobj, &gobj[1], &statususer, cu + 
		    8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
	    if (lvltim >= 2) {
		s1time_(&c_n5, &c__0, &iw[1], leniw, &rw[1], lenrw);
	    }
	    ++iw[194];
	    if (*modefg > 0) {
		++iw[195];
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Scale  x and the derivatives. */
/* ----------------------------------------------------------------- */
    if (scaled) {
	dcopy_(nnl, &rw[lxscal], &c__1, &x[1], &c__1);
	if (*getcon) {
	    dddiv_(nncon, &rw[lascal + *n], &c__1, &fcon[1], &c__1);
	    if (*modefg > 0 && iw[184] > 0) {
		s8sclj_(nncon, nnjac, negcon, n, &rw[lascal], ne, nlocj, &
			locj[1], &indj[1], &gcon[1], &rw[1], lenrw);
	    }
	}
	if (*getobj && modec >= 0) {
	    if (*modefg > 0 && iw[184] > 0) {
		s8sclg_(ngobj, &rw[lascal], &gobj[1], &rw[1], lenrw);
	    }
	}
    }
    if (modec < 0 || modef < 0) {
/* -------------------------------------------------------------- */
/* The user may be saying the function is undefined (mode = -1) */
/* or may just want to stop                         (mode < -1). */
/* -------------------------------------------------------------- */
	if (modec == -1 || modef == -1) {
	    *iexit = -1;
	} else {
	    if (modec < 0) {
		*iexit = 72;
	    } else {
		*iexit = 73;
	    }
	}
    }
/* ================================================================= */
/* Do some housekeeping on the first snOptB entry. */
/* ================================================================= */
    if (status == 1 && *iexit == 0) {
	if (*getcon) {
/* ----------------------------------------------------------- */
/* Count how many Jacobian elements are provided. */
/* ----------------------------------------------------------- */
	    nnglin = 0;
	    ngrad = 0;
	    i__1 = *negcon;
	    for (l = 1; l <= i__1; ++l) {
		if (gcon[l] != gdummy) {
		    ++ngrad;
		}
	    }
	    s_wsfi(&io___21);
	    do_fio(&c__1, (char *)&ngrad, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*negcon), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	    if (ngrad < *negcon) {
/* Some Jacobian elements are missing. */
		if (iw[70] >= 2) {
/* ----------------------------------------------------- */
/* All the Jacobian is known.  Any undefined elements */
/* are assumed constant, and are restored from gConu. */
/* ----------------------------------------------------- */
		    snprnt_(&c__3, " ==>  Some constraint derivatives are mi"
			    "ssing,  assumed constant.", &iw[1], leniw, (
			    ftnlen)65);
		    snprnt_(&c__3, " ", &iw[1], leniw, (ftnlen)1);
		    lg = lgconu;
		    i__1 = *negcon;
		    for (l = 1; l <= i__1; ++l) {
			if (gcon[l] == gdummy) {
			    gcon[l] = rw[lg];
			    ++nnglin;
			}
			++lg;
		    }
		} else {
/* ----------------------------------------------------- */
/* Save a permanent copy of gCon in gConu so that we kn */
/* which derivatives must be estimated. */
/* ----------------------------------------------------- */
		    dcopy_(negcon, &gcon[1], &c__1, &rw[lgconu], &c__1);
		}
	    }
/* ngrad < negCon */
	    if (ngrad + nnglin < *negcon) {
		iw[183] = 1;
	    }
	    if (ngrad > 0) {
		iw[184] = 1;
	    }
	    if (nnglin > 0) {
		iw[185] = 1;
	    }
	}
	if (*getobj) {
/* ----------------------------------------------------------- */
/* Count how many working gradient elements are known. */
/* (These may be the gradients of the FP objective.) */
/* ----------------------------------------------------------- */
	    ngrad = 0;
	    i__1 = *ngobj;
	    for (l = 1; l <= i__1; ++l) {
		if (gobj[l] != gdummy) {
		    ++ngrad;
		}
	    }
	    if (fponly) {
		s_wsfi(&io___23);
		do_fio(&c__1, (char *)&(*ngobj), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	    } else {
		s_wsfi(&io___24);
		do_fio(&c__1, (char *)&ngrad, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nnobj, (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	    }
	    if (ngrad < *ngobj) {
/* Some objective gradients are missing. */
		if (iw[70] == 1 || iw[70] == 3) {
/* ----------------------------------------------------- */
/* The objective gradient was meant to be known. */
/* ----------------------------------------------------- */
		    --iw[70];
		    s_wsfi(&io___25);
		    do_fio(&c__1, (char *)&iw[70], (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
		}
/* -------------------------------------------------------- */
/* Copy gObj into gObju. */
/* -------------------------------------------------------- */
		dcopy_(ngobj, &gobj[1], &c__1, &rw[lgobju], &c__1);
	    }
	    if (ngrad < *ngobj) {
		iw[183] = 1;
	    }
	    if (ngrad > 0) {
		iw[184] = 1;
	    }
	}
    }
    return 0;
} /* s0fgb_ */

