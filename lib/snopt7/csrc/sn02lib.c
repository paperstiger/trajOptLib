/* ../snopt7/src/sn02lib.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__3 = 3;
static integer c__130 = 130;
static integer c__0 = 0;
static integer c__23 = 23;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__5 = 5;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn02lib.f */

/*     snTitl   snInit   snSpec   snchkA   snJac */
/*     snMem    snMemA   snMemA0  snMemB */
/*     snLog    snLog2   snSTOP */
/*     snSet    snSeti   snSetr */
/*     snGet    snGetc   snGeti   snGetr */
/*     snRetH */

/*     09 Mar 2004: snSolF implemented. */
/*     17 Jun 2004: snSolF always flags infeasible jbInf1 as I. */
/*     13 Dec 2013: Moved snEXIT, snWRAP, snSolF to sn04wrap.f. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int sntitl_(char *title, ftnlen title_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

/*     ================================================================== */
/*     snTitl sets the title for snopt. */
/*     ================================================================== */
    s_copy(title, "S N O P T  7.2-12.2 (Jul 2013)", (ftnlen)30, (ftnlen)30);
/* ---------------123456789|123456789|123456789|-------------------------- */
    return 0;
} /* sntitl_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snTitl */
/* Subroutine */ int sninit_(integer *iprint, integer *isumm, char *cw, 
	integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, ftnlen cw_len)
{
    /* Initialized data */

    static char dashes[30] = "==============================";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    char ch__1[39], ch__2[31];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int s3unsetall_(char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static char str[80], str2[80], title[30];
    static integer istdo;
    extern /* Subroutine */ int s1init_(char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    static integer ispecs, inform__;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), sntitl_(char *, 
	    ftnlen), snprnt_(integer *, char *, integer *, integer *, ftnlen);
    extern integer s1outpt_(void);

/*     ================================================================== */
/*     snInit  is called by the user to do the following: */
/*     1. Open default files (Print, Summary). */
/*     2. Initialize title. */
/*     3. Set options to default values. */

/*     15 Nov 1991: First version. */
/*     14 Jul 1997: Thread-safe version. */
/*     02 Oct 1997: Character workspace added. */
/*     15 Oct 2003: snEXIT and snPRNT added. */
/*     18 Jun 2007: First 500 elements of cw, iw and rw are initialized */
/*     18 Jun 2008: Global call-status values added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* start of SNOPT part of rw */
/* end   of SNOPT part of rw */
/* start of SNOPT part of iw */
/* end   of SNOPT part of iw */
/* start of SNOPT part of cw */
/* end   of SNOPT part of cw */
/* # nonlinear Jac, variables */
/* # variables in gObj */
/* # of nonlinear constraints */
/* nonlinear vars */
/* Timing level */
/* QP user-routine call-status */
/* NP user-routine call-status */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    s_copy(solver, "SNINIT", (ftnlen)6, (ftnlen)6);
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/* -------------------------------------------------------------- */
/* Not enough workspace to do ANYTHING! */
/* Print and exit without accessing the work arrays. */
/* -------------------------------------------------------------- */
	inform__ = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(&inform__, solver, str, str2, &iw[1], leniw, (ftnlen)6, (
		ftnlen)80, (ftnlen)80);
	goto L999;
    }
/* ----------------------------------------------------------------- */
/* Initialize cw, iw, rw so that they may be copied safely. */

/* This also sets the options to a specific "undefined" state. */
/* snopt  will check the options later and maybe print them. */
/* ----------------------------------------------------------------- */
    s3unsetall_(cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
/* ----------------------------------------------------------------- */
/* Initialize some default values. */
/* ----------------------------------------------------------------- */
    ispecs = 0;
    istdo = s1outpt_();
    iw[10] = istdo;
/* Standard Output */
    iw[11] = ispecs;
/* Specs file (default) */
    iw[12] = *iprint;
/* Print file */
    iw[13] = *isumm;
/* Summary file */
    iw[6] = 500;
    iw[4] = 500;
    iw[2] = 500;
    iw[7] = *lencw;
    iw[5] = *leniw;
    iw[3] = *lenrw;
/* ----------------------------------------------------------------- */
/* These dimensions need to be initialized for an MPS run. */
/* ----------------------------------------------------------------- */
    iw[23] = 0;
    iw[21] = 0;
    iw[22] = 0;
    iw[24] = 0;
    sntitl_(title, (ftnlen)30);
    s1init_(title, &iw[1], leniw, &rw[1], lenrw, (ftnlen)30);
/* Writing concatenation */
    i__1[0] = 9, a__1[0] = "         ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)39);
    snprnt_(&c__11, ch__1, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
    i__1[0] = 9, a__1[0] = "         ";
    i__1[1] = 30, a__1[1] = title;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)39);
    snprnt_(&c__1, ch__1, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
    i__1[0] = 9, a__1[0] = "         ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)39);
    snprnt_(&c__1, ch__1, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)31);
    snprnt_(&c__12, ch__2, &iw[1], leniw, (ftnlen)31);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = 30, a__1[1] = title;
    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)31);
    snprnt_(&c__2, ch__2, &iw[1], leniw, (ftnlen)31);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = 30, a__1[1] = dashes;
    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)31);
    snprnt_(&c__2, ch__2, &iw[1], leniw, (ftnlen)31);
/* ----------------------------------------------------------------- */
/* Initialize some global values. */
/* ----------------------------------------------------------------- */
    iw[235] = -11111;
    iw[236] = -11111;
    iw[182] = 3;
L999:
    return 0;
} /* sninit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snInit */
/* Subroutine */ int snspec_(integer *ispecs, integer *iexit, char *cw, 
	integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char str[80], str2[80];
    extern /* Subroutine */ int s3opt_();
    static integer calls, isumm;
    extern /* Subroutine */ int s3file_(integer *, integer *, integer *, U_fp,
	     char *, integer *, integer *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static integer iprint;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer errors;

/*     ================================================================== */
/*     snSpec  may be called by the user to read a Specs file. */

/*     07 Feb 1998: First version of snSpec. */
/*     01 Aug 2003: s3file now has a "title" parameter.  Use ' '. */
/*     13 Jul 2005: Included error count in value of iExit. */
/*     22 Apr 2007: Exit 107 added (some keywords unrecognized) */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNSPEC", (ftnlen)6, (ftnlen)6);
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
    if (*ispecs <= 0 || *ispecs > 99) {
	*iexit = 131;
/* iSPECS out of range */
	goto L800;
    }
    iw[11] = *ispecs;
/* Specs (options) file */
    iprint = iw[12];
/* Print file */
    isumm = iw[13];
/* Summary file */
    calls = 1;
/*     ------------------------------------------------------------------ */
/*     Read the Specs file. */
/*     snopt  will check the options later and maybe print them. */
/*     ------------------------------------------------------------------ */
    s3file_(iexit, &calls, ispecs, (U_fp)s3opt_, " ", &iprint, &isumm, &
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)1, (
	    ftnlen)8);
L800:
    if (*iexit == 0) {
	if (errors == 0) {
	    *iexit = 101;
/* SPECS file read successfully */
	} else {
	    *iexit = 107;
/* some SPECS keywords not recognized */
	}
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* snspec_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snSpec */
/* Subroutine */ int snchka_(integer *iexit, integer *nf, integer *n, integer 
	*lvlchk, U_fp userfg, integer *igfun, integer *jgvar, integer *leng, 
	integer *neg, doublereal *x, integer *mincw, integer *miniw, integer *
	minrw, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* Format strings */
    static char fmt_9010[] = "(\002 Total integer   workspace  should be sig"
	    "nificantly\002,\002 more than\002,i8)";
    static char fmt_9020[] = "(\002 Total real      workspace  should be sig"
	    "nificantly\002,\002 more than\002,i8)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer lf, lg, lw, ly, lz, lf1, lg1, lw1, lx1;
    static char str[80], str2[80];
    static integer lset, lelem;
    static doublereal egmax;
    extern /* Subroutine */ int s2mem0_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, ftnlen);
    static integer igmax, maxcw, maxiw, maxrw;
    extern /* Subroutine */ int s7chka_(integer *, integer *, U_fp, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, char *, integer *, integer *,
	     integer *, doublereal *, integer *, ftnlen);
    static integer nextcw;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nextiw;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer nextrw;

    /* Fortran I/O blocks */
    static icilist io___35 = { 0, str, 0, fmt_9010, 80, 1 };
    static icilist io___36 = { 0, str, 0, fmt_9020, 80, 1 };


/*     ================================================================== */
/*     snchkA  is a stand-alone derivative checker for problems in */
/*     snOptA format. */

/*     snchkA is different in two ways from the built-in derivative */
/*     checker in snopta: */
/*       1. The derivatives are checked before the variables and */
/*          constraints are reordered to conform with snoptb format. */

/*       2. The derivatives are checked at the point defined by the */
/*          input argument x.  A feasible point is NOT computed before */
/*          the check is done. */


/*     lvlChk has the following meaning: */

/*       -1         do not perform any check. */
/*        0         do the cheap test only. */
/*       >0         do both cheap and full test on problem derivatives. */

/*     ------------------------------------------------------------------ */
/*     NOTE: Before calling snchkA, there MUST be a call to the */
/*     initialization routine, i.e., */
/*     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw ) */
/*     This sets the default values of the optional parameters. You can */
/*     also alter the default values of iPrint and iSumm before snchkA */
/*     is used.  iPrint = 0, etc, is OK. */
/*     ------------------------------------------------------------------ */

/*     EXIT snchkA 100 -- completed successfully */
/*     EXIT INFO   105 -- user-supplied derivatives appear to be correct */
/*     EXIT INFO   106 -- no derivatives were checked */

/*     EXIT snchkA  50 -- error in the user-supplied functions */
/*     EXIT INFO    55 -- incorrect derivatives */

/*     EXIT snchkA  70 -- user requested termination */
/*     EXIT INFO    71 -- terminated during function evaluation */

/*     EXIT snchkA  80 -- insufficient storage allocated */
/*     EXIT INFO    81 -- work arrays must have at least 500 elements */
/*     EXIT INFO    83 -- not enough integer storage */

/*     SNOPT package maintained by Philip E. Gill, */
/*     Dept of Mathematics, University of California, San Diego. */

/*     01 Jun 2006: First version of snchkA */
/*     22 Apr 2007: INFO 106 added */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --jgvar;
    --igfun;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNCHKA", (ftnlen)6, (ftnlen)6);
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Check memory limits and fetch the workspace starting positions. */
/*     ------------------------------------------------------------------ */
    s2mem0_(iexit, solver, lencw, leniw, lenrw, &iw[1], mincw, miniw, minrw, &
	    maxcw, &maxiw, &maxrw, &nextcw, &nextiw, &nextrw, (ftnlen)6);
    if (*iexit != 0) {
	goto L999;
    }
    lelem = nextiw;
    lset = lelem + *nf;
    *miniw = lset + *n - 1;
    lf = nextrw;
    lg = lf + *nf;
    lx1 = lg + *leng;
    lf1 = lx1 + *n;
    lg1 = lf1 + *nf;
    lw = lg1 + *leng;
    lw1 = lw + *nf;
    ly = lw1 + *nf;
    lz = ly + *n;
    *minrw = lz + *n - 1;
    if (*miniw > maxiw || *minrw > maxrw) {
/*        --------------------------------------------------------------- */
/*        Not enough space to check the derivatives. */
/*        Exit with an (over) estimate of the additional space needed. */
/*        --------------------------------------------------------------- */
	if (*miniw > maxiw) {
	    s_wsfi(&io___35);
	    do_fio(&c__1, (char *)&(*miniw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    *iexit = 83;
	}
	if (*minrw > maxrw) {
	    s_wsfi(&io___36);
	    do_fio(&c__1, (char *)&(*minrw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    *iexit = 84;
	}
	goto L800;
    }
/*     ------------------------------------------------------------------ */
/*     Go for it. */
/*     ------------------------------------------------------------------ */
    s7chka_(iexit, lvlchk, (U_fp)userfg, nf, n, &igmax, &egmax, &iw[lelem], &
	    iw[lset], &igfun[1], &jgvar[1], leng, neg, &x[1], &rw[lf], &rw[lg]
	    , &rw[lx1], &rw[lf1], &rw[lg1], &rw[lw], &rw[lw1], &rw[ly], &iw[1]
	    , leniw, &rw[1], lenrw, cu + 8, lencu, &iu[1], leniu, &ru[1], 
	    lenru, (ftnlen)8);
    if (*iexit != 0) {
	goto L800;
    }
/*     Print the exit conditions. */
L800:
    if (*iexit == 0) {
	*iexit = 105;
/* all derivatives appear to be correct */
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
/* L2000: */
} /* snchka_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snchkA */
/* Subroutine */ int snjac_(integer *iexit, integer *nf, integer *n, U_fp 
	userfg, integer *iafun, integer *javar, integer *lena, integer *nea, 
	doublereal *a, integer *igfun, integer *jgvar, integer *leng, integer 
	*neg, doublereal *x, doublereal *xlow, doublereal *xupp, integer *
	mincw, integer *miniw, integer *minrw, char *cu, integer *lencu, 
	integer *iu, integer *leniu, doublereal *ru, integer *lenru, char *cw,
	 integer *lencw, integer *iw, integer *leniw, doublereal *rw, integer 
	*lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_9010[] = "(\002 Total integer   workspace  should be sig"
	    "nificantly\002,\002 more than\002,i8)";
    static char fmt_9020[] = "(\002 Total real      workspace  should be sig"
	    "nificantly\002,\002 more than\002,i8)";
    static char fmt_1000[] = "(\002 Nonzero derivs  Jij  \002,i8)";
    static char fmt_1010[] = "(\002 Non-constant    Jij's\002,i8,5x,\002Cons"
	    "tant Jij's\002,4x,i8)";
    static char fmt_2000[] = "(\002 -->  largest error in the estimated Jaco"
	    "bian is\002,1p,e12.2,\002  in row\002,i6)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer lf, lg, lw, lx, ly, lz, lfw, lfy, lfz;
    static char str[80], str2[80];
    static integer isum, iprt;
    extern /* Subroutine */ int s7jac_(integer *, U_fp, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static doublereal emaxj;
    extern /* Subroutine */ int s2mem0_(integer *, char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, ftnlen);
    static integer imaxj, maxcw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxiw, maxrw, lxlow, lxupp;
    static doublereal infbnd;
    static integer lgcolw, lgcoly, lgcolz, deropt;
    static logical prtall;
    extern /* Subroutine */ int snseti_(char *, integer *, integer *, integer 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);
    static integer nextcw;
    static char solver[6];
    static integer coltyp;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nextiw;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer nextrw, rowtyp;

    /* Fortran I/O blocks */
    static icilist io___64 = { 0, str, 0, fmt_9010, 80, 1 };
    static icilist io___65 = { 0, str, 0, fmt_9020, 80, 1 };
    static icilist io___73 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___74 = { 0, str, 0, fmt_1010, 80, 1 };
    static icilist io___75 = { 0, str, 0, fmt_2000, 80, 1 };


/*     ================================================================== */
/*     snJac  computes the coordinates of the Jacobian. */

/*     All calculations are based on a point defined by moving the input */
/*     x  inside its upper and lower bounds.  snJac is terminated if */
/*     the problem functions are undefined at this point. */

/*     ------------------------------------------------------------------ */
/*     NOTE: Before calling snJac, your calling program MUST call the */
/*     initialization routine using the call: */
/*     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw ) */
/*     This sets the default values of the optional parameters. You can */
/*     also alter the default values of iPrint and iSumm before snJac */
/*     is used.  iPrint = 0, etc, is OK. */
/*     ------------------------------------------------------------------ */

/*     EXIT snJac  100 -- completed successfully (for auxiliary routines) */
/*     EXIT INFO   102 -- Jacobian structure estimated */

/*     EXIT snJac  12x -- Errors while estimating Jacobian structure */
/*     EXIT INFO   121 -- cannot estimate Jacobian structure at given point */
/*     EXIT INFO    71 -- terminated during function evaluation */
/*     EXIT INFO    81 -- work arrays must have at least 500 elements */
/*     EXIT INFO    83 -- not enough integer storage */

/*     SNOPT package maintained by Philip E. Gill, */
/*     Dept of Mathematics, University of California, San Diego. */

/*     26 Oct 2002: First version of snJac. */
/*     27 Sep 2003: More thorough checks for feasibility */
/*     15 Jun 2008: Call-status implemented correctly. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xupp;
    --xlow;
    --x;
    --a;
    --javar;
    --iafun;
    --jgvar;
    --igfun;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNJAC ", (ftnlen)6, (ftnlen)6);
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Check memory limits and fetch the workspace starting positions. */
/*     ------------------------------------------------------------------ */
    s2mem0_(iexit, solver, lencw, leniw, lenrw, &iw[1], mincw, miniw, minrw, &
	    maxcw, &maxiw, &maxrw, &nextcw, &nextiw, &nextrw, (ftnlen)6);
    if (*iexit != 0) {
	goto L999;
    }
    rowtyp = nextiw;
    coltyp = rowtyp + *nf;
    *miniw = coltyp + *n - 1;
    lf = nextrw;
    lg = lf + *nf;
    lx = lg + *leng;
    lxlow = lx + *n;
    lxupp = lxlow + *n;
    lw = lxupp + *n;
    ly = lw + *n;
    lz = ly + *n;
    lfw = lz + *n;
    lfy = lfw + *nf;
    lfz = lfy + *nf;
    lgcolw = lfz + *nf;
    lgcoly = lgcolw + *nf;
    lgcolz = lgcoly + *nf;
    *minrw = lgcolz + *nf - 1;
    if (*miniw > maxiw || *minrw > maxrw) {
/*        --------------------------------------------------------------- */
/*        Not enough space to build the Jacobian. */
/*        Provide the user an (over) estimate of what is needed. */
/*        --------------------------------------------------------------- */
	if (*miniw > maxiw) {
	    s_wsfi(&io___64);
	    do_fio(&c__1, (char *)&(*miniw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    *iexit = 83;
	}
	if (*minrw > maxrw) {
	    s_wsfi(&io___65);
	    do_fio(&c__1, (char *)&(*minrw), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	    *iexit = 84;
	}
	goto L800;
    }
/*     ------------------------------------------------------------------ */
/*     Go for it. */
/*     ------------------------------------------------------------------ */
    prtall = TRUE_;
/* ignore fixed variables for diffs */
    infbnd = rw[70];
/* definition of an infinite bound */
    if (infbnd < 0.) {
	infbnd = 1e20;
/* User hasn't assigned it */
    }
/*     ------------------------------------------------------------------ */
/*     Make sure that snOptA does not let userf compute derivatives. */
/*     The parameters iPrt and iSum may refer to the Print and Summary */
/*     file respectively.  Setting them to 0 suppresses printing. */
/*     ------------------------------------------------------------------ */
    deropt = 0;
    iprt = 0;
    isum = 0;
    snseti_("Derivative option", &deropt, &iprt, &isum, iexit, cw + 8, lencw, 
	    &iw[1], leniw, &rw[1], lenrw, (ftnlen)17, (ftnlen)8);
    dcopy_(n, &x[1], &c__1, &rw[lx], &c__1);
    dcopy_(n, &xlow[1], &c__1, &rw[lxlow], &c__1);
    dcopy_(n, &xupp[1], &c__1, &rw[lxupp], &c__1);
    s7jac_(iexit, (U_fp)userfg, &prtall, nf, n, &infbnd, &imaxj, &emaxj, &
	    iafun[1], &javar[1], lena, nea, &a[1], &igfun[1], &jgvar[1], leng,
	     neg, &iw[rowtyp], &iw[coltyp], &rw[lx], &rw[lxlow], &rw[lxupp], &
	    rw[lf], &rw[lg], &rw[lw], &rw[ly], &rw[lz], &rw[lfw], &rw[lfy], &
	    rw[lfz], &rw[lgcolw], &rw[lgcoly], &rw[lgcolz], &iw[1], leniw, cu 
	    + 8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
    if (*iexit != 0) {
	goto L800;
    }
    s_wsfi(&io___73);
    i__1 = *nea + *neg;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    s_wsfi(&io___74);
    do_fio(&c__1, (char *)&(*neg), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nea), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
    if (emaxj > .1) {
	*iexit = 121;
/* unable to estimate Jacobian structure */
	s_wsfi(&io___75);
	do_fio(&c__1, (char *)&emaxj, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&imaxj, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
    }
/*     Print the exit conditions. */
L800:
    if (*iexit == 0) {
	*iexit = 102;
/* Jacobian structure estimated */
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* snjac_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snJac */
/* Subroutine */ int snmem_(integer *iexit, integer *m, integer *n, integer *
	ne, integer *negcon, integer *nncon, integer *nnjac, integer *nnobj, 
	integer *mincw, integer *miniw, integer *minrw, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cw_len)
{
    extern /* Subroutine */ int snmemb_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);

/*     ================================================================== */
/*     snMem   estimates the memory requirements for snOptB, */
/*     using the values: */
/*        m    , n    , ne    negCon, */
/*        nnObj, nnCon, nnJac */

/*     These values are used to compute the minimum required storage: */
/*     mincw, miniw, minrw. */

/*     Note: */
/*     1. The initialization routine snInit MUST be called before snMem. */

/*     2. All default parameters must be set before calling snMem, */
/*        since some values affect the amount of memory required. */

/*     3. The arrays rw and iw hold  constants and work-space addresses. */
/*        They must have dimension at least 500. */

/*     4. This version of snMem does not allow user accessible */
/*        partitions of cw, iw and rw. */

/*     Exit messages: */

/*     SNMEMB EXIT  80 -- insufficient storage allocated */
/*     SNMEMB INFO  81 -- work arrays must have at least 500 elements */

/*     SNMEMB EXIT 100 -- finished successfully */
/*     SNMEMB INFO 104 -- requirements estimated */

/*     29 Mar 1998: First version. */
/*     15 Oct 2003: iExit added as an argument. */
/*     15 Oct 2003: Current version of snMem. */
/*     ================================================================== */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    snmemb_(iexit, m, n, ne, negcon, nncon, nnjac, nnobj, mincw, miniw, minrw,
	     cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
    return 0;
} /* snmem_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snMem */
/* Subroutine */ int snmema_(integer *iexit, integer *nf, integer *n, integer 
	*nxname, integer *nfname, integer *nea, integer *neg, integer *mincw, 
	integer *miniw, integer *minrw, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    static integer isumm, iprint;
    extern /* Subroutine */ int snmema0_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);

/*     ================================================================== */
/*     snMemA   estimates the memory requirements for snOptA, */
/*     using the values: */

/*        nF, n, nxname, nFname, neA, neG */

/*     These values are used to compute the minimum required storage: */
/*     mincw, miniw, minrw. */

/*     Note: */
/*     1. The initialization routine snInit MUST be called before snMemA. */

/*     2. Some optional parameter settings affect the amount of memory */
/*        needed, so weird things may happen if some optional parameters */
/*        are set after the call to snMemA. */

/*     3. The arrays rw and iw hold constants and work-space addresses. */
/*        They must have dimension at least 500. */

/*     4. This version of snMemA does not allow user accessible */
/*        partitions of cw, iw and rw. */

/*     Exit messages: */

/*     SNMEMA EXIT  80 -- insufficient storage allocated */
/*     SNMEMA INFO  81 -- work arrays must have at least 500 elements */

/*     SNMEMA EXIT 100 -- finished successfully */
/*     SNMEMA INFO 104 -- requirements estimated */

/*     01 Aug 2002: First version based on snMem. */
/*     31 Jul 2003: snEXIT and snPRNT adopted. */
/*     15 Oct 2003: iExit added as an argument. */
/*     09 Nov 2004: Optional printing added. */
/*     09 Nov 2004: Current version of snMemA. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    iprint = iw[12];
    isumm = iw[13];
    snmema0_(iexit, &iprint, &isumm, nf, n, nxname, nfname, nea, neg, mincw, 
	    miniw, minrw, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)8);
    return 0;
} /* snmema_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snMemA */
/* Subroutine */ int snmema0_(integer *iexit, integer *lprint, integer *lsumm,
	 integer *nf, integer *n, integer *nxname, integer *nfname, integer *
	nea, integer *neg, integer *mincw, integer *miniw, integer *minrw, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer m, ne, nkx;
    static char str[80], str2[80];
    static integer lenr, maxr, maxs;
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s8map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer nnjac, nname, nnobj, nncon, maxcw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxiw, isumm, maxrw;
    extern /* Subroutine */ int s3mapa_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), s2bmap_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     s8dflt_(integer *, integer *, integer *, integer *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer negcon, llencw;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer inform__, lleniw, mqnmod, lvlhes, llenrw;
    static logical prtmem;
    static integer iprint, liwest;
    static char usercw[8*130], solver[6];
    static integer nextcw;
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nextiw, useriw[130], lrwest, nextrw;
    static doublereal userrw[130];

/*     ================================================================== */
/*     snMemA0   does the work for snMemA. */

/*     09 Nov 2004: First version */
/*     09 Nov 2004: Current version of  snMemA0. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    iprint = iw[12];
    isumm = iw[13];
    iw[12] = *lprint;
/* Print (Log) file */
    iw[13] = *lsumm;
/* Specs (options) file */
    s_copy(solver, "SNMEMA", (ftnlen)6, (ftnlen)6);
    *iexit = 0;
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
/*     Save the user's option choices  (weird choices get overwritten). */
    chcopy_(&c__130, cw + 408, &c__1, usercw, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__130, &rw[51], &c__1, userrw, &c__1);
/*     Assign fake values for lencw, leniw, lenrw. */
/*     This will force s2Mem to estimate the memory requirements. */
    llenrw = 500;
    lleniw = 500;
    llencw = 500;
/*     An obligatory call to snInit has `undefined' all options. */
/*     Check the user-defined values and assign undefined values. */
/*     s8dflt needs various problem dimensions in iw. */
/*     Allocate temporary work arrays for s3size. */
    nkx = *n + *nf;
/*     Provide the user an (over) estimate of what is needed. */
    ne = *nea + *neg;
    m = *nf;
    if (*nxname == 1 && *nfname == 1) {
	nname = 1;
    } else {
	nname = *n + m;
    }
    nncon = m;
    nnjac = *n;
    nnobj = *n;
    negcon = ne;
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = m;
/* copy of the number of rows */
    iw[17] = ne;
/* copy of the number of nonzeros in Jcol */
    iw[21] = nnjac;
/* # nonlinear Jacobian variables */
    iw[22] = nnobj;
/* # variables in gObj */
    iw[23] = nncon;
/* # of nonlinear constraints */
    s8dflt_(&m, n, &nncon, &nnjac, &nnobj, cw + 8, &llencw, &iw[1], &lleniw, &
	    rw[1], &llenrw, (ftnlen)8);
    nextcw = 501;
    nextiw = 501;
    nextrw = 501;
    maxcw = *lencw;
    maxiw = *leniw;
    maxrw = *lenrw;
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    nkx = *n + m;
    s8map_(&m, n, &negcon, &nkx, &nncon, &nnjac, &nnobj, &lenr, &maxr, &maxs, 
	    &mqnmod, &lvlhes, &nextcw, &nextiw, &nextrw, &iw[1], leniw);
    s3mapa_(&m, n, &ne, nf, neg, &negcon, &nkx, &nnjac, &nname, &nextcw, &
	    nextiw, &nextrw, &iw[1], leniw);
    s2bmap_(&m, n, &ne, &maxs, &nextiw, &nextrw, &maxiw, &maxrw, &liwest, &
	    lrwest, &iw[1], leniw);
    prtmem = FALSE_;
/* Suppress messages from s2Mem */
    s2mem_(&inform__, &prtmem, &liwest, &lrwest, &nextcw, &nextiw, &nextrw, &
	    maxcw, &maxiw, &maxrw, &llencw, &lleniw, &llenrw, mincw, miniw, 
	    minrw, &iw[1]);
/*     mincw = mincw */
    *miniw = liwest;
    *minrw = lrwest;
/*     Restore the user's choices of options. */
    chcopy_(&c__130, usercw, &c__1, cw + 408, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, useriw, &c__1, &iw[51], &c__1);
    dcopy_(&c__130, userrw, &c__1, &rw[51], &c__1);
/*     Print the exit conditions. */
    if (*iexit == 0) {
	*iexit = 104;
/* memory requirements estimated */
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    iw[12] = iprint;
    iw[13] = isumm;
    return 0;
} /* snmema0_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snMemA0 */
/* Subroutine */ int snmemb_(integer *iexit, integer *m, integer *n, integer *
	ne, integer *negcon, integer *nncon, integer *nnjac, integer *nnobj, 
	integer *mincw, integer *miniw, integer *minrw, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer nkx;
    static char str[80], str2[80];
    static integer lenr, maxr, maxs;
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s8map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxcw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxiw, maxrw;
    extern /* Subroutine */ int s2bmap_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), s8dflt_(integer *, integer *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static integer llencw;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer inform__, lleniw, mqnmod, lvlhes, llenrw;
    static logical prtmem;
    static integer liwest;
    static char usercw[8*130];
    static integer nextcw;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nextiw, useriw[130], lrwest;
    static doublereal userrw[130];
    static integer nextrw;

/*     ================================================================== */
/*     snMemB   estimates the memory requirements for snoptB, */
/*     using the values: */
/*        m    , n    , ne    negCon, */
/*        nnObj, nnCon, nnJac */

/*     These values are used to compute the minimum required storage: */
/*     mincw, miniw, minrw. */

/*     Note: */
/*     1. The initialization routine snInit MUST be called before snMemB. */

/*     2. All default parameters must be set before calling snMemB, */
/*        since some values affect the amount of memory required. */

/*     3. The arrays rw and iw hold  constants and work-space addresses. */
/*        They must have dimension at least 500. */

/*     4. This version of snMemB does not allow user-accessible */
/*        partitions of cw, iw and rw. */

/*     Exit messages: */

/*     SNMEMB EXIT  80 -- insufficient storage allocated */
/*     SNMEMB INFO  81 -- work arrays must have at least 500 elements */

/*     SNMEMB EXIT 100 -- finished successfully */
/*     SNMEMB INFO 104 -- requirements estimated */

/*     29 Mar 1998: First version. */
/*     31 Jul 2003: snEXIT and snPRNT adopted. */
/*     15 Oct 2003: iExit added as an argument. */
/*     19 Feb 2004: Current version of snMemB. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNMEMB", (ftnlen)6, (ftnlen)6);
    *iexit = 0;
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
/*     Save the user's option choices  (weird choices get overwritten). */
    chcopy_(&c__130, cw + 408, &c__1, usercw, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__130, &rw[51], &c__1, userrw, &c__1);
/*     Assign fake values for lencw, leniw, lenrw. */
/*     This will force s2Mem to estimate the memory requirements. */
    llenrw = 500;
    lleniw = 500;
    llencw = 500;
/*     An obligatory call to snInit has `undefined' all options. */
/*     Check the user-defined values and assign undefined values. */
/*     s8dflt needs various problem dimensions in iw. */
    iw[15] = *n;
/* copy of the number of columns */
    iw[16] = *m;
/* copy of the number of rows */
    iw[17] = *ne;
/* copy of the number of nonzeros in Jcol */
    iw[21] = *nnjac;
/* # nonlinear Jacobian variables */
    iw[22] = *nnobj;
/* # variables in gObj */
    iw[23] = *nncon;
/* # of nonlinear constraints */
    s8dflt_(m, n, nncon, nnjac, nnobj, cw + 8, &llencw, &iw[1], &lleniw, &rw[
	    1], &llenrw, (ftnlen)8);
    nextcw = 501;
    nextiw = 501;
    nextrw = 501;
    maxcw = *lencw;
    maxiw = *leniw;
    maxrw = *lenrw;
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    nkx = *n + *m;
    s8map_(m, n, negcon, &nkx, nncon, nnjac, nnobj, &lenr, &maxr, &maxs, &
	    mqnmod, &lvlhes, &nextcw, &nextiw, &nextrw, &iw[1], leniw);
    s2bmap_(m, n, ne, &maxs, &nextiw, &nextrw, &maxiw, &maxrw, &liwest, &
	    lrwest, &iw[1], leniw);
    prtmem = FALSE_;
/* Suppress messages from s2Mem */
    s2mem_(&inform__, &prtmem, &liwest, &lrwest, &nextcw, &nextiw, &nextrw, &
	    maxcw, &maxiw, &maxrw, &llencw, &lleniw, &llenrw, mincw, miniw, 
	    minrw, &iw[1]);
/*     mincw = mincw */
    *miniw = liwest;
    *minrw = lrwest;
/*     Restore the user's choices of options. */
    chcopy_(&c__130, usercw, &c__1, cw + 408, &c__1, (ftnlen)8, (ftnlen)8);
    icopy_(&c__130, useriw, &c__1, &iw[51], &c__1);
    dcopy_(&c__130, userrw, &c__1, &rw[51], &c__1);
/*     Print the exit conditions. */
    if (*iexit == 0) {
	*iexit = 104;
/* memory requirements estimated */
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* snmemb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snMemB */
/* Subroutine */ int snlog_(integer *iabort, integer *info, integer *hqntype, 
	logical *ktcond, integer *mjrprt, integer *minimz, integer *n, 
	integer *nb, integer *nncon0, integer *ns, integer *itn, integer *
	nmajor, integer *nminor, integer *nswap, doublereal *condhz, integer *
	iobj, doublereal *sclobj, doublereal *objadd, doublereal *fmrt, 
	doublereal *pennrm, doublereal *step, doublereal *prinf, doublereal *
	duinf, doublereal *vimax, doublereal *virel, integer *hs, integer *ne,
	 integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	doublereal *ascale, doublereal *bl, doublereal *bu, doublereal *fcon, 
	doublereal *ycon, doublereal *x, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* Initialized data */

    static char flag__[1*2] = " " " ";
    static char key[2*6] = "fr" "lo" "up" "sb" "bs" "fx";

    /* Format strings */
    static char fmt_3000[] = "(i7,i6,i7,1p,e8.1,i7,1x,2(a,e7.1,a),e14.7,i8,i"
	    "6,i7,e8.1,e8.1,a1,a)";
    static char fmt_3100[] = "(i7,i6,i7,1p,e8.1,i7,10x,1(a,e7.1,a),e14.7,i8,"
	    "i6,i7,e8.1,a1,a)";
    static char fmt_3200[] = "(i7,i6,i7,1p,e8.1,i7,10x,1(a,e7.1,a),e14.7,i8,"
	    "i6,i7,a1,a)";
    static char fmt_5000[] = "(i6,i7,1p,e9.1,i7,1x,2(a,e7.1,a),e14.7,i6,e8.1"
	    ",a)";
    static char fmt_5100[] = "(i6,i7,1p,e9.1,i7,10x,1(a,e7.1,a),e14.7,i6,a)";
    static char fmt_7600[] = "(\002 Maximum constraint violation    =\002,1p"
	    ",e12.4,4x,\002 ( =\002,e11.4,\002 normalized)\002)";
    static char fmt_7410[] = "(\002 x(\002,i6,\002)\002,1p,e13.5,1x,a2)";
    static char fmt_7420[] = "(\002   \002,i6,\002 \002,1p,e13.5)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3, i__4;
    char ch__1[106], ch__2[98], ch__3[90];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k, l, k1, k2, ir, lu;
    static char str[80];
    static integer lenl;
    static char form[40];
    static integer lenu;
    static logical prtf;
    static integer itns, mjrs;
    static logical prtj;
    static integer mnrs;
    static logical prtl, prtx, summ1, prnt1;
    static char cflag[1];
    static logical phead;
    static integer nnjac;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dddiv_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer nfobj[4];
    static char buffp[115];
    static integer nfcon[4];
    static char buffs[84];
    static integer nline, nncon, nnobj;
    static doublereal merit;
    static logical prntc;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    s2scla_(integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static logical major0;
    extern integer s2coln_(integer *, integer *, integer *), s2rown_(integer *
	    , integer *, integer *);
    static logical scaled;
    static doublereal infbnd;
    static char ktflag[1*4];
    static logical nlnobj, nlncon, prthdg;
    static doublereal sgnobj;
    static char mjrmsg[8];
    static integer lvlscl;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___171 = { 0, buffp, 0, fmt_3000, 115, 1 };
    static icilist io___172 = { 0, buffp, 0, fmt_3100, 115, 1 };
    static icilist io___173 = { 0, buffp, 0, fmt_3200, 115, 1 };
    static icilist io___174 = { 0, buffs, 0, fmt_5000, 84, 1 };
    static icilist io___175 = { 0, buffs, 0, fmt_5100, 84, 1 };
    static icilist io___176 = { 0, buffs, 0, fmt_5100, 84, 1 };
    static icilist io___186 = { 0, str, 0, form, 80, 1 };
    static icilist io___187 = { 0, str, 0, form, 80, 1 };
    static icilist io___188 = { 0, str, 0, form, 80, 1 };
    static icilist io___189 = { 0, str, 0, fmt_7600, 80, 1 };
    static icilist io___190 = { 0, str, 0, fmt_7410, 80, 1 };
    static icilist io___194 = { 0, str, 0, fmt_7420, 80, 1 };


/*     ================================================================== */
/*     snLog  prints the major iteration log. */

/*     The end-of-line summary is as follows: */

/*     Position    | Possible Entries */
/*     ------------|------------------- */
/*     1 = Blank   | */
/*     2 = Update  | n     s   - */
/*     3 = Modify  |    M  m   - */
/*     4 = HQNType |    R  r   - */
/*     5 = Steps   |    d  l   - */
/*     6 = QPfea   |    i  -   - */
/*     7 = QPerr   |    t  T   u   w  z */
/*     8 = FDiff   |    c  -   -   -  - */

/*     15 Nov 1991: First version. */
/*     19 Jul 1997: Thread-safe version. */
/*     02 Dec 2000: Reordered and sparsified. */
/*     28 Dec 2000: Row and column permutations added. */
/*     31 Jul 2003: snPRNT adopted. */
/*     31 Jul 2003: Current version of snLog. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* >0 => Major heading for iPrint */
/*     ------------------------------------------------------------------ */
/* >0 => Major heading for iSumm */
    /* Parameter adjustments */
    --info;
    --ktcond;
    --x;
    --bu;
    --bl;
    --ascale;
    --hs;
    --ycon;
    --fcon;
    --jcol;
    --indj;
    --locj;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    nncon = iw[23];
/* # of nonlinear constraints */
    nnjac = iw[21];
/* # nonlinear Jacobian variables */
    nnobj = iw[22];
/* # variables in gObj */
    lvlscl = iw[75];
/* scale option */
    lenl = iw[173];
/* size of current  L */
    lenu = iw[174];
/* size of current  U */
    nfcon[0] = iw[189];
/* number of calls of fCon */
    nfcon[1] = iw[190];
/* number of calls of fCon */
    nfcon[2] = iw[191];
/* number of calls of fCon */
    nfcon[3] = iw[192];
/* number of calls of fCon */
    nfobj[0] = iw[194];
/* number of calls of fObj */
    nfobj[1] = iw[195];
/* number of calls of fObj */
    nfobj[2] = iw[196];
/* number of calls of fObj */
    nfobj[3] = iw[197];
/* number of calls of fObj */
    infbnd = rw[70];
/* definition of an infinite bound */
    *iabort = 0;
    lu = lenl + lenu;
    nlncon = nncon > 0;
    nlnobj = nnobj > 0;
    prnt1 = *mjrprt >= 1;
    prntc = *mjrprt >= 100;
    summ1 = *mjrprt >= 1;
    major0 = *nmajor == 0;
    if (major0 && prnt1) {
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    }
    itns = *itn % 1000000;
    mnrs = *nminor % 1000000;
    mjrs = *nmajor % 1000000;
    s_copy(mjrmsg, "_", (ftnlen)8, (ftnlen)1);
/* Used to detect major iteration summary */
    s_copy(buffp, " ", (ftnlen)115, (ftnlen)1);
    s_copy(buffs, " ", (ftnlen)84, (ftnlen)1);
/* Add the alphabet soup. */
    if (info[1] == 0) {
	*(unsigned char *)&mjrmsg[1] = 'n';
/* No update could be made */
    } else if (info[1] == 2) {
	*(unsigned char *)&mjrmsg[1] = 's';
/* Scaled BFGS */
    }
    if (info[2] == 1) {
	*(unsigned char *)&mjrmsg[2] = 'M';
/* BFGS + qN Hessian mod. 1 */
    } else if (info[2] == 2) {
	*(unsigned char *)&mjrmsg[2] = 'm';
/* BFGS + qN Hessian mods. 1 + 2 */
    }
    if (*hqntype == 1) {
	*(unsigned char *)&mjrmsg[3] = 'R';
/* H set to a diagonal */
    } else if (*hqntype == 2) {
	*(unsigned char *)&mjrmsg[3] = 'r';
/* H set to the identity */
    }
    if (info[3] == 1) {
	*(unsigned char *)&mjrmsg[4] = 'd';
/* Violation limited via step */
    } else if (info[3] == 2) {
	*(unsigned char *)&mjrmsg[4] = 'D';
/* User-defined limit */
    } else if (info[3] == 3) {
	*(unsigned char *)&mjrmsg[4] = 'l';
/* Vars limited via step */
    }
    if (info[4] == 1) {
	*(unsigned char *)&mjrmsg[5] = 'i';
/* QP infeasible */
    }
    if (info[5] == 1) {
	*(unsigned char *)&mjrmsg[6] = 'T';
/* terminated by WSlimit */
    } else if (info[5] == 2) {
	*(unsigned char *)&mjrmsg[6] = 't';
/* terminated by itMax */
    } else if (info[5] == 3) {
	*(unsigned char *)&mjrmsg[6] = 'u';
/* QP unbounded */
    } else if (info[5] == 4) {
	*(unsigned char *)&mjrmsg[6] = 'w';
/* Weak QP solutions */
    } else if (info[5] == 5) {
	*(unsigned char *)&mjrmsg[6] = 'z';
/* superbasic limit reached */
    }
    if (info[6] == 1) {
	*(unsigned char *)&mjrmsg[7] = 'c';
/* Central differences */
    }
/* Put ( ) around small primal and dual infeasibilities. */
    for (k = 1; k <= 4; ++k) {
	*(unsigned char *)&ktflag[k - 1] = ' ';
    }
    k = 1;
    for (j = 1; j <= 2; ++j) {
	if (ktcond[j]) {
	    *(unsigned char *)&ktflag[k - 1] = '(';
	    *(unsigned char *)&ktflag[k] = ')';
	}
	k += 2;
    }
    sgnobj = (doublereal) (*minimz);
    merit = sgnobj * (*objadd + *fmrt);
    if (prnt1) {
/*        ------------------------------------------ */
/*        Terse line for the Print file. */
/*        ------------------------------------------ */
	prthdg = iw[224] > 0;
	nline = mjrs % 20;
	phead = nline == 0 || prthdg;
	if (phead) {
	    *(unsigned char *)cflag = *(unsigned char *)&flag__[min(nline,1)];
	    iw[224] = 0;
	}
	if (nlncon) {
	    if (phead) {
/* Writing concatenation */
		i__1[0] = 105, a__1[0] = "   Itns Major Minors    Step   nCo"
			"n Feasible  Optimal  MeritFunction     L+U BSwap    "
			" nS  condHz Penalty";
		i__1[1] = 1, a__1[1] = cflag;
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)106);
		snprnt_(&c__11, ch__1, &iw[1], leniw, (ftnlen)106);
	    }
	    s_wsfi(&io___171);
	    do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mjrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nfcon[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ktflag, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*prinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 1, (ftnlen)1);
	    do_fio(&c__1, ktflag + 2, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 3, (ftnlen)1);
	    do_fio(&c__1, (char *)&merit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nswap), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*condhz), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*pennrm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, flag__, (ftnlen)1);
	    do_fio(&c__1, mjrmsg, (ftnlen)8);
	    e_wsfi();
	    if (*condhz == 0.) {
		s_copy(buffp + 89, " ", (ftnlen)8, (ftnlen)1);
	    }
	    if (*pennrm == 0.) {
		s_copy(buffp + 97, " ", (ftnlen)8, (ftnlen)1);
	    }
	} else if (nlnobj) {
	    if (phead) {
/* Writing concatenation */
		i__1[0] = 97, a__1[0] = "   Itns Major Minors    Step   nObj"
			" Feasible  Optimal      Objective     L+U BSwap     "
			"nS  condHz";
		i__1[1] = 1, a__1[1] = cflag;
		s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)98);
		snprnt_(&c__11, ch__2, &iw[1], leniw, (ftnlen)98);
	    }
	    s_wsfi(&io___172);
	    do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mjrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nfobj[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ktflag + 2, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 3, (ftnlen)1);
	    do_fio(&c__1, (char *)&merit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nswap), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*condhz), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, flag__, (ftnlen)1);
	    do_fio(&c__1, mjrmsg, (ftnlen)8);
	    e_wsfi();
	    if (*condhz == 0.) {
		s_copy(buffp + 89, " ", (ftnlen)8, (ftnlen)1);
	    }
	} else {
	    if (phead) {
		if (*ns > 0) {
/* Writing concatenation */
		    i__1[0] = 89, a__1[0] = "   Itns Major Minors    Step   "
			    "     Feasible  Optimal    LPobjective     L+U BS"
			    "wap     nS";
		    i__1[1] = 1, a__1[1] = cflag;
		    s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)90);
		    snprnt_(&c__11, ch__3, &iw[1], leniw, (ftnlen)90);
		} else {
/* Writing concatenation */
		    i__1[0] = 89, a__1[0] = "   Itns Major Minors    Step   "
			    "     Feasible  Optimal    LPobjective     L+U   "
			    "          ";
		    i__1[1] = 1, a__1[1] = cflag;
		    s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)90);
		    snprnt_(&c__11, ch__3, &iw[1], leniw, (ftnlen)90);
		}
	    }
	    s_wsfi(&io___173);
	    do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mjrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nfobj[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ktflag + 2, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 3, (ftnlen)1);
	    do_fio(&c__1, (char *)&merit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nswap), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, flag__, (ftnlen)1);
	    do_fio(&c__1, mjrmsg, (ftnlen)8);
	    e_wsfi();
	    s_copy(buffp + 28, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*step == 0.) {
	    s_copy(buffp + 20, " ", (ftnlen)8, (ftnlen)1);
	}
	if (*nswap == 0) {
	    s_copy(buffp + 76, " ", (ftnlen)6, (ftnlen)1);
	}
	if (*ns == 0) {
	    s_copy(buffp + 82, " ", (ftnlen)7, (ftnlen)1);
	}
	snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)115);
    }
    if (summ1) {
/*        -------------------------------------------- */
/*        Terse line for the Summary file. */
/*        -------------------------------------------- */
	*(unsigned char *)mjrmsg = ' ';
	prthdg = iw[226] > 0;
	phead = mjrs % 10 == 0 || prthdg || major0;
	if (phead) {
	    iw[226] = 0;
	}
	if (nlncon) {
	    if (phead) {
		snprnt_(&c__12, " Major Minors     Step   nCon Feasible  Opt"
			"imal  MeritFunction    nS Penalty", &iw[1], leniw, (
			ftnlen)76);
	    }
	    s_wsfi(&io___174);
	    do_fio(&c__1, (char *)&mjrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nfcon[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ktflag, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*prinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 1, (ftnlen)1);
	    do_fio(&c__1, ktflag + 2, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 3, (ftnlen)1);
	    do_fio(&c__1, (char *)&merit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*pennrm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, mjrmsg, (ftnlen)8);
	    e_wsfi();
	    if (*pennrm == 0.) {
		s_copy(buffs + 68, " ", (ftnlen)8, (ftnlen)1);
	    }
	} else if (nlnobj) {
	    if (phead) {
		snprnt_(&c__12, " Major Minors     Step   nObj Feasible  Opt"
			"imal      Objective    nS", &iw[1], leniw, (ftnlen)68)
			;
	    }
	    s_wsfi(&io___175);
	    do_fio(&c__1, (char *)&mjrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nfobj[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ktflag + 2, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 3, (ftnlen)1);
	    do_fio(&c__1, (char *)&merit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, mjrmsg, (ftnlen)8);
	    e_wsfi();
	} else {
	    if (phead) {
		if (*ns > 0) {
		    snprnt_(&c__12, " Major Minors     Step        Feasible "
			    " Optimal    LPobjective    nS", &iw[1], leniw, (
			    ftnlen)68);
		} else {
/* Once zero, nS remains zero. */
		    snprnt_(&c__12, " Major Minors     Step        Feasible "
			    " Optimal    LPobjective", &iw[1], leniw, (ftnlen)
			    62);
		}
	    }
	    s_wsfi(&io___176);
	    do_fio(&c__1, (char *)&mjrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nfobj[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ktflag + 2, (ftnlen)1);
	    do_fio(&c__1, (char *)&(*duinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ktflag + 3, (ftnlen)1);
	    do_fio(&c__1, (char *)&merit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, mjrmsg, (ftnlen)8);
	    e_wsfi();
	    s_copy(buffs + 28, " ", (ftnlen)7, (ftnlen)1);
/* Zap nObj for LPs. */
	}
	if (*step == 0.) {
	    s_copy(buffs + 13, " ", (ftnlen)9, (ftnlen)1);
	}
	if (*ns == 0) {
	    s_copy(buffs + 62, " ", (ftnlen)6, (ftnlen)1);
	}
	snprnt_(&c__2, buffs, &iw[1], leniw, (ftnlen)84);
    }
    if (prntc && nncon > 0) {
/*        --------------------------------------------------------------- */
/*        Output heading for detailed log. */
/*        --------------------------------------------------------------- */
	s1page_(&c__0, &iw[1], leniw);
	if (major0) {
	    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
	}
/*        Unscale everything if necessary. */
	scaled = lvlscl >= 2;
	if (scaled) {
	    s2scla_(&c__1, &nncon, n, nb, iobj, &infbnd, sclobj, ne, nlocj, &
		    locj[1], &indj[1], &jcol[1], &ascale[1], &bl[1], &bu[1], &
		    ycon[1], &x[1]);
	    ddscl_(&nncon, &ascale[*n + 1], &c__1, &fcon[1], &c__1);
	}
	i__ = 0;
/* Keeps ftnchek happy */
	l = *mjrprt / 100;
	prtx = l % 10 > 0;
	l /= 10;
	prtl = l % 10 > 0 && mjrs > 0;
	l /= 10;
	prtf = l % 10 > 0;
	l /= 10;
	prtj = l % 10 > 0;
	s_copy(form, "(1p, i9, e13.5)", (ftnlen)40, (ftnlen)15);
	if (prtx) {
	    snprnt_(&c__11, " Jacobian variables", &iw[1], leniw, (ftnlen)19);
	    snprnt_(&c__1, " ------------------", &iw[1], leniw, (ftnlen)19);
	    i__2 = nnjac;
	    for (j = 1; j <= i__2; ++j) {
		s_wsfi(&io___186);
		i__3 = s2coln_(&j, leniw, &iw[1]);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
	if (prtl) {
	    snprnt_(&c__11, " Multiplier estimates", &iw[1], leniw, (ftnlen)
		    21);
	    snprnt_(&c__1, " --------------------", &iw[1], leniw, (ftnlen)21)
		    ;
	    i__2 = nncon;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s_wsfi(&io___187);
		i__3 = s2rown_(&i__, leniw, &iw[1]);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ycon[i__], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
	if (prtf) {
	    snprnt_(&c__11, " Constraint functions", &iw[1], leniw, (ftnlen)
		    21);
	    snprnt_(&c__1, " --------------------", &iw[1], leniw, (ftnlen)21)
		    ;
	    i__2 = nncon;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s_wsfi(&io___188);
		i__3 = s2rown_(&i__, leniw, &iw[1]);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fcon[i__], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	    }
	    s_wsfi(&io___189);
	    do_fio(&c__1, (char *)&(*vimax), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*virel), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
	}
	if (prtj) {
	    snprnt_(&c__11, " x  and  Jacobian", &iw[1], leniw, (ftnlen)17);
	    snprnt_(&c__1, " ----------------", &iw[1], leniw, (ftnlen)17);
	    i__2 = nnjac;
	    for (j = 1; j <= i__2; ++j) {
		l = hs[j];
		s_wsfi(&io___190);
		i__3 = s2coln_(&j, leniw, &iw[1]);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, key + (l + 1 << 1), (ftnlen)2);
		e_wsfi();
		snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
		k1 = locj[j];
		k2 = locj[j + 1] - 1;
		i__3 = k2;
		for (k = k1; k <= i__3; ++k) {
		    ir = indj[k];
		    if (ir > nncon) {
			goto L160;
		    }
		    s_wsfi(&io___194);
		    i__4 = s2rown_(&indj[k], leniw, &iw[1]);
		    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&jcol[k], (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
		}
L160:
		;
	    }
	}
/*        Scale again if necessary. */
	if (scaled) {
	    s2scla_(&c__0, &nncon, n, nb, iobj, &infbnd, sclobj, ne, nlocj, &
		    locj[1], &indj[1], &jcol[1], &ascale[1], &bl[1], &bu[1], &
		    ycon[1], &x[1]);
	    dddiv_(&nncon, &ascale[*n + 1], &c__1, &fcon[1], &c__1);
	}
    }
    return 0;
/*     Major log,  Print file. */
/* nlnCo */
/* nlnOb */
/*     Major log,  Summary file. */
/* LP wi */
} /* snlog_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snLog */
/* Subroutine */ int snlog2_(integer *prob, char *probtag, logical *elastc, 
	logical *gotr, logical *jstfea, logical *feasbl, integer *m, integer *
	mbs, integer *nnh, integer *ns, integer *jsq, integer *jbr, integer *
	jsr, integer *linesp, integer *liness, integer *itn, integer *itqp, 
	integer *kprc, integer *lvlinf, doublereal *pivot, doublereal *step, 
	integer *ninf, doublereal *sinf, doublereal *wtinf, doublereal *
	objprt, doublereal *condhz, doublereal *djqprt, doublereal *rgnorm, 
	integer *kbs, doublereal *xbs, integer *iw, integer *leniw, ftnlen 
	probtag_len)
{
    /* Format strings */
    static char fmt_8010[] = "(\002 Itn\002,i7,\002: Feasible \002,a)";
    static char fmt_8030[] = "(\002 Itn\002,i7,\002: Elastic Phase 2 -- mini"
	    "mizing\002,\002 elastic variables\002)";
    static char fmt_8040[] = "(\002 Itn\002,i7,\002: Elastic Phase 2 -- mini"
	    "mizing\002,\002 obj + weighted elastics\002)";
    static char fmt_3000[] = "(1p,i7,4x,e9.1,e8.1,i7,2e9.1,e15.7,3i7,e9.1,i8"
	    ",i4,i6,e8.1,i7)";
    static char fmt_5000[] = "(1p,i13,e9.1,i7,e9.1,e9.1,e15.7,i6,i7)";
    static char fmt_6000[] = "(i7,g17.8)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer k, ncp;
    static char str[80];
    static integer jbrn, lenl, lenu, itns, jsqn, jsrn, mnrs;
    static logical phead;
    static char buffp[138];
    static integer cgitn;
    static char buffs[75];
    static integer width;
    extern integer s2varn_(integer *, integer *, integer *);
    static integer lprdbg, qpmode;
    static logical prthdp, prthds, prtlog, newset;
    static doublereal mxwdth;
    static integer printp, prints;
    static doublereal rmxint;
    static integer mnrprt;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static logical prtsum;
    static integer lvlsys;
    extern integer s1intmx_(void);

    /* Fortran I/O blocks */
    static icilist io___215 = { 0, str, 0, fmt_8010, 80, 1 };
    static icilist io___216 = { 0, str, 0, fmt_8030, 80, 1 };
    static icilist io___217 = { 0, str, 0, fmt_8040, 80, 1 };
    static icilist io___225 = { 0, buffp, 0, fmt_3000, 138, 1 };
    static icilist io___226 = { 0, buffp, 0, fmt_3000, 138, 1 };
    static icilist io___227 = { 0, buffs, 0, fmt_5000, 75, 1 };
    static icilist io___228 = { 0, buffs, 0, fmt_5000, 75, 1 };
    static icilist io___230 = { 0, buffp, 0, fmt_6000, 138, 1 };
    static icilist io___231 = { 0, buffp, 0, fmt_6000, 138, 1 };


/*     ================================================================== */
/*     snLog2  prints the minor iteration log. */

/*     mBS = m + maxS  .ge.  m + nS. */

/*     MnrHdP is  1  if a new heading is required for some reason other */
/*     than frequency (e.g., after a basis factorization). */

/*     The output consists of a number of ``sections'' of one line */
/*     summaries, with each section preceded by a header message. */
/*     linesP and linesS count the number of lines remaining to be */
/*     printed in each section of the print and summary files */
/*     respectively.   They too may force a new heading. */

/*     01 Dec 1991: First version based on Minos routine m5log. */
/*     17 Nov 2000: First version for SQP minor itns. */
/*     02 Dec 2000: Reordered and sparsified. */
/*     28 Dec 2000: Row and column permutations added. */
/*     01 Aug 2003: cgItn added to Print log. */
/*     04 Jul 2005: Current version of snLog2. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* >0 => Minor heading for Log */
/* >0 => Major heading for Log */
/* >0 => Minor heading for Summary */
/*     ------------------------------------------------------------------ */
/* >0 => Major heading for Summary */
    /* Parameter adjustments */
    --xbs;
    --kbs;
    --iw;

    /* Function Body */
    lvlsys = iw[71];
/* > 0   => print system info */
    lprdbg = iw[85];
/* > 0   => private debug print */
    mnrprt = iw[93];
/* Minor print level */
    ncp = iw[176];
/* no. of LU compressions */
    lenl = iw[173];
/* size of current  L */
    lenu = iw[174];
/* size of current  U */
    qpmode = iw[208];
/* Current QP solver   (based on QPslvr) */
    printp = iw[218];
/* (on/off) current log     file output status */
    prints = iw[219];
/* (on/off) current summary file output status */
    cgitn = iw[387];
/* symmlq itns for the last QP minor itn */
    prtlog = printp == 1;
    prtsum = prints == 1;
    mxwdth = 1e7;
/* Integers printed i7 */
    rmxint = (doublereal) s1intmx_();
/* Largest integer without overflow */
    if (mxwdth > rmxint) {
	mxwdth = rmxint;
    }
    width = (integer) mxwdth;
    itns = *itn % width;
    mnrs = *itqp % width;
    s_copy(buffp, " ", (ftnlen)138, (ftnlen)1);
    s_copy(buffs, " ", (ftnlen)75, (ftnlen)1);
/* If  newly feasible, print something. */
    if (*jstfea && mnrprt >= 10) {
	if (! (*elastc)) {
/* Constraints feasible in Normal mode. */
/* Print a message. */
/* ProbTag is one of the following: */
/* ProbTag = 'QP problem' */
/* ProbTag = 'LP problem' */
/* ProbTag = 'QP subproblem' */
/* ProbTag = 'Linear constraints' */
	    if (*prob != 4 && *prob != 0 && *prob != 3) {
		s_wsfi(&io___215);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, probtag, (ftnlen)20);
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    }
	} else {
/* Elastic mode */
/* Elastic Phase 1 has completed. */
	    if (*lvlinf == 2) {
/* Infinite weight on sumInf. */
/* Minimize the infeasible elastics. */
		s_wsfi(&io___216);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    } else if (*lvlinf == 1) {
/* Finite nonzero weight on sumInf */
/* Minimize a weighted objective. */
		s_wsfi(&io___217);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
	if (lvlsys > 0) {
	    iw[223] = 1;
/* Print the header to the print   file */
	    iw[225] = 1;
/* Print the header to the summary file */
	}
    }
    prthdp = iw[223] > 0;
    prthds = iw[225] > 0;
    if (prtlog) {
/*        -------------------------------------- */
/*        Terse line for the Print file. */
/*        -------------------------------------- */
	newset = *linesp == 0;
	phead = prthdp || newset;
	if (phead) {
	    iw[223] = 0;
	    *linesp = 40;
	}
	iw[224] = 1;
	--(*linesp);
	jsqn = s2varn_(jsq, leniw, &iw[1]);
	jsrn = s2varn_(jsr, leniw, &iw[1]);
	jbrn = s2varn_(jbr, leniw, &iw[1]);
	if (*nnh > 0) {
	    if (phead) {
		s_copy(buffp, "    Itn       QPmult  QPstep   nInf   SumInf "
			"  rgNorm    QPobjective   +SBS   -SBS    -BS    Pivo"
			"t     L+U ncp    nS  condHz", (ftnlen)138, (ftnlen)
			124);
		if (*elastc) {
		    s_copy(buffp + 55, "Elastic QPobj", (ftnlen)13, (ftnlen)
			    13);
		}
		if (qpmode == 1) {
		    s_copy(buffp + 125, "cgItns", (ftnlen)6, (ftnlen)6);
		}
		snprnt_(&c__11, buffp, &iw[1], leniw, (ftnlen)138);
	    }
	    s_wsfi(&io___225);
	    do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&jsqn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jsrn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jbrn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
	    i__1 = lenl + lenu;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*condhz), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&cgitn, (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
/* nnH == 0 */
	    if (phead) {
		s_copy(buffp, "    Itn       LPmult  LPstep   nInf   SumInf "
			"            LPobjective   +SBS   -SBS    -BS    Pivo"
			"t     L+U ncp", (ftnlen)138, (ftnlen)110);
		if (*elastc) {
		    s_copy(buffp + 55, "Elastic LPobj", (ftnlen)13, (ftnlen)
			    13);
		}
		if (*ns > 0) {
		    s_copy(buffp + 114, "nS", (ftnlen)2, (ftnlen)2);
		}
		snprnt_(&c__11, buffp, &iw[1], leniw, (ftnlen)138);
	    }
	    s_wsfi(&io___226);
	    do_fio(&c__1, (char *)&itns, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*step), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&jsqn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jsrn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jbrn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
	    i__1 = lenl + lenu;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	if (*djqprt == 0.) {
	    s_copy(buffp + 12, " ", (ftnlen)8, (ftnlen)1);
	}
	if (*step == 0.) {
	    s_copy(buffp + 20, " ", (ftnlen)8, (ftnlen)1);
	}
	if (*ninf == 0) {
	    s_copy(buffp + 28, " ", (ftnlen)16, (ftnlen)1);
	}
/* nInf and sInf */
	if (*rgnorm == 0.) {
	    s_copy(buffp + 44, " ", (ftnlen)9, (ftnlen)1);
	}
	if (! (*feasbl)) {
	    s_copy(buffp + 53, " ", (ftnlen)15, (ftnlen)1);
	}
	if (*jsq == 0) {
	    s_copy(buffp + 68, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*jsr == 0) {
	    s_copy(buffp + 75, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*jbr == 0) {
	    s_copy(buffp + 82, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*pivot == 0.) {
	    s_copy(buffp + 89, " ", (ftnlen)9, (ftnlen)1);
	}
	if (ncp == 0) {
	    s_copy(buffp + 106, " ", (ftnlen)4, (ftnlen)1);
	}
	if (*ns == 0) {
	    s_copy(buffp + 110, " ", (ftnlen)6, (ftnlen)1);
	}
	if (*condhz == 0.) {
	    s_copy(buffp + 116, " ", (ftnlen)8, (ftnlen)1);
	}
/* condHz */
	if (cgitn == 0) {
	    s_copy(buffp + 124, " ", (ftnlen)7, (ftnlen)1);
	}
	snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)138);
    }
    if (prtsum) {
/*        -------------------------------- */
/*        Terse line for the Summary file. */
/*        -------------------------------- */
	newset = *liness == 0;
	phead = prthds || newset;
	if (phead) {
	    iw[225] = 0;
	    *liness = 10;
	}
	iw[226] = 1;
	--(*liness);
	if (*nnh > 0) {
	    if (phead) {
		s_copy(buffs, "        Minor   QPmult   nInf   SumInf   rgNo"
			"rm    QPobjective", (ftnlen)75, (ftnlen)62);
		if (*elastc) {
		    s_copy(buffs + 49, "Elastic QPobj", (ftnlen)13, (ftnlen)
			    13);
		}
		if (*ns > 0) {
		    s_copy(buffs + 66, "nS", (ftnlen)2, (ftnlen)2);
		}
		if (qpmode == 1) {
		    s_copy(buffs + 69, "cgItns", (ftnlen)6, (ftnlen)6);
		}
		snprnt_(&c__12, buffs, &iw[1], leniw, (ftnlen)75);
	    }
	    s_wsfi(&io___227);
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&cgitn, (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
/* nnH == 0 */
	    if (phead) {
		s_copy(buffs, "        Minor   LPmult   nInf   SumInf   rgNo"
			"rm    LPobjective", (ftnlen)75, (ftnlen)62);
		if (*elastc) {
		    s_copy(buffs + 49, "Elastic LPobj", (ftnlen)13, (ftnlen)
			    13);
		}
		if (*ns > 0) {
		    s_copy(buffs + 66, "nS", (ftnlen)2, (ftnlen)2);
		}
		snprnt_(&c__12, buffs, &iw[1], leniw, (ftnlen)75);
	    }
	    s_wsfi(&io___228);
	    do_fio(&c__1, (char *)&mnrs, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*djqprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*objprt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	if (*djqprt == 0.) {
	    s_copy(buffs + 13, " ", (ftnlen)9, (ftnlen)1);
	}
	if (*ninf == 0) {
	    s_copy(buffs + 22, " ", (ftnlen)7, (ftnlen)1);
	}
	if (*sinf == 0.) {
	    s_copy(buffs + 29, " ", (ftnlen)9, (ftnlen)1);
	}
	if (*rgnorm == 0.) {
	    s_copy(buffs + 40, " ", (ftnlen)7, (ftnlen)1);
	}
	if (! (*feasbl)) {
	    s_copy(buffs + 47, " ", (ftnlen)15, (ftnlen)1);
	}
	if (*ns == 0) {
	    s_copy(buffs + 62, " ", (ftnlen)6, (ftnlen)1);
	}
	if (cgitn == 0) {
	    s_copy(buffs + 68, " ", (ftnlen)7, (ftnlen)1);
	}
	snprnt_(&c__2, buffs, &iw[1], leniw, (ftnlen)75);
    }
/*     ------------------------------------------------------------------ */
/*     Debug output. */
/*     ------------------------------------------------------------------ */
    if (lprdbg == 100) {
	snprnt_(&c__11, " BS values...", &iw[1], leniw, (ftnlen)13);
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    s_wsfi(&io___230);
	    i__2 = s2varn_(&kbs[k], leniw, &iw[1]);
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xbs[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)138);
	}
	snprnt_(&c__11, " SB values...", &iw[1], leniw, (ftnlen)13);
	i__1 = *m + *ns;
	for (k = *m + 1; k <= i__1; ++k) {
	    s_wsfi(&io___231);
	    i__2 = s2varn_(&kbs[k], leniw, &iw[1]);
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xbs[k], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__1, buffp, &iw[1], leniw, (ftnlen)138);
	}
    }
    return 0;
/* Minor log,  Print   file. */
/* Minor log,  Summary file. */
} /* snlog2_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snLog2 */
/* Subroutine */ int snstop_(integer *iabort, integer *info, integer *hqntype,
	 logical *ktcond, integer *mjrprt, integer *minimz, integer *m, 
	integer *maxs, integer *n, integer *nb, integer *nncon0, integer *
	nncon, integer *nnobj0, integer *nnobj, integer *ns, integer *itn, 
	integer *nmajor, integer *nminor, integer *nswap, doublereal *condhz, 
	integer *iobj, doublereal *sclobj, doublereal *objadd, doublereal *
	fmrt, doublereal *pennrm, doublereal *step, doublereal *prinf, 
	doublereal *duinf, doublereal *vimax, doublereal *virel, integer *hs, 
	integer *ne, integer *nlocj, integer *locj, integer *indj, doublereal 
	*jcol, integer *negcon, doublereal *ascale, doublereal *bl, 
	doublereal *bu, doublereal *fcon, doublereal *gcon, doublereal *gobj, 
	doublereal *ycon, doublereal *pi, doublereal *rc, doublereal *rg, 
	doublereal *x, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
/*     ================================================================== */
/*     snSTOP is called every major iteration. */
/*     If iAbort > 0 on exit, the run is terminated. */
/*     By specifying a custom version of snSTOP, the user can arrange for */
/*     snopt to be terminated at any given major iteration. */

/*     14 Oct 2004: First version of   snSTOP. */
/*     29 Aug 2007: Parameter list extended. */
/*     ================================================================== */
    /* Parameter adjustments */
    --info;
    --ktcond;
    --pi;
    --rg;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --ycon;
    --fcon;
    --gobj;
    --jcol;
    --indj;
    --locj;
    --gcon;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    *iabort = 0;
/*     Relax */
    return 0;
} /* snstop_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snSTOP */
/* Subroutine */ int snset_(char *buffer, integer *iprint, integer *isumm, 
	integer *errors, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen buffer_len, ftnlen 
	cw_len)
{
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalue;

/*     ================================================================== */
/*     snSet  decodes the option contained in  buffer. */

/*     The buffer is output to file iPrint, minus trailing blanks. */
/*     Error messages are output to files iPrint and iSumm. */
/*     Buffer is echoed to iPrint but normally not to iSumm. */
/*     It is echoed to iSumm before any error msg. */

/*     On entry, */
/*     iPrint is the print   file.  no output occurs if iPrint .le 0. */
/*     iSumm  is the Summary file.  no output occurs if iSumm  .le 0. */
/*     Errors is the number of errors so far. */

/*     On exit, */
/*     Errors is the number of errors so far. */

/*     27 Nov 1991: first version of snSet. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_true, buffer, key, cvalue, &ivalue, &rvalue, iprint, isumm, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snset_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snSet */
/* Subroutine */ int snseti_(char *buffer, integer *ivalue, integer *iprint, 
	integer *isumm, integer *errors, char *cw, integer *lencw, integer *
	iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen buffer_len,
	 ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char buff72[72];
    static integer lenbuf;
    static char cvalue[8];
    static doublereal rvalue;
    static integer ivalxx;

    /* Fortran I/O blocks */
    static icilist io___237 = { 0, key, 0, "(i16)", 16, 1 };


/*     ================================================================== */
/*     snSeti decodes the option contained in  buffer // ivalue. */
/*     The parameters other than ivalue are as in snSet. */

/*     27 Nov 1991: first version of snSeti. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_wsfi(&io___237);
    do_fio(&c__1, (char *)&(*ivalue), (ftnlen)sizeof(integer));
    e_wsfi();
    lenbuf = i_len(buffer, buffer_len);
    s_copy(buff72, buffer, (ftnlen)72, buffer_len);
    i__1 = lenbuf;
    s_copy(buff72 + i__1, key, lenbuf + 16 - i__1, (ftnlen)16);
    ivalxx = *ivalue;
    s3opt_(&c_true, buff72, key, cvalue, &ivalxx, &rvalue, iprint, isumm, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)72, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snseti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snSeti */
/* Subroutine */ int snsetr_(char *buffer, doublereal *rvalue, integer *
	iprint, integer *isumm, integer *errors, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	buffer_len, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char buff72[72];
    static integer lenbuf;
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalxx;

    /* Fortran I/O blocks */
    static icilist io___244 = { 0, key, 0, "(1p, e16.8)", 16, 1 };


/*     ================================================================== */
/*     snSetr decodes the option contained in  buffer // rvalue. */
/*     The parameters other than rvalue are as in snSet. */

/*     27 Nov 1991: first version of snSetr. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_wsfi(&io___244);
    do_fio(&c__1, (char *)&(*rvalue), (ftnlen)sizeof(doublereal));
    e_wsfi();
    lenbuf = i_len(buffer, buffer_len);
    s_copy(buff72, buffer, (ftnlen)72, buffer_len);
    i__1 = lenbuf;
    s_copy(buff72 + i__1, key, lenbuf + 16 - i__1, (ftnlen)16);
    rvalxx = *rvalue;
    s3opt_(&c_true, buff72, key, cvalue, &ivalue, &rvalxx, iprint, isumm, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)72, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* snsetr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snSetr */
integer snget_(char *buffer, integer *errors, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	buffer_len, ftnlen cw_len)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalue;

/*     ================================================================== */
/*     snGet  decodes the option contained in  buffer */
/*     and returns 1 if the option has previously been set, else 0. */
/*     For example, */
/*     i = snGet ( 'Maximize', Errors, cw, lencw, iw, leniw, rw, lenrw ) */

/*     01 Aug 2003: First version of snGet.  Needed because */
/*                  snGetc, snGeti, snGetr were not well defined */
/*                  for strings that had no numerical value. */
/*     01 Aug 2003: Current version of snGet. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ivalue, &rvalue, &c__0, &c__0, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    ret_val = ivalue;
    return ret_val;
} /* snget_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function snGet */
/* Subroutine */ int sngetc_(char *buffer, char *cvalue, integer *errors, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen buffer_len, ftnlen cvalue_len, ftnlen cw_len)
{
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static integer ivalue;
    static doublereal rvalue;

/*     ================================================================== */
/*     snGetc gets the value of the option contained in  buffer. */
/*     The parameters other than cvalue are as in snSet. */

/*     17 May 1998: first version of snGetc. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ivalue, &rvalue, &c__0, &c__0, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* sngetc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snGetc */
/* Subroutine */ int sngeti_(char *buffer, integer *ivalue, integer *errors, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen buffer_len, ftnlen cw_len)
{
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static doublereal rvalue;

/*     ================================================================== */
/*     snGeti gets the value of the option contained in  buffer. */
/*     The parameters other than ivalue are as in snSet. */

/*     17 May 1998: first version of snGeti. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, ivalue, &rvalue, &c__0, &c__0, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* sngeti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snGeti */
/* Subroutine */ int sngetr_(char *buffer, doublereal *rvalue, integer *
	errors, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen buffer_len, ftnlen cw_len)
{
    static char key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;

/*     ================================================================== */
/*     snGetr gets the value of the option contained in  buffer. */
/*     The parameters other than rvalue are as in snSet. */

/*     17 May 1998: first version of snGetr. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ivalue, rvalue, &c__0, &c__0, 
	    errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* sngetr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snGetr */
/* Subroutine */ int snreth_(integer *errors, integer *lenh, doublereal *h__, 
	integer *nnh, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_9990[] = "(\002 XXX  lenH too small: needs to be at le"
	    "ast\002,i6)";
    static char fmt_9991[] = "(\002 XXX  Full-memory Hessian not requeste"
	    "d\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer lu, ly, ly1;
    static char str[80], str2[80];
    static integer lenu, iexit;
    extern /* Subroutine */ int s4chkp_(integer *, char *, integer *, integer 
	    *, integer *, ftnlen), s8geth_(integer *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    static integer lvlhes;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), snprnt_(integer *, 
	    char *, integer *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___272 = { 0, str, 0, fmt_9990, 80, 1 };
    static icilist io___273 = { 0, str, 0, fmt_9991, 80, 1 };


/*     ================================================================== */
/*     snRetH  retrieves the SNOPT approximate Hessian from memory. */

/*     snRetH must be called immediately after a call to npOpt, snoptB */
/*     or snoptC  with the option "Hessian full memory" selected. */

/*     On entry: */
/*        lenH      (ge 1) is the length the array H, i.e., H(lenH). */

/*     On exit */
/*        Errors    contains the number of errors. */
/*                  if Errors gt 0, then error messages are printed on */
/*                  the standard output */
/*        H(lenH)   contains the upper triangular part of the approximate */
/*                  Hessian, stored by columns. */
/*        nnH       contains the number of columns of H. */


/*     03 Sep 2006: First version of snRetH */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --h__;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "SNRETH", (ftnlen)6, (ftnlen)6);
    *errors = 0;
    if (*lencw < 500 || *leniw < 500 || *lenrw < 500) {
/*        --------------------------------------------------------------- */
/*        Not enough workspace to do ANYTHING! */
/*        Print and exit without accessing the work arrays. */
/*        --------------------------------------------------------------- */
	iexit = 81;
/* Work arrays must have at least 500 elements */
	++(*errors);
	snwrap_(&iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
/*     ------------------------------------------------------------------ */
/*     Retrieve the approximate Hessian. */
/*     ------------------------------------------------------------------ */
    lvlhes = iw[72];
/* LM, FM or Exact Hessian */
    *nnh = iw[24];
/* max( nnObj, nnJac ) */
    ly = iw[311];
/* y (nb)      =  real work vector */
    ly1 = iw[312];
/* y1(nb)      =  real work vector */
    lu = iw[391];
/* U(lenU), BFGS Hessian H = U'U */
    lenu = iw[392];
/*     Check pointers, lengths, etc,  retrieved from memory. */

    s4chkp_(errors, "lvlHes", &lvlhes, &iw[1], leniw, (ftnlen)6);
    s4chkp_(errors, "    ly", &ly, &iw[1], leniw, (ftnlen)6);
    s4chkp_(errors, "   ly1", &ly1, &iw[1], leniw, (ftnlen)6);
    s4chkp_(errors, "    lU", &lu, &iw[1], leniw, (ftnlen)6);
    if (*lenh < lenu) {
	++(*errors);
	s_wsfi(&io___272);
	do_fio(&c__1, (char *)&lenu, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*lenh), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__5, str, &iw[1], leniw, (ftnlen)80);
    }
    if (lvlhes != 1) {
	++(*errors);
	s_wsfi(&io___273);
	e_wsfi();
	snprnt_(&c__5, str, &iw[1], leniw, (ftnlen)80);
    }
    if (*errors == 0) {
	s8geth_(nnh, lenh, &rw[lu], &h__[1], &rw[ly], &rw[ly1]);
    }
L999:
    return 0;
} /* snreth_ */

