/* ./src/np02lib.f -- translated by f2c (version 20100827).
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
static integer c__11 = 11;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__12 = 12;
static integer c__130 = 130;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__0 = 0;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  np02lib.f */

/*     npTitl   npInit   npSpec   npMem */
/*     npSet    npSeti   npSetr */
/*     npGet    npGeti   npGetr */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int nptitl_(char *title, ftnlen title_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

/*     ================================================================== */
/*     npTitl sets the title. */
/*     ================================================================== */
    s_copy(title, "N P O P T  7.2-12.2 (Jul 2013)", (ftnlen)30, (ftnlen)30);
/* ---------------123456789|123456789|123456789|-------------------------- */
    return 0;
} /* nptitl_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npTitl */
/* Subroutine */ int npinit_(integer *iprint, integer *isumm, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
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
    static char cw[8*500];
    extern /* Subroutine */ int s3unsetall_(char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static char str[80], str2[80], title[30];
    static integer istdo;
    extern /* Subroutine */ int s1init_(char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    static integer ispecs, inform__;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), nptitl_(char *, 
	    ftnlen), snprnt_(integer *, char *, integer *, integer *, ftnlen);
    extern integer s1outpt_(void);

/*     ================================================================== */
/*     npInit  is called by the user to do the following: */
/*     1. Open default files (Print, Summary). */
/*     2. Initialize title. */
/*     3. Set options to default values. */

/*     15 Nov 1991: First version. */
/*     14 Jul 1997: Thread-safe version. */
/*     21 Mar 1997: First version based on snopt routine snInit */
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
/* Timing level */
/* QP user-routine call-status */
/* NP user-routine call-status */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    s_copy(solver, "NPINIT", (ftnlen)6, (ftnlen)6);
    if (*leniw < 500 || *lenrw < 500) {
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
    s3unsetall_(cw, &c__500, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
/* ----------------------------------------------------------------- */
/* Initialize some default values. */
/* ----------------------------------------------------------------- */
    ispecs = 0;
    istdo = s1outpt_();
    iw[10] = istdo;
/* Standard Output */
    iw[11] = ispecs;
    iw[12] = *iprint;
/* Print file */
    iw[13] = *isumm;
/* Summary file */
    iw[6] = 500;
    iw[4] = 500;
    iw[2] = 500;
    iw[7] = 500;
    iw[5] = *leniw;
    iw[3] = *lenrw;
    nptitl_(title, (ftnlen)30);
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
} /* npinit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npInit */
/* Subroutine */ int npspec_(integer *ispecs, integer *iexit, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char cw[8*500], str[80], str2[80];
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
/*     npSpec  is called by the user to read the Specs file. */

/*     07 Feb 1998: First version. */
/*     01 Aug 2003: s3file now has a "title" parameter.  Use ' '. */
/*     27 Oct 2003: Current version of npSpec. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "NPSPEC", (ftnlen)6, (ftnlen)6);
    if (*leniw < 500 || *lenrw < 500) {
/* -------------------------------------------------------------- */
/* Not enough workspace to do ANYTHING! */
/* Print and exit without accessing the work arrays. */
/* -------------------------------------------------------------- */
	*iexit = 81;
/* Work arrays must have at least 500 elements */
	snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)
		80, (ftnlen)80);
	goto L999;
    }
    if (*ispecs <= 0) {
	*iexit = 131;
	goto L800;
    }
    iw[11] = *ispecs;
/* Specs (options) file */
    iprint = iw[12];
/* Print file */
    isumm = iw[13];
/* Summary file */
    *iexit = 0;
    calls = 1;
/*     ------------------------------------------------------------------ */
/*     Read the Specs file. */
/*     npopt  will check the options later and maybe print them. */
/*     ------------------------------------------------------------------ */
    s3file_(iexit, &calls, ispecs, (U_fp)s3opt_, " ", &iprint, &isumm, &
	    errors, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, (ftnlen)1, (
	    ftnlen)8);
L800:
    if (*iexit == 0) {
	*iexit = 101;
/* SPECS file read successfully */
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* npspec_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npSpec */
/* Subroutine */ int npmem_(integer *iexit, integer *n, integer *nclin, 
	integer *ncnln, integer *mincw, integer *miniw, integer *minrw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, m, ne;
    static char cw[8*500];
    static integer nkx;
    static char str[80], str2[80];
    static integer ncon, lenr, maxr, maxs;
    extern /* Subroutine */ int s2mem_(integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), s8map_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer nnjac, nnobj, nncol, nncon, maxcw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxiw, maxrw;
    extern /* Subroutine */ int s2bmap_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), s8dflt_(integer *, integer *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static integer negcon, llencw, inform__, lleniw, mqnmod, lvlhes, llenrw;
    static logical prtmem;
    static integer liwest, nextcw;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer useriw[130], nextiw, lrwest;
    static doublereal userrw[130];
    static integer nextrw;

/*     ================================================================== */
/*     npMem   estimates the memory requirements for npopt, */
/*     using the values: */
/*     n      the number of variables (dimension of  x), */

/*     nclin  the number of linear constraints (rows of the matrix  A), */

/*     ncnln  the number of nonlinear constraints (dimension of  c(x)), */

/*     These values are used to compute the minimum required storage: */
/*     miniw, minrw. */

/*     Note: */
/*     1. All default parameters must be set before calling npMem, */
/*        since some values affect the amount of memory required. */

/*     2. The arrays rw and iw hold  constants and work-space addresses. */
/*        They must have dimension at least 500. */

/*     3. This version of npMem does not allow user accessible */
/*        partitions of iw and rw. */

/*     01 May 1998: First version. */
/*     19 Feb 2004: Current version of npMem. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s_copy(solver, "NPMEM ", (ftnlen)6, (ftnlen)6);
    *iexit = 0;
    if (*leniw < 500 || *lenrw < 500) {
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
    icopy_(&c__130, &iw[51], &c__1, useriw, &c__1);
    dcopy_(&c__130, &rw[51], &c__1, userrw, &c__1);
/*     Assign fake values for lencw, leniw, lenrw. */
/*     This will force s2Mem to estimate the memory requirements. */
    llenrw = 500;
    lleniw = 500;
    llencw = 500;
/*     Compute the problem dimensions. */
    ncon = *nclin + *ncnln;
    if (ncon == 0) {
/*        The problem is unconstrained. */
/*        A dummy row of zeros will be included. */
	nncol = 0;
	m = 1;
	ne = 1;
    } else {
	nncol = *n;
	m = ncon;
	ne = m * *n;
    }
    negcon = *ncnln * *n;
    nncon = *ncnln;
    nnjac = nncol;
    nnobj = *n;
/*     An obligatory call to snInit has `undefined' all options. */
/*     However, it could not undefine the char*8 options.  Do it now. */
/*     Check the user-defined values and assign undefined values. */
/*     s8dflt needs various problem dimensions in iw. */
    for (i__ = 51; i__ <= 180; ++i__) {
	s_copy(cw + (i__ - 1 << 3), "-1111111", (ftnlen)8, (ftnlen)8);
    }
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
    s8dflt_(&m, n, &nncon, &nnjac, &nnobj, cw, &llencw, &iw[1], &lleniw, &rw[
	    1], &llenrw, (ftnlen)8);
    nextcw = 501;
    nextiw = 501;
    nextrw = 501;
    maxcw = 500;
    maxiw = *leniw;
    maxrw = *lenrw;
    nkx = *n + m;
    maxr = iw[52];
/* max columns of R. */
    maxs = iw[53];
/* max # of superbasics */
    mqnmod = iw[54];
/* (ge 0) max # of BFGS updates */
    lvlhes = iw[72];
/* 0,1,2  => LM, FM, Exact Hessian */
    lenr = maxr * (maxr + 1) / 2 + (maxs - maxr);
    s8map_(&m, n, &negcon, &nkx, &nncon, &nnjac, &nnobj, &lenr, &maxr, &maxs, 
	    &mqnmod, &lvlhes, &nextcw, &nextiw, &nextrw, &iw[1], leniw);
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
} /* npmem_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npMem */
/* Subroutine */ int npset_(char *buffer, integer *iprint, integer *isumm, 
	integer *iexit, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, ftnlen buffer_len)
{
    static char cw[8*500], key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalue;

/*     ================================================================== */
/*     npSet  decodes the option contained in  buffer. */

/*     The buffer is output to file iPrint, minus trailing blanks. */
/*     Error messages are output to files iPrint and iSumm. */
/*     Buffer is echoed to iPrint but normally not to iSumm. */
/*     It is echoed to iSumm before any error msg. */

/*     On entry, */
/*     iPrint is the print   file.  no output occurs if iPrint .le 0. */
/*     iSumm  is the Summary file.  no output occurs if iSumm  .le 0. */
/*     iExit  is the number of errors so far. */

/*     On exit, */
/*     iExit  is the number of errors so far. */

/*     27 Nov 1991: first version of npSet. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_true, buffer, key, cvalue, &ivalue, &rvalue, iprint, isumm, 
	    iexit, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* npset_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npSet */
/* Subroutine */ int npseti_(char *buffer, integer *ivalue, integer *iprint, 
	integer *isumm, integer *iexit, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen buffer_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char cw[8*500], key[16];
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
    static icilist io___58 = { 0, key, 0, "(i16)", 16, 1 };


/*     ================================================================== */
/*     npSeti decodes the option contained in  buffer // ivalue. */
/*     The parameters other than ivalue are as in npSet. */

/*     27 Nov 1991: first version of npSeti. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s_wsfi(&io___58);
    do_fio(&c__1, (char *)&(*ivalue), (ftnlen)sizeof(integer));
    e_wsfi();
    lenbuf = i_len(buffer, buffer_len);
    s_copy(buff72, buffer, (ftnlen)72, buffer_len);
    i__1 = lenbuf;
    s_copy(buff72 + i__1, key, lenbuf + 16 - i__1, (ftnlen)16);
    ivalxx = *ivalue;
    s3opt_(&c_true, buff72, key, cvalue, &ivalxx, &rvalue, iprint, isumm, 
	    iexit, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, (ftnlen)72, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* npseti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npSeti */
/* Subroutine */ int npsetr_(char *buffer, doublereal *rvalue, integer *
	iprint, integer *isumm, integer *iexit, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen buffer_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char cw[8*500], key[16];
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
    static icilist io___66 = { 0, key, 0, "(1p, e16.8)", 16, 1 };


/*     ================================================================== */
/*     npSetr decodes the option contained in  buffer // rvalue. */
/*     The parameters other than rvalue are as in npSet. */

/*     27 Nov 1991: first version of npSetr. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s_wsfi(&io___66);
    do_fio(&c__1, (char *)&(*rvalue), (ftnlen)sizeof(doublereal));
    e_wsfi();
    lenbuf = i_len(buffer, buffer_len);
    s_copy(buff72, buffer, (ftnlen)72, buffer_len);
    i__1 = lenbuf;
    s_copy(buff72 + i__1, key, lenbuf + 16 - i__1, (ftnlen)16);
    rvalxx = *rvalue;
    s3opt_(&c_true, buff72, key, cvalue, &ivalue, &rvalxx, iprint, isumm, 
	    iexit, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, (ftnlen)72, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* npsetr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npSetr */
integer npget_(char *buffer, integer *iexit, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen buffer_len)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static char cw[8*500], key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalue;

/*     ================================================================== */
/*     npGet  decodes the option contained in  buffer */
/*     and returns 1 if the option has previously been set, else 0. */
/*     For example, */
/*     i = npGet ( 'Maximize', iExit, iw, leniw, rw, lenrw ) */

/*     01 Aug 2003: First version of npGet.  Needed because */
/*                  npGetc, npGeti, npGetr were not well defined */
/*                  for strings that had no numerical value. */
/*     18 Feb 2003: Current version of npGet. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ivalue, &rvalue, &c__0, &c__0, 
	    iexit, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    ret_val = ivalue;
    return ret_val;
} /* npget_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* integer function npGet */
/* Subroutine */ int npgeti_(char *buffer, integer *ivalue, integer *iexit, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	buffer_len)
{
    static char cw[8*500], key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static doublereal rvalue;

/*     ================================================================== */
/*     npGeti gets the value of the option contained in  buffer. */
/*     The parameters other than ivalue are as in npSet. */

/*     17 May 1998: first version of npGeti. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, ivalue, &rvalue, &c__0, &c__0, 
	    iexit, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* npgeti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine npGeti */
/* Subroutine */ int npgetr_(char *buffer, doublereal *rvalue, integer *iexit,
	 integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	buffer_len)
{
    static char cw[8*500], key[16];
    extern /* Subroutine */ int s3opt_(logical *, char *, char *, char *, 
	    integer *, doublereal *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char cvalue[8];
    static integer ivalue;

/*     ================================================================== */
/*     npGetr gets the value of the option contained in  buffer. */
/*     The parameters other than rvalue are as in npSet. */

/*     17 May 1998: first version of npGetr. */
/*     03 Nov 2000: current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    s3opt_(&c_false, buffer, key, cvalue, &ivalue, rvalue, &c__0, &c__0, 
	    iexit, cw, &c__500, &iw[1], leniw, &rw[1], lenrw, buffer_len, (
	    ftnlen)16, (ftnlen)8, (ftnlen)8);
    return 0;
} /* npgetr_ */

