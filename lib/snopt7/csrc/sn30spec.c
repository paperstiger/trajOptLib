/* ./src/sn30spec.f -- translated by f2c (version 20100827).
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

static integer c__72 = 72;
static integer c__1 = 1;
static integer c__14 = 14;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__13 = 13;
static integer c__11 = 11;
static integer c__12 = 12;
static logical c_true = TRUE_;
static integer c__3 = 3;
static integer c__89 = 89;
static integer c__10 = 10;
static integer c__70 = 70;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  sn30spec.f */

/*     s3optc   s3opti   s3optr   s3optl */
/*     s3file   s3key    s3opt    s3tie */
/*     s3unsetAll        s3unsetPrm */
/*     oplook   opnumb   opscan   optokn   opuppr */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s3optc_(logical *set, char *cwork, char *cvalue, ftnlen 
	cwork_len, ftnlen cvalue_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

/*     ================================================================== */
/*     s3optc  sets cwork to cvalue or vice versa depending on the value */
/*     of set. */

/*     17 May 1998: First version of s3optc. */
/*     17 May 1998: Current version. */
/*     ================================================================== */
    if (*set) {
	s_copy(cwork, cvalue, (ftnlen)8, (ftnlen)8);
    } else {
	s_copy(cvalue, cwork, (ftnlen)8, (ftnlen)8);
    }
    return 0;
} /* s3optc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3optc */
/* Subroutine */ int s3opti_(logical *set, integer *iwork, integer *ivalue)
{
/*     ================================================================== */
/*     s3opti  sets iwork to ivalue or vice versa depending on the value */
/*     of set. */

/*     17 May 1998: First version of s3opti. */
/*     17 May 1998: Current version. */
/*     ================================================================== */
    if (*set) {
	*iwork = *ivalue;
    } else {
	*ivalue = *iwork;
    }
    return 0;
} /* s3opti_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3opti */
/* Subroutine */ int s3optr_(logical *set, doublereal *rwork, doublereal *
	rvalue)
{
/*     ================================================================== */
/*     s3optr  sets rwork to rvalue or vice versa depending on the value */
/*     of set. */

/*     17 May 1998: First version of s3optr. */
/*     17 May 1998: Current version. */
/*     ================================================================== */
    if (*set) {
	*rwork = *rvalue;
    } else {
	*rvalue = *rwork;
    }
    return 0;
} /* s3optr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3optr */
/* Subroutine */ int s3optl_(logical *set, integer *iwork, integer *ivalue, 
	integer *l)
{
/*     ================================================================== */
/*     If set=true, s3optz sets iwork = ivalue and ignores l. */
/*     Otherwise  , s3optz sets l     = 1 if iwork = ivalue, else l = 0. */

/*     01 Aug 2003: First version of s3optl. */
/*     01 Aug 2003: Current version of s3optl. */
/*     ================================================================== */
    if (*set) {
	*iwork = *ivalue;
    } else {
	*l = 0;
	if (*ivalue == *iwork) {
	    *l = 1;
	}
    }
    return 0;
} /* s3optl_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3optl */
/* Subroutine */ int s3file_(integer *iexit, integer *ncalls, integer *ispecs,
	 S_fp opset, char *title, integer *iprint, integer *isumm, integer *
	errors, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen title_len, ftnlen cw_len)
{
    /* Initialized data */

    static char dashes[30] = "==============================";

    /* Format strings */
    static char fmt_2000[] = "(\002 XXX  Error while looking for a SPECS fil"
	    "e on unit\002,i6)";
    static char fmt_2200[] = "(\002 XXX  End-of-file encountered while looki"
	    "ng for\002,\002 a BEGIN file on unit\002,i6)";
    static char fmt_2300[] = "(\002 XXX  End-of-file encountered while proce"
	    "ssing\002,\002 a SPECS file on unit\002,i6)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    char ch__1[78], ch__2[39], ch__3[89], ch__4[31], ch__5[81], ch__6[73];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfi(icilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfi(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static char key[16], str[80];
    static integer nkey, nread;
    static char token[16*1];
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    s1trim_(char *, integer *, ftnlen);
    static char buffer[72];
    static integer lenbuf;
    extern /* Subroutine */ int snread_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static char cvalue[8];
    static integer ivalue;
    static doublereal rvalue;
    extern /* Subroutine */ int optokn_(char *, integer *, integer *, char *, 
	    ftnlen, ftnlen), snprnt_(integer *, char *, integer *, integer *, 
	    ftnlen);
    static integer endfile;

    /* Fortran I/O blocks */
    static icilist io___9 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___14 = { 0, str, 0, fmt_2200, 80, 1 };
    static icilist io___15 = { 0, str, 0, fmt_2300, 80, 1 };


/*     ================================================================== */
/*     s3file  reads the specs file from unit  iSpecs  and loads the */
/*     relevant optional parameters, using opset to process each line. */

/*     On exit, Errors says how many errors were encountered. */

/*     15 Nov 1991: First version based on Minos/Npsol routine s3file. */
/*     31 Jul 2003: snPRNT adopted.  iPrint, iSumm now used only by */
/*                  s3opt.  Beware -- they may get changed there. */
/*                  New input parameter "title" so we can remove s3fils. */
/*     27 Oct 2003: Errors counted separately from iExit. */
/*     06 Aug 2005: Fredrik Hellman found that str*72 is too short */
/*                  for format 2300.  Increased it to 80. */
/*     19 Dec 2005: iExit values clarified and reordered. */
/*     30 Apr 2006: snREAD adopted. */
/*     ================================================================== */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    *iexit = 0;
    *errors = 0;
/*     Return if the unit number is out of range. */
    if (*ispecs < 0 || *ispecs > 99) {
	*iexit = 131;
/* iSpecs out of range */
	return 0;
    }
/*     ------------------------------------------------------------------ */
/*     Look for  Begin, Endrun  or  Skip. */
/*     ------------------------------------------------------------------ */
    nread = 0;
L50:
    snread_(ispecs, buffer, &c__72, &endfile, (ftnlen)72);
/*        read (iSpecs, '(a72)', end = 920) buffer */
    if (endfile > 0) {
	goto L920;
    }
    ++nread;
    optokn_(buffer, &c__1, &nkey, token, (ftnlen)72, (ftnlen)16);
    s_copy(key, token, (ftnlen)16, (ftnlen)16);
    if (s_cmp(key, "ENDRUN", (ftnlen)16, (ftnlen)6) == 0) {
	goto L940;
    }
    if (s_cmp(key, "BEGIN", (ftnlen)16, (ftnlen)5) != 0) {
	if (nread == 1 && s_cmp(key, "SKIP", (ftnlen)16, (ftnlen)4) != 0) {
	    ++(*errors);
	    s_wsfi(&io___9);
	    do_fio(&c__1, (char *)&(*ispecs), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__14, str, &iw[1], leniw, (ftnlen)80);
	    snprnt_(&c__4, " XXX  The file should start with Begin, Skip or "
		    "Endrun", &iw[1], leniw, (ftnlen)54);
	    snprnt_(&c__4, " XXX  but the first record found was the followi"
		    "ng:", &iw[1], leniw, (ftnlen)51);
/* Writing concatenation */
	    i__1[0] = 6, a__1[0] = " ---->";
	    i__1[1] = 72, a__1[1] = buffer;
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)78);
	    snprnt_(&c__14, ch__1, &iw[1], leniw, (ftnlen)78);
	    snprnt_(&c__14, " XXX  Continuing to look for SPECS file...", &iw[
		    1], leniw, (ftnlen)42);
	}
	goto L50;
    }
/*     ------------------------------------------------------------------ */
/*     Begin found. */
/*     This is taken to be the first line of a SPECS file. */
/*     Print the title first if it's not blank. */
/*     ------------------------------------------------------------------ */
    s1page_(&c__1, &iw[1], leniw);
    s1trim_(buffer, &lenbuf, (ftnlen)72);
    if (s_cmp(title, " ", title_len, (ftnlen)1) != 0) {
	s_copy(str, title, (ftnlen)80, title_len);
	snprnt_(&c__13, " ", &iw[1], leniw, (ftnlen)1);
/* Writing concatenation */
	i__1[0] = 9, a__1[0] = "         ";
	i__1[1] = 30, a__1[1] = dashes;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)39);
	snprnt_(&c__1, ch__2, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
	i__1[0] = 9, a__1[0] = "         ";
	i__1[1] = 80, a__1[1] = str;
	s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)89);
	snprnt_(&c__1, ch__3, &iw[1], leniw, (ftnlen)89);
/* Writing concatenation */
	i__1[0] = 9, a__1[0] = "         ";
	i__1[1] = 30, a__1[1] = dashes;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)39);
	snprnt_(&c__1, ch__2, &iw[1], leniw, (ftnlen)39);
/* Writing concatenation */
	i__1[0] = 1, a__1[0] = " ";
	i__1[1] = 30, a__1[1] = dashes;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)31);
	snprnt_(&c__2, ch__4, &iw[1], leniw, (ftnlen)31);
/* Writing concatenation */
	i__1[0] = 1, a__1[0] = " ";
	i__1[1] = 80, a__1[1] = str;
	s_cat(ch__5, a__1, i__1, &c__2, (ftnlen)81);
	snprnt_(&c__2, ch__5, &iw[1], leniw, (ftnlen)81);
/* Writing concatenation */
	i__1[0] = 1, a__1[0] = " ";
	i__1[1] = 30, a__1[1] = dashes;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)31);
	snprnt_(&c__2, ch__4, &iw[1], leniw, (ftnlen)31);
    }
    snprnt_(&c__11, " SPECS file", &iw[1], leniw, (ftnlen)11);
    snprnt_(&c__1, " ----------", &iw[1], leniw, (ftnlen)11);
/* Writing concatenation */
    i__1[0] = 6, a__1[0] = "      ";
    i__1[1] = lenbuf, a__1[1] = buffer;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)78);
    snprnt_(&c__11, ch__1, &iw[1], leniw, lenbuf + 6);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = lenbuf, a__1[1] = buffer;
    s_cat(ch__6, a__1, i__1, &c__2, (ftnlen)73);
    snprnt_(&c__12, ch__6, &iw[1], leniw, lenbuf + 1);
/*     ------------------------------------------------------------------ */
/*     Read the rest of the file. */
/*     ------------------------------------------------------------------ */
/* +    while (key .ne. 'END') do */
L100:
    if (s_cmp(key, "END", (ftnlen)16, (ftnlen)3) != 0) {
	snread_(ispecs, buffer, &c__72, &endfile, (ftnlen)72);
	if (endfile > 0) {
	    goto L930;
	}
	(*opset)(&c_true, buffer, key, cvalue, &ivalue, &rvalue, iprint, 
		isumm, errors, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)72, (ftnlen)16, (ftnlen)8, (ftnlen)8);
	goto L100;
    }
/* +    end do */
    return 0;
L920:
    if (*ncalls <= 1) {
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*ispecs), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__14, str, &iw[1], leniw, (ftnlen)80);
    } else {
	snprnt_(&c__3, " Endrun", &iw[1], leniw, (ftnlen)7);
    }
    *iexit = 132;
/* End-of-file found while looking for BEGIN */
    return 0;
L930:
    s_wsfi(&io___15);
    do_fio(&c__1, (char *)&(*ispecs), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__14, str, &iw[1], leniw, (ftnlen)80);
    *iexit = 133;
/* End-of-file found while reading specs */
    return 0;
L940:
/* Writing concatenation */
    i__1[0] = 6, a__1[0] = "      ";
    i__1[1] = lenbuf, a__1[1] = buffer;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)78);
    snprnt_(&c__11, ch__1, &iw[1], leniw, lenbuf + 6);
/* Writing concatenation */
    i__1[0] = 1, a__1[0] = " ";
    i__1[1] = lenbuf, a__1[1] = buffer;
    s_cat(ch__6, a__1, i__1, &c__2, (ftnlen)73);
    snprnt_(&c__12, ch__6, &iw[1], leniw, lenbuf + 1);
    *iexit = 134;
/* Endrun found for empty SPECS file */
    return 0;
} /* s3file_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3file */
/* Subroutine */ int s3key_(char *key, integer *loc, ftnlen key_len)
{
    /* Initialized data */

    static char keys[16*89] = "AIJ             " "BACKUP          " "BOUNDS "
	    "         " "CALL            " "CENTRAL         " "CG            "
	    "  " "CHECK           " "COEFFICIENTS    " "COLD            " 
	    "COLUMNS         " "CRASH           " "CYCLE           " "DEBUG "
	    "          " "DEFAULTS        " "DERIVATIVE      " "DIFFERENCE   "
	    "   " "DUMP            " "ELASTIC         " "ELEMENTS        " 
	    "ERROR           " "EXPAND          " "FACTORIZATION   " "FEASIB"
	    "ILITY     " "FEASIBLE        " "FUNCTION        " "HESSIAN      "
	    "   " "HOT             " "INFEASIBLE      " "INFINITE        " 
	    "INSERT          " "ITERATIONS      " "IW              " "JACOBI"
	    "AN        " "LINESEARCH      " "LIST            " "LOAD         "
	    "   " "LOG             " "LOWER           " "LP              " 
	    "LU              " "MAJOR           " "MAXIMIZE        " "MINIMI"
	    "ZE        " "MINOR           " "MPS             " "NEW          "
	    "   " "NO              " "NON             " "NONDERIVATIVE   " 
	    "NONLINEAR       " "OBJECTIVE       " "OLD             " "OPTIMA"
	    "LITY      " "PARTIAL         " "PENALTY         " "PIVOT        "
	    "   " "PRINT           " "PROBLEM         " "PROXIMAL        " 
	    "PUNCH           " "QP              " "QPSOLVER        " "REDUCE"
	    "D         " "RANGES          " "REPORT          " "RHS          "
	    "   " "ROWS            " "RW              " "SAVE            " 
	    "SCALE           " "SOLUTION        " "START           " "STICKY"
	    "          " "STOP            " "SUBSPACE        " "SUMMARY      "
	    "   " "SUPERBASICS     " "SUPPRESS        " "SYSTEM          " 
	    "TIMING          " "TOTAL           " "UNBOUNDED       " "UPPER "
	    "          " "USER            " "VERIFY          " "VIOLATION    "
	    "   " "WARM            " "WORKING         " "WORKSPACE       ";

    extern /* Subroutine */ int oplook_(integer *, char *, logical *, char *, 
	    integer *, ftnlen, ftnlen);

/*     ================================================================== */
/*     s3key  sets key to be the standard form for the first keyword */
/*     on each line of a SPECS file. */

/*     17 May 1998: First version. */
/*     13 Mar 2004: Hot start option added. */
/*     22 Jun 2004: System info option added. */
/*     01 Sep 2007: Sticky parameters added. */
/*     19 Jun 2008: Call status added. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    oplook_(&c__89, keys, &c_true, key, loc, (ftnlen)16, (ftnlen)16);
    return 0;
} /* s3key_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3key */
/* Subroutine */ int s3opt_(logical *s, char *buffer, char *key, char *c__, 
	integer *i__, doublereal *r__, integer *lprnt, integer *lsumm, 
	integer *errors, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen buffer_len, ftnlen 
	key_len, ftnlen c_len, ftnlen cw_len)
{
    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsfi(icilist *), e_rsfi(
	    void);

    /* Local variables */
    static integer i0, i1, i2, i3, m1;
    extern /* Subroutine */ int s3unsetprm_(char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static char str[132];
    static integer loc1, loc2;
    static char key2[16], key3[16], str1[132];
    static integer lenb;
    static logical more;
    extern /* Subroutine */ int s3tie_(char *, integer *, ftnlen), s3key_(
	    char *, integer *, ftnlen);
    static char value[16], token[16*10];
    extern /* Subroutine */ int s3optc_(logical *, char *, char *, ftnlen, 
	    ftnlen), s1trim_(char *, integer *, ftnlen), s3opti_(logical *, 
	    integer *, integer *), s3optl_(logical *, integer *, integer *, 
	    integer *), s3optr_(logical *, doublereal *, doublereal *);
    static integer lenbuf, maxbuf;
    static logical number;
    static integer ntoken;
    extern logical opnumb_(char *, ftnlen);
    static integer maxint;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), optokn_(char *, integer *, integer *, char *, ftnlen, 
	    ftnlen);
    extern integer s1intmx_(void);

    /* Fortran I/O blocks */
    static icilist io___21 = { 0, str, 0, "(6x,a)", 132, 1 };
    static icilist io___31 = { 0, value, 0, "(bn, e16.0)", 16, 1 };
    static icilist io___38 = { 0, key2, 0, "(bn, i16)", 16, 1 };
    static icilist io___39 = { 0, key2, 0, "(bn, e16.0)", 16, 1 };
    static icilist io___40 = { 0, str, 0, "(2a)", 132, 1 };
    static icilist io___41 = { 0, str, 0, "(2a)", 132, 1 };
    static icilist io___42 = { 0, str, 0, "(2a)", 132, 1 };
    static icilist io___43 = { 0, str, 0, "(a,i8)", 132, 1 };
    static icilist io___45 = { 0, str1, 0, "(1x,a)", 132, 1 };


/*     ================================================================== */
/*     s3opt  decodes the option contained in  buffer  in order to */
/*     set or get a parameter value in the relevant array iw or rw. */

/*     The buffer is output to file iPrint, minus trailing blanks. */
/*     Error messages are output to files iPrint and iSumm. */
/*     buffer is echoed to iPrint but normally not to iSumm. */
/*     It is echoed to iSumm before any error msg. */

/*     On entry, */
/*     buffer contains the option string to be processed. */
/*     s      is true if an option is to be extracted from buffer. */
/*            Otherwise, c, i and r are to be assigned the value of the */
/*            option defined in the option string. */
/*     lPrnt  is iPrint as given to s3file. */
/*     lSumm  is iSumm  as given to s3file. */
/*     Errors is the number of errors so far. */

/*     On exit, */
/*     key    is the first keyword contained in buffer. */
/*     If s is true, c, i and r may be ignored.  (They are usually */
/*            option values that have been saved in cw, iw, rw.) */
/*     If s is false, */
/*     c      is the OBJECTIVE, RHS, RANGE or BOUND name if key is */
/*            one of those words. */
/*     r      is the first numerical value found in buffer (or zero). */
/*     i      is int(r) if  abs(r) < maxint. */
/*     Errors is the number of errors so far. */


/*     s3opt  uses opnumb and the subprograms */
/*                 lookup, scannr, tokens, upcase */
/*     (now called oplook, opscan, optokn, opuppr) */
/*     supplied by Sterling Software, Palo Alto, California. */

/*     15 Nov 1991: First version based on s3key/opkey. */
/*     10 Dec 2002: Recognize LU Diagonal and LU Rook Pivoting */
/*     02 Jul 2003: Options for QN QP solver added. */
/*     01 Aug 2003: snPRNT adopted.  If lPrnt, lSumm are positive, */
/*                  snPRNT outputs to global iPrint, iSumm, iStdo. */
/*                  Sometimes lPrnt, lSumm are zero (for "get" routines). */
/*                  snPRNT( 4, ... ) then outputs error msgs */
/*     21 Dec 2003: Hot start added. */
/*     22 Jun 2004: lvlSys added */
/*     04 Dec 2004: kosher maxint added */
/*     22 Apr 2007: Allowed pure comments to start in any column */
/*     16 Jun 2007: "iw 20" or "rw 20" now grabs iw(20) or rw(20) */
/*     01 Sep 2007: Sticky parameters implemented */
/*     25 Nov 2007: Hessian parameters implemented */
/*     18 Jun 2008: Call status option added */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Minor Phase 1 Opt tol */
/* Minor Phase 2 Opt tol */
/* Major Optimality tolerance */
/* cg tolerance */
/* Minor feasibility tolerance */
/* Major feasibility tolerance */
/* excludes small pivot elems */
/* crash tolerance */
/* LU swap tolerance */
/* LU factor tolerance */
/* LU update tolerance */
/* definition of plus infinity */
/* unbounded objective */
/* unbounded step */
/* relative function precision */
/* forward difference interval */
/* central difference interval */
/* Step limit */
/* violation limit */
/* Quasi-Newton QP rg tolerance */
/* line search tolerance */
/* infeasibility weight */
/* initial penalty parameter */
/* Elastic weightmax */
/* scale tolerance */
/* zero Aij tolerance */
/* default lower bound on x */
/* default upper bound on x */
/* abs tol for small diag of U */
/* rel tol for small diag of U */
/* switch to dense LU */
/* Start of SNOPT part of rw */
/* End   of SNOPT part of rw */
/* Start of SNOPT part of iw */
/* End   of SNOPT part of iw */
/* Start of SNOPT part of cw */
/* End   of SNOPT part of cw */
/* Print   file */
/* Summary file */
/* # nonlinear Jacobian variabl */
/* # variables in gObj */
/* nonlinear constraints */
/*   max( nnObj, nnJac ) */
/* max columns of R */
/* max # of superbasics */
/* (ge 0) max # of BFGS updates */
/* 0(1) => QP(QN) QP solver */
/* >0    => use elastic mode */
/* check (row) frequency */
/* factorization frequency */
/* save basis map */
/* log/print frequency */
/* Summary print frequency */
/* max. expansions of featol */
/* Hessian frequency */
/* Hessian flush */
/* = 0(1) => cold(warm) start */
/* derivative level */
/* > 0   => print system info */
/* 0,1,2 => LM, FM, Exact H */
/* Elastic option */
/* scale option */
/* >0 => deriv. line search */
/* >0 => QN preconditioned CG */
/* Verify level */
/* Proximal Point method for x0 */
/* 0(1 2 3) LU pivoting */
/* > 0  =>  parms are printed */
/* line search debug start itn */
/* > 0  => print the scales */
/* > 0  =>  print the solution */
/* > 0  => private debug print */
/* 1, -1  => MIN, FP, MAX */
/* Crash option */
/* limit on total iterations */
/* limit on major iterations */
/* limit on minor iterations */
/* Major print level */
/* Minor print level */
/* # partial pricing sections */
/* maximum # of new SB */
/* CG iteration limit */
/* Objective row */
/* 0, 1, 2 => derivative option */
/* 1(2) => dense(sparse) deriv. */
/* maximum # errors in MPS data */
/* maximum # lines  of MPS data */
/* problem number */
/* start g derivative checking */
/* stop  g derivative checking */
/* start J derivative checking */
/* stop  J derivative checking */
/* start H derivative checking */
/* stop  H derivative checking */
/* > 0 => sticky parameters */
/* backup file */
/* dump file */
/* load file */
/* MPS file */
/* new basis file */
/* insert file */
/* old basis file */
/* punch file */
/* Report file */
/* Solution file */
/* Row    estimate */
/* Column estimate */
/* Estimated element count */
/* Timing level */
/* QP user-routine call-status */
/* NP user-routine call-status */
/* Objective name */
/* Right-hand side name */
/* Range name */
/*     ------------------------------------------------------------------ */
/* Bnd section name */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    maxint = s1intmx_();
/*     ------------------------------------------------------------------ */
/*     Trim trailing blanks and echo to the Print file. */
/*     ------------------------------------------------------------------ */
    s1trim_(buffer, &lenbuf, buffer_len);
    maxbuf = min(120,lenbuf);
    if (*lprnt > 0) {
	s_wsfi(&io___21);
	do_fio(&c__1, buffer, maxbuf);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
    }
/*     Set lenb = length of buffer without trailing comments. */
/*     Eliminate comments and empty lines. */
/*     A '*' appearing anywhere in buffer terminates the string. */
    *i__ = i_indx(buffer, "*", lenbuf, (ftnlen)1);
    if (*i__ == 0) {
	lenb = lenbuf;
    } else {
	lenb = *i__ - 1;
    }
    if (lenb <= 0) {
	s_copy(key, "*", key_len, (ftnlen)1);
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Extract up to maxtok tokens from the record. */
/*     ntoken returns how many were actually found. */
/*     key, key2, are the first tokens if any, otherwise blank. */
/*     For some values of key (bounds, objective, ranges, rhs) */
/*     we have to save key2 before s3tie (and oplook) alter it. */
/*     For example, if the data is     objective = obj */
/*     oplook will change obj to objective. */
/*     ------------------------------------------------------------------ */
    optokn_(buffer, &c__10, &ntoken, token, lenbuf, (ftnlen)16);
    s_copy(key, token, key_len, (ftnlen)16);
    s_copy(key2, token + 16, (ftnlen)16, (ftnlen)16);
    s_copy(key3, token + 32, (ftnlen)16, (ftnlen)16);
    s_copy(c__, key2, (ftnlen)8, (ftnlen)8);
/*     Certain keywords require no action. */
    if (s_cmp(key, "   ", key_len, (ftnlen)3) == 0) {
	goto L900;
    }
/* blank line */
    if (s_cmp(key, "*  ", key_len, (ftnlen)3) == 0) {
	goto L900;
    }
/* comment starting in column no. > 1 */
    if (s_cmp(key, "END", key_len, (ftnlen)3) == 0) {
	goto L900;
    }
/*     Convert the keywords to their most fundamental form */
/*     (upper case, no abbreviations). */
/*     loci   says where the keywords are in the dictionaries. */
/*     loci = 0 signals that the keyword wasn't there. */
    s3key_(key, &loc1, key_len);
    s3tie_(key2, &loc2, (ftnlen)16);
/*     Most keywords will have an associated integer or real value, */
/*     so look for it no matter what the keyword. */
    s_copy(c__, key2, (ftnlen)8, (ftnlen)8);
    *i__ = 1;
    number = FALSE_;
/* +    while (i .lt. ntoken  .and.  .not. number) loop */
L50:
    if (*i__ < ntoken && ! number) {
	++(*i__);
	s_copy(value, token + (*i__ - 1 << 4), (ftnlen)16, (ftnlen)16);
	number = opnumb_(value, (ftnlen)16);
	goto L50;
    }
/* +    end while */
    *i__ = 0;
    *r__ = 0.;
    if (number) {
	s_rsfi(&io___31);
	do_fio(&c__1, (char *)&(*r__), (ftnlen)sizeof(doublereal));
	e_rsfi();
	*i__ = maxint;
	if (abs(*r__) < (doublereal) maxint) {
	    *i__ = (integer) (*r__);
	}
    }
/*     ------------------------------------------------------------------ */
/*     Decide what to do about each keyword. */
/*     The second keyword (if any) might be needed to break ties. */
/*     Some seemingly redundant testing of more is used */
/*     to avoid compiler limits on the number of consecutive else ifs. */
/*     ------------------------------------------------------------------ */
    m1 = -1;
    i0 = 0;
    i1 = 1;
    i2 = 2;
    i3 = 3;
    more = TRUE_;
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "BACKUP      ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[120], i__);
	} else if (s_cmp(key, "CALL        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[235], i__);
	    s3opti_(s, &iw[236], i__);
	} else if (s_cmp(key, "CENTRAL     ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[77], r__);
	} else if (s_cmp(key, "CG          ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "TOLERANCE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[54], r__);
	    }
	    if (s_cmp(key2, "PRECONDITIONING", (ftnlen)16, (ftnlen)15) == 0) {
		s3opti_(s, &iw[77], i__);
	    }
	    if (s_cmp(key2, "ITERATIONS  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[97], i__);
	    }
	} else if (s_cmp(key, "CHECK       ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[58], i__);
	} else if (s_cmp(key, "COLD        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[69], &i0);
	} else if (s_cmp(key, "CRASH       ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "OPTION      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[88], i__);
	    }
	    if (s_cmp(key2, "TOLERANCE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[62], r__);
	    }
	} else if (s_cmp(key, "DEBUG       ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[85], i__);
	} else if (s_cmp(key, "DEFAULTS    ", key_len, (ftnlen)12) == 0) {
	    s3unsetprm_(cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)
		    8);
	} else if (s_cmp(key, "DERIVATIVE  ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "LEVEL       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[70], i__);
	    }
	    if (s_cmp(key2, "LINESEARCH  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[76], &i1, i__);
	    }
	    if (s_cmp(key2, "OPTION      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[104], i__);
	    }
	} else if (s_cmp(key, "DIFFERENCE  ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[76], r__);
	} else if (s_cmp(key, "DUMP        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[121], i__);
	} else if (s_cmp(key, "ELASTIC     ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "OBJECTIVE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[73], i__);
	    }
	    if (s_cmp(key2, "MODE        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[56], i__);
	    }
	    if (s_cmp(key2, "WEIGHT      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[88], r__);
	    }
	    if (s_cmp(key2, "WEIGHTMAX   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[90], r__);
	    }
	} else if (s_cmp(key, "EXPAND      ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[63], i__);
	} else if (s_cmp(key, "FACTORIZATION", key_len, (ftnlen)13) == 0) {
	    s3opti_(s, &iw[59], i__);
	} else if (s_cmp(key, "FEASIBLE    ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "POINT       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[87], &i0, i__);
	    }
	    if (s_cmp(key2, "EXIT        ", (ftnlen)16, (ftnlen)12) == 0) {
		goto L890;
	    }
	} else if (s_cmp(key, "FEASIBILITY ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "TOLERANCE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[56], r__);
	    }
	} else if (s_cmp(key, "FUNCTION    ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[73], r__);
	} else if (s_cmp(key, "HESSIAN     ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "COLUMNS     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[24], i__);
	    }
	    if (s_cmp(key2, "DIMENSION   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[52], i__);
	    }
	    if (s_cmp(key2, "FREQUENCY   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[64], i__);
	    }
	    if (s_cmp(key2, "FLUSH       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[66], i__);
	    }
	    if (s_cmp(key2, "UPDATES     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[54], i__);
	    }
	    if (s_cmp(key2, "LIMITED     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[72], &i0, i__);
	    }
	    if (s_cmp(key2, "FULL        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[72], &i1, i__);
	    }
	    if (s_cmp(key2, "PRECONDITIONING", (ftnlen)16, (ftnlen)15) == 0) {
		s3opti_(s, &iw[77], i__);
	    }
	} else if (s_cmp(key, "HOT         ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[69], &i3);
	} else if (s_cmp(key, "INFINITE    ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[70], r__);
	} else if (s_cmp(key, "INSERT      ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[125], i__);
	} else if (s_cmp(key, "ITERATIONS  ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[89], i__);
	} else {
	    more = TRUE_;
	}
    }
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "IW          ", key_len, (ftnlen)12) == 0) {
	    if (*i__ < 1 || *i__ > 500) {
		goto L880;
	    } else if (*s) {
/*              Allow things like  iw 21 = 100  to set iw(21) = 100 */
		s_copy(key2, token + 32, (ftnlen)16, (ftnlen)16);
		s_rsfi(&io___38);
		do_fio(&c__1, (char *)&iw[*i__], (ftnlen)sizeof(integer));
		e_rsfi();
	    } else {
/*              Grab the contents of  iw(i) */
		*i__ = iw[*i__];
		*r__ = (doublereal) (*i__);
	    }
	} else if (s_cmp(key, "LINESEARCH  ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "TOLERANCE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[84], r__);
	    }
	    if (s_cmp(key2, "DEBUG       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[82], i__);
	    }
	} else if (s_cmp(key, "LOAD        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[122], i__);
	} else if (s_cmp(key, "LOG         ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[61], i__);
	} else if (s_cmp(key, "LP          ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "FEASIBILITY ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[56], r__);
	    }
	    if (s_cmp(key2, "OPTIMALITY  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[52], r__);
	    }
	} else if (s_cmp(key, "LU          ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "PARTIAL     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[80], &i0, i__);
	    }
	    if (s_cmp(key2, "COMPLETE    ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[80], &i2, i__);
	    }
	    if (s_cmp(key2, "DIAGONAL    ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[80], &i3, i__);
	    }
	    if (s_cmp(key2, "FACTORIZATION", (ftnlen)16, (ftnlen)13) == 0) {
		s3optr_(s, &rw[66], r__);
	    }
	    if (s_cmp(key2, "ROOK        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[80], &i1, i__);
	    }
	    if (s_cmp(key2, "UPDATES     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[67], r__);
	    }
	    if (s_cmp(key2, "DENSITY     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[158], r__);
	    }
	    if (s_cmp(key2, "SINGULARITY ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[154], r__);
		s3optr_(s, &rw[155], r__);
	    }
	    if (s_cmp(key2, "SWAP        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[65], r__);
	    }
/*              if (key2.eq. 'DEFAULTS    ') then */
/*                 if (loc3.eq.  0           ) go to 820 */
/*                 if (key3.eq.'TPP         ') call s3optr(s,rw(tolDpp),r) */
/*                 if (key3.eq.'TCP         ') call s3optr(s,rw(tolDcp),r) */
/*                 if (key3.eq.'UPDATES     ') call s3optr(s,rw(tolDup),r) */
/*              end if */
	} else {
	    more = TRUE_;
	}
    }
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "MAJOR       ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "FEASIBILITY ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[57], r__);
	    }
	    if (s_cmp(key2, "ITERATIONS  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[90], i__);
	    }
	    if (s_cmp(key2, "OPTIMALITY  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[53], r__);
	    }
	    if (s_cmp(key2, "PRINT       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[92], i__);
	    }
	    if (s_cmp(key2, "STEP        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[80], r__);
	    }
	} else if (s_cmp(key, "MAXIMIZE    ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[87], &m1);
	} else if (s_cmp(key, "MINIMIZE    ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[87], &i1);
	} else if (s_cmp(key, "MINOR       ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "ITERATIONS  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[91], i__);
	    }
	    if (s_cmp(key2, "FEASIBILITY ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[56], r__);
	    }
	    if (s_cmp(key2, "OPTIMALITY  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[52], r__);
	    }
	    if (s_cmp(key2, "PHASE1      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[51], r__);
	    }
	    if (s_cmp(key2, "PHASE2      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[52], r__);
	    }
	    if (s_cmp(key2, "PRINT       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[93], i__);
	    }
	    if (s_cmp(key2, "SUPERBASICS ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[95], i__);
	    }
	} else if (s_cmp(key, "NEW         ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "BASIS       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[124], i__);
	    }
	    if (s_cmp(key2, "SUPERBASICS ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[95], i__);
	    }
	} else if (s_cmp(key, "NO          ", key_len, (ftnlen)12) == 0 || 
		s_cmp(key, "NONDERIVATIVE", key_len, (ftnlen)13) == 0 || 
		s_cmp(key, "NON         ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[76], &i0);
	} else if (s_cmp(key, "OBJECTIVE   ", key_len, (ftnlen)12) == 0) {
	    if (s_cmp(key2, "ROW         ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[103], i__);
	    } else {
		s3optc_(s, cw + 416, c__, (ftnlen)8, (ftnlen)8);
	    }
	} else if (s_cmp(key, "OLD         ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[126], i__);
	} else if (s_cmp(key, "OPTIMALITY  ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[53], r__);
	} else if (s_cmp(key, "PARTIAL     ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[94], i__);
	} else if (s_cmp(key, "PENALTY     ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[89], r__);
	} else if (s_cmp(key, "PIVOT       ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[60], r__);
	} else if (s_cmp(key, "PRINT       ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "FILE        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[12], i__);
	    }
	    if (s_cmp(key2, "FREQUENCY   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[61], i__);
	    }
	    if (s_cmp(key2, "LEVEL       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[93], i__);
	    }
	} else if (s_cmp(key, "PROXIMAL    ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[79], i__);
	} else if (s_cmp(key, "PUNCH       ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[127], i__);
	} else {
	    more = TRUE_;
	}
    }
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "QP          ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "COLUMNS     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[24], i__);
	    }
	    if (s_cmp(key2, "FEASIBILITY ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[56], r__);
	    }
	    if (s_cmp(key2, "OPTIMALITY  ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[52], r__);
	    }
	} else if (s_cmp(key, "QPSOLVER    ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "CHOLESKY    ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[55], &i0, i__);
	    }
	    if (s_cmp(key2, "CG          ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[55], &i1, i__);
	    }
	    if (s_cmp(key2, "QN          ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[55], &i2, i__);
	    }
	} else if (s_cmp(key, "REDUCED     ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[52], i__);
	} else if (s_cmp(key, "REPORT      ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[130], i__);
	} else if (s_cmp(key, "ROWS        ", key_len, (ftnlen)12) == 0) {
/*             gams should recognize row tolerance */
/*             but not just          rows */
/*             This is a relic from MINOS */
	    if (s_cmp(key2, "TOLERANCE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[57], r__);
	    } else {
		s3opti_(s, &iw[133], i__);
	    }
	} else if (s_cmp(key, "RW          ", key_len, (ftnlen)12) == 0) {
	    if (*i__ < 1 || *i__ > 500) {
		goto L880;
	    } else if (*s) {
/*              Allow things like rw 21 = 2  to set rw(21) = 2.0 */
		s_copy(key2, token + 32, (ftnlen)16, (ftnlen)16);
		s_rsfi(&io___39);
		do_fio(&c__1, (char *)&rw[*i__], (ftnlen)sizeof(doublereal));
		e_rsfi();
	    } else {
/*              Grab the contents of  rw(i) */
		*r__ = rw[*i__];
	    }
	} else if (s_cmp(key, "SAVE        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[60], i__);
	} else if (s_cmp(key, "SCALE       ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "OPTION      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[75], i__);
	    }
	    if (s_cmp(key2, "TOLERANCE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[92], r__);
	    }
	    if (s_cmp(key2, "PRINT       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[83], &i1, i__);
	    }
	} else {
	    more = TRUE_;
	}
    }
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "SOLUTION    ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "FILE        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[131], i__);
	    }
	    if (s_cmp(key2, "YES         ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[84], &i2, i__);
	    }
	    if (s_cmp(key2, "NO          ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[84], &i0, i__);
	    }
	} else if (s_cmp(key, "START       ", key_len, (ftnlen)12) == 0) {
	    if (s_cmp(key2, "OBJECTIVE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[110], i__);
	    }
	    if (s_cmp(key2, "CONSTRAINTS ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[112], i__);
	    }
	    if (s_cmp(key2, "HESSIAN     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[114], i__);
	    }
	    if (loc2 == 0) {
		goto L820;
	    }
	} else if (s_cmp(key, "STOP        ", key_len, (ftnlen)12) == 0) {
	    if (s_cmp(key2, "OBJECTIVE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[111], i__);
	    }
	    if (s_cmp(key2, "CONSTRAINTS ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[113], i__);
	    }
	    if (s_cmp(key2, "HESSIAN     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[115], i__);
	    }
	    if (loc2 == 0) {
		goto L820;
	    }
	} else {
	    more = TRUE_;
	}
    }
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "STICKY      ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key3, "YES         ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[116], &i1, i__);
	    }
	    if (s_cmp(key3, "NO          ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[116], &i0, i__);
	    }
	} else if (s_cmp(key, "SUBSPACE    ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[83], r__);
	} else if (s_cmp(key, "SUPERBASICS ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[53], i__);
	} else if (s_cmp(key, "SUMMARY     ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "FILE        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[13], i__);
	    }
	    if (s_cmp(key2, "FREQUENCY   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[62], i__);
	    }
	} else if (s_cmp(key, "SUPPRESS    ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[81], i__);
	} else if (s_cmp(key, "TIMING      ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[182], i__);
	} else if (s_cmp(key, "SYSTEM      ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key3, "YES         ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[71], &i1, i__);
	    }
	    if (s_cmp(key3, "NO          ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[71], &i0, i__);
	    }
	} else if (s_cmp(key, "TOTAL       ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "INTEGER     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[5], i__);
	    }
	    if (s_cmp(key2, "REAL        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[3], i__);
	    }
	    if (s_cmp(key2, "CHARACTER   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[7], i__);
	    }
	} else if (s_cmp(key, "UNBOUNDED   ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "OBJECTIVE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[71], r__);
	    }
	    if (s_cmp(key2, "STEP        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optr_(s, &rw[72], r__);
	    }
	} else if (s_cmp(key, "USER        ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "INTEGER     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[4], i__);
	    }
	    if (s_cmp(key2, "REAL        ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[2], i__);
	    }
	    if (s_cmp(key2, "CHARACTER   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[6], i__);
	    }
	} else if (s_cmp(key, "VERIFY      ", key_len, (ftnlen)12) == 0) {
	    if (s_cmp(key2, "            ", (ftnlen)16, (ftnlen)12) == 0) {
		loc2 = 1;
		*i__ = 3;
	    }
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "OBJECTIVE   ", (ftnlen)16, (ftnlen)12) == 0) {
		*i__ = 1;
	    }
	    if (s_cmp(key2, "CONSTRAINTS ", (ftnlen)16, (ftnlen)12) == 0) {
		*i__ = 2;
	    }
	    if (s_cmp(key2, "GRADIENTS   ", (ftnlen)16, (ftnlen)12) == 0) {
		*i__ = 3;
	    }
	    if (s_cmp(key2, "YES         ", (ftnlen)16, (ftnlen)12) == 0) {
		*i__ = 3;
	    }
	    if (s_cmp(key2, "NO          ", (ftnlen)16, (ftnlen)12) == 0) {
		*i__ = 0;
	    }
	    if (s_cmp(key2, "LEVEL       ", (ftnlen)16, (ftnlen)12) == 0) {
		*i__ = *i__;
	    }
	    s3opti_(s, &iw[78], i__);
	} else if (s_cmp(key, "VIOLATION   ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[81], r__);
	} else if (s_cmp(key, "WARM        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[69], &i2);
	} else if (s_cmp(key, "WORKSPACE   ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "(USER)      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[2], i__);
	    }
	    if (s_cmp(key2, "(TOTAL)     ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[3], i__);
	    }
	} else {
	    more = TRUE_;
	}
    }
    if (! more) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Keywords for MPS files. */
/*     ------------------------------------------------------------------ */
    if (more) {
	more = FALSE_;
	if (s_cmp(key, "AIJ         ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[95], r__);
	} else if (s_cmp(key, "BOUNDS      ", key_len, (ftnlen)12) == 0) {
	    s3optc_(s, cw + 440, c__, (ftnlen)8, (ftnlen)8);
	} else if (s_cmp(key, "COEFFICIENTS", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[135], i__);
	} else if (s_cmp(key, "COLUMNS     ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[134], i__);
	} else if (s_cmp(key, "ELEMENTS    ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[135], i__);
	} else if (s_cmp(key, "ERROR       ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[106], i__);
	} else if (s_cmp(key, "INFINITE    ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[70], r__);
	} else if (s_cmp(key, "JACOBIAN    ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "DENSE       ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[105], &i1, i__);
	    }
	    if (s_cmp(key2, "SPARSE      ", (ftnlen)16, (ftnlen)12) == 0) {
		s3optl_(s, &iw[105], &i2, i__);
	    }
	} else if (s_cmp(key, "LIST        ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[107], i__);
	} else if (s_cmp(key, "LOWER       ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[96], r__);
	} else if (s_cmp(key, "MPS         ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[123], i__);
	} else if (s_cmp(key, "NONLINEAR   ", key_len, (ftnlen)12) == 0) {
	    if (loc2 == 0) {
		goto L820;
	    }
	    if (s_cmp(key2, "CONSTRAINTS ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[23], i__);
	    }
	    if (s_cmp(key2, "OBJECTIVE   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[22], i__);
	    }
	    if (s_cmp(key2, "JACOBIAN    ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[21], i__);
	    }
	    if (s_cmp(key2, "VARIABLES   ", (ftnlen)16, (ftnlen)12) == 0) {
		s3opti_(s, &iw[22], i__);
		s3opti_(s, &iw[21], i__);
	    }
	} else if (s_cmp(key, "OBJECTIVE   ", key_len, (ftnlen)12) == 0) {
	    s3optc_(s, cw + 416, c__, (ftnlen)8, (ftnlen)8);
	} else if (s_cmp(key, "PROBLEM     ", key_len, (ftnlen)12) == 0) {
	    s3opti_(s, &iw[108], i__);
	} else if (s_cmp(key, "RANGES      ", key_len, (ftnlen)12) == 0) {
	    s3optc_(s, cw + 432, c__, (ftnlen)8, (ftnlen)8);
	} else if (s_cmp(key, "RHS         ", key_len, (ftnlen)12) == 0) {
	    s3optc_(s, cw + 424, c__, (ftnlen)8, (ftnlen)8);
	} else if (s_cmp(key, "UPPER       ", key_len, (ftnlen)12) == 0) {
	    s3optr_(s, &rw[97], r__);
	} else {
	    more = TRUE_;
	}
    }
    if (! more) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Error messages. */
/*     This is the only way we can think of to concatenate strings */
/*     when one of them is of indeterminate length. */
/*     ------------------------------------------------------------------ */
    s_wsfi(&io___40);
    do_fio(&c__1, " XXX  Keyword not recognized:         ", (ftnlen)38);
    do_fio(&c__1, key, key_len);
    e_wsfi();
    goto L895;
L820:
    s_wsfi(&io___41);
    do_fio(&c__1, " XXX  Second keyword not recognized:  ", (ftnlen)38);
    do_fio(&c__1, key2, (ftnlen)16);
    e_wsfi();
    goto L895;
/* L840: */
    s_wsfi(&io___42);
    do_fio(&c__1, " XXX  Fourth keyword not recognized:  ", (ftnlen)38);
    do_fio(&c__1, key2, (ftnlen)16);
    e_wsfi();
    goto L895;
L880:
    s_wsfi(&io___43);
    do_fio(&c__1, " XXX  parm subscript out of range:    ", (ftnlen)38);
    do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
    e_wsfi();
    goto L895;
L890:
    s_copy(str, " XXX  Obsolete option", (ftnlen)132, (ftnlen)21);
    goto L895;
/*     The buffer should have been output already to the Print file. */
/*     First output it to the Summary file. */
/*     Then print the error message. */
L895:
    ++(*errors);
    if (*lsumm > 0) {
	s_wsfi(&io___45);
	do_fio(&c__1, buffer, maxbuf);
	e_wsfi();
	snprnt_(&c__2, str1, &iw[1], leniw, (ftnlen)132);
    }
    snprnt_(&c__4, str, &iw[1], leniw, (ftnlen)132);
L900:
    return 0;
} /* s3opt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3opt */
/* Subroutine */ int s3tie_(char *tie, integer *loc, ftnlen tie_len)
{
    /* Initialized data */

    static char ties[16*70] = "(TOTAL)         " "(USER)          " "ALL    "
	    "         " "BASIC           " "BASIS           " "BOUND         "
	    "  " "CG              " "CHARACTER       " "CHOLESKY        " 
	    "COLUMNS         " "COMPLETE        " "CONSTRAINTS     " "DAMPIN"
	    "G         " "DEBUG           " "DENSE           " "DENSITY      "
	    "   " "DERIVATIVE      " "DIAGONAL        " "DIFFERENCES     " 
	    "DIMENSION       " "ELEMENTS        " "EXIT            " "FACTOR"
	    "IZATION   " "FEASIBILITY     " "FILE            " "FLUSH        "
	    "   " "FREQUENCY       " "FULL            " "GRADIENTS       " 
	    "INFORMATION     " "INTEGER         " "ITERATIONS      " "JACOBI"
	    "AN        " "LEVEL           " "LIMITED         " "LINEAR       "
	    "   " "LINESEARCH      " "LOG             " "MODE            " 
	    "NEWTON          " "NO              " "NONLINEAR       " "OBJECT"
	    "IVE       " "OPTIMALITY      " "OPTION          " "PARAMETERS   "
	    "   " "PARTIAL         " "PHASE1          " "PHASE2          " 
	    "POINT           " "PRECONDITIONING " "PRINT           " "QN    "
	    "          " "REAL            " "ROOK            " "ROW          "
	    "   " "SINGULARITY     " "SOLVER          " "SPARSE          " 
	    "START           " "STEP            " "STOP            " "SUPERB"
	    "ASICS     " "SWAP            " "TOLERANCE       " "UPDATES      "
	    "   " "VARIABLES       " "WEIGHT          " "WEIGHTMAX       " 
	    "YES             ";

    extern /* Subroutine */ int oplook_(integer *, char *, logical *, char *, 
	    integer *, ftnlen, ftnlen);

/*     ================================================================== */
/*     s3tie  sets key to be the standard form for the second keyword */
/*     on each line of a SPECS file. */

/*     21 May 1998: First version of s3tie. */
/*     10 Dec 2002: New ties: 'Diagonal' and 'Rook' */
/*     17 Jan 2003: New tie: 'Solver' */
/*     01 Aug 2003: New ties for QPsolver: 'Cholesky', 'CG', 'QN' */
/*     22 Jun 2004: New tie for 'System': 'Information' */
/*     02 Sep 2007: New tie for 'Sticky': 'Parameters' */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    oplook_(&c__70, ties, &c_true, tie, loc, (ftnlen)16, (ftnlen)16);
    return 0;
} /* s3tie_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3tie */
/* Subroutine */ int s3unsetall_(char *cw, integer *lencw, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;

/*     ================================================================== */
/*     s3unsetAll sets all optional parameters as undefined. */

/*     07 Feb 1998: First version (with name s3undf). */
/*     10 Jun 2007: Initialize everything, not just the options. */
/*     18 Jun 2007: Switched to an unrolled loop */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    for (i__ = 1; i__ <= 500; i__ += 4) {
	s_copy(cw + (i__ << 3), "-1111111", (ftnlen)8, (ftnlen)8);
	s_copy(cw + (i__ + 1 << 3), "-1111111", (ftnlen)8, (ftnlen)8);
	s_copy(cw + (i__ + 2 << 3), "-1111111", (ftnlen)8, (ftnlen)8);
	s_copy(cw + (i__ + 3 << 3), "-1111111", (ftnlen)8, (ftnlen)8);
	iw[i__] = -11111;
	iw[i__ + 1] = -11111;
	iw[i__ + 2] = -11111;
	iw[i__ + 3] = -11111;
	rw[i__] = -11111.;
	rw[i__ + 1] = -11111.;
	rw[i__ + 2] = -11111.;
	rw[i__ + 3] = -11111.;
    }
    return 0;
} /* s3unsetall_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3unsetAll */
/* Subroutine */ int s3unsetprm_(char *cw, integer *lencw, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;

/*     ================================================================== */
/*     s3unsetPrm sets all optional parameters as undefined. */
/*     Called by the MPS main programs. */

/*     07 Feb 1998: First version of s3undf. */
/*     18 Jun 2007: Renamed s3unsetOpt. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    for (i__ = 51; i__ <= 180; ++i__) {
	s_copy(cw + (i__ << 3), "-1111111", (ftnlen)8, (ftnlen)8);
	iw[i__] = -11111;
	rw[i__] = -11111.;
    }
    return 0;
} /* s3unsetprm_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* subroutine s3unsetPrm */
/* Subroutine */ int oplook_(integer *ndict, char *dictry, logical *alpha, 
	char *key, integer *entry__, ftnlen dictry_len, ftnlen key_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    logical l_lt(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static char flag__[16];
    static integer mark, last, first, length;
    extern /* Subroutine */ int opscan_(char *, integer *, integer *, integer 
	    *, ftnlen);
    static char target[16];



/* Description and usage: */

/*       Performs dictionary lookups.  A pointer is returned if a */
/*    match is found between the input key and the corresponding */
/*    initial characters of one of the elements of the dictionary. */
/*    If a "synonym" has been provided for an entry, the search is */
/*    continued until a match to a primary dictionary entry is found. */
/*    Cases of no match, or multiple matches, are also provided for. */

/*     Dictionary entries must be left-justified, and may be alphabetized */
/*    for faster searches.  Secondary entries, if any, are composed of */
/*    two words separated by one or more characters such as blank, tab, */
/*    comma, colon, or equal sign which are treated as non-significant */
/*    by opscan.  The first entry of each such pair serves as a synonym */
/*    for the second, more fundamental keyword. */

/*       The ordered search stops after the section of the dictionary */
/*    having the same first letters as the key has been checked, or */
/*    after a specified number of entries have been examined.  A special */
/*    dictionary entry, the vertical bar '|', will also terminate the */
/*    search.  This will speed things up if an appropriate dictionary */
/*    length parameter cannot be determined.  Both types of search are */
/*    sequential.  See "Notes" below for some suggestions if efficiency */
/*    is an issue. */


/* Parameters: */

/*    Name    Dimension  Type  I/O/S  Description */
/*    NDICT               I    I      Number of dictionary entries to be */
/*                                    examined. */
/*    DICTRY  NDICT       C    I      Array of dictionary entries, */
/*                                    left-justified in their fields. */
/*                                    May be alphabetized for efficiency, */
/*                                    in which case ALPHA should be */
/*                                    .TRUE.  Entries with synonyms are */
/*                                    of the form */
/*                                    'ENTRY : SYNONYM', where 'SYNONYM' */
/*                                    is a more fundamental entry in the */
/*                                    same dictionary.  NOTE: Don't build */
/*                                    "circular" dictionaries! */
/*    ALPHA               L    I      Indicates whether the dictionary */
/*                                    is in alphabetical order, in which */
/*                                    case the search can be terminated */
/*                                    sooner. */
/*    KEY                 C    I/O    String to be compared against the */
/*                                    dictionary.  Abbreviations are OK */
/*                                    if they correspond to a unique */
/*                                    entry in the dictionary.  KEY is */
/*                                    replaced on termination by its most */
/*                                    fundamental equivalent dictionary */
/*                                    entry (uppercase, left-justified) */
/*                                    if a match was found. */
/*    ENTRY               I      O    Dictionary pointer.  If > 0, it */
/*                                    indicates which entry matched KEY. */
/*                                    In case of trouble, a negative */
/*                                    value means that a UNIQUE match */
/*                                    was not found - the absolute value */
/*                                    of ENTRY points to the second */
/*                                    dictionary entry that matched KEY. */
/*                                    Zero means that NO match could be */
/*                                    found.  ENTRY always refers to the */
/*                                    last search performed - */
/*                                    in searching a chain of synonyms, */
/*                                    a non-positive value will be */
/*                                    returned if there is any break, */
/*                                    even if the original input key */
/*                                    was found. */


/* External references: */

/*    Name    Description */
/*    opscan  Finds first and last significant characters. */


/* Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77). */
/*               Appears to satisfy the ANSI Fortran 77 standard. */


/* Notes: */

/*    (1)  IMPLICIT NONE is non-standard. */

/*    (2)  We have assumed that the dictionary is not too big.  If */
/*         many searches are to be done or if the dictionary has more */
/*         than a dozen or so entries, it may be advantageous to build */
/*         an index array of pointers to the beginning of the section */
/*         of the dictionary containing each letter, then pass in the */
/*         portion of the dictionary beginning with DICTRY (INDEX). */
/*         (This won't generally work for dictionaries with synonyms.) */
/*         For very large problems, a completely different approach may */
/*         be advisable, e.g. a binary search for ordered dictionaries. */

/*    (3)  oplook is case sensitive.  In most applications it will be */
/*         necessary to use an uppercase dictionary, and to convert the */
/*         input key to uppercase before calling oplook.  Companion */
/*         routines optokn and PAIRS, available from the author, already */
/*         take care of this. */

/*    (4)  The key need not be left-justified.  Any leading (or */
/*         trailing) characters which are "non-significant" to opscan */
/*         will be ignored.  These include blanks, horizontal tabs, */
/*         commas, colons, and equal signs.  See opscan for details. */

/*    (5)  The ASCII collating sequence for character data is assumed. */
/*         (N.B. This means the numerals precede the alphabet, unlike */
/*         common practice!)  This should not cause trouble on EBCDIC */
/*         machines if DICTRY just contains alphabetic keywords. */
/*         Otherwise it may be necessary to use the FORTRAN lexical */
/*         library routines to force use of the ASCII sequence. */

/*    (6)  Parameter NUMSIG sets a limit on the length of significant */
/*         dictionary entries.  Special applications may require that */
/*         this be increased.  (It is 16 in the present version.) */

/*    (7)  No protection against "circular" dictionaries is provided: */
/*         don't claim that A is B, and that B is A.  All synonym chains */
/*         must terminate!  Other potential errors not checked for */
/*         include duplicate or mis-ordered entries. */

/*    (8)  The handling of ambiguities introduces some ambiguity: */

/*            ALPHA = .TRUE.  A potential problem, when one entry */
/*                            looks like an abbreviation for another */
/*                            (eg. does 'A' match 'A' or 'AB'?) was */
/*                            resolved by dropping out of the search */
/*                            immediately when an "exact" match is found. */

/*            ALPHA = .FALSE. The programmer must ensure that the above */
/*                            situation does not arise: each dictionary */
/*                            entry must be recognizable, at least when */
/*                            specified to full length.  Otherwise, the */
/*                            result of a search will depend on the */
/*                            order of entries. */


/* Author:  Robert Kennelly, Informatics General Corporation. */


/* Development history: */

/*    24 Feb. 1984  RAK/DAS  Initial design and coding. */
/*    25 Feb. 1984    RAK    Combined the two searches by suitable */
/*                           choice of terminator FLAG. */
/*    28 Feb. 1984    RAK    Optional synonyms in dictionary, no */
/*                           longer update KEY. */
/*    29 Mar. 1984    RAK    Put back replacement of KEY by its */
/*                           corresponding entry. */
/*    21 June 1984    RAK    Corrected bug in error handling for cases */
/*                           where no match was found. */
/*    23 Apr. 1985    RAK    Introduced test for exact matches, which */
/*                           permits use of dictionary entries which */
/*                           would appear to be ambiguous (for ordered */
/*                           case).  Return -I to point to the entry */
/*                           which appeared ambiguous (had been -1). */
/*                           Repaired loop termination - had to use */
/*                           equal length strings or risk quitting too */
/*                           soon when one entry is an abbreviation */
/*                           for another.  Eliminated HIT, reduced */
/*                           NUMSIG to 16. */
/*    15 Nov. 1985    MAS    Loop 20 now tests .LT. FLAG, not .LE. FLAG. */
/*                           If ALPHA is false, FLAG is now '|', not '{'. */
/*    26 Jan. 1986    PEG    Declaration of FLAG and TARGET modified to */
/*                           conform to ANSI-77 standard. */
/*    05 Dec  2004    PEG    Used intrinsics LGE, LLE */
/* ----------------------------------------------------------------------- */
/*     Variable declarations. */
/*     ---------------------- */
/*     Parameters. */
/*     Variables. */
/*     CHARACTER */
/*    &   DICTRY (NDICT) * (*), FLAG * (NUMSIG), */
/*    &   KEY * (*), TARGET * (NUMSIG) */
/*     Procedures. */
/*     Executable statements. */
/*     ---------------------- */
    /* Parameter adjustments */
    dictry -= dictry_len;

    /* Function Body */
    *entry__ = 0;
/*     Isolate the significant portion of the input key (if any). */
    first = 1;
/* Computing MIN */
    i__1 = i_len(key, key_len);
    last = min(i__1,16);
    opscan_(key, &first, &last, &mark, key_len);
    if (mark > 0) {
	s_copy(target, key + (first - 1), (ftnlen)16, mark - (first - 1));
/*        Look up TARGET in the dictionary. */
L10:
	length = mark - first + 1;
/*           Select search strategy by cunning choice of termination test */
/*           flag.  The vertical bar is just about last in both the */
/*           ASCII and EBCDIC collating sequences. */
	if (*alpha) {
	    s_copy(flag__, target, (ftnlen)16, (ftnlen)16);
	} else {
	    s_copy(flag__, "|", (ftnlen)16, (ftnlen)1);
	}
/*           Perform search. */
/*           --------------- */
	i__ = 0;
L20:
	++i__;
	if (s_cmp(target, dictry + i__ * dictry_len, length, length) == 0) {
	    if (*entry__ == 0) {
/*                    First "hit" - must still guard against ambiguities */
/*                    by searching until we've gone beyond the key */
/*                    (ordered dictionary) or until the end-of-dictionary */
/*                    mark is reached (exhaustive search). */
		*entry__ = i__;
/*                    Special handling if match is exact - terminate */
/*                    search.  We thus avoid confusion if one dictionary */
/*                    entry looks like an abbreviation of another. */
/*                    This fix won't generally work for un-ordered */
/*                    dictionaries! */
		first = 1;
		last = 16;
		opscan_(dictry + *entry__ * dictry_len, &first, &last, &mark, 
			dictry_len);
		if (mark == length) {
		    i__ = *ndict;
		}
	    } else {
/*                    Oops - two hits!  Abnormal termination. */
/*                    --------------------------------------- */
		*entry__ = -i__;
		return 0;
	    }
	}
/*           Check whether we've gone past the appropriate section of the */
/*           dictionary.  The test on the index provides insurance and an */
/*           optional means for limiting the extent of the search. */
	if (l_lt(dictry + i__ * dictry_len, flag__, length, (ftnlen)16) && 
		i__ < *ndict) {
	    goto L20;
	}
/*           Check for a synonym. */
/*           -------------------- */
	if (*entry__ > 0) {
/*              Look for a second entry "behind" the first entry.  FIRST */
/*              and MARK were determined above when the hit was detected. */
	    first = mark + 2;
	    opscan_(dictry + *entry__ * dictry_len, &first, &last, &mark, 
		    dictry_len);
	    if (mark > 0) {
/*                 Re-set target and dictionary pointer, then repeat the */
/*                 search for the synonym instead of the original key. */
		s_copy(target, dictry + (*entry__ * dictry_len + (first - 1)),
			 (ftnlen)16, mark - (first - 1));
		*entry__ = 0;
		goto L10;
	    }
	}
    }
    if (*entry__ > 0) {
	s_copy(key, dictry + *entry__ * dictry_len, key_len, dictry_len);
    }
/*     Normal termination. */
/*     ------------------- */
    return 0;
/*     End of oplook */
} /* oplook_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
logical opnumb_(char *string, ftnlen string_len)
{
    /* System generated locals */
    logical ret_val;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    logical l_ge(char *, char *, ftnlen, ftnlen), l_le(char *, char *, ftnlen,
	     ftnlen);

    /* Local variables */
    static integer j;
    static char atom[1];
    static integer nexp, nplus, ndigit, length;
    static logical number;
    static integer npoint, nminus;

/* *********************************************************************** */
/*     Description and usage: */

/*        A simple(-minded) test for numeric data is implemented by */
/*        searching an input string for legitimate characters: */
/*                digits 0 to 9, D, E, -, + and . */
/*        Insurance is provided by requiring that a numeric string */
/*        have at least one digit, at most one D, E or . */
/*        and at most two -s or +s.  Note that a few ambiguities remain: */

/*           (a)  A string might have the form of numeric data but be */
/*                intended as text.  No general test can hope to detect */
/*                such cases. */

/*           (b)  There is no check for correctness of the data format. */
/*                For example a meaningless string such as 'E1.+2-' */
/*                will be accepted as numeric. */

/*        Despite these weaknesses, the method should work in the */
/*        majority of cases. */


/*     Parameters: */

/*        Name    Dimension  Type  I/O/S  Description */
/*        opnumb              L      O    Set .TRUE. if STRING appears */
/*                                        to be numerical data. */
/*        STRING              C    I      Input data to be tested. */


/*     Environment:  ANSI FORTRAN 77. */


/*     Notes: */

/*        (1)  It is assumed that STRING is a token extracted by */
/*             optokn, which will have converted any lower-case */
/*             characters to upper-case. */

/*        (2)  optokn pads STRING with blanks, so that a genuine */
/*             number is of the form  '1234        '. */
/*             Hence, the scan of STRING stops at the first blank. */

/*        (3)  COMPLEX data with parentheses will not look numeric. */


/*     Systems Optimization Laboratory, Stanford University. */
/*     12 Nov  1985    Initial design and coding, starting from the */
/*                     routine ALPHA from Informatics General, Inc. */
/*     05 Dec  2004    Used intrinsics LGE, LLE */
/* *********************************************************************** */
    ndigit = 0;
    nexp = 0;
    nminus = 0;
    nplus = 0;
    npoint = 0;
    number = TRUE_;
    length = i_len(string, string_len);
    j = 0;
L10:
    ++j;
    *(unsigned char *)atom = *(unsigned char *)&string[j - 1];
    if (l_ge(atom, "0", (ftnlen)1, (ftnlen)1) && l_le(atom, "9", (ftnlen)1, (
	    ftnlen)1)) {
	++ndigit;
    } else if (*(unsigned char *)atom == 'D' || *(unsigned char *)atom == 'E')
	     {
	++nexp;
    } else if (*(unsigned char *)atom == '-') {
	++nminus;
    } else if (*(unsigned char *)atom == '+') {
	++nplus;
    } else if (*(unsigned char *)atom == '.') {
	++npoint;
    } else if (*(unsigned char *)atom == ' ') {
	j = length;
    } else {
	number = FALSE_;
    }
    if (number && j < length) {
	goto L10;
    }
    ret_val = number && ndigit >= 1 && nexp <= 1 && nminus <= 2 && nplus <= 2 
	    && npoint <= 1;
    return ret_val;
/*     End of opnumb */
} /* opnumb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* Subroutine */ int opscan_(char *string, integer *first, integer *last, 
	integer *mark, ftnlen string_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char c__[1], ht[1];
    static integer end, begin, length;



/* Description and usage: */

/*       Looks for non-blank fields ("tokens") in a string, where the */
/*    fields are of arbitrary length, separated by blanks, tabs, commas, */
/*    colons, or equal signs.  The position of the end of the 1st token */
/*    is also returned, so this routine may be conveniently used within */
/*    a loop to process an entire line of text. */

/*       The procedure examines a substring, STRING (FIRST : LAST), which */
/*    may of course be the entire string (in which case just call opscan */
/*    with FIRST <= 1 and LAST >= LEN (STRING) ).  The indices returned */
/*    are relative to STRING itself, not the substring. */


/* Parameters: */

/*    Name    Dimension  Type  I/O/S  Description */
/*    STRING              C    I      Text string containing data to be */
/*                                    scanned. */
/*    FIRST               I    I/O    Index of beginning of substring. */
/*                                    If <= 1, the search begins with 1. */
/*                                    Output is index of beginning of */
/*                                    first non-blank field, or 0 if no */
/*                                    token was found. */
/*    LAST                I    I/O    Index of end of substring. */
/*                                    If >= LEN (STRING), the search */
/*                                    begins with LEN (STRING).  Output */
/*                                    is index of end of last non-blank */
/*                                    field, or 0 if no token was found. */
/*    MARK                I      O    Points to end of first non-blank */
/*                                    field in the specified substring. */
/*                                    Set to 0 if no token was found. */


/* Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77). */
/*               ANSI Fortran 77, except for the tab character HT. */

/* Notes: */

/*    (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined */
/*         in a non-standard way:  the CHAR function is not permitted */
/*         in a PARAMETER declaration (OK on VAX, though).  For Absoft */
/*         FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it */
/*         may be best to declare HT as a variable and assign */
/*         HT = CHAR(9) on ASCII machines, or CHAR(5) for EBCDIC. */

/*    (2)  The pseudo-recursive structure was chosen for fun.  It is */
/*         equivalent to three DO loops with embedded GO TOs in sequence. */

/*    (3)  The variety of separators recognized limits the usefulness of */
/*         this routine somewhat.  The intent is to facilitate handling */
/*         such tokens as keywords or numerical values.  In other */
/*         applications, it may be necessary for ALL printing characters */
/*         to be significant.  A simple modification to statement */
/*         function SOLID will do the trick. */


/* Author:  Robert Kennelly, Informatics General Corporation. */


/* Development history: */

/*    29 Dec. 1984    RAK    Initial design and coding, (very) loosely */
/*                           based on SCAN_STRING by Ralph Carmichael. */
/*    25 Feb. 1984    RAK    Added ':' and '=' to list of separators. */
/*    16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY */
/*                           (previous re-use of STRING was ambiguous). */

/* ----------------------------------------------------------------------- */
/*     Variable declarations. */
/*     ---------------------- */
/*     Parameters. */
/*     Variables. */
/*     LOGICAL */
/*    &   SOLID */
/*     Statement functions. */

/*     SOLID (DUMMY) = (DUMMY .NE. BLANK) .AND. */
/*    &                (DUMMY .NE. COLON) .AND. */
/*    &                (DUMMY .NE. COMMA) .AND. */
/*    &                (DUMMY .NE. EQUAL) .AND. */
/*    &                (DUMMY .NE. HT) */
/*     Executable statements. */
/*     ---------------------- */
/* ***  HT     = CHAR(9) for ASCII machines, CHAR(5) for EBCDIC. */
    *(unsigned char *)ht = '\t';
    *mark = 0;
    length = i_len(string, string_len);
    begin = max(*first,1);
    end = min(length,*last);
/*     Find the first significant character ... */
    i__1 = end;
    for (*first = begin; *first <= i__1; ++(*first)) {
	*(unsigned char *)c__ = *(unsigned char *)&string[*first - 1];
/* IF ( SOLID(c) ) THEN */
	if (*(unsigned char *)c__ != ' ' && *(unsigned char *)c__ != ':' && *(
		unsigned char *)c__ != ',' && *(unsigned char *)c__ != '=' && 
		*(unsigned char *)c__ != *(unsigned char *)ht) {
/*           ... then the end of the first token ... */
	    i__2 = end - 1;
	    for (*mark = *first; *mark <= i__2; ++(*mark)) {
		i__3 = *mark;
		s_copy(c__, string + i__3, (ftnlen)1, *mark + 1 - i__3);
/* IF (.NOT.SOLID(c) ) THEN */
		if (*(unsigned char *)c__ != ' ' && *(unsigned char *)c__ != 
			':' && *(unsigned char *)c__ != ',' && *(unsigned 
			char *)c__ != '=' && *(unsigned char *)c__ != *(
			unsigned char *)ht) {
/* relax */
		} else {
/*                 ... and finally the last significant character. */
		    i__3 = *mark;
		    for (*last = end; *last >= i__3; --(*last)) {
			*(unsigned char *)c__ = *(unsigned char *)&string[*
				last - 1];
/* IF ( SOLID(c) ) THEN */
			if (*(unsigned char *)c__ != ' ' && *(unsigned char *)
				c__ != ':' && *(unsigned char *)c__ != ',' && 
				*(unsigned char *)c__ != '=' && *(unsigned 
				char *)c__ != *(unsigned char *)ht) {
			    return 0;
			}
/* L10: */
		    }
/*                 Everything past the first token was a separator. */
		    ++(*last);
		    return 0;
		}
/* L20: */
	    }
/*           There was nothing past the first token. */
	    *last = *mark;
	    return 0;
	}
/* L30: */
    }
/*     Whoops - the entire substring STRING (BEGIN : END) was composed of */
/*     separators ! */
    *first = 0;
    *mark = 0;
    *last = 0;
    return 0;
/*     End of opscan */
} /* opscan_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* Subroutine */ int optokn_(char *string, integer *numin, integer *numout, 
	char *list, ftnlen string_len, ftnlen list_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, mark, last, first, count;
    extern /* Subroutine */ int opscan_(char *, integer *, integer *, integer 
	    *, ftnlen), opuppr_(char *, ftnlen);



/* Description and usage: */

/*       An aid to parsing input data.  The individual "tokens" in a */
/*    character string are isolated, converted to uppercase, and stored */
/*    in an array.  Here, a token is a group of significant, contiguous */
/*    characters.  The following are NON-significant, and hence may */
/*    serve as separators:  blanks, horizontal tabs, commas, colons, */
/*    and equal signs.  See opscan for details.  Processing continues */
/*    until the requested number of tokens have been found or the end */
/*    of the input string is reached. */


/* Parameters: */

/*    Name    Dimension  Type  I/O/S  Description */
/*    STRING              C    I      Input string to be analyzed. */
/*    NUMIN               I    I/O    Number of tokens requested (input) */
/*                                    and found (output). */
/*    LIST    NUMIN       C      O    Array of tokens, changed to upper */
/*                                    case. */


/* External references: */

/*    Name    Description */
/*    opscan  Finds positions of first and last significant characters. */
/*    opuppr  Converts a string to uppercase. */


/* Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77). */
/*               Appears to satisfy the ANSI Fortran 77 standard. */


/* Notes: */

/*    (1)  IMPLICIT NONE is non-standard. */


/* Author:  Robert Kennelly, Informatics General Corporation. */


/* Development history: */

/*    16 Jan. 1984    RAK    Initial design and coding. */
/*    16 Mar. 1984    RAK    Revised header to reflect full list of */
/*                           separators, repaired faulty WHILE clause */
/*                           in "10" loop. */
/*    18 Sep. 1984    RAK    Change elements of LIST to uppercase one */
/*                           at a time, leaving STRING unchanged. */
/*    05 Dec. 2004    PEG    Replaced by NUMIN, NUMOUT */

/* ----------------------------------------------------------------------- */
/*     Variable declarations. */
/*     ---------------------- */
/*     Parameters. */
/*     Variables. */
/*     Procedures. */
/*     Executable statements. */
/*     ---------------------- */
/*     WHILE there are tokens to find, loop UNTIL enough have been found. */
    /* Parameter adjustments */
    list -= list_len;

    /* Function Body */
    first = 1;
    last = i_len(string, string_len);
    count = 0;
L10:
/*        Get delimiting indices of next token, if any. */
    opscan_(string, &first, &last, &mark, string_len);
    if (last > 0) {
	++count;
/*           Pass token to output string array, then change case. */
	s_copy(list + count * list_len, string + (first - 1), list_len, mark 
		- (first - 1));
	opuppr_(list + count * list_len, list_len);
	first = mark + 2;
	if (count < *numin) {
	    goto L10;
	}
    }
/*     Fill the rest of LIST with blanks and set NUMOUT for output. */
    i__1 = *numin;
    for (i__ = count + 1; i__ <= i__1; ++i__) {
	s_copy(list + i__ * list_len, " ", list_len, (ftnlen)1);
/* L20: */
    }
    *numout = count;
/*     Termination. */
/*     ------------ */
    return 0;
/*     End of optokn */
} /* optokn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* Subroutine */ int opuppr_(char *string, ftnlen string_len)
{
    /* Initialized data */

    static char low[26] = "abcdefghijklmnopqrstuvwxyz";
    static char upp[26] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    logical l_ge(char *, char *, ftnlen, ftnlen), l_le(char *, char *, ftnlen,
	     ftnlen);
    integer i_indx(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char c__[1];
    static integer i__, j;


/* ACRONYM:  UPper CASE */

/* PURPOSE:  This subroutine changes all lower case letters in the */
/*           character string to upper case. */

/* METHOD:   Each character in STRING is treated in turn.  The intrinsic */
/*           function INDEX effectively allows a table lookup, with */
/*           the local strings LOW and UPP acting as two tables. */
/*           This method avoids the use of CHAR and ICHAR, which appear */
/*           be different on ASCII and EBCDIC machines. */

/* ARGUMENTS */
/*    ARG       DIM     TYPE I/O/S DESCRIPTION */
/*  STRING       !       C   I/O   Character string possibly containing */
/*                                 some lower-case letters on input; */
/*                                 strictly upper-case letters on output */
/*                                 with no change to any non-alphabetic */
/*                                 characters. */

/* EXTERNAL REFERENCES: */
/*  LEN    - Returns the declared length of a CHARACTER variable. */
/*  INDEX  - Returns the position of second string within first. */

/* ENVIRONMENT:  ANSI FORTRAN 77 */

/* DEVELOPMENT HISTORY: */
/*     DATE  INITIALS  DESCRIPTION */
/*   06/28/83   CLH    Initial design. */
/*   01/03/84   RAK    Eliminated NCHAR input. */
/*   06/14/84   RAK    Used integer PARAMETERs in comparison. */
/*   04/21/85   RAK    Eliminated DO/END DO in favor of standard code. */
/*   09/10/85   MAS    Eliminated CHAR,ICHAR in favor of LOW, UPP, INDEX. */
/*   12/04/04   PEG    Used intrinsics LGE, LLE */

/* AUTHOR: Charles Hooper, Informatics General, Palo Alto, CA. */

/* ----------------------------------------------------------------------- */
    i__1 = i_len(string, string_len);
    for (j = 1; j <= i__1; ++j) {
	*(unsigned char *)c__ = *(unsigned char *)&string[j - 1];
	if (l_ge(c__, "a", (ftnlen)1, (ftnlen)1) && l_le(c__, "z", (ftnlen)1, 
		(ftnlen)1)) {
	    i__ = i_indx(low, c__, (ftnlen)26, (ftnlen)1);
	    if (i__ > 0) {
		*(unsigned char *)&string[j - 1] = *(unsigned char *)&upp[i__ 
			- 1];
	    }
	}
/* L10: */
    }
    return 0;
/*     End of opuppr */
} /* opuppr_ */

