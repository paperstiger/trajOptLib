/* ./src/sn03prnt.f -- translated by f2c (version 20100827).
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

/*     File  sn03prnt.f  (default versions of snPRNT and snREAD) */

/*     snPRNT   snREAD */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snprnt_(integer *mode, char *string, integer *iw, 
	integer *leniw, ftnlen string_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer m;
    static char form[4];
    static integer istdo, isumm;
    extern /* Subroutine */ int s1trim_(char *, integer *, ftnlen);
    static integer length, iprint, lvlsys;
    extern integer s1outpt_(void);

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, form, 0 };
    static cilist io___9 = { 0, 0, 0, form, 0 };
    static cilist io___10 = { 0, 0, 0, form, 0 };
    static cilist io___11 = { 0, 0, 0, form, 0 };


/*     ================================================================== */
/*     snPRNT  prints a trimmed form of "string" on various files. */
/*     If mode = 0,      nothing is output. */
/*     If mode = 1,      string is output to iPrint. */
/*     If mode = 2,      string is output to iSumm. */
/*     If mode = 3 or 4, string is output to iPrint and iSumm. */
/*     If mode = 4,      string is output to iStdo (standard output) */
/*                       if iPrint and iSumm are both zero.  This mode */
/*                       is intended for error messages. */
/*     If mode = 5,      string is output to iStdo (standard output) */
/*                       This mode is to be used when the elements of */
/*                       the integer work array iw cannot be trusted. */

/*     mode 11-15 are the same as mode 1-5 with blank line before output. */

/*     If mode > 15 then nothing is printed unless  lvlSys > 0. */
/*     mode 21-25 are the same as mode 1-5 */
/*     mode 31-35 are the same as mode 11-15 */

/*     25 Sep 2002: First version of snPRNT. */
/*     31 Jul 2003: mode 11-14 added.  form introduced. */
/*     27 Dec 2003: mode 5 added to allow printing before iw is set. */
/*     12 Mar 2004: s1trim called to trim the string. */
/*     22 Jun 2004: System printing option added. */
/*     22 Jun 2004: Current version of snPRNT. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    lvlsys = iw[71];
/* > 0   => print system info */
    m = 0;
    if (*mode <= 0) {
/*        Relax */
    } else if (*mode < 10) {
	m = *mode;
	s_copy(form, "( a)", (ftnlen)4, (ftnlen)4);
    } else if (*mode < 20) {
/* Blank line first */
	m = *mode - 10;
	s_copy(form, "(/a)", (ftnlen)4, (ftnlen)4);
    } else if (lvlsys > 0) {
/* Print system Info */
	if (*mode < 30) {
	    m = *mode - 20;
	    s_copy(form, "( a)", (ftnlen)4, (ftnlen)4);
	} else {
	    m = *mode - 30;
	    s_copy(form, "(/a)", (ftnlen)4, (ftnlen)4);
	}
    }
    if (m > 0) {
/*        length = len_trim(string)     ! An F90 intrinsic */
	s1trim_(string, &length, string_len);
/* The F77 equivalent */
	if (m == 5) {
	    istdo = s1outpt_();
	    if (istdo > 0) {
		io___6.ciunit = istdo;
		s_wsfe(&io___6);
		do_fio(&c__1, string, length);
		e_wsfe();
	    }
	} else {
	    istdo = iw[10];
/* Standard output */
	    iprint = iw[12];
/* Print file */
	    isumm = iw[13];
/* Summary file */
	    if (m == 1 || m >= 3) {
		if (iprint > 0) {
		    io___9.ciunit = iprint;
		    s_wsfe(&io___9);
		    do_fio(&c__1, string, length);
		    e_wsfe();
		}
	    }
	    if (m == 2 || m >= 3) {
		if (isumm > 0) {
		    io___10.ciunit = isumm;
		    s_wsfe(&io___10);
		    do_fio(&c__1, string, length);
		    e_wsfe();
		}
	    }
	    if (m == 4) {
		if (iprint <= 0 && isumm <= 0) {
		    if (istdo > 0) {
			io___11.ciunit = istdo;
			s_wsfe(&io___11);
			do_fio(&c__1, string, length);
			e_wsfe();
		    }
		}
	    }
	}
    }
    return 0;
} /* snprnt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snPRNT */
/* Subroutine */ int snread_(integer *unitno, char *string, integer *nchar, 
	integer *endfile, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_rsfe(cilist *), e_rsfe(void);

    /* Local variables */
    static char frmt[6];

    /* Fortran I/O blocks */
    static icilist io___13 = { 0, frmt, 0, "(a2,i1,a1)", 6, 1 };
    static icilist io___14 = { 0, frmt, 0, "(a2,i2,a1)", 6, 1 };
    static icilist io___15 = { 0, frmt, 0, "(a2,i3,a1)", 6, 1 };
    static cilist io___16 = { 0, 0, 1, frmt, 0 };


/*     ================================================================== */
/*     snREAD reads a string of length nchar (including trailing blanks) */
/*     from the file with logical unit number  unitno. */

/*     Restriction: 0 < nChar < 1000 */

/*     On exit: */
/*       endfile = 0 means that the string was read successfully. */
/*       endfile = 1 means that an end-of-file was encountered or nchar */
/*                   lies outside the range  0 < nChar < 1000. */

/*     30 Apr 2006: First version of snREAD. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    s_copy(frmt, "      ", (ftnlen)6, (ftnlen)6);
    if (*nchar >= 1 && *nchar <= 999) {
	if (*nchar < 10) {
	    s_wsfi(&io___13);
	    do_fio(&c__1, "(a", (ftnlen)2);
	    do_fio(&c__1, (char *)&(*nchar), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ")", (ftnlen)1);
	    e_wsfi();
	} else if (*nchar < 100) {
	    s_wsfi(&io___14);
	    do_fio(&c__1, "(a", (ftnlen)2);
	    do_fio(&c__1, (char *)&(*nchar), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ")", (ftnlen)1);
	    e_wsfi();
	} else {
	    s_wsfi(&io___15);
	    do_fio(&c__1, "(a", (ftnlen)2);
	    do_fio(&c__1, (char *)&(*nchar), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ")", (ftnlen)1);
	    e_wsfi();
	}
	*endfile = 0;
	io___16.ciunit = *unitno;
	i__1 = s_rsfe(&io___16);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = do_fio(&c__1, string, string_len);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L100;
	}
	return 0;
    }
L100:
    *endfile = 1;
    return 0;
} /* snread_ */

