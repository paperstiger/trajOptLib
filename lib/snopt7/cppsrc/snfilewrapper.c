/* cppsrc/snfilewrapper.f -- translated by f2c (version 20100827).
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

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     snFilewrapper   snOpenappend   snClose   snOpenread */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snfilewrapper_(char *name__, integer *ispec, integer *
	inform__, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen name_len, ftnlen cw_len)
{
    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int snspec_(integer *, integer *, char *, integer 
	    *, integer *, integer *, doublereal *, integer *, ftnlen);
    static integer iostat;

/*     ================================================================== */
/*     Read options for snopt from the file named name. inform .eq 0 if */
/*     successful. */

/*     09 Jan 2000: First version of snFilewrapper. */
/*     ================================================================== */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    o__1.oerr = 1;
    o__1.ounit = *ispec;
    o__1.ofnmlen = name_len;
    o__1.ofnm = name__;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    iostat = f_open(&o__1);
    if (iostat != 0) {
	*inform__ = iostat + 2;
    } else {
	snspec_(ispec, inform__, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, 
		(ftnlen)8);
	cl__1.cerr = 0;
	cl__1.cunit = *ispec;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* snfilewrapper_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snFilewrapper */
/* Subroutine */ int snopenappend_(integer *iunit, char *name__, integer *
	inform__, ftnlen name_len)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);

/*     ================================================================== */
/*     Open file named name to Fortran unit iunit. inform .eq. 0 if */
/*     sucessful.  Although opening for appending is not in the f77 */
/*     standard, it is understood by f2c. */

/*     09 Jan 2000: First version of snOpenappend */
/*     ================================================================== */
    o__1.oerr = 1;
    o__1.ounit = *iunit;
    o__1.ofnmlen = name_len;
    o__1.ofnm = name__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    *inform__ = f_open(&o__1);
    return 0;
} /* snopenappend_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snOpenappend */
/* Subroutine */ int snclose_(integer *iunit)
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    integer f_clos(cllist *);

/*     ================================================================== */
/*     Close unit iunit. */

/*     09 Jan 2000: First version of snClose */
/*     ================================================================== */
    cl__1.cerr = 0;
    cl__1.cunit = *iunit;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* snclose_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snClose */
/* Subroutine */ int snopenread_(integer *iunit, char *name__, integer *
	inform__, ftnlen name_len)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);

/*     ================================================================== */
/*     Open file called name to Fortran unit iunit. inform .eq. 0 if */
/*     sucessful.  Although opening for appending is not in the f77 */
/*     standard, it is understood by f2c. */

/*     09 Jan 2000: First version of snOpenread */
/*     ================================================================== */
    o__1.oerr = 1;
    o__1.ounit = *iunit;
    o__1.ofnmlen = name_len;
    o__1.ofnm = name__;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = "FORMATTED";
    o__1.oblnk = 0;
    *inform__ = f_open(&o__1);
    return 0;
} /* snopenread_ */

