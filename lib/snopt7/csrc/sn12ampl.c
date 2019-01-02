/* ./src/sn12ampl.f -- translated by f2c (version 20100827).
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

/*     File  sn10ampl.f                    Machine dependent routine */

/*     getfnm */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int getfnm_(integer *lun, char *filnam, integer *last, 
	ftnlen filnam_len)
{
    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

/*     ------------------------------------------------------------------ */
/*     getfnm  is a machine-dependent routine. */
/*     In principal it opens a file with logical unit number lun */
/*     and positions it at the beginning. */

/*     lun     (input) is the unit number. */

/*     16 May 2001: First version, follows some of the advice offered */
/*                  by David Gay, Bell Laboratories. */
/*     ------------------------------------------------------------------ */
    if (*lun <= 9) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 100;
	ici__1.iciunit = filnam;
	ici__1.icifmt = "(a,i1)";
	s_wsfi(&ici__1);
	do_fio(&c__1, "fort.", (ftnlen)5);
	do_fio(&c__1, (char *)&(*lun), (ftnlen)sizeof(integer));
	e_wsfi();
	*last = 6;
    } else {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 100;
	ici__1.iciunit = filnam;
	ici__1.icifmt = "(a,i2)";
	s_wsfi(&ici__1);
	do_fio(&c__1, "fort.", (ftnlen)5);
	do_fio(&c__1, (char *)&(*lun), (ftnlen)sizeof(integer));
	e_wsfi();
	*last = 7;
    }
    return 0;
} /* getfnm_ */

