/* ../snopt7/src/sn40bfil.f -- translated by f2c (version 20100827).
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
static integer c__0 = 0;
static integer c__13 = 13;
static integer c__5 = 5;
static integer c__11 = 11;
static integer c__3 = 3;
static integer c__20 = 20;
static doublereal c_b298 = -1.;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__12 = 12;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn40bfil.f */

/*     s4getB   s4id     s4ksav   s4name   s4inst   s4load   s4oldB */
/*     s4chek   s4chkP   s4dump   s4newB   s4pnch   s4rept   s4savB */
/*     s4soln   s4solp   s4stat */

/*     09 Mar 2004: snSolF implemented and called by s4soln. */
/*     17 Jun 2004: s4savb saves biggest primal and dual infeasibilities */
/*                  before and after scaling, so GAMS can use them. */
/*     27 Dec 2013: s4savB unscales gObj even for Scale option 1. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s4getb_(integer *iexit, integer *m, integer *n, integer *
	nb, integer *nname, integer *ns, integer *iobj, integer *hs, 
	doublereal *bl, doublereal *bu, doublereal *x, char *names, integer *
	iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen names_len)
{
    static integer ioldb;
    extern /* Subroutine */ int s4load_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), s4oldb_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s4inst_(integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, doublereal *, integer *, ftnlen);
    static integer iloadb, iinsrt;

/*     ================================================================== */
/*     s4getb loads one of the basis files. */

/*     15 Nov 1991: First version based on Minos routine m4getb. */
/*     20 Apr 1999: Current version of s4getb. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --hs;
    names -= 8;
    --iw;
    --rw;

    /* Function Body */
    iloadb = iw[122];
/* load file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/*     Load a basis file if one exists and istart = 0 (Cold start). */
/* old basis file */
    if (ioldb > 0) {
	s4oldb_(iexit, m, n, nb, ns, &hs[1], &bl[1], &bu[1], &x[1], &iw[1], 
		leniw, &rw[1], lenrw);
    } else if (iinsrt > 0) {
	s4inst_(n, nb, ns, iobj, &hs[1], &bl[1], &bu[1], &x[1], names + 8, &
		iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
    } else if (iloadb > 0) {
	s4load_(n, nb, ns, iobj, &hs[1], &x[1], names + 8, &iw[1], leniw, &rw[
		1], lenrw, (ftnlen)8);
    }
    return 0;
} /* s4getb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4getB */
/* Subroutine */ int s4id_(integer *j, integer *n, integer *nb, integer *
	nname, char *names, char *id, ftnlen names_len, ftnlen id_len)
{
    /* Initialized data */

    static char colnm[1] = "x";
    static char rownm[1] = "r";

    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__;

/*     ================================================================== */
/*     s4id   returns a name id for the j-th variable. */
/*     If nName = nb, the name is already in Names. */
/*     Otherwise nName = 1. Some generic column or row name is cooked up */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4id. */
/*     16 Sep 1997: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    names -= 8;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    if (*nname == *nb) {
	s_copy(id, names + (*j << 3), (ftnlen)8, (ftnlen)8);
    } else if (*j <= *n) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 8;
	ici__1.iciunit = id;
	ici__1.icifmt = "(a1,i7)";
	s_wsfi(&ici__1);
	do_fio(&c__1, colnm, (ftnlen)1);
	do_fio(&c__1, (char *)&(*j), (ftnlen)sizeof(integer));
	e_wsfi();
    } else {
	i__ = *j - *n;
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 8;
	ici__1.iciunit = id;
	ici__1.icifmt = "(a1,i7)";
	s_wsfi(&ici__1);
	do_fio(&c__1, rownm, (ftnlen)1);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfi();
    }
    return 0;
} /* s4id_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4id */
/* Subroutine */ int s4ksav_(integer *minimz, integer *m, integer *n, integer 
	*nb, integer *ns, integer *mbs, integer *itn, integer *ninf, 
	doublereal *sinf, doublereal *f, integer *kbs, integer *hs, 
	doublereal *ascale, doublereal *bl, doublereal *bu, doublereal *x, 
	doublereal *xbs, char *cw, integer *lencw, integer *iw, integer *
	leniw, ftnlen cw_len)
{
    static integer k, iback, inewb;
    extern /* Subroutine */ int s4newb_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, char *, integer *, integer *, integer *, 
	    ftnlen, ftnlen), s4stat_(integer *, char *, ftnlen);
    static char istate[4*3];
    static integer itnlim;

/*     ================================================================== */
/*     s4ksav  saves various quantities as determined by the frequency */
/*     control ksav. */

/*     15 Nov 1991: First version. */
/*     20 Apr 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --ascale;
    --hs;
    --xbs;
    --kbs;
    cw -= 8;
    --iw;

    /* Function Body */
    iback = iw[120];
/* backup file */
    inewb = iw[124];
/* new basis file */
    itnlim = iw[89];
/* limit on total iterations */
    if (inewb > 0 && *itn < itnlim) {
	k = 0;
	s4stat_(&k, istate, (ftnlen)4);
	s4newb_(&c__0, &inewb, minimz, m, n, nb, ns, mbs, itn, ninf, sinf, f, 
		&kbs[1], &hs[1], &ascale[1], &bl[1], &bu[1], &x[1], &xbs[1], 
		istate, cw + 8, lencw, &iw[1], leniw, (ftnlen)4, (ftnlen)8);
	if (iback > 0) {
	    s4newb_(&c__0, &iback, minimz, m, n, nb, ns, mbs, itn, ninf, sinf,
		     f, &kbs[1], &hs[1], &ascale[1], &bl[1], &bu[1], &x[1], &
		    xbs[1], istate, cw + 8, lencw, &iw[1], leniw, (ftnlen)4, (
		    ftnlen)8);
	}
    }
    return 0;
} /* s4ksav_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4ksav */
/* Subroutine */ int s4name_(integer *n, char *names, char *id, integer *
	ncard, integer *notfnd, integer *maxmsg, integer *j1, integer *j2, 
	integer *jmark, integer *jfound, integer *iw, integer *leniw, ftnlen 
	names_len, ftnlen id_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 XXX  Line\002,i6,\002  --  name not foun"
	    "d:\002,8x,a8)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfi(icilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfi(void);

    /* Local variables */
    static integer j;
    static char str[60];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___14 = { 0, str, 0, fmt_1000, 60, 1 };


/*     ================================================================== */
/*     s4name searches for names in the array  Names(j), j = j1, j2. */
/*     jmark  will probably speed the search on the next entry. */
/*     Used by subroutines s3mpsc, s4inst, s4load. */

/*     Left-justified alphanumeric data is being tested for a match. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4name. */
/*     01 Aug 2003: snPRNT adopted. */
/*     08 Oct 2003: Current version of s4name. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    names -= 8;
    --iw;

    /* Function Body */
    i__1 = *j2;
    for (j = *jmark; j <= i__1; ++j) {
	if (s_cmp(id, names + (j << 3), (ftnlen)8, (ftnlen)8) == 0) {
	    goto L100;
	}
    }
    i__1 = *jmark;
    for (j = *j1; j <= i__1; ++j) {
	if (s_cmp(id, names + (j << 3), (ftnlen)8, (ftnlen)8) == 0) {
	    goto L100;
	}
    }
/*     Not found. */
    *jfound = 0;
    *jmark = *j1;
    ++(*notfnd);
    if (*notfnd <= *maxmsg) {
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*ncard), (ftnlen)sizeof(integer));
	do_fio(&c__1, id, (ftnlen)8);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)60);
    }
    return 0;
/*     Got it. */
L100:
    *jfound = j;
    *jmark = j;
    return 0;
} /* s4name_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4name */
/* Subroutine */ int s4inst_(integer *n, integer *nb, integer *ns, integer *
	iobj, integer *hs, doublereal *bl, doublereal *bu, doublereal *x, 
	char *names, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw, ftnlen names_len)
{
    /* Initialized data */

    static char lll[4] = " LL ";
    static char lul[4] = " UL ";
    static char lxl[4] = " XL ";
    static char lxu[4] = " XU ";
    static char lsb[4] = " SB ";
    static char lend[4] = "ENDA";

    /* Format strings */
    static char fmt_1999[] = "(\002 INSERT file to be input from file\002,i4)"
	    ;
    static char fmt_1000[] = "(14x,2a4,2x,3a4)";
    static char fmt_2000[] = "(\002 NAME\002,10x,2a4,2x,3a4)";
    static char fmt_1020[] = "(a4,a8,2x,a8,2x,e12.5)";
    static char fmt_2010[] = "(\002 XXX  Line\002,i6,\002  ignored:\002,8x,a"
	    "4,a8,2x,a8)";
    static char fmt_2050[] = "(\002 No. of lines read      \002,i6,\002  Lin"
	    "es ignored\002,i6)";
    static char fmt_2060[] = "(\002 No. of basics specified\002,i6,\002  Sup"
	    "erbasics  \002,i6)";

    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_rsfe(cilist *), e_rsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), f_rew(alist *);

    /* Local variables */
    static integer j, l, l1;
    static char id[4*5];
    static doublereal xj;
    static integer nbs;
    static char key[4], str[80];
    static integer jobj, ndum;
    static char name1[8], name2[8];
    static integer ncard;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer jmark, lmark, istdi;
    static doublereal bplus;
    static integer nloop;
    extern /* Subroutine */ int s4name_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal infbnd;
    static integer ignord, notfnd, iinsrt, mpserr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___27 = { 0, str, 0, fmt_1999, 80, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_1000, 0 };
    static icilist io___30 = { 0, str, 0, fmt_2000, 80, 1 };
    static cilist io___41 = { 0, 0, 0, fmt_1020, 0 };
    static icilist io___48 = { 0, str, 0, fmt_2010, 80, 1 };
    static icilist io___49 = { 0, str, 0, fmt_2050, 80, 1 };
    static icilist io___50 = { 0, str, 0, fmt_2060, 80, 1 };


/*     ================================================================== */
/*     This impression of INSERT reads a file produced by  s4pnch. */
/*     It is intended to read files similar to those produced by */
/*     standard MPS systems.  It recognizes SB as an additional key. */
/*     Also, values are extracted from columns 25--36. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4inst. */
/*     01 Aug 2003: snPRNT adopted. */
/*     08 Oct 2003: Current version of s4inst. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    names -= 8;
    --x;
    --bu;
    --bl;
    --hs;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    infbnd = rw[70];
/* definition of an infinite bound */
    istdi = iw[9];
/* Standard Input */
    iinsrt = iw[125];
/* insert file */
    mpserr = iw[106];
/* maximum # errors in MPS data */
    bplus = infbnd * .9;
    s_wsfi(&io___27);
    do_fio(&c__1, (char *)&iinsrt, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    io___28.ciunit = iinsrt;
    s_rsfe(&io___28);
    do_fio(&c__5, id, (ftnlen)4);
    e_rsfe();
    s_wsfi(&io___30);
    do_fio(&c__5, id, (ftnlen)4);
    e_wsfi();
    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
    l1 = *n + 1;
/*     Make logicals basic. */
    i__1 = *nb - *n;
    iload_(&i__1, &c__3, &hs[*n + 1], &c__1);
    ignord = 0;
    nbs = 0;
    *ns = 0;
    notfnd = 0;
    ncard = 0;
    jmark = 1;
    lmark = l1;
    ndum = *n + 100000;
    if (*iobj == 0) {
	jobj = 0;
    } else {
	jobj = *n + *iobj;
    }
/*     Read names until ENDATA */
    i__1 = ndum;
    for (nloop = 1; nloop <= i__1; ++nloop) {
	io___41.ciunit = iinsrt;
	s_rsfe(&io___41);
	do_fio(&c__1, key, (ftnlen)4);
	do_fio(&c__1, name1, (ftnlen)8);
	do_fio(&c__1, name2, (ftnlen)8);
	do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
	e_rsfe();
	if (s_cmp(key, lend, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L310;
	}
/*        Look for  Name1.  It may be a column or a row, */
/*        since a superbasic variable could be either. */
	ncard = nloop;
	s4name_(nb, names + 8, name1, &ncard, &notfnd, &mpserr, &c__1, nb, &
		jmark, &j, &iw[1], leniw, (ftnlen)8, (ftnlen)8);
	if (j <= 0) {
	    goto L300;
	}
	if (hs[j] > 1) {
	    goto L290;
	}
	if (s_cmp(key, lxl, (ftnlen)4, (ftnlen)4) != 0 && s_cmp(key, lxu, (
		ftnlen)4, (ftnlen)4) != 0) {
	    goto L70;
	}
/*        Look for  Name2.  It has to be a row. */
	s4name_(nb, names + 8, name2, &ncard, &notfnd, &mpserr, &l1, nb, &
		lmark, &l, &iw[1], leniw, (ftnlen)8, (ftnlen)8);
	if (l <= 0) {
	    goto L300;
	}
/*        XL, XU (exchange card)  --  make col j basic,  row l nonbasic. */
	if (l == jobj) {
	    goto L290;
	}
	if (hs[l] != 3) {
	    goto L290;
	}
	++nbs;
	hs[j] = 3;
	if (s_cmp(key, lxu, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L50;
	}
	hs[l] = 0;
	if (bl[l] > -bplus) {
	    x[l] = bl[l];
	}
	goto L250;
L50:
	hs[l] = 1;
	if (bu[l] < bplus) {
	    x[l] = bu[l];
	}
	goto L250;
/*        LL, UL, SB  --  only  j  and  xj  are relevant. */
L70:
	if (s_cmp(key, lll, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L100;
	}
	if (s_cmp(key, lul, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L150;
	}
	if (s_cmp(key, lsb, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L200;
	}
	goto L290;
/*        LO or UP */
L100:
	hs[j] = 0;
	goto L250;
L150:
	hs[j] = 1;
	goto L250;
/*        Make superbasic. */
L200:
	hs[j] = 2;
	++(*ns);
/*        Save  x  values. */
L250:
	if (abs(xj) < bplus) {
	    x[j] = xj;
	}
	goto L300;
/*        Card ignored. */
L290:
	++ignord;
	if (ignord <= mpserr) {
	    s_wsfi(&io___48);
	    do_fio(&c__1, (char *)&ncard, (ftnlen)sizeof(integer));
	    do_fio(&c__1, key, (ftnlen)4);
	    do_fio(&c__1, name1, (ftnlen)8);
	    do_fio(&c__1, name2, (ftnlen)8);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
L300:
	;
    }
L310:
    ignord += notfnd;
    s_wsfi(&io___49);
    do_fio(&c__1, (char *)&ncard, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ignord, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    s_wsfi(&io___50);
    do_fio(&c__1, (char *)&nbs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
    if (iinsrt != istdi) {
	al__1.aerr = 0;
	al__1.aunit = iinsrt;
	f_rew(&al__1);
    }
    return 0;
} /* s4inst_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4inst */
/* Subroutine */ int s4load_(integer *n, integer *nb, integer *ns, integer *
	iobj, integer *hs, doublereal *x, char *names, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen names_len)
{
    /* Initialized data */

    static char lbs[4] = " BS ";
    static char lll[4] = " LL ";
    static char lul[4] = " UL ";
    static char lsb[4] = " SB ";
    static char lend[4] = "ENDA";

    /* Format strings */
    static char fmt_1999[] = "(\002 LOAD file to be input from file\002,i4)";
    static char fmt_1000[] = "(14x,2a4,2x,3a4)";
    static char fmt_2000[] = "(\002 NAME\002,10x,2a4,2x,3a4)";
    static char fmt_1020[] = "(a4,a8,12x,e12.5)";
    static char fmt_2010[] = "(\002 XXX  Line\002,i6,\002  ignored:\002,8x,a"
	    "4,a8,2x,a8)";
    static char fmt_2050[] = "(\002 No. of lines read      \002,i6,\002  Lin"
	    "es ignored\002,i6)";
    static char fmt_2060[] = "(\002 No. of basics specified\002,i6,\002  Sup"
	    "erbasics  \002,i6)";

    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_rsfe(cilist *), e_rsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), f_rew(alist *);

    /* Local variables */
    static integer j;
    static char id[4*5];
    static doublereal xj;
    static integer nbs;
    static char key[4], str[80], name__[8];
    static integer jobj, ndum, ncard, jmark, istdi;
    static doublereal bplus;
    static integer nloop;
    extern /* Subroutine */ int s4name_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer iloadb;
    static doublereal infbnd;
    static integer ignord, notfnd, mpserr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___62 = { 0, str, 0, fmt_1999, 80, 1 };
    static cilist io___63 = { 0, 0, 0, fmt_1000, 0 };
    static icilist io___65 = { 0, str, 0, fmt_2000, 80, 1 };
    static cilist io___74 = { 0, 0, 0, fmt_1020, 0 };
    static icilist io___79 = { 0, str, 0, fmt_2010, 80, 1 };
    static icilist io___80 = { 0, str, 0, fmt_2050, 80, 1 };
    static icilist io___81 = { 0, str, 0, fmt_2060, 80, 1 };


/*     ================================================================== */
/*     s4load  inputs a load file, which may contain a full or partial */
/*     list of row and column names and their states and values. */
/*     Valid keys are   BS, LL, UL, SB. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4load. */
/*     01 Aug 2003: snPRNT adopted. */
/*     01 Aug 2003: Current version of s4load. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    names -= 8;
    --x;
    --hs;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    infbnd = rw[70];
/* definition of an infinite bound */
    istdi = iw[9];
/* Standard Input */
    iloadb = iw[122];
/* load file */
    mpserr = iw[106];
/* maximum # errors in MPS data */
    bplus = infbnd * .9;
    s_wsfi(&io___62);
    do_fio(&c__1, (char *)&iloadb, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    io___63.ciunit = iloadb;
    s_rsfe(&io___63);
    do_fio(&c__5, id, (ftnlen)4);
    e_rsfe();
    s_wsfi(&io___65);
    do_fio(&c__5, id, (ftnlen)4);
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    ignord = 0;
    nbs = 0;
    *ns = 0;
    notfnd = 0;
    ncard = 0;
    jmark = 1;
    ndum = *n + 100000;
    if (*iobj == 0) {
	jobj = 0;
    } else {
	jobj = *n + *iobj;
    }
/*     Read names until ENDATA is found. */
    i__1 = ndum;
    for (nloop = 1; nloop <= i__1; ++nloop) {
	io___74.ciunit = iloadb;
	s_rsfe(&io___74);
	do_fio(&c__1, key, (ftnlen)4);
	do_fio(&c__1, name__, (ftnlen)8);
	do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
	e_rsfe();
	if (s_cmp(key, lend, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L310;
	}
	ncard = nloop;
	s4name_(nb, names + 8, name__, &ncard, &notfnd, &mpserr, &c__1, nb, &
		jmark, &j, &iw[1], leniw, (ftnlen)8, (ftnlen)8);
	if (j <= 0) {
	    goto L300;
	}
/*        The name Name belongs to the j-th variable. */
	if (hs[j] > 1) {
	    goto L290;
	}
	if (j == jobj) {
	    goto L90;
	}
	if (s_cmp(key, lbs, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L90;
	}
	if (s_cmp(key, lll, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L100;
	}
	if (s_cmp(key, lul, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L150;
	}
	if (s_cmp(key, lsb, (ftnlen)4, (ftnlen)4) == 0) {
	    goto L200;
	}
	goto L290;
/*        Make basic. */
L90:
	++nbs;
	hs[j] = 3;
	goto L250;
/*        LO or UP. */
L100:
	hs[j] = 0;
	goto L250;
L150:
	hs[j] = 1;
	goto L250;
/*        Make superbasic. */
L200:
	++(*ns);
	hs[j] = 2;
/*        Save  x  values. */
L250:
	if (abs(xj) < bplus) {
	    x[j] = xj;
	}
	goto L300;
/*        Card ignored. */
L290:
	++ignord;
	if (ignord <= mpserr) {
	    s_wsfi(&io___79);
	    do_fio(&c__1, (char *)&ncard, (ftnlen)sizeof(integer));
	    do_fio(&c__1, key, (ftnlen)4);
	    do_fio(&c__1, name__, (ftnlen)8);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
L300:
	;
    }
L310:
    ignord += notfnd;
    s_wsfi(&io___80);
    do_fio(&c__1, (char *)&ncard, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ignord, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    s_wsfi(&io___81);
    do_fio(&c__1, (char *)&nbs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
/*     Make sure the linear Objective is basic. */
    if (*iobj > 0) {
	if (hs[jobj] != 3) {
	    hs[jobj] = 3;
/*           Swap Obj with last basic variable. */
	    for (j = *nb; j >= 1; --j) {
		if (hs[j] == 3) {
		    goto L860;
		}
	    }
L860:
	    hs[j] = 0;
	}
    }
    if (iloadb != istdi) {
	al__1.aerr = 0;
	al__1.aunit = iloadb;
	f_rew(&al__1);
    }
    return 0;
} /* s4load_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4load */
/* Subroutine */ int s4oldb_(integer *iexit, integer *m, integer *n, integer *
	nb, integer *ns, integer *hs, doublereal *bl, doublereal *bu, 
	doublereal *x, integer *iw, integer *leniw, doublereal *rw, integer *
	lenrw)
{
    /* Format strings */
    static char fmt_1999[] = "(\002 OLD BASIS file to be input from file\002"
	    ",i4)";
    static char fmt_1000[] = "(20a4)";
    static char fmt_2000[] = "(1x,20a4)";
    static char fmt_1005[] = "(13a4,2x,i7,3x,i7,4x,i5)";
    static char fmt_2005[] = "(1x,13a4,\002m=\002,i7,\002 n=\002,i7,\002 sb"
	    "=\002,i5)";
    static char fmt_1010[] = "(80i1)";
    static char fmt_1020[] = "(i8,e24.14)";
    static char fmt_2010[] = "(\002 No. of superbasics loaded\002,i7)";

    /* System generated locals */
    integer hs_dim1, i__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_rsfe(cilist *), e_rsfe(void), f_rew(alist *);

    /* Local variables */
    static integer i__, j;
    static char id[4*20];
    static integer js;
    static doublereal xj;
    static char str[81];
    static integer newm, newn, ioldb, istdi;
    static doublereal bplus, infbnd;
    static integer idummy, ndummy;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___87 = { 0, str, 0, fmt_1999, 81, 1 };
    static cilist io___88 = { 0, 0, 0, fmt_1000, 0 };
    static icilist io___90 = { 0, str, 0, fmt_2000, 81, 1 };
    static cilist io___91 = { 0, 0, 0, fmt_1005, 0 };
    static icilist io___95 = { 0, str, 0, fmt_2005, 81, 1 };
    static cilist io___96 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_1020, 0 };
    static icilist io___103 = { 0, str, 0, fmt_2010, 81, 1 };


/*     ================================================================== */
/*     s4oldB  inputs a compact basis file from file  iOldB. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4OldB. */
/*     01 Aug 2003: snExit and snPRNT adopted. */
/*     04 Dec 2004: Current version of s4oldB. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    hs_dim1 = *nb;
    --hs;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    istdi = iw[9];
/* Standard Input */
    ioldb = iw[126];
/* old basis file */
    bplus = infbnd * .9;
    s_wsfi(&io___87);
    do_fio(&c__1, (char *)&ioldb, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)81);
    io___88.ciunit = ioldb;
    s_rsfe(&io___88);
    do_fio(&c__20, id, (ftnlen)4);
    e_rsfe();
    s_wsfi(&io___90);
    do_fio(&c__20, id, (ftnlen)4);
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)81);
    io___91.ciunit = ioldb;
    s_rsfe(&io___91);
    for (i__ = 1; i__ <= 13; ++i__) {
	do_fio(&c__1, id + (i__ - 1 << 2), (ftnlen)4);
    }
    do_fio(&c__1, (char *)&newm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&newn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_rsfe();
    s_wsfi(&io___95);
    for (i__ = 1; i__ <= 13; ++i__) {
	do_fio(&c__1, id + (i__ - 1 << 2), (ftnlen)4);
    }
    do_fio(&c__1, (char *)&newm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&newn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)81);
    if (newm != *m || newn != *n) {
	goto L900;
    }
    io___96.ciunit = ioldb;
    s_rsfe(&io___96);
    i__1 = 1 * hs_dim1;
    do_fio(&i__1, (char *)&hs[1], (ftnlen)sizeof(integer));
    e_rsfe();
/*     Set values for nonbasic variables. */
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	js = hs[j];
	if (js <= 1) {
	    if (js == 0) {
		xj = bl[j];
	    }
	    if (js == 1) {
		xj = bu[j];
	    }
	    if (abs(xj) < bplus) {
		x[j] = xj;
	    }
	}
    }
/*     Load superbasics. */
    *ns = 0;
    ndummy = *m + *n + 10000;
    i__1 = ndummy;
    for (idummy = 1; idummy <= i__1; ++idummy) {
	io___102.ciunit = ioldb;
	s_rsfe(&io___102);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
	e_rsfe();
	if (j <= 0) {
	    goto L310;
	}
	if (j <= *nb) {
	    x[j] = xj;
	    if (hs[j] == 2) {
		++(*ns);
	    }
	}
    }
L310:
    s_wsfi(&io___103);
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)81);
    goto L990;
/*     Error exits. */
/*     ------------------------------------------- */
/*     Incompatible basis file dimensions. */
/*     ------------------------------------------- */
L900:
    *iexit = 92;
/* Basis file dimensions do not match this problem */
L990:
    if (ioldb != istdi) {
	al__1.aerr = 0;
	al__1.aunit = ioldb;
	f_rew(&al__1);
    }
    return 0;
} /* s4oldb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4oldB */
/* Subroutine */ int s4chek_(integer *m, integer *maxs, integer *mbs, integer 
	*n, integer *nb, logical *needb, integer *gothes, integer *ns, 
	integer *iobj, integer *hs, integer *kbs, doublereal *bl, doublereal *
	bu, doublereal *x, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 WARNING:\002,i6,\002 superbasics in hs(*"
	    ");\002,\002 previously nS =\002,i6,\002.  Hessian not saved\002)";
    static char fmt_1100[] = "(\002 WARNING:\002,i7,\002 basics specified"
	    ";\002,\002 preferably should have been\002,i7)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j;
    static doublereal b1, b2;
    static integer jj, js;
    static doublereal xj;
    static char str[80];
    static doublereal bplus;
    static logical setks;
    static integer nbasic;
    static doublereal infbnd;
    static integer nssave;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___112 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___113 = { 0, str, 0, fmt_1100, 80, 1 };


/*     ================================================================== */
/*     s4chek  takes hs and x and checks they contain reasonable values. */
/*     The entries hs(j) = 2 are used to set  nS  and possibly */
/*     the list of superbasic variables kBS(m+1) thru kBS(m+nS). */
/*     Scaling, if any, has taken place by this stage. */

/*     15 Nov 1991: First version based on Minos routine m4chek. */
/*     01 Aug 2003: snPRNT adopted. */
/*     08 Mar 2004: needB, gotHes, setkS allow for Hot start. */
/*     08 Mar 2004: Current version of s4chek. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --kbs;
    --x;
    --bu;
    --bl;
    --hs;
    --iw;
    --rw;

    /* Function Body */
    infbnd = rw[70];
/* definition of an infinite bound */
    setks = *needb || *gothes <= 0;
/*     Make sure hs(j) = 0, 1, 2 or 3 only. */
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	js = hs[j];
	if (js < 0) {
	    hs[j] = 0;
	}
	if (js >= 4) {
	    hs[j] = js - 4;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Make sure the Objective is basic and free. */
/*     Then count the basics and superbasics, making sure they don't */
/*     exceed m and maxS respectively.  Also, set nS and possibly */
/*     kBS(m+1) thru kBS(m+ns) to define the list of superbasics. */
/*     Mar 1988: Loop 100 now goes backwards to make sure we grab Obj. */
/*     Apr 1992: Backwards seems a bit silly in the documentation. */
/*               We now go forward through the slacks, */
/*               then forward through the columns. */
/*     ------------------------------------------------------------------ */
L100:
    nbasic = 0;
    nssave = *ns;
    *ns = 0;
    if (*iobj > 0) {
	hs[*n + *iobj] = 3;
	bl[*n + *iobj] = -infbnd;
	bu[*n + *iobj] = infbnd;
    }
/*     If too many basics or superbasics, make them nonbasic. */
/*     Do slacks first to make sure we grab the objective slack. */
    j = *n;
    i__1 = *nb;
    for (jj = 1; jj <= i__1; ++jj) {
	++j;
	if (j > *nb) {
	    j = 1;
	}
	js = hs[j];
	if (js == 2) {
	    ++(*ns);
	    if (*ns <= *maxs) {
		if (setks) {
		    kbs[*m + *ns] = j;
		}
	    } else {
		hs[j] = 0;
	    }
	} else if (js == 3) {
	    ++nbasic;
	    if (nbasic > *m) {
		hs[j] = 0;
	    }
	}
    }
/*     Proceed if the superbasic kBS were reset, */
/*     or if nS seems to agree with the superbasics in hs(*). */
/*     Otherwise, give up trying to save the reduced Hessian, */
/*     and reset the superbasic kBS after all. */
    if (setks) {
/*        ok */
    } else if (*ns != nssave) {
	setks = TRUE_;
	*gothes = 0;
	s_wsfi(&io___112);
	do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nssave, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
	goto L100;
    }
/*     Check the number of basics. */
    *ns = min(*ns,*maxs);
    if (nbasic != *m) {
	s_wsfi(&io___113);
	do_fio(&c__1, (char *)&nbasic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    }
/*     ------------------------------------------------------------------ */
/*     Set each nonbasic x(j) to be exactly on its */
/*     nearest bound if it is within tolb of that bound. */
/*     ------------------------------------------------------------------ */
    bplus = infbnd * .1;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	xj = x[j];
	if (abs(xj) >= bplus) {
	    xj = 0.;
	}
	if (hs[j] <= 1) {
	    b1 = bl[j];
	    b2 = bu[j];
	    xj = max(xj,b1);
	    xj = min(xj,b2);
	    if (xj - b1 > b2 - xj) {
		b1 = b2;
	    }
	    if ((d__1 = xj - b1, abs(d__1)) <= 1e-4) {
		xj = b1;
	    }
	    hs[j] = 0;
	    if (xj > bl[j]) {
		hs[j] = 1;
	    }
	}
	x[j] = xj;
    }
    return 0;
} /* s4chek_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4chek */
/* Subroutine */ int s4chkp_(integer *errors, char *cpointr, integer *ipointr,
	 integer *iw, integer *leniw, ftnlen cpointr_len)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 XXX  Pointer out of range:  \002,a6,\002"
	    " = \002,i6)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___119 = { 0, str, 0, fmt_9999, 80, 1 };


/*     ================================================================== */
/*     s4chkP checks the value of a pointer retrieved from workspace */
/*     Error messages are listed on the standard output */

/*     03 Sep 2006: First version of s4chkP. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (*ipointr < 0 || *ipointr > *leniw) {
	++(*errors);
	s_wsfi(&io___119);
	do_fio(&c__1, cpointr, (ftnlen)6);
	do_fio(&c__1, (char *)&(*ipointr), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__5, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s4chkp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4chkP */
/* Subroutine */ int s4dump_(integer *idump, integer *m, integer *n, integer *
	nb, integer *nname, integer *hs, doublereal *x, char *names, char *cw,
	 integer *lencw, integer *iw, integer *leniw, ftnlen names_len, 
	ftnlen cw_len)
{
    /* Initialized data */

    static char key[4*4] = " LL " " UL " " SB " " BS ";

    /* Format strings */
    static char fmt_2000[] = "(\002NAME\002,10x,a8,2x,\002   DUMP/LOAD\002)";
    static char fmt_2100[] = "(a4,a8,12x,1p,e12.5)";
    static char fmt_2200[] = "(\002ENDATA\002)";
    static char fmt_3000[] = "(\002 DUMP file saved on file\002,i4)";

    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsfi(icilist *), e_wsfi(void), f_rew(alist *);

    /* Local variables */
    static integer j, k;
    static char id[8], str[40];
    extern /* Subroutine */ int s4id_(integer *, integer *, integer *, 
	    integer *, char *, char *, ftnlen, ftnlen);
    static char mprob[8];
    static integer iprint;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___123 = { 0, 0, 0, fmt_2000, 0 };
    static cilist io___127 = { 0, 0, 0, fmt_2100, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_2200, 0 };
    static icilist io___130 = { 0, str, 0, fmt_3000, 40, 1 };


/*     ================================================================== */
/*     s4dump outputs basis names in a format compatible with s4load. */
/*     This file is normally easier to modify than a punch file. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4dump. */
/*     01 Aug 2003: snPRNT adopted. */
/*     17 Jul 2005: Dump default names when nName = 1 */
/*     ================================================================== */
    /* Parameter adjustments */
    --x;
    --hs;
    names -= 8;
    cw -= 8;
    --iw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    iprint = iw[12];
/* Print file */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Problem name */
    io___123.ciunit = *idump;
    s_wsfe(&io___123);
    do_fio(&c__1, mprob, (ftnlen)8);
    e_wsfe();
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	s4id_(&j, n, nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	k = hs[j] + 1;
	io___127.ciunit = *idump;
	s_wsfe(&io___127);
	do_fio(&c__1, key + (k - 1 << 2), (ftnlen)4);
	do_fio(&c__1, id, (ftnlen)8);
	do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    io___128.ciunit = *idump;
    s_wsfe(&io___128);
    e_wsfe();
    s_wsfi(&io___130);
    do_fio(&c__1, (char *)&(*idump), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)40);
    if (*idump != iprint) {
	al__1.aerr = 0;
	al__1.aunit = *idump;
	f_rew(&al__1);
    }
    return 0;
} /* s4dump_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4dump */
/* Subroutine */ int s4newb_(integer *job, integer *inewb, integer *minimz, 
	integer *m, integer *n, integer *nb, integer *ns, integer *mbs, 
	integer *itn, integer *ninf, doublereal *sinf, doublereal *f, integer 
	*kbs, integer *hs, doublereal *ascale, doublereal *bl, doublereal *bu,
	 doublereal *x, doublereal *xbs, char *istate, char *cw, integer *
	lencw, integer *iw, integer *leniw, ftnlen istate_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(a8,\002  ITN\002,i8,4x,3a4,\002  NINF\002,i7"
	    ",\002      OBJ\002,1p,e21.12)";
    static char fmt_1005[] = "(\002OBJ=\002,a8,\002 RHS=\002,a8,\002 RNG="
	    "\002,a8,\002 BND=\002,a8,\002 M=\002,i7,\002 N=\002,i7,\002 SB"
	    "=\002,i5)";
    static char fmt_1020[] = "(i8,1p,e24.14)";
    static char fmt_1030[] = "(\002 NEW BASIS file saved on file\002,i4,\002"
	    "    itn =\002,i7)";

    /* System generated locals */
    integer i__1, i__2;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe(void), f_rew(alist *), s_wsfi(
	    icilist *), e_wsfi(void);

    /* Local variables */
    static integer i__, j, k, j1, j2, js;
    static doublereal xj, obj;
    static char str[60], mbnd[8], mobj[8], mrng[8], mrhs[8];
    static integer nnjac;
    static char mprob[8];
    static logical scaled;
    static integer buffer[80], lvlscl, iprint;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___141 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___142 = { 0, 0, 0, fmt_1005, 0 };
    static cilist io___150 = { 0, 0, 0, "(80i1)", 0 };
    static cilist io___152 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___153 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___154 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_1020, 0 };
    static icilist io___157 = { 0, str, 0, fmt_1030, 60, 1 };


/*     ================================================================== */
/*     s4newB  saves a compact basis on file iNewB.  Called from S5QP. */
/*     job = Freq, the save is a periodic one due to the save frequency. */
/*     job = Wrap, S5solv has just finished the current problem. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4newb. */
/*     16 Jan 2003: First line had "or" instead of "and". */
/*     01 Aug 2003: snPRNT adopted. */
/*     01 Apr 2005: Negative hs values converted in-situ. */
/*     17 Jul 2005: Blanks printed for undefined MPS names. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --ascale;
    --hs;
    --xbs;
    --kbs;
    istate -= 4;
    cw -= 8;
    --iw;

    /* Function Body */
    if (*job != 0 && *job != 1) {
	return 0;
    }
    iprint = iw[12];
/* Print file */
    nnjac = iw[21];
/* # nonlinear Jacobian variables */
    lvlscl = iw[75];
/* scale option */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Problem name */
    s_copy(mobj, cw + 416, (ftnlen)8, (ftnlen)8);
/* Objective name */
    s_copy(mrhs, cw + 424, (ftnlen)8, (ftnlen)8);
/* Right-hand side name */
    s_copy(mrng, cw + 432, (ftnlen)8, (ftnlen)8);
/* Range name */
    s_copy(mbnd, cw + 440, (ftnlen)8, (ftnlen)8);
/* Bnd section name */
    scaled = lvlscl > 0;
    if (s_cmp(mprob, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mprob, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mobj, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mobj, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mrhs, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mrhs, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mrng, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mrng, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(mbnd, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(mbnd, "        ", (ftnlen)8, (ftnlen)8);
    }
    if (*ninf == 0) {
	obj = *minimz * *f;
    } else {
	obj = *sinf;
    }
/*     Output header cards and the state vector. */
    io___141.ciunit = *inewb;
    s_wsfe(&io___141);
    do_fio(&c__1, mprob, (ftnlen)8);
    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
    do_fio(&c__3, istate + 4, (ftnlen)4);
    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&obj, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___142.ciunit = *inewb;
    s_wsfe(&io___142);
    do_fio(&c__1, mobj, (ftnlen)8);
    do_fio(&c__1, mrhs, (ftnlen)8);
    do_fio(&c__1, mrng, (ftnlen)8);
    do_fio(&c__1, mbnd, (ftnlen)8);
    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_wsfe();
/*     write(iNewB, 1010) hs */
    j2 = 0;
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; i__ += 80) {
	j1 = j2 + 1;
	j2 += 80;
	j2 = min(j2,*nb);
	k = 0;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    js = hs[j];
	    if (js >= 4 || js == -1) {
		js = 0;
	    }
	    ++k;
	    buffer[k - 1] = js;
	}
	io___150.ciunit = *inewb;
	s_wsfe(&io___150);
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&buffer[j - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
/*     Output the superbasic variables. */
    i__1 = *m + *ns;
    for (k = *m + 1; k <= i__1; ++k) {
	j = kbs[k];
	xj = xbs[k];
	if (scaled) {
	    xj *= ascale[j];
	}
	io___152.ciunit = *inewb;
	s_wsfe(&io___152);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     Output the values of all other (non-SB) nonlinear variables. */
    i__1 = nnjac;
    for (j = 1; j <= i__1; ++j) {
	if (hs[j] != 2) {
	    xj = x[j];
	    if (scaled) {
		xj *= ascale[j];
	    }
	    io___153.ciunit = *inewb;
	    s_wsfe(&io___153);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/*     Output nonbasic variables that are not at a bound. */
    i__1 = *nb;
    for (j = nnjac + 1; j <= i__1; ++j) {
	js = hs[j];
	if (js <= 1 || js >= 4 || js == -1) {
	    xj = x[j];
	    if (xj != bl[j]) {
		if (xj != bu[j]) {
		    if (scaled) {
			xj *= ascale[j];
		    }
		    io___154.ciunit = *inewb;
		    s_wsfe(&io___154);
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }
/*     Terminate the list with a zero. */
    j = 0;
    io___155.ciunit = *inewb;
    s_wsfe(&io___155);
    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    e_wsfe();
    if (*inewb != iprint) {
	al__1.aerr = 0;
	al__1.aunit = *inewb;
	f_rew(&al__1);
    }
    s_wsfi(&io___157);
    do_fio(&c__1, (char *)&(*inewb), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)60);
    return 0;
/* L1010: */
} /* s4newb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4newB */
/* Subroutine */ int s4pnch_(integer *ipnch, integer *m, integer *n, integer *
	nb, integer *nname, integer *hs, doublereal *bl, doublereal *bu, 
	doublereal *x, char *names, char *cw, integer *lencw, integer *iw, 
	integer *leniw, ftnlen names_len, ftnlen cw_len)
{
    /* Initialized data */

    static char key[4*5] = " LL " " UL " " SB " " XL " " XU ";
    static char lblank[8] = "        ";

    /* Format strings */
    static char fmt_2000[] = "(\002NAME\002,10x,a8,2x,\002PUNCH/INSERT\002)";
    static char fmt_2100[] = "(a4,a8,2x,a8,2x,1p,e12.5)";
    static char fmt_2200[] = "(\002ENDATA\002)";
    static char fmt_3000[] = "(\002 PUNCH file saved on file\002,i4)";

    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsfi(icilist *), e_wsfi(void), f_rew(alist *);

    /* Local variables */
    static integer j, k, jn;
    static char str[60];
    extern /* Subroutine */ int s4id_(integer *, integer *, integer *, 
	    integer *, char *, char *, ftnlen, ftnlen);
    static char name__[8];
    static integer irow;
    static char colnm[8], mprob[8];
    static integer iprint;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___162 = { 0, 0, 0, fmt_2000, 0 };
    static cilist io___169 = { 0, 0, 0, fmt_2100, 0 };
    static cilist io___170 = { 0, 0, 0, fmt_2100, 0 };
    static cilist io___171 = { 0, 0, 0, fmt_2100, 0 };
    static cilist io___172 = { 0, 0, 0, fmt_2200, 0 };
    static icilist io___174 = { 0, str, 0, fmt_3000, 60, 1 };


/*     ================================================================== */
/*     s4pnch  outputs a PUNCH file (list of basis names, states and */
/*     values) in a format that is compatible with MPS/360. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4pnch. */
/*     16 Jan 2003: Allow for missing names. */
/*     01 Aug 2003: snPRNT adopted. */
/*     01 Aug 2003: Current version of s4pnch. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --hs;
    names -= 8;
    cw -= 8;
    --iw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    iprint = iw[12];
/* Print file */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Problem name */
    io___162.ciunit = *ipnch;
    s_wsfe(&io___162);
    do_fio(&c__1, mprob, (ftnlen)8);
    e_wsfe();
    irow = *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jn = j;
	s4id_(&jn, n, nb, nname, names + 8, colnm, (ftnlen)8, (ftnlen)8);
	k = hs[j];
	if (k == 3) {
/*           Basics -- find the next row that isn't basic. */
L300:
	    ++irow;
	    if (irow <= *nb) {
		k = hs[irow];
		if (k == 3) {
		    goto L300;
		}
		s4id_(&irow, n, nb, nname, names + 8, name__, (ftnlen)8, (
			ftnlen)8);
		if (k == 2) {
		    k = 0;
		}
		io___169.ciunit = *ipnch;
		s_wsfe(&io___169);
		do_fio(&c__1, key + (k + 3 << 2), (ftnlen)4);
		do_fio(&c__1, colnm, (ftnlen)8);
		do_fio(&c__1, name__, (ftnlen)8);
		do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else {
/*           Skip nonbasic variables with zero lower bounds. */
	    if (k <= 1) {
		if (bl[j] == 0. && x[j] == 0.) {
		    goto L500;
		}
	    }
	    io___170.ciunit = *ipnch;
	    s_wsfe(&io___170);
	    do_fio(&c__1, key + (k << 2), (ftnlen)4);
	    do_fio(&c__1, colnm, (ftnlen)8);
	    do_fio(&c__1, lblank, (ftnlen)8);
	    do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
L500:
	;
    }
/*     Output superbasic slacks. */
    i__1 = *nb;
    for (j = *n + 1; j <= i__1; ++j) {
	if (hs[j] == 2) {
	    jn = j;
	    s4id_(&jn, n, nb, nname, names + 8, name__, (ftnlen)8, (ftnlen)8);
	    io___171.ciunit = *ipnch;
	    s_wsfe(&io___171);
	    do_fio(&c__1, key + 8, (ftnlen)4);
	    do_fio(&c__1, name__, (ftnlen)8);
	    do_fio(&c__1, lblank, (ftnlen)8);
	    do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    io___172.ciunit = *ipnch;
    s_wsfe(&io___172);
    e_wsfe();
    s_wsfi(&io___174);
    do_fio(&c__1, (char *)&(*ipnch), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)60);
    if (*ipnch != iprint) {
	al__1.aerr = 0;
	al__1.aunit = *ipnch;
	f_rew(&al__1);
    }
    return 0;
} /* s4pnch_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4pnch */
/* Subroutine */ int s4rept_(logical *ondisk, integer *m, integer *n, integer 
	*nb, integer *nname, integer *nncon0, integer *nncon, integer *nnobj0,
	 integer *nnobj, integer *ns, integer *ne, integer *nloca, integer *
	loca, integer *inda, doublereal *acol, integer *hs, doublereal *
	ascale, doublereal *bl, doublereal *bu, doublereal *gobj, doublereal *
	pi, doublereal *x, doublereal *fx, char *names, char *istate, integer 
	*iw, integer *leniw, ftnlen names_len, ftnlen istate_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 XXX Report file requested.  s4rept does "
	    "nothing.\002)";

    /* Builtin functions */
    integer s_wsfi(icilist *), e_wsfi(void);

    /* Local variables */
    static char str[80];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___176 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s4rept  has the same parameter list as s4soln, the routine that */
/*     prints the solution.  It will be called if the SPECS file */
/*     specifies  REPORT file  n  for some positive value of  n. */

/*     pi contains the unscaled dual solution. */
/*     x contains the unscaled primal solution.  There are n + m = nb */
/*        values (n structural variables and m slacks, in that order). */
/*     y  contains the true slack values for nonlinear constraints */
/*        in its first nnCon components (computed by s8nslk). */

/*     This version of s4rept does nothing.    Added for PILOT, Oct 1985. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4rept. */
/*     26 Mar 2000: Updated for SNOPT 6.1. */
/*     01 Aug 2003: snPRNT adopted. */
/*     01 Aug 2003: Current version of s4rept. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --x;
    --bu;
    --bl;
    --ascale;
    --hs;
    names -= 8;
    --fx;
    --gobj;
    --acol;
    --inda;
    --loca;
    istate -= 4;
    --iw;

    /* Function Body */
    s_wsfi(&io___176);
    e_wsfi();
    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    return 0;
} /* s4rept_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4rept */
/* Subroutine */ int s4savb_(integer *iexit, integer *task, integer *minimz, 
	integer *m, integer *n, integer *nb, integer *nkx, integer *nncon0, 
	integer *nncon, integer *nnobj0, integer *nnobj, integer *nname, 
	integer *ns, integer *itn, integer *ninf, doublereal *sinf, 
	doublereal *wtinf, doublereal *vimax, integer *iobj, doublereal *
	sclobj, doublereal *objtru, doublereal *pnorm1, doublereal *pnorm2, 
	doublereal *pinorm, doublereal *xnorm, integer *ne, integer *nlocj, 
	integer *locj, integer *indj, doublereal *jcol, integer *kx, integer *
	hestat, integer *hs, doublereal *ascale, doublereal *bl, doublereal *
	bu, doublereal *fx, doublereal *gobj, char *names, doublereal *pi, 
	doublereal *rc, doublereal *x, char *cw, integer *lencw, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw, ftnlen names_len, 
	ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1010[] = "(\002 Max x       (scaled)\002,i9,1p,e8.1,2x"
	    ",\002 Max pi      (scaled)\002,i9,e8.1)";
    static char fmt_1020[] = "(\002 Max x               \002,i9,1p,e8.1,2x"
	    ",\002 Max pi              \002,i9,e8.1)";
    static char fmt_1030[] = "(\002 Max Prim inf(scaled)\002,i9,1p,e8.1,2x"
	    ",\002 Max Dual inf(scaled)\002,i9,e8.1)";
    static char fmt_1040[] = "(\002 Max Primal infeas   \002,i9,1p,e8.1,2x"
	    ",\002 Max Dual infeas     \002,i9,e8.1)";
    static char fmt_1080[] = "(\002 Nonlinear constraint violn\002,1p,e11.1)";
    static char fmt_1200[] = "(\002 Solution printed on file\002,i4)";
    static char fmt_1300[] = "(\002 Solution not printed\002)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer k;
    static char str[80];
    static doublereal eps0, binf, dinf;
    static integer maxp, maxx;
    static logical prnt;
    static doublereal tolx, binf1, dinf1;
    extern /* Subroutine */ int s2rca_(logical *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer maxp1, maxx1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer jbinf;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jdinf;
    extern /* Subroutine */ int dddiv_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ipnch, idump;
    static doublereal pimax;
    static integer isoln, maxvi, jbinf1, jdinf1;
    extern /* Subroutine */ int s2binf_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *), s2dinf_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), s2scla_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal pimax1;
    extern /* Subroutine */ int s4pnch_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, char *, integer *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern integer s2varn_(integer *, integer *, integer *);
    extern /* Subroutine */ int s4dump_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, char *, char *, 
	    integer *, integer *, integer *, ftnlen, ftnlen), s2xmat_(integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *), s2vmax_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    , s4rept_(logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, char *, integer 
	    *, integer *, ftnlen, ftnlen), s4stat_(integer *, char *, ftnlen),
	     s4soln_(logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, char *, integer 
	    *, integer *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static doublereal xnorm1;
    static logical feasbl;
    static doublereal infbnd;
    extern integer idamax_(integer *, doublereal *, integer *);
    static char istate[4*3];
    static integer lvlscl, iprint, ireprt, lprsol;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___206 = { 0, str, 0, fmt_1010, 80, 1 };
    static icilist io___207 = { 0, str, 0, fmt_1020, 80, 1 };
    static icilist io___208 = { 0, str, 0, fmt_1030, 80, 1 };
    static icilist io___209 = { 0, str, 0, fmt_1040, 80, 1 };
    static icilist io___211 = { 0, str, 0, fmt_1080, 80, 1 };
    static icilist io___213 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___214 = { 0, str, 0, fmt_1300, 80, 1 };


/*     ================================================================== */
/*     s4savB  saves basis files  and/or  prints the solution. */
/*     It is called twice at the end of s5solv */
/*              and twice at the end of s8solv. */

/*     If Task = SaveB, the problem is first unscaled, then from 0 to 4 */
/*     files are saved (PUNCH file, DUMP file, SOLUTION file, */
/*     REPORT file, in that order). */
/*     A NEW BASIS file, if any, will already have been saved by s8SQP. */
/*     A call with Task = SaveB must precede a call with Task = PrintS. */

/*     If Task = PrintS, the solution is printed under the control of */
/*     lprSol (which is set by the Solution keyword in the SPECS file). */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4savb. */
/*     19 Feb 1994: Use s4rc to compute reduced costs. */
/*     05 Apr 1996: s2rcA called to get the reduced costs (as in */
/*                  Minos 5.5). Maximum primal and dual infeasibilities */
/*                  computed and printed here. */
/*     14 Jul 1997: Thread-safe version. */
/*     26 Mar 2000: Updated for SNOPT 6.1. */
/*     16 Jan 2003: Reinstated calls to s4dump, s4pnch. */
/*     27 Jul 2003: Print max elements of x and pi (and which ones). */
/*                  This is more helpful than piNorm >= 1. */
/*     01 Aug 2003: snPRNT adopted. */
/*     09 Mar 2004: s4soln now deals with the UNSCALED solution, */
/*                  so we no longer have to save the scaled pinorm. */
/*     17 Jun 2004: Save biggest primal and dual infeasibilities */
/*                  before and after scaling, so GAMS can use them. */
/*     29 Apr 2011: Fixed unscaled primal and dual infeasibilities. */
/*     27 Dec 2013: s4savB unscales gObj even for Scale option 1. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --kx;
    --fx;
    --gobj;
    names -= 8;
    --jcol;
    --indj;
    --locj;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    iprint = iw[12];
/* Print file */
    idump = iw[121];
/* dump file */
    ipnch = iw[127];
/* punch file */
    ireprt = iw[130];
/* Report file */
    isoln = iw[131];
/* Solution file */
    lvlscl = iw[75];
/* scale option */
    lprsol = iw[84];
/* > 0    =>  print the solution */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    tolx = rw[56];
/* Minor feasibility tolerance. */
    infbnd = rw[70];
/* definition of an infinite bound */
    feasbl = *ninf == 0;
    k = min(9,*iexit) / 10 + 1;
    s4stat_(&k, istate, (ftnlen)4);
    if (*task == 0) {
/*        --------------------------------------------------------------- */
/*        Compute rc and unscale Jcol, bl, bu, g, pi, x, xNorm */
/*        and piNorm (but s4soln uses scaled piNorm, so save it). */
/*        Then save basis files. */
/*        --------------------------------------------------------------- */
/*        Compute reduced costs rc(*) for all columns and rows. */
/*        Find the maximum bound and dual infeasibilities. */
	s2rca_(&feasbl, &tolx, iobj, minimz, wtinf, m, n, nb, nnobj0, nnobj, 
		ne, nlocj, &locj[1], &indj[1], &jcol[1], &hestat[1], &hs[1], &
		bl[1], &bu[1], &gobj[1], &pi[1], &rc[1], &x[1]);
	s2binf_(nb, &bl[1], &bu[1], &x[1], &binf, &jbinf);
	s2dinf_(n, nb, iobj, &eps0, &bl[1], &bu[1], &rc[1], &x[1], &dinf, &
		jdinf);
	jbinf = s2varn_(&jbinf, leniw, &iw[1]);
	jdinf = s2varn_(&jdinf, leniw, &iw[1]);
	binf1 = binf;
	dinf1 = dinf;
	jbinf1 = jbinf;
	jdinf1 = jdinf;
	rw[427] = binf;
/* Largest bound infeasibility (  scaled) */
	rw[428] = dinf;
/* Largest dual  infeasibility (  scaled) */
	iw[427] = jbinf;
	iw[428] = jdinf;
	maxx = idamax_(n, &x[1], &c__1);
	maxp = idamax_(m, &pi[1], &c__1);
	*xnorm = (d__1 = x[maxx], abs(d__1));
	pimax = (d__1 = pi[maxp], abs(d__1));
	*pinorm = max(pimax,1.);
	maxx1 = maxx;
	maxp1 = maxp;
	xnorm1 = *xnorm;
	pimax1 = pimax;
	*pnorm1 = *pinorm;
/*        Unscale a, bl, bu, pi, x, rc, Fx, gObj and xNorm, piNorm. */
/*        (Previously, s4soln used the scaled piNorm, but no more.) */
	if (lvlscl > 0) {
	    s2scla_(&c__1, m, n, nb, iobj, &infbnd, sclobj, ne, nlocj, &locj[
		    1], &indj[1], &jcol[1], &ascale[1], &bl[1], &bu[1], &pi[1]
		    , &x[1]);
	    dddiv_(nb, &ascale[1], &c__1, &rc[1], &c__1);
/* 27 Dec 2013: We should unscale even for Scale option 1 */
/* if (lvlScl .eq. 2) then */
	    if (*nncon > 0) {
		ddscl_(nncon, &ascale[*n + 1], &c__1, &fx[1], &c__1);
	    }
	    if (*nnobj > 0) {
		dddiv_(nnobj, &ascale[1], &c__1, &gobj[1], &c__1);
	    }
/* end if */
	    maxx = idamax_(n, &x[1], &c__1);
	    maxp = idamax_(m, &pi[1], &c__1);
	    *xnorm = (d__1 = x[maxx], abs(d__1));
	    pimax = (d__1 = pi[maxp], abs(d__1));
	    *pinorm = max(pimax,1.);
	    s2binf_(nb, &bl[1], &bu[1], &x[1], &binf, &jbinf);
	    s2dinf_(n, nb, iobj, &eps0, &bl[1], &bu[1], &rc[1], &x[1], &dinf, 
		    &jdinf);
	    jbinf = s2varn_(&jbinf, leniw, &iw[1]);
	    jdinf = s2varn_(&jdinf, leniw, &iw[1]);
	}
	*pnorm2 = *pinorm;
	rw[422] = *pinorm;
/* Lagrange multiplier norm */
	rw[429] = binf;
/* Largest bound infeasibility (unscaled) */
	rw[430] = dinf;
/* Largest dual  infeasibility (unscaled) */
	iw[429] = jbinf;
	iw[430] = jdinf;
/*        --------------------------------------------------------------- */
/*        Print various scaled and unscaled norms. */
/*        xNorm1, pNorm1  are  scaled (if scaling was used) */
/*                pNorm2  is unscaled */
/*                piNorm  is unscaled */
/*        --------------------------------------------------------------- */
	if (lvlscl > 0) {
	    s_wsfi(&io___206);
	    do_fio(&c__1, (char *)&maxx1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xnorm1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&maxp1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pimax1, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	}
	s_wsfi(&io___207);
	do_fio(&c__1, (char *)&maxx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*xnorm), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&maxp, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pimax, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	if (lvlscl > 0) {
	    s_wsfi(&io___208);
	    do_fio(&c__1, (char *)&jbinf1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&binf1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&jdinf1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dinf1, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	}
	s_wsfi(&io___209);
	do_fio(&c__1, (char *)&jbinf, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&binf, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&jdinf, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dinf, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
/*        Change the sign of pi and rc if feasible and maximizing. */
	if (*ninf == 0 && *minimz < 0) {
	    dscal_(m, &c_b298, &pi[1], &c__1);
	    dscal_(nb, &c_b298, &rc[1], &c__1);
	}
/*        Compute nonlinear constraint infeasibilities (violations). */
	if (*nncon > 0) {
	    s2vmax_(n, nncon, &maxvi, vimax, &bl[1], &bu[1], &fx[1]);
	    s_wsfi(&io___211);
	    do_fio(&c__1, (char *)&(*vimax), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)80);
	}
/*        Output Punch, Dump, Solution and/or Report files. */
	if (ipnch > 0) {
	    s4pnch_(&ipnch, m, n, nb, nname, &hs[1], &bl[1], &bu[1], &x[1], 
		    names + 8, cw + 8, lencw, &iw[1], leniw, (ftnlen)8, (
		    ftnlen)8);
	}
	if (idump > 0) {
	    s4dump_(&idump, m, n, nb, nname, &hs[1], &x[1], names + 8, cw + 8,
		     lencw, &iw[1], leniw, (ftnlen)8, (ftnlen)8);
	}
/* 09 Mar 2004: No longer worry about scaled solution. */
/* piNorm = pNorm1 */
	if (isoln > 0) {
	    s4soln_(&c_true, minimz, m, n, nb, nkx, nname, nncon0, nncon, 
		    nnobj0, nnobj, ns, iobj, itn, ninf, sinf, objtru, pinorm, 
		    names + 8, ne, nlocj, &locj[1], &indj[1], &jcol[1], &kx[1]
		    , &hs[1], &ascale[1], &bl[1], &bu[1], &gobj[1], &pi[1], &
		    rc[1], &x[1], &fx[1], istate, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)4, (ftnlen)8);
	}
/* 09 Oct 2000: Export A, B or (B S) to Report file 91, 92 or 93 */
	if (ireprt >= 91 && ireprt <= 93) {
	    s2xmat_(&ireprt, n, nb, ne, nlocj, &locj[1], &indj[1], &jcol[1], &
		    hs[1]);
	} else if (ireprt > 0) {
	    s4rept_(&c_true, m, n, nb, nname, nncon0, nncon, nnobj0, nnobj, 
		    ns, ne, nlocj, &locj[1], &indj[1], &jcol[1], &hs[1], &
		    ascale[1], &bl[1], &bu[1], &gobj[1], &pi[1], &x[1], &fx[1]
		    , names + 8, istate, &iw[1], leniw, (ftnlen)8, (ftnlen)4);
	}
/* 09 Mar 2004: No longer worry about scaled solution. */
/* piNorm = pNorm2 */
    } else if (*task == 1) {
/*        --------------------------------------------------------------- */
/*        Print solution if requested. */

/*        lprSol = 0   means   no */
/*               = 1   means   if optimal, infeasible or unbounded */
/*               = 2   means   yes */
/*               = 3   means   if error condition */
/*        --------------------------------------------------------------- */
	prnt = iprint > 0 && lprsol > 0;
	if (lprsol == 1 && *iexit > 2 || lprsol == 3 && *iexit <= 2) {
	    prnt = FALSE_;
	}
	if (prnt) {
/* 09 Mar 2004: No longer worry about scaled solution. */
/* piNorm = pNorm1 */
	    s4soln_(&c_false, minimz, m, n, nb, nkx, nname, nncon0, nncon, 
		    nnobj0, nnobj, ns, iobj, itn, ninf, sinf, objtru, pinorm, 
		    names + 8, ne, nlocj, &locj[1], &indj[1], &jcol[1], &kx[1]
		    , &hs[1], &ascale[1], &bl[1], &bu[1], &gobj[1], &pi[1], &
		    rc[1], &x[1], &fx[1], istate, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)4, (ftnlen)8);
/* 09 Mar 2004: No longer worry about scaled solution. */
/* piNorm = pNorm2 */
	    s_wsfi(&io___213);
	    do_fio(&c__1, (char *)&iprint, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__12, str, &iw[1], leniw, (ftnlen)80);
	} else {
	    s_wsfi(&io___214);
	    e_wsfi();
	    snprnt_(&c__12, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    return 0;
} /* s4savb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4savB */
/* Subroutine */ int s4soln_(logical *ondisk, integer *minimz, integer *m, 
	integer *n, integer *nb, integer *nkx, integer *nname, integer *
	nncon0, integer *nncon, integer *nnobj0, integer *nnobj, integer *ns, 
	integer *iobj, integer *itn, integer *ninf, doublereal *sinf, 
	doublereal *objtru, doublereal *pinorm, char *names, integer *ne, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *kx, integer *hs, doublereal *ascale, doublereal *bl, 
	doublereal *bu, doublereal *gobj, doublereal *pi, doublereal *rc, 
	doublereal *x, doublereal *fx, char *istate, char *cw, integer *lencw,
	 integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	names_len, ftnlen istate_len, ftnlen cw_len)
{
    /* Initialized data */

    static char objtyp[8*3] = "Max     " "Feas    " "Min     ";

    /* Format strings */
    static char fmt_1002[] = "(\002 Name\002,11x,a8,16x,\002 Objective Valu"
	    "e\002,1p,e22.10)";
    static char fmt_1000[] = "(\002 Name\002,11x,a8,16x,\002 Infeasibilitie"
	    "s\002,i7,1p,e16.4)";
    static char fmt_1004[] = "(\002 Status\002,9x,3a4,12x,\002 Iteration\002"
	    ",i7,\002    Superbasics\002,i6)";
    static char fmt_1005[] = "(\002 Objective\002,6x,a8,\002 (\002,a3,\002"
	    ")\002)";
    static char fmt_1006[] = "(\002 RHS\002,12x,a8)";
    static char fmt_1007[] = "(\002 Ranges\002,9x,a8)";
    static char fmt_1008[] = "(\002 Bounds\002,9x,a8)";
    static char fmt_1010[] = "(/\002 Section 1 - Rows\002//\002  Number  ..."
	    "Row.. State  ...Activity...  Slack Activity\002,\002  ..Lower Li"
	    "mit.  ..Upper Limit.  .Dual Activity    ..i\002/)";
    static char fmt_1020[] = "(\002 Section 2 - Columns\002//\002  Number  ."
	    "Column. State  ...Activity...  .Obj Gradient.\002,\002  ..Lower "
	    "Limit.  ..Upper Limit.  Reduced Gradnt    m+j\002/)";
    static char fmt_1400[] = "(\002 SOLUTION file saved on file\002,i4)";

    /* System generated locals */
    integer i__1, i__2;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_wsfe(cilist *), e_wsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), f_rew(alist *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal b1, b2, d1, d2;
    static char id[8];
    static doublereal cj, dj;
    static integer in, ir, jn;
    static doublereal xj, py, slk, row;
    static char str[132];
    extern /* Subroutine */ int s4id_(integer *, integer *, integer *, 
	    integer *, char *, char *, ftnlen, ftnlen);
    static char mbnd[8], line[111], mobj[8];
    static integer jkey;
    static char mrng[8], mrhs[8], mprob[8];
    static integer iloop, jloop, isoln;
    static doublereal bplus;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    s1trim_(char *, integer *, ftnlen), s4solp_(logical *, char *, 
	    doublereal *, integer *, integer *, integer *, char *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical feasbl;
    static doublereal infbnd;
    static integer length, jstate, iprint;
    extern /* Subroutine */ int snsolf_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     doublereal *, integer *), snprnt_(integer *, char *, integer *, 
	    integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___227 = { 0, str, 0, fmt_1002, 132, 1 };
    static cilist io___229 = { 0, 0, 0, "( a)", 0 };
    static icilist io___230 = { 0, str, 0, fmt_1000, 132, 1 };
    static cilist io___231 = { 0, 0, 0, "( a)", 0 };
    static icilist io___232 = { 0, str, 0, fmt_1004, 132, 1 };
    static cilist io___233 = { 0, 0, 0, "(/a)", 0 };
    static icilist io___234 = { 0, str, 0, fmt_1005, 132, 1 };
    static cilist io___235 = { 0, 0, 0, "(/a)", 0 };
    static icilist io___236 = { 0, str, 0, fmt_1006, 132, 1 };
    static cilist io___237 = { 0, 0, 0, "( a)", 0 };
    static icilist io___238 = { 0, str, 0, fmt_1007, 132, 1 };
    static cilist io___239 = { 0, 0, 0, "( a)", 0 };
    static icilist io___240 = { 0, str, 0, fmt_1008, 132, 1 };
    static cilist io___241 = { 0, 0, 0, "( a)", 0 };
    static cilist io___242 = { 0, 0, 0, fmt_1010, 0 };
    static icilist io___243 = { 0, str, 0, "(a)", 132, 1 };
    static icilist io___244 = { 0, str, 0, "(a)", 132, 1 };
    static cilist io___263 = { 0, 0, 0, "(a)", 0 };
    static cilist io___264 = { 0, 0, 0, fmt_1020, 0 };
    static icilist io___265 = { 0, str, 0, "(a)", 132, 1 };
    static icilist io___266 = { 0, str, 0, "(a)", 132, 1 };
    static cilist io___271 = { 0, 0, 0, "(a)", 0 };
    static icilist io___272 = { 0, str, 0, fmt_1400, 132, 1 };


/*     ================================================================== */
/*     s4soln  is the standard output routine for printing the solution. */

/*     On entry, */
/*     pi    contains the dual solution. */
/*     x     contains the primal solution.  There are n + m = nb values */
/*           (n structural variables and m slacks, in that order). */
/*     Fx    contains the true slack values for nonlinear constraints. */

/*     All quantities a, bl, bu, pi, x, Fx, g are unscaled, */
/*     and adjusted in sign if maximizing. */

/*     If ondisk is true, the solution is output to the solution file. */
/*     Otherwise, it is output to the printer. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4soln. */
/*     26 Jul 1996: Slacks modified. */
/*     26 Mar 2000: Updated for SNOPT 6.1. */
/*     10 Oct 2003: snPRNT adopted. */
/*     09 Mar 2004: Solution flags are now set by snSolF. */
/*                  They are defined by the UNSCALED solution! */
/*                  Ascale(*) is no longer referenced. */
/*     05 May 2006: printing to the print file done using snPRNT. */
/*     12 Sep 2007: Solution file was missing each line of Rows and Cols. */
/*                  Added  if (ondisk), write(iSoln, '(a)') line */
/*                  to both sections. */
/*     06 Apr 2008  Updated to reflect changes to kx in snOptA. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --kx;
    names -= 8;
    --fx;
    --gobj;
    --acol;
    --inda;
    --loca;
    istate -= 4;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    iprint = iw[12];
/* Print file */
    isoln = iw[131];
/* Solution file */
    infbnd = rw[70];
/*     lvlScl    = iw( 75) ! scale option */
/* definition of an infinite bound */
    s_copy(mprob, cw + 408, (ftnlen)8, (ftnlen)8);
/* Problem name */
    s_copy(mobj, cw + 416, (ftnlen)8, (ftnlen)8);
/* Objective name */
    s_copy(mrhs, cw + 424, (ftnlen)8, (ftnlen)8);
/* Right-hand side name */
    s_copy(mrng, cw + 432, (ftnlen)8, (ftnlen)8);
/* Range name */
    s_copy(mbnd, cw + 440, (ftnlen)8, (ftnlen)8);
/* Bnd section name */
    bplus = infbnd * .1;
/* !    scale     = one */
    feasbl = *ninf == 0;
/*     maximz    = minimz .lt. 0 */
/* !    scaled    = lvlScl .gt. 0 */
    s1page_(&c__1, &iw[1], leniw);
    if (feasbl) {
	s_wsfi(&io___227);
	do_fio(&c__1, mprob, (ftnlen)8);
	do_fio(&c__1, (char *)&(*objtru), (ftnlen)sizeof(doublereal));
	e_wsfi();
	if (*ondisk) {
	    s1trim_(str, &length, (ftnlen)132);
	    io___229.ciunit = isoln;
	    s_wsfe(&io___229);
	    do_fio(&c__1, str, length);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
    } else {
	s_wsfi(&io___230);
	do_fio(&c__1, mprob, (ftnlen)8);
	do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
	e_wsfi();
	if (*ondisk) {
	    s1trim_(str, &length, (ftnlen)132);
	    io___231.ciunit = isoln;
	    s_wsfe(&io___231);
	    do_fio(&c__1, str, length);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    s_wsfi(&io___232);
    do_fio(&c__3, istate + 4, (ftnlen)4);
    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    e_wsfi();
    if (*ondisk) {
	s1trim_(str, &length, (ftnlen)132);
	io___233.ciunit = isoln;
	s_wsfe(&io___233);
	do_fio(&c__1, str, length);
	e_wsfe();
    } else {
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
    }
    if (s_cmp(mobj, "-1111111", (ftnlen)8, (ftnlen)8) != 0) {
	s_wsfi(&io___234);
	do_fio(&c__1, mobj, (ftnlen)8);
	do_fio(&c__1, objtyp + (*minimz + 1 << 3), (ftnlen)3);
	e_wsfi();
	if (*ondisk) {
	    s1trim_(str, &length, (ftnlen)132);
	    io___235.ciunit = isoln;
	    s_wsfe(&io___235);
	    do_fio(&c__1, str, length);
	    e_wsfe();
	} else {
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
	}
	s_wsfi(&io___236);
	do_fio(&c__1, mrhs, (ftnlen)8);
	e_wsfi();
	if (*ondisk) {
	    s1trim_(str, &length, (ftnlen)132);
	    io___237.ciunit = isoln;
	    s_wsfe(&io___237);
	    do_fio(&c__1, str, length);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
	s_wsfi(&io___238);
	do_fio(&c__1, mrng, (ftnlen)8);
	e_wsfi();
	if (*ondisk) {
	    s1trim_(str, &length, (ftnlen)132);
	    io___239.ciunit = isoln;
	    s_wsfe(&io___239);
	    do_fio(&c__1, str, length);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
	s_wsfi(&io___240);
	do_fio(&c__1, mbnd, (ftnlen)8);
	e_wsfi();
	if (*ondisk) {
	    s1trim_(str, &length, (ftnlen)132);
	    io___241.ciunit = isoln;
	    s_wsfe(&io___241);
	    do_fio(&c__1, str, length);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    if (*ondisk) {
	io___242.ciunit = isoln;
	s_wsfe(&io___242);
	e_wsfe();
    } else {
	s_wsfi(&io___243);
	do_fio(&c__1, " Section 1 - Rows", (ftnlen)17);
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___244);
	do_fio(&c__1, "  Number  ...Row.. State  ...Activity...  Slack Activ"
		"ity  ..Lower Limit.  ..Upper Limit.  .Dual Activity    ..i", (
		ftnlen)111);
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    }
/*     ------------------------------------------------------------------ */
/*     Output the ROWS section. */
/*     ------------------------------------------------------------------ */
    i__1 = *m;
    for (iloop = 1; iloop <= i__1; ++iloop) {
	in = iloop;
	jn = *n + in;
	i__ = kx[jn];
	j = *n + i__;
/* ! if (scaled) scale = Ascale(j) */
	b1 = bl[j];
	b2 = bu[j];
	xj = x[j];
	py = pi[i__];
	dj = rc[j];
/*        Define the row value and slack activities. */
/*        The slack activity is the distance of the row value to its */
/*        nearest bound. (For a free row, it is just the row value). */
	if (i__ <= *nncon) {
	    xj = fx[i__];
	}
	row = xj;
	d1 = b1 - xj;
	d2 = xj - b2;
	slk = -d1;
	if (abs(d1) > abs(d2)) {
	    slk = d2;
	}
	if (abs(slk) >= bplus) {
	    slk = row;
	}
/* ! SCALING NOW IGNORED. */
/* ! d1     =   d1 / scale */
/* ! d2     =   d2 / scale */
	s4id_(&jn, n, nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	snsolf_(m, n, nb, ninf, &j, &jkey, &jstate, &hs[1], &bl[1], &bu[1], &
		rc[1], &x[1], &iw[1], leniw, &rw[1], lenrw);
	s4solp_(ondisk, line, &bplus, &jkey, &jstate, &jn, id, &row, &slk, &
		b1, &b2, &py, &in, (ftnlen)111, (ftnlen)8);
	if (*ondisk) {
	    io___263.ciunit = isoln;
	    s_wsfe(&io___263);
	    do_fio(&c__1, line, (ftnlen)111);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, line, &iw[1], leniw, (ftnlen)111);
	}
    }
/*     ------------------------------------------------------------------ */
/*     Output the COLUMNS section. */
/*     ------------------------------------------------------------------ */
    s1page_(&c__1, &iw[1], leniw);
    if (*ondisk) {
	io___264.ciunit = isoln;
	s_wsfe(&io___264);
	e_wsfe();
    } else {
	s_wsfi(&io___265);
	do_fio(&c__1, " Section 2 - Columns", (ftnlen)20);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___266);
	do_fio(&c__1, "  Number  .Column. State  ...Activity...  .Obj Gradie"
		"nt.  ..Lower Limit.  ..Upper Limit.  Reduced Gradnt    m+j", (
		ftnlen)111);
	e_wsfi();
	snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
	snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    }
    i__1 = *n;
    for (jloop = 1; jloop <= i__1; ++jloop) {
	jn = jloop;
	j = kx[jn];
/* ! if (scaled) scale = Ascale(j) */
	b1 = bl[j];
	b2 = bu[j];
	xj = x[j];
	cj = 0.;
	dj = rc[j];
	i__2 = loca[j + 1] - 1;
	for (k = loca[j]; k <= i__2; ++k) {
	    ir = inda[k];
	    if (ir == *iobj) {
		cj = acol[k];
	    }
	}
	d1 = b1 - xj;
/* ! / scale */
	d2 = xj - b2;
/* ! / scale */
	if (feasbl) {
	    if (j <= *nnobj) {
		cj += gobj[j];
	    }
	}
	s4id_(&jn, n, nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	snsolf_(m, n, nb, ninf, &j, &jkey, &jstate, &hs[1], &bl[1], &bu[1], &
		rc[1], &x[1], &iw[1], leniw, &rw[1], lenrw);
	i__2 = *m + jn;
	s4solp_(ondisk, line, &bplus, &jkey, &jstate, &jn, id, &xj, &cj, &b1, 
		&b2, &dj, &i__2, (ftnlen)111, (ftnlen)8);
	if (*ondisk) {
	    io___271.ciunit = isoln;
	    s_wsfe(&io___271);
	    do_fio(&c__1, line, (ftnlen)111);
	    e_wsfe();
	} else {
	    snprnt_(&c__1, line, &iw[1], leniw, (ftnlen)111);
	}
    }
    if (*ondisk) {
	if (isoln != iprint) {
	    al__1.aerr = 0;
	    al__1.aunit = isoln;
	    f_rew(&al__1);
	}
	s_wsfi(&io___272);
	do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)132);
    }
    return 0;
} /* s4soln_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4soln */
/* Subroutine */ int s4solp_(logical *ondisk, char *line, doublereal *bplus, 
	integer *jkey, integer *jstate, integer *j, char *id, doublereal *xj, 
	doublereal *cj, doublereal *b1, doublereal *b2, doublereal *dj, 
	integer *k, ftnlen line_len, ftnlen id_len)
{
    /* Initialized data */

    static char ckey[1*5] = " " "A" "D" "I" "N";
    static char cstate[4*6] = " LL " " UL " "SBS " " BS " " EQ " " FR ";
    static char lzero[16] = "          .     ";
    static char lone[16] = "         1.0    ";
    static char lmone[16] = "        -1.0    ";
    static char none[16] = "           None ";
    static char e[10] = " 1p,e16.6,";

    /* Format strings */
    static char fmt_1000[] = "(i8,2x,a8,1x,a1,1x,a3,1p,5e16.6,i7)";

    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char form[82];

/*     ================================================================== */
/*     s4solp  prints one line of the Solution file. */

/*     The following conditions are marked by key: */

/*        D  degenerate basic or superbasic variable. */
/*        I  infeasible basic or superbasic variable. */
/*        A  alternative optimum      (degenerate nonbasic dual). */
/*        N  nonoptimal superbasic or nonbasic (infeasible dual). */

/*     Prior to 09 Mar 2004, */
/*     tests for these conditions were performed on scaled quantities */
/*     d1, d2, djtest, */
/*     since the correct indication was then more likely to be given. */
/*     On badly scaled problems, the unscaled solution could then appear */
/*     to be flagged incorrectly, but it would be just an "illusion". */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4solp. */
/*     18 Oct 1993: Replaced by modified Minos 5.4 routine m4solp. */
/*                  Infinite bounds and certain other values treated */
/*                  specially. */
/*     10 Oct 2003: snEXIT and snPRNT adopted. */
/*     09 Mar 2004: Now use jkey and jstate from snSolF. */
/*                  d1, d2, djtest are no longer used in here. */
/*     09 Mar 2004: Large xj, cj, dj (as well as b1, b2) use e format. */
/*     04 Dec 2004: Current version of s4solp. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* Select format for printing. */
    if (*ondisk) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = line_len;
	ici__1.iciunit = line;
	ici__1.icifmt = fmt_1000;
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&(*j), (ftnlen)sizeof(integer));
	do_fio(&c__1, id, (ftnlen)8);
	do_fio(&c__1, ckey + *jkey, (ftnlen)1);
	do_fio(&c__1, cstate + (*jstate << 2), (ftnlen)4);
	do_fio(&c__1, (char *)&(*xj), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*cj), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*b1), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*b2), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*dj), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	e_wsfi();
    } else {
	s_copy(form, "(i8, 2x,  a8, 1x, a1, 1x,a3, 0p,f16.5, 0p,f16.5, 0p,f1"
		"6.5, 0p,f16.5, 0p,f16.5, i7)", (ftnlen)82, (ftnlen)82);
	if (abs(*xj) >= 1e9) {
	    s_copy(form + 28, e, (ftnlen)10, (ftnlen)10);
	}
	if (abs(*cj) >= 1e9) {
	    s_copy(form + 38, e, (ftnlen)10, (ftnlen)10);
	}
	if (abs(*b1) >= 1e9) {
	    s_copy(form + 48, e, (ftnlen)10, (ftnlen)10);
	}
	if (abs(*b2) >= 1e9) {
	    s_copy(form + 58, e, (ftnlen)10, (ftnlen)10);
	}
	if (abs(*dj) >= 1e9) {
	    s_copy(form + 68, e, (ftnlen)10, (ftnlen)10);
	}
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = line_len;
	ici__1.iciunit = line;
	ici__1.icifmt = form;
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&(*j), (ftnlen)sizeof(integer));
	do_fio(&c__1, id, (ftnlen)8);
	do_fio(&c__1, ckey + *jkey, (ftnlen)1);
	do_fio(&c__1, cstate + (*jstate << 2), (ftnlen)4);
	do_fio(&c__1, (char *)&(*xj), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*cj), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*b1), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*b2), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*dj), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	e_wsfi();
    }
/*     Test for 0.0, 1.0 and -1.0 */
    if (*xj == 0.) {
	s_copy(line + 24, lzero, (ftnlen)16, (ftnlen)16);
    }
    if (*xj == 1.) {
	s_copy(line + 24, lone, (ftnlen)16, (ftnlen)16);
    }
    if (*xj == -1.) {
	s_copy(line + 24, lmone, (ftnlen)16, (ftnlen)16);
    }
    if (*cj == 0.) {
	s_copy(line + 40, lzero, (ftnlen)16, (ftnlen)16);
    }
    if (*cj == 1.) {
	s_copy(line + 40, lone, (ftnlen)16, (ftnlen)16);
    }
    if (*cj == -1.) {
	s_copy(line + 40, lmone, (ftnlen)16, (ftnlen)16);
    }
    if (*b1 == 0.) {
	s_copy(line + 56, lzero, (ftnlen)16, (ftnlen)16);
    }
    if (*b1 == 1.) {
	s_copy(line + 56, lone, (ftnlen)16, (ftnlen)16);
    }
    if (*b1 == -1.) {
	s_copy(line + 56, lmone, (ftnlen)16, (ftnlen)16);
    }
    if (*b2 == 0.) {
	s_copy(line + 72, lzero, (ftnlen)16, (ftnlen)16);
    }
    if (*b2 == 1.) {
	s_copy(line + 72, lone, (ftnlen)16, (ftnlen)16);
    }
    if (*b2 == -1.) {
	s_copy(line + 72, lmone, (ftnlen)16, (ftnlen)16);
    }
    if (*dj == 0.) {
	s_copy(line + 88, lzero, (ftnlen)16, (ftnlen)16);
    }
    if (*dj == 1.) {
	s_copy(line + 88, lone, (ftnlen)16, (ftnlen)16);
    }
    if (*dj == -1.) {
	s_copy(line + 88, lmone, (ftnlen)16, (ftnlen)16);
    }
    if (*ondisk) {
/* Relax */
    } else {
	if (*b1 < -(*bplus)) {
	    s_copy(line + 56, none, (ftnlen)16, (ftnlen)16);
	}
	if (*b2 > *bplus) {
	    s_copy(line + 72, none, (ftnlen)16, (ftnlen)16);
	}
    }
    return 0;
/* L2000: */
} /* s4solp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s4solp */
/* Subroutine */ int s4stat_(integer *k, char *istate, ftnlen istate_len)
{
    /* Initialized data */

    static char c__[4*18] = "Proc" "eedi" "ng  " "Opti" "mal " "Soln" "Infe" 
	    "asib" "le  " "Unbo" "unde" "d   " "Exce" "ss i" "tns " "Erro" 
	    "r co" "ndn ";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;

/*     ================================================================== */
/*     s4stat loads istate(*) with words describing the current state. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m4stat. */
/*     20 Apr 1999: Current version. */
/*     ================================================================== */
    /* Parameter adjustments */
    istate -= 4;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    j = 3 * min(*k,5);
    for (i__ = 1; i__ <= 3; ++i__) {
	s_copy(istate + (i__ << 2), c__ + (j + i__ - 1 << 2), (ftnlen)4, (
		ftnlen)4);
    }
    return 0;
} /* s4stat_ */

