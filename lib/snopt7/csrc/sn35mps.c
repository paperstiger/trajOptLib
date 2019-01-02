/* ./src/sn35mps.f -- translated by f2c (version 20090411).
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
static integer c__6 = 6;
static integer c__11 = 11;
static integer c__20 = 20;
static integer c__0 = 0;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__2 = 2;
static doublereal c_b257 = 0.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn35mps.f */

/*     MPSinp   MPSout */
/*     s3dflt   s3getp   s3hash   s3inpt   s3MapM   s3mps */
/*     s3mpsa   s3mpsb   s3mpsc   s3read */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int mpsinp_(integer *imps, integer *maxm, integer *maxn, 
	integer *maxne, integer *nncon, integer *nnjac, integer *nnobj, 
	integer *m, integer *n, integer *ne, integer *iobj, doublereal *
	objadd, char *prbnms, doublereal *acol, integer *inda, integer *loca, 
	doublereal *bl, doublereal *bu, char *names, integer *hs, doublereal *
	x, doublereal *pi, integer *iexit, integer *mincw, integer *miniw, 
	integer *minrw, integer *ns, char *cw, integer *lencw, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw, ftnlen prbnms_len, 
	ftnlen names_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_9830[] = "(\002 Total integer   workspace should be sign"
	    "ificantly\002,\002 more than\002,i8)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer nnl;
    static char key[4], str[80], str2[80];
    static integer lenh, lkbs;
    extern /* Subroutine */ int s3mps_(integer *, char *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), s2mem0_(
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, ftnlen);
    static integer maxcw, maxiw, maxrw;
    extern /* Subroutine */ int s3dflt_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen), s3getp_(integer *, integer *)
	    ;
    static integer negcon, lkeynm, nextcw;
    static char solver[6];
    extern /* Subroutine */ int snwrap_(integer *, char *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nextiw, lhrtyp;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer nextrw;

    /* Fortran I/O blocks */
    static icilist io___14 = { 0, str, 0, fmt_9830, 80, 1 };


/*     ------------------------------------------------------------------ */
/*     MPSinp  inputs constraint data for a linear or nonlinear program */
/*     in MPS format, consisting of NAME, ROWS, COLUMNS, RHS, RANGES and */
/*     BOUNDS sections in that order.  The RANGES and BOUNDS sections are */
/*     optional. */

/*     ------------------------------------------------------------------ */
/*     NOTE: Before calling MPSinp, your calling program MUST call the */
/*     initialization routine using the call: */
/*     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw ) */
/*     This sets the default values of the optional parameters. You can */
/*     also alter the default values of iPrint and iSumm before MPSinp */
/*     is used.  iPrint = 0, etc, is OK. */
/*     ------------------------------------------------------------------ */

/*     In the LP case, MPS format defines a set of constraints of the */
/*     form */
/*              l <= x <= u,      b1 <=  Ax  <= b2, */
/*     where l and u are specified by the BOUNDS section, and b1 and b2 */
/*     are defined somewhat indirectly by the ROWS, RHS and RANGES */
/*     sections.  snmps converts these constraints into the equivalent */
/*     form */
/*              Ax - s = 0,       bl <= ( x ) <= bu, */
/*                                      ( s ) */
/*     where s is a set of slack variables.  This is the way SNOPT deals */
/*     with the data.  The first n components of bl and bu are the same */
/*     as l and u.  The last m components are b1 and b2. */

/*     MPS format gives 8-character names to the rows and columns of A. */
/*     One of the rows of A may be regarded as a linear objective row. */
/*     This will be row iObj, where iObj = 0 means there is no such row. */

/*     The data defines a linear program if nnCon = nnJac = nnObj = 0. */
/*     The nonlinear case is the same except for a few details. */
/*     1. If nnCon = nnJac = 0 but nnObj > 0, the first nnObj columns */
/*        are associated with a nonlinear objective function. */
/*     2. If nnCon > 0, then nnJac > 0 and nnObj may be zero or positive. */
/*        The first nnCon rows and the first nnJac columns are associated */
/*        with a set of nonlinear constraints. */
/*     3. Let nnL = max( nnJac, nnObj ).  The first nnL columns */
/*        correspond to "nonlinear variables". */
/*     4. If an objective row is specified (iObj > 0), then it must be */
/*        such that iObj > nnCon. */
/*     5. "Small" elements (below the Aij tolerance) are ignored only if */
/*        they lie outside the nnCon by nnJac Jacobian, i.e. outside */
/*        the top-left corner of A. */
/*     6. No warning is given if some of the first nnL columns are empty. */


/*     ON ENTRY */
/*     ======== */
/*     iMPS   is the unit containing the MPS file.  On some systems, it */
/*            may be necessary to open file iMPS before calling snmps. */

/*     maxm   is an overestimate of the number of rows in the ROWS */
/*            section of the MPS file. */

/*     maxn   is an overestimate of the number of columns in the COLUMNS */
/*            section of the MPS file. */

/*     maxne  is an overestimate of the number of elements (matrix */
/*            coefficients) in the COLUMNS section. */

/*     nnCon  is the no. of nonlinear constraints in the problem. */
/*            These must be the FIRST rows in the ROWS section. */

/*     nnJac  is the no. of nonlinear Jacobian variables in the problem. */
/*            These must be the FIRST columns in the COLUMNS section. */

/*     nnObj  is the no. of nonlinear objective variables in the problem. */
/*            These must be the FIRST columns in the COLUMNS section, */
/*            overlapping where necessary with the Jacobian variables. */

/*     PrbNms is an array of five 8-character names. */
/*            PrbNms(1) need not be specified... it will be changed to */
/*            the name on the NAME card of the MPS file. */
/*            PrbNms(2) is the name of the objective row to be selected */
/*            from the ROWS section, or blank if snmps should select */
/*            the first type N row encountered. */
/*            Similarly, */
/*            PrbNms(3), PrbNms(4) and PrbNms(5) are the names of the */
/*            RHS, RANGES and BOUNDS to be selected from the */
/*            RHS, RANGES and BOUNDS sections respectively, or blank */
/*            if snmps should select the first ones encountered. */

/*     iw(*)  is a workspace array of length leniw.  It is needed to */
/*            hold the row-name hash table and a few other things. */

/*     leniw  is the length of iw(*).  It should be at least 4*maxm. */

/*     rw(*)  is a workspace array of length lenrw.  It is needed to */
/*            hold the row-name hash table and a few other things. */

/*     lenrw  is the length of rw(*).  It should be at least 4*maxm. */


/*     ON EXIT */
/*     ======= */
/*     m       is the number of rows in the ROWS section. */

/*     n       is the number of columns in the COLUMNS section. */

/*     ne      is the number of matrix coefficients in COLUMNS section. */

/*     iObj    is the row number of the specified objective row, */
/*             or zero if no such row was found. */

/*     ObjAdd  is a real constant extracted from row iObj of the RHS. */
/*             It is zero if the RHS contained no objective entry. */
/*             SNOPT adds ObjAdd to the objective function. */

/*     PrbNms(1)-PrbNms(5) contain the names of the */
/*             Problem, Objective row, RHS, RANGES and BOUNDS */
/*             respectively. */

/*     Acol(*) contains the ne entries for each column of the matrix */
/*             specified in the COLUMNS section. */

/*     indA(*) contains the corresponding row indices. */

/*     locA(j) (j = 1:n) points to the beginning of column j */
/*             in the parallel arrays Acol(*), indA(*). */
/*     locA(n+1) = ne+1. */

/*     bl(*)  contains n+m lower bounds for the columns and slacks. */
/*            If there is no lower bound on x(j), then bl(j) = - 1.0d+20. */

/*     bu(*)  contains n+m lower bounds for the columns and slacks. */
/*            If there is no upper bound on x(j), then bu(j) = + 1.0d+20. */

/*     Names(*) contains n+m column and row names in character*8 format. */
/*            The j-th column name is stored in Names(j). */
/*            The i-th row    name is stored in Names(k), */
/*            where k = n + i. */

/*     hs(*)  contains an initial state for each column and slack. */

/*     x(*)   contains an initial value for each column and slack. */

/*            If there is no INITIAL bounds set, */
/*               x(j) = 0 if that value lies between bl(j) and bu(j), */
/*                     = the bound closest to zero otherwise, */
/*               hs(j) = 0 if x(j) < bu(j), */
/*                     = 1 if x(j) = bu(j). */

/*            If there is an INITIAL bounds set, x(j) and hs(j) are */
/*            set as follows.  Suppose the j-th variable has the name Xj, */
/*            and suppose any numerical value specified happens to be 3. */
/*                                                   x(j)    hs(j) */
/*             FR INITIAL   Xj         3.0           3.0       -1 */
/*             FX INITIAL   Xj         3.0           3.0        2 */
/*             LO INITIAL   Xj                       bl(j)      4 */
/*             UP INITIAL   Xj                       bu(j)      5 */
/*             MI INITIAL   Xj         3.0           3.0        4 */
/*             PL INITIAL   Xj         3.0           3.0        5 */

/*     pi(*)  contains a vector defined by a special RHS called LAGRANGE. */
/*            If the MPS file contains no such RHS, pi(i) = 0.0, i=1:m. */

/*     iExit  =   0  if no fatal errors were encountered. */
/*            =  83  if there is not enough integer workspace. */
/*            = 111  if no MPS file was specified. */
/*            = 112  if maxm, maxn or maxne were too small. */
/*            = 113  if the ROWS or COLUMNS sections were empty */
/*                   or iObj > 0 but iObj <= nnCon, */

/*     nS     is the no. of FX INITIAL entries in the INITIAL bounds set. */


/*     09 Jul 1997: Original version, derived from mimps. */
/*     27 Oct 2003: Current version of MPSinp. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Problem name */
/* Objective name */
/* rhs name */
/* range name */
/*     ------------------------------------------------------------------ */
/* bounds name */
    /* Parameter adjustments */
    --pi;
    --x;
    --hs;
    names -= 8;
    --bu;
    --bl;
    --loca;
    --inda;
    --acol;
    prbnms -= 8;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    *iexit = 0;
    s_copy(solver, "MPSinp", (ftnlen)6, (ftnlen)6);
/*     ------------------------------------------------------------------ */
/*     Check memory limits and fetch the workspace starting positions. */
/*     ------------------------------------------------------------------ */
    s2mem0_(iexit, solver, lencw, leniw, lenrw, &iw[1], mincw, miniw, minrw, &
	    maxcw, &maxiw, &maxrw, &nextcw, &nextiw, &nextrw, (ftnlen)6);
    if (*iexit != 0) {
	goto L999;
    }
/*     Set undefined MPS options to default values. */
/*     Quit if no MPS file is specified. */
    s3dflt_(cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8);
    iw[21] = *nnjac;
/* # nonlinear Jacobian variables */
    iw[22] = *nnobj;
/* # variables in gObj */
    iw[23] = *nncon;
/* # of nonlinear constraints */
    s_copy(cw + 408, prbnms + 8, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 416, prbnms + 16, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 424, prbnms + 24, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 432, prbnms + 32, (ftnlen)8, (ftnlen)8);
    s_copy(cw + 440, prbnms + 40, (ftnlen)8, (ftnlen)8);
    iw[123] = *imps;
    if (*imps <= 0) {
	*iexit = 111;
	goto L800;
    }
    nnl = max(*nnjac,*nnobj);
/*     ------------------------------------------------------------------ */
/*     Allocate workspace for the MPS input routines. */
/*     s3getp finds a prime number for the length of the row hash table. */
/*     ------------------------------------------------------------------ */
    s3getp_(maxm, &lenh);
    lhrtyp = nextiw;
    lkbs = lhrtyp + 1 + *maxm;
    lkeynm = lkbs + 1 + *maxm;
    nextiw = lkeynm + lenh;
    *miniw = nextiw - 1;
    if (*miniw > maxiw) {
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*miniw), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
	*iexit = 83;
	goto L800;
    }
    s3mps_(iexit, key, maxm, maxn, maxne, &lenh, &nnl, nncon, nnjac, m, n, ne,
	     &negcon, ns, iobj, objadd, &loca[1], &inda[1], &acol[1], &bl[1], 
	    &bu[1], names + 8, &hs[1], &iw[lhrtyp], &iw[lkbs], &iw[lkeynm], &
	    x[1], &pi[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)4, (ftnlen)8, (ftnlen)8);
L800:
    if (*iexit == 0) {
	*iexit = 103;
/* MPS file read successfully */
	s_copy(prbnms + 8, cw + 408, (ftnlen)8, (ftnlen)8);
	s_copy(prbnms + 16, cw + 416, (ftnlen)8, (ftnlen)8);
	s_copy(prbnms + 24, cw + 424, (ftnlen)8, (ftnlen)8);
	s_copy(prbnms + 32, cw + 432, (ftnlen)8, (ftnlen)8);
	s_copy(prbnms + 40, cw + 440, (ftnlen)8, (ftnlen)8);
    }
    snwrap_(iexit, solver, str, str2, &iw[1], leniw, (ftnlen)6, (ftnlen)80, (
	    ftnlen)80);
L999:
    return 0;
} /* mpsinp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine MPSinp */
/* Subroutine */ int mpsout_(integer *imps, integer *m, integer *n, integer *
	ne, integer *nname, char *prbnms, doublereal *acol, integer *inda, 
	integer *loca, doublereal *bl, doublereal *bu, char *names, ftnlen 
	prbnms_len, ftnlen names_len)
{
    /* Initialized data */

    static char form[29*6] = "(4x,  a8, 2x,  a8, f8.0     )" "(4x,  a8, 2x, "
	    " a8, 1p, e14.5)" "(4x,  a8, 2x,  a8, f8.0     )" "(4x,  a8, 2x, "
	    " a8, 1p, e14.5)" "(a4,  a8, 2x,  a8, f8.0     )" "(a4,  a8, 2x, "
	    " a8, 1p, e14.5)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    cilist ci__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal b1, b2, ai;
    static char ic[8], id[8];
    static integer nb;
    static doublereal bnd, rng;
    extern /* Subroutine */ int s4id_(integer *, integer *, integer *, 
	    integer *, char *, char *, ftnlen, ftnlen);
    static logical value;
    static integer kform;
    static char bndtyp[4], rowtyp[1];

    /* Fortran I/O blocks */
    static cilist io___20 = { 0, 0, 0, "(a, 10x, a)", 0 };
    static cilist io___21 = { 0, 0, 0, "(a)", 0 };
    static cilist io___28 = { 0, 0, 0, "(1x, a1, 2x, a8)", 0 };
    static cilist io___29 = { 0, 0, 0, "(a)", 0 };
    static cilist io___34 = { 0, 0, 0, "(a)", 0 };
    static cilist io___36 = { 0, 0, 0, "(a)", 0 };
    static cilist io___38 = { 0, 0, 0, "(a)", 0 };
    static cilist io___41 = { 0, 0, 0, "(a)", 0 };


/*     ================================================================== */
/*     mpsout  outputs an MPS file to file number iMPS. */
/*     All parameters are the same as for snOptB in SNOPT 6 onwards. */
/*     They are all input parameters. */

/*     16 Aug 1998: Original version, derived from MINOS 5.4 mpsout. */
/*     27 Oct 2003: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --bu;
    --bl;
    --loca;
    --inda;
    --acol;
    names -= 8;
    prbnms -= 8;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    if (*imps <= 0) {
	return 0;
    }
    nb = *n + *m;
/*     ------------------------------------------------------------------ */
/*     ROWS section. */
/*     Note: b1 and b2 are bounds on ROWS, not slacks. */
/*           The objective row gets its name from name1(*), name2(*). */
/*           name(2) is ignored. */
/*     ------------------------------------------------------------------ */
    io___20.ciunit = *imps;
    s_wsfe(&io___20);
    do_fio(&c__1, "NAME", (ftnlen)4);
    do_fio(&c__1, prbnms + 8, (ftnlen)8);
    e_wsfe();
    io___21.ciunit = *imps;
    s_wsfe(&io___21);
    do_fio(&c__1, "ROWS", (ftnlen)4);
    e_wsfe();
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	b1 = -bu[j];
	b2 = -bl[j];
	if (b1 == b2) {
	    *(unsigned char *)rowtyp = 'E';
	} else if (b1 > -1e19) {
	    *(unsigned char *)rowtyp = 'G';
	} else if (b2 < 1e19) {
	    *(unsigned char *)rowtyp = 'L';
	} else {
	    *(unsigned char *)rowtyp = 'N';
	}
	s4id_(&j, n, &nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	io___28.ciunit = *imps;
	s_wsfe(&io___28);
	do_fio(&c__1, rowtyp, (ftnlen)1);
	do_fio(&c__1, id, (ftnlen)8);
	e_wsfe();
    }
/*     ------------------------------------------------------------------ */
/*     COLUMNS section. */
/*     Note: Objective entries get their name from name1(*), name2(*). */
/*           name(2) is ignored. */
/*     ------------------------------------------------------------------ */
    io___29.ciunit = *imps;
    s_wsfe(&io___29);
    do_fio(&c__1, "COLUMNS", (ftnlen)7);
    e_wsfe();
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s4id_(&j, n, &nb, nname, names + 8, ic, (ftnlen)8, (ftnlen)8);
	i__2 = loca[j + 1] - 1;
	for (k = loca[j]; k <= i__2; ++k) {
	    i__ = inda[k];
	    i__3 = *n + i__;
	    s4id_(&i__3, n, &nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	    ai = acol[k];
	    kform = 2;
	    if (ai == 0. || abs(ai) == 1.) {
		kform = 1;
	    }
	    ci__1.cierr = 0;
	    ci__1.ciunit = *imps;
	    ci__1.cifmt = form + (kform - 1) * 29;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, ic, (ftnlen)8);
	    do_fio(&c__1, id, (ftnlen)8);
	    do_fio(&c__1, (char *)&ai, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/*     ------------------------------------------------------------------ */
/*     RHS section. */
/*     Note: b1 and b2 are bounds on ROWS, not slacks. */
/*     ------------------------------------------------------------------ */
    io___34.ciunit = *imps;
    s_wsfe(&io___34);
    do_fio(&c__1, "RHS", (ftnlen)3);
    e_wsfe();
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	b1 = -bu[j];
	b2 = -bl[j];
	bnd = 0.;
	if (b1 == b2) {
	    bnd = b1;
	} else if (b1 > -1e19) {
	    bnd = b1;
	} else if (b2 < 1e19) {
	    bnd = b2;
	}
	if (bnd != 0.) {
	    s4id_(&j, n, &nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	    kform = 4;
	    if (abs(bnd) == 1.) {
		kform = 3;
	    }
	    ci__1.cierr = 0;
	    ci__1.ciunit = *imps;
	    ci__1.cifmt = form + (kform - 1) * 29;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, prbnms + 24, (ftnlen)8);
	    do_fio(&c__1, id, (ftnlen)8);
	    do_fio(&c__1, (char *)&bnd, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/*     ------------------------------------------------------------------ */
/*     RANGES section. */
/*     ------------------------------------------------------------------ */
    io___36.ciunit = *imps;
    s_wsfe(&io___36);
    do_fio(&c__1, "RANGES", (ftnlen)6);
    e_wsfe();
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	b1 = -bu[j];
	b2 = -bl[j];
	if (b1 < b2 && b1 > -1e19 && b2 < 1e19) {
	    rng = b2 - b1;
	    s4id_(&j, n, &nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	    kform = 4;
	    if (abs(rng) == 1.) {
		kform = 3;
	    }
	    ci__1.cierr = 0;
	    ci__1.ciunit = *imps;
	    ci__1.cifmt = form + (kform - 1) * 29;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, prbnms + 32, (ftnlen)8);
	    do_fio(&c__1, id, (ftnlen)8);
	    do_fio(&c__1, (char *)&rng, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/*     ------------------------------------------------------------------ */
/*     BOUNDS section. */
/*     ------------------------------------------------------------------ */
    io___38.ciunit = *imps;
    s_wsfe(&io___38);
    do_fio(&c__1, "BOUNDS", (ftnlen)6);
    e_wsfe();
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	b1 = bl[j];
	b2 = bu[j];
	s_copy(bndtyp, "    ", (ftnlen)4, (ftnlen)4);
	value = FALSE_;
/*        Output lower bound, except for vanilla variables. */
	if (b1 == b2) {
	    s_copy(bndtyp, " FX ", (ftnlen)4, (ftnlen)4);
	    value = TRUE_;
	} else if (b1 > -1e19) {
	    if (b1 != 0.) {
		s_copy(bndtyp, " LO ", (ftnlen)4, (ftnlen)4);
		value = TRUE_;
	    }
	} else if (b2 < 1e19) {
	    s_copy(bndtyp, " MI ", (ftnlen)4, (ftnlen)4);
	} else {
	    s_copy(bndtyp, " FR ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(bndtyp, "    ", (ftnlen)4, (ftnlen)4) != 0) {
	    s4id_(&j, n, &nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	    if (value) {
		kform = 6;
		if (b1 == 0. || abs(b1) == 1.) {
		    kform = 5;
		}
		ci__1.cierr = 0;
		ci__1.ciunit = *imps;
		ci__1.cifmt = form + (kform - 1) * 29;
		s_wsfe(&ci__1);
		do_fio(&c__1, bndtyp, (ftnlen)4);
		do_fio(&c__1, prbnms + 40, (ftnlen)8);
		do_fio(&c__1, id, (ftnlen)8);
		do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		ci__1.cierr = 0;
		ci__1.ciunit = *imps;
		ci__1.cifmt = form + (kform - 1) * 29;
		s_wsfe(&ci__1);
		do_fio(&c__1, bndtyp, (ftnlen)4);
		do_fio(&c__1, prbnms + 40, (ftnlen)8);
		do_fio(&c__1, id, (ftnlen)8);
		e_wsfe();
	    }
	}
/*        Output second bound if necessary. */
	s_copy(bndtyp, "    ", (ftnlen)4, (ftnlen)4);
	value = FALSE_;
	if (b1 == b2) {
/*           do nothing */
	} else if (b1 > -1e19) {
	    if (b2 < 1e19) {
		s_copy(bndtyp, " UP ", (ftnlen)4, (ftnlen)4);
		value = TRUE_;
	    }
	} else if (b2 < 1e19) {
	    if (b2 != 0.) {
		s_copy(bndtyp, " UP ", (ftnlen)4, (ftnlen)4);
		value = TRUE_;
	    }
	}
	if (s_cmp(bndtyp, "    ", (ftnlen)4, (ftnlen)4) != 0) {
	    s4id_(&j, n, &nb, nname, names + 8, id, (ftnlen)8, (ftnlen)8);
	    if (value) {
		kform = 6;
		if (b2 == 0. || abs(b2) == 1.) {
		    kform = 5;
		}
		ci__1.cierr = 0;
		ci__1.ciunit = *imps;
		ci__1.cifmt = form + (kform - 1) * 29;
		s_wsfe(&ci__1);
		do_fio(&c__1, bndtyp, (ftnlen)4);
		do_fio(&c__1, prbnms + 40, (ftnlen)8);
		do_fio(&c__1, id, (ftnlen)8);
		do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		ci__1.cierr = 0;
		ci__1.ciunit = *imps;
		ci__1.cifmt = form + (kform - 1) * 29;
		s_wsfe(&ci__1);
		do_fio(&c__1, bndtyp, (ftnlen)4);
		do_fio(&c__1, prbnms + 40, (ftnlen)8);
		do_fio(&c__1, id, (ftnlen)8);
		e_wsfe();
	    }
	}
    }
    io___41.ciunit = *imps;
    s_wsfe(&io___41);
    do_fio(&c__1, "ENDATA", (ftnlen)6);
    e_wsfe();
    return 0;
} /* mpsout_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine MPSout */
/* Subroutine */ int s3dflt_(char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cw_len)
{
    /* Initialized data */

    static char cblank[8] = "        ";
    static char id[4*2] = " Den" "Spar";

    /* Format strings */
    static char fmt_2000[] = "(\002 MPS file ..............\002,i10)";
    static char fmt_2010[] = "(\002 Row limit..............\002,i10,6x,\002 "
	    "Problem Number.........\002,i10,6x,\002 Lower bound default..."
	    ".\002,1p,e10.2)";
    static char fmt_2020[] = "(\002 Column limit...........\002,i10,6x,\002 "
	    "List limit.............\002,i10,6x,\002 Upper bound default..."
	    ".\002,1p,e10.2)";
    static char fmt_2030[] = "(\002 Elements limit ........\002,i10,6x,\002 "
	    "Error message limit....\002,i10,6x,\002 Aij tolerance........."
	    ".\002,1p,e10.2)";
    static char fmt_2100[] = "(\002 Jacobian...............\002,4x,a4,\002s"
	    "e\002)";

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal t;
    static char str[132];
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    snprnt_(integer *, char *, integer *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___46 = { 0, str, 0, fmt_2000, 132, 1 };
    static icilist io___47 = { 0, str, 0, fmt_2010, 132, 1 };
    static icilist io___48 = { 0, str, 0, fmt_2020, 132, 1 };
    static icilist io___49 = { 0, str, 0, fmt_2030, 132, 1 };
    static icilist io___50 = { 0, str, 0, fmt_2100, 132, 1 };


/*     ================================================================== */
/*     s3dflt  sets undefined MPS options to their default values and */
/*     prints them. */

/*     15 Feb 1998: First version. */
/*     05 Feb 2001: Current version of s3dflt. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* def. of an infinite bound */
/* def. of an 'unset' value */
/* zero Aij tolerance. */
/* default lower bound on x */
/* default upper bound on x */
/* > 0    => parms are printed */
/* 1, 0, -1  => MIN, FP, MAX */
/* 1(2)    => dense(sparse) J */
/* maximum # errors in MPS data */
/* maximum # lines  of MPS data */
/* problem number */
/* MPS file */
/* # of nonlinear constraints */
/* Estimate of the # of rows */
/* Estimate of the # of columns */
/* Estimate of the # of element */
/* Timing level */
/* Problem name */
/* Objective name */
/* rhs name */
/* range name */
/* bounds name */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    if (s_cmp(cw + 408, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(cw + 408, cblank, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(cw + 416, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(cw + 416, cblank, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(cw + 424, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(cw + 424, cblank, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(cw + 432, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(cw + 432, cblank, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(cw + 440, "-1111111", (ftnlen)8, (ftnlen)8) == 0) {
	s_copy(cw + 440, cblank, (ftnlen)8, (ftnlen)8);
    }
    if (iw[133] <= 0) {
	iw[133] = 100;
    }
    if (iw[134] <= 0) {
	iw[134] = iw[133] * 3;
    }
    if (iw[135] <= 0) {
	iw[135] = iw[134] * 5;
    }
    if (iw[81] < 0) {
	iw[81] = 1;
    }
    if (iw[123] == -11111) {
	iw[123] = 0;
    }
    if (iw[182] < 0) {
	iw[182] = 3;
    }
    if (iw[108] < 0) {
	iw[108] = 0;
    }
    if (iw[105] < 0) {
	iw[105] = 1;
    }
    if (iw[106] < 0) {
	iw[106] = 10;
    }
    if (iw[107] < 0) {
	iw[107] = 0;
    }
    if (iw[87] == -11111) {
	iw[87] = 1;
    }
    if (rw[70] < 0.) {
	rw[70] = 1e20;
    }
    if (rw[95] < 0.) {
	rw[95] = 1e-10;
    }
    if (rw[96] == rw[69]) {
	rw[96] = 0.;
    }
    if (rw[97] == rw[69]) {
	rw[97] = rw[70];
    }
    if (rw[96] > rw[97]) {
	t = rw[96];
	rw[96] = rw[97];
	rw[97] = t;
    }
/*     ------------------------------------------------------------------ */
/*     Print parameters unless SUPPRESS PARAMETERS was specified. */
/*     ------------------------------------------------------------------ */
    if (iw[81] > 0) {
	s1page_(&c__1, &iw[1], leniw);
	snprnt_(&c__1, " MPS Input Data", &iw[1], leniw, (ftnlen)15);
	snprnt_(&c__1, " ==============", &iw[1], leniw, (ftnlen)15);
	s_wsfi(&io___46);
	do_fio(&c__1, (char *)&iw[123], (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___47);
	do_fio(&c__1, (char *)&iw[133], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iw[108], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rw[96], (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___48);
	do_fio(&c__1, (char *)&iw[134], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iw[107], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rw[97], (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	s_wsfi(&io___49);
	do_fio(&c__1, (char *)&iw[135], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iw[106], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rw[95], (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	if (iw[23] > 0) {
	    s_wsfi(&io___50);
	    do_fio(&c__1, id + (iw[105] - 1 << 2), (ftnlen)4);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)132);
	}
    }
    return 0;
} /* s3dflt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3dflt */
/* Subroutine */ int s3getp_(integer *maxm, integer *lenh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;

/*     ------------------------------------------------------------------ */
/*     s3getp finds a prime number lenh suitably larger than maxm. */
/*     It is used as the length of the hash table for the MPS row names. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    *lenh = *maxm << 1;
    *lenh = max(*lenh,100);
    *lenh = (*lenh / 2 << 1) - 1;
    k = *lenh / 20 + 6;
L100:
    ++k;
    *lenh += 2;
    i__1 = k;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	if (*lenh % i__ == 0) {
	    goto L100;
	}
    }
    return 0;
} /* s3getp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3getp */
/* Subroutine */ int s3hash_(integer *len, integer *nen, integer *ncoll, char 
	*key, integer *mode, integer *keytab, char *names, integer *ka, 
	logical *found, ftnlen key_len, ftnlen names_len)
{
    /* System generated locals */
    char ch__1[8];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer ia, ic, iq, ir, is, jq, jr, kt, ix, len2, ikey;

/*     ================================================================== */
/*     s3hash looks up and/or inserts keys in a hash table. */
/*     Reference:  R.P. Brent, CACM 16,2 (Feb 1973), pp. 105-109. */
/*     This version is simplified for the case where no entries are */
/*     deleted. */

/*     keytab is used as an index into a consecutive list of unique */
/*     identifiers Names(*), to find if "key" is in the list. */

/*     On entry, */
/*     len     is the (fixed) length of the hash table. */
/*             It must be significantly more than nen. */
/*     nen     is the number of unique identifiers. */
/*     ncoll   is the number of collisions so far. */
/*     key     is the identifier to be looked up or inserted. */
/*     mode    = 1 if key is just being looked up, */
/*             = 2 if key should be inserted into the hash table. */
/*     keytab  is the integer hash table. */
/*     Names   is the list of identifiers. */

/*     On exit, */
/*     ncoll   is (new) number of collisions so far. */
/*     ka      satisfies Names(ka) = key if key is in the hash table. */
/*     found   is true if key was already in the hash table on entry. */

/*     Approx 1975: First f66 version written for MINOS, using integer */
/*                  names. */
/*     16 Sep 1997: key and Names(*) are now character*8. */
/*                  Statement function "hash" introduced.  It should be */
/*                  the only thing different between f77 and f90. */
/*     ================================================================== */
/*     The following statement function "hash" returns an integer from */
/*     its character*8 argument.  f77 doesn't permit great efficiency */
/*     (or we would have done it this way long ago!).  f90 should be */
/*     better with the transfer function. */
/* -->  f77 */
/*     16 Sep 1997: Try looking at just a few important characters. */
/*     02 Oct 1997: Adding ichar of p(1,2,5,8) gave over a million */
/*                  collisions on pilots.mps. */
/*                  Try left shifts and more chars. */
/* -->  f90 */
/*     16 Sep 1997: Simulate the function used in MINOS (with Name(*) */
/*                  treated as two integer arrays).  It seems reasonably */
/*                  efficient in terms of the number of collisions. */
/*     hash(p) =   transfer( p(1:4), mode ) */
/*    &          - transfer( p(5:8), mode ) */
/*     ================================================================== */
    /* Parameter adjustments */
    --keytab;
    names -= 8;

    /* Function Body */
    len2 = *len - 2;
    ic = -1;
/*     Compute address of first probe (ir) and increment (iq). */
    s_copy(ch__1, key, (ftnlen)8, (ftnlen)8);
    ikey = 16 * (16 * (16 * (16 * (16 * *(unsigned char *)&ch__1[0] + *(
	    unsigned char *)&ch__1[0]) + *(unsigned char *)&ch__1[0]) + *(
	    unsigned char *)&ch__1[0]) + *(unsigned char *)&ch__1[0]) + *(
	    unsigned char *)&ch__1[0];
    iq = ikey % len2 + 1;
    ir = ikey % *len + 1;
    *ka = ir;
/*     Look in the table. */
L20:
    kt = keytab[*ka];
/*     Check for an empty space or a match. */
    if (kt == 0) {
	goto L30;
    }
    if (s_cmp(key, names + (kt << 3), (ftnlen)8, (ftnlen)8) == 0) {
	goto L60;
    }
    ++ic;
    ++(*ncoll);
/*     Compute address of next probe. */
    *ka += iq;
    if (*ka > *len) {
	*ka -= *len;
    }
/*     See if whole table has been searched. */
    if (*ka != ir) {
	goto L20;
    }
/*     The key is not in the table. */
L30:
    *found = FALSE_;
/*     Return with ka = 0 unless an entry has to be made. */
    if (*mode == 2 && ic <= len2) {
	goto L70;
    }
    *ka = 0;
    return 0;
L60:
    *found = TRUE_;
    return 0;
/*     Look for the best way to make an entry. */
L70:
    if (ic <= 0) {
	return 0;
    }
    ia = *ka;
    is = 0;
/*     Compute the maximum length to search along current chain. */
L80:
    ix = ic - is;
    kt = keytab[ir];
/*     Compute increment jq for current chain. */
    s_copy(ch__1, names + (kt << 3), (ftnlen)8, (ftnlen)8);
    ikey = 16 * (16 * (16 * (16 * (16 * *(unsigned char *)&ch__1[0] + *(
	    unsigned char *)&ch__1[0]) + *(unsigned char *)&ch__1[0]) + *(
	    unsigned char *)&ch__1[0]) + *(unsigned char *)&ch__1[0]) + *(
	    unsigned char *)&ch__1[0];
    jq = ikey % len2 + 1;
    jr = ir;
/*     Look along the chain. */
L90:
    jr += jq;
    if (jr > *len) {
	jr -= *len;
    }
/*     Check for a hole. */
    if (keytab[jr] == 0) {
	goto L100;
    }
    --ix;
    if (ix > 0) {
	goto L90;
    }
    goto L110;
/*     Save location of hole. */
L100:
    ia = jr;
    *ka = ir;
    ic -= ix;
/*     Move down to the next chain. */
L110:
    ++is;
    ir += iq;
    if (ir > *len) {
	ir -= *len;
    }
/*     Go back if a better hole might still be found. */
    if (ic > is) {
	goto L80;
    }
/*     If necessary move an old entry. */
    if (ia != *ka) {
	keytab[ia] = keytab[*ka];
    }
    return 0;
} /* s3hash_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3hash */
/* Subroutine */ int s3inpt_(integer *iexit, integer *maxm, integer *maxn, 
	integer *maxne, integer *nncon, integer *nnjac, integer *nnobj, 
	integer *m, integer *n, integer *ne, integer *iobj, doublereal *
	objadd, integer *nextcw, integer *nextiw, integer *nextrw, integer *
	mincw, integer *miniw, integer *minrw, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	cw_len)
{
    /* Format strings */
    static char fmt_9820[] = "(\002 Total character workspace should be sign"
	    "ificantly\002,\002 more than\002,i8)";
    static char fmt_9830[] = "(\002 Total integer   workspace should be sign"
	    "ificantly\002,\002 more than\002,i8)";
    static char fmt_9840[] = "(\002 Total real      workspace should be sign"
	    "ificantly\002,\002 more than\002,i8)";

    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_rew(alist *), s_cmp(char *, char *, ftnlen, ftnlen), s_rsfe(
	    cilist *), e_rsfe(void);

    /* Local variables */
    static integer ns, kx, lx, kbl, lbl, kbu, lbu, kpi, lpi, khs, lhs, nnl;
    static char key[4], str[80];
    static integer lenh, lkbs, imps;
    extern /* Subroutine */ int s3mps_(integer *, char *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static integer linda, kloca, lacol, lloca, ldenj, istdi;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer newcw;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer newiw, newrw;
    extern /* Subroutine */ int s3mapm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     s3getp_(integer *, integer *);
    static integer knames, lnames;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static integer ispecs, lkeynm, idummy, lhrtyp;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___95 = { 0, str, 0, fmt_9820, 80, 1 };
    static icilist io___96 = { 0, str, 0, fmt_9830, 80, 1 };
    static icilist io___97 = { 0, str, 0, fmt_9840, 80, 1 };
    static cilist io___99 = { 0, 0, 0, "(a4)", 0 };


/*     ================================================================== */
/*     s3inpt  inputs constraint data in MPS format, and sets up */
/*     various quantities as follows: */

/*     ObjAdd     (output) is minus the coefficient in row iObj of the */
/*                RHS section (zero by default).  SNOPT adds it to the */
/*                objective function. */

/*     m, n, ne   are the number of rows, columns and elements in A. */

/*     iObj       is the row number for the linear objective (if any). */
/*                It must come after any nonlinear rows. */
/*                iObj = 0 if there is no linear objective. */

/*     locA, indA, Acol */
/*                is the matrix A, */
/*                stored in locations llocA, lindA, lAcol */
/*     bl, bu     are the bounds,  stored at locations lbl, lbu. */
/*     hs, x      are states and values,   stored at   lhs, lx. */

/*     hs(j)      is set to  0, 1  to indicate a plausible initial state */
/*                (at lo or up bnd) for each variable j  (j = 1 to n+m). */
/*                If crash is to be used, i.e., crash option gt 0 and */
/*                if no basis file will be supplied, the initial bounds */
/*                set may initialize hs(j) as follows to assist crash: */

/*     -1      if column or row j is likely to be in the optimal basis, */
/*      4      if column j is likely to be nonbasic at its lower bound, */
/*      5      if column j is likely to be nonbasic at its upper bound, */
/*      2      if column or row j should initially be superbasic, */
/*      0 or 1 otherwise. */

/*     x(j)       is a corresponding set of initial values. */
/*                Safeguards are applied later by SNOPT, so the */
/*                values of hs and x are not desperately critical. */

/*     18 Nov 1991: First version based on Minos routine m3inpt. */
/*     09 Nov 2000: miniw, minrw, mincw set elsewhere. */
/*     27 Oct 2003: Current version of s3inpt. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    istdi = iw[9];
/* Standard Input */
    ispecs = iw[11];
/* Specs (options) file */
    ldenj = iw[105];
/* 1(2)    => dense(sparse) Jacobian */
    imps = iw[123];
/* MPS file */
    nnl = max(*nnjac,*nnobj);
    newcw = *nextcw;
    newiw = *nextiw;
    newrw = *nextrw;
/*     Start.  We may come back here to try again with more workspace. */
/*     key retains the first 4 characters of the NAME, ROWS, COLUMNS */
/*     RHS, RANGES and BOUNDS cards. */
/*     s3getp finds a prime number for the length of the row hash table. */
L100:
    *iexit = 0;
    s3getp_(maxm, &lenh);
/*     Allocate temporary addresses for the problem data */
/*     based on maxm, maxn and maxne. */
    s3mapm_(iexit, maxm, maxn, maxne, &lenh, nextcw, nextiw, nextrw, mincw, 
	    miniw, minrw, lencw, leniw, lenrw, &iw[1]);
    if (*iexit != 0) {
	goto L700;
    }
/*     Get the data from an MPS file. */
/*     Allocate space for the arguments of sqopt and s5solv. */
/*     These are for the data, */
/*        locA, indA, Acol, bl, bu, Names */
/*     and for the solution */
/*        hs, x, pi, rc, hs. */
/*     The true values of m, n and ne are also found. */
    lacol = iw[256];
/* Jcol(ne)    = Constraint Jacobian by columns */
    linda = iw[258];
/* indJ(ne) holds the row indices for Jij */
    lloca = iw[257];
/* locJ(n+1)   = column pointers for indJ */
    lbl = iw[271];
/* bl(nb)      = lower bounds */
    lbu = iw[272];
/* bu(nb)      = upper bounds */
    lx = iw[299];
/* x(nb)       = the solution (x,s) */
    lpi = iw[279];
/* pi(m)       = the pi-vector */
    lhs = iw[282];
/* the column state vector */
    lnames = iw[359];
/*     These are used as work arrays */
/* Names(nName) */
    lhrtyp = iw[284];
/* hhrtyp(mBS), feasibility types */
    lkbs = iw[292];
/* kBS(mBS)    = ( B  S ) list */
    lkeynm = iw[360];
/* keynm(lenh) = hash table keys */
    s3mps_(iexit, key, maxm, maxn, maxne, &lenh, &nnl, nncon, nnjac, m, n, ne,
	     &iw[20], &ns, iobj, objadd, &iw[lloca], &iw[linda], &rw[lacol], &
	    rw[lbl], &rw[lbu], cw + (lnames << 3), &iw[lhs], &iw[lhrtyp], &iw[
	    lkbs], &iw[lkeynm], &rw[lx], &rw[lpi], cw + 8, lencw, &iw[1], 
	    leniw, &rw[1], lenrw, (ftnlen)4, (ftnlen)8, (ftnlen)8);
    if (*iexit != 0) {
	goto L700;
    }
/*     negCon  counted the actual Jacobian entries in the MPS file, */
/*     but for Jacobian = Dense, we have to reset it. */
    if (ldenj == 1) {
	iw[20] = *nncon * *nnjac;
    }
/*     ------------------------------------------------------------------ */
/*     Compress storage, now that we know the size of everything. */
/*     ------------------------------------------------------------------ */
/*     Save current positions of  bl, bu, etc. */
    knames = lnames;
    kloca = lloca;
    khs = lhs;
    kbl = lbl;
    kbu = lbu;
    kx = lx;
    kpi = lpi;
    *nextcw = newcw;
    *nextiw = newiw;
    *nextrw = newrw;
    s3mapm_(iexit, m, n, ne, &lenh, nextcw, nextiw, nextrw, mincw, miniw, 
	    minrw, lencw, leniw, lenrw, &iw[1]);
    if (*iexit != 0) {
	goto L700;
    }
    lacol = iw[256];
/* Jcol(ne)    = Constraint Jacobian by columns */
    linda = iw[258];
/* indJ(ne) holds the row indices for Jij */
    lloca = iw[257];
/* locJ(n+1)   = column pointers for indJ */
    lbl = iw[271];
/* bl(nb)      = lower bounds */
    lbu = iw[272];
/* bu(nb)      = upper bounds */
    lpi = iw[279];
/* pi(m)       = the pi-vector */
    lhs = iw[282];
/* the column state vector */
    lhrtyp = iw[284];
/* hhrtyp(mBS), feasibility types */
    lkbs = iw[292];
/* kBS(mBS)    = ( B  S ) list */
    lx = iw[299];
/* x(nb)       = the solution (x,s) */
    lnames = iw[359];
/* Names(nName) */
    lkeynm = iw[360];
/*     Move bl, bu, etc. into their final positions. */
/* keynm(lenh) = hash table keys */
    i__1 = *n + 1;
    icopy_(&i__1, &iw[kloca], &c__1, &iw[lloca], &c__1);
    i__1 = *n + *m;
    chcopy_(&i__1, cw + (knames << 3), &c__1, cw + (lnames << 3), &c__1, (
	    ftnlen)8, (ftnlen)8);
    i__1 = *n + *m;
    icopy_(&i__1, &iw[khs], &c__1, &iw[lhs], &c__1);
    i__1 = *n + *m;
    dcopy_(&i__1, &rw[kbl], &c__1, &rw[lbl], &c__1);
    i__1 = *n + *m;
    dcopy_(&i__1, &rw[kbu], &c__1, &rw[lbu], &c__1);
    i__1 = *n + *m;
    dcopy_(&i__1, &rw[kx], &c__1, &rw[lx], &c__1);
    dcopy_(m, &rw[kpi], &c__1, &rw[lpi], &c__1);
/*     ------------------------------------------------------------------ */
/*     Check for error exits. */
/*     ------------------------------------------------------------------ */
L700:
    if (*iexit == 82) {
	s_wsfi(&io___95);
	do_fio(&c__1, (char *)&(*mincw), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    } else if (*iexit == 83) {
	s_wsfi(&io___96);
	do_fio(&c__1, (char *)&(*miniw), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    } else if (*iexit == 84) {
	s_wsfi(&io___97);
	do_fio(&c__1, (char *)&(*minrw), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)80);
    } else if (*iexit == 112) {
/* Problem estimates too small */
	if (*m >= *maxm) {
	    *maxm = *m;
	} else {
	    *maxn = *n;
	    *maxne = *ne;
	}
	if (imps != istdi && imps != ispecs) {
	    al__1.aerr = 0;
	    al__1.aunit = imps;
	    f_rew(&al__1);
	    goto L100;
	}
    }
    if (*iexit == 112 || *iexit == 113) {
/*        ------------------------------------ */
/*        Flush MPS file to the ENDATA card */
/*        if it is the same as the SPECS file. */
/*        ------------------------------------ */
	if (imps == ispecs) {
	    for (idummy = 1; idummy <= 100000; ++idummy) {
		if (s_cmp(key, "ENDA", (ftnlen)4, (ftnlen)4) == 0) {
		    goto L900;
		}
		io___99.ciunit = imps;
		s_rsfe(&io___99);
		do_fio(&c__1, key, (ftnlen)4);
		e_rsfe();
	    }
	}
    }
/*     Exit. */
L900:
    if (imps != istdi && imps != ispecs) {
	al__1.aerr = 0;
	al__1.aunit = imps;
	f_rew(&al__1);
    }
    return 0;
} /* s3inpt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3inpt */
/* Subroutine */ int s3mapm_(integer *iexit, integer *m, integer *n, integer *
	ne, integer *lenh, integer *nextcw, integer *nextiw, integer *nextrw, 
	integer *mincw, integer *miniw, integer *minrw, integer *lencw, 
	integer *leniw, integer *lenrw, integer *iw)
{
    static integer nb, lx, lbl, lrc, lbu, lpi, lhs, lkbs, linda, lloca, lacol,
	     lhfeas, lnames, lkeynm, lhetyp;

/*     ================================================================== */
/*     s3MapM   is called from s3inpt to find the minimum storage needed */
/*     for mps input. */

/*     Storage is estimated using the dimensions: */
/*        m    , n    , ne */

/*     18 Feb 1998: First version. */
/*     27 Oct 2003: Current version of s3MapM. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     Allocate space for the data, */
/*              locA, indA, Acol, hEtype, bl, bu, Names */
/*     and for the solution */
/*              hs, x, pi, rc, hs. */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    *iexit = 0;
    nb = *n + *m;
    lnames = *nextcw;
    *nextcw = lnames + nb;
    linda = *nextiw;
    lloca = linda + *ne;
    lhetyp = lloca + *n + 1;
    lhs = lhetyp + nb;
    *nextiw = lhs + nb;
    lacol = *nextrw;
    lbl = lacol + *ne;
    lbu = lbl + nb;
    lx = lbu + nb;
    lpi = lx + nb;
    lrc = lpi + *m;
    *nextrw = lrc + nb;
/*     Allocate arrays needed during and after MPS input. */
    lhfeas = *nextiw;
    lkbs = lhfeas + nb;
    lkeynm = lkbs + nb;
    *nextiw = lkeynm + *lenh;
    iw[257] = lloca;
    iw[258] = linda;
    iw[282] = lhs;
    iw[283] = lhetyp;
    iw[359] = lnames;
    iw[256] = lacol;
    iw[271] = lbl;
    iw[272] = lbu;
    iw[279] = lpi;
    iw[280] = lrc;
    iw[299] = lx;
    iw[284] = lhfeas;
    iw[292] = lkbs;
    iw[360] = lkeynm;
    *mincw = *nextcw - 1;
    *miniw = *nextiw - 1;
    *minrw = *nextrw - 1;
    if (*mincw > *lencw) {
	*iexit = 82;
    } else if (*miniw > *leniw) {
	*iexit = 83;
    } else if (*minrw > *lenrw) {
	*iexit = 84;
    }
    return 0;
} /* s3mapm_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3MapM */
/* Subroutine */ int s3mps_(integer *iexit, char *key, integer *maxm, integer 
	*maxn, integer *maxne, integer *lenh, integer *nnl, integer *nncon, 
	integer *nnjac, integer *m, integer *n, integer *ne, integer *negcon, 
	integer *ns, integer *iobj, doublereal *objadd, integer *loca, 
	integer *inda, doublereal *acol, doublereal *bl, doublereal *bu, char 
	*names, integer *hs, integer *hrtype, integer *kbs, integer *keynam, 
	doublereal *x, doublereal *pi, char *cw, integer *lencw, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw, ftnlen key_len, 
	ftnlen names_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1300[] = "(\002 Length of row-name hash table  \002,i12)";
    static char fmt_1310[] = "(\002 Collisions during table lookup \002,i12)";
    static char fmt_1320[] = "(\002 No. of rejected coefficients   \002,i12)";
    static char fmt_1350[] = "(\002 No. of Jacobian entries specified\002,i1"
	    "0)";
    static char fmt_1400[] = "(\002 No. of INITIAL  bounds  specified\002,i1"
	    "0)";
    static char fmt_1410[] = "(\002 No. of superbasics specified   \002,i12)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer na0, ier[20];
    static char str[80];
    static integer line, lrow, ncard[6];
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer ncoll;
    extern /* Subroutine */ int s3mpsa_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, char *, integer *, 
	    integer *, integer *, char *, char *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), s3mpsb_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, char *, 
	    integer *, integer *, integer *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     char *, integer *, integer *, integer *, doublereal *, integer *,
	     ftnlen, ftnlen, ftnlen), s3mpsc_(integer *, integer *, integer *,
	     integer *, integer *, integer *, char *, integer *, char *, 
	    doublereal *, doublereal *, integer *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static integer lennam;
    extern /* Subroutine */ int chcopy_(integer *, char *, integer *, char *, 
	    integer *, ftnlen, ftnlen), snprnt_(integer *, char *, integer *, 
	    integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___123 = { 0, str, 0, fmt_1300, 80, 1 };
    static icilist io___124 = { 0, str, 0, fmt_1310, 80, 1 };
    static icilist io___125 = { 0, str, 0, fmt_1320, 80, 1 };
    static icilist io___126 = { 0, str, 0, fmt_1350, 80, 1 };
    static icilist io___127 = { 0, str, 0, fmt_1400, 80, 1 };
    static icilist io___128 = { 0, str, 0, fmt_1410, 80, 1 };


/*     ================================================================== */
/*     s3mps  does the work for MPSinp and s3inpt. */

/*     iExit  status */
/*     -----  ------ */
/*      112   too many rows, columns or elements */
/*      113   no rows or columns, or iObj le nnCon */

/*     18 Feb 1998: First version. */
/*     03 Aug 2003: snPRNT adopted. */
/*     26 Oct 2003: Current version of s3mps. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ncard  counts the number of data records in each section. */
    /* Parameter adjustments */
    --pi;
    --hrtype;
    --x;
    --kbs;
    --hs;
    names -= 8;
    --bu;
    --bl;
    --loca;
    --acol;
    --inda;
    --keynam;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    *iexit = 0;
    ncoll = 0;
    s_copy(key, "    ", (ftnlen)4, (ftnlen)4);
    iload_(&c__6, &ncoll, ncard, &c__1);
/*     ------------------------------------------------------------------ */
/*     Input ROWS. */
/*     lrow   is the location of the first rowname in Names. */
/*            The column names will later go at the front. */
/*     lenNam is the initial length of Names, */
/*            i.e. the maximum no. of names allowed for. */
/*     ------------------------------------------------------------------ */
    lrow = *maxn + 1;
    lennam = *maxn + *maxm;
    s3mpsa_(iexit, ier, &line, maxm, maxn, iobj, &ncoll, m, &lrow, &lennam, 
	    lenh, nnl, nncon, key, ncard, &hrtype[1], &keynam[1], names + 8, 
	    cw + 8, lencw, &iw[1], leniw, (ftnlen)4, (ftnlen)8, (ftnlen)8);
    if (*iexit != 0) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     m  is now known. */
/*     Input COLUMNS, RHS, RANGES. */
/*     ------------------------------------------------------------------ */
    s3mpsb_(iexit, ier, &line, maxn, maxne, &lrow, &lennam, lenh, &ncoll, 
	    iobj, objadd, m, n, ne, nnl, nncon, nnjac, negcon, &na0, key, 
	    ncard, &hrtype[1], &keynam[1], names + 8, &loca[1], &inda[1], &
	    acol[1], &bl[1], &bu[1], &kbs[1], &pi[1], cw + 8, lencw, &iw[1], 
	    leniw, &rw[1], lenrw, (ftnlen)4, (ftnlen)8, (ftnlen)8);
    if (*iexit != 0) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     n  and  ne  are now known. */
/*     Move the row names to be contiguous with the column names. */
/*     Input BOUNDS. */
/*     ------------------------------------------------------------------ */
    if (lrow > *n + 1) {
	chcopy_(m, names + (lrow << 3), &c__1, names + (*n + 1 << 3), &c__1, (
		ftnlen)8, (ftnlen)8);
    }
    s3mpsc_(ier, &line, m, n, ns, &lennam, key, ncard, names + 8, &bl[1], &bu[
	    1], &hs[1], &x[1], cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
	    ftnlen)4, (ftnlen)8, (ftnlen)8);
    s_wsfi(&io___123);
    do_fio(&c__1, (char *)&(*lenh), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)80);
    s_wsfi(&io___124);
    do_fio(&c__1, (char *)&ncoll, (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    if (na0 > 0) {
	s_wsfi(&io___125);
	do_fio(&c__1, (char *)&na0, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    }
    if (*nncon > 0) {
	s_wsfi(&io___126);
	do_fio(&c__1, (char *)&(*negcon), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    }
    if (*nnl > 0 || ncard[5] > 0) {
	s_wsfi(&io___127);
	do_fio(&c__1, (char *)&ncard[5], (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	s_wsfi(&io___128);
	do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
    }
L900:
    return 0;
} /* s3mps_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3mps */
/* Subroutine */ int s3mpsa_(integer *iexit, integer *ier, integer *line, 
	integer *maxm, integer *maxn, integer *iobj, integer *ncoll, integer *
	m, integer *lrow, integer *lennm, integer *lenh, integer *nnl, 
	integer *nncon, char *key, integer *ncard, integer *hrtype, integer *
	keynam, char *names, char *cw, integer *lencw, integer *iw, integer *
	leniw, ftnlen key_len, ftnlen names_len, ftnlen cw_len)
{
    /* Initialized data */

    static char cblank[8] = "        ";
    static char lname[4] = "NAME";
    static char lrows[4] = "ROWS";
    static char lcolu[4] = "COLU";
    static char lex[4] = " E  ";
    static char lgx[4] = " G  ";
    static char llx[4] = " L  ";
    static char lnx[4] = " N  ";
    static char lxe[4] = "  E ";
    static char lxg[4] = "  G ";
    static char lxl[4] = "  L ";
    static char lxn[4] = "  N ";

    /* Format strings */
    static char fmt_5000[] = "(\002 Name   \002,a8)";
    static char fmt_1170[] = "(\002 ===>  Note --- row  \002,a8,\002  select"
	    "ed as linear part of objective.\002)";
    static char fmt_1160[] = "(\002 XXXX  Illegal row type at line\002,i7"
	    ",\002... \002,a4,a8)";
    static char fmt_1200[] = "(\002 XXXX  Duplicate row name --\002,a8,\002 "
	    "-- ignored\002)";
    static char fmt_3030[] = "(\002 XXXX  Too many rows.  Limit was\002,i8,4"
	    "x,\002  Actual number is\002,i8)";
    static char fmt_5100[] = "(\002 Rows   \002,i8)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    char ch__1[48];

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer ia;
    static char id[8*3];
    static integer it;
    static char str[100];
    static integer imps, mlst, jrow;
    static doublereal aelem[2];
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static logical found, gotnm;
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *), 
	    s3read_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, char *, char *, doublereal *, integer *, ftnlen, 
	    ftnlen), s3hash_(integer *, integer *, integer *, char *, integer 
	    *, integer *, char *, integer *, logical *, ftnlen, ftnlen);
    static integer inform__, mpserr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___149 = { 0, str, 0, fmt_5000, 100, 1 };
    static icilist io___151 = { 0, str, 0, fmt_1170, 100, 1 };
    static icilist io___152 = { 0, str, 0, fmt_1160, 100, 1 };
    static icilist io___155 = { 0, str, 0, fmt_1200, 100, 1 };
    static icilist io___157 = { 0, str, 0, fmt_3030, 100, 1 };
    static icilist io___158 = { 0, str, 0, fmt_5100, 100, 1 };


/*     ================================================================== */
/*     s3mpsa  inputs the name and rows sections of an MPS file. */

/*     Original version written by Keith Morris, Wellington, 1973. */

/*     15 Nov 1991: First version based on Minos routine m3mpsa. */
/*     19 Jul 1997: Thread-safe version. */
/*     24 Sep 1997: key and Names(*) are now character*8. */
/*     03 Aug 2003: snPRNT adopted. */
/*     03 Aug 2003: Current version of s3mpsa. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ier;
    --hrtype;
    names -= 8;
    --keynam;
    --ncard;
    cw -= 8;
    --iw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    imps = iw[123];
/* MPS file */
    mpserr = iw[106];
/* maximum # errors in MPS data */
    mlst = iw[107];
/* maximum # lines  of MPS data */
    s1page_(&c__1, &iw[1], leniw);
    snprnt_(&c__1, " MPS file", &iw[1], leniw, (ftnlen)9);
    snprnt_(&c__1, " --------", &iw[1], leniw, (ftnlen)9);
    inform__ = 0;
    *iobj = 0;
    *line = 0;
    *m = 0;
    gotnm = s_cmp(cw + 416, cblank, (ftnlen)8, (ftnlen)8) != 0;
    iload_(&c__20, &c__0, &ier[1], &c__1);
    iload_(lenh, &c__0, &keynam[1], &c__1);
/*     Look for the NAME card. */
L10:
    s3read_(&c__1, &imps, &iw[1], leniw, line, &c__5, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (s_cmp(key, lname, (ftnlen)4, (ftnlen)4) != 0) {
	if (ier[1] == 0) {
	    ier[1] = 1;
	    snprnt_(&c__3, " XXXX  Garbage before NAME card", &iw[1], leniw, (
		    ftnlen)31);
	}
	goto L10;
    }
    s_copy(cw + 408, id + 8, (ftnlen)8, (ftnlen)8);
    s_wsfi(&io___149);
    do_fio(&c__1, cw + 408, (ftnlen)8);
    e_wsfi();
    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)100);
/*     Look for the ROWS card. */
    s3read_(&c__1, &imps, &iw[1], leniw, line, &c__5, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    inform__ = 0;
    if (s_cmp(key, lrows, (ftnlen)4, (ftnlen)4) != 0) {
	++ier[1];
	snprnt_(&c__3, " XXXX  ROWS card not found", &iw[1], leniw, (ftnlen)
		26);
	goto L35;
    }
/*     ================================================================== */
/*     Read the row names and check if the relationals are valid. */
/*     ================================================================== */
L30:
    s3read_(&c__1, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L110;
    }
L35:
    if (s_cmp(key, lgx, (ftnlen)4, (ftnlen)4) == 0 || s_cmp(key, lxg, (ftnlen)
	    4, (ftnlen)4) == 0) {
	it = -1;
    } else if (s_cmp(key, lex, (ftnlen)4, (ftnlen)4) == 0 || s_cmp(key, lxe, (
	    ftnlen)4, (ftnlen)4) == 0) {
	it = 0;
    } else if (s_cmp(key, llx, (ftnlen)4, (ftnlen)4) == 0 || s_cmp(key, lxl, (
	    ftnlen)4, (ftnlen)4) == 0) {
	it = 1;
    } else if (s_cmp(key, lnx, (ftnlen)4, (ftnlen)4) == 0 || s_cmp(key, lxn, (
	    ftnlen)4, (ftnlen)4) == 0) {
	it = 2;
/*        Record objective name if we don't already have one. */
	if (*iobj == 0) {
	    if (! gotnm) {
		s_copy(cw + 416, id, (ftnlen)8, (ftnlen)8);
		if (*nnl > 0) {
		    s_wsfi(&io___151);
		    do_fio(&c__1, cw + 416, (ftnlen)8);
		    e_wsfi();
		    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)100);
		}
	    }
	    if (s_cmp(id, cw + 416, (ftnlen)8, (ftnlen)8) == 0) {
		*iobj = *m + 1;
		++ncard[1];
	    }
	}
    } else {
	++ier[3];
	if (ier[3] <= mpserr) {
	    s_wsfi(&io___152);
	    do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
	    do_fio(&c__1, key, (ftnlen)4);
	    do_fio(&c__3, id, (ftnlen)8);
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)100);
	}
	goto L30;
    }
/*     ------------------------------------------------------------------ */
/*     Look up the row name  id(1)  in the hash table. */
/*     ------------------------------------------------------------------ */
    s3hash_(lenh, maxm, ncoll, id, &c__2, &keynam[1], names + (*lrow << 3), &
	    ia, &found, (ftnlen)8, (ftnlen)8);
/*     Error if the row name was already there. */
/*     Otherwise, enter the new name into the hash table. */
    if (found) {
	++ier[4];
	if (ier[4] <= mpserr) {
	    s_wsfi(&io___155);
	    do_fio(&c__3, id, (ftnlen)8);
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)100);
	}
    } else {
	++(*m);
	if (*m <= *maxm) {
	    jrow = *maxn + *m;
	    keynam[ia] = *m;
	    s_copy(names + (jrow << 3), id, (ftnlen)8, (ftnlen)8);
	    hrtype[*m] = it;
	}
    }
    goto L30;
/*     ================================================================== */
/*     Should be COLUMNS card. */
/*     ================================================================== */
L110:
    if (s_cmp(key, lcolu, (ftnlen)4, (ftnlen)4) != 0) {
	++ier[1];
	snprnt_(&c__3, " XXXX  COLUMNS card not found", &iw[1], leniw, (
		ftnlen)29);
    }
/*     Error if no rows or too many rows. */
    if (*m <= 0) {
	snprnt_(&c__3, " XXXX  No rows specified", &iw[1], leniw, (ftnlen)24);
	++ier[1];
	*iexit = 113;
	goto L900;
    } else if (*m > *maxm) {
	s_wsfi(&io___157);
	do_fio(&c__1, (char *)&(*maxm), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)100);
	++ier[1];
	*iexit = 112;
	goto L900;
    }
/*     Warning if no objective row found. */
/*     Error if linear objective is ahead of nonlinear rows. */
    if (*iobj == 0) {
	snprnt_(&c__3, " ===>  Warning - no linear objective selected", &iw[1]
		, leniw, (ftnlen)45);
    } else if (*iobj <= *nncon) {
/* Writing concatenation */
	i__1[0] = 40, a__1[0] = " XXXX  The linear objective card      N ";
	i__1[1] = 8, a__1[1] = cw + 416;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)48);
	snprnt_(&c__13, ch__1, &iw[1], leniw, (ftnlen)48);
	snprnt_(&c__3, " XXXX  is out of place.    Nonlinear constraints", &
		iw[1], leniw, (ftnlen)48);
	snprnt_(&c__3, " XXXX  must be listed first in the ROWS section.", &
		iw[1], leniw, (ftnlen)48);
	*iexit = 113;
	goto L900;
    }
    s_wsfi(&io___158);
    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)100);
/*     ------------------------------------------------------------------ */
/*     Exit */
/*     ------------------------------------------------------------------ */
L900:
    return 0;
/* L1000: */
} /* s3mpsa_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3mpsa */
/* Subroutine */ int s3mpsb_(integer *iexit, integer *ier, integer *line, 
	integer *maxn, integer *maxne, integer *lrow, integer *lennm, integer 
	*lenh, integer *ncoll, integer *iobj, doublereal *objadd, integer *m, 
	integer *n, integer *ne, integer *nnl, integer *nncon, integer *nnjac,
	 integer *negcon, integer *na0, char *key, integer *ncard, integer *
	hrtype, integer *keynam, char *names, integer *loca, integer *inda, 
	doublereal *acol, doublereal *bl, doublereal *bu, integer *kbs, 
	doublereal *pi, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen key_len, ftnlen names_len, 
	ftnlen cw_len)
{
    /* Initialized data */

    static char cblank[8] = "        ";
    static char lrhs[4] = "RHS ";
    static char lrhsx[4] = "RHS'";
    static char lrang[4] = "RANG";
    static char llagr[8] = "LAGRANGE";

    /* Format strings */
    static char fmt_1420[] = "(\002 XXXX  Column  \002,a8,\002  has more tha"
	    "n one entry\002,\002 in row  \002,a8/\002 XXXX  Coefficient\002,"
	    "1p,e15.5,\002  ignored in line\002,i10)";
    static char fmt_1400[] = "(\002 XXXX  Non-existent row    specified --"
	    " \002,a8,\002 -- entry ignored in line\002,i7)";
    static char fmt_1500[] = "(\002 XXXX  No valid row entries in column "
	    " \002,a8)";
    static char fmt_3040[] = "(\002 XXXX  Too many columns.   The limit wa"
	    "s\002,i8,4x,\002  Actual number is\002,i8)";
    static char fmt_3050[] = "(\002 XXXX  Too many elements.  The limit wa"
	    "s\002,i8,4x,\002  Actual number is\002,i8)";
    static char fmt_5200[] = "(\002 Columns\002,i8)";
    static char fmt_5210[] = "(\002 Elements\002,i7)";
    static char fmt_1720[] = "(\002 ===>  Warning - first RHS is LAGRANGE"
	    ".\002,\002   Other RHS's will be ignored.\002)";
    static char fmt_1630[] = "(\002 ===>  Note:  constant\002,1p,e15.7,\002 "
	    " is added to the objective.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfi(icilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfi(void);

    /* Local variables */
    static integer i__, k, ia;
    static char id[8*3];
    static integer ne1;
    static doublereal aij, bnd;
    static char str[100];
    static integer ljac;
    static doublereal arng, brng;
    static integer imps, mlst, irow;
    static doublereal aelem[2];
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), iload_(integer *, integer *, integer *, integer *);
    static integer ldenj;
    static logical dense;
    static char colnm[8];
    static logical found, gotnm;
    static char rownm[8];
    extern /* Subroutine */ int s3read_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, char *, doublereal *, 
	    integer *, ftnlen, ftnlen), s3hash_(integer *, integer *, integer 
	    *, char *, integer *, integer *, char *, integer *, logical *, 
	    ftnlen, ftnlen);
    static doublereal infbnd;
    static integer jslack;
    static doublereal aijtol, biglow, bigupp;
    static integer inform__;
    static doublereal bstruc[2];
    static integer mpserr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___187 = { 0, str, 0, fmt_1420, 100, 1 };
    static icilist io___188 = { 0, str, 0, fmt_1400, 100, 1 };
    static icilist io___189 = { 0, str, 0, fmt_1500, 100, 1 };
    static icilist io___190 = { 0, str, 0, fmt_3040, 100, 1 };
    static icilist io___191 = { 0, str, 0, fmt_3050, 100, 1 };
    static icilist io___192 = { 0, str, 0, fmt_5200, 100, 1 };
    static icilist io___193 = { 0, str, 0, fmt_5210, 100, 1 };
    static icilist io___198 = { 0, str, 0, fmt_1400, 100, 1 };
    static icilist io___199 = { 0, str, 0, fmt_1720, 100, 1 };
    static icilist io___200 = { 0, str, 0, fmt_1400, 100, 1 };
    static icilist io___201 = { 0, str, 0, fmt_1630, 100, 1 };
    static icilist io___204 = { 0, str, 0, fmt_1400, 100, 1 };


/*     ================================================================== */
/*     s3mpsb inputs the COLUMNS, RHS and RANGES sections of an MPS file. */

/*     Original version written by Keith Morris, Wellington, 1973. */

/*     19 Jul 1997: Thread-safe version. */
/*     24 Sep 1997: key and Names(*) are now character*8. */
/*     03 Aug 2003: snPRNT adopted. */
/*     27 Oct 2003: Current version of s3mpsb. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* rhs name */
/*     ------------------------------------------------------------------ */
/* range name */
    /* Parameter adjustments */
    --ier;
    --loca;
    --acol;
    --inda;
    names -= 8;
    --keynam;
    --pi;
    --kbs;
    --bu;
    --bl;
    --hrtype;
    --ncard;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    imps = iw[123];
/* MPS file */
    ldenj = iw[105];
/* 1(2)    => dense(sparse) Jacobian */
    mpserr = iw[106];
/* maximum # errors in MPS data */
    mlst = iw[107];
/* maximum # lines  of MPS data */
    infbnd = rw[70];
/* definition of an infinite bound */
    aijtol = rw[95];
/* zero Aij tolerance. */
    bstruc[0] = rw[96];
/* default lower bound on x */
    bstruc[1] = rw[97];
/* default upper bound on x */
    *objadd = 0.;
    dense = ldenj == 1;
    bigupp = infbnd;
    biglow = -bigupp;
    s_copy(colnm, "12345678", (ftnlen)8, (ftnlen)8);
    *n = 0;
    *na0 = 0;
    *ne = 0;
    ne1 = -1;
    *negcon = 0;
    inform__ = 0;
    iload_(m, &c__0, &kbs[1], &c__1);
    dload_(m, &c_b257, &pi[1], &c__1);
/*     ================================================================== */
/*     Read the next columns card. */
/*     ================================================================== */
L210:
    s3read_(&c__2, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L310;
    }
L220:
    if (s_cmp(id, colnm, (ftnlen)8, (ftnlen)8) != 0) {
/*        Start a new column. */
	if (*ne <= ne1) {
	    goto L310;
	}
	++(*n);
	ne1 = *ne;
	s_copy(colnm, id, (ftnlen)8, (ftnlen)8);
	if (*n <= *maxn) {
	    loca[*n] = *ne + 1;
	    s_copy(names + (*n << 3), colnm, (ftnlen)8, (ftnlen)8);
/*           Make room for a dense Jacobian column. */
	    if (*nncon > 0) {
		ljac = *ne;
		if (dense && *n <= *nnjac) {
		    *ne += *nncon;
		    if (*ne <= *maxne) {
			*ne -= *nncon;
			i__1 = *nncon;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    ++(*ne);
			    inda[*ne] = i__;
			    acol[*ne] = 0.;
			}
		    }
		}
	    }
	}
    }
/*     Process two row names and values. */
    for (i__ = 1; i__ <= 2; ++i__) {
/*        Check for only one on the card. */
	s_copy(rownm, id + (i__ << 3), (ftnlen)8, (ftnlen)8);
	if (s_cmp(rownm, cblank, (ftnlen)8, (ftnlen)8) != 0) {
/*           Look up the row name. */
	    s3hash_(lenh, m, ncoll, rownm, &c__1, &keynam[1], names + (*lrow 
		    << 3), &ia, &found, (ftnlen)8, (ftnlen)8);
	    if (found) {
		aij = aelem[i__ - 1];
		irow = keynam[ia];
/*              Test for a duplicate entry. */
		if (kbs[irow] == *n) {
		    ++ier[8];
		    if (ier[8] <= mpserr) {
			s_wsfi(&io___187);
			do_fio(&c__1, colnm, (ftnlen)8);
			do_fio(&c__1, rownm, (ftnlen)8);
			do_fio(&c__1, (char *)&aij, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(
				integer));
			e_wsfi();
			snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
		    }
		    goto L260;
		}
		kbs[irow] = *n;
		if (irow <= *nncon && *n <= *nnjac) {
/*                 Deal with Jacobian elements. */
		    ++(*negcon);
		    if (dense) {
			acol[ne1 + irow] = aij;
			goto L260;
		    }
/*                 Sparse Jacobian -- make sure the new element is */
/*                 squeezed in ahead of any linear-constraint elements. */
		    ++ljac;
		    if (ljac <= *ne) {
			aij = acol[ljac];
			irow = inda[ljac];
			acol[ljac] = aelem[i__ - 1];
			inda[ljac] = keynam[ia];
		    }
		} else if (abs(aij) < aijtol) {
/*                 Ignore small aijs. */
		    ++(*na0);
		    goto L260;
		}
/*              Pack the nonzero. */
		++(*ne);
		if (*ne <= *maxne) {
		    inda[*ne] = irow;
		    acol[*ne] = aij;
		}
	    } else {
		++ier[5];
		if (ier[5] <= mpserr) {
		    s_wsfi(&io___188);
		    do_fio(&c__1, rownm, (ftnlen)8);
		    do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
		}
	    }
	}
/* RowNm ne Blank */
L260:
	;
    }
    goto L210;
/*     Test for an empty column. */
L310:
    if (*ne <= ne1) {
/*        Column with no rows.   Warning unless variable is nonlinear. */
/*        Insert dummy column with zero in first row. */
	if (*n > *nnl) {
	    ++ier[6];
	    if (ier[6] <= mpserr) {
		s_wsfi(&io___189);
		do_fio(&c__1, colnm, (ftnlen)8);
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
	    }
	}
	++(*ne);
	if (*ne <= *maxne) {
	    inda[*ne] = 1;
	    acol[*ne] = 0.;
	}
	if (inform__ == 0) {
	    goto L220;
	}
    }
/*     ================================================================== */
/*     See if we have hit the RHS. */
/*     ================================================================== */
    if (s_cmp(key, lrhs, (ftnlen)4, (ftnlen)4) != 0 && s_cmp(key, lrhsx, (
	    ftnlen)4, (ftnlen)4) != 0) {
/*        Nope sumpins rong. */
/*        Terminate the COLUMNS section anyway. */
	++ier[7];
	snprnt_(&c__3, " XXXX  RHS card not found", &iw[1], leniw, (ftnlen)25)
		;
    }
/*     Are there any columns at all? */
/*     Or too many columns or elements? */
    if (*n <= 0) {
	snprnt_(&c__3, " XXXX  No columns specified", &iw[1], leniw, (ftnlen)
		27);
	++ier[2];
	*iexit = 113;
	return 0;
    } else if (*n > *maxn) {
	s_wsfi(&io___190);
	do_fio(&c__1, (char *)&(*maxn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)100);
	++ier[2];
	*iexit = 112;
	return 0;
    } else if (*ne > *maxne) {
	s_wsfi(&io___191);
	do_fio(&c__1, (char *)&(*maxne), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)100);
	++ier[2];
	*iexit = 112;
	return 0;
    }
/*     ------------------------------------------------------------------ */
/*     Input the RHS. */
/*     ------------------------------------------------------------------ */
    s_wsfi(&io___192);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)100);
    s_wsfi(&io___193);
    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__2, str, &iw[1], leniw, (ftnlen)100);
/*     We finally know how big the problem is. */
    loca[*n + 1] = *ne + 1;
/*     Set bounds to default values. */
    dload_(n, bstruc, &bl[1], &c__1);
    dload_(n, &bstruc[1], &bu[1], &c__1);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = hrtype[i__];
	jslack = *n + i__;
	if (k <= 0) {
	    bl[jslack] = 0.;
	}
	if (k < 0) {
	    bu[jslack] = bigupp;
	}
	if (k > 0) {
	    bl[jslack] = biglow;
	}
	if (k >= 0) {
	    bu[jslack] = 0.;
	}
	if (k == 2) {
	    bu[jslack] = bigupp;
	}
    }
/*     Check for no RHS. */
    if (s_cmp(key, lrhs, (ftnlen)4, (ftnlen)4) != 0 && s_cmp(key, lrhsx, (
	    ftnlen)4, (ftnlen)4) != 0) {
	goto L600;
    }
    gotnm = s_cmp(cw + 424, cblank, (ftnlen)8, (ftnlen)8) != 0;
    inform__ = 0;
/*     ================================================================== */
/*     Read next RHS card and see if it is the one we want. */
/*     ================================================================== */
L420:
    s3read_(&c__2, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L600;
    }
/*     A normal RHS is terminated if LAGRANGE is found. */
    if (s_cmp(id, llagr, (ftnlen)8, (ftnlen)8) == 0) {
	goto L490;
    }
    if (! gotnm) {
	gotnm = TRUE_;
	s_copy(cw + 424, id, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(id, cw + 424, (ftnlen)8, (ftnlen)8) == 0) {
/*        Look at both halves of the record. */
	for (i__ = 1; i__ <= 2; ++i__) {
	    s_copy(rownm, id + (i__ << 3), (ftnlen)8, (ftnlen)8);
	    if (s_cmp(rownm, cblank, (ftnlen)8, (ftnlen)8) != 0) {
		s3hash_(lenh, m, ncoll, rownm, &c__1, &keynam[1], names + (*
			lrow << 3), &ia, &found, (ftnlen)8, (ftnlen)8);
		if (found) {
		    ++ncard[2];
		    bnd = aelem[i__ - 1];
		    irow = keynam[ia];
		    jslack = *n + irow;
		    k = hrtype[irow];
		    if (irow == *iobj) {
			*objadd = bnd;
		    } else if (k != 2) {
			if (k <= 0) {
			    bl[jslack] = bnd;
			}
			if (k >= 0) {
			    bu[jslack] = bnd;
			}
		    }
		} else {
		    ++ier[5];
		    if (ier[5] <= mpserr) {
			s_wsfi(&io___198);
			do_fio(&c__1, rownm, (ftnlen)8);
			do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(
				integer));
			e_wsfi();
			snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
		    }
		}
	    }
/* if RowNm ne Blank */
	}
    }
    goto L420;
/*     LAGRANGE RHS found. */
L490:
    if (ncard[2] == 0) {
	s_copy(cw + 424, cblank, (ftnlen)8, (ftnlen)8);
	s_wsfi(&io___199);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    }
    goto L520;
/*     ================================================================== */
/*     Read next RHS card and see if it is a LAGRANGE one. */
/*     ================================================================== */
L510:
    s3read_(&c__2, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L600;
    }
    if (s_cmp(id, llagr, (ftnlen)8, (ftnlen)8) != 0) {
	goto L510;
    }
/*     Find which row. */
/*     Look at both halves of the record. */
L520:
    for (i__ = 1; i__ <= 2; ++i__) {
	s_copy(rownm, id + (i__ << 3), (ftnlen)8, (ftnlen)8);
	if (s_cmp(rownm, cblank, (ftnlen)8, (ftnlen)8) != 0) {
	    s3hash_(lenh, m, ncoll, rownm, &c__1, &keynam[1], names + (*lrow 
		    << 3), &ia, &found, (ftnlen)8, (ftnlen)8);
	    if (found) {
		++ncard[5];
		irow = keynam[ia];
		pi[irow] = aelem[i__ - 1];
	    } else {
		++ier[5];
		if (ier[5] <= mpserr) {
		    s_wsfi(&io___200);
		    do_fio(&c__1, rownm, (ftnlen)8);
		    do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
		}
	    }
	}
/* RowNm ne Blank */
    }
    goto L510;
/*     ------------------------------------------------------------------ */
/*     RHS has been input. */
/*     ------------------------------------------------------------------ */
L600:
    if (ncard[2] == 0) {
	snprnt_(&c__3, " ===>  Warning - the RHS is zero", &iw[1], leniw, (
		ftnlen)32);
    }
    if (*objadd != 0.) {
	s_wsfi(&io___201);
	do_fio(&c__1, (char *)&(*objadd), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    }
/*     ------------------------------------------------------------------ */
/*     Input RANGES. */
/*     ------------------------------------------------------------------ */
/*     Check for no RANGES. */
    if (s_cmp(key, lrang, (ftnlen)4, (ftnlen)4) != 0) {
	goto L800;
    }
    gotnm = s_cmp(cw + 432, cblank, (ftnlen)8, (ftnlen)8) != 0;
    inform__ = 0;
/*     ================================================================== */
/*     Read card and see if it is the range we want. */
/*     ================================================================== */
L610:
    s3read_(&c__2, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L800;
    }
    if (! gotnm) {
	gotnm = TRUE_;
	s_copy(cw + 432, id, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(id, cw + 432, (ftnlen)8, (ftnlen)8) == 0) {
/*        Look at both halves of the record. */
	for (i__ = 1; i__ <= 2; ++i__) {
	    s_copy(rownm, id + (i__ << 3), (ftnlen)8, (ftnlen)8);
	    if (s_cmp(rownm, cblank, (ftnlen)8, (ftnlen)8) != 0) {
		s3hash_(lenh, m, ncoll, rownm, &c__1, &keynam[1], names + (*
			lrow << 3), &ia, &found, (ftnlen)8, (ftnlen)8);
		if (found) {
		    ++ncard[3];
		    brng = aelem[i__ - 1];
		    arng = abs(brng);
		    irow = keynam[ia];
		    jslack = *n + irow;
		    k = hrtype[irow];
		    if (k == 0) {
			if (brng > 0.) {
			    bu[jslack] = bl[jslack] + arng;
			}
			if (brng < 0.) {
			    bl[jslack] = bu[jslack] - arng;
			}
		    } else if (k == 2) {
/*                    Relax */
		    } else if (k < 0) {
			bu[jslack] = bl[jslack] + arng;
		    } else if (k > 0) {
			bl[jslack] = bu[jslack] - arng;
		    }
		} else {
		    ++ier[5];
		    if (ier[5] <= mpserr) {
			s_wsfi(&io___204);
			do_fio(&c__1, rownm, (ftnlen)8);
			do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(
				integer));
			e_wsfi();
			snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
		    }
		}
	    }
/* RowNm ne Blank */
/* L640: */
	}
    }
    goto L610;
/*     RANGES have been input. */
L800:
    return 0;
} /* s3mpsb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3mpsb */
/* Subroutine */ int s3mpsc_(integer *ier, integer *line, integer *m, integer 
	*n, integer *ns, integer *lennm, char *key, integer *ncard, char *
	names, doublereal *bl, doublereal *bu, integer *hs, doublereal *x, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen key_len, ftnlen names_len, ftnlen cw_len)
{
    /* Initialized data */

    static char cblank[8] = "        ";
    static char linit[8] = "INITIAL ";
    static char lboun[4] = "BOUN";
    static char lenda[4] = "ENDA";
    static char lfr[4] = " FR ";
    static char lfx[4] = " FX ";
    static char llo[4] = " LO ";
    static char lmi[4] = " MI ";
    static char lpl[4] = " PL ";
    static char lup[4] = " UP ";
    static char objtyp[8*3] = "Max     " "Feas    " "Min     ";

    /* Format strings */
    static char fmt_1400[] = "(\002 XXXX  Non-existent column specified --"
	    " \002,a8,\002 -- entry ignored in line\002,i7)";
    static char fmt_1700[] = "(\002 XXXX  Illegal bound type at line\002,i7"
	    ",\002... \002,a4,a8,2x,a8)";
    static char fmt_1720[] = "(\002 ===>  Warning - first bounds set is  INI"
	    "TIAL .\002,\002   Other bounds will be ignored.\002)";
    static char fmt_1740[] = "(\002 XXXX  Bounds back to front on column\002"
	    ",i6,\002 :\002,1p,2e15.5)";
    static char fmt_1900[] = "(\002 XXXX  Total no. of errors in MPS file"
	    "\002,i6)";
    static char fmt_2100[] = "(\002 Objective\002,6x,a8,\002 (\002,a3,\002"
	    ")\002,i8)";
    static char fmt_2110[] = "(\002 RHS      \002,6x,a8,i14)";
    static char fmt_2120[] = "(\002 RANGES   \002,6x,a8,i14)";
    static char fmt_2130[] = "(\002 BOUNDS   \002,6x,a8,i14)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k;
    static doublereal b1, b2;
    static char id[8*3];
    static integer js;
    static doublereal bnd;
    static char str[100];
    static integer imps, mlst;
    static doublereal aelem[2];
    static integer ioldb, jmark;
    static logical gotnm;
    extern /* Subroutine */ int s3read_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, char *, char *, doublereal *, 
	    integer *, ftnlen, ftnlen), s4name_(integer *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static integer iloadb;
    static doublereal infbnd;
    static logical ignore;
    static doublereal biglow, bigupp;
    static integer minmax, inform__, iinsrt, mpserr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___234 = { 0, str, 0, fmt_1400, 100, 1 };
    static icilist io___235 = { 0, str, 0, fmt_1700, 100, 1 };
    static icilist io___237 = { 0, str, 0, fmt_1720, 100, 1 };
    static icilist io___239 = { 0, str, 0, fmt_1400, 100, 1 };
    static icilist io___241 = { 0, str, 0, fmt_1700, 100, 1 };
    static icilist io___244 = { 0, str, 0, fmt_1740, 100, 1 };
    static icilist io___246 = { 0, str, 0, fmt_1900, 100, 1 };
    static icilist io___247 = { 0, str, 0, fmt_2100, 100, 1 };
    static icilist io___248 = { 0, str, 0, fmt_2110, 100, 1 };
    static icilist io___249 = { 0, str, 0, fmt_2120, 100, 1 };
    static icilist io___250 = { 0, str, 0, fmt_2130, 100, 1 };


/*     ================================================================== */
/*     s3mpsc inputs the BOUNDS section of an MPS file. */

/*     Original version written by Keith Morris, Wellington, 1973. */

/*     19 Jul 1997: Thread-safe version. */
/*     24 Sep 1997: key and Names(*) are now character*8. */
/*     01 Aug 2003: s4name now needs iw, leniw for snPRNT. */
/*     03 Aug 2003: snPRNT adopted. */
/*     08 Oct 2003: Current version of s3mpsc. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ier;
    --x;
    --hs;
    --bu;
    --bl;
    names -= 8;
    --ncard;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    minmax = iw[87];
/* 1, 0, -1  => MIN, FP, MAX */
    mpserr = iw[106];
/* maximum # errors in MPS data */
    mlst = iw[107];
/* maximum # lines  of MPS data */
    iloadb = iw[122];
/* load file */
    imps = iw[123];
/* MPS file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    infbnd = rw[70];
/* definition of an infinite bound */
    inform__ = 1;
    bigupp = infbnd;
    biglow = -bigupp;
/*     Check for no BOUNDS. */
    if (s_cmp(key, lboun, (ftnlen)4, (ftnlen)4) != 0) {
	goto L700;
    }
    gotnm = s_cmp(cw + 440, cblank, (ftnlen)8, (ftnlen)8) != 0;
    inform__ = 0;
    jmark = 1;
/*     ================================================================== */
/*     Read and check BOUNDS cards.  Notice the double plural. */
/*     ================================================================== */
L610:
    s3read_(&c__3, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L700;
    }
/*     A normal bounds set is terminated if INITIAL is found. */
    bnd = aelem[0];
    if (s_cmp(id, linit, (ftnlen)8, (ftnlen)8) == 0) {
	goto L690;
    }
    if (! gotnm) {
	gotnm = TRUE_;
	s_copy(cw + 440, id, (ftnlen)8, (ftnlen)8);
    }
    if (s_cmp(id, cw + 440, (ftnlen)8, (ftnlen)8) == 0) {
/*        Find which column. */
	s4name_(n, names + 8, id + 8, line, &ier[10], &c__0, &c__1, n, &jmark,
		 &j, &iw[1], leniw, (ftnlen)8, (ftnlen)8);
	if (j <= 0) {
	    if (ier[10] <= mpserr) {
		s_wsfi(&io___234);
		do_fio(&c__1, id + 8, (ftnlen)8);
		do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
	    }
	} else {
/*           Select bound type for column j. */
	    ++ncard[4];
	    if (s_cmp(key, lup, (ftnlen)4, (ftnlen)4) == 0) {
		bu[j] = bnd;
	    } else if (s_cmp(key, llo, (ftnlen)4, (ftnlen)4) == 0) {
		bl[j] = bnd;
	    } else if (s_cmp(key, lfx, (ftnlen)4, (ftnlen)4) == 0) {
		bu[j] = bnd;
		bl[j] = bnd;
	    } else if (s_cmp(key, lfr, (ftnlen)4, (ftnlen)4) == 0) {
		bu[j] = bigupp;
		bl[j] = biglow;
	    } else if (s_cmp(key, lmi, (ftnlen)4, (ftnlen)4) == 0) {
		if (bu[j] >= bigupp) {
		    bu[j] = 0.;
		}
		bl[j] = biglow;
	    } else if (s_cmp(key, lpl, (ftnlen)4, (ftnlen)4) == 0) {
		bu[j] = bigupp;
	    } else {
/*              This lad didn't even make it to Form 1. */
		++ier[11];
		if (ier[11] <= mpserr) {
		    s_wsfi(&io___235);
		    do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
		    do_fio(&c__1, key, (ftnlen)4);
		    for (i__ = 1; i__ <= 2; ++i__) {
			do_fio(&c__1, id + (i__ - 1 << 3), (ftnlen)8);
		    }
		    e_wsfi();
		    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
		}
	    }
	}
    }
    goto L610;
/*     INITIAL bounds set found. */
L690:
    if (ncard[4] == 0) {
	s_copy(cw + 440, cblank, (ftnlen)8, (ftnlen)8);
	s_wsfi(&io___237);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    }
/*     ------------------------------------------------------------------ */
/*     End of normal bounds. */
/*     ------------------------------------------------------------------ */
L700:
    *ns = 0;
    bigupp *= .9;
    biglow = -bigupp;
/*     Set variables to be nonbasic at zero (as long as that's feasible). */
/*     All variables will be eligible for the initial basis. */
    i__1 = *n + *m;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = 0., d__2 = bl[j];
	x[j] = max(d__1,d__2);
/* Computing MIN */
	d__1 = x[j], d__2 = bu[j];
	x[j] = min(d__1,d__2);
	hs[j] = 0;
	if (x[j] == bu[j]) {
	    hs[j] = 1;
	}
    }
/*     Ignore INITIAL bounds if a basis will be loaded. */
    if (inform__ != 0) {
	goto L790;
    }
    ignore = ioldb > 0 || iinsrt > 0 || iloadb > 0;
    if (! ignore) {
	jmark = 1;
	goto L720;
    }
/*     ================================================================== */
/*     Read INITIAL bounds set. */
/*     ================================================================== */
L710:
    s3read_(&c__3, &imps, &iw[1], leniw, line, &mlst, key, id, aelem, &
	    inform__, (ftnlen)4, (ftnlen)8);
    if (inform__ != 0) {
	goto L790;
    }
    bnd = aelem[0];
    if (ignore || s_cmp(id, linit, (ftnlen)8, (ftnlen)8) != 0) {
	goto L710;
    }
/*     Find which column. */
L720:
    s4name_(n, names + 8, id + 8, line, &ier[12], &c__0, &c__1, n, &jmark, &j,
	     &iw[1], leniw, (ftnlen)8, (ftnlen)8);
    if (j <= 0) {
	if (ier[12] <= mpserr) {
	    s_wsfi(&io___239);
	    do_fio(&c__1, id + 8, (ftnlen)8);
	    do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
	}
    } else {
/*        Select bound type for column j. */
	++ncard[6];
	if (s_cmp(key, lfr, (ftnlen)4, (ftnlen)4) == 0) {
	    js = -1;
	} else if (s_cmp(key, lfx, (ftnlen)4, (ftnlen)4) == 0) {
	    js = 2;
	    ++(*ns);
	} else if (s_cmp(key, llo, (ftnlen)4, (ftnlen)4) == 0) {
	    js = 4;
	    bnd = bl[j];
	} else if (s_cmp(key, lup, (ftnlen)4, (ftnlen)4) == 0) {
	    js = 5;
	    bnd = bu[j];
	} else if (s_cmp(key, lmi, (ftnlen)4, (ftnlen)4) == 0) {
	    js = 4;
	} else if (s_cmp(key, lpl, (ftnlen)4, (ftnlen)4) == 0) {
	    js = 5;
	} else {
	    ++ier[13];
	    if (ier[13] <= mpserr) {
		s_wsfi(&io___241);
		do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
		do_fio(&c__1, key, (ftnlen)4);
		for (i__ = 1; i__ <= 2; ++i__) {
		    do_fio(&c__1, id + (i__ - 1 << 3), (ftnlen)8);
		}
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
	    }
	    goto L710;
	}
    }
    if (abs(bnd) >= bigupp) {
	bnd = 0.;
    }
    x[j] = bnd;
    hs[j] = js;
    goto L710;
/*     Should be ENDATA card. */
L790:
    if (s_cmp(key, lenda, (ftnlen)4, (ftnlen)4) != 0) {
	ier[14] = 1;
	snprnt_(&c__3, " XXXX  ENDATA card not found", &iw[1], leniw, (ftnlen)
		28);
    }
/*     ------------------------------------------------------------------ */
/*     Pass the Buck - not got to Truman yet. */
/*     ------------------------------------------------------------------ */
/*     Check that  bl .le. bu */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	b1 = bl[j];
	b2 = bu[j];
	if (b1 > b2) {
	    ++ier[20];
	    if (ier[20] <= mpserr) {
		s_wsfi(&io___244);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&b1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&b2, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
	    }
	    bl[j] = b2;
	    bu[j] = b1;
	}
    }
/*     Count the errors. */
    k = 0;
    for (i__ = 1; i__ <= 20; ++i__) {
	k += ier[i__];
    }
    if (k > 0) {
	s_wsfi(&io___246);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)100);
    }
    snprnt_(&c__11, " Names selected", &iw[1], leniw, (ftnlen)15);
    snprnt_(&c__1, " --------------", &iw[1], leniw, (ftnlen)15);
    s_wsfi(&io___247);
    do_fio(&c__1, cw + 416, (ftnlen)8);
    do_fio(&c__1, objtyp + (minmax + 1 << 3), (ftnlen)3);
    do_fio(&c__1, (char *)&ncard[1], (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    s_wsfi(&io___248);
    do_fio(&c__1, cw + 424, (ftnlen)8);
    do_fio(&c__1, (char *)&ncard[2], (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    s_wsfi(&io___249);
    do_fio(&c__1, cw + 432, (ftnlen)8);
    do_fio(&c__1, (char *)&ncard[3], (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    s_wsfi(&io___250);
    do_fio(&c__1, cw + 440, (ftnlen)8);
    do_fio(&c__1, (char *)&ncard[4], (ftnlen)sizeof(integer));
    e_wsfi();
    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    return 0;
} /* s3mpsc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s3mpsc */
/* Subroutine */ int s3read_(integer *mode, integer *imps, integer *iw, 
	integer *leniw, integer *line, integer *mxlist, char *key, char *id, 
	doublereal *aelem, integer *inform__, ftnlen key_len, ftnlen id_len)
{
    /* Initialized data */

    static char lblank[1] = " ";
    static char lstar[1] = "*";

    /* Format strings */
    static char fmt_1000[] = "(a61)";
    static char fmt_2000[] = "(i7,4x,a)";
    static char fmt_1100[] = "(a4,a8,2x,a8,2x,bn,e12.0,3x,a8,2x,e12.0)";

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_wsfi(icilist *), e_wsfi(void), s_rsfi(icilist *), e_rsfi(void);

    /* Local variables */
    static char str[100];
    static integer last;
    static char buff1[1], buffer[61];
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___253 = { 0, 0, 0, fmt_1000, 0 };
    static icilist io___258 = { 0, str, 0, fmt_2000, 100, 1 };
    static icilist io___259 = { 0, buffer, 0, fmt_1100, 61, 1 };
    static icilist io___260 = { 0, buffer, 0, fmt_1100, 61, 1 };
    static icilist io___261 = { 0, buffer, 0, fmt_1100, 61, 1 };
    static icilist io___262 = { 0, buffer, 0, fmt_1100, 61, 1 };


/*     ================================================================== */
/*     s3read  reads data from file iMPS and prints a listing on file */
/*     iPrint.  The data is assumed to be in MPS format, with items of */
/*     interest in the following six fields... */

/*     Field:     1         2         3         4         5         6 */

/*     Columns: 01-04     05-12     15-22     25-36     40-47     50-61 */

/*     Format:    a4        a8        a8      e12.0       a8      e12.0 */

/*     Data:     key      id(1)     id(2)   Aelem(1)    id(3)   Aelem(2) */


/*     Comments may contain a * in column 1 and anything in columns 2-61. */
/*     They are listed and then ignored. */


/*     On entry,  mode    specifies which fields are to be processed. */
/*     On exit ,  inform  is set to 1 if column 1 is not blank. */

/*     15 Nov 1991: First version based on Minos routine m3read. */
/*     24 Sep 1997: key and id(*) are now character*8. */
/*     03 Aug 2003: snPRNT adopted. */
/*     03 Aug 2003: Current version of s3read. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    id -= 8;
    --aelem;

    /* Function Body */
/*     ------------------------------------------------------------------ */
/*     Read a data card and look for keywords and comments. */
/*     ------------------------------------------------------------------ */
L10:
    io___253.ciunit = *imps;
    s_rsfe(&io___253);
    do_fio(&c__1, buffer, (ftnlen)61);
    e_rsfe();
    *(unsigned char *)buff1 = *(unsigned char *)buffer;
    ++(*line);
/*     Print the buffer if column 1 is nonblank */
/*     or if a listing is wanted. */
    if (*(unsigned char *)buff1 != *(unsigned char *)&lblank[0] || *line <= *
	    mxlist) {
/*        Find the last nonblank character. */
	for (last = 61; last >= 2; --last) {
	    if (*(unsigned char *)&buffer[last - 1] != *(unsigned char *)&
		    lblank[0]) {
		goto L30;
	    }
/* L20: */
	}
	last = 1;
L30:
	s_wsfi(&io___258);
	do_fio(&c__1, (char *)&(*line), (ftnlen)sizeof(integer));
	do_fio(&c__1, buffer, last);
	e_wsfi();
	snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)100);
    }
/*     Ignore comments. */
    if (*(unsigned char *)buff1 == *(unsigned char *)&lstar[0]) {
	goto L10;
    }
/*     If column 1 is nonblank, load key and exit. */
/*     The NAME card is unusual in having some data in field 3. */
/*     We have to load it into id(2). */
    if (*(unsigned char *)buff1 != *(unsigned char *)&lblank[0]) {
	s_rsfi(&io___259);
	do_fio(&c__1, key, (ftnlen)4);
	do_fio(&c__1, id + 8, (ftnlen)8);
	do_fio(&c__1, id + 16, (ftnlen)8);
	e_rsfi();
	*inform__ = 1;
	return 0;
    }
/*     ------------------------------------------------------------------ */
/*     Process normal data cards. */
/*     ------------------------------------------------------------------ */
    if (*mode == 1) {
/*        NAME or ROWS sections. */
	s_rsfi(&io___260);
	do_fio(&c__1, key, (ftnlen)4);
	do_fio(&c__1, id + 8, (ftnlen)8);
	do_fio(&c__1, id + 16, (ftnlen)8);
	e_rsfi();
    } else if (*mode == 2) {
/*        COLUMNS, RHS or RANGES sections. */
	s_rsfi(&io___261);
	do_fio(&c__1, key, (ftnlen)4);
	do_fio(&c__1, id + 8, (ftnlen)8);
	do_fio(&c__1, id + 16, (ftnlen)8);
	do_fio(&c__1, (char *)&aelem[1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, id + 24, (ftnlen)8);
	do_fio(&c__1, (char *)&aelem[2], (ftnlen)sizeof(doublereal));
	e_rsfi();
    } else {
/*        BOUNDS section. */
	s_rsfi(&io___262);
	do_fio(&c__1, key, (ftnlen)4);
	do_fio(&c__1, id + 8, (ftnlen)8);
	do_fio(&c__1, id + 16, (ftnlen)8);
	do_fio(&c__1, (char *)&aelem[1], (ftnlen)sizeof(doublereal));
	e_rsfi();
    }
    return 0;
} /* s3read_ */

