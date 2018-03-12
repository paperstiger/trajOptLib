/* ../snopt7/src/sn04wrap.f -- translated by f2c (version 20100827).
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
static integer c__15 = 15;
static integer c__5 = 5;
static integer c__2 = 2;
static integer c__3 = 3;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn04wrap.f */

/*     snEXIT   snWRAP   snSolF */

/*     13 Dec 2013: Three subroutines separated from sn02lib to build */
/*                  standalone SQOPT library. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int snexit_(integer *iexit, char *solver, char *string, char 
	*string2, ftnlen solver_len, ftnlen string_len, ftnlen string2_len)
{
    /* Initialized data */

    static integer indc[15] = { 0,5,11,14,18,23,29,33,38,43,47,55,59,61,66 };
    static char c__[80*70] = "finished successfully                         "
	    "                                  " "optimality conditions satis"
	    "fied                                                 " "feasible"
	    " point found                                                    "
	    "        " "requested accuracy could not be achieved             "
	    "                           " "weak QP minimizer                 "
	    "                                              " "the problem app"
	    "ears to be infeasible                                            "
	     "infeasible linear constraints                                 "
	    "                  " "infeasible linear equalities               "
	    "                                     " "nonlinear infeasibilitie"
	    "s minimized                                             " "infea"
	    "sibilities minimized                                            "
	    "           " "infeasible linear constraints in QP subproblem    "
	    "                              " "the problem appears to be unbou"
	    "nded                                             " "unbounded ob"
	    "jective                                                         "
	    "    " "constraint violation limit reached                       "
	    "                       " "resource limit error                  "
	    "                                          " "iteration limit rea"
	    "ched                                                         " 
	    "major iteration limit reached                                  "
	    "                 " "the superbasics limit is too small          "
	    "                                    " "terminated after numerica"
	    "l difficulties                                         " "curren"
	    "t point cannot be improved                                      "
	    "          " "singular basis                                     "
	    "                             " "cannot satisfy the general const"
	    "raints                                          " "ill-condition"
	    "ed null-space basis                                             "
	    "   " "error in the user-supplied functions                      "
	    "                      " "incorrect objective  derivatives       "
	    "                                         " "incorrect constraint"
	    " derivatives                                                " 
	    "the QP Hessian is indefinite                                   "
	    "                 " "incorrect second derivatives                "
	    "                                    " "incorrect derivatives    "
	    "                                                       " "undefi"
	    "ned user-supplied functions                                     "
	    "          " "undefined function at the first feasible point     "
	    "                             " "undefined function at the initia"
	    "l point                                         " "unable to pro"
	    "ceed into undefined region                                      "
	    "   " "user requested termination                                "
	    "                      " "terminated during function evaluation  "
	    "                                         " "terminated during co"
	    "nstraint evaluation                                         " 
	    "terminated during objective evaluation                         "
	    "                 " "terminated from monitor routine             "
	    "                                    " "insufficient storage allo"
	    "cated                                                  " "work a"
	    "rrays must have at least 500 elements                           "
	    "          " "not enough character storage                       "
	    "                             " "not enough integer storage      "
	    "                                                " "not enough re"
	    "al storage                                                      "
	    "   " "input arguments out of range                              "
	    "                      " "invalid input argument                 "
	    "                                         " "basis file dimension"
	    "s do not match this problem                                 " 
	    "the QP Hessian is indefinite                                   "
	    "                 " "finished successfully                       "
	    "                                    " "SPECS file read          "
	    "                                                       " "Jacobi"
	    "an structure estimated                                          "
	    "          " "MPS file read                                      "
	    "                             " "memory requirements estimated   "
	    "                                                " "user-supplied"
	    " derivatives appear to be correct                               "
	    "   " "no derivatives were checked                               "
	    "                      " "some SPECS keywords were not recognized"
	    "                                         " "errors while process"
	    "ing MPS data                                                " 
	    "no MPS file specified                                          "
	    "                 " "problem-size estimates too small            "
	    "                                    " "fatal error in the MPS fi"
	    "le                                                     " "errors"
	    " while estimating Jacobian structure                            "
	    "          " "cannot find Jacobian structure at given point      "
	    "                             " "fatal errors while reading the S"
	    "PECS                                            " "no SPECS file"
	    " (iSpecs le 0 or iSpecs gt 99)                                  "
	    "   " "End-of-file while looking for a BEGIN                     "
	    "                      " "End-of-file while reading SPECS file   "
	    "                                         " "ENDRUN found before "
	    "any valid SPECS                                             " 
	    "system error                                                   "
	    "                 " "wrong no of basic variables                 "
	    "                                    " "error in basis package   "
	    "                                                       " "Proble"
	    "m dimensions are too large                                      "
	    "          ";

    /* System generated locals */
    integer i__1;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i1, i2, mjr, mnr;
    extern /* Subroutine */ int s1trim_(char *, integer *, ftnlen);
    static integer length;

/*     ================================================================== */
/*     snEXIT  returns the strings associated with EXIT condition iExit. */

/*     On exit, string1 and string2 are trimmed of trailing blanks and */
/*              each have a maximum length of 76 chars. */

/*     25 Sep 2002: First version of snEXIT. */
/*     12 Mar 2004: Trimmed message strings to avoid trouble on SGI's. */
/*     07 May 2006: Second derivative exits added (Exit 54) */
/*     01 Jun 2006: Stand-alone derivative checker added (Exits 55,105) */
/*     01 Jun 2006: Student edition exit added (Exit 66) */
/*     22 Apr 2007: Unrecognized options flagged (Exit 135) */
/*     22 Apr 2007: No derivatives checked (Exit 106) */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Find the "major" and "minor" iExit modes */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
/* EXI */
    mjr = *iexit / 10;
    mnr = *iexit - mjr * 10;
    i1 = indc[mjr];
    i2 = i1 + mnr;
    s1trim_(c__ + i1 * 80, &length, (ftnlen)80);
    ici__1.icierr = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = string_len;
    ici__1.iciunit = string;
    ici__1.icifmt = "(1x,2a,i4,a,(a))";
    s_wsfi(&ici__1);
    do_fio(&c__1, solver, (ftnlen)6);
    do_fio(&c__1, " EXIT", (ftnlen)5);
    i__1 = mjr * 10;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    do_fio(&c__1, " -- ", (ftnlen)4);
    do_fio(&c__1, c__ + i1 * 80, length);
    e_wsfi();
    s1trim_(c__ + i2 * 80, &length, (ftnlen)80);
    ici__1.icierr = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = string2_len;
    ici__1.iciunit = string2;
    ici__1.icifmt = "(1x,2a,i4,a,(a))";
    s_wsfi(&ici__1);
    do_fio(&c__1, solver, (ftnlen)6);
    do_fio(&c__1, " INFO", (ftnlen)5);
    do_fio(&c__1, (char *)&(*iexit), (ftnlen)sizeof(integer));
    do_fio(&c__1, " -- ", (ftnlen)4);
    do_fio(&c__1, c__ + i2 * 80, length);
    e_wsfi();
    return 0;
} /* snexit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snEXIT */
/* Subroutine */ int snwrap_(integer *iexit, char *solver, char *string, char 
	*string2, integer *iw, integer *leniw, ftnlen solver_len, ftnlen 
	string_len, ftnlen string2_len)
{
    extern /* Subroutine */ int s1page_(integer *, integer *, integer *);
    static integer iexit0;
    extern /* Subroutine */ int snexit_(integer *, char *, char *, char *, 
	    ftnlen, ftnlen, ftnlen), snprnt_(integer *, char *, integer *, 
	    integer *, ftnlen);

/*     ================================================================== */
/*     snWRAP  it's a wrap! */

/*     18 Oct 2003: First version of snWRAP. */
/*     21 Jun 2004: iExit saved for use in snSolF. */
/*     21 Jun 2004: Current version of snWRAP. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    iw[424] = *iexit;
/* INFO code from all solvers */
    snexit_(iexit, solver, string, string2, (ftnlen)6, string_len, 
	    string2_len);
    if (*iexit == 81) {
/* Print without using accessing iw, etc. */
	snprnt_(&c__15, string, &iw[1], leniw, string_len);
	snprnt_(&c__5, string2, &iw[1], leniw, string2_len);
    } else {
	iexit0 = *iexit / 10;
	if (iexit0 == 0) {
/* Normal exit */
	    s1page_(&c__1, &iw[1], leniw);
	} else {
	    s1page_(&c__2, &iw[1], leniw);
	}
	snprnt_(&c__3, string, &iw[1], leniw, string_len);
	snprnt_(&c__3, string2, &iw[1], leniw, string2_len);
    }
    return 0;
} /* snwrap_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine snWRAP */
/* Subroutine */ int snsolf_(integer *m, integer *n, integer *nb, integer *
	ninf, integer *j, integer *jkey, integer *jstate, integer *hs, 
	doublereal *bl, doublereal *bu, doublereal *rc, doublereal *xs, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    static doublereal b1, b2, d1, d2, dj;
    static integer js;
    static doublereal xj, tolx;
    static integer iexit, jbinf1, jdinf1;
    static logical feasbl;
    static doublereal tolfea, djtest;
    static integer minimz;
    static logical maximz;
    static doublereal pinorm, tolnlp, tolopt;

/*     ================================================================== */
/*     snSolF sets the solution flags for the j-th variable: */
/*                      ' ' A  D  I  N    and   LL  UL SBS  BS  EQ  FR */
/*     by returning jkey=0  1  2  3  4,  jstate= 0   1   2   3   4   5 */

/*     snSolF is called by SNOPT from s4soln. */
/*     snSolF may also be called externally (e.g. by GAMS) */
/*     following a normal call of SNOPT. */
/*     At this stage the solution will be UNSCALED!! */
/*     Hence, SNOPT (via m4soln) now outputs flags for the UNSCALED soln. */

/*     Input parameters m, n, nb, hs, bl, bu, rc, xs */
/*     are the same as for SNOPT. */

/*     j      (input ) is column j if j <= n;  otherwise row i = j - n. */
/*     jkey   (output) is one of 0 1 2 3 4. */
/*     jstate (output) is one of 0 1 2 3 4 5. */

/*     09 Mar 2004: First version of snSolF, derived from misolf. */
/*     18 Jun 2004: If the scaled problem was infeasible */
/*                  (with max inf at j = jbInf1), always flag that j */
/*                  as infeasible in the unscaled solution. */
/*     21 Jun 2004: Similarly, if the scaled problem wasn't optimal, */
/*                  (with max dual inf at j = jdInf1), always flag that j */
/*                  as nonoptimal in the unscaled solution. */
/*     23 Jun 2004: Suppress nonoptimal flag if iExit <= 2 (not 0). */
/*     ================================================================== */
    /* Parameter adjustments */
    --xs;
    --rc;
    --bu;
    --bl;
    --hs;
    --iw;
    --rw;

    /* Function Body */
    minimz = iw[199];
/* (-1)(+1)    => (max)(min) */
    iexit = iw[424];
/* INFO code from all solvers */
    jbinf1 = iw[427];
/* Largest bound infeasibility (  scaled) */
    jdinf1 = iw[428];
/* Largest dual  infeasibility (  scaled) */
    tolnlp = rw[53];
/* Major Optimality tolerance */
    tolx = rw[56];
/* Minor feasibility tolerance */
    pinorm = rw[422];
/* Lagrange multiplier norm */
    tolfea = tolx;
    tolopt = tolnlp * pinorm;
    feasbl = *ninf == 0;
    maximz = minimz < 0;
    js = hs[*j];
    b1 = bl[*j];
    b2 = bu[*j];
    xj = xs[*j];
    dj = rc[*j];
    d1 = b1 - xj;
    d2 = xj - b2;
    djtest = -dj;
    if (feasbl) {
	if (maximz) {
	    djtest = -djtest;
	}
	jbinf1 = 0;
    }
    if (iexit <= 2) {
	jdinf1 = 0;
    }
/* Set keys and states. */
    *jkey = 0;
/* blank */
    *jstate = js;
/* 0, 1, 2, 3 */
    if (js <= 1) {
/* Nonbasic variables. */
	if (b1 == b2) {
	    *jstate = 4;
	}
	if (-d1 > tolfea && -d2 > tolfea) {
	    *jstate = 5;
	}
	if (*jstate == 1) {
	    djtest = -djtest;
	}
	if (*jstate >= 4) {
	    djtest = abs(djtest);
	}
	if (abs(djtest) <= tolopt) {
	    *jkey = 1;
	}
/* A */
	if (*jstate != 4 && djtest > tolopt) {
	    *jkey = 4;
	}
/* N */
    } else {
/* Basic and superbasic variables. */
	if (abs(d1) <= tolfea || abs(d2) <= tolfea) {
	    *jkey = 2;
	}
/* D */
	if (*jstate == 2 && abs(djtest) > tolopt) {
	    *jkey = 4;
	}
/* N */
	if (d1 > tolfea || d2 > tolfea) {
	    *jkey = 3;
	}
/* I */
	if (*j == jbinf1) {
	    *jkey = 3;
	}
/* I */
    }
    if (*j == jdinf1) {
	*jkey = 4;
    }
/* N */
    return 0;
} /* snsolf_ */

