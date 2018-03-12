/* ./src/sn60srch.f -- translated by f2c (version 20100827).
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

static integer c__3 = 3;
static integer c__1 = 1;
static doublereal c_b36 = -1.;
static integer c__0 = 0;
static doublereal c_b40 = 1.;
static doublereal c_b59 = 0.;
static integer c__23 = 23;
static integer c__21 = 21;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  sn60srch.f */

/*     s6fdG    s6fdG1   s6srch   s6tols   s6usrf */
/*     lsrchc   lsrchq */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s6fdg_(integer *iexit, integer *n, integer *negcon, 
	integer *nncon0, integer *nncon, integer *nnjac, integer *nnl, 
	integer *nnobj0, integer *nnobj, U_fp fgwrap, U_fp fgcon, U_fp fgobj, 
	doublereal *bl, doublereal *bu, doublereal *x, integer *ne, integer *
	nlocj, integer *locj, integer *indj, doublereal *fcon, doublereal *
	fobj, doublereal *gcon, doublereal *gobj, doublereal *y, char *cu, 
	integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer *
	lenru, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    extern /* Subroutine */ int s6fdg1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, U_fp, U_fp,
	     U_fp, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static integer llocg, nlocg, gotfd, lfcon2, lgobju, lgconu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

/*     ================================================================== */
/*     s6fdG   computes any missing objective and constraint gradients */
/*     for the current value of the variables in  x. */

/*     NOTE --  s6dcon overwrites the first  nnCon  elements of  y */
/*     if central differences are needed. */

/*     30 Dec 1991: First version based on Minos 5.4 routine m6grd. */
/*     17 Jul 1997: First thread-safe version. */
/*     11 Oct 1998: s6dcon and s6dobj merged. */
/*     24 Oct 2000: Updated for SNOPT 6.1 */
/*     12 Oct 2003: snEXIT and snPRNT added. */
/*     14 Oct 2003: Enforced feasible perturbation. */
/*     01 Apr 2005: Current version of s6fdG. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --gcon;
    --y;
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
    gotfd = iw[183];
/* > 0 => some differences needed */
    llocg = iw[260];
/* locG(nlocG) = column pointers for indG */
    lfcon2 = iw[318];
/* fCon2(nnCon) work vector */
    lgconu = iw[319];
/* record of unknown derivatives and constants */
    lgobju = iw[323];
/* record of unknown derivatives */
    *iexit = 0;
    nlocg = *nnjac + 1;
    if (gotfd > 0) {
	s6fdg1_(iexit, n, nncon0, nncon, nnjac, nnl, nnobj0, nnobj, (U_fp)
		fgwrap, (U_fp)fgcon, (U_fp)fgobj, &bl[1], &bu[1], &x[1], ne, 
		nlocj, &locj[1], &indj[1], negcon, &nlocg, &iw[llocg], fobj, &
		gobj[1], &fcon[1], &gcon[1], &rw[lgobju], &rw[lfcon2], &rw[
		lgconu], &y[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
		cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (
		ftnlen)8);
	if (*iexit == 63) {
/* The user didn't like some x's */
	    snprnt_(&c__3, " XXX  Unable to apply reversion when differencing"
		    , &iw[1], leniw, (ftnlen)49);
	}
    }
    return 0;
} /* s6fdg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6fdG */
/* Subroutine */ int s6fdg1_(integer *iexit, integer *n, integer *nncon0, 
	integer *nncon, integer *nnjac, integer *nnl, integer *nnobj0, 
	integer *nnobj, S_fp fgwrap, U_fp fgcon, U_fp fgobj, doublereal *bl, 
	doublereal *bu, doublereal *x, integer *ne, integer *nlocj, integer *
	locj, integer *indj, integer *negcon, integer *nlocg, integer *locg, 
	doublereal *fobj, doublereal *gobj, doublereal *fcon, doublereal *
	gcon, doublereal *gobju, doublereal *fconu, doublereal *gconu, 
	doublereal *y, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, l, ir;
    static doublereal xj, buj, dxj;
    static logical done;
    static integer kmax, numf;
    static doublereal tolx, fback;
    static logical needg, needj;
    static doublereal delta, fdint[2];
    static logical found, someg, somej;
    static doublereal infbnd;
    static integer modefg, lvldif;
    static doublereal fforwd;
    static integer lvlder;
    static logical centrl;
    static integer inform__;
    static doublereal gdummy;

/*     ================================================================== */
/*     s6fdG1  estimates missing elements in the objective gradient and */
/*     the columns of the Jacobian using finite differences of the */
/*     problem functions fObj and fCon. */

/*     The arrays y, fConu, gConu  and gObju are used as workspace. */
/*     Dummy elements of gObju and gConu define the unknown derivatives. */

/*     11 Oct 1998: First version based on combining s6dobj and s6dcon. */
/*     24 Oct 2000: Updated for SNOPT 6.1 */
/*     12 Oct 2003: snEXIT and snPRNT added. */
/*     14 Oct 2003: Implemented feasible perturbation. */
/*     16 Jun 2008: Call-status implemented correctly. */
/*     15 Nov 2010: Call-status removed from argument list. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* number of calls of fCon */
/*     ------------------------------------------------------------------ */
/* number of calls of fObj */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --y;
    --fconu;
    --fcon;
    --gobju;
    --gobj;
    --indj;
    --locj;
    --gconu;
    --gcon;
    --locg;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lvlder = iw[70];
/* = 0, 1, 2 or 3, the derivative level */
    lvldif = iw[181];
/*    =1 (2) for forwd (cntrl) diffs */
    tolx = rw[56];
/* Minor feasibility tolerance */
    gdummy = rw[69];
/* definition of an 'unset' value */
    infbnd = rw[70];
/* definition of an infinite bound */
    fdint[0] = rw[76];
/* (1) forwrd diff. interval */
    fdint[1] = rw[77];
/* (2) cntrl  diff. interval */
    *iexit = 0;
/*     The problem functions are called to provide functions only. */
    modefg = 0;
    someg = lvlder == 0 || lvlder == 2;
    somej = lvlder == 0 || lvlder == 1;
    centrl = lvldif == 2;
    delta = fdint[lvldif - 1];
    numf = 0;
    i__1 = *nnl;
    for (j = 1; j <= i__1; ++j) {
/*        Look for the first missing element in this column. */
	found = FALSE_;
	needj = j <= *nnjac && somej;
	needg = j <= *nnobj && someg;
	if (needg) {
	    if (gobju[j] == gdummy) {
		found = TRUE_;
	    }
	}
	if (needj) {
	    l = locg[j];
	    k = locj[j];
	    kmax = locj[j + 1] - 1;
	    done = FALSE_;
/* +          ------------------------------------------------------------ */
/* +          while (k .le. kmax  .and. .not.(found  .or.  done)) do */
L120:
	    if (k <= kmax && ! (found || done)) {
		ir = indj[k];
		if (ir > *nncon) {
		    done = TRUE_;
		} else {
		    if (gconu[l] == gdummy) {
			found = TRUE_;
		    }
		    ++l;
		}
		++k;
		goto L120;
/* +          end while */
/* +          ------------------------------------------------------------ */
	    }
	}
	if (found) {
/*           ------------------------------------------------------------ */
/*           Some missing derivatives for this variable. */
/*           A finite difference is needed. */
/*           ------------------------------------------------------------ */
	    xj = x[j];
	    dxj = delta * (abs(xj) + 1.);
	    buj = bu[j];
	    if (buj < infbnd && xj + dxj + dxj > buj + tolx) {
		dxj = -dxj;
	    }
	    x[j] = xj + dxj;
	    ++numf;
	    (*fgwrap)(&inform__, &modefg, &needj, &needg, n, negcon, nncon0, 
		    nncon, nnjac, nnl, nnobj0, nnobj, (U_fp)fgcon, (U_fp)
		    fgobj, &x[1], ne, nlocj, &locj[1], &indj[1], &fconu[1], &
		    fforwd, &gconu[1], &gobju[1], cu + 8, lencu, &iu[1], 
		    leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1]
		    , lenrw, (ftnlen)8, (ftnlen)8);
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
		} else {
		    *iexit = 63;
/* unable to move into undefined region */
		}
		goto L999;
	    }
	    if (centrl) {
		dxj += dxj;
		x[j] = xj + dxj;
		++numf;
		(*fgwrap)(&inform__, &modefg, &needj, &needg, n, negcon, 
			nncon0, nncon, nnjac, nnl, nnobj0, nnobj, (U_fp)fgcon,
			 (U_fp)fgobj, &x[1], ne, nlocj, &locj[1], &indj[1], &
			y[1], &fback, &gconu[1], &gobju[1], cu + 8, lencu, &
			iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
			leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
		    } else {
			*iexit = 63;
/* unable to move into undefined region */
		    }
		    goto L999;
		}
	    }
	    if (needg) {
		if (gobju[j] == gdummy) {
		    if (centrl) {
			gobj[j] = (fforwd * 4. - *fobj * 3. - fback) / dxj;
		    } else {
			gobj[j] = (fforwd - *fobj) / dxj;
		    }
		}
	    }
	    if (needj) {
		l = locg[j];
		k = locj[j];
		done = FALSE_;
/* +             --------------------------------------------------------- */
/* +             while (k .le. kmax  .and.  .not. done) do */
L140:
		if (k <= kmax && ! done) {
		    ir = indj[k];
		    if (ir > *nncon) {
			done = TRUE_;
		    } else {
			if (gconu[l] == gdummy) {
			    if (centrl) {
				gcon[l] = (fconu[ir] * 4. - fcon[ir] * 3. - y[
					ir]) / dxj;
			    } else {
				gcon[l] = (fconu[ir] - fcon[ir]) / dxj;
			    }
			}
			++l;
		    }
		    ++k;
		    goto L140;
/* +             end while */
/* +             --------------------------------------------------------- */
		}
	    }
/* j <= nnJac */
	    x[j] = xj;
	}
/* found */
    }
/*     ------------------------------------------------------------------ */
/*     The missing derivatives have been estimated. */
/*     Finish up with some housekeeping. */
/*     ------------------------------------------------------------------ */
/* L900: */
    l = lvldif + 1;
    iw[l + 189] += numf;
    iw[l + 194] += numf;
L999:
    return 0;
} /* s6fdg1_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6fdG1 */
/* Subroutine */ int s6srch_(integer *iexit, S_fp fgwrap, U_fp fgcon, U_fp 
	fgobj, logical *debug, logical *elastc, logical *fonly, logical *
	prfeas, integer *iobj, doublereal *sclobj, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, integer *nnjac, integer *nnl, 
	integer *nnobj0, integer *nnobj, integer *itn, doublereal *wolfeg, 
	doublereal *sgnobj, doublereal *step, doublereal *stepmn, doublereal *
	stepmx, doublereal *pnorm, doublereal *xnorm, doublereal *fmrt, 
	doublereal *fmrt1, doublereal *gmrt, doublereal *gmrt1, doublereal *
	sinf, doublereal *sinf1, doublereal *sinf2, doublereal *wtinf, 
	integer *ne, integer *nlocj, integer *locj, integer *indj, doublereal 
	*jcol, integer *negcon, integer *nlocg, integer *locg, doublereal *
	fobj1, doublereal *fcon1, doublereal *gcon1, doublereal *gobj1, 
	doublereal *fobj2, doublereal *fcon2, doublereal *gcon2, doublereal *
	gobj2, doublereal *dx, doublereal *dycon, doublereal *x, doublereal *
	x1, doublereal *x2, doublereal *ycon, doublereal *ycon1, doublereal *
	ycon2, doublereal *xpen, doublereal *y, doublereal *y1, doublereal *
	y2, char *cu, integer *lencu, integer *iu, integer *leniu, doublereal 
	*ru, integer *lenru, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(//\002 --------------------------------------"
	    "------\002/\002 Output from s6srch following iteration\002,i7,5x,"
	    "\002 Norm p =\002,1p,e11.2)";
    static char fmt_1700[] = "(\002 XXX  The line search has evaluated the f"
	    "unctions\002,i5,\002  times\002)";
    static char fmt_1600[] = "(\002 XXX  The search direction is uphill.  gM"
	    "rt  =\002,1p,e9.1)";
    static char fmt_1800[] = "(\002 stepmx =\002,1p,e11.2,\002    pNorm ="
	    "\002,e11.2,\002   gMrt =\002,e11.2,\002    numf  =\002,i3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsfi(icilist *), e_wsfi(void);

    /* Local variables */
    static doublereal g0, fa, fv, fw, gw, xv, xw;
    static char str[80];
    static integer nnj1;
    static doublereal eps0;
    static logical done;
    static integer jobj;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer maxf;
    static doublereal oldf, oldg;
    static integer nlin;
    static doublereal alow;
    static integer numf;
    static doublereal bupp;
    static logical vset, wset;
    static doublereal ftry;
    static integer nout;
    static doublereal gtry, xtry, fmrt2, gmrt2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), ddscl_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static doublereal epsaf, bigfx, dsinf, fbest, gbest;
    static logical moved;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sbest;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical first;
    extern /* Subroutine */ int s2aprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    s8gprd_(integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *), s6tols_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer modefg, nsamea, nsameb;
    static logical crampd, braktd;
    static doublereal factor;
    extern /* Subroutine */ int lsrchc_(integer *, logical *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     logical *, logical *, logical *, logical *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static logical nlnobj;
    static doublereal tolabs;
    static logical nlncon;
    static doublereal targtg;
    static integer inform__;
    extern /* Subroutine */ int lsrchq_(integer *, logical *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, logical *, 
	    logical *, logical *, logical *, logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static logical imprvd;
    static doublereal tolrel;
    static logical gknown, extrap;
    static integer iprint;
    static doublereal tolmax, stpmax;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static doublereal toltny;

    /* Fortran I/O blocks */
    static cilist io___65 = { 0, 0, 0, fmt_1000, 0 };
    static icilist io___89 = { 0, str, 0, fmt_1700, 80, 1 };
    static icilist io___90 = { 0, str, 0, fmt_1600, 80, 1 };
    static icilist io___91 = { 0, str, 0, fmt_1800, 80, 1 };


/*     ================================================================== */
/*     s6srch  finds a step along the search direction  p,  such that */
/*     the function  fMrt  is sufficiently reduced, i.e., */
/*               fMrt(x + step*p)  <  fMrt(x). */

/*     On entry, */
/*     step     Initial estimate of the step. */
/*     stepmx   Maximum allowable step. */
/*     wolfeG   Line search accuracy parameter in the range (0,1). */
/*              0.001 means very accurate.   0.99 means quite sloppy. */
/*     fonly    true if function-only search should be used. */
/*     gMrt1    The directional derivative of the merit function. */
/*     fMrt1    Current value of fMrt(x). */
/*     x(nb)    The base point for the search. */
/*     x2(nb)   x2 = x + dx,  the QP solution. */
/*     yCon     The base point for the multipliers. */
/*     p(nb)    The search direction. */
/*     yCon2    the QP multipliers. */
/*     dx(nb)   The search direction, x2 - x1. */

/*     Output parameters... */

/*     step     The final step. */
/*     fMrt1    The final value of fMrt. */
/*     x1       Final value of the variables. */
/*     fCon1,gCon1,gObj1,yCon1 are defined at the new point x1. */

/*     fCon2,gCon2,gObj2,yCon2 are work arrays. */

/*     iExit    Result */
/*     -----    ------ */
/*      >0      Fatal error */
/*       0      Repeat the search with smaller stpmax. */
/*      -1      The search is successful and step < stpmax. */
/*      -2      The search is successful and step = stpmax. */
/*      -3      A better point was found but too many functions */
/*              were needed (not sufficient decrease). */
/*      -4      stpmax < tolabs (too small to do a search). */
/*      -5      step   < stepmn (lsrchq only -- maybe want to switch */
/*              to central differences to get a better direction). */
/*      -6      No useful step. */
/*              The interval of uncertainty is less than 2*tolabs. */
/*              The minimizer is very close to step = zero */
/*              or the gradients are not sufficiently accurate. */
/*      -7      Too many function calls. */
/*      -8      Bad input parameters */
/*              (stpmax le toltny  or  oldg ge 0). */

/*     30 Dec 1991: First version based on NPSOL 4.6 routine npsrch. */
/*     28 Sep 1993: Allow functions to say "undefined at this point". */
/*                  (Back up and try again.) */
/*     18 Feb 1994: Back up in a similar way if vimax increases a lot. */
/*                  Deleted after first limited-memory version. */
/*     29 Dec 1994: Merit function calculations included explicitly. */
/*     06 Apr 1996: Special coding for the unit step. */
/*                  On entry, x2 is the QP solution. */
/*     18 Oct 1996: First Min-sum version. */
/*     17 Jul 1997: First thread-safe version. */
/*     11 Oct 1998: Facility to combine funobj and funcon added. */
/*     11 Jun 2000: Tolerances computed in a subroutine. */
/*     24 Oct 2000: Updated for SNOPT 6.1 */
/*     04 Aug 2003: snEXIT and snPRNT adopted. */
/*     15 Jun 2008: Call-status implemented correctly. */
/*     15 Nov 2010: Call-status removed from argument list. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y2;
    --y1;
    --y;
    --x2;
    --x1;
    --x;
    --dx;
    --xpen;
    --ycon2;
    --ycon1;
    --ycon;
    --dycon;
    --fcon2;
    --fcon1;
    --gobj2;
    --gobj1;
    --jcol;
    --indj;
    --locj;
    --gcon2;
    --gcon1;
    --locg;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    iprint = iw[12];
/* Print file */
    eps0 = rw[2];
/* eps**(4/5) */
    bigfx = rw[71];
/* unbounded objective. */
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Set the input parameters for lsrchc or lsrchq. */

/*     stepmn  is used by lsrchq.  If  step  would be less than  stepmn, */
/*             the search will be terminated early. */
/*             If p was found using forward or backward differences, */
/*             stepmn  should be positive (and related to the difference */
/*             interval used). */
/*             If p was found using central differences (lvlDif = 2) */
/*             stepmn  should be zero. */

/*     epsaf   is the absolute function precision. If f(x1) and f(x2) are */
/*             as close as  epsaf,  we cannot safely conclude from the */
/*             function values alone which of x1 or x2 is a better point. */

/*     tolabs  is an estimate of the absolute spacing between points */
/*             along  p.  This step should produce a perturbation */
/*             of  epsaf  in the merit function. */

/*     tolrel  is an estimate of the relative spacing between points */
/*             along  p. */

/*     toltny  is the minimum allowable absolute spacing between points */
/*             along  p. */
/*     ------------------------------------------------------------------ */
    stpmax = *stepmx;
    gknown = ! (*fonly);
    if (*fonly) {
	maxf = 15;
	modefg = 0;
    } else {
	maxf = 10;
	modefg = 2;
    }
    nlin = *n - *nnjac;
    nnj1 = *nnjac + 1;
    nout = iprint;
    jobj = *n + *iobj;
    nlncon = *nncon > 0;
    nlnobj = *nnobj > 0;
/*     Define the line search tolerances. */
    s6tols_(nb, &epsaf, stepmx, &tolabs, &tolrel, &toltny, pnorm, xnorm, fmrt,
	     &dx[1], &x[1], &rw[1], lenrw);
    oldf = *fmrt;
    oldg = *gmrt;
    *fmrt1 = *fmrt;
    *gmrt1 = *gmrt;
    dcopy_(nb, &x[1], &c__1, &x1[1], &c__1);
    if (nlncon) {
	dcopy_(nncon, &ycon[1], &c__1, &ycon1[1], &c__1);
    }
    if (*elastc) {
	*sinf1 = *sinf;
	dsinf = *sinf2 - *sinf1;
    }
    *fobj2 = 0.;
/* keeps ftnchek quiet */
    fmrt2 = 0.;
/* keeps ftnchek quiet */
    ftry = 0.;
/* keeps ftnchek quiet */
    gtry = 0.;
/* keeps ftnchek quiet */
    fv = 0.;
/* keeps ftnchek quiet */
    fw = 0.;
/* keeps ftnchek quiet */
    first = TRUE_;
    sbest = 0.;
    fbest = 0.;
    gbest = oldg * .99990000000000001;
    targtg = (1e-4 - *wolfeg) * oldg;
    g0 = gbest;
    if (*debug) {
	io___65.ciunit = nout;
	s_wsfe(&io___65);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*pnorm), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     ------------------------------------------------------------------ */
/*     Commence main loop, entering lsrchc or lsrchq two or more times. */
/*     first = true for the first entry,  false for subsequent entries. */
/*     done  = true indicates termination, in which case inform gives */
/*     the result of the search (with inform = iExit as above). */
/*     ------------------------------------------------------------------ */
/* +    repeat */
L200:
    if (gknown) {
	lsrchc_(&inform__, &first, debug, &done, &imprvd, &maxf, &numf, &nout,
		 &stpmax, &epsaf, &g0, &targtg, &ftry, &gtry, &tolabs, &
		tolrel, &toltny, step, &sbest, &fbest, &gbest, &braktd, &
		crampd, &extrap, &moved, &wset, &nsamea, &nsameb, &alow, &
		bupp, &factor, &xtry, &xw, &fw, &gw, &tolmax);
    } else {
	lsrchq_(&inform__, &first, debug, &done, &imprvd, &maxf, &numf, &nout,
		 &stpmax, stepmn, &epsaf, &g0, &targtg, &ftry, &tolabs, &
		tolrel, &toltny, step, &sbest, &fbest, &braktd, &crampd, &
		extrap, &moved, &vset, &wset, &nsamea, &nsameb, &alow, &bupp, 
		&fa, &factor, &xtry, &xw, &fw, &xv, &fv, &tolmax);
    }
    if (imprvd) {
	*fmrt1 = fmrt2;
	if (nlncon) {
	    dcopy_(nncon, &fcon2[1], &c__1, &fcon1[1], &c__1);
	    if (gknown) {
		dcopy_(negcon, &gcon2[1], &c__1, &gcon1[1], &c__1);
	    }
	}
	if (nlnobj) {
	    *fobj1 = *fobj2;
	    if (gknown) {
		dcopy_(nnobj, &gobj2[1], &c__1, &gobj1[1], &c__1);
	    }
/*              Terminate if the objective is unbounded below in the */
/*              feasible region. */
	    if (*sgnobj * *fobj1 < -bigfx && *prfeas) {
		*iexit = 21;
		goto L900;
	    }
	}
    }
/*        --------------------------------------------------------------- */
/*           done = false  first time through. */
/*        If done = false, the functions must be computed for the next */
/*                  entry to lsrchc or lsrchq. */
/*        If done = true,  this is the last time through and inform ge 1. */
/*        --------------------------------------------------------------- */
/* imprvd */
    if (! done) {
	if (*step != 1.) {
	    dcopy_(nb, &x1[1], &c__1, &x2[1], &c__1);
	    daxpy_(nb, step, &dx[1], &c__1, &x2[1], &c__1);
	    if (*elastc) {
		*sinf2 = *sinf1 + *step * dsinf;
	    }
	}
	if (*nnl > 0) {
	    (*fgwrap)(&inform__, &modefg, &nlncon, &nlnobj, n, negcon, nncon0,
		     nncon, nnjac, nnl, nnobj0, nnobj, (U_fp)fgcon, (U_fp)
		    fgobj, &x2[1], ne, nlocj, &locj[1], &indj[1], &fcon2[1], 
		    fobj2, &gcon2[1], &gobj2[1], cu + 8, lencu, &iu[1], leniu,
		     &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		    lenrw, (ftnlen)8, (ftnlen)8);
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
		} else {
		    inform__ = 0;
/* Redo the search */
		}
		goto L900;
	    }
	}
	if (*iobj == 0) {
	    fmrt2 = 0.;
	} else {
	    fmrt2 = *sgnobj * x2[jobj] * *sclobj;
	}
	if (nlnobj) {
	    fmrt2 += *sgnobj * *fobj2;
	}
	if (nlncon) {
/*              --------------------------------------------------------- */
/*              Compute w and y, the constraint violations and the */
/*              negative of the gradient of the merit function with */
/*              respect to the nonlinear slacks. These quantities define */
/*              the directional derivative of the merit function. */
/*              Finally, add the nonlinear terms to the merit function. */
/*              --------------------------------------------------------- */
	    if (*step <= 1.) {
		dcopy_(nncon, &ycon1[1], &c__1, &ycon2[1], &c__1);
		daxpy_(nncon, step, &dycon[1], &c__1, &ycon2[1], &c__1);
	    }
/*              Compute the constraint violations and aux. multipliers: */
/*              y1 = fCon  + A(linear) x - nonlinear slacks. */
/*              y  = yCon2 - xPen*viol. */
	    dcopy_(nncon, &fcon2[1], &c__1, &y1[1], &c__1);
	    daxpy_(nncon, &c_b36, &x2[*n + 1], &c__1, &y1[1], &c__1);
	    if (nlin > 0) {
		i__1 = nlin + 1;
		s2aprd_(&c__0, &eps0, ne, &i__1, &locj[nnj1], &indj[1], &jcol[
			1], &c_b40, &x2[nnj1], &nlin, &c_b40, &y1[1], nncon);
	    }
	    dcopy_(nncon, &y1[1], &c__1, &y[1], &c__1);
	    ddscl_(nncon, &xpen[1], &c__1, &y[1], &c__1);
	    fmrt2 = fmrt2 - ddot_(nncon, &ycon2[1], &c__1, &y1[1], &c__1) + 
		    ddot_(nncon, &y[1], &c__1, &y1[1], &c__1) * .5;
/*              If we are in elastic mode, include the contribution from */
/*              the violations of the elastic variables. */
	    if (*elastc) {
		fmrt2 += *wtinf * *sinf2;
	    }
	    if (gknown) {
		daxpy_(nncon, &c_b36, &ycon2[1], &c__1, &y[1], &c__1);
		dscal_(nncon, &c_b36, &y[1], &c__1);
	    }
	}
/* nlnCon */
	ftry = fmrt2 - oldf - oldg * 1e-4 * *step;
	if (gknown) {
/*              --------------------------------------------------------- */
/*              A gradient search is requested. */
/*              Compute the directional derivative gtry. */
/*              --------------------------------------------------------- */
	    if (*iobj == 0) {
		gmrt2 = 0.;
	    } else {
		gmrt2 = *sgnobj * dx[jobj] * *sclobj;
	    }
	    if (nlnobj) {
		gmrt2 += *sgnobj * ddot_(nnobj, &gobj2[1], &c__1, &dx[1], &
			c__1);
	    }
	    if (nlncon) {
/*                 ------------------------------------------------------ */
/*                 Form  J*dx (including linear columns). */
/*                 Set  y2 = J*dx + A(linear) dx - dx(slacks). */
/*                 ------------------------------------------------------ */
		s8gprd_(&c__0, &eps0, ne, nlocj, &locj[1], &indj[1], negcon, 
			nlocg, &locg[1], &gcon2[1], &c_b40, &dx[1], nnjac, &
			c_b59, &y2[1], nncon);
		if (nlin > 0) {
		    i__1 = nlin + 1;
		    s2aprd_(&c__0, &eps0, ne, &i__1, &locj[nnj1], &indj[1], &
			    jcol[1], &c_b40, &dx[nnj1], &nlin, &c_b40, &y2[1],
			     nncon);
		}
		daxpy_(nncon, &c_b36, &dx[*n + 1], &c__1, &y2[1], &c__1);
		gmrt2 = gmrt2 - ddot_(nncon, &y[1], &c__1, &y2[1], &c__1) - 
			ddot_(nncon, &y1[1], &c__1, &dycon[1], &c__1);
/*                 If we are elastic mode, include the contribution from */
/*                 the violations of the elastic variables. */
		if (*elastc) {
		    gmrt2 += *wtinf * dsinf;
		}
	    }
/* nlnCon */
	    gtry = gmrt2 - oldg * 1e-4;
	}
/* Gknown */
    }
/* +    until (      done) */
/* not done */
    if (! done) {
	goto L200;
    }
/*     ================================================================== */
/*     The search is done. */
/*     Finish with  x1 = the best point found so far. */
/*     ================================================================== */
    *step = sbest;
    if (imprvd) {
	dcopy_(nb, &x2[1], &c__1, &x1[1], &c__1);
	if (nlncon) {
	    dcopy_(nncon, &ycon2[1], &c__1, &ycon1[1], &c__1);
	}
	if (*elastc) {
	    *sinf1 = *sinf2;
	}
    } else if (*step > 0.) {
	daxpy_(nb, step, &dx[1], &c__1, &x1[1], &c__1);
	if (nlncon) {
	    daxpy_(nncon, step, &dycon[1], &c__1, &ycon1[1], &c__1);
	}
	if (*elastc) {
	    *sinf1 += *step * dsinf;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Print any warning messages. */
/*     ------------------------------------------------------------------ */
    if (inform__ == 7) {
	s_wsfi(&io___89);
	do_fio(&c__1, (char *)&numf, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
    } else if (inform__ == 8) {
	if (oldg >= 0.) {
	    s_wsfi(&io___90);
	    do_fio(&c__1, (char *)&oldg, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	} else {
	    s_wsfi(&io___91);
	    do_fio(&c__1, (char *)&(*stepmx), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*pnorm), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&oldg, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&numf, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    *iexit = -inform__;
L900:
    return 0;
} /* s6srch_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6srch */
/* Subroutine */ int s6tols_(integer *nb, doublereal *epsaf, doublereal *
	stepmx, doublereal *tolabs, doublereal *tolrel, doublereal *toltny, 
	doublereal *xdnorm, doublereal *xnorm, doublereal *f, doublereal *dx, 
	doublereal *x, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal q, s, t, eps, eps0, epsrf, tolax, tolrx;

/*     ================================================================== */
/*     s6tols defines various tolerances for the line search. */

/*     11 Jun 2000: First   version of s6tols. */
/*     11 Jun 2000: Current version of s6tols. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --dx;
    --rw;

    /* Function Body */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5) */
    epsrf = rw[73];
/*     ------------------------------------------------------------------ */
/* relative function precision. */
    *epsaf = max(epsrf,eps) * (abs(*f) + 1.);
    tolax = eps0;
    tolrx = eps0;
    t = *xnorm * tolrx + tolax;
    if (t < *xdnorm * *stepmx) {
	*tolabs = t / *xdnorm;
    } else {
	*tolabs = *stepmx;
    }
    *tolrel = max(tolrx,eps);
    t = 0.;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	s = (d__1 = dx[j], abs(d__1));
	q = (d__1 = x[j], abs(d__1)) * tolrx + tolax;
	if (s > q * t) {
	    t = s / q;
	}
    }
    if (t * *tolabs > 1.) {
	*toltny = 1. / t;
    } else {
	*toltny = *tolabs;
    }
    return 0;
} /* s6tols_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6tols */
/* Subroutine */ int s6usrf_(integer *mode, integer *n, doublereal *x, 
	doublereal *xpert, doublereal *damper, S_fp userfg, integer *needf, 
	integer *nf, doublereal *fpert, integer *needg, integer *leng, 
	doublereal *gpert, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, integer *iw, integer *leniw, 
	ftnlen cu_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, tries;
    extern /* Subroutine */ int s8stat_(integer *, integer *, integer *);
    static integer status;

/*     ================================================================== */
/*     s6usrf  attempts to compute the snoptA problem functions at a */
/*     point xPert.  If the functions are undefined at xPert, then the */
/*     evaluation takes place at a point on the ray joining xPert and */
/*     a point  x  at which the functions are known to be well-defined. */

/*     On entry: */
/*        xPert  is the point at which the functions are required. */

/*        x      is a point at which the problem functions have been */
/*               computed successfully. */

/*     On exit: */
/*        mode   is nonnegative if the problem function were evaluated */
/*               successfully.   Otherwise, five attempts to evaluate */
/*               the functions at points closer to x were unsuccessful. */

/*        If mode ge 0,  then the output values of damper and xPert are */
/*        defined as follows: */

/*        damper (0 lt damper le 1) is 1 if the functions were */
/*               evaluated successfuly at the input value of xPert. */
/*               Otherwise the problem functions were evaluated */
/*               successfuly at xPert + damper*(xPert - x). */

/*        xPert  is  xPert(in) + damper*(xPert(in) - x), */

/*     26 Oct 2002: First version. */
/*     22 Apr 2007: damper added as an output argument. */
/*     15 Jun 2008: Call-status implemented correctly. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xpert;
    --x;
    --fpert;
    --gpert;
    cu -= 8;
    --iu;
    --ru;
    --iw;

    /* Function Body */
    *mode = 0;
/* Determine the status of this call. */
    s8stat_(&status, &iw[1], leniw);
    *damper = 1.;
    tries = 0;
/*     ================================================================== */
/*     Repeat                       (until problem functions are defined) */
L100:
    ++tries;
    (*userfg)(&status, n, &xpert[1], needf, nf, &fpert[1], needg, leng, &
	    gpert[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, (ftnlen)8);
    *mode = status;
    status = 0;
    if (*mode == -1) {
	*damper = *damper * *damper / 10.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    xpert[j] = *damper * xpert[j] + (1. - *damper) * x[j];
	}
    }
/* +    until (.not. (mode .eq. -1  .and.  tries .lt. 5)) */
    if (*mode == -1 && tries < 5) {
	goto L100;
    }
/*     ================================================================== */
    return 0;
} /* s6usrf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s6usrf */
/* Subroutine */ int lsrchc_(integer *iexit, logical *first, logical *debug, 
	logical *done, logical *imprvd, integer *maxf, integer *numf, integer 
	*nout, doublereal *alfmax, doublereal *epsaf, doublereal *g0, 
	doublereal *targtg, doublereal *ftry, doublereal *gtry, doublereal *
	tolabs, doublereal *tolrel, doublereal *toltny, doublereal *alfa, 
	doublereal *alfbst, doublereal *fbest, doublereal *gbest, logical *
	braktd, logical *crampd, logical *extrap, logical *moved, logical *
	wset, integer *nsamea, integer *nsameb, doublereal *a, doublereal *b, 
	doublereal *factor, doublereal *xtry, doublereal *xw, doublereal *fw, 
	doublereal *gw, doublereal *tolmax)
{
    /* Format strings */
    static char fmt_1000[] = "(/\002     g0  tolabs  alfmax        \002,1p,2"
	    "e22.14,e16.8/\002 targtg  tolrel   epsaf        \002,1p,2e22.14,"
	    "e16.8/\002 crampd                        \002,l3)";
    static char fmt_1100[] = "(/\002 alfa    ftry    gtry          \002,1p,2"
	    "e22.14,e16.8)";
    static char fmt_1200[] = "(/\002 a       b       b - a   tol   \002,1p,2"
	    "e22.14,2e16.8/\002 nsamea  nsameb  numf          \002,3i3/\002 b"
	    "raktd  extrap  closef  imprvd\002,4l3/\002 found   quitI        "
	    "         \002,2l3/\002 alfbst  fbest   gbest         \002,1p,3e2"
	    "2.14/\002 alfaw   fw      gw            \002,1p,3e22.14)";
    static char fmt_2200[] = "(\002 Parabola.\002)";
    static char fmt_2100[] = "(\002 Cubic.   \002)";
    static char fmt_2400[] = "(\002 Geo. bisection. xtry,daux,dtry\002,1p,3e"
	    "22.14)";
    static char fmt_2300[] = "(\002 Bisection.              xmidpt\002,1p,e2"
	    "2.14)";
    static char fmt_2500[] = "(\002 Polynomial fit accepted.  xtry\002,1p,e2"
	    "2.14)";
    static char fmt_3000[] = "(\002 ----------------------------------------"
	    "------------\002/)";

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal q, r__, s, tol, absr, daux, dtry, scale;
    static logical found, fitok;
    static doublereal truea, trueb;
    static logical quitf, quiti, setxw, badfun;
    static doublereal artifa, artifb;
    static logical closef;
    static doublereal xmidpt;

    /* Fortran I/O blocks */
    static cilist io___108 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___109 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___119 = { 0, 0, 0, fmt_2200, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_2100, 0 };
    static cilist io___129 = { 0, 0, 0, fmt_2400, 0 };
    static cilist io___130 = { 0, 0, 0, fmt_2300, 0 };
    static cilist io___131 = { 0, 0, 0, fmt_2500, 0 };
    static cilist io___132 = { 0, 0, 0, fmt_3000, 0 };


/*     ================================================================== */
/*     lsrchc  finds a sequence of improving estimates of a minimizer of */
/*     the univariate function f(alpha) in the interval (0,alfmax]. */
/*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0. */
/*     lsrchc requires both  f(alpha)  and  f'(alpha)  to be evaluated at */
/*     points in the interval.  Estimates of the minimizer are computed */
/*     using safeguarded cubic interpolation. */

/*     Reverse communication is used to allow the calling program to */
/*     evaluate f and f'.  Some of the parameters must be set or tested */
/*     by the calling program.  The remainder would ordinarily be local */
/*     variables. */

/*     Input parameters (relevant to the calling program) */
/*     -------------------------------------------------- */

/*     first         must be true on the first entry. It is subsequently */
/*                   altered by lsrchc. */

/*     debug         specifies whether detailed output is wanted. */

/*     maxf          is an upper limit on the number of times lsrchc is */
/*                   to be entered consecutively with done = false */
/*                   (following an initial entry with first = true). */

/*     alfa          is the first estimate of a minimizer.  alfa is */
/*                   subsequently altered by lsrchc (see below). */

/*     alfmax        is the upper limit of the interval to be searched. */

/*     epsaf         is an estimate of the absolute precision in the */
/*                   computed value of f(0). */

/*     ftry, gtry    are the values of f, f'  at the new point */
/*                   alfa = alfbst + xtry. */

/*     g0            is the value of f'(0).  g0 must be negative. */

/*     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs */
/*                   such that if f has already been evaluated at alfa, */
/*                   it will not be evaluated closer than tol(alfa). */
/*                   These values may be reduced by lsrchc. */

/*     targtg        is the target value of abs(f'(alfa)). The search */
/*                   is terminated when */
/*                    abs(f'(alfa)) le targtg and f(alfa) lt 0. */

/*     toltny        is the smallest value that tolabs is allowed to be */
/*                   reduced to. */

/*     Output parameters (relevant to the calling program) */
/*     --------------------------------------------------- */

/*     imprvd        is true if the previous alfa was the best point so */
/*                   far.  Any related quantities should be saved by the */
/*                   calling program (e.g., gradient arrays) before */
/*                   paying attention to the variable done. */

/*     done = false  means the calling program should evaluate */
/*                      ftry = f(alfa),  gtry = f'(alfa) */
/*                   for the new trial alfa, and re-enter lsrchc. */

/*     done = true   means that no new alfa was calculated.  The value */
/*                   of iExit gives the result of the search as follows */

/*                   iExit = 1 means the search has terminated */
/*                             successfully with alfbst < alfmax. */

/*                   iExit = 2 means the search has terminated */
/*                             successfully with alfbst = alfmax. */

/*                   iExit = 3 means that the search failed to find a */
/*                             point of sufficient decrease. */
/*                             The function is either decreasing at */
/*                             alfmax or maxf function evaluations */
/*                             have been exceeded. */

/*                   iExit = 4 means alfmax is so small that a search */
/*                             should not have been attempted. */

/*                   iExit = 5 is never set by lsrchc. */

/*                   iExit = 6 means the search has failed to find a */
/*                             useful step.  The interval of uncertainty */
/*                             is [0,b] with b < 2*tolabs. A minimizer */
/*                             lies very close to alfa = 0, or f'(0) is */
/*                             not sufficiently accurate. */

/*                   iExit = 7 if no better point could be found after */
/*                             maxf  function calls. */

/*                   iExit = 8 means the input parameters were bad. */
/*                             alfmax le toltny  or g0 ge zero. */
/*                             No function evaluations were made. */

/*     numf          counts the number of times lsrchc has been entered */
/*                   consecutively with done = false (i.e., with a new */
/*                   function value ftry). */

/*     alfa          is the point at which the next function ftry and */
/*                   derivative gtry must be computed. */

/*     alfbst        should be accepted by the calling program as the */
/*                   approximate minimizer, whenever lsrchc returns */
/*                   iExit = 1 or 2 (and possibly 3). */

/*     fbest, gbest  will be the corresponding values of f, f'. */


/*     The following parameters retain information between entries */
/*     ----------------------------------------------------------- */

/*     braktd        is false if f and f' have not been evaluated at */
/*                   the far end of the interval of uncertainty.  In this */
/*                   case, the point b will be at alfmax + tol(alfmax). */

/*     crampd        is true if alfmax is very small (le tolabs).  If the */
/*                   search fails, this indicates that a zero step should */
/*                   be taken. */

/*     extrap        is true if xw lies outside the interval of */
/*                   uncertainty.  In this case, extra safeguards are */
/*                   applied to allow for instability in the polynomial */
/*                   fit. */

/*     moved         is true if a better point has been found, i.e., */
/*                   alfbst gt 0. */

/*     wset          records whether a second-best point has been */
/*                   determined it will always be true when convergence */
/*                   is tested. */

/*     nsamea        is the number of consecutive times that the */
/*                   left-hand end point of the interval of uncertainty */
/*                   has remained the same. */

/*     nsameb        similarly for the right-hand end. */

/*     a, b, alfbst  define the current interval of uncertainty. */
/*                   A minimizer lies somewhere in the interval */
/*                   [alfbst + a, alfbst + b]. */

/*     alfbst        is the best point so far.  It is always at one end */
/*                   of the interval of uncertainty.  hence we have */
/*                   either  a lt 0,  b = 0  or  a = 0,  b gt 0. */

/*     fbest, gbest  are the values of f, f' at the point alfbst. */

/*     factor        controls the rate at which extrapolated estimates */
/*                   of alfa may expand into the interval of uncertainty. */
/*                   factor is not used if a minimizer has been bracketed */
/*                   (i.e., when the variable braktd is true). */

/*     fw, gw        are the values of f, f' at the point alfbst + xw. */
/*                   they are not defined until wset is true. */

/*     xtry          is the trial point in the shifted interval (a, b). */

/*     xw            is such that  alfbst + xw  is the second-best point. */
/*                   it is not defined until  wset  is true. */
/*                   in some cases,  xw  will replace a previous  xw */
/*                   that has a lower function but has just been excluded */
/*                   from the interval of uncertainty. */


/*     Systems Optimization Laboratory, Stanford University, California. */
/*     Original version February 1982.  Rev. May 1983. */
/*     Original f77 version 22-August-1985. */
/*     14 Sep 1992: Introduced quitI, quitF, etc. */
/*     22 Nov 1995: Altered criterion for reducing the step below tolabs. */
/*     17 Jul 1997: Removed saved variables for thread-safe version. */
/*     19 Apr 2000: QuitF only allowed after a move. */
/*     19 Apr 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Local variables */
/*     =============== */

/*     closef     is true if the new function ftry is within epsaf of */
/*                fbest (up or down). */

/*     found      is true if the sufficient decrease conditions hold at */
/*                alfbst. */

/*     quitF      is true when  maxf  function calls have been made. */

/*     quitI      is true when the interval of uncertainty is less than */
/*                2*tol. */
/*  --------------------------------------------------------------------- */
    badfun = FALSE_;
    quitf = FALSE_;
    quiti = FALSE_;
    *imprvd = FALSE_;
    if (*first) {
/*        --------------------------------------------------------------- */
/*        First entry.  Initialize various quantities, check input data */
/*        and prepare to evaluate the function at the initial alfa. */
/*        --------------------------------------------------------------- */
	*first = FALSE_;
	*numf = 0;
	*alfbst = 0.;
	badfun = *alfmax <= *toltny || *g0 >= 0.;
	*done = badfun;
	*moved = FALSE_;
	if (! (*done)) {
	    *braktd = FALSE_;
	    *crampd = *alfmax <= *tolabs;
	    *extrap = FALSE_;
	    *wset = FALSE_;
	    *nsamea = 0;
	    *nsameb = 0;
	    *tolmax = *tolabs + *tolrel * *alfmax;
	    *a = 0.;
	    *b = *alfmax + *tolmax;
	    *factor = 5.;
	    tol = *tolabs;
	    *xtry = *alfa;
	    if (*debug) {
		io___108.ciunit = *nout;
		s_wsfe(&io___108);
		do_fio(&c__1, (char *)&(*g0), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*tolabs), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*alfmax), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*targtg), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*tolrel), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*epsaf), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*crampd), (ftnlen)sizeof(logical));
		e_wsfe();
	    }
	}
    } else {
/*        --------------------------------------------------------------- */
/*        Subsequent entries. The function has just been evaluated at */
/*        alfa = alfbst + xtry,  giving ftry and gtry. */
/*        --------------------------------------------------------------- */
	if (*debug) {
	    io___109.ciunit = *nout;
	    s_wsfe(&io___109);
	    do_fio(&c__1, (char *)&(*alfa), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ftry), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*gtry), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	++(*numf);
	++(*nsamea);
	++(*nsameb);
	if (! (*braktd)) {
	    *tolmax = *tolabs + *tolrel * *alfmax;
	    *b = *alfmax - *alfbst + *tolmax;
	}
/*        See if the new step is better.  If alfa is large enough that */
/*        ftry can be distinguished numerically from zero,  the function */
/*        is required to be sufficiently negative. */
	closef = (d__1 = *ftry - *fbest, abs(d__1)) <= *epsaf;
	if (closef) {
	    *imprvd = abs(*gtry) <= abs(*gbest);
	} else {
	    *imprvd = *ftry < *fbest;
	}
	if (*imprvd) {
/*           We seem to have an improvement.  The new point becomes the */
/*           origin and other points are shifted accordingly. */
	    *fw = *fbest;
	    *fbest = *ftry;
	    *gw = *gbest;
	    *gbest = *gtry;
	    *alfbst = *alfa;
	    *moved = TRUE_;
	    *a -= *xtry;
	    *b -= *xtry;
	    *xw = 0. - *xtry;
	    *wset = TRUE_;
	    *extrap = *xw < 0. && *gbest < 0. || *xw > 0. && *gbest > 0.;
/*           Decrease the length of the interval of uncertainty. */
	    if (*gtry <= 0.) {
		*a = 0.;
		*nsamea = 0;
	    } else {
		*b = 0.;
		*nsameb = 0;
		*braktd = TRUE_;
	    }
	} else {
/*           The new function value is not better than the best point so */
/*           far.  The origin remains unchanged but the new point may */
/*           qualify as xw.  xtry must be a new bound on the best point. */
	    if (*xtry <= 0.) {
		*a = *xtry;
		*nsamea = 0;
	    } else {
		*b = *xtry;
		*nsameb = 0;
		*braktd = TRUE_;
	    }
/*           If xw has not been set or ftry is better than fw, update the */
/*           points accordingly. */
	    if (*wset) {
		setxw = *ftry < *fw || ! (*extrap);
	    } else {
		setxw = TRUE_;
	    }
	    if (setxw) {
		*xw = *xtry;
		*fw = *ftry;
		*gw = *gtry;
		*wset = TRUE_;
		*extrap = FALSE_;
	    }
	}
/*        --------------------------------------------------------------- */
/*        Check the termination criteria.  wset will always be true. */
/*        --------------------------------------------------------------- */
	tol = *tolabs + *tolrel * *alfbst;
	truea = *alfbst + *a;
	trueb = *alfbst + *b;
	found = abs(*gbest) <= *targtg;
	quitf = *numf >= *maxf && *moved;
	quiti = *b - *a <= tol + tol;
	if (quiti && ! (*moved)) {
/*           The interval of uncertainty appears to be small enough, */
/*           but no better point has been found.  Check that changing */
/*           alfa by b-a changes f by less than epsaf. */
	    tol /= 10.;
	    *tolabs = tol;
	    quiti = tol <= *toltny || abs(*fw) <= *epsaf && *gw <= *epsaf;
	}
	*done = quitf || quiti || found;
	if (*debug) {
	    io___115.ciunit = *nout;
	    s_wsfe(&io___115);
	    do_fio(&c__1, (char *)&truea, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&trueb, (ftnlen)sizeof(doublereal));
	    d__1 = *b - *a;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*nsamea), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nsameb), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*numf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*braktd), (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&(*extrap), (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&closef, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&(*imprvd), (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&found, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&quiti, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&(*alfbst), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*fbest), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*gbest), (ftnlen)sizeof(doublereal));
	    d__2 = *alfbst + *xw;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*fw), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*gw), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/*        --------------------------------------------------------------- */
/*        Proceed with the computation of a trial steplength. */
/*        The choices are... */
/*        1. Parabolic fit using derivatives only, if the f values are */
/*           close. */
/*        2. Cubic fit for a minimizer, using both f and f'. */
/*        3. Damped cubic or parabolic fit if the regular fit appears to */
/*           be consistently overestimating the distance to a minimizer. */
/*        4. Bisection, geometric bisection, or a step of  tol  if */
/*           choices 2 or 3 are unsatisfactory. */
/*        --------------------------------------------------------------- */
	if (! (*done)) {
	    xmidpt = (*a + *b) * .5;
	    s = 0.;
	    q = 0.;
	    if (closef) {
/*              --------------------------------------------------------- */
/*              Fit a parabola to the two best gradient values. */
/*              --------------------------------------------------------- */
		s = *gbest;
		q = *gbest - *gw;
		if (*debug) {
		    io___119.ciunit = *nout;
		    s_wsfe(&io___119);
		    e_wsfe();
		}
	    } else {
/*              --------------------------------------------------------- */
/*              Fit cubic through  fbest  and  fw. */
/*              --------------------------------------------------------- */
		if (*debug) {
		    io___120.ciunit = *nout;
		    s_wsfe(&io___120);
		    e_wsfe();
		}
		fitok = TRUE_;
		r__ = (*fbest - *fw) * 3. / *xw + *gbest + *gw;
		absr = abs(r__);
		s = sqrt((abs(*gbest))) * sqrt((abs(*gw)));
/*              Compute  q =  the square root of  r*r - gbest*gw. */
/*              The method avoids unnecessary underflow and overflow. */
		if (*gw < 0. && *gbest > 0. || *gw > 0. && *gbest < 0.) {
		    scale = absr + s;
		    if (scale == 0.) {
			q = 0.;
		    } else {
/* Computing 2nd power */
			d__1 = absr / scale;
/* Computing 2nd power */
			d__2 = s / scale;
			q = scale * sqrt(d__1 * d__1 + d__2 * d__2);
		    }
		} else if (absr >= s) {
		    q = sqrt(absr + s) * sqrt(absr - s);
		} else {
		    fitok = FALSE_;
		}
		if (fitok) {
/*                 Compute a minimizer of the fitted cubic. */
		    if (*xw < 0.) {
			q = -q;
		    }
		    s = *gbest - r__ - q;
		    q = *gbest - *gw - q - q;
		}
	    }
/*           ------------------------------------------------------------ */
/*           Construct an artificial interval  (artifa, artifb)  in which */
/*           the new estimate of a minimizer must lie.  Set a default */
/*           value of xtry that will be used if the polynomial fit fails. */
/*           ------------------------------------------------------------ */
	    artifa = *a;
	    artifb = *b;
	    if (! (*braktd)) {
/*              A minimizer has not been bracketed.  Set an artificial */
/*              upper bound by expanding the interval  xw  by a suitable */
/*              factor. */
		*xtry = -(*factor) * *xw;
		artifb = *xtry;
		if (*alfbst + *xtry < *alfmax) {
		    *factor *= 5.;
		}
	    } else if (*extrap) {
/*              The points are configured for an extrapolation. */
/*              Set a default value of  xtry  in the interval  (a, b) */
/*              that will be used if the polynomial fit is rejected.  In */
/*              the following,  dtry  and  daux  denote the lengths of */
/*              the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if */
/*              appropriate).  The value of  xtry is the point at which */
/*              the exponents of  dtry  and  daux  are approximately */
/*              bisected. */
		daux = abs(*xw);
		dtry = *b - *a;
		if (daux >= dtry) {
		    *xtry = dtry * 5. * (dtry / daux + .1) / 11.;
		} else {
		    *xtry = sqrt(daux) * .5 * sqrt(dtry);
		}
		if (*xw > 0.) {
		    *xtry = -(*xtry);
		}
		if (*debug) {
		    io___129.ciunit = *nout;
		    s_wsfe(&io___129);
		    do_fio(&c__1, (char *)&(*xtry), (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&daux, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&dtry, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
/*              Reset the artificial bounds.  If the point computed by */
/*              extrapolation is rejected,  xtry will remain at the */
/*              relevant artificial bound. */
		if (*xtry <= 0.) {
		    artifa = *xtry;
		}
		if (*xtry > 0.) {
		    artifb = *xtry;
		}
	    } else {
/*              The points are configured for an interpolation.  The */
/*              default value xtry bisects the interval of uncertainty. */
/*              the artificial interval is just (a, b). */
		*xtry = xmidpt;
		if (*debug) {
		    io___130.ciunit = *nout;
		    s_wsfe(&io___130);
		    do_fio(&c__1, (char *)&(*xtry), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfe();
		}
		if (*nsamea >= 3 || *nsameb >= 3) {
/*                 If the interpolation appears to be overestimating the */
/*                 distance to a minimizer,  damp the interpolation. */
		    *factor /= 5.;
		    s = *factor * s;
		} else {
		    *factor = 1.;
		}
	    }
/*           ------------------------------------------------------------ */
/*           The polynomial fits give  (s/q)*xw  as the new step. */
/*           Reject this step if it lies outside  (artifa, artifb). */
/*           ------------------------------------------------------------ */
	    if (q != 0.) {
		if (q < 0.) {
		    s = -s;
		}
		if (q < 0.) {
		    q = -q;
		}
		if (s * *xw >= q * artifa && s * *xw <= q * artifb) {
/*                 Accept the polynomial fit. */
		    if ((d__1 = s * *xw, abs(d__1)) >= q * tol) {
			*xtry = s / q * *xw;
		    } else {
			*xtry = 0.;
		    }
		    if (*debug) {
			io___131.ciunit = *nout;
			s_wsfe(&io___131);
			do_fio(&c__1, (char *)&(*xtry), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }
/*     ================================================================== */
    if (! (*done)) {
	*alfa = *alfbst + *xtry;
	if (*braktd || *alfa < *alfmax - *tolmax) {
/*           The function must not be evaluated too close to a or b. */
/*           (It has already been evaluated at both those points.) */
	    if (*xtry <= *a + tol || *xtry >= *b - tol) {
		if ((*a + *b) * .5 <= 0.) {
		    *xtry = -tol;
		} else {
		    *xtry = tol;
		}
		*alfa = *alfbst + *xtry;
	    }
	} else {
/*           The step is close to, or larger than alfmax, replace it by */
/*           alfmax to force evaluation of  f  at the boundary. */
	    *braktd = TRUE_;
	    *xtry = *alfmax - *alfbst;
	    *alfa = *alfmax;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Exit. */
/*     ------------------------------------------------------------------ */
    if (*done) {
	if (badfun) {
	    *iexit = 8;
/* bad arguments */
	} else if (found) {
	    if (*alfbst < *alfmax) {
		*iexit = 1;
/* Sufficient decrease */
	    } else {
		*iexit = 2;
/* Suff. Decrease on the boundary */
	    }
	} else if (*moved) {
	    *iexit = 3;
/* Decr at boundary or max funs */
	} else if (quitf) {
	    *iexit = 7;
/* No new point after max funs */
	} else if (*crampd) {
	    *iexit = 4;
/* alfmax too mall */
	} else {
	    *iexit = 6;
/* [a,b] too small */
	}
    }
    if (*debug) {
	io___132.ciunit = *nout;
	s_wsfe(&io___132);
	e_wsfe();
    }
    return 0;
} /* lsrchc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lsrchc */
/* Subroutine */ int lsrchq_(integer *iexit, logical *first, logical *debug, 
	logical *done, logical *imprvd, integer *maxf, integer *numf, integer 
	*nout, doublereal *alfmax, doublereal *alfsml, doublereal *epsaf, 
	doublereal *g0, doublereal *targtg, doublereal *ftry, doublereal *
	tolabs, doublereal *tolrel, doublereal *toltny, doublereal *alfa, 
	doublereal *alfbst, doublereal *fbest, logical *braktd, logical *
	crampd, logical *extrap, logical *moved, logical *vset, logical *wset,
	 integer *nsamea, integer *nsameb, doublereal *a, doublereal *b, 
	doublereal *fa, doublereal *factor, doublereal *xtry, doublereal *xw, 
	doublereal *fw, doublereal *xv, doublereal *fv, doublereal *tolmax)
{
    /* Format strings */
    static char fmt_1000[] = "(/\002     g0  tolabs  alfmax        \002,1p,2"
	    "e22.14,e16.8/\002 targtg  tolrel   epsaf        \002,1p,2e22.14,"
	    "e16.8/\002 crampd                        \002,l3)";
    static char fmt_1100[] = "(/\002 alfa    ftry                  \002,1p,2"
	    "e22.14)";
    static char fmt_1200[] = "(/\002 a       b       b - a   tol   \002,1p,2"
	    "e22.14,2e16.8/\002 nsamea  nsameb  numf          \002,3i3/\002 b"
	    "raktd  extrap  closef  imprvd\002,4l3/\002 found   quitI   quitF"
	    "Z  quitS \002,4l3/\002 alfbst  fbest                 \002,1p,2e2"
	    "2.14/\002 alfaw   fw                    \002,1p,2e22.14)";
    static char fmt_1300[] = "(\002 alfav   fv                    \002,1p,2e"
	    "22.14/)";
    static char fmt_2200[] = "(\002 Parabolic fit,  three points. \002)";
    static char fmt_2100[] = "(\002 Parabolic fit,    two points. \002)";
    static char fmt_2500[] = "(\002 Geo. bisection. xtry,daux,dtry\002,1p,3e"
	    "22.14)";
    static char fmt_2400[] = "(\002 Exponent reduced.  Trial point\002,1p,e2"
	    "2.14)";
    static char fmt_2600[] = "(\002 Polynomial fit accepted.  xtry\002,1p,e2"
	    "2.14)";
    static char fmt_3000[] = "(\002 ----------------------------------------"
	    "------------\002/)";

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal q, s, gw, gv, tol, daux, dtry;
    static logical found;
    static doublereal truea, trueb;
    static logical quitf, quiti, quits, setxv, xinxw, badfun;
    static doublereal artifa, artifb;
    static logical closef;
    static doublereal endpnt, xmidpt;
    static logical quitfz;

    /* Fortran I/O blocks */
    static cilist io___139 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___140 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___147 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___148 = { 0, 0, 0, fmt_1300, 0 };
    static cilist io___154 = { 0, 0, 0, fmt_2200, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_2100, 0 };
    static cilist io___161 = { 0, 0, 0, fmt_2500, 0 };
    static cilist io___162 = { 0, 0, 0, fmt_2400, 0 };
    static cilist io___163 = { 0, 0, 0, fmt_2600, 0 };
    static cilist io___164 = { 0, 0, 0, fmt_3000, 0 };


/*     ================================================================== */
/*     lsrchq  finds a sequence of improving estimates of a minimizer of */
/*     the univariate function f(alpha) in the interval (0,alfmax]. */
/*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0. */
/*     lsrchq  requires  f(alpha) (but not f'(alpha)) to be evaluated */
/*     in the interval.  New estimates of a minimizer are computed using */
/*     safeguarded quadratic interpolation. */

/*     Reverse communication is used to allow the calling program to */
/*     evaluate f.  Some of the parameters must be set or tested by the */
/*     calling program.  The remainder would ordinarily be local */
/*     variables. */

/*     Input parameters (relevant to the calling program) */
/*     -------------------------------------------------- */

/*     first         must be true on the first entry.  It is subsequently */
/*                   altered by lsrchq. */

/*     debug         specifies whether detailed output is wanted. */

/*     maxf          is an upper limit on the number of times lsrchq is */
/*                   to be entered consecutively with done = false */
/*                   (following an initial entry with first = true). */

/*     alfa          is the first estimate of a minimizer.  alfa is */
/*                   subsequently altered by lsrchq (see below). */

/*     alfmax        is the upper limit of the interval to be searched. */

/*     alfsml        is intended to prevent inefficiency when a minimizer */
/*                   is very small, for cases where the calling program */
/*                   would prefer to redefine f'(alfa).  alfsml is */
/*                   allowed to be zero.  Early termination will occur if */
/*                   lsrchq determines that a minimizer lies somewhere in */
/*                   the interval [0, alfsml) (but not if alfmax is */
/*                   smaller that alfsml). */

/*     epsaf         is an estimate of the absolute precision in the */
/*                   computed value of f(0). */

/*     ftry          the value of f at the new point */
/*                   alfa = alfbst + xtry. */

/*     g0            is the value of f'(0).  g0 must be negative. */

/*     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs */
/*                   such that if f has already been evaluated at alfa, */
/*                   it will not be evaluated closer than tol(alfa). */
/*                   These values may be reduced by lsrchq. */

/*     targtg        is the target value of abs(f'(alfa)). The search */
/*                   is terminated when */
/*                    abs(f'(alfa)) le targtg and f(alfa) lt 0. */

/*     toltny        is the smallest value that tolabs is allowed to be */
/*                   reduced to. */

/*     Output parameters (relevant to the calling program) */
/*     --------------------------------------------------- */

/*     imprvd        is true if the previous alfa was the best point so */
/*                   far.  Any related quantities should be saved by the */
/*                   calling program (e.g., arrays) before paying */
/*                   attention to the variable done. */

/*     done = false  means the calling program should evaluate ftry */
/*                   for the new trial step alfa, and reenter lsrchq. */

/*     done = true   means that no new alfa was calculated.  The value */
/*                   of iExit gives the result of the search as follows */

/*                   iExit = 1 means the search has terminated */
/*                             successfully with alfbst < alfmax. */

/*                   iExit = 2 means the search has terminated */
/*                             successfully with alfbst = alfmax. */

/*                   iExit = 3 means that the search failed to find a */
/*                             point of sufficient decrease in maxf */
/*                             functions, but a lower point was found. */

/*                   iExit = 4 means alfmax is so small that a search */
/*                             should not have been attempted. */

/*                   iExit = 5 means that the search was terminated */
/*                             because of alfsml (see above). */

/*                   iExit = 6 means the search has failed to find a */
/*                             useful step.  The interval of uncertainty */
/*                             is [0,b] with b < 2*tolabs. A minimizer */
/*                             lies very close to alfa = 0, or f'(0) is */
/*                             not sufficiently accurate. */

/*                   iExit = 7 if no better point could be found after */
/*                             maxf  function calls. */

/*                   iExit = 8 means the input parameters were bad. */
/*                             alfmax le toltny  or  g0 ge zero. */
/*                             No function evaluations were made. */

/*     numf          counts the number of times lsrchq has been entered */
/*                   consecutively with done = false (i.e., with a new */
/*                   function value ftry). */

/*     alfa          is the point at which the next function ftry must */
/*                   be computed. */

/*     alfbst        should be accepted by the calling program as the */
/*                   approximate minimizer, whenever lsrchq returns */
/*                   iExit = 1, 2 or 3. */

/*     fbest         will be the corresponding value of f. */

/*     The following parameters retain information between entries */
/*     ----------------------------------------------------------- */

/*     braktd        is false if f has not been evaluated at the far end */
/*                   of the interval of uncertainty.  In this case, the */
/*                   point b will be at alfmax + tol(alfmax). */

/*     crampd        is true if alfmax is very small (le tolabs).  If the */
/*                   search fails, this indicates that a zero step should */
/*                   be taken. */

/*     extrap        is true if alfbst has moved at least once and xv */
/*                   lies outside the interval of uncertainty.  In this */
/*                   case, extra safeguards are applied to allow for */
/*                   instability in the polynomial fit. */

/*     moved         is true if a better point has been found, i.e., */
/*                   alfbst gt 0. */

/*     vset          records whether a third-best point has been defined. */

/*     wset          records whether a second-best point has been */
/*                   defined.  It will always be true by the time the */
/*                   convergence test is applied. */

/*     nsamea        is the number of consecutive times that the */
/*                   left-hand end point of the interval of uncertainty */
/*                   has remained the same. */

/*     nsameb        similarly for the right-hand end. */

/*     a, b, alfbst  define the current interval of uncertainty. */
/*                   A minimizer lies somewhere in the  interval */
/*                   [alfbst + a, alfbst + b]. */

/*     alfbst        is the best point so far.  It lies strictly within */
/*                   [atrue,btrue]  (except when alfbst has not been */
/*                   moved, in which case it lies at the left-hand end */
/*                   point).  Hence we have a .le. 0 and b .gt. 0. */

/*     fbest         is the value of f at the point alfbst. */

/*     fa            is the value of f at the point alfbst + a. */

/*     factor        controls the rate at which extrapolated estimates of */
/*                   alfa  may expand into the interval of uncertainty. */
/*                   Factor is not used if a minimizer has been bracketed */
/*                   (i.e., when the variable braktd is true). */

/*     fv, fw        are the values of f at the points alfbst + xv  and */
/*                   alfbst + xw.  They are not defined until  vset  or */
/*                   wset  are true. */

/*     xtry          is the trial point within the shifted interval */
/*                   (a, b).  The new trial function value must be */
/*                   computed at the point alfa = alfbst + xtry. */

/*     xv            is such that alfbst + xv is the third-best point. */
/*                   It is not defined until vset is true. */

/*     xw            is such that alfbst + xw is the second-best point. */
/*                   It is not defined until wset is true.  In some */
/*                   cases,  xw will replace a previous xw that has a */
/*                   lower function but has just been excluded from */
/*                   (a,b). */

/*     Systems Optimization Laboratory, Stanford University, California. */
/*     Original version February 1982.  Rev. May 1983. */
/*     Original F77 version 22-August-1985. */
/*     17 Jul 1997: Removed saved variables for thread-safe version. */
/*     31 Jul 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Local variables */
/*     =============== */

/*     closef     is true if the worst function fv is within epsaf of */
/*                fbest (up or down). */

/*     found      is true if the sufficient decrease conditions holds at */
/*                alfbst. */

/*     quitF      is true when  maxf  function calls have been made. */

/*     quitFZ     is true when the three best function values are within */
/*                epsaf of each other, and the new point satisfies */
/*                fbest le ftry le fbest+epsaf. */

/*     quitI      is true when the interval of uncertainty is less than */
/*                2*tol. */

/*     quitS      is true as soon as alfa is too small to be useful; */
/*                i.e., btrue le alfsml. */

/*     xinxw      is true if xtry is in (xw,0) or (0,xw). */
/*     ------------------------------------------------------------------ */
    *imprvd = FALSE_;
    badfun = FALSE_;
    quitf = FALSE_;
    quitfz = FALSE_;
    quits = FALSE_;
    quiti = FALSE_;
    if (*first) {
/*        --------------------------------------------------------------- */
/*        First entry.  Initialize various quantities, check input data */
/*        and prepare to evaluate the function at the initial step alfa. */
/*        --------------------------------------------------------------- */
	*first = FALSE_;
	*numf = 0;
	*alfbst = 0.;
	badfun = *alfmax <= *toltny || *g0 >= 0.;
	*done = badfun;
	*moved = FALSE_;
	if (! (*done)) {
	    *braktd = FALSE_;
	    *crampd = *alfmax <= *tolabs;
	    *extrap = FALSE_;
	    *vset = FALSE_;
	    *wset = FALSE_;
	    *nsamea = 0;
	    *nsameb = 0;
	    *tolmax = *tolrel * *alfmax + *tolabs;
	    *a = 0.;
	    *b = *alfmax + *tolmax;
	    *fa = 0.;
	    *factor = 5.;
	    tol = *tolabs;
	    *xtry = *alfa;
	    if (*debug) {
		io___139.ciunit = *nout;
		s_wsfe(&io___139);
		do_fio(&c__1, (char *)&(*g0), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*tolabs), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*alfmax), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*targtg), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*tolrel), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*epsaf), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*crampd), (ftnlen)sizeof(logical));
		e_wsfe();
	    }
	}
    } else {
/*        --------------------------------------------------------------- */
/*        Subsequent entries.  The function has just been evaluated at */
/*        alfa = alfbst + xtry,  giving ftry. */
/*        --------------------------------------------------------------- */
	if (*debug) {
	    io___140.ciunit = *nout;
	    s_wsfe(&io___140);
	    do_fio(&c__1, (char *)&(*alfa), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ftry), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	++(*numf);
	++(*nsamea);
	++(*nsameb);
	if (! (*braktd)) {
	    *tolmax = *tolabs + *tolrel * *alfmax;
	    *b = *alfmax - *alfbst + *tolmax;
	}
/*        Check if xtry is in the interval (xw,0) or (0,xw). */
	if (*wset) {
	    xinxw = 0. < *xtry && *xtry <= *xw || *xw <= *xtry && *xtry < 0.;
	} else {
	    xinxw = FALSE_;
	}
	*imprvd = *ftry < *fbest;
	if (*vset) {
	    closef = (d__1 = *fbest - *fv, abs(d__1)) <= *epsaf;
	} else {
	    closef = FALSE_;
	}
	if (*imprvd) {
/*           We seem to have an improvement.  The new point becomes the */
/*           origin and other points are shifted accordingly. */
	    if (*wset) {
		*xv = *xw - *xtry;
		*fv = *fw;
		*vset = TRUE_;
	    }
	    *xw = 0. - *xtry;
	    *fw = *fbest;
	    *wset = TRUE_;
	    *fbest = *ftry;
	    *alfbst = *alfa;
	    *moved = TRUE_;
	    *a -= *xtry;
	    *b -= *xtry;
	    *extrap = ! xinxw;
/*           Decrease the length of (a,b). */
	    if (*xtry >= 0.) {
		*a = *xw;
		*fa = *fw;
		*nsamea = 0;
	    } else {
		*b = *xw;
		*nsameb = 0;
		*braktd = TRUE_;
	    }
	} else if (closef && *ftry - *fbest < *epsaf) {
/*           Quit if there has been no progress and ftry, fbest, fw */
/*           and fv are all within epsaf of each other. */
	    quitfz = TRUE_;
	} else {
/*           The new function value is no better than the current best */
/*           point.  xtry must an end point of the new (a,b). */
	    if (*xtry < 0.) {
		*a = *xtry;
		*fa = *ftry;
		*nsamea = 0;
	    } else {
		*b = *xtry;
		*nsameb = 0;
		*braktd = TRUE_;
	    }
/*           The origin remains unchanged but xtry may qualify as xw. */
	    if (*wset) {
		if (*ftry < *fw) {
		    *xv = *xw;
		    *fv = *fw;
		    *vset = TRUE_;
		    *xw = *xtry;
		    *fw = *ftry;
		    if (*moved) {
			*extrap = xinxw;
		    }
		} else if (*moved) {
		    if (*vset) {
			setxv = *ftry < *fv || ! (*extrap);
		    } else {
			setxv = TRUE_;
		    }
		    if (setxv) {
			if (*vset && xinxw) {
			    *xw = *xv;
			    *fw = *fv;
			}
			*xv = *xtry;
			*fv = *ftry;
			*vset = TRUE_;
		    }
		} else {
		    *xw = *xtry;
		    *fw = *ftry;
		}
	    } else {
		*xw = *xtry;
		*fw = *ftry;
		*wset = TRUE_;
	    }
	}
/*        --------------------------------------------------------------- */
/*        Check the termination criteria. */
/*        --------------------------------------------------------------- */
	tol = *tolabs + *tolrel * *alfbst;
	truea = *alfbst + *a;
	trueb = *alfbst + *b;
	found = *moved && (d__1 = *fa - *fbest, abs(d__1)) <= -(*a) * *targtg;
	quitf = *numf >= *maxf;
	quiti = *b - *a <= tol + tol;
	quits = trueb <= *alfsml;
	if (quiti && ! (*moved)) {
/*           The interval of uncertainty appears to be small enough, */
/*           but no better point has been found.  Check that changing */
/*           alfa by b-a changes f by less than epsaf. */
	    tol /= 10.;
	    *tolabs = tol;
	    quiti = abs(*fw) <= *epsaf || tol <= *toltny;
	}
	*done = quitf || quitfz || quits || quiti || found;
	if (*debug) {
	    io___147.ciunit = *nout;
	    s_wsfe(&io___147);
	    do_fio(&c__1, (char *)&truea, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&trueb, (ftnlen)sizeof(doublereal));
	    d__1 = *b - *a;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*nsamea), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nsameb), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*numf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*braktd), (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&(*extrap), (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&closef, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&(*imprvd), (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&found, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&quiti, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&quitfz, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&quits, (ftnlen)sizeof(logical));
	    do_fio(&c__1, (char *)&(*alfbst), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*fbest), (ftnlen)sizeof(doublereal));
	    d__2 = *alfbst + *xw;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*fw), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    if (*vset) {
		io___148.ciunit = *nout;
		s_wsfe(&io___148);
		d__1 = *alfbst + *xv;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*fv), (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
/*        --------------------------------------------------------------- */
/*        Proceed with the computation of an estimate of a minimizer. */
/*        The choices are... */
/*        1. Parabolic fit using function values only. */
/*        2. Damped parabolic fit if the regular fit appears to be */
/*           consistently overestimating the distance to a minimizer. */
/*        3. Bisection, geometric bisection, or a step of tol if the */
/*           parabolic fit is unsatisfactory. */
/*        --------------------------------------------------------------- */
	if (! (*done)) {
	    xmidpt = (*a + *b) * .5;
	    s = 0.;
	    q = 0.;
/*           ============================================================ */
/*           Fit a parabola. */
/*           ============================================================ */
/*           See if there are two or three points for the parabolic fit. */
	    gw = (*fw - *fbest) / *xw;
	    if (*vset && *moved) {
/*              Three points available.  Use fbest, fw and fv. */
		gv = (*fv - *fbest) / *xv;
		s = gv - *xv / *xw * gw;
		q = (gv - gw) * 2.;
		if (*debug) {
		    io___154.ciunit = *nout;
		    s_wsfe(&io___154);
		    e_wsfe();
		}
	    } else {
/*              Only two points available.  Use fbest, fw and g0. */
		if (*moved) {
		    s = *g0 - gw * 2.;
		} else {
		    s = *g0;
		}
		q = (*g0 - gw) * 2.;
		if (*debug) {
		    io___155.ciunit = *nout;
		    s_wsfe(&io___155);
		    e_wsfe();
		}
	    }
/*           ------------------------------------------------------------ */
/*           Construct an artificial interval (artifa, artifb) in which */
/*           the new estimate of the steplength must lie.  Set a default */
/*           value of  xtry  that will be used if the polynomial fit is */
/*           rejected. In the following, the interval (a,b) is considered */
/*           the sum of two intervals of lengths  dtry  and  daux, with */
/*           common end point the best point (zero).  dtry is the length */
/*           of the interval into which the default xtry will be placed */
/*           and endpnt denotes its non-zero end point.  The magnitude of */
/*           xtry is computed so that the exponents of dtry and daux are */
/*           approximately bisected. */
/*           ------------------------------------------------------------ */
	    artifa = *a;
	    artifb = *b;
	    if (! (*braktd)) {
/*              A minimizer has not yet been bracketed. */
/*              Set an artificial upper bound by expanding the interval */
/*              xw  by a suitable factor. */
		*xtry = -(*factor) * *xw;
		artifb = *xtry;
		if (*alfbst + *xtry < *alfmax) {
		    *factor *= 5.;
		}
	    } else if (*vset && *moved) {
/*              Three points exist in the interval of uncertainty. */
/*              Check if the points are configured for an extrapolation */
/*              or an interpolation. */
		if (*extrap) {
/*                 The points are configured for an extrapolation. */
		    if (*xw < 0.) {
			endpnt = *b;
		    }
		    if (*xw > 0.) {
			endpnt = *a;
		    }
		} else {
/*                 If the interpolation appears to be overestimating the */
/*                 distance to a minimizer,  damp the interpolation step. */
		    if (*nsamea >= 3 || *nsameb >= 3) {
			*factor /= 5.;
			s = *factor * s;
		    } else {
			*factor = 1.;
		    }
/*                 The points are configured for an interpolation.  The */
/*                 artificial interval will be just (a,b).  Set endpnt so */
/*                 that xtry lies in the larger of the intervals (a,b) */
/*                 and  (0,b). */
		    if (xmidpt > 0.) {
			endpnt = *b;
		    } else {
			endpnt = *a;
		    }
/*                 If a bound has remained the same for three iterations, */
/*                 set endpnt so that  xtry  is likely to replace the */
/*                 offending bound. */
		    if (*nsamea >= 3) {
			endpnt = *a;
		    }
		    if (*nsameb >= 3) {
			endpnt = *b;
		    }
		}
/*              Compute the default value of  xtry. */
		dtry = abs(endpnt);
		daux = *b - *a - dtry;
		if (daux >= dtry) {
		    *xtry = dtry * 5. * (dtry / daux + .1) / 11.;
		} else {
		    *xtry = sqrt(daux) * .5 * sqrt(dtry);
		}
		if (endpnt < 0.) {
		    *xtry = -(*xtry);
		}
		if (*debug) {
		    io___161.ciunit = *nout;
		    s_wsfe(&io___161);
		    do_fio(&c__1, (char *)&(*xtry), (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&daux, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&dtry, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
/*              If the points are configured for an extrapolation set the */
/*              artificial bounds so that the artificial interval lies */
/*              within (a,b).  If the polynomial fit is rejected,  xtry */
/*              will remain at the relevant artificial bound. */
		if (*extrap) {
		    if (*xtry <= 0.) {
			artifa = *xtry;
		    } else {
			artifb = *xtry;
		    }
		}
	    } else {
/*              The gradient at the origin is being used for the */
/*              polynomial fit.  Set the default xtry to one tenth xw. */
		if (*extrap) {
		    *xtry = -(*xw);
		} else {
		    *xtry = *xw / 10.;
		}
		if (*debug) {
		    io___162.ciunit = *nout;
		    s_wsfe(&io___162);
		    do_fio(&c__1, (char *)&(*xtry), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfe();
		}
	    }
/*           ------------------------------------------------------------ */
/*           The polynomial fits give (s/q)*xw as the new step.  Reject */
/*           this step if it lies outside (artifa, artifb). */
/*           ------------------------------------------------------------ */
	    if (q != 0.) {
		if (q < 0.) {
		    s = -s;
		}
		if (q < 0.) {
		    q = -q;
		}
		if (s * *xw >= q * artifa && s * *xw <= q * artifb) {
/*                 Accept the polynomial fit. */
		    if ((d__1 = s * *xw, abs(d__1)) >= q * tol) {
			*xtry = s / q * *xw;
		    } else {
			*xtry = 0.;
		    }
		    if (*debug) {
			io___163.ciunit = *nout;
			s_wsfe(&io___163);
			do_fio(&c__1, (char *)&(*xtry), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }
/*     ================================================================== */
    if (! (*done)) {
	*alfa = *alfbst + *xtry;
	if (*braktd || *alfa < *alfmax - *tolmax) {
/*           The function must not be evaluated too close to a or b. */
/*           (It has already been evaluated at both those points.) */
	    xmidpt = (*a + *b) * .5;
	    if (*xtry <= *a + tol || *xtry >= *b - tol) {
		if (xmidpt <= 0.) {
		    *xtry = -tol;
		} else {
		    *xtry = tol;
		}
	    }
	    if (abs(*xtry) < tol) {
		if (xmidpt <= 0.) {
		    *xtry = -tol;
		} else {
		    *xtry = tol;
		}
	    }
	    *alfa = *alfbst + *xtry;
	} else {
/*           The step is close to or larger than alfmax, replace it by */
/*           alfmax to force evaluation of the function at the boundary. */
	    *braktd = TRUE_;
	    *xtry = *alfmax - *alfbst;
	    *alfa = *alfmax;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Exit. */
/*     ------------------------------------------------------------------ */
    if (*done) {
	if (badfun) {
	    *iexit = 8;
	} else if (quits) {
	    *iexit = 5;
	} else if (found) {
	    if (*alfbst < *alfmax) {
		*iexit = 1;
	    } else {
		*iexit = 2;
	    }
	} else if (*moved) {
	    *iexit = 3;
	} else if (quitf) {
	    *iexit = 7;
	} else if (*crampd) {
	    *iexit = 4;
	} else {
	    *iexit = 6;
	}
    }
    if (*debug) {
	io___164.ciunit = *nout;
	s_wsfe(&io___164);
	e_wsfe();
    }
    return 0;
} /* lsrchq_ */

