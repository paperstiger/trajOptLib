/* ./src/sn80ncon.f -- translated by f2c (version 20100827).
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
static doublereal c_b5 = 1.;
static doublereal c_b6 = -1.;
static integer c__3 = 3;
static integer c__31 = 31;
static integer c__22 = 22;
static integer c__23 = 23;
static integer c__33 = 33;
static integer c__6 = 6;
static integer c__26 = 26;
static integer c__25 = 25;
static integer c__4 = 4;
static integer c__21 = 21;
static integer c__5 = 5;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn80ncon.f */

/*     s8feas   s8FD     s8Fv     s8Fx     s8Gcpy   s8getR   s8Gloc */
/*     s8Gprd   s8Infs   s8iQN    s8iQP    s8iQP2   s8mrt    s8PPHx */
/*     s8qpHx   s8rand   s8rc     s8sclJ   s8sInf   s8step   s8sOpt */
/*     s8wInf */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s8feas_(integer *iexit, U_fp mnrlog, integer *lenr, 
	integer *m, integer *maxs, integer *mbs, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, integer *nnl0, integer *nnl, integer 
	*ndegen, integer *ns, integer *numlc, integer *numliq, integer *itn, 
	integer *itnlim, integer *itqp, integer *mnrprt, doublereal *sclobj, 
	doublereal *tolqp, doublereal *tolx, integer *ninf, doublereal *sinf, 
	doublereal *wtinf, doublereal *pinorm, doublereal *rgnorm, integer *
	ne, integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	integer *hetype, integer *hestat, integer *hfeas, integer *hs, 
	integer *kbs, doublereal *ascale, doublereal *bl, doublereal *bu, 
	doublereal *blsav, doublereal *busav, doublereal *blbs, doublereal *
	bubs, doublereal *gbs, doublereal *gqp, doublereal *hdx, doublereal *
	pbs, doublereal *pi, doublereal *r__, doublereal *rc, doublereal *rg, 
	doublereal *qprhs, doublereal *x0, doublereal *x, doublereal *xbs, 
	integer *iy, integer *iy1, doublereal *y, doublereal *y1, doublereal *
	y2, char *cw, integer *lencw, integer *iw, integer *leniw, doublereal 
	*rw, integer *lenrw, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_8000[] = "(\002 Itn\002,i7,\002: Feasible linear rows"
	    "\002)";
    static char fmt_8100[] = "(\002 Itn\002,i7,\002: PP\002,i1,\002.  Minimi"
	    "zing  Norm(x-x0)\002)";
    static char fmt_8200[] = "(\002 Itn\002,i7,\002: PP\002,i1,\002.  Norm(x"
	    "-x0) approximately minimized  (\002,1p,e8.2,\002)\002)";
    static char fmt_8300[] = "(\002 Itn\002,i7,\002: PP1.  Making nonlinear "
	    "variables feasible\002)";
    static char fmt_8400[] = "(\002 Itn\002,i7,\002: PP1. \002,i7,\002 nonli"
	    "near variables made feasible\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, itqptargt;
    static doublereal x0j, blj, buj;
    static char str[80];
    static doublereal eps0, eps2;
    extern /* Subroutine */ int s5lp_(integer *, integer *, char *, logical *,
	     integer *, U_fp, logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s5qp_(integer *, integer *, char *, logical *, integer *, U_fp, 
	    U_fp, U_fp, logical *, logical *, integer *, logical *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal obja;
    static logical gotr;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static logical needx;
    static doublereal objpp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal tolfp;
    static integer nviol, nobjp0;
    extern /* Subroutine */ int s2aprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static doublereal hcndbd;
    extern /* Subroutine */ int s8pphx_(), s8qphx_();
    static doublereal zcndbd;
    static integer lemode;
    static logical elastc, needlu;
    static integer iobjpp, inform__, lvlinf, msbsav, nobjpp, minimz, lvlppm, 
	    subopt;
    static doublereal tolqpp;
    static integer typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer hvcalls;
    static char probtag[20];
    static integer itqpmax;

    /* Fortran I/O blocks */
    static icilist io___19 = { 0, str, 0, fmt_8000, 80, 1 };
    static icilist io___20 = { 0, str, 0, fmt_8100, 80, 1 };
    static icilist io___26 = { 0, str, 0, fmt_8200, 80, 1 };
    static icilist io___27 = { 0, str, 0, fmt_8300, 80, 1 };
    static icilist io___28 = { 0, str, 0, fmt_8400, 80, 1 };
    static icilist io___39 = { 0, str, 0, fmt_8200, 80, 1 };


/*     ================================================================== */
/*     s8feas   finds a feasible point for a set of linear constraints. */
/*     A basis is assumed to be specified by nS, hs(*), x(*) and */
/*     kBS(m+1:m+nS).  In particular, there must be nS values hs(j) = 2, */
/*     and the corresponding j's must be listed in kBS(m+1:m+nS). */
/*     The ordering in kBS matches the reduced Hessian R (if any). */

/*     On entry, blSav and blSav contain copies of the true (possibly */
/*     scaled) upper and bounds set in s5getB. */

/*      iExit       Result */
/*      -----       ------ */
/*       >0         Fatal error */
/*        0         Feasible point found */

/*     11 May 1994: First version of s8feas. */
/*     19 Aug 1996: First minsum version. */
/*     05 Feb 1998: Proximal point norm changed to one-norm. */
/*     23 Dec 1999: Optional Proximal Point methods 0 and  2 added. */
/*     03 Aug 2003: snPRNT and snEXIT adopted. */
/*     16 May 2006: Explicit target itQP added. */
/*     18 Jun 2008: Hdx, pBS and rg added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* max # of new superbasics */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --pi;
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --y2;
    --y1;
    --y;
    --iy1;
    --iy;
    --x;
    --x0;
    --rc;
    --pbs;
    --busav;
    --blsav;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --qprhs;
    --hdx;
    --gqp;
    --jcol;
    --indj;
    --locj;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    lvlppm = iw[79];
/* 1(2)-norm proximal point method for x0 */
    eps0 = rw[2];
/* eps**(4/5) */
    eps2 = rw[4];
/* eps**(1/2) */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    obja = 0.;
    iobjpp = 0;
    minimz = 1;
/* Local value */
    lemode = 1;
/* Enter elastic mode if infeasible */
    lvlinf = 2;
/*   W1 = 0  W2 = 1   W1*true obj + W2*sInfE */
/* In elastic mode, use: */
    elastc = FALSE_;
    needlu = TRUE_;
    needx = needlu;
/*     Set the LP rhs to make x satisfy the (relaxed) nonlinear rows. */
/*     The array  QPrhs  contains the rhs. */
/*     Use a fairly tight optimality tolerance for phase 1. */
    if (*nncon > 0) {
	dcopy_(nncon, &x[*n + 1], &c__1, &qprhs[1], &c__1);
	s2aprd_(&c__0, &eps0, ne, nlocj, &locj[1], &indj[1], &jcol[1], &c_b5, 
		&x[1], n, &c_b6, &qprhs[1], nncon);
    }
    if (*numliq > 0 || *ninf > 0) {
/*        --------------------------------------------------------------- */
/*        Find a feasible point for the linear constraints. */
/*        If none exists, minimize the sum of infeasibilities of the */
/*        linear rows, subject to the column bounds. */
/*        --------------------------------------------------------------- */
	iload_(numlc, &c__3, &hetype[*n + *nncon + 1], &c__1);
	s_copy(probtag, "linear rows", (ftnlen)20, (ftnlen)11);
	subopt = -1;
	tolfp = eps2;
	s5lp_(&inform__, &c__0, probtag, &elastc, &subopt, (U_fp)mnrlog, &
		needlu, &needx, m, n, nb, ndegen, itqp, itnlim, itn, &lemode, 
		&lvlinf, mnrprt, &minimz, &iobjpp, sclobj, &obja, &tolfp, 
		tolqp, tolx, ninf, sinf, wtinf, pinorm, rgnorm, ne, nlocj, &
		locj[1], &indj[1], &jcol[1], &hetype[1], &hestat[1], &hfeas[1]
		, &hs[1], &kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1], &
		bubs[1], &gbs[1], &pi[1], &rc[1], nncon0, nncon, &qprhs[1], &
		x[1], &xbs[1], &x0[1], &iy[1], &iy1[1], &y[1], &y1[1], cw + 8,
		 lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)20, (ftnlen)8);
/*        Check for trouble in s5LP. */
/*        iExit        Status */
/*        -----        ------ */
/*         -3          Too many iterations */
/*         -2          Phase 1 is unbounded */
/*         -1          infeasible nonelastics */
/*          0          infeasibilities minimized */
/*         >0          Fatal error */

/*        If the linear constraints are infeasible, the sum of */
/*        infeasibilities will have been minimized. */
	if (inform__ != 0 || *ninf > 0) {
	    if (inform__ > 0) {
		*iexit = inform__;
/* Fatal error */
	    } else if (inform__ == -3) {
		*iexit = 31;
/* iterations limit */
	    } else if (*ninf > 0) {
		*iexit = 11;
/* infeasible linear constraints */
	    }
	    if (*iexit != 0) {
		goto L800;
	    }
	}
/*        Now the linear rows are feasible, they are never allowed */
/*        to be infeasible again. */
	iload_(numlc, &c__0, &hetype[*n + *nncon + 1], &c__1);
/*        Print something brief if s5LP didn't already do so. */
	if (*mnrprt >= 1) {
	    s_wsfi(&io___19);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
	    snprnt_(&c__22, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    if (lvlppm > 0 && *nnl > 0) {
/*        =============================================================== */
/*        x  is feasible for the linear constraints. */
/*        Find a feasible point closest to x0. */
/*        Minimize norm(x - x0). */
/*        =============================================================== */
	if (*mnrprt >= 1) {
	    s_wsfi(&io___20);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	}
	if (lvlppm == 1) {
/*           ------------------------------------------------------------ */
/*           Minimize the one-norm of (x-x0) by fixing the nonlinear */
/*           variables so that bl = x0 = bu.  Any bl or bu that is moved */
/*           to  x0  is made elastic. */
/*           ------------------------------------------------------------ */
	    i__1 = *nnl;
	    for (j = 1; j <= i__1; ++j) {
		blj = bl[j];
		buj = bu[j];
		if (blj == buj) {
/*                 Relax */
		} else {
		    x0j = x0[j];
		    bl[j] = x0j;
		    bu[j] = x0j;
		    hetype[j] = 3;
		    if (hs[j] <= 1) {
			x[j] = x0j;
		    }
		}
	    }
	    s_copy(probtag, "norm(x-x0) problem  ", (ftnlen)20, (ftnlen)20);
	    iw[223] = 1;
/* New LP print   header */
	    iw[225] = 1;
/* New LP summary header */
	    needx = TRUE_;
	    subopt = -1;
	    tolfp = .01;
/* Sloppy phase 1 optimality tol for PP. */
	    s5lp_(&inform__, &c__0, probtag, &elastc, &subopt, (U_fp)mnrlog, &
		    needlu, &needx, m, n, nb, ndegen, itqp, itnlim, itn, &
		    lemode, &lvlinf, mnrprt, &minimz, &iobjpp, sclobj, &obja, 
		    &tolfp, tolqp, tolx, ninf, sinf, wtinf, pinorm, rgnorm, 
		    ne, nlocj, &locj[1], &indj[1], &jcol[1], &hetype[1], &
		    hestat[1], &hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1],
		     &bu[1], &blbs[1], &bubs[1], &gbs[1], &pi[1], &rc[1], 
		    nncon0, nncon, &qprhs[1], &x[1], &xbs[1], &x0[1], &iy[1], 
		    &iy1[1], &y[1], &y1[1], cw + 8, lencw, &iw[1], leniw, &rw[
		    1], lenrw, (ftnlen)20, (ftnlen)8);
/*           Some elastic variables may have moved outside their bounds. */
/*           Count them.  Reset the true bounds. */
/*           If necessary,  get feasible again with the normal tolQP. */
	    nviol = 0;
	    i__1 = *nnl;
	    for (j = 1; j <= i__1; ++j) {
		bl[j] = blsav[j];
		bu[j] = busav[j];
		if (x[j] < bl[j] - *tolx || x[j] > bu[j] + *tolx) {
		    ++nviol;
		}
	    }
/*           Check for errors in s5LP. */
/*           inform values are = -3,-2,-1, 0, >0 */
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal error */
		} else if (inform__ == -3) {
		    *iexit = 31;
/* iterations limit */
		}
		if (*iexit != 0) {
		    goto L800;
		}
	    }
	    if (inform__ == 0 && *mnrprt >= 1) {
		s_wsfi(&io___26);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__33, str, &iw[1], leniw, (ftnlen)80);
		if (nviol > 0) {
		    s_wsfi(&io___27);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	    if (nviol > 0) {
		s_copy(probtag, "linear rows again   ", (ftnlen)20, (ftnlen)
			20);
		elastc = FALSE_;
		needx = TRUE_;
		subopt = -1;
		tolfp = eps2;
/* Revert to accurate phase 1 opt tol */
		if (inform__ != 0) {
		    needlu = TRUE_;
		}
		s5lp_(&inform__, &c__0, probtag, &elastc, &subopt, (U_fp)
			mnrlog, &needlu, &needx, m, n, nb, ndegen, itqp, 
			itnlim, itn, &lemode, &lvlinf, mnrprt, &minimz, &
			iobjpp, sclobj, &obja, &tolfp, tolqp, tolx, ninf, 
			sinf, wtinf, pinorm, rgnorm, ne, nlocj, &locj[1], &
			indj[1], &jcol[1], &hetype[1], &hestat[1], &hfeas[1], 
			&hs[1], &kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1],
			 &bubs[1], &gbs[1], &pi[1], &rc[1], nncon0, nncon, &
			qprhs[1], &x[1], &xbs[1], &x0[1], &iy[1], &iy1[1], &y[
			1], &y1[1], cw + 8, lencw, &iw[1], leniw, &rw[1], 
			lenrw, (ftnlen)20, (ftnlen)8);
/*              Possible inform values are = -3,-2,-1, 0, >0 */
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
/* Fatal error */
		    } else if (inform__ == -3) {
			*iexit = 31;
/* iterations limit */
		    } else if (*ninf > 0) {
			*iexit = 11;
/* infeasible (should not happen here) */
		    }
		    if (*iexit != 0) {
			goto L800;
		    }
		}
		if (inform__ == 0 && *mnrprt >= 1) {
		    s_wsfi(&io___28);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&nviol, (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	    *ninf = 0;
	    *sinf = 0.;
/*           Now the nonlinear variables are feasible, they are never */
/*           allowed to be infeasible again. */
	    iload_(nnl, &c__0, &hetype[1], &c__1);
	} else if (lvlppm == 2) {
/*           ------------------------------------------------------------ */
/*           Minimize the two-norm of (x-x0). */
/*           ------------------------------------------------------------ */
/*           Now the linear rows are feasible, they are never allowed */
/*           to be infeasible again. */
	    iload_(numlc, &c__0, &hetype[*n + *nncon + 1], &c__1);
	    nobjpp = 0;
/* No explicit gradient in proximal point */
	    nobjp0 = 1;
	    gotr = FALSE_;
	    needlu = FALSE_;
	    typelu = 3;
	    hvcalls = 0;
	    s_copy(probtag, "norm(x-x0) problem  ", (ftnlen)20, (ftnlen)20);
	    iw[223] = 1;
/* Switch to QP print   heading */
	    iw[225] = 1;
/* Switch to QP summary heading */
	    needx = FALSE_;
	    itqpmax = 100;
/* Limit the number of minor iterations */
	    itqptargt = 100;
	    msbsav = iw[95];
	    iw[95] = 100;
/* and the number of new superbasics */
	    subopt = 0;
	    tolfp = eps2;
	    tolqpp = .01;
/* Sloppy phase 2 opt tol */
	    s5qp_(&inform__, &c__6, probtag, &elastc, &subopt, (U_fp)s8pphx_, 
		    (U_fp)s8qphx_, (U_fp)mnrlog, &gotr, &needlu, &typelu, &
		    needx, lenr, m, maxs, mbs, n, nb, ndegen, &hvcalls, nnl0, 
		    nnl, &nobjp0, &nobjpp, nnl0, nnl, ns, itqp, &itqpmax, &
		    itqptargt, itn, &lemode, &lvlinf, mnrprt, &minimz, &
		    iobjpp, sclobj, &obja, &objpp, &hcndbd, &zcndbd, &tolfp, &
		    tolqpp, tolx, ninf, sinf, wtinf, pinorm, rgnorm, ne, 
		    nlocj, &locj[1], &indj[1], &jcol[1], &hetype[1], &hestat[
		    1], &hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1], &bu[1]
		    , &blbs[1], &bubs[1], &gbs[1], &gqp[1], &gqp[1], &hdx[1], 
		    &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], nncon0, nncon, &
		    qprhs[1], nnl0, nnl, &x0[1], &x[1], &xbs[1], &x0[1], &iy[
		    1], &iy1[1], &y[1], &y1[1], &y2[1], cw + 8, lencw, &iw[1],
		     leniw, &rw[1], lenrw, cw + 8, lencw, &iw[1], leniw, &rw[
		    1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
	    iw[95] = msbsav;
/*           Check for trouble. */
/*           Possible inform values are = -8(-1)-1, 0, >0 */
	    if (inform__ != 0 || *ninf > 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal LU error */
		} else if (inform__ == -3) {
		    *iexit = 31;
/* iterations limit */
		} else if (inform__ == -1 || *ninf > 0) {
		    *iexit = 11;
/* infeasible (should not happen here) */
		}
		if (*iexit != 0) {
		    goto L800;
		}
	    }
/*           Note: ObjQP is an updated quantity that may be slightly */
/*           negative. */
	    if (*mnrprt >= 1) {
		s_wsfi(&io___39);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&lvlppm, (ftnlen)sizeof(integer));
		d__1 = abs(objpp);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
		snprnt_(&c__22, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
/* Proximal Point method 1 */
    }
/* nnL > 0 */
L800:
    return 0;
} /* s8feas_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8feas */
/* Subroutine */ int s8fd_(integer *nncon0, integer *nncon, integer *nnobj, 
	integer *itn, integer *cditns, logical *centrl, logical *goodg, 
	logical *newg, logical *usefd, integer *info, doublereal *duinf, 
	doublereal *fcon, doublereal *fobj, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002 -- Central differences i"
	    "nvoked.\002,\002  Small reduced gradient.\002)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    static doublereal epsrf, cnorm, fdint1;
    extern doublereal dnrm1s_(integer *, doublereal *, integer *);
    static doublereal objsiz, rgnorm, rgtest;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___47 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s8FD   controls the switch from forward to central differences and */
/*     vice versa. */

/*     If the forward-difference estimate of the reduced gradient of the */
/*     Lagrangian is small,  a switch is made to central differences. */
/*     In this case, the derivatives are recomputed and the QP is solved */
/*     again. */

/*     On the other hand, if central differences have produced a large */
/*     reduced-gradient norm, switch back to forward differences. */

/*     31 Mar 2000: First version of s8FD written for SNOPT 6.1. */
/*     03 Aug 2003: snPRNT adopted. */
/*     03 Aug 2003: Current version of s8FD. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* forwd diffs or cntrl diffs */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fcon;
    --info;
    --iw;
    --rw;

    /* Function Body */
    epsrf = rw[73];
/* relative function precision. */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    cnorm = 0.;
    if (*nncon > 0) {
	cnorm = dnrm1s_(nncon, &fcon[1], &c__1);
    }
    objsiz = 0.;
    if (*nnobj > 0) {
	objsiz = abs(*fobj);
    }
    *goodg = TRUE_;
    rgtest = (objsiz + 1. + cnorm) * epsrf / fdint1;
    rgnorm = *duinf;
    if (*centrl) {
	if (rgnorm > rgtest * 10. && *cditns > 0) {
	    iw[181] = 1;
	    *centrl = FALSE_;
	    if (*usefd) {
		info[6] = 0;
	    }
	}
    } else {
	if (rgnorm <= rgtest) {
	    *cditns = 0;
	    iw[181] = 2;
	    if (*usefd) {
		*goodg = FALSE_;
		*newg = TRUE_;
		info[6] = 1;
		s_wsfi(&io___47);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
    return 0;
} /* s8fd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8FD */
/* Subroutine */ int s8fv_(logical *elastc, integer *n, integer *nncon, 
	doublereal *tolz, doublereal *wtinf, doublereal *bl, doublereal *bu, 
	doublereal *fv, doublereal *x, doublereal *ycon, doublereal *fx)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, blj, buj, fvi, fxi, fvl, fvu, yconi, yconv, yconw;

/*     ================================================================== */
/*     s8Fv  computes the vector of nonlinear constraint violations: */
/*        Fv = fCon + A(linear)*x - (nonlinear slacks) */

/*     If the Lagrange multiplier is zero, the violation can be set to */
/*     any value without changing the merit function.  In this case we */
/*     try and set the slack so that Fv is zero (subject to the slack */
/*     being feasible). */

/*     In elastic mode we implicitly adjust the variables v and w such */
/*     that   c - s(feas) + v - w = 0,  with  v >= 0  and  w >= 0. */

/*     On entry, */
/*        x   =  the current x. */
/*        Fx  =  fCon + A(linear)*x,   defined in s8Fx. */

/*     On exit, */
/*        x   =  x containing the modified slacks. */
/*        Fv  =  fCon + A(linear)*x -  slacks. */
/*        Fx  =  unaltered. */

/*     19 Apr 2001: First version based on s8sOpt */
/*     19 Apr 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --fx;
    --ycon;
    --x;
    --fv;
    --bu;
    --bl;

    /* Function Body */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	xj = x[j];
	fxi = fx[i__];
	fvi = fxi - xj;
	yconi = ycon[i__];
	blj = bl[j];
	buj = bu[j];
	fvu = fxi - buj;
	fvl = fxi - blj;
	yconv = (d__1 = *wtinf - yconi, abs(d__1));
/* Multiplier for v in elastic mode */
	yconw = (d__1 = *wtinf + yconi, abs(d__1));
/* Multiplier for w in elastic mode */
	if (*elastc && xj <= blj && yconv <= *tolz) {
	    if (fvi > 0.) {
		fvi = max(0.,fvl);
	    } else {
		fvi = 0.;
	    }
	} else if (*elastc && xj >= buj && yconw <= *tolz) {
	    if (fvi < 0.) {
		fvi = min(0.,fvu);
	    } else {
		fvi = 0.;
	    }
	} else {
	    if (yconi <= *tolz && fvi > 0.) {
		fvi = max(0.,fvu);
	    } else if (yconi >= -(*tolz) && fvi < 0.) {
		fvi = min(0.,fvl);
	    }
	}
	xj = fxi - fvi;
	fv[i__] = fvi;
	x[j] = xj;
    }
    return 0;
} /* s8fv_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Fv */
/* Subroutine */ int s8fx_(integer *n, integer *nncon, integer *nnjac, 
	doublereal *tolz, integer *ne, integer *nlocj, integer *locj, integer 
	*indj, doublereal *jcol, doublereal *fcon, doublereal *x, doublereal *
	fx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer nlin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2aprd_(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *);

/*     ================================================================== */
/*     s8Fx  defines the nonlinear constraint values */
/*       Fx  =  true nonlinear slack = fCon + A(linear)*x, */

/*     09 Jan 1992: First version based on Minos routine m8viol. */
/*     16 Nov 1998: Norm x changed to include only columns. */
/*     21 Oct 2000: Made compatible with SNOPT 6.1 */
/*     21 Oct 2000: Current version of s8Fx */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Compute the nonlinear constraint value. */
/*     Set  Fx  =  fCon + (linear A)*x,   excluding slacks. */
    /* Parameter adjustments */
    --fx;
    --x;
    --fcon;
    --jcol;
    --indj;
    --locj;

    /* Function Body */
    dcopy_(nncon, &fcon[1], &c__1, &fx[1], &c__1);
    nlin = *n - *nnjac;
    if (nlin > 0) {
	i__1 = nlin + 1;
	s2aprd_(&c__0, tolz, ne, &i__1, &locj[*nnjac + 1], &indj[1], &jcol[1],
		 &c_b5, &x[*nnjac + 1], &nlin, &c_b5, &fx[1], nncon);
    }
    return 0;
} /* s8fx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Fx */
/* Subroutine */ int s8gcpy_(integer *nncon, integer *nnjac, integer *ne, 
	integer *nlocj, integer *locj, integer *indj, integer *neg1, integer *
	nlocg1, integer *locg1, doublereal *g1, integer *neg2, integer *
	nlocg2, integer *locg2, doublereal *g2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l1, l2, ir;

/*     ================================================================== */
/*     s8Gcpy  copies G1 into G2 when either  G1 or  G2 */
/*     is stored in the upper-left hand corner of J. */

/*     16 Sep 1993: First version. */
/*     26 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --indj;
    --locj;
    --g1;
    --locg1;
    --g2;
    --locg2;

    /* Function Body */
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	l1 = locg1[j];
	l2 = locg2[j];
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
		goto L100;
	    }
	    g2[l2] = g1[l1];
	    ++l1;
	    ++l2;
	}
L100:
	;
    }
    return 0;
} /* s8gcpy_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gcpy */
/* Subroutine */ int s8getr_(integer *iexit, U_fp hprod, U_fp hprod1, integer 
	*hqntype, integer *hvcalls, logical *gotr, integer *typelu, integer *
	lureq, integer *itn, integer *lenr, integer *m, integer *mbs, integer 
	*n, integer *nb, integer *nncon0, integer *nncon, integer *nnh, 
	integer *ns, integer *mjrprt, integer *minimz, integer *iobj, 
	doublereal *u0ii, doublereal *targth, doublereal *targtz, integer *ne,
	 integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	integer *hs, integer *kbs, doublereal *bl, doublereal *bu, doublereal 
	*blbs, doublereal *bubs, doublereal *r__, doublereal *qprhs, 
	doublereal *xqp, doublereal *xbs, integer *iy, integer *iy1, 
	doublereal *y, doublereal *y1, doublereal *y2, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Reduced Hessian rese"
	    "t\002)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Indefinite reduced Hess"
	    "ian\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Ill-conditioned QP null"
	    "-space basis.\002,\002 Cond = \002,1p,e8.1)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int s8h0_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static doublereal eps;
    static char str[80];
    extern /* Subroutine */ int s5hz_(integer *, U_fp, U_fp, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , char *, integer *, integer *, integer *, doublereal *, integer *
	    , char *, integer *, integer *, integer *, doublereal *, integer *
	    , ftnlen, ftnlen);
    static integer eigh;
    static logical newb;
    static integer maxr;
    static logical luok;
    static doublereal hdmax, flmax;
    static integer nswap;
    static logical newlu;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), s5hfac_(integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static doublereal condz0;
    extern /* Subroutine */ int s5rchk_(integer *, U_fp, U_fp, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static logical rcheck, needlu;
    static integer inform__, lvlhes;
    static doublereal plinfy;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2trylu_(integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, integer *, doublereal *, integer 
	    *);

    /* Fortran I/O blocks */
    static icilist io___81 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___82 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___84 = { 0, str, 0, fmt_1300, 80, 1 };


/*     ================================================================== */
/*     s8getR   computes the Cholesky factor of the reduced Hessian. */

/*     On entry, the LU factorization is assumed to be known. */
/*       gotR = .false. */

/*     iExit      Result */
/*     -----      ------ */
/*       0        the reduced Hessian has been computed successfully */
/*      >0        fatal error */

/*     LUreq =  1  Frequency */
/*     LUreq =  2  LU nonzeros increased */
/*     LUreq =  3 */
/*     LUreq =  4 */
/*     LUreq =  5  Singular after LU mod */
/*     LUreq =  6  Unstable LU mod (growth in new column of U) */
/*     LUreq =  7  Not enough memory */
/*     LUreq =  8 */
/*     LUreq =  9 */
/*     LUreq = 10  Row error in setx */
/*     LUreq = 11  Big  dx   in setx */

/*     LUreq = 20 */
/*     LUreq = 21  Iterative refinement failed in QP */
/*     LUreq = 22  Unbounded QP */
/*     LUreq = 23 */
/*     LUreq = 24  Small directional derivative in QP */
/*     LUreq = 25  Ill-conditioned Z in QP */
/*     LUreq = 26  Indefinite Z'HZ in QP */
/*     LUreq = 27  R singular after bound swap in QP */

/*     On output, */
/*     QPerr points to ' ', 't', 'u' or 'w'. */
/*     QPfea points to ' '  or 'i'. */

/*     14 Mar 2001: First version. */
/*     03 Aug 2003: snPRNT adopted. */
/*     29 Jun 2005: Current version of s8getR. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* condition estimate of Z */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --xbs;
    --bubs;
    --blbs;
    --kbs;
    --y2;
    --y1;
    --y;
    --iy1;
    --iy;
    --xqp;
    --bu;
    --bl;
    --hs;
    --qprhs;
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
    maxr = iw[52];
/* max columns of R */
    lvlhes = iw[72];
/* 0,1,2 => LM, FM, Exact Hessian */
    eigh = iw[200];
/* =1(0) for definite QP Hessian */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    flmax = rw[8];
/* est. of the largest pos. real */
    plinfy = flmax;
    condz0 = plinfy;
/* save an initial estimate of cond(Z) */
    luok = TRUE_;
/*     ================================================================== */
/* +    while (LUok  .and. .not. gotR) do */
L100:
    if (luok && ! (*gotr)) {
/*     ------------------------------------------------------------------ */
	inform__ = 0;
	needlu = *lureq > 0;
	if (needlu) {
	    s2bfac_(iexit, typelu, &needlu, &newlu, &newb, iobj, itn, mjrprt, 
		    lureq, m, mbs, n, nb, nnh, ns, &nswap, ne, nlocj, &locj[1]
		    , &indj[1], &jcol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &
		    blbs[1], &bubs[1], nncon0, nncon, &qprhs[1], &xqp[1], &
		    xbs[1], &iy[1], &iy1[1], &y[1], &y1[1], &iw[1], leniw, &
		    rw[1], lenrw);
	    if (*iexit != 0) {
		goto L900;
	    }
	}
	if (*ns > 0) {
/* ----------------------------------------------------------- */
/* Compute and factorize  Z'HZ. */
/* ----------------------------------------------------------- */
	    s5hz_(&inform__, (U_fp)hprod, (U_fp)hprod1, &maxr, lenr, minimz, 
		    m, mbs, n, nb, nnh, ns, hvcalls, ne, nlocj, &locj[1], &
		    indj[1], &jcol[1], &hdmax, &rw[192], targtz, &kbs[1], &
		    r__[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], 
		    leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1]
		    , lenrw, (ftnlen)8, (ftnlen)8);
/* Check for trouble in s5Hz.  Possible exits are: */
/* inform  Status */
/* ------  ------ */
/*  -1     Ill-conditioned null-space basis */
/*   0     reduced Hessian computed successfully */
/*  >0     Fatal error in LU solve */
	    if (inform__ == 0) {
		s5hfac_(&inform__, &eigh, itn, lenr, m, &maxr, mbs, nb, ns, 
			targth, &hdmax, &hs[1], &kbs[1], &iy[1], &bl[1], &bu[
			1], &blbs[1], &bubs[1], &xqp[1], &xbs[1], &r__[1], &
			iw[1], leniw, &rw[1], lenrw);
/* Check for trouble in s5Hfac. */
/* inform    Status */
/* ------    ------ */
/*  -2      H singular (but should be positive definite) */
/*  -1      H indefinite */
/*   0      normal exit */
		if (inform__ != 0) {
/* The reduced Hessian is not positive definite. */
/* Reset the H = I for QN Hessian if it has not been */
/* done already.  Otherwise refactorize B,  possibly */
/* with tighter tols. */
		    if (*hqntype != 2 && lvlhes != 2) {
/* -------------------------------------------------- */
/* Set unit Hessian. */
/* Z'HZ must be computed again. */
/* -------------------------------------------------- */
			s_wsfi(&io___81);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
			*hqntype = 2;
			s8h0_(hqntype, nnh, u0ii, &iw[1], leniw, &rw[1], 
				lenrw);
		    } else {
/* -------------------------------------------------- */
/* H = I, or is exact, so Z must be ill-conditioned. */
/* Refactorize B with tighter factor tol. */
/* -------------------------------------------------- */
			s_wsfi(&io___82);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
			*targth = 1. / (eps * eps);
			s2trylu_(itn, &c__26, ns, lureq, &luok, typelu, &iw[1]
				, leniw, &rw[1], lenrw);
		    }
		}
		rcheck = FALSE_;
		if (rcheck) {
		    s5rchk_(iexit, (U_fp)hprod, (U_fp)hprod1, itn, minimz, &
			    maxr, lenr, m, mbs, n, nb, hvcalls, nnh, ns, ne, 
			    nlocj, &locj[1], &indj[1], &jcol[1], &kbs[1], &
			    r__[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[
			    1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
			    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		    if (*iexit != 0) {
			goto L900;
		    }
		}
	    } else if (inform__ == -1) {
/* Ill-conditioned Z in s5Hz. */
/* Refactorize B, possibly with a reduced factor tol. */
/* If factor tol is already tight, accept Z, however bad. */
/* To avoid repeated factorizations, accept Z if condZ */
/* wasn't even reduced by the last factorize. */
		s_wsfi(&io___84);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rw[192], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__25, ns, lureq, &luok, typelu, &iw[1], leniw,
			 &rw[1], lenrw);
		if (! luok) {
		    *targtz = plinfy;
		    luok = TRUE_;
		} else if (rw[192] < condz0) {
		    condz0 = rw[192];
		} else {
		    *targtz = plinfy;
		}
	    } else if (inform__ > 0) {
/* LU error in s5Hz */
		*iexit = inform__;
		goto L900;
	    }
	}
	*gotr = inform__ == 0;
	goto L100;
    }
/* +    end while */
/*     ------------------------------------------------------------------ */
    if (! (*gotr)) {
	*iexit = 44;
    }
/* unable to factor Z'Z */
L900:
    return 0;
} /* s8getr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8getR */
/* Subroutine */ int s8gloc_(integer *nncon, integer *nnjac, integer *ne, 
	integer *nlocj, integer *locj, integer *indj, integer *negcon, 
	integer *nlocg, integer *locg)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, ir, neg;

/*     ================================================================== */
/*     s8Gloc  counts the number of nonlinear Jacobian elements and */
/*     assembles their column pointers in locG. */

/*     29 Oct 2000: First version of s8Gloc. */
/*     31 Aug 2008: Local variable used for negCon. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --indj;
    --locj;
    --locg;

    /* Function Body */
    neg = 0;
    locg[1] = 1;
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
		goto L100;
	    }
	    ++neg;
	}
L100:
	locg[j + 1] = neg + 1;
    }
    return 0;
} /* s8gloc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gloc */
/* Subroutine */ int s8gprd_(integer *task, doublereal *tolz, integer *ne, 
	integer *nlocj, integer *locj, integer *indj, integer *negcon, 
	integer *nlocg, integer *locg, doublereal *gcon, doublereal *alpha, 
	doublereal *x, integer *lenx, doublereal *beta, doublereal *y, 
	integer *leny)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ig, ij, ir;
    static doublereal xj, sum, alphxj;

/*     ================================================================== */
/*     s8Gprd computes matrix-vector products involving J and x.  The */
/*     variable task specifies the operation to be performed as follows: */
/*       task = 'N' (normal)          y := alpha*J *x + beta*y, */
/*       task = 'T' (transpose)       y := alpha*J'*x + beta*y, */
/*     where alpha and beta are scalars, x and y are vectors and J is a */
/*     sparse matrix whose columns are in natural order. */

/*     26 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --indj;
    --locj;
    --gcon;
    --locg;
    --x;
    --y;

    /* Function Body */
    if (*alpha == 0. && *beta == 1.) {
	return 0;
    }
/*     First form  y := beta*y. */
    if (*beta != 1.) {
	if (*beta == 0.) {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = 0.;
	    }
	} else {
	    i__1 = *leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		y[i__] = *beta * y[i__];
	    }
	}
    }
    if (*alpha == 0.) {
/*        Relax */
    } else if (*alpha == -1.) {
	if (*task == 0) {
	    i__1 = *lenx;
	    for (j = 1; j <= i__1; ++j) {
		xj = x[j];
		if (abs(xj) > *tolz) {
		    ig = locg[j];
		    i__2 = locj[j + 1] - 1;
		    for (ij = locj[j]; ij <= i__2; ++ij) {
			ir = indj[ij];
			if (ir > *leny) {
			    goto L100;
			}
			y[ir] -= gcon[ig] * xj;
			++ig;
		    }
		}
L100:
		;
	    }
	} else if (*task == 1) {
	    i__1 = *leny;
	    for (j = 1; j <= i__1; ++j) {
		sum = y[j];
		ig = locg[j];
		i__2 = locj[j + 1] - 1;
		for (ij = locj[j]; ij <= i__2; ++ij) {
		    ir = indj[ij];
		    if (ir > *lenx) {
			goto L200;
		    }
		    sum -= gcon[ig] * x[ir];
		    ++ig;
		}
L200:
		y[j] = sum;
	    }
	}
    } else {
/* General alpha */
	if (*task == 0) {
	    i__1 = *lenx;
	    for (j = 1; j <= i__1; ++j) {
		alphxj = *alpha * x[j];
		if (abs(alphxj) > *tolz) {
		    ig = locg[j];
		    i__2 = locj[j + 1] - 1;
		    for (ij = locj[j]; ij <= i__2; ++ij) {
			ir = indj[ij];
			if (ir > *leny) {
			    goto L300;
			}
			y[ir] += gcon[ig] * alphxj;
			++ig;
		    }
		}
L300:
		;
	    }
	} else if (*task == 1) {
	    i__1 = *leny;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.;
		ig = locg[j];
		i__2 = locj[j + 1] - 1;
		for (ij = locj[j]; ij <= i__2; ++ij) {
		    ir = indj[ij];
		    if (ir > *lenx) {
			goto L400;
		    }
		    sum += gcon[ig] * x[ir];
		    ++ig;
		}
L400:
		y[j] += *alpha * sum;
	    }
	}
/* task .eq. Normal */
    }
/* general alpha */
    return 0;
} /* s8gprd_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gprd */
/* Subroutine */ int s8gsiz_(integer *m, integer *nncon, integer *nnjac, 
	integer *ne, integer *nlocj, integer *locj, integer *indj, integer *
	negcon)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, ir, neg, last, nlocg;

/* ================================================================= */
/* s8Gsiz  counts the number of nonlinear Jacobian elements. */

/* 04 Nov 2000: First version of s8Gsiz */
/* 31 Aug 2008: Local variable used for negCon. */
/* ================================================================= */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --indj;
    --locj;

    /* Function Body */
    neg = 0;
    nlocg = *nnjac + 1;
    if (*nncon > 0) {
	last = locj[nlocg] - 1;
	if (*nncon == *m) {
	    neg = last;
	} else {
	    i__1 = last;
	    for (k = 1; k <= i__1; ++k) {
		ir = indj[k];
		if (ir <= *nncon) {
		    ++neg;
		}
	    }
	}
    }
    *negcon = max(1,neg);
    return 0;
} /* s8gsiz_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Gsiz */
/* Subroutine */ int s8infs_(logical *elastc, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, doublereal *tolx, doublereal *wtinf, 
	doublereal *prinf, doublereal *duinf, integer *jprinf, integer *
	jduinf, doublereal *bl, doublereal *bu, doublereal *fx, doublereal *
	rc, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal v, w, dj, tol, viol, slack;

/*     ================================================================== */
/*     s8Infs computes the maximum primal and dual infeasibilities, */
/*     using bl, bu, rc, x and the true nonlinear slacks Fxslk. */
/*     The linear constraints and bounds are assumed to be satisfied. */
/*     The primal infeasibility is therefore the maximum violation of */
/*     the nonlinear constraints. */
/*     The dual infeasibility is the maximum complementarity gap */
/*     for the bound constraints (with bounds assumed to be no further */
/*     than 1.0 from each x(j)). */

/*     prInf, duInf   return the max primal and dual infeas. */

/*     20 Feb 1994: First version based on Minos 5.5 routine m8infs. */
/*     25 Oct 1996: Elastic mode added. */
/*     29 Apr 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --rc;
    --bu;
    --bl;
    --fx;

    /* Function Body */
    *jprinf = 0;
    *prinf = 0.;
    tol = *tolx;
/*     See how much  Fx  violates the bounds on the nonlinear slacks. */
/*     prInf is the maximum violation. */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	slack = fx[i__];
/* Computing MAX */
	d__1 = 0., d__2 = bl[j] - slack, d__1 = max(d__1,d__2), d__2 = slack 
		- bu[j];
	viol = max(d__1,d__2);
	if (*prinf < viol) {
	    *prinf = viol;
	    *jprinf = j;
	}
    }
/*     ------------------------------------------------------------------ */
/*     + rc(j)  is the multiplier for lower bound constraints. */
/*     - rc(j)  is the multiplier for upper bound constraints. */
/*     duInf is the maximum component-wise complementarity measure. */
/*     ------------------------------------------------------------------ */
    *jduinf = 0;
    *duinf = 0.;
    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {
	dj = rc[j];
	if (dj != 0.) {
	    if (dj > 0.) {
/* Computing MIN */
		d__1 = x[j] - bl[j];
		dj *= min(d__1,1.);
	    } else if (dj < 0.) {
/* Computing MIN */
		d__1 = bu[j] - x[j];
		dj = -dj * min(d__1,1.);
	    }
	    if (*duinf < dj) {
		*duinf = dj;
		*jduinf = j;
	    }
	}
/* dj nonzero */
    }
/*     ------------------------------------------------------------------ */
/*     Include contributions from the elastic variables. */
/*     ------------------------------------------------------------------ */
    if (*elastc) {
	i__1 = *n + *nncon;
	for (j = *n + 1; j <= i__1; ++j) {
	    dj = rc[j];
	    v = bl[j] - x[j];
	    w = x[j] - bu[j];
	    if (v > tol) {
		dj = (d__1 = *wtinf - dj, abs(d__1)) * min(v,1.);
	    } else if (w > tol) {
		dj = (d__1 = *wtinf + dj, abs(d__1)) * min(w,1.);
	    } else {
		dj = 0.;
	    }
	    if (*duinf < dj) {
		*duinf = dj;
		*jduinf = j;
	    }
	}
    }
    return 0;
} /* s8infs_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8Infs */
/* Subroutine */ int s8iqn_(integer *iexit, integer *info, U_fp hprod, U_fp 
	hprod1, U_fp mnrlog, integer *hqntype, integer *hvcalls, logical *
	elastc, logical *gotr, integer *itn, integer *itqp, integer *lenr, 
	integer *m, integer *maxs, integer *mbs, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, integer *nnobj0, integer *nnobj, 
	integer *nnl0, integer *nnl, integer *ns, integer *ndegen, integer *
	mjrprt, integer *mnrprt, integer *minimz, integer *iobj, doublereal *
	condhz, doublereal *sclobj, doublereal *objadd, doublereal *objqp, 
	doublereal *tolfp, doublereal *tolqpk, doublereal *tolx, integer *
	ninf, doublereal *sinf, doublereal *wtinf, doublereal *u0ii, 
	doublereal *pinorm, integer *ne, integer *nlocj, integer *locj, 
	integer *indj, doublereal *jcol, integer *hetype, integer *hestat, 
	integer *hfeas, integer *hs, integer *kbs, doublereal *ascale, 
	doublereal *bl, doublereal *bu, doublereal *blbs, doublereal *bubs, 
	doublereal *gbs, doublereal *gqp, doublereal *gobj, doublereal *hdx, 
	doublereal *pbs, doublereal *pi, doublereal *r__, doublereal *rc, 
	doublereal *rg, doublereal *rg2, doublereal *qprhs, doublereal *x, 
	doublereal *xbs, doublereal *xqp0, doublereal *xqp, integer *iy, 
	integer *iy1, doublereal *y, doublereal *y1, doublereal *y2, char *cu,
	 integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer 
	*lenru, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1500[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP feasibility\002,\002 phase\002)";
    static char fmt_1800[] = "(\002 Itn\002,i7,\002: Ill-conditioned CG null"
	    "-space basis.\002,\002 Cond = \002,1p,e8.1)";
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Infeasible subproblem"
	    ".\002,\002 Elastic mode started with weight = \002,1p,e8.1)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Feasible QP non-elast"
	    "ics\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Feasible QP subproblem"
	    " \002)";
    static char fmt_1400[] = "(\002 Itn\002,i7,\002: Large multipliers.\002"
	    ",\002 Elastic mode started with weight = \002,1p,e8.1)";
    static char fmt_1510[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP optimality\002,\002 phase\002)";
    static char fmt_1600[] = "(\002 Itn\002,i7,\002: Unbounded QP subproble"
	    "m\002)";
    static char fmt_1950[] = "(\002 Itn\002,i7,\002: Hessian reset\002)";
    static char fmt_1900[] = "(\002 Itn\002,i7,\002: Too many subspace itera"
	    "tions\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer itqptargt;
    extern /* Subroutine */ int s8h0_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static integer nnh;
    static char str[80];
    static integer nnh0;
    extern /* Subroutine */ int s5qn_(integer *, integer *, char *, logical *,
	     integer *, U_fp, U_fp, U_fp, logical *, logical *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical done, newb;
    static integer ngqp;
    static logical luok;
    static integer ngqp0;
    static doublereal objfp;
    static logical needx;
    static doublereal flmax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lureq;
    static logical newlu;
    static integer nswap, zngqp;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static logical feasbl;
    static doublereal zcndbd;
    static integer lemode;
    static logical needlu;
    extern doublereal dnormj_(integer *, doublereal *, integer *);
    static integer inform__, lvlinf;
    static logical solved;
    static integer itnlim, mminor;
    static logical normin;
    static doublereal pinnln, plinfy, targtz;
    static integer subopt, typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2trylu_(integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, integer *, doublereal *, integer 
	    *);
    static char probtag[20];
    static integer itqpmax;

    /* Fortran I/O blocks */
    static icilist io___141 = { 0, str, 0, fmt_1500, 80, 1 };
    static icilist io___142 = { 0, str, 0, fmt_1800, 80, 1 };
    static icilist io___143 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___144 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___145 = { 0, str, 0, fmt_1300, 80, 1 };
    static icilist io___148 = { 0, str, 0, fmt_1400, 80, 1 };
    static icilist io___149 = { 0, str, 0, fmt_1510, 80, 1 };
    static icilist io___150 = { 0, str, 0, fmt_1600, 80, 1 };
    static icilist io___151 = { 0, str, 0, fmt_1950, 80, 1 };
    static icilist io___152 = { 0, str, 0, fmt_1800, 80, 1 };
    static icilist io___153 = { 0, str, 0, fmt_1900, 80, 1 };
    static icilist io___154 = { 0, str, 0, fmt_1950, 80, 1 };


/*     ================================================================== */
/*     s8iQN   computes  xQP, the solution of the QP subproblem. */
/*     By construction, the problem has  nnL  nonlinear variables, */

/*     The SQP base point  x  is not altered. */

/*     On entry, the LU factorization is assumed to be known. */
/*     The arrays  xBS, blBS and buBS are defined. */

/*     iExit     Status */
/*     -----     ------ */
/*      >0         Fatal error */
/*       0         QP solution found */
/*      -1         Too many iterations */
/*      -2         Too many superbasics */

/*     LUreq =  1  Frequency */
/*     LUreq =  2  LU nonzeros increased */
/*     LUreq =  3 */
/*     LUreq =  4 */
/*     LUreq =  5  Singular after LU mod */
/*     LUreq =  6  Unstable LU mod (growth in new column of U) */
/*     LUreq =  7  Not enough memory */
/*     LUreq =  8 */
/*     LUreq =  9 */
/*     LUreq = 10  Row error in setx */
/*     LUreq = 11  Big  dx   in setx */

/*     LUreq = 20 */
/*     LUreq = 21  Iterative refinement failed in QP */
/*     LUreq = 22  Unbounded QP */
/*     LUreq = 23 */
/*     LUreq = 24  Small directional derivative in QP */
/*     LUreq = 25  Ill-conditioned Z */
/*     LUreq = 26  Indefinite Z'HZ in QP */
/*     LUreq = 27  R singular after bound swap in QP */

/*     On output, */
/*     QPerr points to ' ', 't', 'u' or 'w'. */
/*     QPfea points to ' '  or 'i'. */

/*     15 Jun 2001: First version of s8iQN based on s8iQP. */
/*     31 Jul 2003: dnormj used for norm of the nonlinear pis. */
/*     03 Aug 2003: snEXIT and snPRNT adopted. */
/*     26 Dec 2003: calls to s2tryLU added. */
/*     19 Jun 2008: Hprod, Hprod1 added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* condition estimate of Z */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --info;
    --r__;
    --pi;
    --rg2;
    --rg;
    --xbs;
    --pbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --y2;
    --y1;
    --y;
    --iy1;
    --iy;
    --xqp;
    --xqp0;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --qprhs;
    --gobj;
    --hdx;
    --gqp;
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
    itnlim = iw[89];
/* limit on total iterations */
    mminor = iw[91];
/* limit on minor iterations */
    flmax = rw[8];
/* est. of the largest pos. real */
    zcndbd = rw[86];
/* bound on the condition of Z */
    plinfy = flmax;
    *iexit = 0;
    targtz = zcndbd;
    *condhz = 0.;
    feasbl = FALSE_;
/* Status of the non-elastic variables */
    normin = ! (*elastc);
    s_copy(probtag, "QP subproblem", (ftnlen)20, (ftnlen)13);
    *itqp = 0;
    iw[220] = 0;
    iw[221] = 0;
    iw[223] = 1;
    iw[225] = 1;
    info[5] = 0;
    info[4] = 0;
    nnh = *nnl;
    nnh0 = *nnl0;
    ngqp = *nnl;
    ngqp0 = max(ngqp,1);
/*     Set lEmode to switch to Elastic mode on infeasibility. */
/*     When in elastic mode, set lvlInf to use the composite objective: */
/*     w1*Obj + w2*sInf,  with W1 = 0, W2 = wtInf. This minimizes the */
/*     sum of the infeasibilities of the elastic constraints subject to */
/*     the nonelastic constraints. */
    lemode = 1;
    lvlinf = 2;
    typelu = 3;
    lureq = 0;
/*     ================================================================== */
/*     Find a feasible point for this linearization. */
/*     If the constraints are linear, x is already feasible. */
/*     ================================================================== */
    if (*nncon > 0) {
/*        --------------------------------------------------------------- */
/*        Find a feasible point. */
/*        If the reduced Hessian is defined, then it is updated. */
/*        If the constraints are infeasible, minimize the sum of the */
/*        elastic variables, subject to keeping the non-elastic variables */
/*        feasible.  Elastic variables can move outside their bounds. */
/*        --------------------------------------------------------------- */
	zngqp = 0;
/* No objective term */
	itqpmax = itnlim;
	itqptargt = itnlim;
	subopt = -1;
	luok = TRUE_;
	done = FALSE_;
/*        =============================================================== */
/* +       while (.not. done  .and.  LUok) do */
L500:
	if (! done && luok) {
	    needlu = lureq > 0;
	    if (needlu) {
		s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, 
			mjrprt, &lureq, m, mbs, n, nb, &nnh, ns, &nswap, ne, 
			nlocj, &locj[1], &indj[1], &jcol[1], &kbs[1], &hs[1], 
			&bl[1], &bu[1], &blbs[1], &bubs[1], nncon0, nncon, &
			qprhs[1], &xqp[1], &xbs[1], &iy[1], &iy1[1], &y[1], &
			y1[1], &iw[1], leniw, &rw[1], lenrw);
		if (*iexit != 0) {
		    goto L900;
		}
		if (nswap > 0) {
		    *gotr = FALSE_;
		}
		lureq = 0;
	    }
	    needx = needlu;
	    s5qn_(&inform__, &c__4, probtag, elastc, &subopt, (U_fp)hprod, (
		    U_fp)hprod1, (U_fp)mnrlog, gotr, &needlu, &typelu, &needx,
		     lenr, m, maxs, mbs, n, nb, ndegen, hvcalls, &ngqp0, &
		    zngqp, nnobj0, nnobj, &nnh0, &nnh, ns, itqp, &itqpmax, &
		    itqptargt, itn, &lemode, &lvlinf, mnrprt, minimz, iobj, 
		    sclobj, objadd, &objfp, condhz, &targtz, tolfp, tolqpk, 
		    tolx, ninf, sinf, wtinf, pinorm, ne, nlocj, &locj[1], &
		    indj[1], &jcol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[
		    1], &kbs[1], &ascale[1], &bl[1], &bu[1], &blbs[1], &bubs[
		    1], &gbs[1], &gobj[1], &gqp[1], &hdx[1], &pbs[1], &pi[1], 
		    &r__[1], &rc[1], &rg[1], &rg2[1], nncon0, nncon, &qprhs[1]
		    , nnl0, nnl, &x[1], &xqp[1], &xbs[1], &x[1], &iy[1], &iy1[
		    1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], leniu, &
		    ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
		     (ftnlen)20, (ftnlen)8, (ftnlen)8);
/*           Check for trouble.  Here are the possibilities: */
/*           inform      Result */
/*           ------      ------ */
/*            >0         Fatal LU error */
/*             0         Found a feasible point for the nonelastics */
/*            -1         The nonelastics are infeasible */
/*            -2         Phase 1 is unbounded */
/*            -3         Too many iterations */
/*            -5         Superbasic limit exceeded */
/*            -6         Void */
/*            -7         Void */
/*            -8         Ill-conditioned Z */
/*            -9         Too many subspace iterations (should not happen) */
	    if (inform__ > 0) {
		*iexit = inform__;
/* Fatal LU error */
		goto L900;
	    }
	    done = inform__ == 0 || inform__ == -3 || inform__ == -5;
	    if (! done) {
/*              ========================================================= */
/*              Trouble. */
/*              inform = -2 means that the phase 1 was unbounded, which */
/*                          can only occur if a bad basis gives a large */
/*                          search direction */
/*              inform = -1 means that the nonelastics are infeasible, */
/*                          which should not happen since we already */
/*                          know a feasible point for the nonelastics. */
/*              inform = -8 means that  R  is being updated but a crude */
/*                          estimate of condZ is bigger than targtZ. */
/*                          Refactorize B, possibly with a reduced factor */
/*                          tol. If the factor tol is already tight, */
/*                          accept Z, however bad. */
/*              ========================================================= */
		if (inform__ == -1 || inform__ == -2) {
/*                 Treat both cases as infeasible. Repeatedly refactorize */
/*                 with tighter tols before declaring LC infeasibility. */
		    inform__ = -1;
		    s_wsfi(&io___141);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		} else if (inform__ == -8) {
		    s_wsfi(&io___142);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&rw[192], (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    s2trylu_(itn, &c__25, ns, &lureq, &luok, &typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		    if (! luok) {
			targtz = plinfy;
			luok = TRUE_;
		    }
		}
	    }
	    goto L500;
	}
/* +       end while */
/*        --------------------------------------------------------------- */
	if (inform__ < 0) {
	    goto L800;
	}
/* Itns or infeasible linear constr */
	if (*elastc && normin) {
/*           ------------------------------------------------------------ */
/*           The QP switched to elastic mode. */
/*           The linearized constraints are infeasible. */
/*           ------------------------------------------------------------ */
	    if (*mjrprt >= 1 || *mnrprt >= 10) {
		s_wsfi(&io___143);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    }
	    *gotr = FALSE_;
	    *hqntype = 2;
	    s8h0_(hqntype, &nnh, u0ii, &iw[1], leniw, &rw[1], lenrw);
	} else if (*mjrprt > 10 && *mnrprt > 10) {
/*           No change in mode. */
	    if (*elastc) {
		s_wsfi(&io___144);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    } else {
		s_wsfi(&io___145);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     The inelastic variables (x's and linear slacks) are now feasible. */
/*     Save them in xQP0 for use with the BFGS update. */

/*     Solve the QP subproblem. */
/*     Loop back sometimes if we need a BS factorize. */
/*     ------------------------------------------------------------------ */
/* nlnCon */
    dcopy_(nb, &xqp[1], &c__1, &xqp0[1], &c__1);
    feasbl = TRUE_;
/* the nonelastics are feasible */
    itqptargt = *itqp + mminor;
    itqpmax = itnlim;
    lureq = 0;
    typelu = 3;
    luok = TRUE_;
    done = FALSE_;
    solved = FALSE_;
/*     ================================================================== */
/* +    while (.not. (solved  .or.  done)  .and.  LUok) do */
L600:
    if (! (solved || done) && luok) {
/*        --------------------------------------------------------------- */
/*        Refactorize the basis if necessary. */
/*        --------------------------------------------------------------- */
	needlu = lureq > 0;
	if (needlu) {
	    s2bfac_(iexit, &typelu, &needlu, &newlu, &newb, iobj, itn, mjrprt,
		     &lureq, m, mbs, n, nb, &nnh, ns, &nswap, ne, nlocj, &
		    locj[1], &indj[1], &jcol[1], &kbs[1], &hs[1], &bl[1], &bu[
		    1], &blbs[1], &bubs[1], nncon0, nncon, &qprhs[1], &xqp[1],
		     &xbs[1], &iy[1], &iy1[1], &y[1], &y1[1], &iw[1], leniw, &
		    rw[1], lenrw);
	    if (*iexit != 0) {
		goto L900;
	    }
	    if (nswap > 0) {
		*gotr = FALSE_;
	    }
	    lureq = 0;
	}
/*        How do we update R if the superbasics change? */
/*        --------------------------------------------------------------- */
/*        Solve the QP subproblem using a quasi-Newton method. */
/*        --------------------------------------------------------------- */
	if (*mnrprt >= 10) {
	    iw[223] = 1;
/* QN print   header */
	    iw[225] = 1;
/* QN summary header */
	}
/*        Set lEmode to switch to Elastic mode on infeasibility. */
/*        Set lvlInf to use the composite objective  Obj + wtInf*sInf */
/*        after any switch to elastic mode. */
	lvlinf = 1;
	if (*nnl > 0) {
	    subopt = 0;
	} else {
	    subopt = -1;
	}
	needx = needlu;
	s5qn_(&inform__, &c__5, probtag, elastc, &subopt, (U_fp)hprod, (U_fp)
		hprod1, (U_fp)mnrlog, gotr, &needlu, &typelu, &needx, lenr, m,
		 maxs, mbs, n, nb, ndegen, hvcalls, &ngqp0, &ngqp, nnobj0, 
		nnobj, &nnh0, &nnh, ns, itqp, &itqpmax, &itqptargt, itn, &
		lemode, &lvlinf, mnrprt, minimz, iobj, sclobj, objadd, objqp, 
		condhz, &targtz, tolfp, tolqpk, tolx, ninf, sinf, wtinf, 
		pinorm, ne, nlocj, &locj[1], &indj[1], &jcol[1], &hetype[1], &
		hestat[1], &hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1], &
		bu[1], &blbs[1], &bubs[1], &gbs[1], &gobj[1], &gqp[1], &hdx[1]
		, &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], &rg2[1], nncon0, 
		nncon, &qprhs[1], nnl0, nnl, &x[1], &xqp[1], &xbs[1], &x[1], &
		iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], 
		leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/*        iExit       Result */
/*        -----       ------ */
/*         >0         Fatal LU error */
/*          0         QP solution found */
/*         -1         The nonelastics are infeasible */
/*         -2         The QP subproblem is unbounded */
/*         -3         Too many iterations */
/*         -4         Void */
/*         -5         Too many superbasics */
/*         -6         Void */
/*         -7         Void */
/*         -8         Ill-conditioned Z */
/*         -9         too many subspace iterations */
	if (inform__ > 0) {
	    *iexit = inform__;
	    goto L900;
	}
	solved = inform__ == 0 || inform__ == -9;
	done = inform__ == -3 || inform__ == -5;
	if (done) {
/*           ============================================================ */
/*           Relax */
/*           ============================================================ */
	} else if (solved) {
/*           ============================================================ */
/*           Finish if there are no large nonlinear pi's. */
/*           Otherwise, re-solve the QP in elastic mode */
/*           ============================================================ */
	    if (! (*elastc)) {
		pinnln = dnormj_(nncon, &pi[1], &c__1);
		if (pinnln > *wtinf) {
		    *elastc = TRUE_;
		    solved = FALSE_;
		    s_wsfi(&io___148);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(
			    doublereal));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	} else {
/*           ============================================================ */
/*           Trouble. */
/*           ============================================================ */
	    if (inform__ == -1) {
/*              --------------------------------------------------------- */
/*              The nonelastics are infeasible. This should not happen. */
/*              Phase 1 has already found a feasible point for the */
/*              nonelastics, so the basis must be ill-conditioned. */
/*              Refactorize with tighter tols and restart at the known */
/*              feasible point.  Reduce the feasibility tol to try and */
/*              prevent repeats. */
/*              --------------------------------------------------------- */
		s_wsfi(&io___149);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		*elastc = FALSE_;
		dcopy_(nb, &xqp0[1], &c__1, &xqp[1], &c__1);
	    } else if (inform__ == -2) {
/*              --------------------------------------------------------- */
/*              The QP is unbounded. */
/*              As the Hessian is positive definite, this is probably */
/*              because of an ill-conditioned reduced Hessian. */
/*              Reset both the full and reduced Hessian. */
/*              --------------------------------------------------------- */
		s_wsfi(&io___150);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		if (*hqntype != 2) {
		    *gotr = FALSE_;
		    *hqntype = 2;
		    s8h0_(hqntype, &nnh, u0ii, &iw[1], leniw, &rw[1], lenrw);
		    s_wsfi(&io___151);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
		s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
	    } else if (inform__ == -8) {
/*              --------------------------------------------------------- */
/*              condZ > targtZ  while computing the search direction. */
/*              Refactorize B, possibly with a reduced factor tol. If */
/*              the factor tol is already tight, accept Z, however bad. */
/*              --------------------------------------------------------- */
		s_wsfi(&io___152);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rw[192], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__25, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    targtz = plinfy;
		    luok = TRUE_;
		}
	    } else if (inform__ == -9) {
/*              --------------------------------------------------------- */
/*              Too many CG subspace iterations. */
/*              --------------------------------------------------------- */
		s_wsfi(&io___153);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		if (*hqntype != 2) {
		    *gotr = FALSE_;
		    *hqntype = 2;
		    s8h0_(hqntype, &nnh, u0ii, &iw[1], leniw, &rw[1], lenrw);
		    s_wsfi(&io___154);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	}
	goto L600;
    }
/* +    end while */
/*     ------------------------------------------------------------------ */
L800:
    if (*ninf > 0) {
	info[4] = 1;
    }
    if (inform__ == 0) {
	info[5] = max(subopt,0);
    } else if (inform__ == -1) {
	*iexit = 15;
/* infeasible nonelastics */
    } else if (inform__ == -2) {
	info[5] = 3;
/* unbounded subproblem */
    } else if (inform__ == -3) {
	info[5] = 2;
	*iexit = -1;
/* too many iterations */
    } else if (inform__ == -5 && feasbl) {
	info[5] = 5;
	*iexit = -2;
/* superbasic limit */
    } else if (inform__ == -5) {
	*iexit = 33;
/* superbasic limit */
    } else if (inform__ == -9) {
/*        Relax and hope for the best */
    } else {
	*iexit = 44;
/* ill-conditioned null-space basis */
    }
L900:
    return 0;
} /* s8iqn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8iQN */
/* Subroutine */ int s8iqp_(integer *iexit, integer *info, U_fp hprod, U_fp 
	hprod1, U_fp mnrlog, integer *hqntype, integer *hvcalls, logical *
	elastc, logical *gotr, integer *itn, integer *itqp, integer *lenr, 
	integer *m, integer *maxs, integer *mbs, integer *n, integer *nb, 
	integer *nncon0, integer *nncon, integer *nnobj0, integer *nnobj, 
	integer *nnl0, integer *nnl, integer *ns, integer *ndegen, integer *
	mjrprt, integer *mnrprt, integer *minimz, integer *iobj, doublereal *
	sclobj, doublereal *objadd, doublereal *objqp, doublereal *tolfp, 
	doublereal *tolqpk, doublereal *tolx, integer *ninf, doublereal *sinf,
	 doublereal *wtinf, doublereal *u0ii, doublereal *pinorm, integer *ne,
	 integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	integer *hetype, integer *hestat, integer *hfeas, integer *hs, 
	integer *kbs, doublereal *ascale, doublereal *bl, doublereal *bu, 
	doublereal *blbs, doublereal *bubs, doublereal *gbs, doublereal *gqp, 
	doublereal *gobj, doublereal *hdx, doublereal *pbs, doublereal *pi, 
	doublereal *r__, doublereal *rc, doublereal *rg, doublereal *qprhs, 
	doublereal *x, doublereal *xbs, doublereal *xqp0, doublereal *xqp, 
	integer *iy, integer *iy1, doublereal *y, doublereal *y1, doublereal *
	y2, char *cu, integer *lencu, integer *iu, integer *leniu, doublereal 
	*ru, integer *lenru, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1500[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP feasibility\002,\002 phase\002)";
    static char fmt_1100[] = "(\002 Itn\002,i7,\002: Infeasible subproblem"
	    ".\002,\002 Elastic mode started with weight = \002,1p,e8.1)";
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Feasible QP non-elast"
	    "ics\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Feasible QP subproblem"
	    " \002)";
    static char fmt_1400[] = "(\002 Itn\002,i7,\002: Large multipliers.\002"
	    ",\002 Elastic mode started with weight = \002,1p,e8.1)";
    static char fmt_1510[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics "
	    "in QP optimality\002,\002 phase\002)";
    static char fmt_1600[] = "(\002 Itn\002,i7,\002: Indefinite QP reduced H"
	    "essian\002)";
    static char fmt_1700[] = "(\002 Itn\002,i7,\002: Large QP reduced gradie"
	    "nt\002)";
    static char fmt_1900[] = "(\002 Itn\002,i7,\002: Reduced Hessian rese"
	    "t\002)";
    static char fmt_1800[] = "(\002 Itn\002,i7,\002: Ill-conditioned QP null"
	    "-space basis.\002,\002 Cond = \002,1p,e8.1)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer itqptargt, nnh;
    static doublereal eps;
    extern /* Subroutine */ int s8h0_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static char str[80];
    static integer nnh0;
    extern /* Subroutine */ int s5qp_(integer *, integer *, char *, logical *,
	     integer *, U_fp, U_fp, U_fp, logical *, logical *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static logical done;
    static integer ngqp;
    static logical luok;
    static integer znnh, ngqp0;
    static doublereal objfp;
    static logical needx;
    static doublereal flmax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lureq, zngqp;
    extern /* Subroutine */ int s8getr_(integer *, U_fp, U_fp, integer *, 
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal hcndbd;
    static logical feasbl;
    static doublereal zcndbd;
    static integer lemode;
    static logical needlu;
    extern doublereal dnormj_(integer *, doublereal *, integer *);
    static integer inform__, lvlinf;
    static logical solved;
    static integer itnlim, mminor;
    static logical normin;
    static doublereal pinnln, plinfy, rgnorm, targth;
    static integer subopt;
    static doublereal targtz;
    static integer typelu;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2trylu_(integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, integer *, doublereal *, integer 
	    *);
    static char probtag[20];
    static integer itqpmax;

    /* Fortran I/O blocks */
    static icilist io___188 = { 0, str, 0, fmt_1500, 80, 1 };
    static icilist io___189 = { 0, str, 0, fmt_1100, 80, 1 };
    static icilist io___190 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___191 = { 0, str, 0, fmt_1300, 80, 1 };
    static icilist io___194 = { 0, str, 0, fmt_1400, 80, 1 };
    static icilist io___195 = { 0, str, 0, fmt_1510, 80, 1 };
    static icilist io___196 = { 0, str, 0, fmt_1600, 80, 1 };
    static icilist io___197 = { 0, str, 0, fmt_1700, 80, 1 };
    static icilist io___198 = { 0, str, 0, fmt_1900, 80, 1 };
    static icilist io___199 = { 0, str, 0, fmt_1800, 80, 1 };


/*     ================================================================== */
/*     s8iQP  computes  xQP, the solution of the QP subproblem. */
/*     By construction, the problem has  nnL  nonlinear variables, */

/*     The SQP base point  x  is not altered. */

/*     On entry, the LU factorization is assumed to be known. */
/*     The arrays  xBS, blBS and buBS are defined. */

/*     iExit     Status */
/*     -----     ------ */
/*      >0         Fatal error */
/*       0         QP solution found */
/*      -1         Too many iterations */
/*      -2         Too many superbasics */

/*     LUreq =  1  Frequency */
/*     LUreq =  2  LU nonzeros increased */
/*     LUreq =  3 */
/*     LUreq =  4 */
/*     LUreq =  5  Singular after LU mod */
/*     LUreq =  6  Unstable LU mod (growth in new column of U) */
/*     LUreq =  7  Not enough memory */
/*     LUreq =  8 */
/*     LUreq =  9 */
/*     LUreq = 10  Row error in setx */
/*     LUreq = 11  Big  dx   in setx */

/*     LUreq = 20 */
/*     LUreq = 21  Iterative refinement failed in QP */
/*     LUreq = 22  Unbounded QP */
/*     LUreq = 23 */
/*     LUreq = 24  Small directional derivative in QP */
/*     LUreq = 25  Ill-conditioned Z */
/*     LUreq = 26  Indefinite Z'HZ in QP */
/*     LUreq = 27  R singular after bound swap in QP */

/*     On output, */
/*     QPerr points to ' ', 't', 'u' or 'w'. */
/*     QPfea points to ' '  or 'i'. */

/*     30 Dec 1991: First version of s8iQP. */
/*     19 Jul 1997: Thread-safe version. */
/*     31 Jul 2003: dnormj used for norm of the nonlinear pis. */
/*     03 Aug 2003: snEXIT and snPRNT  adopted. */
/*     19 Jun 2008: Hprod, Hprod1 added as arguments. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* condition estimate of Z */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --info;
    --r__;
    --pi;
    --rg;
    --xbs;
    --pbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --y2;
    --y1;
    --y;
    --iy1;
    --iy;
    --xqp;
    --xqp0;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --qprhs;
    --gobj;
    --hdx;
    --gqp;
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
    itnlim = iw[89];
/* limit on total iterations */
    mminor = iw[91];
/* limit on minor iterations */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    flmax = rw[8];
/* est. of the largest pos. real */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    zcndbd = rw[86];
/* bound on the condition of Z */
    plinfy = flmax;
    *iexit = 0;
    *itqp = 0;
    targtz = zcndbd;
    targth = hcndbd;
    feasbl = FALSE_;
/* Status of the non-elastic variables */
    normin = ! (*elastc);
    s_copy(probtag, "QP subproblem", (ftnlen)20, (ftnlen)13);
    iw[220] = 0;
    iw[221] = 0;
    iw[223] = 1;
    iw[225] = 1;
    info[5] = 0;
    info[4] = 0;
    nnh = *nnl;
    nnh0 = *nnl0;
    ngqp = *nnl;
    ngqp0 = max(ngqp,1);
/* Set lEmode to switch to Elastic mode on infeasibility. */
/* When in elastic mode, set lvlInf to use the composite objective: */
/* w1*Obj + w2*sInf,  with w1 = 0, w2 = wtInf. This minimizes the */
/* sum of the infeasibilities of the elastic constraints subject to */
/* the nonelastic constraints. */
    lemode = 1;
    lvlinf = 2;
    typelu = 3;
    lureq = 0;
/* ----------------------------------------------------------------- */
/* Find a feasible point. */
/* If the constraints are linear, x is already feasible. */
/* Otherwise, find a feasible x for this linearization. */
/* Minimize the sum of the elastic variables */
/* subject to keeping the non-elastic variables feasible. */
/* Elastic variables can move outside their bounds. */
/* ----------------------------------------------------------------- */
    if (*nncon > 0) {
	zngqp = 0;
/* No objective term in phase 1 */
	znnh = 0;
/* No Hessian either */
	itqpmax = itnlim;
	itqptargt = itnlim;
	subopt = -1;
	*gotr = FALSE_;
	luok = TRUE_;
	done = FALSE_;
/* while (.not. done  .and.  LUok) do */
L500:
	if (! done && luok) {
	    needlu = lureq > 0;
	    needx = needlu;
	    s5qp_(&inform__, &c__4, probtag, elastc, &subopt, (U_fp)hprod, (
		    U_fp)hprod1, (U_fp)mnrlog, gotr, &needlu, &typelu, &needx,
		     lenr, m, maxs, mbs, n, nb, ndegen, hvcalls, &ngqp0, &
		    zngqp, nnobj0, nnobj, &nnh0, &znnh, ns, itqp, &itqpmax, &
		    itqptargt, itn, &lemode, &lvlinf, mnrprt, minimz, iobj, 
		    sclobj, objadd, &objfp, &targth, &targtz, tolfp, tolqpk, 
		    tolx, ninf, sinf, wtinf, pinorm, &rgnorm, ne, nlocj, &
		    locj[1], &indj[1], &jcol[1], &hetype[1], &hestat[1], &
		    hfeas[1], &hs[1], &kbs[1], &ascale[1], &bl[1], &bu[1], &
		    blbs[1], &bubs[1], &gbs[1], &gobj[1], &gqp[1], &hdx[1], &
		    pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], nncon0, nncon, &
		    qprhs[1], nnl0, nnl, &x[1], &xqp[1], &xbs[1], &x[1], &iy[
		    1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1],
		     leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[
		    1], lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* Check for trouble.  Here are the possibilities: */

/* inform      Result */
/* ------      ------ */
/*  >0         Fatal LU error */
/*   0         Found a feasible point for the nonelastics */
/*  -1         The nonelastics are infeasible */
/*  -2         Phase 1 is unbounded */
/*  -3         Too many iterations */
/*  -4         Weak minimizer (after starting elastic mode) */
/*  -5         Superbasic limit exceeded */
	    if (inform__ > 0) {
		*iexit = inform__;
/* Fatal LU error */
		goto L900;
	    }
	    done = inform__ == 0 || inform__ == -3 || inform__ == -4 || 
		    inform__ == -5;
	    if (! done) {
/* ======================================================== */
/* Trouble. */
/* inform = -2 implies that the phase 1 was unbounded, */
/*             which can only occur if a bad basis gives */
/*             a large search direction */
/* inform = -1 means that the nonelastics are infeasible, */
/*             which should not happen since we already */
/*             know a feasible point for the nonelastics. */
/* ======================================================== */
/* Treat both cases as infeasible. Repeatedly refactorize */
/* with tighter tols before declaring LC infeasibility. */
		inform__ = -1;
		s_wsfi(&io___188);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
	    }
	    goto L500;
	}
/* end while */
	if (inform__ < 0) {
	    goto L800;
	}
/* Itns, infeas inelastcs, or nZlim */
	if (*elastc && normin) {
/* The QP switched to elastic mode. */
/* The linearized constraints are infeasible. */
	    if (*mjrprt >= 1 || *mnrprt >= 10) {
		s_wsfi(&io___189);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    }
	    *gotr = FALSE_;
	    *hqntype = 2;
	    s8h0_(hqntype, &nnh, u0ii, &iw[1], leniw, &rw[1], lenrw);
	} else if (*mjrprt > 10 && *mnrprt > 10) {
/* No change in mode. */
	    if (*elastc) {
		s_wsfi(&io___190);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    } else {
		s_wsfi(&io___191);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)80);
	    }
	}
    }
/* The inelastic variables (x and linear slacks) are now feasible. */
/* Save them in xQP0 for use with the BFGS update. */
/* nlnCon */
    dcopy_(nb, &xqp[1], &c__1, &xqp0[1], &c__1);
/* ----------------------------------------------------------------- */
/* Solve the QP subproblem. */
/* Loop back if we need better LU factors. */
/* ----------------------------------------------------------------- */
    feasbl = TRUE_;
/* the nonelastics are feasible */
    itqpmax = itnlim;
    itqptargt = *itqp + mminor;
    lureq = 0;
    typelu = 3;
    luok = TRUE_;
    done = FALSE_;
    solved = FALSE_;
/* while (.not. (solved  .or.  done)  .and.  LUok) do */
L600:
    if (! (solved || done) && luok) {
	if (! (*gotr) || lureq > 0) {
/* Compute and factorize the initial Z'HZ. */
/* The basis is refactorized if necessary. */
	    s8getr_(iexit, (U_fp)hprod, (U_fp)hprod1, hqntype, hvcalls, gotr, 
		    &typelu, &lureq, itn, lenr, m, mbs, n, nb, nncon0, nncon, 
		    &nnh, ns, mjrprt, minimz, iobj, u0ii, &targth, &targtz, 
		    ne, nlocj, &locj[1], &indj[1], &jcol[1], &hs[1], &kbs[1], 
		    &bl[1], &bu[1], &blbs[1], &bubs[1], &r__[1], &qprhs[1], &
		    xqp[1], &xbs[1], &iy[1], &iy1[1], &y[1], &y1[1], &y2[1], 
		    cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	    if (*iexit != 0) {
		goto L900;
	    }
	    lureq = 0;
	}
/* -------------------------------------------------------------- */
/* Solve the QP subproblem. */
/* -------------------------------------------------------------- */
	if (*mnrprt >= 10) {
	    iw[223] = 1;
/* QP print   header */
	    iw[225] = 1;
/* QP summary header */
	}
/* Set lvlInf to use the composite objective  Obj + wtInf*sInf */
/* after any switch to elastic mode. */
	if (*nncon > 0) {
/* Nonlinear constraints. */
/* In elastic mode, use the composite objective: */
/* w1*Obj + w2*sInf,  with w1 = 1, w2 = wtInf. */
	    lvlinf = 1;
	} else {
/* Linear constraints.  In theory, the subproblem should be */
/* feasible. If it is not, do not switch to Elastic mode. */
	    lemode = 0;
	    lvlinf = 0;
	}
	if (*nnl > 0) {
	    subopt = 0;
	} else {
	    subopt = -1;
	}
	needlu = lureq > 0;
	needx = needlu;
	s5qp_(&inform__, &c__5, probtag, elastc, &subopt, (U_fp)hprod, (U_fp)
		hprod1, (U_fp)mnrlog, gotr, &needlu, &typelu, &needx, lenr, m,
		 maxs, mbs, n, nb, ndegen, hvcalls, &ngqp0, &ngqp, nnobj0, 
		nnobj, &nnh0, &nnh, ns, itqp, &itqpmax, &itqptargt, itn, &
		lemode, &lvlinf, mnrprt, minimz, iobj, sclobj, objadd, objqp, 
		&targth, &targtz, tolfp, tolqpk, tolx, ninf, sinf, wtinf, 
		pinorm, &rgnorm, ne, nlocj, &locj[1], &indj[1], &jcol[1], &
		hetype[1], &hestat[1], &hfeas[1], &hs[1], &kbs[1], &ascale[1],
		 &bl[1], &bu[1], &blbs[1], &bubs[1], &gbs[1], &gobj[1], &gqp[
		1], &hdx[1], &pbs[1], &pi[1], &r__[1], &rc[1], &rg[1], nncon0,
		 nncon, &qprhs[1], nnl0, nnl, &x[1], &xqp[1], &xbs[1], &x[1], 
		&iy[1], &iy1[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1],
		 leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
		lenrw, (ftnlen)20, (ftnlen)8, (ftnlen)8);
/* inform      Result */
/* ------      ------ */
/*  >0         Fatal LU error */
/*   0         QP solution found */
/*  -1         The nonelastics are infeasible */
/*  -2         The QP subproblem is unbounded */
/*  -3         Too many iterations */
/*  -4         The QP subproblem has a weak minimizer */
/*  -5         Too many superbasics */
/*  -6         QP Hessian not positive semidefinite */
/*  -7         Z'g could not be made sufficiently small */
/*  -8         Ill-conditioned QP null-space basis */
	if (inform__ > 0) {
	    *iexit = inform__;
	    goto L900;
	}
	solved = inform__ == 0 || inform__ == -4;
	done = inform__ == -2 || inform__ == -3 || inform__ == -5;
	if (done) {
/* Relax */
	} else if (solved) {
/* Finish if there are no large nonlinear pi's. */
/* Otherwise, re-solve the QP in elastic mode */
	    if (! (*elastc)) {
		pinnln = dnormj_(nncon, &pi[1], &c__1);
		if (pinnln > *wtinf) {
		    *elastc = TRUE_;
		    solved = FALSE_;
		    s_wsfi(&io___194);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(
			    doublereal));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		}
	    }
	} else {
/* trouble while solving the QP. */
	    *gotr = FALSE_;
	    if (inform__ == -1) {
/* The nonelastics are infeasible. This should not happen. */
/* Phase 1 has already found a feasible point for the */
/* nonelastics, so the basis must be ill-conditioned. */
/* Refactorize with tighter tols and restart at the known */
/* feasible point.  Reduce the feasibility tol to try and */
/* prevent repeats. */
		s_wsfi(&io___195);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__22, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		*elastc = FALSE_;
		dcopy_(nb, &xqp0[1], &c__1, &xqp[1], &c__1);
	    } else if (inform__ == -6 || inform__ == -7) {
/* Indefinite Z'HZ  or large Z'g. */
/* Most likely an ill-conditioned Z'HZ. */
/* Try to reset the Hessian to I. */
		if (inform__ == -6) {
		    s_wsfi(&io___196);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		} else {
		    s_wsfi(&io___197);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		}
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		if (*hqntype != 2) {
		    *hqntype = 2;
		    s8h0_(hqntype, &nnh, u0ii, &iw[1], leniw, &rw[1], lenrw);
		    s_wsfi(&io___198);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		} else {
		    if (inform__ == -6) {
			targth = 1. / (eps * eps);
			lureq = 26;
		    } else if (inform__ == -7) {
			lureq = 21;
		    }
		    s2trylu_(itn, &lureq, ns, &lureq, &luok, &typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		}
	    } else if (inform__ == -8) {
/* condZ > targtZ  while forming Z'HZ for a freq. check. */
/* Refactorize B, possibly with a reduced factor tol. If */
/* the factor tol is already tight, accept Z, however bad. */
		s_wsfi(&io___199);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rw[192], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__25, ns, &lureq, &luok, &typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    targtz = plinfy;
		    luok = TRUE_;
		}
	    }
	}
	goto L600;
    }
/* end while */
L800:
    if (*ninf > 0) {
	info[4] = 1;
    }
    if (inform__ == 0) {
	info[5] = max(subopt,0);
    } else if (inform__ == -1) {
	*iexit = 15;
/* infeasible nonelastics */
    } else if (inform__ == -2) {
	info[5] = 3;
/* unbounded subproblem */
    } else if (inform__ == -3) {
	info[5] = 2;
	*iexit = -1;
/* too many iterations */
    } else if (inform__ == -4) {
	info[5] = 4;
/* weak QP solution */
    } else if (inform__ == -5 && feasbl) {
	info[5] = 5;
	*iexit = -2;
/* superbasic limit */
    } else if (inform__ == -5) {
	*iexit = 33;
/* superbasic limit */
    } else {
	*iexit = 44;
/* ill-conditioned null-space basis */
    }
L900:
    return 0;
} /* s8iqp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8iQP */
/* Subroutine */ int s8mrt_(integer *nncon, doublereal *fmrt, doublereal *
	gmrt, doublereal *hmrt, logical *incrun, doublereal *pendmp, 
	doublereal *penmax, doublereal *pennrm, doublereal *fv, doublereal *
	xpen, doublereal *y, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal eps0;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *), ddot_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    dnrm2_(integer *, doublereal *, integer *);
    static doublereal xpen0;
    extern /* Subroutine */ int ddscl_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal ppscl, xpeni;
    static logical boost;
    static doublereal ynorm, penold, penmin, pennew;
    static logical overfl;
    static doublereal rtundf, penlty;

/*     ================================================================== */
/*     s8mrt  computes the contributions to the merit function and its */
/*     directional derivative from the nonlinear constraints. */
/*     The penalty parameters  xPen(j)  are increased if */
/*     the directional derivative is not sufficiently negative. */

/*     On entry: */
/*         Fv     is the violation c(x) + A(linear)x - s,  where */
/*                s  minimizes the merit function with respect to the */
/*                nonlinear slacks only. */

/*     30 Dec 1991: First version based on Npsol 4.0 routine npmrt. */
/*     02 Nov 1996: Multipliers no longer updated here. */
/*     19 Jul 1997: Thread-safe version. */
/*     21 Oct 2000: Made compatible with SNOPT 6.1 */
/*     21 Oct 2000: Current version of s8mrt. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y;
    --xpen;
    --fv;
    --rw;

    /* Function Body */
    eps0 = rw[2];
    rtundf = rw[10];
    xpen0 = rw[89];
    overfl = FALSE_;
/*     Find the quantities that define  penMin, the vector of minimum */
/*     two-norm such that the directional derivative is one half of */
/*     approximate curvature   - (p)'H(p). */
/*     The factor  rtUndf  tends to keep  xPen  sparse. */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = fv[i__], abs(d__1)) <= rtundf) {
	    y[i__] = 0.;
	} else {
/* Computing 2nd power */
	    d__1 = fv[i__];
	    y[i__] = d__1 * d__1;
	}
    }
    ynorm = dnrm2_(nncon, &y[1], &c__1);
    d__1 = *gmrt + *hmrt * .5;
    ppscl = ddiv_(&d__1, &ynorm, &overfl);
    if (abs(ppscl) <= *penmax && ! overfl) {
/*        --------------------------------------------------------------- */
/*        Bounded  penMin  found.  The final value of  xPen(i)  will */
/*        never be less than  penMin(i).  A trial value  penNew  is */
/*        computed that is equal to the geometric mean of the previous */
/*        xPen  and a damped value of penMin.  The new  xPen  is defined */
/*        as  penNew  if it is less than half the previous  xPen  and */
/*        greater than  penMin. */
/*        --------------------------------------------------------------- */
	i__1 = *nncon;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = y[i__] / ynorm * ppscl;
	    penmin = max(d__1,0.);
	    xpeni = xpen[i__];
	    pennew = sqrt(xpeni * (*pendmp + penmin));
	    if (pennew < xpeni * .5) {
		xpeni = pennew;
	    }
	    xpeni = max(xpeni,penmin);
	    xpen[i__] = max(xpeni,xpen0);
	}
	penold = *pennrm;
	*pennrm = dnrm2_(nncon, &xpen[1], &c__1);
/*        --------------------------------------------------------------- */
/*        If  IncRun = true,  there has been a run of iterations in */
/*        which the norm of  xPen  has not decreased.  Conversely, */
/*        IncRun = false  implies that there has been a run of */
/*        iterations in which the norm of xPen has not increased.  If */
/*        IncRun changes during this iteration the damping parameter */
/*        PenDmp is increased by a factor of two.  This ensures that */
/*        xPen(j) will oscillate only a finite number of times. */
/*        --------------------------------------------------------------- */
	boost = FALSE_;
	if (*incrun && *pennrm < penold) {
	    boost = TRUE_;
	}
	if (! (*incrun) && *pennrm > penold) {
	    boost = TRUE_;
	}
	if (boost) {
/* Computing MIN */
	    d__1 = 1 / eps0, d__2 = *pendmp * 2.;
	    *pendmp = min(d__1,d__2);
	    *incrun = ! (*incrun);
	}
    }
/*     ------------------------------------------------------------------ */
/*     Compute the new value and directional derivative of the */
/*     merit function. */
/*     ------------------------------------------------------------------ */
    dcopy_(nncon, &fv[1], &c__1, &y[1], &c__1);
    ddscl_(nncon, &xpen[1], &c__1, &y[1], &c__1);
    penlty = ddot_(nncon, &y[1], &c__1, &fv[1], &c__1);
    *fmrt += penlty * .5;
    *gmrt -= penlty;
    return 0;
} /* s8mrt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine  s8mrt */
/* Subroutine */ int s8pphx_(U_fp hprod, integer *nnh, doublereal *x, 
	doublereal *hx, integer *status, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);

/*     ================================================================== */
/*     s8PPHx  defines the product  H*x  for the proximal-point QP */
/*     subproblem of snopt. */

/*     On exit,    Hx   = x. */

/*     23 Oct 1993: First version of s8PPHx. */
/*     02 Aug 2000: Current version. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hx;
    --x;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    dcopy_(nnh, &x[1], &c__1, &hx[1], &c__1);
    return 0;
} /* s8pphx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8PPHx */
/* Subroutine */ int s8qphx_(integer *nnh, doublereal *x, doublereal *hx, 
	integer *status, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, ftnlen cu_len)
{
/*     ================================================================== */
/*     s8qpHx is the argument qpHx for s5solv when s5solv is called from */
/*     one of the snOpt wrappers. */

/*     04 Dec 2004: First version of s8qpHx. */
/*     04 Dec 2004: Current version of s8qpHx. */
/*     ================================================================== */
/* Relax */
    /* Parameter adjustments */
    --hx;
    --x;
    cu -= 8;
    --iu;
    --ru;

    /* Function Body */
    return 0;
} /* s8qphx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8qpHx */
/* Subroutine */ int s8rand_(integer *leng, integer *neg, doublereal *g)
{
    static integer seeds[3];
    extern /* Subroutine */ int ddrand_(integer *, doublereal *, integer *, 
	    integer *);

/*     ================================================================== */
/*     s8rand  fills the array g with random numbers. */

/*     15 Nov 1991: First version of s8rand in s8aux. */
/*     30 Jun 1999: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --g;

    /* Function Body */
    if (*neg <= 0) {
	return 0;
    }
    seeds[0] = 1547;
    seeds[1] = 2671;
    seeds[2] = 3770;
    ddrand_(neg, &g[1], &c__1, seeds);
    return 0;
} /* s8rand_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8rand */
/* Subroutine */ int s8rc_(doublereal *sclobj, integer *minimz, integer *iobj,
	 integer *m, integer *n, integer *nb, integer *nnobj0, integer *nnobj,
	 integer *nncon, integer *nnjac, integer *negcon, integer *ne, 
	integer *nlocj, integer *locj, integer *indj, doublereal *jcol, 
	doublereal *gobj, doublereal *gcon, doublereal *pi, doublereal *rc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal dj;
    static integer ir;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal sgnobj;

/*     ================================================================== */
/*     s8rc   computes reduced costs rc = gObj - ( A  -I )'*pi, */
/*     using  gCon  as the top left-hand corner of A. */
/*     gCon, gObj and pi are assumed to exist. */

/*     s8rc   is called by s8SQP. */

/*     28 Sep 1993: First version, derived from m4rc. */
/*     31 Oct 1996: Min sum option added. */
/*     30 Oct 2000: Current version of s8rc. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --rc;
    --gobj;
    --gcon;
    --jcol;
    --indj;
    --locj;

    /* Function Body */
    l = 0;
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	dj = 0.;
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir <= *nncon) {
		++l;
		dj += pi[ir] * gcon[l];
	    } else {
		dj += pi[ir] * jcol[k];
	    }
	}
	rc[j] = -dj;
    }
    i__1 = *n;
    for (j = *nnjac + 1; j <= i__1; ++j) {
	dj = 0.;
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    dj += pi[ir] * jcol[k];
	}
	rc[j] = -dj;
    }
    dcopy_(m, &pi[1], &c__1, &rc[*n + 1], &c__1);
/*     Include the nonlinear objective gradient. */
    sgnobj = (doublereal) (*minimz);
    if (*nnobj > 0) {
	daxpy_(nnobj, &sgnobj, &gobj[1], &c__1, &rc[1], &c__1);
    }
    if (*iobj > 0) {
	rc[*n + *iobj] += sgnobj * *sclobj;
    }
    return 0;
} /* s8rc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8rc */
/* Subroutine */ int s8sclg_(integer *nnobj, doublereal *ascale, doublereal *
	gobj, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal grad, gdummy;

/*     ================================================================== */
/*     s8sclg  scales the objective gradient. */
/*     s8sclg is called by fgwrap only if modefg = 2. */
/*     Hence, it is used to scale known gradient elements (if any), */
/*     but is not called when missing gradients are being estimated */
/*     by s6dobj. */

/*     17 Feb 1992: First version. */
/*     16 Jul 1997: Thread-safe version. */
/*     02 Jan 2001: Current version of s8sclg. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --gobj;
    --ascale;
    --rw;

    /* Function Body */
    gdummy = rw[69];
/* definition of 'unset' value */
    i__1 = *nnobj;
    for (j = 1; j <= i__1; ++j) {
	grad = gobj[j];
	if (grad != gdummy) {
	    gobj[j] = grad * ascale[j];
	}
    }
    return 0;
} /* s8sclg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8sclg */
/* Subroutine */ int s8sclj_(integer *nncon, integer *nnjac, integer *negcon, 
	integer *n, doublereal *ascale, integer *ne, integer *nlocj, integer *
	locj, integer *indj, doublereal *gcon, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l, ir;
    static doublereal grad, cscale, gdummy;

/*     ================================================================== */
/*     s8sclJ  scales the Jacobian. */
/*     s8sclJ is called by fgwrap only if modefg = 2. */
/*     Hence, it is used to scale known gradient elements (if any), */
/*     but is not called when missing gradients are being estimated */
/*     by s6dcon. */

/*     17 Feb 1992: First version based on Minos routine m8sclj. */
/*     16 Jul 1997: Thread-safe version. */
/*     02 Dec 2001: Current version of s8sclJ. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --gcon;
    --ascale;
    --indj;
    --locj;
    --rw;

    /* Function Body */
    gdummy = rw[69];
/* definition of 'unset' value */
    l = 0;
    i__1 = *nnjac;
    for (j = 1; j <= i__1; ++j) {
	cscale = ascale[j];
	i__2 = locj[j + 1] - 1;
	for (k = locj[j]; k <= i__2; ++k) {
	    ir = indj[k];
	    if (ir > *nncon) {
		goto L300;
	    }
	    ++l;
	    grad = gcon[l];
	    if (grad != gdummy) {
		gcon[l] = grad * cscale / ascale[*n + ir];
	    }
	}
L300:
	;
    }
    return 0;
} /* s8sclj_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8sclJ */
/* Subroutine */ int s8sinf_(integer *n, integer *nb, integer *nncon, 
	doublereal *tolx, integer *ninf, doublereal *sinf, doublereal *bl, 
	doublereal *bu, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal tol, slack, violl, violu;

/*     ================================================================== */
/*     s8sInf computes the sum of infeasibilities of the nonlinear slacks */
/*     using bl, bu and x. */

/*     10 Jan 1997: First version of s8sInf. */
/*     30 Oct 2000: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;

    /* Function Body */
    *ninf = 0;
    *sinf = 0.;
    tol = *tolx;
/*     See how much  x(n+1:n+nnCon) violates its bounds. */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	slack = x[j];
	violl = bl[j] - slack;
	violu = slack - bu[j];
	if (violl > tol || violu > tol) {
	    ++(*ninf);
	    *sinf += max(violl,violu);
	}
    }
    return 0;
} /* s8sinf_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8sInf */
/* Subroutine */ int s8step_(logical *centrl, logical *usefls, integer *nb, 
	integer *nncon, integer *nnobj, integer *nmajor, integer *nskip, 
	doublereal *step, doublereal *stepmn, doublereal *steplm, doublereal *
	stepmx, doublereal *tolz, doublereal *xdnorm, doublereal *xnorm, 
	doublereal *bl, doublereal *bu, doublereal *x, doublereal *dx, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal res;
    extern doublereal ddiv_(doublereal *, doublereal *, logical *);
    static doublereal tolp, bigdx;
    static integer gotfd;
    static doublereal xdlim, pivot, fdint1, pivabs;
    static logical overfl, switch__;
    static doublereal stepqp, tolpiv;
    static integer hqntype;

/*     ================================================================== */
/*     s8step  finds the maximum, minimum and initial value for the */
/*     linesearch step. */

/*     For problems with nonlinear constraints, the maximum step stepmx */
/*     is one.  If there are only linear constraints the maximum step is */
/*     the largest step such that x + step*dx  reaches one of its bounds. */

/*     All step sizes are subject to the user-specified limit  steplm. */

/*     04 Dec 1992: First version of s8step based on npsol routine npalf. */
/*     31 Mar 2000: Updated for SNOPT 6.1. */
/*     19 Mar 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --dx;
    --x;
    --bu;
    --bl;
    --iw;
    --rw;

    /* Function Body */
    tolpiv = rw[60];
/* excludes small elements of y. */
    bigdx = rw[72];
/* unbounded step. */
    fdint1 = rw[76];
/* (1) forwrd diff. interval */
    xdlim = rw[80];
/* Step limit */
    gotfd = iw[183];
/* > 0 => some differences needed */
    hqntype = iw[202];
/* Current approximate Hessian type */
    overfl = FALSE_;
/*     ================================================================== */
/*     switch  indicates if there is an option to switch to */
/*             central differences to get a better search direction. */
/*     stepQP  is the step predicted by the QP subproblem (usually 1). */
/*     stepmx  is the largest feasible steplength subject to a */
/*             user-defined limit, bigdx, on the change in  x. */
/*     step    is initialized subject to a user-defined limit, xdlim. */
/*     ================================================================== */
    if (*nncon == 0 && *nnobj == 0) {
/* LP !! */
	*step = 1.;
	*stepmn = 1.;
	*steplm = 1.;
	*stepmx = 1.;
    } else {
	switch__ = gotfd > 0 && ! (*centrl);
	*stepmn = 0.;
	if (*usefls && switch__) {
	    *stepmn = fdint1 * (*xnorm + 1.) / *xdnorm;
	}
	stepqp = 1.;
	if (*nncon > 0 && (*nskip == 0 || hqntype != 2)) {
	    *stepmx = 1.;
	} else {
	    tolp = tolpiv * *xdnorm;
	    *stepmx = ddiv_(&bigdx, xdnorm, &overfl);
	    *step = *stepmx;
	    j = 1;
/* +          while (j .le. nb  .and.  step .gt. stepQP) do */
L100:
	    if (j <= *nb && *step > stepqp) {
		pivot = dx[j];
		pivabs = abs(pivot);
		if (pivabs > tolp) {
		    if (pivot <= 0.) {
			res = x[j] - bl[j];
			if (*step * pivabs > res) {
			    *step = res / pivabs;
			}
		    } else {
			res = bu[j] - x[j];
			if (*step * pivabs > res) {
			    *step = res / pivabs;
			}
		    }
		}
		++j;
		goto L100;
/* +          end while */
	    }
	    *step = max(*step,stepqp);
	    if (*step < stepqp + *tolz) {
		*step = stepqp;
	    }
	    *stepmx = *step;
	}
	d__1 = (*xnorm + 1.) * xdlim;
	*steplm = ddiv_(&d__1, xdnorm, &overfl);
	if (*nmajor <= 1) {
/* Computing MIN */
	    d__1 = *steplm, d__2 = ddiv_(&c_b5, xdnorm, &overfl);
	    *steplm = min(d__1,d__2);
	}
	*stepmx = min(*steplm,*stepmx);
	*step = min(*steplm,1.);
    }
    return 0;
} /* s8step_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8step */
/* Subroutine */ int s8sopt_(integer *n, integer *nb, integer *nncon, integer 
	*hestat, doublereal *pinorm, doublereal *tolz, doublereal *wtinf, 
	doublereal *bl, doublereal *bu, doublereal *fv, doublereal *x, 
	doublereal *ycon, doublereal *xpen, doublereal *fx)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, blj, con, buj;
    static integer jes;
    static doublereal fvi, tol, fvbl, fvbu, yconi, xpeni, fvlow, fvupp, 
	    dfvmax;

/*     ================================================================== */
/*     s8sOpt computes the vector of nonlinear constraint violations: */
/*        Fv = fCon + A(linear)*x - (optimal nonlinear slacks) */

/*     The optimal nonlinear slacks are computed as follows: */
/*     (1) Feasible  nonlinear slacks are adjusted so that they minimize */
/*         the merit function subject to  x  and  yCon  being held */
/*         constant. */
/*     (2) Infeasible slacks are compared with the true nonlinear slacks, */
/*         and, if necessary, they are adjusted so that the sum of */
/*         infeasibilities is reduced. */

/*     If yCon is zero, the violation can be set to any value without */
/*     changing the merit function.  In this case we choose the slack to */
/*     so that  the violation is zero (subject to the constraints above). */

/*     On entry, */
/*        x   =  the current x. */
/*        Fx  =  fCon + A(linear)*x,   defined in s8Fx. */

/*     On exit, */
/*        x   =  x containing the optimal slacks. */
/*        Fv  =  fCon + A(linear)*x - optimal slacks. */
/*        Fx  =  unaltered. */

/*     09 Jan 1992: First version based on Npsol routine npslk. */
/*     09 Oct 1996: First infeasible slack version. */
/*     28 Jul 2003: Test hEstat for slacks that are allowed to move. */
/*     29 Apr 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --bu;
    --bl;
    --hestat;
    --fx;
    --xpen;
    --ycon;
    --fv;

    /* Function Body */
    tol = *tolz * *pinorm;
/* piNorm ge 1 */
    i__1 = *nncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + i__;
	con = fx[i__];
	xj = x[j];
	fvi = con - xj;
	xpeni = xpen[i__];
	yconi = ycon[i__];
	blj = bl[j];
	buj = bu[j];
	fvbu = con - buj;
	fvbl = con - blj;
/* -------------------------------------------------------------- */
/* Redefine  xj  so that it minimizes the merit function */
/* subject to upper and lower bounds determined by the current */
/* multipliers. For convenience (but not clarity),  instead of */
/* checking that  xj  is within these bounds, the violation */
/* Fvi = c - xj  is checked against  FvLow  and  FvUpp, the */
/* violations at the upper and lower bounds on xj. */
/* -------------------------------------------------------------- */
/* Impose   FvLow <= Fv(opt) <= FvUpp */
/* First, define default bounds (tbl, tbu). */
	dfvmax = (abs(fvi) + 1.) * 10.;
	fvlow = fvi - dfvmax;
	fvupp = fvi + dfvmax;
	jes = hestat[j];
	if (jes == 1 && xj <= blj) {
/* ----------------------------------------------------------- */
/* This slack is at or below its lower bound in elastic mode. */
/* ----------------------------------------------------------- */
	    if (yconi < 0.) {
/* xj is eligible to increase. */
/* Require                  bl <=  xj <= min( bu,tbu ). */
		fvlow = max(fvbu,fvlow);
		fvupp = fvbl;
	    } else if (yconi > 0.) {
/* xj is eligible to decrease and violate its lower bound. */
/* Require              -infty <=  xj <= bl */
		yconi -= *wtinf;
		fvlow = fvbl;
	    } else {
/* xj can either increase or decrease. */
/* Require              -infty <=  xj <= min( bu,tbu ). */
		fvlow = max(fvbu,fvlow);
	    }
	} else if (jes == 2 && xj >= buj) {
/* ----------------------------------------------------------- */
/* This slack is at or above its upper bound in elastic mode. */
/* ----------------------------------------------------------- */
	    if (yconi > 0.) {
/* xj is eligible to decrease. */
/* Require      max( bl, tbl ) <=  xj <= bu. */
		fvlow = fvbu;
		fvupp = min(fvbl,fvupp);
	    } else if (yconi < 0.) {
/* xj is eligible to increase and violate its upper bound. */
/* Require                  bu <=  xj <= +infty */
		yconi += *wtinf;
		fvupp = fvbu;
	    } else {
/* xj can either increase or decrease. */
/* Require      max( bl, tbl ) <=  xj <= +infty */
		fvupp = min(fvbl,fvupp);
	    }
	} else {
/* ----------------------------------------------------------- */
/* Feasible slack.  xj can move either way. */
/* ----------------------------------------------------------- */
/* Require      max( bl, tbl ) <=  xj <= min( bu,tbu ). */
	    fvlow = max(fvbu,fvlow);
	    fvupp = min(fvbl,fvupp);
	}
	if (abs(yconi) <= tol) {
/* Computing MIN */
	    d__1 = max(0.,fvlow);
	    fvi = min(d__1,fvupp);
	} else if (xpeni >= *tolz) {
	    if (yconi >= xpeni * fvupp) {
		fvi = fvupp;
	    } else if (yconi <= xpeni * fvlow) {
		fvi = fvlow;
	    } else {
		fvi = yconi / xpeni;
	    }
	}
	xj = con - fvi;
	fv[i__] = fvi;
	x[j] = xj;
    }
    return 0;
} /* s8sopt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s8sOpt */
/* Subroutine */ int s8winf_(integer *job, logical *boostd, integer *itn, 
	doublereal *gnorm, doublereal *wtinf0, doublereal *wtinf, doublereal *
	wtmax, doublereal *weight, doublereal *wtfac, doublereal *wtscal, 
	integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Elastic weight increase"
	    "d to \002,1p,e11.3)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    static doublereal newwt;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___268 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s8wInf  initializes or updates the elastic weight  wtInf. */
/*     The elastic weight is given by  wtInf = wtScal*weight, */
/*     where wtScal is some scale-dependent quantity (fObj here). */
/*     wtInf is increased by redefining weight as weight*wtFac, where */
/*     wtFac is a constant factor. */

/*     weight, wtFac and wtScal are 'saved' local variables. */

/*     20 Feb 1997: First version of s8wInf. */
/*     27 Apr 2001: wtMax introduced as parameter instead of local. */
/*     03 Aug 2003: snPRNT adopted. */
/*     03 Aug 2003: Current version of s8wInf. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (*job == 0) {
/*        Set the weight. */
/*        weight is the ``unscaled'' weight on the infeasibilities. */
/*        wtScal is a scale factor based on the current gradient. */
	*wtscal = *gnorm;
	*wtfac = 10.;
	*weight = *wtinf0;
	*wtinf = *wtscal * *weight;
    } else if (*job == 1) {
/*        If possible, boost the weight. */
/* Computing MIN */
	d__1 = *wtfac * *weight;
	newwt = min(d__1,*wtmax);
	*boostd = newwt > *weight;
	if (*boostd) {
	    *weight = newwt;
	    *wtinf = *weight * *wtscal;
	    *wtfac *= 10.;
	    s_wsfi(&io___268);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*wtinf), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    return 0;
} /* s8winf_ */

