/* ./src/sn55qp.f -- translated by f2c (version 20090411).
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

static integer c__0 = 0;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__23 = 23;
static integer c__21 = 21;
static integer c__11 = 11;
static integer c__27 = 27;
static integer c__2 = 2;
static doublereal c_b76 = -1.;
static integer c__31 = 31;
static integer c__26 = 26;
static integer c__25 = 25;
static doublereal c_b146 = 1.;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  sn55qp.f */

/*     s5QP */
/*     s5chkp   s5chzq   s5getp   s5getR   s5Hfac   s5Hz     s5QPfg */
/*     s5QPit   s5Rchk   s5Rcol   s5rg     s5Rsng   s5Sdel   s5Zp */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5qp_(integer *iexit, integer *prob, char *probtag, 
	logical *elastc, integer *subopt, U_fp hprod, U_fp hprod1, S_fp qplog,
	 logical *gotr, logical *needlu, integer *typelu, logical *needx, 
	integer *lenr, integer *m, integer *maxs, integer *mbs, integer *n, 
	integer *nb, integer *ndegen, integer *hvcalls, integer *ngqp0, 
	integer *ngqp, integer *ngobj0, integer *ngobj, integer *nnh0, 
	integer *nnh, integer *ns, integer *itqp, integer *itqpmax, integer *
	itqptargt, integer *itn, integer *lemode, integer *lvlinf, integer *
	prtlvl, integer *minimz, integer *iobj, doublereal *sclobj, 
	doublereal *objadd, doublereal *objqp, doublereal *targth, doublereal 
	*targtz, doublereal *tolfp, doublereal *tolqp, doublereal *tolx, 
	integer *ninf, doublereal *sinf, doublereal *wtinf, doublereal *
	pinorm, doublereal *rgnorm, integer *ne, integer *nloca, integer *
	loca, integer *inda, doublereal *acol, integer *hetype, integer *
	hestat, integer *hfeas, integer *hs, integer *kbs, doublereal *ascale,
	 doublereal *bl, doublereal *bu, doublereal *blbs, doublereal *bubs, 
	doublereal *gbs, doublereal *gobj, doublereal *gqp, doublereal *hdx, 
	doublereal *pbs, doublereal *pi, doublereal *r__, doublereal *rc, 
	doublereal *rg, integer *nrhs0, integer *nrhs, doublereal *rhs, 
	integer *lenx0, integer *nx0, doublereal *x0, doublereal *x, 
	doublereal *xbs, doublereal *xfreez, integer *iy, integer *iy1, 
	doublereal *y, doublereal *y1, doublereal *y2, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen probtag_len, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1030[] = "(\002 Itn\002,i7,\002: Infeasible nonelastics."
	    "  Num =\002,i5,1p,\002  Sum of Infeasibilities =\002,e8.1)";
    static char fmt_1610[] = "(\002 Itn\002,i7,\002: Suboptimize: \002,i7"
	    ",\002 new superbasics\002)";
    static char fmt_1620[] = "(\002 Itn\002,i7,\002: Suboptimize: \002,i7"
	    ",\002 minor iterations\002)";
    static char fmt_1040[] = "(\002 Itn\002,i7,\002: Infinite pi-vector\002)";
    static char fmt_1600[] = "(\002 Itn\002,i7,\002: Singularity after a "
	    "\002,\002bound swap.  Basis refactorized\002)";
    static char fmt_8050[] = "(\002 Itn\002,i7,\002: Infeasible \002,a)";
    static char fmt_8060[] = "(\002 Itn\002,i7,\002: Elastic Phase 1 -- maki"
	    "ng \002,\002non-elastic variables feasible\002)";
    static char fmt_1010[] = "(\002 Biggest dj =\002,1p,e11.3,\002 (variabl"
	    "e\002,i7,\002)\002,\002    norm rg =\002,e11.3,\002   norm pi "
	    "=\002,e11.3)";
    static char fmt_1020[] = "(\002 Norm rg =\002,1p,e11.3,\002   norm pi "
	    "=\002,e11.3)";
    static char fmt_1000[] = "(\002 ==> LU file has increased by a factor o"
	    "f\002,f6.1)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal c6;
    static integer jq, kp;
    static doublereal obj;
    static integer jbq, jbr;
    static doublereal djq;
    static integer nbs, jsq, jsr, lrs;
    static char str[132];
    static doublereal djq0, eps0, eps2;
    extern /* Subroutine */ int s5rg_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), s5hs_(integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *), s5hz_(integer *, U_fp, 
	    U_fp, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    static integer kfac, nfac, eigh, kchk;
    static doublereal bold;
    static logical newb;
    static integer klog;
    static logical gote, goth;
    static integer kprc, ksav;
    static logical prt10, luok;
    static integer maxr, nfix[2];
    static doublereal drsq;
    static logical newx;
    static doublereal step;
    extern /* Subroutine */ int s5inf_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tolx0;
    static logical needf;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical chkpi;
    static integer ninfe;
    static logical needv;
    static doublereal hdmax, sinfe, drmin;
    static logical newsb;
    static integer lumax, lureq;
    static logical newlu;
    static integer ksumm, nsmax, nswap;
    static doublereal anorm, drmax, normg, pivot, rgtol[2];
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
	     integer *, doublereal *, integer *), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), s5dgen_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), s5einf_(integer *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), s5egrd_(integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *), 
	    s6rcnd_(integer *, integer *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *), s2bsol_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s5qpfg_(U_fp, U_fp, integer *, integer *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen), s5pric_(logical *, 
	    logical *, logical *, logical *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s5eset_(integer *, integer *, integer *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static logical jstph1;
    extern /* Subroutine */ int s4ksav_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, ftnlen);
    static logical chkfea;
    extern /* Subroutine */ int s5rsng_(integer *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *), 
	    s5setp_(integer *, integer *, logical *, doublereal *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *), 
	    s2unpk_(integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static integer lusiz0;
    extern /* Subroutine */ int s5qpit_(integer *, U_fp, U_fp, logical *, 
	    logical *, logical *, logical *, logical *, logical *, logical *, 
	    logical *, logical *, logical *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    s5setx_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *);
    static logical feasbl;
    static integer kdegen;
    static logical deadpt;
    static doublereal infbnd;
    static logical needlm, checkx, needpi, unbndd, jstfea;
    static doublereal featol;
    static logical posdef, maxref, incres;
    static doublereal objslk, condhz;
    static logical qpdone;
    static integer nelast;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__;
    static logical bndswp, singlr;
    static integer itnlim;
    static logical gotgqp;
    static integer itnfix;
    static logical optiml, prtlog, usegqp;
    static integer jqsave;
    static logical statpt;
    static integer lvltol, kprprt, mnewsb, nfmove, nfreez, nonopt;
    static logical prtsum;
    static integer nuncon, sqstat;
    static doublereal bgrwth, djqprt, objprt, rgtest, rowerr, sgnobj, tolinc, 
	    weight;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2gathr_(integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *), s5einit_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *), s5ewrap_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), s2trylu_(integer *, integer *, 
	    integer *, integer *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___79 = { 0, str, 0, fmt_1030, 132, 1 };
    static icilist io___83 = { 0, str, 0, fmt_1610, 132, 1 };
    static icilist io___84 = { 0, str, 0, fmt_1620, 132, 1 };
    static icilist io___85 = { 0, str, 0, fmt_1040, 132, 1 };
    static icilist io___94 = { 0, str, 0, fmt_1600, 132, 1 };
    static icilist io___102 = { 0, str, 0, fmt_8050, 132, 1 };
    static icilist io___103 = { 0, str, 0, fmt_8060, 132, 1 };
    static icilist io___104 = { 0, str, 0, fmt_1010, 132, 1 };
    static icilist io___105 = { 0, str, 0, fmt_1020, 132, 1 };
    static icilist io___110 = { 0, str, 0, fmt_1000, 132, 1 };


/* ================================================================= */
/* s5QP   solves a linear or quadratic program. */
/* The problem type can be: */
/*   Prob = 0 FP   feasible point only */
/*   Prob = 1 LP   LP problem */
/*   Prob = 2 QP   QP problem */
/*   Prob = 3 FPE  feasible point for equalities only */
/*   Prob = 4 FPS  feasible point for QP subProblem */
/*   Prob = 5 QPS  QP subproblem */
/*   Prob = 6 QPP  FP subproblem with proximal point objective */

/* ngQP = max( nnH, ngObj ) */

/* The optimization can pass through the following phases: */

/*   Phase 1               find a feasible point for all variables */

/*   Elastic Phase 1       make the non-elastic variables feasible */
/*                         while allowing infeasible elastics */

/*   Phase 2               minimize the objective */

/*   Elastic Phase 2       minimize a composite objective while */
/*                         keeping the non-elastics feasible */

/*                         In this phase, lvlInf means the followin */

/*             lvlInf = 0  zero     weight on the infeasibilities */
/*                                  (infeasibillities are ignored) */
/*                      1  finite   weight on the infeasibilities */
/*                      2  infinite weight on the infeasibilities */
/*                                  (the objective is ignored) */

/* The array kBS is a permutation on the column indices. */
/* kBS(1  :m )    holds the col. indices of the basic variables. */
/* kBS(m+1:m+nS)  holds the col. indices of the superbasic variable */
/*                These nS columns indices must have hs(j) = 2. */

/*  iExit       Result */
/*  -----       ------ */
/*   >0         Fatal LU error */
/*    0         QP solution found */
/*   -1         QP is infeasible */
/*   -2         QP is unbounded */
/*   -3         Too many iterations */
/*   -4         Weak QP minimizer */
/*   -5         Too many superbasics */
/*   -6         QP Hessian not positive semidefinite after pricing */
/*   -7         Z'g could not be made sufficiently small */
/*   -8         Ill-conditioned Z */

/* 30 Sep 1991: First version of s5QP  based on Qpsol routine qpcor */
/* 29 Oct 1993: QP objective computed separately. */
/* 19 May 1995: Bordered Hessian updated. */
/* 30 Jul 1995: Border updates removed. */
/* 04 Jan 1996: Positive semi-definite H treated correctly. */
/* 20 Jul 1996: Slacks changed to be the row value. */
/* 09 Aug 1996: First Min Sum version. */
/* 15 Jul 1997: Thread-safe version. */
/* 02 Feb 1998: Piecewise linear line search added. */
/* 07 Nov 1998: Explicit Hessian option added. */
/* 24 Dec 1999: Sub-optimization option added. */
/* 25 Jul 2001: Exit on nS > maxR activated. */
/* 30 Jul 2003: Superbasic slacks allowed. */
/* 02 Aug 2003: snEXIT and snPRNT adopted. */
/* 24 Dec 2003: pi checked for NaN and Inf entries. */
/* 16 May 2006: Explicit target itQP added */
/* 26 May 2013: infBnd used to identify infinite bounds. */
/* ================================================================= */
/*     ------------------------------------------------------------------ */
/* condition estimate of Z */
/* size of L0 */
/* size of initial  U */
/* size of current  L */
/* size of current  U */
/* phase 1 dj tol for p.p. */
/* phase 2 dj tol for p.p. */
/* current optimality tol */
/* xBS(kObj) is the obj. slack */
/* itns since last factorize */
/* number of LU mods */
/* (on/off) log     status */
/* (on/off) summary status */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor heading in log file */
/* >0 => Minor heading for iSumm */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
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
    --xfreez;
    --x;
    --rc;
    --bu;
    --bl;
    --ascale;
    --hs;
    --hestat;
    --hetype;
    --gqp;
    --hdx;
    --acol;
    --inda;
    --loca;
    --gobj;
    --rhs;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    eps2 = rw[4];
/* eps**(1/2)       IEEE DP  1.49e-08 */
    infbnd = rw[70];
/* definition of an infinite bound */
    maxr = iw[52];
/* max columns of R. */
    kchk = iw[58];
/* check (row) frequency */
    kfac = iw[59];
/* factorization frequency */
    ksav = iw[60];
/* save basis map */
    klog = iw[61];
/* log/print frequency */
    ksumm = iw[62];
/* Summary print frequency */
    kdegen = iw[63];
/* max. expansions of featol */
    itnlim = iw[89];
/* limit on total iterations */
    mnewsb = iw[95];
/* max # of new superbasics */
    eigh = iw[200];
/* -1,0,1 for indef, psd and pdef H */
    nfac = iw[210];
/* # of LU factorizations */
    *iexit = 0;
    sqstat = 0;
    c6 = max(1e-6,eps2);
    if (nfac > 0) {
	lusiz0 = iw[171] + iw[172];
	lumax = lusiz0 << 1;
    }
    prt10 = *prtlvl >= 10;
    prtlog = *prtlvl >= 1 && (*itqp % klog == 0 && *itqp != 0 || klog == 1);
    prtsum = *prtlvl >= 1 && (*itqp % ksumm == 0 && *itqp != 0 || ksumm == 1);
    iw[218] = 0;
    iw[219] = 0;
    if (prtlog) {
	iw[218] = 1;
    }
    if (prtsum) {
	iw[219] = 1;
    }
/*     ------------------------------------------------------------------ */
/*     s5QP operates in either ``Normal'' or ``Elastic'' mode. */
/*     Everything is normal unless a weighted sum is being minimized or */
/*     the constraints are infeasible. */
/*     The logical feasbl refers to the feasibility of the nonelastics. */
/*     wtInf  is the optional parameter Infeasibility Weight. */
/*     ------------------------------------------------------------------ */
    feasbl = FALSE_;
    gote = FALSE_;
    goth = *nnh > 0;
    gotgqp = *ngqp > 0;
/*     jstPh1 = stop at the end of phase 1 (either regular or elastic) */
    jstph1 = *prob == 0 || *prob == 3 || *prob == 4;
/*     The phase 2 objective is F1 + wtInf*F2. */
    if (*elastc) {
	needf = *lvlinf != 2;
/* F1 required in phase 2 */
	needv = *lvlinf != 0;
/* F2 required in phase 2 */
    } else {
	needf = TRUE_;
	needv = FALSE_;
    }
    objslk = 0.;
    *objqp = 0.;
    obj = 0.;
    pivot = 0.;
    step = 0.;
    ninfe = 0;
    jq = 0;
    djq = 0.;
    djq0 = 0.;
    djqprt = 0.;
    jbq = 0;
/* x(jBq) is the incoming  BS */
    jbr = 0;
/* x(jBr) is the outgoing  BS */
    jsq = 0;
/* x(jSq) is the incoming SBS */
    jsr = 0;
/* x(jSr) is the outgoing SBS */
    jqsave = 0;
    kprprt = 0;
    sgnobj = (doublereal) (*minimz);
    rgtol[0] = min(*tolqp,c6);
/* relaxed   rgTol */
    rgtol[1] = eps0;
/* stringent rgTol */
    lvltol = 2;
/* working   rgTol */
    rw[184] = 100.;
    rw[185] = 100.;
    kprc = 0;
/* last sec scanned in part. prc */
    lureq = 0;
    bndswp = FALSE_;
    chkfea = TRUE_;
    chkpi = TRUE_;
    deadpt = FALSE_;
    needpi = TRUE_;
    newlu = TRUE_;
    newx = FALSE_;
    unbndd = FALSE_;
    posdef = FALSE_;
    qpdone = FALSE_;
/*     nUncon  counts the number of unconstrained (i.e., Newton) steps. */
/*             If the test for a minimizer were scale-independent, */
/*             Uncon would never be larger than 1. */
/*     nfmove  counts the number of times that the QP obj is decreased, */
    nfmove = 0;
    nuncon = 0;
/*     subopt nonzero implies that optimization occurs with a subset of */
/*     the variables frozen at their initial values. */
/*     During suboptimization, nFreez is the number of frozen variables. */
    nfreez = 0;
    nsmax = *ns + mnewsb;
    s5hs_(&c__0, nb, &bl[1], &bu[1], &hs[1], &x[1]);
    s5dgen_(&inform__, &c__0, prtlvl, nb, ninf, itn, &featol, tolx, &tolinc, &
	    hs[1], &bl[1], &bu[1], &x[1], &itnfix, nfix, &tolx0, &iw[1], 
	    leniw, &rw[1], lenrw);
/*     ======================Start of main loop========================== */
/* +    do while (.not. QPdone  .and.  iExit .eq. 0) */
L100:
    if (! qpdone && *iexit == 0) {
/* ============================================================== */
/* Check the initial  x  and move it onto  ( A  -I )*x = b. */
/* If needLU is true, this will require a basis factorization. */
/* ============================================================== */
/* If necessary,  factorize the basis  ( B = LU ) and set x. */
/* If needLU is false on entry to s5QP, the first call to s2Bfac */
/* will try to use existing factors. */
/* If needLU is true on entry to s5QP, an LU factorization of */
/* type typeLU is computed. */

/* The reason for the requested LU is as follows. */

/* LUreq =  0  First LU for a given subproblem */
/* LUreq =  1  Frequency */
/* LUreq =  2  LU nonzeros increased */
/* LUreq =  3 */
/* LUreq =  4 */
/* LUreq =  5  Singular after LU mod */
/* LUreq =  6  Unstable LU mod (growth in new column of U) */
/* LUreq =  7  Not enough memory */
/* LUreq =  8 */
/* LUreq =  9 */
/* LUreq = 10  Row error in setx */
/* LUreq = 11  Big  dx   in setx */

/* LUreq = 20 */
/* LUreq = 21  Iterative refinement failed in QP */
/* LUreq = 22  Unbounded QP */
/* LUreq = 23  Infeasibility after refactorization */
/* LUreq = 24  Small directional derivative in QP */
/* LUreq = 25  Ill-conditioned Z in QP */
/* LUreq = 26  Indefinite Z'HZ in QP */
/* LUreq = 27  R singular after bound swap in QP */
/* -------------------------------------------------------------- */
	jstfea = FALSE_;
	if (lureq > 0) {
	    *needlu = TRUE_;
	}
	if (*needx || *needlu) {
	    s2bfac_(iexit, typelu, needlu, &newlu, &newb, iobj, itn, prtlvl, &
		    lureq, m, mbs, n, nb, nnh, ns, &nswap, ne, nloca, &loca[1]
		    , &inda[1], &acol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &
		    blbs[1], &bubs[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], &
		    iy[1], &iy1[1], &y[1], &y1[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    if (newlu) {
		lusiz0 = iw[171] + iw[172];
		lumax = lusiz0 << 1;
		*gotr = FALSE_;
/* Reset R. */
		if (prt10) {
		    iw[223] = 1;
		}
	    }
	    if (*iexit != 0) {
		goto L100;
	    }
	    gote = FALSE_;
/* Check hEstat in elastic mode. */
	    needpi = TRUE_;
/* Recalculate the pi's. */
	    *needx = FALSE_;
	    newx = TRUE_;
	    chkfea = TRUE_;
	    chkpi = TRUE_;
	    pivot = 0.;
	    jqsave = 0;
	    nuncon = 0;
	}
	nbs = *m + *ns;
	newsb = FALSE_;
	*ninf = 0;
	*sinf = 0.;
	optiml = FALSE_;
	dload_(&nbs, &c_b5, &gbs[1], &c__1);
	normg = 1.;
	if (*elastc && ! gote) {
/*           ------------------------------------------------------------ */
/*           Reset blBS and buBS for any violated elastics. */
/*           These values are used in s5step. */
/*           Strictly feasible elastics are returned to normality. */
/*           ------------------------------------------------------------ */
	    s5eset_(&nbs, nb, &nelast, &featol, &infbnd, &hetype[1], &hestat[
		    1], &kbs[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &xbs[1]);
	    gote = TRUE_;
	}
	if (chkfea) {
/*           In Phase 1 or just after a factorize, check the feasibility */
/*           of the basic and superbasic non-elastics. */
/*           jstFea  indicates that we have just become feasible. */
/*           jstFea is turned off once a step is taken. */
	    s5inf_(&nbs, &featol, &infbnd, ninf, sinf, &hfeas[1], &blbs[1], &
		    bubs[1], &gbs[1], &xbs[1]);
	    if (*ninf > 0) {
/*              Non-elastics are infeasible. */
/*              If necessary, switch back to the feasibility phase, after */
/*              refactorization (possibly with tighter tols). */
/*              Print something if the basis has just been refactorized. */
		if (feasbl) {
		    s2trylu_(itn, &c__23, ns, &lureq, &luok, typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 11;
		    }
		    feasbl = FALSE_;
		    goto L100;
		}
		*gotr = FALSE_;
		if (prt10 && iw[215] == 0) {
		    s_wsfi(&io___79);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
		}
	    }
/*           feasbl = true means that the nonelastics are feasible. */
/*                    This defines the start of Phase 2. */
	    if (! feasbl) {
		jstfea = *ninf == 0;
	    }
	    feasbl = *ninf == 0;
	    chkfea = *ninf > 0;
	}
/* if chkFea */
	if (*elastc) {
/*           ------------------------------------------------------------ */
/*           Compute the sum of infeasibilities of the elastic variables. */
/*           ------------------------------------------------------------ */
	    s5einf_(nb, &nbs, &hestat[1], &kbs[1], &ninfe, &sinfe, &bl[1], &
		    bu[1], &x[1]);
	    *ninf += ninfe;
	    *sinf += sinfe;
	}
	if (feasbl && jstph1) {
/* The non-elastics are feasible.  Exit. */
	    condhz = 0.;
	    djqprt = 0.;
	    *rgnorm = 0.;
	    dload_(m, &c_b5, &pi[1], &c__1);
	    *pinorm = 1.;
/* pinorm = max(norm(pi), 1.0) */
	    deadpt = FALSE_;
	    optiml = TRUE_;
	} else {
	    if (feasbl) {
/*              --------------------------------------------------------- */
/*              Feasible for the nonelastics. */
/*              (Elastc = false means no elastics.) */
/*              --------------------------------------------------------- */
/*              If just feasible, compute the QP objective (and gradient) */
/*              and R. */
		objslk = 0.;
		if (*iobj != 0) {
		    objslk = xbs[iw[205]] * *sclobj;
		}
		obj = sgnobj * objslk;
		if (jstfea || newx) {
		    if (needf) {
/*                    =================================================== */
/*                    Initialize the QP objective and gradient. */
/*                    ObjQP is the linear plus quadratic term of the */
/*                    objective (not scaled by sgnObj).   It is updated */
/*                    after each QP step. */
/*                    =================================================== */
			if (gotgqp) {
			    if (*hvcalls == 0) {
				sqstat = 1;
			    }
			    s5qpfg_((U_fp)hprod, (U_fp)hprod1, ngqp, ngobj0, 
				    ngobj, nnh, &sqstat, hvcalls, objqp, &
				    gobj[1], &gqp[1], lenx0, nx0, &x0[1], &x[
				    1], &y[1], cu + 8, lencu, &iu[1], leniu, &
				    ru[1], lenru, cw + 8, lencw, &iw[1], 
				    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)
				    8);
			    obj += sgnobj * *objqp;
			    sqstat = 0;
			}
			if (goth && ! (*gotr)) {
/*                       ------------------------------------------------ */
/*                       Load and factor the reduced Hessian. */
/*                       This happens after every LU factorize. */
/*                       If the reduced Hessian is not positive definite, */
/*                       reduce the LU factor tolerances to get a better */
/*                       conditioned Z. */
/*                       ------------------------------------------------ */
			    if (*ns > 0) {
				s5hz_(&inform__, (U_fp)hprod, (U_fp)hprod1, &
					maxr, lenr, minimz, m, mbs, n, nb, 
					nnh, ns, hvcalls, ne, nloca, &loca[1],
					 &inda[1], &acol[1], &hdmax, &rw[192],
					 targtz, &kbs[1], &r__[1], &y[1], &y1[
					1], &y2[1], cu + 8, lencu, &iu[1], 
					leniu, &ru[1], lenru, cw + 8, lencw, &
					iw[1], leniw, &rw[1], lenrw, (ftnlen)
					8, (ftnlen)8);
/*                          inform = -1, 0, >0 */
				if (inform__ != 0) {
				    if (inform__ == -1) {
					*iexit = -8;
/* Ill-conditioned Z */
				    } else {
					*iexit = inform__;
/* Fatal error in LU */
				    }
				    goto L100;
				}
				s5hfac_(&inform__, &c__1, itn, lenr, m, &maxr,
					 mbs, nb, ns, targth, &hdmax, &hs[1], 
					&kbs[1], &iy[1], &bl[1], &bu[1], &
					blbs[1], &bubs[1], &x[1], &xbs[1], &
					r__[1], &iw[1], leniw, &rw[1], lenrw);
/*                          inform = -2, -1, 0 */
				if (inform__ != 0) {
				    *iexit = -6;
/* Z'HZ not positive definite */
				    goto L100;
				}
			    }
/* nS > 0 */
			    *gotr = TRUE_;
			}
/* gotH and not gotR */
		    }
/* needf */
		    nbs = *m + *ns;
		    posdef = *gotr || *ns == 0;
		}
/*              --------------------------------------------------------- */
/*              Gather the QP gradient in BS order. */
/*              Assign the nonzero components of gBS. */
/*              --------------------------------------------------------- */
/* jstFea .or. newLU */
		if (needf) {
		    if (gotgqp) {
			s2gathr_(ngqp, &nbs, &kbs[1], &sgnobj, &gqp[1], &gbs[
				1]);
		    }
		    if (*iobj > 0) {
			gbs[iw[205]] = sgnobj * *sclobj;
		    }
		}
		if (*elastc && ninfe > 0 && needv) {
		    s5egrd_(nb, &nbs, wtinf, &hestat[1], &kbs[1], &gbs[1]);
		}
		normg = dnormi_(&nbs, &gbs[1], &c__1);
/*              --------------------------------------------------------- */
/*              See if it's time to suboptimize. */
/*              NOTE: We must not suboptimize if all steps have been */
/*              degenerate. */
/*              --------------------------------------------------------- */
		if (*subopt != 0 || nfmove == 0) {
/*                 Relax */
		} else {
		    if (*ns >= nsmax) {
			*subopt = 1;
			if (prt10) {
			    s_wsfi(&io___83);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&mnewsb, (ftnlen)sizeof(
				    integer));
			    e_wsfi();
			    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
			}
		    } else if (*itqp >= *itqptargt) {
			*subopt = 2;
			if (prt10) {
			    s_wsfi(&io___84);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&(*itqptargt), (ftnlen)
				    sizeof(integer));
			    e_wsfi();
			    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
			}
		    }
		}
	    }
/* feasible */
	    if (needpi) {
		dcopy_(m, &gbs[1], &c__1, &y[1], &c__1);
		s5setp_(&inform__, m, &chkpi, pinorm, &y[1], &pi[1], &iw[1], 
			leniw, &rw[1], lenrw);
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
		    } else {
/* pi is infinite or contains a NaN/Inf. */
			s_wsfi(&io___85);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
			s2trylu_(itn, &c__11, ns, &lureq, &luok, typelu, &iw[
				1], leniw, &rw[1], lenrw);
			if (! luok) {
			    *iexit = 43;
			}
		    }
		    goto L100;
		}
		needpi = FALSE_;
	    }
	    *rgnorm = 0.;
	    if (*ns > 0) {
		s5rg_(m, &nbs, n, ns, &eps0, ne, nloca, &loca[1], &inda[1], &
			acol[1], &gbs[1], &pi[1], &rg[1], rgnorm, &kbs[1]);
	    }
/*           ============================================================ */
/*           Determine if the reduced Hessian is positive definite. */
/*           ============================================================ */
	    condhz = 0.;
	    if (*gotr) {
		s6rcnd_(&maxr, ns, lenr, &r__[1], &drmax, &drmin, &condhz);
	    }
	    if (! posdef) {
		if (feasbl && *gotr) {
		    if (*ns > 0) {
			lrs = (*ns - 1) * maxr + (3 - *ns) * *ns / 2;
/* Magic formul */
/* Computing 2nd power */
			d__1 = r__[lrs];
			drsq = d__1 * d__1;
		    }
		    s5rsng_(&eigh, &posdef, &singlr, itn, &maxr, lenr, ns, &
			    drsq, &r__[1], &iw[1], leniw, &rw[1], lenrw);
		} else {
		    posdef = *ns == 0;
		}
	    }
/*           ============================================================ */
/*           Check for a stationary point.  Use a stringent rgTol after */
/*           a constrained step to help avoid false stationary points. */
/*           In theory, the reduced gradient is zero and the reduced */
/*           Hessian is positive definite after a bound swap. */

/*           If x is a minimizer,  reduced costs are calculated. */
/*           ============================================================ */
	    if (feasbl) {
		rw[186] = *tolqp;
	    } else {
		rw[186] = *tolfp;
	    }
	    rgtest = max(*pinorm,normg);
	    if (! feasbl) {
		statpt = *rgnorm <= rgtol[0] * rgtest;
	    } else if (nuncon >= 1) {
		statpt = *rgnorm <= rgtol[lvltol - 1] * rgtest;
	    } else {
		statpt = *rgnorm <= rgtol[1] * rgtest;
	    }
	    if (feasbl) {
		maxref = nuncon > 1;
		if ((maxref || bndswp) && ! statpt) {
/*                 If this point should be stationary but isn't. */
/*                 If possible, relax the reduced-gradient tolerance. */
		    if (lvltol == 2) {
			lvltol = 1;
			statpt = *rgnorm <= rgtol[lvltol - 1] * rgtest;
		    }
		    if (! statpt) {
			s2trylu_(itn, &c__21, ns, &lureq, &luok, typelu, &iw[
				1], leniw, &rw[1], lenrw);
			if (! luok) {
			    *iexit = -7;
/* Large Z'g */
			    goto L100;
			}
		    }
		}
		deadpt = statpt && needf && ! posdef;
	    }
	    if (statpt || bndswp) {
		jqsave = 0;
		nuncon = 0;
		bndswp = FALSE_;
		if (*gotr && ! posdef) {
		    s_wsfi(&io___94);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
		    s2trylu_(itn, &c__27, ns, &lureq, &luok, typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 44;
		    }
/* Ill-conditioned Z */
		    goto L100;
		}
	    }
	    needlm = statpt;
	    kprprt = kprc;
	    jq = 0;
	    if (needlm) {
/*              --------------------------------------------------------- */
/*              Compute Lagrange multipliers. */
/*              --------------------------------------------------------- */
		djq0 = djq;
/* save djq in case of bad statpt */
		djq = 0.;
		nuncon = 0;
		usegqp = feasbl && needf && gotgqp;
		weight = 0.;
		if (*elastc && feasbl) {
		    weight = *wtinf;
		}
		s5pric_(elastc, &feasbl, &incres, &usegqp, subopt, itn, m, n, 
			nb, ngqp0, ngqp, nnh, ns, &nfreez, &nonopt, &weight, &
			sgnobj, pinorm, &jq, &djq, &kprc, &rw[184], ne, nloca,
			 &loca[1], &inda[1], &acol[1], &hetype[1], &hs[1], &
			bl[1], &bu[1], &gqp[1], &pi[1], &rc[1], &x[1], &
			xfreez[1], &iw[1], leniw, &rw[1], lenrw);
		optiml = nonopt == 0;
		newsb = nonopt > 0;
	    }
/* needLM */
	}
/* jstPh1 */
	qpdone = optiml || deadpt || unbndd;
	if (qpdone) {
/*           ------------------------------------------------------------ */
/*           Apparently we are finished. */
/*           See if any nonbasics have to be set back on their bounds. */
/*           ------------------------------------------------------------ */
	    s5dgen_(&inform__, &c__1, prtlvl, nb, ninf, itn, &featol, tolx, &
		    tolinc, &hs[1], &bl[1], &bu[1], &x[1], &itnfix, nfix, &
		    tolx0, &iw[1], leniw, &rw[1], lenrw);
	    qpdone = inform__ == 0;
	    if (qpdone) {
/*              --------------------------------------------------------- */
/*              So far so good.  Now check the row residuals. */
/*              --------------------------------------------------------- */
		if (iw[215] > 0) {
		    s5setx_(&inform__, &c__1, itn, m, n, nb, &nbs, &rowerr, 
			    ne, nloca, &loca[1], &inda[1], &acol[1], &kbs[1], 
			    &xbs[1], nrhs0, nrhs, &rhs[1], &x[1], &y[1], &y1[
			    1], &iw[1], leniw, &rw[1], lenrw);
		    qpdone = inform__ == 0;
		    lureq = inform__;
		    if (lureq > 0) {
			*typelu = 2;
		    }
		}
	    }
	    if (qpdone) {
		if (unbndd) {
		    *iexit = -2;
		}
		if (deadpt) {
		    *iexit = -4;
		}
	    } else {
		*needx = TRUE_;
		unbndd = FALSE_;
		needpi = TRUE_;
		goto L100;
	    }
	}
/* done */
	if (jstph1 && optiml) {
/*           Relax, we are about to exit without printing anything. */
	} else {
/*           ============================================================ */
/*           Print the details of this iteration. */
/*           ============================================================ */
	    objprt = 0.;
	    if (feasbl) {
		if (needf) {
		    objprt = *objadd + objslk + *objqp;
		}
		if (needv) {
		    objprt += sgnobj * *wtinf * *sinf;
		}
	    }
	    (*qplog)(prob, probtag, elastc, gotr, &jstfea, &feasbl, m, mbs, 
		    nnh, ns, &jsq, &jbr, &jsr, &iw[220], &iw[221], itn, itqp, 
		    &kprprt, lvlinf, &pivot, &step, ninf, sinf, wtinf, &
		    objprt, &condhz, &djqprt, rgnorm, &kbs[1], &xbs[1], &iw[1]
		    , leniw, (ftnlen)20);
	}
	jbq = 0;
	jbr = 0;
	jsq = 0;
	jsr = 0;
	kprprt = 0;
	djqprt = 0.;
	if (qpdone) {
/*           ------------------------------------------------------------ */
/*           Convergence. */
/*           ------------------------------------------------------------ */
	    if (*ninf > 0) {
/* No feasible point. */
/* Stop or continue in elastic mode, depending on the */
/* specified level of infeasibility. */
		if (*lemode > 0) {
/* Enter elastic mode */
		    if (*elastc) {
/* Already in elastic mode, so we are done. */
			if (feasbl) {
/* Phase 2 elastic mode */
/* Find the final sumInf for the elastics */
			    s5ewrap_(&nbs, nb, ninf, sinf, &featol, &hestat[1]
				    , &kbs[1], &bl[1], &bu[1], &xbs[1]);
			} else {
/* Infeasible (this should not happen) */
/* The nonelastics cannot be satisfied */
/* by relaxing the elastic variables.  Exit. */
			    *iexit = -1;
/* Infeasible nonelastics */
			}
		    } else {
/* The constraints are infeasible in Normal mode. */
/* Print a message and start Elastic Phase 1. */
/* if .not. Elastc */
			if (prt10) {
			    s_wsfi(&io___102);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, probtag, (ftnlen)20);
			    e_wsfi();
			    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
			    s_wsfi(&io___103);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    e_wsfi();
			    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)132);
			    iw[223] = 1;
			    iw[225] = 1;
			}
			*elastc = TRUE_;
			qpdone = FALSE_;
			needf = *lvlinf != 2;
/* Need F1 in phase 2 */
			needv = *lvlinf != 0;
/* Need F2 in phase 2 */
			needpi = TRUE_;
			djq = 0.;
			step = 0.;
			s5einit_(nb, &nbs, &nelast, &featol, &infbnd, &hetype[
				1], &hestat[1], &kbs[1], &bl[1], &bu[1], &
				blbs[1], &bubs[1], &xbs[1]);
			gote = TRUE_;
		    }
		    goto L100;
		}
	    }
	    if (prt10 && ! jstph1) {
		if (jq != 0) {
		    djqprt = sgnobj * djq;
		    if (prt10 && klog == 1) {
			s_wsfi(&io___104);
			do_fio(&c__1, (char *)&djq, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&jq, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
		    }
		} else {
		    if (prt10 && klog == 1) {
			s_wsfi(&io___105);
			do_fio(&c__1, (char *)&(*rgnorm), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)132);
		    }
		}
	    }
	} else {
/* ---------------------------------------------------------- */
/* A nonbasic has been selected to become superbasic. */
/* Compute the vector y such that B y = column jq. */
/* ---------------------------------------------------------- */
/* not done */
	    if (newsb) {
/* -------------------------------------------------------- */
/* The price has selected a nonbasic to become superbasic. */
/* -------------------------------------------------------- */
		if (*ns + 1 > maxr) {
		    *iexit = -5;
		    goto L100;
		}
		lvltol = 2;
		djqprt = djq;
/*              --------------------------------------------------------- */
/*              Compute the vector pBS such that B pB = column jq. */
/*              pBS is a multiple of part of the new column of  Z  and */
/*              is used to define the QP search direction and update R. */
/*              --------------------------------------------------------- */
/*              Unpack column jq into  y1  and solve  B*y = y1. */
/*              The solve computes  y1  such that  L*y1 = ajq. */
/*              It is used below to modify L and U in s5QPit. */
		s2unpk_(&jq, m, n, ne, &anorm, nloca, &loca[1], &inda[1], &
			acol[1], &y1[1]);
		s2bsol_(iexit, &c__1, m, &y1[1], &pbs[1], &iw[1], leniw, &rw[
			1], lenrw);
		if (*iexit != 0) {
		    return 0;
		}
	    }
/*           ============================================================ */
/*           Take a step. */
/*           ============================================================ */
	    if (*itn >= itnlim || *itqp >= *itqpmax) {
		*iexit = -3;
		goto L100;
	    }
	    ++(*itqp);
	    ++(*itn);
/*           Decide if we want to print something this iteration. */
	    prtlog = *prtlvl >= 1 && *itqp % klog == 0;
	    prtsum = *prtlvl >= 1 && *itqp % ksumm == 0;
	    iw[218] = 0;
	    iw[219] = 0;
	    if (prtlog) {
		iw[218] = 1;
	    }
	    if (prtsum) {
		iw[219] = 1;
	    }
/*           ------------------------------------------------------------ */
/*           Take a reduced gradient step. */
/*           The new  x  will either minimize the objective on the */
/*           working set or lie on the boundary of a new constraint. */
/*           ------------------------------------------------------------ */
	    s5qpit_(&inform__, (U_fp)hprod, (U_fp)hprod1, &bndswp, elastc, &
		    feasbl, &gotgqp, &goth, gotr, &incres, &needf, &needv, &
		    needpi, &newsb, &posdef, itn, lenr, m, mbs, &maxr, maxs, 
		    n, nb, hvcalls, nnh0, nnh, ns, ngqp0, ngqp, ndegen, &
		    lureq, &kp, &jbq, &jsq, &jbr, &jsr, &jq, &jqsave, &nfmove,
		     &nuncon, &djq0, &djq, minimz, &obj, objqp, &featol, &
		    pivot, &step, &tolinc, wtinf, ne, nloca, &loca[1], &inda[
		    1], &acol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[1], &
		    kbs[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &gbs[1], &gqp[
		    1], &hdx[1], &pbs[1], &rg[1], &r__[1], &x[1], &xbs[1], &y[
		    1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
		    lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		    ftnlen)8, (ftnlen)8);
/*           Check for trouble in s5QPit. */
/*           inform values are -2, -1, 0 >0 */
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal LU error */
		} else if (inform__ == -1) {
		    *iexit = -2;
/* unbounded */
		} else if (inform__ == -2) {
		    *iexit = -6;
/* Hz not positive definite */
		}
		goto L100;
	    }
	    if (lureq > 0) {
		s2trylu_(itn, &lureq, ns, &lureq, &luok, typelu, &iw[1], 
			leniw, &rw[1], lenrw);
/* $$$               call s2tryLU2 */
/* $$$     &            ( itn, .false., LUreq, nS, LUreq, LUok, typeLU, */
/* $$$     &              iw, leniw, rw, lenrw ) */
		if (! luok) {
		    *iexit = 43;
		    goto L100;
		}
	    }
	    ++iw[215];
	    newlu = FALSE_;
	    newx = FALSE_;
	    chkpi = FALSE_;
/*           Increment featol every iteration. */
	    featol += tolinc;
/*           ============================================================ */
/*           Test for error condition and/or frequency interrupts. */
/*           ============================================================ */
/*           (1) Save a basis map (frequency controlled). */
/*           (2) Every kdegen iterations, reset featol and move nonbasic */
/*               variables onto their bounds if they are very close. */
/*           (3) Refactorize the basis if it has been modified too many */
/*               times. */
/*           (4) Update the LU factors of the basis if requested. */
/*           (5) Check row error (frequency controlled). */
	    if (*itn % ksav == 0) {
		s4ksav_(minimz, m, n, nb, ns, mbs, itn, ninf, sinf, objqp, &
			kbs[1], &hs[1], &ascale[1], &bl[1], &bu[1], &x[1], &
			xbs[1], cw + 8, lencw, &iw[1], leniw, (ftnlen)8);
	    }
	    if (*itn % kdegen == 0) {
		s5dgen_(&inform__, &c__2, prtlvl, nb, ninf, itn, &featol, 
			tolx, &tolinc, &hs[1], &bl[1], &bu[1], &x[1], &itnfix,
			 nfix, &tolx0, &iw[1], leniw, &rw[1], lenrw);
		*needx = inform__ > 0;
	    }
	    if (lureq == 0) {
		if (iw[216] >= kfac - 1) {
		    lureq = 1;
		} else if (iw[216] >= 20 && iw[173] + iw[174] > lumax) {
		    bgrwth = (doublereal) (iw[173] + iw[174]);
		    bold = (doublereal) lusiz0;
		    bgrwth /= bold;
		    if (prt10) {
			s_wsfi(&io___110);
			do_fio(&c__1, (char *)&bgrwth, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)132);
		    }
		    lureq = 2;
		} else {
		    checkx = iw[215] % kchk == 0;
		    if (checkx && ! (*needx)) {
			s5setx_(&inform__, &c__1, itn, m, n, nb, &nbs, &
				rowerr, ne, nloca, &loca[1], &inda[1], &acol[
				1], &kbs[1], &xbs[1], nrhs0, nrhs, &rhs[1], &
				x[1], &y[1], &y1[1], &iw[1], leniw, &rw[1], 
				lenrw);
			lureq = inform__;
		    }
		}
		if (lureq > 0) {
		    *typelu = 3;
		}
	    }
	}
/* not optiml */
	goto L100;
/* +    end while */
    }
/*     ======================end of main loop============================ */

    s5hs_(&c__1, nb, &bl[1], &bu[1], &hs[1], &x[1]);
    if (*subopt > 0) {
	if (nfreez > 0) {
/*           Relax */
	} else {
	    *subopt = 0;
	}
    }
    return 0;
/* L1500: */
} /* s5qp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QP */
/* Subroutine */ int s5chkp_(integer *iexit, integer *itn, integer *nbs, 
	integer *jqsave, integer *kbs, doublereal *gtp, doublereal *pbs, 
	integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Bad directional derivat"
	    "ive \002,1p,e9.1)";
    static char fmt_9000[] = "(\002 XXX  s5chkp.  kSave not found. jqSave ="
	    " \002,i5)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, jq;
    static char str[80];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ksave;
    static doublereal psave;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___118 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___119 = { 0, str, 0, fmt_9000, 80, 1 };


/*     ================================================================== */
/*     s5chkp  makes  pBS  a feasible direction. */

/*     16 Jun 1995: First version of s5chkp. */
/*     02 Aug 2003: snPRNT adopted. */
/*     02 Aug 2003: Current version of s5chkp. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pbs;
    --kbs;
    --iw;

    /* Function Body */
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Find the element of  pBS  corresponding to the most recently freed */
/*     variable. Usually, it will be pBS(nBS). */
/*     ------------------------------------------------------------------ */
    jq = abs(*jqsave);
    ksave = 0;
    for (k = *nbs; k >= 1; --k) {
	j = kbs[k];
	if (j == jq) {
	    ksave = k;
	    goto L100;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Choose the sign of  pBS  so that the most recently freed */
/*     variable continues to increase or decrease. */
/*     ------------------------------------------------------------------ */
L100:
    if (ksave > 0) {
	psave = pbs[ksave];
	if (*jqsave < 0 && psave > 0. || *jqsave > 0 && psave < 0.) {
	    dscal_(nbs, &c_b76, &pbs[1], &c__1);
	    *gtp = -(*gtp);
	}
	if (*gtp > 0.) {
/*           ------------------------------------------------------------ */
/*           Looks as though the sign of gtp cannot be relied upon. */
/*           In later versions we'll fix this variable. */
/*           For now, we just print a warning and stop. */
/*           ------------------------------------------------------------ */
	    s_wsfi(&io___118);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*gtp), (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
	    *iexit = 1;
/* Bad directional derivative */
	}
    } else {
/*        --------------------------------------------------------------- */
/*        Couldn't find the index of the most recently freed variable. */
/*        This should never happen! */
/*        --------------------------------------------------------------- */
	s_wsfi(&io___119);
	do_fio(&c__1, (char *)&(*jqsave), (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
    }
    return 0;
} /* s5chkp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5chkp */
/* Subroutine */ int s5chzq_(integer *m, integer *mbs, integer *n, integer *
	nb, integer *ns, integer *kbsq, doublereal *pivot, doublereal *tolpiv,
	 integer *ne, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *kbs, doublereal *bl, doublereal *bu, 
	doublereal *xbs, doublereal *y, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 XXX  s5chzq.  Max pivot is too small:"
	    "\002,1p,e11.1)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k;
    static doublereal d1, d2;
    static integer m1;
    static doublereal xj, tol;
    static char str[80];
    static doublereal eps0, dpiv;
    extern /* Subroutine */ int s2bprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___123 = { 0, str, 0, fmt_1000, 80, 1 };


/*     ================================================================== */
/*     s5chzq  selects a superbasic to replace the kp-th basic variable. */
/*     On entry,  y  contains the kp-th row of B(inverse). */
/*     On exit, pivot and  y(m+1), ..., y(m+nS) define the S-part of */
/*     the modifying vector w. */

/*     01 Dec 1991: First version based on Minos routine m7chzq. */
/*     02 Aug 2003: snPRNT adopted. */
/*     30 Jun 2005: Current version of s5chzq. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y;
    --xbs;
    --kbs;
    --bu;
    --bl;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/*     Set yS = 0 -  S'*y. */
/* eps**(4/5) */
    m1 = *m + 1;
    s2bprd_(&c__1, &eps0, n, ns, &kbs[m1], ne, nloca, &loca[1], &inda[1], &
	    acol[1], &c_b76, &y[1], m, &c_b5, &y[m1], ns);
    *kbsq = *m + idamax_(ns, &y[m1], &c__1);
    *pivot = (d__1 = y[*kbsq], abs(d__1));
/*     Exit if the pivot is too small. */
    if (*pivot < *tolpiv) {
	s_wsfi(&io___123);
	do_fio(&c__1, (char *)&(*pivot), (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)80);
	*kbsq = -(*m + *ns);
    } else {
/*        Choose one away from its bounds if possible. */
	tol = *pivot * .1;
	dpiv = -1.;
	i__1 = *m + *ns;
	for (k = m1; k <= i__1; ++k) {
	    if ((d__1 = y[k], abs(d__1)) >= tol) {
		j = kbs[k];
		xj = xbs[k];
		d1 = xj - bl[j];
		d2 = bu[j] - xj;
/* Computing MIN */
		d__1 = abs(d1), d__2 = abs(d2);
		d1 = min(d__1,d__2);
		if (dpiv <= d1) {
		    dpiv = d1;
		    *kbsq = k;
		}
	    }
	}
	*pivot = -y[*kbsq];
    }
/* pivot .ge. tolpiv */
    return 0;
} /* s5chzq_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5chzq */
/* Subroutine */ int s5getp_(logical *feasbl, logical *gotr, logical *newsb, 
	logical *posdef, integer *maxr, integer *lenr, integer *n, doublereal 
	*djq, doublereal *r__, doublereal *g, doublereal *p, doublereal *gp, 
	doublereal *php)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal diag;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ldiag;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), s6rsol_(integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal dirctn;

/*     ================================================================== */
/*     s5getp  computes a search direction  p  for the superbasic */
/*     variables, using the current reduced gradient  g. */

/*     29 Mar 2001: R stored by rows. */
/*     20 May 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --p;
    --g;

    /* Function Body */
    if (*feasbl && *gotr) {
	if (*posdef) {
/* Compute the Newton direction. */
	    dcopy_(n, &g[1], &c__1, &p[1], &c__1);
	    s6rsol_(&c__1, maxr, n, lenr, &r__[1], &p[1]);
	    *php = ddot_(n, &p[1], &c__1, &p[1], &c__1);
	    s6rsol_(&c__0, maxr, n, lenr, &r__[1], &p[1]);
	} else {
/* Compute a direction of zero or negative curvature. */
	    ldiag = (*n - 1) * *maxr + (3 - *n) * *n / 2;
/* Magic formula! */
	    diag = r__[ldiag];
/* Save last diag of R. */
	    r__[ldiag] = 1.;
/* Computing 2nd power */
	    d__1 = diag;
	    *php = d__1 * d__1;
	    dload_(n, &c_b5, &p[1], &c__1);
	    if (*newsb) {
		dirctn = g[*n];
	    } else {
		dirctn = *djq;
	    }
	    if (dirctn >= 0.) {
		p[*n] = 1.;
	    } else {
		p[*n] = -1.;
	    }
	    s6rsol_(&c__0, maxr, n, lenr, &r__[1], &p[1]);
	    r__[ldiag] = diag;
/* Restore diag of R. */
	}
    } else {
/* -------------------------------------------------------------- */
/* Direction of steepest-descent. */
/* -------------------------------------------------------------- */
	dcopy_(n, &g[1], &c__1, &p[1], &c__1);
	*php = 0.;
    }
/* ----------------------------------------------------------------- */
/* Fix the sign of p. */
/* ---------------------------------------------------------------- */
/* feasbl and gotR */
    dscal_(n, &c_b76, &p[1], &c__1);
    *gp = ddot_(n, &g[1], &c__1, &p[1], &c__1);
    return 0;
} /* s5getp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5getp */
/* Subroutine */ int s5getr_(integer *iexit, U_fp hprod, U_fp hprod1, integer 
	*hvcalls, logical *gotr, integer *typelu, integer *lureq, integer *
	itn, integer *lenr, integer *m, integer *mbs, integer *n, integer *nb,
	 integer *nnh, integer *ns, integer *prtlvl, integer *minimz, integer 
	*iobj, doublereal *targth, doublereal *targtz, integer *ne, integer *
	nloca, integer *loca, integer *inda, doublereal *acol, integer *hs, 
	integer *kbs, doublereal *bl, doublereal *bu, doublereal *blbs, 
	doublereal *bubs, doublereal *r__, integer *nrhs0, integer *nrhs, 
	doublereal *rhs, doublereal *x, doublereal *xbs, integer *iy, integer 
	*iy1, doublereal *y, doublereal *y1, doublereal *y2, char *cu, 
	integer *lencu, integer *iu, integer *leniu, doublereal *ru, integer *
	lenru, char *cw, integer *lencw, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1200[] = "(\002 Itn\002,i7,\002: Indefinite reduced Hess"
	    "ian\002)";
    static char fmt_1300[] = "(\002 Itn\002,i7,\002: Ill-conditioned Z.  Con"
	    "d(Z) = \002,1p,e8.2)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
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
	     integer *, doublereal *, integer *), s5rchk_(integer *, U_fp, 
	    U_fp, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, char *,
	     integer *, integer *, integer *, doublereal *, integer *, char *,
	     integer *, integer *, integer *, doublereal *, integer *, ftnlen,
	     ftnlen);
    static logical rcheck, needlu;
    static integer inform__;
    static doublereal plinfy;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2trylu_(integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, integer *, doublereal *, integer 
	    *);

    /* Fortran I/O blocks */
    static icilist io___147 = { 0, str, 0, fmt_1200, 80, 1 };
    static icilist io___149 = { 0, str, 0, fmt_1300, 80, 1 };


/*     ================================================================== */
/*     s5getR   computes the Cholesky factor of the reduced Hessian. */

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
/*     LUreq = 25  Ill-conditioned null-space basis in QP */
/*     LUreq = 26  Indefinite Z'HZ in QP */
/*     LUreq = 27  R singular after bound swap in QP */

/*     25 Oct 2003: First version of s5getR based on s8getR. */
/*     26 Dec 2003: s2newLU added. */
/*     09 Dec 2004: Changed to column-packed format for H. */
/*     09 Dec 2004: Current version of s5getR. */
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
    --x;
    --bu;
    --bl;
    --hs;
    --acol;
    --inda;
    --loca;
    --rhs;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    maxr = iw[52];
/* max columns of R */
    eigh = iw[200];
/* =1(0) for definite QP Hessian */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    flmax = rw[8];
/* est. of the largest pos. real */
    *iexit = 0;
    plinfy = flmax;
    luok = TRUE_;
/*     ================================================================== */
/* +    while (LUok  .and. .not. gotR) do */
L100:
    if (luok && ! (*gotr)) {
/*     ------------------------------------------------------------------ */
	inform__ = 0;
	needlu = *lureq > 0;
	if (needlu) {
	    s2bfac_(iexit, typelu, &needlu, &newlu, &newb, iobj, itn, prtlvl, 
		    lureq, m, mbs, n, nb, nnh, ns, &nswap, ne, nloca, &loca[1]
		    , &inda[1], &acol[1], &kbs[1], &hs[1], &bl[1], &bu[1], &
		    blbs[1], &bubs[1], nrhs0, nrhs, &rhs[1], &x[1], &xbs[1], &
		    iy[1], &iy1[1], &y[1], &y2[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    if (*iexit != 0) {
		goto L900;
	    }
	}
	if (*ns > 0) {
/*           ------------------------------------------------------------ */
/*           Compute and factorize  Z'HZ. */
/*           ------------------------------------------------------------ */
	    s5hz_(&inform__, (U_fp)hprod, (U_fp)hprod1, &maxr, lenr, minimz, 
		    m, mbs, n, nb, nnh, ns, hvcalls, ne, nloca, &loca[1], &
		    inda[1], &acol[1], &hdmax, &rw[192], targtz, &kbs[1], &
		    r__[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[1], 
		    leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1]
		    , lenrw, (ftnlen)8, (ftnlen)8);
/*           Check for trouble in s5Hz.  Possible exits are: */
/*           inform  Status */
/*           ------  ------ */
/*            -1     Ill-conditioned null-space basis */
/*             0     reduced Hessian computed successfully */
/*            >0     Fatal error in LU solve */
	    if (inform__ == 0) {
		s5hfac_(&inform__, &eigh, itn, lenr, m, &maxr, mbs, nb, ns, 
			targth, &hdmax, &hs[1], &kbs[1], &iy[1], &bl[1], &bu[
			1], &blbs[1], &bubs[1], &x[1], &xbs[1], &r__[1], &iw[
			1], leniw, &rw[1], lenrw);
/*              Check for trouble in s5Hfac. */
/*              inform    Status */
/*              ------    ------ */
/*                -2      H singular (but should be positive definite) */
/*                -1      H indefinite */
/*                 0      normal exit */
		if (inform__ != 0) {
/*                 The reduced Hessian appears to be indefinite. */
/*                 Refactorize B with reduced factor tol. */
/*                 If the factor tol is already tight, give up. */
		    s_wsfi(&io___147);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    s2trylu_(itn, &c__26, ns, lureq, &luok, typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		}
		rcheck = FALSE_;
		if (rcheck) {
		    s5rchk_(iexit, (U_fp)hprod, (U_fp)hprod1, itn, minimz, &
			    maxr, lenr, m, mbs, n, nb, hvcalls, nnh, ns, ne, 
			    nloca, &loca[1], &inda[1], &acol[1], &kbs[1], &
			    r__[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, &iu[
			    1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
			    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		    if (*iexit != 0) {
			goto L900;
		    }
		}
	    } else if (inform__ == -1) {
/*              Ill-conditioned Z from s5Hz. */
/*              Refactorize B, possibly with a reduced factor tol. */
/*              If factor tol is already tight, accept Z, however bad. */
		s_wsfi(&io___149);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rw[192], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		s2trylu_(itn, &c__25, ns, lureq, &luok, typelu, &iw[1], leniw,
			 &rw[1], lenrw);
		if (! luok) {
		    *targtz = plinfy;
		    luok = TRUE_;
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
/* indefinite reduced Hessian */
	*iexit = 94;
    }
L900:
    return 0;
/* L1100: */
} /* s5getr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5getR */
/* Subroutine */ int s5hfac_(integer *iexit, integer *eigh, integer *itn, 
	integer *lenr, integer *m, integer *maxr, integer *mbs, integer *nb, 
	integer *ns, doublereal *hcndbd, doublereal *hdmax, integer *hs, 
	integer *kbs, integer *perm, doublereal *bl, doublereal *bu, 
	doublereal *blbs, doublereal *bubs, doublereal *x, doublereal *xbs, 
	doublereal *r__, integer *iw, integer *leniw, doublereal *rw, integer 
	*lenrw)
{
    /* Format strings */
    static char fmt_9060[] = "(\002 Itn\002,i7,\002: Reduced Hessian appears"
	    " to be indefinite.\002,\002 dpiv, Hdmin = \002,1p,e9.2,\002,\002"
	    ",e9.2)";
    static char fmt_9000[] = "(\002 Itn\002,i7,\002: Reduced Hessian appears"
	    " to have \002,i6,\002 small eigenvalues.  PD tol = \002,1p,e9.2)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static integer js;
    static doublereal eps;
    static char str[90];
    static integer jmax, kmax;
    static doublereal dpiv, hdmin;
    static integer ksave, pivot;
    extern /* Subroutine */ int s6chol_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);
    static integer inform__, rankhz, nssave;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___163 = { 0, str, 0, fmt_9060, 90, 1 };
    static icilist io___164 = { 0, str, 0, fmt_9000, 90, 1 };


/*     ================================================================== */
/*     s5Hfac  factorizes the reduced Hessian Z'HZ. */

/*      iExit       Result */
/*      -----       ------ */
/*       -2         H  singular   (but should be positive definite) */
/*       -1         H  indefinite */
/*        0         positive-definite, factors computed successfully */

/*     13 Oct 1992: First version based on Qpsol routine Qpcrsh. */
/*     15 Oct 1994: Dependent columns fixed at their current value. */
/*     02 Aug 2003: snPRNT adopted. */
/*     01 Sep 2005: Current version of s5Hfac. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --perm;
    --xbs;
    --bubs;
    --blbs;
    --kbs;
    --x;
    --bu;
    --bl;
    --hs;
    --iw;
    --rw;

    /* Function Body */
    *iexit = 0;
    eps = max(*hdmax,1.) / *hcndbd;
/* Computing MAX */
    d__1 = *hdmax / *hcndbd;
    hdmin = max(d__1,eps);
    if (*eigh == 1) {
	pivot = 0;
    } else if (*eigh == -1 || *eigh == 0) {
	pivot = 1;
    } else {
/* Relax -- there are no other options. */
    }
    s6chol_(&inform__, &pivot, maxr, ns, lenr, &r__[1], &hdmin, &dpiv, &
	    rankhz, &perm[1]);
/*     inform > 0 implies rankHz < nS. */
    if (pivot == 1) {
/*        ----------------------- */
/*        Apply any interchanges. */
/*        ----------------------- */
	i__1 = min(rankhz,*ns);
	for (j = 1; j <= i__1; ++j) {
	    jmax = perm[j];
	    if (jmax > j) {
		kmax = *m + jmax;
		k = *m + j;
		ksave = kbs[kmax];
		kbs[kmax] = kbs[k];
		kbs[k] = ksave;
		s = xbs[kmax];
		xbs[kmax] = xbs[k];
		xbs[k] = s;
		s = blbs[kmax];
		blbs[kmax] = blbs[k];
		blbs[k] = s;
		s = bubs[kmax];
		bubs[kmax] = bubs[k];
		bubs[k] = s;
	    }
	}
    }
    if (dpiv < -hdmin) {
/*        --------------------------------------- */
/*        H  appears to be indefinite. */
/*        --------------------------------------- */
	*iexit = -1;
/* Indefinite H */
    } else if (dpiv < hdmin) {
/*        --------------------------------------- */
/*        H  appears to be positive semidefinite. */
/*        rankHz < nS */
/*        --------------------------------------- */
	if (*eigh == 1) {
/* H should be PD */
	    *iexit = -2;
	    s_wsfi(&io___163);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dpiv, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&hdmin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)90);
	} else {
/* singular H allowed */
	    s_wsfi(&io___164);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    i__1 = *ns - rankhz;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&hdmin, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)90);
	    nssave = *ns;
	    i__1 = nssave;
	    for (js = rankhz + 1; js <= i__1; ++js) {
		k = *m + js;
		j = kbs[k];
/*              Make variable  j  nonbasic (it is already feasible). */
/*              hs(j) = -1 means x(j) is strictly between its bounds. */
		if (x[j] <= bl[j]) {
		    x[j] = bl[j];
		    hs[j] = 0;
		} else if (x[j] >= bu[j]) {
		    x[j] = bu[j];
		    hs[j] = 1;
		} else {
		    hs[j] = -1;
		}
		if (bl[j] == bu[j]) {
		    hs[j] = 4;
		}
		--(*ns);
	    }
	    *ns = min(*ns,rankhz);
	}
/* eigH == POSDEF */
    }
    return 0;
} /* s5hfac_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Hfac */
/* Subroutine */ int s5hz_(integer *iexit, S_fp hprod, U_fp hprod1, integer *
	maxr, integer *lenr, integer *minimz, integer *m, integer *mbs, 
	integer *n, integer *nb, integer *nnh, integer *ns, integer *hvcalls, 
	integer *ne, integer *nloca, integer *loca, integer *inda, doublereal 
	*acol, doublereal *hdmax, doublereal *condz, doublereal *targtz, 
	integer *kbs, doublereal *r__, doublereal *v, doublereal *w, 
	doublereal *y, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer jq, js, nbs;
    static doublereal eps0, diag;
    static integer ldiag;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm;
    extern /* Subroutine */ int s2bprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), s2bsol_(integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s2unpk_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), s6rrow_(integer *, integer *, integer *, integer *,
	     doublereal *, doublereal *, integer *);
    static doublereal sgnobj;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer sqstat;
    extern /* Subroutine */ int s2gathr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), s2scatr_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Hz    computes the reduced Hessian and loads it by columns into */
/*     the upper triangle R. */

/*      iExit       Status */
/*      -----       ------ */
/*       -1         Ill-conditioned null-space basis */
/*        0         reduced Hessian computed successfully */
/*       >0         Fatal error in LU solve */

/*     13 Oct 1992: First version based on QPSOL routine Qpcrsh. */
/*     15 Oct 1994: Dependent columns fixed at their current value. */
/*     04 Dec 2000: R converted to row-wise storage. */
/*     25 Mar 2005: Current version */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --kbs;
    --y;
    --w;
    --v;
    --acol;
    --inda;
    --loca;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    if (*ns == 0) {
	return 0;
    }
    *iexit = 0;
    nbs = *m + *ns;
    *hdmax = 0.;
    sgnobj = (doublereal) (*minimz);
    *condz = 1.;
    sqstat = 0;
/*     ------------------------------------------------------------------ */
/*     Main loop to find a column of Z'HZ. */
/*     ------------------------------------------------------------------ */
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
/* -------------------------------------------------------------- */
/* Get the nonlinear elements of the column of Z. */
/* Find y such that B y = column jq. */
/* Scatter the nonlinear part of y into w. */
/* -------------------------------------------------------------- */
	jq = kbs[*m + js];
	s2unpk_(&jq, m, n, ne, &anorm, nloca, &loca[1], &inda[1], &acol[1], &
		w[1]);
	s2bsol_(iexit, &c__1, m, &w[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit > 0) {
	    return 0;
	}
/* Computing MAX */
	d__1 = dnormi_(m, &y[1], &c__1) / anorm;
	*condz = max(d__1,*condz);
	if (*condz >= *targtz) {
	    *iexit = -1;
	    return 0;
	}
	s2scatr_(nnh, m, &kbs[1], &c_b76, &y[1], &w[1]);
	if (jq <= *nnh) {
	    w[jq] = 1.;
	}
/* Set v = H w. */
	if (*nnh > 0) {
	    (*hprod)((U_fp)hprod1, nnh, &w[1], &v[1], &sqstat, cu + 8, lencu, 
		    &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	    ++(*hvcalls);
	    if (*minimz < 0) {
		dscal_(nnh, &sgnobj, &v[1], &c__1);
	    }
	}
/* -------------------------------------------------------------- */
/* Gather w = vBS and compute v = Z' w. */
/* Solve  B' vB = wB  and  form  wS = wS - S' vB. */
/* -------------------------------------------------------------- */
	s2gathr_(nnh, &nbs, &kbs[1], &c_b146, &v[1], &w[1]);
	s2bsol_(iexit, &c__2, m, &w[1], &v[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit > 0) {
	    return 0;
	}
	s2bprd_(&c__1, &eps0, n, ns, &kbs[*m + 1], ne, nloca, &loca[1], &inda[
		1], &acol[1], &c_b76, &v[1], m, &c_b146, &w[*m + 1], ns);
/* -------------------------------------------------------------- */
/* Store w(1:nS) in the jS-th row of R. */
/* R is NO LONGER SYMMETRIZED. */
/* -------------------------------------------------------------- */
	s6rrow_(&js, maxr, ns, lenr, &r__[1], &w[*m + 1], &ldiag);
	diag = r__[ldiag];
/* Computing MAX */
	d__1 = *hdmax, d__2 = abs(diag);
	*hdmax = max(d__1,d__2);
    }
    return 0;
} /* s5hz_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Hz */
/* Subroutine */ int s5qpfg_(S_fp hprod, U_fp hprod1, integer *ngqp, integer *
	ngobj0, integer *ngobj, integer *nnh, integer *sqstat, integer *
	hvcalls, doublereal *fqp, doublereal *gobj, doublereal *gqp, integer *
	lenx0, integer *nx0, doublereal *x0, doublereal *x, doublereal *dx, 
	char *cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru,
	 integer *lenru, char *cw, integer *lencw, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer nzero;

/*     ================================================================== */
/*     s5QPfg  computes various quantities associated with the LP/QP. */

/*       1.  fQP =  gObj'*(x-x0)  + half*(x - x0)'*H*(x - x0) */
/*       2.  gQP =  gradient of fQP */

/*     On entry, */
/*     ngQP         is max( ngObj, nnH ) */
/*     x(ngQP)      are the nonlinear variables */
/*     x0(ngQP)     is the base point x0 */
/*     gObj(ngObj)  defines the explicit QP linear term */

/*     On exit, */
/*     fQP          is the QP quadratic term (1) above */
/*     gQP(ngQP)    is the gradient of fQP */
/*     dx(ngQP)     is  x-x0 */

/*     02 May 1992: First version of s5QPfg. */
/*     23 Oct 1993: Hx added as an argument. */
/*     29 Oct 1993: Modified to compute only the QP objective. */
/*     07 Oct 1994: gQP added as an argument. */
/*     09 Dec 2004: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --dx;
    --x;
    --gqp;
    --gobj;
    --x0;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    if (*ngqp <= 0) {
	return 0;
    }
    dcopy_(ngqp, &x[1], &c__1, &dx[1], &c__1);
    if (*nx0 > 0) {
	daxpy_(ngqp, &c_b76, &x0[1], &c__1, &dx[1], &c__1);
    }
    *fqp = 0.;
    if (*nnh > 0) {
	(*hprod)((U_fp)hprod1, nnh, &dx[1], &gqp[1], sqstat, cu + 8, lencu, &
		iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &
		rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	++(*hvcalls);
	*fqp = ddot_(nnh, &dx[1], &c__1, &gqp[1], &c__1) * .5;
    }
    nzero = *ngqp - *nnh;
    if (nzero > 0) {
	dload_(&nzero, &c_b5, &gqp[*nnh + 1], &c__1);
    }
    if (*ngobj > 0) {
	*fqp += ddot_(ngobj, &gobj[1], &c__1, &dx[1], &c__1);
	daxpy_(ngobj, &c_b146, &gobj[1], &c__1, &gqp[1], &c__1);
    }
    return 0;
} /* s5qpfg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QPfg */
/* Subroutine */ int s5qpit_(integer *iexit, S_fp hprod, U_fp hprod1, logical 
	*bndswp, logical *elastc, logical *feasbl, logical *gotgqp, logical *
	goth, logical *gotr, logical *incres, logical *needf, logical *needv, 
	logical *needpi, logical *newsb, logical *posdef, integer *itn, 
	integer *lenr, integer *m, integer *mbs, integer *maxr, integer *maxs,
	 integer *n, integer *nb, integer *hvcalls, integer *nnh0, integer *
	nnh, integer *ns, integer *ngqp0, integer *ngqp, integer *ndegen, 
	integer *lureq, integer *kp, integer *jbq, integer *jsq, integer *jbr,
	 integer *jsr, integer *jq, integer *jqsave, integer *nfmove, integer 
	*nuncon, doublereal *djq0, doublereal *djq, integer *minimz, 
	doublereal *obj, doublereal *objqp, doublereal *featol, doublereal *
	pivot, doublereal *step, doublereal *tolinc, doublereal *wtinf, 
	integer *ne, integer *nloca, integer *loca, integer *inda, doublereal 
	*acol, integer *hetype, integer *hestat, integer *hfeas, integer *hs, 
	integer *kbs, doublereal *bl, doublereal *bu, doublereal *blbs, 
	doublereal *bubs, doublereal *gbs, doublereal *gqp, doublereal *hdx, 
	doublereal *pbs, doublereal *rg, doublereal *r__, doublereal *x, 
	doublereal *xbs, doublereal *y, doublereal *y1, doublereal *y2, char *
	cu, integer *lencu, integer *iu, integer *leniu, doublereal *ru, 
	integer *lenru, char *cw, integer *lencw, integer *iw, integer *leniw,
	 doublereal *rw, integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Bad direction after add"
	    "ing a superbasic.\002)";
    static char fmt_9999[] = "(\002 Itn\002,i7,\002: Chzq failed in s5QPit!"
	    "!\002)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer jr, js;
    static doublereal gp;
    static integer ns1, nbs;
    static doublereal eps, php;
    static integer ksq;
    static char str[80];
    static integer nbs1;
    static doublereal pbs1, eps0;
    extern /* Subroutine */ int s5zp_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *);
    static integer eigh;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jesq, kbsq;
    static logical move;
    static doublereal drsq, gpqp, tolp;
    static integer ntry;
    static doublereal tolp0;
    extern /* Subroutine */ int s5bsx_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static integer ldiag;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *);
    static doublereal bigdx;
    static logical onbnd;
    static doublereal exact, bound, anorm;
    static logical uncon;
    static doublereal stepb;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal phpqp, pnorm, stepp;
    extern /* Subroutine */ int s5chkp_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), s5rchk_(integer *, S_fp, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen), s5sdel_(integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s6rdel_(
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), s2bsol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s5rcol_(integer *, S_fp, U_fp, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     char *, integer *, integer *, integer *, doublereal *, integer *,
	     char *, integer *, integer *, integer *, doublereal *, integer *,
	     ftnlen, ftnlen), s5getp_(logical *, logical *, logical *, 
	    logical *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s5chzq_(integer *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), s5rsng_(integer *, logical *
	    , logical *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s2unpk_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), s5step_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , logical *, logical *, logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal objchg;
    static logical rcheck;
    static doublereal infbnd;
    extern /* Subroutine */ int s6rswp_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , s2bmod2_(integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, doublereal *, integer *);
    static logical unbndd;
    static doublereal sgnobj;
    static logical hitcon, hposdf, hsingr, singlr;
    static integer inform__, infpiv;
    static doublereal sclpiv;
    static logical hitlow;
    static integer jqstat, jrstat;
    static doublereal tolpiv;
    static integer sqstat;
    static doublereal stepmx;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2scatr_(integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static icilist io___203 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___226 = { 0, str, 0, fmt_9999, 80, 1 };


/*     ================================================================== */
/*     s5QPit performs a QP step. */

/*     On entry, */
/*        newSB = true implies that variable jq just went superbasic. */
/*                In this case: */
/*                pBS  satisfies B pBS = a(jq). */
/*                y1   satisfies L  y1 = a(jq). */

/*     On exit, */
/*        pBS contains the most recent QP search direction. */

/*      iExit       Result */
/*      -----       ------ */
/*       -2         reduced Hessian is not positive semidefinite */
/*       -1         unbounded */
/*        0         normal exit */
/*       >0         Fatal LU error */

/*     25 Nov 1991: First version of s5QPit. */
/*     05 Jan 1996: Positive semidefinite R treated correctly. */
/*     29 Aug 1996: First min sum version added. */
/*     27 Jul 1997: Thread-safe version. */
/*     02 Feb 1998: Piecewise linear line search added. */
/*     23 Mar 2000: gQP  and  H  scaled. */
/*     16 Oct 2000: Reverted to non-bordered version of s5QPit. */
/*     04 Dec 2000: R converted to row-wise storage. */
/*     02 Aug 2003: snPRNT adopted. */
/*     07 May 2006: s5Zp added to compute Z*p. */
/*     08 Apr 2008: hEstat accessed in elastic mode only. */
/*     04 Jul 2008: Modify both bl and bu in elastic phase 1. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* number of LU mods */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --xbs;
    --pbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --rg;
    --y2;
    --y1;
    --y;
    --x;
    --bu;
    --bl;
    --hs;
    --hestat;
    --hetype;
    --hdx;
    --gqp;
    --acol;
    --inda;
    --loca;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps = rw[1];
/* machine precision.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5)          IEEE DP  3.00e-13 */
    tolpiv = rw[60];
/* excludes small elements of pBS. */
    infbnd = rw[70];
/* definition of an infinite bound */
    bigdx = rw[72];
/* unbounded step. */
    eigh = iw[200];
/* -1,0,1 for indef, psd and pdef QP Hessian */
    *iexit = 0;
    sqstat = 0;
    unbndd = FALSE_;
    sgnobj = (doublereal) (*minimz);
    nbs = *m + *ns;
    if (*newsb) {
/*        --------------------------------------------------------------- */
/*        New superbasic. */
/*        PosDef must be true if there is a new superbasic. */
/*        --------------------------------------------------------------- */
	ns1 = *ns + 1;
	nbs1 = nbs + 1;
	kbs[nbs1] = *jq;
	xbs[nbs1] = x[*jq];
	blbs[nbs1] = bl[*jq];
	bubs[nbs1] = bu[*jq];
	jqstat = hs[*jq];
	*posdef = FALSE_;
	if (*gotr) {
/*           ------------------------------------------------------------ */
/*           Add the new column to R at position nS+1. */
/*           Check for a singular or indefinite reduced Hessian. */
/*           ------------------------------------------------------------ */
	    kbs[nbs1] = *jq;
	    s5rcol_(iexit, (S_fp)hprod, (U_fp)hprod1, minimz, jq, &ns1, &drsq,
		     &ldiag, maxr, lenr, m, mbs, n, nb, hvcalls, nnh, &ns1, 
		    ne, nloca, &loca[1], &inda[1], &acol[1], &kbs[1], &r__[1],
		     &y[1], &y2[1], &pbs[1], cu + 8, lencu, &iu[1], leniu, &
		    ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
		     (ftnlen)8, (ftnlen)8);
	    if (*iexit != 0) {
		goto L900;
	    }
	    s5rsng_(&eigh, &hposdf, &hsingr, itn, maxr, lenr, &ns1, &drsq, &
		    r__[1], &iw[1], leniw, &rw[1], lenrw);
	    if (*feasbl) {
		*posdef = hposdf;
		singlr = hsingr;
		if (! (*posdef || singlr)) {
		    *iexit = -2;
		    goto L900;
		}
	    } else {
		*gotr = hposdf;
	    }
	    if (*gotr) {
		r__[ldiag] = sqrt(drsq);
	    }
	}
/* -->     R can be checked here. */
/* gotR */
	rcheck = FALSE_;
	if (rcheck && *gotr) {
	    s5rchk_(iexit, (S_fp)hprod, (U_fp)hprod1, itn, minimz, maxr, lenr,
		     m, mbs, n, nb, hvcalls, nnh, &ns1, ne, nloca, &loca[1], &
		    inda[1], &acol[1], &kbs[1], &r__[1], &y[1], &y2[1], &pbs[
		    1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, 
		    lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
		    ;
	    if (*iexit != 0) {
		goto L900;
	    }
	}
	if (*incres) {
	    *jqsave = *jq;
	} else {
	    *jqsave = -(*jq);
	}
	*ns = ns1;
	nbs = nbs1;
	hfeas[nbs] = 0;
	if (*feasbl && *needf && *gotgqp && *jq <= *ngqp) {
	    gbs[nbs] = sgnobj * gqp[*jq];
	} else {
	    gbs[nbs] = 0.;
	}
/*        =============================================================== */
/*        Set hEstat(jq) and the elastic parts of blBS and buBS. */
/*        =============================================================== */
	if (*elastc) {
/*           If the new superbasic is an elastic variable */
/*           and it wants to move infeasible, set its elastic state. */
	    jesq = hestat[*jq];
	    if (hetype[*jq] > 0) {
		if (*incres) {
		    if (jqstat == 1 || jqstat == 4) {
			hestat[*jq] = 2;
			blbs[nbs] = bu[*jq];
			bubs[nbs] = infbnd;
			if (*feasbl) {
			    gbs[nbs] += *wtinf;
			}
		    }
		} else {
		    if (jqstat == 0 || jqstat == 4) {
			hestat[*jq] = 1;
			blbs[nbs] = -infbnd;
			bubs[nbs] = bl[*jq];
			if (*feasbl) {
			    gbs[nbs] -= *wtinf;
			}
		    }
		}
	    }
	}
/*        --------------------------------------------------------------- */
/*        In phase 1, or phase 2 for an LP, price can select nonbasics */
/*        floating free between their bounds with zero reduced cost. */
/*        We have to check that dqj is not zero. */
/*        --------------------------------------------------------------- */
/* Elastc */
	rg[*ns] = *djq;
	if (! (*feasbl) || *needf && ! (*goth)) {
	    if (hs[*jq] == -1) {
		if (*incres) {
		    rg[*ns] = -1.;
		} else {
		    rg[*ns] = 1.;
		}
	    }
	}
	*jsq = *jq;
	hs[*jq] = 2;
    }
/*     ------------------------------------------------------------------ */
/*     Store the free components of the search direction in pBS(1:nBS). */
/*     First, find the search direction pS for the superbasics, store it */
/*     in  pBS(m+1:nBS), and find its norm.  Put the search direction for */
/*     the basic variables in pBS(1)  ,...,pBS(m). */
/*     ------------------------------------------------------------------ */
/* newSB */
L100:
    singlr = ! (*posdef);
    s5getp_(feasbl, gotr, newsb, posdef, maxr, lenr, ns, djq, &r__[1], &rg[1],
	     &pbs[*m + 1], &gp, &php);
    pbs1 = pbs[*m + *ns];
    s5zp_(iexit, m, mbs, n, nb, ns, &eps0, &pnorm, ne, nloca, &loca[1], &inda[
	    1], &acol[1], &kbs[1], &pbs[1], &y2[1], &iw[1], leniw, &rw[1], 
	    lenrw);
    if (*iexit != 0) {
	goto L900;
    }
    if (*feasbl) {
/* -------------------------------------------------------------- */
/* If R is singular, ensure that pBS is a feasible direction. */
/* A nonzero exit value of inform implies that the directional */
/* derivative is too small to be relied upon. */
/* -------------------------------------------------------------- */
	if (*gotr && singlr) {
	    s5chkp_(&inform__, itn, &nbs, jqsave, &kbs[1], &gp, &pbs[1], &iw[
		    1], leniw);
	    if (inform__ > 0) {
		*lureq = 24;
		goto L900;
	    }
	}
	if (*newsb && *gotr) {
/* Check for a feasible direction. */
/* A large  rgTol  may give a pBS(nBS) with the wrong sign. */
/* If so, continue minimizing with the old superbasic set. */
	    if (*djq * pbs1 > 0.) {
		s_wsfi(&io___203);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		--(*ns);
		--nbs;
		hs[*jq] = jqstat;
		if (*elastc) {
		    hestat[*jq] = jesq;
		}
		*jq = 0;
		*djq = *djq0;
		*jsq = -(*jsq);
		*jqsave = 0;
		*posdef = TRUE_;
		*newsb = FALSE_;
		goto L100;
	    }
	    *bndswp = FALSE_;
	}
/*        --------------------------------------------------------------- */
/*        Compute y = pBS(scattered) and Hdx(scattered). */
/*        The vector Hdx is used to update the objective and gradient of */
/*        the QP.  Form  gpQP  and  pHpQP  for the quadratic. */
/*        gp = gpQP - pBS(kObj) + terms from the elastic gradient. */
/*        --------------------------------------------------------------- */
	if (*needf && (*gotgqp || *goth)) {
	    s2scatr_(ngqp, &nbs, &kbs[1], &c_b146, &pbs[1], &y[1]);
	    if (*gotgqp) {
		gpqp = ddot_(ngqp, &gqp[1], &c__1, &y[1], &c__1);
	    }
	    if (*goth) {
		phpqp = 0.;
		(*hprod)((U_fp)hprod1, nnh, &y[1], &hdx[1], &sqstat, cu + 8, 
			lencu, &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &
			iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
		++(*hvcalls);
		phpqp += ddot_(nnh, &y[1], &c__1, &hdx[1], &c__1);
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*     Find the nearest constraint in direction  x + step*pBS (step > 0). */
/*     Exact  is the step that takes xBS(kp) exactly onto bound. */
/*     It may be positive or slightly negative. (Not defined if Unbndd.) */

/*     If onbnd  is true, step is a step that reaches a bound exactly. */
/*     xBS(kp) reaches the value bound.  If we take a constrained step, */
/*     bound is used to put the new nonbasic variable x(jr) exactly on */
/*     its bound. */

/*     If Unbndd is true, step = stepmx. */
/*     ------------------------------------------------------------------ */
/* feasbl */
    stepmx = bigdx / pnorm;
    sclpiv = 1.;
    tolp0 = tolpiv;
    tolp = tolpiv * pnorm;
    ntry = 0;
/* +    Repeat */
L200:
    tolp /= sclpiv;
    tolp0 /= sclpiv;
    s5step_(&nbs, ndegen, featol, &infbnd, &stepmx, tolinc, &tolp, &hfeas[1], 
	    &blbs[1], &bubs[1], &xbs[1], &pbs[1], &hitlow, &move, &onbnd, &
	    unbndd, &infpiv, kp, &bound, &exact, &stepb, &stepp);
/* Find if the step is constrained or unconstrained. */
/* If R has been flagged as singular, we double check by trying */
/* to compute the QP minimizer along pBS.  If the minimizer */
/* exists,  the singularity tolerance must be too large. */
    if (*feasbl) {
	if (*posdef) {
	    uncon = stepp > 1.;
	} else {
	    uncon = stepp * php > -gp;
	}
	unbndd = unbndd && ! uncon || stepmx <= 1.;
    } else {
	uncon = FALSE_;
    }
    sclpiv = 10.;
    ++ntry;
/* +    until    (infpiv .eq. 0 .and. (.not. Unbndd .or. feasbl) .or. */
/* +                ntry .ge. mtry) */
    if (! (infpiv == 0 && (! unbndd || *feasbl) || ntry >= 6)) {
	goto L200;
    }
    if (unbndd) {
	*iexit = -1;
	goto L900;
    }
    hitcon = ! uncon;
    *needpi = TRUE_;
    if (hitcon) {
	*nuncon = 0;
	*step = stepb;
    } else {
	++(*nuncon);
	*pivot = 0.;
	if (*posdef) {
	    *step = 1.;
	} else {
	    *step = -gp / php;
	    *posdef = TRUE_;
	}
    }
/* ----------------------------------------------------------------- */
/* Compute ObjChg, the change in ObjQP (minimized or maximized). */
/* Note: Obj = sgnObj*ObjQP */
/* pHp = sgnObj*pHpQP */
/* ----------------------------------------------------------------- */
    if (*feasbl) {
/* Computing 2nd power */
	d__1 = *step;
	objchg = *step * gp + php * .5 * (d__1 * d__1);
	*obj += objchg;
	if (*needf) {
	    if (*gotgqp) {
		*objqp += *step * gpqp;
	    }
	    if (*goth) {
/* Computing 2nd power */
		d__1 = *step;
		*objqp += phpqp * .5 * (d__1 * d__1);
		if (*step > 0.) {
		    daxpy_(nnh, step, &hdx[1], &c__1, &gqp[1], &c__1);
		}
	    }
	}
    }
    if (*feasbl && move) {
	++(*nfmove);
    }
/*     ------------------------------------------------------------------ */
/*     Update the basic variables xBS. */
/*     ------------------------------------------------------------------ */
    daxpy_(&nbs, step, &pbs[1], &c__1, &xbs[1], &c__1);
    s5bsx_(&c__1, &nbs, nb, &kbs[1], &x[1], &xbs[1]);
    if (hitcon) {
/*        =============================================================== */
/*        There is a blocking variable. */
/*        It could be a fixed variable, whose new state must be 4. */
/*        =============================================================== */
	*pivot = -pbs[*kp];
	jr = kbs[*kp];
	*bndswp = jr == abs(*jqsave);
/* $$$!        10 Mar 2004: Care is needed to prevent the */
/* $$$!        new nonbasic variable jr from ending up slightly inside */
/* $$$!        its bound.  EXPAND normally ensures that x(jr) will be */
/* $$$!        ON or slightly OUTSIDE its bound, but now we realise that */
/* $$$!        rounding error might make it slightly INSIDE. */
/* $$$ */
/* $$$         if (onbnd) then */
/* $$$            x(jr) = bound */
/* $$$         else if (hitlow) then */
/* $$$            x(jr) = min( x(jr), bl(jr) ) */
/* $$$         else */
/* $$$            x(jr) = max( x(jr), bu(jr) ) */
/* $$$         end if */
	if (onbnd) {
	    x[jr] = bound;
	}
	js = 0;
	if (*elastc) {
	    js = hestat[jr];
	    hestat[jr] = 0;
	}
	if (js == 0) {
	    if (blbs[*kp] == bubs[*kp]) {
		jrstat = 4;
	    } else if (hitlow) {
		jrstat = 0;
	    } else {
		jrstat = 1;
	    }
	} else if (js == 1) {
	    if (bl[jr] == bu[jr]) {
		jrstat = 4;
	    } else if (onbnd) {
		jrstat = 0;
	    } else if (x[jr] < bu[jr]) {
		jrstat = -1;
	    } else {
		jrstat = 1;
	    }
	} else {
/*   js .eq. 2 */
	    if (bl[jr] == bu[jr]) {
		jrstat = 4;
	    } else if (onbnd) {
		jrstat = 1;
	    } else if (x[jr] > bl[jr]) {
		jrstat = -1;
	    } else {
		jrstat = 0;
	    }
	}
	if (*kp <= *m) {
/*           ============================================================ */
/*           A variable in B hit a bound. */
/*           Find column kSq = kBSq-m  of S to replace column kp of B. */
/*           If nS = 1 there is no choice. */
/*           ============================================================ */
	    if (*ns == 1) {
		kbsq = nbs;
		*pivot /= pbs1;
	    } else {
		dload_(m, &c_b5, &y2[1], &c__1);
		y2[*kp] = 1.;
		s2bsol_(iexit, &c__2, m, &y2[1], &y[1], &iw[1], leniw, &rw[1],
			 lenrw);
		if (*iexit != 0) {
		    return 0;
		}
		s5chzq_(m, mbs, n, nb, ns, &kbsq, pivot, &tolp0, ne, nloca, &
			loca[1], &inda[1], &acol[1], &kbs[1], &bl[1], &bu[1], 
			&xbs[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
		if (kbsq <= 0) {
		    s_wsfi(&io___226);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    e_wsfi();
		    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)80);
		    kbsq = nbs;
		}
	    }
	    ksq = kbsq - *m;
	    hs[jr] = jrstat;
	    *jbr = jr;
/* Outgoing basic */
	    *jsr = kbs[kbsq];
/* Outgoing superbasic */
	    kbs[kbsq] = *jbr;
	    *jbq = *jsr;
/* Incoming basic */
	    kbs[*kp] = *jsr;
	    blbs[*kp] = blbs[kbsq];
	    bubs[*kp] = bubs[kbsq];
	    xbs[*kp] = xbs[kbsq];
	    hs[*jbq] = 3;
	    if (*ns > 1 && *gotr) {
/*              Finish computing y(m+1), ..., y(m+nS). */
		y[kbsq] = -(*pivot + 1.);
		d__1 = 1. / *pivot;
		dscal_(ns, &d__1, &y[*m + 1], &c__1);
		s6rswp_(maxr, ns, lenr, &r__[1], &y2[1], &y[*m + 1], &ksq, &
			eps0);
	    }
/*           ------------------------------------------------------------ */
/*           Get a new  y1, used to modify L and U.  If the outgoing */
/*           superbasic just came in, we already have it. */
/*           ------------------------------------------------------------ */
	    if (*jsr != *jq) {
		s2unpk_(jbq, m, n, ne, &anorm, nloca, &loca[1], &inda[1], &
			acol[1], &y1[1]);
		s2bsol_(iexit, &c__0, m, &y1[1], &y[1], &iw[1], leniw, &rw[1],
			 lenrw);
		if (*iexit != 0) {
		    return 0;
		}
	    }
/*           Update the LU factors. */
	    ++iw[216];
	    s2bmod2_(&inform__, kp, m, &y1[1], &iw[1], leniw, &rw[1], lenrw);
	    if (inform__ == -1) {
		*lureq = 5;
	    }
/* Singular after LU mod */
	    if (inform__ == 2) {
		*lureq = 6;
	    }
/* Unstable LU mod */
	    if (inform__ == 7) {
		*lureq = 7;
	    }
/* Insufficient free memory */
	} else {
/*           ============================================================ */
/*           A variable in S hit a bound. */
/*           ============================================================ */
	    hs[jr] = jrstat;
	    *jsr = jr;
	    kbsq = *kp;
	    ksq = kbsq - *m;
	}
/*        Delete the kSq-th superbasic and adjust all arrays in BS order. */
	s5sdel_(&ksq, m, ns, &nbs, &kbs[1], &blbs[1], &bubs[1], &gbs[1], &rg[
		1], &xbs[1]);
	if (*gotr) {
/*           ------------------------------------------------------------ */
/*           Cyclically demote column kSq of R to position nS. */
/*           ------------------------------------------------------------ */
	    if (ksq < *ns) {
		s6rdel_(&ksq, maxr, ns, lenr, &r__[1], &eps);
	    }
	}
/* feasbl and gotH */
	--(*ns);
	--nbs;
/* -->     R can be checked here. */
    }
/* hitcon */
L900:
    return 0;
} /* s5qpit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QPit */
/* Subroutine */ int s5rchk_(integer *iexit, S_fp hprod, U_fp hprod1, integer 
	*itn, integer *minimz, integer *maxr, integer *lenr, integer *m, 
	integer *mbs, integer *n, integer *nb, integer *hvcalls, integer *nnh,
	 integer *ns, integer *ne, integer *nloca, integer *loca, integer *
	inda, doublereal *acol, integer *kbs, doublereal *r__, doublereal *v, 
	doublereal *w, doublereal *y, char *cu, integer *lencu, integer *iu, 
	integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer jq, js, nbs;
    static doublereal whw, eps0;
    static integer eigh;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal drsq;
    static integer ldiag;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm;
    extern /* Subroutine */ int s2bprd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), s2bsol_(integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s6rcol_(integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *), s5rsng_(integer *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s2unpk_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), s6rsol_(integer *, integer *, integer *, integer *,
	     doublereal *, doublereal *);
    static integer lencol;
    static logical posdef;
    static doublereal sgnobj;
    static logical singlr;
    static doublereal rnrmsq;
    static integer status;
    extern /* Subroutine */ int s2gathr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), s2scatr_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Rchk  computes the Cholesky factor R such that */
/*     R'R = Z'HZ.  The update corresponds to the addition of a new */
/*     column to Z. */

/*     On entry, */
/*        R     holds the columns of the factor associated with the */
/*              first jRadd-1 columns of Q. */

/*        nS    is the number of columns in R. */

/*     14 Mar 2001: First version (s6Rchk) based on SNOPT routine s5Rcol. */
/*     20 Jun 2008: Renamed s5Rchk. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --y;
    --kbs;
    --w;
    --v;
    --acol;
    --inda;
    --loca;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
    eigh = iw[200];
/* -1,0,1 for indef, psd and pdef QP Hessian */
    nbs = *m + *ns;
    sgnobj = (doublereal) (*minimz);
    *iexit = 0;
/*     ------------------------------------------------------------------ */
/*     Main loop to find a column of Z'HZ. */
/*     ------------------------------------------------------------------ */
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
/* Computing MIN */
	i__2 = js - 1;
	lencol = min(i__2,*nnh);
/*        --------------------------------------------------------------- */
/*        Get the nonlinear elements of the column of Z. */
/*        Find y such that B y = column jq. */
/*        Scatter the nonlinear part of y into w. */
/*        --------------------------------------------------------------- */
	jq = kbs[*m + js];
	s2unpk_(&jq, m, n, ne, &anorm, nloca, &loca[1], &inda[1], &acol[1], &
		w[1]);
	s2bsol_(iexit, &c__1, m, &w[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    return 0;
	}
	s2scatr_(nnh, m, &kbs[1], &c_b76, &y[1], &w[1]);
	if (jq <= *nnh) {
	    w[jq] = 1.;
	}
/*        --------------------------------------------------------------- */
/*        Compute  H*w  and  w'*H*w. */
/*        --------------------------------------------------------------- */
	whw = 0.;
	if (*nnh > 0) {
	    status = 0;
	    (*hprod)((U_fp)hprod1, nnh, &w[1], &v[1], &status, cu + 8, lencu, 
		    &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
	    ++(*hvcalls);
	    whw += ddot_(nnh, &w[1], &c__1, &v[1], &c__1);
	    if (*minimz < 0) {
		dscal_(nnh, &sgnobj, &v[1], &c__1);
		whw = sgnobj * whw;
	    }
	}
	rnrmsq = 0.;
	if (js > 1) {
/*           ------------------------------------------------------------ */
/*           Gather the nonlinear elements of v in w (= vBS). */
/*           Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB). */
/*           ------------------------------------------------------------ */
	    s2gathr_(nnh, &nbs, &kbs[1], &c_b146, &v[1], &w[1]);
	    s2bsol_(iexit, &c__2, m, &w[1], &v[1], &iw[1], leniw, &rw[1], 
		    lenrw);
	    if (*iexit != 0) {
		return 0;
	    }
	    if (*ns > 0) {
		s2bprd_(&c__1, &eps0, n, ns, &kbs[*m + 1], ne, nloca, &loca[1]
			, &inda[1], &acol[1], &c_b76, &v[1], m, &c_b146, &w[*
			m + 1], ns);
	    }
/*           ------------------------------------------------------------ */
/*           Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jS). */
/*           ------------------------------------------------------------ */
	    s6rsol_(&c__1, maxr, &lencol, lenr, &r__[1], &w[*m + 1]);
	    rnrmsq = ddot_(&lencol, &w[*m + 1], &c__1, &w[*m + 1], &c__1);
	}
	if (js <= *nnh) {
	    drsq = whw - rnrmsq;
/* New diagonal of R. */
	} else {
	    drsq = 0.;
	}
	w[*m + js] = drsq;
/*        Insert w(m+1:m+jS) as column jS of R. */
	s6rcol_(&js, maxr, &js, lenr, &r__[1], &w[*m + 1], &ldiag);
	s5rsng_(&eigh, &posdef, &singlr, itn, maxr, lenr, ns, &drsq, &r__[1], 
		&iw[1], leniw, &rw[1], lenrw);
	if (! (posdef || singlr)) {
	    *iexit = 6;
	    return 0;
	}
	r__[ldiag] = sqrt(drsq);
    }
    return 0;
} /* s5rchk_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Rchk */
/* Subroutine */ int s5rcol_(integer *iexit, S_fp hprod, U_fp hprod1, integer 
	*minimz, integer *jq, integer *jradd, doublereal *drsq, integer *
	ldiag, integer *maxr, integer *lenr, integer *m, integer *mbs, 
	integer *n, integer *nb, integer *hvcalls, integer *nnh, integer *ns, 
	integer *ne, integer *nloca, integer *loca, integer *inda, doublereal 
	*acol, integer *kbs, doublereal *r__, doublereal *v, doublereal *w, 
	doublereal *y, char *cu, integer *lencu, integer *iu, integer *leniu, 
	doublereal *ru, integer *lenru, char *cw, integer *lencw, integer *iw,
	 integer *leniw, doublereal *rw, integer *lenrw, ftnlen cu_len, 
	ftnlen cw_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer nbs;
    static doublereal whw, eps0;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s2bprd_(integer *, doublereal *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), s2bsol_(integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s6rcol_(integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *), s6rsol_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static integer lencol;
    static doublereal sgnobj;
    static integer sqstat;
    static doublereal rnrmsq;
    extern /* Subroutine */ int s2gathr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), s2scatr_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Rcol  computes column jRadd of the Cholesky factor R such that */
/*     R'R = Z'HZ.  The update corresponds to the addition of a new */
/*     column to Z. */

/*     On entry, */
/*        R     holds the columns of the factor associated with the */
/*              first jRadd-1 columns of Q. */

/*        y     is the vector such that B y = a(jq). */

/*        nS    is the number of columns in R. */

/*     11 Dec 1991: First version based on Qpsol routine Qpcolr. */
/*     24 Apr 1994: Columns of Nx no longer in Q. */
/*     27 Oct 2000: Previous version of s5Rcol. */
/*     04 Dec 2000: R converted to row-wise storage. */
/*     09 Dec 2004: Current version of s5Rcol. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --y;
    --kbs;
    --w;
    --v;
    --acol;
    --inda;
    --loca;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
    *iexit = 0;
    nbs = *m + *ns;
    sgnobj = (doublereal) (*minimz);
/* Computing MIN */
    i__1 = *jradd - 1;
    lencol = min(i__1,*nnh);
/*     ------------------------------------------------------------------ */
/*     Get w, the vector of nonlinear components of the new column of Z. */
/*     ------------------------------------------------------------------ */
/*     The input vector y satisfies B y = column jq. */
/*     Scatter the nonlinear components of y into w. */
    s2scatr_(nnh, m, &kbs[1], &c_b76, &y[1], &w[1]);
    if (*jq <= *nnh) {
	w[*jq] = 1.;
    }
/*     ------------------------------------------------------------------ */
/*     Compute  H*w  and  w'*H*w. */
/*     ------------------------------------------------------------------ */
    whw = 0.;
    if (*nnh > 0) {
	sqstat = 0;
	(*hprod)((U_fp)hprod1, nnh, &w[1], &v[1], &sqstat, cu + 8, lencu, &iu[
		1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1]
		, lenrw, (ftnlen)8, (ftnlen)8);
	++(*hvcalls);
	whw += ddot_(nnh, &w[1], &c__1, &v[1], &c__1);
	if (*minimz < 0) {
	    dscal_(nnh, &sgnobj, &v[1], &c__1);
	    whw = sgnobj * whw;
	}
    }
    rnrmsq = 0.;
    if (*jradd > 1) {
/* ------------------------------------------------------------- */
/* Gather the nonlinear elements of v in w (= vBS). */
/* Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB). */
/* ------------------------------------------------------------- */
	s2gathr_(nnh, &nbs, &kbs[1], &c_b146, &v[1], &w[1]);
	s2bsol_(iexit, &c__2, m, &w[1], &v[1], &iw[1], leniw, &rw[1], lenrw);
	if (*iexit != 0) {
	    return 0;
	}
	if (*ns > 0) {
	    s2bprd_(&c__1, &eps0, n, ns, &kbs[*m + 1], ne, nloca, &loca[1], &
		    inda[1], &acol[1], &c_b76, &v[1], m, &c_b146, &w[*m + 1], 
		    ns);
	}
/* ------------------------------------------------------------- */
/* Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jRadd). */
/* ------------------------------------------------------------- */
	s6rsol_(&c__1, maxr, &lencol, lenr, &r__[1], &w[*m + 1]);
	rnrmsq = ddot_(&lencol, &w[*m + 1], &c__1, &w[*m + 1], &c__1);
    }
    if (*jradd <= *nnh) {
	*drsq = whw - rnrmsq;
/* Square of the new diagonal of R. */
    } else {
	*drsq = 0.;
    }
    w[*m + *jradd] = *drsq;
/* Insert w(m+1:m+jRadd) as column jRadd of R. */
    s6rcol_(jradd, maxr, jradd, lenr, &r__[1], &w[*m + 1], ldiag);
    return 0;
} /* s5rcol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Rcol */
/* Subroutine */ int s5rg_(integer *m, integer *nbs, integer *n, integer *ns, 
	doublereal *tolz, integer *ne, integer *nloca, integer *loca, integer 
	*inda, doublereal *acol, doublereal *gbs, doublereal *pi, doublereal *
	rg, doublereal *rgnorm, integer *kbs)
{
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2bprd_(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    extern doublereal dnormi_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s5rg    calculates the reduced gradient  rg = gS - S'*pi. */

/*     23 Nov 1991: First version based on Minos routine m7rg. */
/*     16 Nov 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pi;
    --kbs;
    --gbs;
    --rg;
    --acol;
    --inda;
    --loca;

    /* Function Body */
    dcopy_(ns, &gbs[*m + 1], &c__1, &rg[1], &c__1);
    s2bprd_(&c__1, tolz, n, ns, &kbs[*m + 1], ne, nloca, &loca[1], &inda[1], &
	    acol[1], &c_b76, &pi[1], m, &c_b146, &rg[1], ns);
    *rgnorm = dnormi_(ns, &rg[1], &c__1);
    return 0;
} /* s5rg_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5rg */
/* Subroutine */ int s5rsng_(integer *eigh, logical *posdef, logical *singlr, 
	integer *itn, integer *maxr, integer *lenr, integer *ns, doublereal *
	drsq, doublereal *r__, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: Reduced Hessian is semi"
	    "definite.\002,\002 Square of diag, min diag = \002,1p,2e9.1)";
    static char fmt_9000[] = "(\002 Itn\002,i7,\002: Reduced Hessian is inde"
	    "finite.\002,\002 Square of diag, min diag = \002,1p,2e9.1)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[110];
    static doublereal condh, drmin, drmax;
    extern /* Subroutine */ int s6rcnd_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal hcndbd, drsqmn;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___257 = { 0, str, 0, fmt_1000, 110, 1 };
    static icilist io___258 = { 0, str, 0, fmt_9000, 110, 1 };


/*     ================================================================== */
/*     s5Rsng  estimates the inertia of the current reduced Hessian. */

/*     15 Jul 1995: First version of s5Rsng. */
/*     02 Aug 2003: snPRNT adopted. */
/*     17 Jun 2004: Current version of s5Rsng. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --iw;
    --rw;

    /* Function Body */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    if (*ns == 0) {
/*        --------------------------------------------------------------- */
/*        Vertices are positive definite by definition. */
/*        --------------------------------------------------------------- */
	*posdef = TRUE_;
	*singlr = FALSE_;
    } else {
/*        --------------------------------------------------------------- */
/*        Compute dRsqmn, the square of the smallest possible diagonal */
/*        of a positive-definite reduced Hessian. */
/*        --------------------------------------------------------------- */
	i__1 = *ns - 1;
	s6rcnd_(maxr, &i__1, lenr, &r__[1], &drmax, &drmin, &condh);
	drsqmn = drmax * (drmax / hcndbd);
	*posdef = *drsq >= drsqmn;
	*singlr = abs(*drsq) < drsqmn;
	if (*singlr || *posdef) {
	    if (*drsq < 0.) {
		if (*eigh == 1) {
		    s_wsfi(&io___257);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*drsq), (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&drsqmn, (ftnlen)sizeof(doublereal))
			    ;
		    e_wsfi();
		    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)110);
		}
		*drsq = max(0.,*drsq);
	    }
	} else if (*eigh == 0 || *eigh == 1) {
	    s_wsfi(&io___258);
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*drsq), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&drsqmn, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)110);
	}
    }
    return 0;
} /* s5rsng_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Rsng */
/* Subroutine */ int s5sdel_(integer *ksq, integer *m, integer *ns, integer *
	nbs, integer *kbs, doublereal *blbs, doublereal *bubs, doublereal *
	gbs, doublereal *rg, doublereal *xbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;

/*     ================================================================== */
/*     s5Sdel  deletes the kSqth superbasic variable from the arrays */
/*     kBS, blBS, blBS, gBS, rg and xBS. */

/*     16 Jun 2001: First version of s5Bswp. */
/*     16 Jun 2001: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     Shift all the arrays one place to the left. */
    /* Parameter adjustments */
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;

    /* Function Body */
    i__1 = *ns - 1;
    for (j = *ksq; j <= i__1; ++j) {
	k = *m + j;
	kbs[k] = kbs[k + 1];
	blbs[k] = blbs[k + 1];
	bubs[k] = bubs[k + 1];
	gbs[k] = gbs[k + 1];
	xbs[k] = xbs[k + 1];
	rg[j] = rg[j + 1];
    }
    return 0;
} /* s5sdel_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Sdel */
/* Subroutine */ int s5zp_(integer *iexit, integer *m, integer *mbs, integer *
	n, integer *nb, integer *ns, doublereal *eps0, doublereal *pnorm, 
	integer *ne, integer *nloca, integer *loca, integer *inda, doublereal 
	*acol, integer *kbs, doublereal *pbs, doublereal *y, integer *iw, 
	integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer nbs;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), s2bprd_(integer *, doublereal *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), s2bsol_(integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    extern doublereal dnormi_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s5Zp computes the free components of the search direction */
/*     p = Z pS, where pS is the search direction for the superbasics, */
/*     stored in  pBS(m+1:nBS) */

/*     On exit, the  free components of the search direction are stored */
/*     in pBS(1:nBS). The search direction for the basic variables is */
/*     stored in pBS(1),...,pBS(m). */

/*     20 Dec 2005: First version of s5Zp. */
/*     01 Dec 2012: Added pNorm <= 0.0 test. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --pbs;
    --kbs;
    --y;
    --acol;
    --inda;
    --loca;
    --iw;
    --rw;

    /* Function Body */
    *iexit = 0;
    nbs = *m + *ns;
    *pnorm = dnormi_(ns, &pbs[*m + 1], &c__1);
    if (*pnorm <= 0.) {
	*pnorm = 1.;
    }
/* First, compute  y = - S*pS and prepare to solve  B*pB = y */
/* for pB, the search direction for the basic variables. */
/* We first normalize y so the LU solver won't ignore */
/* too many "small" elements while computing pB. */
    d__1 = 1. / *pnorm;
    dscal_(ns, &d__1, &pbs[*m + 1], &c__1);
    s2bprd_(&c__0, eps0, n, ns, &kbs[*m + 1], ne, nloca, &loca[1], &inda[1], &
	    acol[1], &c_b76, &pbs[*m + 1], ns, &c_b5, &y[1], m);
/* Solve  B*pBS = y  and unnormalize all of pBS. */
    s2bsol_(iexit, &c__1, m, &y[1], &pbs[1], &iw[1], leniw, &rw[1], lenrw);
    if (*iexit != 0) {
	return 0;
    }
    dscal_(&nbs, pnorm, &pbs[1], &c__1);
    *pnorm = dnormi_(&nbs, &pbs[1], &c__1);
    return 0;
} /* s5zp_ */

