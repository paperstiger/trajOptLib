/* ./src/sn56qncg.f -- translated by f2c (version 20090411).
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
static integer c__31 = 31;
static integer c__2 = 2;
static doublereal c_b85 = -1.;
static doublereal c_b102 = 1.;
static doublereal c_b184 = .33333;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     file  sn56qncg.f */

/*     s5QN */
/*     s5QNdir  s5QNit   s5Sswp   s5Hzx   s5Hzx1   s5Zswp */
/*     SYMMLQ   s5Msolv */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5qn_(integer *iexit, integer *prob, char *probtag, 
	logical *elastc, integer *subopt, U_fp hprod, U_fp hprod1, S_fp qplog,
	 logical *gotr, logical *needlu, integer *typelu, logical *needx, 
	integer *lenr, integer *m, integer *maxs, integer *mbs, integer *n, 
	integer *nb, integer *ndegen, integer *hvcalls, integer *ngqp0, 
	integer *ngqp, integer *ngobj0, integer *ngobj, integer *nnh0, 
	integer *nnh, integer *ns, integer *itqp, integer *itqpmax, integer *
	itqptargt, integer *itn, integer *lemode, integer *lvlinf, integer *
	prtlvl, integer *minimz, integer *iobj, doublereal *sclobj, 
	doublereal *objadd, doublereal *objqp, doublereal *condhz, doublereal 
	*targtz, doublereal *tolfp, doublereal *tolqp, doublereal *tolx, 
	integer *ninf, doublereal *sinf, doublereal *wtinf, doublereal *
	pinorm, integer *ne, integer *nloca, integer *loca, integer *inda, 
	doublereal *acol, integer *hetype, integer *hestat, integer *hfeas, 
	integer *hs, integer *kbs, doublereal *ascale, doublereal *bl, 
	doublereal *bu, doublereal *blbs, doublereal *bubs, doublereal *gbs, 
	doublereal *gobj, doublereal *gqp, doublereal *hdx, doublereal *pbs, 
	doublereal *pi, doublereal *r__, doublereal *rc, doublereal *rg, 
	doublereal *rg2, integer *nrhs0, integer *nrhs, doublereal *rhs, 
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
    static char fmt_1700[] = "(\002 Itn\002,i7,\002: Subspace tolerance ="
	    " \002,1p,e11.3)";
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
    static integer jq, kp, jbq, jbr;
    static doublereal djq;
    static integer nbs, jsq, jsr;
    static char str[120];
    static doublereal djq0, eps0;
    extern /* Subroutine */ int s5rg_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), s5hs_(integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *);
    static integer kfac, nfac, kchk;
    static doublereal bold;
    static logical done;
    static doublereal fobj;
    static logical newb;
    static integer klog;
    static logical gote, goth;
    static integer kprc, ksav;
    static logical conv[3], prt10, luok;
    static integer maxr, nfix[2];
    static logical newx;
    static doublereal ftol[2], step, xtol[2], fobj0;
    extern /* Subroutine */ int s5inf_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static logical badcg;
    static doublereal tolx0;
    static logical needf;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical needr, chkpi;
    static integer ninfe;
    static logical needv;
    static doublereal etarg, sinfe, drmin;
    static logical newsb;
    static integer lumax, lureq;
    static logical newlu;
    static integer isumm, ksumm, nsmax, nswap;
    static doublereal anorm, drmax, ftoli, pivot, rgtol[2], tolrg, xtoli;
    extern /* Subroutine */ int s2bfac_(integer *, integer *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), icopy_(integer 
	    *, integer *, integer *, integer *, integer *), s5dgen_(integer *,
	     integer *, integer *, integer *, integer *, integer *, 
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
    static doublereal psnrm1;
    static logical chkfea;
    static doublereal hcndbd;
    extern /* Subroutine */ int s5setp_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s2unpk_(integer *, integer *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static integer lusiz0;
    extern /* Subroutine */ int s6rset_(logical *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    static doublereal xsnrm1;
    extern /* Subroutine */ int s5qnit_(integer *, U_fp, U_fp, logical *, 
	    logical *, logical *, logical *, logical *, logical *, logical *, 
	    logical *, logical *, logical *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen), s5setx_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    static logical feasbl;
    static integer kdegen;
    static doublereal deltaf, infbnd;
    static logical checkx;
    extern /* Subroutine */ int s2bmod2_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *);
    static logical unbndd, jstfea;
    static doublereal djqmod, featol, deltax;
    static logical parhes, cnvrgd, hitcon, incres, qpdone;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__, itnlim;
    static logical gotgqp;
    static integer itnfix;
    static logical optiml, prtlog, usegqp;
    static integer iprint, lvltol, mnewsb;
    static logical statpt;
    static integer muncon, nelast, kprprt, nfmove, nfreez, nonopt, nuncon, 
	    precon;
    static logical prtsum;
    static integer qpmode, sqstat;
    static doublereal bgrwth, djqprt, gsnorm, objprt, psnorm, rgnorm, rowerr, 
	    sgnobj, tolinc, weight, objslk, xsnorm;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2gathr_(integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *), s5zswp_(logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , s5einit_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), s5ewrap_(
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *), 
	    s2trylu_(integer *, integer *, integer *, integer *, logical *, 
	    integer *, integer *, integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static icilist io___82 = { 0, str, 0, fmt_1030, 120, 1 };
    static icilist io___85 = { 0, str, 0, fmt_1610, 120, 1 };
    static icilist io___86 = { 0, str, 0, fmt_1620, 120, 1 };
    static icilist io___95 = { 0, str, 0, fmt_1700, 120, 1 };
    static icilist io___98 = { 0, str, 0, fmt_8050, 120, 1 };
    static icilist io___99 = { 0, str, 0, fmt_8060, 120, 1 };
    static icilist io___100 = { 0, str, 0, fmt_1010, 120, 1 };
    static icilist io___101 = { 0, str, 0, fmt_1020, 120, 1 };
    static icilist io___121 = { 0, str, 0, fmt_1000, 120, 1 };


/*     ================================================================== */
/*     s5QN   solves a linear or quadratic program.  If the objective */
/*     is quadratic, a BFGS quasi-Newton method is used. */

/*     The problem type can be: */
/*       Prob = 0 FP   feasible point only */
/*       Prob = 1 LP   LP problem */
/*       Prob = 2 QP   QP problem */
/*       Prob = 3 FPE  feasible point for equalities only */
/*       Prob = 4 FPS  feasible point for QP subProblem */
/*       Prob = 5 QPS  QP subproblem */
/*       Prob = 6 QPP  FP subproblem with proximal point objective */

/*     ngQP = max( nnH, ngObj ) */

/*     The optimization can pass through the following phases: */

/*       Phase 1               find a feasible point for all variables */

/*       Elastic Phase 1       make the non-elastic variables feasible */
/*                             while allowing infeasible elastics */

/*       Phase 2               minimize the objective */

/*       Elastic Phase 2       minimize a composite objective while */
/*                             keeping the non-elastics feasible */

/*                             In this phase, lvlInf means the following: */

/*                 lvlInf = 0  zero     weight on the infeasibilities */
/*                                      (infeasibillities are ignored) */
/*                          1  finite   weight on the infeasibilities */
/*                          2  infinite weight on the infeasibilities */
/*                                      (the objective is ignored) */

/*     The array kBS is a permutation on the column indices. */
/*     kBS(1  :m )    holds the col. indices of the basic variables. */
/*     kBS(m+1:m+nS)  holds the col. indices of the superbasic variables. */
/*                    These nS columns indices must have hs(j) = 2. */

/*     On exit: */

/*      iExit       Result */
/*      -----       ------ */
/*       >0         Fatal error */
/*        0         QP solution found */
/*       -1         QP is infeasible */
/*       -2         QP is unbounded */
/*       -3         Too many iterations */
/*       -4         Void */
/*       -5         Too many superbasics */
/*       -6         Void */
/*       -7         Void */
/*       -8         Ill-conditioned CG null-space basis */
/*       -9         Too many subspace iterations */

/*     22 Jun 2001: First version of s5QN  based on Qpsol routine qpcore. */
/*     29 Jul 2003: Make sure R is never too ill-conditioned. */
/*     30 Jul 2003: Superbasic slacks allowed. */
/*     02 Aug 2003: snEXIT and snPRNT adopted. */
/*     03 Mar 2013: pinorm initialized to one instead of zero. */
/*     26 May 2013: infBnd used to identify infinite bounds. */
/*     ================================================================== */
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
/* > 0  => modify the LU factors */
/* (on/off) log     status */
/* (on/off) summary status */
/* # lines in log     file */
/* # lines in summary file */
/* >0 => Minor heading for iPrint */
/* >0 => Minor heading for iSumm */
/*     ------------------------------------------------------------------ */
/* symmlq itns for last minor */
    /* Parameter adjustments */
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
    infbnd = rw[70];
/* definition of an infinite bound */
    etarg = rw[83];
/* rgNorm tolerance */
    hcndbd = rw[85];
/* bound on the condition of Hz */
    iprint = iw[12];
/* Print file */
    isumm = iw[13];
/* Summary file */
    maxr = iw[52];
/* max columns of R */
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
    nfac = iw[210];
/* # of LU factorizations */
    qpmode = iw[208];
/* Current QP solver */
    precon = iw[209];
/* Current precon mode (based on QPslvr) */
    *iexit = 0;
    muncon = 20;
    sqstat = 0;
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
    kprc = 0;
/* last sec scanned in part. prc */
    lureq = 0;
/*     ------------------------------------------------------------------ */
/*     s5QN operates in either ``Normal'' or ``Elastic'' mode. */
/*     Everything is normal unless a weighted sum is being minimized or */
/*     the constraints are infeasible. */
/*     The logical feasbl refers to the non-elastic variables. */
/*     Note that feasbl can be false while in elastic mode. */
/*     wtInf  is the optional parameter Infeasibility Weight. */
/*     ------------------------------------------------------------------ */
    feasbl = FALSE_;
    gote = FALSE_;
    goth = *nnh > 0;
    gotgqp = *ngqp > 0;
/*     jstPh1 = stop at the end of phase 1 (in normal or elastic mode) */
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
/*     needR  = true if R needs to be set or reset. */
    needr = ! jstph1 && (qpmode == 2 || qpmode == 1 && precon == 1);
    parhes = *ns > maxr && *gotr;
/*     call dcopy ( 2, rw( ftol1), 1,  ftol, 1 ) */
/*     call dcopy ( 2, rw( xtol1), 1,  xtol, 1 ) */
/*     call dcopy ( 2, rw(rgtol1), 1, rgtol, 1 ) */
    xtol[0] = .1;
    xtol[1] = 1e-6;
    ftol[0] = xtol[0] * .1;
/* Computing 2nd power */
    d__1 = xtol[1];
    ftol[1] = d__1 * d__1;
    rgtol[0] = .001;
    rgtol[1] = 1e-7;
    *objqp = 0.;
    pivot = 0.;
    step = 0.;
    *sinf = 0.;
    tolrg = 0.;
    ninfe = 0;
    jq = 0;
    djq = 0.;
    djq0 = 0.;
    djqprt = 0.;
    objslk = 0.;
    jbq = 0;
/* x(jBq) is the incoming   BS */
    jbr = 0;
/* x(jBr) is the outgoing   BS */
    jsq = 0;
/* x(jSq) is the incoming SBS */
    jsr = 0;
/* x(jSr) is the outgoing SBS */
    lvltol = 1;
    kprprt = 0;
    sgnobj = (doublereal) (*minimz);
    iw[387] = 0;
    rw[184] = 100.;
    rw[185] = 100.;
    chkfea = TRUE_;
    chkpi = TRUE_;
    cnvrgd = FALSE_;
    newlu = TRUE_;
    newx = FALSE_;
    unbndd = FALSE_;
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
/* !    ======================Start of main loop========================== */
/* +    do while (.not. QPdone  .and.  iExit .eq. 0) */
L100:
    if (! qpdone && *iexit == 0) {
/*        =============================================================== */
/*        Check the initial  x  and move it onto  ( A  -I )*x = b. */
/*        If needLU is true, this will require a basis factorization. */
/*        =============================================================== */
/*        If necessary,  factorize the basis  ( B = LU ) and set x. */
/*        If needLU is false on entry to s5QN, the first call to s2Bfac */
/*        will try to use existing factors. */
/*        If needLU is true on entry to s5QN, an LU factorization of */
/*        type typeLU is computed. */

/*        The reason for the requested LU is as follows. */

/*        LUreq =  0  First LU for a given subproblem */
/*        LUreq =  1  Frequency */
/*        LUreq =  2  LU nonzeros increased */
/*        LUreq =  3 */
/*        LUreq =  4 */
/*        LUreq =  5  Singular after LU mod */
/*        LUreq =  6  Unstable LU mod (growth in new column of U) */
/*        LUreq =  7  Not enough memory */
/*        LUreq =  8 */
/*        LUreq =  9 */
/*        LUreq = 10  Row error in setx */
/*        LUreq = 11  Big  dx   in setx */

/*        LUreq = 20 */
/*        LUreq = 21  Iterative refinement failed in QP */
/*        LUreq = 22  Unbounded QP */
/*        LUreq = 23  Infeasibility after refactorization */
/*        LUreq = 24  Small directional derivative in QP */
/*        LUreq = 25  Ill-conditioned Z in QP */
/*        LUreq = 26  Indefinite Z'HZ in QP */
/*        LUreq = 27  R singular after bound swap in QP */
/*        --------------------------------------------------------------- */
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
		if (nswap > 0) {
		    *gotr = FALSE_;
		}
/* Reset R. */
		if (prt10) {
		    iw[223] = 1;
		}
	    }
	    if (*iexit != 0) {
		goto L100;
	    }
	    parhes = *ns > maxr && *gotr;
	    cnvrgd = FALSE_;
	    gote = FALSE_;
/* Check hEstat in elastic mode. */
	    *needx = FALSE_;
	    newx = TRUE_;
	    chkfea = TRUE_;
	    chkpi = TRUE_;
	    pivot = 0.;
	    nuncon = 0;
	}
	nbs = *m + *ns;
	newsb = FALSE_;
	*ninf = 0;
	*sinf = 0.;
	optiml = FALSE_;
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
	    dload_(&nbs, &c_b5, &gbs[1], &c__1);
	    s5inf_(&nbs, &featol, &infbnd, ninf, sinf, &hfeas[1], &blbs[1], &
		    bubs[1], &gbs[1], &xbs[1]);
	    if (*ninf > 0) {
/*              Non-elastics are infeasible. */
/*              If necessary, switch back to the feasibility phase, after */
/*              refactorization. */
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
		if (prt10 && iw[215] == 0) {
		    s_wsfi(&io___82);
		    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*ninf), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*sinf), (ftnlen)sizeof(doublereal)
			    );
		    e_wsfi();
		    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)120);
		}
	    }
/*           Feasbl = true  means that the non-elastics are feasible. */
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
/*           ------------------------------------------------------------ */
/*           The non-elastic variables just became feasible. Exit. */
/*           ------------------------------------------------------------ */
	    *condhz = 0.;
	    djqprt = 0.;
	    rgnorm = 0.;
	    dload_(m, &c_b5, &pi[1], &c__1);
	    *pinorm = 1.;
/* pinorm = max(norm(pi), 1.0) */
	    optiml = TRUE_;
	} else {
	    if (feasbl) {
/*              --------------------------------------------------------- */
/*              The nonelastics are feasible. */
/*              (Elastc = false means no elastics.) */
/*              --------------------------------------------------------- */
/*              The objective row (if defined) is computed in all phases. */
		objslk = 0.;
		if (*iobj != 0) {
		    objslk = xbs[iw[205]] * *sclobj;
		}
		if (needf && (jstfea || newx)) {
/*                 ====================================================== */
/*                 If just feasible or x has changed, initialize */
/*                 the QP objective, the QP gradient and R. */
/*                 ObjQP = quadratic term to be minimized (or maximized). */
/*                 The objective row and ObjAdd do not appear in ObjQP. */
/*                 ObjQP is updated after each QP step. */
/*                 ====================================================== */
		    if (gotgqp) {
			if (*hvcalls == 0) {
			    sqstat = 1;
			}
			s5qpfg_((U_fp)hprod, (U_fp)hprod1, ngqp, ngobj0, 
				ngobj, nnh, &sqstat, hvcalls, objqp, &gobj[1],
				 &gqp[1], lenx0, nx0, &x0[1], &x[1], &y[1], 
				cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
				cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
				ftnlen)8, (ftnlen)8);
			sqstat = 0;
		    }
		    if (goth && needr && ! (*gotr)) {
			s6rset_(gotr, &maxr, ns, lenr, &r__[1], &y[1], condhz)
				;
/*                    gotR = .true. ! done by s6Rset */
			parhes = *ns > maxr;
		    }
/*                 ------------------------------------------------------ */
/*                 Gather the QP gradient in BS order. */
/*                 Assign the nonzero components of gBS. */
/*                 ------------------------------------------------------ */
		    if (gotgqp) {
			s2gathr_(ngqp, &nbs, &kbs[1], &sgnobj, &gqp[1], &gbs[
				1]);
		    }
		    if (*iobj > 0) {
			gbs[iw[205]] = sgnobj * *sclobj;
		    }
		    if (*elastc && ninfe > 0) {
			s5egrd_(nb, &nbs, wtinf, &hestat[1], &kbs[1], &gbs[1])
				;
		    }
		}
/*              --------------------------------------------------------- */
/*              See if it's time to suboptimize. */
/*              NOTE: We must not suboptimize if all steps have been */
/*              degenerate. */
/*              --------------------------------------------------------- */
/* jstFea .or. newx */
		if (*subopt != 0 || nfmove == 0) {
/*                 Relax */
		} else {
		    if (*ns >= nsmax) {
			*subopt = 1;
			if (prt10) {
			    s_wsfi(&io___85);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&mnewsb, (ftnlen)sizeof(
				    integer));
			    e_wsfi();
			    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)120);
			}
		    } else if (*itqp >= *itqptargt) {
			*subopt = 2;
			if (prt10) {
			    s_wsfi(&io___86);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&(*itqptargt), (ftnlen)
				    sizeof(integer));
			    e_wsfi();
			    snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)120);
			}
		    }
		}
	    }
/*           ============================================================ */
/*           Check for an approximate stationary point. */
/*           If the gradient has changed, compute the reduced  gradient */
/*           ============================================================ */
/* feasible */
	    if (! feasbl || jstfea || newx) {
		dcopy_(m, &gbs[1], &c__1, &y[1], &c__1);
		s5setp_(&inform__, m, &chkpi, pinorm, &y[1], &pi[1], &iw[1], 
			leniw, &rw[1], lenrw);
		if (inform__ != 0) {
		    if (inform__ > 0) {
			*iexit = inform__;
		    } else {
/* pi is infinite or contains a NaN/Inf. */
			s2trylu_(itn, &c__11, ns, &lureq, &luok, typelu, &iw[
				1], leniw, &rw[1], lenrw);
			if (! luok) {
			    *iexit = 43;
			}
		    }
		    goto L100;
		}
		rgnorm = 0.;
		if (*ns > 0) {
		    s5rg_(m, &nbs, n, ns, &eps0, ne, nloca, &loca[1], &inda[1]
			    , &acol[1], &gbs[1], &pi[1], &rg[1], &rgnorm, &
			    kbs[1]);
		}
	    }
/* .not. feasbl  .or.  jstFea  .or.  newx */
	    if (tolrg == 0.) {
		tolrg = etarg * rgnorm;
	    }
	    if (feasbl) {
		rw[186] = *tolqp;
	    } else {
		rw[186] = *tolfp;
	    }
	    if (parhes) {
		gsnorm = dnormi_(&maxr, &rg[1], &c__1);
	    } else {
		gsnorm = rgnorm;
	    }
	    statpt = gsnorm <= tolrg * .1 || gsnorm <= rgtol[1] * *pinorm || 
		    cnvrgd;
	    if (parhes && statpt) {
/*              Swap the largest reduced gradient in Z2 into the front of */
/*              Z2 and see if it is significantly large. */
/*              Reset statpt if necessary. */
		s5zswp_(&statpt, m, &maxr, maxs, &nbs, lenr, ns, &tolrg, &kbs[
			1], &blbs[1], &bubs[1], &gbs[1], &r__[1], &rg[1], &
			xbs[1], &rw[1], lenrw);
	    }
	    if (*gotr) {
		s6rcnd_(&maxr, ns, lenr, &r__[1], &drmax, &drmin, condhz);
		if (*condhz > hcndbd) {
		    *gotr = FALSE_;
		    *condhz = 0.;
		    if (needr) {
			s6rset_(gotr, &maxr, ns, lenr, &r__[1], &y[1], condhz)
				;
/*                    gotR = .true. ! done by s6Rset */
		    }
		}
	    }
	    kprprt = kprc;
	    jq = 0;
	    if (statpt) {
/*              --------------------------------------------------------- */
/*              Compute Lagrange multipliers. */
/*              --------------------------------------------------------- */
		djq0 = djq;
/* save djq in case of bad statpt */
		djq = 0.;
		nuncon = 0;
		usegqp = feasbl && gotgqp;
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
		if (lvltol == 1 && (*ns >= *maxs || optiml)) {
		    lvltol = 2;
		    tolrg = rw[186] * *pinorm;
		    if (prt10) {
			s_wsfi(&io___95);
			do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer)
				);
			do_fio(&c__1, (char *)&tolrg, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)120);
		    }
		    if (rgnorm > tolrg) {
			optiml = FALSE_;
			newsb = FALSE_;
		    }
		}
	    }
/* statpt */
	}
/* jstPh1 */
	qpdone = optiml || unbndd;
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
			    &xbs[1], nrhs0, nrhs, &rhs[1], &x[1], &y[1], &y2[
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
	    } else {
		*needx = TRUE_;
		unbndd = FALSE_;
		cnvrgd = FALSE_;
		goto L100;
	    }
	    if (jstph1) {
/*              Relax, we are about to exit */
	    } else {
/*              ========================================================= */
/*              Print the details of the final iteration. */
/*              ========================================================= */
		objprt = 0.;
		if (feasbl) {
		    if (needf) {
			objprt = *objadd + objslk + *objqp;
		    }
		    if (needv) {
			objprt += sgnobj * *wtinf * *sinf;
		    }
		}
		(*qplog)(prob, probtag, elastc, gotr, &jstfea, &feasbl, m, 
			mbs, nnh, ns, &jsq, &jbr, &jsr, &iw[220], &iw[221], 
			itn, itqp, &kprprt, lvlinf, &pivot, &step, ninf, sinf,
			 wtinf, &objprt, condhz, &djqprt, &rgnorm, &kbs[1], &
			xbs[1], &iw[1], leniw, (ftnlen)20);
	    }
	    jbq = 0;
	    jbr = 0;
	    jsq = 0;
	    jsr = 0;
	    kprprt = 0;
	    iw[387] = 0;
	    djqprt = 0.;
	    *condhz = 0.;
/*           ------------------------------------------------------------ */
/*           Convergence. */
/*           ------------------------------------------------------------ */
	    if (*ninf > 0) {
/*              No feasible point. */
/*              Stop or continue in elastic mode, depending on the */
/*              specified level of infeasibility. */
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
/* infeasible (should not happen) */
			}
		    } else {
/*                    The constraints are infeasible in Normal mode. */
/*                    Print a message and start Elastic Phase 1. */
/* if .not. Elastc */
			if (prt10) {
			    s_wsfi(&io___98);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, probtag, (ftnlen)20);
			    e_wsfi();
			    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
			    s_wsfi(&io___99);
			    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(
				    integer));
			    e_wsfi();
			    snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
			    iw[223] = 1;
			    iw[225] = 1;
			}
			*elastc = TRUE_;
			needf = *lvlinf != 2;
/* F1 required in phase 2 */
			needv = *lvlinf != 0;
/* F2 required in phase 2 */
			qpdone = FALSE_;
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
			s_wsfi(&io___100);
			do_fio(&c__1, (char *)&djq, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&jq, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&rgnorm, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)120);
		    }
		} else {
		    if (prt10 && klog == 1) {
			s_wsfi(&io___101);
			do_fio(&c__1, (char *)&rgnorm, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&(*pinorm), (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__31, str, &iw[1], leniw, (ftnlen)120);
		    }
		}
	    }
	} else {
/*           ============================================================ */
/*           Take a series of reduced gradient steps until the new point */
/*           approximately minimizes the objective on the working set, */
/*           or lies on the boundary of a new constraint. */
/* +          ============================================================ */
/*           Repeat                   (until  hitcon or cnvrgd or statpt) */
/* not done */
L500:
	    cnvrgd = FALSE_;
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
		    objprt, condhz, &djqprt, &rgnorm, &kbs[1], &xbs[1], &iw[1]
		    , leniw, (ftnlen)20);
	    jbq = 0;
	    jbr = 0;
	    jsq = 0;
	    jsr = 0;
	    kprprt = 0;
	    iw[387] = 0;
	    djqprt = 0.;
	    *condhz = 0.;
	    if (newsb) {
/*                 ------------------------------------------------------ */
/*                 A nonbasic has been selected to become superbasic. */
/*                 Compute the vector y such that B y = column jq. */
/*                 ------------------------------------------------------ */
/*                 Set the level to which rgNorm must be reduced in the */
/*                 current subspace before we consider moving off another */
/*                 constraint. */
		djqmod = abs(djq);
		rgnorm = max(rgnorm,djqmod);
		tolrg = etarg * djqmod;
		if (*ns + 1 > *maxs) {
		    *iexit = -5;
		    goto L100;
		}
		djqprt = djq;
/*                 ------------------------------------------------------ */
/*                 Compute the vector pBS such that B pB = column jq. */
/*                 pBS is a multiple of part of the new column of  Z  and */
/*                 is used to define the QN search direction. */
/*                 ------------------------------------------------------ */
/*                 Unpack column jq into  y1  and solve  B*pB = y1. */
/*                 The altered  y1  satisfies  L*y1 = ajq. */
/*                 It is used later in s5QNit to modify L and U. */
		s2unpk_(&jq, m, n, ne, &anorm, nloca, &loca[1], &inda[1], &
			acol[1], &y1[1]);
		s2bsol_(iexit, &c__1, m, &y1[1], &pbs[1], &iw[1], leniw, &rw[
			1], lenrw);
		if (*iexit != 0) {
		    return 0;
		}
		rw[192] = dnormi_(m, &pbs[1], &c__1) / anorm;
/* $$$                  if (rw(condZ) .ge. targtZ  .and.  nS .gt. 0) then */
/* $$$                     iExit = -8 */
/* $$$                     go to 100 */
/* $$$                  end if */
	    }
/* newSB */
	    if (*itn >= itnlim || *itqp >= *itqpmax) {
		*iexit = -3;
		goto L100;
	    }
	    ++(*itqp);
	    ++(*itn);
/*              Decide if we want to print something this iteration. */
	    prtlog = iprint > 0 && *prtlvl >= 1 && *itqp % klog == 0;
	    prtsum = isumm > 0 && *prtlvl >= 1 && *itqp % ksumm == 0;
	    iw[218] = 0;
	    iw[219] = 0;
	    if (prtlog) {
		iw[218] = 1;
	    }
	    if (prtsum) {
		iw[219] = 1;
	    }
	    if (feasbl) {
		fobj0 = sgnobj * (objslk + *objqp) + *wtinf * *sinf;
	    }
/*              --------------------------------------------------------- */
/*              Take a ``reduced gradient'' step. */
/*              The new  x  will either minimize the objective on the */
/*              working set or lie on the boundary of a new constraint. */
/*              --------------------------------------------------------- */
	    s5qnit_(&inform__, (U_fp)hprod, (U_fp)hprod1, elastc, &feasbl, &
		    gotgqp, &goth, gotr, &hitcon, &incres, &needf, &needv, &
		    newsb, itn, lenr, m, mbs, &maxr, maxs, n, nb, hvcalls, 
		    nnh0, nnh, ns, ngqp0, ngqp, ndegen, &lureq, &kp, &jbq, &
		    jsq, &jbr, &jsr, &jq, &lvltol, &nfmove, &nuncon, &djq0, &
		    djq, minimz, iobj, sclobj, objqp, condhz, &featol, &pivot,
		     pinorm, &rgnorm, &step, &tolinc, ninf, sinf, wtinf, &
		    psnorm, &psnrm1, &xsnorm, &xsnrm1, ne, nloca, &loca[1], &
		    inda[1], &acol[1], &hetype[1], &hestat[1], &hfeas[1], &hs[
		    1], &kbs[1], &bl[1], &bu[1], &blbs[1], &bubs[1], &gbs[1], 
		    &gqp[1], &hdx[1], &pbs[1], &pi[1], &rg[1], &rg2[1], &r__[
		    1], &x[1], &xbs[1], &y[1], &y1[1], &y2[1], cu + 8, lencu, 
		    &iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], 
		    leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
/*              Check for trouble. */
/*              inform values are  -1, 0 >0 */
	    if (inform__ != 0) {
		if (inform__ > 0) {
		    *iexit = inform__;
/* Fatal LU error */
		} else if (inform__ == -1) {
		    *iexit = -2;
/* unbounded */
		}
		goto L100;
	    }
	    nbs = *m + *ns;
/*              --------------------------------------------------------- */
/*              Stop if there are errors. */
/*              --------------------------------------------------------- */
	    done = lureq != 0;
	    badcg = nuncon >= muncon && lvltol == 2;
	    if (! done) {
/*                 ------------------------------------------------------ */
/*                 Reset R if it is too ill-conditioned. */
/*                 ------------------------------------------------------ */
		if (*gotr) {
		    s6rcnd_(&maxr, ns, lenr, &r__[1], &drmax, &drmin, condhz);
		    if (*condhz > hcndbd) {
			*gotr = FALSE_;
			*condhz = 0.;
			if (needr) {
			    s6rset_(gotr, &maxr, ns, lenr, &r__[1], &y[1], 
				    condhz);
/*                          gotR = .true. ! done by s6Rset */
			}
		    }
		}
		if (! hitcon) {
/*                    --------------------------------------------------- */
/*                    The QP step was unconstrained. */
/*                    Test for convergence on the current working set. */
/*                    --------------------------------------------------- */
		    parhes = *ns > maxr && *gotr;
		    gsnorm = rgnorm;
		    if (parhes) {
			psnorm = psnrm1;
			xsnorm = xsnrm1;
			gsnorm = dnormi_(&maxr, &rg[1], &c__1);
		    }
		    if (feasbl && needf) {
			if (*iobj > 0) {
			    objslk = xbs[iw[205]] * *sclobj;
			}
			fobj = sgnobj * (objslk + *objqp) + *wtinf * *sinf;
			ftoli = ftol[lvltol - 1];
			xtoli = xtol[lvltol - 1];
			deltax = step * psnorm;
			deltaf = fobj0 - fobj;
			conv[0] = deltax <= xtoli * (xsnorm + 1.);
			conv[1] = deltaf <= ftoli * (abs(fobj0) + 1.);
			conv[2] = gsnorm <= tolrg;
			cnvrgd = conv[0] && conv[1] && conv[2];
		    }
		    statpt = gsnorm <= tolrg * .1 || gsnorm <= rgtol[1] * *
			    pinorm;
		}
/* .not. Hitcon */
	    }
/* +          until    (Hitcon .or. statpt .or. cnvrgd.or. done.or.badCG) */
	    if (! (hitcon || statpt || cnvrgd || done || badcg)) {
		goto L500;
	    }
/*           ------------------------------------------------------------ */
	    if (badcg) {
		*iexit = -9;
/* too many subspace iterations */
		goto L100;
	    }
	    if (lureq > 0) {
		s2trylu_(itn, &lureq, ns, &lureq, &luok, typelu, &iw[1], 
			leniw, &rw[1], lenrw);
		if (! luok) {
		    *iexit = 43;
		    goto L100;
		}
	    }
	    ++iw[215];
/* itns since the last LU */
	    newlu = FALSE_;
	    newx = FALSE_;
	    chkpi = FALSE_;
/*           ============================================================ */
/*           Test for error condition and/or frequency interrupts. */
/*           ============================================================ */
	    if (*itn % ksav == 0) {
		icopy_(nb, &hs[1], &c__1, &iy1[1], &c__1);
		s5hs_(&c__1, nb, &bl[1], &bu[1], &iy1[1], &x[1]);
		s4ksav_(minimz, m, n, nb, ns, mbs, itn, ninf, sinf, objqp, &
			kbs[1], &iy1[1], &ascale[1], &bl[1], &bu[1], &x[1], &
			xbs[1], cw + 8, lencw, &iw[1], leniw, (ftnlen)8);
	    }
/*           Increment featol every iteration. */
	    featol += tolinc;
/*           Every kdegen iterations, reset featol and move nonbasic */
/*           variables onto their bounds if they are very close. */
	    if (*itn % kdegen == 0) {
		s5dgen_(&inform__, &c__2, prtlvl, nb, ninf, itn, &featol, 
			tolx, &tolinc, &hs[1], &bl[1], &bu[1], &x[1], &itnfix,
			 nfix, &tolx0, &iw[1], leniw, &rw[1], lenrw);
		*needx = inform__ > 0;
	    }
/*           Refactorize the basis if it has been modified too much. */
	    if (lureq == 0) {
		if (iw[216] >= kfac - 1) {
		    lureq = 1;
		} else if (iw[216] >= 20 && iw[173] + iw[174] > lumax) {
		    bgrwth = (doublereal) (iw[173] + iw[174]);
		    bold = (doublereal) lusiz0;
		    bgrwth /= bold;
		    if (prt10) {
			s_wsfi(&io___121);
			do_fio(&c__1, (char *)&bgrwth, (ftnlen)sizeof(
				doublereal));
			e_wsfi();
			snprnt_(&c__21, str, &iw[1], leniw, (ftnlen)120);
		    }
		    lureq = 2;
		}
		if (lureq > 0) {
		    *typelu = 3;
		}
	    }
/*           Update the LU factors. */
	    if (lureq == 0 && iw[217] > 0) {
		++iw[216];
		s2bmod2_(&inform__, &kp, m, &y1[1], &iw[1], leniw, &rw[1], 
			lenrw);
		if (inform__ == -1) {
		    lureq = 5;
		}
/* Singular after LU mod */
		if (inform__ == 2) {
		    lureq = 6;
		}
/* Unstable LU mod */
		if (inform__ == 7) {
		    lureq = 7;
		}
/* Insufficient free memory */
		if (lureq > 0) {
		    s2trylu_(itn, &lureq, ns, &lureq, &luok, typelu, &iw[1], 
			    leniw, &rw[1], lenrw);
		    if (! luok) {
			*iexit = 43;
			goto L100;
		    }
		}
	    }
	    if (lureq == 0) {
		checkx = iw[215] % kchk == 0;
		if (checkx && ! (*needx)) {
		    s5setx_(&inform__, &c__1, itn, m, n, nb, &nbs, &rowerr, 
			    ne, nloca, &loca[1], &inda[1], &acol[1], &kbs[1], 
			    &xbs[1], nrhs0, nrhs, &rhs[1], &x[1], &y[1], &y2[
			    1], &iw[1], leniw, &rw[1], lenrw);
		    lureq = inform__;
		    cnvrgd = FALSE_;
		}
		if (lureq > 0) {
		    *typelu = 3;
		}
	    }
	}
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
} /* s5qn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QN */
/* Subroutine */ int s5qndir_(U_fp hprod, U_fp hprod1, logical *feasbl, 
	logical *exactp, integer *itn, logical *qpterm, doublereal *condhz, 
	doublereal *rgnorm, integer *maxr, integer *ns, integer *lenr, 
	doublereal *r__, doublereal *g, doublereal *p, doublereal *w, 
	doublereal *gp, integer *ne, integer *nloca, integer *loca, integer *
	inda, doublereal *acol, char *cu, integer *lencu, integer *iu, 
	integer *leniu, doublereal *ru, integer *lenru, char *cw, integer *
	lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: CG gives  gp = \002,1p,"
	    "e8.1,\002 rtol = \002,e8.1)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer lr1, lr2, ls1, ls2, ls3;
    static char str[120];
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal rtol;
    static integer nout;
    extern /* Subroutine */ int s5hzx_();
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical goodb;
    static doublereal tolcg, shift;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal rnorm;
    static integer istop;
    static doublereal ynorm;
    static logical checka;
    extern /* Subroutine */ int s6rsol_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static integer lprdbg, precon, qpmode, cgitmx;
    static doublereal hznorm;
    static logical sympre;
    extern /* Subroutine */ int symmlq_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, U_fp, U_fp, logical *, logical *, logical *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     U_fp, U_fp, integer *, integer *, integer *, integer *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen), snprnt_(integer *, char *, integer 
	    *, integer *, ftnlen);
    extern /* Subroutine */ int s5msolv_();

    /* Fortran I/O blocks */
    static icilist io___144 = { 0, str, 0, fmt_1000, 120, 1 };


/*     ================================================================== */
/*     s5QNdir  computes a search direction  p  for the superbasic */
/*     variables, using the current reduced gradient  g. */

/*     16 Nov 2001: First version of s5QNdir. */
/*     26 Jul 2003: cgItn added. */
/*     14 Mar 2004: Argument QPterm added for lvlInf = 2. */
/*     04 Dec 2004: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* Total symmlq iterations */
/* symmlq itns for last minor */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --w;
    --p;
    --g;
    --r__;
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
    lprdbg = iw[85];
/* >0    => private debug print */
    cgitmx = iw[97];
/* CG iteration limit */
    nout = iw[151];
/* unit # for printed messages */
    qpmode = iw[208];
/* Current QP solver   (based on QPslvr) */
    precon = iw[209];
/* Current precon mode (based on QPslvr) */
    lr1 = iw[353];
/* r1(maxS) SYMMLQ real work vector */
    lr2 = iw[354];
/* r1(maxS) SYMMLQ real work vector */
    ls1 = iw[355];
/* s1(maxS) SYMMLQ real work vector */
    ls2 = iw[356];
/* s2(maxS) SYMMLQ real work vector */
    ls3 = iw[357];
/* s3(maxS) SYMMLQ real work vector */
    tolcg = rw[54];
/* cg tolerance */
    if (lprdbg == 0) {
	nout = 0;
    }
    checka = FALSE_;
    *exactp = TRUE_;
    goodb = FALSE_;
    sympre = FALSE_;
/* No preconditioning done by SYMMLQ */
    iw[387] = 0;
    shift = -1e-6;
/* symmlq solves with  (A - shift*I) ! */
    rtol = tolcg * min(1.,*rgnorm);
/*     Steepest-descent direction. */
    dcopy_(ns, &g[1], &c__1, &p[1], &c__1);
    if (*feasbl && *qpterm) {
/*        Preconditioned CG.  Save  w  such that  R'*w = rg. */
	if (qpmode == 1 && precon == 1 || qpmode == 2) {
	    s6rsol_(&c__1, maxr, ns, lenr, &r__[1], &p[1]);
	}
	dcopy_(ns, &p[1], &c__1, &w[1], &c__1);
	if (qpmode == 2 && *ns > *maxr || qpmode == 1) {
	    *exactp = FALSE_;
	    symmlq_(ns, &w[1], &rw[lr1], &rw[lr2], &rw[ls1], &rw[ls2], &p[1], 
		    &rw[ls3], (U_fp)s5hzx_, (U_fp)s5msolv_, &checka, &goodb, &
		    sympre, &shift, &nout, &cgitmx, &rtol, &istop, &iw[387], &
		    hznorm, condhz, &rnorm, &ynorm, (U_fp)hprod, (U_fp)hprod1,
		     ne, nloca, &loca[1], &inda[1], &acol[1], cu + 8, lencu, &
		    iu[1], leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw,
		     &rw[1], lenrw, (ftnlen)8, (ftnlen)8);
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
	    if (istop == -1) {
		dcopy_(ns, &w[1], &c__1, &p[1], &c__1);
	    }
	}
	dcopy_(ns, &p[1], &c__1, &w[1], &c__1);
	if (qpmode == 1 && precon == 1 || qpmode == 2) {
	    s6rsol_(&c__0, maxr, ns, lenr, &r__[1], &p[1]);
	}
    }
/*     Change the sign of p. */
/* feasbl */
    dscal_(ns, &c_b85, &p[1], &c__1);
    *gp = ddot_(ns, &g[1], &c__1, &p[1], &c__1);
    if (*gp >= 0.) {
	s_wsfi(&io___144);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*gp), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rtol, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
    }
    iw[386] += iw[387];
    return 0;
} /* s5qndir_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QNdir */
/* Subroutine */ int s5qnit_(integer *iexit, S_fp hprod, U_fp hprod1, logical 
	*elastc, logical *feasbl, logical *gotgqp, logical *goth, logical *
	gotr, logical *hitcon, logical *incres, logical *needf, logical *
	needv, logical *newsb, integer *itn, integer *lenr, integer *m, 
	integer *mbs, integer *maxr, integer *maxs, integer *n, integer *nb, 
	integer *hvcalls, integer *nnh0, integer *nnh, integer *ns, integer *
	ngqp0, integer *ngqp, integer *ndegen, integer *lureq, integer *kp, 
	integer *jbq, integer *jsq, integer *jbr, integer *jsr, integer *jq, 
	integer *lvltol, integer *nfmove, integer *nuncon, doublereal *djq0, 
	doublereal *djq, integer *minimz, integer *iobj, doublereal *sclobj, 
	doublereal *objqp, doublereal *condhz, doublereal *featol, doublereal 
	*pivot, doublereal *pinorm, doublereal *rgnorm, doublereal *step, 
	doublereal *tolinc, integer *ninf, doublereal *sinf, doublereal *
	wtinf, doublereal *psnorm, doublereal *psnrm1, doublereal *xsnorm, 
	doublereal *xsnrm1, integer *ne, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, integer *hetype, integer *hestat, 
	integer *hfeas, integer *hs, integer *kbs, doublereal *bl, doublereal 
	*bu, doublereal *blbs, doublereal *bubs, doublereal *gbs, doublereal *
	gqp, doublereal *hdx, doublereal *pbs, doublereal *pi, doublereal *rg,
	 doublereal *rg2, doublereal *r__, doublereal *x, doublereal *xbs, 
	doublereal *y, doublereal *y1, doublereal *y2, char *cu, integer *
	lencu, integer *iu, integer *leniu, doublereal *ru, integer *lenru, 
	char *cw, integer *lencw, integer *iw, integer *leniw, doublereal *rw,
	 integer *lenrw, ftnlen cu_len, ftnlen cw_len)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: CG gives  gp = \002,1p,"
	    "e8.1)";
    static char fmt_9999[] = "(\002 Itn\002,i7,\002: Chzq failed in s5QNit!"
	    "!\002)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static doublereal t;
    static integer m1, jr, js;
    static doublereal gp;
    static integer nz1, nz2, nbs;
    static doublereal eps, php;
    static integer ksq;
    static char str[120];
    static doublereal pbs1, eps0, eps2;
    extern /* Subroutine */ int s5rg_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), s5zp_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *);
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jesq, kbsq;
    static logical move;
    static doublereal gpqp, tolp;
    static integer ntry, modr1;
    static doublereal tolp0;
    extern /* Subroutine */ int s5bsx_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), s6rqn_(integer *, logical 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dload_(integer *, 
	    doublereal *, doublereal *, integer *), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal bigdx;
    static logical chkpi;
    static doublereal rgdel;
    static logical onbnd;
    static doublereal exact, bound, anorm;
    static logical uncon;
    static doublereal stepb, phpqp, pnorm, stepp;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dcopy_(integer *, doublereal 
	    *, integer *, doublereal *, integer *), s6radd_(integer *, 
	    integer *, integer *, doublereal *), s5einf_(integer *, integer *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), s5egrd_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *), s6rdel_(
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *), s5sdel_(integer *, integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *), s2bsol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), s5chzq_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static doublereal psnrm2;
    extern /* Subroutine */ int s2unpk_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), s5step_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, logical *, logical *, logical *, logical *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), s5setp_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), s6rset_(logical *, integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *);
    static doublereal xsnrm2, objchg, infbnd;
    extern /* Subroutine */ int s6rswp_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , s5sswp_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static logical unbndd, parhes, exactp;
    static doublereal sgnobj;
    extern doublereal dnormi_(integer *, doublereal *, integer *);
    static integer inform__, infpiv;
    static doublereal sclpiv;
    static logical hitlow;
    static integer jqstat, jrstat;
    static doublereal tolpiv;
    static integer sqstat;
    static doublereal stepmx;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen), s2gathr_(integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *), s2scatr_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *), s5qndir_(
	    S_fp, U_fp, logical *, logical *, integer *, logical *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, char *, integer *, integer *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___173 = { 0, str, 0, fmt_1000, 120, 1 };
    static icilist io___194 = { 0, str, 0, fmt_9999, 120, 1 };


/*     ================================================================== */
/*     s5QNit performs a QP step. */

/*     On entry, */
/*        newSB = true implies that variable jq just went superbasic. */
/*                In this case: */
/*                pBS(1:m)  satisfies B pBS = a(jq). */
/*                 y1(1:m)  satisfies L  y1 = a(jq). */

/*      iExit       Result */
/*      -----       ------ */
/*       -1         unbounded */
/*        0         normal exit */
/*       >0         Fatal LU error */
/*     On exit, */
/*        pBS contains the most recent QP search direction. */

/*     incres says if the new variable should increase or not. */
/*            It is used only in linear mode when s5Pric is moving free */
/*            nonbasics toward their closest bound (and djq = zero). */

/*     gotR   says if a useful quasi-Newton preconditioner  R  exists. */
/*            It will be false if no preconditioning is being used. */
/*            If  preconditioned CG (QN) is being used, then gotR may */
/*            be true even if the current itn is infeasible. */
/*            This allows  R  to be updated in during a temporary loss */
/*            of feasibility. */

/*     parHes is true if  R  is not big enough to store a full */
/*            preconditioner for the reduced Hessian. */
/*            The null space is then  Z = ( Z1  Z2 ),  and  R  estimates */
/*            the Hessian only in the subspace  Z1.  A diagonal estimate */
/*            is used (also in  R) for the subspace  Z2. */

/*     12 Jun 2001: First   version of s5QNit based on s5Qpit */
/*                  and MINOS routine m7rgit. */
/*     14 Jul 2001: Different gp needed for parHes. */
/*     02 Aug 2003: snPRNT adopted. */
/*     19 Jun 2008: hEstat accessed in elastic mode only. */
/*     19 Jun 2006: s5Zp added to compute Z*p. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* xBS(kObj) is the obj. slack */
/* > 0  => modify the LU factors */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --pi;
    --xbs;
    --pbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --hfeas;
    --rg2;
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
/* unit round-off.  IEEE DP  2.22e-16 */
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    eps2 = rw[4];
/* eps**(1/2)       IEEE DP  1.49e-08 */
    tolpiv = rw[60];
/* excludes small elements of pBS. */
    infbnd = rw[70];
/* definition of an infinite bound */
    bigdx = rw[72];
/* unbounded step. */
    *iexit = 0;
    inform__ = 0;
    sqstat = 0;
    chkpi = FALSE_;
    unbndd = FALSE_;
    sgnobj = (doublereal) (*minimz);
    m1 = *m + 1;
    nbs = *m + *ns;
    if (*newsb) {
/*        --------------------------------------------------------------- */
/*        New superbasic. */
/*        --------------------------------------------------------------- */
	++(*ns);
	++nbs;
	kbs[nbs] = *jq;
	xbs[nbs] = x[*jq];
	blbs[nbs] = bl[*jq];
	bubs[nbs] = bu[*jq];
	jqstat = hs[*jq];
	if (*gotr) {
/* Add a unit column to R at position nS. */
	    s6radd_(maxr, lenr, ns, &r__[1]);
	}
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
		++(*ninf);
		if (*incres) {
		    if (jqstat == 1 || jqstat == 4) {
			hestat[*jq] = 2;
			bubs[nbs] = infbnd;
			if (*feasbl) {
			    gbs[nbs] += *wtinf;
			    blbs[nbs] = bu[*jq];
			}
		    }
		} else {
		    if (jqstat == 0 || jqstat == 4) {
			hestat[*jq] = 1;
			blbs[nbs] = -infbnd;
			if (*feasbl) {
			    gbs[nbs] -= *wtinf;
			    bubs[nbs] = bl[*jq];
			}
		    }
		}
	    }
	}
/*        --------------------------------------------------------------- */
/*        In phase 1, or in phase 2 for an LP, price can select nonbasics */
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
/*     Compute the search direction for the superbasics. */
/*     ------------------------------------------------------------------ */
/* newSB */
    parhes = *ns > *maxr && *gotr;
    if (parhes) {
	nz1 = *maxr;
    } else {
	nz1 = *ns;
    }
    nz2 = *ns - *maxr;
/*     ------------------------------------------------------------------ */
/*     Store the free components of the search direction in pBS(1:nBS). */
/*     First, find the search direction pS for the superbasics. */
/*     The vector v such that R*v = rg is held in y2 for the BFGS update. */
/*     ------------------------------------------------------------------ */
L100:
    s5qndir_((S_fp)hprod, (U_fp)hprod1, feasbl, &exactp, itn, needf, condhz, 
	    rgnorm, maxr, ns, lenr, &r__[1], &rg[1], &pbs[m1], &y2[1], &gp, 
	    ne, nloca, &loca[1], &inda[1], &acol[1], cu + 8, lencu, &iu[1], 
	    leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw,
	     (ftnlen)8, (ftnlen)8);
    if (gp >= 0.) {
	*gotr = FALSE_;
	*lureq = 24;
	goto L900;
    }
    pbs1 = pbs[*m + *ns];
    if (*newsb && *feasbl) {
/*        Check that pBS(nBS) is feasible wrt blBS(nBS) and buBS(nBS). */
/*        A large rgtol may give a pBS(nBS) with the wrong sign. */
	if (*djq * pbs1 > 0.) {
	    --(*ns);
	    --nbs;
	    hs[*jq] = jqstat;
	    if (*elastc) {
		hestat[*jq] = jesq;
	    }
	    *jq = 0;
	    *djq = *djq0;
	    *jsq = -(*jsq);
	    *newsb = FALSE_;
	    if (*lvltol == 2 && *gotr) {
		*gotr = FALSE_;
		s6rset_(gotr, maxr, ns, lenr, &r__[1], &y[1], condhz);
/*              gotR   = .true.  ! done by s6Rset */
	    } else {
		*lvltol = 2;
	    }
	    goto L100;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Find the norm of  pS.  Put the search direction for the basic */
/*     variables in pBS(1)  ,...,pBS(m). */
/*     ------------------------------------------------------------------ */
    *xsnrm1 = dnormi_(&nz1, &xbs[m1], &c__1);
    *psnrm1 = dnormi_(&nz1, &pbs[m1], &c__1);
    xsnrm2 = 0.;
    psnrm2 = 0.;
    if (parhes) {
	xsnrm2 = dnormi_(&nz2, &xbs[m1 + *maxr], &c__1);
	psnrm2 = dnormi_(&nz2, &pbs[m1 + *maxr], &c__1);
    }
    *psnorm = max(*psnrm1,psnrm2);
    *xsnorm = max(*xsnrm1,xsnrm2);
    s5zp_(iexit, m, mbs, n, nb, ns, &eps0, &pnorm, ne, nloca, &loca[1], &inda[
	    1], &acol[1], &kbs[1], &pbs[1], &y[1], &iw[1], leniw, &rw[1], 
	    lenrw);
    if (*iexit != 0) {
	goto L900;
    }
    if (*feasbl) {
/*        --------------------------------------------------------------- */
/*        Compute y = pBS(scattered) and Hdx(scattered). */
/*        The vector Hdx is used to update the objective and gradient of */
/*        the QP.  Form  gpQP  and  pHpQP  for the quadratic term. */
/*        gp = gpQP - pBS(kObj) + terms from the elastic gradient. */
/*        --------------------------------------------------------------- */
	php = 0.;
	if (*needf && (*gotgqp || *goth)) {
	    s2scatr_(ngqp, &nbs, &kbs[1], &c_b102, &pbs[1], &y[1]);
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
		php = sgnobj * phpqp;
	    }
	}
    }
/* feasbl */
    if (gp >= 0.) {
	s_wsfi(&io___173);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&gp, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
	*lureq = 24;
	goto L900;
    }
    *newsb = FALSE_;
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
/*        Find if the step is constrained or unconstrained. */
/*        If R has been flagged as singular, we double check by trying */
/*        to compute the QP minimizer along pBS.  If the minimizer */
/*        exists,  the singularity tolerance must be too large. */
    if (*feasbl) {
	uncon = stepp * php > -gp;
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
    *hitcon = ! uncon;
    if (*hitcon) {
	*nuncon = 0;
	*step = stepb;
	*pivot = -pbs[*kp];
    } else {
	++(*nuncon);
	*step = -gp / php;
	*pivot = 0.;
    }
/*     ------------------------------------------------------------------ */
/*     Compute ObjChg, the predicted change in fObj. */
/*     Note: pHpQP and pHp are defined before and after scaling by */
/*     sgnObj. */
/*     Update the value and gradient of the quadratic obj term. */
/*     ------------------------------------------------------------------ */
    if (*feasbl) {
	if (*step > 0.) {
/* Computing 2nd power */
	    d__1 = *step;
	    objchg = *step * gp + php * .5 * (d__1 * d__1);
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
    }
    if (*feasbl && move) {
	++(*nfmove);
    }
/*     ------------------------------------------------------------------ */
/*     Update the basic variables xBS. */
/*     ------------------------------------------------------------------ */
    daxpy_(&nbs, step, &pbs[1], &c__1, &xbs[1], &c__1);
    s5bsx_(&c__1, &nbs, nb, &kbs[1], &x[1], &xbs[1]);
    if (*hitcon) {
    }
    modr1 = 0;
    if (*feasbl && *step > 0.) {
/*        =============================================================== */
/*        gQP has been changed. */
/*        Compute the new gBS and reduced gradient rg2. */
/*        Update the reduced Hessian. */
/*        =============================================================== */
	dload_(&nbs, &c_b5, &gbs[1], &c__1);
	if (*needf) {
	    if (*gotgqp) {
		s2gathr_(ngqp, &nbs, &kbs[1], &sgnobj, &gqp[1], &gbs[1]);
	    }
	    if (*iobj > 0) {
		gbs[iw[205]] = sgnobj * *sclobj;
	    }
	}
	if (*elastc && *needv) {
	    s5einf_(nb, &nbs, &hestat[1], &kbs[1], ninf, sinf, &bl[1], &bu[1],
		     &x[1]);
	    s5egrd_(nb, &nbs, wtinf, &hestat[1], &kbs[1], &gbs[1]);
	}
	dcopy_(m, &gbs[1], &c__1, &y[1], &c__1);
	s5setp_(iexit, m, &chkpi, pinorm, &y[1], &pi[1], &iw[1], leniw, &rw[1]
		, lenrw);
	if (*iexit != 0) {
	    goto L900;
	}
	s5rg_(m, &nbs, n, ns, &eps0, ne, nloca, &loca[1], &inda[1], &acol[1], 
		&gbs[1], &pi[1], &rg2[1], rgnorm, &kbs[1]);
	if (*gotr) {
	    s6rqn_(&modr1, &exactp, maxr, ns, lenr, &r__[1], &gp, &rg[1], &
		    rg2[1], &pbs[m1], &y2[1], step, &eps2, &eps0);
	}
    } else {
	dcopy_(ns, &rg[1], &c__1, &rg2[1], &c__1);
    }
    if (uncon) {
/*        =============================================================== */
/*        The step is unconstrained. */
/*        =============================================================== */
	dcopy_(ns, &rg2[1], &c__1, &rg[1], &c__1);
	iw[217] = 0;
/* No LU update */
    } else {
/*        =============================================================== */
/*        There is a blocking variable. */
/*        It could be a fixed variable, whose new state must be 4. */
/*        =============================================================== */
/* hit constraint */
	jr = kbs[*kp];
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
/*           If nS = 1, it must be the entering SB column. */
/*           ============================================================ */
/*           if (nS .eq. 1) then */
/*              kBSq  = nBS */
/*              pivot = pivot/pBS1 */
/*           else */
	    dload_(m, &c_b5, &y2[1], &c__1);
	    y2[*kp] = 1.;
/* Solve  B'yB = ep */
/* Set      y2 = ep */
	    s2bsol_(&inform__, &c__2, m, &y2[1], &y[1], &iw[1], leniw, &rw[1],
		     lenrw);
	    s5chzq_(m, mbs, n, nb, ns, &kbsq, pivot, &tolp0, ne, nloca, &loca[
		    1], &inda[1], &acol[1], &kbs[1], &bl[1], &bu[1], &xbs[1], 
		    &y[1], &iw[1], leniw, &rw[1], lenrw);
	    if (kbsq <= 0) {
		s_wsfi(&io___194);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)120);
		kbsq = nbs;
	    }
/*           end if */
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
	    gbs[*kp] = gbs[kbsq];
	    hs[*jbq] = 3;
/*           Finish computing yS = (y(m+1), ..., y(m+nS)). */
	    y[kbsq] = -(*pivot + 1.);
	    d__1 = 1. / *pivot;
	    dscal_(ns, &d__1, &y[m1], &c__1);
	    if (*gotr && ksq <= *maxr) {
		s6rswp_(maxr, &nz1, lenr, &r__[1], &y2[1], &y[m1], &ksq, &
			eps0);
	    }
/*           Modify  pi  using  y  where  B' y = e(kp). */
	    t = rg2[ksq] / *pivot;
	    daxpy_(m, &t, &y[1], &c__1, &pi[1], &c__1);
/* Computing MAX */
	    d__1 = dnormi_(m, &pi[1], &c__1);
	    *pinorm = max(d__1,1.);
/*           ------------------------------------------------------------ */
/*           Get a new  y1, used to modify L and U.  If the outgoing */
/*           superbasic just came in, we already have it. */
/*           ------------------------------------------------------------ */
/*           if (jSr .ne. jq) then */
	    s2unpk_(jbq, m, n, ne, &anorm, nloca, &loca[1], &inda[1], &acol[1]
		    , &y1[1]);
	    s2bsol_(&inform__, &c__0, m, &y1[1], &y[1], &iw[1], leniw, &rw[1],
		     lenrw);
	    iw[217] = 1;
/*           end if */
/* Update the LU */
	} else {
/*           ============================================================ */
/*           A variable in S hit a bound. */
/*           ============================================================ */
	    hs[jr] = jrstat;
	    *jsr = jr;
	    kbsq = *kp;
	    ksq = kbsq - *m;
	    iw[217] = 0;
/* No LU update */
	}
	if (*feasbl) {
	    s5rg_(m, &nbs, n, ns, &eps0, ne, nloca, &loca[1], &inda[1], &acol[
		    1], &gbs[1], &pi[1], &rg[1], rgnorm, &kbs[1]);
	}
/*        =============================================================== */
/*        If necessary, swap the largest reduced-gradient in  Z2  into */
/*        the front of  Z2,  so it will end up at the end of  Z1. */
/*        =============================================================== */
	rgdel = (d__1 = rg[ksq], abs(d__1));
	if (parhes && ksq <= *maxr) {
	    s5sswp_(m, maxr, lenr, ns, &nbs, &kbs[1], &blbs[1], &bubs[1], &
		    gbs[1], &r__[1], &rg[1], &xbs[1]);
	}
/*        --------------------------------------------------------------- */
/*        Delete the kSq-th superbasic, updating R if it exists and */
/*        adjusting all arrays in BS order. */
/*        --------------------------------------------------------------- */
	if (*gotr && ksq < *ns) {
	    s6rdel_(&ksq, maxr, ns, lenr, &r__[1], &eps);
	}
	s5sdel_(&ksq, m, ns, &nbs, &kbs[1], &blbs[1], &bubs[1], &gbs[1], &rg[
		1], &xbs[1]);
	--(*ns);
	nbs = *m + *ns;
	if (*rgnorm <= rgdel) {
	    *rgnorm = 0.;
	    if (*ns > 0) {
		*rgnorm = dnormi_(ns, &rg[1], &c__1);
	    }
	}
    }
/* if uncon */
L900:
    return 0;
} /* s5qnit_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5QNit */
/* Subroutine */ int s5sswp_(integer *m, integer *maxr, integer *lenr, 
	integer *ns, integer *nbs, integer *kbs, doublereal *blbs, doublereal 
	*bubs, doublereal *gbs, doublereal *r__, doublereal *rg, doublereal *
	xbs)
{
    static integer j, k, j1, k1, k2;
    static doublereal bl1, bu1, rg1;
    static integer nz2;
    static doublereal gbs1, xbs1;
    static integer lastr, ldiag1, ldiag2;
    static doublereal rdiag1;
    extern integer idamax_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s5Sswp  (superbasic swap)  finds the largest reduced gradient in */
/*     the range  rg(maxR+1), ..., rg(nS)  and swaps it into position */
/*     maxR + 1  (so we know where it is). */

/*     16 Jun 2001: First version based on MINOS routine m6swap. */
/*     16 Jun 2001: Current version of s5Sswp. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --r__;
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;

    /* Function Body */
    k1 = *maxr + 1;
    if (*ns > k1) {
	nz2 = *ns - *maxr;
	k2 = *maxr + idamax_(&nz2, &rg[k1], &c__1);
	if (k2 > k1) {
	    j = *m + k1;
	    k = *m + k2;
	    lastr = *maxr * k1 / 2;
	    ldiag1 = lastr + 1;
	    ldiag2 = lastr + (k2 - *maxr);
	    rdiag1 = r__[ldiag1];
	    rg1 = rg[k1];
	    j1 = kbs[j];
	    bl1 = blbs[j];
	    bu1 = bubs[j];
	    gbs1 = gbs[j];
	    xbs1 = xbs[j];
	    r__[ldiag1] = r__[ldiag2];
	    rg[k1] = rg[k2];
	    kbs[j] = kbs[k];
	    blbs[j] = blbs[k];
	    bubs[j] = bubs[k];
	    gbs[j] = gbs[k];
	    xbs[j] = xbs[k];
	    r__[ldiag2] = rdiag1;
	    rg[k2] = rg1;
	    kbs[k] = j1;
	    blbs[k] = bl1;
	    bubs[k] = bu1;
	    gbs[k] = gbs1;
	    xbs[k] = xbs1;
	}
    }
    return 0;
} /* s5sswp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Sswp */
/* Subroutine */ int s5hzx_(integer *ns, doublereal *x, doublereal *zhzx, 
	U_fp hprod, U_fp hprod1, integer *ne, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    static integer m, n, nb, lr, ly, ly3, mbs, nnh, lkbs, lenr, maxr, maxs;
    extern /* Subroutine */ int s5hzx1_(U_fp, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);

/*     ================================================================== */
/*     s5Hzx  is called by SYMMLQ to compute the vector */
/*       ZHZx = Z'HZ x  or  R^(-T) Z'HZ R^(-1) x. */
/*     (Preconditioning with R'R is done explicitly, */
/*     not via SYMMLQ's Msolve.) */

/*     16 Nov 2001: First version of s5Hzx. */
/*     17 Aug 2004: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* number of Hx products */
    /* Parameter adjustments */
    --zhzx;
    --x;
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
    lenr = iw[28];
/* R(lenR) is the reduced Hessian factor */
    n = iw[15];
/* copy of the number of columns */
    m = iw[16];
/* copy of the number of rows */
    nnh = iw[24];
/*    max( nnObj, nnJac ) */
    maxr = iw[52];
/* max columns of R */
    maxs = iw[53];
/* max # of superbasics */
    lkbs = iw[292];
/* kBS(mBS)    = ( B  S ) list */
    lr = iw[295];
/* R(lenR)     = factor of Z'HZ */
    ly = iw[311];
/* y (nb)      =  real work vector */
    ly3 = iw[314];
/* y3(nb)      =  real work vector */
    nb = n + m;
    mbs = m + maxs;
    s5hzx1_((U_fp)hprod, (U_fp)hprod1, &iw[188], &nnh, &m, &mbs, &n, &nb, ns, 
	    &maxr, ne, nloca, &loca[1], &inda[1], &acol[1], &iw[lkbs], &x[1], 
	    &zhzx[1], &rw[ly], &rw[ly3], &lenr, &rw[lr], cu + 8, lencu, &iu[1]
	    , leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
	    lenrw, (ftnlen)8, (ftnlen)8);
    return 0;
} /* s5hzx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Hzx */
/* Subroutine */ int s5hzx1_(S_fp hprod, U_fp hprod1, integer *hvcalls, 
	integer *nnh, integer *m, integer *mbs, integer *n, integer *nb, 
	integer *ns, integer *maxr, integer *ne, integer *nloca, integer *
	loca, integer *inda, doublereal *acol, integer *kbs, doublereal *x, 
	doublereal *zhzx, doublereal *y, doublereal *y1, integer *lenr, 
	doublereal *r__, char *cu, integer *lencu, integer *iu, integer *
	leniu, doublereal *ru, integer *lenru, char *cw, integer *lencw, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw, ftnlen 
	cu_len, ftnlen cw_len)
{
    static integer m1, nbs;
    static doublereal eps0;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), s2bprd_(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), s2bsol_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *), s6rsol_(integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *);
    static doublereal sgnobj;
    static integer precon, inform__, minimz, sqstat;
    extern /* Subroutine */ int s2gathr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), s2scatr_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Hzx1  does the work for s5Hzx. */
/*     Computes the vector ZHZx = Z'HZ * x  or  R^(-T) Z'HZ R^(-1) x. */

/*     16 Nov 2001: First version of s5Hzx1. */
/*     17 Aug 2004: Current version. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --kbs;
    --y1;
    --y;
    --zhzx;
    --x;
    --acol;
    --inda;
    --loca;
    --r__;
    cu -= 8;
    --iu;
    --ru;
    cw -= 8;
    --iw;
    --rw;

    /* Function Body */
    eps0 = rw[2];
/* eps**(4/5)       IEEE DP  3.00e-13 */
    minimz = iw[199];
/* (-1)(+1)    => (max)(min) */
    precon = iw[209];
/* Current precon mode (based on QPslvr) */
    nbs = *m + *ns;
    m1 = *m + 1;
    sgnobj = (doublereal) minimz;
    sqstat = 0;
/*     ------------------------------------------------------------------ */
/*     Set y = Z*x. */
/*     ------------------------------------------------------------------ */
    dcopy_(ns, &x[1], &c__1, &y[m1], &c__1);
    if (precon == 1) {
	s6rsol_(&c__0, maxr, ns, lenr, &r__[1], &y[m1]);
    }
/*     Compute  y1 = - S*yS. */
    s2bprd_(&c__0, &eps0, n, ns, &kbs[m1], ne, nloca, &loca[1], &inda[1], &
	    acol[1], &c_b85, &y[m1], ns, &c_b5, &y1[1], m);
/*     Solve  B*y = y1. */
    s2bsol_(&inform__, &c__1, m, &y1[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
/*     ------------------------------------------------------------------ */
/*     Set  y = H*y = H*Z*x. */
/*     ------------------------------------------------------------------ */
    s2scatr_(nnh, &nbs, &kbs[1], &c_b102, &y[1], &y1[1]);
/*     call s8Hx */
/*    &   ( Hvcalls, nnH, y1, y, cw, lencw, iw, leniw, rw, lenrw ) */
    (*hprod)((U_fp)hprod1, nnh, &y1[1], &y[1], &sqstat, cu + 8, lencu, &iu[1],
	     leniu, &ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], 
	    lenrw, (ftnlen)8, (ftnlen)8);
/* y1 = dx,  y = Hdx */
    ++(*hvcalls);
/*     ------------------------------------------------------------------ */
/*     Compute  ZHZx = Z'y. */
/*     ------------------------------------------------------------------ */
/*     Gather y1 = yBS and solve  B'y = y1B. */
    s2gathr_(nnh, &nbs, &kbs[1], &sgnobj, &y[1], &y1[1]);
    s2bsol_(&inform__, &c__2, m, &y1[1], &y[1], &iw[1], leniw, &rw[1], lenrw);
/*     Set ZHZx = y1S - S'y. */
    dcopy_(ns, &y1[m1], &c__1, &zhzx[1], &c__1);
    s2bprd_(&c__1, &eps0, n, ns, &kbs[m1], ne, nloca, &loca[1], &inda[1], &
	    acol[1], &c_b85, &y[1], m, &c_b102, &zhzx[1], ns);
    if (precon == 1) {
	s6rsol_(&c__1, maxr, ns, lenr, &r__[1], &zhzx[1]);
    }
    return 0;
} /* s5hzx1_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Hz1 */
/* Subroutine */ int s5zswp_(logical *statpt, integer *m, integer *maxr, 
	integer *maxs, integer *nbs, integer *lenr, integer *ns, doublereal *
	tolrg, integer *kbs, doublereal *blbs, doublereal *bubs, doublereal *
	gbs, doublereal *r__, doublereal *rg, doublereal *xbs, doublereal *rw,
	 integer *lenrw)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k, k1, k2, lr;
    static doublereal rg1;
    static integer jz1;
    static doublereal eps, gbs1, xbs1, blbs1, bubs1;
    static integer lastr, ldiag1;
    static doublereal rdiag1, rgmin1;
    extern /* Subroutine */ int s6rdel_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), s5sdel_(integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal rgnrm2;
    extern /* Subroutine */ int s5sswp_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/*     ================================================================== */
/*     s5Zswp  (subspace convergence)  decides whether or not */
/*     optimization should continue on the current subspace. */
/*     On exit,  statpt = .false.  means it should, */
/*               statpt = .true    means it should not. */

/*     R  is a partial Hessian for the */
/*     first  maxR  superbasic variables.  The superbasics are then */
/*     in two sets  Z1  (containing   nZ1 = maxR       variables) */
/*             and  Z2  (containing   nZ2 = nS - maxR  variables). */

/*     The null-space matrix is similarly partitioned as  Z = ( Z1  Z2 ). */

/*     The normal convergence test is first applied to  Z1.  If it looks */
/*     like we are near enough to an optimum in that restricted */
/*     subspace, we find the largest reduced gradient in  Z2  (index k2). */
/*     If this is less than  tolrg  we exit with  nxtphs = 3  (to ask for */
/*     price).  Otherwise we make room in  Z1  for the corresponding */
/*     variable by the moving the superbasic with the smallest reduced */
/*     gradient in  Z1  (index  k1)  to the end of  Z2. */

/*     16 Jun 2001: First version based on MINOS routine m7sscv. */
/*     19 Jul 2001: Current version of s5Sswp. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --rg;
    --xbs;
    --gbs;
    --bubs;
    --blbs;
    --kbs;
    --r__;
    --rw;

    /* Function Body */
    eps = rw[1];
/*     Swap the largest reduced gradient in  Z2  to the front of  Z2 */
/*     and see if it is significantly large. */
/* unit round-off.  IEEE DP  2.22e-16 */
    s5sswp_(m, maxr, lenr, ns, nbs, &kbs[1], &blbs[1], &bubs[1], &gbs[1], &
	    r__[1], &rg[1], &xbs[1]);
    k2 = *maxr + 1;
    rgnrm2 = (d__1 = rg[k2], abs(d__1));
    if (rgnrm2 > *tolrg) {
/*        Find the smallest component of  Z1'g. */
	rgmin1 = abs(rg[1]);
	k1 = 1;
	i__1 = *maxr;
	for (k = 1; k <= i__1; ++k) {
	    if (rgmin1 >= (d__1 = rg[k], abs(d__1))) {
		rgmin1 = (d__1 = rg[k], abs(d__1));
		k1 = k;
	    }
	}
	if (rgmin1 < rgnrm2) {
/*           Save the relevant values. */
	    *statpt = FALSE_;
	    lastr = *maxr * (*maxr + 1) / 2;
	    ldiag1 = (k1 - 1) * *maxr + (3 - k1) * k1 / 2;
/* Magic formula! */
	    rdiag1 = r__[ldiag1];
	    rg1 = rg[k1];
	    k = *m + k1;
	    jz1 = kbs[k];
	    blbs1 = blbs[k];
	    bubs1 = bubs[k];
	    gbs1 = gbs[k];
	    xbs1 = xbs[k];
/*           Delete the k1-th variable from  Z1,  and shift the remaining */
/*           superbasics in  Z1  and  Z2  one place to the left. */
	    s6rdel_(&k1, maxr, ns, lenr, &r__[1], &eps);
	    s5sdel_(&k1, m, ns, nbs, &kbs[1], &blbs[1], &bubs[1], &gbs[1], &
		    rg[1], &xbs[1]);
/*           Put the old k1-th superbasic in at the very end. */
	    lr = lastr + (*ns - *maxr);
	    r__[lr] = rdiag1;
	    rg[*ns] = rg1;
	    kbs[*nbs] = jz1;
	    blbs[*nbs] = blbs1;
	    bubs[*nbs] = bubs1;
	    gbs[*nbs] = gbs1;
	    xbs[*nbs] = xbs1;
	}
    }
    return 0;
} /* s5zswp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s5Zswp */
/* Subroutine */ int symmlq_(integer *n, doublereal *b, doublereal *r1, 
	doublereal *r2, doublereal *v, doublereal *w, doublereal *x, 
	doublereal *y, S_fp aprod, S_fp msolve, logical *checka, logical *
	goodb, logical *precon, doublereal *shift, integer *nout, integer *
	itnlim, doublereal *rtol, integer *istop, integer *itn, doublereal *
	anorm, doublereal *acond, doublereal *rnorm, doublereal *ynorm, U_fp 
	hprod, U_fp hprod1, integer *ne, integer *nloca, integer *loca, 
	integer *inda, doublereal *acol, char *cu, integer *lencu, integer *
	iu, integer *leniu, doublereal *ru, integer *lenru, char *cw, integer 
	*lencw, integer *iw, integer *leniw, doublereal *rw, integer *lenrw, 
	ftnlen cu_len, ftnlen cw_len)
{
    /* Initialized data */

    static char enter[16] = " Enter SYMMLQ.  ";
    static char exit[16] = " Exit  SYMMLQ.  ";
    static char msg[52*11] = "beta2 = 0.  If M = I, b and x are eigenvectors"
	    " of A " "beta1 = 0.  The exact solution is  x = 0            " 
	    "Requested accuracy achieved, as determined by rtol  " "Reasonab"
	    "le accuracy achieved, given eps             " "x has converged t"
	    "o an eigenvector                   " "Acond has exceeded 0.1/eps"
	    "                          " "The iteration limit was reached    "
	    "                 " "Aprod  does not define a symmetric matrix   "
	    "        " "Msolve does not define a symmetric matrix           " 
	    "Msolve does not define a pos-def preconditioner     " "x is not"
	    " a descent direction                        ";

    /* Format strings */
    static char fmt_1000[] = "(//1p,a,5x,\002Solution of symmetric   Ax = "
	    "b\002,7x,\002bnorm  =\002,e10.2/\002 n      =\002,i7,5x,\002chec"
	    "kA =\002,l4,12x,\002goodb  =\002,l4,7x,\002precon =\002,l4/\002 "
	    "itnlim =\002,i7,5x,\002rtol   =\002,e11.2,5x,\002shift  =\002,e2"
	    "3.14)";
    static char fmt_1100[] = "(/1p,\002 beta1  =\002,e10.2,3x,\002alpha1 "
	    "=\002,e10.2/\002 (v1,v2) before and after \002,e14.2/\002 local "
	    "reorthogonalization\002,e14.2)";
    static char fmt_1200[] = "(//5x,\002itn\002,7x,\002x1(cg)\002,10x,\002no"
	    "rm(r)\002,5x,\002bstep\002,7x,\002norm(A)\002,3x,\002cond(A)\002)"
	    ;
    static char fmt_1300[] = "(1p,i8,e19.10,e11.2,e14.5,2e10.2)";
    static char fmt_1500[] = "(1x)";
    static char fmt_2000[] = "(/1p,a,6x,\002istop =\002,i3,15x,\002itn   "
	    "=\002,i8/a,6x,\002Anorm =\002,e12.4,6x,\002Acond =\002,e12.4/a,6"
	    "x,\002rnorm =\002,e12.4,6x,\002ynorm =\002,e12.4)";
    static char fmt_3000[] = "(a,6x,a)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s, t, z__, b1, cs, sn, eps, x1cg, rhs1, rhs2, alfa, 
	    diag, dbar, beta, gbar, oldb, epsa;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal gmin, gmax, zbar, epsr, epsx, beta1;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal gamma;
    extern /* Subroutine */ int dload_(integer *, doublereal *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *);
    static doublereal delta, denom, bnorm, bstep, epsln;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal tnorm, ynorm2, cgnorm, snprod, lqnorm, qrnorm;

    /* Fortran I/O blocks */
    static cilist io___256 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___266 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___285 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___286 = { 0, 0, 0, fmt_1300, 0 };
    static cilist io___287 = { 0, 0, 0, fmt_1500, 0 };
    static cilist io___294 = { 0, 0, 0, fmt_2000, 0 };
    static cilist io___295 = { 0, 0, 0, fmt_3000, 0 };


/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/*     ------------------------------------------------------------------ */

/*     SYMMLQ  is designed to solve the system of linear equations */

/*                Ax = b */

/*     where A is an n by n symmetric matrix and b is a given vector. */
/*     The matrix A is not required to be positive definite. */
/*     (If A is known to be definite, the method of conjugate gradients */
/*     might be preferred, since it will require about the same number of */
/*     iterations as SYMMLQ but slightly less work per iteration.) */


/*     The matrix A is intended to be large and sparse.  It is accessed */
/*     by means of a subroutine call of the form */

/*                call Aprod ( n, x, y ) */

/*     which must return the product y = Ax for any given vector x. */


/*     More generally, SYMMLQ is designed to solve the system */

/*                (A - shift*I) x = b */

/*     where  shift  is a specified scalar value.  If  shift  and  b */
/*     are suitably chosen, the computed vector x may approximate an */
/*     (unnormalized) eigenvector of A, as in the methods of */
/*     inverse iteration and/or Rayleigh-quotient iteration. */
/*     Again, the matrix (A - shift*I) need not be positive definite. */
/*     The work per iteration is very slightly less if  shift = 0. */


/*     A further option is that of preconditioning, which may reduce */
/*     the number of iterations required.  If M = C C' is a positive */
/*     definite matrix that is known to approximate  (A - shift*I) */
/*     in some sense, and if systems of the form  My = x  can be */
/*     solved efficiently, the parameters precon and Msolve may be */
/*     used (see below).  When  precon = .true., SYMMLQ will */
/*     implicitly solve the system of equations */

/*             P (A - shift*I) P' xbar  =  P b, */

/*     i.e.                  Abar xbar  =  bbar */
/*     where                         P  =  C**(-1), */
/*                                Abar  =  P (A - shift*I) P', */
/*                                bbar  =  P b, */

/*     and return the solution       x  =  P' xbar. */
/*     The associated residual is rbar  =  bbar - Abar xbar */
/*                                      =  P (b - (A - shift*I)x) */
/*                                      =  P r. */

/*     In the discussion below, eps refers to the machine precision. */
/*     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16 */
/*     for IBM mainframes and IEEE double-precision arithmetic. */

/*     Parameters */
/*     ---------- */

/*     n       input      The dimension of the matrix A. */

/*     b(n)    input      The rhs vector b. */

/*     r1(n)   workspace */
/*     r2(n)   workspace */
/*     v(n)    workspace */
/*     w(n)    workspace */

/*     x(n)    output     Returns the computed solution  x. */

/*     y(n)    workspace */

/*     Aprod   external   A subroutine defining the matrix A. */
/*                        For a given vector x, the statement */

/*                              call Aprod ( n, x, y ) */

/*                        must return the product y = Ax */
/*                        without altering the vector x. */

/*     Msolve  external   An optional subroutine defining a */
/*                        preconditioning matrix M, which should */
/*                        approximate (A - shift*I) in some sense. */
/*                        M must be positive definite. */
/*                        For a given vector x, the statement */

/*                              call Msolve( n, x, y ) */

/*                        must solve the linear system My = x */
/*                        without altering the vector x. */

/*                        In general, M should be chosen so that Abar has */
/*                        clustered eigenvalues.  For example, */
/*                        if A is positive definite, Abar would ideally */
/*                        be close to a multiple of I. */
/*                        If A or A - shift*I is indefinite, Abar might */
/*                        be close to a multiple of diag( I  -I ). */

/*                        NOTE.  The program calling SYMMLQ must declare */
/*                        Aprod and Msolve to be external. */

/*     checkA  input      If checkA = .true., an extra call of Aprod will */
/*                        be used to check if A is symmetric.  Also, */
/*                        if precon = .true., an extra call of Msolve */
/*                        will be used to check if M is symmetric. */

/*     goodb   input      Usually, goodb should be .false. */
/*                        If x is expected to contain a large multiple of */
/*                        b (as in Rayleigh-quotient iteration), */
/*                        better precision may result if goodb = .true. */
/*                        See Lewis (1977) below. */
/*                        When goodb = .true., an extra call to Msolve */
/*                        is required. */

/*     precon  input      If precon = .true., preconditioning will */
/*                        be invoked.  Otherwise, subroutine Msolve */
/*                        will not be referenced; in this case the */
/*                        actual parameter corresponding to Msolve may */
/*                        be the same as that corresponding to Aprod. */

/*     shift   input      Should be zero if the system Ax = b is to be */
/*                        solved.  Otherwise, it could be an */
/*                        approximation to an eigenvalue of A, such as */
/*                        the Rayleigh quotient b'Ab / (b'b) */
/*                        corresponding to the vector b. */
/*                        If b is sufficiently like an eigenvector */
/*                        corresponding to an eigenvalue near shift, */
/*                        then the computed x may have very large */
/*                        components.  When normalized, x may be */
/*                        closer to an eigenvector than b. */

/*     nout    input      A file number. */
/*                        If nout .gt. 0, a summary of the iterations */
/*                        will be printed on unit nout. */

/*     itnlim  input      An upper limit on the number of iterations. */

/*     rtol    input      A user-specified tolerance.  SYMMLQ terminates */
/*                        if it appears that norm(rbar) is smaller than */
/*                              rtol * norm(Abar) * norm(xbar), */
/*                        where rbar is the transformed residual vector, */
/*                              rbar = bbar - Abar xbar. */

/*                        If shift = 0 and precon = .false., SYMMLQ */
/*                        terminates if norm(b - A*x) is smaller than */
/*                              rtol * norm(A) * norm(x). */

/*     iStop   output     An integer giving the reason for termination... */

/*              -1        beta2 = 0 in the Lanczos iteration; i.e. the */
/*                        second Lanczos vector is zero.  This means the */
/*                        rhs is very special. */
/*                        If there is no preconditioner, b is an */
/*                        eigenvector of A. */
/*                        Otherwise (if precon is true), let My = b. */
/*                        If shift is zero, y is a solution of the */
/*                        generalized eigenvalue problem Ay = lambda My, */
/*                        with lambda = alpha1 from the Lanczos vectors. */

/*                        In general, (A - shift*I)x = b */
/*                        has the solution         x = (1/alpha1) y */
/*                        where My = b. */

/*               0        b = 0, so the exact solution is x = 0. */
/*                        No iterations were performed. */

/*               1        Norm(rbar) appears to be less than */
/*                        the value  rtol * norm(Abar) * norm(xbar). */
/*                        The solution in  x  should be acceptable. */

/*               2        Norm(rbar) appears to be less than */
/*                        the value  eps * norm(Abar) * norm(xbar). */
/*                        This means that the residual is as small as */
/*                        seems reasonable on this machine. */

/*               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps, */
/*                        which should indicate that x has essentially */
/*                        converged to an eigenvector of A */
/*                        corresponding to the eigenvalue shift. */

/*               4        Acond (see below) has exceeded 0.1/eps, so */
/*                        the matrix Abar must be very ill-conditioned. */
/*                        x may not contain an acceptable solution. */

/*               5        The iteration limit was reached before any of */
/*                        the previous criteria were satisfied. */

/*               6        The matrix defined by Aprod does not appear */
/*                        to be symmetric. */
/*                        For certain vectors y = Av and r = Ay, the */
/*                        products y'y and r'v differ significantly. */

/*               7        The matrix defined by Msolve does not appear */
/*                        to be symmetric. */
/*                        For vectors satisfying My = v and Mr = y, the */
/*                        products y'y and r'v differ significantly. */

/*               8        An inner product of the form  x' M**(-1) x */
/*                        was not positive, so the preconditioning matrix */
/*                        M does not appear to be positive definite. */

/*                        If iStop .ge. 5, the final x may not be an */
/*                        acceptable solution. */

/*     itn     output     The number of iterations performed. */

/*     Anorm   output     An estimate of the norm of the matrix operator */
/*                        Abar = P (A - shift*I) P',   where P = C**(-1). */

/*     Acond   output     An estimate of the condition of Abar above. */
/*                        This will usually be a substantial */
/*                        under-estimate of the true condition. */

/*     rnorm   output     An estimate of the norm of the final */
/*                        transformed residual vector, */
/*                           P (b  -  (A - shift*I) x). */

/*     ynorm   output     An estimate of the norm of xbar. */
/*                        This is sqrt( x'Mx ).  If precon is false, */
/*                        ynorm is an estimate of norm(x). */



/*     To change precision */
/*     ------------------- */

/*     Alter the words */
/*            double precision, */
/*            daxpy, dcopy, ddot, dnrm2 */
/*     to their single or double equivalents. */
/*     ------------------------------------------------------------------ */


/*     This routine is an implementation of the algorithm described in */
/*     the following references: */

/*     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite */
/*          Systems of Linear Equations, */
/*          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629. */

/*     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems, */
/*          Report STAN-CS-77-595, Computer Science Department, */
/*          Stanford University, Stanford, California, March 1977. */

/*     Applications of SYMMLQ and the theory of preconditioning */
/*     are described in the following references: */

/*     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient */
/*          Type Methods to Eigenvalue Calculations, */
/*          in R. Vichnevetsky and R.S. Steplman (editors), */
/*          Advances in Computer Methods for Partial Differential */
/*          Equations -- III, IMACS, 1979, 167-173. */

/*     D.B. Szyld,  A Two-level Iterative Method for Large Sparse */
/*          Generalized Eigenvalue Calculations, */
/*          Ph. D. dissertation, Department of Mathematics, */
/*          New York University, New York, October 1983. */

/*     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders, */
/*          Preconditioners for indefinite systems arising in */
/*          optimization, SIMAX 13, 1, 292--311, January 1992. */
/*          (SIAM J. on Matrix Analysis and Applications) */
/*     ------------------------------------------------------------------ */


/*     SYMMLQ development: */
/*            1972: First version. */
/*            1975: John Lewis recommended modifications to help with */
/*                  inverse iteration: */
/*                  1. Reorthogonalize v1 and v2. */
/*                  2. Regard the solution as x = x1  +  bstep * b, */
/*                     with x1 and bstep accumulated separately */
/*                     and bstep * b added at the end. */
/*                     (In inverse iteration, b might be close to the */
/*                     required x already, so x1 may be a lot smaller */
/*                     than the multiple of b.) */
/*            1978: Daniel Szyld and Olof Widlund implemented the first */
/*                  form of preconditioning. */
/*                  This required both a solve and a multiply with M. */
/*            1979: Implemented present method for preconditioning. */
/*                  This requires only a solve with M. */
/*            1984: Sven Hammarling noted corrections to tnorm and x1lq. */
/*                  SYMMLQ added to NAG Fortran Library. */
/*     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib. */
/*     16 Feb 1989: First F77 version. */

/*     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure) */
/*                  if Abar = const*I.  iStop = -1 added for this case. */

/*     01 Mar 1989: Hans Mittelmann observed premature termination on */
/*                  ( 1  1  1 )     (   )                   ( 1  1    ) */
/*                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ). */
/*                  ( 1     1 )     (   )                   (    1  1 ) */
/*                  T2 is exactly singular, so estimating cond(A) from */
/*                  the diagonals of Lbar is unsafe.  We now use */
/*                  L       or  Lbar         depending on whether */
/*                  lqnorm  or  cgnorm       is least. */

/*     03 Mar 1989: eps computed internally instead of coming in as a */
/*                  parameter. */
/*     07 Jun 1989: ncheck added as a parameter to say if A and M */
/*                  should be checked for symmetry. */
/*                  Later changed to checkA (see below). */
/*     20 Nov 1990: goodb added as a parameter to make Lewis's changes */
/*                  an option.  Usually b is NOT much like x.  Setting */
/*                  goodb = .false. saves a call to Msolve at the end. */
/*     20 Nov 1990: Residual not computed exactly at end, to save time */
/*                  when only one or two iterations are required */
/*                  (e.g. if the preconditioner is very good). */
/*                  Beware, if precon is true, rnorm estimates the */
/*                  residual of the preconditioned system, not Ax = b. */
/*     04 Sep 1991: Parameter list changed and reordered. */
/*                  integer ncheck is now logical checkA. */
/*     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti */
/*                  showed that beta2 = 0 (iStop = -1) means that */
/*                  b is an eigenvector when M = I. */
/*                  More complicated if there is a preconditioner; */
/*                  not clear yet how to describe it. */
/*     20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero. */
/*                  Need to estimate Anorm from column of Tk. */
/*     18 Nov 2001: bnorm printed in first line, since rnorm at itn = 0 */
/*                  isn't bnorm as we might expect -- it's cgnorm after */
/*                  the proverbial step to the CG point! */
/*     02 Aug 2003: eps grabbed from rw (special to SNOPT). */
/*     23 Dec 2003: Normal backward error exit (iStop = 1) terminates */
/*                  too early when rtol = 0.01 say (rather large). */
/*                  Keep it with   rtol = eps  (iStop = 2) */
/*                  but replace normal test with ||r|| < rtol*||b|| */
/*                  like most other CG solvers.  This fits better with */
/*                  the inexact Newton context.  Note that it is correct */
/*                  only with precon = .false. */

/*     08 Apr 2004: Always move to the CG point. */

/*     Michael A. Saunders                   na.msaunders@na-net.ornl.gov */
/*     Systems Optimization Laboratory       saunders@stanford.edu */
/*     Dept of Management Science and Engineering */
/*     Stanford University */
/*     Stanford, CA 94305-4026                             (650) 723-1875 */
/*     ------------------------------------------------------------------ */


/*     Subroutines and functions */

/*     USER       Aprod, Msolve */
/*     BLAS       daxpy, dcopy, ddot , dnrm2 */


/*     Intrinsics and local variables */
/* Added for SNOPT */
/* Added for SNOPT */
    /* Parameter adjustments */
    --y;
    --x;
    --w;
    --v;
    --r2;
    --r1;
    --b;
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
/*     ------------------------------------------------------------------ */
    eps = rw[1];
/*     Print heading and initialize. */
/* unit round-off.  IEEE DP  2.22e-16 */
    bnorm = dnrm2_(n, &b[1], &c__1);
    beta1 = bnorm;
    if (*nout > 0) {
	io___256.ciunit = *nout;
	s_wsfe(&io___256);
	do_fio(&c__1, enter, (ftnlen)16);
	do_fio(&c__1, (char *)&beta1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*checka), (ftnlen)sizeof(logical));
	do_fio(&c__1, (char *)&(*goodb), (ftnlen)sizeof(logical));
	do_fio(&c__1, (char *)&(*precon), (ftnlen)sizeof(logical));
	do_fio(&c__1, (char *)&(*itnlim), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*rtol), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*shift), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    *istop = 0;
    *itn = 0;
    *anorm = 0.;
    *acond = 0.;
    *rnorm = 0.;
    *ynorm = 0.;
    dload_(n, &c_b5, &x[1], &c__1);
/*     Set up y for the first Lanczos vector v1. */
/*     y is really beta1 * P * v1  where  P = C**(-1). */
/*     y and beta1 will be zero if b = 0. */
    dcopy_(n, &b[1], &c__1, &y[1], &c__1);
    dcopy_(n, &b[1], &c__1, &r1[1], &c__1);
    if (*precon) {
	(*msolve)(n, &r1[1], &y[1]);
    }
    if (*goodb) {
	b1 = y[1];
    } else {
	b1 = 0.;
    }
    beta1 = ddot_(n, &r1[1], &c__1, &y[1], &c__1);
/*     See if Msolve is symmetric. */
    if (*checka && *precon) {
	(*msolve)(n, &y[1], &r2[1]);
	s = ddot_(n, &y[1], &c__1, &y[1], &c__1);
	t = ddot_(n, &r1[1], &c__1, &r2[1], &c__1);
	z__ = (d__1 = s - t, abs(d__1));
	epsa = (s + eps) * pow_dd(&eps, &c_b184);
	if (z__ > epsa) {
	    *istop = 7;
	    goto L900;
	}
    }
/*     Test for an indefinite preconditioner. */
    if (beta1 < 0.) {
	*istop = 8;
	goto L900;
    }
/*     If b = 0 exactly, stop with x = 0. */
    if (beta1 == 0.) {
	goto L900;
    }
/*     Here and later, v is really P * (the Lanczos v). */
    beta1 = sqrt(beta1);
    s = 1. / beta1;
    dcopy_(n, &y[1], &c__1, &v[1], &c__1);
    dscal_(n, &s, &v[1], &c__1);
    (*aprod)(n, &v[1], &y[1], (U_fp)hprod, (U_fp)hprod1, ne, nloca, &loca[1], 
	    &inda[1], &acol[1], cu + 8, lencu, &iu[1], leniu, &ru[1], lenru, 
	    cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8, (ftnlen)8)
	    ;
/*     See if Aprod  is symmetric. */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
    if (*checka) {
	(*aprod)(n, &y[1], &r2[1], (U_fp)hprod, (U_fp)hprod1, ne, nloca, &
		loca[1], &inda[1], &acol[1], cu + 8, lencu, &iu[1], leniu, &
		ru[1], lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (
		ftnlen)8, (ftnlen)8);
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
	s = ddot_(n, &y[1], &c__1, &y[1], &c__1);
	t = ddot_(n, &v[1], &c__1, &r2[1], &c__1);
	z__ = (d__1 = s - t, abs(d__1));
	epsa = (s + eps) * pow_dd(&eps, &c_b184);
	if (z__ > epsa) {
	    *istop = 6;
	    goto L900;
	}
    }
/*     Set up y for the second Lanczos vector. */
/*     Again, y is beta * P * v2  where  P = C**(-1). */
/*     y and beta will be zero or very small if b is an eigenvector. */
    d__1 = -(*shift);
    daxpy_(n, &d__1, &v[1], &c__1, &y[1], &c__1);
    alfa = ddot_(n, &v[1], &c__1, &y[1], &c__1);
    d__1 = -alfa / beta1;
    daxpy_(n, &d__1, &r1[1], &c__1, &y[1], &c__1);
/*     Make sure  r2  will be orthogonal to the first  v. */
    z__ = ddot_(n, &v[1], &c__1, &y[1], &c__1);
    s = ddot_(n, &v[1], &c__1, &v[1], &c__1);
    d__1 = -z__ / s;
    daxpy_(n, &d__1, &v[1], &c__1, &y[1], &c__1);
    dcopy_(n, &y[1], &c__1, &r2[1], &c__1);
    if (*precon) {
	(*msolve)(n, &r2[1], &y[1]);
    }
    oldb = beta1;
    beta = ddot_(n, &r2[1], &c__1, &y[1], &c__1);
    if (beta < 0.) {
	*istop = 8;
	goto L900;
    }
/*     Cause termination (later) if beta is essentially zero. */
    beta = sqrt(beta);
    if (beta <= eps) {
	*istop = -1;
    }
/*     See if the local reorthogonalization achieved anything. */
    denom = sqrt(s) * dnrm2_(n, &r2[1], &c__1) + eps;
    s = z__ / denom;
    t = ddot_(n, &v[1], &c__1, &r2[1], &c__1) / denom;
    if (*nout > 0 && *goodb) {
	io___266.ciunit = *nout;
	s_wsfe(&io___266);
	do_fio(&c__1, (char *)&beta1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&alfa, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     Initialize other quantities. */
    cgnorm = beta1;
    gbar = alfa;
    dbar = beta;
    rhs1 = beta1;
    rhs2 = 0.;
    bstep = 0.;
    snprod = 1.;
/* Computing 2nd power */
    d__1 = alfa;
/* Computing 2nd power */
    d__2 = beta;
    tnorm = d__1 * d__1 + d__2 * d__2;
    ynorm2 = 0.;
    gmax = abs(alfa) + eps;
    gmin = gmax;
    if (*goodb) {
	dload_(n, &c_b5, &w[1], &c__1);
    } else {
	dcopy_(n, &v[1], &c__1, &w[1], &c__1);
    }
/*     ------------------------------------------------------------------ */
/*     Main iteration loop. */
/*     ------------------------------------------------------------------ */
/* +    Repeat                                          Until iStop .ne. 0 */
/*        Estimate various norms and test for convergence. */
L100:
    *anorm = sqrt(tnorm);
    *ynorm = sqrt(ynorm2);
    epsa = *anorm * eps;
    epsx = *anorm * *ynorm * eps;
    epsr = *anorm * *ynorm * *rtol;
    diag = gbar;
    if (diag == 0.) {
	diag = epsa;
    }
/* Computing 2nd power */
    d__1 = rhs1;
/* Computing 2nd power */
    d__2 = rhs2;
    lqnorm = sqrt(d__1 * d__1 + d__2 * d__2);
    qrnorm = snprod * beta1;
    cgnorm = qrnorm * beta / abs(diag);
/*        Estimate  cond(A). */
/*        In this version we look at the diagonals of  L  in the */
/*        factorization of the tridiagonal matrix,  T = L*Q. */
/*        Sometimes, T(k) can be misleadingly ill-conditioned when */
/*        T(k+1) is not, so we must be careful not to overestimate Acond. */
    if (lqnorm <= cgnorm) {
	*acond = gmax / gmin;
    } else {
/* Computing MIN */
	d__1 = gmin, d__2 = abs(diag);
	denom = min(d__1,d__2);
	*acond = gmax / denom;
    }
/*        See if any of the stopping criteria are satisfied. */
/*        In rare cases, iStop is already -1 from above (Abar = const * I). */
    if (*istop == 0) {
	if (*itn >= *itnlim) {
	    *istop = 5;
	}
	if (*acond >= .1 / eps) {
	    *istop = 4;
	}
	if (epsx >= beta1) {
	    *istop = 3;
	}
	if (cgnorm <= epsx) {
	    *istop = 2;
	}
	if (cgnorm <= bnorm * *rtol) {
	    *istop = 1;
	}
/*           if (cgnorm .le. epsr      ) iStop = 1 ! for inexact Newton. */
/* Replaces next line */
    }
/*        =============================================================== */
/*        See if it is time to print something. */
    if ((*n <= 40 || *itn <= 10 || *itn >= *itnlim - 10 || *itn % 10 == 0 || 
	    cgnorm <= epsx * 10. || cgnorm <= epsr * 10. || *acond >= .01 / 
	    eps || *istop != 0) && *nout > 0) {
/*           Print a line for this iteration. */
	zbar = rhs1 / diag;
	z__ = (snprod * zbar + bstep) / beta1;
/*           x1lq   = x(1)  +  b1 * bstep / beta1 */
	x1cg = x[1] + w[1] * zbar + b1 * z__;
	if (*itn == 0) {
	    io___285.ciunit = *nout;
	    s_wsfe(&io___285);
	    e_wsfe();
	}
	io___286.ciunit = *nout;
	s_wsfe(&io___286);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x1cg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cgnorm, (ftnlen)sizeof(doublereal));
	d__1 = bstep / beta1;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*anorm), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*acond), (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (*itn % 10 == 0) {
	    io___287.ciunit = *nout;
	    s_wsfe(&io___287);
	    e_wsfe();
	}
    }
/*        =============================================================== */
/*        Obtain the current Lanczos vector  v = (1 / beta)*y */
/*        and set up  y  for the next iteration. */
    if (*istop == 0) {
	s = 1. / beta;
	dcopy_(n, &y[1], &c__1, &v[1], &c__1);
	dscal_(n, &s, &v[1], &c__1);
	(*aprod)(n, &v[1], &y[1], (U_fp)hprod, (U_fp)hprod1, ne, nloca, &loca[
		1], &inda[1], &acol[1], cu + 8, lencu, &iu[1], leniu, &ru[1], 
		lenru, cw + 8, lencw, &iw[1], leniw, &rw[1], lenrw, (ftnlen)8,
		 (ftnlen)8);
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
/* Added for SNOPT */
	d__1 = -(*shift);
	daxpy_(n, &d__1, &v[1], &c__1, &y[1], &c__1);
	d__1 = -beta / oldb;
	daxpy_(n, &d__1, &r1[1], &c__1, &y[1], &c__1);
	alfa = ddot_(n, &v[1], &c__1, &y[1], &c__1);
	d__1 = -alfa / beta;
	daxpy_(n, &d__1, &r2[1], &c__1, &y[1], &c__1);
	dcopy_(n, &r2[1], &c__1, &r1[1], &c__1);
	dcopy_(n, &y[1], &c__1, &r2[1], &c__1);
	if (*precon) {
	    (*msolve)(n, &r2[1], &y[1]);
	}
	oldb = beta;
	beta = ddot_(n, &r2[1], &c__1, &y[1], &c__1);
	if (beta < 0.) {
	    *istop = 6;
	    goto L800;
	}
	beta = sqrt(beta);
/* Computing 2nd power */
	d__1 = alfa;
/* Computing 2nd power */
	d__2 = oldb;
/* Computing 2nd power */
	d__3 = beta;
	tnorm = tnorm + d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/*           Compute the next plane rotation for  Q. */
/* Computing 2nd power */
	d__1 = gbar;
/* Computing 2nd power */
	d__2 = oldb;
	gamma = sqrt(d__1 * d__1 + d__2 * d__2);
	cs = gbar / gamma;
	sn = oldb / gamma;
	delta = cs * dbar + sn * alfa;
	gbar = sn * dbar - cs * alfa;
	epsln = sn * beta;
	dbar = -cs * beta;
/*           Update  x. */
	z__ = rhs1 / gamma;
	s = z__ * cs;
	t = z__ * sn;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = w[i__] * s + v[i__] * t + x[i__];
	    w[i__] = w[i__] * sn - v[i__] * cs;
	}
/*           Accumulate the step along the direction  b, */
/*           and go round again. */
	bstep = snprod * cs * z__ + bstep;
	snprod *= sn;
	gmax = max(gmax,gamma);
	gmin = min(gmin,gamma);
/* Computing 2nd power */
	d__1 = z__;
	ynorm2 = d__1 * d__1 + ynorm2;
	rhs1 = rhs2 - delta * z__;
	rhs2 = -epsln * z__;
	++(*itn);
    }
/* +    until    (iStop .ne. 0) */
    if (! (*istop != 0)) {
	goto L100;
    }
/*     ------------------------------------------------------------------ */
/*     End of main iteration loop. */
/*     ------------------------------------------------------------------ */
/*     Move to the CG point if it seems better. */
/*     In this version of SYMMLQ, the convergence tests involve */
/*     only cgnorm, so we're unlikely to stop at an LQ point, */
/*     EXCEPT if the iteration limit interferes. */

/*     April 8, 2004. Always move to the CG point for SNOPT application */
L800:
    zbar = rhs1 / diag;
    bstep = snprod * zbar + bstep;
/* Computing 2nd power */
    d__1 = zbar;
    *ynorm = sqrt(ynorm2 + d__1 * d__1);
    *rnorm = cgnorm;
    daxpy_(n, &zbar, &w[1], &c__1, &x[1], &c__1);
    if (*goodb) {
/*        Add the step along  b. */
	bstep /= beta1;
	dcopy_(n, &b[1], &c__1, &y[1], &c__1);
	if (*precon) {
	    (*msolve)(n, &b[1], &y[1]);
	}
	daxpy_(n, &bstep, &y[1], &c__1, &x[1], &c__1);
    }
/*     ================================================================== */
/*     Display final status. */
/*     ================================================================== */
L900:
    if (*nout > 0) {
	io___294.ciunit = *nout;
	s_wsfe(&io___294);
	do_fio(&c__1, exit, (ftnlen)16);
	do_fio(&c__1, (char *)&(*istop), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, exit, (ftnlen)16);
	do_fio(&c__1, (char *)&(*anorm), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*acond), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, exit, (ftnlen)16);
	do_fio(&c__1, (char *)&(*rnorm), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*ynorm), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___295.ciunit = *nout;
	s_wsfe(&io___295);
	do_fio(&c__1, exit, (ftnlen)16);
	do_fio(&c__1, msg + (*istop + 1) * 52, (ftnlen)52);
	e_wsfe();
    }
    return 0;
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     end of SYMMLQ */
} /* symmlq_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s5msolv_(integer *n, doublereal *x, doublereal *y)
{
/*     ================================================================== */
/*     s5Msolv  dummy Msolve for SYMMLQ. */

/*     04 Dec 2004: First version of s5Msolv. */
/*     04 Dec 2004: Current version. */
/*     ================================================================== */
/*     Relax */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    return 0;
} /* s5msolv_ */

