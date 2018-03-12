/* ./src/sn25bfac.f -- translated by f2c (version 20100827).
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
static integer c__11 = 11;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__23 = 23;
static integer c__40 = 40;
static integer c__14 = 14;
static integer c__13 = 13;
static integer c__3 = 3;
static integer c__6 = 6;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn25bfac.f */

/*     s2Bfac   s2Bkbs   s2Bmap   s2newB   s2BLU     s2Bmod   s2Bmod2 */
/*     s2Bsol   s2sb     s2sing   s2tols   s2tryLU */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s2bfac_(integer *iexit, integer *typelu, logical *needlu,
	 logical *newlu, logical *newb, integer *iobj, integer *itn, integer *
	lprint, integer *lureq, integer *m, integer *mbs, integer *n, integer 
	*nb, integer *nnl, integer *ns, integer *nswap, integer *ne, integer *
	nloca, integer *loca, integer *inda, doublereal *acol, integer *kbs, 
	integer *hs, doublereal *bl, doublereal *bu, doublereal *blbs, 
	doublereal *bubs, integer *nrhs0, integer *nrhs, doublereal *rhs, 
	doublereal *x, doublereal *xbs, integer *iy, integer *iy1, doublereal 
	*y, doublereal *y1, integer *iw, integer *leniw, doublereal *rw, 
	integer *lenrw)
{
    /* Initialized data */

    static char lutype[2*4] = "B " "BR" "BS" "BT";

    /* Format strings */
    static char fmt_1005[] = "(\002 Factor\002,i7,\002  Demand\002,i7,\002  "
	    "Itn\002,i11)";
    static char fmt_2000[] = "(\002 Itn\002,i7,\002: \002,a,\002 factoriz"
	    "e\002)";
    static char fmt_3000[] = "(\002 BS Factorize.   nSwap = \002,i6)";
    static char fmt_1100[] = "(\002 Nonlin\002,i7,\002  Linear\002,i7,\002  "
	    "Slacks\002,i8,2x,a,\002 factorize\002)";
    static char fmt_7000[] = "(i7,g17.8)";
    static char fmt_9032[] = "(\002 XXX Singular basis after \002,i5,\002 fa"
	    "ctorization attempts\002)";
    static char fmt_9080[] = "(24x,\002        Current    Recommended\002)";
    static char fmt_9081[] = "(\002 Total integer workspace\002,2i15)";
    static char fmt_9082[] = "(\002 Total real    workspace\002,2i15)";
    static char fmt_9101[] = "(\002 XXX Wrong no. of basic variables:\002,i8)"
	    ;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer j, k, ip, iq, lua, lin, nbs;
    static char str[72];
    static doublereal eps2;
    extern /* Subroutine */ int s2sb_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *);
    static integer indc;
    static logical badx, bigu;
    static integer indr, locr, cols, more, newi, newr;
    static doublereal utol;
    static integer ntry;
    extern /* Subroutine */ int s2blu_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, integer *);
    static logical brfac, bsfac, btfac;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer markc;
    static logical needx;
    static integer btask, markr, ntask, ptask;
    static logical prnt10;
    static integer maxiw, maxrw;
    extern /* Subroutine */ int s2bkbs_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), s2dmat_(integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *), s2newb_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *), s2sing_(integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), s2tols_(integer *, logical *, integer *, integer *, 
	    integer *, doublereal *, integer *), s5setx_(integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static integer nbasic;
    static logical brdone, bsdone;
    static integer lprdbg, nslack, matfil, lenalu, inform__, nonlin;
    static logical singlr;
    static integer maxlui;
    static logical newtol;
    static integer maxlur;
    static doublereal rowerr;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___28 = { 0, str, 0, fmt_1005, 72, 1 };
    static icilist io___33 = { 0, str, 0, fmt_2000, 72, 1 };
    static icilist io___35 = { 0, str, 0, fmt_3000, 72, 1 };
    static icilist io___42 = { 0, str, 0, fmt_1100, 72, 1 };
    static icilist io___43 = { 0, str, 0, fmt_2000, 72, 1 };
    static icilist io___50 = { 0, str, 0, fmt_7000, 72, 1 };
    static icilist io___51 = { 0, str, 0, fmt_9032, 72, 1 };
    static icilist io___55 = { 0, str, 0, fmt_9080, 72, 1 };
    static icilist io___56 = { 0, str, 0, fmt_9081, 72, 1 };
    static icilist io___57 = { 0, str, 0, fmt_9082, 72, 1 };
    static icilist io___58 = { 0, str, 0, fmt_9101, 72, 1 };


/*     ================================================================= */
/*     s2Bfac  computes an acceptable x such that  ( A  -I )*x = b. */
/*     The LU factorization of the basis is computed if necessary. */

/*     If typeLU = B , the usual B = LU is computed. */
/*     If typeLU = BR, we want to check rank(B) with a special B = LU */
/*                     (TRP, tight tols) before computing a normal LU. */
/*     If typeLU = BS, there are some superbasics and we want to */
/*                     choose a good basis from the columns of (B S). */
/*                     We first factorize (B S)' to obtain a new B. */
/*                     Then B = LU is computed as usual. */
/*     If typeLU = BT, we should TRY 'B ' first and go back to 'BS' */
/*                     only if B seems ill-conditioned. */

/*     15 Nov 1991: First version based on Minos routine m2bfac. */
/*     29 Oct 1993: typeLU options implemented. */
/*                  nSwap returns the number of (B S) changes. */
/*     22 Apr 1994: Retry with reduced LU Factor tol */
/*                  if s2BLU says there was large growth in U. */
/*                  'BT' option implemented to save R more often each */
/*                  major iteration. */
/*     02 Apr 1996: kObj added to mark position of Obj slack in B. */
/*     14 Jul 1997: Thread-safe version. */
/*     10 Mar 2001: BR option implemented to help with CUTE prob TRAINH. */
/*                  BTfac, BSfac or Bfac tried first. */
/*                  BRfac used only if growth detected in U or x. */
/*     17 Jun 2001: If BSfac is singular, ask for BRfac. */
/*     18 Jun 2001: If B fac is singular, ask for BRfac. */
/*     18 Nov 2001: lPrint added as parameter to s2sing. */
/*     31 Jul 2003: snEXIT and snPRNT adopted. */
/*     11 Apr 2004: kBS re-used to reduce the value of nswap */
/*     20 Dec 2004: If B continues to be singular after BRdone, */
/*                  tighten the LU tols each time. */
/*     21 Dec 2004: Insert slacks before tightening tols (for next LU). */
/*     02 Jul 2005: Calls to s2kbs and s2bs added. */
/*     04 Jul 2010: Bigger values of miniw and minrw saved in workspace. */
/*     ================================================================= */
/*     ------------------------------------------------------------------ */
/*     LUSOL arguments. */
/* # of singularities in w(*) */
/* minimum recommended lenaLU */

/* minimum diagonal in  U. */
/* saved smallest U diagonal */
/* xBS(kObj) is the obj. slack */
/* # of LU factorizations */
/* # consecutive `B' facts */
/* itns since last factorize */
/* number of LU mods */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xbs;
    --bubs;
    --blbs;
    --kbs;
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
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    eps2 = rw[4];
/* eps**(1/2)       IEEE DP  1.49e-08 */
    maxrw = iw[3];
/* end of SNOPT part of rw */
    maxiw = iw[5];
/* end of SNOPT part of iw */
    lprdbg = iw[85];
/* > 0    => private debug print */
    lenalu = iw[213];
/* space allotted for LU factors */
    maxlui = iw[361];
/* max LU nonzeros in iw(*) */
    maxlur = iw[362];
/* max LU nonzeros in rw(*) */
    ip = iw[363];

    iq = iw[364];

    locr = iw[368];

    lua = iw[371];

    indc = iw[373];

    indr = iw[374];

    cols = iw[375];
/* 3 new work arrays for lu1mxr */
    markc = iw[376];
/* passed to lu1fac and lu1fad */
    markr = iw[377];
/* from s2BLU */
    *newb = FALSE_;
    brdone = FALSE_;
    bsdone = FALSE_;
    brfac = FALSE_;
    bsfac = FALSE_;
    btfac = FALSE_;
    nbs = *m + *ns;
    *nswap = 0;
    ntask = 0;
    ptask = *typelu;
    prnt10 = *lprint >= 10;
/*     Initialize Umin and nBFac on first entry. */
/*     nBFac  counts consecutive B factorizations (reset if BS is done). */
/*     Umin   is the smallest diagonal of U after last BS factor. */
    if (iw[210] == 0) {
	rw[190] = 0.;
	iw[211] = 0;
    }
    if (*needlu) {
	if (iw[210] >= 1) {
	    ++iw[210];
	} else {
	    ++iw[210];
	}
	++iw[211];
	if (prnt10) {
	    s_wsfi(&io___28);
	    do_fio(&c__1, (char *)&iw[210], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*lureq), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__11, str, &iw[1], leniw, (ftnlen)72);
	}
	iw[215] = 0;
	iw[216] = 0;
	*lureq = 0;
/*        --------------------------------------------------------------- */
/*        Set local logicals to select the required type of LU. */
/*        We come back to 100 if a BT factorize looks doubtful. */
/*        If BT was requested but we haven't done BS yet, */
/*        might as well do BS now. */
/*        --------------------------------------------------------------- */
	btfac = *typelu == 3 && *ns > 0;
	bsfac = *typelu == 2 && *ns > 0 || btfac && rw[190] == 0.;
	brfac = *typelu == 2 && *ns == 0;
	if (! (btfac || bsfac)) {
	    ptask = 0;
	}
    }
/*     ------------------------------------------------------------------ */
/*     We may come back here to do a BSfac after all. */
/*     ------------------------------------------------------------------ */
/* needLU */
L100:
    *iexit = 0;
    singlr = FALSE_;
    if (bsfac) {
/*        --------------------------------------------------------------- */
/*        Repartition (B S) to get a better B. */
/*        --------------------------------------------------------------- */
	btfac = FALSE_;
	iw[211] = 1;
	ptask = 2;
/* Load the basics,  x first, into kBS. */
	s2bkbs_(&c__0, iobj, m, mbs, n, nb, nnl, &hs[1], &kbs[1], &iw[205], &
		nbasic, &nonlin, &nslack);
	if (nbasic == *m) {
/*           ------------------------------------------------------------ */
/*           We have the right number of basics. */
/*           1. Factorize (B S)'. */
/*           2. Apply the resulting row permutation to the columns */
/*              of (B S). */
/*           ------------------------------------------------------------ */
	    ++ntask;
	    if (prnt10) {
		if (ntask > 1) {
		    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
		}
		s_wsfi(&io___33);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, lutype + (ptask << 1), (ftnlen)2);
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)72);
	    }
	    s2blu_(&inform__, &c__2, lprint, m, n, nb, &nbs, ne, nloca, &loca[
		    1], &inda[1], &acol[1], &kbs[1], &iw[ip], &rw[lua], &iw[
		    indc], &iw[indr], &lenalu, &iy[1], &iy1[1], &y[1], &iw[1],
		     leniw, &rw[1], lenrw);
	    if (inform__ >= 7) {
		*iexit = 80;
/* insufficient storage for the LU. */
		goto L400;
	    } else if (inform__ >= 3) {
		*iexit = 142;
/* Error in basis package */
		goto L400;
	    }
	    singlr = iw[161] > 0;
	    bsdone = TRUE_;
/* Once only! */
	    ptask = 0;
	    s2newb_(m, mbs, nb, ns, &hs[1], &iw[ip], &kbs[1], &iy[1], &iw[
		    locr], nswap);
	    if (*nswap > 0) {
		*newb = TRUE_;
	    }
	    if (prnt10) {
		s_wsfi(&io___35);
		do_fio(&c__1, (char *)&(*nswap), (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)72);
	    }
	    if (singlr && ! brdone) {
		brfac = TRUE_;
		ptask = 1;
	    }
	}
/* m basics */
    }
/* BS */
    needx = TRUE_;
    newtol = TRUE_;
/* +    while (needx  .and.  NewTol) do */
L150:
    if (needx && newtol) {
/*        =============================================================== */
/*        Main loop to obtain a good x. */
/*        typeLU is not used.  (We are always factoring just B.) */
/*        =============================================================== */
	ntry = 0;
	bigu = FALSE_;
/* +       while (needLU  .and.  ntry .le. mtry  .and.  .not. bigU) do */
L200:
	if (*needlu && ntry <= 40 && ! bigu) {
/* ----------------------------------------------------------- */
/* Main loop to find a good  B = LU. */
/* ----------------------------------------------------------- */
/* Load and count the basics in kBS, with the slacks loaded */
/* first.  kObj keeps track of the linear objective. */
	    s2bkbs_(&c__1, iobj, m, mbs, n, nb, nnl, &hs[1], &kbs[1], &iw[205]
		    , &nbasic, &nonlin, &nslack);
	    if (nbasic > *m) {
		*iexit = 141;
/* Too many basic variables */
		goto L400;
	    }
	    if (nbasic < *m) {
/* Too few basics. */
/* Set remaining kBS(k) = 0 for s2sing. */
		i__1 = *m - nbasic;
		iload_(&i__1, &c__0, &kbs[nbasic + 1], &c__1);
	    }
	    lin = nbasic - nslack - nonlin;
	    if (lin < 0) {
		lin = 0;
	    }
/* ----------------------------------------------------------- */
/* Load the basis matrix into the LU arrays and factorize it. */
/* ----------------------------------------------------------- */
	    ++ntry;
	    ++ntask;
	    if (brfac) {
		btask = 1;
	    } else {
		btask = 0;
	    }
	    if (prnt10) {
		if (ntask > 1) {
		    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
		}
		s_wsfi(&io___42);
		do_fio(&c__1, (char *)&nonlin, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&lin, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nslack, (ftnlen)sizeof(integer));
		do_fio(&c__1, lutype + (ptask << 1), (ftnlen)2);
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)72);
	    }
	    if (ptask == 1) {
		s_wsfi(&io___43);
		do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
		do_fio(&c__1, lutype + (ptask << 1), (ftnlen)2);
		e_wsfi();
		snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)72);
	    }
/* !!!! DEBUG */
/*     Report file 91, 92, 93 activates this call */
	    matfil = iw[130];
/* Report file n */
	    if (matfil > 0) {
		s2dmat_(&matfil, n, nb, ne, nloca, &acol[1], &inda[1], &loca[
			1], &hs[1]);
	    }
	    s2blu_(&inform__, &btask, lprint, m, n, nb, &nbs, ne, nloca, &
		    loca[1], &inda[1], &acol[1], &kbs[1], &iw[ip], &rw[lua], &
		    iw[indc], &iw[indr], &lenalu, &iy[1], &iy1[1], &y[1], &iw[
		    1], leniw, &rw[1], lenrw);
	    *needlu = inform__ > 0 || brfac;
	    bigu = inform__ == 2;
	    singlr = iw[161] > 0;
	    if (inform__ >= 7) {
		*iexit = 80;
/* Insufficient storage for the LU. */
	    } else if (inform__ >= 3) {
		*iexit = 142;
/* error in basis package */
	    } else if (inform__ > 0 && ntry > 40) {
		*iexit = 42;
/* Singular basis */
	    }
	    if (*iexit > 0) {
		goto L400;
	    }
	    if (brfac) {
		brfac = FALSE_;
		brdone = TRUE_;
		ptask = 0;
	    }
	    if (bsfac) {
/*              -------------------------------------------------------- */
/*              We started with a BS factorize this time. */
/*              Save the smallest diag of U. */
/*              -------------------------------------------------------- */
		rw[190] = rw[164];
	    } else if (btfac) {
/*              -------------------------------------------------------- */
/*              (We come here only once.) */
/*              See if we should have done a BS factorize after all. */
/*              In this version we do it if any of the following hold: */
/*              1. dUmin (the smallest diagonal of U) is noticeably */
/*                 smaller than Umin (its value at the last BS factor). */
/*              2. dUmin is pretty small anyway. */
/*              3. B was singular. */
/*              nBFac  makes BS increasingly likely the longer we */
/*              keep doing B and not BS. */
/*              -------------------------------------------------------- */
		btfac = FALSE_;
		utol = rw[190] * .1 * iw[211];
		bsfac = rw[164] <= utol || rw[164] <= eps2 || singlr;
		if (bsfac) {
		    *needlu = TRUE_;
		    ptask = 2;
		    goto L100;
		}
	    } else if (*ns == 0) {
		rw[190] = rw[164];
	    }
/* ---------------------------------------------------------- */
/* Deal with singularity. */
/* ---------------------------------------------------------- */
/* BS */
	    if (singlr) {
		if (! brdone) {
		    brfac = TRUE_;
		    ptask = 1;
		    goto L150;
		}
/* Suspect columns are indicated by y(j) <= 0. */
/* Replace them by suitable slacks. */
/* Then check if any superbasic slacks were made basic. */
		s2sing_(lprint, mbs, m, n, nb, &y[1], &iw[ip], &iw[iq], &bl[1]
			, &bu[1], &hs[1], &kbs[1], &x[1], &iw[1], leniw);
		if (*ns > 0) {
		    s2sb_(m, mbs, nb, ns, &nbs, &hs[1], &kbs[1], &xbs[1]);
		}
/* Tighten tols for next LU if possible. */
		s2tols_(&c__1, &newtol, itn, &iw[1], leniw, &rw[1], lenrw);
		*newb = TRUE_;
		ptask = 0;
	    }
/* singlr */
	    goto L200;
	}
/* +       end while needLU */
/*        --------------------------------------------------------------- */
/*        We have a nonsingular B such that B = LU. */
/*        Compute the basic variables and check that  (A -I)*x = b. */
/*        s5setx also loads the basic/superbasic variables in xBS. */
/*        If the row check fails (or U was big earlier), request BRfac. */
/*        --------------------------------------------------------------- */
	s5setx_(&inform__, &c__0, itn, m, n, nb, &nbs, &rowerr, ne, nloca, &
		loca[1], &inda[1], &acol[1], &kbs[1], &xbs[1], nrhs0, nrhs, &
		rhs[1], &x[1], &y[1], &y1[1], &iw[1], leniw, &rw[1], lenrw);
	badx = inform__ > 0;
	needx = bigu || badx;
	if (needx) {
	    *needlu = TRUE_;
	    if (! bsdone) {
		bsfac = *ns > 0;
		if (bsfac) {
		    ptask = 2;
		    goto L100;
		}
	    }
	    if (! brdone) {
		brfac = TRUE_;
		ptask = 1;
	    } else {
		s2tols_(&c__1, &newtol, itn, &iw[1], leniw, &rw[1], lenrw);
		brfac = TRUE_;
		brdone = FALSE_;
		ptask = 1;
	    }
	}
	goto L150;
    }
/* +    end while needx and NewTol */
    if (needx) {
	*iexit = 43;
    }
/*     ================================================================== */
/*     Tidy up */
/*     ------------------------------------------------------------------ */
/* Cannot satisfy the linearized constraints */
L400:
    if (*iexit == 0) {
/*        -------------------------------------------------------------- */
/*        Normal exit. */
/*        Load the basic/superbasic bounds into blBS, buBS. */
/*        -------------------------------------------------------------- */
	i__1 = nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    blbs[k] = bl[j];
	    bubs[k] = bu[j];
	}
	*newlu = ntry > 0;
	if (lprdbg == 100) {
	    snprnt_(&c__11, " BS and SB values:", &iw[1], leniw, (ftnlen)18);
	    i__1 = nbs;
	    for (k = 1; k <= i__1; ++k) {
		s_wsfi(&io___50);
		do_fio(&c__1, (char *)&kbs[k], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&xbs[k], (ftnlen)sizeof(doublereal));
		e_wsfi();
		snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)72);
	    }
	}
    } else {
/*        -------------------------------------------------------------- */
/*        Error exits (all fatal) */
/*        Some values of iExit invoke additional output. */
/*        -------------------------------------------------------------- */
	if (*iexit == 42) {
/*           --------------------------------------------- */
/*           The basis is singular after mtry tries. */
/*           Time to give up. */
/*           --------------------------------------------- */
	    s_wsfi(&io___51);
	    do_fio(&c__1, (char *)&c__40, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__14, str, &iw[1], leniw, (ftnlen)72);
	} else if (*iexit == 80) {
/* --------------------------------------------- */
/* Insufficient storage to factorize B. */
/* --------------------------------------------- */
	    more = iw[163] - lenalu;
	    newi = maxiw + (more << 1);
	    newr = maxrw + more;
	    iw[48] = newi;
/* minimum length of iw */
	    iw[49] = newr;
/* minimum length of rw */
	    if (maxlui < maxlur) {
		*iexit = 83;
/* not enough integer storage for LU */
	    } else {
		*iexit = 84;
/* not enough real    storage for LU */
	    }
	    s_wsfi(&io___55);
	    e_wsfi();
	    snprnt_(&c__13, str, &iw[1], leniw, (ftnlen)72);
	    s_wsfi(&io___56);
	    do_fio(&c__1, (char *)&maxiw, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&newi, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)72);
	    s_wsfi(&io___57);
	    do_fio(&c__1, (char *)&maxrw, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&newr, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)72);
	} else if (*iexit == 141) {
/*           --------------------------------------------- */
/*           Wrong number of basics. */
/*           --------------------------------------------- */
	    s_wsfi(&io___58);
	    do_fio(&c__1, (char *)&nbasic, (ftnlen)sizeof(integer));
	    e_wsfi();
	    snprnt_(&c__14, str, &iw[1], leniw, (ftnlen)72);
	}
    }
    return 0;
} /* s2bfac_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bfac */
/* Subroutine */ int s2bkbs_(integer *start, integer *iobj, integer *m, 
	integer *mbs, integer *n, integer *nb, integer *nnl, integer *hs, 
	integer *kbs, integer *kobj, integer *nbasic, integer *nonlin, 
	integer *nslack)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, jobj;

/*     ================================================================= */
/*     s2Bkbs  loads the basic variables into  kBS. */

/*     11 Apr 2004: First version of   s2Bkbs. */
/*     01 Jul 2005: Added option to load basics from the front. */
/*     03 Jul 2005: Current version of s2Bkbs. */
/*     ================================================================= */
/*     ----------------------------------------------------------------- */
/*     ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --kbs;
    --hs;

    /* Function Body */
    k = 0;
/* Counts the number of nonbasics */
    if (*start == 0) {
/* Load basics,  x first. */
	i__1 = *nb;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] == 3) {
		++k;
		kbs[k] = j;
	    }
	}
    } else if (*start == 1) {
/* Load basics,  slacks first. */
/* kObj points to the linear objective */
	jobj = *n + *iobj;
	*nonlin = 0;
	*kobj = 0;
	i__1 = *nb;
	for (j = *n + 1; j <= i__1; ++j) {
	    if (hs[j] == 3) {
		++k;
		kbs[k] = j;
		if (j == jobj) {
		    *kobj = k;
		}
	    }
	}
	*nslack = k;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (hs[j] == 3) {
		++k;
		if (k <= *m) {
		    kbs[k] = j;
		    if (j <= *nnl) {
			++(*nonlin);
		    }
		}
	    }
	}
    }
    *nbasic = k;
    return 0;
} /* s2bkbs_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bkbs */
/* Subroutine */ int s2bmap_(integer *m, integer *n, integer *ne, integer *
	maxs, integer *nextiw, integer *nextrw, integer *maxiw, integer *
	maxrw, integer *liwest, integer *lrwest, integer *iw, integer *leniw)
{
    static integer ip, iq, mbs, lua, mlu, nlu, indc, locc, lenc, mina, indr, 
	    locr, cols, markc, iploc, iqloc, lenri, markr, necola, lenalu, 
	    maxlui, lastiw, maxlur, lastrw;

/*     ================================================================== */
/*     s2Bmap sets up the core allocation for the basis factors. */

/*     Normally the storage is for B = LU. */
/*     For nonlinear problems, we may also need to factorize (B S)' = LU, */
/*     where S has at most maxS columns. */

/*     On entry: */
/*     nextiw, nextrw  say where the LU arrays can start in iw(*), rw(*). */
/*     maxiw, maxrw  say where the LU arrays must end  in iw(*), rw(*). */

/*     On exit: */
/*     liwEst, lrwEst  estimate the minimum length of iw(*), rw(*). */
/*     The LU routines will have some room to maneuver if the arrays */
/*     are at least that long.  (Later LU messages may ask for more.) */
/*     ------------------------------------------------------------------ */

/*     15 Nov 1991: First version based on Minos 5.4 routine m2bmap. */
/*     08 Nov 1993: Generalized to allow room for (B S)'. */
/*     11 Nov 1994: rw(*) replaced by iw(*) and rw(*). */
/*     14 Jul 1997: Thread-safe version. */
/*     17 Mar 1998: nextiw now points to start of integer LU workspace. */
/*     07 Mar 2002: Fixed storage allocation (minA) to allow for BS factors. */
/*     11 Sep 2006: Anders Goran reports overflow on QP problem with */
/*                  m = 20736, n = 28, ne = 242752. */
/*                  Changed   minA   = 6 * mBS * necolA */
/*                  to        minA   = 6 * min(mBS,n ) * necolA */
/*                  and ensured maxLUi, maxLUr are nonnegative. */
/*     18 May 2007: Fixed bug that saved different values for maxLUi */
/*                  and  maxLUr. */
/*     04 Jul 2010: maxLUi and maxLUr saved in workspace. */
/*     07 Jul 2013: cols, markc, markr needed for LUSOL's lu1mxr. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/*     Allocate arrays for an  mLU x nLU  matrix. */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    mbs = *m + *maxs;
    mlu = mbs;
    nlu = *m;
/*     LU integer workspace is  iw(iP:maxiw). */
/*     nextiw points to the start of indc(*), indr(*). */
/*     indc and indr are made as long as possible. */
    ip = *nextiw;
    iq = ip + mlu;
    lenc = iq + nlu;
    lenri = lenc + nlu;
    locc = lenri + mlu;
    locr = locc + nlu;
    iploc = locr + mlu;
    iqloc = iploc + nlu;
    cols = iqloc + nlu;
    markc = cols + nlu;
    markr = markc + nlu;
    lastiw = markr + mlu;
    *nextiw = lastiw;
    maxlui = (*maxiw - lastiw) / 2;
    maxlui = max(maxlui,0);
    indc = lastiw;
    indr = indc + maxlui;
/*     LU real workspace is  rw(LUa:maxrw) */
/*     nextrw points to the start of A(*). */
    lua = *nextrw;
    lastrw = *nextrw;
    maxlur = *maxrw - lastrw;
    maxlur = max(maxlur,0);
/*     LUSOL thinks indc(*), indr(*) and A(*) are all of length lenaLU. */
    lenalu = min(maxlui,maxlur);
/*     Estimate the number of nonzeros in the basis factorization. */
/*     necolA = estimate of nonzeros per column of  A. */
/*     We guess that the density of the basis factorization is */
/*     5 times as great, and then allow 1 more such lot for elbow room. */
    necola = *ne / *n;
    necola = max(necola,5);
/* !!   minA   = 6 * min( m, n ) * necolA  ! Too little for BSfac */
    mina = min(mbs,*n) * 6 * necola;
/* mBS is biggest LU dimension */
    mina += 10000;
/* So tiny problems have plenty */
    *liwest = lastiw + (mina << 1);
    *lrwest = lastrw + mina;
    iw[213] = lenalu;
/* space allotted for LU factors */
    iw[361] = maxlui;
/* max LU nonzeros in iw(*) */
    iw[362] = maxlur;
/* max LU nonzeros in rw(*) */
    iw[363] = ip;

    iw[364] = iq;

    iw[365] = lenc;

    iw[366] = lenri;

    iw[367] = locc;

    iw[368] = locr;

    iw[369] = iploc;

    iw[370] = iqloc;

    iw[371] = lua;

    iw[373] = indc;

    iw[374] = indr;

    iw[375] = cols;
/* 3 new work arrays for lu1mxr */
    iw[376] = markc;

    iw[377] = markr;

    return 0;
} /* s2bmap_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bmap */
/* Subroutine */ int s2newb_(integer *m, integer *mbs, integer *nb, integer *
	ns, integer *hs, integer *ip, integer *kbs, integer *kbsold, integer *
	locr, integer *nswap)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, m1, nbs;
    extern /* Subroutine */ int icopy_(integer *, integer *, integer *, 
	    integer *, integer *);

/*     ================================================================== */
/*     s2newB  permutes kBS(*) to reflect the permutation (B S)P, */
/*     where P is in iP(*).  It updates hs(*) accordingly. */
/*     kBSold(*) and locr(*) are needed for workspace. */

/*     30 Oct 1993: First version. */
/*     04 Nov 1993: kBSold, nSwap used to save old R if there's no */
/*                  change in the set of superbasics. */
/*     16 Sep 2000: Superbasic slacks counted. */
/*     02 Aug 2003: Superbasic slacks allowed. */
/*     02 Aug 2003: Current version of s2newB. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --kbsold;
    --kbs;
    --ip;
    --hs;

    /* Function Body */
    *nswap = 0;
    m1 = *m + 1;
    nbs = *m + *ns;
    icopy_(&nbs, &kbs[1], &c__1, &locr[1], &c__1);
    icopy_(ns, &kbs[m1], &c__1, &kbsold[m1], &c__1);
    i__1 = nbs;
    for (k = 1; k <= i__1; ++k) {
	i__ = ip[k];
	j = locr[i__];
	kbs[k] = j;
	if (k <= *m) {
	    hs[j] = 3;
	} else {
	    if (hs[j] != 2) {
		++(*nswap);
	    }
	    hs[j] = 2;
	}
    }
/*     Restore the old S ordering if S contains the same variables. */
    if (*nswap == 0) {
	icopy_(ns, &kbsold[m1], &c__1, &kbs[m1], &c__1);
    }
    return 0;
} /* s2newb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2newB */
/* Subroutine */ int s2blu_(integer *iexit, integer *task, integer *lprint, 
	integer *m, integer *n, integer *nb, integer *nbs, integer *ne, 
	integer *nloca, integer *loca, integer *inda, doublereal *acol, 
	integer *kbs, integer *ip, doublereal *alu, integer *indc, integer *
	indr, integer *lenalu, integer *iy, integer *iy1, doublereal *y, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, iq, ir, nz, locc, lenc, locr, cols;
    static doublereal oldl1;
    extern /* Subroutine */ int iload_(integer *, integer *, integer *, 
	    integer *);
    static integer markc, iploc, iqloc, lenri, markr, oldtp;
    extern /* Subroutine */ int lu1fac_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    static doublereal brtol1, bstol1;
    static integer lprdbg;
    static doublereal growth;

/*     ================================================================== */
/*     s2BLU  factorizes the basis. */

/*     Task = B   Extract basis from the constraint matrix */
/*                and factorize B = L U with current tols. */

/*     Task = BR  Factorize B = L U to check its rank. */
/*                Use TCP with tight tols and don't save L and U. */

/*     Task = BS  Factorize transpose of (B S), so that  (B') = L U, */
/*                                                       (S') */
/*                to get a new partition of (B S). */
/*                Use TCP with tight tols and don't save L and U. */

/*     The following tolerances are used... */

/*     luparm(3) = maxcol   lu1fac: Maximum number of columns */
/*                          searched allowed in a Markowitz-type */
/*                          search for the next pivot element. */
/*     luparm(6) = TPivot   TPivot = 0 means threshold partial  pivoting. */
/*                                   0 means threshold complete pivoting. */
/*     luparm(8) = keepLU   keepLU = 1 means keep L and U, */
/*                                   0 means discard them. */
/*     parmlu(1) = Lmax1  = Maximum multiplier allowed in  L  during */
/*                          refactorization. */
/*     parmlu(2) = Lmax2  = Maximum multiplier allowed during updates. */
/*     parmlu(3) = small  = Minimum element kept in  B  or in */
/*                          transformed matrix during elimination. */
/*     parmlu(4) = Utol1  = Abs tol for flagging small diagonals of  U. */
/*     parmlu(5) = Utol2  = Rel tol for flagging small diagonals of  U. */
/*     parmlu(6) = Uspace = Factor allowing waste space in row/col lists. */
/*     parmlu(7) = dens1    The density at which the Markowitz */
/*                          strategy should search maxcol columns */
/*                          and no rows. */
/*     parmlu(8) = dens2    The density at which the Markowitz */
/*                          strategy should search only 1 column. */
/*                          (In one version of lu1fac, the remaining */
/*                          matrix is treated as dense if there is */
/*                          sufficient storage.) */

/*     On exit, */
/*     iExit = 2  if there was excessive growth in U. */
/*                Other iExit values may be set by lu1fac: */
/*     iExit = 0  if the LU factors were computed. */
/*           = 1  if there are singularities (nSing gt 0). */
/*           = 2  if there was large growth in U. */
/*           = 3  if the matrix B has an illegal row or column index. */
/*           = 4  if an entry of B has the same indices as an earlier entry. */
/*           = 7  if insufficient storage for the LU. */
/*                minlen is an estimate of the necessary value of  lenaLU. */
/*           = 8  if there is a fatal error in lu1fac. */

/*     20 Oct 1990  Initial version based on Minos routine m2bsol. */
/*     07 Nov 1993: Add option to factorize (B S)' */
/*     06 Mar 1994: Include all rows of (B S), even if B contains slacks. */
/*     22 Apr 1994: Test for excessive growth in U. */
/*     14 Jul 1997: Thread-safe version. */
/*     23 Sep 2000: LU statistics now printed by LUSOL. */
/*     10 Mar 2001: Task = BR implemented. */
/*     09 May 2001: For BRfac, BSfac, force Lmax1 <= BRtol1, BStol1 */
/*                  (with TCP). */
/*     20 Dec 2004: Reduce BStol1 from 3.99 to 2.5. */
/*     20 Dec 2004: Current version of s2BLU. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* max L-multiplier in factor */

/* 0(1) TPP(TRP) */

/*     ------------------------------------------------------------------ */
/* minimum recommended lenaLU */
    /* Parameter adjustments */
    --y;
    --iy1;
    --iy;
    --ip;
    --kbs;
    --acol;
    --inda;
    --loca;
    --indr;
    --indc;
    --alu;
    --iw;
    --rw;

    /* Function Body */
    lprdbg = iw[85];
/* > 0    => private debug print */
    iq = iw[364];

    lenc = iw[365];

    lenri = iw[366];

    locc = iw[367];

    locr = iw[368];

    iploc = iw[369];

    iqloc = iw[370];

    cols = iw[375];
/* 3 new work arrays for lu1mxr */
    markc = iw[376];
/* passed to lu1fac and lu1fad */
    markr = iw[377];
/* from s2BLU */
    oldl1 = rw[151];
/* Save Lmax1 */
    oldtp = iw[156];
/* Save the TP state */
    brtol1 = 2.5;
/* !    BRtol1    = 3.99d+0 !!!!!! TEST FOR SIGEST PAPER (DRCAVTY2) */
/* Max Lmax1 for BRfac */
    bstol1 = 2.5;
/* Max Lmax1 for BSfac */
    *iexit = 0;
/*     -------------------------------------------------- */
/*     Set Print level for lu1fac. */
/*        iw(LUprnt) = 1  for errors, */
/*                   = 10 for statistics */
/*                   = 50 for debug info */
/*     -------------------------------------------------- */
    iw[152] = min(*lprint,10);
    if (lprdbg == 51) {
	iw[152] = 50;
    }
/*     -------------------------------------------------- */
/*     Set other parameters for lu1fac. */
/*     -------------------------------------------------- */
    if (*task == 0) {
	iw[158] = 1;
    } else if (*task == 1) {
/* Force TRP */
	if (oldtp == 0) {
	    iw[156] = 1;
/* ! iw(TPivot) = TPP !!!!! TEST FOR SIGEST PAPER (DRCAVTY2) */
	}
	iw[158] = 0;
	rw[151] = min(oldl1,brtol1);
/* Make sure Lmax1 is reasona */
    } else if (*task == 2) {
/* Stay with current TPP or TRP */
	iw[158] = 0;
	rw[151] = min(oldl1,bstol1);
/* Make sure Lmax1 is reasona */
    }
    if (*task == 0 || *task == 1) {
/*        ------------------------------------- */
/*        Estimate the number of nonzeros in B. */
/*        ------------------------------------- */
	nz = 0;
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j == 0) {
/*              -------------------------- */
/*              Relax, just a zero column. */
/*              -------------------------- */
	    } else if (j <= *n) {
/*              ---------------------------- */
/*              Basic column from A. */
/*              ---------------------------- */
		nz = nz + loca[j + 1] - loca[j];
	    } else {
/*              --------------------- */
/*              Basic slack. */
/*              --------------------- */
		++nz;
	    }
	}
	iw[163] = nz * 5 / 4;
	if (iw[163] > *lenalu) {
	    *iexit = 7;
	    goto L900;
	}
/*        --------------------------------------------------------------- */
/*        Load B into LUSOL. */
/*        --------------------------------------------------------------- */
	nz = 0;
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j == 0) {
/*              -------------------------- */
/*              Relax, just a zero column. */
/*              -------------------------- */
	    } else if (j <= *n) {
/*              --------------------- */
/*              Basic column from A. */
/*              --------------------- */
		i__2 = loca[j + 1] - 1;
		for (i__ = loca[j]; i__ <= i__2; ++i__) {
		    ir = inda[i__];
		    ++nz;
		    alu[nz] = acol[i__];
		    indc[nz] = ir;
		    indr[nz] = k;
		}
	    } else {
/*              --------------------- */
/*              Basic slacks. */
/*              --------------------- */
		++nz;
		alu[nz] = -1.;
		indc[nz] = j - *n;
		indr[nz] = k;
	    }
	}
/*        iy and iy1 are work vectors */
/*        y is an output parameter, used by s2sing. */
	lu1fac_(m, m, &nz, lenalu, &iw[151], &rw[151], &alu[1], &indc[1], &
		indr[1], &ip[1], &iw[iq], &iw[lenc], &iw[lenri], &iw[locc], &
		iw[locr], &iw[iploc], &iw[iqloc], &iy[1], &iy1[1], &iw[cols], 
		&iw[markc], &iw[markr], &y[1], iexit);
/*        Test for excessive growth in U. */
	growth = rw[166];
/* TPP: Umax/Amax    TRP: Akmax/Amax */
	if (*iexit == 0 && growth >= 1e8) {
	    *iexit = 2;
	}
    } else if (*task == 2) {
/*        --------------------------------------------------------------- */
/*        Factorize (B S)' = LU without keeping L and U. */
/*        --------------------------------------------------------------- */
/*        Extract (B S)'. */
/*        iP(1:m) is needed for workspace. */
/*        iP(i) = 0 except for rows with a basic or superbasic slack. */
/*        We can ignore all of these rows except for the slack itself. */
/*        06 Mar 1994: Keep all rows.  (Made a difference in MINOS) */
/*        22 Apr 1994: Make sure the objective slack remains basic!! */
/*        11 Nov 1994: Go back to ignoring rows with a slack in B. */
/*                   This means we don't have to worry about the Obj. */
/*        --------------------------------------------------------------- */
	iload_(m, &c__0, &ip[1], &c__1);
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j > *n) {
		ip[j - *n] = 1;
	    }
	}
/*        Count the number of nonzeros in ( B S ). */
	nz = 0;
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *n) {
		i__2 = loca[j + 1] - 1;
		for (i__ = loca[j]; i__ <= i__2; ++i__) {
		    ir = inda[i__];
		    if (ip[ir] == 0) {
			++nz;
		    }
		}
	    } else {
		++nz;
	    }
	}
	iw[163] = nz * 5 / 4;
	if (iw[163] > *lenalu) {
	    *iexit = 7;
	    goto L900;
	}
	nz = 0;
	i__1 = *nbs;
	for (k = 1; k <= i__1; ++k) {
	    j = kbs[k];
	    if (j <= *n) {
		i__2 = loca[j + 1] - 1;
		for (i__ = loca[j]; i__ <= i__2; ++i__) {
		    ir = inda[i__];
		    if (ip[ir] == 0) {
			++nz;
			alu[nz] = acol[i__];
			indc[nz] = k;
			indr[nz] = ir;
		    }
		}
	    } else {
		++nz;
		alu[nz] = -1.;
		indc[nz] = k;
		indr[nz] = j - *n;
	    }
	}
	lu1fac_(nbs, m, &nz, lenalu, &iw[151], &rw[151], &alu[1], &indc[1], &
		indr[1], &ip[1], &iw[iq], &iw[lenc], &iw[lenri], &iw[locc], &
		iw[locr], &iw[iploc], &iw[iqloc], &iy[1], &iy1[1], &iw[cols], 
		&iw[markc], &iw[markr], &y[1], iexit);
    }
/*     -------------------------------------------------- */
/*     Restore Lmax1 etc. for BR and BS. */
/*     -------------------------------------------------- */
L900:
    rw[151] = oldl1;
    iw[156] = oldtp;
    return 0;
} /* s2blu_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2BLU */
/* Subroutine */ int s2bmod_(integer *iexit, integer *jrep, integer *m, 
	doublereal *z__, integer *iw, integer *leniw, doublereal *rw, integer 
	*lenrw)
{
    static integer ip, iq, lua;
    static doublereal diag;
    static integer indc, locc, lenc, indr, locr, lenri;
    static doublereal znorm;
    extern /* Subroutine */ int lu8rpc_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static integer lenalu;

/*     ================================================================== */
/*     s2Bmod  updates the LU factors of B when column "jrep" is replaced */
/*             by a vector  v.  On entry,   z  must satisfy  L z = v. */
/*             It is overwritten. */

/*     20 Oct 1990  Initial version. */
/*     14 Jul 1997: Thread-safe version. */
/*     01 Jun 1999: Current version of s2Bmod. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --z__;
    --iw;
    --rw;

    /* Function Body */
    lenalu = iw[213];
/* space allotted for LU factors */
    ip = iw[363];

    iq = iw[364];

    lenc = iw[365];

    lenri = iw[366];

    locc = iw[367];

    locr = iw[368];

    lua = iw[371];

    indc = iw[373];

    indr = iw[374];

    lu8rpc_(&c__1, &c__2, m, m, jrep, &z__[1], &z__[1], &lenalu, &iw[151], &
	    rw[151], &rw[lua], &iw[indc], &iw[indr], &iw[ip], &iw[iq], &iw[
	    lenc], &iw[lenri], &iw[locc], &iw[locr], iexit, &diag, &znorm);
    return 0;
} /* s2bmod_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bmod */
/* Subroutine */ int s2bmod2_(integer *iexit, integer *jrep, integer *m, 
	doublereal *z__, integer *iw, integer *leniw, doublereal *rw, integer 
	*lenrw)
{
    extern /* Subroutine */ int s2bmod_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *);
    static doublereal utol1s, utol2s;

/*     ================================================================== */
/*     s2Bmod2  updates the LU factors of B with tight singularity tol. */

/*     22 Jun 2004  Initial version. */
/*     22 Jun 2004: Current version of s2Bmod2. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* abs tol for small diag of U */
/* rel tol for small diag of U */
/* abs tol for small diag of U. */
/*     ------------------------------------------------------------------ */
/* rel tol for small diag of U. */
    /* Parameter adjustments */
    --z__;
    --iw;
    --rw;

    /* Function Body */
    utol1s = rw[154];
    utol2s = rw[155];
    rw[154] = rw[63];
    rw[155] = rw[64];
    s2bmod_(iexit, jrep, m, &z__[1], &iw[1], leniw, &rw[1], lenrw);
    rw[154] = utol1s;
    rw[155] = utol2s;
    return 0;
} /* s2bmod2_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bmod2 */
/* Subroutine */ int s2bsol_(integer *iexit, integer *task, integer *m, 
	doublereal *z__, doublereal *y, integer *iw, integer *leniw, 
	doublereal *rw, integer *lenrw)
{
    static integer ip, iq, lua, indc, locc, lenc, indr, locr, lenri;
    static doublereal small0;
    extern /* Subroutine */ int lu6sol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *);
    static integer lenalu;
    extern doublereal dnormi_(integer *, doublereal *, integer *);

/*     ================================================================== */
/*     s2Bsol  solves various systems with the LU factors of B. */
/*     Task  selects one of the following: */
/*      Task          Action */
/*      ----          ------ */
/*      with L    Solve  L*z = z(input).    y  is not touched. */
/*      with B    Solve  L*z = z(input)  and solve  B*y = z(input). */
/*      with Bt   Solve  B(transpose)*y = z.  Note that  z  is destroyed. */

/*     20 Oct 1990: Initial version. */
/*     16 Nov 2001: dnormi added. */
/*     06 Aug 2003: adopted snEXIT. */
/*     06 Aug 2003: Current version of s2Bsol */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
/* defn of small real */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --y;
    --z__;
    --iw;
    --rw;

    /* Function Body */
    lenalu = iw[213];
/* space allotted for LU factors */
    ip = iw[363];

    iq = iw[364];

    lenc = iw[365];

    lenri = iw[366];

    locc = iw[367];

    locr = iw[368];

    lua = iw[371];

    indc = iw[373];

    indr = iw[374];

    if (*task == 0 || *task == 1) {
/*        --------------------------------------------------------------- */
/*        Solve   L*z = z(input). */
/*        When LU*y = z is being solved in SNOPT, norm(z) will sometimes */
/*        be small (e.g. after periodic refactorization).  Hence for */
/*        solves with L we scale parmlu(3) to alter what lu6sol thinks */
/*        is small. */
/*        --------------------------------------------------------------- */
	small0 = rw[153];
	if (*task == 1) {
	    rw[153] = small0 * dnormi_(m, &z__[1], &c__1);
	}
	lu6sol_(&c__1, m, m, &z__[1], &y[1], &lenalu, &iw[151], &rw[151], &rw[
		lua], &iw[indc], &iw[indr], &iw[ip], &iw[iq], &iw[lenc], &iw[
		lenri], &iw[locc], &iw[locr], iexit);
	rw[153] = small0;
	if (*task == 1) {
/*           ------------------------------------------------------------ */
/*           Task = solve with B.   Solve  U*y = z. */
/*           ------------------------------------------------------------ */
	    lu6sol_(&c__3, m, m, &z__[1], &y[1], &lenalu, &iw[151], &rw[151], 
		    &rw[lua], &iw[indc], &iw[indr], &iw[ip], &iw[iq], &iw[
		    lenc], &iw[lenri], &iw[locc], &iw[locr], iexit);
	}
    } else if (*task == 2) {
/*        --------------------------------------------------------------- */
/*        Task = solve with B transpose.  Solve  B'*y = z. */
/*        --------------------------------------------------------------- */
	lu6sol_(&c__6, m, m, &y[1], &z__[1], &lenalu, &iw[151], &rw[151], &rw[
		lua], &iw[indc], &iw[indr], &iw[ip], &iw[iq], &iw[lenc], &iw[
		lenri], &iw[locc], &iw[locr], iexit);
    }
    if (*iexit != 0) {
	*iexit = 142;
    }
/* error in basis package */
    return 0;
} /* s2bsol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2Bsol */
/* Subroutine */ int s2sb_(integer *m, integer *mbs, integer *nb, integer *ns,
	 integer *nbs, integer *hs, integer *kbs, doublereal *xbs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, jq, ns0;

/*     ================================================================== */
/*     s2SB  finds if any superbasic slacks were made basic by s2sing. */
/*     If any are found, nS and nBS are updated and the corresponding */
/*     entries are removed from  kBS and xBS. */

/*     02 Jul 2005: First version of s2sb. */
/*     02 Jul 2005: Current version of s2sb. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --xbs;
    --kbs;
    --hs;

    /* Function Body */
    if (*ns == 0) {
	return 0;
    }
    ns0 = *ns;
    for (jq = ns0; jq >= 1; --jq) {
	j = kbs[*m + jq];
	if (hs[j] == 3) {
	    --(*ns);
	    *nbs = *m + *ns;
	    i__1 = *nbs;
	    for (k = *m + jq; k <= i__1; ++k) {
		kbs[k] = kbs[k + 1];
		xbs[k] = xbs[k + 1];
	    }
	}
    }
    return 0;
} /* s2sb_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2sb */
/* Subroutine */ int s2sing_(integer *lprint, integer *mbs, integer *m, 
	integer *n, integer *nb, doublereal *z__, integer *ip, integer *iq, 
	doublereal *bl, doublereal *bu, integer *hs, integer *kbs, doublereal 
	*x, integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Column\002,i7,\002  replaced by slack"
	    "\002,i7)";
    static char fmt_1100[] = "(\002 and so on.  Total slacks inserted =\002,"
	    "i6)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static integer i__, j, k;
    static char str[72];
    static integer nsing;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___146 = { 0, str, 0, fmt_1000, 72, 1 };
    static icilist io___147 = { 0, str, 0, fmt_1100, 72, 1 };


/*     ================================================================= */
/*     s2sing  is called if the LU factorization of the basis appears */
/*     to be singular.   If  z(j)  is not positive, the  jth  basic */
/*     variable  kBS(j)  is replaced by the appropriate slack. */
/*     If any kBS(j) = 0, only a partial basis was supplied. */

/*     30 Sep 1991: First version based on minos routine m2sing. */
/*     29 May 1995: Optional swapping of slack and basic. */
/*     12 Jul 1997: Thread-safe version. */
/*     18 Nov 2001: lPrint added as parameter. */
/*     31 Jul 2003: snPRNT adopted. */
/*     ================================================================= */
/*     ----------------------------------------------------------------- */
/*     ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --kbs;
    --iq;
    --ip;
    --z__;
    --x;
    --hs;
    --bu;
    --bl;
    --iw;

    /* Function Body */
    nsing = 0;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	j = iq[k];
	if (z__[j] <= 0.) {
	    j = kbs[j];
	    if (j > 0) {
/*              Make variable  j  nonbasic (and feasible). */
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
	    }
/*           Make the appropriate slack basic. */
	    i__ = ip[k];
	    hs[*n + i__] = 3;
	    ++nsing;
	    if (*lprint >= 10 && nsing <= 5) {
		s_wsfi(&io___146);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfi();
		snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)72);
	    }
	}
    }
    if (*lprint >= 10 && nsing > 5) {
	s_wsfi(&io___147);
	do_fio(&c__1, (char *)&nsing, (ftnlen)sizeof(integer));
	e_wsfi();
	snprnt_(&c__3, str, &iw[1], leniw, (ftnlen)72);
    }
    return 0;
} /* s2sing_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2sing */
/* Subroutine */ int s2tols_(integer *mode, logical *newtol, integer *itn, 
	integer *iw, integer *leniw, doublereal *rw, integer *lenrw)
{
    /* Initialized data */

    static char tp[8*3] = "partial " "rook    " "complete";

    /* Format strings */
    static char fmt_1000[] = "(\002 Itn\002,i7,\002: LU \002,a,\002 pivoting"
	    " tols \002,2f10.2)";

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[72];
    static doublereal oldl1, oldl2, lmax1, lmax2, tolfac, toldcp, toldpp, 
	    toldrp, toldup, tolupd;
    static integer lvlpiv, tpivot;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___162 = { 0, str, 0, fmt_1000, 72, 1 };


/*     ================================================================= */
/*     s2tols  sets the LU Factor and Update tolerances. */

/*     mode   (input) says what should be done: */
/*     mode = DefTol  Set default tols for current TP (TPP, TRP or TCP). */
/*     mode = RstTol  Set default tols for TP specified by user. */
/*     mode = RedTol  Reduce tols or switch from TPP to TRP. */
/*     mode = MinTol  Set minimum tols for current TP. */

/*     NewTol (output) says if the LU tols were changed */
/*                     or we switched from TPP to TRP. */

/*     30 Aug 2000: First version of s2tols. */
/*     10 Dec 2002: Rook and diagonal pivoting options added. */
/*     31 Jul 2003: snPRNT adopted. */
/*     21 Dec 2004: Current version of s2tols. */
/*     ================================================================= */
/*     ------------------------------------------------------------------ */
/* Minimum LU tols allowed */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
/*     ------------------------------------------------------------------ */
    tolfac = rw[66];
/* LU factor tolerance = user-defined Lmax1 */
    tolupd = rw[67];
/* LU update tolerance = user-defined Lmax2 */
    lvlpiv = iw[80];
/* 0(1 2 3) LU Partial (Rook Complete Diag) piv */
    tpivot = iw[156];
/* Current LU pivot option */
    lmax1 = rw[151];
/* max allowable L-multiplier in factor */
    lmax2 = rw[152];
/* max allowable L-multiplier in update */
    toldpp = rw[181];
/* default Lmax1 for TPP */
    toldrp = rw[187];
/* default Lmax1 for TRP */
    toldcp = rw[182];
/*    tolDdp    = rw(188) ! default Lmax1 for TDP */
/* default Lmax1 for TCP */
    toldup = rw[183];
/* default Lmax2 */
    oldl1 = lmax1;
    oldl2 = lmax2;
    *newtol = FALSE_;
    if (*mode == 0) {
/*        --------------------------------------------------------------- */
/*        Set the default LU Factor tol and LU Update tol. */
/*        --------------------------------------------------------------- */
	if (tpivot == 0) {
	    lmax1 = toldpp;
	} else if (tpivot == 1) {
	    lmax1 = toldrp;
	} else if (tpivot == 2) {
	    lmax1 = toldcp;
	}
	lmax2 = toldup;
    } else if (*mode == 3) {
/*        --------------------------------------------------------------- */
/*        Reset the user LU Factor tol and LU Update tol. */
/*        --------------------------------------------------------------- */
	tpivot = lvlpiv;
	lmax1 = tolfac;
	lmax2 = tolupd;
    } else if (*mode == 1) {
/*        --------------------------------------------------------------- */
/*        Reduce LU Factor tol and LU Update tol, */
/*        perhaps switching to TRP first. */
/*        13 Aug 2007: Note: TCP is not activated here -- only TRP. */
/*        TCP must be requested via the LU Complete Pivoting option. */
/*        It is then used for all factorizations. */
/*        --------------------------------------------------------------- */
	if (tpivot == 0 && lmax1 <= 1.011) {
/* Switch */
	    tpivot = 1;
	    lmax1 = toldrp;
	    lmax2 = toldup;
/* else if (TPivot .eq. TRP  .and.  Lmax1 .le. tolTRP) then */
/*  TPivot = TCP */
/*  Lmax1  = tolDcp   !!! Comment out for SIGEST test (DRCAVTY2) */
/*  Lmax2  = tolDup   !!! Comment out for SIGEST test (DRCAVTY2) */
	} else {
	    if (lmax1 > 4.) {
		lmax1 *= .5;
	    } else {
		lmax1 = sqrt(lmax1);
	    }
	    if (lmax2 > 4.) {
		lmax2 *= .5;
	    } else {
		lmax2 = sqrt(lmax2);
	    }
	}
    } else if (*mode == 2) {
/*        --------------------------------------------------------------- */
/*        Set LU Factor tol and LU Update tol to their minimum values. */
/*        --------------------------------------------------------------- */
	if (tpivot == 0) {
	    lmax1 = 1.011;
	} else if (tpivot == 1) {
	    lmax1 = 1.011;
	} else if (tpivot == 2) {
	    lmax1 = 1.011;
	}
	lmax2 = 1.011;
    }
/*     ------------------------------------------------------------------ */
/*     Make sure the tols aren't too small or big. */
/*     ------------------------------------------------------------------ */
    if (tpivot == 0) {
	lmax1 = max(lmax1,1.011);
    } else if (tpivot == 1) {
	lmax1 = max(lmax1,1.011);
    } else if (tpivot == 2) {
	lmax1 = max(lmax1,1.011);
    }
    lmax2 = max(lmax2,1.011);
    lmax2 = min(lmax1,lmax2);
/*     ------------------------------------------------------------------ */
/*     Print significant changes. */
/*     ------------------------------------------------------------------ */
/* Update tol never more than Factol */
    *newtol = lmax1 != oldl1 || lmax2 != oldl2;
/* !    TPivot = 0 !!!!!! TEST FOR SIGEST PAPER (DRCAVTY2) */
    if (*newtol) {
	s_wsfi(&io___162);
	do_fio(&c__1, (char *)&(*itn), (ftnlen)sizeof(integer));
	do_fio(&c__1, tp + (tpivot << 3), (ftnlen)8);
	do_fio(&c__1, (char *)&lmax1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&lmax2, (ftnlen)sizeof(doublereal));
	e_wsfi();
	snprnt_(&c__23, str, &iw[1], leniw, (ftnlen)72);
    }
    iw[156] = tpivot;
    rw[151] = lmax1;
    rw[152] = lmax2;
    return 0;
} /* s2tols_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s2tols */
/* Subroutine */ int s2trylu_(integer *itn, integer *lureq0, integer *ns, 
	integer *lureq, logical *luok, integer *typelu, integer *iw, integer *
	leniw, doublereal *rw, integer *lenrw)
{
    static integer nbfac, lumod;
    extern /* Subroutine */ int s2tols_(integer *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    static logical newtol;

/*     ================================================================= */
/*     s2tryLU  is called when numerical difficulties imply that a new */
/*     factorize is needed. */

/*     On exit: */
/*        LUok    says if a new LU is possible. */
/*        typeLU  defines the type of factorization to be done. */
/*        LUreq   is used for printing and indicates the reason for */
/*                the LU request. */

/*     26 Dec 2003: First version of s2tryLU. */
/*     09 Apr 2004: Current version of s2tryLU. */
/*     ================================================================= */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --rw;

    /* Function Body */
    nbfac = iw[211];
/* number of consecutive `B' factorizes */
    lumod = iw[216];
/* number of LU mods since the last factorize */
    *lureq = *lureq0;
    *luok = TRUE_;
    if (nbfac > 1 && *ns > 1) {
/*        The LU has been computed since the last BS factorize. */
/*        Try to find a better basis. */
	*typelu = 2;
    } else if (lumod > 0 && (lumod != 1 || *lureq != 5)) {
/*        The LU has been modified since the last factorize. */
/*        Refactorize with the current LU tolerances. */
/*        A BS factorize will be used if needed. */
	*typelu = 3;
    } else {
/*        We have the LU at the current point. */
/*        Refactorize with smaller LU tolerances. */
	s2tols_(&c__1, &newtol, itn, &iw[1], leniw, &rw[1], lenrw);
	if (newtol) {
	    *luok = TRUE_;
	    if (*ns == 0) {
		*typelu = 0;
	    } else {
		*typelu = 2;
	    }
	} else {
	    *luok = FALSE_;
	    *lureq = 0;
	}
    }
    return 0;
} /* s2trylu_ */

