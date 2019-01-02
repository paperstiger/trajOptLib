/* ./src/sn27lu77.f -- translated by f2c (version 20100827).
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
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__2 = 2;

/* *********************************************************************** */
/*     File mi27lu.f and sn27lu.f. */

/*     This file contains most of the sparse LU package LUSOL */
/*     (the parts needed by MINOS, SQOPT and SNOPT). */
/*     The parts included are */

/*        27HEAD.f    These comments. */
/*        lusol1.f    Factor a given matrix A from scratch (lu1fac). */
/*        lusol2.f    Heap-management routines for lu1fac. */
/*        lusol6a.f   Solve with the current LU factors. */
/*        lusol7a.f   Utilities for all update routines. */
/*        lusol8a.f   Replace a column (Bartels-Golub update). */

/* *********************************************************************** */

/*     File lusol1.f */

/*     lu1fac   lu1fad   lu1gau   lu1mar   lu1mRP   lu1mCP   lu1mSP */
/*     lu1pen   lu1mxc   lu1mxr   lu1or1   lu1or2   lu1or3   lu1or4 */
/*     lu1pq1   lu1pq2   lu1pq3   lu1rec   lu1slk */
/*     lu1ful   lu1DPP   lu1DCP */

/* 26 Apr 2002: TCP implemented using heap data structure. */
/* 01 May 2002: lu1DCP implemented. */
/* 07 May 2002: lu1mxc must put 0.0 at top of empty columns. */
/* 09 May 2002: lu1mCP implements Markowitz with cols searched */
/*              in heap order. */
/*              Often faster (searching 20 or 40 cols) but more dense. */
/* 11 Jun 2002: TRP implemented. */
/*              lu1mRP implements Markowitz with Threshold Rook Pivoting. */
/*              lu1mxc maintains max col elements.  (Previously lu1max.) */
/*              lu1mxr maintains max row elements. */
/* 12 Jun 2002: lu1mCP seems too slow on big problems (e.g. memplus). */
/*              Disabled it for the moment.  (Use lu1mar + TCP.) */
/* 14 Dec 2002: TSP implemented. */
/*              lu1mSP implements Markowitz with */
/*              Threshold Symmetric Pivoting. */
/* 07 Mar 2003: character*1, character*2 changed to f90 form. */
/*              Comments changed from * in column to ! in column 1. */
/*              Comments kept within column 72 to avoid compiler warning. */
/* 19 Dec 2004: Hdelete(...) has new input argument Hlenin. */
/* 21 Dec 2004: Print Ltol and Lmax with e10.2 instead of e10.1. */
/* 26 Mar 2006: lu1fad: Ignore nsing from lu1ful. */
/*              lu1DPP: nsing redefined (but not used by lu1fad). */
/*              lu1DCP: nsing redefined (but not used by lu1fad). */
/* 30 Mar 2011: In lu1pen, loop 620 revised following advice from Ralf Östermark, */
/*              School of Business and Economics at Åbo Akademi University. */
/*              The previous version produced a core dump on the Cray XT. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int lu1fac_(integer *m, integer *n, integer *nelem, integer *
	lena, integer *luparm, doublereal *parmlu, doublereal *a, integer *
	indc, integer *indr, integer *ip, integer *iq, integer *lenc, integer 
	*lenr, integer *locc, integer *locr, integer *iploc, integer *iqloc, 
	integer *ipinv, integer *iqinv, doublereal *w, integer *inform__)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 m\002,i12,\002 \002,a,\002n\002,i12,\002"
	    "  Elems\002,i9,\002  Amax\002,1p,e10.1,\002  Density\002,0p,f7.2)"
	    ;
    static char fmt_1300[] = "(/\002 lu1fac  error...  entry  a(\002,i8,\002"
	    ")  has an illegal\002,\002 row or column index\002//\002 indc, i"
	    "ndr =\002,2i8)";
    static char fmt_1400[] = "(/\002 lu1fac  error...  entry  a(\002,i8,\002"
	    ")  has the same\002,\002 indices as an earlier entry\002//\002 i"
	    "ndc, indr =\002,2i8)";
    static char fmt_1700[] = "(/\002 lu1fac  error...  insufficient storag"
	    "e\002//\002 Increase  lena  from\002,i10,\002  to at least\002,i"
	    "10)";
    static char fmt_1800[] = "(/\002 lu1fac  error...  fatal bug\002,\002   "
	    "(sorry --- this should never happen)\002)";
    static char fmt_1900[] = "(/\002 lu1fac  error...  TSP used but\002,\002"
	    " diagonal pivot could not be found\002)";
    static char fmt_1100[] = "(\002 Merit\002,0p,f8.1,\002  lenL\002,i9,\002"
	    "  L+U\002,i11,\002  Cmpressns\002,i5,\002  Incres\002,0p,f8.2"
	    "/\002 Utri\002,i9,\002  lenU\002,i9,\002  Ltol\002,1p,e10.2,\002"
	    "  Umax\002,e10.1,\002  Ugrwth\002,e8.1/\002 Ltri\002,i9,\002  de"
	    "nse1\002,i7,\002  Lmax\002,e10.2)";
    static char fmt_1120[] = "(\002 Mer\002,a2,0p,f8.1,\002  lenL\002,i9,"
	    "\002  L+U\002,i11,\002  Cmpressns\002,i5,\002  Incres\002,0p,f8."
	    "2/\002 Utri\002,i9,\002  lenU\002,i9,\002  Ltol\002,1p,e10.2,"
	    "\002  Umax\002,e10.1,\002  Ugrwth\002,e8.1/\002 Ltri\002,i9,\002"
	    "  dense1\002,i7,\002  Lmax\002,e10.2,\002  Akmax\002,e9.1,\002  "
	    "Agrwth\002,e8.1)";
    static char fmt_1200[] = "(\002 bump\002,i9,\002  dense2\002,i7,\002  DU"
	    "max\002,1p,e9.1,\002  DUmin\002,e9.1,\002  condU\002,e9.1)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, l, l2;
    static doublereal dm, dn;
    static integer ll, lm, lu, ncp;
    static logical tcp, tpp, trp, tsp;
    static integer loch;
    static doublereal amax;
    static integer lenh, lenl;
    static doublereal lmax;
    static integer lenu, lerr;
    static char kpiv[2*4];
    static integer lpiv;
    static doublereal ltol, umax;
    static integer lrow, nout, lena2, numl0;
    static doublereal delem, dincr, akmax;
    static integer lenlk;
    static doublereal small;
    static integer nrank, jsing;
    static doublereal condu, dumin;
    static integer nsing;
    static doublereal dumax;
    static integer lenuk, nbump, jumin;
    static char mnkey[1];
    extern /* Subroutine */ int lu1or1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *), lu1or2_(integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *);
    static integer nltri;
    extern /* Subroutine */ int lu1or3_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     lu1or4_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), lu1pq1_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    static integer lmaxr, ltopl, nmove, nutri, numnz;
    extern /* Subroutine */ int lu1fad_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer nelem0, ndens1, ndens2;
    extern /* Subroutine */ int lu6chk_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static doublereal avgmer;
    static integer minlen;
    static logical keeplu;
    static integer llsave;
    static doublereal agrwth;
    static integer idummy;
    static doublereal densty;
    static integer lprint, mersum;
    static doublereal growth, ugrwth;

    /* Fortran I/O blocks */
    static cilist io___40 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_1300, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_1400, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_1700, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_1800, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_1900, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_1120, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_1200, 0 };


/*     ------------------------------------------------------------------ */
/*     lu1fac computes a factorization A = L*U, where A is a sparse */
/*     matrix with m rows and n columns, P*L*P' is lower triangular */
/*     and P*U*Q is upper triangular for certain permutations P, Q */
/*     (which are returned in the arrays ip, iq). */
/*     Stability is ensured by limiting the size of the elements of L. */

/*     The nonzeros of A are input via the parallel arrays a, indc, indr, */
/*     which should contain nelem entries of the form    aij,    i,    j */
/*     in any order.  There should be no duplicate pairs         i,    j. */

/*     ****************************************************************** */
/*     *        Beware !!!   The row indices i must be in indc,         * */
/*     *              and the column indices j must be in indr.         * */
/*     *              (Not the other way round!)                        * */
/*     ****************************************************************** */

/*     It does not matter if some of the entries in a(*) are zero. */
/*     Entries satisfying  abs( a(i) ) .le. parmlu(3)  are ignored. */
/*     Other parameters in luparm and parmlu are described below. */

/*     The matrix A may be singular.  On exit, nsing = luparm(11) gives */
/*     the number of apparent singularities.  This is the number of */
/*     "small" diagonals of the permuted factor U, as judged by */
/*     the input tolerances Utol1 = parmlu(4) and  Utol2 = parmlu(5). */
/*     The diagonal element diagj associated with column j of A is */
/*     "small" if */
/*                 abs( diagj ) .le. Utol1 */
/*     or */
/*                 abs( diagj ) .le. Utol2 * max( uj ), */

/*     where max( uj ) is the maximum element in the j-th column of U. */
/*     The position of such elements is returned in w(*).  In general, */
/*     w(j) = + max( uj ),  but if column j is a singularity, */
/*     w(j) = - max( uj ).  Thus, w(j) .le. 0 if column j appears to be */
/*     dependent on the other columns of A. */

/*     NOTE: lu1fac (like certain other sparse LU packages) does not */
/*     treat dense columns efficiently.  This means it will be slow */
/*     on "arrow matrices" of the form */
/*                  A = (x       a) */
/*                      (  x     b) */
/*                      (    x   c) */
/*                      (      x d) */
/*                      (x x x x e) */
/*     if the numerical values in the dense column allow it to be */
/*     chosen LATE in the pivot order. */

/*     With TPP (Threshold Partial Pivoting), the dense column is */
/*     likely to be chosen late. */

/*     With TCP (Threshold Complete Pivoting), if any of a,b,c,d */
/*     is significantly larger than other elements of A, it will */
/*     be chosen as the first pivot and the dense column will be */
/*     eliminated, giving reasonably sparse factors. */
/*     However, if element e is so big that TCP chooses it, the factors */
/*     will become dense.  (It's hard to win on these examples!) */
/*     ================================================================== */


/*     Notes on the array names */
/*     ------------------------ */

/*     During the LU factorization, the sparsity pattern of the matrix */
/*     being factored is stored twice: in a column list and a row list. */

/*     The column list is ( a, indc, locc, lenc ) */
/*     where */
/*           a(*)    holds the nonzeros, */
/*           indc(*) holds the indices for the column list, */
/*           locc(j) points to the start of column j in a(*) and indc(*), */
/*           lenc(j) is the number of nonzeros in column j. */

/*     The row list is    (    indr, locr, lenr ) */
/*     where */
/*           indr(*) holds the indices for the row list, */
/*           locr(i) points to the start of row i in indr(*), */
/*           lenr(i) is the number of nonzeros in row i. */


/*     At all stages of the LU factorization, ip contains a complete */
/*     row permutation.  At the start of stage k,  ip(1), ..., ip(k-1) */
/*     are the first k-1 rows of the final row permutation P. */
/*     The remaining rows are stored in an ordered list */
/*                          ( ip, iploc, ipinv ) */
/*     where */
/*           iploc(nz) points to the start in ip(*) of the set of rows */
/*                     that currently contain nz nonzeros, */
/*           ipinv(i)  points to the position of row i in ip(*). */

/*     For example, */
/*           iploc(1) = k   (and this is where rows of length 1 begin), */
/*           iploc(2) = k+p  if there are p rows of length 1 */
/*                          (and this is where rows of length 2 begin). */

/*     Similarly for iq, iqloc, iqinv. */
/*     ================================================================== */


/*     00 Jun 1983  Original version. */
/*     00 Jul 1987  nrank  saved in luparm(16). */
/*     12 Apr 1989  ipinv, iqinv added as workspace. */
/*     26 Apr 1989  maxtie replaced by maxcol in Markowitz search. */
/*     16 Mar 1992  jumin  saved in luparm(19). */
/*     10 Jun 1992  lu1fad has to move empty rows and cols to the bottom */
/*                  (via lu1pq3) before doing the dense LU. */
/*     12 Jun 1992  Deleted dense LU (lu1ful, lu1vlu). */
/*     25 Oct 1993  keepLU implemented. */
/*     07 Feb 1994  Added new dense LU (lu1ful, lu1den). */
/*     21 Dec 1994  Bugs fixed in lu1fad (nrank) and lu1ful (ipvt). */
/*     08 Aug 1995  Use ip instead of w as parameter to lu1or3 (for F90). */
/*     13 Sep 2000  TPP and TCP options implemented. */
/*     17 Oct 2000  Fixed troubles due to A = empty matrix (Todd Munson). */
/*     01 Dec 2000  Save Lmax, Umax, etc. after both lu1fad and lu6chk. */
/*                  lu1fad sets them when keepLU = false. */
/*                  lu6chk sets them otherwise, and includes items */
/*                  from the dense LU. */
/*     11 Mar 2001  lu6chk now looks at diag(U) when keepLU = false. */
/*     26 Apr 2002  New TCP implementation using heap routines to */
/*                  store largest element in each column. */
/*                  New workspace arrays Ha, Hj, Hk required. */
/*                  For compatibility, borrow space from a, indc, indr */
/*                  rather than adding new input parameters. */
/*     01 May 2002  lu1den changed to lu1DPP (dense partial  pivoting). */
/*                  lu1DCP implemented       (dense complete pivoting). */
/*                  Both TPP and TCP now switch to dense mode and end. */

/*     Systems Optimization Laboratory, Stanford University. */
/*  --------------------------------------------------------------------- */


/*  INPUT PARAMETERS */

/*  m      (not altered) is the number of rows in A. */
/*  n      (not altered) is the number of columns in A. */
/*  nelem  (not altered) is the number of matrix entries given in */
/*         the arrays a, indc, indr. */
/*  lena   (not altered) is the dimension of  a, indc, indr. */
/*         This should be significantly larger than nelem. */
/*         Typically one should have */
/*            lena > max( 2*nelem, 10*m, 10*n, 10000 ) */
/*         but some applications may need more. */
/*         On machines with virtual memory it is safe to have */
/*         lena "far bigger than necessary", since not all of the */
/*         arrays will be used. */
/*  a      (overwritten) contains entries   Aij  in   a(1:nelem). */
/*  indc   (overwritten) contains the indices i in indc(1:nelem). */
/*  indr   (overwritten) contains the indices j in indr(1:nelem). */

/*  luparm input parameters:                                Typical value */

/*  luparm( 1) = nout     File number for printed messages.         6 */

/*  luparm( 2) = lprint   Print level.                              0 */
/*                   <  0 suppresses output. */
/*                   =  0 gives error messages. */
/*                  >= 10 gives statistics about the LU factors. */
/*                  >= 50 gives debug output from lu1fac */
/*                        (the pivot row and column and the */
/*                        no. of rows and columns involved at */
/*                        each elimination step). */

/*  luparm( 3) = maxcol   lu1fac: maximum number of columns         5 */
/*                        searched allowed in a Markowitz-type */
/*                        search for the next pivot element. */
/*                        For some of the factorization, the */
/*                        number of rows searched is */
/*                        maxrow = maxcol - 1. */

/*  luparm( 6) = 0    =>  TPP: Threshold Partial   Pivoting.        0 */
/*             = 1    =>  TRP: Threshold Rook      Pivoting. */
/*             = 2    =>  TCP: Threshold Complete  Pivoting. */
/*             = 3    =>  TSP: Threshold Symmetric Pivoting. */
/*             = 4    =>  TDP: Threshold Diagonal  Pivoting. */
/*                             (TDP not yet implemented). */
/*                        TRP and TCP are more expensive than TPP but */
/*                        more stable and better at revealing rank. */
/*                        Take care with setting parmlu(1), especially */
/*                        with TCP. */
/*                        NOTE: TSP and TDP are for symmetric matrices */
/*                        that are either definite or quasi-definite. */
/*                        TSP is effectively TRP for symmetric matrices. */
/*                        TDP is effectively TCP for symmetric matrices. */

/*  luparm( 8) = keepLU   lu1fac: keepLU = 1 means the numerical    1 */
/*                        factors will be computed if possible. */
/*                        keepLU = 0 means L and U will be discarded */
/*                        but other information such as the row and */
/*                        column permutations will be returned. */
/*                        The latter option requires less storage. */

/*  parmlu input parameters:                                Typical value */

/*  parmlu( 1) = Ltol1    Max Lij allowed during Factor. */
/*                                                  TPP     10.0 or 100.0 */
/*                                                  TRP      4.0 or  10.0 */
/*                                                  TCP      5.0 or  10.0 */
/*                                                  TSP      4.0 or  10.0 */
/*                        With TRP and TCP (Rook and Complete Pivoting), */
/*                        values less than 25.0 may be expensive */
/*                        on badly scaled data.  However, */
/*                        values less than 10.0 may be needed */
/*                        to obtain a reliable rank-revealing */
/*                        factorization. */
/*  parmlu( 2) = Ltol2    Max Lij allowed during Updates.            10.0 */
/*                        during updates. */
/*  parmlu( 3) = small    Absolute tolerance for       eps**0.8 = 3.0d-13 */
/*                        treating reals as zero. */
/*  parmlu( 4) = Utol1    Absolute tol for flagging    eps**0.67= 3.7d-11 */
/*                        small diagonals of U. */
/*  parmlu( 5) = Utol2    Relative tol for flagging    eps**0.67= 3.7d-11 */
/*                        small diagonals of U. */
/*                        (eps = machine precision) */
/*  parmlu( 6) = Uspace   Factor limiting waste space in  U.      3.0 */
/*                        In lu1fac, the row or column lists */
/*                        are compressed if their length */
/*                        exceeds Uspace times the length of */
/*                        either file after the last compression. */
/*  parmlu( 7) = dens1    The density at which the Markowitz      0.3 */
/*                        pivot strategy should search maxcol */
/*                        columns and no rows. */
/*                        (Use 0.3 unless you are experimenting */
/*                        with the pivot strategy.) */
/*  parmlu( 8) = dens2    the density at which the Markowitz      0.5 */
/*                        strategy should search only 1 column, */
/*                        or (if storage is available) */
/*                        the density at which all remaining */
/*                        rows and columns will be processed */
/*                        by a dense LU code. */
/*                        For example, if dens2 = 0.1 and lena is */
/*                        large enough, a dense LU will be used */
/*                        once more than 10 per cent of the */
/*                        remaining matrix is nonzero. */


/*  OUTPUT PARAMETERS */

/*  a, indc, indr     contain the nonzero entries in the LU factors of A. */
/*         If keepLU = 1, they are in a form suitable for use */
/*         by other parts of the LUSOL package, such as lu6sol. */
/*         U is stored by rows at the start of a, indr. */
/*         L is stored by cols at the end   of a, indc. */
/*         If keepLU = 0, only the diagonals of U are stored, at the */
/*         end of a. */
/*  ip, iq    are the row and column permutations defining the */
/*         pivot order.  For example, row ip(1) and column iq(1) */
/*         defines the first diagonal of U. */
/*  lenc(1:numl0) contains the number of entries in nontrivial */
/*         columns of L (in pivot order). */
/*  lenr(1:m) contains the number of entries in each row of U */
/*         (in original order). */
/*  locc(1:n) = 0 (ready for the LU update routines). */
/*  locr(1:m) points to the beginning of the rows of U in a, indr. */
/*  iploc, iqloc, ipinv, iqinv  are undefined. */
/*  w      indicates singularity as described above. */
/*  inform = 0 if the LU factors were obtained successfully. */
/*         = 1 if U appears to be singular, as judged by lu6chk. */
/*         = 3 if some index pair indc(l), indr(l) lies outside */
/*             the matrix dimensions 1:m , 1:n. */
/*         = 4 if some index pair indc(l), indr(l) duplicates */
/*             another such pair. */
/*         = 7 if the arrays a, indc, indr were not large enough. */
/*             Their length "lena" should be increase to at least */
/*             the value "minlen" given in luparm(13). */
/*         = 8 if there was some other fatal error.  (Shouldn't happen!) */
/*         = 9 if no diagonal pivot could be found with TSP or TDP. */
/*             The matrix must not be sufficiently definite */
/*             or quasi-definite. */

/*  luparm output parameters: */

/*  luparm(10) = inform   Return code from last call to any LU routine. */
/*  luparm(11) = nsing    No. of singularities marked in the */
/*                        output array w(*). */
/*  luparm(12) = jsing    Column index of last singularity. */
/*  luparm(13) = minlen   Minimum recommended value for  lena. */
/*  luparm(14) = maxlen   ? */
/*  luparm(15) = nupdat   No. of updates performed by the lu8 routines. */
/*  luparm(16) = nrank    No. of nonempty rows of U. */
/*  luparm(17) = ndens1   No. of columns remaining when the density of */
/*                        the matrix being factorized reached dens1. */
/*  luparm(18) = ndens2   No. of columns remaining when the density of */
/*                        the matrix being factorized reached dens2. */
/*  luparm(19) = jumin    The column index associated with DUmin. */
/*  luparm(20) = numL0    No. of columns in initial  L. */
/*  luparm(21) = lenL0    Size of initial  L  (no. of nonzeros). */
/*  luparm(22) = lenU0    Size of initial  U. */
/*  luparm(23) = lenL     Size of current  L. */
/*  luparm(24) = lenU     Size of current  U. */
/*  luparm(25) = lrow     Length of row file. */
/*  luparm(26) = ncp      No. of compressions of LU data structures. */
/*  luparm(27) = mersum   lu1fac: sum of Markowitz merit counts. */
/*  luparm(28) = nUtri    lu1fac: triangular rows in U. */
/*  luparm(29) = nLtri    lu1fac: triangular rows in L. */
/*  luparm(30) = */



/*  parmlu output parameters: */

/*  parmlu(10) = Amax     Maximum element in  A. */
/*  parmlu(11) = Lmax     Maximum multiplier in current  L. */
/*  parmlu(12) = Umax     Maximum element in current  U. */
/*  parmlu(13) = DUmax    Maximum diagonal in  U. */
/*  parmlu(14) = DUmin    Minimum diagonal in  U. */
/*  parmlu(15) = Akmax    Maximum element generated at any stage */
/*                        during TCP factorization. */
/*  parmlu(16) = growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax */
/*  parmlu(17) = */
/*  parmlu(18) = */
/*  parmlu(19) = */
/*  parmlu(20) = resid    lu6sol: residual after solve with U or U'. */
/*  ... */
/*  parmlu(30) = */
/*  --------------------------------------------------------------------- */
/*     Grab relevant input parameters. */
    /* Parameter adjustments */
    --ipinv;
    --iqloc;
    --locr;
    --lenr;
    --ip;
    --w;
    --iqinv;
    --iploc;
    --locc;
    --lenc;
    --iq;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    nelem0 = *nelem;
    nout = luparm[1];
    lprint = luparm[2];
    lpiv = luparm[6];
    keeplu = luparm[8] != 0;
    ltol = parmlu[1];
/* Limit on size of Lij */
    small = parmlu[3];
/* Drop tolerance */
    tpp = lpiv == 0;
/* Threshold Partial   Pivoting (normal). */
    trp = lpiv == 1;
/* Threshold Rook      Pivoting */
    tcp = lpiv == 2;
/* Threshold Complete  Pivoting. */
    tsp = lpiv == 3;
/* Threshold Symmetric Pivoting. */
    s_copy(kpiv, "PP", (ftnlen)2, (ftnlen)2);
    s_copy(kpiv + 2, "RP", (ftnlen)2, (ftnlen)2);
    s_copy(kpiv + 4, "CP", (ftnlen)2, (ftnlen)2);
    s_copy(kpiv + 6, "SP", (ftnlen)2, (ftnlen)2);
/*     Initialize output parameters. */
    *inform__ = 0;
    minlen = *nelem + (*m + *n << 1);
    numl0 = 0;
    lenl = 0;
    lenu = 0;
    lrow = 0;
    mersum = 0;
    nutri = *m;
    nltri = 0;
    ndens1 = 0;
    ndens2 = 0;
    nrank = 0;
    nsing = 0;
    jsing = 0;
    jumin = 0;
    amax = 0.;
    lmax = 0.;
    umax = 0.;
    dumax = 0.;
    dumin = 0.;
    akmax = 0.;
    if (*m > *n) {
	*(unsigned char *)mnkey = '>';
    } else if (*m == *n) {
	*(unsigned char *)mnkey = '=';
    } else {
	*(unsigned char *)mnkey = '<';
    }
/*     Float version of dimensions. */
    dm = (doublereal) (*m);
    dn = (doublereal) (*n);
    delem = (doublereal) (*nelem);
/*     Initialize workspace parameters. */
    luparm[26] = 0;
/* ncp */
    if (*lena < minlen) {
	goto L970;
    }
/*     ------------------------------------------------------------------ */
/*     Organize the  aij's  in  a, indc, indr. */
/*     lu1or1  deletes small entries, tests for illegal  i,j's, */
/*             and counts the nonzeros in each row and column. */
/*     lu1or2  reorders the elements of  A  by columns. */
/*     lu1or3  uses the column list to test for duplicate entries */
/*             (same indices  i,j). */
/*     lu1or4  constructs a row list from the column list. */
/*     ------------------------------------------------------------------ */
    lu1or1_(m, n, nelem, lena, &small, &a[1], &indc[1], &indr[1], &lenc[1], &
	    lenr[1], &amax, &numnz, &lerr, inform__);
    if (nout > 0 && lprint >= 10) {
	densty = delem * 100. / (dm * dn);
	io___40.ciunit = nout;
	s_wsfe(&io___40);
	do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	do_fio(&c__1, mnkey, (ftnlen)1);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nelem), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&amax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&densty, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*inform__ != 0) {
	goto L930;
    }
    *nelem = numnz;
    lu1or2_(n, nelem, lena, &a[1], &indc[1], &indr[1], &lenc[1], &locc[1]);
    lu1or3_(m, n, lena, &indc[1], &lenc[1], &locc[1], &ip[1], &lerr, inform__)
	    ;
    if (*inform__ != 0) {
	goto L940;
    }
    lu1or4_(m, n, nelem, lena, &indc[1], &indr[1], &lenc[1], &lenr[1], &locc[
	    1], &locr[1]);
/*     ------------------------------------------------------------------ */
/*     Set up lists of rows and columns with equal numbers of nonzeros, */
/*     using  indc(*)  as workspace. */
/*     ------------------------------------------------------------------ */
    lu1pq1_(m, n, &lenr[1], &ip[1], &iploc[1], &ipinv[1], &indc[*nelem + 1]);
    lu1pq1_(n, m, &lenc[1], &iq[1], &iqloc[1], &iqinv[1], &indc[*nelem + 1]);
/*     ------------------------------------------------------------------ */
/*     For TCP, allocate Ha, Hj, Hk at the end of a, indc, indr. */
/*     Then compute the factorization  A = L*U. */
/*     ------------------------------------------------------------------ */
    if (tpp || tsp) {
	lenh = 1;
	lena2 = *lena;
	loch = *lena;
	lmaxr = 1;
    } else if (trp) {
	lenh = 1;
/* Dummy */
	lena2 = *lena - *m;
/* Reduced length of      a */
	loch = *lena;
/* Dummy */
	lmaxr = lena2 + 1;
/* Start of Amaxr      in a */
    } else if (tcp) {
	lenh = *n;
/* Length of heap */
	lena2 = *lena - lenh;
/* Reduced length of      a, indc, indr */
	loch = lena2 + 1;
/* Start of Ha, Hj, Hk in a, indc, indr */
	lmaxr = 1;
/* Dummy */
    }
    lu1fad_(m, n, nelem, &lena2, &luparm[1], &parmlu[1], &a[1], &indc[1], &
	    indr[1], &ip[1], &iq[1], &lenc[1], &lenr[1], &locc[1], &locr[1], &
	    iploc[1], &iqloc[1], &ipinv[1], &iqinv[1], &w[1], &lenh, &a[loch],
	     &indc[loch], &indr[loch], &a[lmaxr], inform__, &lenl, &lenu, &
	    minlen, &mersum, &nutri, &nltri, &ndens1, &ndens2, &nrank, &lmax, 
	    &umax, &dumax, &dumin, &akmax);
    luparm[16] = nrank;
    luparm[23] = lenl;
    if (*inform__ == 7) {
	goto L970;
    }
    if (*inform__ == 9) {
	goto L985;
    }
    if (*inform__ > 0) {
	goto L980;
    }
    if (keeplu) {
/*        --------------------------------------------------------------- */
/*        The LU factors are at the top of  a, indc, indr, */
/*        with the columns of  L  and the rows of  U  in the order */

/*        ( free )   ... ( u3 ) ( l3 ) ( u2 ) ( l2 ) ( u1 ) ( l1 ). */

/*        Starting with ( l1 ) and ( u1 ), move the rows of  U  to the */
/*        left and the columns of  L  to the right, giving */

/*        ( u1 ) ( u2 ) ( u3 ) ...   ( free )   ... ( l3 ) ( l2 ) ( l1 ). */

/*        Also, set  numl0 = the number of nonempty columns of L. */
/*        --------------------------------------------------------------- */
	lu = 0;
	ll = *lena + 1;
	lm = lena2 + 1;
	ltopl = ll - lenl - lenu;
	lrow = lenu;
	i__1 = nrank;
	for (k = 1; k <= i__1; ++k) {
	    i__ = ip[k];
	    lenuk = -lenr[i__];
	    lenr[i__] = lenuk;
	    j = iq[k];
	    lenlk = -lenc[j] - 1;
	    if (lenlk > 0) {
		++numl0;
		iqloc[numl0] = lenlk;
	    }
	    if (lu + lenuk < ltopl) {
/*              ========================================================= */
/*              There is room to move ( uk ).  Just right-shift ( lk ). */
/*              ========================================================= */
		i__2 = lenlk;
		for (idummy = 1; idummy <= i__2; ++idummy) {
		    --ll;
		    --lm;
		    a[ll] = a[lm];
		    indc[ll] = indc[lm];
		    indr[ll] = indr[lm];
		}
	    } else {
/*              ========================================================= */
/*              There is no room for ( uk ) yet.  We have to */
/*              right-shift the whole of the remaining LU file. */
/*              Note that ( lk ) ends up in the correct place. */
/*              ========================================================= */
		llsave = ll - lenlk;
		nmove = lm - ltopl;
		i__2 = nmove;
		for (idummy = 1; idummy <= i__2; ++idummy) {
		    --ll;
		    --lm;
		    a[ll] = a[lm];
		    indc[ll] = indc[lm];
		    indr[ll] = indr[lm];
		}
		ltopl = ll;
		ll = llsave;
		lm = ll;
	    }
/*           ====================================================== */
/*           Left-shift ( uk ). */
/*           ====================================================== */
	    locr[i__] = lu + 1;
	    l2 = lm - 1;
	    lm -= lenuk;
	    i__2 = l2;
	    for (l = lm; l <= i__2; ++l) {
		++lu;
		a[lu] = a[l];
		indr[lu] = indr[l];
	    }
	}
/*        --------------------------------------------------------------- */
/*        Save the lengths of the nonempty columns of  L, */
/*        and initialize  locc(j)  for the LU update routines. */
/*        --------------------------------------------------------------- */
	i__1 = numl0;
	for (k = 1; k <= i__1; ++k) {
	    lenc[k] = iqloc[k];
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    locc[j] = 0;
	}
/*        --------------------------------------------------------------- */
/*        Test for singularity. */
/*        lu6chk  sets  nsing, jsing, jumin, Lmax, Umax, DUmax, DUmin */
/*        (including entries from the dense LU). */
/*        inform = 1  if there are singularities (nsing gt 0). */
/*        --------------------------------------------------------------- */
	lu6chk_(&c__1, m, n, &w[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &ip[1], &iq[1], &lenc[1], &lenr[1], &locc[
		1], &locr[1], inform__);
	nsing = luparm[11];
	jsing = luparm[12];
	jumin = luparm[19];
	lmax = parmlu[11];
	umax = parmlu[12];
	dumax = parmlu[13];
	dumin = parmlu[14];
    } else {
/*        --------------------------------------------------------------- */
/*        keepLU = 0.  L and U were not kept, just the diagonals of U. */
/*        lu1fac will probably be called again soon with keepLU = .true. */
/*        11 Mar 2001: lu6chk revised.  We can call it with keepLU = 0, */
/*                     but we want to keep Lmax, Umax from lu1fad. */
/*        05 May 2002: Allow for TCP with new lu1DCP.  Diag(U) starts */
/*                     below lena2, not lena.  Need lena2 in next line. */
/*        --------------------------------------------------------------- */
	lu6chk_(&c__1, m, n, &w[1], &lena2, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &ip[1], &iq[1], &lenc[1], &lenr[1], &locc[
		1], &locr[1], inform__);
	nsing = luparm[11];
	jsing = luparm[12];
	jumin = luparm[19];
	dumax = parmlu[13];
	dumin = parmlu[14];
    }
    goto L990;
/*     ------------ */
/*     Error exits. */
/*     ------------ */
L930:
    *inform__ = 3;
    if (lprint >= 0) {
	io___59.ciunit = nout;
	s_wsfe(&io___59);
	do_fio(&c__1, (char *)&lerr, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&indc[lerr], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&indr[lerr], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L990;
L940:
    *inform__ = 4;
    if (lprint >= 0) {
	io___60.ciunit = nout;
	s_wsfe(&io___60);
	do_fio(&c__1, (char *)&lerr, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&indc[lerr], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&indr[lerr], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L990;
L970:
    *inform__ = 7;
    if (lprint >= 0) {
	io___61.ciunit = nout;
	s_wsfe(&io___61);
	do_fio(&c__1, (char *)&(*lena), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&minlen, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L990;
L980:
    *inform__ = 8;
    if (lprint >= 0) {
	io___62.ciunit = nout;
	s_wsfe(&io___62);
	e_wsfe();
    }
    goto L990;
L985:
    *inform__ = 9;
    if (lprint >= 0) {
	io___63.ciunit = nout;
	s_wsfe(&io___63);
	e_wsfe();
    }
/*     Store output parameters. */
L990:
    *nelem = nelem0;
    luparm[10] = *inform__;
    luparm[11] = nsing;
    luparm[12] = jsing;
    luparm[13] = minlen;
    luparm[15] = 0;
    luparm[16] = nrank;
    luparm[17] = ndens1;
    luparm[18] = ndens2;
    luparm[19] = jumin;
    luparm[20] = numl0;
    luparm[21] = lenl;
    luparm[22] = lenu;
    luparm[23] = lenl;
    luparm[24] = lenu;
    luparm[25] = lrow;
    luparm[27] = mersum;
    luparm[28] = nutri;
    luparm[29] = nltri;
    parmlu[10] = amax;
    parmlu[11] = lmax;
    parmlu[12] = umax;
    parmlu[13] = dumax;
    parmlu[14] = dumin;
    parmlu[15] = akmax;
    agrwth = akmax / (amax + 1e-20);
    ugrwth = umax / (amax + 1e-20);
    if (tpp) {
	growth = ugrwth;
    } else {
/* TRP or TCP or TSP */
	growth = agrwth;
    }
    parmlu[16] = growth;
/*     ------------------------------------------------------------------ */
/*     Print statistics for the LU factors. */
/*     ------------------------------------------------------------------ */
    ncp = luparm[26];
    condu = dumax / max(dumin,1e-20);
    dincr = (doublereal) (lenl + lenu - *nelem);
    dincr = dincr * 100. / max(delem,1.);
    avgmer = (doublereal) mersum;
    avgmer /= dm;
    nbump = *m - nutri - nltri;
    if (nout > 0 && lprint >= 10) {
	if (tpp) {
	    io___72.ciunit = nout;
	    s_wsfe(&io___72);
	    do_fio(&c__1, (char *)&avgmer, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lenl, (ftnlen)sizeof(integer));
	    i__1 = lenl + lenu;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dincr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nutri, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lenu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ltol, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&umax, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ugrwth, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nltri, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndens1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lmax, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___73.ciunit = nout;
	    s_wsfe(&io___73);
	    do_fio(&c__1, kpiv + (lpiv << 1), (ftnlen)2);
	    do_fio(&c__1, (char *)&avgmer, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lenl, (ftnlen)sizeof(integer));
	    i__1 = lenl + lenu;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncp, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dincr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nutri, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lenu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ltol, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&umax, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ugrwth, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nltri, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndens1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lmax, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&akmax, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&agrwth, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___74.ciunit = nout;
	s_wsfe(&io___74);
	do_fio(&c__1, (char *)&nbump, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ndens2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dumax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dumin, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&condu, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* lu1fac_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1fac */
/* Subroutine */ int lu1fad_(integer *m, integer *n, integer *nelem, integer *
	lena, integer *luparm, doublereal *parmlu, doublereal *a, integer *
	indc, integer *indr, integer *ip, integer *iq, integer *lenc, integer 
	*lenr, integer *locc, integer *locr, integer *iploc, integer *iqloc, 
	integer *ipinv, integer *iqinv, doublereal *w, integer *lenh, 
	doublereal *ha, integer *hj, integer *hk, doublereal *amaxr, integer *
	inform__, integer *lenl, integer *lenu, integer *minlen, integer *
	mersum, integer *nutri, integer *nltri, integer *ndens1, integer *
	ndens2, integer *nrank, doublereal *lmax, doublereal *umax, 
	doublereal *dumax, doublereal *dumin, doublereal *akmax)
{
    /* Format strings */
    static char fmt_1100[] = "(/1x,a)";
    static char fmt_1200[] = "(\002 nrowu\002,i7,\002   i,jbest\002,2i7,\002"
	    "   nrowd,ncold\002,2i6,\002   i,jmax\002,2i7,\002   aijmax\002,1"
	    "p,e10.2)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer h__, i__, j, k, l;
    static doublereal v;
    static integer lc, ld, kk, ll, lq, lr, ls, lu, lc1, ll1, lq1, lq2, lr1, 
	    lu1;
    static doublereal lij;
    static logical tcp, tpp, trp, tsp;
    static doublereal diag;
    static integer lend, hlen;
    static doublereal amax;
    static integer lcol, lenj, leni, imax, jmax, last, hops;
    static logical ltri;
    static doublereal ltol;
    static integer lpiv;
    static logical utri;
    static integer lrow, nout;
    static doublereal dens1, dens2;
    static integer lfile, lfree;
    static logical dense;
    static integer nfree, ncold, melim, nelim;
    static doublereal abest, small;
    static integer mleft, ilast, jlast, minmn, maxmn, nleft, jbest, ibest, 
	    mbest, nrowd, lsave, limit, lpivr, lpivc, kbest, nfill;
    extern /* Subroutine */ int lu1pq2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), lu1pq3_(
	    integer *, integer *, integer *, integer *, integer *);
    static integer mrank, nsing, nrowu;
    extern /* Subroutine */ int lu1rec_(integer *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *);
    static integer lpivc1, lpivc2;
    extern /* Subroutine */ int lu1gau_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, integer *), lu1mar_(
	    integer *, integer *, integer *, integer *, logical *, doublereal 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *), lu1pen_(integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *), lu1ful_(integer *, integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static logical spars1, spars2;
    extern /* Subroutine */ int lu1mxc_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), lu1slk_(integer *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *);
    static integer lpivr1, lpivr2;
    extern /* Subroutine */ int lu1mrp_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *), lu1msp_(integer *, integer *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *), lu1mxr_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static integer ldiagu;
    extern /* Subroutine */ int hbuild_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *);
    static doublereal aijmax;
    static integer hlenin;
    static doublereal uspace;
    static integer minfre, maxcol;
    static doublereal aijtol;
    static logical keeplu;
    static integer nzchng, nspare;
    static logical denslu;
    static integer nzleft, lfirst, lprint, maxrow;
    extern /* Subroutine */ int hchange_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *), hdelete_(doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___125 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___127 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___137 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___140 = { 0, 0, 0, fmt_1200, 0 };


/*     ------------------------------------------------------------------ */
/*     lu1fad  is a driver for the numerical phase of lu1fac. */
/*     At each stage it computes a column of  L  and a row of  U, */
/*     using a Markowitz criterion to select the pivot element, */
/*     subject to a stability criterion that bounds the elements of  L. */

/*     00 Jan 1986  Version documented in LUSOL paper: */
/*                  Gill, Murray, Saunders and Wright (1987), */
/*                  Maintaining LU factors of a general sparse matrix, */
/*                  Linear algebra and its applications 88/89, 239-270. */

/*     02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest */
/*                  element in each column is now kept at the start of */
/*                  the column, i.e. in position locc(j) of a and indc. */
/*                  This should speed up the Markowitz searches. */
/*                  To save time on highly triangular matrices, we wait */
/*                  until there are no further columns of length 1 */
/*                  before setting and maintaining that property. */

/*     12 Apr 1989  ipinv and iqinv added (inverses of ip and iq) */
/*                  to save searching ip and iq for rows and columns */
/*                  altered in each elimination step.  (Used in lu1pq2) */

/*     19 Apr 1989  Code segmented to reduce its size. */
/*                  lu1gau does most of the Gaussian elimination work. */
/*                  lu1mar does just the Markowitz search. */
/*                  lu1mxc moves biggest elements to top of columns. */
/*                  lu1pen deals with pending fill-in in the row list. */
/*                  lu1pq2 updates the row and column permutations. */

/*     26 Apr 1989  maxtie replaced by maxcol, maxrow in the Markowitz */
/*                  search.  maxcol, maxrow change as density increases. */

/*     25 Oct 1993  keepLU implemented. */

/*     07 Feb 1994  Exit main loop early to finish off with a dense LU. */
/*                  densLU tells lu1fad whether to do it. */
/*     21 Dec 1994  Bug fixed.  nrank was wrong after the call to lu1ful. */
/*     12 Nov 1999  A parallel version of dcopy gave trouble in lu1ful */
/*                  during left-shift of dense matrix D within a(*). */
/*                  Fixed this unexpected problem here in lu1fad */
/*                  by making sure the first and second D don't overlap. */

/*     13 Sep 2000  TCP (Threshold Complete Pivoting) implemented. */
/*                  lu2max added */
/*                  (finds aijmax from biggest elems in each col). */
/*                  Utri, Ltri and Spars1 phases apply. */
/*                  No switch to Dense CP yet.  (Only TPP switches.) */
/*     14 Sep 2000  imax needed to remember row containing aijmax. */
/*     22 Sep 2000  For simplicity, lu1mxc always fixes */
/*                  all modified cols. */
/*                  (TPP spars2 used to fix just the first maxcol cols.) */
/*     08 Nov 2000: Speed up search for aijmax. */
/*                  Don't need to search all columns if the elimination */
/*                  didn't alter the col containing the current aijmax. */
/*     21 Nov 2000: lu1slk implemented for Utri phase with TCP */
/*                  to guard against deceptive triangular matrices. */
/*                  (Utri used to have aijtol >= 0.9999 to include */
/*                  slacks, but this allows other 1s to be accepted.) */
/*                  Utri now accepts slacks, but applies normal aijtol */
/*                  test to other pivots. */
/*     28 Nov 2000: TCP with empty cols must call lu1mxc and lu2max */
/*                  with ( lq1, n, ... ), not just ( 1, n, ... ). */
/*     23 Mar 2001: lu1fad bug with TCP. */
/*                  A col of length 1 might not be accepted as a pivot. */
/*                  Later it appears in a pivot row and temporarily */
/*                  has length 0 (when pivot row is removed */
/*                  but before the column is filled in).  If it is the */
/*                  last column in storage, the preceding col also thinks */
/*                  it is "last".  Trouble arises when the preceding col */
/*                  needs fill-in -- it overlaps the real "last" column. */
/*                  (Very rarely, same trouble might have happened if */
/*                  the drop tolerance caused columns to have length 0.) */

/*                  Introduced ilast to record the last row in row file, */
/*                             jlast to record the last col in col file. */
/*                  lu1rec returns ilast = indr(lrow + 1) */
/*                              or jlast = indc(lcol + 1). */
/*                  (Should be an output parameter, but didn't want to */
/*                  alter lu1rec's parameter list.) */
/*                  lu1rec also treats empty rows or cols safely. */
/*                  (Doesn't eliminate them!) */

/*     26 Apr 2002: Heap routines added for TCP. */
/*                  lu2max no longer needed. */
/*                  imax, jmax used only for printing. */
/*     01 May 2002: lu1DCP implemented (dense complete pivoting). */
/*                  Both TPP and TCP now switch to dense LU */
/*                  when density exceeds dens2. */
/*     06 May 2002: In dense mode, store diag(U) in natural order. */
/*     09 May 2002: lu1mCP implemented (Markowitz TCP via heap). */
/*     11 Jun 2002: lu1mRP implemented (Markowitz TRP). */
/*     28 Jun 2002: Fixed call to lu1mxr. */
/*     14 Dec 2002: lu1mSP implemented (Markowitz TSP). */
/*     15 Dec 2002: Both TPP and TSP can grab cols of length 1 */
/*                  during Utri. */
/*     19 Dec 2004: Hdelete(...) has new input argument Hlenin. */
/*     26 Mar 2006: lu1fad returns nrank  = min( mrank, nrank ) */
/*                  and ignores nsing from lu1ful */
/*     26 Mar 2006: Allow for empty columns before calling Hbuild. */

/*     Systems Optimization Laboratory, Stanford University. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Local variables */
/*     --------------- */

/*     lcol   is the length of the column file.  It points to the last */
/*            nonzero in the column list. */
/*     lrow   is the analogous quantity for the row file. */
/*     lfile  is the file length (lcol or lrow) after the most recent */
/*            compression of the column list or row list. */
/*     nrowd  and  ncold  are the number of rows and columns in the */
/*            matrix defined by the pivot column and row.  They are the */
/*            dimensions of the submatrix D being altered at this stage. */
/*     melim  and  nelim  are the number of rows and columns in the */
/*            same matrix D, excluding the pivot column and row. */
/*     mleft  and  nleft  are the number of rows and columns */
/*            still left to be factored. */
/*     nzchng is the increase in nonzeros in the matrix that remains */
/*            to be factored after the current elimination */
/*            (usually negative). */
/*     nzleft is the number of nonzeros still left to be factored. */
/*     nspare is the space we leave at the end of the last row or */
/*            column whenever a row or column is being moved to the end */
/*            of its file.  nspare = 1 or 2 might help reduce the */
/*            number of file compressions when storage is tight. */

/*     The row and column ordering permutes A into the form */

/*                        ------------------------ */
/*                         \                     | */
/*                          \         U1         | */
/*                           \                   | */
/*                            -------------------- */
/*                            |\ */
/*                            | \ */
/*                            |  \ */
/*            P A Q   =       |   \ */
/*                            |    \ */
/*                            |     -------------- */
/*                            |     |            | */
/*                            |     |            | */
/*                            | L1  |     A2     | */
/*                            |     |            | */
/*                            |     |            | */
/*                            -------------------- */

/*     where the block A2 is factored as  A2 = L2 U2. */
/*     The phases of the factorization are as follows. */

/*     Utri   is true when U1 is being determined. */
/*            Any column of length 1 is accepted immediately (if TPP). */

/*     Ltri   is true when L1 is being determined. */
/*            lu1mar exits as soon as an acceptable pivot is found */
/*            in a row of length 1. */

/*     spars1 is true while the density of the (modified) A2 is less */
/*            than the parameter dens1 = parmlu(7) = 0.3 say. */
/*            lu1mar searches maxcol columns and maxrow rows, */
/*            where  maxcol = luparm(3),  maxrow = maxcol - 1. */
/*            lu1mxc is used to keep the biggest element at the top */
/*            of all remaining columns. */

/*     spars2 is true while the density of the modified A2 is less */
/*            than the parameter dens2 = parmlu(8) = 0.6 say. */
/*            lu1mar searches maxcol columns and no rows. */
/*            lu1mxc could fix up only the first maxcol cols (with TPP). */
/*            22 Sep 2000:  For simplicity, lu1mxc fixes all */
/*                          modified cols. */

/*     dense  is true once the density of A2 reaches dens2. */
/*            lu1mar searches only 1 column (the shortest). */
/*            lu1mxc could fix up only the first column (with TPP). */
/*            22 Sep 2000:  For simplicity, lu1mxc fixes all */
/*                          modified cols. */

/*     ------------------------------------------------------------------ */
/*     To eliminate timings, comment out all lines containing "time". */
/*     ------------------------------------------------------------------ */
/*     integer            eltime, mktime */
/*     call timer ( 'start ', 3 ) */
/*     ntime  = (n / 4.0) */
    /* Parameter adjustments */
    --amaxr;
    --ipinv;
    --iqloc;
    --locr;
    --lenr;
    --ip;
    --w;
    --iqinv;
    --iploc;
    --locc;
    --lenc;
    --iq;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;
    --hk;
    --hj;
    --ha;

    /* Function Body */
    nout = luparm[1];
    lprint = luparm[2];
    maxcol = luparm[3];
    lpiv = luparm[6];
    keeplu = luparm[8] != 0;
    tpp = lpiv == 0;
/* Threshold Partial   Pivoting (normal). */
    trp = lpiv == 1;
/* Threshold Rook      Pivoting */
    tcp = lpiv == 2;
/* Threshold Complete  Pivoting. */
    tsp = lpiv == 3;
/* Threshold Symmetric Pivoting. */
    denslu = FALSE_;
    maxrow = maxcol - 1;
    ilast = *m;
/* Assume row m is last in the row file. */
    jlast = *n;
/* Assume col n is last in the col file. */
    lfile = *nelem;
    lrow = *nelem;
    lcol = *nelem;
    minmn = min(*m,*n);
    maxmn = max(*m,*n);
    nzleft = *nelem;
    nspare = 1;
    if (keeplu) {
	lu1 = *lena + 1;
    } else {
/*        Store only the diagonals of U in the top of memory. */
	ldiagu = *lena - *n;
	lu1 = ldiagu + 1;
    }
    ltol = parmlu[1];
    small = parmlu[3];
    uspace = parmlu[6];
    dens1 = parmlu[7];
    dens2 = parmlu[8];
    utri = TRUE_;
    ltri = FALSE_;
    spars1 = FALSE_;
    spars2 = FALSE_;
    dense = FALSE_;
/*     Check parameters. */
    ltol = max(ltol,1.0001);
    dens1 = min(dens1,dens2);
/*     Initialize output parameters. */
/*     lenL, lenU, minlen, mersum, nUtri, nLtri, ndens1, ndens2, nrank */
/*     are already initialized by lu1fac. */
    *lmax = 0.;
    *umax = 0.;
    *dumax = 0.;
    *dumin = 1e20;
    if (*nelem == 0) {
	*dumin = 0.;
    }
    *akmax = 0.;
    hops = 0;
/* More initialization. */
    if (tpp || tsp) {
/* Don't worry yet about lu1mxc. */
	aijmax = 0.;
	aijtol = 0.;
	hlen = 1;
    } else {
/* Move biggest element to top of each column. */
/* Set w(*) to mark slack columns (unit vectors). */
/* TRP or TCP */
	lu1mxc_(&c__1, n, &iq[1], &a[1], &indc[1], &lenc[1], &locc[1]);
	lu1slk_(m, n, lena, &iq[1], &iqloc[1], &a[1], &locc[1], &w[1]);
    }
    if (trp) {
/* Find biggest element in each row. */
	lu1mxr_(&c__1, m, &ip[1], &amaxr[1], &a[1], &indc[1], &lenc[1], &locc[
		1], &indr[1], &lenr[1], &locr[1]);
    }
    if (tcp) {
/* Set Ha(1:Hlen) = biggest element in each column, */
/*     Hj(1:Hlen) = corresponding column indices. */
	hlen = 0;
	i__1 = *n;
	for (kk = 1; kk <= i__1; ++kk) {
	    ++hlen;
	    j = iq[kk];
	    if (lenc[j] > 0) {
		lc = locc[j];
		amax = (d__1 = a[lc], abs(d__1));
	    } else {
		amax = 0.;
	    }
	    ha[hlen] = amax;
	    hj[hlen] = j;
	    hk[j] = hlen;
	}
/* Build the heap, creating new Ha, Hj and setting Hk(1:Hlen). */
	hbuild_(&ha[1], &hj[1], &hk[1], &hlen, &hlen, &hops);
    }
/*     ------------------------------------------------------------------ */
/*     Start of main loop. */
/*     ------------------------------------------------------------------ */
    mleft = *m + 1;
    nleft = *n + 1;
    i__1 = minmn;
    for (nrowu = 1; nrowu <= i__1; ++nrowu) {
/*        mktime = (nrowu / ntime) + 4 */
/*        eltime = (nrowu / ntime) + 9 */
	--mleft;
	--nleft;
/*        Bail out if there are no nonzero rows left. */
	if (iploc[1] > *m) {
	    goto L900;
	}
/* For TCP, the largest Aij is at the top of the heap. */
	if (tcp) {
	    aijmax = ha[1];
/* Marvelously easy ! */
	    *akmax = max(*akmax,aijmax);
	    aijtol = aijmax / ltol;
	}
/*        =============================================================== */
/*        Find a suitable pivot element. */
/*        =============================================================== */
	if (utri) {
/*           ------------------------------------------------------------ */
/*           So far all columns have had length 1. */
/*           We are still looking for the (backward) triangular part of A */
/*           that forms the first rows and columns of U. */
/*           ------------------------------------------------------------ */
	    lq1 = iqloc[1];
	    lq2 = *n;
	    if (*m > 1) {
		lq2 = iqloc[2] - 1;
	    }
	    if (lq1 <= lq2) {
/* There are more cols of length 1. */
		if (tpp || tsp) {
		    jbest = iq[lq1];
/* Grab the first one. */
		} else {
/* TRP or TCP    ! Scan all columns of length 1. */
		    jbest = 0;
		    i__2 = lq2;
		    for (lq = lq1; lq <= i__2; ++lq) {
			j = iq[lq];
			if (w[j] > 0.) {
/* Accept a slack */
			    jbest = j;
			    goto L250;
			}
			lc = locc[j];
			amax = (d__1 = a[lc], abs(d__1));
			if (trp) {
			    i__ = indc[lc];
			    aijtol = amaxr[i__] / ltol;
			}
			if (amax >= aijtol) {
			    jbest = j;
			    goto L250;
			}
		    }
		}
L250:
		if (jbest > 0) {
		    lc = locc[jbest];
		    ibest = indc[lc];
		    mbest = 0;
		    goto L300;
		}
	    }
/*           This is the end of the U triangle. */
/*           We will not return to this part of the code. */
/*           TPP and TSP call lu1mxc for the first time */
/*           (to move biggest element to top of each column). */
	    if (lprint >= 50) {
		io___125.ciunit = nout;
		s_wsfe(&io___125);
		do_fio(&c__1, "Utri ended.  spars1 = true", (ftnlen)26);
		e_wsfe();
	    }
	    utri = FALSE_;
	    ltri = TRUE_;
	    spars1 = TRUE_;
	    *nutri = nrowu - 1;
	    if (tpp || tsp) {
		lu1mxc_(&lq1, n, &iq[1], &a[1], &indc[1], &lenc[1], &locc[1]);
	    }
	}
	if (spars1) {
/*           ------------------------------------------------------------ */
/*           Perform a Markowitz search. */
/*           Search cols of length 1, then rows of length 1, */
/*           then   cols of length 2, then rows of length 2, etc. */
/*           ------------------------------------------------------------ */
/*           call timer ( 'start ', mktime ) */
/* if (TPP) then ! 12 Jun 2002: Next line disables lu1mCP below */
	    if (tpp || tcp) {
		lu1mar_(m, n, lena, &maxmn, &tcp, &aijtol, &ltol, &maxcol, &
			maxrow, &ibest, &jbest, &mbest, &a[1], &indc[1], &
			indr[1], &ip[1], &iq[1], &lenc[1], &lenr[1], &locc[1],
			 &locr[1], &iploc[1], &iqloc[1]);
	    } else if (trp) {
		lu1mrp_(m, n, lena, &maxmn, &ltol, &maxcol, &maxrow, &ibest, &
			jbest, &mbest, &a[1], &indc[1], &indr[1], &ip[1], &iq[
			1], &lenc[1], &lenr[1], &locc[1], &locr[1], &iploc[1],
			 &iqloc[1], &amaxr[1]);
/*           else if (TCP) then ! Disabled by test above */
/*              call lu1mCP( m    , n     , lena  , aijtol, */
/*    &              ibest, jbest , mbest , */
/*    &              a    , indc  , indr  , */
/*    &              lenc , lenr  , locc  , */
/*    &              Hlen , Ha    , Hj    ) */
	    } else if (tsp) {
		lu1msp_(m, n, lena, &maxmn, &ltol, &maxcol, &ibest, &jbest, &
			mbest, &a[1], &indc[1], &iq[1], &locc[1], &iqloc[1]);
		if (ibest == 0) {
		    goto L990;
		}
	    }
/*           call timer ( 'finish', mktime ) */
	    if (ltri) {
/*              So far all rows have had length 1. */
/*              We are still looking for the (forward) triangle of A */
/*              that forms the first rows and columns of L. */
		if (mbest > 0) {
		    ltri = FALSE_;
		    *nltri = nrowu - 1 - *nutri;
		    if (lprint >= 50) {
			io___126.ciunit = nout;
			s_wsfe(&io___126);
			do_fio(&c__1, "Ltri ended.", (ftnlen)11);
			e_wsfe();
		    }
		}
	    } else {
/*              See if what's left is as dense as dens1. */
		if ((doublereal) nzleft >= dens1 * mleft * nleft) {
		    spars1 = FALSE_;
		    spars2 = TRUE_;
		    *ndens1 = nleft;
		    maxrow = 0;
		    if (lprint >= 50) {
			io___127.ciunit = nout;
			s_wsfe(&io___127);
			do_fio(&c__1, "spars1 ended.  spars2 = true", (ftnlen)
				28);
			e_wsfe();
		    }
		}
	    }
	} else if (spars2 || dense) {
/*           ------------------------------------------------------------ */
/*           Perform a restricted Markowitz search, */
/*           looking at only the first maxcol columns.  (maxrow = 0.) */
/*           ------------------------------------------------------------ */
/*           call timer ( 'start ', mktime ) */
/* if (TPP) then ! 12 Jun 2002: Next line disables lu1mCP below */
	    if (tpp || tcp) {
		lu1mar_(m, n, lena, &maxmn, &tcp, &aijtol, &ltol, &maxcol, &
			maxrow, &ibest, &jbest, &mbest, &a[1], &indc[1], &
			indr[1], &ip[1], &iq[1], &lenc[1], &lenr[1], &locc[1],
			 &locr[1], &iploc[1], &iqloc[1]);
	    } else if (trp) {
		lu1mrp_(m, n, lena, &maxmn, &ltol, &maxcol, &maxrow, &ibest, &
			jbest, &mbest, &a[1], &indc[1], &indr[1], &ip[1], &iq[
			1], &lenc[1], &lenr[1], &locc[1], &locr[1], &iploc[1],
			 &iqloc[1], &amaxr[1]);
/*           else if (TCP) then ! Disabled by test above */
/*              call lu1mCP( m    , n     , lena  , aijtol, */
/*    &              ibest, jbest , mbest , */
/*    &              a    , indc  , indr  , */
/*    &              lenc , lenr  , locc  , */
/*    &              Hlen , Ha    , Hj    ) */
	    } else if (tsp) {
		lu1msp_(m, n, lena, &maxmn, &ltol, &maxcol, &ibest, &jbest, &
			mbest, &a[1], &indc[1], &iq[1], &locc[1], &iqloc[1]);
		if (ibest == 0) {
		    goto L985;
		}
	    }
/*           call timer ( 'finish', mktime ) */
/*           See if what's left is as dense as dens2. */
	    if (spars2) {
		if ((doublereal) nzleft >= dens2 * mleft * nleft) {
		    spars2 = FALSE_;
		    dense = TRUE_;
		    *ndens2 = nleft;
		    maxcol = 1;
		    if (lprint >= 50) {
			io___128.ciunit = nout;
			s_wsfe(&io___128);
			do_fio(&c__1, "spars2 ended.  dense = true", (ftnlen)
				27);
			e_wsfe();
		    }
		}
	    }
	}
/*        --------------------------------------------------------------- */
/*        See if we can finish quickly. */
/*        --------------------------------------------------------------- */
	if (dense) {
	    lend = mleft * nleft;
	    nfree = lu1 - 1;
	    if (nfree >= lend << 1) {
/*              There is room to treat the remaining matrix as */
/*              a dense matrix D. */
/*              We may have to compress the column file first. */
/*              12 Nov 1999: D used to be put at the */
/*                           beginning of free storage (lD = lcol + 1). */
/*                           Now put it at the end     (lD = lu1 - lenD) */
/*                           so the left-shift in lu1ful will not */
/*                           involve overlapping storage */
/*                           (fatal with parallel dcopy). */

		denslu = TRUE_;
		*ndens2 = nleft;
		ld = lu1 - lend;
		if (lcol >= ld) {
		    lu1rec_(n, &c_true, &luparm[1], &lcol, lena, &a[1], &indc[
			    1], &lenc[1], &locc[1]);
		    lfile = lcol;
		    jlast = indc[lcol + 1];
		}
		goto L900;
	    }
	}
/*        =============================================================== */
/*        The best  aij  has been found. */
/*        The pivot row  ibest  and the pivot column  jbest */
/*        Define a dense matrix  D  of size  nrowd  by  ncold. */
/*        =============================================================== */
L300:
	ncold = lenr[ibest];
	nrowd = lenc[jbest];
	melim = nrowd - 1;
	nelim = ncold - 1;
	*mersum += mbest;
	*lenl += melim;
	*lenu += ncold;
	if (lprint >= 50) {
	    if (nrowu == 1) {
		io___136.ciunit = nout;
		s_wsfe(&io___136);
		do_fio(&c__1, "lu1fad debug:", (ftnlen)13);
		e_wsfe();
	    }
	    if (tpp || trp || tsp) {
		io___137.ciunit = nout;
		s_wsfe(&io___137);
		do_fio(&c__1, (char *)&nrowu, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ibest, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jbest, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nrowd, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncold, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
/* TCP */
		jmax = hj[1];
		imax = indc[locc[jmax]];
		io___140.ciunit = nout;
		s_wsfe(&io___140);
		do_fio(&c__1, (char *)&nrowu, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ibest, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jbest, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nrowd, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncold, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&imax, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jmax, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&aijmax, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
/*        =============================================================== */
/*        Allocate storage for the next column of  L  and next row of  U. */
/*        Initially the top of a, indc, indr are used as follows: */

/*                   ncold       melim       ncold        melim */

/*        a      |...........|...........|ujbest..ujn|li1......lim| */

/*        indc   |...........|  lenr(i)  |  lenc(j)  |  markl(i)  | */

/*        indr   |...........| iqloc(i)  |  jfill(j) |  ifill(i)  | */

/*              ^           ^             ^           ^            ^ */
/*              lfree   lsave             lu1         ll1          oldlu1 */

/*        Later the correct indices are inserted: */

/*        indc   |           |           |           |i1........im| */

/*        indr   |           |           |jbest....jn|ibest..ibest| */

/*        =============================================================== */
	if (keeplu) {
/*           relax */
	} else {
/*           Always point to the top spot. */
/*           Only the current column of L and row of U will */
/*           take up space, overwriting the previous ones. */
	    lu1 = ldiagu + 1;
	}
	ll1 = lu1 - melim;
	lu1 = ll1 - ncold;
	lsave = lu1 - nrowd;
	lfree = lsave - ncold;
/*        Make sure the column file has room. */
/*        Also force a compression if its length exceeds a certain limit. */
	limit = (integer) (uspace * lfile + *m + *n + 1000);
	minfre = ncold + melim;
	nfree = lfree - lcol;
	if (nfree < minfre || lcol > limit) {
	    lu1rec_(n, &c_true, &luparm[1], &lcol, lena, &a[1], &indc[1], &
		    lenc[1], &locc[1]);
	    lfile = lcol;
	    jlast = indc[lcol + 1];
	    nfree = lfree - lcol;
	    if (nfree < minfre) {
		goto L970;
	    }
	}
/*        Make sure the row file has room. */
	minfre = melim + ncold;
	nfree = lfree - lrow;
	if (nfree < minfre || lrow > limit) {
	    lu1rec_(m, &c_false, &luparm[1], &lrow, lena, &a[1], &indr[1], &
		    lenr[1], &locr[1]);
	    lfile = lrow;
	    ilast = indr[lrow + 1];
	    nfree = lfree - lrow;
	    if (nfree < minfre) {
		goto L970;
	    }
	}
/*        =============================================================== */
/*        Move the pivot element to the front of its row */
/*        and to the top of its column. */
/*        =============================================================== */
	lpivr = locr[ibest];
	lpivr1 = lpivr + 1;
	lpivr2 = lpivr + nelim;
	i__2 = lpivr2;
	for (l = lpivr; l <= i__2; ++l) {
	    if (indr[l] == jbest) {
		goto L335;
	    }
/* L330: */
	}
L335:
	indr[l] = indr[lpivr];
	indr[lpivr] = jbest;
	lpivc = locc[jbest];
	lpivc1 = lpivc + 1;
	lpivc2 = lpivc + melim;
	i__2 = lpivc2;
	for (l = lpivc; l <= i__2; ++l) {
	    if (indc[l] == ibest) {
		goto L345;
	    }
/* L340: */
	}
L345:
	indc[l] = indc[lpivc];
	indc[lpivc] = ibest;
	abest = a[l];
	a[l] = a[lpivc];
	a[lpivc] = abest;
	if (keeplu) {
/*           relax */
	} else {
/*           Store just the diagonal of U, in natural order. */
/* !!         a(ldiagU + nrowu) = abest ! This was in pivot order. */
	    a[ldiagu + jbest] = abest;
	}
/* ============================================================== */
/* Delete pivot col from heap. */
/* Hk tells us where it is in the heap. */
/* ============================================================== */
	if (tcp) {
	    kbest = hk[jbest];
	    hlenin = hlen;
	    hdelete_(&ha[1], &hj[1], &hk[1], &hlenin, &hlen, n, &kbest, &h__);
	    hops += h__;
	}
/*        =============================================================== */
/*        Delete the pivot row from the column file */
/*        and store it as the next row of  U. */
/*        set  indr(lu) = 0     to initialize jfill ptrs on columns of D, */
/*             indc(lu) = lenj  to save the original column lengths. */
/*        =============================================================== */
	a[lu1] = abest;
	indr[lu1] = jbest;
	indc[lu1] = nrowd;
	lu = lu1;
	diag = abs(abest);
	*umax = max(*umax,diag);
	*dumax = max(*dumax,diag);
	*dumin = min(*dumin,diag);
	i__2 = lpivr2;
	for (lr = lpivr1; lr <= i__2; ++lr) {
	    ++lu;
	    j = indr[lr];
	    lenj = lenc[j];
	    lenc[j] = lenj - 1;
	    lc1 = locc[j];
	    last = lc1 + lenc[j];
	    i__3 = last;
	    for (l = lc1; l <= i__3; ++l) {
		if (indc[l] == ibest) {
		    goto L355;
		}
/* L350: */
	    }
L355:
	    a[lu] = a[l];
	    indr[lu] = 0;
	    indc[lu] = lenj;
/* Computing MAX */
	    d__2 = *umax, d__3 = (d__1 = a[lu], abs(d__1));
	    *umax = max(d__2,d__3);
	    a[l] = a[last];
	    indc[l] = indc[last];
	    indc[last] = 0;
/* ???        if (j .eq. jlast) lcol = lcol - 1 */
/* Free entry */
/* L360: */
	}
/*        =============================================================== */
/*        Delete the pivot column from the row file */
/*        and store the nonzeros of the next column of  L. */
/*        Set  indc(ll) = 0     to initialize markl(*) markers, */
/*             indr(ll) = 0     to initialize ifill(*) row fill-in cntrs, */
/*             indc(ls) = leni  to save the original row lengths, */
/*             indr(ls) = iqloc(i)    to save parts of  iqloc(*), */
/*             iqloc(i) = lsave - ls  to point to the nonzeros of  L */
/*                      = -1, -2, -3, ... in mark(*). */
/*        =============================================================== */
	indc[lsave] = ncold;
	if (melim == 0) {
	    goto L700;
	}
	ll = ll1 - 1;
	ls = lsave;
	abest = 1. / abest;
	i__2 = lpivc2;
	for (lc = lpivc1; lc <= i__2; ++lc) {
	    ++ll;
	    ++ls;
	    i__ = indc[lc];
	    leni = lenr[i__];
	    lenr[i__] = leni - 1;
	    lr1 = locr[i__];
	    last = lr1 + lenr[i__];
	    i__3 = last;
	    for (l = lr1; l <= i__3; ++l) {
		if (indr[l] == jbest) {
		    goto L385;
		}
/* L380: */
	    }
L385:
	    indr[l] = indr[last];
	    indr[last] = 0;
/* ???        if (i .eq. ilast) lrow = lrow - 1 */
/* Free entry */
	    a[ll] = -a[lc] * abest;
	    lij = (d__1 = a[ll], abs(d__1));
	    *lmax = max(*lmax,lij);
/* !!!! DEBUG */
/*           if (Lij .gt. Ltol) then */
/*              write( *  ,*) ' Big Lij!!!', nrowu */
/*              write(nout,*) ' Big Lij!!!', nrowu */
/*           end if */
	    indc[ll] = 0;
	    indr[ll] = 0;
	    indc[ls] = leni;
	    indr[ls] = iqloc[i__];
	    iqloc[i__] = lsave - ls;
/* L390: */
	}
/*        =============================================================== */
/*        Do the Gaussian elimination. */
/*        This involves adding a multiple of the pivot column */
/*        to all other columns in the pivot row. */

/*        Sometimes more than one call to lu1gau is needed to allow */
/*        compression of the column file. */
/*        lfirst  says which column the elimination should start with. */
/*        minfre  is a bound on the storage needed for any one column. */
/*        lu      points to off-diagonals of u. */
/*        nfill   keeps track of pending fill-in in the row file. */
/*        =============================================================== */
	if (nelim == 0) {
	    goto L700;
	}
	lfirst = lpivr1;
	minfre = mleft + nspare;
	lu = 1;
	nfill = 0;
/* 400    call timer ( 'start ', eltime ) */
L400:
	lu1gau_(m, &melim, &ncold, &nspare, &small, &lpivc1, &lpivc2, &lfirst,
		 &lpivr2, &lfree, &minfre, &ilast, &jlast, &lrow, &lcol, &lu, 
		&nfill, &a[1], &indc[1], &indr[1], &lenc[1], &lenr[1], &locc[
		1], &locr[1], &iqloc[1], &a[ll1], &indc[ll1], &a[lu1], &indr[
		ll1], &indr[lu1]);
/*        call timer ( 'finish', eltime ) */
	if (lfirst > 0) {
/*           The elimination was interrupted. */
/*           Compress the column file and try again. */
/*           lfirst, lu and nfill have appropriate new values. */
	    lu1rec_(n, &c_true, &luparm[1], &lcol, lena, &a[1], &indc[1], &
		    lenc[1], &locc[1]);
	    lfile = lcol;
	    jlast = indc[lcol + 1];
	    lpivc = locc[jbest];
	    lpivc1 = lpivc + 1;
	    lpivc2 = lpivc + melim;
	    nfree = lfree - lcol;
	    if (nfree < minfre) {
		goto L970;
	    }
	    goto L400;
	}
/*        =============================================================== */
/*        The column file has been fully updated. */
/*        Deal with any pending fill-in in the row file. */
/*        =============================================================== */
	if (nfill > 0) {
/*           Compress the row file if necessary. */
/*           lu1gau has set nfill to be the number of pending fill-ins */
/*           plus the current length of any rows that need to be moved. */
	    minfre = nfill;
	    nfree = lfree - lrow;
	    if (nfree < minfre) {
		lu1rec_(m, &c_false, &luparm[1], &lrow, lena, &a[1], &indr[1],
			 &lenr[1], &locr[1]);
		lfile = lrow;
		ilast = indr[lrow + 1];
		lpivr = locr[ibest];
		lpivr1 = lpivr + 1;
		lpivr2 = lpivr + nelim;
		nfree = lfree - lrow;
		if (nfree < minfre) {
		    goto L970;
		}
	    }
/*           Move rows that have pending fill-in to end of the row file. */
/*           Then insert the fill-in. */
	    lu1pen_(m, &melim, &ncold, &nspare, &ilast, &lpivc1, &lpivc2, &
		    lpivr1, &lpivr2, &lrow, &lenc[1], &lenr[1], &locc[1], &
		    locr[1], &indc[1], &indr[1], &indr[ll1], &indr[lu1]);
	}
/*        =============================================================== */
/*        Restore the saved values of  iqloc. */
/*        Insert the correct indices for the col of L and the row of U. */
/*        =============================================================== */
L700:
	lenr[ibest] = 0;
	lenc[jbest] = 0;
	ll = ll1 - 1;
	ls = lsave;
	i__2 = lpivc2;
	for (lc = lpivc1; lc <= i__2; ++lc) {
	    ++ll;
	    ++ls;
	    i__ = indc[lc];
	    iqloc[i__] = indr[ls];
	    indc[ll] = i__;
	    indr[ll] = ibest;
/* L710: */
	}
	lu = lu1 - 1;
	i__2 = lpivr2;
	for (lr = lpivr; lr <= i__2; ++lr) {
	    ++lu;
	    indr[lu] = indr[lr];
/* L720: */
	}
/*        =============================================================== */
/*        Free the space occupied by the pivot row */
/*        and update the column permutation. */
/*        Then free the space occupied by the pivot column */
/*        and update the row permutation. */

/*        nzchng is found in both calls to lu1pq2, but we use it only */
/*        after the second. */
/*        =============================================================== */
	lu1pq2_(&ncold, &nzchng, &indr[lpivr], &indc[lu1], &lenc[1], &iqloc[1]
		, &iq[1], &iqinv[1]);
	lu1pq2_(&nrowd, &nzchng, &indc[lpivc], &indc[lsave], &lenr[1], &iploc[
		1], &ip[1], &ipinv[1]);
	nzleft += nzchng;
/*        =============================================================== */
/*        lu1mxr resets Amaxr(i) in each modified row i. */
/*        lu1mxc moves the largest aij to the top of each modified col j. */
/*        28 Jun 2002: Note that cols of L have an implicit diag of 1.0, */
/*                     so lu1mxr is called with ll1, not ll1+1, whereas */
/*                        lu1mxc is called with          lu1+1. */
/*        =============================================================== */
	if (utri && tpp) {
/* Relax -- we're not keeping big elements at the top yet. */
	} else {
	    if (trp && melim > 0) {
		lu1mxr_(&ll1, &ll, &indc[1], &amaxr[1], &a[1], &indc[1], &
			lenc[1], &locc[1], &indr[1], &lenr[1], &locr[1]);
	    }
	    if (nelim > 0) {
		i__2 = lu1 + 1;
		lu1mxc_(&i__2, &lu, &indr[1], &a[1], &indc[1], &lenc[1], &
			locc[1]);
		if (tcp) {
/* Update modified columns in heap */
		    i__2 = lu;
		    for (kk = lu1 + 1; kk <= i__2; ++kk) {
			j = indr[kk];
			k = hk[j];
			v = (d__1 = a[locc[j]], abs(d__1));
/* Biggest aij in column j */
			hchange_(&ha[1], &hj[1], &hk[1], &hlen, n, &k, &v, &j,
				 &h__);
			hops += h__;
		    }
		}
	    }
	}
/*        =============================================================== */
/*        Negate lengths of pivot row and column so they will be */
/*        eliminated during compressions. */
/*        =============================================================== */
	lenr[ibest] = -ncold;
	lenc[jbest] = -nrowd;
/*        Test for fatal bug: row or column lists overwriting L and U. */
	if (lrow > lsave) {
	    goto L980;
	}
	if (lcol > lsave) {
	    goto L980;
	}
/*        Reset the file lengths if pivot row or col was at the end. */
	if (ibest == ilast) {
	    lrow = locr[ibest];
	}
	if (jbest == jlast) {
	    lcol = locc[jbest];
	}
/* L800: */
    }
/*     ------------------------------------------------------------------ */
/*     End of main loop. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Normal exit. */
/*     Move empty rows and cols to the end of ip, iq. */
/*     Then finish with a dense LU if necessary. */
/*     ------------------------------------------------------------------ */
L900:
    *inform__ = 0;
    lu1pq3_(m, &lenr[1], &ip[1], &ipinv[1], &mrank);
    lu1pq3_(n, &lenc[1], &iq[1], &iqinv[1], nrank);
    *nrank = min(mrank,*nrank);
    if (denslu) {
/*        call timer ( 'start ', 17 ) */
	lu1ful_(m, n, lena, &lend, &lu1, &tpp, &mleft, &nleft, nrank, &nrowu, 
		lenl, lenu, &nsing, &keeplu, &small, &a[1], &a[ld], &indc[1], 
		&indr[1], &ip[1], &iq[1], &lenc[1], &lenr[1], &locc[1], &
		ipinv[1], &locr[1]);
/* ***     21 Dec 1994: Bug in next line. */
/* ***     nrank  = nrank - nsing.  Changed to next line: */
/* ***     nrank  = minmn - nsing */
/* ***     26 Mar 2006: Previous line caused bug with m<n and nsing>0. */
/*        Don't mess with nrank any more.  Let end of lu1fac handle it. */
/*        call timer ( 'finish', 17 ) */
    }
    *minlen = *lenl + *lenu + (*m + *n << 1);
    goto L990;
/*     Not enough space free after a compress. */
/*     Set  minlen  to an estimate of the necessary value of  lena. */
L970:
    *inform__ = 7;
    *minlen = *lena + lfile + (*m + *n << 1);
    goto L990;
/*     Fatal error.  This will never happen! */
/*    (Famous last words.) */
L980:
    *inform__ = 8;
    goto L990;
/*     Fatal error with TSP.  Diagonal pivot not found. */
L985:
    *inform__ = 9;
/*     Exit. */
L990:
/*     call timer ( 'finish', 3 ) */
    return 0;
} /* lu1fad_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1fad */
/* Subroutine */ int lu1gau_(integer *m, integer *melim, integer *ncold, 
	integer *nspare, doublereal *small, integer *lpivc1, integer *lpivc2, 
	integer *lfirst, integer *lpivr2, integer *lfree, integer *minfre, 
	integer *ilast, integer *jlast, integer *lrow, integer *lcol, integer 
	*lu, integer *nfill, doublereal *a, integer *indc, integer *indr, 
	integer *lenc, integer *lenr, integer *locc, integer *locr, integer *
	mark, doublereal *al, integer *markl, doublereal *au, integer *ifill, 
	integer *jfill)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, l1, l2, lc, ll, lr;
    static doublereal uj;
    static integer lc1, lc2, lr1;
    static doublereal aij;
    static integer leni, lenj, lrep, last;
    static logical atend;
    static integer nfree, ndone, ndrop;

/*     ------------------------------------------------------------------ */
/*     lu1gau does most of the work for each step of */
/*     Gaussian elimination. */
/*     A multiple of the pivot column is added to each other column j */
/*     in the pivot row.  The column list is fully updated. */
/*     The row list is updated if there is room, but some fill-ins may */
/*     remain, as indicated by ifill and jfill. */


/*  Input: */
/*     ilast    is the row    at the end of the row    list. */
/*     jlast    is the column at the end of the column list. */
/*     lfirst   is the first column to be processed. */
/*     lu + 1   is the corresponding element of U in au(*). */
/*     nfill    keeps track of pending fill-in. */
/*     a(*)     contains the nonzeros for each column j. */
/*     indc(*)  contains the row indices for each column j. */
/*     al(*)    contains the new column of L.  A multiple of it is */
/*              used to modify each column. */
/*     mark(*)  has been set to -1, -2, -3, ... in the rows */
/*              corresponding to nonzero 1, 2, 3, ... of the col of L. */
/*     au(*)    contains the new row of U.  Each nonzero gives the */
/*              required multiple of the column of L. */

/*  Workspace: */
/*     markl(*) marks the nonzeros of L actually used. */
/*              (A different mark, namely j, is used for each column.) */

/*  Output: */
/*     ilast     New last row    in the row    list. */
/*     jlast     New last column in the column list. */
/*     lfirst    = 0 if all columns were completed, */
/*               > 0 otherwise. */
/*     lu        returns the position of the last nonzero of U */
/*               actually used, in case we come back in again. */
/*     nfill     keeps track of the total extra space needed in the */
/*               row file. */
/*     ifill(ll) counts pending fill-in for rows involved in the new */
/*               column of L. */
/*     jfill(lu) marks the first pending fill-in stored in columns */
/*               involved in the new row of U. */

/*     16 Apr 1989: First version of lu1gau. */
/*     23 Apr 1989: lfirst, lu, nfill are now input and output */
/*                  to allow re-entry if elimination is interrupted. */
/*     23 Mar 2001: Introduced ilast, jlast. */
/*     27 Mar 2001: Allow fill-in "in situ" if there is already room */
/*                  up to but NOT INCLUDING the end of the */
/*                  row or column file. */
/*                  Seems safe way to avoid overwriting empty rows/cols */
/*                  at the end.  (May not be needed though, now that we */
/*                  have ilast and jlast.) */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ifill;
    --markl;
    --al;
    --jfill;
    --au;
    --a;
    --indc;
    --indr;
    --lenc;
    --lenr;
    --locc;
    --locr;
    --mark;

    /* Function Body */
    i__1 = *lpivr2;
    for (lr = *lfirst; lr <= i__1; ++lr) {
	j = indr[lr];
	lenj = lenc[j];
	nfree = *lfree - *lcol;
	if (nfree < *minfre) {
	    goto L900;
	}
/*        --------------------------------------------------------------- */
/*        Inner loop to modify existing nonzeros in column  j. */
/*        Loop 440 performs most of the arithmetic involved in the */
/*        whole LU factorization. */
/*        ndone counts how many multipliers were used. */
/*        ndrop counts how many modified nonzeros are negligibly small. */
/*        --------------------------------------------------------------- */
	++(*lu);
	uj = au[*lu];
	lc1 = locc[j];
	lc2 = lc1 + lenj - 1;
	atend = j == *jlast;
	ndone = 0;
	if (lenj == 0) {
	    goto L500;
	}
	ndrop = 0;
	i__2 = lc2;
	for (l = lc1; l <= i__2; ++l) {
	    i__ = indc[l];
	    ll = -mark[i__];
	    if (ll > 0) {
		++ndone;
		markl[ll] = j;
		a[l] += al[ll] * uj;
		if ((d__1 = a[l], abs(d__1)) <= *small) {
		    ++ndrop;
		}
	    }
/* L440: */
	}
/*        --------------------------------------------------------------- */
/*        Remove any negligible modified nonzeros from both */
/*        the column file and the row file. */
/*        --------------------------------------------------------------- */
	if (ndrop == 0) {
	    goto L500;
	}
	k = lc1;
	i__2 = lc2;
	for (l = lc1; l <= i__2; ++l) {
	    i__ = indc[l];
	    if ((d__1 = a[l], abs(d__1)) <= *small) {
		goto L460;
	    }
	    a[k] = a[l];
	    indc[k] = i__;
	    ++k;
	    goto L480;
/*           Delete the nonzero from the row file. */
L460:
	    --lenj;
	    --lenr[i__];
	    lr1 = locr[i__];
	    last = lr1 + lenr[i__];
	    i__3 = last;
	    for (lrep = lr1; lrep <= i__3; ++lrep) {
		if (indr[lrep] == j) {
		    goto L475;
		}
/* L470: */
	    }
L475:
	    indr[lrep] = indr[last];
	    indr[last] = 0;
	    if (i__ == *ilast) {
		--(*lrow);
	    }
L480:
	    ;
	}
/*        Free the deleted elements from the column file. */
	i__2 = lc2;
	for (l = k; l <= i__2; ++l) {
	    indc[l] = 0;
/* L490: */
	}
	if (atend) {
	    *lcol = k - 1;
	}
/*        --------------------------------------------------------------- */
/*        Deal with the fill-in in column j. */
/*        --------------------------------------------------------------- */
L500:
	if (ndone == *melim) {
	    goto L590;
	}
/*        See if column j already has room for the fill-in. */
	if (atend) {
	    goto L540;
	}
	last = lc1 + lenj - 1;
	l1 = last + 1;
	l2 = last + (*melim - ndone);
/* 27 Mar 2001: Be sure it's not at or past end of the col file. */
	if (l2 >= *lcol) {
	    goto L520;
	}
	i__2 = l2;
	for (l = l1; l <= i__2; ++l) {
	    if (indc[l] != 0) {
		goto L520;
	    }
/* L510: */
	}
	goto L540;
/*        We must move column j to the end of the column file. */
/*        First, leave some spare room at the end of the */
/*        current last column. */
L520:
	i__2 = *lcol + *nspare;
	for (l = *lcol + 1; l <= i__2; ++l) {
	    *lcol = l;
	    indc[l] = 0;
/* Spare space is free. */
/* L522: */
	}
	atend = TRUE_;
	*jlast = j;
	l1 = lc1;
	lc1 = *lcol + 1;
	locc[j] = lc1;
	i__2 = last;
	for (l = l1; l <= i__2; ++l) {
	    ++(*lcol);
	    a[*lcol] = a[l];
	    indc[*lcol] = indc[l];
	    indc[l] = 0;
/* Free space. */
/* L525: */
	}
/*        --------------------------------------------------------------- */
/*        Inner loop for the fill-in in column j. */
/*        This is usually not very expensive. */
/*        --------------------------------------------------------------- */
L540:
	last = lc1 + lenj - 1;
	ll = 0;
	i__2 = *lpivc2;
	for (lc = *lpivc1; lc <= i__2; ++lc) {
	    ++ll;
	    if (markl[ll] == j) {
		goto L560;
	    }
	    aij = al[ll] * uj;
	    if (abs(aij) <= *small) {
		goto L560;
	    }
	    ++lenj;
	    ++last;
	    a[last] = aij;
	    i__ = indc[lc];
	    indc[last] = i__;
	    leni = lenr[i__];
/*           Add 1 fill-in to row i if there is already room. */
/*           27 Mar 2001: Be sure it's not at or past the end */
/*                        of the row file. */
	    l = locr[i__] + leni;
	    if (l >= *lrow) {
		goto L550;
	    }
	    if (indr[l] > 0) {
		goto L550;
	    }
	    indr[l] = j;
	    lenr[i__] = leni + 1;
	    goto L560;
/*           Row i does not have room for the fill-in. */
/*           Increment ifill(ll) to count how often this has */
/*           happened to row i.  Also, add m to the row index */
/*           indc(last) in column j to mark it as a fill-in that is */
/*           still pending. */

/*           If this is the first pending fill-in for row i, */
/*           nfill includes the current length of row i */
/*           (since the whole row has to be moved later). */

/*           If this is the first pending fill-in for column j, */
/*           jfill(lu) records the current length of column j */
/*           (to shorten the search for pending fill-ins later). */
L550:
	    if (ifill[ll] == 0) {
		*nfill = *nfill + leni + *nspare;
	    }
	    if (jfill[*lu] == 0) {
		jfill[*lu] = lenj;
	    }
	    ++(*nfill);
	    ++ifill[ll];
	    indc[last] = *m + i__;
L560:
	    ;
	}
	if (atend) {
	    *lcol = last;
	}
/*        End loop for column  j.  Store its final length. */
L590:
	lenc[j] = lenj;
/* L600: */
    }
/*     Successful completion. */
    *lfirst = 0;
    return 0;
/*     Interruption.  We have to come back in after the */
/*     column file is compressed.  Give lfirst a new value. */
/*     lu and nfill will retain their current values. */
L900:
    *lfirst = lr;
    return 0;
} /* lu1gau_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1gau */
/* Subroutine */ int lu1mar_(integer *m, integer *n, integer *lena, integer *
	maxmn, logical *tcp, doublereal *aijtol, doublereal *ltol, integer *
	maxcol, integer *maxrow, integer *ibest, integer *jbest, integer *
	mbest, doublereal *a, integer *indc, integer *indr, integer *ip, 
	integer *iq, integer *lenc, integer *lenr, integer *locc, integer *
	locr, integer *iploc, integer *iqloc)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, lc, lp, lq, lr, nz, lc1, lc2, lp1, lq1, lq2, lp2, 
	    lr1, lr2, nz1;
    static doublereal aij;
    static integer len1;
    static doublereal amax, cmax;
    static integer ncol, nrow;
    static doublereal abest, lbest;
    static integer kbest, merit;

/*     ------------------------------------------------------------------ */
/*     lu1mar  uses a Markowitz criterion to select a pivot element */
/*     for the next stage of a sparse LU factorization, */
/*     subject to a Threshold Partial Pivoting stability criterion (TPP) */
/*     that bounds the elements of L. */

/*     00 Jan 1986  Version documented in LUSOL paper: */
/*                  Gill, Murray, Saunders and Wright (1987), */
/*                  Maintaining LU factors of a general sparse matrix, */
/*                  Linear algebra and its applications 88/89, 239-270. */

/*     02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest */
/*                  element in each column is now kept at the start of */
/*                  the column, i.e. in position locc(j) of a and indc. */
/*                  This should speed up the Markowitz searches. */

/*     26 Apr 1989  Both columns and rows searched during spars1 phase. */
/*                  Only columns searched during spars2 phase. */
/*                  maxtie replaced by maxcol and maxrow. */
/*     05 Nov 1993  Initializing  "mbest = m * n"  wasn't big enough when */
/*                  m = 10, n = 3, and last column had 7 nonzeros. */
/*     09 Feb 1994  Realised that "mbest = maxmn * maxmn" might overflow. */
/*                  Changed to    "mbest = maxmn * 1000". */
/*     27 Apr 2000  On large example from Todd Munson, */
/*                  that allowed  "if (mbest .le. nz1**2) go to 900" */
/*                  to exit before any pivot had been found. */
/*                  Introduced kbest = mbest / nz1. */
/*                  Most pivots can be rejected with no integer multiply. */
/*                  True merit is evaluated only if it's as good as the */
/*                  best so far (or better).  There should be no danger */
/*                  of integer overflow unless A is incredibly */
/*                  large and dense. */

/*     10 Sep 2000  TCP, aijtol added for Threshold Complete Pivoting. */

/*     Systems Optimization Laboratory, Stanford University. */
/*     ------------------------------------------------------------------ */
/*     gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper. */
/*     ------------------------------------------------------------------ */
/*     Search cols of length nz = 1, then rows of length nz = 1, */
/*     then   cols of length nz = 2, then rows of length nz = 2, etc. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iqloc;
    --locr;
    --lenr;
    --ip;
    --iploc;
    --locc;
    --lenc;
    --iq;
    --indr;
    --indc;
    --a;

    /* Function Body */
    abest = 0.;
    lbest = 0.;
    *ibest = 0;
    kbest = *maxmn + 1;
    *mbest = -1;
    ncol = 0;
    nrow = 0;
    nz1 = 0;
    i__1 = *maxmn;
    for (nz = 1; nz <= i__1; ++nz) {
/*        nz1    = nz - 1 */
/*        if (mbest .le. nz1**2) go to 900 */
	if (kbest <= nz1) {
	    goto L900;
	}
	if (*ibest > 0) {
	    if (ncol >= *maxcol) {
		goto L200;
	    }
	}
	if (nz > *m) {
	    goto L200;
	}
/*        --------------------------------------------------------------- */
/*        Search the set of columns of length  nz. */
/*        --------------------------------------------------------------- */
	lq1 = iqloc[nz];
	lq2 = *n;
	if (nz < *m) {
	    lq2 = iqloc[nz + 1] - 1;
	}
	i__2 = lq2;
	for (lq = lq1; lq <= i__2; ++lq) {
	    ++ncol;
	    j = iq[lq];
	    lc1 = locc[j];
	    lc2 = lc1 + nz1;
	    amax = (d__1 = a[lc1], abs(d__1));
/*           Test all aijs in this column. */
/*           amax is the largest element (the first in the column). */
/*           cmax is the largest multiplier if aij becomes pivot. */
	    if (*tcp) {
		if (amax < *aijtol) {
		    goto L180;
		}
/* Nothing in whole column */
	    }
	    i__3 = lc2;
	    for (lc = lc1; lc <= i__3; ++lc) {
		i__ = indc[lc];
		len1 = lenr[i__] - 1;
/*              merit  = nz1 * len1 */
/*              if (merit .gt. mbest) go to 160 */
		if (len1 > kbest) {
		    goto L160;
		}
/*              aij  has a promising merit. */
/*              Apply the stability test. */
/*              We require  aij  to be sufficiently large compared to */
/*              all other nonzeros in column  j.  This is equivalent */
/*              to requiring cmax to be bounded by Ltol. */
		if (lc == lc1) {
/*                 This is the maximum element, amax. */
/*                 Find the biggest element in the rest of the column */
/*                 and hence get cmax.  We know cmax .le. 1, but */
/*                 we still want it exactly in order to break ties. */
/*                 27 Apr 2002: Settle for cmax = 1. */
		    aij = amax;
		    cmax = 1.;
/*                 cmax   = zero */
/*                 do 140 l = lc1 + 1, lc2 */
/*                    cmax  = max( cmax, abs( a(l) ) ) */
/*  140            continue */
/*                 cmax   = cmax / amax */
		} else {
/*                 aij is not the biggest element, so cmax .ge. 1. */
/*                 Bail out if cmax will be too big. */
		    aij = (d__1 = a[lc], abs(d__1));
		    if (*tcp) {
/* Absolute test for Complete Pivoting */
			if (aij < *aijtol) {
			    goto L160;
			}
		    } else {
/* !! TPP */
			if (aij * *ltol < amax) {
			    goto L160;
			}
		    }
		    cmax = amax / aij;
		}
/*              aij  is big enough.  Its maximum multiplier is cmax. */
		merit = nz1 * len1;
		if (merit == *mbest) {
/*                 Break ties. */
/*                 (Initializing mbest < 0 prevents getting here if */
/*                 nothing has been found yet.) */
/*                 In this version we minimize cmax */
/*                 but if it is already small we maximize the pivot. */
		    if (lbest <= 2. && cmax <= 2.) {
			if (abest >= aij) {
			    goto L160;
			}
		    } else {
			if (lbest <= cmax) {
			    goto L160;
			}
		    }
		}
/*              aij  is the best pivot so far. */
		*ibest = i__;
		*jbest = j;
		kbest = len1;
		*mbest = merit;
		abest = aij;
		lbest = cmax;
		if (nz == 1) {
		    goto L900;
		}
L160:
		;
	    }
/*           Finished with that column. */
	    if (*ibest > 0) {
		if (ncol >= *maxcol) {
		    goto L200;
		}
	    }
L180:
	    ;
	}
/*        --------------------------------------------------------------- */
/*        Search the set of rows of length  nz. */
/*        --------------------------------------------------------------- */
/* 200    if (mbest .le. nz*nz1) go to 900 */
L200:
	if (kbest <= nz) {
	    goto L900;
	}
	if (*ibest > 0) {
	    if (nrow >= *maxrow) {
		goto L290;
	    }
	}
	if (nz > *n) {
	    goto L290;
	}
	lp1 = iploc[nz];
	lp2 = *m;
	if (nz < *n) {
	    lp2 = iploc[nz + 1] - 1;
	}
	i__2 = lp2;
	for (lp = lp1; lp <= i__2; ++lp) {
	    ++nrow;
	    i__ = ip[lp];
	    lr1 = locr[i__];
	    lr2 = lr1 + nz1;
	    i__3 = lr2;
	    for (lr = lr1; lr <= i__3; ++lr) {
		j = indr[lr];
		len1 = lenc[j] - 1;
/*              merit  = nz1 * len1 */
/*              if (merit .gt. mbest) go to 260 */
		if (len1 > kbest) {
		    goto L260;
		}
/*              aij  has a promising merit. */
/*              Find where  aij  is in column  j. */
		lc1 = locc[j];
		lc2 = lc1 + len1;
		amax = (d__1 = a[lc1], abs(d__1));
		i__4 = lc2;
		for (lc = lc1; lc <= i__4; ++lc) {
		    if (indc[lc] == i__) {
			goto L230;
		    }
/* L220: */
		}
/*              Apply the same stability test as above. */
L230:
		aij = (d__1 = a[lc], abs(d__1));
		if (*tcp) {
/* !! Absolute test for Complete Pivoting */
		    if (aij < *aijtol) {
			goto L260;
		    }
		}
		if (lc == lc1) {
/*                 This is the maximum element, amax. */
/*                 Find the biggest element in the rest of the column */
/*                 and hence get cmax.  We know cmax .le. 1, but */
/*                 we still want it exactly in order to break ties. */
/*                 27 Apr 2002: Settle for cmax = 1. */
		    cmax = 1.;
/*                 cmax   = zero */
/*                 do 240 l = lc1 + 1, lc2 */
/*                    cmax  = max( cmax, abs( a(l) ) ) */
/* 240             continue */
/*                 cmax   = cmax / amax */
		} else {
/*                 aij is not the biggest element, so cmax .ge. 1. */
/*                 Bail out if cmax will be too big. */
		    if (*tcp) {
/* relax */
		    } else {
			if (aij * *ltol < amax) {
			    goto L260;
			}
		    }
		    cmax = amax / aij;
		}
/*              aij  is big enough.  Its maximum multiplier is cmax. */
		merit = nz1 * len1;
		if (merit == *mbest) {
/*                 Break ties as before. */
/*                 (Initializing mbest < 0 prevents getting here if */
/*                 nothing has been found yet.) */
		    if (lbest <= 2. && cmax <= 2.) {
			if (abest >= aij) {
			    goto L260;
			}
		    } else {
			if (lbest <= cmax) {
			    goto L260;
			}
		    }
		}
/*              aij  is the best pivot so far. */
		*ibest = i__;
		*jbest = j;
		kbest = len1;
		*mbest = merit;
		abest = aij;
		lbest = cmax;
		if (nz == 1) {
		    goto L900;
		}
L260:
		;
	    }
/*           Finished with that row. */
	    if (*ibest > 0) {
		if (nrow >= *maxrow) {
		    goto L290;
		}
	    }
/* L280: */
	}
/*        See if it's time to quit. */
L290:
	if (*ibest > 0) {
	    if (nrow >= *maxrow && ncol >= *maxcol) {
		goto L900;
	    }
	}
/*        Press on with next nz. */
	nz1 = nz;
	if (*ibest > 0) {
	    kbest = *mbest / nz1;
	}
    }
L900:
    return 0;
} /* lu1mar_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1mar */
/* Subroutine */ int lu1mrp_(integer *m, integer *n, integer *lena, integer *
	maxmn, doublereal *ltol, integer *maxcol, integer *maxrow, integer *
	ibest, integer *jbest, integer *mbest, doublereal *a, integer *indc, 
	integer *indr, integer *ip, integer *iq, integer *lenc, integer *lenr,
	 integer *locc, integer *locr, integer *iploc, integer *iqloc, 
	doublereal *amaxr)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, lc, lp, lq, lr, nz, lc1, lc2, lp1, lp2, lq1, lq2, 
	    lr1, lr2, nz1;
    static doublereal aij;
    static integer len1;
    static doublereal amax;
    static integer ncol, nrow;
    static doublereal abest;
    static integer kbest;
    static doublereal atoli, atolj;
    static integer merit;

/*     ------------------------------------------------------------------ */
/*     lu1mRP  uses a Markowitz criterion to select a pivot element */
/*     for the next stage of a sparse LU factorization, */
/*     subject to a Threshold Rook Pivoting stability criterion (TRP) */
/*     that bounds the elements of L and U. */

/*     11 Jun 2002: First version of lu1mRP derived from lu1mar. */
/*     11 Jun 2002: Current version of lu1mRP. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Search cols of length nz = 1, then rows of length nz = 1, */
/*     then   cols of length nz = 2, then rows of length nz = 2, etc. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --amaxr;
    --iqloc;
    --locr;
    --lenr;
    --ip;
    --iploc;
    --locc;
    --lenc;
    --iq;
    --indr;
    --indc;
    --a;

    /* Function Body */
    abest = 0.;
    *ibest = 0;
    kbest = *maxmn + 1;
    *mbest = -1;
    ncol = 0;
    nrow = 0;
    nz1 = 0;
    i__1 = *maxmn;
    for (nz = 1; nz <= i__1; ++nz) {
/*        nz1    = nz - 1 */
/*        if (mbest .le. nz1**2) go to 900 */
	if (kbest <= nz1) {
	    goto L900;
	}
	if (*ibest > 0) {
	    if (ncol >= *maxcol) {
		goto L200;
	    }
	}
	if (nz > *m) {
	    goto L200;
	}
/*        --------------------------------------------------------------- */
/*        Search the set of columns of length  nz. */
/*        --------------------------------------------------------------- */
	lq1 = iqloc[nz];
	lq2 = *n;
	if (nz < *m) {
	    lq2 = iqloc[nz + 1] - 1;
	}
	i__2 = lq2;
	for (lq = lq1; lq <= i__2; ++lq) {
	    ++ncol;
	    j = iq[lq];
	    lc1 = locc[j];
	    lc2 = lc1 + nz1;
	    amax = (d__1 = a[lc1], abs(d__1));
	    atolj = amax / *ltol;
/*           Test all aijs in this column. */
/* Min size of pivots in col j */
	    i__3 = lc2;
	    for (lc = lc1; lc <= i__3; ++lc) {
		i__ = indc[lc];
		len1 = lenr[i__] - 1;
/*              merit  = nz1 * len1 */
/*              if (merit .gt. mbest) go to 160 */
		if (len1 > kbest) {
		    goto L160;
		}
/*              aij  has a promising merit. */
/*              Apply the Threshold Rook Pivoting stability test. */
/*              First we require aij to be sufficiently large */
/*              compared to other nonzeros in column j. */
/*              Then  we require aij to be sufficiently large */
/*              compared to other nonzeros in row    i. */
		aij = (d__1 = a[lc], abs(d__1));
		if (aij < atolj) {
		    goto L160;
		}
		if (aij * *ltol < amaxr[i__]) {
		    goto L160;
		}
/*              aij  is big enough. */
		merit = nz1 * len1;
		if (merit == *mbest) {
/*                 Break ties. */
/*                 (Initializing mbest < 0 prevents getting here if */
/*                 nothing has been found yet.) */
		    if (abest >= aij) {
			goto L160;
		    }
		}
/*              aij  is the best pivot so far. */
		*ibest = i__;
		*jbest = j;
		kbest = len1;
		*mbest = merit;
		abest = aij;
		if (nz == 1) {
		    goto L900;
		}
L160:
		;
	    }
/*           Finished with that column. */
	    if (*ibest > 0) {
		if (ncol >= *maxcol) {
		    goto L200;
		}
	    }
/* L180: */
	}
/*        --------------------------------------------------------------- */
/*        Search the set of rows of length  nz. */
/*        --------------------------------------------------------------- */
/* 200    if (mbest .le. nz*nz1) go to 900 */
L200:
	if (kbest <= nz) {
	    goto L900;
	}
	if (*ibest > 0) {
	    if (nrow >= *maxrow) {
		goto L290;
	    }
	}
	if (nz > *n) {
	    goto L290;
	}
	lp1 = iploc[nz];
	lp2 = *m;
	if (nz < *n) {
	    lp2 = iploc[nz + 1] - 1;
	}
	i__2 = lp2;
	for (lp = lp1; lp <= i__2; ++lp) {
	    ++nrow;
	    i__ = ip[lp];
	    lr1 = locr[i__];
	    lr2 = lr1 + nz1;
	    atoli = amaxr[i__] / *ltol;
/* Min size of pivots in row i */
	    i__3 = lr2;
	    for (lr = lr1; lr <= i__3; ++lr) {
		j = indr[lr];
		len1 = lenc[j] - 1;
/*              merit  = nz1 * len1 */
/*              if (merit .gt. mbest) go to 260 */
		if (len1 > kbest) {
		    goto L260;
		}
/*              aij  has a promising merit. */
/*              Find where  aij  is in column j. */
		lc1 = locc[j];
		lc2 = lc1 + len1;
		amax = (d__1 = a[lc1], abs(d__1));
		i__4 = lc2;
		for (lc = lc1; lc <= i__4; ++lc) {
		    if (indc[lc] == i__) {
			goto L230;
		    }
		}
/*              Apply the Threshold Rook Pivoting stability test. */
/*              First we require aij to be sufficiently large */
/*              compared to other nonzeros in row    i. */
/*              Then  we require aij to be sufficiently large */
/*              compared to other nonzeros in column j. */
L230:
		aij = (d__1 = a[lc], abs(d__1));
		if (aij < atoli) {
		    goto L260;
		}
		if (aij * *ltol < amax) {
		    goto L260;
		}
/*              aij  is big enough. */
		merit = nz1 * len1;
		if (merit == *mbest) {
/*                 Break ties as before. */
/*                 (Initializing mbest < 0 prevents getting here if */
/*                 nothing has been found yet.) */
		    if (abest >= aij) {
			goto L260;
		    }
		}
/*              aij  is the best pivot so far. */
		*ibest = i__;
		*jbest = j;
		kbest = len1;
		*mbest = merit;
		abest = aij;
		if (nz == 1) {
		    goto L900;
		}
L260:
		;
	    }
/*           Finished with that row. */
	    if (*ibest > 0) {
		if (nrow >= *maxrow) {
		    goto L290;
		}
	    }
/* L280: */
	}
/*        See if it's time to quit. */
L290:
	if (*ibest > 0) {
	    if (nrow >= *maxrow && ncol >= *maxcol) {
		goto L900;
	    }
	}
/*        Press on with next nz. */
	nz1 = nz;
	if (*ibest > 0) {
	    kbest = *mbest / nz1;
	}
    }
L900:
    return 0;
} /* lu1mrp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1mRP */
/* Subroutine */ int lu1mcp_(integer *m, integer *n, integer *lena, 
	doublereal *aijtol, integer *ibest, integer *jbest, integer *mbest, 
	doublereal *a, integer *indc, integer *indr, integer *lenc, integer *
	lenr, integer *locc, integer *hlen, doublereal *ha, integer *hj)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, lc, lc1, lc2, nz1;
    static doublereal aij;
    static integer len1;
    static doublereal amax;
    static integer lenj;
    static doublereal cmax;
    static integer ncol, kheap;
    static doublereal abest, lbest;
    static integer merit, maxcol;

/*     ------------------------------------------------------------------ */
/*     lu1mCP  uses a Markowitz criterion to select a pivot element */
/*     for the next stage of a sparse LU factorization, */
/*     subject to a Threshold Complete Pivoting stability criterion (TCP) */
/*     that bounds the elements of L and U. */

/*     09 May 2002: First version of lu1mCP. */
/*                  It searches columns only, using the heap that */
/*                  holds the largest element in each column. */
/*     09 May 2002: Current version of lu1mCP. */
/*     ------------------------------------------------------------------ */
/*     gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper. */
/*     ------------------------------------------------------------------ */
/*     Search up to maxcol columns stored at the top of the heap. */
/*     The very top column helps initialize mbest. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --lenr;
    --locc;
    --lenc;
    --indr;
    --indc;
    --a;
    --hj;
    --ha;

    /* Function Body */
    abest = 0.;
    lbest = 0.;
    *ibest = 0;
    *jbest = hj[1];
/* Column at the top of the heap */
    lenj = lenc[*jbest];
    *mbest = lenj * *hlen;
/* Bigger than any possible merit */
    maxcol = 40;
/* ??? Big question */
    ncol = 0;
/* No. of columns searched */
    i__1 = *hlen;
    for (kheap = 1; kheap <= i__1; ++kheap) {
	amax = ha[kheap];
	if (amax < *aijtol) {
	    goto L500;
	}
	++ncol;
	j = hj[kheap];
/*        --------------------------------------------------------------- */
/*        This column has at least one entry big enough (the top one). */
/*        Search the column for other possibilities. */
/*        --------------------------------------------------------------- */
	lenj = lenc[j];
	nz1 = lenj - 1;
	lc1 = locc[j];
	lc2 = lc1 + nz1;
/* --      amax   = abs( a(lc1) ) */
/*        Test all aijs in this column. */
/*        amax is the largest element (the first in the column). */
/*        cmax is the largest multiplier if aij becomes pivot. */
	i__2 = lc2;
	for (lc = lc1; lc <= i__2; ++lc) {
	    i__ = indc[lc];
	    len1 = lenr[i__] - 1;
	    merit = nz1 * len1;
	    if (merit > *mbest) {
		goto L160;
	    }
/*           aij  has a promising merit. */
	    if (lc == lc1) {
/*              This is the maximum element, amax. */
/*              Find the biggest element in the rest of the column */
/*              and hence get cmax.  We know cmax .le. 1, but */
/*              we still want it exactly in order to break ties. */
/*              27 Apr 2002: Settle for cmax = 1. */
		aij = amax;
		cmax = 1.;
/*              cmax   = zero */
/*              do 140 l = lc1 + 1, lc2 */
/*                 cmax  = max( cmax, abs( a(l) ) ) */
/*  140         continue */
/*              cmax   = cmax / amax */
	    } else {
/*              aij is not the biggest element, so cmax .ge. 1. */
/*              Bail out if cmax will be too big. */
		aij = (d__1 = a[lc], abs(d__1));
		if (aij < *aijtol) {
		    goto L160;
		}
		cmax = amax / aij;
	    }
/*           aij  is big enough.  Its maximum multiplier is cmax. */
	    if (merit == *mbest) {
/*              Break ties. */
/*              (Initializing mbest "too big" prevents getting here if */
/*              nothing has been found yet.) */
/*              In this version we minimize cmax */
/*              but if it is already small we maximize the pivot. */
		if (lbest <= 2. && cmax <= 2.) {
		    if (abest >= aij) {
			goto L160;
		    }
		} else {
		    if (lbest <= cmax) {
			goto L160;
		    }
		}
	    }
/*           aij  is the best pivot so far. */
	    *ibest = i__;
	    *jbest = j;
	    *mbest = merit;
	    abest = aij;
	    lbest = cmax;
	    if (merit == 0) {
		goto L900;
	    }
/* Col or row of length 1 */
L160:
	    ;
	}
	if (ncol >= maxcol) {
	    goto L900;
	}
L500:
	;
    }
L900:
    return 0;
} /* lu1mcp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1mCP */
/* Subroutine */ int lu1msp_(integer *m, integer *n, integer *lena, integer *
	maxmn, doublereal *ltol, integer *maxcol, integer *ibest, integer *
	jbest, integer *mbest, doublereal *a, integer *indc, integer *iq, 
	integer *locc, integer *iqloc)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, lc, lq, nz, lc1, lc2, lq1, lq2, nz1;
    static doublereal aij, amax;
    static integer ncol;
    static doublereal abest;
    static integer kbest;
    static doublereal atolj;
    static integer merit;

/*     ------------------------------------------------------------------ */
/*     lu1mSP  is intended for symmetric matrices that are either */
/*     definite or quasi-definite. */
/*     lu1mSP  uses a Markowitz criterion to select a pivot element for */
/*     the next stage of a sparse LU factorization of a symmetric matrix, */
/*     subject to a Threshold Symmetric Pivoting stability criterion */
/*     (TSP) restricted to diagonal elements to preserve symmetry. */
/*     This bounds the elements of L and U and should have rank-revealing */
/*     properties analogous to Threshold Rook Pivoting for unsymmetric */
/*     matrices. */

/*     14 Dec 2002: First version of lu1mSP derived from lu1mRP. */
/*                  There is no safeguard to ensure that A is symmetric. */
/*     14 Dec 2002: Current version of lu1mSP. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     Search cols of length nz = 1, then cols of length nz = 2, etc. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iqloc;
    --locc;
    --iq;
    --indc;
    --a;

    /* Function Body */
    abest = 0.;
    *ibest = 0;
    kbest = *maxmn + 1;
    *mbest = -1;
    ncol = 0;
    nz1 = 0;
    i__1 = *maxmn;
    for (nz = 1; nz <= i__1; ++nz) {
/*        nz1    = nz - 1 */
/*        if (mbest .le. nz1**2) go to 900 */
	if (kbest <= nz1) {
	    goto L900;
	}
	if (*ibest > 0) {
	    if (ncol >= *maxcol) {
		goto L200;
	    }
	}
	if (nz > *m) {
	    goto L200;
	}
/*        --------------------------------------------------------------- */
/*        Search the set of columns of length  nz. */
/*        --------------------------------------------------------------- */
	lq1 = iqloc[nz];
	lq2 = *n;
	if (nz < *m) {
	    lq2 = iqloc[nz + 1] - 1;
	}
	i__2 = lq2;
	for (lq = lq1; lq <= i__2; ++lq) {
	    ++ncol;
	    j = iq[lq];
	    lc1 = locc[j];
	    lc2 = lc1 + nz1;
	    amax = (d__1 = a[lc1], abs(d__1));
	    atolj = amax / *ltol;
/*           Test all aijs in this column. */
/*           Ignore everything except the diagonal. */
/* Min size of pivots in col j */
	    i__3 = lc2;
	    for (lc = lc1; lc <= i__3; ++lc) {
		i__ = indc[lc];
		if (i__ != j) {
		    goto L160;
		}
/*              merit  = nz1 * nz1 */
/*              if (merit .gt. mbest) go to 160 */
/* Skip off-diagonals. */
		if (nz1 > kbest) {
		    goto L160;
		}
/*              aij  has a promising merit. */
/*              Apply the Threshold Partial Pivoting stability test */
/*              (which is equivalent to Threshold Rook Pivoting for */
/*              symmetric matrices). */
/*              We require aij to be sufficiently large */
/*              compared to other nonzeros in column j. */
		aij = (d__1 = a[lc], abs(d__1));
		if (aij < atolj) {
		    goto L160;
		}
/*              aij  is big enough. */
		merit = nz1 * nz1;
		if (merit == *mbest) {
/*                 Break ties. */
/*                 (Initializing mbest < 0 prevents getting here if */
/*                 nothing has been found yet.) */
		    if (abest >= aij) {
			goto L160;
		    }
		}
/*              aij  is the best pivot so far. */
		*ibest = i__;
		*jbest = j;
		kbest = nz1;
		*mbest = merit;
		abest = aij;
		if (nz == 1) {
		    goto L900;
		}
L160:
		;
	    }
/*           Finished with that column. */
	    if (*ibest > 0) {
		if (ncol >= *maxcol) {
		    goto L200;
		}
	    }
/* L180: */
	}
/*        See if it's time to quit. */
L200:
	if (*ibest > 0) {
	    if (ncol >= *maxcol) {
		goto L900;
	    }
	}
/*        Press on with next nz. */
	nz1 = nz;
	if (*ibest > 0) {
	    kbest = *mbest / nz1;
	}
    }
L900:
    return 0;
} /* lu1msp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1mSP */
/* Subroutine */ int lu1pen_(integer *m, integer *melim, integer *ncold, 
	integer *nspare, integer *ilast, integer *lpivc1, integer *lpivc2, 
	integer *lpivr1, integer *lpivr2, integer *lrow, integer *lenc, 
	integer *lenr, integer *locc, integer *locr, integer *indc, integer *
	indr, integer *ifill, integer *jfill)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, lc, ll, lr, lu, lc1, lc2, lr1, lr2, last;

/*     ------------------------------------------------------------------ */
/*     lu1pen deals with pending fill-in in the row file. */
/*     ifill(ll) says if a row involved in the new column of L */
/*               has to be updated.  If positive, it is the total */
/*               length of the final updated row. */
/*     jfill(lu) says if a column involved in the new row of U */
/*               contains any pending fill-ins.  If positive, it points */
/*               to the first fill-in in the column that has yet to be */
/*               added to the row file. */

/*     16 Apr 1989: First version of lu1pen. */
/*     23 Mar 2001: ilast used and updated. */
/*     30 Mar 2011: Loop 620 revised following advice from Ralf Östermark, */
/*                  School of Business and Economics at Åbo Akademi University */
/*                  after consulting Martti Louhivuori at the */
/*                  Centre of Scientific Computing in Helsinki. */
/*                  The previous version produced a core dump on the Cray XT. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ifill;
    --jfill;
    --lenc;
    --lenr;
    --locc;
    --locr;
    --indc;
    --indr;

    /* Function Body */
    ll = 0;
    i__1 = *lpivc2;
    for (lc = *lpivc1; lc <= i__1; ++lc) {
	++ll;
	if (ifill[ll] == 0) {
	    goto L650;
	}
/* Another row has pending fill. */
/* First, add some spare space at the end */
/* of the current last row. */
/* 30 Mar 2011: Better not to alter loop limits inside the lo */
/*              Some compilers might not freeze the limits */
/*              before entering the loop (even though the */
/*              Fortran standards say they should). */
	i__2 = *lrow + *nspare;
	for (l = *lrow + 1; l <= i__2; ++l) {
/*  lrow    = l */
	    indr[l] = 0;
/* L620: */
	}
	*lrow += *nspare;
/* Now move row i to the end of the row file. */
	i__ = indc[lc];
	*ilast = i__;
	lr1 = locr[i__];
	lr2 = lr1 + lenr[i__] - 1;
	locr[i__] = *lrow + 1;
	i__2 = lr2;
	for (lr = lr1; lr <= i__2; ++lr) {
	    ++(*lrow);
	    indr[*lrow] = indr[lr];
	    indr[lr] = 0;
/* L630: */
	}
	*lrow += ifill[ll];
L650:
	;
    }
/*        Scan all columns of  D  and insert the pending fill-in */
/*        into the row file. */
    lu = 1;
    i__1 = *lpivr2;
    for (lr = *lpivr1; lr <= i__1; ++lr) {
	++lu;
	if (jfill[lu] == 0) {
	    goto L680;
	}
	j = indr[lr];
	lc1 = locc[j] + jfill[lu] - 1;
	lc2 = locc[j] + lenc[j] - 1;
	i__2 = lc2;
	for (lc = lc1; lc <= i__2; ++lc) {
	    i__ = indc[lc] - *m;
	    if (i__ > 0) {
		indc[lc] = i__;
		last = locr[i__] + lenr[i__];
		indr[last] = j;
		++lenr[i__];
	    }
/* L670: */
	}
L680:
	;
    }
    return 0;
} /* lu1pen_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1pen */
/* Subroutine */ int lu1mxc_(integer *k1, integer *k2, integer *iq, 
	doublereal *a, integer *indc, integer *lenc, integer *locc)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, lc, lc1, lc2;
    static doublereal amax;
    static integer lenj;

/*     ------------------------------------------------------------------ */
/*     lu1mxc  moves the largest element in each of columns iq(k1:k2) */
/*     to the top of its column. */
/*     If k1 > k2, nothing happens. */

/*     06 May 2002: (and earlier) */
/*                  All columns k1:k2 must have one or more elements. */
/*     07 May 2002: Allow for empty columns.  The heap routines need to */
/*                  find 0.0 as the "largest element". */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iq;
    --a;
    --indc;
    --lenc;
    --locc;

    /* Function Body */
    i__1 = *k2;
    for (k = *k1; k <= i__1; ++k) {
	j = iq[k];
	lc1 = locc[j];
	lenj = lenc[j];
	if (lenj == 0) {
	    a[lc1] = 0.;
	} else {
/* The next 10 lines are equivalent to */
/* l      = idamax( lenc(j), a(lc1), 1 )  +  lc1 - 1 */
/* >>>>>>>> */
	    lc2 = lc1 + lenc[j] - 1;
	    amax = (d__1 = a[lc1], abs(d__1));
	    l = lc1;
	    i__2 = lc2;
	    for (lc = lc1 + 1; lc <= i__2; ++lc) {
		if (amax < (d__1 = a[lc], abs(d__1))) {
		    amax = (d__1 = a[lc], abs(d__1));
		    l = lc;
		}
	    }
/* >>>>>>>> */
	    if (l > lc1) {
		amax = a[l];
		a[l] = a[lc1];
		a[lc1] = amax;
		i__ = indc[l];
		indc[l] = indc[lc1];
		indc[lc1] = i__;
	    }
	}
    }
    return 0;
} /* lu1mxc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1mxc */
/* Subroutine */ int lu1mxr_(integer *k1, integer *k2, integer *ip, 
	doublereal *amaxr, doublereal *a, integer *indc, integer *lenc, 
	integer *locc, integer *indr, integer *lenr, integer *locr)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, lc, lr, lc1, lc2, lr1, lr2;
    static doublereal amax;

/*     ------------------------------------------------------------------ */
/*     lu1mxr  finds the largest element in each of row ip(k1:k2) */
/*     and stores it in Amaxr(*).  The nonzeros are stored column-wise */
/*     in (a,indc,lenc,locc) and their structure is row-wise */
/*     in (  indr,lenr,locr). */
/*     If k1 > k2, nothing happens. */

/*     11 Jun 2002: First version of lu1mxr. */
/*                  Allow for empty columns. */
/*     19 Dec 2004: Declare Amaxr(*) instead of Amaxr(k2) */
/*                  to stop grumbles from the NAG compiler. */
/*                  (Mentioned by Mick Pont.) */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ip;
    --amaxr;
    --a;
    --indc;
    --lenc;
    --locc;
    --indr;
    --lenr;
    --locr;

    /* Function Body */
    i__1 = *k2;
    for (k = *k1; k <= i__1; ++k) {
	amax = 0.;
	i__ = ip[k];
/* Find largest element in row i. */
	lr1 = locr[i__];
	lr2 = lr1 + lenr[i__] - 1;
	i__2 = lr2;
	for (lr = lr1; lr <= i__2; ++lr) {
	    j = indr[lr];
/* Find where  aij  is in column  j. */
	    lc1 = locc[j];
	    lc2 = lc1 + lenc[j] - 1;
	    i__3 = lc2;
	    for (lc = lc1; lc <= i__3; ++lc) {
		if (indc[lc] == i__) {
		    goto L230;
		}
	    }
L230:
/* Computing MAX */
	    d__2 = amax, d__3 = (d__1 = a[lc], abs(d__1));
	    amax = max(d__2,d__3);
	}
	amaxr[i__] = amax;
    }
    return 0;
} /* lu1mxr_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1mxr */
/* Subroutine */ int lu1or1_(integer *m, integer *n, integer *nelem, integer *
	lena, doublereal *small, doublereal *a, integer *indc, integer *indr, 
	integer *lenc, integer *lenr, doublereal *amax, integer *numnz, 
	integer *lerr, integer *inform__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, l, ldummy;

/*     ------------------------------------------------------------------ */
/*     lu1or1  organizes the elements of an  m by n  matrix  A  as */
/*     follows.  On entry, the parallel arrays   a, indc, indr, */
/*     contain  nelem  entries of the form     aij,    i,    j, */
/*     in any order.  nelem  must be positive. */

/*     Entries not larger than the input parameter  small  are treated as */
/*     zero and removed from   a, indc, indr.  The remaining entries are */
/*     defined to be nonzero.  numnz  returns the number of such nonzeros */
/*     and  Amax  returns the magnitude of the largest nonzero. */
/*     The arrays  lenc, lenr  return the number of nonzeros in each */
/*     column and row of  A. */

/*     inform = 0  on exit, except  inform = 1  if any of the indices in */
/*     indc, indr  imply that the element  aij  lies outside the  m by n */
/*     dimensions of  A. */

/*     xx Feb 1985: Original version. */
/*     17 Oct 2000: a, indc, indr now have size lena to allow nelem = 0. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --lenr;
    --lenc;
    --indr;
    --indc;
    --a;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lenr[i__] = 0;
/* L10: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	lenc[j] = 0;
/* L20: */
    }
    *amax = 0.;
    *numnz = *nelem;
    l = *nelem + 1;
    i__1 = *nelem;
    for (ldummy = 1; ldummy <= i__1; ++ldummy) {
	--l;
	if ((d__1 = a[l], abs(d__1)) > *small) {
	    i__ = indc[l];
	    j = indr[l];
/* Computing MAX */
	    d__2 = *amax, d__3 = (d__1 = a[l], abs(d__1));
	    *amax = max(d__2,d__3);
	    if (i__ < 1 || i__ > *m) {
		goto L910;
	    }
	    if (j < 1 || j > *n) {
		goto L910;
	    }
	    ++lenr[i__];
	    ++lenc[j];
	} else {
/*           Replace a negligible element by last element.  Since */
/*           we are going backwards, we know the last element is ok. */
	    a[l] = a[*numnz];
	    indc[l] = indc[*numnz];
	    indr[l] = indr[*numnz];
	    --(*numnz);
	}
/* L100: */
    }
    *inform__ = 0;
    return 0;
L910:
    *lerr = l;
    *inform__ = 1;
    return 0;
} /* lu1or1_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1or1 */
/* Subroutine */ int lu1or2_(integer *n, integer *numa, integer *lena, 
	doublereal *a, integer *inum, integer *jnum, integer *len, integer *
	loc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, ja, jb;
    static doublereal ace;
    static integer ice, jce;
    static doublereal acep;
    static integer icep, jcep;

/*     ------------------------------------------------------------------ */
/*     lu1or2  sorts a list of matrix elements  a(i,j)  into column */
/*     order, given  numa  entries  a(i,j),  i,  j  in the parallel */
/*     arrays  a, inum, jnum  respectively.  The matrix is assumed */
/*     to have  n  columns and an arbitrary number of rows. */

/*     On entry,  len(*)  must contain the length of each column. */

/*     On exit,  a(*) and inum(*)  are sorted,  jnum(*) = 0,  and */
/*     loc(j)  points to the start of column j. */

/*     lu1or2  is derived from  mc20ad,  a routine in the Harwell */
/*     Subroutine Library, author J. K. Reid. */

/*     xx Feb 1985: Original version. */
/*     17 Oct 2000: a, inum, jnum now have size lena to allow nelem = 0. */
/*     ------------------------------------------------------------------ */
/*     Set  loc(j)  to point to the beginning of column  j. */
    /* Parameter adjustments */
    --loc;
    --len;
    --jnum;
    --inum;
    --a;

    /* Function Body */
    l = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	loc[j] = l;
	l += len[j];
/* L150: */
    }
/*     Sort the elements into column order. */
/*     The algorithm is an in-place sort and is of order  numa. */
    i__1 = *numa;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        Establish the current entry. */
	jce = jnum[i__];
	if (jce == 0) {
	    goto L230;
	}
	ace = a[i__];
	ice = inum[i__];
	jnum[i__] = 0;
/*        Chain from current entry. */
	i__2 = *numa;
	for (j = 1; j <= i__2; ++j) {
/*           The current entry is not in the correct position. */
/*           Determine where to store it. */
	    l = loc[jce];
	    ++loc[jce];
/*           Save the contents of that location. */
	    acep = a[l];
	    icep = inum[l];
	    jcep = jnum[l];
/*           Store current entry. */
	    a[l] = ace;
	    inum[l] = ice;
	    jnum[l] = 0;
/*           If next current entry needs to be processed, */
/*           copy it into current entry. */
	    if (jcep == 0) {
		goto L230;
	    }
	    ace = acep;
	    ice = icep;
	    jce = jcep;
/* L200: */
	}
L230:
	;
    }
/*     Reset loc(j) to point to the start of column j. */
    ja = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jb = loc[j];
	loc[j] = ja;
	ja = jb;
/* L250: */
    }
    return 0;
} /* lu1or2_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1or2 */
/* Subroutine */ int lu1or3_(integer *m, integer *n, integer *lena, integer *
	indc, integer *lenc, integer *locc, integer *iw, integer *lerr, 
	integer *inform__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, l1, l2;

/*     ------------------------------------------------------------------ */
/*     lu1or3  looks for duplicate elements in an  m by n  matrix  A */
/*     defined by the column list  indc, lenc, locc. */
/*     iw  is used as a work vector of length  m. */

/*     xx Feb 1985: Original version. */
/*     17 Oct 2000: indc, indr now have size lena to allow nelem = 0. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --locc;
    --lenc;
    --indc;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iw[i__] = 0;
/* L100: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (lenc[j] > 0) {
	    l1 = locc[j];
	    l2 = l1 + lenc[j] - 1;
	    i__2 = l2;
	    for (l = l1; l <= i__2; ++l) {
		i__ = indc[l];
		if (iw[i__] == j) {
		    goto L910;
		}
		iw[i__] = j;
/* L150: */
	    }
	}
/* L200: */
    }
    *inform__ = 0;
    return 0;
L910:
    *lerr = l;
    *inform__ = 1;
    return 0;
} /* lu1or3_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1or3 */
/* Subroutine */ int lu1or4_(integer *m, integer *n, integer *nelem, integer *
	lena, integer *indc, integer *indr, integer *lenc, integer *lenr, 
	integer *locc, integer *locr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, l1, l2, lr, jdummy;

/*     ------------------------------------------------------------------ */
/*     lu1or4     constructs a row list  indr, locr */
/*     from a corresponding column list  indc, locc, */
/*     given the lengths of both columns and rows in  lenc, lenr. */

/*     xx Feb 1985: Original version. */
/*     17 Oct 2000: indc, indr now have size lena to allow nelem = 0. */
/*     ------------------------------------------------------------------ */
/*     Initialize  locr(i)  to point just beyond where the */
/*     last component of row  i  will be stored. */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --locc;
    --lenc;
    --indr;
    --indc;

    /* Function Body */
    l = 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l += lenr[i__];
	locr[i__] = l;
/* L10: */
    }
/*     By processing the columns backwards and decreasing  locr(i) */
/*     each time it is accessed, it will end up pointing to the */
/*     beginning of row  i  as required. */
    l2 = *nelem;
    j = *n + 1;
    i__1 = *n;
    for (jdummy = 1; jdummy <= i__1; ++jdummy) {
	--j;
	if (lenc[j] > 0) {
	    l1 = locc[j];
	    i__2 = l2;
	    for (l = l1; l <= i__2; ++l) {
		i__ = indc[l];
		lr = locr[i__] - 1;
		locr[i__] = lr;
		indr[lr] = j;
/* L30: */
	    }
	    l2 = l1 - 1;
	}
/* L40: */
    }
    return 0;
} /* lu1or4_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1or4 */
/* Subroutine */ int lu1pq1_(integer *m, integer *n, integer *len, integer *
	iperm, integer *loc, integer *inv, integer *num)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, l, nz, nzero;

/*     ------------------------------------------------------------------ */
/*     lu1pq1  constructs a permutation  iperm  from the array  len. */

/*     On entry: */
/*     len(i)  holds the number of nonzeros in the i-th row (say) */
/*             of an m by n matrix. */
/*     num(*)  can be anything (workspace). */

/*     On exit: */
/*     iperm   contains a list of row numbers in the order */
/*             rows of length 0,  rows of length 1,..., rows of length n. */
/*     loc(nz) points to the first row containing  nz  nonzeros, */
/*             nz = 1, n. */
/*     inv(i)  points to the position of row i within iperm(*). */
/*     ------------------------------------------------------------------ */
/*     Count the number of rows of each length. */
    /* Parameter adjustments */
    --inv;
    --iperm;
    --len;
    --num;
    --loc;

    /* Function Body */
    nzero = 0;
    i__1 = *n;
    for (nz = 1; nz <= i__1; ++nz) {
	num[nz] = 0;
	loc[nz] = 0;
/* L10: */
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nz = len[i__];
	if (nz == 0) {
	    ++nzero;
	} else {
	    ++num[nz];
	}
/* L20: */
    }
/*     Set starting locations for each length. */
    l = nzero + 1;
    i__1 = *n;
    for (nz = 1; nz <= i__1; ++nz) {
	loc[nz] = l;
	l += num[nz];
	num[nz] = 0;
/* L60: */
    }
/*     Form the list. */
    nzero = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nz = len[i__];
	if (nz == 0) {
	    ++nzero;
	    iperm[nzero] = i__;
	} else {
	    l = loc[nz] + num[nz];
	    iperm[l] = i__;
	    ++num[nz];
	}
/* L100: */
    }
/*     Define the inverse of iperm. */
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
	i__ = iperm[l];
	inv[i__] = l;
/* L120: */
    }
    return 0;
} /* lu1pq1_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1pq1 */
/* Subroutine */ int lu1pq2_(integer *nzpiv, integer *nzchng, integer *indr, 
	integer *lenold, integer *lennew, integer *iqloc, integer *iq, 
	integer *iqinv)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, l, lr, nz, jnew, lnew, next, nznew;

/*     =============================================================== */
/*     lu1pq2 frees the space occupied by the pivot row, */
/*     and updates the column permutation iq. */

/*     Also used to free the pivot column and update the row perm ip. */

/*     nzpiv   (input)    is the length of the pivot row (or column). */
/*     nzchng  (output)   is the net change in total nonzeros. */

/*     14 Apr 1989  First version. */
/*     =============================================================== */
    /* Parameter adjustments */
    --lenold;
    --indr;
    --lennew;
    --iqloc;
    --iq;
    --iqinv;

    /* Function Body */
    *nzchng = 0;
    i__1 = *nzpiv;
    for (lr = 1; lr <= i__1; ++lr) {
	j = indr[lr];
	indr[lr] = 0;
	nz = lenold[lr];
	nznew = lennew[j];
	if (nz != nznew) {
	    l = iqinv[j];
	    *nzchng += nznew - nz;
/*           l above is the position of column j in iq  (so j = iq(l)). */
	    if (nz < nznew) {
/*              Column  j  has to move towards the end of  iq. */
L110:
		next = nz + 1;
		lnew = iqloc[next] - 1;
		if (lnew != l) {
		    jnew = iq[lnew];
		    iq[l] = jnew;
		    iqinv[jnew] = l;
		}
		l = lnew;
		iqloc[next] = lnew;
		nz = next;
		if (nz < nznew) {
		    goto L110;
		}
	    } else {
/*              Column  j  has to move towards the front of  iq. */
L120:
		lnew = iqloc[nz];
		if (lnew != l) {
		    jnew = iq[lnew];
		    iq[l] = jnew;
		    iqinv[jnew] = l;
		}
		l = lnew;
		iqloc[nz] = lnew + 1;
		--nz;
		if (nz > nznew) {
		    goto L120;
		}
	    }
	    iq[lnew] = j;
	    iqinv[j] = lnew;
	}
/* L200: */
    }
    return 0;
} /* lu1pq2_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1pq2 */
/* Subroutine */ int lu1pq3_(integer *n, integer *len, integer *iperm, 
	integer *iw, integer *nrank)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, nzero;

/*     ------------------------------------------------------------------ */
/*     lu1pq3  looks at the permutation  iperm(*)  and moves any entries */
/*     to the end whose corresponding length  len(*)  is zero. */

/*     09 Feb 1994: Added work array iw(*) to improve efficiency. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;
    --iperm;
    --len;

    /* Function Body */
    *nrank = 0;
    nzero = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__ = iperm[k];
	if (len[i__] == 0) {
	    ++nzero;
	    iw[nzero] = i__;
	} else {
	    ++(*nrank);
	    iperm[*nrank] = i__;
	}
/* L10: */
    }
    i__1 = nzero;
    for (k = 1; k <= i__1; ++k) {
	iperm[*nrank + k] = iw[k];
/* L20: */
    }
    return 0;
} /* lu1pq3_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1pq3 */
/* Subroutine */ int lu1rec_(integer *n, logical *reals, integer *luparm, 
	integer *ltop, integer *lena, doublereal *a, integer *ind, integer *
	len, integer *loc)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 lu1rec.  File compressed from\002,i10"
	    ",\002   to\002,i10,l3,\002  nempty =\002,i8)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, k, l, leni, nout, ilast, klast, lprint, nempty;

    /* Fortran I/O blocks */
    static cilist io___365 = { 0, 0, 0, fmt_1000, 0 };


/*     ------------------------------------------------------------------ */
/*     00 Jun 1983: Original version of lu1rec followed John Reid's */
/*                  compression routine in LA05.  It recovered */
/*                  space in ind(*) and optionally a(*) */
/*                  by eliminating entries with ind(l) = 0. */
/*                  The elements of ind(*) could not be negative. */
/*                  If len(i) was positive, entry i contained */
/*                  that many elements, starting at  loc(i). */
/*                  Otherwise, entry i was eliminated. */

/*     23 Mar 2001: Realised we could have len(i) = 0 in rare cases! */
/*                  (Mostly during TCP when the pivot row contains */
/*                  a column of length 1 that couldn't be a pivot.) */
/*                  Revised storage scheme to */
/*                     keep        entries with       ind(l) >  0, */
/*                     squeeze out entries with -n <= ind(l) <= 0, */
/*                  and to allow len(i) = 0. */
/*                  Empty items are moved to the end of the compressed */
/*                  ind(*) and/or a(*) arrays are given one empty space. */
/*                  Items with len(i) < 0 are still eliminated. */

/*     27 Mar 2001: Decided to use only ind(l) > 0 and = 0 in lu1fad. */
/*                  Still have to keep entries with len(i) = 0. */

/*     On exit: */
/*     ltop         is the length of useful entries in ind(*), a(*). */
/*     ind(ltop+1)  is "i" such that len(i), loc(i) belong to the last */
/*                  item in ind(*), a(*). */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --loc;
    --len;
    --luparm;
    --ind;
    --a;

    /* Function Body */
    nempty = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	leni = len[i__];
	if (leni > 0) {
	    l = loc[i__] + leni - 1;
	    len[i__] = ind[l];
	    ind[l] = -(*n + i__);
	} else if (leni == 0) {
	    ++nempty;
	}
/* L10: */
    }
    k = 0;
    klast = 0;
/* Previous k */
    ilast = 0;
/* Last entry moved. */
    i__1 = *ltop;
    for (l = 1; l <= i__1; ++l) {
	i__ = ind[l];
	if (i__ > 0) {
	    ++k;
	    ind[k] = i__;
	    if (*reals) {
		a[k] = a[l];
	    }
	} else if (i__ < -(*n)) {
/*           This is the end of entry  i. */
	    i__ = -(i__ + *n);
	    ilast = i__;
	    ++k;
	    ind[k] = len[i__];
	    if (*reals) {
		a[k] = a[l];
	    }
	    loc[i__] = klast + 1;
	    len[i__] = k - klast;
	    klast = k;
	}
/* L20: */
    }
/*     Move any empty items to the end, adding 1 free entry for each. */
    if (nempty > 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (len[i__] == 0) {
		++k;
		loc[i__] = k;
		ind[k] = 0;
		ilast = i__;
	    }
	}
    }
    nout = luparm[1];
    lprint = luparm[2];
    if (lprint >= 50) {
	io___365.ciunit = nout;
	s_wsfe(&io___365);
	do_fio(&c__1, (char *)&(*ltop), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*reals), (ftnlen)sizeof(logical));
	do_fio(&c__1, (char *)&nempty, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    ++luparm[26];
/*     Return ilast in ind(ltop + 1). */
/* ncp */
    *ltop = k;
    ind[*ltop + 1] = ilast;
    return 0;
} /* lu1rec_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1rec */
/* Subroutine */ int lu1slk_(integer *m, integer *n, integer *lena, integer *
	iq, integer *iqloc, doublereal *a, integer *locc, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer j, lq, lc1, lq1, lq2;

/*     ------------------------------------------------------------------ */
/*     lu1slk  sets w(j) > 0 if column j is a unit vector. */

/*     21 Nov 2000: First version.  lu1fad needs it for TCP. */
/*                  Note that w(*) is nominally an integer array, */
/*                  but the only spare space is the double array w(*). */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iqloc;
    --w;
    --locc;
    --iq;
    --a;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	w[j] = 0.;
    }
    lq1 = iqloc[1];
    lq2 = *n;
    if (*m > 1) {
	lq2 = iqloc[2] - 1;
    }
    i__1 = lq2;
    for (lq = lq1; lq <= i__1; ++lq) {
	j = iq[lq];
	lc1 = locc[j];
	if ((d__1 = a[lc1], abs(d__1)) == 1.) {
	    w[j] = 1.;
	}
    }
    return 0;
} /* lu1slk_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1slk */
/* Subroutine */ int lu1ful_(integer *m, integer *n, integer *lena, integer *
	lend, integer *lu1, logical *tpp, integer *mleft, integer *nleft, 
	integer *nrank, integer *nrowu, integer *lenl, integer *lenu, integer 
	*nsing, logical *keeplu, doublereal *small, doublereal *a, doublereal 
	*d__, integer *indc, integer *indr, integer *ip, integer *iq, integer 
	*lenc, integer *lenr, integer *locc, integer *ipinv, integer *ipvt)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, l1, l2;
    static doublereal ai, aj;
    static integer la, lc, ld, ll, lq, lu, lc1, lc2, lkk, lkn, ncold, ibest, 
	    jbest;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nrowd;
    extern /* Subroutine */ int lu1dcp_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    lu1dpp_(doublereal *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *);
    static integer ldbase, ipbase, ldiagu;

/*     ------------------------------------------------------------------ */
/*     lu1ful computes a dense (full) LU factorization of the */
/*     mleft by nleft matrix that remains to be factored at the */
/*     beginning of the nrowu-th pass through the main loop of lu1fad. */

/*     02 May 1989: First version. */
/*     05 Feb 1994: Column interchanges added to lu1DPP. */
/*     08 Feb 1994: ipinv reconstructed, since lu1pq3 may alter ip. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     If lu1pq3 moved any empty rows, reset ipinv = inverse of ip. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ipvt;
    --ipinv;
    --lenr;
    --ip;
    --locc;
    --lenc;
    --iq;
    --indr;
    --indc;
    --a;
    --d__;

    /* Function Body */
    if (*nrank < *m) {
	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
	    i__ = ip[l];
	    ipinv[i__] = l;
/* L100: */
	}
    }
/*     ------------------------------------------------------------------ */
/*     Copy the remaining matrix into the dense matrix D. */
/*     ------------------------------------------------------------------ */
/*     call dload ( lenD, zero, d, 1 ) */
    i__1 = *lend;
    for (j = 1; j <= i__1; ++j) {
	d__[j] = 0.;
    }
    ipbase = *nrowu - 1;
    ldbase = 1 - *nrowu;
    i__1 = *n;
    for (lq = *nrowu; lq <= i__1; ++lq) {
	j = iq[lq];
	lc1 = locc[j];
	lc2 = lc1 + lenc[j] - 1;
	i__2 = lc2;
	for (lc = lc1; lc <= i__2; ++lc) {
	    i__ = indc[lc];
	    ld = ldbase + ipinv[i__];
	    d__[ld] = a[lc];
/* L150: */
	}
	ldbase += *mleft;
/* L200: */
    }
/*     ------------------------------------------------------------------ */
/*     Call our favorite dense LU factorizer. */
/*     ------------------------------------------------------------------ */
    if (*tpp) {
	lu1dpp_(&d__[1], mleft, mleft, nleft, small, nsing, &ipvt[1], &iq[*
		nrowu]);
    } else {
	lu1dcp_(&d__[1], mleft, mleft, nleft, small, nsing, &ipvt[1], &iq[*
		nrowu]);
    }
/*     ------------------------------------------------------------------ */
/*     Move D to the beginning of A, */
/*     and pack L and U at the top of a, indc, indr. */
/*     In the process, apply the row permutation to ip. */
/*     lkk points to the diagonal of U. */
/*     ------------------------------------------------------------------ */
    dcopy_(lend, &d__[1], &c__1, &a[1], &c__1);
    ldiagu = *lena - *n;
    lkk = 1;
    lkn = *lend - *mleft + 1;
    lu = *lu1;
    i__1 = min(*mleft,*nleft);
    for (k = 1; k <= i__1; ++k) {
	l1 = ipbase + k;
	l2 = ipbase + ipvt[k];
	if (l1 != l2) {
	    i__ = ip[l1];
	    ip[l1] = ip[l2];
	    ip[l2] = i__;
	}
	ibest = ip[l1];
	jbest = iq[l1];
	if (*keeplu) {
/*           =========================================================== */
/*           Pack the next column of L. */
/*           =========================================================== */
	    la = lkk;
	    ll = lu;
	    nrowd = 1;
	    i__2 = *mleft;
	    for (i__ = k + 1; i__ <= i__2; ++i__) {
		++la;
		ai = a[la];
		if (abs(ai) > *small) {
		    ++nrowd;
		    --ll;
		    a[ll] = ai;
		    indc[ll] = ip[ipbase + i__];
		    indr[ll] = ibest;
		}
/* L410: */
	    }
/*           =========================================================== */
/*           Pack the next row of U. */
/*           We go backwards through the row of D */
/*           so the diagonal ends up at the front of the row of  U. */
/*           Beware -- the diagonal may be zero. */
/*           =========================================================== */
	    la = lkn + *mleft;
	    lu = ll;
	    ncold = 0;
	    i__2 = k;
	    for (j = *nleft; j >= i__2; --j) {
		la -= *mleft;
		aj = a[la];
		if (abs(aj) > *small || j == k) {
		    ++ncold;
		    --lu;
		    a[lu] = aj;
		    indr[lu] = iq[ipbase + j];
		}
/* L420: */
	    }
	    lenr[ibest] = -ncold;
	    lenc[jbest] = -nrowd;
	    *lenl = *lenl + nrowd - 1;
	    *lenu += ncold;
	    ++lkn;
	} else {
/*           =========================================================== */
/*           Store just the diagonal of U, in natural order. */
/*           =========================================================== */
	    a[ldiagu + jbest] = a[lkk];
	}
	lkk = lkk + *mleft + 1;
/* L450: */
    }
    return 0;
} /* lu1ful_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1ful */
/* Subroutine */ int lu1dpp_(doublereal *a, integer *lda, integer *m, integer 
	*n, doublereal *small, integer *nsing, integer *ipvt, integer *iq)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer kp1, last;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ranku;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer lencol;

/*     ------------------------------------------------------------------ */
/*     lu1DPP factors a dense m x n matrix A by Gaussian elimination, */
/*     using row interchanges for stability, as in dgefa from LINPACK. */
/*     This version also uses column interchanges if all elements in a */
/*     pivot column are smaller than (or equal to) "small".  Such columns */
/*     are changed to zero and permuted to the right-hand end. */

/*     As in LINPACK, ipvt(*) keeps track of pivot rows. */
/*     Rows of U are interchanged, but we don't have to physically */
/*     permute rows of L.  In contrast, column interchanges are applied */
/*     directly to the columns of both L and U, and to the column */
/*     permutation vector iq(*). */

/*     02 May 1989: First version derived from dgefa */
/*                  in LINPACK (version dated 08/14/78). */
/*     05 Feb 1994: Generalized to treat rectangular matrices */
/*                  and use column interchanges when necessary. */
/*                  ipvt is retained, but column permutations are applied */
/*                  directly to iq(*). */
/*     21 Dec 1994: Bug found via example from Steve Dirkse. */
/*                  Loop 100 added to set ipvt(*) for singular rows. */
/*     26 Mar 2006: nsing redefined (see below). */
/*                  Changed to implicit none. */
/*     ------------------------------------------------------------------ */

/*     On entry: */

/*        a       Array holding the matrix A to be factored. */

/*        lda     The leading dimension of the array  a. */

/*        m       The number of rows    in  A. */

/*        n       The number of columns in  A. */

/*        small   A drop tolerance.  Must be zero or positive. */

/*     On exit: */

/*        a       An upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        nsing   Number of singularities detected. */
/*                26 Mar 2006: nsing redefined to be more meaningful. */
/*                Users may define rankU = n - nsing and regard */
/*                U as upper-trapezoidal, with the first rankU columns */
/*                being triangular and the rest trapezoidal. */
/*                It would be better to return rankU, but we still */
/*                return nsing for compatibility (even though lu1fad */
/*                no longer uses it). */

/*        ipvt    Records the pivot rows. */

/*        iq      A vector to which column interchanges are applied. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ipvt;
    --iq;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ranku = 0;
    k = 1;
    last = *n;
/*     ------------------------------------------------------------------ */
/*     Start of elimination loop. */
/*     ------------------------------------------------------------------ */
L10:
    kp1 = k + 1;
    lencol = *m - k + 1;
/* Find l, the pivot row. */
    l = idamax_(&lencol, &a[k + k * a_dim1], &c__1) + k - 1;
    ipvt[k] = l;
    if ((d__1 = a[l + k * a_dim1], abs(d__1)) <= *small) {
/* ============================================================== */
/* Do column interchange, changing old pivot column to zero. */
/* Reduce "last" and try again with same k. */
/* ============================================================== */
	j = iq[last];
	iq[last] = iq[k];
	iq[k] = j;
	i__1 = k - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = a[i__ + last * a_dim1];
	    a[i__ + last * a_dim1] = a[i__ + k * a_dim1];
	    a[i__ + k * a_dim1] = t;
	}
	i__1 = *m;
	for (i__ = k; i__ <= i__1; ++i__) {
	    t = a[i__ + last * a_dim1];
	    a[i__ + last * a_dim1] = 0.;
	    a[i__ + k * a_dim1] = t;
	}
	--last;
	if (k <= last) {
	    goto L10;
	}
    } else {
	++ranku;
	if (k < *m) {
/* =========================================================== */
/* Do row interchange if necessary. */
/* =========================================================== */
	    if (l != k) {
		t = a[l + k * a_dim1];
		a[l + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = t;
	    }
/* =========================================================== */
/* Compute multipliers. */
/* Do row elimination with column indexing. */
/* =========================================================== */
	    t = -1. / a[k + k * a_dim1];
	    i__1 = *m - k;
	    dscal_(&i__1, &t, &a[kp1 + k * a_dim1], &c__1);
	    i__1 = last;
	    for (j = kp1; j <= i__1; ++j) {
		t = a[l + j * a_dim1];
		if (l != k) {
		    a[l + j * a_dim1] = a[k + j * a_dim1];
		    a[k + j * a_dim1] = t;
		}
		i__2 = *m - k;
		daxpy_(&i__2, &t, &a[kp1 + k * a_dim1], &c__1, &a[kp1 + j * 
			a_dim1], &c__1);
	    }
	    ++k;
	    if (k <= last) {
		goto L10;
	    }
	}
    }
/* Set ipvt(*) for singular rows. */
    i__1 = *m;
    for (k = last + 1; k <= i__1; ++k) {
	ipvt[k] = k;
/* L100: */
    }
    *nsing = *n - ranku;
    return 0;
} /* lu1dpp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1DPP */
/* Subroutine */ int lu1dcp_(doublereal *a, integer *lda, integer *m, integer 
	*n, doublereal *small, integer *nsing, integer *ipvt, integer *iq)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer kp1, imax, jmax, jnew, last;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal ajmax;
    static integer jlast, ranku;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal aijmax;
    static integer lencol;

/*     ------------------------------------------------------------------ */
/*     lu1DCP factors a dense m x n matrix A by Gaussian elimination, */
/*     using Complete Pivoting (row and column interchanges) for */
/*     stability. */
/*     This version also uses column interchanges if all elements in a */
/*     pivot column are smaller than (or equal to) "small".  Such columns */
/*     are changed to zero and permuted to the right-hand end. */

/*     As in LINPACK's dgefa, ipvt(*) keeps track of pivot rows. */
/*     Rows of U are interchanged, but we don't have to physically */
/*     permute rows of L.  In contrast, column interchanges are applied */
/*     directly to the columns of both L and U, and to the column */
/*     permutation vector iq(*). */

/*     01 May 2002: First dense Complete Pivoting, derived from lu1DPP. */
/*     07 May 2002: Another break needed at end of first loop. */
/*     26 Mar 2006: Cosmetic mods while looking for "nsing" bug when m<n. */
/*                  nsing redefined (see below). */
/*                  Changed to implicit none. */
/*     ------------------------------------------------------------------ */

/*     On entry: */

/*        a       Array holding the matrix A to be factored. */

/*        lda     The leading dimension of the array  a. */

/*        m       The number of rows    in  A. */

/*        n       The number of columns in  A. */

/*        small   A drop tolerance.  Must be zero or positive. */

/*     On exit: */

/*        a       An upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        nsing   Number of singularities detected. */
/*                26 Mar 2006: nsing redefined to be more meaningful. */
/*                Users may define rankU = n - nsing and regard */
/*                U as upper-trapezoidal, with the first rankU columns */
/*                being triangular and the rest trapezoidal. */
/*                It would be better to return rankU, but we still */
/*                return nsing for compatibility (even though lu1fad */
/*                no longer uses it). */

/*        ipvt    Records the pivot rows. */

/*        iq      A vector to which column interchanges are applied. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ipvt;
    --iq;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ranku = 0;
    lencol = *m + 1;
    last = *n;
/* ----------------------------------------------------------------- */
/* Start of elimination loop. */
/* ----------------------------------------------------------------- */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	--lencol;
/* Find the biggest aij in row imax and column jmax. */
	aijmax = 0.;
	imax = k;
	jmax = k;
	jlast = last;
	i__2 = jlast;
	for (j = k; j <= i__2; ++j) {
L10:
	    l = idamax_(&lencol, &a[k + j * a_dim1], &c__1) + k - 1;
	    ajmax = (d__1 = a[l + j * a_dim1], abs(d__1));
	    if (ajmax <= *small) {
/* ======================================================== */
/* Do column interchange, changing old column to zero. */
/* Reduce  "last"  and try again with same j. */
/* ======================================================== */
		jnew = iq[last];
		iq[last] = iq[j];
		iq[j] = jnew;
		i__3 = k - 1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    t = a[i__ + last * a_dim1];
		    a[i__ + last * a_dim1] = a[i__ + j * a_dim1];
		    a[i__ + j * a_dim1] = t;
		}
		i__3 = *m;
		for (i__ = k; i__ <= i__3; ++i__) {
		    t = a[i__ + last * a_dim1];
		    a[i__ + last * a_dim1] = 0.;
		    a[i__ + j * a_dim1] = t;
		}
		--last;
		if (j <= last) {
		    goto L10;
		}
/* repeat */
		goto L200;
/* break */
	    }
/* Check if this column has biggest aij so far. */
	    if (aijmax < ajmax) {
		aijmax = ajmax;
		imax = l;
		jmax = j;
	    }
	    if (j >= last) {
		goto L200;
	    }
/* break */
	}
L200:
	ipvt[k] = imax;
	if (jmax != k) {
/* ========================================================== */
/* Do column interchange (k and jmax). */
/* ========================================================== */
	    jnew = iq[jmax];
	    iq[jmax] = iq[k];
	    iq[k] = jnew;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		t = a[i__ + jmax * a_dim1];
		a[i__ + jmax * a_dim1] = a[i__ + k * a_dim1];
		a[i__ + k * a_dim1] = t;
	    }
	}
	if (k < *m) {
/* =========================================================== */
/* Do row interchange if necessary. */
/* =========================================================== */
	    t = a[imax + k * a_dim1];
	    if (imax != k) {
		a[imax + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = t;
	    }
/* =========================================================== */
/* Compute multipliers. */
/* Do row elimination with column indexing. */
/* =========================================================== */
	    t = -1. / t;
	    i__2 = *m - k;
	    dscal_(&i__2, &t, &a[kp1 + k * a_dim1], &c__1);
	    i__2 = last;
	    for (j = kp1; j <= i__2; ++j) {
		t = a[imax + j * a_dim1];
		if (imax != k) {
		    a[imax + j * a_dim1] = a[k + j * a_dim1];
		    a[k + j * a_dim1] = t;
		}
		i__3 = *m - k;
		daxpy_(&i__3, &t, &a[kp1 + k * a_dim1], &c__1, &a[kp1 + j * 
			a_dim1], &c__1);
	    }
	} else {
	    goto L500;
/* break */
	}
	if (k >= last) {
	    goto L500;
	}
/* break */
    }
/* Set ipvt(*) for singular rows. */
L500:
    i__1 = *m;
    for (k = last + 1; k <= i__1; ++k) {
	ipvt[k] = k;
    }
    *nsing = *n - ranku;
    return 0;
} /* lu1dcp_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  lusol2.f */

/*     Hbuild   Hchange  Hdelete  Hdown    Hinsert  Hup */

/*     Heap-management routines for LUSOL's lu1fac. */
/*     May be useful for other applications. */

/* 11 Feb 2002: MATLAB version derived from "Algorithms" by R. Sedgewick. */
/* 03 Mar 2002: F77    version derived from MATLAB version. */
/* 07 May 2002: Safeguard input parameters k, N, Nk. */
/*              We don't want them to be output! */
/* 19 Dec 2004: Hdelete: Nin is new input parameter for length of Hj, Ha. */
/* 19 Dec 2004: Current version of lusol2.f. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     For LUSOL, the heap structure involves three arrays of length N. */
/*     N        is the current number of entries in the heap. */
/*     Ha(1:N)  contains the values that the heap is partially sorting. */
/*              For LUSOL they are double precision values -- the largest */
/*              element in each remaining column of the updated matrix. */
/*              The biggest entry is in Ha(1), the top of the heap. */
/*     Hj(1:N)  contains column numbers j. */
/*              Ha(k) is the biggest entry in column j = Hj(k). */
/*     Hk(1:N)  contains indices within the heap.  It is the */
/*              inverse of Hj(1:N), so  k = Hk(j)  <=>  j = Hj(k). */
/*              Column j is entry k in the heap. */
/*     hops     is the number of heap operations, */
/*              i.e., the number of times an entry is moved */
/*              (the number of "hops" up or down the heap). */
/*     Together, Hj and Hk let us find values inside the heap */
/*     whenever we want to change one of the values in Ha. */
/*     For other applications, Ha may need to be some other data type, */
/*     like the keys that sort routines operate on. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu1DCP */
/* Subroutine */ int hbuild_(doublereal *ha, integer *hj, integer *hk, 
	integer *n, integer *nk, integer *hops)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer h__, k;
    static doublereal v;
    static integer kk, jv, nkk;
    extern /* Subroutine */ int hinsert_(doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *);

/*     ================================================================== */
/*     Hbuild initializes the heap by inserting each element of Ha. */
/*     Input:  Ha, Hj. */
/*     Output: Ha, Hj, Hk, hops. */

/*     01 May 2002: Use k for new length of heap, not k-1 for old length. */
/*     05 May 2002: Use kk in call to stop loop variable k being altered. */
/*                  (Actually Hinsert no longer alters that parameter.) */
/*     07 May 2002: ftnchek wants us to protect Nk, Ha(k), Hj(k) too. */
/*     07 May 2002: Current version of Hbuild. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hj;
    --ha;
    --hk;

    /* Function Body */
    nkk = *nk;
    *hops = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = k;
	v = ha[k];
	jv = hj[k];
	hinsert_(&ha[1], &hj[1], &hk[1], &kk, &nkk, &v, &jv, &h__);
	*hops += h__;
    }
    return 0;
} /* hbuild_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine Hbuild */
/* Subroutine */ int hchange_(doublereal *ha, integer *hj, integer *hk, 
	integer *n, integer *nk, integer *k, doublereal *v, integer *jv, 
	integer *hops)
{
    static doublereal v1;
    static integer kx, nx;
    extern /* Subroutine */ int hup_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer nkx;
    extern /* Subroutine */ int hdown_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

/*     ================================================================== */
/*     Hchange changes Ha(k) to v in heap of length N. */

/*     01 May 2002: Need Nk for length of Hk. */
/*     07 May 2002: Protect input parameters N, Nk, k. */
/*     07 May 2002: Current version of Hchange. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hj;
    --ha;
    --hk;

    /* Function Body */
    nx = *n;
    nkx = *nk;
    kx = *k;
    v1 = ha[*k];
    ha[*k] = *v;
    hj[*k] = *jv;
    hk[*jv] = *k;
    if (v1 < *v) {
	hup_(&ha[1], &hj[1], &hk[1], &nx, &nkx, &kx, hops);
    } else {
	hdown_(&ha[1], &hj[1], &hk[1], &nx, &nkx, &kx, hops);
    }
    return 0;
} /* hchange_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine Hchange */
/* Subroutine */ int hdelete_(doublereal *ha, integer *hj, integer *hk, 
	integer *nin, integer *n, integer *nk, integer *k, integer *hops)
{
    static doublereal v;
    static integer jv, kx, nx, nkx;
    extern /* Subroutine */ int hchange_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *);

/*     ================================================================== */
/*     Hdelete deletes Ha(k) from heap of length N. */

/*     03 Apr 2002: Current version of Hdelete. */
/*     01 May 2002: Need Nk for length of Hk. */
/*     07 May 2002: Protect input parameters N, Nk, k. */
/*     19 Dec 2004: Nin is new input parameter for length of Hj, Ha. */
/*     19 Dec 2004: Current version of Hdelete. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hj;
    --ha;
    --hk;

    /* Function Body */
    kx = *k;
    nkx = *nk;
    nx = *n;
    v = ha[*n];
    jv = hj[*n];
    --(*n);
    *hops = 0;
    if (*k <= *n) {
	hchange_(&ha[1], &hj[1], &hk[1], &nx, &nkx, &kx, &v, &jv, hops);
    }
    return 0;
} /* hdelete_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine Hdelete */
/* Subroutine */ int hdown_(doublereal *ha, integer *hj, integer *hk, integer 
	*n, integer *nk, integer *kk, integer *hops)
{
    static integer j, k;
    static doublereal v;
    static integer n2, jj, jv;

/*     ================================================================== */
/*     Hdown  updates heap by moving down tree from node k. */

/*     01 May 2002: Need Nk for length of Hk. */
/*     05 May 2002: Change input paramter k to kk to stop k being output. */
/*     05 May 2002: Current version of Hdown. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hj;
    --ha;
    --hk;

    /* Function Body */
    k = *kk;
    *hops = 0;
    v = ha[k];
    jv = hj[k];
    n2 = *n / 2;
/*     while 1 */
L100:
    if (k > n2) {
	goto L200;
    }
/* break */
    ++(*hops);
    j = k + k;
    if (j < *n) {
	if (ha[j] < ha[j + 1]) {
	    ++j;
	}
    }
    if (v >= ha[j]) {
	goto L200;
    }
/* break */
    ha[k] = ha[j];
    jj = hj[j];
    hj[k] = jj;
    hk[jj] = k;
    k = j;
    goto L100;
/*     end while */
L200:
    ha[k] = v;
    hj[k] = jv;
    hk[jv] = k;
    return 0;
} /* hdown_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine Hdown */
/* Subroutine */ int hinsert_(doublereal *ha, integer *hj, integer *hk, 
	integer *n, integer *nk, doublereal *v, integer *jv, integer *hops)
{
    static integer kk, nkk;
    extern /* Subroutine */ int hup_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer nnew;

/*     ================================================================== */
/*     Hinsert inserts (v,jv) into heap of length N-1 */
/*     to make heap of length N. */

/*     03 Apr 2002: First version of Hinsert. */
/*     01 May 2002: Require N to be final length, not old length. */
/*                  Need Nk for length of Hk. */
/*     07 May 2002: Protect input parameters N, Nk. */
/*     07 May 2002: Current version of Hinsert. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hj;
    --ha;
    --hk;

    /* Function Body */
    nnew = *n;
    nkk = *nk;
    kk = nnew;
    ha[nnew] = *v;
    hj[nnew] = *jv;
    hk[*jv] = nnew;
    hup_(&ha[1], &hj[1], &hk[1], &nnew, &nkk, &kk, hops);
    return 0;
} /* hinsert_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine Hinsert */
/* Subroutine */ int hup_(doublereal *ha, integer *hj, integer *hk, integer *
	n, integer *nk, integer *kk, integer *hops)
{
    static integer j, k;
    static doublereal v;
    static integer k2, jv;

/*     ================================================================== */
/*     Hup updates heap by moving up tree from node k. */

/*     01 May 2002: Need Nk for length of Hk. */
/*     05 May 2002: Change input paramter k to kk to stop k being output. */
/*     05 May 2002: Current version of Hup. */
/*     ================================================================== */
    /* Parameter adjustments */
    --hj;
    --ha;
    --hk;

    /* Function Body */
    k = *kk;
    *hops = 0;
    v = ha[k];
    jv = hj[k];
/*     while 1 */
L100:
    if (k < 2) {
	goto L200;
    }
/* break */
    k2 = k / 2;
    if (v < ha[k2]) {
	goto L200;
    }
/* break */
    ++(*hops);
    ha[k] = ha[k2];
    j = hj[k2];
    hj[k] = j;
    hk[j] = k;
    k = k2;
    goto L100;
/*     end while */
L200:
    ha[k] = v;
    hj[k] = jv;
    hk[jv] = k;
    return 0;
} /* hup_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  lusol6a.f */

/*     lu6sol   lu6L     lu6Lt     lu6U     Lu6Ut   lu6LD */
/*     lu6chk */

/* 26 Apr 2002: lu6 routines put into a separate file. */
/* 15 Dec 2002: lu6sol modularized via lu6L, lu6Lt, lu6U, lu6Ut. */
/*              lu6LD implemented to allow solves with LDL' or L|D|L'. */
/* 23 Apr 2004: lu6chk modified.  TRP can judge singularity better */
/*              by comparing all diagonals to DUmax. */
/* 27 Jun 2004: lu6chk.  Allow write only if nout .gt. 0 . */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine Hup */
/* Subroutine */ int lu6sol_(integer *mode, integer *m, integer *n, 
	doublereal *v, doublereal *w, integer *lena, integer *luparm, 
	doublereal *parmlu, doublereal *a, integer *indc, integer *indr, 
	integer *ip, integer *iq, integer *lenc, integer *lenr, integer *locc,
	 integer *locr, integer *inform__)
{
    extern /* Subroutine */ int lu6l_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *), lu6u_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *), lu6ld_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *), lu6lt_(
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, integer *), 
	    lu6ut_(integer *, integer *, integer *, doublereal *, doublereal *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/*     lu6sol  uses the factorization  A = L U  as follows: */

/*     mode */
/*      1    v  solves   L v = v(input).   w  is not touched. */
/*      2    v  solves   L'v = v(input).   w  is not touched. */
/*      3    w  solves   U w = v.          v  is not altered. */
/*      4    v  solves   U'v = w.          w  is destroyed. */
/*      5    w  solves   A w = v.          v  is altered as in 1. */
/*      6    v  solves   A'v = w.          w  is destroyed. */

/*     If mode = 3,4,5,6, v and w must not be the same arrays. */

/*     If lu1fac has just been used to factorize a symmetric matrix A */
/*     (which must be definite or quasi-definite), the factors A = L U */
/*     may be regarded as A = LDL', where D = diag(U).  In such cases, */

/*     mode */
/*      7    v  solves   A v = L D L'v = v(input).   w  is not touched. */
/*      8    v  solves       L |D| L'v = v(input).   w  is not touched. */

/*     ip(*), iq(*)      hold row and column numbers in pivotal order. */
/*     lenc(k)           is the length of the k-th column of initial L. */
/*     lenr(i)           is the length of the i-th row of U. */
/*     locc(*)           is not used. */
/*     locr(i)           is the start  of the i-th row of U. */

/*     U is assumed to be in upper-trapezoidal form (nrank by n). */
/*     The first entry for each row is the diagonal element */
/*     (according to the permutations  ip, iq).  It is stored at */
/*     location locr(i) in a(*), indr(*). */

/*     On exit, inform = 0 except as follows. */
/*     If mode = 3,4,5,6 and if U (and hence A) is singular, then */
/*     inform = 1 if there is a nonzero residual in solving the system */
/*     involving U.  parmlu(20) returns the norm of the residual. */

/*       July 1987: Early version. */
/*     09 May 1988: f77 version. */
/*     27 Apr 2000: Abolished the dreaded "computed go to". */
/*                  But hard to change other "go to"s to "if then else". */
/*     15 Dec 2002: lu6L, lu6Lt, lu6U, lu6Ut added to modularize lu6sol. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --v;
    --locc;
    --lenc;
    --iq;
    --w;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    if (*mode == 1) {
/* Solve  L v(new) = v. */
	lu6l_(inform__, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &lenc[1]);
    } else if (*mode == 2) {
/* Solve  L'v(new) = v. */
	lu6lt_(inform__, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &lenc[1]);
    } else if (*mode == 3) {
/* Solve  U w = v. */
	lu6u_(inform__, m, n, &v[1], &w[1], lena, &luparm[1], &parmlu[1], &a[
		1], &indr[1], &ip[1], &iq[1], &lenr[1], &locr[1]);
    } else if (*mode == 4) {
/* Solve  U'v = w. */
	lu6ut_(inform__, m, n, &v[1], &w[1], lena, &luparm[1], &parmlu[1], &a[
		1], &indr[1], &ip[1], &iq[1], &lenr[1], &locr[1]);
    } else if (*mode == 5) {
/* Solve  A w      = v */
	lu6l_(inform__, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &lenc[1]);
/* via    L v(new) = v */
	lu6u_(inform__, m, n, &v[1], &w[1], lena, &luparm[1], &parmlu[1], &a[
		1], &indr[1], &ip[1], &iq[1], &lenr[1], &locr[1]);
/* and    U w = v(new). */
    } else if (*mode == 6) {
/* Solve  A'v = w */
	lu6ut_(inform__, m, n, &v[1], &w[1], lena, &luparm[1], &parmlu[1], &a[
		1], &indr[1], &ip[1], &iq[1], &lenr[1], &locr[1]);
/* via    U'v = w */
	lu6lt_(inform__, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &lenc[1]);
/* and    L'v(new) = v. */
    } else if (*mode == 7) {
	lu6ld_(inform__, &c__1, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[
		1], &indc[1], &indr[1], &lenc[1], &locr[1]);
/* Solve  LDv(bar) = v */
	lu6lt_(inform__, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &lenc[1]);
/* and    L'v(new) = v(bar). */
    } else if (*mode == 8) {
	lu6ld_(inform__, &c__2, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[
		1], &indc[1], &indr[1], &lenc[1], &locr[1]);
/* Solve  L|D|v(bar) = v */
	lu6lt_(inform__, m, n, &v[1], lena, &luparm[1], &parmlu[1], &a[1], &
		indc[1], &indr[1], &lenc[1]);
/* and    L'v(new) = v(bar). */
    }
    return 0;
} /* lu6sol_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6sol */
/* Subroutine */ int lu6l_(integer *inform__, integer *m, integer *n, 
	doublereal *v, integer *lena, integer *luparm, doublereal *parmlu, 
	doublereal *a, integer *indc, integer *indr, integer *lenc)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, l1, len, lenl, ipiv, numl;
    static doublereal vpiv;
    static integer lenl0, numl0;
    static doublereal small;
    static integer ldummy;

/*     ------------------------------------------------------------------ */
/*     lu6L   solves   L v = v(input). */

/*     15 Dec 2002: First version derived from lu6sol. */
/*     15 Dec 2002: Current version. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --v;
    --lenc;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    numl0 = luparm[20];
    lenl0 = luparm[21];
    lenl = luparm[23];
    small = parmlu[3];
    *inform__ = 0;
    l1 = *lena + 1;
    i__1 = numl0;
    for (k = 1; k <= i__1; ++k) {
	len = lenc[k];
	l = l1;
	l1 -= len;
	ipiv = indr[l1];
	vpiv = v[ipiv];
	if (abs(vpiv) > small) {
/* ***** This loop could be coded specially. */
	    i__2 = len;
	    for (ldummy = 1; ldummy <= i__2; ++ldummy) {
		--l;
		j = indc[l];
		v[j] += a[l] * vpiv;
	    }
	}
    }
    l = *lena - lenl0 + 1;
    numl = lenl - lenl0;
/* ***** This loop could be coded specially. */
    i__1 = numl;
    for (ldummy = 1; ldummy <= i__1; ++ldummy) {
	--l;
	i__ = indr[l];
	if ((d__1 = v[i__], abs(d__1)) > small) {
	    j = indc[l];
	    v[j] += a[l] * v[i__];
	}
    }
/*     Exit. */
    luparm[10] = *inform__;
    return 0;
} /* lu6l_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6L */
/* Subroutine */ int lu6lt_(integer *inform__, integer *m, integer *n, 
	doublereal *v, integer *lena, integer *luparm, doublereal *parmlu, 
	doublereal *a, integer *indc, integer *indr, integer *lenc)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, l1, l2, len;
    static doublereal sum;
    static integer lenl, ipiv, lenl0, numl0;
    static doublereal small;

/*     ------------------------------------------------------------------ */
/*     lu6Lt  solves   L'v = v(input). */

/*     15 Dec 2002: First version derived from lu6sol. */
/*     15 Dec 2002: Current version. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --v;
    --lenc;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    numl0 = luparm[20];
    lenl0 = luparm[21];
    lenl = luparm[23];
    small = parmlu[3];
    *inform__ = 0;
    l1 = *lena - lenl + 1;
    l2 = *lena - lenl0;
/* ***** This loop could be coded specially. */
    i__1 = l2;
    for (l = l1; l <= i__1; ++l) {
	j = indc[l];
	if ((d__1 = v[j], abs(d__1)) > small) {
	    i__ = indr[l];
	    v[i__] += a[l] * v[j];
	}
    }
    for (k = numl0; k >= 1; --k) {
	len = lenc[k];
	sum = 0.;
	l1 = l2 + 1;
	l2 += len;
/* ***** This loop could be coded specially. */
	i__1 = l2;
	for (l = l1; l <= i__1; ++l) {
	    j = indc[l];
	    sum += a[l] * v[j];
	}
	ipiv = indr[l1];
	v[ipiv] += sum;
    }
/*     Exit. */
    luparm[10] = *inform__;
    return 0;
} /* lu6lt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6Lt */
/* Subroutine */ int lu6u_(integer *inform__, integer *m, integer *n, 
	doublereal *v, doublereal *w, integer *lena, integer *luparm, 
	doublereal *parmlu, doublereal *a, integer *indr, integer *ip, 
	integer *iq, integer *lenr, integer *locr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer l1, l2, l3;
    static doublereal resid, small;
    static integer nrank, klast, nrank1;

/*     ------------------------------------------------------------------ */
/*     lu6U   solves   U w = v.          v  is not altered. */

/*     15 Dec 2002: First version derived from lu6sol. */
/*     15 Dec 2002: Current version. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --v;
    --iq;
    --w;
    --indr;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    nrank = luparm[16];
    small = parmlu[3];
    *inform__ = 0;
    nrank1 = nrank + 1;
    resid = 0.;
/*     Find the first nonzero in v(1:nrank), counting backwards. */
    for (klast = nrank; klast >= 1; --klast) {
	i__ = ip[klast];
	if ((d__1 = v[i__], abs(d__1)) > small) {
	    goto L320;
	}
    }
L320:
    i__1 = *n;
    for (k = klast + 1; k <= i__1; ++k) {
	j = iq[k];
	w[j] = 0.;
    }
/*     Do the back-substitution, using rows 1:klast of U. */
    for (k = klast; k >= 1; --k) {
	i__ = ip[k];
	t = v[i__];
	l1 = locr[i__];
	l2 = l1 + 1;
	l3 = l1 + lenr[i__] - 1;
/* ***** This loop could be coded specially. */
	i__1 = l3;
	for (l = l2; l <= i__1; ++l) {
	    j = indr[l];
	    t -= a[l] * w[j];
	}
	j = iq[k];
	if (abs(t) <= small) {
	    w[j] = 0.;
	} else {
	    w[j] = t / a[l1];
	}
    }
/*     Compute residual for overdetermined systems. */
    i__1 = *m;
    for (k = nrank1; k <= i__1; ++k) {
	i__ = ip[k];
	resid += (d__1 = v[i__], abs(d__1));
    }
/*     Exit. */
    if (resid > 0.) {
	*inform__ = 1;
    }
    luparm[10] = *inform__;
    parmlu[20] = resid;
    return 0;
} /* lu6u_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6U */
/* Subroutine */ int lu6ut_(integer *inform__, integer *m, integer *n, 
	doublereal *v, doublereal *w, integer *lena, integer *luparm, 
	doublereal *parmlu, doublereal *a, integer *indr, integer *ip, 
	integer *iq, integer *lenr, integer *locr)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer l1, l2;
    static doublereal resid, small;
    static integer nrank, nrank1;

/*     ------------------------------------------------------------------ */
/*     lu6Ut  solves   U'v = w.          w  is destroyed. */

/*     15 Dec 2002: First version derived from lu6sol. */
/*     15 Dec 2002: Current version. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --v;
    --iq;
    --w;
    --indr;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    nrank = luparm[16];
    small = parmlu[3];
    *inform__ = 0;
    nrank1 = nrank + 1;
    resid = 0.;
    i__1 = *m;
    for (k = nrank1; k <= i__1; ++k) {
	i__ = ip[k];
	v[i__] = 0.;
    }
/*     Do the forward-substitution, skipping columns of U(transpose) */
/*     when the associated element of w(*) is negligible. */
    i__1 = nrank;
    for (k = 1; k <= i__1; ++k) {
	i__ = ip[k];
	j = iq[k];
	t = w[j];
	if (abs(t) <= small) {
	    v[i__] = 0.;
	    goto L480;
	}
	l1 = locr[i__];
	t /= a[l1];
	v[i__] = t;
	l2 = l1 + lenr[i__] - 1;
	++l1;
/* ***** This loop could be coded specially. */
	i__2 = l2;
	for (l = l1; l <= i__2; ++l) {
	    j = indr[l];
	    w[j] -= t * a[l];
	}
L480:
	;
    }
/*     Compute residual for overdetermined systems. */
    i__1 = *n;
    for (k = nrank1; k <= i__1; ++k) {
	j = iq[k];
	resid += (d__1 = w[j], abs(d__1));
    }
/*     Exit. */
    if (resid > 0.) {
	*inform__ = 1;
    }
    luparm[10] = *inform__;
    parmlu[20] = resid;
    return 0;
} /* lu6ut_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6Ut */
/* Subroutine */ int lu6ld_(integer *inform__, integer *mode, integer *m, 
	integer *n, doublereal *v, integer *lena, integer *luparm, doublereal 
	*parmlu, doublereal *a, integer *indc, integer *indr, integer *lenc, 
	integer *locr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, l, l1, len;
    static doublereal diag;
    static integer ipiv;
    static doublereal vpiv;
    static integer numl0;
    static doublereal small;
    static integer ldummy;

/* ----------------------------------------------------------------------- */
/*     lu6LD  assumes lu1fac has computed factors A = LU of a */
/*     symmetric definite or quasi-definite matrix A, */
/*     using Threshold Symmetric Pivoting (TSP),   luparm(6) = 3, */
/*     or    Threshold Diagonal  Pivoting (TDP),   luparm(6) = 4. */
/*     It also assumes that no updates have been performed. */
/*     In such cases,  U = D L', where D = diag(U). */
/*     lu6LDL returns v as follows: */

/*     mode */
/*      1    v  solves   L D v = v(input). */
/*      2    v  solves   L|D|v = v(input). */

/*     15 Dec 2002: First version of lu6LD. */
/*     15 Dec 2002: Current version. */
/* ----------------------------------------------------------------------- */
/* Solve L D v(new) = v  or  L|D|v(new) = v, depending on mode. */
/* The code for L is the same as in lu6L, */
/* but when a nonzero entry of v arises, we divide by */
/* the corresponding entry of D or |D|. */
    /* Parameter adjustments */
    --locr;
    --v;
    --lenc;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    numl0 = luparm[20];
    small = parmlu[3];
    *inform__ = 0;
    l1 = *lena + 1;
    i__1 = numl0;
    for (k = 1; k <= i__1; ++k) {
	len = lenc[k];
	l = l1;
	l1 -= len;
	ipiv = indr[l1];
	vpiv = v[ipiv];
	if (abs(vpiv) > small) {
/* ***** This loop could be coded specially. */
	    i__2 = len;
	    for (ldummy = 1; ldummy <= i__2; ++ldummy) {
		--l;
		j = indc[l];
		v[j] += a[l] * vpiv;
	    }
/* Find diag = U(ipiv,ipiv) and divide by diag or |diag|. */
	    l = locr[ipiv];
	    diag = a[l];
	    if (*mode == 2) {
		diag = abs(diag);
	    }
	    v[ipiv] = vpiv / diag;
	}
    }
    return 0;
} /* lu6ld_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6LD */
/* Subroutine */ int lu6chk_(integer *mode, integer *m, integer *n, 
	doublereal *w, integer *lena, integer *luparm, doublereal *parmlu, 
	doublereal *a, integer *indc, integer *indr, integer *ip, integer *iq,
	 integer *lenc, integer *lenr, integer *locc, integer *locr, integer *
	inform__)
{
    /* Format strings */
    static char fmt_1100[] = "(\002 Singular(m\002,a,\002n)\002,\002  ran"
	    "k\002,i9,\002  n-rank\002,i8,\002  nsing\002,i9)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, l, l1, l2;
    static doublereal aij;
    static logical trp;
    static doublereal diag;
    static integer lenl;
    static doublereal lmax, umax;
    static integer nout;
    static doublereal utol1, utol2;
    static integer nrank, jsing;
    static doublereal dumin;
    static integer nsing;
    static doublereal dumax;
    static integer jumin;
    static char mnkey[1];
    static integer ndefic, ldiagu;
    static logical keeplu;
    static integer lprint;

    /* Fortran I/O blocks */
    static cilist io___537 = { 0, 0, 0, fmt_1100, 0 };


/*     ------------------------------------------------------------------ */
/*     lu6chk  looks at the LU factorization  A = L*U. */

/*     If mode = 1, lu6chk is being called by lu1fac. */
/*     (Other modes not yet implemented.) */
/*     The important input parameters are */

/*                    lprint = luparm(2) */
/*                             luparm(6) = 1 if TRP */
/*                    keepLU = luparm(8) */
/*                    Utol1  = parmlu(4) */
/*                    Utol2  = parmlu(5) */

/*     and the significant output parameters are */

/*                    inform = luparm(10) */
/*                    nsing  = luparm(11) */
/*                    jsing  = luparm(12) */
/*                    jumin  = luparm(19) */
/*                    Lmax   = parmlu(11) */
/*                    Umax   = parmlu(12) */
/*                    DUmax  = parmlu(13) */
/*                    DUmin  = parmlu(14) */
/*                    and      w(*). */

/*     Lmax  and Umax  return the largest elements in L and U. */
/*     DUmax and DUmin return the largest and smallest diagonals of U */
/*                     (excluding diagonals that are exactly zero). */

/*     In general, w(j) is set to the maximum absolute element in */
/*     the j-th column of U.  However, if the corresponding diagonal */
/*     of U is small in absolute terms or relative to w(j) */
/*     (as judged by the parameters Utol1, Utol2 respectively), */
/*     then w(j) is changed to - w(j). */

/*     Thus, if w(j) is not positive, the j-th column of A */
/*     appears to be dependent on the other columns of A. */
/*     The number of such columns, and the position of the last one, */
/*     are returned as nsing and jsing. */

/*     Note that nrank is assumed to be set already, and is not altered. */
/*     Typically, nsing will satisfy      nrank + nsing = n,  but if */
/*     Utol1 and Utol2 are rather large,  nsing > n - nrank   may occur. */

/*     If keepLU = 0, */
/*     Lmax  and Umax  are already set by lu1fac. */
/*     The diagonals of U are in the top of A. */
/*     Only Utol1 is used in the singularity test to set w(*). */

/*     inform = 0  if  A  appears to have full column rank  (nsing = 0). */
/*     inform = 1  otherwise  (nsing .gt. 0). */

/*     00 Jul 1987: Early version. */
/*     09 May 1988: f77 version. */
/*     11 Mar 2001: Allow for keepLU = 0. */
/*     17 Nov 2001: Briefer output for singular factors. */
/*     05 May 2002: Comma needed in format 1100 (via Kenneth Holmstrom). */
/*     06 May 2002: With keepLU = 0, diags of U are in natural order. */
/*                  They were not being extracted correctly. */
/*     23 Apr 2004: TRP can judge singularity better by comparing */
/*                  all diagonals to DUmax. */
/*     27 Jun 2004: (PEG) Allow write only if nout .gt. 0 . */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --locc;
    --lenc;
    --iq;
    --w;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    nout = luparm[1];
    lprint = luparm[2];
    trp = luparm[6] == 1;
/* Threshold Rook Pivoting */
    keeplu = luparm[8] != 0;
    nrank = luparm[16];
    lenl = luparm[23];
    utol1 = parmlu[4];
    utol2 = parmlu[5];
    *inform__ = 0;
    lmax = 0.;
    umax = 0.;
    nsing = 0;
    jsing = 0;
    jumin = 0;
    dumax = 0.;
    dumin = 1e30;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	w[j] = 0.;
    }
    if (keeplu) {
/* -------------------------------------------------------------- */
/* Find  Lmax. */
/* -------------------------------------------------------------- */
	i__1 = *lena;
	for (l = *lena + 1 - lenl; l <= i__1; ++l) {
/* Computing MAX */
	    d__2 = lmax, d__3 = (d__1 = a[l], abs(d__1));
	    lmax = max(d__2,d__3);
	}
/* -------------------------------------------------------------- */
/* Find Umax and set w(j) = maximum element in j-th column of U. */
/* -------------------------------------------------------------- */
	i__1 = nrank;
	for (k = 1; k <= i__1; ++k) {
	    i__ = ip[k];
	    l1 = locr[i__];
	    l2 = l1 + lenr[i__] - 1;
	    i__2 = l2;
	    for (l = l1; l <= i__2; ++l) {
		j = indr[l];
		aij = (d__1 = a[l], abs(d__1));
/* Computing MAX */
		d__1 = w[j];
		w[j] = max(d__1,aij);
		umax = max(umax,aij);
	    }
	}
	parmlu[11] = lmax;
	parmlu[12] = umax;
/* -------------------------------------------------------------- */
/* Find DUmax and DUmin, the extreme diagonals of U. */
/* -------------------------------------------------------------- */
	i__1 = nrank;
	for (k = 1; k <= i__1; ++k) {
	    j = iq[k];
	    i__ = ip[k];
	    l1 = locr[i__];
	    diag = (d__1 = a[l1], abs(d__1));
	    dumax = max(dumax,diag);
	    if (dumin > diag) {
		dumin = diag;
		jumin = j;
	    }
	}
    } else {
/* -------------------------------------------------------------- */
/* keepLU = 0. */
/* Only diag(U) is stored.  Set w(*) accordingly. */
/* Find DUmax and DUmin, the extreme diagonals of U. */
/* -------------------------------------------------------------- */
	ldiagu = *lena - *n;
	i__1 = nrank;
	for (k = 1; k <= i__1; ++k) {
	    j = iq[k];
/* !diag   = abs( a(ldiagU + k) ) ! 06 May 2002: Diags */
	    diag = (d__1 = a[ldiagu + j], abs(d__1));
/* are in natural order */
	    w[j] = diag;
	    dumax = max(dumax,diag);
	    if (dumin > diag) {
		dumin = diag;
		jumin = j;
	    }
	}
    }
/* -------------------------------------------------------------- */
/* Negate w(j) if the corresponding diagonal of U is */
/* too small in absolute terms or relative to the other elements */
/* in the same column of  U. */

/* 23 Apr 2004: TRP ensures that diags are NOT small relative to */
/*              other elements in their own column. */
/*              Much better, we can compare all diags to DUmax. */
/* -------------------------------------------------------------- */
    if (*mode == 1 && trp) {
/* Computing MAX */
	d__1 = utol1, d__2 = utol2 * dumax;
	utol1 = max(d__1,d__2);
    }
    if (keeplu) {
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    j = iq[k];
	    if (k > nrank) {
		diag = 0.;
	    } else {
		i__ = ip[k];
		l1 = locr[i__];
		diag = (d__1 = a[l1], abs(d__1));
	    }
	    if (diag <= utol1 || diag <= utol2 * w[j]) {
		++nsing;
		jsing = j;
		w[j] = -w[j];
	    }
	}
    } else {
/* keepLU = 0 */
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    j = iq[k];
	    diag = w[j];
	    if (diag <= utol1) {
		++nsing;
		jsing = j;
		w[j] = -w[j];
	    }
	}
    }
/* ----------------------------------------------------------------- */
/* Set output parameters. */
/* ----------------------------------------------------------------- */
    if (jumin == 0) {
	dumin = 0.;
    }
    luparm[11] = nsing;
    luparm[12] = jsing;
    luparm[19] = jumin;
    parmlu[13] = dumax;
    parmlu[14] = dumin;
    if (nsing > 0) {
/* The matrix has been judged singular. */
	*inform__ = 1;
	ndefic = *n - nrank;
	if (nout > 0 && lprint >= 0) {
	    if (*m > *n) {
		*(unsigned char *)mnkey = '>';
	    } else if (*m == *n) {
		*(unsigned char *)mnkey = '=';
	    } else {
		*(unsigned char *)mnkey = '<';
	    }
	    io___537.ciunit = nout;
	    s_wsfe(&io___537);
	    do_fio(&c__1, mnkey, (ftnlen)1);
	    do_fio(&c__1, (char *)&nrank, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndefic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nsing, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
/*     Exit. */
    luparm[10] = *inform__;
    return 0;
} /* lu6chk_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* *********************************************************************** */

/*     File  lusol7a.f */

/*     lu7add   lu7cyc   lu7elm   lu7for   lu7rnk   lu7zap */

/*     Utilities for LUSOL's update routines. */
/*     lu7for is the most important -- the forward sweep. */

/* 01 May 2002: Derived from LUSOL's original lu7a.f file. */
/* 01 May 2002: Current version of lusol7a.f. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu6chk */
/* Subroutine */ int lu7add_(integer *m, integer *n, integer *jadd, 
	doublereal *v, integer *lena, integer *luparm, doublereal *parmlu, 
	integer *lenl, integer *lenu, integer *lrow, integer *nrank, 
	doublereal *a, integer *indr, integer *ip, integer *lenr, integer *
	locr, integer *inform__, integer *klast, doublereal *vnorm)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, lr1, lr2, leni, nfree;
    static doublereal small;
    extern /* Subroutine */ int lu1rec_(integer *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *);
    static integer minfre;

/*     ------------------------------------------------------------------ */
/*     lu7add  inserts the first nrank elements of the vector v(*) */
/*     as column  jadd  of  U.  We assume that  U  does not yet have any */
/*     entries in this column. */
/*     Elements no larger than  parmlu(3)  are treated as zero. */
/*     klast  will be set so that the last row to be affected */
/*     (in pivotal order) is row  ip(klast). */

/*     09 May 1988: First f77 version. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --v;
    --indr;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    small = parmlu[3];
    *vnorm = 0.;
    *klast = 0;
    i__1 = *nrank;
    for (k = 1; k <= i__1; ++k) {
	i__ = ip[k];
	if ((d__1 = v[i__], abs(d__1)) <= small) {
	    goto L200;
	}
	*klast = k;
	*vnorm += (d__1 = v[i__], abs(d__1));
	leni = lenr[i__];
/*        Compress row file if necessary. */
	minfre = leni + 1;
	nfree = *lena - *lenl - *lrow;
	if (nfree < minfre) {
	    lu1rec_(m, &c_true, &luparm[1], lrow, lena, &a[1], &indr[1], &
		    lenr[1], &locr[1]);
	    nfree = *lena - *lenl - *lrow;
	    if (nfree < minfre) {
		goto L970;
	    }
	}
/*        Move row  i  to the end of the row file, */
/*        unless it is already there. */
/*        No need to move if there is a gap already. */
	if (leni == 0) {
	    locr[i__] = *lrow + 1;
	}
	lr1 = locr[i__];
	lr2 = lr1 + leni - 1;
	if (lr2 == *lrow) {
	    goto L150;
	}
	if (indr[lr2 + 1] == 0) {
	    goto L180;
	}
	locr[i__] = *lrow + 1;
	i__2 = lr2;
	for (l = lr1; l <= i__2; ++l) {
	    ++(*lrow);
	    a[*lrow] = a[l];
	    j = indr[l];
	    indr[l] = 0;
	    indr[*lrow] = j;
/* L140: */
	}
L150:
	lr2 = *lrow;
	++(*lrow);
/*        Add the element of  v. */
L180:
	++lr2;
	a[lr2] = v[i__];
	indr[lr2] = *jadd;
	lenr[i__] = leni + 1;
	++(*lenu);
L200:
	;
    }
/*     Normal exit. */
    *inform__ = 0;
    goto L990;
/*     Not enough storage. */
L970:
    *inform__ = 7;
L990:
    return 0;
} /* lu7add_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu7add */
/* Subroutine */ int lu7cyc_(integer *kfirst, integer *klast, integer *ip)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, ifirst;

/*     ------------------------------------------------------------------ */
/*     lu7cyc performs a cyclic permutation on the row or column ordering */
/*     stored in ip, moving entry kfirst down to klast. */
/*     If kfirst .ge. klast, lu7cyc should not be called. */
/*     Sometimes klast = 0 and nothing should happen. */

/*     09 May 1988: First f77 version. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ip;

    /* Function Body */
    if (*kfirst < *klast) {
	ifirst = ip[*kfirst];
	i__1 = *klast - 1;
	for (k = *kfirst; k <= i__1; ++k) {
	    ip[k] = ip[k + 1];
/* L100: */
	}
	ip[*klast] = ifirst;
    }
    return 0;
} /* lu7cyc_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu7cyc */
/* Subroutine */ int lu7elm_(integer *m, integer *n, integer *jelm, 
	doublereal *v, integer *lena, integer *luparm, doublereal *parmlu, 
	integer *lenl, integer *lenu, integer *lrow, integer *nrank, 
	doublereal *a, integer *indc, integer *indr, integer *ip, integer *iq,
	 integer *lenr, integer *locc, integer *locr, integer *inform__, 
	doublereal *diag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l, l1, l2;
    static doublereal vi;
    static integer imax, kmax, lmax;
    static doublereal vmax;
    static integer nfree;
    static doublereal small;
    static integer nrank1;
    extern /* Subroutine */ int lu1rec_(integer *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *);
    static integer minfre;

/*     ------------------------------------------------------------------ */
/*     lu7elm  eliminates the subdiagonal elements of a vector  v(*), */
/*     where  L*v = y  for some vector y. */
/*     If  jelm > 0,  y  has just become column  jelm  of the matrix  A. */
/*     lu7elm  should not be called unless  m  is greater than  nrank. */

/*     inform = 0 if y contained no subdiagonal nonzeros to eliminate. */
/*     inform = 1 if y contained at least one nontrivial subdiagonal. */
/*     inform = 7 if there is insufficient storage. */

/*     09 May 1988: First f77 version. */
/*                  No longer calls lu7for at end.  lu8rpc, lu8mod do so. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --v;
    --locc;
    --iq;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    small = parmlu[3];
    nrank1 = *nrank + 1;
    *diag = 0.;
/*     Compress row file if necessary. */
    minfre = *m - *nrank;
    nfree = *lena - *lenl - *lrow;
    if (nfree >= minfre) {
	goto L100;
    }
    lu1rec_(m, &c_true, &luparm[1], lrow, lena, &a[1], &indr[1], &lenr[1], &
	    locr[1]);
    nfree = *lena - *lenl - *lrow;
    if (nfree < minfre) {
	goto L970;
    }
/*     Pack the subdiagonals of  v  into  L,  and find the largest. */
L100:
    vmax = 0.;
    kmax = 0;
    l = *lena - *lenl + 1;
    i__1 = *m;
    for (k = nrank1; k <= i__1; ++k) {
	i__ = ip[k];
	vi = (d__1 = v[i__], abs(d__1));
	if (vi <= small) {
	    goto L200;
	}
	--l;
	a[l] = v[i__];
	indc[l] = i__;
	if (vmax >= vi) {
	    goto L200;
	}
	vmax = vi;
	kmax = k;
	lmax = l;
L200:
	;
    }
    if (kmax == 0) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Remove  vmax  by overwriting it with the last packed  v(i). */
/*     Then set the multipliers in  L  for the other elements. */
/*     ------------------------------------------------------------------ */
    imax = ip[kmax];
    vmax = a[lmax];
    a[lmax] = a[l];
    indc[lmax] = indc[l];
    l1 = l + 1;
    l2 = *lena - *lenl;
    *lenl += l2 - l;
    i__1 = l2;
    for (l = l1; l <= i__1; ++l) {
	a[l] = -a[l] / vmax;
	indr[l] = imax;
/* L300: */
    }
/*     Move the row containing vmax to pivotal position nrank + 1. */
    ip[kmax] = ip[nrank1];
    ip[nrank1] = imax;
    *diag = vmax;
/*     ------------------------------------------------------------------ */
/*     If jelm is positive, insert  vmax  into a new row of  U. */
/*     This is now the only subdiagonal element. */
/*     ------------------------------------------------------------------ */
    if (*jelm > 0) {
	++(*lrow);
	locr[imax] = *lrow;
	lenr[imax] = 1;
	a[*lrow] = vmax;
	indr[*lrow] = *jelm;
    }
    *inform__ = 1;
    goto L990;
/*     No elements to eliminate. */
L900:
    *inform__ = 0;
    goto L990;
/*     Not enough storage. */
L970:
    *inform__ = 7;
L990:
    return 0;
} /* lu7elm_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu7elm */
/* Subroutine */ int lu7for_(integer *m, integer *n, integer *kfirst, integer 
	*klast, integer *lena, integer *luparm, doublereal *parmlu, integer *
	lenl, integer *lenu, integer *lrow, doublereal *a, integer *indc, 
	integer *indr, integer *ip, integer *iq, integer *lenr, integer *locc,
	 integer *locr, integer *inform__, doublereal *diag)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer j, k, l, iv, iw;
    static doublereal wj, vj;
    static integer lv, jv, lw, lv1, lw1, lw2, lv2, lv3, lenv, lenw;
    static doublereal ltol;
    static integer ldiag, nfree;
    static doublereal small;
    static integer jlast, limit;
    static doublereal amult;
    static integer kstop;
    extern /* Subroutine */ int lu1rec_(integer *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *);
    static integer kbegin;
    static doublereal uspace;
    static integer minfre;
    static logical swappd;
    static integer jfirst, lfirst, kstart;

/*     ------------------------------------------------------------------ */
/*     lu7for  (forward sweep) updates the LU factorization  A = L*U */
/*     when row  iw = ip(klast)  of  U  is eliminated by a forward */
/*     sweep of stabilized row operations, leaving  ip * U * iq  upper */
/*     triangular. */

/*     The row permutation  ip  is updated to preserve stability and/or */
/*     sparsity.  The column permutation  iq  is not altered. */

/*     kfirst  is such that row  ip(kfirst)  is the first row involved */
/*     in eliminating row  iw.  (Hence,  kfirst  marks the first nonzero */
/*     in row  iw  in pivotal order.)  If  kfirst  is unknown it may be */
/*     input as  1. */

/*     klast   is such that row  ip(klast)  is the row being eliminated. */
/*     klast   is not altered. */

/*     lu7for  should be called only if  kfirst .le. klast. */
/*     If  kfirst = klast,  there are no nonzeros to eliminate, but the */
/*     diagonal element of row  ip(klast)  may need to be moved to the */
/*     front of the row. */

/*     On entry,  locc(*)  must be zero. */

/*     On exit: */
/*     inform = 0  if row iw has a nonzero diagonal (could be small). */
/*     inform = 1  if row iw has no diagonal. */
/*     inform = 7  if there is not enough storage to finish the update. */

/*     On a successful exit (inform le 1),  locc(*)  will again be zero. */

/*        Jan 1985: Final f66 version. */
/*     09 May 1988: First f77 version. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --locc;
    --iq;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    ltol = parmlu[2];
    small = parmlu[3];
    uspace = parmlu[6];
    kbegin = *kfirst;
    swappd = FALSE_;
/*     We come back here from below if a row interchange is performed. */
L100:
    iw = ip[*klast];
    lenw = lenr[iw];
    if (lenw == 0) {
	goto L910;
    }
    lw1 = locr[iw];
    lw2 = lw1 + lenw - 1;
    jfirst = iq[kbegin];
    if (kbegin >= *klast) {
	goto L700;
    }
/*     Make sure there is room at the end of the row file */
/*     in case row  iw  is moved there and fills in completely. */
    minfre = *n + 1;
    nfree = *lena - *lenl - *lrow;
    if (nfree < minfre) {
	lu1rec_(m, &c_true, &luparm[1], lrow, lena, &a[1], &indr[1], &lenr[1],
		 &locr[1]);
	lw1 = locr[iw];
	lw2 = lw1 + lenw - 1;
	nfree = *lena - *lenl - *lrow;
	if (nfree < minfre) {
	    goto L970;
	}
    }
/*     Set markers on row  iw. */
    i__1 = lw2;
    for (l = lw1; l <= i__1; ++l) {
	j = indr[l];
	locc[j] = l;
/* L120: */
    }
/*     ================================================================== */
/*     Main elimination loop. */
/*     ================================================================== */
    kstart = kbegin;
    kstop = min(*klast,*n);
    i__1 = kstop;
    for (k = kstart; k <= i__1; ++k) {
	jfirst = iq[k];
	lfirst = locc[jfirst];
	if (lfirst == 0) {
	    goto L490;
	}
/*        Row  iw  has its first element in column  jfirst. */
	wj = a[lfirst];
	if (k == *klast) {
	    goto L490;
	}
/*        --------------------------------------------------------------- */
/*        We are about to use the first element of row  iv */
/*               to eliminate the first element of row  iw. */
/*        However, we may wish to interchange the rows instead, */
/*        to preserve stability and/or sparsity. */
/*        --------------------------------------------------------------- */
	iv = ip[k];
	lenv = lenr[iv];
	lv1 = locr[iv];
	vj = 0.;
	if (lenv == 0) {
	    goto L150;
	}
	if (indr[lv1] != jfirst) {
	    goto L150;
	}
	vj = a[lv1];
	if (swappd) {
	    goto L200;
	}
	if (ltol * abs(wj) < abs(vj)) {
	    goto L200;
	}
	if (ltol * abs(vj) < abs(wj)) {
	    goto L150;
	}
	if (lenv <= lenw) {
	    goto L200;
	}
/*        --------------------------------------------------------------- */
/*        Interchange rows  iv  and  iw. */
/*        --------------------------------------------------------------- */
L150:
	ip[*klast] = iv;
	ip[k] = iw;
	kbegin = k;
	swappd = TRUE_;
	goto L600;
/*        --------------------------------------------------------------- */
/*        Delete the eliminated element from row  iw */
/*        by overwriting it with the last element. */
/*        --------------------------------------------------------------- */
L200:
	a[lfirst] = a[lw2];
	jlast = indr[lw2];
	indr[lfirst] = jlast;
	indr[lw2] = 0;
	locc[jlast] = lfirst;
	locc[jfirst] = 0;
	--lenw;
	--(*lenu);
	if (*lrow == lw2) {
	    --(*lrow);
	}
	--lw2;
/*        --------------------------------------------------------------- */
/*        Form the multiplier and store it in the  L  file. */
/*        --------------------------------------------------------------- */
	if (abs(wj) <= small) {
	    goto L490;
	}
	amult = -wj / vj;
	l = *lena - *lenl;
	a[l] = amult;
	indr[l] = iv;
	indc[l] = iw;
	++(*lenl);
/*        --------------------------------------------------------------- */
/*        Add the appropriate multiple of row  iv  to row  iw. */
/*        We use two different inner loops.  The first one is for the */
/*        case where row  iw  is not at the end of storage. */
/*        --------------------------------------------------------------- */
	if (lenv == 1) {
	    goto L490;
	}
	lv2 = lv1 + 1;
	lv3 = lv1 + lenv - 1;
	if (lw2 == *lrow) {
	    goto L400;
	}
/*        ............................................................... */
/*        This inner loop will be interrupted only if */
/*        fill-in occurs enough to bump into the next row. */
/*        ............................................................... */
	i__2 = lv3;
	for (lv = lv2; lv <= i__2; ++lv) {
	    jv = indr[lv];
	    lw = locc[jv];
	    if (lw > 0) {
/*              No fill-in. */
		a[lw] += amult * a[lv];
		if ((d__1 = a[lw], abs(d__1)) <= small) {
/*                 Delete small computed element. */
		    a[lw] = a[lw2];
		    j = indr[lw2];
		    indr[lw] = j;
		    indr[lw2] = 0;
		    locc[j] = lw;
		    locc[jv] = 0;
		    --(*lenu);
		    --lenw;
		    --lw2;
		}
	    } else {
/*              Row  iw  doesn't have an element in column  jv  yet */
/*              so there is a fill-in. */
		if (indr[lw2 + 1] != 0) {
		    goto L360;
		}
		++(*lenu);
		++lenw;
		++lw2;
		a[lw2] = amult * a[lv];
		indr[lw2] = jv;
		locc[jv] = lw2;
	    }
/* L350: */
	}
	goto L490;
/*        Fill-in interrupted the previous loop. */
/*        Move row  iw  to the end of the row file. */
L360:
	lv2 = lv;
	locr[iw] = *lrow + 1;
	i__2 = lw2;
	for (l = lw1; l <= i__2; ++l) {
	    ++(*lrow);
	    a[*lrow] = a[l];
	    j = indr[l];
	    indr[l] = 0;
	    indr[*lrow] = j;
	    locc[j] = *lrow;
/* L370: */
	}
	lw1 = locr[iw];
	lw2 = *lrow;
/*        ............................................................... */
/*        Inner loop with row  iw  at the end of storage. */
/*        ............................................................... */
L400:
	i__2 = lv3;
	for (lv = lv2; lv <= i__2; ++lv) {
	    jv = indr[lv];
	    lw = locc[jv];
	    if (lw > 0) {
/*              No fill-in. */
		a[lw] += amult * a[lv];
		if ((d__1 = a[lw], abs(d__1)) <= small) {
/*                 Delete small computed element. */
		    a[lw] = a[lw2];
		    j = indr[lw2];
		    indr[lw] = j;
		    indr[lw2] = 0;
		    locc[j] = lw;
		    locc[jv] = 0;
		    --(*lenu);
		    --lenw;
		    --lw2;
		}
	    } else {
/*              Row  iw  doesn't have an element in column  jv  yet */
/*              so there is a fill-in. */
		++(*lenu);
		++lenw;
		++lw2;
		a[lw2] = amult * a[lv];
		indr[lw2] = jv;
		locc[jv] = lw2;
	    }
/* L450: */
	}
	*lrow = lw2;
/*        The  k-th  element of row  iw  has been processed. */
/*        Reset  swappd  before looking at the next element. */
L490:
	swappd = FALSE_;
/* L500: */
    }
/*     ================================================================== */
/*     End of main elimination loop. */
/*     ================================================================== */
/*     Cancel markers on row  iw. */
L600:
    lenr[iw] = lenw;
    if (lenw == 0) {
	goto L910;
    }
    i__1 = lw2;
    for (l = lw1; l <= i__1; ++l) {
	j = indr[l];
	locc[j] = 0;
/* L620: */
    }
/*     Move the diagonal element to the front of row  iw. */
/*     At this stage,  lenw gt 0  and  klast le n. */
L700:
    i__1 = lw2;
    for (l = lw1; l <= i__1; ++l) {
	ldiag = l;
	if (indr[l] == jfirst) {
	    goto L730;
	}
/* L720: */
    }
    goto L910;
L730:
    *diag = a[ldiag];
    a[ldiag] = a[lw1];
    a[lw1] = *diag;
    indr[ldiag] = indr[lw1];
    indr[lw1] = jfirst;
/*     If an interchange is needed, repeat from the beginning with the */
/*     new row  iw,  knowing that the opposite interchange cannot occur. */
    if (swappd) {
	goto L100;
    }
    *inform__ = 0;
    goto L950;
/*     Singular. */
L910:
    *diag = 0.;
    *inform__ = 1;
/*     Force a compression if the file for  U  is much longer than the */
/*     no. of nonzeros in  U  (i.e. if  lrow  is much bigger than  lenU). */
/*     This should prevent memory fragmentation when there is far more */
/*     memory than necessary  (i.e. when  lena  is huge). */
L950:
    limit = (integer) (uspace * *lenu + *m + *n + 1000);
    if (*lrow > limit) {
	lu1rec_(m, &c_true, &luparm[1], lrow, lena, &a[1], &indr[1], &lenr[1],
		 &locr[1]);
    }
    goto L990;
/*     Not enough storage. */
L970:
    *inform__ = 7;
/*     Exit. */
L990:
    return 0;
} /* lu7for_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu7for */
/* Subroutine */ int lu7rnk_(integer *m, integer *n, integer *jsing, integer *
	lena, integer *luparm, doublereal *parmlu, integer *lenl, integer *
	lenu, integer *lrow, integer *nrank, doublereal *a, integer *indc, 
	integer *indr, integer *ip, integer *iq, integer *lenr, integer *locc,
	 integer *locr, integer *inform__, doublereal *diag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer l, l1, l2, iw, lmax, jmax, kmax, lenw;
    static doublereal umax, utol1;

/*     ------------------------------------------------------------------ */
/*     lu7rnk (check rank) assumes U is currently nrank by n */
/*     and determines if row nrank contains an acceptable pivot. */
/*     If not, the row is deleted and nrank is decreased by 1. */

/*     jsing is an input parameter (not altered).  If jsing is positive, */
/*     column jsing has already been judged dependent.  A substitute */
/*     (if any) must be some other column. */

/*     -- Jul 1987: First version. */
/*     09 May 1988: First f77 version. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --locc;
    --iq;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    utol1 = parmlu[4];
    *diag = 0.;
/*     Find Umax, the largest element in row nrank. */
    iw = ip[*nrank];
    lenw = lenr[iw];
    if (lenw == 0) {
	goto L400;
    }
    l1 = locr[iw];
    l2 = l1 + lenw - 1;
    umax = 0.;
    lmax = l1;
    i__1 = l2;
    for (l = l1; l <= i__1; ++l) {
	if (umax < (d__1 = a[l], abs(d__1))) {
	    umax = (d__1 = a[l], abs(d__1));
	    lmax = l;
	}
/* L100: */
    }
/*     Find which column that guy is in (in pivotal order). */
/*     Interchange him with column nrank, then move him to be */
/*     the new diagonal at the front of row nrank. */
    *diag = a[lmax];
    jmax = indr[lmax];
    i__1 = *n;
    for (kmax = *nrank; kmax <= i__1; ++kmax) {
	if (iq[kmax] == jmax) {
	    goto L320;
	}
/* L300: */
    }
L320:
    iq[kmax] = iq[*nrank];
    iq[*nrank] = jmax;
    a[lmax] = a[l1];
    a[l1] = *diag;
    indr[lmax] = indr[l1];
    indr[l1] = jmax;
/*     See if the new diagonal is big enough. */
    if (umax <= utol1) {
	goto L400;
    }
    if (jmax == *jsing) {
	goto L400;
    }
/*     ------------------------------------------------------------------ */
/*     The rank stays the same. */
/*     ------------------------------------------------------------------ */
    *inform__ = 0;
    return 0;
/*     ------------------------------------------------------------------ */
/*     The rank decreases by one. */
/*     ------------------------------------------------------------------ */
L400:
    *inform__ = -1;
    --(*nrank);
    if (lenw > 0) {
/*        Delete row nrank from U. */
	*lenu -= lenw;
	lenr[iw] = 0;
	i__1 = l2;
	for (l = l1; l <= i__1; ++l) {
	    indr[l] = 0;
/* L420: */
	}
	if (l2 == *lrow) {
/*           This row was at the end of the data structure. */
/*           We have to reset lrow. */
/*           Preceding rows might already have been deleted, so we */
/*           have to be prepared to go all the way back to 1. */
	    i__1 = l2;
	    for (l = 1; l <= i__1; ++l) {
		if (indr[*lrow] > 0) {
		    goto L900;
		}
		--(*lrow);
/* L450: */
	    }
	}
    }
L900:
    return 0;
} /* lu7rnk_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu7rnk */
/* Subroutine */ int lu7zap_(integer *m, integer *n, integer *jzap, integer *
	kzap, integer *lena, integer *lenu, integer *lrow, integer *nrank, 
	doublereal *a, integer *indr, integer *ip, integer *iq, integer *lenr,
	 integer *locr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, l, lr1, lr2, leni;

/*     ------------------------------------------------------------------ */
/*     lu7zap  eliminates all nonzeros in column  jzap  of  U. */
/*     It also sets  kzap  to the position of  jzap  in pivotal order. */
/*     Thus, on exit we have  iq(kzap) = jzap. */

/*     -- Jul 1987: nrank added. */
/*     10 May 1988: First f77 version. */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --iq;
    --indr;
    --a;

    /* Function Body */
    i__1 = *nrank;
    for (k = 1; k <= i__1; ++k) {
	i__ = ip[k];
	leni = lenr[i__];
	if (leni == 0) {
	    goto L90;
	}
	lr1 = locr[i__];
	lr2 = lr1 + leni - 1;
	i__2 = lr2;
	for (l = lr1; l <= i__2; ++l) {
	    if (indr[l] == *jzap) {
		goto L60;
	    }
/* L50: */
	}
	goto L90;
/*        Delete the old element. */
L60:
	a[l] = a[lr2];
	indr[l] = indr[lr2];
	indr[lr2] = 0;
	lenr[i__] = leni - 1;
	--(*lenu);
/*        Stop if we know there are no more rows containing  jzap. */
L90:
	*kzap = k;
	if (iq[k] == *jzap) {
	    goto L800;
	}
/* L100: */
    }
/*     nrank must be smaller than n because we haven't found kzap yet. */
    i__1 = *n;
    for (k = *nrank + 1; k <= i__1; ++k) {
	*kzap = k;
	if (iq[k] == *jzap) {
	    goto L800;
	}
/* L200: */
    }
/*     See if we zapped the last element in the file. */
L800:
    if (*lrow > 0) {
	if (indr[*lrow] == 0) {
	    --(*lrow);
	}
    }
    return 0;
} /* lu7zap_ */

/* *********************************************************************** */

/*     File  lusol8a.f */

/*     lu8rpc */

/*     Sparse LU update: Replace Column */
/*     LUSOL's sparse implementation of the Bartels-Golub update. */

/* 01 May 2002: Derived from LUSOL's original lu8a.f file. */
/* 01 May 2002: Current version of lusol8a.f. */
/* 15 Sep 2004: Test nout. gt. 0 to protect write statements. */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine lu7zap */
/* Subroutine */ int lu8rpc_(integer *mode1, integer *mode2, integer *m, 
	integer *n, integer *jrep, doublereal *v, doublereal *w, integer *
	lena, integer *luparm, doublereal *parmlu, doublereal *a, integer *
	indc, integer *indr, integer *ip, integer *iq, integer *lenc, integer 
	*lenr, integer *locc, integer *locr, integer *inform__, doublereal *
	diag, doublereal *vnorm)
{
    /* Format strings */
    static char fmt_1100[] = "(/\002 lu8rpc  warning.  Singularity after rep"
	    "lacing column.\002,\002    jrep =\002,i8,\002    diag =\002,1p,e"
	    "12.2)";
    static char fmt_1200[] = "(/\002 lu8rpc  warning.  Instability after rep"
	    "lacing column.\002,\002    jrep =\002,i8,\002    diag =\002,1p,e"
	    "12.2)";
    static char fmt_1700[] = "(/\002 lu8rpc  error...  Insufficient storage"
	    ".\002,\002    lena =\002,i8)";
    static char fmt_1800[] = "(/\002 lu8rpc  error...  jrep  is out of ran"
	    "ge.\002,\002    m =\002,i8,\002    n =\002,i8,\002    jrep =\002"
	    ",i8)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer j1, l1, iw, lenl, krep, lenu, lrow, nout;
    static doublereal utol1, utol2;
    static integer nrank, jsing, klast;
    extern /* Subroutine */ int lu7add_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *);
    static integer nrank0;
    extern /* Subroutine */ int lu7elm_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *), lu7cyc_(integer *, integer *, integer *), 
	    lu7for_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *), lu7zap_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *), lu7rnk_(integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *), lu6sol_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *);
    static logical singlr;
    static integer lprint;

    /* Fortran I/O blocks */
    static cilist io___628 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___629 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___630 = { 0, 0, 0, fmt_1700, 0 };
    static cilist io___631 = { 0, 0, 0, fmt_1800, 0 };


/*     ------------------------------------------------------------------ */
/*     lu8rpc  updates the LU factorization  A = L*U  when column  jrep */
/*     is replaced by some vector  a(new). */

/*     lu8rpc  is an implementation of the Bartels-Golub update, */
/*     designed for the case where A is rectangular and/or singular. */
/*     L is a product of stabilized eliminations (m x m, nonsingular). */
/*     P U Q is upper trapezoidal (m x n, rank nrank). */

/*     If  mode1 = 0,  the old column is taken to be zero */
/*                     (so it does not have to be removed from  U). */

/*     If  mode1 = 1,  the old column need not have been zero. */

/*     If  mode2 = 0,  the new column is taken to be zero. */
/*                     v(*)  is not used or altered. */

/*     If  mode2 = 1,  v(*)  must contain the new column  a(new). */
/*                     On exit,  v(*)  will satisfy  L*v = a(new). */

/*     If  mode2 = 2,  v(*)  must satisfy  L*v = a(new). */

/*     The array  w(*)  is not used or altered. */

/*     On entry, all elements of  locc  are assumed to be zero. */
/*     On a successful exit (inform ne 7), this will again be true. */

/*     On exit: */
/*     inform = -1  if the rank of U decreased by 1. */
/*     inform =  0  if the rank of U stayed the same. */
/*     inform =  1  if the rank of U increased by 1. */
/*     inform =  2  if the update seemed to be unstable */
/*                  (diag much bigger than vnorm). */
/*     inform =  7  if the update was not completed (lack of storage). */
/*     inform =  8  if jrep is not between 1 and n. */

/*     -- Jan 1985: Original F66 version. */
/*     -- Jul 1987: Modified to maintain U in trapezoidal form. */
/*     10 May 1988: First f77 version. */
/*     16 Oct 2000: Added test for instability (inform = 2). */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --locr;
    --lenr;
    --ip;
    --v;
    --locc;
    --lenc;
    --iq;
    --w;
    --indr;
    --indc;
    --a;
    --luparm;
    --parmlu;

    /* Function Body */
    nout = luparm[1];
    lprint = luparm[2];
    nrank = luparm[16];
    lenl = luparm[23];
    lenu = luparm[24];
    lrow = luparm[25];
    utol1 = parmlu[4];
    utol2 = parmlu[5];
    nrank0 = nrank;
    *diag = 0.;
    *vnorm = 0.;
    if (*jrep < 1) {
	goto L980;
    }
    if (*jrep > *n) {
	goto L980;
    }
/*     ------------------------------------------------------------------ */
/*     If mode1 = 0, there are no elements to be removed from  U */
/*     but we still have to set  krep  (using a backward loop). */
/*     Otherwise, use lu7zap to remove column  jrep  from  U */
/*     and set  krep  at the same time. */
/*     ------------------------------------------------------------------ */
    if (*mode1 == 0) {
	krep = *n + 1;
L10:
	--krep;
	if (iq[krep] != *jrep) {
	    goto L10;
	}
    } else {
	lu7zap_(m, n, jrep, &krep, lena, &lenu, &lrow, &nrank, &a[1], &indr[1]
		, &ip[1], &iq[1], &lenr[1], &locr[1]);
    }
/*     ------------------------------------------------------------------ */
/*     Insert a new column of u and find klast. */
/*     ------------------------------------------------------------------ */
    if (*mode2 == 0) {
	klast = 0;
    } else {
	if (*mode2 == 1) {
/*           Transform v = a(new) to satisfy  L*v = a(new). */
	    lu6sol_(&c__1, m, n, &v[1], &w[1], lena, &luparm[1], &parmlu[1], &
		    a[1], &indc[1], &indr[1], &ip[1], &iq[1], &lenc[1], &lenr[
		    1], &locc[1], &locr[1], inform__);
	}
/*        Insert into  U  any nonzeros in the top of  v. */
/*        row  ip(klast)  will contain the last nonzero in pivotal order. */
/*        Note that  klast  will be in the range  ( 0, nrank ). */
	lu7add_(m, n, jrep, &v[1], lena, &luparm[1], &parmlu[1], &lenl, &lenu,
		 &lrow, &nrank, &a[1], &indr[1], &ip[1], &lenr[1], &locr[1], 
		inform__, &klast, vnorm);
	if (*inform__ == 7) {
	    goto L970;
	}
    }
/*     ------------------------------------------------------------------ */
/*     In general, the new column causes U to look like this: */

/*                 krep        n                 krep  n */

/*                ....a.........          ..........a... */
/*                 .  a        .           .        a  . */
/*                  . a        .            .       a  . */
/*                   .a        .             .      a  . */
/*        P U Q =     a        .    or        .     a  . */
/*                    b.       .               .    a  . */
/*                    b .      .                .   a  . */
/*                    b  .     .                 .  a  . */
/*                    b   ......                  ..a...  nrank */
/*                    c                             c */
/*                    c                             c */
/*                    c                             c     m */

/*     klast points to the last nonzero "a" or "b". */
/*     klast = 0 means all "a" and "b" entries are zero. */
/*     ------------------------------------------------------------------ */
    if (*mode2 == 0) {
	if (krep > nrank) {
	    goto L900;
	}
    } else if (nrank < *m) {
/*        Eliminate any "c"s (in either case). */
/*        Row nrank + 1 may end up containing one nonzero. */
	lu7elm_(m, n, jrep, &v[1], lena, &luparm[1], &parmlu[1], &lenl, &lenu,
		 &lrow, &nrank, &a[1], &indc[1], &indr[1], &ip[1], &iq[1], &
		lenr[1], &locc[1], &locr[1], inform__, diag);
	if (*inform__ == 7) {
	    goto L970;
	}
	if (*inform__ == 1) {
/*           The nonzero is apparently significant. */
/*           Increase nrank by 1 and make klast point to the bottom. */
	    ++nrank;
	    klast = nrank;
	}
    }
    if (nrank < *n) {
/*        The column rank is low. */

/*        In the first case, we want the new column to end up in */
/*        position nrank, so the trapezoidal columns will have a chance */
/*        later on (in lu7rnk) to pivot in that position. */

/*        Otherwise the new column is not part of the triangle.  We */
/*        swap it into position nrank so we can judge it for singularity. */
/*        lu7rnk might choose some other trapezoidal column later. */
	if (krep < nrank) {
	    klast = nrank;
	} else {
	    iq[krep] = iq[nrank];
	    iq[nrank] = *jrep;
	    krep = nrank;
	}
    }
/*     ------------------------------------------------------------------ */
/*     If krep .lt. klast, there are some "b"s to eliminate: */

/*                  krep */

/*                ....a......... */
/*                 .  a        . */
/*                  . a        . */
/*                   .a        . */
/*        P U Q =     a        .  krep */
/*                    b.       . */
/*                    b .      . */
/*                    b  .     . */
/*                    b   ......  nrank */

/*     If krep .eq. klast, there are no "b"s, but the last "a" still */
/*     has to be moved to the front of row krep (by lu7for). */
/*     ------------------------------------------------------------------ */
    if (krep <= klast) {
/*        Perform a cyclic permutation on the current pivotal order, */
/*        and eliminate the resulting row spike.  krep becomes klast. */
/*        The final diagonal (if any) will be correctly positioned at */
/*        the front of the new krep-th row.  nrank stays the same. */
	lu7cyc_(&krep, &klast, &ip[1]);
	lu7cyc_(&krep, &klast, &iq[1]);
	lu7for_(m, n, &krep, &klast, lena, &luparm[1], &parmlu[1], &lenl, &
		lenu, &lrow, &a[1], &indc[1], &indr[1], &ip[1], &iq[1], &lenr[
		1], &locc[1], &locr[1], inform__, diag);
	if (*inform__ == 7) {
	    goto L970;
	}
	krep = klast;
/*        Test for instability (diag much bigger than vnorm). */
	singlr = *vnorm < utol2 * abs(*diag);
	if (singlr) {
	    goto L920;
	}
    }
/*     ------------------------------------------------------------------ */
/*     Test for singularity in column krep (where krep .le. nrank). */
/*     ------------------------------------------------------------------ */
    *diag = 0.;
    iw = ip[krep];
    singlr = lenr[iw] == 0;
    if (! singlr) {
	l1 = locr[iw];
	j1 = indr[l1];
	singlr = j1 != *jrep;
	if (! singlr) {
	    *diag = a[l1];
	    singlr = abs(*diag) <= utol1 || abs(*diag) <= utol2 * *vnorm;
	}
    }
    if (singlr && krep < nrank) {
/*        Perform cyclic permutations to move column jrep to the end. */
/*        Move the corresponding row to position nrank */
/*        then eliminate the resulting row spike. */
	lu7cyc_(&krep, &nrank, &ip[1]);
	lu7cyc_(&krep, n, &iq[1]);
	lu7for_(m, n, &krep, &nrank, lena, &luparm[1], &parmlu[1], &lenl, &
		lenu, &lrow, &a[1], &indc[1], &indr[1], &ip[1], &iq[1], &lenr[
		1], &locc[1], &locr[1], inform__, diag);
	if (*inform__ == 7) {
	    goto L970;
	}
    }
/*     Find the best column to be in position nrank. */
/*     If singlr, it can't be the new column, jrep. */
/*     If nothing satisfactory exists, nrank will be decreased. */
    if (singlr || nrank < *n) {
	jsing = 0;
	if (singlr) {
	    jsing = *jrep;
	}
	lu7rnk_(m, n, &jsing, lena, &luparm[1], &parmlu[1], &lenl, &lenu, &
		lrow, &nrank, &a[1], &indc[1], &indr[1], &ip[1], &iq[1], &
		lenr[1], &locc[1], &locr[1], inform__, diag);
    }
/*     ------------------------------------------------------------------ */
/*     Set inform for exit. */
/*     ------------------------------------------------------------------ */
L900:
    if (nrank == nrank0) {
	*inform__ = 0;
    } else if (nrank < nrank0) {
	*inform__ = -1;
	if (nrank0 == *n) {
	    if (nout > 0 && lprint >= 0) {
		io___628.ciunit = nout;
		s_wsfe(&io___628);
		do_fio(&c__1, (char *)&(*jrep), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*diag), (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    } else {
	*inform__ = 1;
    }
    goto L990;
/*     Instability. */
L920:
    *inform__ = 2;
    if (nout > 0 && lprint >= 0) {
	io___629.ciunit = nout;
	s_wsfe(&io___629);
	do_fio(&c__1, (char *)&(*jrep), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*diag), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    goto L990;
/*     Not enough storage. */
L970:
    *inform__ = 7;
    if (nout > 0 && lprint >= 0) {
	io___630.ciunit = nout;
	s_wsfe(&io___630);
	do_fio(&c__1, (char *)&(*lena), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L990;
/*     jrep  is out of range. */
L980:
    *inform__ = 8;
    if (nout > 0 && lprint >= 0) {
	io___631.ciunit = nout;
	s_wsfe(&io___631);
	do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*jrep), (ftnlen)sizeof(integer));
	e_wsfe();
    }
/*     Exit. */
L990:
    luparm[10] = *inform__;
    ++luparm[15];
    luparm[16] = nrank;
    luparm[23] = lenl;
    luparm[24] = lenu;
    luparm[25] = lrow;
    return 0;
} /* lu8rpc_ */

