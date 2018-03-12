/* ./src/sn10mach.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     File  sn10mach.f                    Machine dependent routines */

/*     s1cpu                             =>  Timing routine */
/*     s1eps    s1flmx   s1flmn          =>  Floating-point arithmetic */
/*     s1intmx                           =>  largest integer */
/*     s1inpt   s1outpt                  =>  standard input/output */
/*     s1file                            =>  Default File types */
/*     s1clos   s1envt   s1open   s1page =>  Bibs n bobs */

/*     22 Jun 2004: s1cpu: Added GAMS clock gfclck(). */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Subroutine */ int s1cpu_(integer *mode, real *time)
{
    extern real etime_(real *);
    static real tarray[2];

/*     ------------------------------------------------------------------ */
/*     s1cpu is a machine-dependent routine to return time = cpu time */
/*     in seconds, so that 2 consecutive calls will indicate the */
/*     time difference of operations between the 2 calls. */
/*     The parameter 'mode' indicates what function should be done */
/*     to the timer.  This allows necessary initialization for certain */
/*     machines. */
/*     mode =  1  indicates initialization, */
/*     mode =  0  indicates normal use, */
/*     mode = -1  indicates stop the timer. */

/*     1988:  Used Irv Lustig's approach... */
/*     On DEC VAX/VMS systems we need to call the correct library */
/*     routine to get the timer statistics.  These statistics are */
/*     found by using the times() function in the VAX C Runtime library. */
/*     To use this version of s1cpu, one must create an options file */
/*     called  vmsc.opt  with the line */
/*        SYS$LIBRARY:VAXCRTL/SHARE */
/*     in it.   Then link using the usual command and append ,vmsc/opt */
/*     to the end of the line.  The name vmsc can be anything. */

/*     02 Apr 1993: Went back to VMS Fortran routines to avoid linking */
/*                  to the C library.  (On DEC AXP, the C runtime lib */
/*                  appears to be translated from the VAX executable, */
/*                  and therefore requires linking with /NONATIVE, */
/*                  which possibly adds a small overhead to all */
/*                  subroutine calls. */
/*     21 Oct 1999: Timer for WinNT with DEC F90 compiler.  Code provided */
/*                  by Thomas Kronseder. */
/*     22 Jun 2004: (At GAMS) Added GAMS clock gfclck(). */
/*                  1. Uncomment two separate lines: */
/*                        external gfclck */
/*                           ... */
/*                        time   = gfclck() */
/*                  2. Following Unix (Sun, SGI, etc), comment out */
/*                        time   = etime ( tarray ) */
/*     13 Apr 2008: Intrinsic cpu_time works with g95, gfortran */
/*     ------------------------------------------------------------------ */
/* -->  WinNT with DEC F90 */
/* -->  USE DFPORT, ONLY: RTC */
/* -->  REAL*8             dectime, decinit */
/* -->  SAVE               decinit */
/* -->  DEC OpenVMS with Fortran runtime library. */
/* -->  external        lib$init_timer, lib$stat_timer, lib$free_timer */
/* -->  integer         itimad, istatu, idata */
/* -->  save            itimad */
/* -->  DEC VAX/VMS with C runtime library. */
/* -->  integer            itimad(4) */
/* -->  PC Lahey Fortran */
/* -->  integer            itimad(4) */
/* -->  AIX */
/* -->  integer          mclock */
/* -->  intrinsic        real */
/* -->  Linux, Unix (SGI, Sun, DECstation), gfortran, g95 */
/* -->  external           etime */
/*     13 Apr 2008: Apparently it's best not to declare etime "external" */
/*                  because some compilers think of it as an intrinsic! */
/*                  http://icl.cs.utk.edu/lapack-forum/archives/lapack/msg00193.html */
/* -->  GAMS clock */
/*     external           gfclck */
/*     double precision   gfclck */
    if (*mode == 1) {
/*        --------------------------------------------------------------- */
/*        Initialize. */
/*        --------------------------------------------------------------- */
	*time = 0.f;
/* -->     DEC OpenVMS with Fortran library. */
/* -->     istatu = lib$init_timer( itimad ) */
/* -->     if (.not. istatu) call lib$signal( %val(istatu) ) */
/* -->     WinNT with DEC F90 */
/* -->     decinit = RTC() */
    } else if (*mode == 0) {
/*        --------------------------------------------------------------- */
/*        Normal call. */
/*        Return current timer value here. */
/*        --------------------------------------------------------------- */
/* -->     DEC OpenVMS with Fortran library. */
/* -->     istatu returns the number of  centiseconds. */
/* -->     istatu = lib$stat_timer( 2, idata, itimad ) */
/* -->     if (.not. istatu) call lib$signal( %val(istatu) ) */
/* -->     time   = idata */
/* -->     time   = time * 0.01d+0 */
/* -->     DEC VAX/VMS with C library. */
/* -->     itimad(1) returns the number of  centiseconds. */
/* -->     call times ( itimad ) */
/* -->     time   = itimad(1) */
/* -->     time   = time * 0.01d+0 */
/* -->     PC Lahey Fortran, itimad(1) returns the number of  centiseconds. */
/* -->     call timer ( itimad ) */
/* -->     time   = itimad(1) */
/* -->     time   = time * 0.01d+0 */
/* -->     On AIX, mclock returns hundredths of a second */
/* -->     time = real( mclock( ) ) / 100.0 */
/* -->     On Unix (SGI Irix, Sun Solaris), etime returns seconds. */
/* -->     Linux g77, Absoft f77, g95, gfortran, NagWare f95 */
/* -->     etime must be declared intrinsic for some f90 compilers */
	*time = etime_(tarray);
/* -->     g95, gfortran, NagWare f95, use f95 intrinsic. */
/* -->     call cpu_time(time) */
/* -->     On UNIX (NeXTstation M68040), using routine in ftime.c */
/* -->     call ftime(time) */
/* -->     WinNT with DEC F90 */
/* -->     Time in secs since 00:00:00 GMT Jan 1st, 1970. */
/* -->     Bias must be subtracted to give adequate precision for real*4 */
/* -->     dectime = RTC() - decinit */
/* -->     time    = real(dectime) */
/* -->     GAMS clock */
/*        time   = gfclck() */
/* -->     On other machines, to forget about timing, just say */
/* -->     time   = -1.0 */
    } else if (*mode == -1) {
/*        --------------------------------------------------------------- */
/*        Stop the clock. */
/*        --------------------------------------------------------------- */
	*time = 0.f;
/* -->     DEC OpenVMS with Fortran library. */
/* -->     istatu = lib$free_timer( itimad ) */
/* -->     if (.not. istatu) call lib$signal( %val(istatu) ) */
    }
    return 0;
} /* s1cpu_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1cpu */
doublereal s1eps_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal u, base;
    static integer nbase, ndigit;


/*     Compute the machine precision. */
/*     IEEE Floating point double precision. */

/*     ----------------------------------------------------------------- */
    nbase = 2;
    ndigit = 53;
    base = (doublereal) nbase;
    i__1 = -ndigit;
    u = pow_di(&base, &i__1);
    ret_val = u * 2.;
    return ret_val;
} /* s1eps_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* function s1eps */
doublereal s1flmx_(void)
{
    /* System generated locals */
    doublereal ret_val;


/*     IEEE Floating point double precision. */

    ret_val = 1.7977e307;
    return ret_val;
} /* s1flmx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* function s1flmx */
doublereal s1flmn_(void)
{
    /* System generated locals */
    doublereal ret_val;


/*     IEEE Floating point double precision. */

    ret_val = 2.2251e-308;
    return ret_val;
} /* s1flmn_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* function s1flmn */
integer s1intmx_(void)
{
    /* System generated locals */
    integer ret_val;


/*     2**31-1, the largest positive 32-bit integer */

    ret_val = 2147483647;
    return ret_val;
} /* s1intmx_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* function s1intmx */
integer s1inpt_(void)
{
    /* System generated locals */
    integer ret_val;


/*     Fortran standard input */

    ret_val = 5;
    return ret_val;
} /* s1inpt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* function s1inpt */
integer s1outpt_(void)
{
    /* System generated locals */
    integer ret_val;


/*     Fortran standard output */

    ret_val = 6;
    return ret_val;
} /* s1outpt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* function s1outpt */
/* Subroutine */ int s1file_(integer *task, integer *iw, integer *leniw)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 ===>  Warning: the Specs file and \002"
	    ",a,\002 file are on the same unit\002)";
    static char fmt_2000[] = "(\002 ===>  Warning: the  MPS  file and \002"
	    ",a,\002 file are on the same unit\002)";

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    static char str[80];
    static integer iback, ioldb, ipnch, inewb, idump, isoln, isumm;
    extern /* Subroutine */ int s1open_(integer *, integer *, char *, ftnlen);
    extern integer s1inpt_(void);
    static integer iloadb, ispecs, iprint, ireprt, iinsrt;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___19 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___20 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___21 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___22 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___23 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___24 = { 0, str, 0, fmt_1000, 80, 1 };
    static icilist io___25 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___26 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___27 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___28 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___29 = { 0, str, 0, fmt_2000, 80, 1 };
    static icilist io___30 = { 0, str, 0, fmt_2000, 80, 1 };


/*     ------------------------------------------------------------------ */
/*     s1file  is a machine-dependent routine for opening various files. */
/*     It calls s1open (which is also machine-dependent). */

/*     SNOPT uses sequential files only */
/*     and does not need to read and write to the same file. */

/*     iPrint, iSumm  are defined in snInit */
/*     iSpecs         is  defined in snSpec */
/*     iStdi          is defined in s1inpt, but loaded into memory here. */

/*     iStdi and iPrint have the following use: */
/*        Input  files (MPS, Old Basis, Insert, Load) */
/*        are rewound after being read, */
/*        but not if they are the same as  iStdi. */
/*        Output files (Backup, New Basis, Punch, Dump, Solution, Report) */
/*        are rewound after being written, */
/*        but not if they are the same as  iPrint. */

/*     iStdi  = (conceptually) the Keyboard that can't be rewound. */
/*              SNOPT does not use this file, so there is no 'open'. */
/*     iSumm  = the SUMMARY file.  Sometimes this is the Terminal. */
/*              If so, it may not need to be opened. */
/*     iSpecs = the SPECS file, containing one or more problem specs. */
/*              This file is not rewound after use, because it may */
/*              contain another SPECS file. */

/*     Here are all the files used by SNOPT. */
/*     The associated Index is passed to s1open */
/*     and must match the list of names in s1open, if that routine */
/*     uses method = 1. */

/*        Unit    Index    Description       Status */
/*        iSpecs     1     Specs     file    In */
/*        iPrint     2     Print     file    Out */
/*        iSumm      3     Summary   file    Out */
/*        iMPS       4     MPS       file    In */
/*        iOldB      5     Old Basis file    In */
/*        iInsrt     6     Insert    file    In */
/*        iLoadB     7     Load      file    In */
/*        iBack      8     Backup    file    Out */
/*        iNewB      9     New Basis file    Out */
/*        iPnch     10     Punch     file    Out */
/*        iDump     11     Dump      file    Out */
/*        iSoln     12     Solution  file    Out */
/*        iReprt    13     Report    file    Out */
/*        iStdi            Not opened, but used as described above */

/*     15 Nov 1991: First version based on Minos 5.4 routine mifile. */
/*     03 Oct 1997: Modified to comply with Minos 5.5. */
/*     31 Jul 2003: snPRNT adopted. */
/*     22 Nov 2003: Current version. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* Standard input */
/* MPS file */
/* Global value of iPrint */
/*     ------------------------------------------------------------------ */
/* Global value of iSumm */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ispecs = iw[11];
/* Specs (options) file */
    iprint = iw[12];
/* Print file */
    isumm = iw[13];
/* Summary file */
    iback = iw[120];
/* backup file */
    idump = iw[121];
/* dump file */
    iloadb = iw[122];
/* load file */
    inewb = iw[124];
/* new basis file */
    iinsrt = iw[125];
/* insert file */
    ioldb = iw[126];
/* old basis file */
    ipnch = iw[127];
/* punch file */
    ireprt = iw[130];
/* Report file */
    isoln = iw[131];
/*     ------------------------------------------------------------------ */
/* -->  Machine dependency. */
/*     Set iStdi = some input unit number that should not be rewound. */
/*     Set iStdi = 0 if this is irrelevant. */
/*     ------------------------------------------------------------------ */
/* Solution file */
    iw[9] = s1inpt_();
    if (*task == 2) {
/*        Relax, do nothing */
    } else if (*task == 0) {
/*        --------------------------------------------------------------- */
/*        Task = Default: Open the Specs, Print and Summary files. */
/*        iSpecs          remains the same throughout the run. */
/*        iPrint, iSumm   may be altered by the SPECS file.  They may */
/*                        need to be opened by both Tasks Deflt and Open. */
/*        --------------------------------------------------------------- */
	iw[228] = iprint;
	iw[229] = isumm;
	s1open_(&ispecs, &c__1, "IN ", (ftnlen)3);
	s1open_(&iprint, &c__2, "OUT", (ftnlen)3);
	s1open_(&isumm, &c__3, "OUT", (ftnlen)3);
    } else if (*task == 1) {
/*        --------------------------------------------------------------- */
/*        Task = OpenF: Open files mentioned in the SPECS file just read. */
/*        Input files are opened first.  Only one basis file is needed. */
/*        --------------------------------------------------------------- */
	if (iw[123] <= 0) {
	    iw[123] = ispecs;
	}
	if (iw[123] != ispecs) {
	    s1open_(&iw[123], &c__4, "IN ", (ftnlen)3);
	}
	if (ioldb > 0) {
	    s1open_(&ioldb, &c__5, "IN ", (ftnlen)3);
	} else if (iinsrt > 0) {
	    s1open_(&iinsrt, &c__6, "IN ", (ftnlen)3);
	} else if (iloadb > 0) {
	    s1open_(&iloadb, &c__7, "IN ", (ftnlen)3);
	}
	s1open_(&iback, &c__8, "OUT", (ftnlen)3);
	s1open_(&inewb, &c__9, "OUT", (ftnlen)3);
	s1open_(&ipnch, &c__10, "OUT", (ftnlen)3);
	s1open_(&idump, &c__11, "OUT", (ftnlen)3);
	s1open_(&isoln, &c__12, "OUT", (ftnlen)3);
	s1open_(&ireprt, &c__13, "OUT", (ftnlen)3);
/*        Open new Print or Summary files if they were altered */
/*        by the Specs file. */
	if (iprint != iw[228]) {
	    s1open_(&iprint, &c__2, "OUT", (ftnlen)3);
	}
	if (isumm != iw[229]) {
	    s1open_(&isumm, &c__3, "OUT", (ftnlen)3);
	}
    }
/*     Check that output files are different from Specs or MPS. */
    if (ispecs > 0) {
	if (iback == ispecs) {
	    s_wsfi(&io___19);
	    do_fio(&c__1, "Backup", (ftnlen)6);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (inewb == ispecs) {
	    s_wsfi(&io___20);
	    do_fio(&c__1, "New Basis", (ftnlen)9);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (ipnch == ispecs) {
	    s_wsfi(&io___21);
	    do_fio(&c__1, "Punch", (ftnlen)5);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (idump == ispecs) {
	    s_wsfi(&io___22);
	    do_fio(&c__1, "Dump", (ftnlen)4);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (isoln == ispecs) {
	    s_wsfi(&io___23);
	    do_fio(&c__1, "Solution", (ftnlen)8);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (ireprt == ispecs) {
	    s_wsfi(&io___24);
	    do_fio(&c__1, "Report", (ftnlen)6);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    if (iw[123] > 0) {
	if (iback == iw[123]) {
	    s_wsfi(&io___25);
	    do_fio(&c__1, "Backup", (ftnlen)6);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (inewb == iw[123]) {
	    s_wsfi(&io___26);
	    do_fio(&c__1, "New Basis", (ftnlen)9);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (ipnch == iw[123]) {
	    s_wsfi(&io___27);
	    do_fio(&c__1, "Punch", (ftnlen)5);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (idump == iw[123]) {
	    s_wsfi(&io___28);
	    do_fio(&c__1, "Dump", (ftnlen)4);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (isoln == iw[123]) {
	    s_wsfi(&io___29);
	    do_fio(&c__1, "Solution", (ftnlen)8);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
	if (ireprt == iw[123]) {
	    s_wsfi(&io___30);
	    do_fio(&c__1, "Report", (ftnlen)6);
	    e_wsfi();
	    snprnt_(&c__1, str, &iw[1], leniw, (ftnlen)80);
	}
    }
    return 0;
} /* s1file_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1file */
/* Subroutine */ int s1clos_(integer *lun)
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    integer f_clos(cllist *);

/*     ================================================================== */
/*     s1clos  closes the file with logical unit number lun. */
/*     This version is trivial and so far is not even used by SNOPT. */
/*     Perhaps some implementations will need something fancier. */
/*     ================================================================== */
    cl__1.cerr = 0;
    cl__1.cunit = *lun;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* s1clos_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1clos */
/* Subroutine */ int s1envt_(integer *mode, integer *iw, integer *leniw)
{
/*     ================================================================== */
/*     s1envt specifies the environment within which SNOPT is being used. */

/*     When mode = 0, information about the environment should be */
/*     initialized. */

/*     iPage1 says whether new pages are ever wanted on file iPrint. */
/*     iPage2 says whether new pages are ever wanted on file iSumm. */

/*     When mode is in the range 1 to 99, each environment does its */
/*     own thing. */

/*     The only environment at present is: */
/*     ALONE: */
/*     This means SNOPT is in stand-alone mode---the normal case. */
/*     Nothing special is done. */

/*     16 Sep 1987. */
/*     ================================================================== */
/* > 0    =>  stand-alone */
/* > 0    =>  Page 1 */
/*     ------------------------------------------------------------------ */
/* > 0    =>  Page 2 */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    if (*mode <= 0) {
/*        --------------------------------------------------------------- */
/*        mode = 0.    Initialize. */
/*        --------------------------------------------------------------- */
	iw[238] = 1;
	iw[241] = 1;
	iw[242] = 0;
    } else {
/*        Relax, do nothing in this version. */
    }
    return 0;
} /* s1envt_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1envt */
/* Subroutine */ int s1open_(integer *lun, integer *index, char *state, 
	ftnlen state_len)
{
    /* Initialized data */

    static char names[9*13] = "snopt.spc" "snopt.prn" "snopt.sum" "snopt.mps" 
	    "snopt.olb" "snopt.ins" "snopt.lod" "snopt.bak" "snopt.nwb" "sno"
	    "pt.pun" "snopt.dmp" "snopt.sol" "snopt.rpt";

    /* System generated locals */
    integer i__1;
    olist o__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_cmp(char *, char *, ftnlen, ftnlen), f_open(
	    olist *), s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), 
	    e_wsfi(void), f_rew(alist *);

    /* Local variables */
    static integer last;
    static logical uopen, input;
    static char filnam[100];
    static integer screen, method;
    extern /* Subroutine */ int getfnm_(integer *, char *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___37 = { 0, filnam, 0, "(a,i1)", 100, 1 };
    static icilist io___38 = { 0, filnam, 0, "(a,i2)", 100, 1 };
    static icilist io___39 = { 0, filnam, 0, "(a,i1)", 100, 1 };
    static icilist io___40 = { 0, filnam, 0, "(a,i2)", 100, 1 };


/*     ------------------------------------------------------------------ */
/*     s1open  is a machine-dependent routine. */
/*     In principal it opens a file with logical unit number lun */
/*     and positions it at the beginning. */

/*     Input files are treated that way. */
/*     An input file is opened with status='OLD' */
/*     (and F77 will terminate with an error if the file doesn't exist). */

/*     Output files are more machine-dependent. */
/*     With status='NEW', F77 would terminate if the file DID exist. */
/*     With status='UNKNOWN', existing files would be overwritten. */
/*     This is normal with Unix, but on systems that have file */
/*     version numbers (e.g. DEC OpenVMS), a new version is preferable. */
/*     It is then better not to open the file at all, but let a new */
/*     version be created by the first "write". */

/*     Nothing happens if */
/*     1. lun <= 0.  Saves us from testing lun before calling s1open. */

/*     2. Unit lun is already open.  Helps applications that call */
/*        snopt -- they can open files themselves if they want to. */

/*     3. lun = screen, where  screen  is a local machine-dependent */
/*        variable, typically 6.  It seems inadvisable to do an OPEN */
/*        on a file that is predefined to be an interactive screen. */
/*        With Unix on DECstations, open(6, file='snopt.sum') sends */
/*        output to file snopt.sum, not to the screen. */


/*     lun     (input) is the unit number. */

/*     index   (input) points to one of the hardwired names below. */
/*             Used only if method = 1. */

/*     state   (input) is 'IN ' or 'OUT', indicating whether the file */
/*             is to be input or output. */

/*     15 Jul 1989: First version, follows some of the advice offered */
/*                  by David Gay, Bell Laboratories. */
/*     -- --- 1990: Added parameter "state". */
/*     03 Feb 1994: Added parameter "index" to help when method = 1. */
/*                  Local variable "method" must be set here to select */
/*                  various methods for naming files. */
/*                  Chris Jaensch, IFR Stuttgart, recommends not opening */
/*                  input files if they are already open.  This version */
/*                  ignores all open files (input or output), assuming */
/*                  that they are taken care of by the calling program. */
/*     13 May 1994: Local variable "screen" used to avoid opening screen. */
/*     20 May 2001: Method 4 is now the default. The default getfnm uses */
/*                  method 2. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* -->  Machine dependency. */
/*     names(*) is needed if method = 1 below. */
/*     Make sure "character*n" sets n big enough below. */
/*     It doesn't matter if it is bigger than necessary, since */
/*     "open( lun, file=name )" allows name to have trailing blanks. */
/*     ------------------------------------------------------------------ */
/*     ------------------------------------------------------------------ */
/* -->  Machine-dependency. */
/*     Set "method" to suit your operating system. */
/*     It determines how a file name is assigned to unit number "lun". */
/*     Typically, */
/*     method = 1 for fixed file names (e.g. PCs under DOS). */
/*     method = 2 for names like fort.7, fort.15 (e.g. Unix on Sun, SGI). */
/*     method = 3 for names like FTN07 , FTN15   (e.g. Unix on HP). */
/*     method = 4 if names are assigned by an external routine getfnm */
/*                (default getfnm sets names as in method 2). */
/*     method = 5 if explicit file names are not needed (e.g. OpenVMS). */
/*     See more comments below. */

/*     Set "screen" to a unit number that never needs to be opened. */
/*     (Typically, screen = 6.  If unknown, set screen = 0.) */
/*     ------------------------------------------------------------------ */
    method = 4;
/* Default value for SNOPT */
    screen = 6;
/*     ------------------------------------------------------------------ */
/*     Quit if lun<=0 or lun = iscreen or unit lun is already open. */
/*     ------------------------------------------------------------------ */
    if (*lun <= 0) {
	goto L900;
    }
    if (*lun == screen) {
	goto L900;
    }
    ioin__1.inerr = 0;
    ioin__1.inunit = *lun;
    ioin__1.infile = 0;
    ioin__1.inex = 0;
    ioin__1.inopen = &uopen;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (uopen) {
	goto L900;
    }
/*     ------------------------------------------------------------------ */
/*     Open file lun by a specified machine-dependent method. */
/*     Note that 'UNKNOWN' is equivalent to trying first with 'OLD', */
/*     and then with 'NEW' if the file doesn't exist. */
/*     ------------------------------------------------------------------ */
    input = s_cmp(state, "IN ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(state, 
	    "in ", (ftnlen)3, (ftnlen)3) == 0;
    if (method == 1) {
/*        --------------------------------------------------------------- */
/*        Hardwired filename. */
/*        We use "index" to get it from names(*) above. */
/*        Typical machines: IBM PC under DOS. */
/*        --------------------------------------------------------------- */
	if (input) {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = 9;
	    o__1.ofnm = names + (*index - 1) * 9;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	} else {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = 9;
	    o__1.ofnm = names + (*index - 1) * 9;
	    o__1.orl = 0;
	    o__1.osta = "UNKNOWN";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
    } else if (method == 2) {
/*        --------------------------------------------------------------- */
/*        Construct a name like fort.7 or fort.15. */
/*        (This approach suggested by Chris Jaensch, IFR Stuttgart.) */
/*        Typical machines:  Unix on Sun, SGI, DEC. */
/*        --------------------------------------------------------------- */
	if (*lun <= 9) {
	    s_wsfi(&io___37);
	    do_fio(&c__1, "fort.", (ftnlen)5);
	    do_fio(&c__1, (char *)&(*lun), (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
	    s_wsfi(&io___38);
	    do_fio(&c__1, "fort.", (ftnlen)5);
	    do_fio(&c__1, (char *)&(*lun), (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	if (input) {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = 100;
	    o__1.ofnm = filnam;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	} else {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = 100;
	    o__1.ofnm = filnam;
	    o__1.orl = 0;
	    o__1.osta = "UNKNOWN";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
    } else if (method == 3) {
/*        --------------------------------------------------------------- */
/*        Construct a name like FTN07 or FTN15. */
/*        Typical machines:  Unix on HP. */
/*        --------------------------------------------------------------- */
	if (*lun <= 9) {
	    s_wsfi(&io___39);
	    do_fio(&c__1, "FTN0", (ftnlen)4);
	    do_fio(&c__1, (char *)&(*lun), (ftnlen)sizeof(integer));
	    e_wsfi();
	} else {
	    s_wsfi(&io___40);
	    do_fio(&c__1, "FTN", (ftnlen)3);
	    do_fio(&c__1, (char *)&(*lun), (ftnlen)sizeof(integer));
	    e_wsfi();
	}
	if (input) {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = 100;
	    o__1.ofnm = filnam;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	} else {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = 100;
	    o__1.ofnm = filnam;
	    o__1.orl = 0;
	    o__1.osta = "UNKNOWN";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
    } else if (method == 4) {
/*        --------------------------------------------------------------- */
/*        Assume some routine "getfnm" will provide a name at run-time. */
/*        (This approach is used by David Gay, Bell Labs.) */
/*        The default getfnm.f provided in sn12ampl.f constructs names */
/*        as in method 2. */
/*        --------------------------------------------------------------- */
	getfnm_(lun, filnam, &last, (ftnlen)100);
	if (input) {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = last;
	    o__1.ofnm = filnam;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	} else {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnmlen = last;
	    o__1.ofnm = filnam;
	    o__1.orl = 0;
	    o__1.osta = "UNKNOWN";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
    } else if (method == 5) {
/*        --------------------------------------------------------------- */
/*        Explicit file names are not needed for the open statement. */
/*        The operating system uses a default name */
/*        (e.g. fort.7 or for007) */
/*        or has already assigned a name to this unit */
/*        (e.g. via a script or command file). */
/*        Typical machines:  Unix, */
/*                           DEC OpenVMS, */
/*                           IBM Mainframes. */
/*        --------------------------------------------------------------- */
	if (input) {
	    o__1.oerr = 0;
	    o__1.ounit = *lun;
	    o__1.ofnm = 0;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	} else {
/*           Let the first "write" do it. */
	}
    }
/*     ------------------------------------------------------------------ */
/*     Rewind input files. */
/*     (Some systems position existing files at the end */
/*     rather than the beginning.) */
/*     err=900 covers files that have not yet been opened. */
/*     ------------------------------------------------------------------ */
    if (input) {
	al__1.aerr = 1;
	al__1.aunit = *lun;
	i__1 = f_rew(&al__1);
	if (i__1 != 0) {
	    goto L900;
	}
    }
L900:
    return 0;
} /* s1open_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* subroutine s1open */
/* Subroutine */ int s1page_(integer *mode, integer *iw, integer *leniw)
{
    static integer ipage1, ipage2;
    extern /* Subroutine */ int snprnt_(integer *, char *, integer *, integer 
	    *, ftnlen);

/*     ------------------------------------------------------------------ */
/*     s1page is an installation-dependent routine.  It is called at */
/*     points where some users might want output to files iPrint or iSumm */
/*     to begin on a new page. */

/*     iPage1 and iPage2 have already been set by s1envt. */
/*     If they are true, a page eject and a blank line are output. */
/*     Otherwise, just a blank line is output. */

/*     If mode = 0  and Summary level = 0, nothing is output to the */
/*                  Summary file.  At present, this is so s8log */
/*                  will print just one line per major iteration, with */
/*                  no blank line in between. */
/*     If mode = 1, just the page control is relevant. */
/*     If mode = 2, SNOPT has encountered an error condition. */
/*                  At the moment, this case is treated the same as */
/*                  mode = 1. */

/*     15 Nov 1991: First version based on Minos 5.4 routine m1page. */
/*     31 Jul 2003: snPRNT adopted. */
/*     31 Jul 2003: Current version of s1page. */
/*     ================================================================== */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ipage1 = iw[241];
/* > 0    =>  Page 1 */
    ipage2 = iw[242];
/* > 0    =>  Page 2 */
    if (ipage1 > 0) {
	snprnt_(&c__1, "1", &iw[1], leniw, (ftnlen)1);
    }
    snprnt_(&c__1, " ", &iw[1], leniw, (ftnlen)1);
    if (ipage2 > 0) {
	snprnt_(&c__2, "1", &iw[1], leniw, (ftnlen)1);
    }
    if (*mode != 0) {
	snprnt_(&c__2, " ", &iw[1], leniw, (ftnlen)1);
    }
    return 0;
} /* s1page_ */

