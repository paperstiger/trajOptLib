*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     snFilewrapper   snOpenappend   snClose   snOpenread
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFilewrapper
     &   ( name, ispec, inform, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     snSpec
      character*(*)
     &     name
      integer
     &     inform, lencw, leniw, lenrw, iw(leniw), ispec
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     Read options for snopt from the file named name. inform .eq 0 if
*     successful.
*
*     09 Jan 2000: First version of snFilewrapper.
*     ==================================================================
      integer
     &     iostat

      open( ispec, iostat=iostat, file=name, status='old' )
      if (iostat .ne. 0) then
         inform = 2 + iostat
      else
         call snSpec( ispec, inform, cw, lencw, iw, leniw, rw, lenrw )
         close( ispec )
      end if

      end ! subroutine snFilewrapper

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snOpenappend
     &   ( iunit, name, inform )

      implicit
     &     none
      integer
     &     iunit
      character*(*)
     &     name
      integer
     &     inform

*     ==================================================================
*     Open file named name to Fortran unit iunit. inform .eq. 0 if
*     sucessful.  Although opening for appending is not in the f77
*     standard, it is understood by f2c.
*
*     09 Jan 2000: First version of snOpenappend
*     ==================================================================
      open( iunit, iostat=inform, file=name, access='append' )

      end ! subroutine snOpenappend

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snClose
     &   ( iunit )

      implicit
     &     none
      integer
     &     iunit

*     ==================================================================
*     Close unit iunit.
*
*     09 Jan 2000: First version of snClose
*     ==================================================================
      close( iunit )

      end ! subroutine snClose

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snOpenread
     &   ( iunit, name, inform )

      implicit
     &     none
      integer
     &     iunit
      character*(*)
     &     name
      integer
     &     inform

*     ==================================================================
*     Open file called name to Fortran unit iunit. inform .eq. 0 if
*     sucessful.  Although opening for appending is not in the f77
*     standard, it is understood by f2c.
*
*     09 Jan 2000: First version of snOpenread
*     ==================================================================
      open ( iunit, iostat=inform, file=name, form = 'FORMATTED',
     *     status = 'OLD' )

      end ! subroutine snOpenread
