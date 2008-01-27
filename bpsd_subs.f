c     $Id$
c=======================================================================
      module bpsd_subs
c
      contains
c
c-----------------------------------------------------------------------
      subroutine spl1D_bpsd(data1D,nd,ierr)
c-----------------------------------------------------------------------
c
!      use libspl_mod
      use bpsd_types_internal
      implicit none
      type(bpsd_1ddatax_type) :: data1D
      integer :: nd     ! position of dependent variable
      integer :: ierr    ! error indicator
      real(8), dimension(:), allocatable :: deriv
c
      allocate(deriv(data1D%nrmax))
      call spl1D(data1D%s,data1D%data(1,nd),deriv,data1D%spline(1,1,nd),
     &           data1D%nrmax,0,ierr)
      if(ierr.ne.0) 
     &     write(6,*) 'XX spl1D_bpsd : '//data1D%kid(nd)//
     &                ': ierr=',ierr
      deallocate(deriv)
      return
      end subroutine spl1D_bpsd
c
c-----------------------------------------------------------------------
      subroutine spl1DF_bpsd(pos,val,data1D,nd,ierr)
c-----------------------------------------------------------------------
c
!      use libspl_mod
      use bpsd_types_internal
      implicit none
      real(8) :: pos     ! value of independent variable
      real(8) :: val     ! value of dependent variable
      type(bpsd_1ddatax_type) :: data1D
      integer :: nd      ! position of dependent variable
      integer :: ierr    ! error indicator

      call spl1DF(pos,val,
     &            data1D%s,data1D%spline(1,1,nd),data1D%nrmax,ierr)
      if(ierr.ne.0) then
         write(6,*) 'XX spl1DF_bpsd : '//data1D%kid(nd)//
     &                ': ierr=',ierr
         write(6,'(1P3E12.4)')  
     &        pos,data1D%s(1),data1D%s(data1D%nrmax)
      endif
      end subroutine spl1DF_bpsd
c
      end module bpsd_subs
