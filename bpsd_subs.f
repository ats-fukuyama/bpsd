c     $Id$
c=======================================================================
      module bpsd_subs
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_spl1D(data1D,nd,ierr)
c-----------------------------------------------------------------------
c
!      use libspl_mod
      use bpsd_types_internal
      implicit none
      type(bpsd_data1Dx_type) :: data1D
      integer :: nd     ! position of dependent variable
      integer :: ierr    ! error indicator
      real(8), dimension(:), pointer :: deriv
c
      allocate(deriv(data1D%nrmax))
      call spl1D(data1D%s,data1D%data(1,nd),deriv,data1D%spline(1,1,nd),
     &           data1D%nrmax,0,ierr)
      if(ierr.ne.0) 
     &     write(6,*) 'XX spl1D_bpsd : '//data1D%kid(nd)//
     &                ': ierr=',ierr
      deallocate(deriv)
      return
      end subroutine bpsd_spl1D
c
c-----------------------------------------------------------------------
      subroutine bpsd_spl1DF(pos,val,data1D,nd,ierr)
c-----------------------------------------------------------------------
c
!      use libspl_mod
      use bpsd_types_internal
      implicit none
      real(8) :: pos     ! value of independent variable
      real(8) :: val     ! value of dependent variable
      type(bpsd_data1Dx_type) :: data1D
      integer :: nd      ! position of dependent variable
      integer :: ierr    ! error indicator

      call spl1DF(pos**2,val,
     &            data1D%s,data1D%spline(1,1,nd),data1D%nrmax,ierr)
      if(ierr.ne.0) then
         write(6,*) 'XX spl1DF_bpsd : '//data1D%kid(nd)//
     &                ': ierr=',ierr
         write(6,'(1P3E12.4)')  
     &        pos**2,data1D%s(1),data1D%s(data1D%nrmax)
      endif
      end subroutine bpsd_spl1DF
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_shotx(fid,datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_shotx_type),intent(in):: datax
      integer,intent(out):: ierr
c
      write(fid,IOSTAT=ierr,ERR=8) 'data:str'
      write(fid,IOSTAT=ierr,ERR=8) datax%dataName
      write(fid,IOSTAT=ierr,ERR=8) datax%deviceID
      write(fid,IOSTAT=ierr,ERR=8) datax%shotID,datax%modelID

      write(6,*) '+++ ',datax%dataName,': saved'
      ierr=0
      return

    8 continue
      return
      end subroutine bpsd_save_shotx
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_data0Dx(fid,datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data0Dx_type),intent(in):: datax
      integer,intent(out):: ierr
c
      write(fid,IOSTAT=ierr,ERR=8) 'data::0D'
      write(fid,IOSTAT=ierr,ERR=8) datax%dataName
      write(fid,IOSTAT=ierr,ERR=8) datax%time
      write(fid,IOSTAT=ierr,ERR=8) datax%ndmax
      write(fid,IOSTAT=ierr,ERR=8) datax%created_date,
     &                             datax%created_time,
     &                             datax%created_timezone
      write(fid,IOSTAT=ierr,ERR=8) datax%kid
      write(fid,IOSTAT=ierr,ERR=8) datax%kunit
      write(fid,IOSTAT=ierr,ERR=8) datax%data

      write(6,*) '+++ ',datax%dataName,': saved'
      ierr=0
      return

    8 continue
      write(6,*) 'XXX', ierr
      return
      end subroutine bpsd_save_data0Dx
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_data1Dx(fid,datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data1Dx_type),intent(in):: datax
      integer,intent(out):: ierr
c
      write(fid,IOSTAT=ierr,ERR=8) 'data::1D'
      write(fid,IOSTAT=ierr,ERR=8) datax%dataName
      write(fid,IOSTAT=ierr,ERR=8) datax%time
      write(fid,IOSTAT=ierr,ERR=8) datax%nrmax,datax%ndmax
      write(fid,IOSTAT=ierr,ERR=8) datax%created_date,
     &                             datax%created_time,
     &                             datax%created_timezone
      write(fid,IOSTAT=ierr,ERR=8) datax%kid
      write(fid,IOSTAT=ierr,ERR=8) datax%kunit
      write(fid,IOSTAT=ierr,ERR=8) datax%rho
      write(fid,IOSTAT=ierr,ERR=8) datax%data

      write(6,*) '+++ ',datax%dataName,': saved'
      ierr=0
      return

    8 continue
      return
      end subroutine bpsd_save_data1Dx
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_data2Dx(fid,datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data2Dx_type),intent(in):: datax
      integer,intent(out):: ierr
c
      write(fid,IOSTAT=ierr,ERR=8) 'data::2D'
      write(fid,IOSTAT=ierr,ERR=8) datax%dataName
      write(fid,IOSTAT=ierr,ERR=8) datax%time
      write(fid,IOSTAT=ierr,ERR=8) datax%nrmax,datax%nthmax,datax%ndmax
      write(fid,IOSTAT=ierr,ERR=8) datax%created_date,
     &                             datax%created_time,
     &                             datax%created_timezone
      write(fid,IOSTAT=ierr,ERR=8) datax%kid
      write(fid,IOSTAT=ierr,ERR=8) datax%rho
      write(fid,IOSTAT=ierr,ERR=8) datax%th
      write(fid,IOSTAT=ierr,ERR=8) datax%data

      write(6,*) '+++ ',datax%dataName,': saved'
      ierr=0
      return

    8 continue
      return
      end subroutine bpsd_save_data2Dx
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_data3Dx(fid,datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data3Dx_type),intent(in):: datax
      integer,intent(out):: ierr
c
      write(fid,IOSTAT=ierr,ERR=8) 'data::3D'
      write(fid,IOSTAT=ierr,ERR=8) datax%dataName
      write(fid,IOSTAT=ierr,ERR=8) datax%time
      write(fid,IOSTAT=ierr,ERR=8) 
     &     datax%nrmax,datax%nthmax,datax%nphmax,datax%ndmax
      write(fid,IOSTAT=ierr,ERR=8) datax%created_date,
     &                             datax%created_time,
     &                             datax%created_timezone
      write(fid,IOSTAT=ierr,ERR=8) datax%kid
      write(fid,IOSTAT=ierr,ERR=8) datax%rho
      write(fid,IOSTAT=ierr,ERR=8) datax%th
      write(fid,IOSTAT=ierr,ERR=8) datax%ph
      write(fid,IOSTAT=ierr,ERR=8) datax%data

      write(6,*) '+++ ',datax%dataName,': saved'
      ierr=0
      return

    8 continue
      return
      end subroutine bpsd_save_data3Dx
      end module bpsd_subs
