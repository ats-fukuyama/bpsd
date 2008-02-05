c     $Id$
c=======================================================================
      module bpsd_metric1D
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_metric1D,bpsd_get_metric1D,
     &       bpsd_save_metric1D,bpsd_load_metric1D
      private
c
      logical, save :: bpsd_metric1Dx_init_flag = .TRUE.
      type(bpsd_data1Dx_type), save :: metric1Dx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_metric1Dx_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      metric1Dx%status=0
      metric1Dx%dataName='metric1D'
      metric1Dx%ndmax=15
      allocate(metric1Dx%kid(15))
      metric1Dx%kid( 1)='metric1D%pvol'
      metric1Dx%kid( 2)='metric1D%psur'
      metric1Dx%kid( 3)='metric1D%dvpsit'
      metric1Dx%kid( 4)='metric1D%dvpsip'
      metric1Dx%kid( 5)='metric1D%aver2'
      metric1Dx%kid( 6)='metric1D%aver2i'
      metric1Dx%kid( 7)='metric1D%aveb2'
      metric1Dx%kid( 8)='metric1D%aveb2i'
      metric1Dx%kid( 9)='metric1D%avegv2'
      metric1Dx%kid(10)='metric1D%avegvr2'
      metric1Dx%kid(11)='metric1D%avegpp2'
      metric1Dx%kid(12)='metric1D%rr'
      metric1Dx%kid(13)='metric1D%rs'
      metric1Dx%kid(14)='metric1D%elip'
      metric1Dx%kid(15)='metric1D%trig'
c
      bpsd_metric1Dx_init_flag = .FALSE.
c
      return
      end subroutine bpsd_metric1Dx_init
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_metric1D(metric1D_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_metric1D_type):: metric1D_in
      integer :: ierr
      integer :: nr,nd
c
      if(bpsd_metric1Dx_init_flag) call bpsd_metric1Dx_init
c
      if(metric1Dx%status.ne.0) then
         if(metric1D_in%nrmax.ne.metric1Dx%nrmax) then
            if(metric1Dx%status.ge.3) deallocate(metric1Dx%spline)
            deallocate(metric1Dx%data)
            deallocate(metric1Dx%s)
            metric1Dx%status=0
         endif
      endif
c
      if(metric1Dx%status.eq.0) then
         metric1Dx%nrmax=metric1D_in%nrmax
         allocate(metric1Dx%s(metric1Dx%nrmax))
         allocate(metric1Dx%data(metric1Dx%nrmax,metric1Dx%ndmax))
         metric1Dx%status=1
      endif
c
      metric1Dx%time = metric1D_in%time
      do nr=1,metric1D_in%nrmax
         metric1Dx%s(nr) = metric1D_in%s(nr)
         metric1Dx%data(nr, 1) = metric1D_in%data(nr)%pvol
         metric1Dx%data(nr, 2) = metric1D_in%data(nr)%psur
         metric1Dx%data(nr, 3) = metric1D_in%data(nr)%dvpsit
         metric1Dx%data(nr, 4) = metric1D_in%data(nr)%dvpsip
         metric1Dx%data(nr, 5) = metric1D_in%data(nr)%aver2
         metric1Dx%data(nr, 6) = metric1D_in%data(nr)%aver2i
         metric1Dx%data(nr, 7) = metric1D_in%data(nr)%aveb2
         metric1Dx%data(nr, 8) = metric1D_in%data(nr)%aveb2i
         metric1Dx%data(nr, 9) = metric1D_in%data(nr)%avegv2
         metric1Dx%data(nr,10) = metric1D_in%data(nr)%avegvr2
         metric1Dx%data(nr,11) = metric1D_in%data(nr)%avegpp2
         metric1Dx%data(nr,12) = metric1D_in%data(nr)%rr
         metric1Dx%data(nr,13) = metric1D_in%data(nr)%rs
         metric1Dx%data(nr,14) = metric1D_in%data(nr)%elip
         metric1Dx%data(nr,15) = metric1D_in%data(nr)%trig
      enddo
      if(metric1Dx%status.ge.3) then
         metric1Dx%status=3
      else
         metric1Dx%status=2
      endif
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_set_metric1D'
         write(6,*) '---- metric1Dx%s'
         write(6,'(1P5E12.4)') 
     &        (metric1Dx%s(nr),nr=1,metric1Dx%nrmax)
         do nd=1,metric1Dx%ndmax
            write(6,*) '---- ',metric1Dx%kid(nd)
            write(6,'(1P5E12.4)') 
     &           (metric1Dx%data(nr,nd),nr=1,metric1Dx%nrmax)
         enddo
      endif
      return
      end subroutine bpsd_set_metric1D
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_metric1D(metric1D_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_metric1D_type) :: metric1D_out
      integer :: ierr
      real(8), dimension(:,:), allocatable ::  temp
      real(8), dimension(:), allocatable ::  deriv
      integer :: nr, nd
      real(8) :: s
      real(8), dimension(15) :: v
c
      if(bpsd_metric1Dx_init_flag) call bpsd_metric1Dx_init
c
      if(metric1Dx%status.eq.0) then
         write(6,*) 
     &   'XX bpsd_get_metric1D: no space allocated to metric1Dx%data'
         ierr=1
         return
      endif
c
      if(metric1Dx%status.eq.1) then
         write(6,*) 'XX bpsd_get_metric1D: no data in metric1Dx%data'
         ierr=2
         return
      endif
c
      if(metric1D_out%nrmax.eq.0) then
         if(allocated(metric1D_out%data)) then
            if(metric1Dx%nrmax.ne.size(metric1D_out%data,1)) then
               deallocate(metric1D_out%data)
               metric1D_out%nrmax = metric1Dx%nrmax
               allocate(metric1D_out%data(metric1D_out%nrmax))
            endif
         else
            metric1D_out%nrmax = metric1Dx%nrmax
            allocate(metric1D_out%data(metric1D_out%nrmax))
         endif
      endif
c
      if(allocated(metric1D_out%data)) then
         if(metric1Dx%nrmax.le.size(metric1D_out%data,1)) then
            metric1D_out%time  = metric1Dx%time
            metric1D_out%nrmax = metric1Dx%nrmax
            do nr=1,metric1Dx%nrmax
               metric1D_out%data(nr)%pvol     = metric1Dx%data(nr, 1)
               metric1D_out%data(nr)%psur     = metric1Dx%data(nr, 2)
               metric1D_out%data(nr)%dvpsit   = metric1Dx%data(nr, 3)
               metric1D_out%data(nr)%dvpsip   = metric1Dx%data(nr, 4)
               metric1D_out%data(nr)%aver2    = metric1Dx%data(nr, 5)
               metric1D_out%data(nr)%aver2i   = metric1Dx%data(nr, 6)
               metric1D_out%data(nr)%aveb2    = metric1Dx%data(nr, 7)
               metric1D_out%data(nr)%aveb2i   = metric1Dx%data(nr, 8)
               metric1D_out%data(nr)%avegv2   = metric1Dx%data(nr, 9)
               metric1D_out%data(nr)%avegvr2  = metric1Dx%data(nr,10)
               metric1D_out%data(nr)%avegpp2  = metric1Dx%data(nr,11)
               metric1D_out%data(nr)%rr       = metric1Dx%data(nr,12)
               metric1D_out%data(nr)%rs       = metric1Dx%data(nr,13)
               metric1D_out%data(nr)%elip     = metric1Dx%data(nr,14)
               metric1D_out%data(nr)%trig     = metric1Dx%data(nr,15)
            enddo
            ierr=0
            return
         endif
      else
         ierr=3
         return
      endif
c
      if(metric1Dx%status.eq.2) then
         allocate(metric1Dx%spline(4,metric1Dx%nrmax,metric1Dx%ndmax))
         metric1Dx%status=3
      endif
c
      if(metric1Dx%status.eq.3) then
         do nd=1,metric1Dx%ndmax
            call bpsd_spl1D(metric1Dx,nd,ierr)
         enddo
         metric1Dx%status=4
      endif
c
      do nr=1,metric1D_out%nrmax
         s = metric1D_out%s(nr)
         do nd=1,metric1Dx%ndmax
            call bpsd_spl1DF(s,v(nd),metric1Dx,nd,ierr)
         enddo
         metric1D_out%data(nr)%pvol     = v( 1)
         metric1D_out%data(nr)%psur     = v( 2)
         metric1D_out%data(nr)%dvpsit   = v( 3)
         metric1D_out%data(nr)%dvpsip   = v( 4)
         metric1D_out%data(nr)%aver2    = v( 5)
         metric1D_out%data(nr)%aver2i   = v( 6)
         metric1D_out%data(nr)%aveb2    = v( 7)
         metric1D_out%data(nr)%aveb2i   = v( 8)
         metric1D_out%data(nr)%avegv2   = v( 9)
         metric1D_out%data(nr)%avegvr2  = v(10)
         metric1D_out%data(nr)%avegpp2  = v(11)
         metric1D_out%data(nr)%rr       = v(12)
         metric1D_out%data(nr)%rs       = v(13)
         metric1D_out%data(nr)%elip     = v(14)
         metric1D_out%data(nr)%trig     = v(15)
      enddo
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_get_metric1D'
         write(6,*) '---- metric1Dx%s'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%s(nr),nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%pvol'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%pvol,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%psur'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%psur,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%dvpssit'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%dvpsit,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%dvpssit'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%dvpsit,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%aver2'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%aver2,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%aver2i'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%aver2i,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%aveb2'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%aveb2,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%aveb2i'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%aveb2i,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%avegv2'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%avegv2,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%avegvr2'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%avegvr2,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%avegpp2'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%avegpp2,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%rr'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%rr,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%rs'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%rs,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%elip'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%elip,nr=1,metric1D_out%nrmax)
         write(6,*) '---- metric1Dx%trig'
         write(6,'(1P5E12.4)') 
     &        (metric1D_out%data(nr)%trig,nr=1,metric1D_out%nrmax)
      endif
      return
      end subroutine bpsd_get_metric1D
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_metric1D(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      integer,intent(in) :: fid
      integer,intent(out) :: ierr
c
      if(bpsd_metric1Dx_init_flag) call bpsd_metric1Dx_init
c
      if(metric1Dx%status.gt.1) 
     &     call bpsd_save_data1Dx(fid,metric1Dx,ierr)
      return
c
      end subroutine bpsd_save_metric1D
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_metric1D(datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_data1Dx_type),intent(in) :: datax
      integer,intent(out) :: ierr
      integer:: ns,nr,nd
c
      if(bpsd_metric1Dx_init_flag) call bpsd_metric1Dx_init
c
      if(metric1Dx%status.ne.0) then
         if(datax%nrmax.ne.metric1Dx%nrmax) then
            if(metric1Dx%status.ge.3) deallocate(metric1Dx%spline)
            deallocate(metric1Dx%data)
            deallocate(metric1Dx%s)
            metric1Dx%status=0
         endif
      endif
c
      if(metric1Dx%status.eq.0) then
         metric1Dx%nrmax=datax%nrmax
         allocate(metric1Dx%s(metric1Dx%nrmax))
         allocate(metric1Dx%data(metric1Dx%nrmax,metric1Dx%ndmax))
         metric1Dx%status=1
      endif
c
      metric1Dx%time = datax%time
      do nr=1,metric1Dx%nrmax
         metric1Dx%s(nr) = datax%s(nr)
         do nd=1,metric1Dx%ndmax
            metric1Dx%data(nr,nd) = datax%data(nr,nd)
         enddo
      enddo
      metric1Dx%status=2
      ierr=0
      return
c
      end subroutine bpsd_load_metric1D
c
      end module bpsd_metric1D
