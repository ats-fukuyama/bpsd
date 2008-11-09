c     $Id$
c=======================================================================
      module bpsd_plasmaf
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_plasmaf,bpsd_get_plasmaf,
     &       bpsd_save_plasmaf,bpsd_load_plasmaf
      private
c
      logical, save :: bpsd_plasmafx_init_flag = .TRUE.
      type(bpsd_data1Dx_type), save :: plasmafx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_plasmafx_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      plasmafx%status=0
      plasmafx%dataName='plasmaf'
      plasmafx%ndmax=0
c
      bpsd_plasmafx_init_flag = .FALSE.
c
      return
      end subroutine bpsd_plasmafx_init
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_plasmaf(plasmaf_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_plasmaf_type),intent(in):: plasmaf_in
      integer,intent(out) :: ierr
      integer :: ns,nr,nd
c
      if(bpsd_plasmafx_init_flag) call bpsd_plasmafx_init
c
      if(plasmafx%status.ne.0) then
         if((plasmaf_in%nrmax.ne.plasmafx%nrmax) .or.
     &      (plasmaf_in%nsmax*5.ne.plasmafx%ndmax))  then
            if(plasmafx%status.ge.3) deallocate(plasmafx%spline)
            deallocate(plasmafx%data)
            deallocate(plasmafx%rho)
            deallocate(plasmafx%kid)
            plasmafx%status=0
         endif
      endif
c
      if(plasmafx%status.eq.0) then
         plasmafx%ndmax=plasmaf_in%nsmax*5+1
         plasmafx%nrmax=plasmaf_in%nrmax
         allocate(plasmafx%kid(plasmafx%ndmax))
         allocate(plasmafx%kunit(plasmafx%ndmax))
         allocate(plasmafx%rho(plasmafx%nrmax))
         allocate(plasmafx%data(plasmafx%nrmax,plasmafx%ndmax))
         do ns=1,plasmaf_in%nsmax
            nd=5*(ns-1)
            plasmafx%kid(nd+1)='plasmaf%pn'
            plasmafx%kid(nd+2)='plasmaf%pt'
            plasmafx%kid(nd+3)='plasmaf%ptpr'
            plasmafx%kid(nd+4)='plasmaf%ptpp'
            plasmafx%kid(nd+5)='plasmaf%pu'
            plasmafx%kunit(nd+1)='10^20/m^3'
            plasmafx%kunit(nd+2)='eV'
            plasmafx%kunit(nd+3)='eV'
            plasmafx%kunit(nd+4)='eV'
            plasmafx%kunit(nd+5)='m/s'
         enddo
         plasmafx%kid(plasmafx%ndmax)='plasmaf%qinv'
         plasmafx%kunit(plasmafx%ndmax)=' '
         plasmafx%status=1
      endif
c
      plasmafx%time = plasmaf_in%time
      do nr=1,plasmaf_in%nrmax
         plasmafx%rho(nr) = plasmaf_in%rho(nr)
         do ns=1,plasmaf_in%nsmax
            nd=5*(ns-1)
            plasmafx%data(nr,nd+1) = plasmaf_in%data(nr,ns)%pn
            plasmafx%data(nr,nd+2) = plasmaf_in%data(nr,ns)%pt
            plasmafx%data(nr,nd+3) = plasmaf_in%data(nr,ns)%ptpr
            plasmafx%data(nr,nd+4) = plasmaf_in%data(nr,ns)%ptpp
            plasmafx%data(nr,nd+5) = plasmaf_in%data(nr,ns)%pu
         enddo
         plasmafx%data(nr,plasmafx%ndmax) = plasmaf_in%qinv(nr)
      enddo
      CALL DATE_AND_TIME(plasmafx%created_date,
     &                   plasmafx%created_time,
     &                   plasmafx%created_timezone)

      if(plasmafx%status.ge.3) then 
         plasmafx%status=3
      else
         plasmafx%status=2
      endif

      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_set_plasmaf'
         write(6,*) '---- plasmafx%rho'
         write(6,'(1P5E12.4)') 
     &        (plasmafx%rho(nr),nr=1,plasmafx%nrmax)
         do nd=1,plasmafx%ndmax
            write(6,*) '---- ',plasmafx%kid(nd)
            write(6,'(1P5E12.4)') 
     &           (plasmafx%data(nr,nd),nr=1,plasmafx%nrmax)
         enddo
      endif
      return
      end subroutine bpsd_set_plasmaf
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_plasmaf(plasmaf_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_plasmaf_type),intent(out) :: plasmaf_out
      integer,intent(out) :: ierr
      integer :: nr, nd, ns, mode
      real(8) :: s
      real(8), dimension(:), pointer :: v

      if(bpsd_plasmafx_init_flag) call bpsd_plasmafx_init

      if(plasmafx%status.eq.0) then
         write(6,*) 
     &   'XX bpsd_get_plasmaf: no space allocated to plasmafx%data'
         ierr=1
         return
      endif

      if(plasmafx%status.eq.1) then
         write(6,*) 'XX bpsd_get_plasmaf: no data in plasmafx%data'
         ierr=2
         return
      endif

      if(plasmaf_out%nrmax.eq.0) then
         mode=0
         plasmaf_out%nrmax = plasmafx%nrmax
      else
         mode=1
      endif
      plasmaf_out%nsmax = (plasmafx%ndmax-1)/5

      if(associated(plasmaf_out%rho)) then
         if(plasmaf_out%nrmax.ne.size(plasmaf_out%rho)) then
            deallocate(plasmaf_out%rho)
            allocate(plasmaf_out%rho(plasmaf_out%nrmax))
         endif
      else
         allocate(plasmaf_out%rho(plasmaf_out%nrmax))
      endif
      if(associated(plasmaf_out%qinv)) then
         if(plasmaf_out%nrmax.ne.size(plasmaf_out%qinv)) then
            deallocate(plasmaf_out%qinv)
            allocate(plasmaf_out%qinv(plasmaf_out%nrmax))
         endif
      else
         allocate(plasmaf_out%qinv(plasmaf_out%nrmax))
      endif
      if(associated(plasmaf_out%data)) then
         if(plasmaf_out%nrmax.ne.size(plasmaf_out%data,1).or.
     &      plasmaf_out%nsmax.ne.size(plasmaf_out%data,2)) then
            deallocate(plasmaf_out%data)
            allocate(plasmaf_out%data(plasmaf_out%nrmax,
     &                                plasmaf_out%nsmax))
         endif
      else
         allocate(plasmaf_out%data(plasmaf_out%nrmax,
     &                             plasmaf_out%nsmax))
      endif

      if(mode.eq.0) then
         plasmaf_out%time  = plasmafx%time
         do nr=1,plasmafx%nrmax
            plasmaf_out%rho(nr)=plasmafx%rho(nr)
            do ns=1,plasmaf_out%nsmax
               nd=5*(ns-1)
               plasmaf_out%data(nr,ns)%pn  =plasmafx%data(nr,nd+1)
               plasmaf_out%data(nr,ns)%pt  =plasmafx%data(nr,nd+2)
               plasmaf_out%data(nr,ns)%ptpr=plasmafx%data(nr,nd+3)
               plasmaf_out%data(nr,ns)%ptpp=plasmafx%data(nr,nd+4)
               plasmaf_out%data(nr,ns)%pu  =plasmafx%data(nr,nd+5)
            enddo
            plasmaf_out%qinv(nr)=plasmafx%data(nr,plasmafx%ndmax)
         enddo
         ierr=0
         return
      endif
c
      if(plasmafx%status.eq.2) then
         allocate(plasmafx%spline(4,plasmafx%nrmax,plasmafx%ndmax))
         plasmafx%status=3
      endif
c
      if(plasmafx%status.eq.3) then
         if(associated(plasmafx%s)) then
            if(plasmafx%nrmax.ne.size(plasmafx%s,1)) then
               deallocate(plasmafx%s)
               allocate(plasmafx%s(plasmafx%nrmax))
            endif
         else
            allocate(plasmafx%s(plasmafx%nrmax))
         endif
         do nr=1,plasmafx%nrmax
            plasmafx%s(nr)=plasmafx%rho(nr)**2
         enddo
         do nd=1,plasmafx%ndmax
            call bpsd_spl1D(plasmafx,nd,ierr)
         enddo
         plasmafx%status=4
      endif
!
      allocate(v(plasmafx%ndmax))
      do nr=1,plasmaf_out%nrmax
         s = plasmaf_out%rho(nr)**2
         do nd=1,plasmafx%ndmax
            call bpsd_spl1DF(s,v(nd),plasmafx,nd,ierr)
         enddo
         do ns=1,plasmaf_out%nsmax
            nd=5*(ns-1)
            plasmaf_out%data(nr,ns)%pn   = v(nd+1)
            plasmaf_out%data(nr,ns)%pt   = v(nd+2)
            plasmaf_out%data(nr,ns)%ptpr = v(nd+3)
            plasmaf_out%data(nr,ns)%ptpp = v(nd+4)
            plasmaf_out%data(nr,ns)%pu   = v(nd+5)
         enddo
         plasmaf_out%qinv(nr)  = v(plasmafx%ndmax)
      enddo
      deallocate(v)
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_get_plasmaf'
         write(6,*) '---- plasmafx%rho'
         write(6,'(1P5E12.4)') 
     &        (plasmaf_out%rho(nr),nr=1,plasmaf_out%nrmax)
         do ns=1,plasmaf_out%nsmax
            write(6,*) '---- plasmafx%pn(',ns,')'
            write(6,'(1P5E12.4)') 
     &           (plasmaf_out%data(nr,ns)%pn,nr=1,plasmaf_out%nrmax)
            write(6,*) '---- plasmafx%pt(',ns,')'
            write(6,'(1P5E12.4)') 
     &           (plasmaf_out%data(nr,ns)%pt,nr=1,plasmaf_out%nrmax)
            write(6,*) '---- plasmafx%ptpr(',ns,')'
            write(6,'(1P5E12.4)') 
     &           (plasmaf_out%data(nr,ns)%ptpr,nr=1,plasmaf_out%nrmax)
            write(6,*) '---- plasmafx%ptpp(',ns,')'
            write(6,'(1P5E12.4)') 
     &           (plasmaf_out%data(nr,ns)%ptpp,nr=1,plasmaf_out%nrmax)
            write(6,*) '---- plasmafx%pu(',ns,')'
            write(6,'(1P5E12.4)') 
     &           (plasmaf_out%data(nr,ns)%pu,nr=1,plasmaf_out%nrmax)
         enddo
         write(6,*) '---- plasmafx%qinv'
         write(6,'(1P5E12.4)') 
     &        (plasmaf_out%qinv(nr),nr=1,plasmaf_out%nrmax)
      endif
      return
      end subroutine bpsd_get_plasmaf
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_plasmaf(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      integer,intent(in) :: fid
      integer,intent(out) :: ierr
c
      if(bpsd_plasmafx_init_flag) call bpsd_plasmafx_init
c
      if(plasmafx%status.gt.1) 
     &     call bpsd_save_data1Dx(fid,plasmafx,ierr)
      return
c
      end subroutine bpsd_save_plasmaf
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_plasmaf(datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_data1Dx_type),intent(in) :: datax
      integer,intent(out) :: ierr
      integer:: ns,nr,nd
c
      if(bpsd_plasmafx_init_flag) call bpsd_plasmafx_init
c
      if(plasmafx%status.ne.0) then
         if((datax%nrmax.ne.plasmafx%nrmax) .or.
     &      (datax%ndmax.ne.plasmafx%ndmax))  then
            if(plasmafx%status.ge.3) deallocate(plasmafx%spline)
            deallocate(plasmafx%data)
            deallocate(plasmafx%rho)
            deallocate(plasmafx%kunit)
            deallocate(plasmafx%kid)
            plasmafx%status=0
         endif
      endif
c
      if(plasmafx%status.eq.0) then
         allocate(plasmafx%kid(datax%ndmax))
         allocate(plasmafx%kunit(datax%ndmax))
         allocate(plasmafx%rho(datax%nrmax))
         allocate(plasmafx%data(datax%nrmax,datax%ndmax))
         plasmafx%status=1
      endif
c
      plasmafx%dataName=datax%dataName
      plasmafx%ndmax=datax%ndmax
      plasmafx%nrmax=datax%nrmax
      plasmafx%time = datax%time
      do nr=1,plasmafx%nrmax
         plasmafx%rho(nr) = datax%rho(nr)
         do nd=1,plasmafx%ndmax
            plasmafx%data(nr,nd) = datax%data(nr,nd)
         enddo
      enddo
      do nd=1,plasmafx%ndmax
         plasmafx%kid(nd)=datax%kid(nd)
         plasmafx%kunit(nd)=datax%kunit(nd)
      enddo
      plasmafx%created_date = datax%created_date
      plasmafx%created_time = datax%created_time
      plasmafx%created_timezone = datax%created_timezone
      plasmafx%status=2
      ierr=0
      return
c
      end subroutine bpsd_load_plasmaf
c
      end module bpsd_plasmaf
