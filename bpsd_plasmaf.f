c     $Id$
c=======================================================================
      module bpsd_plasmaf
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
c
      logical, private, save :: bpsd_plasmaf_init_flag = .TRUE.
      type(bpsd_1ddatax_type), private, save :: plasmafx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_plasmaf_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      plasmafx%status=0
      plasmafx%dataName='plasmaf'
      plasmafx%ndmax=0
c
      bpsd_plasmaf_init_flag = .FALSE.
c
      return
      end subroutine bpsd_plasmaf_init
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_plasmaf(plasmaf_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_plasmaf_type):: plasmaf_in
      integer :: ierr
      integer :: ns,nr,nd
c
      if(bpsd_plasmaf_init_flag) call bpsd_plasmaf_init
c
      if(plasmafx%status.ne.0) then
         if((plasmaf_in%nrmax.ne.plasmafx%nrmax) .or.
     &      (plasmaf_in%nsmax*5.ne.plasmafx%ndmax))  then
            if(plasmafx%status.ge.3) deallocate(plasmafx%spline)
            deallocate(plasmafx%data)
            deallocate(plasmafx%s)
            deallocate(plasmafx%kid)
            plasmafx%status=0
         endif
      endif
c
      if(plasmafx%status.eq.0) then
         plasmafx%ndmax=plasmaf_in%nsmax*5+1
         plasmafx%nrmax=plasmaf_in%nrmax
         allocate(plasmafx%kid(plasmafx%ndmax))
         allocate(plasmafx%s(plasmafx%nrmax))
         allocate(plasmafx%data(plasmafx%nrmax,plasmafx%ndmax))
         do ns=1,plasmaf_in%nsmax
            nd=5*(ns-1)
            plasmafx%kid(nd+1)='plasmaf%pn'
            plasmafx%kid(nd+2)='plasmaf%pt'
            plasmafx%kid(nd+3)='plasmaf%ptpr'
            plasmafx%kid(nd+4)='plasmaf%ptpp'
            plasmafx%kid(nd+5)='plasmaf%pu'
         enddo
         plasmafx%kid(plasmafx%ndmax)='plasmaf%qinv'
         plasmafx%status=1
      endif
c
      plasmafx%time = plasmaf_in%time
      do nr=1,plasmaf_in%nrmax
         plasmafx%s(nr) = plasmaf_in%s(nr)
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
      if(plasmafx%status.ge.3) then 
         plasmafx%status=3
      else
         plasmafx%status=2
      endif

      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_set_plasmaf'
         write(6,*) '---- plasmafx%s'
         write(6,'(1P5E12.4)') 
     &        (plasmafx%s(nr),nr=1,plasmafx%nrmax)
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
      type(bpsd_plasmaf_type) :: plasmaf_out
      integer :: ierr
      integer :: nr, nd, ns
      real(8) :: s
      real(8), dimension(:), allocatable :: v
c
      if(bpsd_plasmaf_init_flag) call bpsd_plasmaf_init
c
      if(plasmafx%status.eq.0) then
         write(6,*) 
     &   'XX bpsd_get_plasmaf: no space allocated to plasmafx%data'
         ierr=1
         return
      endif
c
      if(plasmafx%status.eq.1) then
         write(6,*) 'XX bpsd_get_plasmaf: no data in plasmafx%data'
         ierr=1
         return
      endif
c
      if(plasmaf_out%nrmax.eq.0) then
         if(allocated(plasmaf_out%data)) then
            if(plasmafx%nrmax.ne.size(plasmaf_out%data,1)) then
               deallocate(plasmaf_out%qinv)
               deallocate(plasmaf_out%data)
               deallocate(plasmaf_out%s)
               plasmaf_out%nrmax = plasmafx%nrmax
               plasmaf_out%nsmax = (plasmafx%ndmax-1)/5
               allocate(plasmaf_out%s(plasmaf_out%nrmax))
               allocate(plasmaf_out%data(plasmaf_out%nrmax,
     &                                   plasmaf_out%nsmax))
               allocate(plasmaf_out%qinv(plasmaf_out%nrmax))
            endif
         else
            plasmaf_out%nrmax = plasmafx%nrmax
            plasmaf_out%nsmax = (plasmafx%ndmax-1)/5
            allocate(plasmaf_out%s(plasmaf_out%nrmax))
            allocate(plasmaf_out%data(plasmaf_out%nrmax,
     &                                plasmaf_out%nsmax))
            allocate(plasmaf_out%qinv(plasmaf_out%nrmax))
         endif
      endif
c
      if(allocated(plasmaf_out%data)) then
         if(plasmafx%nrmax.le.size(plasmaf_out%data,1)) then
            plasmaf_out%time  = plasmafx%time
            plasmaf_out%nrmax = plasmafx%nrmax
            plasmaf_out%nsmax = (plasmafx%ndmax-1)/5
            do nr=1,plasmafx%nrmax
               plasmaf_out%s(nr)=plasmafx%s(nr)
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
      else
         ierr=3
         return
      endif
c
      if(plasmafx%status.eq.2) then
         allocate(plasmafx%spline(4,plasmafx%nrmax,plasmafx%ndmax))
         plasmafx%status=3
      endif
c
      if(plasmafx%status.eq.3) then
         do nd=1,plasmafx%ndmax
            call spl1D_bpsd(plasmafx,nd,ierr)
         enddo
         plasmafx%status=4
      endif
c
      allocate(v(plasmafx%ndmax))
      do nr=1,plasmaf_out%nrmax
         s = plasmaf_out%s(nr)
         do nd=1,plasmafx%ndmax
            call spl1DF_bpsd(s,v(nd),plasmafx,nd,ierr)
         enddo
         do ns=1,plasmaf_out%nsmax
            nd=5*(ns-1)
            plasmaf_out%data(nr,ns)%pn   = v(nd+1)
            plasmaf_out%data(nr,ns)%pt   = v(nd+2)
            plasmaf_out%data(nr,ns)%ptpr = v(nd+3)
            plasmaf_out%data(nr,ns)%ptpp = v(nd+4)
            plasmaf_out%data(nr,ns)%pu   = v(nd+5)
         enddo
         plasmaf_out%qinv(nr)  = plasmafx%data(nr,plasmafx%ndmax)
      enddo
      deallocate(v)
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_get_plasmaf'
         write(6,*) '---- plasmafx%s'
         write(6,'(1P5E12.4)') 
     &        (plasmaf_out%s(nr),nr=1,plasmaf_out%nrmax)
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
      end module bpsd_plasmaf
