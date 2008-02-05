c     $Id$
c=======================================================================
      module bpsd_equ1d
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_equ1d,bpsd_get_equ1d
      private
c
      logical, save :: bpsd_equ1d_init_flag = .TRUE.
      type(bpsd_1ddatax_type), save :: equ1dx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_equ1d_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      equ1Dx%status=0
      equ1Dx%ndmax=6
      equ1Dx%dataName='equ1D'
      allocate(equ1Dx%kid(6))
      equ1Dx%kid(1)='equ1D%psit'
      equ1Dx%kid(2)='equ1D%psip'
      equ1Dx%kid(3)='equ1D%ppp'
      equ1Dx%kid(4)='equ1D%piq'
      equ1Dx%kid(5)='equ1D%pip'
      equ1Dx%kid(6)='equ1D%pit'
c
      bpsd_equ1d_init_flag = .FALSE.
c
      return
      end subroutine bpsd_equ1d_init
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_equ1D(equ1D_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_equ1D_type):: equ1D_in
      integer :: ierr
      integer :: nr,nd
c
      if(bpsd_equ1d_init_flag) call bpsd_equ1d_init
c
      if(equ1Dx%status.ne.0) then
         if(equ1D_in%nrmax.ne.equ1Dx%nrmax) then
            if(equ1Dx%status.ge.3) deallocate(equ1Dx%spline)
            deallocate(equ1Dx%data)
            deallocate(equ1Dx%s)
            equ1Dx%status=0
         endif
      endif
c
      if(equ1Dx%status.eq.0) then
         equ1Dx%nrmax=equ1D_in%nrmax
         allocate(equ1Dx%s(equ1Dx%nrmax))
         allocate(equ1Dx%data(equ1Dx%nrmax,6))
         equ1Dx%status=1
      endif
c
      equ1Dx%time=equ1D_in%time
      do nr=1,equ1Dx%nrmax
         equ1Dx%s(nr) = equ1D_in%s(nr)
         equ1Dx%data(nr,1) = equ1D_in%data(nr)%psit
         equ1Dx%data(nr,2) = equ1D_in%data(nr)%psip
         equ1Dx%data(nr,3) = equ1D_in%data(nr)%ppp
         equ1Dx%data(nr,4) = equ1D_in%data(nr)%piq
         equ1Dx%data(nr,5) = equ1D_in%data(nr)%pip
         equ1Dx%data(nr,6) = equ1D_in%data(nr)%pit
      enddo
      if(equ1Dx%status.ge.3) then
         equ1Dx%status=3
      else
         equ1Dx%status=2
      endif
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_set_equ1D'
         write(6,*) '---- equ1Dx%s'
         write(6,'(1P5E12.4)') 
     &        (equ1Dx%s(nr),nr=1,equ1Dx%nrmax)
         do nd=1,equ1Dx%ndmax
            write(6,*) '---- ',equ1Dx%kid(nd)
            write(6,'(1P5E12.4)') 
     &           (equ1Dx%data(nr,nd),nr=1,equ1Dx%nrmax)
         enddo
      endif
      return
      end subroutine bpsd_set_equ1D
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_equ1D(equ1D_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_equ1D_type) :: equ1D_out
      integer :: ierr
      integer :: nr, nd
      real(8) :: s
      real(8), dimension(6) :: v
c
      if(bpsd_equ1d_init_flag) call bpsd_equ1d_init
c
      if(equ1Dx%status.eq.0) then
         write(6,*) 
     &        'XX bpsd_get_equ1D: no space allocated to equ1Dx%data'
         ierr=1
         return
      endif
c
      if(equ1Dx%status.eq.1) then
         write(6,*) 'XX bpsd_get_equ1D: no data in equ1Dx%data'
         ierr=2
         return
      endif
c
      if(equ1D_out%nrmax.eq.0) then
         if(allocated(equ1D_out%data)) then
            if(equ1Dx%nrmax.ne.size(equ1D_out%data,1)) then
               deallocate(equ1D_out%data)
               equ1D_out%nrmax = equ1Dx%nrmax
               allocate(equ1D_out%data(equ1D_out%nrmax))
            endif
         else
            equ1D_out%nrmax = equ1Dx%nrmax
            allocate(equ1D_out%data(equ1D_out%nrmax))
         endif
      endif
c
      if(allocated(equ1D_out%data)) then
         if(equ1Dx%nrmax.le.size(equ1D_out%data,1)) then
            equ1D_out%time  = equ1Dx%time
            equ1D_out%nrmax = equ1Dx%nrmax
            do nr=1,equ1D_out%nrmax
               equ1D_out%s(nr)         = equ1Dx%s(nr)
               equ1D_out%data(nr)%psit = equ1Dx%data(nr,1)
               equ1D_out%data(nr)%psip = equ1Dx%data(nr,2)
               equ1D_out%data(nr)%ppp  = equ1Dx%data(nr,3)
               equ1D_out%data(nr)%piq  = equ1Dx%data(nr,4)
               equ1D_out%data(nr)%pip  = equ1Dx%data(nr,5)
               equ1D_out%data(nr)%pit  = equ1Dx%data(nr,6)
            enddo
            ierr=0
            return
         endif
      else
         ierr=3
         return
      endif
c
      if(equ1Dx%status.eq.2) then
         allocate(equ1Dx%spline(4,equ1Dx%nrmax,6))
         equ1Dx%status=3
      endif
c
      if(equ1Dx%status.eq.3) then
         do nd=1,6
            call spl1D_bpsd(equ1Dx,nd,ierr)
         enddo
         equ1Dx%status=4
      endif
c
      do nr=1,equ1D_out%nrmax
         s = equ1D_out%s(nr)
         do nd=1,6
            call spl1DF_bpsd(s,v(nd),equ1Dx,nd,ierr)
         enddo
         equ1D_out%data(nr)%psit = v(1)
         equ1D_out%data(nr)%psip = v(2)
         equ1D_out%data(nr)%ppp  = v(3)
         equ1D_out%data(nr)%piq  = v(4)
         equ1D_out%data(nr)%pip  = v(5)
         equ1D_out%data(nr)%pit  = v(6)
      enddo
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_get_equ1D'
         write(6,*) '---- equ1Dx%s'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%s(nr),nr=1,equ1D_out%nrmax)
         write(6,*) '---- equ1Dx%psit'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%data(nr)%psit,nr=1,equ1D_out%nrmax)
         write(6,*) '---- equ1Dx%psip'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%data(nr)%psip,nr=1,equ1D_out%nrmax)
         write(6,*) '---- equ1Dx%ppp'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%data(nr)%ppp,nr=1,equ1D_out%nrmax)
         write(6,*) '---- equ1Dx%piq'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%data(nr)%piq,nr=1,equ1D_out%nrmax)
         write(6,*) '---- equ1Dx%pip'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%data(nr)%pip,nr=1,equ1D_out%nrmax)
         write(6,*) '---- equ1Dx%pit'
         write(6,'(1P5E12.4)') 
     &        (equ1D_out%data(nr)%pit,nr=1,equ1D_out%nrmax)
      endif
      return
      end subroutine bpsd_get_equ1D
c
      end module bpsd_equ1d
