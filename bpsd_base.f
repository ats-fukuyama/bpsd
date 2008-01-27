c     $Id$
c=======================================================================
      module bpsd
c
      use bpsd_kind
      use bpsd_constant
      use bpsd_flags
      use bpsd_types
      use bpsd_plasmaf
c
      logical, private, save :: bpsd_init_flag = .TRUE.
      type(bpsd_shotx_type), private, save :: shot
      type(bpsd_0ddatax_type), private, save :: devicex
      type(bpsd_0ddatax_type), private, save :: speciesx
      type(bpsd_1ddatax_type), private, save :: equ1Dx
      type(bpsd_1ddatax_type), private, save :: metric1Dx
c      type(bpsd_1ddatax_type), private, save :: plasmafx
c
      interface bpsd_set_data
         module procedure bpsd_set_shot, 
     &                    bpsd_set_device, 
     &                    bpsd_set_species, 
     &                    bpsd_set_equ1D, 
     &                    bpsd_set_metric1D, 
     &                    bpsd_set_plasmaf
      end interface bpsd_set_data
c
      interface bpsd_get_data
         module procedure bpsd_get_shot, 
     &                    bpsd_get_device, 
     &                    bpsd_get_species, 
     &                    bpsd_get_equ1D, 
     &                    bpsd_get_metric1D, 
     &                    bpsd_get_plasmaf
      end interface bpsd_get_data
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_init
c-----------------------------------------------------------------------
      implicit none
c
      devicex%status=0
      devicex%dataName='device'
      devicex%ndmax=8
      allocate(devicex%kid(8))
      allocate(devicex%data(8))
      devicex%kid(1)='device%rr'
      devicex%kid(2)='device%zz'
      devicex%kid(3)='device%ra'
      devicex%kid(4)='device%rb'
      devicex%kid(5)='device%bb'
      devicex%kid(6)='device%ip'
      devicex%kid(7)='device%elip'
      devicex%kid(8)='device%trig'
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
      speciesx%status=0
      speciesx%dataName='species'
      speciesx%ndmax=0
c      plasmafx%status=0
c      plasmafx%dataName='plasmaf'
c      plasmafx%ndmax=0
c
      bpsd_init_flag = .FALSE.
c
      return
      end subroutine bpsd_init
c-----------------------------------------------------------------------
      subroutine bpsd_set_shot(shot_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_shot_type) :: shot_in
      integer :: ierr
c
      shot%deviceID = shot_in%deviceID
      shot%dataName = 'shot'
      shot%shotID = shot_in%shotID
      shot%modelID = shot_in%modelID
      shot%status = 2
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,'(A/A32,A32/A32,I12/A32,I12)')
     &        '-- bpsd_set_shot',
     &        'shot%deviceID: ',shot%deviceID,
     &        'shot%shotID  : ',shot%shotID,
     &        'shot%modelID : ',shot%modelID
      endif
      return
      end subroutine bpsd_set_shot
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_device(device_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_device_type) :: device_in
      integer :: ierr, nd
c
      if(bpsd_init_flag) call bpsd_init
c
      devicex%dataName = 'device'
      devicex%time = 0.D0
      devicex%data(1) = device_in%rr
      devicex%data(2) = device_in%zz
      devicex%data(3) = device_in%ra
      devicex%data(4) = device_in%rb
      devicex%data(5) = device_in%bb
      devicex%data(6) = device_in%ip
      devicex%data(7) = device_in%elip
      devicex%data(8) = device_in%trig
      devicex%status = 2
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_set_device'
         do nd=1,devicex%ndmax
            write(6,'(A32,1PE12.4)') 
     &           devicex%kid(nd),devicex%data(nd)
         enddo
      endif
      return
      end subroutine bpsd_set_device
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_species(species_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_species_type):: species_in
      integer :: ierr
      integer :: ns, nd
c
      if(bpsd_init_flag) call bpsd_init
c
      if(speciesx%status.ne.0) then
         if(species_in%nsmax*3.ne.speciesx%ndmax) then
            deallocate(speciesx%data)
            deallocate(speciesx%kid)
            speciesx%status=0
         endif
      endif
c
      if(speciesx%status.eq.0) then
         speciesx%ndmax=species_in%nsmax*3
         allocate(speciesx%kid(speciesx%ndmax))
         allocate(speciesx%data(speciesx%ndmax))
         do ns=1,species_in%nsmax
            nd=3*(ns-1)
            speciesx%kid(nd+1)='species%pa'
            speciesx%kid(nd+2)='species%pz'
            speciesx%kid(nd+3)='species%pz0'
         enddo
         speciesx%status=1
      endif
c
      speciesx%time=0.D0
      do ns=1,species_in%nsmax
         nd=3*(ns-1)
         speciesx%data(nd+1)=species_in%data(ns)%pa
         speciesx%data(nd+2)=species_in%data(ns)%pz
         speciesx%data(nd+3)=species_in%data(ns)%pz0
      enddo
      speciesx%status = 2
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_set_species'
         do nd=1,speciesx%ndmax
            write(6,'(A32,1PE12.4)') 
     &           speciesx%kid(nd),speciesx%data(nd)
         enddo
      endif
      return
      end subroutine bpsd_set_species
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
      if(bpsd_init_flag) call bpsd_init
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
      subroutine bpsd_set_metric1D(metric1D_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_metric1D_type):: metric1D_in
      integer :: ierr
      integer :: nr,nd
c
      if(bpsd_init_flag) call bpsd_init
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
         allocate(metric1Dx%data(metric1Dx%nrmax,15))
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
         write(6,*) '-- bpsd_set_metric1d'
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
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$      subroutine bpsd_set_plasmaf(plasmaf_in,ierr)
c$$$c-----------------------------------------------------------------------
c$$$c
c$$$      implicit none
c$$$      type(bpsd_plasmaf_type):: plasmaf_in
c$$$      integer :: ierr
c$$$      integer :: ns,nr,nd
c$$$c
c$$$      if(bpsd_init_flag) call bpsd_init
c$$$c
c$$$      if(plasmafx%status.ne.0) then
c$$$         if((plasmaf_in%nrmax.ne.plasmafx%nrmax) .or.
c$$$     &      (plasmaf_in%nsmax*5.ne.plasmafx%ndmax))  then
c$$$            if(plasmafx%status.ge.3) deallocate(plasmafx%spline)
c$$$            deallocate(plasmafx%data)
c$$$            deallocate(plasmafx%s)
c$$$            deallocate(plasmafx%kid)
c$$$            plasmafx%status=0
c$$$         endif
c$$$      endif
c$$$c
c$$$      if(plasmafx%status.eq.0) then
c$$$         plasmafx%ndmax=plasmaf_in%nsmax*5+1
c$$$         plasmafx%nrmax=plasmaf_in%nrmax
c$$$         allocate(plasmafx%kid(plasmafx%ndmax))
c$$$         allocate(plasmafx%s(plasmafx%nrmax))
c$$$         allocate(plasmafx%data(plasmafx%nrmax,plasmafx%ndmax))
c$$$         do ns=1,plasmaf_in%nsmax
c$$$            nd=5*(ns-1)
c$$$            plasmafx%kid(nd+1)='plasmaf%pn'
c$$$            plasmafx%kid(nd+2)='plasmaf%pt'
c$$$            plasmafx%kid(nd+3)='plasmaf%ptpr'
c$$$            plasmafx%kid(nd+4)='plasmaf%ptpp'
c$$$            plasmafx%kid(nd+5)='plasmaf%pu'
c$$$         enddo
c$$$         plasmafx%kid(plasmafx%ndmax)='plasmaf%qinv'
c$$$         plasmafx%status=1
c$$$      endif
c$$$c
c$$$      plasmafx%time = plasmaf_in%time
c$$$      do nr=1,plasmaf_in%nrmax
c$$$         plasmafx%s(nr) = plasmaf_in%s(nr)
c$$$         do ns=1,plasmaf_in%nsmax
c$$$            nd=5*(ns-1)
c$$$            plasmafx%data(nr,nd+1) = plasmaf_in%data(nr,ns)%pn
c$$$            plasmafx%data(nr,nd+2) = plasmaf_in%data(nr,ns)%pt
c$$$            plasmafx%data(nr,nd+3) = plasmaf_in%data(nr,ns)%ptpr
c$$$            plasmafx%data(nr,nd+4) = plasmaf_in%data(nr,ns)%ptpp
c$$$            plasmafx%data(nr,nd+5) = plasmaf_in%data(nr,ns)%pu
c$$$         enddo
c$$$         plasmafx%data(nr,plasmafx%ndmax) = plasmaf_in%qinv(nr)
c$$$      enddo
c$$$      if(plasmafx%status.ge.3) then 
c$$$         plasmafx%status=3
c$$$      else
c$$$         plasmafx%status=2
c$$$      endif
c$$$
c$$$      ierr = 0
c$$$c
c$$$      if(bpsd_debug_flag) then
c$$$         write(6,*) '-- bpsd_set_plasmaf'
c$$$         write(6,*) '---- plasmafx%s'
c$$$         write(6,'(1p5e12.4)') 
c$$$     &        (plasmafx%s(nr),nr=1,plasmafx%nrmax)
c$$$         do nd=1,plasmafx%ndmax
c$$$            write(6,*) '---- ',plasmafx%kid(nd)
c$$$            write(6,'(1p5e12.4)') 
c$$$     &           (plasmafx%data(nr,nd),nr=1,plasmafx%nrmax)
c$$$         enddo
c$$$      endif
c$$$      return
c$$$      end subroutine bpsd_set_plasmaf
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_shot(shot_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_shot_type) :: shot_out
      integer :: ierr
c
      if(shot%status.eq.1) then
         write(6,*) 'XX bpsd_get_shot: no data in shot'
         ierr=2
         return
      endif
c
      shot_out%deviceID = shot%deviceID
      shot_out%shotID = shot%shotID
      shot_out%modelID = shot%modelID
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,'(A/A32,A32/A32,I12/A32,I12)')
     &        '-- bpsd_get_shot',
     &        'shot%deviceID: ',shot%deviceID,
     &        'shot%shotID  : ',shot%shotID,
     &        'shot%modelID : ',shot%modelID
      endif
      return
      end subroutine bpsd_get_shot
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_device(device_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_device_type) :: device_out
      integer :: ierr, nd
c
      if(devicex%status.eq.1) then
         write(6,*) 'XX bpsd_get_device: no data in device'
         ierr=2
         return
      endif
c
      device_out%rr    = devicex%data(1)
      device_out%zz    = devicex%data(2)
      device_out%ra    = devicex%data(3)
      device_out%rb    = devicex%data(4)
      device_out%bb    = devicex%data(5)
      device_out%ip    = devicex%data(6)
      device_out%elip  = devicex%data(7)
      device_out%trig  = devicex%data(8)
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_get_device'
         do nd=1,devicex%ndmax
            write(6,'(A32,1PE12.4)') 
     &           devicex%kid(nd),devicex%data(nd)
         enddo
      endif
      return
      end subroutine bpsd_get_device
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_species(species_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_species_type) :: species_out
      integer :: ierr
      integer :: nr, nd, ns
      real(8) :: s
      real(8), dimension(:), allocatable :: v
c
      if(speciesx%status.eq.0) then
         write(6,*) 
     &   'XX bpsd_get_species: no space allocated to species%data'
         ierr=1
         return
      endif
c
      if(speciesx%status.eq.1) then
         write(6,*) 'XX bpsd_get_species: no data in species%data'
         ierr=2
         return
      endif
c
      if(species_out%nsmax.eq.0) then
         if(allocated(species_out%data)) then
            if(speciesx%ndmax.ne.size(species_out%data,1)*3) then
               deallocate(species_out%data)
               species_out%nsmax = speciesx%ndmax/3
               allocate(species_out%data(species_out%nsmax))
            endif
         else
            species_out%nsmax = speciesx%ndmax/3
            allocate(species_out%data(species_out%nsmax))
         endif
      endif
c
      if(allocated(species_out%data)) then
         if(speciesx%ndmax.le.size(species_out%data,1)*3) then
            species_out%nsmax = speciesx%ndmax/3
            do ns=1,species_out%nsmax
               nd=3*(ns-1)
               species_out%data(ns)%pa =speciesx%data(nd+1)
               species_out%data(ns)%pz =speciesx%data(nd+2)
               species_out%data(ns)%pz0=speciesx%data(nd+3)
            enddo
         endif
      else
         ierr=3
         return
      endif
c
      if(bpsd_debug_flag) then
         write(6,*) '-- bpsd_get_species'
         do nd=1,speciesx%ndmax
            write(6,'(A32,1PE12.4)') 
     &           speciesx%kid(nd),speciesx%data(nd)
         enddo
      endif
      return
      end subroutine bpsd_get_species
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
c-----------------------------------------------------------------------
      subroutine bpsd_get_metric1D(metric1D_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
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
            call spl1D_bpsd(metric1Dx,nd,ierr)
         enddo
         metric1Dx%status=4
      endif
c
      do nr=1,metric1D_out%nrmax
         s = metric1D_out%s(nr)
         do nd=1,metric1Dx%ndmax
            call spl1DF_bpsd(s,v(nd),metric1Dx,nd,ierr)
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
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$      subroutine bpsd_get_plasmaf(plasmaf_out,ierr)
c$$$c-----------------------------------------------------------------------
c$$$c
c$$$      implicit none
c$$$      type(bpsd_plasmaf_type) :: plasmaf_out
c$$$      integer :: ierr
c$$$      integer :: nr, nd, ns
c$$$      real(8) :: s
c$$$      real(8), dimension(:), allocatable :: v
c$$$c
c$$$      if(plasmafx%status.eq.0) then
c$$$         write(6,*) 
c$$$     &   'XX bpsd_get_plasmaf: no space allocated to plasmafx%data'
c$$$         ierr=1
c$$$         return
c$$$      endif
c$$$c
c$$$      if(plasmafx%status.eq.1) then
c$$$         write(6,*) 'XX bpsd_get_plasmaf: no data in plasmafx%data'
c$$$         ierr=1
c$$$         return
c$$$      endif
c$$$c
c$$$      if(plasmaf_out%nrmax.eq.0) then
c$$$         if(allocated(plasmaf_out%data)) then
c$$$            if(plasmafx%nrmax.ne.size(plasmaf_out%data,1)) then
c$$$               deallocate(plasmaf_out%qinv)
c$$$               deallocate(plasmaf_out%data)
c$$$               deallocate(plasmaf_out%s)
c$$$               plasmaf_out%nrmax = plasmafx%nrmax
c$$$               plasmaf_out%nsmax = (plasmafx%ndmax-1)/5
c$$$               allocate(plasmaf_out%s(plasmaf_out%nrmax))
c$$$               allocate(plasmaf_out%data(plasmaf_out%nrmax,
c$$$     &                                   plasmaf_out%nsmax))
c$$$               allocate(plasmaf_out%qinv(plasmaf_out%nrmax))
c$$$            endif
c$$$         else
c$$$            plasmaf_out%nrmax = plasmafx%nrmax
c$$$            plasmaf_out%nsmax = (plasmafx%ndmax-1)/5
c$$$            allocate(plasmaf_out%s(plasmaf_out%nrmax))
c$$$            allocate(plasmaf_out%data(plasmaf_out%nrmax,
c$$$     &                                plasmaf_out%nsmax))
c$$$            allocate(plasmaf_out%qinv(plasmaf_out%nrmax))
c$$$         endif
c$$$      endif
c$$$c
c$$$      if(allocated(plasmaf_out%data)) then
c$$$         if(plasmafx%nrmax.le.size(plasmaf_out%data,1)) then
c$$$            plasmaf_out%time  = plasmafx%time
c$$$            plasmaf_out%nrmax = plasmafx%nrmax
c$$$            plasmaf_out%nsmax = (plasmafx%ndmax-1)/5
c$$$            do nr=1,plasmafx%nrmax
c$$$               plasmaf_out%s(nr)=plasmafx%s(nr)
c$$$               do ns=1,plasmaf_out%nsmax
c$$$                  nd=5*(ns-1)
c$$$                  plasmaf_out%data(nr,ns)%pn  =plasmafx%data(nr,nd+1)
c$$$                  plasmaf_out%data(nr,ns)%pt  =plasmafx%data(nr,nd+2)
c$$$                  plasmaf_out%data(nr,ns)%ptpr=plasmafx%data(nr,nd+3)
c$$$                  plasmaf_out%data(nr,ns)%ptpp=plasmafx%data(nr,nd+4)
c$$$                  plasmaf_out%data(nr,ns)%pu  =plasmafx%data(nr,nd+5)
c$$$               enddo
c$$$               plasmaf_out%qinv(nr)=plasmafx%data(nr,plasmafx%ndmax)
c$$$            enddo
c$$$            ierr=0
c$$$            return
c$$$         endif
c$$$      else
c$$$         ierr=3
c$$$         return
c$$$      endif
c$$$c
c$$$      if(plasmafx%status.eq.2) then
c$$$         allocate(plasmafx%spline(4,plasmafx%nrmax,plasmafx%ndmax))
c$$$         plasmafx%status=3
c$$$      endif
c$$$c
c$$$      if(plasmafx%status.eq.3) then
c$$$         do nd=1,plasmafx%ndmax
c$$$            call spl1D_bpsd(plasmafx,nd,ierr)
c$$$         enddo
c$$$         plasmafx%status=4
c$$$      endif
c$$$c
c$$$      allocate(v(plasmafx%ndmax))
c$$$      do nr=1,plasmaf_out%nrmax
c$$$         s = plasmaf_out%s(nr)
c$$$         do nd=1,plasmafx%ndmax
c$$$            call spl1DF_bpsd(s,v(nd),plasmafx,nd,ierr)
c$$$         enddo
c$$$         do ns=1,plasmaf_out%nsmax
c$$$            nd=5*(ns-1)
c$$$            plasmaf_out%data(nr,ns)%pn   = v(nd+1)
c$$$            plasmaf_out%data(nr,ns)%pt   = v(nd+2)
c$$$            plasmaf_out%data(nr,ns)%ptpr = v(nd+3)
c$$$            plasmaf_out%data(nr,ns)%ptpp = v(nd+4)
c$$$            plasmaf_out%data(nr,ns)%pu   = v(nd+5)
c$$$         enddo
c$$$         plasmaf_out%qinv(nr)  = plasmafx%data(nr,plasmafx%ndmax)
c$$$      enddo
c$$$      deallocate(v)
c$$$      ierr = 0
c$$$c
c$$$      if(bpsd_debug_flag) then
c$$$         write(6,*) '-- bpsd_get_plasmaf'
c$$$         write(6,*) '---- plasmafx%s'
c$$$         write(6,'(1P5E12.4)') 
c$$$     &        (plasmaf_out%s(nr),nr=1,plasmaf_out%nrmax)
c$$$         do ns=1,plasmaf_out%nsmax
c$$$            write(6,*) '---- plasmafx%pn(',ns,')'
c$$$            write(6,'(1P5E12.4)') 
c$$$     &           (plasmaf_out%data(nr,ns)%pn,nr=1,plasmaf_out%nrmax)
c$$$            write(6,*) '---- plasmafx%pt(',ns,')'
c$$$            write(6,'(1P5E12.4)') 
c$$$     &           (plasmaf_out%data(nr,ns)%pt,nr=1,plasmaf_out%nrmax)
c$$$            write(6,*) '---- plasmafx%ptpr(',ns,')'
c$$$            write(6,'(1P5E12.4)') 
c$$$     &           (plasmaf_out%data(nr,ns)%ptpr,nr=1,plasmaf_out%nrmax)
c$$$            write(6,*) '---- plasmafx%ptpp(',ns,')'
c$$$            write(6,'(1P5E12.4)') 
c$$$     &           (plasmaf_out%data(nr,ns)%ptpp,nr=1,plasmaf_out%nrmax)
c$$$            write(6,*) '---- plasmafx%pu(',ns,')'
c$$$            write(6,'(1P5E12.4)') 
c$$$     &           (plasmaf_out%data(nr,ns)%pu,nr=1,plasmaf_out%nrmax)
c$$$         enddo
c$$$         write(6,*) '---- plasmafx%qinv'
c$$$         write(6,'(1P5E12.4)') 
c$$$     &        (plasmaf_out%qinv(nr),nr=1,plasmaf_out%nrmax)
c$$$      endif
c$$$      return
c$$$      end subroutine bpsd_get_plasmaf
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$      subroutine spl1D_bpsd(data1D,nd,ierr)
c$$$c-----------------------------------------------------------------------
c$$$c
c$$$!      use libspl_mod
c$$$      implicit none
c$$$      type(bpsd_1ddatax_type) :: data1D
c$$$      integer :: nd     ! position of dependent variable
c$$$      integer :: ierr    ! error indicator
c$$$      real(8), dimension(:), allocatable :: deriv
c$$$c
c$$$      allocate(deriv(data1D%nrmax))
c$$$      call spl1D(data1D%s,data1D%data(1,nd),deriv,data1D%spline(1,1,nd),
c$$$     &           data1D%nrmax,0,ierr)
c$$$      if(ierr.ne.0) 
c$$$     &     write(6,*) 'XX spl1D_bpsd : '//data1D%kid(nd)//
c$$$     &                ': ierr=',ierr
c$$$      deallocate(deriv)
c$$$      return
c$$$      end subroutine spl1D_bpsd
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$      subroutine spl1DF_bpsd(pos,val,data1D,nd,ierr)
c$$$c-----------------------------------------------------------------------
c$$$c
c$$$!      use libspl_mod
c$$$      implicit none
c$$$      real(8) :: pos     ! value of independent variable
c$$$      real(8) :: val     ! value of dependent variable
c$$$      type(bpsd_1ddatax_type) :: data1D
c$$$      integer :: nd      ! position of dependent variable
c$$$      integer :: ierr    ! error indicator
c$$$
c$$$      call spl1DF(pos,val,
c$$$     &            data1D%s,data1D%spline(1,1,nd),data1D%nrmax,ierr)
c$$$      if(ierr.ne.0) then
c$$$         write(6,*) 'XX spl1DF_bpsd : '//data1D%kid(nd)//
c$$$     &                ': ierr=',ierr
c$$$         write(6,'(1P3E12.4)')  
c$$$     &        pos,data1D%s(1),data1D%s(data1D%nrmax)
c$$$      endif
c$$$      end subroutine spl1DF_bpsd
c
      end module bpsd
