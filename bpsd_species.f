c     $Id$
c=======================================================================
      module bpsd_species
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_species,bpsd_get_species,
     &       bpsd_save_species,bpsd_load_species
      private
c
      logical, save :: bpsd_speciesx_init_flag = .TRUE.
      type(bpsd_data0Dx_type), save :: speciesx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_speciesx_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      speciesx%status=0
      speciesx%dataName='species'
      speciesx%ndmax=0
c
      bpsd_speciesx_init_flag = .FALSE.
c
      return
      end subroutine bpsd_speciesx_init
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
      if(bpsd_speciesx_init_flag) call bpsd_speciesx_init
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
      if(bpsd_speciesx_init_flag) call bpsd_speciesx_init
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
      subroutine bpsd_save_species(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      integer,intent(in) :: fid
      integer,intent(out) :: ierr
c
      if(bpsd_speciesx_init_flag) call bpsd_speciesx_init
c
      if(speciesx%status.gt.1) 
     &     call bpsd_save_data0Dx(fid,speciesx,ierr)
      return
c
      end subroutine bpsd_save_species
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_species(datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_data0Dx_type),intent(in) :: datax
      integer,intent(out) :: ierr
      integer:: ns,nd
c
      if(bpsd_speciesx_init_flag) call bpsd_speciesx_init
c
      if(speciesx%status.eq.0) then
         speciesx%ndmax=datax%ndmax
         allocate(speciesx%kid(speciesx%ndmax))
         allocate(speciesx%data(speciesx%ndmax))
         do ns=1,speciesx%ndmax/3
            nd=3*(ns-1)
            speciesx%kid(nd+1)='species%pa'
            speciesx%kid(nd+2)='species%pz'
            speciesx%kid(nd+3)='species%pz0'
         enddo
         speciesx%status=1
      endif
c
      speciesx%time = datax%time
      do nd=1,speciesx%ndmax
         speciesx%data(nd) = datax%data(nd)
      enddo
      speciesx%status=2
      ierr=0
      return
c
      end subroutine bpsd_load_species
c
      end module bpsd_species
