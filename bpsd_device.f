c     $Id$
c=======================================================================
      module bpsd_device
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_device,bpsd_get_device
      private
c
      logical, save :: bpsd_device_init_flag = .TRUE.
      type(bpsd_0ddatax_type), save :: devicex
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_device_init
c-----------------------------------------------------------------------
      use bpsd_subs
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
      bpsd_device_init_flag = .FALSE.
c
      return
      end subroutine bpsd_device_init
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
      if(bpsd_device_init_flag) call bpsd_device_init
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
      subroutine bpsd_get_device(device_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_device_type) :: device_out
      integer :: ierr, nd
c
      if(bpsd_device_init_flag) call bpsd_device_init
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
      end module bpsd_device
