c     $Id$
c=======================================================================
      module bpsd_types_internal
c
      use bpsd_kinds
c
      type bpsd_shotx_type
         integer :: status = 1! 1:undef 2:assigned
         integer :: shotID
         integer :: modelID
!         integer :: idum      ! Dummy
         character(len=32) :: dataName
         character(len=32) :: deviceID
      end type bpsd_shotx_type
c
      type bpsd_data0Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: data
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data0Dx_type
c
      type bpsd_data1Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: s 
         real(rkind), dimension(:,:), pointer :: data
         real(rkind), dimension(:,:,:), pointer :: spline
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data1Dx_type
c
      type bpsd_data2Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: nthmax    ! Number of poloidal points
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: s 
         real(rkind), dimension(:), pointer :: th
         real(rkind), dimension(:,:,:), pointer :: data
         real(rkind), dimension(:,:,:,:), pointer :: spline
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data2Dx_type
c
      type bpsd_data3Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: nthmax    ! Number of poloidal points
         integer :: nphmax    ! Number of toroidal points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: s 
         real(rkind), dimension(:), pointer :: th
         real(rkind), dimension(:), pointer :: ph
         real(rkind), dimension(:,:,:,:), pointer :: data
         real(rkind), dimension(:,:,:,:,:,:), pointer :: spline
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data3Dx_type
c
      end module bpsd_types_internal
