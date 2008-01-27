c     $Id$
c=======================================================================
      module bpsd_types_internal
c
      use bpsd_kind
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
      type bpsd_0ddatax_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: data
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_0ddatax_type
c
      type bpsd_1ddatax_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
         real(rkind), dimension(:,:), allocatable :: data
         real(rkind), dimension(:,:,:), allocatable :: spline
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_1ddatax_type
c
      type bpsd_2ddatax_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: nthmax    ! Number of poloidal points
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
         real(rkind), dimension(:), allocatable :: th
         real(rkind), dimension(:,:,:), allocatable :: data
         real(rkind), dimension(:,:,:,:), allocatable :: spline
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_2ddatax_type
c
      type bpsd_3ddatax_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: nthmax    ! Number of poloidal points
         integer :: nphmax    ! Number of toroidal points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
         real(rkind), dimension(:), allocatable :: th
         real(rkind), dimension(:), allocatable :: ph
         real(rkind), dimension(:,:,:,:), allocatable :: data
         real(rkind), dimension(:,:,:,:,:,:), allocatable :: spline
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_3ddatax_type
c
      end module bpsd_types_internal
