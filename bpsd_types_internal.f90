!     $Id$
!=======================================================================
      module bpsd_types_internal

      use bpsd_kinds

      type bpsd_shotx_type
         integer :: status = 1! 1:undef 2:assigned
         integer :: shotID
         integer :: modelID
         character(len=32) :: dataName
         character(len=32) :: deviceID
         character(len=8)  :: created_date
         character(len=10) :: created_time
         character(len=5)  :: created_timezone
         character(len=9)  :: dummy
      end type bpsd_shotx_type

      type bpsd_data0Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: ndmax     ! Number of data
         real(rkind) :: time
         real(rkind), dimension(:), ALLOCATABLE :: data
         character(len=32) :: dataName
         character(len=8)  :: created_date
         character(len=10) :: created_time
         character(len=5)  :: created_timezone
         character(len=9)  :: dummy
         character(len=32), dimension(:), ALLOCATABLE :: kid
         character(len=32), dimension(:), ALLOCATABLE :: kunit
      end type bpsd_data0Dx_type

      type bpsd_data1Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         real(rkind) :: time
         real(rkind), dimension(:), ALLOCATABLE :: rho
         real(rkind), dimension(:,:), ALLOCATABLE :: data
         real(rkind), dimension(:), ALLOCATABLE :: s
         real(rkind), dimension(:,:,:), ALLOCATABLE :: spline
         character(len=32) :: dataName
         character(len=8)  :: created_date
         character(len=10) :: created_time
         character(len=5)  :: created_timezone
         character(len=9)  :: dummy
         character(len=32), dimension(:), ALLOCATABLE :: kid
         character(len=32), dimension(:), ALLOCATABLE :: kunit
      end type bpsd_data1Dx_type

      type bpsd_data2Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: nthmax    ! Number of poloidal points
         integer :: ndmax     ! Number of data
         real(rkind) :: time
         real(rkind), dimension(:), ALLOCATABLE :: rho
         real(rkind), dimension(:), ALLOCATABLE :: th
         real(rkind), dimension(:,:,:), ALLOCATABLE :: data
         real(rkind), dimension(:), ALLOCATABLE :: s
         real(rkind), dimension(:,:,:,:), ALLOCATABLE :: spline
         character(len=32) :: dataName
         character(len=8)  :: created_date
         character(len=10) :: created_time
         character(len=5)  :: created_timezone
         character(len=9)  :: dummy
         character(len=32), dimension(:), ALLOCATABLE :: kid
         character(len=32), dimension(:), ALLOCATABLE :: kunit
      end type bpsd_data2Dx_type

      type bpsd_data3Dx_type
         integer :: status = 0! 0:unalloc 1:undef 2:assigned 
                              ! 3:sp-alloc 4:splined
         integer :: nrmax     ! Number of radial points
         integer :: nthmax    ! Number of poloidal points
         integer :: nphmax    ! Number of toroidal points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         real(rkind) :: time
         real(rkind), dimension(:), ALLOCATABLE :: rho
         real(rkind), dimension(:), ALLOCATABLE :: th
         real(rkind), dimension(:), ALLOCATABLE :: ph
         real(rkind), dimension(:,:,:,:), ALLOCATABLE :: data
         real(rkind), dimension(:), ALLOCATABLE :: s
         real(rkind), dimension(:,:,:,:,:,:), ALLOCATABLE :: spline
         character(len=32) :: dataName
         character(len=8)  :: created_date
         character(len=10) :: created_time
         character(len=5)  :: created_timezone
         character(len=9)  :: dummy
         character(len=32), dimension(:), ALLOCATABLE :: kid
         character(len=32), dimension(:), ALLOCATABLE :: kunit
      end type bpsd_data3Dx_type

      end module bpsd_types_internal
