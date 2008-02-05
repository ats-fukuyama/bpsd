c     $Id$
c=======================================================================
      module bpsd_types
c
      use bpsd_kind
c
      type bpsd_shot_type
         integer :: shotID
         integer :: modelID
         character(len=32) :: deviceID
      end type bpsd_shot_type
c
      type bpsd_device_type
         real(rkind) :: rr     ! Geometrical major radius [m]
         real(rkind) :: zz     ! Geometrical vetical position [m]
         real(rkind) :: ra     ! Typical minor radius (Rmax-Rmin)/2 [m]
         real(rkind) :: rb     ! Typical wall radius [m]
         real(rkind) :: bb     ! Vacuum toroidal magnetic field at rr [T]
         real(rkind) :: ip     ! Typical plasma current [A]
         real(rkind) :: elip   ! Typical ellipticity
         real(rkind) :: trig   ! Typical triangularity
      end type bpsd_device_type
c
      type bpsd_species_data
         real(rkind) :: pa     ! Mass number (n. of protons + n. of neutrons)
         real(rkind) :: pz     ! Charge number (n. of protons - n. of electrons)
         real(rkind) :: pz0    ! Atomic number (n. of protons)
      end type bpsd_species_data

      type bpsd_species_type
         integer :: nsmax     ! Number of particle species
         integer :: idum      ! Dummy
         type(bpsd_species_data), dimension(:), allocatable :: data
      end type bpsd_species_type
c
      type bpsd_equ1D_data
         real(rkind) :: psit   ! Toroidal magnetic flux [Wb] ~pi*r^2*B
         real(rkind) :: psip   ! Poloidal magnetic flux [Wb] ~pi*R*r*Bp
         real(rkind) :: ppp    ! Plasma pressure [Pa]
         real(rkind) :: piq    ! Inverse of safety factor, iota
         real(rkind) :: pip    ! Poloidal current [A] ~2*pi*R*B/mu_0
         real(rkind) :: pit    ! Toroidal current [A] ~2*pi*r*Bp/mu_0
      end type bpsd_equ1D_data
      type bpsd_equ1D_type
         integer :: nrmax     ! Number of radial points
         integer :: idum      ! Dummy
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
                              ! (rho^2) normarized toroidal magnetic flux
         type(bpsd_equ1D_data), dimension(:), allocatable :: data
      end type bpsd_equ1D_type
c
      type bpsd_metric1D_data
         real(rkind) :: pvol     ! Plasma volude [m^3] ~2*pi*R*pi*r^2
         real(rkind) :: psur     ! Plasma surface [m^2] ~pi*r^2
         real(rkind) :: dvpsit   ! dV/dPsit
         real(rkind) :: dvpsip   ! dV/dPsip
         real(rkind) :: aver2    ! <R^2>
         real(rkind) :: aver2i   ! <1/R^2>
         real(rkind) :: aveb2    ! <B^2>
         real(rkind) :: aveb2i   ! <1/B^2>
         real(rkind) :: avegv2   ! <|gradV|^2>
         real(rkind) :: avegvr2  ! <|gradV|^2/R^2>
         real(rkind) :: avegpp2  ! <|gradPsip|^2>
         real(rkind) :: rr       ! R
         real(rkind) :: rs       ! r
         real(rkind) :: elip     ! elipticity
         real(rkind) :: trig     ! triangularity
      end type bpsd_metric1D_data
c
      type bpsd_metric1D_type
         integer :: nrmax       ! Number of radial points
         integer :: idum      ! Dummy
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
                                ! (rho^2) normarized toroidal magnetic flux
         type(bpsd_metric1D_data), dimension(:), allocatable :: data
      end type bpsd_metric1D_type
c
      type bpsd_plasmaf_data
         real(rkind) :: pn     ! Number density [m^-3]
         real(rkind) :: pt     ! Temperature [eV]
         real(rkind) :: ptpr   ! Parallel temperature [eV]
         real(rkind) :: ptpp   ! Perpendicular temperature [eV]
         real(rkind) :: pu     ! Parallel flow velocity [m/s]
      end type bpsd_plasmaf_data

      type bpsd_plasmaf_type
         integer :: nrmax     ! Number of radial points
         integer :: nsmax     ! Number of particle species
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
                              ! (rho^2) : normarized toroidal magnetic flux
         real(rkind), dimension(:), allocatable :: qinv 
                              ! 1/q : inverse of safety factor
         type(bpsd_plasmaf_data), dimension(:,:), allocatable :: data
      end type bpsd_plasmaf_type
c
      type bpsd_data0D_type
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: data
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_data0D_type
c
      type bpsd_data1D_type
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: s 
         real(rkind), dimension(:,:), allocatable :: data
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_data1D_type
c
      type bpsd_data2D_type
         integer :: nthmax    ! Number of poloidal points
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: th
         real(rkind), dimension(:), allocatable :: s 
         real(rkind), dimension(:,:,:), allocatable :: data
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_data2D_type
c
      type bpsd_data3D_type
         integer :: nphmax    ! Number of toroidal points
         integer :: nthmax    ! Number of poloidal points
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), allocatable :: ph
         real(rkind), dimension(:), allocatable :: th
         real(rkind), dimension(:), allocatable :: s 
         real(rkind), dimension(:,:,:,:), allocatable :: data
         character(len=32), dimension(:), allocatable :: kid
      end type bpsd_data3D_type
c
      end module bpsd_types
