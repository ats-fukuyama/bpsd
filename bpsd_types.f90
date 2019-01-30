!     $Id$
!=======================================================================
      module bpsd_types

      use bpsd_kinds

      type bpsd_shot_type
         integer :: shotID
         integer :: modelID
         character(len=32) :: deviceID
      end type bpsd_shot_type

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

      type bpsd_species_data
         real(rkind) :: pa     ! Mass number (n. of protons + n. of neutrons)
         real(rkind) :: pz     ! Charge number (n. of protons - n. of electrons)
         integer ::     npa    ! Atomic number (n. of protons)
      end type bpsd_species_data

      type bpsd_species_type
         integer :: nsmax     ! Number of particle species
         integer :: idum      ! Dummy
         type(bpsd_species_data), dimension(:), pointer :: data
      end type bpsd_species_type

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
         real(rkind), dimension(:), pointer :: rho
                              ! normalized minor radius
         type(bpsd_equ1D_data), dimension(:), pointer :: data
      end type bpsd_equ1D_type

      type bpsd_equ2D_data
         real(rkind) :: psip2d  ! Poloidal magnetic flux [Wb] ~pi*R*r*Bp
      end type bpsd_equ2D_data
      type bpsd_equ2D_type
         integer :: nrrmax      ! Number of major radius points
         integer :: nzzmax      ! Number of vertical points
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: rr
                              ! Major radius
         real(rkind), dimension(:), pointer :: zz
                              ! Vertical position
         real(rkind) :: psip_axis,psip_boundary
         type(bpsd_equ2D_data), dimension(:,:), pointer :: data
      end type bpsd_equ2D_type

      type bpsd_equ3D_data
         real(rkind) :: psip3d  ! Poloidal magnetic flux [Wb] ~pi*R*r*Bp
      end type bpsd_equ3D_data
      type bpsd_equ3D_type
         integer :: nrmax       ! Number of radial radius points
         integer :: nthmax      ! Number of poloidal points
         integer :: nphmax      ! Number of toroidal points
         integer :: idum      ! Dummy
         real(rkind) :: time
         real(rkind) :: psip_axis,psip_boundary
         real(rkind), dimension(:), pointer :: rho
                              ! normalized minor radius
         real(rkind), dimension(:), pointer :: theta
                              ! poloidal angle
         real(rkind), dimension(:), pointer :: phi
                              ! toroidal angle
         type(bpsd_equ3D_data), dimension(:,:,:), pointer :: data
      end type bpsd_equ3D_type

      type bpsd_metric1D_data
         real(rkind) :: pvol     ! Plasma volume [m^3] ~2*pi*R*pi*r^2
         real(rkind) :: psur     ! Plasma surface [m^2] ~4*pi*R*pi*r
         real(rkind) :: dvpsit   ! dV/dPsit
         real(rkind) :: dvpsip   ! dV/dPsip
         real(rkind) :: aver2    ! <R^2>
         real(rkind) :: aver2i   ! <1/R^2>
         real(rkind) :: aveb2    ! <B^2>
         real(rkind) :: aveb2i   ! <1/B^2>
         real(rkind) :: avegv    ! <|gradV|>
         real(rkind) :: avegv2   ! <|gradV|^2>
         real(rkind) :: avegvr2  ! <|gradV|^2/R^2>
         real(rkind) :: avegr    ! <|grad rho|>
         real(rkind) :: avegr2   ! <|grad rho|^2>
         real(rkind) :: avegrr2  ! <|grad rho|^2/R^2>
         real(rkind) :: avegpp2  ! <|gradPsip|^2>
         real(rkind) :: rr       ! R
         real(rkind) :: rs       ! r
         real(rkind) :: elip     ! elipticity
         real(rkind) :: trig     ! triangularity
         real(rkind) :: aveb     ! <B>
      end type bpsd_metric1D_data

      type bpsd_metric1D_type
         integer :: nrmax       ! Number of radial points
         integer :: idum      ! Dummy
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: rho 
                                ! normailized minor radius
         type(bpsd_metric1D_data), dimension(:), pointer :: data
      end type bpsd_metric1D_type

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
         real(rkind), dimension(:), pointer :: rho
                              ! normalized minor radius
         real(rkind), dimension(:), pointer :: qinv 
                              ! 1/q : inverse of safety factor
         type(bpsd_plasmaf_data), dimension(:,:), pointer :: data
      end type bpsd_plasmaf_type

      type bpsd_dielectric_data
         complex(rkind),dimension(3,3) :: tensor     ! permeability
      end type bpsd_dielectric_data
      type bpsd_dielectric_type
         integer :: nrmax     ! Number of radial points
         integer :: nsmax     ! Number of particle species
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: rho
                              ! normalized minor radius
         type(bpsd_dielectric_data), dimension(:,:), pointer :: data
      end type bpsd_dielectric_type

      type bpsd_data0D_type
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: data
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data0D_type

      type bpsd_data1D_type
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: s 
         real(rkind), dimension(:,:), pointer :: data
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data1D_type

      type bpsd_data2D_type
         integer :: nthmax    ! Number of poloidal points
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         integer :: idum      ! Dummy
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: th
         real(rkind), dimension(:), pointer :: s 
         real(rkind), dimension(:,:,:), pointer :: data
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data2D_type

      type bpsd_data3D_type
         integer :: nphmax    ! Number of toroidal points
         integer :: nthmax    ! Number of poloidal points
         integer :: nrmax     ! Number of radial points
         integer :: ndmax     ! Number of data
         character(len=32) :: dataName
         real(rkind) :: time
         real(rkind), dimension(:), pointer :: ph
         real(rkind), dimension(:), pointer :: th
         real(rkind), dimension(:), pointer :: s 
         real(rkind), dimension(:,:,:,:), pointer :: data
         character(len=32), dimension(:), pointer :: kid
      end type bpsd_data3D_type

      end module bpsd_types
