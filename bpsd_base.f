c     $Id$
c=======================================================================
      module bpsd
c
      use bpsd_kind
      use bpsd_constant
      use bpsd_flags
      use bpsd_types
      use bpsd_shot
      use bpsd_device
      use bpsd_species
      use bpsd_equ1d
      use bpsd_metric1d
      use bpsd_plasmaf
c
      logical, private, save :: bpsd_init_flag = .TRUE.
c      type(bpsd_shotx_type), private, save :: shot
c      type(bpsd_0ddatax_type), private, save :: devicex
c      type(bpsd_0ddatax_type), private, save :: speciesx
c      type(bpsd_1ddatax_type), private, save :: equ1Dx
c      type(bpsd_1ddatax_type), private, save :: metric1Dx
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
      end module bpsd
