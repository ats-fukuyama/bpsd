module bpsd_constants
  use bpsd_kinds
  implicit none

  real(rkind),parameter :: ZERO = 0.0_dp
  real(rkind),parameter :: HALF = 0.5_dp
  real(rkind),parameter :: ONE  = 1.0_dp
  real(rkind),parameter :: TWO  = 2.0_dp

  real(rkind),parameter :: PI   = 3.14159265358979323846_dp
  real(rkind),parameter :: TWOPI= PI+PI

  ! ***** Physical constants based on CODATA 2006 *****

  real(rkind),parameter :: AEE  = 1.602176487E-19_dp ! elementary charge
  real(rkind),parameter :: AME  = 9.10938215E-31_dp  ! electron mass
  real(rkind),parameter :: AMP  = 1.672621637E-27_dp ! proton mass
  real(rkind),parameter :: VC   = 2.99792458E8_dp  ! speed of light
  real(rkind),parameter :: RMU0 = 4.E-7_dp*PI      ! permeability
  real(rkind),parameter :: EPS0 = ONE/(VC*VC*RMU0) ! permittivity

  complex(rkind),parameter :: CI   = (0.0_dp,1.0_dp)
end module bpsd_constants
