! bpsd_kinds.f90

module bpsd_kinds
  USE ISO_FORTRAN_ENV
  integer, parameter :: rkind=REAL64
  integer, parameter :: ikind=INT32
  integer, parameter :: dp=REAL64
  integer, parameter :: sp=REAL32
  integer, parameter :: long=INT64
!  integer, parameter :: rkind=selected_real_kind(12,100)
!  integer, parameter :: ikind=selected_int_kind(8)
!  integer, parameter :: dp=rkind
!  integer, parameter :: sp=selected_real_kind(6)
end module bpsd_kinds
