! bpsd_kinds.f90

module bpsd_kinds
  implicit none
  integer, parameter :: rkind=selected_real_kind(12,100)
  integer, parameter :: ikind=selected_int_kind(8)
  integer, parameter :: dp=rkind
end module bpsd_kinds
