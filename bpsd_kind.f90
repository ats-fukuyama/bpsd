module bpsd_kind
  implicit none
  integer, parameter :: rkind=selected_real_kind(12,100)
  integer, parameter :: ikind=selected_int_kind(8)
  integer, parameter :: dp=rkind
end module bpsd_kind
