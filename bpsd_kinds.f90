module bpsd_kinds
  implicit none
  integer, parameter :: rkind=selected_real_kind(12,100)
  integer, parameter :: ikind=selected_int_kind(8)
end module bpsd_kinds
