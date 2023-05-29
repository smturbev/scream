module physics_utils
  use shr_kind_mod,   only: rtype=>shr_kind_r8, itype=>shr_kind_i8
  integer,parameter,public :: btype = kind(.true.) ! native logical
  integer,parameter,public :: rtype8 = selected_real_kind(15, 307) ! 8 byte real, compatible with c type double

end module physics_utils
