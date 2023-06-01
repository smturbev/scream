module physics_utils
  use shr_kind_mod,   only: rtype=>shr_kind_r8, itype=>shr_kind_i8
  integer,parameter,public :: btype = kind(.true.) ! native logical
  integer,parameter,public :: rtype8 = selected_real_kind(15, 307) ! 8 byte real, compatible with c type double
  real,parameter,public ::  NumCirrusINP = 2.e-3 !in units per cm3! BG added some reasonable background upper tropospheric dust/het ice nuclei value for cirrus conditions
  !BG: based on simulated upper tropospheric (200 hPa) values of dust in Tropical Western Pacific for ECHAM-HAM GCM simulations
  real,public :: NumCirrusSulf = 20. !20 will effectively limit the ice nucleation for one event to 20 k /L: better limit, as no way to deplete air mass of nsulf in this code
  !100. !100/cm3 
  !BG max number of sulphate aerosol used for homog freezing -> taken some reasonable value for upper troposphere
  !used this number as Liu and Penner, 2005
end module physics_utils
