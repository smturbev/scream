module cldera_passive_tracers_indices

  implicit none
  save

  integer, public :: ifirst ! global index of first constituent
  integer, public :: ixaoa  ! global index for AOA tracer
  integer, public :: ixbcu  ! global index for BCU tracer
  integer, public :: ixnuc  ! global index for NUC tracer
  integer, public :: ixnucni ! global index for num ice when NUC tracer is active
  integer, public :: ixnucw  ! global index for vert. velocity when NUC tracer is active

end module cldera_passive_tracers_indices
