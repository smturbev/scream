!===============================================================================
! CLDERA passive tracers
! provides dissipation rates and sources for diagnostic constituents
!
! AOA tracer implemented as defined in:
! Gupta, A., Gerber, E.P. and Lauritzen, P.H. (2020) 
! Numerical impacts on tracer transport: A proposed intercomparison test of Atmospheric 
! General Circulation Models, Quarterly Journal of the Royal Meteorological Society,
! doi:10.1002/qj.3881.
!
! E90 tracer implemented as defined in:
! Abalos, M. et al. (2017)
! Using the Artificial Tracer e90 to Examine Present and Future UTLS Tracer Transport in WACCM, 
! Journal of the Atmospheric Sciences,
! doi:10.1175/JAS-D-17-0135.1.
! 
! ST80_25 tracer implemented as defined in
! Eyring, V. et al. (2013)
! Overview of IGAC/SPARC Chemistry-Climate Model Initiative (CCMI) Community Simulations in 
! Support of Upcoming Ozone and Climate Assessments
!===============================================================================

module cldera_passive_tracers

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver
  use constituents,   only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,    only: iulog
  use ref_pres,       only: pref_mid_norm
  use cam_abortutils, only: endrun
  use micro_p3_interface, only: ixcldliq, ixcldice
  use cldera_passive_tracers_indices, only: ifirst, ixaoa, ixbcu, ixnuc, ixnucni, ixnucw

  implicit none
  private
  save

  ! Public interfaces
  public :: cldera_passive_tracers_register         ! register constituents
  public :: cldera_passive_tracers_implements_cnst  ! true if constituent is implemented by this package
  public :: cldera_passive_tracers_init_cnst        ! initialize constituent field
  public :: cldera_passive_tracers_init             ! initialize history fields, datasets
  public :: cldera_passive_tracers_timestep_tend    ! calculate tendencies
  public :: cldera_passive_tracers_readnl           ! read namelist options

  ! ----- Private module data
  integer, parameter :: ncnst=5  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'AOA     ','BCU     ', 'NUC     ','NI_NUC  ', 'W_NUC  '/)

  ! integer :: ifirst ! global index of first constituent
  ! integer :: ixaoa  ! global index for AOA tracer
  ! integer :: ixbcu  ! global index for BCU tracer
  ! integer, public :: ixnuc = -1  ! global index for NUC tracer
  ! integer :: ixe90  ! global index for E90 tracer
  ! integer :: ixst80 ! global index for ST80_25 tracer

  ! Data from namelist variables
  logical :: cldera_passive_tracers_flag  = .true.    ! true => turn on test tracer code, namelist variable
  logical :: cldera_passive_read_from_ic_file = .false. ! true => tracers initialized from IC file

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine cldera_passive_tracers_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cldera_passive_tracers_readnl'


    namelist /cldera_passive_tracers_nl/ cldera_passive_tracers_flag, cldera_passive_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'cldera_passive_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, cldera_passive_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(cldera_passive_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(cldera_passive_read_from_ic_file, 1, mpilog,  0, mpicom)
#endif

  endsubroutine cldera_passive_tracers_readnl

!================================================================================

  subroutine cldera_passive_tracers_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents
    ! 
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    !  initialize ix to be -1
    ixnuc = -1
    ixbcu = -1
    ixaoa = -1
    ixnucni = -1
    ixnucw = -1

    if (.not. cldera_passive_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixaoa,  readiv=cldera_passive_read_from_ic_file, &
                  longname='Age-of-air tracer')
    ifirst = ixaoa
    call cnst_add(c_names(2), 28._r8/1000._r8, cpair, 0._r8, ixbcu,  &
                  readiv=cldera_passive_read_from_ic_file, longname='Buoyant convective updraft tracer', mixtype='dry')
    call cnst_add(c_names(3), 28._r8/1000._r8, cpair, 0._r8, ixnuc, &
                  readiv=cldera_passive_read_from_ic_file, longname='Nucleation tracer', mixtype='dry')
    call cnst_add(c_names(4), 28._r8/1000._r8, cpair, 0._r8, ixnucni,  &
                  readiv=cldera_passive_read_from_ic_file, longname='Num ice at nucleation', mixtype='dry')
    call cnst_add(c_names(5), 28._r8/1000._r8, cpair, 0._r8, ixnucw, &
                  readiv=cldera_passive_read_from_ic_file, longname='W at nucleation', mixtype='dry')
    
  end subroutine cldera_passive_tracers_register

!===============================================================================

  function cldera_passive_tracers_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: cldera_passive_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    cldera_passive_tracers_implements_cnst = .false.

    if (.not. cldera_passive_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          cldera_passive_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function cldera_passive_tracers_implements_cnst

!===============================================================================

  subroutine cldera_passive_tracers_init_cnst(name, q, gcid)

    !----------------------------------------------------------------------- 
    !
    ! Purpose: initialize test tracers mixing ratio fields 
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, q, gcid)
       endif
    end do

  end subroutine cldera_passive_tracers_init_cnst

!===============================================================================

  subroutine cldera_passive_tracers_init

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, horiz_only

    integer :: m, mm
    !-----------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) return

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld (cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
       call add_default (cnst_name(mm), 1, ' ')
    end do
    
  end subroutine cldera_passive_tracers_init

!===============================================================================

  subroutine cldera_passive_tracers_timestep_tend(state, ptend, dt, cflx)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_rlat_all_p , get_lat_all_p
    use physconst,     only: gravit, avogad
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep
    use ref_pres,      only: pref_mid_norm
    use time_manager,  only: get_curr_time
    use tropopause,    only: tropopause_find, TROP_ALG_TWMO, TROP_ALG_STOBIE, NOTFOUND

    ! Arguments
    type(physics_state), intent(inout) :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep
    real(r8),            intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk             ! chunk identifier
    integer :: ncol              ! no. of column in chunk
    integer :: nstep             ! current timestep number
    integer :: trop_level(pcols) ! tropopause level for all columns 
    ! integer :: ixcldliq, ixcldice ! from micro_p3_interface
    logical  :: lq(pcnst)

    integer  :: day,sec          ! date variables
    real(r8) :: t                ! tracer boundary condition
    real(r8) :: aoa_scaling      ! scale AOA1 from nstep to time
    real(r8) :: bcu_scaling      ! scale AOA1 from nstep to time
    real(r8) :: nuc_scaling      ! scale AOA1 from nstep to time

    !------------------------------------------------------------------

    if (.not. cldera_passive_tracers_flag) then
       !Initialize an empty ptend for use with physics_update
       call physics_ptend_init(ptend,state%psetcols,'cldera_passive_trc_ts')
       return
    end if

    lq(:)      = .FALSE.
    lq(ixaoa)  = .TRUE.
    lq(ixbcu)  = .TRUE.
    lq(ixnuc)  = .TRUE.
    lq(ixnucni)  = .TRUE.
    lq(ixnucw) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'cldera_passive_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk
    ncol  = state%ncol

    ! ---- compute nuc and bcu time scaling (1 s in hours)
    nuc_scaling = 1._r8/3600._r8
    bcu_scaling = 1._r8/3600._r8

    ! ---- compute AOA time scaling (1 s in days)
    aoa_scaling = 1._r8/86400._r8

    ! -------------------- TRACER TENDENCIES --------------------
    do k = 1, pver
       do i = 1, ncol

          ! ============ AOA ============
          ! clock tracer with a source of 1 day/day everywhere above ~700hPa
          if (pref_mid_norm(k) <= 0.7_r8) then
              ptend%q(i,k,ixaoa) = 1.0_r8 * aoa_scaling
          else
              ptend%q(i,k,ixaoa) = 0.0_r8
              state%q(i,k,ixaoa) = 0.0_r8
          end if

          ! ============ BCU ============
          ! clock tracer with a source of 1 hour everywhere in 
          ! a cloudy, rising parcel; set ptend
          ! else decay with timescale ~ 1 hour (3600 s)
          if ( (state%omega(i,k) <= -0.1_r8) .and. ((state%q(i,k,ixcldliq)+state%q(i,k,ixcldice)) > 1.e-5_r8 ) ) then
              ptend%q(i,k,ixbcu) = (1.0_r8 - state%q(i,k,ixbcu))/ dt
          else 
              ptend%q(i,k,ixbcu) = -state%q(i,k,ixbcu) * bcu_scaling
          end if

          ! ============ NUC ============
          ! Decay everywhere here but reset tend in mphys (P3 interface)
          ! use this module and indices ixnuc then reset anytime we have fresh nucleation in P3
          ptend%q(i,k,ixnuc) = -state%q(i,k,ixnuc) * nuc_scaling
          ! Decay everywhere for ixnucni and ixnucw at the same rate as ixnuc
          ptend%q(i,k,ixnucni) = -state%q(i,k,ixnucni) * nuc_scaling
          ptend%q(i,k,ixnucw)  = -state%q(i,k,ixnucw)  * nuc_scaling

    ! --------------- TRACER FLUXES --------------------- 
    do i = 1, ncol

       ! ====== AOA ======
       ! no surface flux
       cflx(i,ixaoa) = 0._r8

       ! ====== BCU ======
       ! no surface flux
       cflx(i,ixbcu) = 0._r8
       
       ! ====== NUC ======
       ! no surface flux
       cflx(i,ixnuc) = 0._r8
       cflx(i,ixnucni) = 0._r8
       cflx(i,ixnucw) = 0._r8
    end do
           

  end subroutine cldera_passive_tracers_timestep_tend

!===========================================================================

  subroutine init_cnst_3d(m, q, gcid)

    use dyn_grid,    only : get_horiz_grid_d, get_horiz_grid_dim_d
    use dycore,      only : dycore_is
    use ref_pres,    only : pref_mid_norm

    integer,  intent(in)  :: m       ! global constituent index
    real(r8), intent(out) :: q(:,:)  ! kg tracer/kg dry air (gcol,plev)
    integer,  intent(in)  :: gcid(:) ! global column id

    real(r8), allocatable :: lat(:)
    integer :: plon, plat, ngcols
    integer :: j, k, gsize
    !-----------------------------------------------------------------------
    ! initialization below is on the DYNAMICS grid; length of gcol will 
    ! be number of GLL points, not number of physics columns

    if (masterproc) write(iulog,*) 'CLDERA PASSIVE CONSTITUENTS: INITIALIZING ',cnst_name(m),m
    
    ! ====== AOA ======
    if (m == ixaoa) then
       q(:,:) = 0.0_r8
    
    ! ====== BCU ======
    else if (m == ixbcu) then
       q(:,:) = 0.0_r8
       
    ! ====== NUC ======
    else if (m == ixnuc) then
       q(:,:) = 0.0_r8
    else if (m == ixnucni) then
       q(:,:) = 0.0_r8
    else if (m == ixnucw) then
       q(:,:) = 0.0_r8
    end if

  end subroutine init_cnst_3d

!=====================================================================


end module cldera_passive_tracers
