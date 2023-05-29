module IceNucleationOriginal

  ! get real kind from utils
  use physics_utils, only: rtype,rtype8,btype

  implicit none

  public :: ice_nucleation

  ! Constants from micro_p3_utils
  real(rtype), parameter :: mi0         = 3.77e-14 ! 4._rtype*piov3*900._rtype*1.e-18_rtype
  real(rtype), parameter :: T_zerodegc  = 273.15_rtype
  real(rtype), parameter :: T_icenuc    = 273.15_rtype - 15._rtype

contains

  subroutine ice_nucleation(t_atm,inv_rho,ni,ni_activated,qv_supersat_i,inv_dt,do_predict_nc, do_prescribed_CCN,   &
       qinuc,ni_nucleat_tend)

    !................................................................
    ! deposition/condensation-freezing nucleation
    ! allow ice nucleation if < -15 C and > 5% ice supersaturation
    ! use CELL-AVERAGE values, freezing of vapor

    implicit none

    real(rtype), intent(in) :: t_atm
    real(rtype), intent(in) :: inv_rho
    real(rtype), intent(in) :: ni
    real(rtype), intent(in) :: ni_activated
    real(rtype), intent(in) :: qv_supersat_i
    real(rtype), intent(in) :: inv_dt
    logical(btype), intent(in) :: do_predict_nc, do_prescribed_CCN

    real(rtype), intent(inout) :: qinuc
    real(rtype), intent(inout) :: ni_nucleat_tend


    real(rtype) :: dum, N_nuc, Q_nuc

    ! initialize nucleation tendency
    ni_nucleat_tend = 0.

    if ( t_atm .lt.T_icenuc .and. qv_supersat_i.ge.0.05_rtype) then
      if(.not. do_predict_nc .or. do_prescribed_CCN) then
        !         ! dum = exp(-0.639+0.1296*100.*qv_supersat_i(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
!!        dum = 0.005_rtype*bfb_exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
        dum = 0.005_rtype*exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
        dum = min(dum,100.e3_rtype*inv_rho)
        N_nuc = max(0._rtype,(dum-ni)*inv_dt)
        if (N_nuc.ge.1.e-20_rtype) then
          ni_nucleat_tend = N_nuc
        endif
      else
        ! Ice nucleation predicted by aerosol scheme
        ni_nucleat_tend = max(0._rtype, (ni_activated - ni)*inv_dt)
      endif
    endif

    ! Assume all ice nuclei have mass mi0.
    qinuc = ni_nucleat_tend * mi0

  end subroutine ice_nucleation

end module IceNucleationOriginal
