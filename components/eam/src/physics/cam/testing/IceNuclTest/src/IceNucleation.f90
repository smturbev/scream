module micro_p3_minimal
 
   ! get real kind from utils
   use physics_utils, only: rtype,rtype8,btype

   ! physical and mathematical constants
   use micro_p3_utils, only: mi0,nccnst,T_zerodegc, T_icenuc

   implicit none

   public :: ice_nucleation

     !................................................................
      ! deposition/condensation-freezing nucleation
      call ice_nucleation(t_atm(k),inv_rho(k),&
           ni(k),ni_activated(k),qv_supersat_i(k),inv_dt,do_predict_nc, do_prescribed_CCN, &
           qinuc, ni_nucleat_tend)

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

   if ( t_atm .lt.T_icenuc .and. qv_supersat_i.ge.0.05_rtype) then
      if(.not. do_predict_nc .or. do_prescribed_CCN) then
!         ! dum = exp(-0.639+0.1296*100.*qv_supersat_i(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
         dum = 0.005_rtype*bfb_exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
         dum = min(dum,100.e3_rtype*inv_rho)
         N_nuc = max(0._rtype,(dum-ni)*inv_dt)
         if (N_nuc.ge.1.e-20_rtype) then
            Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt)
            qinuc = Q_nuc
            ni_nucleat_tend = N_nuc
         endif
      else
      ! Ice nucleation predicted by aerosol scheme
         ni_nucleat_tend = max(0._rtype, (ni_activated - ni)*inv_dt)
         qinuc = ni_nucleat_tend * mi0
      endif
   endif

end subroutine

end module micro_p3_minimal
