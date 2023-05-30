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

subroutine ice_nucleation(t_atm, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt, &
   qc, & ! added for new ice freezing from BG
   do_predict_nc, do_prescribed_CCN,   & ! old 
   do_default, do_bg_freezing, no_cirrus_mohler_ice_nucleation, & ! added for new ice freezing from BG
   qinuc, ni_nucleat_tend)

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
   real(rtype), intent(in) :: qc
   logical(btype), intent(in) :: do_predict_nc, do_prescribed_CCN, do_default, do_bg_freezing, no_cirrus_mohler_ice_nucleation

   real(rtype), intent(inout) :: qinuc
   real(rtype), intent(inout) :: ni_nucleat_tend


   real(rtype) :: dum, ndust, N_nuc, Q_nuc ! ST added ndust from BG new code
   
   ! Meyers or Cooper? 
   ! dum = exp(-0.639+0.1296*100.*qv_supersat_i(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
   dum = 0.005_rtype*bfb_exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
   
   ! choose default of bg_freezing 
   if ( do_default ) then
      ! default freezing from SCREAM
      if ( t_atm .lt.T_icenuc .and. qv_supersat_i.ge.0.05_rtype) then
         if(.not. do_predict_nc .or. do_prescribed_CCN) then
            dum = 0.005_rtype*bfb_exp(0.304_rtype*(T_zerodegc-t_atm))*1000.
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
   elseif ( do_bg_freezing ) then
      ! TODO: inseart ice freezing parameterization from Blaz Gasparini (BG)
      ! !................................................................
      ! 1.) deposition/condensation-freezing nucleation for MIXED PHASE
      ! allow ice nucleation if < -15 C and > 5% ice supersaturation
      !BG added a min temp condition (-38 deg C, hom frz threshold)
      !BG added a qc>qsmall condition: based on recent findings, deposition freezing is negligible in mixed phase (e.g. Ansmann et al., 2018)
      !BG therefore in mixed phase we can freeze only by immersion and contact and maybe also condensation freezing 
      !BG (e.g. Hoose and Mohler, 2012, Lohmann et al., 2016)
      !if ( (t(i,k).lt.258.15) .and. (t(i,k) .gt. 236.15) .and. (supi(i,k).ge.0.05) .and. (qc(i,k) .gt. qsmall) ) thenn
      !BG added a min temp condition (-38 deg C, hom frz threshold)
      !BG added a qc>qsmall condition: based on recent findings, deposition freezing is negligible in mixed phase (e.g. Ansmann et al., 2018)
      if ( ( (t_atm.lt.258.15) .and. (t_atm .ge. 236.15) .and. (qv_supersat_i.ge.0.05) .and. (qc .gt. qsmall) ) .or. ( (t_atm.lt.241.15) .and. (t_atm .ge. 236.15) .and. (qv_supersat_i.ge.0.005) .and. (qc .gt. qsmall) ) ) then  !BG added this to prevent too much "hom frz at -37"
         ! freezing in mixed phase
         dum = min(dum,150.e3*inv_rho(i,k)) !BG increased max limit from 100 to 150/L
         N_nuc = max(0._rtype,(dum-ni)*inv_dt)
         if (N_nuc.ge.1.e-20) then
            Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt)
            qinuc = Q_nuc
            ni_nucleat_tend = N_nuc
         endif
      endif
      if (no_cirrus_mohler_ice_nucleation .eq. .false.) then
         ! 2.) deposition/condensation-freezing nucleation for CIRRUS
         ! following Mohler et al., 2006 lab results for dust deposition freezing
         if (t_atm .lt. 220.15) then 
            scrit=0.1 !critical supersaturation for T<-53degC
         else
            scrit=0.2
         endif
         
         if ( (t_atm .lt. 236.15) .and. (qv_supersat_i .ge. scrit) ) then
            ! T < -27degC and supersat > critical value
            ndust = NumCirrusINP ! get value from utils
            dum = ndust *1e6*inv_rho !from /cm3 to kg-1 !assume some small INP/dust concentration, say 2/L which all freeze by deposition freezing
            dum = min(dum,100.e3*inv_rho) !max to 100/liter
            N_nuc =max(0.,(dum-ni)*inv_dt)
            if (N_nuc.ge.1.e-20) then
               Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*odt) 
               qinuc = Q_nuc
               ni_nucleat_tend = N_nuc
            endif
         endif 
      endif ! no_cirrus_mohler_ice_nucleation .eq. .false. 
      if (no_lphom_ice_nucleation .eq. .false.) then 
         ! 3.) HOM nucleation by using Liu and Penner, 2005 => ONLY HOMOG NUCLEATION!!!
         if ( (t_atm.lt.236.15)  .and. (qv_supersat_i .ge. 0.42)) then !added some very conservative supi condition not to always go in that loop ! ST - is supersat_cld the same as qv_supersat_i? 
            call hf(t_atm, uzpl, relhum(i,k), nsulf, nihf)
         endif ! 
      endif ! no_lphom_ice_nucleation .eq. .false.
   endif ! do_default elseif do_bg_freezing

end subroutine

!BG added homo freezing by Liu and Penner, 2005
!===============================================================================
!it assumes a fast and slow growth regime based on RHw and temperature.
!it also needs sulphate aerosol number, which we prescribe to 30/cc.

subroutine hf(T,ww,RH,Na,Ni) !tc, uzpl,relhum, nsulf,nihom
       real, intent(in)  :: T, ww, RH, Na
       real, intent(out)  :: Ni

       real ::    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
       real ::    A2_fast,B2_fast
       real ::    C1_fast,C2_fast,k1_fast,k2_fast
       real ::    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
       real ::    C1_slow,C2_slow,k1_slow,k2_slow
       real ::    regm
       real ::    A,B,C
       real ::    RHw
  !print*, 'we are in homogfrz'
  !print*, T, 'temp in homogfrz'
  !print*, ww,'vertv in homogfrz'
  !print*, RH, 'relh in homogfrz'
  !print*, Na, 'sulf in homogfrz'
 !---------------------------------------------------------------------
 ! parameters

       A1_fast  =0.0231
       A21_fast =-1.6387  !(T>-64 deg)
       A22_fast =-6.045   !(T<=-64 deg)
       B1_fast  =-0.008
       B21_fast =-0.042   !(T>-64 deg)
       B22_fast =-0.112   !(T<=-64 deg)
       C1_fast  =0.0739
       C2_fast  =1.2372

       A1_slow  =-0.3949
       A2_slow  =1.282
       B1_slow  =-0.0156
       B2_slow  =0.0111
       B3_slow  =0.0217
       C1_slow  =0.120
       C2_slow  =2.312

       Ni = 0.0

 !----------------------------
 !RHw parameters
       A = 6.0e-4*log(ww)+6.6e-3
       B = 6.0e-2*log(ww)+1.052
       C = 1.68  *log(ww)+129.35
       RHw=(A*T*T+B*T+C)!*0.01 in %

       if((T.le. -37.0) .and. ((RH).ge.RHw)) then
         !print*,'RHw crit in hf', RHw !BG
         !print*,'RH in hf', RH  !BG
         !print*, 'wind in hf' ,ww !BG
         regm = 6.07*log(ww)-55.0

         if(T.ge.regm) then    ! fast-growth regime

           if(T.gt.-64.0) then
             A2_fast=A21_fast
             B2_fast=B21_fast
           else
             A2_fast=A22_fast
             B2_fast=B22_fast
           endif

           k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
           k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

           Ni = k1_fast*Na**(k2_fast)
           Ni = min(Ni,Na)

         else       ! slow-growth regime

           k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
           k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

           Ni = k1_slow*Na**(k2_slow)
           Ni = min(Ni,Na)

         endif

       end if

end subroutine hf

end module micro_p3_minimal
