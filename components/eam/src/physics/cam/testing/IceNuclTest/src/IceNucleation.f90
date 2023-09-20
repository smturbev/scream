module micro_p3_minimal
 
   ! get real kind from utils
   use physics_utils, only: rtype,rtype8,btype

   ! physical and mathematical constants
!    use micro_p3_utils, only! : mi0, T_zerodegc, T_icenuc, NumCirrusINP, NumCirrusSulf

   implicit none

   public :: ice_nucleation
   
   ! Constants from micro_p3_utils
     real(rtype), parameter :: mi0         = 3.77e-14_rtype ! 4._rtype*piov3*900._rtype*1.e-18_rtype
     real(rtype), parameter :: T_zerodegc  = 273.15_rtype
     real(rtype), parameter :: T_icenuc    = 273.15_rtype - 15._rtype
     real(rtype), parameter :: NumCirrusINP= 2.0e-3_rtype
     real(rtype), parameter :: NumCirrusSulf= 20._rtype

     ! !................................................................
!       ! deposition/condensation-freezing nucleation
!       call ice_nucleation(t_atm(k),inv_rho(k),&
!            ni(k),ni_activated(k),qv_supersat_(k),qv_supersat_i(k),inv_dt, qc(k), uzpl(k), &
!            do_predict_nc, do_prescribed_CCN, &
!            do_new_lp_freezing, no_cirrus_mohler_ice_nucleation, no_lphom_ice_nucleation, &
!            qinuc, ni_nucleat_tend)

contains

subroutine ice_nucleation(t_atm, inv_rho, ni, ni_activated, qv_supersat_l, qv_supersat_i, inv_dt, &
   qc, uzpl, & ! added for new ice freezing from BG
   do_predict_nc, do_prescribed_CCN,   & ! old 
   do_new_lp_freezing, no_cirrus_mohler_ice_nucleation, no_lphom_ice_nucleation,  & ! added for new ice freezing from BG ! TODO: add use_preexisting_ice
   qinuc, ni_nucleat_tend,nnuc0,nnuc1,nnuc2,nnuc3,nnuc4)

   !................................................................
   ! deposition/condensation-freezing nucleation
   ! allow ice nucleation if < -15 C and > 5% ice supersaturation
   ! use CELL-AVERAGE values, freezing of vapor
   ! TODO: ST - make different ice categories for each type of nucleation,
   !       then add together at the end (or make separate output for each)

   implicit none

   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: inv_rho
   real(rtype), intent(in) :: ni
   real(rtype), intent(in) :: ni_activated
   real(rtype), intent(in) :: qv_supersat_i
   real(rtype), intent(in) :: qv_supersat_l
   real(rtype), intent(in) :: inv_dt
   real(rtype), intent(in) :: qc   ! cloud water mixing ratio (kg/kg)
   real(rtype), intent(in) :: uzpl ! vertical velocity (Pa/s)
   logical(btype), intent(in) :: do_predict_nc, do_prescribed_CCN, no_cirrus_mohler_ice_nucleation, no_lphom_ice_nucleation, do_new_lp_freezing

   real(rtype), intent(inout) :: qinuc
   real(rtype), intent(inout) :: ni_nucleat_tend
   real(rtype), intent(out) :: nnuc0, nnuc1, nnuc2, nnuc3, nnuc4 ! separate into categories based on pathway of nucleation

   ! local variables
   real(rtype) :: dum, N_nuc, Q_nuc
   real(rtype) :: ndust, nsulf, qsmall, w, niimm, nidep, nihf, scrit ! for new BG code
   real(rtype) :: wbar1, wbar2, deles, esi, A, B, regm, n1, tc ! work variables
   
   ! convert uzpl to w (Pa/s -> m/s)
   !w = - uzpl * inv_rho * 0.102 ! omega * 1/rho * 1/g  [m/s]
   w = uzpl ! for tests only
   
   ! minium allowable prognostic variables
    qsmall = 1.e-14
   
   ! get value from utils
   nsulf = NumCirrusSulf
   ndust = NumCirrusINP
   
   ! Meyers or Cooper? 
   ! dum = exp(-0.639+0.1296*100.*qv_supersat_i(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
   ! dum = 0.005_rtype**exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho   !Cooper (1986)
   
   ! choose default of bg_freezing 
   if ( do_new_lp_freezing .eq. .false. ) then
      ! default freezing from SCREAM
      if ( t_atm .lt.T_icenuc .and. qv_supersat_i.ge.0.05_rtype) then
         if(.not. do_predict_nc .or. do_prescribed_CCN) then
            dum = 0.005_rtype*exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho ! Cooper 1986
            dum = min(dum,100000._rtype*inv_rho)
            N_nuc = max(0._rtype,(dum-ni)*inv_dt)
            nnuc0=N_nuc
            if (N_nuc.ge.1.e-20_rtype) then
               Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt)
               qinuc = Q_nuc
               ni_nucleat_tend = N_nuc
            endif
         else
            ! Ice nucleation predicted by aerosol scheme
            ni_nucleat_tend = max(0._rtype, (ni_activated - ni)*inv_dt)
            nnuc0=ni_nucleat_tend
            qinuc = ni_nucleat_tend * mi0
         endif
      else
          nnuc0=0._rtype
      endif
   else ! do_new_lp_freezing .eq. .true.
      ! Adding four main components:
      ! 1.) deposition/condensation freezing nucleation in MIXED PHASE
      ! 2.) deposition/condensation-freezing nucleation for CIRRUS
      !       (only if no_cirrus_mohler_ice_nucleation.eq..false.)
      ! 3.) HOM nucleation by using Liu and Penner, 2005 => ONLY HOMOG NUCLEATION!!!
      ! 4.) competition with HOMOG + HETEROG NUCLEATION (+preex ice if set to true)
      !      a.) HET freezing only, Liu et al., 2007, eq 10. 
      !      b.) HOM nucleation only
      !      c.) transition between homogeneous and heterogeneous: interpolate in-between
      ! !................................................................
      ! 1.) deposition/condensation-freezing nucleation for MIXED PHASE
      ! allow ice nucleation if < -15 C and > 5% ice supersaturation
      !BG added a min temp condition (-38 deg C, hom frz threshold)
      !BG added a qc>qsmall condition: based on recent findings, deposition freezing 
      ! is negligible in mixed phase (e.g. Ansmann et al., 2018)
      !BG therefore in mixed phase we can freeze only by immersion and contact and maybe also condensation freezing 
      !BG (e.g. Hoose and Mohler, 2012, Lohmann et al., 2016)
      !if ( (t(i,k).lt.258.15) .and. (t(i,k) .gt. 236.15) .and. (supi(i,k).ge.0.05) .and. (qc(i,k) .gt. qsmall) ) thenn
      !BG added a min temp condition (-38 deg C, hom frz threshold)
      !BG added a qc>qsmall condition: based on recent findings, deposition freezing is negligible in mixed phase (e.g. Ansmann et al., 2018)
      tc= t_atm-T_zerodegc ! convert K to degC
      ! --------------------------------------------------------
      if ( ( (t_atm.lt.258.15_rtype) .and. (t_atm .ge. 236.15_rtype) &
      .and. (qv_supersat_i.ge.0.05_rtype) .and. (qc .gt. qsmall) ) .or. &
      ( (t_atm.lt.241.15_rtype) .and. (t_atm .ge. 236.15_rtype) .and. &
      (qv_supersat_i.ge.0.005_rtype) .and. (qc .gt. qsmall) ) ) then  
         !BG added this ^ to prevent too much "hom frz at -37" freezing in mixed phase
         ! 1.) deposition/condensation-freezing nucleation for MIXED PHASE
         ! allow ice nucleation if < -15 C and > 5% ice supersaturation
         !BG added a min temp condition (-38 deg C, hom frz threshold)
         !BG added a qc>qsmall condition: based on recent findings, deposition freezing 
         ! is negligible in mixed phase (e.g. Ansmann et al., 2018)
         !BG therefore in mixed phase we can freeze only by immersion and contact and maybe also condensation freezing 
         !BG (e.g. Hoose and Mohler, 2012, Lohmann et al., 2016)
         !if ( (t(i,k).lt.258.15) .and. (t(i,k) .gt. 236.15) .and. (supi(i,k).ge.0.05) .and. (qc(i,k) .gt. qsmall) ) then
         !BG added a min temp condition (-38 deg C, hom frz threshold)
         !BG added a qc>qsmall condition: based on recent findings, deposition freezing is negligible in mixed phase (e.g. Ansmann et al., 2018)
         dum = 0.005_rtype*exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho ! Cooper 1986
         dum = min(dum,150000._rtype*inv_rho) !BG increased max limit from 100 to 150/L
         nnuc1 = max(0._rtype,(dum-ni)*inv_dt)
         !print*,"in mixed phase - nnuc1",nnuc1, dum, ni
      else
         nnuc1=0._rtype
      endif
      ! -----------------------------------------------------------
      if (no_cirrus_mohler_ice_nucleation .eq. .false.) then
         ! 2.) deposition/condensation-freezing nucleation for CIRRUS
         ! following Mohler et al., 2006 lab results for dust deposition freezing
         if (t_atm .lt. 220.15_rtype) then 
            scrit=0.1_rtype !critical supersaturation for T<-53degC
         else
            scrit=0.2_rtype
         endif
         
         if ( (t_atm .lt. 236.15_rtype) .and. (qv_supersat_i .ge. scrit) ) then
            ! T < -37degC and supersat > critical value
            dum = ndust*1e6_rtype*inv_rho !from /cm3 to kg-1 !assume some small INP/dust concentration, say 2/L which all freeze by deposition freezing
            dum = min(dum,100.e3_rtype*inv_rho) !max to 100/liter
            nnuc2 =max(0._rtype,(dum-ni)*inv_dt)
            !print*,"in cirrus mohler - nnuc2",nnuc2
         else
            nnuc2=0._rtype
         endif 
      else
          nnuc2=0._rtype
      endif ! no_cirrus_mohler_ice_nucleation .eq. .false. 
      ! -------------------------------------------------------------
      if (no_lphom_ice_nucleation .eq. .false.) then 
         ! 3.) HOM nucleation by using Liu and Penner, 2005 => ONLY HOMOG NUCLEATION!!!
         if ( (t_atm.lt.236.15_rtype)  .and. (qv_supersat_i .ge. 0.42_rtype) ) then !added some very conservative supi condition not to always go in that loop ! 
             !print*,"in lphom",(t_atm.lt.236.15_rtype),(qv_supersat_i .ge. 0.42_rtype), qv_supersat_i
             call hf(tc, w, qv_supersat_l, nsulf, nihf)

             dum = nihf*1.e6_rtype*inv_rho !from cm-3 to m-3 to kg-1 air
             dum = min(dum,80.e6_rtype*inv_rho) !set max to 80000 per L or 80/cc
             nnuc3 =max(0._rtype,(dum-ni)*inv_dt)
             !print*,"in lp hom - nnuc3", nnuc3,dum,nihf
        else
             nnuc3=0._rtype
        endif ! ( (t_atm.lt.236.15)  .and. (qv_supersat_i .ge. 0.42) )  
      else
          nnuc3=0._rtype
      endif ! no_lphom_ice_nucleation .eq. .false.
      
      !---------------------------------------------------------------
      !ST copied from BG copied from e3sm code and nucleate_ice.F90
      ! 4.) competition with HOMOG + HETEROG NUCLEATION (+preex ice if set to true)
      ! temp variables that depend on use_preexisting_ice

      ! Method: The current method is based on Liu & Penner (2005)
      !  It related the ice nucleation with the aerosol number, temperature and the
      !  updraft velocity. It includes homogeneous freezing of sulfate, immersion
      !  freezing of soot, and Meyers et al. (1992) deposition nucleation
      !
      ! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010, and S. Turbeville 2023
      !----------------------------------------------------------------

      wbar1 = w
      wbar2 = w
      ! TODO: add pre-existing ice section from BG
      
      ! new lp freezing
      !temp in celsius tc = tair - 273.15_r8
      ! TODO: ST- check that all temperature has been converted back to Kelvin

      ! initialize
      niimm = 0._rtype
      nidep = 0._rtype
      nihf  = 0._rtype
      deles = 0._rtype
      esi   = 0._rtype
       
      if((t_atm.le.238.15_rtype) .and. (qv_supersat_i.ge.0.2_rtype) .and. (wbar1 .ge. 1.e-6_rtype)  ) then 
      ! 120% is here taken as heterogeneous freezing limit
      !? BG - use higher RHi threshold? 
      !BG added wbar1>1e-6 for num stability (log of small number->problems)
          A = -1.4938_rtype * log(ndust) + 12.884_rtype
          B = -10.41_rtype  * log(ndust) - 67.69_rtype

          regm = A * log(wbar1) + B ! regm from LP2005 eq. 4.5
          
          ! 4a.) heterogeneous nucleation only
          if (tc .gt. regm) then
              !	BG critical soot (in our case dust/some undefined HET ice nucleating particle) 
              ! ndust num conc above which only 
              ! HET freezing occurs, it depends on INP number and w
                if((tc.lt.-40_rtype) .and. (wbar1.gt.1._rtype)) then ! exclude T<-40 & W>1m/s from hetero. nucleation
                    call hf(tc,wbar1,qv_supersat_l,nsulf,nihf)
                    niimm=0._rtype
                    nidep=0._rtype
                    n1=nihf ! TODO: add preexisting ice condition
                    !print*,"in hom only"
                else ! heterogeneous freezing
                    call hetero(tc,wbar2,ndust,niimm,nidep)
                    ! TODO: add preexisting ice condition
                    nihf=0._rtype
                    n1=niimm+nidep ! nidep=0 by definition
                    !print*,"in hetero only"
                endif !(tc.lt.-40. .and. wbar1.gt.1.)
                
          ! 4b.) homogeneous nucleation only
          else if (tc.lt.regm-5._rtype) then
               call hf(tc,wbar1,qv_supersat_l,nsulf,nihf)
               niimm=0._rtype
               nidep=0._rtype
               n1=nihf ! TODO: add preexisting ice condition
          ! 4c.) transition between homogeneous and heterogeneous: interpolate in-between
          else
              if ((tc.lt.-40._rtype) .and. (wbar1.gt.1._rtype)) then ! exclude T<-40 & W>1m/s from hetero. nucleation
                  call hf(tc,wbar1,qv_supersat_l,nsulf,nihf)
                  !print*,"in hom only again"
                  niimm=0._rtype
                  nidep=0._rtype
                  n1=nihf ! TODO: add preexisting ice condition
              else  !this is now the "real" transitional regime with calls for both homogeneous and heterogeneous nucleation
                  call hf(regm-5._rtype,wbar1,qv_supersat_l,nsulf,nihf)
                  call hetero(regm,wbar2,ndust,niimm,nidep) !BG the way I programmed it, nidep is 0 by definition
                  if (nihf .le. (niimm+nidep)) then
                       n1 = nihf 
                  else
                       n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._rtype) ! competition LP2005 eq. 4.10
                       !print*,"true hom vs het comp", n1
                  endif
                  !print*,"in het vs hom comp"
              end if
          end if ! hetero vs homog frz
          nnuc4 = n1 !cm-3 
          !print*,"hom vs hetero frz - nnuc4",nnuc4
     else
          nnuc4=0._rtype
     end if  ! ((t_atm.le.238.15) .and. (qv_supersat_i.ge.0.2) .and. (wbar1 .ge. 1.e-6)  )

     dum = nnuc4*1.e+6_rtype*inv_rho ! change unit from #/cm3 to #/kg
     dum = min(dum,80.e+6_rtype*inv_rho) !set max to 80000 per L or 80/cc
     
 
     !N_nuc =max(0.,(dum-sum(nitot(i,k,:)))*odt)
     !BG I think there is no logic behind the upper statement in case we do "real" freezing. 
     !   It may be reasonable for Meyers/Cooper parameterization, but not here...in my understanding!!!
     !BG changed therefore to: 
     nnuc4 = max(0._rtype, dum-ni)*inv_dt
     N_nuc = nnuc1 + nnuc2 + nnuc3 + nnuc4
     if (N_nuc.ge.1.e-20_rtype) then
       Q_nuc = max(0._rtype,N_nuc*mi0)
       qinuc = Q_nuc
       ni_nucleat_tend = N_nuc
     !else qinuc and ni_nucleat_tend do not change from previous timestep
     endif ! N_nuc.ge.1.e-20
   end if ! do_new_lp_frz 
end subroutine

!BG added homo freezing by Liu and Penner, 2005
! Same as from nucleate_ice.F90
!===============================================================================
!it assumes a fast and slow growth regime based on RHw and temperature.
!it also needs sulphate aerosol number, which we prescribe to 30/cc.
! ST changed input to tk (temp in  Kelvin) and supersat wrt liquid
subroutine hf(T,ww,supersat,Na,Ni) !tc , w, qv_supersat_l, nsulf, nihom
       real(rtype), intent(in)   :: T, ww, supersat, Na
       real(rtype), intent(out)  :: Ni

       real(rtype) ::    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
       real(rtype) ::    A2_fast,B2_fast
       real(rtype) ::    C1_fast,C2_fast,k1_fast,k2_fast
       real(rtype) ::    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
       real(rtype) ::    C1_slow,C2_slow,k1_slow,k2_slow
       real(rtype) ::    regm
       real(rtype) ::    A,B,C
       real(rtype) ::    Sw
  !print*, 'we are in homogfrz'
  !print*, T, 'temp in homogfrz'
  !print*, ww,'vertv in homogfrz'
  !print*, RH, 'relh in homogfrz'
  !print*, Na, 'sulf in homogfrz'
 !---------------------------------------------------------------------
 ! parameters from LP2005 table 1

       A1_fast  =0.0231_rtype
       A21_fast =-1.6387_rtype  !(T>-64 degC)
       A22_fast =-6.045_rtype   !(T<=-64 degC)
       B1_fast  =-0.008_rtype
       B21_fast =-0.042_rtype   !(T>-64 degC)
       B22_fast =-0.112_rtype   !(T<=-64 degC)
       C1_fast  =0.0739_rtype
       C2_fast  =1.2372_rtype

       A1_slow  =-0.3949_rtype
       A2_slow  =1.282_rtype
       B1_slow  =-0.0156_rtype
       B2_slow  =0.0111_rtype
       B3_slow  =0.0217_rtype
       C1_slow  =0.120_rtype
       C2_slow  =2.312_rtype

       Ni = 0.0_rtype

 !----------------------------
 !RHw parameters
       A = 6.0e-4_rtype*log(ww)+6.6e-3_rtype
       B = 6.0e-2_rtype*log(ww)+1.052_rtype
       C = 1.68_rtype  *log(ww)+129.35_rtype
       Sw=((A*T*T+B*T+C)*0.01_rtype)-1._rtype ! LP2005 eq 3.1; in this case around 0.3

       if((T.le.-37.0_rtype) .and. ((supersat).ge.Sw)) then
         !print*,'RHw crit in hf', RHw !BG
         !print*,'RH in hf', RH  !BG
         !print*, 'wind in hf' ,ww !BG
         regm = 6.07_rtype*log(ww)-55.0_rtype ! from LP2005 eq. 3.5

         if (T.ge.regm) then ! fast-growth regime

           if(T.gt.-64.0_rtype) then
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
     !print*, 'in call hf: Sw > Sw_crit? Ni', (supersat.ge.Sw), Ni

end subroutine hf

!-----BG added for ice nucleation-----
subroutine hetero(tc,ww,Ns,Nis,Nid)

    real(rtype), intent(in)  :: tc, ww, Ns
    real(rtype), intent(out) :: Nis, Nid

    real(rtype) :: A11,A12,A21,A22,B11,B12,B21,B22
    real(rtype) :: B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_rtype
      A12 = -0.0185_rtype
      A21 = 2.758_rtype
      A22 = 1.3221_rtype
      B11 = -0.008_rtype
      B12 = -0.0468_rtype
      B21 = -0.2667_rtype
      B22 = -1.4588_rtype

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)
      !print*,"in call hetero"

      Nis = exp(A22) * Ns**B22 * exp(B*tc) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_rtype    ! don't include deposition nucleation for cirrus clouds when T<-37C
      ! BG we assume the current het freezing represents in some sense both het by imm and depo
end subroutine hetero

end module micro_p3_minimal