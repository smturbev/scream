module nucleate_ice_bg

!-------------------------------------------------------------------------------
! Purpose:
!  A parameterization of ice nucleation.
!
!  *** This module is intended to be a "portable" code layer.  Ideally it should
!  *** not contain any use association of modules that belong to the model framework.
!
!
! Method:
!  The current method is based on Liu & Penner (2005; hereafter LP2005) for 
!  homogeneous vs heterogeneous competition. It relates the ice 
!  nucleation with the aerosol number, temperature and the updraft velocity.
!  It includes homogeneous freezing of sulfate & immersion
!  freezing on dust in cirrus clouds and either Cooper (1986) or
!  Meyers et al. (1992) deposition nucleation in mixed-phase clouds.
!  
!  There are 'switches' for cirrus deposition based on Mohler et al. (2006) based
!  on their lab experiements or homogeneous freezing only based on LP2005. If 
!  we turn on homogeneous freezing only, then we turn off HET freezing.
!
!  The effect of preexisting ice crystals on ice nucleation in cirrus clouds is included, 
!  and also consider the sub-grid variability of temperature in cirrus clouds,
!  following X. Shi et al. ACP (2014).
!
!
! Authors:
!  Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010,
!  Xiangjun Shi & Xiaohong Liu, 01/2014,
!  S. Turbeville 2023
!
!  With help from C. C. Chen and B. Eaton (2014),
!     B. Gasparini and P. Blossey (2023)
!-------------------------------------------------------------------------------

use wv_saturation,  only: svp_water, svp_ice
use cam_logfile,    only: iulog

implicit none
private
save

integer, parameter :: r8 = selected_real_kind(12)

public :: nucleati_bg_init, nucleati_bg

logical  :: use_preexisting_ice
logical  :: do_meyers    
logical  :: do_new_bg_lp_frz   
logical  :: do_ci_mohler_dep  
logical  :: do_lphom     
logical  :: no_limits      

real(r8) :: pi
real(r8) :: mincld

! Subgrid scale factor on relative humidity (dimensionless)
real(r8) :: subgrid

real(r8), parameter :: Shet   = 0.2_r8     ! het freezing threshold for supersat of ice
real(r8), parameter :: rhoice = 0.5e3_r8   ! kg/m3, Wpice is not sensitive to rhoice
real(r8), parameter :: minweff= 0.001_r8   ! m/s
real(r8), parameter :: gamma1=1.0_r8 
real(r8), parameter :: gamma2=1.0_r8 
real(r8), parameter :: gamma3=2.0_r8 
real(r8), parameter :: gamma4=6.0_r8 

real(r8) :: ci

!===============================================================================
contains
!===============================================================================

subroutine nucleati_bg_init( &
   use_preexisting_ice_in, do_meyers_in, do_new_bg_lp_frz_in, &
   do_ci_mohler_dep_in, do_lphom_in, no_limits_in, &
   iulog_in, pi_in, mincld_in, subgrid_in)

   logical,  intent(in) :: use_preexisting_ice_in
   logical,  intent(in) :: do_meyers_in
   logical,  intent(in) :: do_new_bg_lp_frz_in
   logical,  intent(in) :: do_ci_mohler_dep_in
   logical,  intent(in) :: do_lphom_in
   logical,  intent(in) :: no_limits_in
   integer,  intent(in) :: iulog_in
   real(r8), intent(in) :: pi_in
   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: subgrid_in

   use_preexisting_ice = use_preexisting_ice_in
   do_meyers           = use_hetfrz_classnuc_in
   do_new_bg_lp_frz    = do_new_bg_lp_frz_in
   do_ci_mohler_dep    = do_ci_mohler_dep_in
   do_lphom            = do_lphom_in
   no_limits           = no_limits_in
   iulog               = iulog_in
   pi                  = pi_in
   mincld              = mincld_in
   subgrid             = subgrid_in

   ci = rhoice*pi/6._r8

end subroutine nucleati_init

!===============================================================================

subroutine nucleati_bg(  &
   wbar, tair, pmid, relhum, supersat_i &
   qc, qi, ni, inv_rho,                 &
   so4_num, dst_num,                    &
   nuci, onihf, oniimm, onidep, onimix, &
   wpice, weff, fhom )

   !---------------------------------------------------------------
   ! Purpose:
   !  The parameterization of ice nucleation.
   !
   ! Method: The current method is based on Liu & Penner (2005)
   !  It related the ice nucleation with the aerosol number, temperature and the
   !  updraft velocity. It includes homogeneous freezing of sulfate, immersion
   !  freezing of dust, and Meyers et al. (1992) deposition nucleation
   !
   ! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010,
   !  modifications by S. Turbeville 2023
   !----------------------------------------------------------------

   ! Input Arguments
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair        ! temperature (K)
   real(r8), intent(in) :: pmid        ! pressure at layer midpoints (pa)
   real(r8), intent(in) :: relhum      ! relative humidity with respect to liquid
   real(r8), intent(in) :: supersat_i  ! supersaturation with respect to ice
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(r8), intent(in) :: qi          ! ice mass mixing ratio (kg/kg)
   real(r8), intent(in) :: ni          ! grid-mean preexisting cloud ice number conc (#/kg) 
   real(r8), intent(in) :: inv_rho     ! inverse rho (1 / air density (kg/m3))
   real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   
   ! Output Arguments
   real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(r8), intent(out) :: onimix     ! nucleated number from deposition nucleation  
                                       !   (meyers: mixed phase)
   real(r8), intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by 
                                       !   preexisting ice (m/s), at Shom
   real(r8), intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s);
                                       !   weff=wbar-wpice
   real(r8), intent(out) :: fhom       ! how much fraction of cloud can reach Shom

   ! Local workspace
   real(r8) :: nihf                    ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                   ! nucleated number from immersion freezing
   real(r8) :: nidep                   ! nucleated number from deposition nucleation
   real(r8) :: nimix                   ! nucleated number from deposition nucleation (mixed phase)
   real(r8) :: nimoh                   ! nucleated number from cirrus deposition (Mohler 2006)
   real(r8) :: nilphf                  ! nucleated number from homogeneous freezing only (LP2005)
   real(r8) :: n1, ni                  ! nucleated number
   real(r8) :: qsmall                  ! min allowable cloud condensate to be considered cloud
   real(r8) :: tc, A, B, regm, dum     ! work variable
   real(r8) :: esl, esi, deles         ! work variable
   real(r8) :: wbar1, wbar2            ! work variable for preexisting ice

   ! used in SUBROUTINE Vpreice
   real(r8) :: Ni_preice        ! cloud ice number conc (1/m3)   
   real(r8) :: lami,Ri_preice   ! mean cloud ice radius (m)
   real(r8) :: Shom             ! initial ice saturation ratio; if <1, use hom threshold Si
   real(r8) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
   real(r8) :: wpicehet         ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet

   real(r8) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet 
   !-------------------------------------------------------------------------------

   ! define work variables
   tc = tair-273.15_r8
   qsmall = 1.e-14_r8
   
   ! initialize ni values
   nimix  = 0._r8
   nimh   = 0._r8
   nilphf = 0._r8
   
   !---MIXED-PHASE--------------------------------------------------------------------------

   if ( ( (tc.lt.-35._r8) .and. (tc.ge.-37._r8) .and. (supersat_i.ge.0.05_r8)  .and. (qc.gt.qsmall) ) .or. &
        ( (tc.lt.-32._r8) .and. (tc.ge.-37._r8) .and. (supersat_i.ge.0.005_r8) .and. (qc.gt.qsmall) ) ) then 
        
        ! dep/cond-frzing for MIXED PHASE
        ! following Cooper 1986 or Meyers 1992
        
        if ( do_meyers .eq. .true. ) then
           ! deposition/condensation nucleation in mixed clouds (Meyers, 1992)
           nimix = exp(-0.639_r8+12.96_r8.*supersat_i)*1.e-3_r8  ! cm-3
        else
           ! deposition/condensation nucleation in mixed clouds (Cooper 1986)
           nimix = 0.005_r8*exp(0.304_r8*tc)*1.e-3_r8 ! cm-3
        endif
          
   endif
        
   !---CIRRUS-MOHLER-----------------------------------------------------------------------
   
   if (do_cirrus_mohler_ice_nucleation .eq. .true.) then
   
        ! deposition/condensation-freezing nucleation for CIRRUS
        ! following Mohler et al., 2006 lab results for dust deposition freezing
        ! this should be false when you allow HET frz in LP competition - WHY?
         
        if (tc.lt.-53._r8) then 
           scrit=0.1_r8 !critical supersaturation for T<-53degC
        else
           scrit=0.2_r8
        endif
        
        if ( (tc.lt.-37_r8) .and. (supersat_i.ge.scrit) ) then
           nimoh = dst_num ! #/cm3
        endif 
   endif ! mohler 
   
   !---LP-HOM----------------------------------------------------------------------------
        
   if (do_lphom_ice_nucleation .eq. .true.) then 
         
         ! HOM freezing only
         ! following Liu & Penner, 2005 (LP2005)
         
         if ( (tc.lt.-37._r8)  .and. (supersat_i.ge.0.42_r8) ) then 
            ! BG added some very conservative supi condition not to always go in that loop
            call hf(tc, w, relhum, so4_num, nilphf)
         endif  
      
   endif ! do_lphom_ice_nucleation
   
   !---LP2005-HOM-HET-----------------------------------------------------------------------
      
   if ( do_new_bg_lp_frz .eq. .true. ) then
      ! temp variables that depend on use_preexisting_ice
      wbar1 = wbar
      wbar2 = wbar

      if (use_preexisting_ice) then

         Ni_preice = ni/inv_rho                    ! (convert from #/kg -> #/m3)
         Ni_preice = Ni_preice / max(mincld,cldn)   ! in-cloud ice number density 

!!== KZ_BUGFIX 
         if (Ni_preice > 10.0_r8 .and. qi > 1.e-10_r8 ) then    ! > 0.01/L = 10/m3   
!!== KZ_BUGFIX 
            Shom = -1.5_r8   ! if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, according to Ren & McKenzie, 2005
            lami = (gamma4*ci*ni/qi)**(1._r8/3._r8)
            Ri_preice = 0.5_r8/lami                  ! radius
            Ri_preice = max(Ri_preice, 1e-8_r8)       ! >0.01micron
            call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shom, wpice)
            call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shet, wpicehet)
         else
            wpice    = 0.0_r8
            wpicehet = 0.0_r8
         endif

         weff     = max(wbar-wpice, minweff)
         wpice    = min(wpice, wbar)
         weffhet  = max(wbar-wpicehet,minweff)
         wpicehet = min(wpicehet, wbar)

         wbar1 = weff
         wbar2 = weffhet

         detaT   = wbar/0.23_r8
         RHimean = 1.0_r8
         call frachom(tair, RHimean, detaT, fhom)

      end if

      ni = 0._r8

      ! initialize
      niimm = 0._r8
      nidep = 0._r8
      nihf  = 0._r8
      deles = 0._r8
      esi   = 0._r8
   
   
      if ( (tc.le.-35._r8) .and. (supersat_i.ge.Shet) .and. (wbar1.ge.1.e-6_r8)  ) then 
       
         ! HET vs HOM competition
         ! following LP2005
         ! BG added wbar1 limit for num stability (log of small num problem) 

         A = -1.4938_r8 * log(dst_num) + 12.884_r8
         B = -10.41_r8  * log(dst_num) - 67.69_r8

         regm = A * log(wbar1) + B

         if ( tc .gt. regm ) then
       
            ! heterogeneous nucleation only

            if ( tc.lt.-40._r8 .and. wbar1.gt.1._r8 ) then 
           
               ! exclude T<-40 & W>1m/s from hetero. nucleation

               call hf(tc,wbar1,relhum,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8

               if (use_preexisting_ice) then
                  if (nihf.gt.1e-3_r8) then ! hom occur,  add preexisting ice
                     niimm=min(dst_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                     nihf=nihf + Ni_preice*1e-6_r8 - niimm
                  endif
                  nihf=nihf*fhom
                  n1=nihf+niimm 
               else
                  n1=nihf
               end if

            else

               call hetero(tc,wbar2,dst_num,niimm,nidep)
             
               if (use_preexisting_ice) then
                  if (niimm .gt. 1e-6_r8) then ! het freezing occur, add preexisting ice
                     niimm = niimm + Ni_preice*1e-6_r8
                     niimm = min(dst_num, niimm)        ! niimm < dst_num 
                  end if
               end if
             
               nihf=0._r8
               n1=niimm+nidep

            endif

         else if (tc.lt.regm-5._r8) then
       
            ! homogeneous nucleation only

            call hf(tc,wbar1,relhum,so4_num,nihf)
            niimm=0._r8
            nidep=0._r8

            if (use_preexisting_ice) then
               if (nihf.gt.1e-3_r8) then !  hom occur,  add preexisting ice
                  niimm=min(dst_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                  nihf=nihf + Ni_preice*1e-6_r8 - niimm
               endif
               nihf=nihf*fhom
               n1=nihf+niimm 
            else
               n1=nihf
            end if

         else
       
            ! transition between homogeneous and heterogeneous: interpolate in-between

            if (tc.lt.-40._r8 .and. wbar1.gt.1._r8) then 
          
               ! exclude T<-40 & W>1m/s from hetero. nucleation

               call hf(tc, wbar1, relhum, so4_num, nihf)
               niimm = 0._r8
               nidep = 0._r8

               if (use_preexisting_ice) then
                  if (nihf .gt. 1e-3_r8) then ! hom occur,  add preexisting ice
                     niimm = min(dst_num, Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                     nihf  = nihf + Ni_preice*1e-6_r8 - niimm
                  endif
                  nihf = nihf*fhom
                  n1   = nihf + niimm
               else
                  n1 = nihf
               end if

            else
          
               ! hetero. and homog. nucleation compete

               call hf(regm-5._r8,wbar1,relhum,so4_num,nihf)
               call hetero(regm,wbar2,dst_num,niimm,nidep)

               if (use_preexisting_ice) then
                  nihf = nihf*fhom
               end if

               if (nihf .le. (niimm+nidep)) then
                  n1 = nihf
               else
                  n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
               endif

               if (use_preexisting_ice) then
                  if (n1 .gt. 1e-3_r8) then   ! add preexisting ice
                     n1    = n1 + Ni_preice*1e-6_r8
                     niimm = min(dst_num, n1)  ! assuming all dst_num freezing earlier than hom  !!
                     nihf  = n1 - niimm
                  else
                     n1    = 0._r8
                     niimm = 0._r8
                     nihf  = 0._r8                         
                  endif
               end if
            end if
         end if ! regimes for hom, het, hom vs het

         ni = n1

      end if ! temp, supersat_i, wbar1 thresholding
   end if ! new_bg_lp_frz

   !---end----------------------------------------------------

   ! impose limits
   if ( no_limits .eq. .true.) then
      nimix = min(nimix,150.e+3_r8) !BG increased max limit from 100 to 150/L
      nimoh = min(nimoh,100.e+3_r8) !max to 100/L
      nilphf = min(nilphf,80.e+6_r8)!max to 80,000/L
      ni = min(ni,80.e+6_r8)        !max to 80,000/L
   endif
   
   nuci=ni+nimix+nimoh+nilphf
   
   if ( (nuci.gt.9999._r8) .or. (nuci.lt.0._r8) ) then
      write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
      nuci=0._r8
   endif
   

   ! change unit from #/cm3 to #/kg
   nuci   = nuci*1.e+6_r8*inv_rho
   onimix = nimix*1.e+6_r8*inv_rho
   onidep = (nidep+nimoh)*1.e+6_r8*inv_rho
   oniimm = niimm*1.e+6_r8*inv_rho
   onihf  = (nihf+nilphf)*1.e+6_r8*inv_rho

end subroutine nucleati_bg

!===============================================================================

subroutine hetero(T,ww,Ns,Nis,Nid)

    real(r8), intent(in)  :: T, ww, Ns
    real(r8), intent(out) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_r8
      A12 = -0.0185_r8
      A21 = 2.758_r8
      A22 = 1.3221_r8
      B11 = -0.008_r8
      B12 = -0.0468_r8
      B21 = -0.2667_r8
      B22 = -1.4588_r8

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

end subroutine hetero

!===============================================================================

subroutine hf(T,ww,RH,Na,Ni)

      real(r8), intent(in)  :: T, ww, RH, Na
      real(r8), intent(out) :: Ni

      real(r8)    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real(r8)    A2_fast,B2_fast
      real(r8)    C1_fast,C2_fast,k1_fast,k2_fast
      real(r8)    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real(r8)    C1_slow,C2_slow,k1_slow,k2_slow
      real(r8)    regm
      real(r8)    A,B,C
      real(r8)    RHw

!---------------------------------------------------------------------
! parameters from LP2005 table 1

      A1_fast  =0.0231_r8
      A21_fast =-1.6387_r8  !(T>-64 deg)
      A22_fast =-6.045_r8   !(T<=-64 deg)
      B1_fast  =-0.008_r8
      B21_fast =-0.042_r8   !(T>-64 deg)
      B22_fast =-0.112_r8   !(T<=-64 deg)
      C1_fast  =0.0739_r8
      C2_fast  =1.2372_r8

      A1_slow  =-0.3949_r8
      A2_slow  =1.282_r8
      B1_slow  =-0.0156_r8
      B2_slow  =0.0111_r8
      B3_slow  =0.0217_r8
      C1_slow  =0.120_r8
      C2_slow  =2.312_r8

      Ni = 0.0_r8

!----------------------------
!RHw parameters
      A = 6.0e-4_r8*log(ww)+6.6e-3_r8
      B = 6.0e-2_r8*log(ww)+1.052_r8
      C = 1.68_r8  *log(ww)+129.35_r8
      RHw=(A*T*T+B*T+C)*0.01_r8 ! LP2005 eq 3.1; in this case around RHw = 1.3

      if((T.le.-37.0_r8) .and. ((RH*subgrid).ge.RHw)) then

        regm = 6.07_r8*log(ww)-55.0_r8

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0_r8) then
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

!===============================================================================

SUBROUTINE Vpreice(P_in, T_in, R_in, C_in, S_in, V_out)

   !  based on  Karcher et al. (2006)
   !  VERTICAL VELOCITY CALCULATED FROM DEPOSITIONAL LOSS TERM

   ! SUBROUTINE arguments
   REAL(r8), INTENT(in)  :: P_in       ! [Pa],INITIAL AIR pressure 
   REAL(r8), INTENT(in)  :: T_in       ! [K] ,INITIAL AIR temperature 
   REAL(r8), INTENT(in)  :: R_in       ! [m],INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS 
   REAL(r8), INTENT(in)  :: C_in       ! [m-3],INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY, [1/cm3]
   REAL(r8), INTENT(in)  :: S_in       ! [-],INITIAL ICE SATURATION RATIO;; if <1, use hom threshold Si 
   REAL(r8), INTENT(out) :: V_out      ! [m/s], VERTICAL VELOCITY REDUCTION (caused by preexisting ice)

   ! SUBROUTINE parameters
   REAL(r8), PARAMETER :: ALPHAc  = 0.5_r8 ! density of ice (g/cm3), !!!V is not related to ALPHAc 
   REAL(r8), PARAMETER :: FA1c    = 0.601272523_r8        
   REAL(r8), PARAMETER :: FA2c    = 0.000342181855_r8
   REAL(r8), PARAMETER :: FA3c    = 1.49236645E-12_r8        
   REAL(r8), PARAMETER :: WVP1c   = 3.6E+10_r8   
   REAL(r8), PARAMETER :: WVP2c   = 6145.0_r8
   REAL(r8), PARAMETER :: FVTHc   = 11713803.0_r8
   REAL(r8), PARAMETER :: THOUBKc = 7.24637701E+18_r8
   REAL(r8), PARAMETER :: SVOLc   = 3.23E-23_r8    ! SVOL=XMW/RHOICE
   REAL(r8), PARAMETER :: FDc     = 249.239822_r8
   REAL(r8), PARAMETER :: FPIVOLc = 3.89051704E+23_r8         
   REAL(r8) :: T,P,S,R,C
   REAL(r8) :: A1,A2,A3,B1,B2
   REAL(r8) :: T_1,PICE,FLUX,ALP4,CISAT,DLOSS,VICE

   T = T_in          ! K  , K
   P = P_in*1e-2_r8  ! Pa , hpa

   IF (S_in.LT.1.0_r8) THEN
      S = 2.349_r8 - (T/259.0_r8) ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
   ELSE
      S = S_in                    ! INPUT ICE SATURATION RATIO, -,  >1
   ENDIF

   R     = R_in*1e2_r8   ! m  => cm
   C     = C_in*1e-6_r8  ! m-3 => cm-3
   T_1   = 1.0_r8/ T
   PICE  = WVP1c * EXP(-(WVP2c*T_1))
   ALP4  = 0.25_r8 * ALPHAc      
   FLUX  = ALP4 * SQRT(FVTHc*T)
   CISAT = THOUBKc * PICE * T_1   
   A1    = ( FA1c * T_1 - FA2c ) * T_1 
   A2    = 1.0_r8/ CISAT      
   A3    = FA3c * T_1 / P
   B1    = FLUX * SVOLc * CISAT * ( S-1.0_r8 ) 
   B2    = FLUX * FDc * P * T_1**1.94_r8 
   DLOSS = FPIVOLc * C * B1 * R**2 / ( 1.0_r8+ B2 * R )         
   VICE  = ( A2 + A3 * S ) * DLOSS / ( A1 * S )  ! 2006,(19)
   V_out = VICE*1e-2_r8  ! cm/s => m/s

END SUBROUTINE Vpreice

subroutine frachom(Tmean,RHimean,detaT,fhom)
   ! How much fraction of cirrus might reach Shom  
   ! base on "A cirrus cloud scheme for general circulation models",
   ! B. Karcher and U. Burkhardt 2008

   real(r8), intent(in)  :: Tmean, RHimean, detaT
   real(r8), intent(out) :: fhom

   real(r8), parameter :: seta = 6132.9_r8  ! K
   integer,  parameter :: Nbin=200          ! (Tmean - 3*detaT, Tmean + 3*detaT)

   real(r8) :: PDF_T(Nbin)    ! temperature PDF;  ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
   real(r8) :: Sbin(Nbin)     ! the fluctuations of Si that are driven by the T variations 
   real(r8) :: Sihom, deta
   integer  :: i

   Sihom = 2.349_r8-Tmean/259.0_r8   ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
   fhom  = 0.0_r8

   do i = Nbin, 1, -1

      deta     = (i - 0.5_r8 - Nbin/2)*6.0_r8/Nbin   ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
      Sbin(i)  = RHimean*exp(deta*detaT*seta/Tmean**2.0_r8)
      PDF_T(i) = exp(-deta**2.0_r8/2.0_r8)*6.0_r8/(sqrt(2.0_r8*Pi)*Nbin)
      

      if (Sbin(i).ge.Sihom) then
         fhom = fhom + PDF_T(i)
      else
         exit
      end if
   end do

   fhom=fhom/0.997_r8   ! accounting for the finite limits (-3 , 3)

end subroutine frachom

end module nucleate_ice

