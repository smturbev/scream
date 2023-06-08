module micro_p3_minimal
 
   ! get real kind from utils
   use physics_utils, only: rtype,rtype8,btype

   ! physical and mathematical constants
!    use micro_p3_utils, only! : mi0, T_zerodegc, T_icenuc, NumCirrusINP, NumCirrusSulf

   implicit none

   public :: ice_nucleation
   
   ! Constants from micro_p3_utils
     real(rtype), parameter :: mi0         = 3.77e-14 ! 4._rtype*piov3*900._rtype*1.e-18_rtype
     real(rtype), parameter :: T_zerodegc  = 273.15_rtype
     real(rtype), parameter :: T_icenuc    = 273.15_rtype - 15._rtype
     real(rtype), parameter :: NumCirrusINP= 2.0e-3
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
   qinuc, ni_nucleat_tend)

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

   ! local variables
   real(rtype) :: dum, N_nuc, Q_nuc 
   real(rtype) :: ndust, nsulf, qsmall, w, niimm, nidep, nihf, scrit ! for new BG code
   real(rtype) :: wbar1, wbar2, deles, esi, A, B, regm, n1 ! work variables
   
   ! convert uzpl to w (Pa/s -> m/s)
   !w = - uzpl * inv_rho * 0.102 ! omega * 1/rho * 1/g  [m/s]
   
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
            dum = 0.005_rtype*exp(0.304_rtype*(T_zerodegc-t_atm))*1000._rtype*inv_rho
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
   else ! do_new_lp_freezing .eq. .true.
	  ! Adding three main components:
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
      if ( ( (t_atm.lt.258.15) .and. (t_atm .ge. 236.15) .and. (qv_supersat_i.ge.0.05) .and. (qc .gt. qsmall) ) .or. ( (t_atm.lt.241.15) .and. (t_atm .ge. 236.15) .and. (qv_supersat_i.ge.0.005) .and. (qc .gt. qsmall) ) ) then  !BG added this to prevent too much "hom frz at -37"
         ! freezing in mixed phase
         dum = min(dum,150.e3*inv_rho) !BG increased max limit from 100 to 150/L
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
            ! T < -37degC and supersat > critical value
            dum = ndust*1e6*inv_rho !from /cm3 to kg-1 !assume some small INP/dust concentration, say 2/L which all freeze by deposition freezing
            dum = min(dum,100.e3*inv_rho) !max to 100/liter
            N_nuc =max(0.,(dum-ni)*inv_dt)
            if (N_nuc.ge.1.e-20) then
               Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt) 
               qinuc = Q_nuc
               ni_nucleat_tend = N_nuc
            endif
         endif 
      endif ! no_cirrus_mohler_ice_nucleation .eq. .false. 
	  
      if (no_lphom_ice_nucleation .eq. .false.) then 
         ! 3.) HOM nucleation by using Liu and Penner, 2005 => ONLY HOMOG NUCLEATION!!!
         if ( (t_atm.lt.236.15)  .and. (qv_supersat_i .ge. 0.42) ) then !added some very conservative supi condition not to always go in that loop ! 
             call hf(t_atm, w, qv_supersat_l, nsulf, nihf)
			 
			 dum = nihf*1.e+6*inv_rho !from cm-3 to m-3 to kg-1 air
        	 dum = min(dum,80000.e+3*inv_rho) !set max to 80000 per L or 80/cc
             N_nuc =max(0.,(dum-ni)*inv_dt)
             if (N_nuc.ge.1.e-20) then
                Q_nuc = max(0._rtype,(dum-ni)*mi0*inv_dt)
                qinuc = Q_nuc
                ni_nucleat_tend = N_nuc
             endif
			 
         endif ! ( (t_atm.lt.236.15)  .and. (qv_supersat_i .ge. 0.42) )
      endif ! no_lphom_ice_nucleation .eq. .false.
	  !ST copied from BG copied from e3sm code
	  ! 4.) competition with HOMOG + HETEROG NUCLEATION (+preex ice if set to true)
	  
	  ! temp variables that depend on use_preexisting_ice
	  wbar1 = w
	  wbar2 = w
	  ! TODO: add pre-existing ice section from BG
  
	  ! lp freezing
      !temp in celsius tc = tair - 273.15_r8
	  ! TODO: ST- check that all temperature has been converted back to Kelvin

      ! initialize
      niimm = 0.
      nidep = 0.
      nihf  = 0.
      deles = 0.
      esi   = 0.
  	      
      if((t_atm.le.238.15) .and. (qv_supersat_i.ge.0.2) .and. (wbar1 .ge. 1.e-6)  ) then 
	  ! 120% is here taken as heterogeneous freezing limit
	  !? BG - use higher RHi threshold? 
      !BG added wbar1>1e-6 for num stability (log of small number->problems)
  		  A = -1.4938 * log(ndust) + 12.884
          B = -10.41  * log(ndust) - 205.46 ! ST convert back to Kelvin -67.69 -> +205.46

          regm = A * log(wbar1) + B ! regm from LP2005 eq. 4.5

		  ! 4a.) heterogeneous nucleation only
          if (t_atm .gt. regm) then
              !	BG critical soot (in our case dust/some undefined HET ice nucleating particle) 
			  ! ndust num conc above which only 
              ! HET freezing occurs, it depends on INP number and w
				if((t_atm.lt.233.15) .and. (wbar1.gt.1.)) then ! exclude T<-40 & W>1m/s from hetero. nucleation
					call hf(t_atm,wbar1,qv_supersat_l,nsulf,nihf)
                    niimm=0.
                    nidep=0.
					n1=nihf ! TODO: add preexisting ice condition
				else ! heterogeneous freezing
					call hetero(t_atm,wbar2,ndust,niimm,nidep)
					! TODO: add preexisting ice condition
					nihf=0.
					n1=niimm+nidep
				endif !(tc.lt.-40. .and. wbar1.gt.1.)
		  ! 4b.) homogeneous nucleation only
		  else if (t_atm.lt.regm-5.) then
               call hf(t_atm,wbar1,qv_supersat_l,nsulf,nihf)
               niimm=0.
               nidep=0.
               n1=nihf ! TODO: add preexisting ice condition
		  ! 4c.) transition between homogeneous and heterogeneous: interpolate in-between
		  else
			  if ((t_atm.lt.233.15) .and. (wbar1.gt.1.)) then ! exclude T<-40 & W>1m/s from hetero. nucleation
			  		call hf(t_atm,wbar1,qv_supersat_l,nsulf,nihf)
                    niimm=0.
                    nidep=0.
					n1=nihf ! TODO: add preexisting ice condition
              else  !this is now the "real" transitional regime with calls for both homogeneous and heterogeneous nucleation
				  call hf(regm-5.,wbar1,qv_supersat_l,nsulf,nihf)
				  call hetero(regm,wbar2,ndust,niimm,nidep) !BG the way I programmed it, nidep is 0 by definition
				  if (nihf .le. (niimm+nidep)) then
                       n1 = nihf
                  else
                       n1=(niimm+nidep)*((niimm+nidep)/nihf)**((t_atm-regm)/5.)
                  endif
			  end if
		  end if ! hetero vs homog frz
		  N_nuc = n1 !cm-3 
	 end if  ! ((t_atm.le.238.15) .and. (qv_supersat_i.ge.0.2) .and. (wbar1 .ge. 1.e-6)  )
	 dum = N_nuc*1.e+6*inv_rho !from cm-3 to m-3 to kg-1 air
	 dum = min(dum,80.e+6*inv_rho) !set max to 80000 per L or 80/cc
 
	 !N_nuc =max(0.,(dum-sum(nitot(i,k,:)))*odt)
	 !BG I think there is no logic behind the upper statement in case we do "real" freezing. 
     !   It may be reasonable for Meyers/Cooper parameterization, but not here...in my understanding!!!
     !BG changed therefore to: 
     N_nuc = dum*inv_dt
 
	 if (N_nuc.ge.1.e-20) then
		 ! BG Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*odt)
	     Q_nuc = dum*mi0*inv_dt !BG changed!!!
		 qinuc = Q_nuc
		 ni_nucleat_tend = N_nuc
	 endif ! no ice nucleation
   end if ! do_new_lp_frz	  
end subroutine

!BG added homo freezing by Liu and Penner, 2005
!===============================================================================
!it assumes a fast and slow growth regime based on RHw and temperature.
!it also needs sulphate aerosol number, which we prescribe to 30/cc.
! ST changed input to tk (temp in  Kelvin) and supersat wrt liquid
subroutine hf(tk,ww,supersat,Na,Ni) !t_atm , w, qv_supersat_l, nsulf, nihom
       real(rtype), intent(in)   :: tk, ww, supersat, Na
       real(rtype), intent(out)  :: Ni

       real(rtype) ::    T,A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
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

       A1_fast  =0.0231
       A21_fast =-1.6387  !(T>-64 degC)
       A22_fast =-6.045   !(T<=-64 degC)
       B1_fast  =-0.008
       B21_fast =-0.042   !(T>-64 degC)
       B22_fast =-0.112   !(T<=-64 degC)
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
	   T = tk-273.15 ! convert to degC

 !----------------------------
 !RHw parameters
       A = 6.0e-4*log(ww)+6.6e-3
       B = 6.0e-2*log(ww)+1.052
       C = 1.68  *log(ww)+129.35
       Sw=((A*T*T+B*T+C)*0.01)-1 

       if((T.le.-37.0) .and. ((supersat).ge.Sw)) then
         !print*,'RHw crit in hf', RHw !BG
         !print*,'RH in hf', RH  !BG
         !print*, 'wind in hf' ,ww !BG
         regm = 6.07*log(ww)-55.0 ! from LP2005 eq. 3.5

         if (T.ge.regm) then ! fast-growth regime

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

!-----BG added for ice nucleation-----
! ST changed input to tk (temp in Kelvin)
subroutine hetero(tk,ww,Ns,Nis,Nid)

    real(rtype), intent(in)  :: tk, ww, Ns
    real(rtype), intent(out) :: Nis, Nid

    real(rtype) :: A11,A12,A21,A22,B11,B12,B21,B22
    real(rtype) :: T,B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263
      A12 = -0.0185
      A21 = 2.758
      A22 = 1.3221
      B11 = -0.008
      B12 = -0.0468
      B21 = -0.2667
      B22 = -1.4588
	  T = tk - 273.15

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0    ! don't include deposition nucleation for cirrus clouds when T<-37C
                   ! BG we assume the current het freezing represents in some sense both het by imm and depo
end subroutine hetero

end module micro_p3_minimal