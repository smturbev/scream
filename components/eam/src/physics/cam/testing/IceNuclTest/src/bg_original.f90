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

subroutine ice_nucleation(t, supi, qc, nitot, iSCF, &
        do_meyers, frzmodes, &
        ) 
        
        dum, dum1, D_new, N_nuc, Q_nuc, nCat, mi0, qsmall, iice_dest
        
        mi0        = 4.*3.14159265/3.*900.*1.e-18
        qsmall     = 1.e-14 ! min allowable prognostic variables
        inv_dt     = 1./dt   ! inverse model time step
        diam_ice   = 0. ! initialize
        iice_dest = -99
        
        diam_ice(i,k,iice) = ((qitot(i,k,iice)*6.)/(nitot(i,k,iice)*dum2*pi))**thrd


  !................................................................
  ! 1.) deposition/condensation-freezing nucleation for MIXED PHASE
  ! allow ice nucleation if < -15 C and > 5% ice supersaturation
  !BG added a min temp condition (-38 deg C, hom frz threshold)
  !BG added a qc>qsmall condition: based on recent findings, deposition freezing is negligible in mixed phase (e.g. Ansmann et al., 2018)
  !BG therefore in mixed phase we can freeze only by immersion and contact and maybe also condensation freezing 
  !BG (e.g. Hoose and Mohler, 2012, Lohmann et al., 2016)
  !if ( (t(i,k).lt.258.15) .and. (t(i,k) .gt. 236.15) .and. (supi(i,k).ge.0.05) .and. (qc(i,k) .gt. qsmall) ) then
  if ( ( (t(i,k).lt.258.15) .and. (t(i,k) .ge. 236.15) .and. (supi(i,k).ge.0.05) .and. (qc(i,k) .gt. qsmall) ) &
  .or. ( (t(i,k).lt.241.15) .and. (t(i,k) .ge. 236.15) .and. (supi(i,k).ge.0.005) .and. (qc(i,k) .gt. qsmall) ) ) then  
  !BG added this to prevent too much "hom frz at -37"
   !print*,'frz in mix phase'!BG
    if (do_meyers) then !BG
     dum = exp(-0.639+0.1296*100.*supi(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992) per m3
   else
     dum = 0.005*exp(0.304*(273.15-t(i,k)))*1000.*inv_rho(i,k)   !Cooper (1986)
   endif !BG Meyers of Cooper?
     dum = min(dum,150.e3*inv_rho(i,k)*SCF(k)) !BG increased max limit from 100 to 150/L
     N_nuc =max(0.,(dum-sum(nitot(i,k,:)))*inv_dt)

     if (N_nuc.ge.1.e-20) then
        Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*inv_dt)
        iice_dest = 1

      !BG
      if (no_ice_nucleation .or. no_het_ice_nucleation) then
        qinuc(iice_dest) = 0.
        ninuc(iice_dest) = 0.
      else
        qinuc(iice_dest) = Q_nuc !BG =0. 1
        ninuc(iice_dest) = N_nuc !BG =0. 2
      end if
      !BG end

     endif

  endif !(t(i,k).lt.258.15 .and. supi(i,k).ge.0.05 .and. (qc(i,k) .gt. qsmall))

  if (no_cirrus_mohler_ice_nucleation .eq. .false.) then
     ! 2.) deposition/condensation-freezing nucleation for CIRRUS
     ! following Mohler et al., 2006 lab results for dust deposition freezing
     if (t(i,k) .lt. 220.15) then
        scrit=0.1 !critical supersaturation
     else
        scrit=0.2
     endif

     if ( (t(i,k).lt.236.15) .and. (supi(i,k).ge. scrit) ) then

        dum = ndust *1e6*inv_rho(i,k) !from /cm3 to kg-1 !assume some small INP/dust concentration, say 2/L which all freeze by deposition freezing
        dum = min(dum,100.e3*inv_rho(i,k)*SCF(k)) !max to 100/liter
        N_nuc =max(0.,(dum-sum(nitot(i,k,:)))*inv_dt)

        if (N_nuc.ge.1.e-20) then
           Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*inv_dt)
           iice_dest = 1
         !BG
         if (no_ice_nucleation .or. no_het_ice_nucleation) then
           qinuc2(iice_dest) = 0.
           ninuc2(iice_dest) = 0.
         else
           qinuc2(iice_dest) = Q_nuc
           ninuc2(iice_dest) = N_nuc
         end if
         !BG end

        endif !N_nuc .ge. 1e-20

    endif !(t(i,k).lt.258.15  .and. gt 236 .and. supi(i,k).ge.0.05)
  end if ! (no_cirrus_mohler_ice_nucleation)

   if (no_lphom_ice_nucleation .eq. .false.) then
     !3.) HOM nucleation by using Liu and Penner, 2005 => ONLY HOMOG NUCLEATION!!!
     if ( (t(i,k).lt.236.15)  .and. (supi_cld.ge. 0.42)) then !added some very conservative supi condition not to always go in that loop
        call hf(tc(i,k),uzpl(i,k),relhum(i,k),nsulf,nihf)

        dum = nihf*1.e+6*inv_rho(i,k) !from cm-3 to m-3 to kg-1 air
        dum = min(dum,80000.e+3*inv_rho(i,k)*SCF(k)) !set max to 80000 per L or 80/cc
        N_nuc =max(0.,(dum-sum(nitot(i,k,:)))*inv_dt)

        if (N_nuc.ge.1.e-20) then
           Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*inv_dt)
           iice_dest = 1

         !BG
         if (no_ice_nucleation) then
           qinuc3(iice_dest) = 0.
           ninuc3(iice_dest) = 0.
         else
           qinuc3(iice_dest) = Q_nuc
           ninuc3(iice_dest) = N_nuc
         end if
         !BG end

        endif

      endif !(t(i,k).lt.236.15)  .and. (supi(i,k).ge. 0.42) !hom freezing in cirrus - liu penner
   endif ! (no_lphom_ice_nucleation)

   !BG copied from e3sm code
   !competition with HOMOG + HETEROG NUCLEATION (+preex ice if set to true)
   if (do_new_lp_frz) then
     ! temp variables that depend on use_preexisting_ice
     wbar1 = uzpl(i,k)
     wbar2 = uzpl(i,k)
     !===========================
    if (use_preexisting_ice) then
    !The pre-existing ice formulation is based on Karcher et al., 2006, and is taken from the code implementation (in CAM5, taken actually from E3SM v1)
    !based on Shi et al., 2015, ACP

       Ni_preice = sum(nitot(i,k,:))*rho(i,k)                    ! (convert from #/kg -> #/m3)
       !BG no need for this here Ni_preice = Ni_preice / max(mincld,cldn)   ! in-cloud ice number density


       if (Ni_preice > 10.0 .and. sum(qitot(i,k,:)) > 1.e-10 ) then    ! > 0.01/L = 10/m3

          Shom = -1.5   ! if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, according to Ren & McKenzie, 2005
          ci = rhoice*pi/6.
          lami = (gamma4*ci*sum(nitot(i,k,:))/sum(qitot(i,k,:)) )**(1./3.) !BG sum because of multiple categories
          Ri_preice = 0.5/lami                   ! radius  BG i think assuming sphere
          Ri_preice = max(Ri_preice, 1e-8)       ! >0.01micron
          call Vpreice(pres(i,k), t(i,k), Ri_preice, Ni_preice, Shom, wpice)  !BG: pres instead of pmid -> i think is ok
          call Vpreice(pres(i,k), t(i,k), Ri_preice, Ni_preice, Shet, wpicehet)
       else
          wpice    = 0.0
          wpicehet = 0.0
       endif

       weff     = max(uzpl(i,k)-wpice, minweff)
       wpice    = min(wpice, uzpl(i,k))
       weffhet  = max(uzpl(i,k)-wpicehet,minweff)
       wpicehet = min(wpicehet, uzpl(i,k))
       !if ((weff .gt. 0.1) .and. (Ni_preice .gt. 50000.)) then
       !    print*,weff,'weff'
       !    print*,wpice,'wpice'
       !    print*,weffhet,'weffhet'
       !    print*,wpicehet,'wpicehet'
       !    print*,uzpl(i,k),'updraft'
       !    print*,sum(nitot(i,k,:))*rho(i,k)*1e-3,'nitot per L'
       !end if
       wbar1 = weff
       wbar2 = weffhet

       !BG added a switch to activate mesoscale variability inside of existing cirrus (with RHimean set for some reason to 1
       !fhom rarely reaches high values close to 1  => assumes big grid boxes, GCM case! 
       !it was first implemented in CAM5, where they were having troubles by generating too many ice crystals
       !so while the freezing itself considers only the portion of the box that is "freezing", the fhom is used to "mix that" throughout the gridbox,
       ! decreasing the number of nucleated ICs
       ! if needed, it can be easiliy modified to act as an additional source of waves/updraft/supersat
       ! currently it should not be used at horizontal resolutions below 10 km or so, but it may make sense if we switch subgrid cloud fraction on
       ! even at smaller horizontal resolutions
       !the logic behind is better described in Shi et al., 2015, ACP section 2.3
       if (do_mesoscale_variab) then
          detaT   = uzpl(i,k)/0.23
          !BG: this is what E3SM used: RHimean = 1.0
          RHimean = 1.0 + supi_cld !BG use calculated supersat value
          call frachom(t(i,k), RHimean, detaT, fhom)
       else
          fhom = 1.0 !for cloud fractions 0/1 this makes most sense
       end if

    end if !preexisting ice
    !========================

   ni = 0.
    !temp in celsius tc = tair - 273.15_r8

    ! initialize
    niimm = 0.
    nidep = 0.
    nihf  = 0.
    deles = 0.
    esi   = 0.

    !? BG why cld cover? if(so4_num >= 1.0e-10_r8 .and. (soot_num+dst3_num) >= 1.0e-10_r8 .and. cldn > 0._r8) then

    !if((tc(i,k).le.-35.0) .and. (supi(i,k).ge. 0.2) ) then ! use higher RHi threshold ?BG?
    !120% is here taken as heterogeneous freezing limit
    if((tc(i,k).le.-35.0) .and. (supi_cld.ge. 0.2) .and. (wbar1 .ge. 1.e-6)  ) then
    !BG added wbar1>1e-6 for num stability (log of small number->problems)

            A = -1.4938 * log(ndust) + 12.884
            B = -10.41  * log(ndust) - 67.69

            regm = A * log(wbar1) + B

            !print*,regm,'regmstart'
            !print*,wbar1,'wbar1start'
            !print*,'log(wbar1start)',log(wbar1)

            ! heterogeneous nucleation only
            if (tc(i,k) .gt. regm) then
            !BG critical soot (in our case dust/some undefined HET ice nucleating particle) num conc above which only 
            !   HET freezing occurs, Liu et al., 2007, eq 10.
            !   it depends on INP number and W!

               if(tc(i,k).lt.-40. .and. wbar1.gt.1.) then ! exclude T<-40 & W>1m/s from hetero. nucleation
                  call hf(tc(i,k),wbar1,relhum(i,k),nsulf,nihf)
                  niimm=0.
                  nidep=0.

                  if (use_preexisting_ice) then
                     if (nihf.gt.1e-3) then ! hom occur,  add preexisting ice
                        niimm=min(ndust,Ni_preice*1e-6)       ! assuming ndust freeze first
                        nihf=nihf + Ni_preice*1e-6 - niimm   !all in cm-3 units
                     endif
                     nihf=nihf*fhom
                     n1=nihf+niimm
                  else
                     n1=nihf
                  end if !preex

               else

                  call hetero(tc(i,k),wbar2,ndust,niimm,nidep)
                  if (use_preexisting_ice) then
                     if (niimm .gt. 1e-6) then ! het freezing occur, add preexisting ice
                        niimm = niimm + Ni_preice*1e-6
                        niimm = min(ndust, niimm)        ! niimm < ndust
                     end if
                  end if !preex
                  nihf=0.
                  n1=niimm+nidep

               endif !(tc(i,k).lt.-40. .and. wbar1.gt.1.)

           ! homogeneous nucleation only
            else if (tc(i,k).lt.regm-5.) then
               !print*,tc(i,k),'tc before homfrz2'
               call hf(tc(i,k),wbar1,relhum(i,k),nsulf,nihf)
               niimm=0.
               nidep=0.

               if (use_preexisting_ice) then
                  if (nihf.gt.1e-3) then !  hom occur,  add preexisting ice
                     niimm=min(ndust,Ni_preice*1e-6)       ! assuming dst_num freeze firstly
                     nihf=nihf + Ni_preice*1e-6 - niimm
                  endif
                  nihf=nihf*fhom
                  n1=nihf+niimm
               else
                  n1=nihf
               end if

            ! transition between homogeneous and heterogeneous: interpolate in-between
            else

               if (tc(i,k).lt.-40. .and. wbar1.gt.1.) then ! exclude T<-40 & W>1m/s from hetero. nucleation
                  call hf(tc(i,k),wbar1,relhum(i,k),nsulf,nihf)
                  niimm = 0.
                  nidep = 0.

                  if (use_preexisting_ice) then
                     if (nihf .gt. 1e-3) then ! hom occur,  add preexisting ice
                        niimm = min(ndust, Ni_preice*1e-6)       ! assuming dst_num freeze firstly
                        nihf  = nihf + Ni_preice*1e-6 - niimm
                     endif
                     nihf = nihf*fhom
                     n1   = nihf + niimm
                  else
                     n1 = nihf
                  end if

               else  !this is now the "real" transitional regime with calls for both homogeneous and heterogeneous nucleation 

                  call hf(regm-5.,wbar1,relhum(i,k),nsulf,nihf)
                  call hetero(regm,wbar2,ndust,niimm,nidep) !BG the way I programmed it, nidep is 0 by definition

                  if (use_preexisting_ice) then
                     nihf = nihf*fhom
                  end if

                  if (nihf .le. (niimm+nidep)) then
                     n1 = nihf
                  else
                     n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc(i,k)-regm)/5.)
                  endif

                  if (use_preexisting_ice) then
                     if (n1 .gt. 1e-3) then   ! add preexisting ice
                        n1    = n1 + Ni_preice*1e-6
                        niimm = min(ndust, n1)  ! assuming all ndust freezing earlier than hom  !!
                        nihf  = n1 - niimm
                     else
                        n1    = 0.
                        niimm = 0.
                        nihf  = 0.
                     endif
                  end if

               end if
            end if

            ni = n1 !cm-3 
              !if (nihf .gt. 0.) then
                   !print*, niimm,'niimm end'
                   !print*, nihf,'nihf end'
                   !print*, fhom, 'fhom end'
                   !print*, supi_cld, 'supi end'
                   !print*, uzpl(i,k), 'w end'
                   !print*, tc(i,k),'tc end'
              !end if
         end if !((tc(i,k).le.-35.0) .and. (supi(i,k).ge. 0.2) .and. (wbar1 .ge. 1.e-6)  )


   dum = ni*1.e+6*inv_rho(i,k) !from cm-3 to m-3 to kg-1 air
   dum = min(dum,80.e+6*inv_rho(i,k)*SCF(k)) !set max to 80000 per L or 80/cc

  !N_nuc =max(0.,(dum-sum(nitot(i,k,:)))*inv_dt)
  !BG I think there is no logic behind the upper statement in case we do "real" freezing. 
  !   It may be reasonable for Meyers/Cooper parameterization, but not here...in my understanding!!!
  !BG changed therefore to: 
   N_nuc = dum*inv_dt

   if (N_nuc.ge.1.e-20) then

      !BG Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*inv_dt)
     Q_nuc = dum*mi0*inv_dt !BG changed!!!
     if (nCat>1) then
             if (frzmodes) then !BG freezing separate category separate source
                iice_dest=5   !LP with preex ice hom+het freezing= cat 5 (TO DO: separate HET FRZ????) 
             else
               !determine destination ice-phase category:
               dum1  = 900.     !density of new ice
               D_new = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
               call icecat_destination(qitot(i,k,:)*iSCF(k),diam_ice(i,k,:),D_new,deltaD_init,iice_dest)

               if (global_status /= STATUS_OK) return
             endif !BG, frzmodes
     else !nCat=1
         iice_dest = 1
     endif

    !BG
     if (no_ice_nucleation) then
       qinuc3(iice_dest) = 0.
       ninuc3(iice_dest) = 0.
     else
       qinuc3(iice_dest) = Q_nuc
       ninuc3(iice_dest) = N_nuc
     end if
     !BG end

   endif !no ice nucleation



  end if !new_lp_frz 
  
  
  end subroutine ice_nucleation
  
  
  !==========================================================================================!
 subroutine icecat_destination(Qi,Di,D_nuc,deltaD_init,iice_dest)

 !--------------------------------------------------------------------------------------!
 ! Returns the index of the destination ice category into which new ice is nucleated.
 !
 ! New ice will be nucleated into the category in which the existing ice is
 ! closest in size to the ice being nucleated.  The exception is that if the
 ! size difference between the nucleated ice and existing ice exceeds a threshold
 ! value for all categories, then ice is initiated into a new category.
 !
 ! D_nuc        = mean diameter of new particles being added to a category
 ! D(i)         = mean diameter of particles in category i
 ! diff(i)      = |D(i) - D_nuc|
 ! deltaD_init  = threshold size difference to consider a new (empty) category
 ! mindiff      = minimum of all diff(i) (for non-empty categories)
 !
 ! POSSIBLE CASES                      DESTINATION CATEGORY
 !---------------                      --------------------
 ! case 1:  all empty                  category 1
 ! case 2:  all full                   category with smallest diff
 ! case 3:  partly full
 !  case 3a:  mindiff <  diff_thrs     category with smallest diff
 !  case 3b:  mindiff >= diff_thrs     first empty category
 !--------------------------------------------------------------------------------------!

 implicit none

! arguments:
 real, intent(in), dimension(:) :: Qi,Di
 real, intent(in)               :: D_nuc,deltaD_init
 integer, intent(out)           :: iice_dest

! local variables:
 logical                        :: all_full,all_empty
 integer                        :: i_firstEmptyCategory,iice,i_mindiff,n_cat
 real                           :: mindiff,diff
 real, parameter                :: qsmall_loc = 1.e-14

 !--------------------------------------------------------------------------------------!

 n_cat     = size(Qi)
 iice_dest = -99

!-- test:
! iice_dest = 1
! return
!==

 if (sum(Qi(:))<qsmall_loc) then

 !case 1:
    iice_dest = 1
    return

 else

    all_full  = .true.
    all_empty = .false.
    mindiff   = 9.e+9
    i_firstEmptyCategory = 0

    do iice = 1,n_cat
       if (Qi(iice) .ge. qsmall_loc) then
          all_empty = .false.
          diff      = abs(Di(iice)-D_nuc)
          if (diff .lt. mindiff) then
             mindiff   = diff
             i_mindiff = iice
          endif
       else
          !print*,iice,'init iice first empty cat'
          all_full = .false.
          if (i_firstEmptyCategory.eq.0) i_firstEmptyCategory = iice
          !print*,i_firstEmptyCategory,'init i_firstEmptyCategory'
       endif
    enddo

    if (all_full) then
 !case 2:
       iice_dest = i_mindiff
       !print*,iice_dest,'2 iice_dest when all_full'
       return
    else
 !BG commented out the standard way of starting a new ice cat, with a fixed theshold
 !   and tried a factor based criterion: where deltaD_init is now a scaling factor, e.g. 3
 !   (i.e. if Dnew>3*Dold, or Dnew<1/3*Dold, make a new category)

       !if (mindiff .lt. deltaD_init) then 
 !case 3a:
       !   iice_dest = i_mindiff
       !   return
       !else
 !case 3b:
       !   iice_dest = i_firstEmptyCategory
       !   return
       !endif
 !      if ( (D_nuc .gt. deltaD_init*Di(iice)) .or. (D_nuc .lt. Di(iice)/deltaD_init) ) then !BG
 !case 3b:
 !         print*,'in new cat'
 !         iice_dest = i_firstEmptyCategory
 !         print*,iice_dest,'3b iice_dest new cat'
 !         print*, D_nuc,'3b d new, in new cat'
 !         print*, Di(iice),'3b d preex in new cat'
 !         return
 !      else
 !case 3a:
 !         print*, 'in old cat'
 !         iice_dest = i_mindiff
 !         print*,iice_dest,'3a iice_dest old cat'
 !         print*, D_nuc,'3a d new, in old cat'
 !         print*, Di(iice),'3a d preex in old cat'
 !         return
 !      endif
 !   endif  !all full

        if ( (D_nuc .gt. (1./deltaD_init)*Di(iice)) .and. (D_nuc .lt. deltaD_init*Di(iice) ) ) then !BG
 !case 3a:
          !print*,'3a in old cat'
          iice_dest = i_mindiff
          !print*,iice_dest,'3a iice_dest old cat'
          !print*, D_nuc,'3a d new, in old cat'
          !print*, Di(iice),'3a d preex in old cat'
          return
       else
 !case 3b:
          !print*, '3b in old cat'
          iice_dest = i_firstEmptyCategory
          !print*,iice_dest,'3b iice_dest new cat'
          !print*, D_nuc,'3b d new, in new cat'
          !print*, Di(iice),'3b d preex in new cat'
          return
       endif
    endif  !all full

 endif ! loop (sum(Qi(:))<qsmall_loc) if

 print*, 'ERROR in s/r icecat_destination -- made it to end'
 global_status = STATUS_ERROR
 return

 end subroutine icecat_destination

