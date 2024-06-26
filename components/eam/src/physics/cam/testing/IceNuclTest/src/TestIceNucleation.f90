program TestIceNucleation

  ! get real kind from utils
  use physics_utils, only: rtype,rtype8,btype
  use micro_p3_minimal, only: ice_nucleation,hf,hetero

  implicit none

  integer :: i, j, k, l, m
  real(rtype) :: inv_rho, ni, ni_activated, inv_dt, t_atm, qv_supersat_l, qv_supersat_i, qc, uzpl
  real(rtype) :: qinuc, ni_nucleat_tend
  logical(btype) :: do_predict_nc, do_prescribed_CCN, do_new_lp_freezing, no_cirrus_mohler_ice_nucleation, no_lphom_ice_nucleation
  real(rtype) :: qq, nnuc
  real(rtype) :: nn0,nn1,nn2,nn3,nn4

  inv_rho = 1.
  ni = 0.
  ni_activated = 0.
  inv_dt = 10.
  ! qv_supersat_l = -0.2 ! RHw = 0.8
  do_predict_nc = .false.
  do_prescribed_CCN = .false.
  do_new_lp_freezing = .false.
  no_cirrus_mohler_ice_nucleation = .false.
  no_lphom_ice_nucleation = .false.
    
  do m = 1,6 ! loop over S_liq
    qv_supersat_l = -0.3 + 0.1*real(m-1)
    do l = 1,4,3 ! loop over vertical velocity, uzpl (m/s)
      uzpl = 0.5*real(l)
      do k = 1,3! 0,0.0001,0.001 ! loop of cloud mixing ratios, qc
        if ( int(k).eq.1 ) then
          qc=0
        else if ( int(k).eq.2 ) then
          qc=0.0001
        else
          qc=0.001
        endif
        do j = 1,18 ! loop over temperatures 
          t_atm = 270. - 5.*real(j-1)
          do i = 1,10 ! loop over supersaturations
            qv_supersat_i = 0.0 + 0.05*real(i-1)
       
            call ice_nucleation(t_atm,inv_rho,&
             ni,ni_activated,qv_supersat_l,& 
             qv_supersat_i,inv_dt,qc,uzpl,&
             do_predict_nc,do_prescribed_CCN,&
             do_new_lp_freezing,no_cirrus_mohler_ice_nucleation,&
             no_lphom_ice_nucleation,qinuc,ni_nucleat_tend,&
             nn0,nn1,nn2,nn3,nn4)

            qq = qinuc/inv_dt
            nnuc = ni_nucleat_tend/inv_dt
            nn0=nn0/inv_dt
            nn1=nn1/inv_dt
            nn2=nn2/inv_dt
            nn3=nn3/inv_dt
            nn4=nn4/inv_dt

            if(i+j+l+k+m.eq.5) write(*,*) 'w  qc  T  Sice  Sliq  Nnuc_old  qnuc_old  nnuc0  nnuc1  nnuc2  nnuc3  nnuc4'
            write(*,992) uzpl, qc, t_atm, qv_supersat_i, qv_supersat_l, nnuc, qq, nn0, nn1, nn2, nn3, nn4
992 format(6f14.3, 1e14.4, 5f14.2)
          end do
        end do
      end do
    end do
  end do
end program TestIceNucleation
