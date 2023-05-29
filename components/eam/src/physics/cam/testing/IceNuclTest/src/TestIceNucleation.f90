program TestIceNucleation

  ! get real kind from utils
  use physics_utils, only: rtype,rtype8,btype
  use IceNucleationOriginal, only: ice_nucleation

  implicit none

  integer :: i, j, k
  real(rtype) :: inv_rho, ni, ni_activated, inv_dt, t_atm, qv_supersat_i
  real(rtype) :: qinuc, ni_nucleat_tend
  logical(btype) :: do_predict_nc, do_prescribed_CCN
  real(rtype) :: qq(10,18), nnuc(10,18)


  inv_rho = 1.
  ni = 0.
  ni_activated = 0.
  inv_dt = 10.
  do_predict_nc = .false.
  do_prescribed_CCN = .true.

  do j = 1,18 ! loop over temperatures 
    t_atm = 270. - 5.*real(j-1)
    do i = 1,10 ! loop over supersaturations
      qv_supersat_i = 0.0 + 0.05*real(i-1)

      
      call ice_nucleation(t_atm,inv_rho,&
       ni,ni_activated,qv_supersat_i,inv_dt,do_predict_nc, do_prescribed_CCN, &
       qinuc, ni_nucleat_tend)

      qq(i,j) = qinuc/inv_dt
      nnuc(i,j) = ni_nucleat_tend/inv_dt

      if(i+j.eq.2) write(*,*) 'T, Sice, Nnuc_old, qnuc_old '
      write(*,992) t_atm, qv_supersat_i, nnuc(i,j), qq(i,j)
992 format(2f8.3,4e12.4)
    end do
  end do

end program TestIceNucleation
