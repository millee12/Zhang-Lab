module run_variables
  implicit none
  save
  real(8) :: dt
  integer :: its,nts,nts_start
  real(8) :: tt,T
  integer :: ntsbout,N
  integer :: restart_freq,restart_klok,restart_unit
  integer,parameter :: restart_u1 = 9611, restart_u2 = 9612
end module run_variables
