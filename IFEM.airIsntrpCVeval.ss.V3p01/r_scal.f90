!     calculation of fp
subroutine r_scalfp(fp,ocpp,i)
  use r_common, only: bpre,cpre,hp
  implicit none 

  real(8) :: fp,ocpp
  integer :: i

  fp=-ocpp*(bpre-cpre)*hp(i)
  return
end subroutine r_scalfp

!     calculation of fu
subroutine r_scalfu(fu,isd,ni)
  use r_common
  use solid_variables, only : nsd_solid
  implicit none

  real(8) :: fu
  integer :: isd,ni
  integer :: ksd

  fu=0.0d0
  do ksd = 1,nsd_solid
    fu = fu + bd(ksd,ni)*PK1str_tens(ksd,isd)
  enddo
  return
end subroutine r_scalfu

!     calculation of fu in current configuration
subroutine r_scalfu_curr(fu,i,ni,cstr_element)
  use r_common
  use solid_variables, only : nsd_solid
  implicit none

  real(8) :: fu
  integer :: i,ni
  real(8),intent(in) :: cstr_element(2*nsd_solid)  !...Cauchy stress

  fu=0.0d0
  if (i == 1) then
     fu=fu+cstr_element(1)*bd_curr(1,ni)
     fu=fu+cstr_element(6)*bd_curr(2,ni)
     fu=fu+cstr_element(5)*bd_curr(3,ni)
  elseif (i == 2) then
     fu=fu+cstr_element(6)*bd_curr(1,ni)
     fu=fu+cstr_element(2)*bd_curr(2,ni)
     fu=fu+cstr_element(4)*bd_curr(3,ni)
  elseif (i == 3) then
     fu=fu+cstr_element(5)*bd_curr(1,ni)
     fu=fu+cstr_element(4)*bd_curr(2,ni)
     fu=fu+cstr_element(3)*bd_curr(3,ni)
  endif

  return
end subroutine r_scalfu_curr

!     calculation of kpp
subroutine r_scalkpp(fkpp,ocpp,k,m)
  use r_common
  use solid_variables, only : nsd_solid
  implicit none

  real(8) :: fkpp,ocpp
  integer :: k,m

  fkpp = ocpp*hp(k)*hp(m)
  return
end subroutine r_scalkpp

!     calculation of kup
subroutine r_scalkup(fkup,ocup,i,k,ni)
  use r_common
  use solid_variables, only : nsd_solid
  implicit none
  
  real(8) :: fkup,ocup(2*nsd_solid)
  integer :: i,k,ni
  integer :: m

  fkup=0.0d0
  do m=1,2*nsd_solid
     fkup=fkup+ocup(m)*dge(m,i,ni)*hp(k)
  enddo
  return
end subroutine r_scalkup


