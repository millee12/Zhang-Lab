subroutine s_dam(vel,KT,fdam)
use solid_variables
implicit none
real(8) :: vel(ndof_solid),fdam(ndof_solid)
real(8) :: KT(ndof_solid,ndof_solid),C(ndof_solid,ndof_solid) 
 C=damp_solid(1)*Mass+damp_solid(2)*KT
 fdam=matmul(vel,C)
end subroutine s_dam