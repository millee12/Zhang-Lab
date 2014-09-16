subroutine s_dam(ndof,damp_solid,vel,M,KT,fdam)
implicit none
integer :: ndof
real(8) :: damp_solid
real(8) :: vel(ndof),fdam(ndof)
real(8) :: M(ndof,ndof),KT(ndof,ndof),C(ndof,ndof) 
 C=damp_solid(1)*M+damp_solid(2)*KT
 fdam=matmul(vel,C)
end subroutine s_dam