subroutine s_dam(ndof,d1,d2,vel,M,KT,fdam)
implicit none
integer :: ndof
real(8) :: d1,d2
real(8) :: vel(ndof),fdam(ndof)
real(8) :: M(ndof,ndof),KT(ndof,ndof),C(ndof,ndof) 
 C=d1*M+d2*KT
 fdam=matmul(vel,C)
end subroutine s_dam