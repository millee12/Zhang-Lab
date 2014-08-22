subroutine s_kin(nee,ndof,acc,M,km,fkin,ka)
implicit none
integer :: ndof,nee
real(8) :: M(ndof,ndof),fkin(ndof),acc(ndof),km(ndof,ndof),ka(nee,nee)
ka(1:ndof,1:ndof)=ka(1:ndof,1:ndof)+km(:,:)
fkin=matmul(M,acc)
end subroutine s_kin