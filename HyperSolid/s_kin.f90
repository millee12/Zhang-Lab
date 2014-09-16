subroutine s_kin(neq_solid,ndof,acc,Mass,km,fkin,ka)
implicit none
integer :: ndof,neq_solid,i,j
real(8),intent(in) :: Mass(ndof,ndof)
real(8) :: fkin(ndof),acc(ndof),km(ndof,ndof),ka(neq_solid,neq_solid)
ka(1:ndof,1:ndof)=ka(1:ndof,1:ndof)+km(:,:)
fkin(:)=0.0d0
do i=1,ndof
	do j=1,ndof
		fkin(i)=fkin(i)+Mass(i,j)*acc(j)
	enddo
enddo
end subroutine s_kin