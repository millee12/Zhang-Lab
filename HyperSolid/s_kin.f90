subroutine s_kin(acc,fkin,ka)
use solid_variables
implicit none
integer :: i,j
real(8) :: fkin(ndof_solid),acc(ndof_solid),ka(neq_solid,neq_solid)
ka(1:ndof_solid,1:ndof_solid)=ka(1:ndof_solid,1:ndof_solid)+km(:,:)
fkin(:)=0.0d0
do i=1,ndof_solid
	do j=1,ndof_solid
		fkin(i)=fkin(i)+Mass(i,j)*acc(j)
	enddo
enddo
end subroutine s_kin