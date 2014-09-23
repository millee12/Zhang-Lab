subroutine s_pen(dis,mlag,rpen,fpen,ka)
use solid_variables
implicit none 
integer :: i,p
real(8) :: fpen(ndof_solid),mlag(nsol_ebc),rpen(nsol_ebc),IQ(ndof_solid,nsol_ebc),dis(ndof_solid)
real(8) :: ka(neq_solid,neq_solid)
IQ(:,:)=0.0d0
rpen(:)=0.0d0
do i=1,nsol_ebc
	p=gdof(i)
	IQ(p,i)=2.0d0*kappa
	rpen(i)=2*kappa*(gx(i)-dis(p))
enddo
fpen(:)=matmul(IQ,mlag)
ka(ndof_solid+1:neq_solid,1:ndof_solid)=transpose(IQ(1:ndof_solid,1:nsol_ebc))
ka(1:ndof_solid,ndof_solid+1:neq_solid)=IQ(1:ndof_solid,1:nsol_ebc)
end subroutine s_pen