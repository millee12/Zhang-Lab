subroutine s_pen(ndof,nee,ng,dis,gdof,gx,kappa,mlag,rpen,fpen,ka)
implicit none 
integer :: gdof(ng),ng,i,ndof,nee,p
real(8) :: fpen(ndof),gx(ng),mlag(ng),rpen(ng),IQ(ndof,ng),dis(ndof),kappa
real(8) :: ka(nee,nee)
IQ(:,:)=0.0d0
rpen(:)=0.0d0
do i=1,ng
	p=gdof(i)
	IQ(p,i)=kappa
	rpen(i)=kappa*(gx(i)-dis(p))
enddo
fpen(:)=matmul(IQ,mlag)
ka(ndof+1:nee,1:ndof)=transpose(IQ(1:ndof,1:ng))
ka(1:ndof,ndof+1:nee)=IQ(1:ndof,1:ng)
end subroutine s_pen