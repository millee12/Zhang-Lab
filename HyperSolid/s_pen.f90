subroutine s_pen(ndof,ng,dg,gdof,gx,kt,kappa,mlag,rpen,fpen,ka)
real :: fpen(ndof),gx(ng),mlag(ng),rpen(ng),IQ(ndof,ng),dg(ndof),kappa
real :: kt(ndof,ndof),ka(ndof+ng,ndof+ng)
integer :: gdof(ng),ng,i
IQ(:,:)=0.0
do i=1,ng
	IQ(gdof(i),i)=kappa
enddo
rpen=kappa*(gx(:)-dg(gdof))
fpen(:)=matmul(IQ,mlag)
ka(:,:)=0.0
ka(1:ndof,1:ndof)=kt(1:ndof,1:ndof)
ka(ndof+1:ndof+ng,1:ndof)=transpose(IQ)
ka(1:ndof,ndof+1:ndof+ng)=IQ
end subroutine s_pen