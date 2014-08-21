subroutine s_pen(ndof,nee,ng,dg,gdof,gx,kappa,mlag,rpen,fpen,ka)
implicit none 
integer :: gdof(ng),ng,i,ndof,nee
real(8) :: fpen(ndof),gx(ng),mlag(ng),rpen(ng),IQ(ndof,ng),dg(ndof),kappa
real(8) :: ka(nee,nee)

IQ(:,:)=0.0d0
do i=1,ng
	IQ(gdof(i),i)=kappa
enddo
rpen=kappa*(gx(:)-dg(gdof))
fpen(:)=matmul(IQ,mlag)
ka(ndof+1:nee,1:ndof)=transpose(IQ)
ka(1:ndof,ndof+1:nee)=IQ
end subroutine s_pen