subroutine defgrad(nen,nsd,xref,dis,nx,ien,F)
implicit none 
integer :: nsd,ien(8,1),dof,nen
real(8) :: xref(8),dis(8),nx(nen,nsd),F(nsd,nsd),detF 
integer :: a,p
F(:,:)=0.0d0
do a=1,nen
   dof=nsd*(a-1)+1
   p=ien(dof,1)
   F(1,1)=F(1,1)+(xref(p)+dis(p))*nx(a,1)
   F(1,2)=F(1,2)+(xref(p)+dis(p))*nx(a,2)
   dof=nsd*(a-1)+2
   p=ien(dof,1)
   F(2,1)=F(2,1)+(xref(p)+dis(p))*nx(a,1)
   F(2,2)=F(2,2)+(xref(p)+dis(p))*nx(a,2)
enddo
detF=F(1,1)*F(2,2)-F(1,2)*F(2,2)
!Check for element inversion
if (detF < 0.0d0) then
write(*,*) 'warning: negative jacobian'
endif
end subroutine