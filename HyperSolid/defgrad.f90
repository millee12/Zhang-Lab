subroutine defgrad(el,ne,ndof,eldof,nen,nsd_solid,xref,dis,nx,ien,F)
implicit none 
integer :: nsd_solid,ien(eldof,ne),dof,nen,ndof,eldof,ne,el
real(8) :: xref(ndof),dis(ndof),nx(nen,nsd_solid),F(nsd_solid,nsd_solid),detF 
integer :: a,p
F(:,:)=0.0d0
do a=1,nen
   dof=nsd_solid*(a-1)+1
   p=ien(dof,el)
   F(1,1)=F(1,1)+(xref(p)+dis(p))*nx(a,1)
   F(1,2)=F(1,2)+(xref(p)+dis(p))*nx(a,2)
   dof=nsd_solid*(a-1)+2
   p=ien(dof,el)
   F(2,1)=F(2,1)+(xref(p)+dis(p))*nx(a,1)
   F(2,2)=F(2,2)+(xref(p)+dis(p))*nx(a,2)
enddo
detF=F(1,1)*F(2,2)-F(1,2)*F(2,2)
!Check for element inversion
if (detF < 0.0d0) then
write(*,*) 'error: negative jacobian'
stop
endif
end subroutine