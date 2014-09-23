subroutine defgrad(el,dis,hx,F)
use solid_variables
implicit none 
integer :: dof,el
real(8) :: dis(ndof_solid),hx(nen_solid,nsd_solid),F(nsd_solid,nsd_solid),detF 
integer :: a,p
F(:,:)=0.0d0
do a=1,nen_solid
   dof=nsd_solid*(a-1)+1
   p=lm_solid(dof,el)
   F(1,1)=F(1,1)+(xref(p)+dis(p))*hx(a,1)
   F(1,2)=F(1,2)+(xref(p)+dis(p))*hx(a,2)
   dof=nsd_solid*(a-1)+2
   p=lm_solid(dof,el)
   F(2,1)=F(2,1)+(xref(p)+dis(p))*hx(a,1)
   F(2,2)=F(2,2)+(xref(p)+dis(p))*hx(a,2)
enddo
detF=F(1,1)*F(2,2)-F(1,2)*F(2,2)
!Check for element inversion
if (detF < 0.0d0) then
write(*,*) 'error: negative jacobian'
stop
endif
end subroutine