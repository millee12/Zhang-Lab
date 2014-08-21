subroutine defgrad(nen,nsd,xref,dg,nx,ien,F)
integer :: nsd,ien(8,1),dof,nen
real :: xref(8),dg(8),nx(nen,nsd),F(nsd,nsd) 
integer :: a,p
F(:,:)=0.0
do a=1,nen
   dof=nsd*(a-1)+1
   p=ien(dof,1)
   F(1,1)=F(1,1)+(xref(p)+dg(p))*nx(a,1)
   F(1,2)=F(1,2)+(xref(p)+dg(p))*nx(a,2)
   dof=nsd*(a-1)+2
   p=ien(dof,1)
   F(2,1)=F(2,1)+(xref(p)+dg(p))*nx(a,1)
   F(2,2)=F(2,2)+(xref(p)+dg(p))*nx(a,2)
enddo
end subroutine