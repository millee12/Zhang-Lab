subroutine s_int(nee,nsd,nnode,nen,nel,ngp,ndof,eldof,ng,rc1,rc2,kappa,xref,dg,nx,detjac,ien,fint,ka)
!For the internal forces (calculated using Mooney-Rivlin)
integer :: nsd,nnode,nen,nel,ngp,ndof,eldof,ng
real(8) :: fint(ndof),sint(ndof)
real(8) :: F(nsd,nsd),pk2(3),rc1,rc2,kappa,CMR(3,3)
real(8) :: xref(ndof),dg(ndof)
real(8) :: detjac(nel)
real(8) :: nx(nen,nsd,ngp,nel)
real(8) :: fpen(ndof)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el,ien(eldof,nel)
real(8) :: ka(nee,nee),kt(ndof,ndof),kl(eldof,eldof),knl(eldof,eldof),ssuper(4,4)
!For the B matrices
real(8) :: l(nsd,nsd)
real(8) :: BL(3,eldof),BNL(4,eldof),dummy(3,eldof)
!Solve for internal nodal forces and tangent stiffness matrix
ka(:,:)=0.0d0
fint(:)=0.0d0
kt(:,:)=0.0d0
do el=1,nel
	kl(:,:)=0.0d0
	knl(:,:)=0.0d0
	sint(:)=0.0d0
	do gp=1,ngp
		!Calculate Deformation Gradient
		call defgrad(nen,nsd,xref,dg,nx(:,:,gp,el),ien,F)
		! Calculate stress and stiffness using the Mooney-Rivlin material model
		call mooney(F,rc1,rc2,kappa,CMR,pk2,ssuper)
		!Make strain-displacement matrices
		call getB(ndof,eldof,nen,nsd,nx(:,:,gp,el),dg,BL,BNL)
		! Calculate stress gradient at gauss points
		sint(:)=sint(:)+matmul(transpose(BL),pk2)*detjac(el)
		kl(:,:)=kl(:,:)+matmul(transpose(BL),matmul(CMR,BL))*detjac(el)
        knl(:,:)=knl(:,:)+matmul(transpose(BNL),matmul(ssuper,BNL))*detjac(el)
	enddo
	do dof1=1,eldof
		p=ien(dof1,el)
		fint(p)=fint(p)+sint(dof1)
		do dof2=1,eldof
			q=ien(dof2,el)
			kt(p,q)=kt(p,q)+kl(dof1,dof2)+knl(dof1,dof2)
		enddo
	enddo
enddo
ka(1:ndof,1:ndof)=kt(1:ndof,1:ndof)
end subroutine s_int