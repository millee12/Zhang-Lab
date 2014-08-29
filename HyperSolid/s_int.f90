subroutine s_int(nee,nsd,nnode,nen,ne,ngp,ndof,eldof,ng,rc1,rc2,kappa,xref,dis,nx,detjac,ien,fint,ka,sel)
implicit none
!For the internal forces (calculated using Mooney-Rivlin)
integer :: nsd,nnode,nen,ne,ngp,ndof,eldof,ng,nee
real(8) :: fint(ndof),sint(ndof)
real(8) :: F(nsd,nsd),pk2(3),rc1,rc2,kappa,CMR(3,3),sel(3,ne)
real(8) :: xref(ndof),dis(ndof)
real(8) :: detjac(ne)
real(8) :: nx(nen,nsd,ngp,ne)
real(8) :: fpen(ndof)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el,ien(eldof,ne)
real(8) :: ka(nee,nee),kt(ndof,ndof),kl(eldof,eldof),knl(eldof,eldof),ssuper(4,4)
!For the B matrices
real(8) :: l(nsd,nsd)
real(8) :: BL(3,eldof),BNL(4,eldof),dummy(3,eldof)
!Solve for internal nodal forces and tangent stiffness matrix
ka(:,:)=0.0d0
fint(:)=0.0d0
kt(:,:)=0.0d0
sel(:,:)=0.0d0
do el=1,ne
	kl(:,:)=0.0d0
	knl(:,:)=0.0d0
	sint(:)=0.0d0
	do gp=1,ngp
		!Calculate Deformation Gradient
		call defgrad(el,ne,ndof,eldof,nen,nsd,xref,dis,nx(:,:,gp,el),ien,F)
		! Calculate stress and stiffness using the Mooney-Rivlin material model
		call mooney(nsd,F,rc1,rc2,kappa,CMR,pk2,ssuper)
		!write(*,*) pk2
		!Make strain-displacement matrices
		call getB(el,ne,ndof,eldof,nen,ien,nsd,nx(:,:,gp,el),dis,BL,BNL)
		! Calculate stress gradient at gauss points
		sint(:)=sint(:)+matmul(transpose(BL),pk2)*detjac(el)
		kl(:,:)=kl(:,:)+matmul(transpose(BL),matmul(CMR,BL))*detjac(el)
        knl(:,:)=knl(:,:)+matmul(transpose(BNL),matmul(ssuper,BNL))*detjac(el)
       	sel(:,el)=	sel(:,el)+(1.0d0/ngp)*pk2(:)
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