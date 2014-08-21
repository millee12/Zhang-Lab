subroutine s_int(nsd,nnode,nen,nel,ngp,ndof,ng,rc1,rc2,kappa,xref,dg,ien,fint,kt)
!For the internal forces (calculated using Mooney-Rivlin)
integer :: nsd,nnode,nen,nel,ngp,ndof,ng
real :: fint(8),sint(8)
real :: F(2,2),pk2(3),rc1,rc2,kappa,CMR(3,3)
real :: xref(8),dg(8)
real, dimension(1) :: detjac
real,dimension(4,2,4,1) :: nx
real :: fpen(8)
integer :: i,j,gp,a,b,p,q,dof1,dof2,gnode(3),el,ien(8,1)
!For the B matrices
real :: l(2,2)
real :: BL(3,8),BNL(4,8)
real :: kt(8,8),kl(8,8),knl(8,8),ssuper(4,4)
!Solve for internal nodal forces and tangent stiffness matrix
fint(:)=0.0
do el=1,nel
	sint(:)=0.0
	do gp=1,ngp
		!Calculate Deformation Gradient
		call defgrad(nen,nsd,xref,dg,nx(:,:,gp,el),ien,F)
		! Calculate stress and stiffness using the Mooney-Rivlin material model
		call mooney(F,rc1,rc2,kappa,CMR,pk2,ssuper)
		!Make strain-displacement matrices
		call getB(nx,dg,BL,BNL)
		! Calculate stress gradient at gauss points
		sint(:)=matmul(transpose(BL),pk2)*detjac(el)
		kl(:,:)=matmul(transpose(BL),matmul(CMR,BL))*detjac(el)
        knl(:,:)=matmul(transpose(BNL),matmul(ssuper,BNL))*detjac(el)
		do dof1=1,ndof
			p=ien(dof1,el)
			fint(p)=fint(p)+sint(dof1)
			do dof2=1,ndof
					q=ien(dof2,el)
					kt(p,q)=kt(p,q)+kl(dof1,dof2)+knl(dof1,dof2)
			enddo
		enddo
	enddo
enddo

end subroutine