subroutine s_int(dnew,fint,kt,ka,sel)
use solid_variables
implicit none
!For the internal forces (calculated using Mooney-Rivlin)
real(8) :: fint(ndof_solid)
real(8) :: sint(ndof_solid)
real(8) :: F(nsd_solid,nsd_solid)
real(8) :: pk2(3)
real(8) :: CMR(3,3)
real(8) :: sel(3,ne_solid)
real(8) :: dnew(ndof_solid)
real(8) :: fpen(ndof_solid)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el
real(8) :: ka(neq_solid,neq_solid)
real(8) :: kt(ndof_solid,ndof_solid)
real(8) :: kl(eldof_solid,eldof_solid)
real(8) :: knl(eldof_solid,eldof_solid)
real(8) :: ssuper(4,4)
!For the B matrices
real(8) :: l(nsd_solid,nsd_solid)
real(8) :: BL(3,eldof_solid)
real(8) :: BNL(4,eldof_solid)
real(8) :: hx(nen_solid,nsd_solid)
!Solve for internal nodal forces and tangent stiffness matrix
ka(:,:)=0.0d0
fint(:)=0.0d0
kt(:,:)=0.0d0
sel(:,:)=0.0d0
do el=1,ne_solid
	kl(:,:)=0.0d0
	knl(:,:)=0.0d0
	sint(:)=0.0d0
	do gp=1,nquad_solid
		do i=1,nen_solid
			do j=1,nsd_solid
				hx(i,j)=nx(i,j,gp,el)
			enddo
		enddo
		!Calculate Deformation Gradient
		call defgrad(el,dnew,hx,F)

		! Calculate stress and stiffness using the Mooney-Rivlin material model
		call mooney(F,CMR,pk2,ssuper)
		!Make strain-displacement matrices
		call getB(el,dnew,hx,BL,BNL)
		! Calculate stress gradient at gauss points
		sint(:)=sint(:)+matmul(transpose(BL),pk2)*detjac(el)	
		kl(:,:)=kl(:,:)+matmul(transpose(BL),matmul(CMR,BL))*detjac(el)
		knl(:,:)=knl(:,:)+matmul(transpose(BNL),matmul(ssuper,BNL))*detjac(el)
       	sel(:,el)=	sel(:,el)+(1.0d0/nquad_solid)*pk2(:)
	enddo
	do dof1=1,eldof_solid
<<<<<<< HEAD
		p=ien(dof1,el)
		fint(p)=fint(p)+sint(dof1)
		do dof2=1,eldof_solid
			q=ien(dof2,el)
=======
		p=ien_solid(dof1,el)
		fint(p)=fint(p)+sint(dof1)
		do dof2=1,eldof_solid
			q=ien_solid(dof2,el)
>>>>>>> 9ca02509ab857d8dbc3512aa29c539835bd53dbc
			kt(p,q)=kt(p,q)+kl(dof1,dof2)+knl(dof1,dof2)
		enddo
	enddo

enddo
ka(1:ndof_solid,1:ndof_solid)=kt(1:ndof_solid,1:ndof_solid)
end subroutine s_int
