!Main File to test Solid Solver Solid Solver
program main
!total: dimensions, nodes, nodes/element, elements,
!gauss points, degrees of freedom,time,essential BCS
integer :: nsd=2,nnode=4,nen=4,nel=1,ngp=4,ndof,tend=50,ng=2
!reference config,displacements,residual AND delta d
real,allocatable :: xref(:),dg(:),rf(:)
!element dof to global dof
integer,allocatable :: ien(:,:)
!counters
integer :: i,j,gp,a,b,p,q,dof1,dof2,el,t
!shape function derivative dN/dX, 
!mapping determinant (not detF), eternal force
real, allocatable :: nx(:,:,:,:),detjac(:),fext(:)
!For lagrange multiplier force:
!penalty force, prescribed displacements,lagrange multiplier, 
!essential BC residual, augmented stiffness matrix
real,allocatable :: fpen(:),gx(:),mlag(:),rpen(:),ka(:,:)
integer,allocatable :: gdof(:)
!For the internal forces (calculated using Mooney-Rivlin):
!internal force,tangent stiffness
real,allocatable :: fint(:),kt(:,:)
real :: rc1,rc2,kappa
!Solver
real :: res
!lapack solver feedback
real, allocatable :: IPIV(:)
!lapack solver feedback, number of unknowns in Ax=b
integer :: INFO,nee
ndof=nnode*nsd
nee=ndof+ng
!allocate variables
allocate(xref(ndof))
allocate(dg(ndof))
allocate(rf(nee))
allocate(ien(ndof,nel))
allocate(nx(nen,nsd,ngp,nel))
allocate(detjac(nel))
allocate(fext(ndof))
allocate(fpen(ndof))
allocate(gx(ng))
allocate(mlag(ng))
allocate(rpen(ng))
allocate(ka(nee,nee))
allocate(gdof(ng))
allocate(fint(ndof))
allocate(kt(ndof,ndof))
allocate(IPIV(nee))
!!!!!!!!!!!!!!!!!!!!
! Mesh Information
ien(:,1)=[7,8,5,6,1,2,3,4]
xref(:)=[0.0,0.0,1.0,0.0,0.0,1.0,1.0,1.0]
fext(:)=0.0
fext(3)=5000.0
fext(7)=5000.0
dg(:)=0.0
!specify essential boundaries
gdof(:)=[1,5]
gx(:)=0.0
mlag(:)=0.0
!MR constants
rc1=8.6207d3
rc2=0.0
kappa=1.6667d5

call svshp(ndof,nen,nsd,ngp,nel,ien,xref,nx,detjac)

do t=1,tend
	!internal forces and tangent stiffness matrix
	call s_int(nsd,nnode,nen,nel,ngp,ndof,ng,rc1,rc2,kappa,xref,dg,nx,detjac,ien,fint,kt)
	call s_pen(ndof,ng,dg,gdof,gx,kt,kappa,mlag,rpen,fpen,ka)
	rf(:)=[fext(:)-fint(:)-fpen(:),rpen]
	call sgesv(10, 1, ka, 10, IPIV, rf(1:10), 10, INFO )
	! SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
	if (INFO .eq. 0) then
		write(*,*) 'stiffness matrix inverted'
		else
		write(*,*) 'warning: linear solver exited with errors'
	endif
	dg(:)=dg(:)+rf(1:ndof)
	mlag(:)=mlag(:)+rf(ndof+1:ndof+ng)
	res=sqrt(sum(rf*rf))
	write(*,*) 'res= ', res
enddo	
end program main