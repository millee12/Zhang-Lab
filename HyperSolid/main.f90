!Main File to test Solid Solver Solid Solver
program main
implicit none 
!total: dimensions, nodes, nodes/element, elements,
!gauss points,time,essential BCS
integer :: nsd=2,nnode=4,nen=4,nel=1,ngp=4,tend=50,ng=2
!timestep
real(8) :: dt=1.0d-3,beta=0.25d0,gamma=0.5d0,tol=1.0d-6
!=====================
!total: global degrees of freedom,element degrees of freedom
integer :: ndof,eldof
!reference config,displacements,residual AND delta d
real(8),allocatable :: xref(:),dis(:),vel(:),acc(:),rf(:)
!element dof to global dof
integer,allocatable :: ien(:,:)
!=====================
!counters
integer :: i,j,gp,a,b,p,q,el,t,w
!=====================
!shape function derivative dN/dX, 
!mapping determinant (not detF), eternal force
real(8), allocatable :: shp(:,:,:),nx(:,:,:,:),detjac(:),fext(:)
!=====================
!For lagrange multiplier force:
!penalty force, prescribed displacements,lagrange multiplier, 
!essential BC residual, augmented stiffness matrix
real(8),allocatable :: fpen(:),gx(:),mlag(:),rpen(:),ka(:,:)
integer,allocatable :: gdof(:)
!=====================
!For the internal forces (calculated using Mooney-Rivlin):
!internal force,tangent stiffness
real(8),allocatable :: fint(:),kt(:,:)
real(8) :: rc1,rc2,kappa
!=====================
!For the kinematic forces:
!kinematic force,kinematic stiffness
real(8),allocatable :: fkin(:),km(:,:),M(:,:)
real(8) :: rho=1.0d0
!=====================
!Solver
real(8) :: res
!lapack solver feedback
real(8), allocatable :: IPIV(:)
!displacement and velocity guesses
real(8), allocatable :: dtil(:),vtil(:)
!lagrange multiplier, displacement, velocity, and acceleration updates
real(8), allocatable :: mlagnew(:),dnew(:),vnew(:),anew(:)
!previous displacement and velocity
real(8), allocatable :: uold(:),mlagold(:)
!lapack solver feedback, number of unknowns in Ax=b
integer :: INFO,nee
!=====================
ndof=nnode*nsd
eldof=nen*nsd
nee=ndof+ng
!allocate variables
allocate(xref(ndof))
allocate(dis(ndof))
allocate(rf(nee))
allocate(ien(eldof,nel))
allocate(shp(0:nsd,nen,ngp))
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
allocate(fkin(ndof))
allocate(km(ndof,ndof))
allocate(M(ndof,ndof))
allocate(vel(ndof))
allocate(acc(ndof))
allocate(mlagnew(ng))
allocate(dnew(ndof))
allocate(vnew(ndof))
allocate(anew(ndof))
allocate(dtil(ndof))
allocate(vtil(ndof))
!!!!!!!!!!!!!!!!!!!!

open(unit = 1, file = 'residual.txt')
open(unit = 2, file = 'ka.txt')
open(unit = 3, file = 'mxyz.txt')
! Mesh Information
ien(:,1)=[7,8,5,6,1,2,3,4]
xref(:)=[0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,1.0d0,1.0d0,1.0d0]
fext(:)=0.0d0
fext(3)=5.0d4
fext(7)=5.0d4
dis(:)=0.0d0
vel(:)=0.0d0
!specify essential boundaries
gdof(:)=[1,5]
gx(:)=0.0d0
mlag(:)=0.0d0
!MR constants
rc1=8.6207d3
rc2=0.0d0
kappa=1.6667d5
!save shape functions and derivatives
call svshp(ndof,eldof,nen,nsd,ngp,nel,ien,xref,nx,detjac,shp)
!save mass matrix and kinematic stiffness
call getM(nsd,dt,beta,gamma,rho,ndof,ngp,eldof,nen,nel,detjac,ien,shp,M,km)

!Time Loop
do t=1,tend
    ! Predict
    dtil=dis+dt*vel+(dt**2/2.0d0)*(1.0d0-2.0d0*beta)*acc
    vtil=vel+(1.0d0-gamma)*dt*acc
    ! Correct
    w=0;
    dnew=dtil
    vnew=vtil
    mlagnew=mlag
	res=huge(1.d0)
	do while ((res .ge. tol) .and. (w .le. 100))
	 	w=w+1
        anew=(1/(beta*dt**2))*(dnew-dtil)
        vnew=vtil+gamma*dt*anew
		!internal forces and tangent stiffness matrix
		call s_int(nee,nsd,nnode,nen,nel,ngp,ndof,eldof,ng,rc1,rc2,kappa,xref,dnew,nx,detjac,ien,fint,ka)
		call s_kin(nee,ndof,anew,M,km,fkin,ka)
		call s_pen(ndof,nee,ng,dnew,gdof,gx,kappa,mlagnew,rpen,fpen,ka)
		!do i=1,ndof
		!write(*,*) km(i,:)
		!enddo
		rf=[fext-fint-fpen-fkin,rpen]
		write(*,*) rf
		call dgesv(nee, 1, ka, nee, IPIV, rf, nee, INFO )
		! SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
		dnew=dnew+rf(1:ndof)
		mlagnew=mlagnew+rf(ndof+1:ndof+ng)
		res=sqrt(sum(rf*rf))
		if (INFO .eq. 0) then
		write(*,*) 'w= ', w
		else
		write(*,*) 'warning: linear solver exited with errors'
		endif
		write(1,*) res
	enddo
	do i=1,nnode
	write(3,*) dnew(nsd*(i-1)+1)+xref(nsd*(i-1)+1), dnew(nsd*(i-1)+2)+xref(nsd*(i-1)+2)
	enddo
	stop
enddo	
!dallocate variables
deallocate(xref)
deallocate(dis)
deallocate(rf)
deallocate(ien)
deallocate(nx)
deallocate(detjac)
deallocate(fext)
deallocate(fpen)
deallocate(gx)
deallocate(mlag)
deallocate(rpen)
deallocate(ka)
deallocate(gdof)
deallocate(fint)
deallocate(kt)
deallocate(IPIV)
end program main