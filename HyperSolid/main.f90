!Main File to test Solid Solver Solid Solver
program main 
use solid_variables
use mesh_convert_variables
use meshgen_solid
implicit none
!total: dimensions, nodes, nodes/element, elements,
!gauss points,time,essential BCS
integer :: nsd=2,nn=9,nen=4,ne=4,ngp=4,tend=50,ng=3
!timestep
real(8) :: dt=1.0d-3,beta=0.25d0,gamma=0.5d0,tol=1.0d-6
!=====================
!total: global degrees of freedom,element degrees of freedom
integer :: ndof,eldof,file
!reference config,displacements,residual AND delta d
real(8),allocatable :: xref(:),dis(:),vel(:),acc(:),rf(:),xyz(:,:)
!element dof to global dof
integer,allocatable :: ien(:,:),solid_con(:,:),mtype(:)
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
integer,allocatable :: gdof(:),tmpgdof(:,:)
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
ndof=nn*nsd
eldof=nen*nsd
nee=ndof+ng
!allocate variables
allocate(xref(ndof))
allocate(dis(ndof))
allocate(rf(nee))
allocate(ien(eldof,ne))
allocate(shp(0:nsd,nen,ngp))
allocate(nx(nen,nsd,ngp,ne))
allocate(detjac(ne))
allocate(fext(ndof))
allocate(fpen(ndof))
allocate(gx(ng))
allocate(mlag(ng))
allocate(rpen(ng))
allocate(tmpgdof(ndof,ndof))
allocate(ka(nee,nee))
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
!new
allocate(xyz(nn,nsd))
allocate(solid_con(ne,nen))
allocate(mtype(ne))
!!!!!!!!!!!!!!!!!!!!
!write(*,*) 'input nrng'
!read(*,*) nrng_solid
nrng_solid=2
allocate(bc_solid(nrng_solid,2))
!do i=1,nrng_solid
!	write(*,*) 'input', i ,'BC'
!	read(*,*) bc_solid(i,2)
!enddo
bc_solid(1,2)=10100
bc_solid(2,2)=10000
call read_abaqus_solid
call readien_solid(solid_con,ne,nen,mtype)
do el=1,ne
	!write(*,*) solid_con(el,:)
	do a=1,nen
		if (a .gt. 2) then
			b=a-2
		else
			b=a+2
		endif
		do i=1,nsd
			p=nsd*(solid_con(el,a)-1)+i
			q=nsd*(b-1)+i
			ien(q,el)=p
		enddo
	enddo
enddo
call readx_solid(xyz,nn,nsd)
do a=1,nn
	p=nsd*(a-1)+1
	q=nsd*(a-1)+2
	xref(p)=xyz(a,1)
	xref(q)=xyz(a,2)
enddo
open(unit = 1, file = 'residual.out')
open(unit = 3, file = 'mxyz.out')
! Mesh Information
!fext(:)=0.0d0
!fext(3)=5.0d2
!fext(7)=5.0d2
!Initial Conditions
dis(:)=0.0d0
vel(:)=0.0d0
fint(:)=0.0d0
!specify essential boundaries
call s_ess(nsd,nen,ne,ndof,eldof,ien,ng,tmpgdof,fext)
allocate(gdof(ng))
call s_gdof(ndof,ng,tmpgdof,gdof,gx,mlag)
!MR constants
rc1=8.6207d3
rc2=0.0d0
kappa=1.6667d5
!save shape functions and derivatives
call svshp(ndof,eldof,nen,nsd,ngp,ne,ien,xref,nx,detjac,shp)
!save mass matrix and kinematic stiffness
call getM(nsd,dt,beta,gamma,rho,ndof,ngp,eldof,nen,ne,detjac,ien,shp,M,km)
!initial acceleration
acc=fext-fint
call dgesv(ndof, 1, M, ndof, IPIV, acc, ndof, INFO )
!Time Loop
do t=1,tend
	write(*,*) 't=', t
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
		call s_int(nee,nsd,nn,nen,ne,ngp,ndof,eldof,ng,rc1,rc2,kappa,xref,dnew,nx,detjac,ien,fint,ka)
		write(*,*) fint
		call s_kin(nee,ndof,anew,M,km,fkin,ka)
		call s_pen(ndof,nee,ng,dnew,gdof,gx,kappa,mlagnew,rpen,fpen,ka)
		rf=[fext-fint-fpen-fkin,rpen]
		call dgesv(nee, 1, ka, nee, IPIV, rf, nee, INFO )
		dnew=dnew+rf(1:ndof)
		mlagnew=mlagnew+rf(ndof+1:ndof+ng)
		res=sqrt(sum(rf*rf))
		write(1,*) res
		if (w .ge. 1) then
		stop
		endif
	enddo

	write(1,*) '-------------'
	dis=dnew
	vel=vnew
	acc=anew
	do i=1,nn
	write(3,*) dnew(nsd*(i-1)+1)+xref(nsd*(i-1)+1), dnew(nsd*(i-1)+2)+xref(nsd*(i-1)+2)
	enddo
	write(3,*) '-------------'
enddo	

!dallocate variables
!deallocate(xref)
!deallocate(dis)
!deallocate(rf)
!deallocate(ien)
!deallocate(nx)
!deallocate(detjac)
!deallocate(fext)
!deallocate(fpen)
!deallocate(gx)
!deallocate(mlag)
!deallocate(rpen)
!deallocate(ka)
!deallocate(gdof)
!deallocate(fint)
!deallocate(kt)
!deallocate(IPIV)
end program main