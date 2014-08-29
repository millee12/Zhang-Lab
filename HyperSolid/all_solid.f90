implicit none
!total: dimensions, nodes, nodes/element, elements,
!gauss points,time,essential BCS
integer :: nsd=2,nn=9,nen=4,ne=4,ngp=4,tend=250,ng=3
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
real(8),allocatable :: fkin(:),km(:,:),Mass(:,:),mtmp(:,:)
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

allocate(xref(ndof))
allocate(dis(ndof))
allocate(rf(nee))
allocate(ien(eldof,ne))
allocate(shp(0:nsd,nen,ngp))
allocate(nx(nen,nsd,ngp,ne))
allocate(detjac(ne))
allocate(fext(ndof))
allocate(fpen(ndof))
allocate(mlag(ng))
allocate(rpen(ng))
allocate(tmpgdof(ndof,ndof))
allocate(ka(nee,nee))
allocate(fint(ndof))
allocate(kt(ndof,ndof))
allocate(IPIV(nee))
allocate(fkin(ndof))
allocate(km(ndof,ndof))
allocate(Mass(ndof,ndof))
allocate(mtmp(ndof,ndof))
allocate(vel(ndof))
allocate(acc(ndof))
allocate(mlagnew(ng))
allocate(dnew(ndof))
allocate(vnew(ndof))
allocate(anew(ndof))
allocate(dtil(ndof))
allocate(vtil(ndof))
allocate(xyz(nn,nsd))
allocate(solid_con(ne,nen))
allocate(mtype(ne))