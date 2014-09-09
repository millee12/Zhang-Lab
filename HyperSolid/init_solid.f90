module init_solid
implicit none
save
!total: dimensions, nodes, nodes/element, elements,
!gauss points,time,essential BCS
integer :: nsd=2,nn=21,nen=4,ne=12,ngp=4,tend=1000,ng
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
real(8),allocatable :: fint(:),kt(:,:),sel(:,:)
real(8) :: rc1,rc2,kappa
!=====================
!For the kinematic forces:
!kinematic force,kinematic stiffness
real(8),allocatable :: fkin(:),km(:,:),Mass(:,:),mtmp(:,:)
real(8) :: rho=1.0d0
!=====================
!For the damping forces:
!kinematic force,kinematic stiffness
real(8),allocatable :: fdam(:)
real(8) :: d1,d2
!=====================
!Solver
real(8) :: res
!lapack solver feedback
real(8), allocatable :: IPIV(:)
!displacement and velocity guesses
real(8), allocatable :: dtil(:),vtil(:),bf(:)
!lagrange multiplier, displacement, velocity, and acceleration updates
real(8), allocatable :: mlagnew(:),dnew(:),vnew(:),anew(:)
!previous displacement and velocity
real(8), allocatable :: uold(:),mlagold(:)
!lapack solver feedback, number of unknowns in Ax=b
integer :: INFO,nee
!=====================
end module init_solid