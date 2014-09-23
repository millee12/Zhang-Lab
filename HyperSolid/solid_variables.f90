module solid_variables
  implicit none
  save
  integer :: nsd_solid
  integer :: nn_solid
  integer :: nen_solid
  integer :: ndof_solid
  integer :: eldof_solid
  integer :: ne_solid
  integer,parameter :: nquad_solid=4
  integer :: neq_solid
  real(8),dimension(1:3) :: solid_scale  !...scale the size of the structure
  integer,parameter :: nsurface = 4
  integer :: nn_solid_1,ne_solid_1
  integer :: n_solid !...number of solids ( x times the solid, which is read in coortable
  integer :: iquad_solid  !...quadratur type, see "quadxdxn.f"

  integer,parameter :: ndfpad_solid=5,nsdpad_solid=3,nenpad_solid=8,nquadpad_solid=8
  real(8) xq_solid(nsdpad_solid,nquadpad_solid),wq_solid(nquadpad_solid)
  integer :: nsol_ebc !...number of nodes with essential BC (displacement)
  real(8),dimension(:,:),allocatable :: solid_ess_BC
  real(8),dimension(:,:),allocatable :: shift
  real(8),allocatable :: mirror(:,:)
  real(8),parameter :: beta=0.25d0
  real(8),parameter  :: gamma=0.5d0
  real(8),parameter  :: tol=1.0d-6

  integer :: nep1,nep2 !number of elements for each solid parts
  integer ne_sbc_1 ! number of edges on solid interface for parameter read in
  integer nn_sbc_1 ! number of nodes on solid interface for parameter read in
  !!!!!Remove!!!!!!
  integer,parameter :: ne_sbc=12 ! number of edges on solid interface
  integer,parameter ::  nn_sbc=14 ! number of nodes on solid interface
  integer node_sfcon ! number of solid nodes on fluid-solid connected boundary
  integer node_sfcon_1

  real(8),dimension(2),parameter :: damp_solid=[0.0d0, 0.0d0] ! solid damping coefficients
  integer :: parallel_solid

 real(8),parameter :: rho=1.0d0 
!MR constants
 real(8),parameter  :: rc1=8.6207d3
 real(8),parameter  :: rc2=0.0d0
 real(8),parameter  :: kappa=1.6667d8

real(8),allocatable :: detjac(:)
real(8),allocatable :: xref_solid(:)
real(8), allocatable :: nx(:,:,:,:)
real(8), allocatable :: shp(:,:,:)
integer, allocatable :: lm_solid(:,:)
integer, allocatable :: ien_solid(:,:)
integer, allocatable :: ien_sbc(:,:)

real(8), allocatable :: IPIV(:)
integer :: INFO
real(8),allocatable :: gx(:)
integer,allocatable :: gdof(:)
real(8),allocatable :: Mass(:,:)
real(8),allocatable :: km(:,:)
real(8),allocatable :: xyz_solid(:,:)

end module solid_variables
