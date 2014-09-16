module solid_variables
  implicit none
  save
  real(8),dimension(1:3) :: solid_scale  !...scale the size of the structure
  integer,parameter :: nsurface = 4
  integer :: nn_solid  ,ne_solid
  integer :: nn_solid_1,ne_solid_1
  integer :: nsd_solid
  integer :: n_solid !...number of solids ( x times the solid, which is read in coortable
  integer :: nen_solid   !...number of nodes per element
  integer :: iquad_solid  !...quadratur type, see "quadxdxn.f"
  integer :: nquad_solid
  integer :: ndof_solid
  integer :: eldof
  integer :: neq_solid
  integer,parameter :: ndfpad_solid=5,nsdpad_solid=3,nenpad_solid=8,nquadpad_solid=8
  real(8) xq_solid(nsdpad_solid,nquadpad_solid),wq_solid(nquadpad_solid)
  integer :: nsol_ebc !...number of nodes with essential BC (displacement)
  real(8),dimension(:,:),allocatable :: solid_ess_BC
  real(8),dimension(:,:),allocatable :: shift
  real(8),allocatable :: mirror(:,:)

  integer :: nep1,nep2 !number of elements for each solid parts
  integer ne_sbc_1 ! number of edges on solid interface for parmeter read in
  integer nn_sbc_1 ! number of nodes on solid interface for parmeter read in
  integer ne_sbc ! number of edges on solid interface
  integer nn_sbc ! number of nodes on solid interface
  integer node_sfcon ! number of solid nodes on fluid-solid connected boundary
  integer node_sfcon_1

  real(8) damp_solid(2) ! solid damping coefficients

  integer :: parallel_solid

end module solid_variables
