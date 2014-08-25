module mesh_convert_variables
  implicit none
  save
  integer,parameter :: DimOfElem(4) = (/ 2, 2, 3, 3 /)
  integer,parameter :: NdPerElem(4) = (/ 3, 4, 4, 8 /)
  integer,parameter :: EgPerElem(4) = (/ 3, 4, 4, 6 /)
  integer,parameter :: NdPerEdge(4) = (/ 2, 2, 3, 4 /)
  integer :: fluid_mesh_type, solid_mesh_type 
              ! 0 = All required files already generated
              ! 1 = Abaqus file called fluid.inp and solid.inp
  integer :: elem_type, elem_type_solid
  ! Abaqus
  integer,allocatable :: nefboundary(:)
  integer :: nrng_solid
  integer,allocatable :: bc_solid(:,:)
end module