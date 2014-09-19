!Main File to test Solid Solver Solid Solver
program main 
use solid_variables
use mesh_convert_variables
use meshgen_solid
implicit none
integer :: tend=1000
real(8), allocatable :: fext(:)
integer,allocatable :: solid_con(:,:),mtype(:)
real(8),allocatable :: mtmp(:,:)
integer :: i,j,gp,a,b,p,q,el,t,w
real(8),allocatable :: xyz(:,:)
real(8),allocatable :: sel(:,:)
real(8), allocatable ::bf(:)
real(8),allocatable :: dis(:),vel(:),acc(:),mlag(:)
! Hard-coded inputs
nrng_solid=1
allocate(bc_solid(nrng_solid,2))
bc_solid(1,2)=10110
!====read====
call read_abaqus_solid
ne_solid=nep1
nn_solid=nn_solid_1
ndof_solid=nn_solid*nsd_solid
eldof_solid=nen_solid*nsd_solid
!===========
!====allocate=====
include 'all_solid.FOR'
allocate(fext(ndof_solid))
fext(:)=0.0d0
allocate(bf(nsd_solid))
bf(:)=[0.0d0,-5.0d0]
allocate(solid_con(ne_solid,nen_solid))
allocate(mtype(ne_solid))
!==================
!====read====
call readien_solid(solid_con,ne_solid,nen_solid,mtype)
allocate(xyz(nn_solid,nsd_solid))
call readx_solid(xyz,nn_solid,nsd_solid)
!===========
!====convert====
allocate(ien_solid(eldof_solid,ne_solid))
open(unit=35,file='solidcon.out')
do i=1,nen_solid
	write(35,*) solid_con(i,:)
enddo 
ien_solid(:,:)=0
do el=1,ne_solid
	do a=1,nen_solid
		if (a .gt. 2) then
			b=a-2
		else
			b=a+2
		endif
		do i=1,nsd_solid
			p=nsd_solid*(solid_con(el,a)-1)+i
			q=nsd_solid*(b-1)+i
			ien_solid(q,el)=p
		enddo
	enddo
enddo
allocate(xref(ndof_solid))
do a=1,nn_solid
	p=nsd_solid*(a-1)+1
	q=nsd_solid*(a-1)+2
	xref(p)=xyz(a,1)
	xref(q)=xyz(a,2)
enddo
!===========
open(unit=34,file='ien_solid.out')
do i=1,eldof_solid
	write(34,*) ien_solid(i,:)
enddo 
!specify essential boundaries
call s_ess(fext)
gx(:)=0.0d0
!save shape functions and derivatives
call svshp
!save mass matrix and kinematic stiffness and add body forces
call getM(bf,fext)
!Initial Conditions
allocate(dis(ndof_solid))
allocate(vel(ndof_solid))
allocate(acc(ndof_solid))
allocate(mlag(nsol_ebc))
dis(:)=0.0d0
vel(:)=0.0d0
acc(:)=0.0d0
mlag(:)=0.0d0
!initial acceleration
acc=fext
allocate(mtmp(ndof_solid,ndof_solid))
mtmp(:,:)=Mass(:,:)
call dgesv(ndof_solid, 1, mtmp,ndof_solid, IPIV, acc, ndof_solid, INFO )
!write initial configuration
		open(unit=28,file='acc.out')
		write(28,*) '--initial--'
		do i=1,ndof_solid
		write(28,*) acc(i)
		enddo
t=0
call paraout(t,dis,vel,solid_con,sel)
!Time Loop
do t=1,tend
	write(*,*) 't=', t
    call solve_solid(fext,dis,vel,acc,mlag,sel)
	call paraout(t,dis,vel,solid_con,sel)
enddo	
end program main
