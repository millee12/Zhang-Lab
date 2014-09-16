!Main File to test Solid Solver Solid Solver
program main 
use init_solid
use solid_variables
use mesh_convert_variables
use meshgen_solid
nsd_solid=2
ndof_solid=nn_solid*nsd_solid
eldof=nen*nsd_solid

!allocate variables
include 'all_solid.FOR'
nrng_solid=1
allocate(bc_solid(nrng_solid,2))
bc_solid(1,2)=10110
!bc_solid(2,2)=10000
bf(:)=[0.0d0,-5.0d0]
call read_abaqus_solid
call readien_solid(solid_con,ne,nen,mtype)
do el=1,ne
	do a=1,nen
		if (a .gt. 2) then
			b=a-2
		else
			b=a+2
		endif
		do i=1,nsd_solid
			p=nsd_solid*(solid_con(el,a)-1)+i
			q=nsd_solid*(b-1)+i
			ien(q,el)=p
		enddo
	enddo
enddo
call readx_solid(xyz,nn,nsd_solid)
do a=1,nn
	p=nsd_solid*(a-1)+1
	q=nsd_solid*(a-1)+2
	xref(p)=xyz(a,1)
	xref(q)=xyz(a,2)
enddo
open(unit = 12, file = 'residual.out')
open(unit=10,file='dg.out')
open(unit=4,file='d.out')
!Initial Conditions
dis(:)=0.0d0
vel(:)=0.0d0
fint(:)=0.0d0
!specify essential boundaries
call s_ess(nsd_solid,nen,ne,ndof,eldof,ien,ng,tmpgdof,fext)
n_solid_ess_BC=ng
neq_solid=ndof+ng
allocate(rf(neq_solid))
allocate(IPIV(neq_solid))
allocate(ka(neq_solid,neq_solid))
call s_gdof(ndof,ng,tmpgdof,gdof,gx,mlag)

!write(*,*) gdof
!stop
!MR constants
rc1=8.6207d3
rc2=0.0d0
kappa=1.6667d5
!Damping Constants
damp_solid(1)=0.25d0
damp_solid(2)=3.65d-4
!save shape functions and derivatives
call svshp(ndof,eldof,nen,nsd_solid,ngp,ne,ien,xref,nx,detjac,shp)
!save mass matrix and kinematic stiffness and add body forces
call getM(nsd_solid,dt,beta,gamma,rho,bf,ndof,ngp,eldof,nen,ne,detjac,ien,shp,Mass,km,fext)
!initial acceleration
acc=fext-fint
mtmp=Mass
call dgesv(ndof, 1, mtmp, ndof, IPIV, acc, ndof, INFO )
!write initial configuration
t=0
call paraout(ne,nen,t,nsd_solid,nn,ndof,dnew,xref,solid_con,sel,vel)
!Time Loop
do t=1,tend
	write(*,*) 't=', t
	write(*,*) ng
    call solve_solid(fext,n_solid_ess_BC,dis,vel,acc,dt,km,Mass)
	call paraout(ne,nen,t,nsd_solid,nn,ndof,dnew,xref,solid_con,sel,vel)
enddo	
end program main
