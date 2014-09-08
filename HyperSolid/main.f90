!Main File to test Solid Solver Solid Solver
program main 
use init_solid
use solid_variables
use mesh_convert_variables
use meshgen_solid
ndof=nn*nsd
eldof=nen*nsd
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
open(unit = 12, file = 'residual.out')
open(unit=10,file='dg.out')
open(unit=4,file='d.out')
!Initial Conditions
dis(:)=0.0d0
vel(:)=0.0d0
fint(:)=0.0d0
!specify essential boundaries
call s_ess(nsd,nen,ne,ndof,eldof,ien,ng,tmpgdof,fext)
allocate(mlagnew(ng))
allocate(mlag(ng))
allocate(rpen(ng))
allocate(gdof(ng))
allocate(gx(ng))
nee=ndof+ng
allocate(rf(nee))
allocate(IPIV(nee))
allocate(ka(nee,nee))
call s_gdof(ndof,ng,tmpgdof,gdof,gx,mlag)
!write(*,*) gdof
!stop
!MR constants
rc1=8.6207d3
rc2=0.0d0
kappa=1.6667d5
!Damping Constants
d1=0.25d0
d2=3.65d-4
!save shape functions and derivatives
call svshp(ndof,eldof,nen,nsd,ngp,ne,ien,xref,nx,detjac,shp)
!save mass matrix and kinematic stiffness and add body forces
call getM(nsd,dt,beta,gamma,rho,bf,ndof,ngp,eldof,nen,ne,detjac,ien,shp,Mass,km,fext)
!initial acceleration
acc=fext-fint
mtmp=Mass
call dgesv(ndof, 1, mtmp, ndof, IPIV, acc, ndof, INFO )
!write initial configuration
t=0
call paraout(ne,nen,t,nsd,nn,ndof,dnew,xref,solid_con,sel,vel)
!Time Loop
do t=1,tend
	write(*,*) 't=', t
	write(*,*) ng
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
		ka(:,:)=0.0d0
	 	w=w+1
		!write(*,*) 'w= ',w
        anew=(1/(beta*dt**2))*(dnew-dtil)
        vnew=vtil+gamma*dt*anew
		!internal forces and tangent stiffness matrix
		call s_int(nee,nsd,nn,nen,ne,ngp,ndof,eldof,ng,rc1,rc2,kappa,xref,dnew,nx,detjac,ien,fint,kt,ka,sel)
		call s_dam(ndof,d1,d2,vnew,Mass,kt,fdam)
		call s_kin(nee,ndof,anew,Mass,km,fkin,ka)
		call s_pen(ndof,nee,ng,dnew,gdof,gx,kappa,mlagnew,rpen,fpen,ka)
		rf(1:ndof)=fext-fint-fpen-fkin-fdam
		rf(ndof+1:nee)=rpen
		res=sqrt(sum(rf**2))/ne
		write(12,*) res
		write(10,*) '----residual------'
		do i=1,nee
			write(10,*) rf(i)
		enddo
		call dgesv(nee, 1, ka, nee, IPIV, rf, nee, INFO )
		if (INFO .ne. 0) then
			write(*,*) 'error: unable to solve for incremental displacements'
			write(*,*) INFO
			stop
		endif
		dnew=dnew+rf(1:ndof)
		mlagnew=mlagnew+rf(ndof+1:ndof+ng)
	enddo
	if (w .gt. 1) then
		write(*,*) 'succesful exit after ', w, 'iterations'
	endif
	dis=dnew
	vel=vnew
	acc=anew
	mlag=mlagnew
	call paraout(ne,nen,t,nsd,nn,ndof,dnew,xref,solid_con,sel,vel)
enddo	
end program main
