subroutine solve_solid(ien_sbc,solid_stress,solid_dis,solid_vel,solid_acc&
,mlag,solid_coor_curr,sel,solid_fint,solid_fkin)
use solid_variables
use run_variables, only: dt,its,ntsbout
use mpi_variables
implicit none
real(8) :: fint(ndof_solid),fkin(ndof_solid)
real(8) :: fdam(ndof_solid),fpen(ndof_solid),fext(ndof_solid)
real(8) :: CMR(3,3),sel(3,ne_solid)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el
real(8) :: ka(neq_solid,neq_solid),kt(ndof_solid,ndof_solid)
real(8) :: kl(eldof_solid,eldof_solid),knl(eldof_solid,eldof_solid),ssuper(4,4)
real(8) :: rf(neq_solid)
real(8) :: dnew(ndof_solid),vnew(ndof_solid),anew(ndof_solid)
real(8) :: dtil(ndof_solid),vtil(ndof_solid)
real(8) :: res
integer :: w
real(8) :: mlagnew(nsol_ebc)
real(8) :: rpen(nsol_ebc)
real(8) :: dis(ndof_solid),vel(ndof_solid),acc(ndof_solid),mlag(nsol_ebc)
real(8) :: solid_dis(nsd_solid,nn_solid)
real(8) :: solid_vel(nsd_solid,nn_solid)
real(8) :: solid_acc(nsd_solid,nn_solid)
!=====
real(8) solid_coor_curr(nsd_solid,nn_solid)
real(8) solid_stress(nsd_solid*2,nn_solid)
real(8) solid_bcforce(nsd_solid,nn_solid)
real(8) solid_fint(nsd_solid,nn_solid)
real(8) solid_fkin(nsd_solid,nn_solid)
integer :: ien_sbc(ne_sbc,nen_solid+2)
!-------------------------------------------------
real(8) sq_solid(0:3,8,8)
real(8) ffsi(nsd_solid,nn_solid)
integer iq
!Convert to single Column
do a=1,nn_solid
	do i=1,nsd_solid
		p=nsd_solid*(a-1)+i 
		dis(p)=solid_dis(i,a)
		vel(p)=solid_vel(i,a)
		acc(p)=solid_acc(i,a)
	enddo
enddo
call apply_2ndbc_solid2d(dis,ien_sbc,ne_sbc,solid_bcforce,solid_stress)
if (myid==0) then
write(*,*) 'FX= ', sum(solid_bcforce(1,:)),'FY= ', sum(solid_bcforce(2,:))
end if
fext(:)=0.0d0
do a=1,nn_solid
	p=nsd_solid*(a-1)+1
	q=nsd_solid*(a-1)+2
	fext(p)=fext(p)-solid_bcforce(1,a)
	fext(q)=fext(q)-solid_bcforce(2,a)
enddo
!write(*,*) 'fext=',fext
! =======Predict=========
    dtil=dis+dt*vel+((dt**2)/2.0d0)*(1.0d0-2.0d0*beta)*acc
    vtil=vel+((1.0d0-gamma)*dt)*acc
    !======Correct========
    w=0;
    dnew(:)=dtil(:)
    vnew(:)=vtil(:)
    mlagnew(:)=mlag(:)
	res=huge(1.d0)
	!=======Newton-Raphson iteration======
	do while ((res .ge. tol) .and. (w .le. 100))
		ka(:,:)=0.0d0
		rf(:)=0.0d0
		kt(:,:)=0.0d0
	 	w=w+1
        anew=(1/(beta*(dt**2)))*(dnew-dtil)
        vnew=vtil+gamma*dt*anew
		!internal forces and tangent stiffness matrix
		call s_int(dnew,fint,kt,ka,sel)
		!Damping Forces
		call s_dam(vnew,kt,fdam)
		!Kinematic Forces
		call s_kin(anew,fkin,ka)
		!Lagrange Multiplier (Reaction) Forces
		call s_pen(dnew,mlagnew,rpen,fpen,ka)
		rf(1:ndof_solid)=fext-fint-fpen-fkin-fdam
		rf(ndof_solid+1:neq_solid)=rpen
		res=sqrt(sum(rf**2))/ne_solid
		!include "printsolid.fi"
		!Linear System Solver
		call dgesv(neq_solid, 1, ka, neq_solid, IPIV, rf, neq_solid, INFO )
		!Add Increments
		dnew(:)=dnew(:)+rf(1:ndof_solid)
		mlagnew(:)=mlagnew(:)+rf(ndof_solid+1:neq_solid)
	enddo
	!Check Solver for Sucessful completion
	if (myid==0) then
		if (w .ge. 100) then
			write(*,*) 'solid solver was unable to converge'
			stop
		else if (w .gt. 1) then
			write(*,*) 'succesful exit after ', w, 'iterations'
		endif
	endif
!Convert to nsd Columns
do a=1,nn_solid
	do i=1,nsd_solid
		p=nsd_solid*(a-1)+i
		solid_coor_curr(i,a)=dnew(p)+xref_solid(p) 
		solid_dis(i,a)=dnew(p)
		solid_vel(i,a)=vnew(p)
		solid_acc(i,a)=anew(p)
		solid_fint(i,a)=fint(p)
		solid_fkin(i,a)=fkin(p)
	enddo
enddo
mlag=mlagnew
end subroutine solve_solid
