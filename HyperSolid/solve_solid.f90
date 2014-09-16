subroutine solve_solid(fext,dis,vel,acc,dt,km,Mass)
use solid_variables
implicit none
real(8) :: dt
real(8) :: fint(ndof_solid),fkin(ndof_solid),fdam(ndof_solid),fpen(ndof_solid)
real(8) :: rc1,rc2,kappa,CMR(3,3),sel(3,ne_solid)
real(8) :: xref(ndof_solid)
real(8) :: detjac(ne_solid)
real(8) :: nx(nen_solid,nsd_solid,ne_solid)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el
real(8) :: ka(neq_solid,neq_solid),kt(ndof_solid,ndof_solid),kl(eldof,eldof),knl(eldof,eldof),ssuper(4,4)
integer :: INFO
real(8) :: rf(neq_solid),fext(ndof_solid)
real(8) :: dis(ndof_solid),vel(ndof_solid),acc(ndof_solid)
real(8) :: dnew(ndof_solid),vnew(ndof_solid),anew(ndof_solid)
real(8) :: dtil(ndof_solid),vtil(ndof_solid)
real(8) :: beta,gamma,tol
real(8) :: ipiv
real(8) :: res
integer :: w
real(8) :: mlagnew(nsol_ebc)
real(8) :: mlag(nsol_ebc)
real(8) :: rpen(nsol_ebc)
real(8) :: gdof(nsol_ebc)
real(8) :: gx(nsol_ebc)
real(8) :: km(ndof_solid,ndof_solid)
real(8) :: Mass(ndof_solid,ndof_solid)
integer :: ien(eldof,ne_solid)
beta=0.25d0
gamma=0.5d0
tol=1.0d-6
neq_solid=ndof_solid+nsol_ebc
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
		call s_int(neq_solid,nsd_solid,nn_solid,nen_solid,ne_solid,nquad_solid,&
			ndof_solid,eldof,nsol_ebc,rc1,rc2,kappa,xref,dnew,nx,detjac,ien,fint,kt,ka,sel)
		call s_dam(ndof_solid,damp_solid,vnew,Mass,kt,fdam)
		call s_kin(neq_solid,ndof_solid,anew,Mass,km,fkin,ka)
		call s_pen(ndof_solid,neq_solid,nsol_ebc,dnew,gdof,gx,kappa,mlagnew,rpen,fpen,ka)
		rf(1:ndof_solid)=fext-fint-fpen-fkin-fdam
		rf(ndof_solid+1:neq_solid)=rpen
		res=sqrt(sum(rf**2))/ne_solid
		write(12,*) res
		write(10,*) '----residual------'
		do i=1,neq_solid
			write(10,*) rf(i)
		enddo
		call dgesv(neq_solid, 1, ka, neq_solid, IPIV, rf, neq_solid, INFO )
		if (INFO .ne. 0) then
			write(*,*) 'error: unable to solve for incremental displacements'
			write(*,*) INFO
			stop
		end if
		dnew=dnew+rf(1:ndof_solid)
		mlagnew=mlagnew+rf(ndof_solid+1:neq_solid)
	enddo
	if (w .gt. 1) then
		write(*,*) 'succesful exit after ', w, 'iterations'
	endif
	dis=dnew
	vel=vnew
	acc=anew
	mlag=mlagnew
	end subroutine solve_solid