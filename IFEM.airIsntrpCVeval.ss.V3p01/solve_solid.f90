subroutine solve_solid(fext,dis,vel,acc,mlag,sel)
use solid_variables
<<<<<<< HEAD
=======
use run_variables, only: dt
>>>>>>> 9ca02509ab857d8dbc3512aa29c539835bd53dbc
implicit none
real(8) :: fint(ndof_solid),fkin(ndof_solid),fdam(ndof_solid),fpen(ndof_solid)
real(8) :: CMR(3,3),sel(3,ne_solid)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el
real(8) :: ka(neq_solid,neq_solid),kt(ndof_solid,ndof_solid)
real(8) :: kl(eldof_solid,eldof_solid),knl(eldof_solid,eldof_solid),ssuper(4,4)
real(8) :: rf(neq_solid),fext(ndof_solid)
real(8) :: dnew(ndof_solid),vnew(ndof_solid),anew(ndof_solid)
real(8) :: dtil(ndof_solid),vtil(ndof_solid)
real(8) :: res
integer :: w
real(8) :: mlagnew(nsol_ebc)
real(8) :: rpen(nsol_ebc)
real(8) :: dis(ndof_solid),vel(ndof_solid),acc(ndof_solid),mlag(nsol_ebc)
!=====
integer ien_sbc(ne_sbc,nen_solid+2)
real(8) solid_stress(nsd_solid*2,nn_solid)
!-------------------------------------------------
real(8) sq_solid(0:3,8,8)
integer iq
! Predict
	!write(*,*) ndof_solid
    dtil=dis+dt*vel+((dt**2)/2.0d0)*(1.0d0-2.0d0*beta)*acc
    vtil=vel+((1.0d0-gamma)*dt)*acc
    	open(unit=29,file='dtil.out')
		do i=1,ndof_solid
		write(29,*) dtil(i)
		enddo
		write(29,*) '---'
		open(unit=30,file='vtil.out')
		do i=1,ndof_solid
		write(30,*) vtil(i)
		enddo
		write(30,*) '---'
    ! Correct
    w=0;
    dnew(:)=dtil(:)
    vnew(:)=vtil(:)
    mlagnew(:)=mlag(:)
	res=huge(1.d0)
	do while ((res .ge. tol) .and. (w .le. 100))
		ka(:,:)=0.0d0
		rf(:)=0.0d0
		kt(:,:)=0.0d0
	 	w=w+1
		!write(*,*) 'w= ',w
        anew=(1/(beta*(dt**2)))*(dnew-dtil)
        vnew=vtil+gamma*dt*anew
		!internal forces and tangent stiffness matrix
		call s_int(dnew,fint,kt,ka,sel)
		call s_dam(vnew,kt,fdam)
		call s_kin(anew,fkin,ka)
		call s_pen(dnew,mlagnew,rpen,fpen,ka)
		rf(1:ndof_solid)=fext-fint-fpen-fkin-fdam
		rf(ndof_solid+1:neq_solid)=rpen
		res=sqrt(sum(rf**2))/ne_solid
		open(unit=23,file='fint.out')
		write(23,*) 'w= ',w
		do i=1,ndof_solid
		write(23,*) fint(i)
		enddo
		open(unit=24,file='fpen.out')
		write(24,*) 'w= ',w
		do i=1,ndof_solid
		write(24,*) fpen(i)
		enddo
		open(unit=25,file='fkin.out')
		write(25,*) 'w= ',w
		do i=1,ndof_solid
		write(25,*) fkin(i)
		enddo
		open(unit=28,file='acc.out')
		write(28,*) 'w= ',w
		do i=1,ndof_solid
		write(28,*) anew(i)
		enddo
		open(unit=31,file='residual.out')
		do i=1,neq_solid
		write(31,*) rf(i)
		enddo
		open(unit=32,file='ka.out')
		do i=1,neq_solid
		write(32,*) ka(i,:)
		enddo
				write(32,*) '==='

		open(unit=35,file='kt.out')
		do i=1,neq_solid
		write(35,*) ka(i,:)
		enddo
				write(35,*) '==='
		call dgesv(neq_solid, 1, ka, neq_solid, IPIV, rf, neq_solid, INFO )
		if (INFO .ne. 0) then
			write(*,*) 'error: unable to solve for incremental displacements'
			write(*,*) INFO
			stop
		end if
		dnew(:)=dnew(:)+rf(1:ndof_solid)
		mlagnew(:)=mlagnew(:)+rf(ndof_solid+1:neq_solid)
	enddo
	if (w .gt. 1) then
		write(*,*) 'succesful exit after ', w, 'iterations'
	endif
	dis=dnew
	vel=vnew
	acc=anew
	mlag=mlagnew
	end subroutine solve_solid
