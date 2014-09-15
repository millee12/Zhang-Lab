subroutine solve_solid_disp(x,x_curr,kid,ien,node_sbc,xpre1,solid_prevel,solid_preacc,solid_premlag,&
                            ien_sbc,solid_stress,solid_bcvel,ng)
                            
! subroutine to solve the solid displacement at every time step
! x  -- initial solid configuration
! x_curr --- current solid configuration
! disp --- solid displacement 
! kid --- 1st type boundary condition nodes id
! ien --- solid connective matrix
! xpre1 --- last 1 time step solid position
! xpre2 --- last 2 time step solid position
! solid_acc --- solid acceleration at n+1
! solid_preacc --- solid acceleration at n
! solid_vel --- solid velocity at n+1
! solid_prevel --- solid velocity at n
! solid_stress --- solid traction B.C.
! solid_bcvel --- solid velocity B.C.

use mpi_variables, only: myid
use solid_variables
use run_variables, only: dt
implicit none
include 'mpif.h'

real(8) x(nsd_solid,nn_solid)
real(8) x_curr(nsd_solid,nn_solid)
real(8) disp(nsd_solid,nn_solid)
integer kid(nsd_solid,nn_solid)
integer ien(ne_solid,nen_solid)
integer node_sbc(nn_sbc)
real(8) xpre1(nsd_solid,nn_solid)
real(8) solid_prevel(nsd_solid,nn_solid)
real(8) solid_preacc(nsd_solid,nn_solid)
real(8) solid_acc(nsd_solid,nn_solid)
real(8) solid_vel(nsd_solid,nn_solid)
real(8) solid_bcvel(nsd_solid,nn_solid)
integer mtype(ne_solid)
!real(8) xpre2(nsd_solid,nn_solid)
!real(8) disp1(nsd_solid,nn_solid)
!real(8) disp2(nsd_solid,nn_solid)
integer ien_sbc(ne_sbc,nen_solid+2)
real(8) solid_stress(nsd_solid*2,nn_solid)
!-------------------------------------------------
real(8) sq_solid(0:3,8,8)
integer iq
!real(8) :: xq_solid(nsdpad_solid,nquadpad_solid),wq(nquadpad_solid)
!-------------------------------------------------
integer inner
integer outer
real(8) tol
parameter (inner = 100) ! solid equation inner 50 should be sufficient
parameter (outer = 10)  ! solid equation outer 5 should be sufficient
parameter (tol = 1.0d-6)  ! solid equation outer 5 should be sufficient
!-------------------------------------------------
real(8) dg(nsd_solid,nn_solid) ! disp correction
real(8) w(nsd_solid,nn_solid)  ! pre-conditioner
real(8) p(nsd_solid,nn_solid)  ! residual vector
real(8) res
real(8) del
integer i,j,a
real(8) time
!----------------------------
        real(8) alpha
        real(8) beta
        real(8) gamma
!----------------------------
!added by Eric Miller
integer :: k,ng,gx(ng),gdof(ng),ndof,nee
real(8) :: ka(nsd_solid*nn_solid,nsd_solid*nn_solid)
real(8) :: dis(nsd_solid*nn_solid)
real(8) :: acc(nsd_solid*nn_solid)
real(8) :: vel(nsd_solid*nn_solid)
real(8) :: dtil(nsd_solid*nn_solid)
real(8) :: vtil(nsd_solid*nn_solid)
real(8) :: atil(nsd_solid*nn_solid)
real(8) :: dnew(nsd_solid*nn_solid)
real(8) :: vnew(nsd_solid*nn_solid)
real(8) :: anew(nsd_solid*nn_solid)
real(8) :: solid_premlag(nsd_solid*nn_solid)
real(8) :: mlagnew(nsd_solid*nn_solid)
real(8) :: fext(nsd_solid*nn_solid)
real(8) :: fint(nsd_solid*nn_solid)
real(8) :: fkin(nsd_solid*nn_solid)
real(8) :: fpen(nsd_solid*nn_solid)
real(8) :: fdam(nsd_solid*nn_solid)
real(8) :: rf(nsd_solid*nn_solid)
real(8) :: rpen(ng)
ndof=nsd_solid*nn_solid
nee=nn_solid*nsd_solid+ng
!Newmark Method
! define the numerical parameters
beta=0.25d0
gamma=0.5d0
!beta = (1.0 - alpha**2 ) * 0.25
!----------------------------------
!Decode
k=1
do i=1,nn_solid
  do j=1,nsd_solid
    a=nsd_solid*(i-1)+j
    dis(a)=xpre1(j,i)
    vel(a)=solid_prevel(j,i)
    acc(a)=solid_preacc(j,i)
    if (kid(j,i)==0) then 
      gdof(k)=a
      k=k+1
    end if
  enddo
enddo
!------------------------------------------------------------------

w(:,:)=0.0d0
p(:,:)=0.0d0
dg(:,:)=0.0d0 
call getnorm(dg,dg,nsd_solid*nn_solid,res)
 res = sqrt(res)
 if (myid == 0) write(*,*) '===Initial error for solid displacement===', res
    ! Predict
    dtil=dis+dt*vel+(dt**2/2.0d0)*(1.0d0-2.0d0*beta)*acc
    vtil=vel+(1.0d0-gamma)*dt*acc
    ! Correct
    k=0;
    dnew=dtil
    vnew=vtil
  res=huge(1.d0)
  do while ((res .ge. tol) .and. (k .le. 100))
    ka(:,:)=0.0d0
    w=w+1
    !write(*,*) 'w= ',w
        anew=(1/(beta*dt**2))*(dnew-dtil)
        vnew=vtil+gamma*dt*anew
    !internal forces and tangent stiffness matrix
    !call s_int(nee,nsd,nn,nen,ne,ngp,ndof,eldof,ng,rc1,rc2,kappa,xref,dnew,nx,detjac,ien,fint,kt,ka,sel)
    !call s_dam(ndof,d1,d2,vnew,Mass,kt,fdam)
    !call s_kin(nee,ndof,anew,Mass,km,fkin,ka)
    !call s_pen(ndof,nee,ng,dnew,gdof,gx,kappa,mlagnew,rpen,fpen,ka)
    rf=fext-fint-fpen-fkin-fdam
    rf(ndof+1:nee)=rpen
    res=sqrt(sum(rf**2))/nen_solid
    write(12,*) res
    write(10,*) '----residual------'
    do i=1,nee
      write(10,*) rf(i)
    enddo
    !call dgesv(nee, 1, ka, nee, IPIV, rf, nee, INFO )
    !if (INFO .ne. 0) then
    !  write(*,*) 'error: unable to solve for incremental displacements'
    !  write(*,*) INFO
    ! stop
    !endif
    dnew=dnew+rf(1:ndof)
    mlagnew=mlagnew+rf(ndof+1:ndof+ng)
  enddo
  if (k .gt. 1) then
    write(*,*) 'succesful exit after ', w, 'iterations'
  endif
  dis=dnew
  vel=vnew
  acc=anew
  mlagnew=solid_premlag

call getnorm(dg,dg,nsd_solid*nn_solid,del)
 del = sqrt(del)
if (myid == 0) write(*,*) '===solid displacement correction norm===', del
stop
!update 


!Encode
do i=1,nn_solid
  do j=1,nsd_solid
    a=nsd_solid*(i-1)+j
    x_curr(j,i)=dnew(a)
    solid_prevel(j,i)=vnew(a)
    solid_preacc(j,i)=anew(a)
  enddo
enddo
solid_premlag=mlagnew
return
end 
