!Main File to test Solid Solver Solid Solver
program main
integer :: nsd=2,nnode=4,nen=4,nel=1,ngp=4,ndof,tend=50,t
real :: xref(8),dg(8),rf(10)
integer :: i,j,gp,a,b,p,q,dof1,dof2,el,ien(8,1)
real,dimension(4,2,4,1) :: nx
real, dimension(1) :: detjac
real :: fext(8)
!For lagrange multiplier force
real :: fpen(8),gx(2),mlag(2),rpen(2),ka(10,10)
integer :: gdof(2),ng
!For the internal forces (calculated using Mooney-Rivlin)
real :: fint(8),kt(8,8)
real :: rc1,rc2,kappa
!Solve
real ::delta(10),res
integer :: INFO,nee,ninc=1,IPIV(10)
ndof=nnode*nsd

! Mesh Information
ien(:,1)=[7,8,5,6,1,2,3,4]
xref(:)=[0.0,0.0,1.0,0.0,0.0,1.0,1.0,1.0]
fext(:)=0.0
fext(3)=5000.0
fext(7)=5000.0
dg(:)=0.0
!specify essential boundaries
ng=2
gdof(:)=[1,5]
gx(:)=0.0
mlag(:)=0.0
nee=ndof+ng
!MR constants
rc1=8.6207d3
rc2=0.0
kappa=1.6667d5
call svshp(ndof,nen,nsd,ngp,nel,ien,xref,nx,detjac)
do t=1,tend
	!internal forces and tangent stiffness matrix
	call s_int(nsd,nnode,nen,nel,ngp,ndof,ng,rc1,rc2,kappa,xref,dg,nx,detjac,ien,fint,kt)
	call s_pen(ndof,ng,dg,gdof,gx,kt,kappa,mlag,rpen,fpen,ka)
	write(*,*) fint
	write(*,*) fpen
	rf(:)=[fext(:)-fint(:)-fpen(:),rpen]
	call sgesv(10, 1, ka, 10, IPIV, rf(1:10), 10, INFO )
	! SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
	if (INFO .eq. 0) then
		write(*,*) 'stiffness matrix inverted'
		else
		write(*,*) 'warning: linear solver exited with errors'
	endif
	dg(:)=dg(:)+rf(1:ndof)
	mlag(:)=mlag(:)+rf(ndof+1:ndof+ng)
	res=sqrt(sum(rf*rf))
	write(*,*) 'res= ', res
enddo	
end program main