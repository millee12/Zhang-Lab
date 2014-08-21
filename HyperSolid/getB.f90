subroutine getB(ndof,nen,nsd,nx,dg,BL,BNL)
real :: dg(ndof),l(nsd,nsd),nx(nen,nsd)
real :: BL(3,ndof),BL0(3,ndof),BL1(3,ndof),BNL(4,ndof)
integer :: i
!Make strain-displacement matrices
l(1,1)=nx(1,1)*dg(1)+nx(2,1)*dg(3)+nx(3,1)*dg(5)+nx(4,1)*dg(7)
l(2,2)=nx(1,2)*dg(2)+nx(2,2)*dg(4)+nx(3,2)*dg(6)+nx(4,2)*dg(8)
l(2,1)=nx(1,1)*dg(2)+nx(2,1)*dg(4)+nx(3,1)*dg(6)+nx(4,1)*dg(8)
l(1,2)=nx(1,2)*dg(1)+nx(2,2)*dg(3)+nx(3,2)*dg(5)+nx(4,2)*dg(7)
! Make BL0
BL0(1,:)=[nx(1,1),0.0, nx(2,1),0.0, nx(3,1), 0.0, nx(4,1), 0.0]
BL0(2,:)=[0.0, nx(1,2), 0.0, nx(2,2), 0.0, nx(3,2), 0.0, nx(4,2)]
BL0(3,:)=[nx(1,2), nx(1,1), nx(2,2), nx(2,1), nx(3,2), nx(3,1), nx(4,2),nx(4,1)]
! Make BL1
BL1(1,:)=[l(1,1)*nx(1,1), l(2,1)*nx(1,1), l(1,1)*nx(2,1), l(2,1)*nx(2,1),&
    l(1,1)*nx(3,1), l(2,1)*nx(3,1), l(1,1)*nx(4,1), l(2,1)*nx(4,1)]
BL1(2,:)=[l(1,2)*nx(1,2), l(2,2)*nx(1,2), l(1,2)*nx(2,2), l(2,2)*nx(2,2),&
    l(1,2)*nx(3,2), l(2,2)*nx(3,2), l(1,2)*nx(4,2), l(2,2)*nx(4,2)]
BL1(3,:)=[l(1,1)*nx(1,2)+l(1,2)*nx(1,1), l(2,1)*nx(1,2)+l(2,2)*nx(1,1),&
    l(1,1)*nx(2,2)+l(1,2)*nx(2,1), l(2,1)*nx(2,2)+l(2,2)*nx(2,1),&
    l(1,1)*nx(3,2)+l(1,2)*nx(3,1), l(2,1)*nx(3,2)+l(2,2)*nx(3,1),&
    l(1,1)*nx(4,2)+l(1,2)*nx(4,1), l(2,1)*nx(4,2)+l(2,2)*nx(4,1)]
!Make BL
BL=BL0+BL1
! Make BNL
BNL(1,:)=[nx(1,1), 0.0, nx(2,1), 0.0, nx(3,1), 0.0, nx(4,1), 0.0]
BNL(2,:)=[nx(1,2), 0.0, nx(2,2), 0.0, nx(3,2), 0.0, nx(4,2), 0.0]
BNL(3,:)=[0.0,nx(1,1), 0.0, nx(2,1), 0.0, nx(3,1), 0.0, nx(4,1)]
BNL(4,:)=[0.0,nx(1,2), 0.0, nx(2,2), 0.0, nx(3,2), 0.0, nx(4,2)]
	!	do i=1,3
	!	write(*,*) BL(i,:)
	!enddo
	!write(*,*) '         '
end subroutine