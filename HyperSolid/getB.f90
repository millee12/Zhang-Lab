subroutine getB(ndof,eldof,nen,nsd,nx,dis,BL,BNL)
integer :: ndof,eldof
real(8) :: dis(ndof),l(nsd,nsd),nx(nen,nsd)
real(8) :: BL(3,eldof),BL0(3,eldof),BL1(3,eldof),BNL(4,eldof)
!Make strain-displacement matrices
l(1,1)=nx(1,1)*dis(1)+nx(2,1)*dis(3)+nx(3,1)*dis(5)+nx(4,1)*dis(7)
l(2,2)=nx(1,2)*dis(2)+nx(2,2)*dis(4)+nx(3,2)*dis(6)+nx(4,2)*dis(8)
l(2,1)=nx(1,1)*dis(2)+nx(2,1)*dis(4)+nx(3,1)*dis(6)+nx(4,1)*dis(8)
l(1,2)=nx(1,2)*dis(1)+nx(2,2)*dis(3)+nx(3,2)*dis(5)+nx(4,2)*dis(7)
! Make BL0
BL0(1,:)=[nx(1,1),0.0d0, nx(2,1),0.0d0, nx(3,1), 0.0d0, nx(4,1), 0.0d0]
BL0(2,:)=[0.0d0, nx(1,2), 0.0d0, nx(2,2), 0.0d0, nx(3,2), 0.0d0, nx(4,2)]
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
BNL(1,:)=[nx(1,1), 0.0d0, nx(2,1), 0.0d0, nx(3,1), 0.0d0, nx(4,1), 0.0d0]
BNL(2,:)=[nx(1,2), 0.0d0, nx(2,2), 0.0d0, nx(3,2), 0.0d0, nx(4,2), 0.0d0]
BNL(3,:)=[0.0d0,nx(1,1), 0.0d0, nx(2,1), 0.0d0, nx(3,1), 0.0d0, nx(4,1)]
BNL(4,:)=[0.0d0,nx(1,2), 0.0d0, nx(2,2), 0.0d0, nx(3,2), 0.0d0, nx(4,2)]
end subroutine