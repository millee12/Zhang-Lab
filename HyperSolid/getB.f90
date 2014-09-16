subroutine getB(el,ne,ndof,eldof,nen,ien,nsd_solid,nx,dis,BL,BNL)
implicit none
integer :: ndof,eldof,nsd_solid,nen,a,p,el,ien(eldof,ne),ne
real(8) :: dis(ndof),uel(eldof),l(nsd_solid,nsd_solid),nx(nen,nsd_solid)
real(8) :: BL(3,eldof),BL0(3,eldof),BL1(3,eldof),BNL(4,eldof)
do a=1,eldof
	p=ien(a,el)
	uel(a)=dis(p)
enddo
!Make strain-displacement matrices
l(1,1)=nx(1,1)*uel(1)+nx(2,1)*uel(3)+nx(3,1)*uel(5)+nx(4,1)*uel(7)
l(2,2)=nx(1,2)*uel(2)+nx(2,2)*uel(4)+nx(3,2)*uel(6)+nx(4,2)*uel(8)
l(2,1)=nx(1,1)*uel(2)+nx(2,1)*uel(4)+nx(3,1)*uel(6)+nx(4,1)*uel(8)
l(1,2)=nx(1,2)*uel(1)+nx(2,2)*uel(3)+nx(3,2)*uel(5)+nx(4,2)*uel(7)
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
end subroutine getB 