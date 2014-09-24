subroutine getB(el,dis,hx,BL,BNL)
use solid_variables
implicit none
integer :: a,p,el
real(8) :: dis(ndof_solid),uel(eldof_solid),l(nsd_solid,nsd_solid),hx(nen_solid,nsd_solid)
real(8) :: BL(3,eldof_solid),BL0(3,eldof_solid),BL1(3,eldof_solid),BNL(4,eldof_solid)
do a=1,eldof_solid
	p=lm_solid(a,el)
	uel(a)=dis(p)
enddo
!Make strain-displacement matrices
l(1,1)=hx(1,1)*uel(1)+hx(2,1)*uel(3)+hx(3,1)*uel(5)+hx(4,1)*uel(7)
l(2,2)=hx(1,2)*uel(2)+hx(2,2)*uel(4)+hx(3,2)*uel(6)+hx(4,2)*uel(8)
l(2,1)=hx(1,1)*uel(2)+hx(2,1)*uel(4)+hx(3,1)*uel(6)+hx(4,1)*uel(8)
l(1,2)=hx(1,2)*uel(1)+hx(2,2)*uel(3)+hx(3,2)*uel(5)+hx(4,2)*uel(7)
! Make BL0
BL0(1,:)=[hx(1,1),0.0d0, hx(2,1),0.0d0, hx(3,1), 0.0d0, hx(4,1), 0.0d0]
BL0(2,:)=[0.0d0, hx(1,2), 0.0d0, hx(2,2), 0.0d0, hx(3,2), 0.0d0, hx(4,2)]
BL0(3,:)=[hx(1,2), hx(1,1), hx(2,2), hx(2,1), hx(3,2), hx(3,1), hx(4,2),hx(4,1)]
! Make BL1
BL1(1,:)=[l(1,1)*hx(1,1), l(2,1)*hx(1,1), l(1,1)*hx(2,1), l(2,1)*hx(2,1),&
    l(1,1)*hx(3,1), l(2,1)*hx(3,1), l(1,1)*hx(4,1), l(2,1)*hx(4,1)]
BL1(2,:)=[l(1,2)*hx(1,2), l(2,2)*hx(1,2), l(1,2)*hx(2,2), l(2,2)*hx(2,2),&
    l(1,2)*hx(3,2), l(2,2)*hx(3,2), l(1,2)*hx(4,2), l(2,2)*hx(4,2)]
BL1(3,:)=[l(1,1)*hx(1,2)+l(1,2)*hx(1,1), l(2,1)*hx(1,2)+l(2,2)*hx(1,1),&
    l(1,1)*hx(2,2)+l(1,2)*hx(2,1), l(2,1)*hx(2,2)+l(2,2)*hx(2,1),&
    l(1,1)*hx(3,2)+l(1,2)*hx(3,1), l(2,1)*hx(3,2)+l(2,2)*hx(3,1),&
    l(1,1)*hx(4,2)+l(1,2)*hx(4,1), l(2,1)*hx(4,2)+l(2,2)*hx(4,1)]
!Make BL
BL=BL0+BL1
! Make BNL
BNL(1,:)=[hx(1,1), 0.0d0, hx(2,1), 0.0d0, hx(3,1), 0.0d0, hx(4,1), 0.0d0]
BNL(2,:)=[hx(1,2), 0.0d0, hx(2,2), 0.0d0, hx(3,2), 0.0d0, hx(4,2), 0.0d0]
BNL(3,:)=[0.0d0,hx(1,1), 0.0d0, hx(2,1), 0.0d0, hx(3,1), 0.0d0, hx(4,1)]
BNL(4,:)=[0.0d0,hx(1,2), 0.0d0, hx(2,2), 0.0d0, hx(3,2), 0.0d0, hx(4,2)]
end subroutine getB 