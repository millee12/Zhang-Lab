subroutine mooney(F,rc1,rc2,kappa,CMR,pk2,ssuper)
real(8) :: t13,t14,t12,t23,t43,t32,t53,t73,t49,t83,t109
real(8) :: F(2,2), RCG(2,2),pk2(3),rc1,rc2,kappa,CMR(3,3),ssuper(4,4)
real(8) :: I1,I2,I3,dI1(3),dI2(3),dI3(3),ddI1(3,3),ddI2(3,3),ddI3(3,3)
real(8) :: J1,J2,J3,dJ1(3),dJ2(3),dJ3(3),ddJ1(3,3),ddJ2(3,3),ddJ3(3,3)
! Calculate stress and stiffness
t13=-1.0d0/3.0
t14=-1.0d0/4.0
t12=-1.0d0/2.0
t23=-2.0/3.0
t43=-4.0/3.0
t32=-3.0/2.0
t53=-5.0/3.0
t73=-7.0/3.0
t49=-4.0/9.0
t83=-8.0/3.0
t109=-10.0d0/9.0
RCG(1,1)=F(1,1)**2+F(2,1)**2
RCG(1,2)=F(1,1)*F(1,2)+F(2,2)*F(2,1)
RCG(2,1)=F(1,2)*F(1,1)+F(2,2)*F(2,1)
RCG(2,2)=F(1,2)**2+F(2,2)**2
!Invariants
I1=RCG(1,1)+RCG(2,2)+1.0d0
I2=0.5*(I1**2-(RCG(1,1)*RCG(1,1)+RCG(2,1)*RCG(2,1)+RCG(1,2)*RCG(1,2)+RCG(2,2)*RCG(2,2))-1.0d0)
I3=RCG(1,1)*RCG(2,2)-RCG(1,2)*RCG(2,1)

!Reduced Invariants 
J1=I1*I3**t13
J2=I2*I3**t23
J3=I3**(0.5)

! (I)*
dI1(1)=2.0
dI1(2)=2.0
dI1(3)=0.0d0
dI2(1)=2.0*RCG(2,2)+2.0
dI2(2)=2.0*RCG(1,1)+2.0
dI2(3)=-1.0d0*(RCG(1,2)+RCG(2,1))
dI3(1)=2.0*RCG(2,2)
dI3(2)=2.0*RCG(1,1)
dI3(3)=-1.0d0*RCG(1,2)-RCG(2,1)
!(J)*
do i=1,3
	dJ1(i)=(I3**t13)*dI1(i)+t13*I1*(I3**t43)*dI3(i)
	dJ2(i)=(I3**t23)*dI2(i)+t23*I2*(I3**t53)*dI3(i)
	dJ3(i)=-1.0d0*t12*(I3**t12)*dI3(i)
enddo

!(I)**
ddI1(:,:)=0.0d0

ddI2(:,:)=0.0d0
ddI2(1,2)=4.0
ddI2(2,1)=4.0
ddI2(3,3)=-8.0


ddI3(:,:)=0.0d0
ddI3(1,2)=4.0
ddI3(2,1)=4.0
ddI3(3,3)=-8.0

! (J)**
do i=1,3
    do j=1,3
        ddJ1(i,j)=(I3**(t13))*ddI1(i,j)&
            +t13*I3**(t43)*(dI1(i)*dI3(j)+dI3(i)*dI1(j)+I1*ddI3(i,j))&
            -t49*(I1*I3**(t73))*dI3(i)*dI3(j)
        ddJ2(i,j)=(I3**t23)*ddI2(i,j)&
            +t23*(I3**(t53))*(dI2(i)*dI3(j)+dI3(i)*dI2(j)+I2*ddI3(i,j))&
            -t109*(I2*(I3**t83))*dI3(i)*dI3(j)
        ddJ3(i,j)=t14*(I3**t32)*(dI3(i)*dI3(j))-t12*(I3**t12)*ddI3(i,j)
    enddo
enddo
!PK2 Stress, Voight Notation
pk2=rc1*dJ1+rc2*dJ2+kappa*(J3-1.0d0)*dJ3
!Super Stress Matrix for assembling nonlinear component
!of tangent stiffness matrix
ssuper(:,:)=0.0d0
ssuper(1,1)=pk2(1)
ssuper(1,2)=pk2(3)
ssuper(2,1)=pk2(3)
ssuper(2,2)=pk2(2)
ssuper(3,3)=pk2(1)
ssuper(3,4)=pk2(3)
ssuper(4,3)=pk2(3)
ssuper(4,4)=pk2(2)
!Calculate Constitutive Matrix
do i=1,3
    do j=1,3
        CMR(i,j)=rc1*ddJ1(i,j)+rc2*ddJ2(i,j) &
            +kappa*(dJ3(i)*dJ3(j)+(J3-1.0d0)*ddJ3(i,j))
    enddo
enddo

end