subroutine svshp(ndof,nen,nsd,ngp,nel,ien,xref,nx,detjac)
integer :: nsd,nen,ngp,ng,nel
real :: xref(ndof)
real :: nx(nen,nsd,ngp,nel)
real :: detjac(nel)
integer :: i,j,gp,a,p,dof,el
integer :: ien(ndof,nel)
real :: jac(nsd,nsd),jacinv(nsd,nsd)
real ::shp(0:nsd,nen)
real :: xigp(ngp,nsd)
real ::xi(nen,nsd)
xi(:,1)=[1.0,-1.0,-1.0,1.0]
xi(:,2)=[1.0,1.0,-1.0,-1.0]
xigp(:,1)=xi(:,1)*(3**(-0.5))
xigp(:,2)=xi(:,2)*(3**(-0.5))
! Find Mesh Scaling
do el=1,nel
    do gp=1,ngp
        jac(:,:)=0.0
        do a=1,nen
        !Make Shape Functions
            shp(0,a)=0.25*(1+xi(a,1)*xigp(gp,1))*(1+xi(a,2)*xigp(gp,2))
            shp(1,a)=0.25*xi(a,1)*(1+xi(a,2)*xigp(gp,2))
            shp(2,a)=0.25*xi(a,2)*(1+xi(a,1)*xigp(gp,1))
            do i=1,nsd
                p=nsd*(a-1)+i
                dof=ien(p,el)
                do j=1,nsd
                   jac(i,j)=jac(i,j)+shp(j,a)*xref(dof)
                enddo
            enddo
        enddo
        detjac(el)=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
        jacinv(1,1)=jac(2,2)/detjac(el)
        jacinv(2,2)=jac(1,1)/detjac(el)
        jacinv(1,2)=-1.0*jac(1,2)/detjac(el)
        jacinv(2,1)=-1.0*jac(2,1)/detjac(el)
        ! Shape function derivatives in the reference domain
        do a=1,nen
            nx(a,1,gp,el)=shp(1,a)*jacinv(1,1)+shp(2,a)*jacinv(2,1);
            nx(a,2,gp,el)=shp(1,a)*jacinv(1,2)+shp(2,a)*jacinv(2,2);      
        enddo
    enddo
enddo
end subroutine svshp