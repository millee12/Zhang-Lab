subroutine svshp(ndof,eldof,nen,nsd,ngp,nel,ien,xref,nx,detjac,shp)
implicit none
integer :: nsd,nen,ngp,ng,nel,ndof,eldof
real(8) :: xref(ndof)
real(8) :: nx(nen,nsd,ngp,nel)
real(8) :: detjac(nel)
integer :: i,j,gp,a,p,dof,el
integer :: ien(eldof,nel)
real(8) :: jac(nsd,nsd),jacinv(nsd,nsd)
real(8) ::shp(0:nsd,nen,ngp)
real(8) :: xigp(ngp,nsd)
real(8) ::xi(nen,nsd)
xi(:,1)=[1.0d0,-1.0d0,-1.0d0,1.0d0]
xi(:,2)=[1.0d0,1.0d0,-1.0d0,-1.0d0]
xigp(:,1)=xi(:,1)*(3.0d0**(-0.5d0))
xigp(:,2)=xi(:,2)*(3.0d0**(-0.5d0))
! Find Mesh Scaling
do el=1,nel
    do gp=1,ngp
        jac(:,:)=0.0d0
        do a=1,nen
        !Make Shape Functions
            shp(0,a,gp)=0.25d0*(1+xi(a,1)*xigp(gp,1))*(1+xi(a,2)*xigp(gp,2))
            shp(1,a,gp)=0.25d0*xi(a,1)*(1+xi(a,2)*xigp(gp,2))
            shp(2,a,gp)=0.25d0*xi(a,2)*(1+xi(a,1)*xigp(gp,1))
            do i=1,nsd
                p=nsd*(a-1)+i
                dof=ien(p,el)
                do j=1,nsd
                   jac(i,j)=jac(i,j)+shp(j,a,gp)*xref(dof)
                enddo
            enddo
        enddo
        detjac(el)=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
        jacinv(1,1)=jac(2,2)/detjac(el)
        jacinv(2,2)=jac(1,1)/detjac(el)
        jacinv(1,2)=-1.0d0*jac(1,2)/detjac(el)
        jacinv(2,1)=-1.0d0*jac(2,1)/detjac(el)
        ! Shape function derivatives in the reference domain
        do a=1,nen
            nx(a,1,gp,el)=shp(1,a,gp)*jacinv(1,1)+shp(2,a,gp)*jacinv(2,1);
            nx(a,2,gp,el)=shp(1,a,gp)*jacinv(1,2)+shp(2,a,gp)*jacinv(2,2);
        enddo
    enddo
enddo
end subroutine svshp