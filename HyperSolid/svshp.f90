subroutine svshp
use solid_variables
implicit none
integer :: i,j,gp,a,p,dof,el
real(8) :: jac(nsd_solid,nsd_solid),jacinv(nsd_solid,nsd_solid)
real(8) :: xigp(nquad_solid,nsd_solid)
real(8) ::xi(nen_solid,nsd_solid)
allocate(shp(0:2,nen_solid,nquad_solid))
allocate(nx(nen_solid,nsd_solid,nquad_solid,ne_solid))
allocate(detjac(ne_solid))
xi(:,1)=[1.0d0,-1.0d0,-1.0d0,1.0d0]
xi(:,2)=[1.0d0,1.0d0,-1.0d0,-1.0d0]
xigp(:,1)=xi(:,1)*(3.0d0**(-0.5d0))
xigp(:,2)=xi(:,2)*(3.0d0**(-0.5d0))
! Find Mesh Scaling
do el=1,ne_solid
    do gp=1,nquad_solid
        jac(:,:)=0.0d0
        do a=1,nen_solid
        !Make Shape Functions
            shp(0,a,gp)=0.25d0*(1.0d0+xi(a,1)*xigp(gp,1))*(1.0d0+xi(a,2)*xigp(gp,2))
            shp(1,a,gp)=0.25d0*xi(a,1)*(1.0d0+xi(a,2)*xigp(gp,2))
            shp(2,a,gp)=0.25d0*xi(a,2)*(1.0d0+xi(a,1)*xigp(gp,1))
            do i=1,nsd_solid
                p=nsd_solid*(a-1)+i
                dof=lm_solid(p,el)
                do j=1,nsd_solid
                   jac(i,j)=jac(i,j)+shp(j,a,gp)*xref_solid(dof)
                enddo
            enddo
        enddo
        detjac(el)=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
        jacinv(1,1)=jac(2,2)/detjac(el)
        jacinv(2,2)=jac(1,1)/detjac(el)
        jacinv(1,2)=-1.0d0*jac(1,2)/detjac(el)
        jacinv(2,1)=-1.0d0*jac(2,1)/detjac(el)
        ! Shape function derivatives in the reference domain
        do a=1,nen_solid
            nx(a,1,gp,el)=shp(1,a,gp)*jacinv(1,1)+shp(2,a,gp)*jacinv(2,1);
            nx(a,2,gp,el)=shp(1,a,gp)*jacinv(1,2)+shp(2,a,gp)*jacinv(2,2);
        enddo
    enddo
enddo
end subroutine svshp