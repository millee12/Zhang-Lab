subroutine getM(nsd_solid,dt,beta,gamma,rho,bf,ndof,ngp,eldof,nen,ne,detjac,ien,shp,Mass,km,fext)
implicit none 
integer :: el,ne,ndof,eldof,gp,ngp,nen,al,bl,pg,qg,i,ien(eldof,ne),pl,ql,nsd_solid 
real(8) :: Mass(ndof,ndof),fkin(ndof),mel(eldof,eldof),km(ndof,ndof),detjac(ne),shp(0:nsd_solid,nen,ngp)
real(8) :: rho,dt,beta,gamma,bf(nsd_solid),fext(ndof)
do el=1,ne
    mel(:,:)=0.0d0
    do gp=1,ngp
        do al=1,nen
            do bl=1,nen
                mel(al,bl)=mel(al,bl)+rho*shp(0,al,gp)*shp(0,bl,gp)*detjac(el)
            enddo
        enddo
    enddo
    do al=1,nen
        do bl=1,nen
            do i=1,nsd_solid
                pl=nsd_solid*(al-1)+i;
                ql=nsd_solid*(bl-1)+i;
                pg=ien(pl,el);
                qg=ien(ql,el);
                Mass(pg,qg)=Mass(pg,qg)+mel(al,bl);
            enddo
        enddo
    enddo
    do al=1,nen
        do i=1,nsd_solid
            pl=nsd_solid*(al-1)+i;
            pg=ien(pl,el);
            fext(pg)=fext(pg)+rho*bf(i)
        enddo
    enddo
enddo
km=(1/(beta*dt**2))*Mass
end subroutine getM