subroutine getM(nsd,dt,beta,gamma,rho,ndof,ngp,eldof,nen,nel,detjac,ien,shp,M,km)
implicit none 
integer :: el,nel,ndof,eldof,gp,ngp,nen,al,bl,pg,qg,i,ien(eldof,nel),pl,ql,nsd 
real(8) :: M(ndof,ndof),fkin(ndof),mel(eldof,eldof),km(ndof,ndof),detjac(nel),shp(0:nsd,nen,ngp)
real(8) :: rho,dt,beta,gamma
do el=1,nel
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
            do i=1,nsd
                pl=nsd*(al-1)+i;
                ql=nsd*(bl-1)+i;
                pg=ien(pl,el);
                qg=ien(ql,el);
                M(pg,qg)=M(pg,qg)+mel(al,bl);
            enddo
        enddo
    enddo
enddo
km=(1/(beta*dt**2))*M
end subroutine getM