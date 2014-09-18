subroutine getM(bf,fext)
use solid_variables
implicit none 
integer :: el,al,bl,pg,qg,i,pl,ql,gp
real(8) :: mel(eldof_solid,eldof_solid)
real(8) :: bf(nsd_solid),fext(ndof_solid)
allocate(Mass(ndof_solid,ndof_solid))
allocate(km(ndof_solid,ndof_solid))
Mass(:,:)=0.0d0
do el=1,ne_solid
    mel(:,:)=0.0d0
    do gp=1,nquad_solid
        do al=1,nen_solid
            do bl=1,nen_solid
                mel(al,bl)=mel(al,bl)+rho*shp(0,al,gp)*shp(0,bl,gp)*detjac(el)
            enddo
        enddo
    enddo
    do al=1,nen_solid
        do bl=1,nen_solid
            do i=1,nsd_solid
                pl=nsd_solid*(al-1)+i;
                ql=nsd_solid*(bl-1)+i;
                pg=ien(pl,el);
                qg=ien(ql,el);
                Mass(pg,qg)=Mass(pg,qg)+mel(al,bl);
            enddo
        enddo
    enddo
    do al=1,nen_solid
        do i=1,nsd_solid
            pl=nsd_solid*(al-1)+i;
            pg=ien(pl,el);
            fext(pg)=fext(pg)+rho*bf(i)
        enddo
    enddo
enddo
km(:,:)=(1/(beta*dt**2))*Mass(:,:)

        open(unit=33,file='km.out')
        open(unit=26,file='fext.out')
        open(unit=27,file='mass.out')
        do i=1,ndof_solid
        write(26,*) fext(i)
        write(27,*) Mass(i,:)
        write(33,*) km(i,:)
        enddo

end subroutine getM