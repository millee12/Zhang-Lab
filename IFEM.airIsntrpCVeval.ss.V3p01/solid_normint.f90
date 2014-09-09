subroutine solid_normint(x_solid,nsd,nn_solid,ien_sbc,ne_sbc,nen_solid,&
			ien_solid,ne_solid,solid_norm)
!---------------------------------
! Calculate solid surface normal integral
use solid_variables, only: wq_solid, nquad_solid, xq_solid
use r_common, only: h
implicit none

real(8) x_solid(nsd,nn_solid)
integer nsd
integer nn_solid
integer ien_sbc(ne_sbc,nen_solid+2)
integer ne_sbc
integer nen_solid
integer ien_solid(ne_solid,nen_solid)
integer ne_solid
real(8) solid_norm(nsd,nn_solid)
!---------------------------------
real(8) x(nsd,nen_solid)
real(8) xj(nsd,nsd)
real(8) xji(nsd,nsd)
real(8) det
real(8) rs(nsd)
integer tmp_ien(nen_solid)
real(8) tmpx(nen_solid,nsd)
integer ibs
integer nos
integer ntem
integer ine
integer iq
integer bcnode1
integer bcnode2
real(8) out_norm(nsd)
integer isd
integer jsd
real(8) w
real(8) tot_len
integer snode


iq=1
tot_len=0.0
solid_norm(:,:)=0.0
out_norm(:)=0.0

do ibs=1,ne_sbc

!write(*,*) 'iensbc', ien_sbc(ibs,:)

        ine=ien_sbc(ibs,1)
		do nos=1,nen_solid
	        ntem=ien_solid(ine,nos) !...connectivity
        	x(1:nsd,nos)   = x_solid(1:nsd,ntem)
		end do

if (nsd == 2) then
        if (nen_solid == 3) then ! triangle case
                if (ien_sbc(ibs,2) == 0) then ! node 2,3 on the edge
                        rs(1)=0.0
                        rs(2)=0.5
                        bcnode1=2
                        bcnode2=3
                end if

                if (ien_sbc(ibs,3) == 0) then ! node 1,3 on the edge
                        rs(1)=0.5
                        rs(2)=0.0
                        bcnode1=1
                        bcnode2=3
                end if

                if (ien_sbc(ibs,4) == 0) then ! node 1,2 on the edge
                        rs(1)=0.5
                        rs(2)=0.5
                        bcnode1=1
                        bcnode2=2
                end if
        end if
        if(nen_solid == 4) then ! quad case
                if (ien_sbc(ibs,2)==1 .and. ien_sbc(ibs,3)==1) then ! node 1,2 on the edge
                        rs(1)=0.0
                        rs(2)=-1.0
                        bcnode1=1
                        bcnode2=2
                end if
        
                if (ien_sbc(ibs,3)==1 .and. ien_sbc(ibs,4)==1) then ! node 2,3 on the edge
                        rs(1)=1.0
                        rs(2)=0.0
                        bcnode1=2
                        bcnode2=3
                end if

                if (ien_sbc(ibs,4)==1 .and. ien_sbc(ibs,5)==1) then ! node 3,4 on the edge
                        rs(1)=0.0
                        rs(2)=1.0
                        bcnode1=3
                        bcnode2=4
                end if

                if (ien_sbc(ibs,5)==1 .and. ien_sbc(ibs,2)==1) then ! node 4,1 on the edge
                        rs(1)=-1.0
                        rs(2)=0.0
                        bcnode1=4
                        bcnode2=1
                end if

        end if

else
write(*,*) '********NOT READY FOR 3-D NOW**************'
end if

        call r_element(rs)
!	call r_jacob(x,xj,xji,det)

        tmp_ien(1:nen_solid)=ien_sbc(ibs,2:nen_solid+1)

        do isd=1,nsd
                do jsd=1,nen_solid
                        tmpx(jsd,isd)=x(isd,jsd)
                end do
        end do

        call outnormal_2d(tmpx,tmp_ien,out_norm,nsd,nen_solid,w)

!write(*,*) 'out_norm', ibs, out_norm(:), w


	tot_len=tot_len+w
	snode=ien_solid(ine,bcnode1)
	do isd=1,nsd
		solid_norm(isd,snode)=solid_norm(isd,snode)+out_norm(isd)*w*h(bcnode1)
	end do

        snode=ien_solid(ine,bcnode2)
        do isd=1,nsd
                solid_norm(isd,snode)=solid_norm(isd,snode)+out_norm(isd)*w*h(bcnode2)
        end do
	
	
end do
!	write(*,*) 'sum', sum(solid_norm(1,:))
!	do ibs=1,nn_solid
!	write(*,*) 'solid norm', ibs, solid_norm(:,ibs)
!	end do

return
end 
