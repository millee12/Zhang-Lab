subroutine apply_2ndbc_solid2d(dis,ien_sbc,ne_sbc,solid_bcforce,solid_stress)
!---------------------------------
! Calculate solid surface normal integral
use solid_variables, only: ien_solid,ne_solid,nsd_solid,nen_solid,nn_solid,xyz_solid,ndof_solid
use r_common, only: h
implicit none

real(8) x_solid(nsd_solid,nn_solid)
integer ien_sbc(ne_sbc,nen_solid+2)
integer ne_sbc
integer a,i,j,p
real(8) solid_bcforce(nsd_solid,nn_solid)
real(8) solid_stress(nsd_solid*2,nn_solid)
!---------------------------------
real(8) x(nsd_solid,nen_solid)
real(8) xj(nsd_solid,nsd_solid)
real(8) xji(nsd_solid,nsd_solid)
real(8) det
real(8) rs(nsd_solid)
integer tmp_ien(nen_solid)
real(8) tmpx(nen_solid,nsd_solid)
integer ibs
integer nos
integer ntem
integer ine
integer iq
integer bcnode1
integer bcnode2
real(8) out_norm(nsd_solid)
integer isd
integer jsd
real(8) w
real(8) tot_len
integer snode
real(8) stress_tmp(nsd_solid,nsd_solid)
real(8) dis(ndof_solid)
!!!! Remove this!!!!
!!!!!!!!!!!!!!!!!!!!
iq=1
tot_len=0.0
solid_bcforce(:,:)=0.0
out_norm(:)=0.0
do a=1,nn_solid
        do j=1,nsd_solid
                p=nsd_solid*(a-1)+j
                x_solid(j,a)=xyz_solid(j,a)+dis(p)
        enddo
enddo
!write(*,*) nn_solid
!write(*,*) 'x_solid= '
!do i=1,nn_solid
!        write(*,*) x_solid(:,i)
!enddo
!write(*,*) '======='
!write(*,*) ne_sbc
do ibs=1,ne_solid
!write(*,*) ien_solid(ibs,:)
enddo
do ibs=1,ne_sbc
ine=ien_sbc(ibs,1)
!	write(*,*) 'iensbc', ien_sbc(ibs,:)
		do nos=1,nen_solid
!	        write(*,*) 'ien_solid@ine=',ine,'and nos=',nos,'is',ien_solid(ine,nos)
		ntem=ien_solid(ine,nos) !...connectivity 
       	x(1:nsd_solid,nos)   = x_solid(1:nsd_solid,ntem)
		end do
if (nsd_solid == 2) then
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

        do isd=1,nsd_solid
                do jsd=1,nen_solid
                        tmpx(jsd,isd)=x(isd,jsd)
                end do
        end do
        !write(*,*) tmp_ien
!       do jsd=1,nen_solid
!               write(*,*) ien_solid(ine,jsd), tmpx(jsd,:)
!       enddo
!write(*,*) tmpx,tmp_ien
call outnormal_2d(tmpx,tmp_ien,out_norm,nsd_solid,nen_solid,w)

!write(*,*) 'out_norm', ibs, out_norm(:), w

	tot_len=tot_len+w

if (ien_sbc(ibs,nen_solid+2) == -999) then

	snode=ien_solid(ine,bcnode1)

                stress_tmp(1,1)=solid_stress(1,snode)
                stress_tmp(2,2)=solid_stress(2,snode)
                stress_tmp(1,2)=solid_stress(3,snode)  
                stress_tmp(2,1)=solid_stress(3,snode)
continue
	do isd=1,nsd_solid	
		do jsd=1,nsd_solid
		solid_bcforce(isd,snode)=solid_bcforce(isd,snode)+&
		stress_tmp(isd,jsd)*out_norm(jsd)*w*h(bcnode1)
		end do
	end do

        snode=ien_solid(ine,bcnode2)

                stress_tmp(1,1)=solid_stress(1,snode)
                stress_tmp(2,2)=solid_stress(2,snode)
                stress_tmp(1,2)=solid_stress(3,snode)         
                stress_tmp(2,1)=solid_stress(3,snode)

        do isd = 1,nsd_solid
		  do jsd = 1,nsd_solid
                solid_bcforce(isd,snode)=solid_bcforce(isd,snode)+&
		          stress_tmp(isd,jsd)*out_norm(jsd)*w*h(bcnode2)
!		write(*,*) 'norm',out_norm(jsd)  
		end do
        end do
end if
	
end do
!write(*,*) solid_bcforce

return
end 
