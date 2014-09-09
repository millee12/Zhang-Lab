subroutine nature_totpre(xloc,dloc,ien,rng,bdel,ne_inflow,bdindex,pin,res_bc,pin_s)
! Apply constant total pressure boundary condition
  use global_constants
  use fluid_variables, only: nn,nsd,nen,ne,ndf,neface,vis_liq,den_liq
! input arguments
 real(8) xloc(nsd,nn) ! coordinates of fluid nodes
 real(8) dloc(ndf,nn) ! [u,v,p] for 2-D [u,v,w,p] for 3d
 integer ien(nen,ne)  ! connectivity matrix
 integer rng(neface,ne) ! all boundary information
 integer bdel(ne_inflow) ! element inndex for those on the inflow boundary
 integer ne_inflow ! number of elements on the inflow boundary
 integer bdindex ! set which edge has inflow boundary condition in the whole fluid domain
 real(8) pin ! pressure difference for inflow
 real(8) res_bc(nsd,nn) ! residuale coming from nature B.C. integration
! local variables
real(8) x(nsd,nen)
integer ie
integer inl
real(8) d(ndf,nen)
integer iface
integer ie_inflow
integer li(2) ! the local index of unknowns only for 2-D case
integer n_edge ! number of nodes per edge
! loop variables
integer ia
real(8) vel(nsd)
real(8) v2
real(8) pin_s
real(8) len
real(8) fl
!=======================================
n_edge=2 ! set number of nodes per edge NOTE: only for 2-D case right now !!!!
res_bc(:,:)=0.0d0
fl=0.0d0
len=0.0d0
pin_s=0.0d0
!===================================================
  do ie_inflow=1,ne_inflow            ! loop over elements on the inflow boundary
     ie=bdel(ie_inflow)               ! global element index
        
	do inl=1,nen
             x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
             d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
        end do
!===========================================================
     do iface=1,neface
        ! decide the knowns and unknows based on iface
        if (rng(iface,ie) == bdindex) then
                if (nsd==2) then ! only for 2-D case right now
                        if (nen==4) then ! Quadratrial case
                                        if (iface == 1) then ! case I local nodes 1,2 on the bounndary
                                                li(1)=1
                                                li(2)=2
                                        end if
                                        if (iface == 2) then ! case II local nodes 2,3 on the boundary
                                                li(1)=2
                                                li(2)=3
                                        end if
                                        if (iface == 3) then ! case III local nodes 3,4 on the boundary
                                                li(1)=3
                                                li(2)=4
                                        end if
                                        if (iface == 4) then ! case IV local nodes 4,1 on the boundary
                                                li(1)=1
                                                li(2)=4
                                        end if
                                ! These four cases are decided based on how we define mrng.in

                        end if
                        if (nen==3) then ! Triangle case
                                        if (iface == 1) then ! case I local nodes 1,2 on the bounndary
                                                li(1)=1
                                                li(2)=2
                                        end if
                                        if (iface == 2) then ! case II local nodes 2,3 on the boundary
                                                li(1)=2
                                                li(2)=3
                                        end if
                                        if (iface == 3) then ! case III local nodes 3,4 on the boundary
                                                li(1)=1
                                                li(2)=3
                                        end if
                        end if
                                ! These three cases are decided based on how we define mrng.in
                end if
        end if
     end do ! for iface
!=========================
	dl=(x(1,li(1))-x(1,li(2)))**2 + (x(2,li(1))-x(2,li(2)))**2
	dl=sqrt(dl) ! length of the line element, Only work for 2-D case !!!!
!=========================
	vel(1:nsd)=(d(1:nsd,li(1))+d(1:nsd,li(2)))*0.5
	v2=sqrt(vel(1)**2+vel(2)**2)
! v^2 at the boundary edge
	fl=fl+v2*dl
	len=len+dl
!=========================
end do
	v2=fl/len
	pin_s=pin-0.5*v2**2*den_liq
return
end
       


