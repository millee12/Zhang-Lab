! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! *.fi files are used to shorten hypo.f (keeping the overview)
! the include command reads these files and replaces the include line
! with the content of these files
subroutine hypo
use global_simulation_parameter
use global_constants
use run_variables
use delta_nonuniform
use solid_variables
use fluid_variables
use r_common, only: ninit, vis_solid, density_solid
use meshgen_fluid
use meshgen_solid
use form
use ensight_output
use mpi_variables ! call mpi variable module
use lumpedmass
implicit none
include 'mpif.h'
!==============================
! Definition of variables
integer :: klok,j
integer infdomain(nn_solid)
real(8) mass_center(2)
real(8) sfxyz(nsd,node_sfcon)
integer inode_sf
! Variables for different fluid density using by implicit form
integer mdata(nn_solid)
integer n_mdata
real(8) res_l0
real(8) del_l0
integer ie, inen
! For output pressure on the solid nodes
real(8) pre_inter(nn_solid)
!============================
! Variables for boudary equations
integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration
real(8) time
real(8) time_com
real(8) pin_s
! integer ng !number of solid BCs (EM 9/15/14)
! real(8) :: solid_mlag(nsd_solid,nn_solid)
real(8), allocatable :: solid_fint(:,:)
real(8), allocatable :: solid_fkin(:,:)
real(8), allocatable :: vsolid_fint(:,:)
real(8), allocatable :: vsolid_fkin(:,:)
real(8), allocatable :: solid_trac(:,:)
real(8), allocatable :: solid_vel(:,:)
real(8), allocatable :: solid_dis(:,:)
real(8), allocatable :: solid_acc(:,:)
real(8), allocatable :: solid_mlag(:)
real(8), allocatable :: solid_stress(:,:)
real(8), allocatable :: vsolid_vel(:,:)
real(8), allocatable :: vsolid_vel_old(:,:)
real(8), allocatable :: vsolid_acc(:,:)
real(8), allocatable :: ffsi_vs(:,:)
!============================
character(len=14) filename
character(len=14) resfilename
character(len=7) fileroot
integer pg,qg,ag,bg,al,bl
!============================
! Define local variables
include "hypo_declaration_solid.fi"
include "hypo_declaration_fluid.fi"
!============================
! Define varibales on each processor
include "hypo_declaration_part.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information
include "hypo_restart_file_check.fi"
include "hypo_prepare_solid.fi"
include "hypo_prepare_fluid.fi"
!===================================
! Prepare for MPI
!call readpartele(partele)
include "hypo_prepare_part.fi"
include "hypo_prepare_com_node.fi"
!=============================
! define the influence domain matrix
! integer infdomain(nn_solid)
call mpi_barrier(mpi_comm_world,ierror)
write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger
!=============================
ndof_solid=nn_solid*nsd_solid
eldof_solid=nen_solid*nsd_solid
!=============================
if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if
!===================================
! save the orignal position of solid nodes at fluid boundary
if (node_sfcon .ne. 0 ) then
do inode_sf=1,node_sfcon
sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
end do
end if
!if (restart == 0) then
!if (myid == 0) then
!dg(:,:)=0.0
!include 'hypo_write_output.fi'
!end if
!else
!include "hypo_restart_read.fi"
!endif
call lumpmassmatrix(x,d,dold,p,hg,ien,f_fluids,ne_intlocal,ien_intlocal,&
node_local,nn_local,fden,fvis,I_fluid,rng)
allocate(lm_solid(eldof_solid,ne_solid))
lm_solid(:,:)=0
do ie=1,ne_solid
        do al=1,nen_solid
                if (al .gt. 2) then
                        bl=al-2
                else
                        bl=al+2
                endif
                do j=1,nsd_solid
                        pg=nsd_solid*(ien_solid(ie,al)-1)+j
                        qg=nsd_solid*(bl-1)+j
                        lm_solid(qg,ie)=pg
                enddo
        enddo
enddo
!save shape functions and derivatives
allocate(xref_solid(ndof_solid))
do ag=1,nn_solid
        pg=nsd_solid*(ag-1)+1
        qg=nsd_solid*(ag-1)+2
        xref_solid(pg)=xyz_solid(1,ag)
        xref_solid(qg)=xyz_solid(2,ag)
enddo
!specify essential boundaries
!write(*,*) 'ne_sbc= ',ne_sbc_1
call s_ess(ien_sbc)
gx(:)=0.0d0
call svshp
!save mass matrix and kinematic stiffness and add body forces
call getM
!=================================================================
allocate(solid_dis(nsd_solid,nn_solid))
allocate(solid_vel(nsd_solid,nn_solid))
allocate(solid_acc(nsd_solid,nn_solid))
allocate(solid_mlag(nsol_ebc))
allocate(vsolid_vel(nsd_solid,nn_solid))
allocate(vsolid_vel_old(nsd_solid,nn_solid))
allocate(vsolid_acc(nsd_solid,nn_solid))
allocate(solid_stress(3,ne_solid))
allocate(solid_trac(nsd_solid,nn_solid))
allocate(solid_fint(nsd_solid,nn_solid))
allocate(solid_fkin(nsd_solid,nn_solid))
allocate(vsolid_fint(nsd_solid,nn_solid))
allocate(vsolid_fkin(nsd_solid,nn_solid))
allocate(ffsi_vs(nsd_solid,nn_solid))
solid_pave(:)=0.0d0
solid_vel(:,:) = 0.0d0
solid_dis(:,:) = 0.0d0
solid_acc(:,:) = 0.0d0
solid_mlag(:)=0.0d0
vsolid_vel(:,:) = 0.0d0
vsolid_vel_old(:,:)=0.0d0
vsolid_acc(:,:)=0.0d0
solid_fint(:,:)=0.0d0
vsolid_fint(:,:)=0.0d0
vsolid_fkin(:,:)=0.0d0
solid_fkin(:,:)=0.0d0
solid_bcvel_old(:,:) = 0.0
ffsi_vs(:,:)=0.0d0
! Write Initial Conditions
if (myid==0) then
		write(*,*) 'Wrote Output Files:'
		call paraout_solid(its,solid_coor_curr,solid_vel,solid_acc,solid_fint,solid_fkin,solid_trac,'pt1',143)
		call paraout_fluid(its,xref,ien,d,I_fluid,f_fluids)
		call paraout_solid(its,solid_coor_curr,vsolid_vel,vsolid_acc,vsolid_fint,vsolid_fkin,ffsi_vs,'vs1',145)
endif
!=================================================================
! time loop
!=================================================================
time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps
call mpi_barrier(mpi_comm_world,ierror)
if (myid ==0) then
write (6,*) ' '
write (6,*) 'TIME STEP = ', its
write (6,*) ' '
write (7,*) ' '
write (7,*) 'TIME STEP = ', its
write (7,*) ' '
!=================================================================
! Write restart information in binary file
!include "hypo_restart_write.fi"
end if
tt = tt + dt !....update real time
klok = klok + 1 !....update counter for output
if (myid ==0) then
write (6,'(" physical time = ",f14.10," s")') tt
write (7,'(" physical time = ",f14.10," s")') tt
end if
!--------------------------------
! Update the solid coor first
solid_coor_pre2(:,:) = solid_coor_pre1(:,:)
solid_coor_pre1(:,:) = solid_coor_curr(:,:)
! solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)
!write(*,*) 'solid_coor_curr',solid_coor_curr
solid_trac(:,:) = 0.0d0
do ie=1,nn_solid
!Temporary!
solid_trac(1:nsd_solid,ie) =0.0d0*solid_pave(ie)
end do
!-------------------------------
! correct the curr solid coor by solving the solid mon equations
!call form_solidid12(id_solidbc,nsd_solid,nn_solid,ien_sbc,ne_sbc,nen_solid,ne_solid,solid_fem_con)
if (myid==0) then
write(*,*) '******SOLVE SOLID HERE*******'
endif
call solve_solid(ien_sbc,solid_trac,solid_dis,solid_vel,solid_acc,solid_mlag,solid_coor_curr,solid_stress,solid_fint,solid_fkin)
!call solve_solid_disp(solid_coor_init,solid_coor_curr,id_solidbc,solid_fem_con,node_sbc, &
!                        solid_coor_pre1,solid_vel,solid_accel,ien_sbc,solid_trac,solid_bcvel,mtyp
!--------------------------------
! Find the fluid nodes overlapping with solid domain
call search_inf_re(solid_coor_curr,x,nn,nn_solid,nsd,ne_solid,nen_solid,solid_fem_con,&
flag_fnode,node_local,nn_local)
! Construction of the dirac deltafunctions at actual solid and fluid node positions
time=mpi_wtime()
call rkpm_nodevolume(x,nsd,nn,ien,ne,nen,ien_intlocal,ne_intlocal,dvolume,sp_radius)
call rkpm_init(solid_coor_curr,nn_solid,x,nsd,nn,dvolume,sp_radius)
time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for initiate RKPM interpolation function---', time
!==================================================================
! Solve Laplace equation get indicatior field
if (myid == 0) then
time=mpi_wtime()
lp_source(:,:)=0.0
call source_laplace(x,nn,nsd,solid_coor_curr,nn_solid,solid_fem_con,ne_solid,nen_solid,&
ien_sbc,ne_sbc,node_sbc,nn_sbc,lp_source)
call solve_laplace(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_fluid,flag_fnode)
! Update density and viscosity field for fluid
!I_fluid(:)=0.0
fden(:)=den_liq+density_solid*I_fluid(:)
fvis(:)=vis_liq+vis_solid*I_fluid(:)
time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for update indicator field---', time
!=================================================================
! Distribution of the solid forces to the fluid domain
! f^fsi(t) -> f(t)
	write(*,*) 'calculating delta'
	call delta_exchange(ffsi_vs,nn_solid,f_fluids,nn,ndelta,dvolume,nsd, &
	delta_exchange_solid_to_fluid,solid_pave,d(ndf,:))
	do ie = 1, nn
		f_fluids(:,ie) = f_fluids(:,ie) * fden(ie)
	end do
endif
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
call mpi_barrier(mpi_comm_world,ierror)
call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
call mpi_bcast(fden(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
call mpi_bcast(fvis(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
call mpi_bcast(solid_pave(1),nn_solid,mpi_double_precision,0,mpi_comm_world,ierror)
call mpi_bcast(I_fluid(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
time=mpi_wtime()
include "hypo_fluid_solver.f90"
!==================
time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for fluid solver---', time
!=================================================================
! Interpolation fluid velocity -> immersed material points
! v^f(t+dt) -> v^s(t+dt)
!if (myid==0) then
call delta_exchange(vsolid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
delta_exchange_fluid_to_solid,solid_pave,d(ndf,:))
!write(*,*) 'passed delta_exchange'
call get_fsi(solid_fint,solid_fkin,solid_dis,vsolid_vel,vsolid_vel_old,vsolid_acc,solid_pave,vsolid_fint,vsolid_fkin,ffsi_vs)
!get_fsi(fint_s,fkin_s,dis_s,vel_vs,vel_vs_old,acc_vs,p_vs,fint_vs,fkin_vs,ffsi)
!endif
vsolid_vel_old(:,:)=vsolid_vel(:,:)
call mpi_barrier(mpi_comm_world,ierror)
!=================================================================
if (myid == 0) then
	if (mod(its,ntsbout) .eq. 0) then
		write(*,*) 'Wrote Output Files:'
		call paraout_solid(its,solid_coor_curr,solid_vel,solid_acc,solid_fint,solid_fkin,solid_trac,'pt1',143)
		call paraout_fluid(its,xref,ien,d,I_fluid,f_fluids)
		call paraout_solid(its,solid_coor_curr,vsolid_vel,vsolid_acc,vsolid_fint,vsolid_fkin,ffsi_vs,'vs1',145)
	endif
endif
call mpi_barrier(mpi_comm_world,ierror)
!TEMPORARY!
ffsi_vs(:,:)=0.0d0
enddo time_loop
end subroutine hypo
