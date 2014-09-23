subroutine read_abaqus_solid

    use solid_variables
    !use mpi_variables
    use mesh_convert_variables

    implicit none
    !include "mpif.h"
    
integer :: file_in_solid,myid=0
    integer :: file_connectsolid, file_mien_solid, file_mxyz_solid
    integer :: file_sbc_solid, file_sbcnode_solid, file_sfcon
    character(len=132) :: line
    logical :: iexist1, iexist2, iexist3, iexist4, iexist5, iexist6
    integer :: error_id1, error_id_open(6), test_error
    integer :: i, j, k, ibc, ie, iedge, inode, ierror
    real(8),allocatable :: tempreal(:)
    integer,allocatable :: connectivity(:,:),tempcon(:,:), tempmrng(:,:)
    integer :: nnsboundary(nrng_solid), nesboundary(nrng_solid)
    integer :: tempint(3)
    integer,allocatable :: nodebc(:,:), elembc(:,:), uniquenodes(:)
    integer :: node, prod_length, elem
    integer :: num_elem_edges, num_edge_nodes, num_unique_nodes
    integer :: npart, tempnesp(2), temp_nrng_solid
    integer :: sfcon_flag(4), nn_sfcon, flag_node_boundary(8)
integer,allocatable :: sfcon_nodes(:)

    file_in_solid = 9
    file_mxyz_solid = 10
    file_mien_solid = 11
    file_connectsolid = 12
    file_sbc_solid = 13
    file_sbcnode_solid = 14
    file_sfcon = 15


    open(file_in_solid,file='solid.inp',status='old',action='read', iostat=error_id1)
    if (error_id1/=0) then 
        if (myid==0) then
            write(*,*) 'No solid.inp file detected'
        endif
        !call mpi_barrier(mpi_comm_world,ierror)
        stop
    endif

    inquire (file='mxyz_solid.in', exist=iexist1)
    inquire (file='mien_solid.in', exist=iexist2)
    inquire (file='connectsolid.in', exist=iexist3)
    inquire (file='sbc_solid.in', exist=iexist4)
    inquire (file='sbcnode_solid.in', exist=iexist5)
    inquire (file='sfcon.in', exist=iexist6)

    if (iexist1) then
        error_id_open(1) = 1
        if (myid==0) then
            write(*,*) 'mxyz_solid.in file detected'
        endif
    else
        error_id_open(1) = 0
        if (myid==0) then
            write(*,*) 'No mxyz_solid.in file detected'
            open(file_mxyz_solid,file='mxyz_solid.in',status='new',action='write', iostat=ierror)
        endif
    endif

    if (iexist2) then
        error_id_open(2) = 1
        if (myid==0) then
            write(*,*) 'mien_solid.in file detected'
        endif
    else
        error_id_open(2) = 0
        if (myid==0) then
            write(*,*) 'No mien_solid.in file detected'
            open(file_mien_solid,file='mien_solid.in',status='new',action='write', iostat=ierror)
        endif
    endif

    if ((iexist3) .and. (iexist2)) then
        error_id_open(3) = 1
        if (myid==0) then
            write(*,*) 'connectsolid.in file detected'
        endif
    else
        error_id_open(3) = 0
        if (myid==0) then
            write(*,*) 'No connectsolid.in file detected'
            open(file_connectsolid,file='connectsolid.in',status='replace',action='write', iostat=ierror)
        endif
    endif

    if (iexist4) then
        error_id_open(4) = 1
        if (myid==0) then
            write(*,*) 'sbc_solid.in file detected'
        endif
    else
        error_id_open(4) = 0
        if (myid==0) then
            write(*,*) 'No sbc_solid.in file detected'
            open(file_sbc_solid,file='sbc_solid.in',status='new',action='write', iostat=ierror)
        endif
    endif

    if (iexist5) then
        error_id_open(5) = 1
        if (myid==0) then
            write(*,*) 'sbcnode_solid.in file detected'
        endif
    else
        error_id_open(5) = 0
        if (myid==0) then
            write(*,*) 'No sbcnode_solid.in file detected'
            open(file_sbcnode_solid,file='sbcnode_solid.in',status='new',action='write', iostat=ierror)
        endif
    endif

    if (iexist6) then
        error_id_open(6) = 1
        if (myid==0) then
            write(*,*) 'sfcon.in file detected'
        endif
    else
        error_id_open(6) = 0
        if (myid==0) then
            write(*,*) 'No sfcon.in file detected'
            open(file_sfcon,file='sfcon.in',status='new',action='write', iostat=ierror)
        endif
    endif


    ! Determine element type, number of spatial dimensions, number of elements,
    ! nodes, boundaries, in the solid
    
    ! Skip header
    do
        read(file_in_solid,'(a132)',iostat=error_id1) line
        if (error_id1/=0) then
            if (myid==0) then
                write(*,*) 'Error reading solid.inp'
            endif
            !call mpi_barrier(mpi_comm_world,ierror)
            stop
        endif

        if (index('*Node',line(1:5))/=0) exit
    enddo

    ! Find nsd_solid

    if (myid==0) then
        write(*,*) 'Finding nsd_solid'
    endif

    read(file_in_solid,'(a132)',iostat=error_id1) line
    if (error_id1/=0) then
        if (myid==0) then
            write(*,*) 'Error reading solid.inp'
        endif
        !call mpi_barrier(mpi_comm_world,ierror)
        stop
    endif
    k = 132
    do while (index('0123456789',line(k:k))==0)
        k = k-1
    enddo

    nsd_solid = 0
    j = index(line(1:k),',')
    i = 1
    do while (i/=0)
        nsd_solid = nsd_solid+1
        i = index(line((j+1):k),',')
        j = j + i
    enddo

    if (myid==0) then
        write(*,'(" nsd_solid = ",i5)') nsd_solid
    endif



    ! Read number of nodes and write mxyz_solid.in
  
    if (myid==0) then
        if (error_id_open(1)==0) then
            write(*,*) 'Writing mxyz_solid.in'
        endif
        write(*,*) 'Finding nn_solid_1 and nen_solid'
    endif

    allocate(tempreal(nsd_solid))

    nn_solid_1 = 0

    do
        nn_solid_1 = nn_solid_1+1

        if (error_id_open(1)==0) then
            j = index(line(1:k),',')
            read(line((j+1):k),*, iostat=test_error) tempreal(1:nsd_solid)
  
            if ((nsd_solid==2) .and. (myid==0)) then
                write(file_mxyz_solid,800) tempreal(1:nsd_solid)
            elseif ((nsd_solid==3) .and. (myid==0)) then
                write(file_mxyz_solid,810) tempreal(1:nsd_solid)
            endif
        endif

        read(file_in_solid,'(a132)',iostat=error_id1) line
        if (error_id1/=0) then
            if (myid==0) then
                write(*,*) 'Error reading solid.inp'
            endif
            !call mpi_barrier(mpi_comm_world,ierror)
            stop
        endif
        k = 132
        do while (index('0123456789',line(k:k))==0)
            k = k-1
        enddo

        if (index(line(1:1),'*')/=0) exit
    enddo

    if (myid==0) then
        close(file_mxyz_solid)
    endif

    deallocate(tempreal)

 800 format (2F18.10)
 810 format (3F18.10)  

    ! Read element type
    if (index(line,'CPS4R')/=0) then
        elem_type_solid = 2
    elseif (index(line,'CPS3')/=0) then
        elem_type_solid = 1
    elseif (index(line,'')/=0) then
        elem_type_solid = 3
    elseif (index(line,'')/=0) then
        elem_type_solid = 4
    else
        if (myid==0) then
            write(*,*) 'Unknown solid element type'
            write(*,*) line
        endif
        !call mpi_barrier(mpi_comm_world,ierror)
        stop
    endif

    nen_solid = NdPerElem(elem_type_solid)

    if (myid==0) then
        write(*,'(" nn_solid_1 = ",i7)') nn_solid_1
        write(*,'(" nen_solid = ",i5)') nen_solid
    endif



    ! Read number of elements and write mien_solid.in

    if (myid==0) then
        if (error_id_open(2)==0) then
            write(*,*) 'Writing mien_solid.in and connectsolid.in'
        endif
        write(*,*) 'Finding ne_solid_1'
    endif

    if ((error_id_open(2)==0) .or. (error_id_open(4)==0) .or. (error_id_open(5)==0) .or. (error_id_open(6)==0)) then
        i = 5*nn_solid_1
        if (allocated(tempcon)) then
            deallocate(tempcon)
        endif
        allocate(tempcon(i,nen_solid+1))
        tempcon(1:i,1:(nen_solid+1)) = 0
    endif

    ne_solid_1 = 0

    read(file_in_solid,'(a132)',iostat=error_id1) line
    k = 132
    do while (index('0123456789',line(k:k))==0)
        k = k-1
    enddo

    do
        ne_solid_1 = ne_solid_1+1

        if ((error_id_open(2)==0) .or. (error_id_open(4)==0) .or. &
            (error_id_open(5)==0) .or. (error_id_open(6)==0)) then
            
            if (ne_solid_1>i) then
          
                allocate(connectivity(i,nen_solid+1))
                connectivity(1:i,1:(nen_solid+1)) = tempcon(1:i,1:(nen_solid+1))
  
                deallocate(tempcon)
                allocate(tempcon(2*i,nen_solid+1))
  
                tempcon(1:i,1:(nen_solid+1)) = connectivity(1:i,1:(nen_solid+1))
                tempcon((i+1):(2*i),1:(nen_solid+1)) = 0
  
                deallocate(connectivity)
  
                i = 2*i
            endif
        
            j = index(line(1:k),',')
            read(line((j+1):k),*, iostat=test_error) tempcon(ne_solid_1,1:nen_solid)

        endif

        read(file_in_solid,'(a132)',iostat=error_id1) line
        if (error_id1/=0) then
            if (myid==0) then
                write(*,*) 'Error reading solid.inp'
            endif
            !call mpi_barrier(mpi_comm_world,ierror)
            stop
        endif
        k = 132
        do while (index('0123456789',line(k:k))==0)
            k = k-1
        enddo

        if (index(line(1:1),'*')/=0) exit
    enddo

    if (myid==0) then
        write(*,'(" ne_solid_1 = ",i7)') ne_solid_1
    endif




    npart = 0
    tempnesp(1:2) = 0

    do while (index('*Nset, nset=sbc',line(1:15))==0)
      
        npart = npart + 1


        do while (index('*Elset',line(1:6))==0)

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif


        enddo

        if (index(line,'generate')/=0) then

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif
            k = 132
            do while (index('0123456789',line(k:k))==0)
                k = k-1
            enddo

            read(line(1:k),*,iostat=test_error) tempint(1:3)
            tempnesp(npart) = 1 + (tempint(2)-tempint(1))/tempint(3)

            if ((error_id_open(2)==0) .or. (error_id_open(4)==0) .or. &
                (error_id_open(5)==0) .or. (error_id_open(6)==0)) then
            
                do i=1,tempnesp(npart)
                    elem = tempint(1)+(i-1)*tempint(3)
                    tempcon(elem,nen_solid+1) = npart
                enddo
            endif

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif

        else

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif   
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif
            k = 132
            do while (index('0123456789',line(k:k))==0)
                k = k-1
            enddo

            do while (index('*Nset',line(1:5))==0)
                
                j = 0
                do 
                    tempnesp(npart) = tempnesp(npart) + 1

                    if ((error_id_open(2)==0) .or. (error_id_open(4)==0) .or. &
                        (error_id_open(5)==0) .or. (error_id_open(6)==0)) then
            
                        read(line((j+1):k),*, iostat=test_error) elem
                        tempcon(elem,nen_solid+1) = npart
                    endif
                    i = index(line((j+1):k),',')
                    if (i==0) exit
                    j = j+i
                enddo
                read(file_in_solid,'(a132)',iostat=error_id1) line
                if (error_id1/=0) then
                    if (myid==0) then
                        write(*,*) 'Error reading solid.inp'
                    endif
                    !call mpi_barrier(mpi_comm_world,ierror)
                    stop
                endif
                k = 132
                do while (index('0123456789',line(k:k))==0)
                    k = k-1
                enddo

            enddo

        endif

    enddo

    nep1 = tempnesp(1)
    nep2 = tempnesp(2)
    if (myid==0) then
        write(*,'(" nep1 = ",i7)') nep1
        write(*,'(" nep2 = ",i7)') nep2
    endif

    if (ne_solid_1/=(tempnesp(1)+tempnesp(2))) then
        if (myid==0) then
            write(*,*) 'Number of elements in each part do not sum to total number of elements'
        endif
        !call mpi_barrier(mpi_comm_world,ierror)
        stop
    endif

    if ((error_id_open(2)==0) .and. (myid==0)) then
        if (elem_type_solid==1) then
            do ie=1,ne_solid_1
                write(file_mien_solid,830) tempcon(ie,1:(nen_solid+1))
                write(file_connectsolid,820) tempcon(ie,1:nen_solid)
            enddo
        elseif (elem_type_solid==2) then
            do ie=1,ne_solid_1
                write(file_mien_solid,840) tempcon(ie,1:(nen_solid+1))
                write(file_connectsolid,830) tempcon(ie,1:nen_solid)
            enddo
        elseif (elem_type_solid==3) then
            do ie=1,ne_solid_1
                write(file_mien_solid,840) tempcon(ie,1:(nen_solid+1))
                write(file_connectsolid,830) tempcon(ie,1:nen_solid)
            enddo
        elseif (elem_type_solid==4) then
            do ie=1,ne_solid_1
                write(file_mien_solid,860) tempcon(ie,1:(nen_solid+1))
                write(file_connectsolid,850) tempcon(ie,1:nen_solid)
            enddo
        endif
    elseif ((error_id_open(2)/=0) .and. (error_id_open(3)==0) .and. (myid==0)) then
        if (elem_type_solid==1) then
            do ie=1,ne_solid_1
                write(file_connectsolid,820) tempcon(ie,1:nen_solid)
            enddo
        elseif (elem_type_solid==2) then
            do ie=1,ne_solid_1
                write(file_connectsolid,830) tempcon(ie,1:nen_solid)
            enddo
        elseif (elem_type_solid==3) then
            do ie=1,ne_solid_1
                write(file_connectsolid,830) tempcon(ie,1:nen_solid)
            enddo
        elseif (elem_type_solid==4) then
            do ie=1,ne_solid_1
                write(file_connectsolid,850) tempcon(ie,1:nen_solid)
            enddo
        endif
    endif

    if (myid==0) then
        close(file_mien_solid)
        close(file_connectsolid)
    endif

    if ((error_id_open(2)==0) .or. (error_id_open(4)==0) .or. (error_id_open(5)==0) .or. (error_id_open(6)==0)) then
        if ((error_id_open(4)==0) .or. (error_id_open(5)==0) .or. (error_id_open(6)==0)) then
            allocate(connectivity(ne_solid_1,nen_solid))
            connectivity(1:ne_solid_1,1:nen_solid) = tempcon(1:ne_solid_1,1:nen_solid)
        endif
        deallocate(tempcon)
    endif


 820 format (3I8)
 830 format (4I8)
 840 format (5I8)
 850 format (8I8)
 860 format (9I8) 

    ! Write mrng.in and save pertinent boundary information

    if (myid==0) then
        if (error_id_open(4)==0) then
            write(*,*) 'Writing sbc_solid.in'
        endif
        if (error_id_open(5)==0) then
            write(*,*) 'Writing sbcnode_solid.in'
        endif
        if (error_id_open(6)==0) then
            write(*,*) 'Writing sfcon.in if needed'
        endif
    endif

    temp_nrng_solid = 0
    nnsboundary(1:nrng_solid) = 0
    nesboundary(1:nrng_solid) = 0
    allocate(nodebc(nrng_solid,nn_solid_1))
    allocate(elembc(nrng_solid,ne_solid_1))
    nodebc(1:nrng_solid,1:nn_solid_1) = 0
    elembc(1:nrng_solid,1:ne_solid_1) = 0


    allocate(uniquenodes(2*nn_solid_1))
    uniquenodes(1:(2*nn_solid_1)) = 0
    num_unique_nodes = 0
    allocate(sfcon_nodes(2*nn_solid_1))
    sfcon_nodes(1:(2*nn_solid_1)) = 0
    nn_sfcon = 0


    do while (index('*End',line(1:4))==0)

        do while (index('*Nset',line(1:5))==0)

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (myid==0) then
                write(*,*) 'I should not be in this Nset loop 2'
            endif
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif

        enddo
      
        temp_nrng_solid = temp_nrng_solid + 1
        if (nrng_solid<temp_nrng_solid) then
            if (myid==0) then
                write(*,*) 'Number of solid boundary conditions is incorrect (not enough)'
            endif
            !call mpi_barrier(mpi_comm_world,ierror)
            stop
        endif

        if (index(line,'generate')/=0) then

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif
            k = 132
            do while (index('0123456789',line(k:k))==0)
                k = k-1
            enddo

            read(line(1:k),*, iostat=test_error) tempint(1:3)
            nnsboundary(temp_nrng_solid) = 1 + (tempint(2)-tempint(1))/tempint(3)
            do i=1,nnsboundary(temp_nrng_solid)
                nodebc(temp_nrng_solid,i) = tempint(1)+(i-1)*tempint(3)
            enddo

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif

        else

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif   
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif
            k = 132
            do while (index('0123456789',line(k:k))==0)
                k = k-1
            enddo

            do while (index('*Elset',line(1:6))==0)
                
                j = 0
                do 
                    nnsboundary(temp_nrng_solid) = nnsboundary(temp_nrng_solid) + 1
                    read(line((j+1):k),*,iostat=test_error) nodebc(temp_nrng_solid,nnsboundary(temp_nrng_solid))
                    i = index(line((j+1):k),',')
                    if (i==0) exit
                    j = j+i
                enddo
                read(file_in_solid,'(a132)',iostat=error_id1) line
                if (error_id1/=0) then
                    if (myid==0) then
                        write(*,*) 'Error reading solid.inp'
                    endif
                    !call mpi_barrier(mpi_comm_world,ierror)
                    stop
                endif
                k = 132
                do while (index('0123456789',line(k:k))==0)
                    k = k-1
                enddo

            enddo

        endif


        do i=1,nnsboundary(temp_nrng_solid)
            if (num_unique_nodes==0) then
                num_unique_nodes = num_unique_nodes+1
                uniquenodes(num_unique_nodes) = nodebc(temp_nrng_solid,i)
            else
                prod_length = 1
                do j=1,num_unique_nodes
                    node = nodebc(temp_nrng_solid,i)
                    if (node==uniquenodes(j)) then
                        prod_length = 0
                    elseif (node>uniquenodes(j)) then
                        k = j+1
                    endif
                enddo
                if (prod_length==1) then
                    num_unique_nodes = num_unique_nodes+1
                    uniquenodes((k+1):num_unique_nodes) = uniquenodes(k:(num_unique_nodes-1))
                    uniquenodes(k) = node
                endif
            endif
        enddo

        


        sfcon_flag(1:3) = (bc_solid(temp_nrng_solid,2)-10000)
        sfcon_flag(1) = int(sfcon_flag(1)/100)
        sfcon_flag(2) = int((sfcon_flag(2)-sfcon_flag(1)*100)/10)
        sfcon_flag(3) = sfcon_flag(3)-sfcon_flag(1)*100-sfcon_flag(2)*10
        if ((sfcon_flag(1)==2) .or. (sfcon_flag(2)==2) .or. (sfcon_flag(3)==2)) then
            sfcon_flag(4) = 1
        else
            sfcon_flag(4) = 0
        endif
        if (sfcon_flag(4)==1) then
            do i=1,nnsboundary(temp_nrng_solid)
                if (nn_sfcon==0) then
                    nn_sfcon = nn_sfcon+1
                    sfcon_nodes(nn_sfcon) = nodebc(temp_nrng_solid,i)
                else
                    prod_length = 1
                    do j=1,nn_sfcon
                        node = nodebc(temp_nrng_solid,i)
                        if (node==sfcon_nodes(j)) then
                            prod_length = 0
                        elseif (node>sfcon_nodes(j)) then
                            k = j+1
                        endif
                    enddo
                    if (prod_length==1) then
                        nn_sfcon = nn_sfcon+1
                        sfcon_nodes((k+1):nn_sfcon) = sfcon_nodes(k:(nn_sfcon-1))
                        sfcon_nodes(k) = node
                    endif
                endif
            enddo
        endif




        do while (index('*Elset',line(1:6))==0)

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (myid==0) then
                write(*,*) 'I should not be in this Elset loop 2'
            endif
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif

        enddo

        if (index(line,'generate')/=0) then

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif
            k = 132
            do while (index('0123456789',line(k:k))==0)
                k = k-1
            enddo

            read(line(1:k),*, iostat=test_error) tempint(1:3)
            nesboundary(temp_nrng_solid) = 1 + (tempint(2)-tempint(1))/tempint(3)
            do i=1,nesboundary(temp_nrng_solid)
                elembc(temp_nrng_solid,i) = tempint(1)+(i-1)*tempint(3)
            enddo

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif

        else

            read(file_in_solid,'(a132)',iostat=error_id1) line
            if (error_id1/=0) then
                if (myid==0) then
                    write(*,*) 'Error reading solid.inp'
                endif   
                !call mpi_barrier(mpi_comm_world,ierror)
                stop
            endif
            k = 132
            do while (index('0123456789',line(k:k))==0)
                k = k-1
            enddo

            do while ((index('*Nset',line(1:5))==0) .and. (index('*End',line(1:4))==0))
                
                j = 0
                do 
                    nesboundary(temp_nrng_solid) = nesboundary(temp_nrng_solid) + 1
                    read(line((j+1):k),*,iostat=test_error) elembc(temp_nrng_solid,nesboundary(temp_nrng_solid))
                    i = index(line((j+1):k),',')
                    if (i==0) exit
                    j = j+i
                enddo
                read(file_in_solid,'(a132)',iostat=error_id1) line
                if (error_id1/=0) then
                    if (myid==0) then
                        write(*,*) 'Error reading solid.inp'
                    endif
                    !call mpi_barrier(mpi_comm_world,ierror)
                    stop
                endif
                k = 132
                do while (index('0123456789',line(k:k))==0)
                    k = k-1
                enddo

            enddo

        endif

    enddo

    if (nrng_solid>temp_nrng_solid) then
        if (myid==0) then
            write(*,*) 'Number of solid boundary conditions is incorrect (too many)'
        endif
        !call mpi_barrier(mpi_comm_world,ierror)
        stop
    endif


    if ((error_id_open(5)==0) .and. (myid==0)) then
        do i=1,num_unique_nodes
            write(file_sbcnode_solid,'(I8)') uniquenodes(i)
        enddo
    endif
    nn_sbc_1 = num_unique_nodes
    ne_sbc_1 = sum(nesboundary(1:nrng_solid))
    deallocate(uniquenodes)
    if (myid==0) then
        write(*,'(" nn_sbc_1 = ",i7)') nn_sbc_1
        write(*,'(" ne_sbc_1 = ",i7)') ne_sbc_1
    endif

    if (myid==0) then
        close(file_sbcnode_solid)
    endif
    
    
    if ((error_id_open(6)==0) .and. (myid==0)) then
        do i=1,nn_sfcon
            write(file_sfcon,'(I8)') sfcon_nodes(i)
        enddo
    endif
    node_sfcon_1 = nn_sfcon
    deallocate(sfcon_nodes)
    if (myid==0) then
        write(*,'(" node_sfcon_1 = ",I7)') node_sfcon_1
    endif

    if (myid==0) then
        close(file_sfcon)
    endif



    if ((error_id_open(4)==0) .and. (myid==0)) then

        do ibc=1,nrng_solid

            do ie=1,nesboundary(ibc)

                elem = elembc(ibc,ie)

                do inode=1,nen_solid

                    node = connectivity(elem,inode)
                    flag_node_boundary(inode) = 0
                    do i=1,nnsboundary(ibc)
                        if (node==nodebc(ibc,i)) then
                            flag_node_boundary(inode) = 1
                        endif
                    enddo
                
                enddo

                if (nen_solid==3) then
                    write(file_sbc_solid,870) elem, flag_node_boundary(1:nen_solid), bc_solid(ibc,2)
                elseif (nen_solid==4) then
                    write(file_sbc_solid,880) elem, flag_node_boundary(1:nen_solid), bc_solid(ibc,2)
                elseif (nen_solid==8) then
                    write(file_sbc_solid,890) elem, flag_node_boundary(1:nen_solid), bc_solid(ibc,2)
                endif

            enddo

        enddo

    endif

    deallocate(elembc)
    deallocate(nodebc)

 870 format (5I8)
 880 format (6I8)
 890 format (10I8) 


    if (myid==0) then
        close(file_sbc_solid)
    endif

    close(file_in_solid)
    !call mpi_barrier(mpi_comm_world,ierror)

    return
end subroutine read_abaqus_solid


                




        












