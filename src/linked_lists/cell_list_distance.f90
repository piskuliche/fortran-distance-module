subroutine cell_list_distance(r1, r2, box, cell_length, rc_sq, dists, atom1, atom2 &
                    , count, same_array, include_vector, dist_components, verbosity)
    ! This subroutine calculates the distance between all pairs of atoms in a system
    ! using a cell linked list. 
    ! 
    ! Right now the MPI parallelization scheme is to build the linked list on each processor for the entire cell, but THEN
    ! to assign parts of the grid to the different MPI processes. This is not the most efficient way to do this, but it is
    ! still likely to be much better than the double loop approach for the same problem.
    !
    ! Inputs:
    !   r1: 3D array of atom positions
    !   r2: 3D array of atom positions
    !   box: 3D array of box lengths
    !   cell_length: Length of the cell
    !   rc_sq: Square of the distance cutoff
    !   same_array: Flag to indicate whether r1 and r2 are the same array (optional)

    !  
    ! Outputs:
    !   dr: Array of distances
    !   id1: Array of atom indices
    !   id2: Array of atom indices
    !
    ! Allocates:
    !  head_r1, head_r2
    !
    ! Deallocates:
    !  head_r1, head_r2
    ! 
    ! Assigns:
    !  head_r1: 3D array of the head of the linked list for r1
    !  list_r1: 1D array of the linked list for r1
    !  head_r2: 3D array of the head of the linked list for r2
    !  list_r2: 1D array of the linked list for r2
    !  nbins: 1D array of number of bins

        use mpi_f08

        implicit none
        
        ! Inputs **************************************************************
        real, intent(in) :: rc_sq, cell_length
        real, dimension(:,:), intent(in) :: r1, r2
        real, dimension(:), intent(in) :: box
        logical, intent(in), optional :: same_array
        logical, intent(in), optional :: include_vector
        integer, intent(in), optional :: verbosity
        ! Outputs *************************************************************
        integer, intent(out) :: count
        REAL, ALLOCATABLE, INTENT(out) :: dists(:)
        INTEGER, ALLOCATABLE, INTENT(out) :: atom1(:), atom2(:)
        REAL, ALLOCATABLE, INTENT(out), optional :: dist_components(:,:)
        
        ! Local Variables *****************************************************
        integer :: i, j, k, l, m, n 
        integer :: ii, jj, kk       
        integer :: di               
        integer :: ihead, jhead     
        integer :: inner_count, rank_count
        integer :: add_offset        
        real :: rsq   
        REAL, ALLOCATABLE :: dr(:)
        INTEGER, ALLOCATABLE :: id1(:), id2(:) 
        REAL, ALLOCATABLE :: dr_components(:,:)           

        integer, dimension(3) :: ir         
        real, dimension(1000) :: dr_tmp ! Temporary array for storing from cell_internal_distance
        integer, dimension(1000) :: id1_tmp, id2_tmp ! Temporary array for storing from cell_internal_distance
        real, dimension(1000, 3) :: dr_comp_tmp
        

        REAL, ALLOCATABLE :: tmpdr(:) ! Temporary array for expanding storage arrays
        INTEGER, ALLOCATABLE :: tmp1(:), tmp2(:) ! Temporary array for expanding storage arrays
        REAL, ALLOCATABLE :: tmp_components(:,:)! Temporary array for expanding storage arrays

        integer, dimension(3,500) :: map

        integer, dimension(3) :: nbins 
        integer, dimension(max_cells) :: list_r1, list_r2 
        integer, dimension(:,:,:), allocatable :: head_r1, head_r2 

        ! MPI Variables ********************************************************
        integer :: rank, nranks, ierror

        integer, dimension(3) :: mpi_nbins_start, mpi_nbins_stop
        INTEGER, ALLOCATABLE :: counts(:), displs(:)
        ! *********************************************************************

        ! MPI Communicator
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ! Allocate the MPI arrays for counts and displacements
        ALLOCATE(counts(nranks), displs(nranks))
        counts = 0; displs=0
        IF (ALLOCATED(dists)) THEN
            dists=0; atom1=0; atom2=0
            IF (present(dist_components)) THEN
                dist_components =0.0
            END IF
        END IF


        
        IF (rank == 0) THEN
            IF (PRESENT(verbosity) .and. verbosity > 1) THEN
                write(*,*) "*******************************************************"
                write(*,*) "Cell-List-Distance called with the following parameters"
                write(*,*) "cell_length: ", cell_length
                write(*,*) "rc_sq: ", rc_sq
                write(*,*) "box: ", box
                write(*,*) "*******************************************************"
            END IF
        END IF
        
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        ! Set up the grid for the distance calculation
        ! Also sets the mpi_nbins_start and stop variables that define which bins the rank is responsible for.
        call setup_cell_grid(cell_length, box, nbins, map, mpi_nbins_start, mpi_nbins_stop)

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

        ! Allocate the head arrays
        allocate(head_r1(nbins(1), nbins(2), nbins(3)))
        allocate(head_r2(nbins(1), nbins(2), nbins(3)))
        IF (rank == 0) THEN
            IF (PRESENT(verbosity) .and. verbosity > 0) THEN
                write(*,*) "nbins:", nbins(1), nbins(2), nbins(3)
            END IF
        END IF
        
        ! Allocates the rank data arrays and resizes them if necessary
        if (.NOT. ALLOCATED (dr) )  THEN
            allocate(dr(1000)); allocate(id1(1000)); allocate(id2(1000))
            allocate(dr_components(1000,3))
        END IF

        ! Initialize the head arrays
        head_r1 = 0; head_r2 = 0
        ! Initialize the Linked List arrays
        list_r1 = 0; list_r2 = 0
        ! Initialze the output arrays
        dr = 0.0; id1 = 0; id2 = 0
        dr_components = 0
        
        IF (rank == 0) THEN
            IF (PRESENT(verbosity) .and. verbosity > 1) THEN
                write(*,*) "Building linked list calculation"
            END IF
        END IF 

        ! Build cell linked list
        call build_linked_list(r1, nbins, box, head_r1, list_r1)

        ! Sets r2 to be r1 IF same_array is true, or calls linked list again. 
        IF (present(same_array) .and. (same_array .eqv. .true.)) THEN
            head_r2 = head_r1
            list_r2 = list_r1
        ELSE
            call build_linked_list(r2, nbins, box, head_r2, list_r2)
        END IF



        IF (rank == 0) THEN
            IF (PRESENT(verbosity) .and. verbosity > 1) THEN
                write(*,*) "Starting Distance Calculation"
            END IF
        END IF 

        IF (PRESENT(verbosity) .and. verbosity > 2) THEN
            DO i=1,3
                write(*,*) "rank goes from ", mpi_nbins_start(i), " to ", mpi_nbins_stop(i)
            END DO 
        END IF

    
        ! **** MPI BARRIER *******************************************************
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)


        ! Distance Calculation ****************************************************
        ! Loop over all cells using basic domain decomposition
        ! Note - there is definitely a better way of doing this for the future
        ! But for now this has to do.

        ! Set the rank_count and count to 0
        count = 0
        rank_count = 0
        
        ! Loop over bins associated with rank
        Do i=mpi_nbins_start(1), mpi_nbins_stop(1)
        Do j=mpi_nbins_start(2), mpi_nbins_stop(2)
        Do k=mpi_nbins_start(3), mpi_nbins_stop(3)
            ! Loop over nearest neighbor cells
            Do l=-2, 2
            Do m=-2, 2
            Do n=-2, 2
                ! Map the cell indices
                ii = map(1,i+l+2)
                jj = map(2,j+m+2)
                kk = map(3,k+n+2)
                ! Head of current cell
                ihead = head_r1(i,  j,  k ) 
                jhead = head_r2(ii, jj, kk)
                inner_count = 0
                
                ! Call the cell internal distance subroutine
                call cell_internal_distance(ihead, jhead, list_r1, list_r2, r1, r2, &
                        box, rc_sq, dr_tmp, id1_tmp, id2_tmp, inner_count &
                        , same_array=same_array, include_vector=include_vector, dr_components=dr_comp_tmp)

                ! If the rank_count + inner_count is greater than the size of the arrays
                IF (rank_count + inner_count > size(dr)) THEN
                    IF (PRESENT(verbosity) .and. verbosity > 2) THEN
                        WRITE(*,*) "Rank ", rank, " is reallocating arrays"
                    END IF
                    ! Note this block reallocates the arrays to bigger size if needed.
                    allocate(tmpdr(size(dr))); allocate(tmp1(size(id1))); allocate(tmp2(size(id2)))
                    allocate(tmp_components(size(dr_components,1), 3))
                    ! Save temporary arrays
                    tmp1 = id1; tmp2 = id2; tmpdr = dr
                    tmp_components = dr_components
                    ! Deallocate old arrays
                    deallocate(id2); deallocate(id1); deallocate(dr)
                    deallocate(dr_components)
                    ! Allocate to new size
                    allocate(dr(size(tmpdr)*2)); allocate(id1(size(tmp1)*2)); allocate(id2(size(tmp2)*2))
                    allocate(dr_components(size(tmp_components,1)*2,3))
                    ! Set to zero
                    dr = 0; id1 = 0; id2 = 0
                    dr_components = 0.0
                    ! Set the values
                    dr(1:size(tmpdr)) = tmpdr
                    id1(1:size(tmp1)) = tmp1
                    id2(1:size(tmp2)) = tmp2
                    dr_components(1:size(tmp_components,1),:) = tmp_components
                    ! Finally deallocate temporary arrays
                    deallocate(tmp2); deallocate(tmp1); deallocate(tmpdr)
                    deallocate(tmp_components)
                END IF

                ! Append the temporary arrays to the main arrays
                dr(rank_count+1:rank_count+inner_count) = dr_tmp(1:inner_count)

                ! If the id arrays are present, append them as well
                id1(rank_count+1:rank_count+inner_count) = id1_tmp(1:inner_count)
                id2(rank_count+1:rank_count+inner_count) = id2_tmp(1:inner_count)

                dr_components(rank_count+1:rank_count+inner_count, :) = dr_comp_tmp(1:inner_count,:)

                ! Increment the rank count by the total counts
                rank_count = rank_count + inner_count  

            EndDo !n
            EndDo !m
            EndDo !l

        EndDo !k
        EndDo !j
        EndDo !i
        IF (PRESENT(verbosity) .and. verbosity > 2) THEN
            write(*,*) "Rank", rank, " has ", rank_count, " counts"
        END IF



        ! **** MPI SECTION ********************************************************
        ! This section calculates the total counts by aggregating the information from each core
        ! It does this through the reduce argument first, and then gathers the counts for a list of 
        ! what came from each of the cores to build the displacements. It then Gathers the data into
        ! the final output array

        ! (1) Sum up the counts from each rank
        CALL MPI_ALLREDUCE(rank_count, count, 1, MPI_INTEGER, MPI_SUM,  MPI_COMM_WORLD, ierror)

        ! (1b) Allocate the output arrays based on this total count
        IF (ALLOCATED(dists)) THEN
            deallocate(dists, atom1, atom2)
            IF (present(dist_components)) THEN
                deallocate(dist_components)
            ENDIF
        END IF
        allocate(dists(count), atom1(count), atom2(count))
        IF (present(dist_components)) THEN
            allocate(dist_components(count,3))
        ENDIF
        
        ! (2) Gather the counts from each rank into a total array of counts
        CALL MPI_GATHER(rank_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        
        ! (2b) Build the displacements array based on the counts
        IF (rank == 0) THEN
            displs(1) = 0
            DO i=2, nranks
                displs(i) = sum(counts(1:i-1))
            END DO
            !write(*,*) counts(:)
            !write(*,*) displs(:)
        END IF

        ! (3) Gather the data from each rank into the final output arrays
        CALL MPI_GATHERV(dr, rank_count, MPI_REAL, dists, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_GATHERV(id1, rank_count, MPI_INTEGER, atom1, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_GATHERV(id2, rank_count, MPI_INTEGER, atom2, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        IF (rank == 0) THEN
            IF (PRESENT(verbosity) .and. verbosity > 0) THEN
                write(*,*) "Cell List Distances: ", count
                write(*,*) "********************************"
            END IF
        END IF
        ! *************************************************************************
        ! Deallocate the head arrays
        deallocate(head_r1)
        deallocate(head_r2)
        deallocate(dr, id1, id2)
        deallocate(dr_components)
        deallocate(counts, displs)
            
    end subroutine cell_list_distance