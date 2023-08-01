module linked_lists
    ! this module builds a linked list and implements its basic operations
    implicit none

    ! Parameters **************************************************************
    integer, parameter :: max_cells = 1000000
    integer, parameter :: max_bins = 500

    ! Variables ****************************************************************

contains

    function assign_to_grid(coord, nbins, box) result(arr)
    ! This function takes a set of coordinates and assigns it to the
    ! grid defined by nbins for the linked list approach. 
    ! This assumes a orthogonal unit cell for the distance calculation
    ! 
    ! Inputs:
    !   coord: 1D array of a single set of coordinates
    !   nbins: 1D array of number of bins
    !   box:   1D array of the box side lengths
    !
    ! Returns:
    !   arr: 1D array of bin assignment indices.
        implicit none

        ! Inputs **************************************************************
        integer, dimension(3), intent(in) :: nbins
        real, dimension(3), intent(in) :: coord, box
        ! Local Variables *****************************************************
        integer :: i
        real, dimension(3) :: arr
        ! *********************************************************************

        arr = 0
        Do i=1,3
            ! Assign bin index, making sure it is within the bounds of the grid
            arr(i) = min(nbins(i), max(1, int(coord(i)/box(i)*nbins(i))))
        EndDo

    end function assign_to_grid

    subroutine setup_cell_grid(cell_length, box, nbins, map, mpi_nbins_start, mpi_nbins_stop)
    ! This subroutine sets up the grid for the linked list approach
    !
    ! Inputs:
    !   cell_length: Length of the cell 
    !   box: 1D array of box lengths
    !
    ! Outputs:
    !   nbins: 1D array of number of bins
    !   map: 2D array of bin indices
    !   bins_per_rank: 1D array of number of bins per rank
    !   mpi_nbins_start: 1D array of starting bin indices for each rank
    !   mpi_nbins_stop: 1D array of ending bin indices for each rank

        use mpi_f08
        implicit none
        
        ! Inputs **************************************************************
        real, intent(in) :: cell_length
        real, dimension(:), intent(in) :: box
        ! Outputs *************************************************************
        integer, dimension(:), intent(out) :: nbins
        integer, dimension(:,:), intent(out) :: map
        ! See the MPI OUTPUT below
        ! Local Variables *****************************************************
        integer :: i, j

        ! MPI Variables
        integer :: rank, nranks, ierror
        integer, dimension(3) :: bins_per_rank
        integer, dimension(3), intent(out) :: mpi_nbins_start, mpi_nbins_stop
        ! *********************************************************************

        ! MPI Communiator
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ! Reset the map and binsize
        map = 0; nbins = 0

        
        do i=1,3
            ! Calculate the number of bins in each direction
            nbins(i) = max(1, int(box(i)/cell_length))

            ! Check to make sure the number of bins is not too large
            IF ( nbins(i) > max_bins ) THEN
                write(*,*) "Error: Number of bins in direction ", i, " is greater than ", max_bins
                stop "Error: Number of bins is too large - increase max_bins"
            END IF

            ! Set up the map array to handle PBCs
            ! Note - assumes that offset is 2 in this case
            ! For five bins, this looks something like:
            ! i:          1 2 3 4 5 6 7 8 9
            ! map(i,:):   4 5 1 2 3 4 5 1 2
            map(i,1) = nbins(i)-1
            map(i,2) = nbins(i)
            do j=1, nbins(i)
                map(i,j+2) = j
            enddo
            map(i,nbins(i)+3) = 1
            map(i,nbins(i)+4) = 2
        enddo

        ! MPI -- assign bins for main loop based on the number of ranks
        ! This is not necessarily load balanced - if you have a high density region in one cell, versus 
        ! a low density region in another cell you will have a load imbalance.
        bins_per_rank = 0
        mpi_nbins_start = 0
        mpi_nbins_stop = 0
        DO i=1, 3
            bins_per_rank(i) = nbins(i) / nranks
            mpi_nbins_start(i) = max(rank * bins_per_rank(i) + 1, 1) ! Make sure we don't go below 1
            mpi_nbins_stop(i) =  min((rank + 1) * bins_per_rank(i),nbins(i)) ! Make sure we don't go above nbins
            write(*,*) "Rank ", rank, " has bins ", mpi_nbins_start(i), " to ", mpi_nbins_stop(i), " for dimension ", i
        END DO

    end subroutine setup_cell_grid

    subroutine build_linked_list(r, nbins, box, head, list)
    ! Subroutine that builds the linked list and the head.
    !
    ! Inputs:
    !  r: 2D array of atom positions
    !  nbins: 1D array of number of bins
    !  box:   1D array of the box side lengths
    ! 
    ! Outputs:
    !  head: 3D array of the head of the linked list
    !  list: 1D array of the linked list

        implicit none

        ! Inputs **************************************************************
        integer, dimension(3), intent(in) :: nbins
        real, dimension(3), intent(in) :: box
        real, dimension(:,:), intent(in) :: r

        ! Outputs *************************************************************
        integer, dimension(:,:,:), intent(out) :: head
        integer, dimension(:), intent(out) :: list

        ! Local Variables *****************************************************
        integer :: i, j
        integer, dimension(3) :: ir
        ! *********************************************************************

        ! Initialize the arrays to zero
        head=0; list=0

        ! Loop over all atoms in r
        Do i=1, size(r,1)
            ! Find cell indices by calling the function assign_bins
            ir(:) = assign_to_grid(r(i,:), nbins, box)
            ! Build linked list
            list(i) = head(ir(1), ir(2), ir(3))
            ! Update Head List
            head(ir(1), ir(2), ir(3)) = i
        EndDo

    end subroutine build_linked_list

    subroutine cell_internal_distance(ihead_init, jhead_init, ll_1, ll_2, r1, r2, &
                        box, rc_sq, dr, id_atom1, id_atom2, inner_count, same_array)
    ! This subroutine calculates the distance between all pairs of atoms in a single cell
    ! TODO: Implement same_array feature.
    
    use distance_module

    implicit none
    ! Arguments
    integer, intent(in) :: ihead_init, jhead_init
    integer, dimension(:), intent(in) :: ll_1, ll_2
    real, dimension(:,:), intent(in) :: r1, r2
    real, dimension(:), intent(in) :: box
    real, intent(in) :: rc_sq
    ! Outputs
    real, dimension(:), intent(out) :: dr
    integer, dimension(:), intent(out) :: id_atom1, id_atom2
    integer, intent(out) :: inner_count
    ! Optional Arguments
    logical, intent(in), optional :: same_array

    ! Local Variables
    real :: rsq

    logical :: same_condition, diff_condition
    integer :: ihead, jhead
    
    inner_count = 0

    ihead   = ihead_init
    Do While (ihead /= 0)
        jhead   = jhead_init

        Do While (jhead /= 0) ! Might need to add somehting here
            ! Calculate the distance between the atoms
            rsq = periodic_distance2(r1(ihead,:), r2(jhead,:), box)

            ! Set up logical conditions
            same_condition = (same_array .and. ihead < jhead)
            diff_condition = (.not. same_array)

            ! If the distance is less than the cutoff, and IF
            ! it meets one of the two conditions, THEN add it to the list
            If (rsq < rc_sq .and. (same_condition .or. diff_condition) ) Then
                inner_count = inner_count + 1
                dr(inner_count) = rsq
                id_atom1(inner_count) = ihead
                id_atom2(inner_count) = jhead
            End If

            ! Move to the next atom in the linked list
            jhead = ll_2(jhead)
        End Do

        ihead = ll_1(ihead)
    End Do

    end subroutine cell_internal_distance

    subroutine cell_list_distance(r1, r2, box, cell_length, rc_sq, dists, atom1, atom2, count, same_array, offset)
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
        integer, intent(in), optional :: offset
        ! Outputs *************************************************************
        integer, intent(out) :: count
        REAL, ALLOCATABLE, INTENT(out) :: dists(:)
        INTEGER, ALLOCATABLE, INTENT(out) :: atom1(:), atom2(:)
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

        integer, dimension(3) :: ir         
        real, dimension(1000) :: dr_tmp ! Temporary array for storing from cell_internal_distance
        integer, dimension(1000) :: id1_tmp, id2_tmp ! Temporary array for storing from cell_internal_distance

        REAL, ALLOCATABLE :: tmpdr(:) ! Temporary array for expanding storage arrays
        INTEGER, ALLOCATABLE :: tmp1(:), tmp2(:) ! Temporary array for expanding storage arrays

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


        ! Set up the grid for the distance calculation
        ! Also sets the mpi_nbins_start and stop variables that define which bins the rank is responsible for.
        call setup_cell_grid(cell_length, box, nbins, map, mpi_nbins_start, mpi_nbins_stop)

        ! Allocate the head arrays
        allocate(head_r1(nbins(1), nbins(2), nbins(3)))
        allocate(head_r2(nbins(1), nbins(2), nbins(3)))
        
        ! Allocates the rank data arrays and resizes them if necessary
        if (.NOT. ALLOCATED (dr) )  THEN
            allocate(dr(1000)); allocate(id1(1000)); allocate(id2(1000))
        END IF

        ! Initialize the head arrays
        head_r1 = 0; head_r2 = 0
        ! Initialize the Linked List arrays
        list_r1 = 0; list_r2 = 0
        ! Initialze the output arrays
        dr = 0.0; id1 = 0; id2 = 0
        
        IF (rank == 0) THEN
            write(*,*) "Building linked list calculation"
        END IF 
        ! Build cell linked list
        call build_linked_list(r1, nbins, box, head_r1, list_r1)
        ! Sets r2 to be r1 IF same_array is true, or calls linked list again. 
        IF (present(same_array) .and. same_array .eqv. .true.) THEN
            head_r2 = head_r1
            list_r2 = list_r1
        ELSE
            call build_linked_list(r2, nbins, box, head_r2, list_r2)
        END IF



        IF (rank == 0) THEN
            write(*,*) "Starting Distance Calculation"
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
                ihead = head_r1(i,j,k) 
                jhead = head_r2(ii, jj, kk)
                inner_count = 0
                
                ! Call the cell internal distance subroutine
                call cell_internal_distance(ihead, jhead, list_r1, list_r2, r1, r2, &
                        box, rc_sq, dr_tmp, id1_tmp, id2_tmp, inner_count, same_array)

                ! If the rank_count + inner_count is greater than the size of the arrays
                IF (rank_count + inner_count > size(dr)) THEN
                    WRITE(*,*) "Rank ", rank, " is reallocating arrays"
                    ! Note this block reallocates the arrays to bigger size if needed.
                    allocate(tmpdr(size(dr))); allocate(tmp1(size(id1))); allocate(tmp2(size(id2)))
                    tmp1 = id1; tmp2 = id2; tmpdr = dr
                    deallocate(id2); deallocate(id1); deallocate(dr)
                    allocate(dr(size(tmpdr)*2)); allocate(id1(size(tmp1)*2)); allocate(id2(size(tmp2)*2))
                    dr = 0; id1 = 0; id2 = 0
                    dr(1:size(tmpdr)) = tmpdr
                    id1(1:size(tmp1)) = tmp1
                    id2(1:size(tmp2)) = tmp2
                    deallocate(tmp2); deallocate(tmp1); deallocate(tmpdr)
                END IF

                ! Append the temporary arrays to the main arrays
                dr(rank_count+1:rank_count+inner_count) = dr_tmp(1:inner_count)

                ! If the id arrays are present, append them as well
                id1(rank_count+1:rank_count+inner_count) = id1_tmp(1:inner_count)
                id2(rank_count+1:rank_count+inner_count) = id2_tmp(1:inner_count)

                ! Increment the rank count by the total counts
                rank_count = rank_count + inner_count  

            EndDo !n
            EndDo !m
            EndDo !l
        EndDo !k
        EndDo !j
        EndDo !i



        ! **** MPI SECTION ********************************************************
        ! This section calculates the total counts by aggregating the information from each core
        ! It does this through the reduce argument first, and then gathers the counts for a list of 
        ! what came from each of the cores to build the displacements. It then Gathers the data into
        ! the final output array

        ! (1) Sum up the counts from each rank
        CALL MPI_REDUCE(rank_count, count, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        ! (1b) Allocate the output arrays based on this total count
        allocate(dists(count), atom1(count), atom2(count))
        
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

        IF (rank == 0) THEN
            write(*,*) "count", count
        END IF
        ! *************************************************************************
        ! Deallocate the head arrays
        deallocate(head_r1)
        deallocate(head_r2)
        deallocate(dr, id1, id2)
        deallocate(counts, displs)
            
    end subroutine cell_list_distance
        


end module linked_lists
