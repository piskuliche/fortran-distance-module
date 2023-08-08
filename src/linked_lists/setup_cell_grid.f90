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
            nbins(i) = max(1, FLOOR(box(i)/cell_length))

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
            bins_per_rank(i) = INT(CEILING(nbins(i) / REAL(nranks)))
            mpi_nbins_start(i) = max(rank * bins_per_rank(i) + 1, 1) ! Make sure we don't go below 1
            mpi_nbins_stop(i) =  min((rank + 1) * bins_per_rank(i),nbins(i)) ! Make sure we don't go above nbins
        END DO

    end subroutine setup_cell_grid