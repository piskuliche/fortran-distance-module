subroutine setup_cell_grid(cell_length, box, nbins, map)
    ! This subroutine sets up the grid for the linked list approach
    !
    ! Inputs:
    !   cell_length: Length of the cell 
    !   box: 1D array of box lengths
    !
    ! Outputs:
    !   nbins: 1D array of number of bins
    !   map: 2D array of bin indices
    !

        implicit none
        
        ! Inputs **************************************************************
        real, intent(in) :: cell_length
        real, dimension(:), intent(in) :: box
        ! Outputs *************************************************************
        integer, dimension(:), intent(out) :: nbins
        integer, dimension(:,:), intent(out) :: map
        ! Local Variables *****************************************************
        integer :: i, j
        ! *********************************************************************

        ! Reset the map and binsize
        map = 0; nbins = 0

        
        do i=1,3
            ! Calculate the number of bins in each direction
            nbins(i) = max(1, int(box(i)/cell_length))

            ! Check to make sure the number of bins is not too large
            if ( nbins(i) > max_bins ) then
                write(*,*) "Error: Number of bins in direction ", i, " is greater than ", max_bins
                stop "Error: Number of bins is too large - increase max_bins"
            endif

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

    end subroutine setup_cell_grid