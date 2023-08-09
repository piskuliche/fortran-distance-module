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
            arr(i) = min(nbins(i), max(1, CEILING(coord(i)/box(i)*nbins(i))))

        EndDo


    end function assign_to_grid