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
        real, dimension(3) :: arr, rtmp
        ! *********************************************************************

        arr = 0
        Do i=1,3
            IF (coord(i) < 0) THEN
                rtmp(i) = coord(i) + box(i)*CEILING(abs(coord(i))/box(i))
            ELSE IF (coord(i) > box(i)) THEN
                rtmp(i) = coord(i) - box(i)*CEILING(abs(coord(i))/box(i))
            ELSE
                rtmp(i) = coord(i)
            END IF
            ! Assign bin index, making sure it is within the bounds of the grid
            arr(i) = min(nbins(i), max(1, CEILING(rtmp(i)/box(i)*real(nbins(i)))))
        EndDo


    end function assign_to_grid