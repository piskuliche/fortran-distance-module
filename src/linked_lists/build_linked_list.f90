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