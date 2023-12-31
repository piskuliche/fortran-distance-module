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
        ! General procedure:
        ! (1) Assign ith coordinate to grid ir(nbins,nbins,nbins)
        ! (2) Assign ith element of list to the value of head for that bin
        ! -> two options here:
        ! -> -> A: No element in that cell yet, aka head(ir(1), ir(2), ir(3)) == 0
        !           This means atom i is first atom in the cell, end of the list
        ! -> -> B: head /= 0: means that i is not the first element in the cell
        !           Point list to previous element.
        ! (3) Set head of that cell to be the atom just placed.
        Do i=1, size(r,1)
            ! Find cell indices by calling the function assign_bins
            ir(:) = assign_to_grid(r(i,:), nbins, box)
            ! Build linked list
            list(i) = head(ir(1), ir(2), ir(3))
            ! Update Head List
            head(ir(1), ir(2), ir(3)) = i
        EndDo

    end subroutine build_linked_list