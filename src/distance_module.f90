module distance_module
    ! This module contains the functions and subroutines for calculating the distance between
    ! all pairs of atoms in a system (within a cutoff) using a cell linked list.
    !


    implicit none

    ! Parameters **************************************************************
    integer, parameter :: max_cells = 1000000
    integer, parameter :: max_bins = 500

    ! Variables ****************************************************************
    integer, dimension(3) :: nbins 
    integer, dimension(max_cells) :: list_r1, list_r2 
    integer, dimension(:,:,:), allocatable :: head_r1, head_r2 

contains

    ! **************************************************************************
    ! Functions ****************************************************************
    ! **************************************************************************

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
    ! 
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

    function periodic_distance2(r1, r2, box) result(dr_sq)
    ! This function calculates the square of the distance between two atoms
    ! in a periodic system.
    !
    ! Inputs:
    !   r1: 1D array of atom positions
    !   r2: 1D array of atom positions
    !   box: 1D array of box lengths
    !
    ! Returns:
    !   dr_sq: Square of the distance between the two atoms
    !
        implicit none

        ! Inputs **************************************************************
        real, dimension(3), intent(in) :: r1, r2, box
        ! Outputs *************************************************************
        real :: dr_sq
        ! Local Variables *****************************************************
        integer :: i
        real, dimension(3) :: dr_tmp
        ! *********************************************************************

        dr_sq = 0.0; dr_tmp = 0.0
        Do i=1,3
            ! Calculate the distance between the two atoms
            dr_tmp(i) = r1(i) - r2(i)
            ! Apply periodic boundary conditions
            dr_tmp(i) = dr_tmp(i) - box(i)*nint(dr_tmp(i)/box(i))
        EndDo
        dr_sq = sum(dr_tmp**2)

    end function periodic_distance2

    function periodic_distance_and_vector(r1, r2, box) result(dr_arr)
    ! This function calculates the square of the distance between two atoms
    ! in a periodic system.
    !
    ! Inputs:
    !   r1: 1D array of atom positions
    !   r2: 1D array of atom positions
    !   box: 1D array of box lengths
    !
    ! Returns:
    !   dr_arr: 1D array that has first three elements as the distance vector
    !           and the last element is the distance.
    !
        implicit none

        ! Inputs **************************************************************
        real, dimension(3), intent(in) :: r1, r2, box
        ! Outputs *************************************************************
        real, dimension(4) :: dr_arr
        ! Local Variables *****************************************************
        integer :: i
        real :: dr
        
        real, dimension(3) :: dr_tmp
        ! *********************************************************************

        dr = 0.0; dr_arr = 0.0; dr_tmp = 0.0
        Do i=1,3
            ! Calculate the distance between the two atoms
            dr_tmp(i) = r1(i) - r2(i)
            ! Apply periodic boundary conditions
            dr_arr(i) = dr_tmp(i) - box(i)*nint(dr_tmp(i)/box(i))
        EndDo
        ! Calculate distance
        dr = sqrt(sum(dr_tmp**2))
        ! Assign values to dr_arr and get unit vector
        dr_arr(1:3) = dr_tmp(1:3)/dr
        ! Assign distance to dr_arr
        dr_arr(4) = dr

    end function periodic_distance_and_vector
    ! **************************************************************************



    ! **************************************************************************
    ! Subroutines **************************************************************
    ! **************************************************************************

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

    subroutine setup_grid(cell_length, box, nbins, map)
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

    end subroutine setup_grid

    subroutine cell_list_distance(r1, r2, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2, same_array, offset)
    ! This subroutine calculates the distance between all pairs of atoms in a system
    ! using a cell linked list. 
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
    !   dr_values: Array of distances
    !   dr_atom1: Array of atom indices
    !   dr_atom2: Array of atom indices
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


        implicit none
        
        ! Inputs **************************************************************
        real, intent(in) :: rc_sq, cell_length
        real, dimension(:,:), intent(in) :: r1, r2
        real, dimension(:), intent(in) :: box
        integer, intent(in), optional :: same_array
        integer, intent(in), optional :: offset
        ! Outputs *************************************************************
        real, dimension(:), intent(out) :: dr_values
        integer, dimension(:), intent(out) :: dr_atom1, dr_atom2
        ! Local Variables *****************************************************
        integer :: i, j, k, l, m, n 
        integer :: ii, jj, kk       
        integer :: di               
        integer :: ihead, jhead     
        integer :: count
        integer :: add_offset        
        real :: rsq               

        integer, dimension(3) :: ir         
        real, dimension(3) :: dr_tmp        
        integer, dimension(3,500) :: map   
        ! *********************************************************************

        ! Set up the grid for the distance calculation
        call setup_grid(cell_length, box, nbins, map)

        ! Allocate the head arrays
        allocate(head_r1(nbins(1), nbins(2), nbins(3)))
        allocate(head_r2(nbins(1), nbins(2), nbins(3)))

        ! Initialize the head arrays
        head_r1 = 0; head_r2 = 0
        ! Initialize the Linked List arrays
        list_r1 = 0; list_r2 = 0
        ! Initialze the output arrays
        dr_values = 0.0; dr_atom1 = 0; dr_atom2 = 0
        
        ! Build cell linked list
        call build_linked_list(r1, nbins, box, head_r1, list_r1)
        call build_linked_list(r2, nbins, box, head_r2, list_r2)

        if (present(offset) .and. offset > 0) then
            add_offset = offset
        else
            add_offset = 0
        endif

        count = 0 ! Set Count to Zero
        ! Distance Calculation ****************************************************
        if ( present(same_array) .and. same_array == 1) then
            ! Array is the Same ***************************************************
            Do i=1, nbins(1)
            Do j=1, nbins(2)
            Do k=1, nbins(3)
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
                    ! Loop over atoms in the current cell
                    Do While (ihead /= 0)
                        ! Loop over atoms in the neighboring cell
                        jhead = head_r2(ii, jj, kk)

                        Do While (jhead /= 0 .and. (jhead /= ihead + add_offset))
                            ! Calculate the distance between the atoms
                            rsq = periodic_distance2(r1(ihead,:), r2(jhead,:), box)

                            ! Distance cutoff
                            If (rsq < rc_sq .and. (ihead+add_offset < jhead)) Then
                                count = count + 1
                                dr_values(count) = rsq
                                dr_atom1(count) = ihead+add_offset
                                dr_atom2(count) = jhead
                            EndIf

                            ! Move to the next atom in neighboring cell
                            jhead = list_r2(jhead)
                        EndDo
                        ! Move to nxt atom in current cell
                        ihead = list_r1(ihead)
                    EndDo
                EndDo !n
                EndDo !m
                EndDo !l
            EndDo !k
            EndDo !j
            EndDo !i
        Else
            ! Array is Different ****************************************************
            Do i=1, nbins(1)
            Do j=1, nbins(2)
            Do k=1, nbins(3)
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
                    ! Loop over atoms in the current cell
                    Do While (ihead /= 0)
                        ! Loop over atoms in the neighboring cell
                        jhead = head_r2(ii, jj, kk)

                        Do While (jhead /= 0)
                            ! Calculate the distance between the atoms
                            rsq = periodic_distance2(r1(ihead,:), r2(jhead,:), box)

                            ! Apply distance cutoff
                            If (rsq < rc_sq .and. rsq > 0) Then
                                count = count + 1
                                dr_values(count) = rsq
                                dr_atom1(count) = ihead
                                dr_atom2(count) = jhead
                            EndIf

                            ! Move to the next atom in neighboring cell
                            jhead = list_r2(jhead)
                        EndDo
                        ! Move to nxt atom in current cell
                        ihead = list_r1(ihead)
                    EndDo
                EndDo !n
                EndDo !m
                EndDo !l
            EndDo !k
            EndDo !j
            EndDo !i
            ! *********************************************************************
        EndIf
        write(*,*) "count", count
        ! *************************************************************************
        ! Deallocate the head arrays
        deallocate(head_r1)
        deallocate(head_r2)
        
    end subroutine cell_list_distance

    subroutine double_loop_distance(r1, r2, box, rc_sq &
        , dr_values, dr_atom1, dr_atom2, cell_assign_1, cell_assign_2, compare_cell, same_array, cell_length)
    ! This subroutine calculates the distance between all pairs of atoms in a system
    ! using a double loop.
    !
    ! Inputs:
    !   r1: 3D array of atom positions
    !   r2: 3D array of atom positions
    !   box: 3D array of box lengths
    !   rc_sq: Square of the distance cutoff
    !   same_array: Flag to indicate whether r1 and r2 are the same array (optional)
    
    ! Outputs:
    !   dr_values: Array of distances
    !   dr_atom1: Array of atom indices
    !   dr_atom2: Array of atom indices
    !   cell_assign_1: Array of cell indices - for comparing with cell list method
    !   cell_assign_2: Array of cell indices - for comparing with cell list method


        implicit none

        ! Inputs **************************************************************
        real, intent(in) :: rc_sq
        real, dimension(:,:), intent(in) :: r1, r2 ! Coordinates
        real, dimension(:), intent(in) :: box ! Box dimensions
        integer, intent(in), optional :: compare_cell ! Flag to compare with cell list method
        integer, intent(in), optional :: same_array   ! Flag to compare with same array (default 0)
        real, intent(in), optional :: cell_length
        ! Outputs *************************************************************
        real, dimension(:), intent(out) :: dr_values
        integer, dimension(:), intent(out) :: dr_atom1, dr_atom2
        integer, dimension(:,:), intent(out) :: cell_assign_1, cell_assign_2
        ! Local Variables *****************************************************
        integer :: i, j, di ! Loop Indices
        integer :: count    ! Number of pairs found
        real :: rsq         ! Temporary distance squared


        real, dimension(3) :: dr_tmp        ! Temporary distance vector
        integer, dimension(3,500) :: map    ! Map of bin indices
        ! ************************************************************************

        ! Set up the grid for the distance calculation if compare_cell is set
        if (present(compare_cell) .and. compare_cell == 1) then
            call setup_grid(cell_length, box, nbins, map)
        end if

        ! Initialize the sparse matrix arrays
        dr_values = 0.0; dr_atom1 = 0; dr_atom2 = 0
        ! Initialize the cell_assign arrays
        cell_assign_1 = 0; cell_assign_2 = 0

        count = 0
        ! Distance Calculation ****************************************************
        if (present(same_array) .and. same_array == 1) then
            ! Array is the Same ***************************************************
            Do i=1, size(r1,1)
                Do j=i+1, size(r2,1)
                    ! Periodic distance calculation
                    rsq = periodic_distance2(r1(i,:), r2(j,:), box)
                    ! Distance cutoff & store non-zero elements
                    If (rsq < rc_sq .and. rsq > 0) Then
                        count = count + 1
                        dr_values(count) = rsq
                        dr_atom1(count) = i
                        dr_atom2(count) = j
                        ! Assign the cell indices if compare_cell == 1
                        If (present(compare_cell) .and. compare_cell == 1) then
                            cell_assign_1(count,:) = assign_to_grid(r1(i,:), nbins, box)
                            cell_assign_2(count,:) = assign_to_grid(r2(j,:), nbins, box)
                        EndIf
                    EndIf
                EndDo
            EndDo
            ! **********************************************************************
        else
            ! Array is Different ***************************************************
            Do i=1, size(r1,1)
                Do j=1, size(r2,1)
                    ! Periodic distance calculation
                    rsq = periodic_distance2(r1(i,:), r2(j,:), box)
                    ! Distance cutoff & store non-zero elements
                    If (rsq < rc_sq) Then
                        count = count + 1
                        dr_values(count) = rsq
                        dr_atom1(count) = i
                        dr_atom2(count) = j
                        ! Assign the cell indices if compare_cell == 1
                        If (present(compare_cell) .and. compare_cell == 1) then
                            cell_assign_1(count,:) = assign_to_grid(r1(i,:), nbins, box)
                            cell_assign_2(count,:) = assign_to_grid(r2(j,:), nbins, box)
                        EndIf
                    EndIf
                EndDo
            EndDo
            ! **********************************************************************
        EndIf
        ! **************************************************************************
        write(*,*) "dcount", count
    end subroutine double_loop_distance

end module distance_module