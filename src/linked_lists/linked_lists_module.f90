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
    real, dimension(1000), intent(out) :: dr
    integer, dimension(1000), intent(out) :: id_atom1, id_atom2
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

            ! If the distance is less than the cutoff, and if
            ! it meets one of the two conditions, then add it to the list
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
        logical, intent(in), optional :: same_array
        integer, intent(in), optional :: offset
        ! Outputs *************************************************************
        real, dimension(:), intent(out) :: dr_values
        integer, dimension(:), intent(out) :: dr_atom1, dr_atom2
        ! Local Variables *****************************************************
        integer :: i, j, k, l, m, n 
        integer :: ii, jj, kk       
        integer :: di               
        integer :: ihead, jhead     
        integer :: count, inner_count
        integer :: add_offset        
        real :: rsq               

        integer, dimension(3) :: ir         
        real, dimension(3) :: dr_tmp        
        integer, dimension(3,500) :: map

        integer, dimension(3) :: nbins 
        integer, dimension(max_cells) :: list_r1, list_r2 
        integer, dimension(:,:,:), allocatable :: head_r1, head_r2 
        ! *********************************************************************

        ! Set up the grid for the distance calculation
        call setup_cell_grid(cell_length, box, nbins, map)

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
        ! Sets r2 to be r1 if same_array is true
        if (present(same_array) .and. same_array .eqv. .true.) then
            head_r2 = head_r1
            list_r2 = list_r1
        else
            call build_linked_list(r2, nbins, box, head_r2, list_r2)
        endif

        if (present(offset)) then
            if (offset > 0) then
                add_offset = offset
            else
                add_offset = 0
            endif
        endif

        count = 0 ! Set Count to Zero
        ! Distance Calculation ****************************************************
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
                jhead = head_r2(ii, jj, kk)
                inner_count = 0
                call cell_internal_distance(ihead, jhead, list_r1, list_r2, r1, r2, &
                        box, rc_sq, dr_values, dr_atom1, dr_atom2, inner_count, same_array)
                count = count + inner_count   
            EndDo !n
            EndDo !m
            EndDo !l
        EndDo !k
        EndDo !j
        EndDo !i
        write(*,*) "count", count
        ! *************************************************************************
        ! Deallocate the head arrays
        deallocate(head_r1)
        deallocate(head_r2)
            
    end subroutine cell_list_distance
        


end module linked_lists
