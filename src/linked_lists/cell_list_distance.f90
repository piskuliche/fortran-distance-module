
subroutine cell_internal_distance(ihead_init, jhead_init, list_r1, list_r2, r1, r2, &
                            box, rc_sq, dr_values, dr_atom1, dr_atom2, inner_count, same_array)
    ! This subroutine calculates the distance between all pairs of atoms in a single cell
    ! TODO: Implement same_array feature.
    
    implicit none
    ! Arguments
    integer, intent(in) :: ihead, jhead, ihead_init, jhead_init
    integer, dimension(:), intent(in) :: list_r1
    real, dimension(:,:), intent(in) :: r1, r2
    real, dimension(:), intent(in) :: box
    real, intent(in) :: rc_sq
    ! Outputs
    real, dimension(1000), intent(out) :: dr
    integer, dimension(1000), intent(out) :: id_atom1, id_atom2
    integer, intent(out) :: inner_count
    ! Local Variables
    real :: rsq
    integer :: same_array
    logical :: same_condition, diff_condition
    
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
                dr_values(inner_count) = rsq
                id_atom1(inner_count) = ihead
                id_atom2(inner_count) = jhead
            End If

            ! Move to the next atom in the linked list
            jhead = list_r2(jhead)
        End Do

        ihead = list_r1(ihead)
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
    ! Sets r2 to be r1 if same_array is true
    if (present(same_array) .and. same_array == 1) then
        head_r2 = head_r1
        list_r2 = list_r1
        r2 = r1
    else
        call build_linked_list(r2, nbins, box, head_r2, list_r2)
    endif

    if (present(offset) .and. offset > 0) then
        add_offset = offset
    else
        add_offset = 0
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
    
