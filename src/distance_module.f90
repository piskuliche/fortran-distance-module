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

    function angle_between_points(r1, r2, r3, box) result(theta)
    ! This function calculates the angle between three points in a periodic system.
    !
    implicit none

    ! Inputs **************************************************************
    real, dimension(3), intent(in) :: r1, r2, r3, box
    ! Outputs *************************************************************
    real :: theta
    ! Local Variables *****************************************************
    real, dimension(4) :: dr12, dr23
    ! *********************************************************************

    dr12 = 0.0; dr23 = 0.0
    dr12 = periodic_distance_and_vector(r2, r1, box)
    dr23 = periodic_distance_and_vector(r2, r3, box)
    theta = acos(dot_product(dr12(1:3), dr23(1:3)))

    end function angle_between_points
    ! **************************************************************************



    ! **************************************************************************
    ! Subroutines **************************************************************
    ! **************************************************************************

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