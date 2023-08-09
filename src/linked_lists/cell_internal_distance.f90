subroutine cell_internal_distance(ihead_init, jhead_init, ll_1, ll_2, r1, r2, &
                        box, rc_sq, dr, id_atom1, id_atom2, inner_count, same_array, include_vector, dr_components)
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
    real, dimension(:), intent(out) :: dr
    integer, dimension(:), intent(out) :: id_atom1, id_atom2
    REAL, dimension(:,:), intent(out), optional :: dr_components
    integer, intent(out) :: inner_count
    ! Optional Arguments
    logical, intent(in), optional :: same_array
    logical, intent(in), optional :: include_vector


    ! Local Variables
    real, dimension(4) :: dr_arr
    real :: rsq

    logical :: same_condition, diff_condition
    integer :: ihead, jhead
    
    inner_count = 0
    ihead   = ihead_init
    Do While (ihead /= 0)
        jhead  = jhead_init

        Do While (jhead /= 0) ! Might need to add somehting here
            ! Calculate the distance between the atoms
            IF (present(include_vector) .and. .NOT. include_vector) THEN
                rsq = periodic_distance2(r1(ihead,:), r2(jhead,:), box)
            ELSE
                dr_arr = periodic_distance_and_vector(r1(ihead,:), r2(jhead,:), box)
                rsq = dr_arr(4)**2.0
            END IF

            ! Set up logical conditions
            same_condition = (same_array .and. ihead < jhead)
            diff_condition = (.not. same_array)

            ! If the distance is less than the cutoff, and IF
            ! it meets one of the two conditions, THEN add it to the list
            If (rsq < rc_sq .and. (same_condition .or. diff_condition) ) Then
                inner_count = inner_count + 1
                dr(inner_count) = rsq
                id_atom1(inner_count) = ihead
                id_atom2(inner_count) = jhead
                IF (present(include_vector) .and. include_vector) THEN
                    dr_components(inner_count,:) = dr_arr(1:3)
                END IF
            End If

            ! Move to the next atom in the linked list
            jhead = ll_2(jhead)
        End Do

        ihead = ll_1(ihead)
    End Do
    end subroutine cell_internal_distance