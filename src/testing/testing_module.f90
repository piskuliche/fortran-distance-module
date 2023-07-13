module testing 
    implicit none

contains

    subroutine coordinate_generator(natoms, box, r)
        implicit none
        ! Input Arguments
        integer, intent(in) :: natoms
        real, dimension(3), intent(in) :: box
        real, dimension(natoms,3), intent(out) :: r
        ! Local Variables
        integer :: i, j
        real, dimension(3) :: rtmp

    
        do i=1, natoms
            rtmp = 0.0
            do j=1, 3
                call random_number(rtmp(j))
            enddo
            r(i,:) = rtmp(:)*box(:)
        end do

    end subroutine coordinate_generator

    subroutine compare_dr_values(ll_count, ll_dr, ll_id1, ll_id2, loop_count, loop_dr, loop_id1, loop_id2)
        real, dimension(:), intent(in) :: ll_dr, loop_dr
        integer, dimension(:), intent(in) :: ll_id1, ll_id2, loop_id1, loop_id2
        integer, intent(in) :: ll_count, loop_count
        integer :: i, j


        do i=1, ll_count
            do j=1, loop_count
                if (ll_id1(i) == loop_id1(j) .and. ll_id2(i) == loop_id2(j)) then
                    if (abs(ll_dr(i) - loop_dr(j)) > 1e-6) then
                        write(*,*) "Error: ", ll_dr(i), loop_dr(j)
                        write(*,*) "Error: ", ll_id1(i), loop_id1(j)
                        write(*,*) "Error: ", ll_id2(i), loop_id2(j)
                    endif
                else if (ll_id1(i) == loop_id2(j) .and. ll_id2(i) == loop_id1(j)) then
                    if (abs(ll_dr(i) - loop_dr(j)) > 1e-6) then
                        write(*,*) "Error: ", ll_dr(i), loop_dr(j)
                        write(*,*) "Error: ", ll_id1(i), loop_id1(j)
                        write(*,*) "Error: ", ll_id2(i), loop_id2(j)
                    endif
                endif
            enddo
        enddo
   
    end subroutine compare_dr_values

    subroutine test_timing_comparison(natoms, rc, elapsed_time)
        use distance_module
        use linked_lists

        implicit none

        integer, parameter :: ntimes = 5

        ! Input variables ******************************************************
        integer, intent(in) :: natoms
        real, intent (in) :: rc

        ! Output Variables *****************************************************
        real, dimension(2), intent(out) :: elapsed_time

        ! Local variables ******************************************************
        integer :: i, j, k 
        integer :: ll_count, dl_count

        real, dimension(3) :: box, rtmp
        real, dimension(natoms,3) :: r

        real :: cell_length, rc_sq

        real, dimension(natoms*500) :: dr_values_naive, dr_values
        integer, dimension(natoms*500) :: dr_atom1, dr_atom2, dr_atom1_naive, dr_atom2_naive
        integer, dimension(natoms*500,3) :: cell_assign_1, cell_assign_2

        real :: start_time, end_time
        ! ************************************************************************
        
        ! Initialize random number generator
        call random_seed()

        ! Initialize box side lengths
        box(1) = 103
        box(2) = 102
        box(3) = 101

        ! Initialize cell length and rc_sq
        cell_length = rc/2.0
        rc_sq = rc**.2

        elapsed_time = 0.0
        dr_values = 0.0; dr_atom1 = 0; dr_atom2 = 0
        dr_values_naive = 0.0; dr_atom1_naive = 0; dr_atom2_naive = 0
        cell_assign_1 = 0; cell_assign_2 = 0
        write(*,*) "Testing timing comparison"
        do i=1, ntimes
            ! Generate Coordinates ****************************************************
            ! Initialize positions within box
            call coordinate_generator(natoms, box, r)
            
            ! Cell-List Approach ******************************************************
            call cpu_time(start_time)

            call cell_list_distance(r, r, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2, ll_count, same_array=.true.)

            call cpu_time(end_time)
            elapsed_time(1) = elapsed_time(1) + end_time - start_time
            write(*,*) "Finished loop cell"
            ! Naive Approach **********************************************************
            call cpu_time(start_time)

            call double_loop_distance(r, r, box, rc_sq, dr_values_naive, dr_atom1_naive, dr_atom2_naive &
                            , cell_assign_1, cell_assign_2, dl_count, same_array=.true., cell_length=cell_length)

            call cpu_time(end_time)
            elapsed_time(2) = elapsed_time(2) + end_time - start_time

            call compare_dr_values(ll_count, dr_values, dr_atom1, dr_atom2, &
                         dl_count, dr_values_naive, dr_atom1_naive, dr_atom2_naive)
        enddo
        write(*,*) "Timing complete"
        elapsed_time = elapsed_time/ntimes

    end subroutine test_timing_comparison

end module testing