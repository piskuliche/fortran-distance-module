module testing 

contains

    subroutine test_timing_comparison(natoms, rc, elapsed_time)
        use distance_module

        implicit none

        integer, parameter :: ntimes = 5

        ! Input variables ******************************************************
        integer, intent(in) :: natoms
        real, intent (in) :: rc

        ! Output Variables *****************************************************
        real, dimension(2), intent(out) :: elapsed_time

        ! Local variables ******************************************************
        integer :: i, j, k 

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

        do i=1, ntimes
            ! Generate Coordinates ****************************************************
            ! Initialize positions within box
            do j=1, natoms
                rtmp = 0.0
                do k=1, 3
                    call random_number(rtmp(k))
                enddo
                r(j,:) = rtmp(:)*box(:)
            end do
            ! Cell-List Approach ******************************************************
            call cpu_time(start_time)
            call cell_list_distance(r, r, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2, same_array=1)
            call cpu_time(end_time)
            elapsed_time(1) = elapsed_time(1) + end_time - start_time
            ! Naive Approach **********************************************************
            call cpu_time(start_time)
            call double_loop_distance(r, r, box, rc_sq, dr_values_naive, dr_atom1_naive, dr_atom2_naive &
                            , cell_assign_1, cell_assign_2, compare_cell=1, same_array=1, cell_length=cell_length)
            call cpu_time(end_time)
            elapsed_time(2) = elapsed_time(2) + end_time - start_time
        enddo
        elapsed_time = elapsed_time/ntimes

        !open(10, file="map_reg.dat")
        !open(11, file="map_naive.dat")
        !i = 1
        !do while ( .True. )
        !    if ( dr_values(i) == 0.0 ) exit
        !    write(10,*) dr_values(i), dr_atom1(i), dr_atom2(i), dr_values_naive(i), dr_atom1_naive(i), dr_atom2_naive(i)
        !    write(11,*) dr_atom1_naive(i), cell_assign_1(i,1), cell_assign_1(i,2), cell_assign_1(i,3)
        !    write(11,*) dr_atom2_naive(i), cell_assign_2(i,1), cell_assign_2(i,2), cell_assign_2(i,3)
        !    i = i + 1
        !enddo
        !close(10)
        !close(11)


    end subroutine test_timing_comparison

end module testing