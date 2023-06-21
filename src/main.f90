program distance_calculation
    use distance_module
    use testing

    implicit none
    integer :: i
    integer, dimension(12) :: natoms
    real :: rc
    real, dimension(2) :: elapsed_time

    rc = 3.5
    
    natoms(1) = 10
    Do i=2, 12
        natoms(i) = natoms(i-1)*2
    EndDo

    !call test_boundary_handling()

    Do i=1,12
        elapsed_time = 0.0
        write(*,*) natoms(i), " atoms"
        call test_timing_comparison(natoms(i), rc, elapsed_time)
        write(*,*) "Elapsed time cell-list: ", elapsed_time(1), " seconds"
        write(*,*) "Elapsed time double loop: ", elapsed_time(2), " seconds"
    EndDo
    

end program distance_calculation