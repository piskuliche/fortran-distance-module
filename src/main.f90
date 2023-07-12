program distance_calculation
    use distance_module
    use testing

    implicit none
    integer :: i
    integer, dimension(20) :: natoms
    real :: rc
    real, dimension(2) :: elapsed_time

 
    
    natoms(1) = 10
    Do i=2, 20
        natoms(i) = natoms(i-1)*2
    EndDo

    !call test_boundary_handling()

    open(15, file="rc_3.0.dat")
    rc = 3.0
    Do i=1,20
        elapsed_time = 0.0
        write(*,*) natoms(i), " atoms"
        call test_timing_comparison(natoms(i), rc, elapsed_time)
        write(*,*) "Elapsed time cell-list: ", elapsed_time(1), " seconds"
        write(*,*) "Elapsed time double loop: ", elapsed_time(2), " seconds"
        write(15,*) natoms(i), elapsed_time(1), elapsed_time(2)
    EndDo
    close(15)


end program distance_calculation