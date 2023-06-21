program distance_calculation
    use distance_module
    use testing

    implicit none
    integer :: i
    integer, dimension(6) :: natoms
    real :: rc

    rc = 3.5
    
    natoms(1) = 10
    Do i=2, 5
        natoms(i) = natoms(i-1)*10
    EndDo

    !call test_boundary_handling()

    Do i=1,6
        write(*,*) natoms(i), " atoms"
        call test_and_time_distance(natoms(i), rc)
    EndDo
    

end program distance_calculation