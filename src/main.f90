program distance_calculation
    use mpi_f08
    use distance_module
    use testing

    implicit none
    integer :: i
    integer, dimension(20) :: natoms
    real :: rc
    real, dimension(2) :: elapsed_time

    integer :: rank, ierror, nranks
    

    
    natoms(1) = 10
    Do i=2, 20
        natoms(i) = natoms(i-1)*2
    EndDo

    ! Initialize MPI
    CALL MPI_INIT(ierror)

    ! Get the number of processors & the rank of this one
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    IF (rank == 0) OPEN(15, file="rc_3.0.dat")
    rc = 2


    elapsed_time = 0.0
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    CALL test_timing_comparison(20, rc, elapsed_time)
    IF ( rank == 0) THEN
        WRITE(*,*) 20, " atoms"
        WRITE(*,*) "Elapsed time cell-list: ", elapsed_time(1), " seconds"
        WRITE(*,*) "Elapsed time double loop: ", elapsed_time(2), " seconds"
        WRITE(15,*) 20, elapsed_time(1), elapsed_time(2)
    ENDIF

    CLOSE(15)

    CALL MPI_FINALIZE(ierror)


end program distance_calculation