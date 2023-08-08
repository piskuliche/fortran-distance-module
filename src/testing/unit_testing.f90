SUBROUTINE unit_testing()

    use MPI_f08
    use coordinates_module
    use distance_module
    use linked_lists

    REAL :: rc, cell_length, rc_sq

    REAL, ALLOCATABLE :: r(:,:), ll_dr(:), loop_dr(:)
    INTEGER, ALLOCATABLE :: ll_id1(:), ll_id2(:), loop_id1(:), loop_id2(:)
    REAL, DIMENSION(3) :: box

    INTEGER :: i, ll_count, loop_count
    integer :: nranks, rank, ierror
    CHARACTER(len=40) :: filename
    LOGICAL :: load_balance


    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    ! First Unit Test - basic coordinates
    IF (rank == 0) THEN
        write(*,*) "Performing first unit test"
    END IF

    box = 30.0
    rc = 1.0
    rc_sq = rc**2.0
    cell_length = rc/2.0
    ALLOCATE(r(8,3))
    r = 0
    DO i=1, 8
        r(i,1) = REAL(i) + 10
    END DO

    load_balance = .true.
    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
    IF (rank == 0) THEN
        write(*,*) "Unit Test One"
        write(*,*) "Unit Test Result", ll_count, loop_count
        write(*,*) "Expected Result", 0
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    rc = 2.0
    rc_sq = rc**2.0
    cell_length = rc/2.0
    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
    IF (rank == 0) THEN
        write(*,*) "Unit Test Two"
        write(*,*) "Unit Test Result", ll_count, loop_count
        write(*,*) "Expected Result", 7
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    rc = 4.0
    rc_sq = rc**2.0
    cell_length = rc/2.0

    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
    IF (rank == 0) THEN
        write(*,*) "Unit Test Three"
        write(*,*) "Unit Test Result", ll_count, loop_count
        write(*,*) "Expected Result", 18
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)


    rc = 5.0
    rc_sq = rc**2.0
    cell_length = rc/2.0

    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
    IF (rank == 0) THEN
        write(*,*) "Unit Test Four"
        write(*,*) "Unit Test Result", ll_count, loop_count
        write(*,*) "Expected Result", 22
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    r = 0
    rc = 5.0
    rc_sq = rc**2.0
    cell_length = rc/2.0

    DO i=1, 4
        r(i,1) = REAL(i) 
    END DO

    DO i=5,8
        r(i,1) = REAL(i) + 10
    END DO 

     CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
    IF (rank == 0) THEN
        write(*,*) "Unit Test Five"
        write(*,*) "Unit Test Result", ll_count, loop_count
        write(*,*) "Expected Result", 12
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    r = 0
    rc = 5.0
    rc_sq = rc**2.0
    cell_length = rc/2.0

    DO i=1, 4
        r(i,1) = REAL(i) 
    END DO

    DO i=5,8
        r(i,1) = REAL(i) + 21
    END DO 

     CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
    IF (rank == 0) THEN
        write(*,*) "Unit Test Five"
        write(*,*) "Unit Test Result", ll_count, loop_count
        write(*,*) "Expected Result", 18
        DO i=1, ll_count
            write(*,*) ll_id1(i), ll_id2(i), ll_dr(i)
        END DO
        write(*,*) "*clear*"
        DO i=1, loop_count
            write(*,*) loop_id1(i), loop_id2(i), loop_dr(i)
        END DO
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)




END SUBROUTINE