SUBROUTINE unit_testing()

    use MPI_f08
    use coordinates_module
    use distance_module
    use linked_lists

    REAL :: rc, cell_length, rc_sq

    REAL, ALLOCATABLE :: r(:,:), ll_dr(:), loop_dr(:)
    REAL, ALLOCATABLE :: loop_drx(:), loop_dry(:), loop_drz(:)
    REAL, ALLOCATABLE :: ll_drx(:), ll_dry(:), ll_drz(:)
    INTEGER, ALLOCATABLE :: ll_id1(:), ll_id2(:), loop_id1(:), loop_id2(:)
    REAL, DIMENSION(3) :: box

    INTEGER :: i, j, ll_count, loop_count
    INTEGER :: natoms
    integer :: nranks, rank, ierror
    CHARACTER(len=40) :: filename
    LOGICAL :: load_balance
    REAL :: rnd


    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    IF (.false.) THEN
    ! First Unit Test - basic coordinates
    IF (rank == 0) THEN
        write(*,*) "Performing Unit Testing"
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
                    , load_balance=load_balance, verbosity=0)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count &
                            , same_array=.true., include_vector=.false.)
    IF (rank == 0) THEN
        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test One"
            write(*,*) "Unit Test Result", ll_count, loop_count
            write(*,*) "Expected Result", 0
        END IF
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    rc = 2.0
    rc_sq = rc**2.0
    cell_length = rc/2.0
    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbosity=0)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count&
                            , same_array=.true., include_vector=.false.)
    IF (rank == 0) THEN
        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test Two"
            write(*,*) "Unit Test Result", ll_count, loop_count
            write(*,*) "Expected Result", 7
        END IF
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    rc = 4.0
    rc_sq = rc**2.0
    cell_length = rc/2.0

    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbosity=0)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count&
                            , same_array=.true., include_vector=.false.)
    IF (rank == 0) THEN
        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test Three"
            write(*,*) "Unit Test Result", ll_count, loop_count
            write(*,*) "Expected Result", 18
        END IF
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)


    rc = 5.0
    rc_sq = rc**2.0
    cell_length = rc/2.0

    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbosity=0)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count&
                            , same_array=.true., include_vector=.false.)
    IF (rank == 0) THEN
        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test Four"
            write(*,*) "Unit Test Result", ll_count, loop_count
            write(*,*) "Expected Result", 22
        END IF
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
                    , load_balance=load_balance, verbosity=0)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count&
                            , same_array=.true., include_vector=.false., verbosity=2)
    IF (rank == 0) THEN
        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test Five"
            write(*,*) "Unit Test Result", ll_count, loop_count
            write(*,*) "Expected Result", 12
        END IF
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
                    , load_balance=load_balance, verbosity=0)
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count&
                            , same_array=.true., include_vector=.false.)
    IF (rank == 0) THEN

        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test Six"
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
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    DEALLOCATE(r)
    ALLOCATE(r(5000,3))
    ! This unit test 

        DO j=1, 5
            IF (rank == 0) THEN
                call random_number(rnd)
                box = 30 + rnd*5.0
            END IF
            CALL MPI_BCAST(box, 3, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
            r = 0
            rc = 3.0
            rc_sq = rc**2.0
            cell_length = rc/2.0
            IF (rank == 0) THEN
                CALL coordinate_generator(5000, box, r)
            END IF

            CALL MPI_BCAST(r, size(r), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
            
            loop_count = 0; ll_count = 0
            CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                            , loop_count, same_array=.true. , cell_length=cell_length &
                            , load_balance=load_balance, verbosity=0)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
            CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count&
                                    , same_array=.true., include_vector=.false., verbosity=3)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
            IF (rank == 0) THEN
                write(*,*) loop_count, ll_count
                IF (ll_count /= loop_count) THEN
                    write(*,*) "Unit Test Six, count:", j

                    write(*,*) "Box", box
                    write(*,*) "Unit Test Result", ll_count, loop_count
                    write(*,*) "Expected Result", 18
                    !DO i=1, ll_count
                    !    write(*,*) ll_id1(i), ll_id2(i), ll_dr(i)
                    !END DO
                    !write(*,*) "*clear*"
                    !DO i=1, loop_count
                    !    write(*,*) loop_id1(i), loop_id2(i), loop_dr(i)
                    !END DO
                END IF
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        END DO
    END IF

    write(*,*) "*************************************"
    write(*,*) "Final Test"
    IF (ALLOCATED(R)) THEN 
        deallocate(r)
    END IF
    if (rank == 0) THEN
        CALL read_xyz("../test_files/test.xyz", r)
        natoms = size(r,1)
    END IF
    CALL MPI_BCAST(natoms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
    IF (.NOT. allocated(r)) THEN
        allocate(r(natoms,3))
    ENDIF
    CALL MPI_BCAST(r, size(r), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

    box = 12.41380
    rc = 3.0
    rc_sq = rc**2.0
    cell_length = rc/2.0
    loop_count = 0; ll_count = 0
    write(*,*) "Test"

    CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , loop_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, include_vector=.true. &
                    , drx=loop_drx, dry=loop_dry, drz=loop_drz, verbosity=0)

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     write(*,*) "Test2"
    CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count &
                            , same_array=.true., include_vector=.true. &
                            , drx=ll_drx, dry=ll_dry, drz=ll_drz, verbosity=3)

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    IF (rank == 0) THEN
        write(*,*) loop_count, ll_count
        IF (ll_count /= loop_count) THEN
            write(*,*) "Unit Test Eight"
            write(*,*) "Unit Test Result", loop_count, ll_count
            write(*,*) "Expected Result", 22
        END IF
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    write(*,*) "*************************************"

    


END SUBROUTINE