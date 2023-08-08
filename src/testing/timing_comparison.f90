subroutine test_timing_comparison(natoms, rc, elapsed_time)
        use MPI_f08
        use coordinates_module
        use distance_module
        use linked_lists


        implicit none

        integer, parameter :: ntimes = 1

        ! Input variables ******************************************************
        integer, intent(in) :: natoms
        real, intent (in) :: rc

        ! Output Variables *****************************************************
        real, dimension(2), intent(out) :: elapsed_time

        ! Local variables ******************************************************
        integer :: i, j, k 
        integer :: ll_count, dl_count
        integer :: nranks, rank, ierror

        real, dimension(3) :: box, rtmp
        real, dimension(natoms,3) :: r

        real :: cell_length, rc_sq

        real, allocatable :: loop_dr (:)
        integer, allocatable :: loop_id1(:), loop_id2(:)

        real, allocatable :: ll_dr(:)
        integer, allocatable :: ll_id1(:), ll_id2(:)

        real :: start_time, end_time

        logical :: load_balance
        ! ************************************************************************

        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ! Initialize random number generator
        IF (rank == 0) THEN
            CALL random_seed()
        ENDIF


        ! Initialize box side lengths
        box(1) = 103
        box(2) = 102
        box(3) = 101

        ! Initialize cell length and rc_sq
        cell_length = rc/2.0
        rc_sq = rc**2.0

        IF (rank == 0) THEN
            write(*,*) cell_length, rc, rc_sq
        END IF

        elapsed_time = 0.0
        write(*,*) "Testing timing comparison"
        do i=1, ntimes
            ! For each times there are a few steps that happen.
            ! (1) Generate coordinates
            ! (2) Initialize MPI
            ! (3) Run cell list calculation 
            ! (4) Run double loop calculation
            ! Generate Coordinates ****************************************************
            ! Initialize positions within box
            ll_count = 0
            dl_count = 0
            write(*,*) "Generating Coordinates"
            IF (rank == 0) THEN
                CALL coordinate_generator(natoms, box, r)
                CALL cpu_time(start_time)
            ENDIF
            write(*,*) "Finished Generating Coordinates"

            ! *** MPI **************************************************************
            ! (2) Intialize MPI
            CALL MPI_BCAST(r, size(r), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

            ! (3) Cell List distance calculation
            if (rank == 0) THEN
                write(*,*) "Starting Cell-List Loop"
            END IF
            CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)

            IF (rank == 0) THEN
                CALL cpu_time(end_time)
                elapsed_time(1) = elapsed_time(1) + end_time - start_time
                write(*,*) "Finished Cell-List Loop"
            ENDIF
            ! ***********************************************************************

            ! Double loop Approach **********************************************************
            IF (rank == 0) THEN
                CALL cpu_time(start_time)
            ENDIF

            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

            ! (4) Double loop distance calculation
            load_balance = .true.
            CALL double_loop_distance(r, r, box, rc_sq, loop_dr, loop_id1, loop_id2 &
                    , dl_count, same_array=.true. , cell_length=cell_length &
                    , load_balance=load_balance, verbose=.true.)

            IF (rank == 0) THEN 
                call cpu_time(end_time)
                elapsed_time(2) = elapsed_time(2) + end_time - start_time
            ENDIF
            
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

            ! WRITE FILES
            IF (rank == 0) THEN
                call CSR_INTEGER_WRITER(21, loop_id1, "loop_id1.dat")
                call CSR_INTEGER_WRITER(21, loop_id2, "loop_id2.dat")
                call CSR_REAL_WRITER(21, loop_dr, "loop_dr.dat")
                call CSR_INTEGER_WRITER(21, ll_id1, "ll_id1.dat")
                call CSR_INTEGER_WRITER(21, ll_id2, "ll_id2.dat")
                call CSR_REAL_WRITER(21, ll_dr, "ll_dr.dat")
            END IF

            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
            ! Cleanup
            IF (ALLOCATED(loop_dr)) DEALLOCATE(loop_dr)
            IF (ALLOCATED(loop_id1)) DEALLOCATE(loop_id1)
            IF (ALLOCATED(loop_id2)) DEALLOCATE(loop_id2)
            IF (ALLOCATED(ll_dr)) DEALLOCATE(ll_dr)
            IF (ALLOCATED(ll_id1)) DEALLOCATE(ll_id1)
            IF (ALLOCATED(ll_id2)) DEALLOCATE(ll_id2)
        enddo
        write(*,*) "Timing complete"
        elapsed_time = elapsed_time/ntimes


    end subroutine test_timing_comparison