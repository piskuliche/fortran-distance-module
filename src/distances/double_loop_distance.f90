subroutine double_loop_distance(r1, r2, box, rc_sq &
        , dists, atom1, atom2, count, same_array, cell_length, load_balance &
        , verbose)
    ! This subroutine calculates the distance between all pairs of atoms in a system
    ! using a double loop.
    !
    ! Inputs:
    !   r1: 3D array of atom positions
    !   r2: 3D array of atom positions
    !   box: 3D array of box lengths
    !   rc_sq: Square of the distance cutoff
    !   same_array: Flag to indicate whether r1 and r2 are the same array (optional)
    
    ! Outputs:
    !   dr_values: Array of distances
    !   dr_atom1: Array of atom indices
    !   dr_atom2: Array of atom indices
    !   cell_assign_1: Array of cell indices - for comparing with cell list method
    !   cell_assign_2: Array of cell indices - for comparing with cell list method

        USE mpi_f08

        IMPLICIT NONE

        ! Inputs **************************************************************
        REAL, INTENT(in) :: rc_sq
        REAL, DIMENSION(:,:), INTENT(in) :: r1, r2 ! Coordinates
        REAL, DIMENSION(:), INTENT(in) :: box ! Box dimensions
        logical, INTENT(in), optional :: same_array   ! Flag to compare with same array (default 0)
        REAL, INTENT(in), optional :: cell_length
        LOGICAL, INTENT(inout), optional :: load_balance
        LOGICAL, INTENT(in), optional :: verbose
        ! Outputs *************************************************************
        REAL, ALLOCATABLE, INTENT(out) :: dists(:)
        INTEGER, ALLOCATABLE, INTENT(out) :: atom1(:), atom2(:)
        INTEGER, INTENT(out) :: count    ! Number of pairs found
        ! Local Variables *****************************************************
        INTEGER :: i, j, di ! Loop Indices
        REAL :: rsq         ! Temporary distance squared
        INTEGER :: jstart
        INTEGER :: ierror, nranks, rank
        INTEGER :: istart, istop,  ncalc
        INTEGER :: dr_count

        INTEGER(kind=8) :: nper, ncoords


        REAL, DIMENSION(3) :: dr_tmp        ! Temporary distance vector
        INTEGER, DIMENSION(3,500) :: map    ! Map of bin indices

        REAL, ALLOCATABLE :: dr(:), tmpdr(:) ! TEMPORARY Distance
        INTEGER, ALLOCATABLE :: id1(:), id2(:), tmp1(:), tmp2(:)
        INTEGER, ALLOCATABLE :: counts(:), displs(:), ifinishes(:), istarts(:)

        INTEGER :: ncalcs, iprev
        ! ************************************************************************
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ALLOCATE(counts(nranks))
        ALLOCATE(displs(nranks))


        ! Allocate initial block for arrays
        IF ( .NOT. ALLOCATED(dr) ) THEN
            allocate(dr(1000)); allocate(id1(1000)); allocate(id2(1000))
        ENDIF
        
        ALLOCATE(ifinishes(nranks))
        ALLOCATE(istarts(nranks))


        ! Get the number of coordinates of r1
        ncoords = size(r1, 1)
        ! Calculate the number of atoms per proc
        nper = CEILING(REAL(ncoords/nranks))
        ! Only load balance if r1 and r2 are the same
        IF ( present(same_array) .and. same_array .eqv. .false.) load_balance=.false.


        ! Calculate the start and stop indices for each rank
        ! ********** Load balancing ************************************************
        ! This does load balancing for the double loop by assuming that 
        ! if r1 and r2 are the same, then only the upper traingle is needed.
        ! For instance, take the following matrix
        !       (1,1) (2,1) (3,1) (4,1)
        !       (1,2) (2,2) (3,2) (4,2)
        !       (1,3) (2,3) (3,3) (4,3)
        !       (1,4) (2,4) (3,4) (4,4)
        !
        !  Now: note that if r1 and r2 are the same, then (1,1) (2,2) are zero, and (2,1) = (1,2)
        !  So, we only need to calculate the upper triangle.
        ! 
        if (load_balance) THEN
            if (rank == 0) THEN
                write(*,*) "Load balancing has been selected"
            write(*,*) nper
            END IF 
            !write(*,*) ncoords*(ncoords-1), ncoords*(ncoords-1)/(2*nranks)
            ! Calculate the number of calculations per rank
            nper = ncoords*(ncoords-1)/(2*nranks)
            ! Write verbose info to the screen
            IF (verbose) THEN 
                write(*,*) "rank load", ncoords, nranks, nper
            END IF
            DO i=1, nranks
                IF (i == 1) THEN
                    iprev = 1
                ELSE
                    iprev = ifinishes(i-1)+1
                END IF
                !write(*,*) rank, iprev, FLOOR(1 + SQRT(1.0 + 4.0*(2*nper + iprev*(iprev-1)))/2.0), nper
                istarts(i) = iprev
                ifinishes(i) = FLOOR(1 + SQRT(1.0 + 4.0*(2*nper + iprev*(iprev-1)))/2.0)
                ! Set the rank start and stop
                IF (rank == i - 1) THEN
                    istart = MAX(1,ncoords-CEILING(REAL(ifinishes(i)))+1)
                    istop = ncoords-istarts(i)+1
                ENDIF
            END DO    
        ELSE
            istart = 1 + rank*nper
            istop = MIN((rank+1)*nper, ncoords)
        END IF
        ! *************************************************************************

        IF (verbose) THEN
            write(*,*) "Rank ", rank, " made it to the double loop"
            write(*,*) "Rank ", rank, " is calculating distances ", istart, " to ", istop
        END IF
        
        ! Distance Calculation ****************************************************
        count = 0
        dr_count = 0
        Do i=istart, istop
            if (present(same_array)) then
                if (same_array) then
                    jstart = i+1
                else
                    jstart = 1
                end if
            else 
                jstart = 1
            end if

            Do j=jstart, size(r2,1)
                ! Periodic distance calculation
                rsq = periodic_distance2(r1(i,:), r2(j,:), box)
                ! Distance cutoff & store non-zero elements
                If (rsq < rc_sq) Then
                    dr_count = dr_count + 1
                    
                    ! This code does dynamics reallocation if needed to increase the size of the arrays
                    ! This could eventually be made a subroutine - but haven't done that yet.
                    IF (dr_count > size(dr)) THEN
                        ! Write out reallocations
                        IF (verbose) THEN
                            WRITE(*,*) "Rank ", rank, " is reallocating arrays"
                        END IF 

                        ! Note this block reallocates the arrays to bigger size if needed.
                        allocate(tmpdr(size(dr))); allocate(tmp1(size(id1))); allocate(tmp2(size(id2)))
                        tmp1 = id1; tmp2 = id2; tmpdr = dr
                        deallocate(id2); deallocate(id1); deallocate(dr)
                        allocate(dr(size(tmpdr)*2)); allocate(id1(size(tmp1)*2)); allocate(id2(size(tmp2)*2))
                        dr = 0; id1 = 0; id2 = 0
                        dr(1:size(tmpdr)) = tmpdr
                        id1(1:size(tmp1)) = tmp1
                        id2(1:size(tmp2)) = tmp2
                        deallocate(tmp2); deallocate(tmp1); deallocate(tmpdr)
                    END IF
                    ! End Dynamic Reallocation

                    dr(dr_count) = rsq
                    id1(dr_count) = i
                    id2(dr_count) = j
                EndIf
            EndDo
        EndDo

        ! Write to the screen
        IF (verbose) THEN
            write(*,*) "Rank ", rank, " has ", dr_count, " distances"
        END IF
        
        ! *** MPI *** *************************************************************
        CALL MPI_REDUCE(dr_count, count, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        allocate(dists(count)); allocate(atom1(count)); allocate(atom2(count))

        CALL MPI_GATHER(dr_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        IF (rank == 0) THEN
            displs(1) = 0
            DO i=2, nranks
                displs(i) = sum(counts(1:i-1))
            END DO
            !write(*,*) counts(:)
            !write(*,*) displs(:)
        END IF

        CALL MPI_GATHERV(dr, dr_count, MPI_REAL, dists, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_GATHERV(id1, dr_count, MPI_INTEGER, atom1, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_GATHERV(id2, dr_count, MPI_INTEGER, atom2, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        ! **************************************************************************
        IF (rank == 0 .AND. verbose) THEN 
            write(*,*) "double loop distances:", count
        ENDIF


        deallocate(dr)
        deallocate(id1)
        deallocate(id2)
        deallocate(ifinishes)
        deallocate(istarts)
        DEALLOCATE(counts, displs)
    end subroutine double_loop_distance