module distance_module
    ! This module contains the functions and subroutines for calculating the distance between
    ! all pairs of atoms in a system (within a cutoff) using a cell linked list.
    !


    implicit none

contains

    ! **************************************************************************
    ! Functions ****************************************************************
    ! **************************************************************************



    function periodic_distance2(r1, r2, box) result(dr_sq)
    ! This function calculates the square of the distance between two atoms
    ! in a periodic system.
    !
    ! Inputs:
    !   r1: 1D array of atom positions
    !   r2: 1D array of atom positions
    !   box: 1D array of box lengths
    !
    ! Returns:
    !   dr_sq: Square of the distance between the two atoms
    !
        implicit none

        ! Inputs **************************************************************
        REAL, DIMENSION(3), INTENT(in) :: r1, r2, box
        ! Outputs *************************************************************
        REAL :: dr_sq
        ! Local Variables *****************************************************
        INTEGER :: i
        REAL, DIMENSION(3) :: dr_tmp
        ! *********************************************************************

        dr_sq = 0.0; dr_tmp = 0.0
        Do i=1,3
            ! Calculate the distance between the two atoms
            dr_tmp(i) = r1(i) - r2(i)
            ! Apply periodic boundary conditions
            dr_tmp(i) = dr_tmp(i) - box(i)*nint(dr_tmp(i)/box(i))
        EndDo
        dr_sq = sum(dr_tmp**2)

    end function periodic_distance2

    function periodic_distance_and_vector(r1, r2, box) result(dr_arr)
    ! This function calculates the square of the distance between two atoms
    ! in a periodic system.
    !
    ! Inputs:
    !   r1: 1D array of atom positions
    !   r2: 1D array of atom positions
    !   box: 1D array of box lengths
    !
    ! Returns:
    !   dr_arr: 1D array that has first three elements as the distance vector
    !           and the last element is the distance.
    !
        implicit none

        ! Inputs **************************************************************
        REAL, DIMENSION(3), INTENT(in) :: r1, r2, box
        ! Outputs *************************************************************
        REAL, DIMENSION(4) :: dr_arr
        ! Local Variables *****************************************************
        INTEGER :: i
        REAL :: dr
        
        REAL, DIMENSION(3) :: dr_tmp
        ! *********************************************************************

        dr = 0.0; dr_arr = 0.0; dr_tmp = 0.0
        Do i=1,3
            ! Calculate the distance between the two atoms
            dr_tmp(i) = r1(i) - r2(i)
            ! Apply periodic boundary conditions
            dr_arr(i) = dr_tmp(i) - box(i)*nint(dr_tmp(i)/box(i))
        EndDo
        ! Calculate distance
        dr = sqrt(sum(dr_tmp**2))
        ! Assign values to dr_arr and get unit vector
        dr_arr(1:3) = dr_tmp(1:3)/dr
        ! Assign distance to dr_arr
        dr_arr(4) = dr

    end function periodic_distance_and_vector

    function angle_between_points(r1, r2, r3, box) result(theta)
    ! This function calculates the angle between three points in a periodic system.
    !
    implicit none

    ! Inputs **************************************************************
    REAL, DIMENSION(3), INTENT(in) :: r1, r2, r3, box
    ! Outputs *************************************************************
    REAL :: theta
    ! Local Variables *****************************************************
    REAL, DIMENSION(4) :: dr12, dr23
    ! *********************************************************************

    dr12 = 0.0; dr23 = 0.0
    dr12 = periodic_distance_and_vector(r2, r1, box)
    dr23 = periodic_distance_and_vector(r2, r3, box)
    theta = acos(dot_product(dr12(1:3), dr23(1:3)))

    end function angle_between_points
    ! **************************************************************************


    ! **************************************************************************
    ! Subroutines **************************************************************
    ! **************************************************************************

    subroutine double_loop_distance(r1, r2, box, rc_sq &
        , dists, atom1, atom2, cell_assign_1, cell_assign_2, count, same_array, cell_length)
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
        ! Outputs *************************************************************
        REAL, ALLOCATABLE, INTENT(out) :: dists(:)
        INTEGER, ALLOCATABLE, INTENT(out) :: atom1(:), atom2(:)
        INTEGER, DIMENSION(:,:), INTENT(out) :: cell_assign_1, cell_assign_2
        INTEGER, INTENT(out) :: count    ! Number of pairs found
        ! Local Variables *****************************************************
        INTEGER :: i, j, di ! Loop Indices
        REAL :: rsq         ! Temporary distance squared
        INTEGER :: jstart
        INTEGER :: ierror, nranks, rank
        INTEGER :: istart, istop, ncoords, nper
        INTEGER :: dr_count


        REAL, DIMENSION(3) :: dr_tmp        ! Temporary distance vector
        INTEGER, DIMENSION(3,500) :: map    ! Map of bin indices

        REAL, ALLOCATABLE :: dr(:), tmpdr(:) ! TEMPORARY Distance
        INTEGER, ALLOCATABLE :: id1(:), id2(:), tmp1(:), tmp2(:)
        INTEGER, ALLOCATABLE :: counts(:), displs(:)
        ! ************************************************************************
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ALLOCATE(counts(nranks))
        ALLOCATE(displs(nranks))

        ! Initialize the cell_assign arrays
        cell_assign_1 = 0; cell_assign_2 = 0

        write(*,*) "Rank ", rank, " made it to the double loop"

        ! Allocate initial block for arrays
        IF ( .NOT. ALLOCATED(dr) ) THEN
            allocate(dr(1000)); allocate(id1(1000)); allocate(id2(1000))
        ENDIF

        ! Get the number of coordinates of r1
        ncoords = size(r1, 1)
        ! Calculate the number of atoms per proc
        nper = CEILING(REAL(ncoords/nranks))
        istart = 1 + rank*nper
        istop = MIN((rank+1)*nper, ncoords)
        write(*,*) "Rank ", rank, " is calculating distances ", istart, " to ", istop

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
            end if
            Do j=jstart, size(r2,1)
                ! Periodic distance calculation
                rsq = periodic_distance2(r1(i,:), r2(j,:), box)
                ! Distance cutoff & store non-zero elements
                If (rsq < rc_sq) Then
                    dr_count = dr_count + 1

                    IF (dr_count > size(dr)) THEN
                        WRITE(*,*) "Rank ", rank, " is reallocating arrays"
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

                    dr(dr_count) = rsq
                    id1(dr_count) = i
                    id2(dr_count) = j
                EndIf
            EndDo
        EndDo
        write(*,*) "Rank ", rank, " has ", dr_count, " distances"

        CALL MPI_REDUCE(dr_count, count, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        allocate(dists(count)); allocate(atom1(count)); allocate(atom2(count))

        CALL MPI_GATHER(dr_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        IF (rank == 0) THEN
            displs(1) = 0
            DO i=2, nranks
                displs(i) = sum(counts(1:i-1))
            END DO
            write(*,*) counts(:)
            write(*,*) displs(:)
        END IF

        CALL MPI_GATHERV(dr, dr_count, MPI_REAL, dists, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierror)

        ! **************************************************************************
        IF (rank == 0) THEN 
            write(*,*) "dcount", count
        ENDIF
        deallocate(dr)
        deallocate(id1)
        deallocate(id2)
        DEALLOCATE(counts, displs)
    end subroutine double_loop_distance

end module distance_module