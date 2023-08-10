SUBROUTINE pick_subset(dr, id1, id2, bonds, osc_indices)
    ! *********************************************************************
    ! 
    !
    ! Input:
    !  
    ! Output:
    ! 
    ! *********************************************************************

    IMPLICIT NONE

    ! Input
    REAL, dimension(:), INTENT(IN) :: dr
    INTEGER, dimension(:), INTENT(IN) :: id1, id2
    INTEGER, dimension(:,:), INTENT(IN) :: bonds

    ! Output
    INTEGER, ALLOCATABLE, INTENT(out) :: osc_indices(:)

    ! Local
    INTEGER :: i, j
    INTEGER :: count

    ! Could do an MPI SCATTER here, but not necessary

    IF (.NOT. ALLOCATED(osc_indices)) THEN
        ALLOCATE(osc_indices(size(bonds,1)))
    END IF

    count = 0
    DO i=1, size(bonds,1)
        DO j=1, size(dr)
            IF (bonds(i,1) == id1(j) .AND. bonds(i,2) == id2(j)) THEN
                count = count + 1
                osc_indices(count) = j
            ELSE IF (bonds(i,1) == id2(j) .AND. bonds(i,2) == id1(j)) THEN
                count = count + 1
                osc_indices(count) = j
            END IF
        END DO
    END DO



END SUBROUTINE