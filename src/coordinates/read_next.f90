SUBROUTINE read_next_xyz(unitnumber, r)

    IMPLICIT NONE

    ! Input Arguments
    INTEGER, INTENT(IN) :: unitnumber

    ! Output Arguments
    REAL, ALLOCATABLE, INTENT(OUT) :: r(:,:)


    ! Local Variables
    INTEGER :: i, j, natoms
    CHARACTER(len=4) :: ctmp

    read(unitnumber,*) natoms
    read(unitnumber,*)

    IF (.NOT. ALLOCATED(r)) THEN
        ALLOCATE(r(natoms,3))
    END IF

    DO i=1,natoms
        read(unitnumber,*) ctmp, r(i,1), r(i,2), r(i,3)
    END DO

END SUBROUTINE read_next_xyz

SUBROUTINE read_next_L(unitnumber, L)

    IMPLICIT NONE

    ! Input Arguments
    INTEGER, INTENT(IN) :: unitnumber

    ! Output Arguments
    REAL, DIMENSION(3), INTENT(OUT) :: L

    read(unitnumber,*) L(1), L(2), L(3)


END SUBROUTINE read_next_L