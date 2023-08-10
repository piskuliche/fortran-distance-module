SUBROUTINE read_osc_bonds(bndname, bonds)
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
    CHARACTER(len=*), INTENT(IN) :: bndname

    ! Output
    INTEGER, ALLOCATABLE, INTENT(OUT) :: bonds(:,:)

    ! Local
    INTEGER :: i
    INTEGER :: nbonds

    ! Read the bond file
    OPEN(12, FILE=trim(bndname), STATUS='old')

        ! Read the number of bonds
        READ(12,*) nbonds
        WRITE(*,*) 'Number of bonds: ', nbonds

        ! Allocate the array if needed
        IF (.NOT. ALLOCATED(bonds)) THEN
            ALLOCATE(bonds(nbonds,2))
        END IF

        ! Initialize
        bonds = 0

        ! Read the bonds into bonds
        DO i=1,nbonds
            READ(12,*) bonds(i,1), bonds(i,2)
        END DO

    CLOSE(12)

END SUBROUTINE read_osc_bonds