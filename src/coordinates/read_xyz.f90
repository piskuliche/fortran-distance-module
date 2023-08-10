SUBROUTINE read_xyz(filename, r)

    IMPLICIT NONE

    ! Input Arguments
    CHARACTER(len=*), INTENT(IN) :: filename

    ! Output Arguments
    REAL, ALLOCATABLE, INTENT(OUT) :: r(:,:)


    ! Local Variables
    INTEGER :: i, j, natoms




    OPEN(11, FILE=TRIM(filename), STATUS='old')

    read(11,*) natoms
    read(11,*)

    ALLOCATE(r(natoms,3))

    DO i=1,natoms
        read(11,*) j, r(i,1), r(i,2), r(i,3)
    END DO

END SUBROUTINE read_xyz