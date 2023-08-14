SUBROUTINE Read_Field_Files(nframes, nfiles, field, dipole_vec)

    IMPLICIT NONE

    ! Input Variables
    INTEGER, INTENT(IN) :: nfiles
    INTEGER, INTENT(IN) :: nframes

    ! Output variables
    REAL, DIMENSION(nfiles, nframes) :: field
    REAL, DIMENSION(nfiles, nframes, 3) :: dipole_vec

    INTEGER :: i
    
    DO i=1, nfiles
        CALL Read_Individual_Field(i+100, nframes, field(i,:), dipole_vec(i,:,:))
    END DO

END SUBROUTINE Read_Field_Files

SUBROUTINE Read_Individual_Field(unit, nframes, field, dipole_vec)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit     ! File unit number
    INTEGER, INTENT(IN) :: nframes  ! Number of frames
    REAL, DIMENSION(nframes), INTENT(OUT) :: field ! Field values at each frame
    REAL, DIMENSION(nframes,3), INTENT(OUT) :: dipole_vec ! Dipole vector at each frame

    ! Local
    CHARACTER(len=40) ext
    INTEGER :: i

    IF (unit < 100) THEN
        WRITE(*,*) "File Unit", unit, "is too small. Must be greater than 100."
        STOP 1
    END IF

    OPEN(unit, FILE='field_files/field."'//ext, action="read")

    DO i=1, nframes
        READ(unit, *) field(i), dipole_vec(i,1), dipole_vec(i,2), dipole_vec(i,3)
    END DO

END SUBROUTINE Read_Individual_Field
