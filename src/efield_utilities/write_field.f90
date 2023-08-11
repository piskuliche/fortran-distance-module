SUBROUTINE Write_Field_Files(field, dipole_vec, initialize)
    ! *********************************************************************
    ! This subroutine takes all of the field and dipole vectors and writes
    ! them to field files using the Write_Individiual_Field subroutine.
    !
    ! Input:
    !   field - the field values (REAL, DIMENSION(:))
    !   dipole_vec - the dipole vector (REAL, DIMENSION(:,:))
    !   initialize - if true, initialize the files (LOGICAL, OPTIONAL)
    !
    ! Output:
    !   none
    ! 
    ! TODO:
    ! 1) Add MPI SUPPORT
    ! 2) Add subdirectory support
    ! 
    ! *********************************************************************

    IMPLICIT NONE

    ! Input
    REAL, DIMENSION(:), INTENT(IN) :: field
    REAL, DIMENSION(:,:), INTENT(IN) :: dipole_vec
    LOGICAL, OPTIONAL, INTENT(IN) :: initialize

    ! Output

    ! Local
    INTEGER :: i
    LOGICAL :: init

    ! MPI
    !INTEGER :: rank, nranks, ierror

    !CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    !CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)

    IF (PRESENT(initialize) .and. initialize) THEN
        init = .TRUE.
    ELSE
        init = .FALSE.
    END IF
    
    DO i=1, size(field)
        CALL Write_Individual_Field(i+100, field(i), dipole_vec(i,:), initialize=init)
    END DO

END SUBROUTINE Write_Field_Files

SUBROUTINE Write_Individual_Field(unit, field, dipole_vec, initialize)

    ! *********************************************************************
    ! This subroutine writes an individual field file (field.unit)
    ! Now - note that if you have lots of molecules, this will be a lot of files
    ! Thus - it is best to use subdirectories to store the files
    !
    ! Input:
    !   unit - the file unit to write to (INTEGER)
    !   field - the field value (REAL)
    !   dipole - the dipole moment (REAL, DIMENSION(3))
    !   initialize - if true, initialize the file (LOGICAL, OPTIONAL)
    !  
    ! Output:
    !   none
    ! 
    ! TODO:
    !  1) Add subdirectories
    ! 
    ! *********************************************************************

    IMPLICIT NONE

    ! Input
    INTEGER, INTENT(IN) :: unit
    REAL, INTENT(IN) :: field
    REAL, DIMENSION(3), INTENT(IN) :: dipole_vec
    LOGICAL, OPTIONAL, INTENT(IN) :: initialize

    ! Local
    CHARACTER(len=40) ext


    IF (unit < 100) THEN
        WRITE(*,*) "File Unit", unit, "is too small. Must be greater than 100."
        STOP 1
    END IF

    write(ext,'(I0.5)') unit

    IF (PRESENT(initialize) .and. initialize) THEN
        OPEN(unit, file="field."//ext, action="write")
    ELSE
        OPEN(unit, file="field."//ext, status="old", position="append", action="write")   
    END IF

    write(unit, '(4F10.5)') field, dipole_vec

    close(unit)

END SUBROUTINE Write_Individual_Field