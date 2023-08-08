SUBROUTINE read_efield_input(inputfile, natoms, nframes, n_osc, charges, oscs)
    ! This subroutine reads an electric field calculation input file for
    ! use to generate field files for a spectroscopy calculation.

    ! Input:
    !   inputfile: input filename for the calculation (Character)
    ! 
    ! Output:
    !   natoms : number of atoms in the simulation (Integer)
    !   nframes : number of frames in the simulation (Integer)
    !   n_osc : Number of oscillators in the simulation (Integer)
    !   charges : Simulation charges (1D Array of size natoms)
    !   oscs : Mapping of atoms to oscillators (1D Array of size natoms)
    ! 
    ! Allocates charges, oscs

    IMPLICIT NONE

    ! Input Variables
    CHARACTER(LEN=*), intent(in) :: inputfile

    ! Output Variables

    INTEGER, INTENT(OUT) :: natoms, nframes, n_osc
    REAL, ALLOCATABLE, INTENT(OUT) :: charges(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: oscs(:)
    ! Local Variables
    INTEGER :: traj_format
    INTEGER :: nframes
    CHARACTER(LEN=40) :: charge_file, osc_file, traj_file

    ! Open the input file
    OPEN(21, FILE=TRIM(inputfile), STATUS="OLD")
    
    
    READ(21,*)  ! Trajectory file name, format of trajectory (1:xtc, 2:xyz)
    READ(21,*) traj_file, traj_format
    READ(21,*)  ! charge_file name, oscillator_name
    READ(21,*) charge_file, osc_file
    READ(21,*) 
    READ(22,*) natoms, nframes

    ALLOCATE(charges(natoms))
    ALLOCATE(oscs(natoms))

    read_charges(chargE_file, natoms, charges)
    read_osc(osc_file, natoms, oscs, n_osc)


    CLOSE(21)
 


END SUBROUTINE read_efield_input


SUBROUTINE read_charges(charge_file, natoms, charges)
    ! This subroutine reads the charge file and stores information
    ! about the charge on each atom within the system.

    IMPLICIT NONE

    ! Input Variables
    CHARACTER(LEN=*), INTENT(IN) :: charge_file
    INTEGER, INTENT(IN) :: natoms
    ! Output Variables
    REAL, DIMENSION(natoms), INTENT(OUT) :: charges(:)
    ! Local Variables
    INTEGER :: i

    ! Initialize to zero charge
    charges = 0.0

    OPEN(22, FILE=TRIM(charge_file), STATUS='old')
    DO i=1, natoms
        read(22,*) charges(i)
    END DO
    CLOSE(22)

END SUBROUTINE read_charges

SUBROUTINE read_osc(osc_file, natoms, contributes, count)
    ! This subroutine reads the oscillator file and stores information
    ! about which atoms to calculate the location of the electric field
    ! and the strength and direction on them. 

    IMPLICIT NONE

    ! Input Variables
    CHARACTER(LEN=*), INTENT(IN) :: osc_file
    INTEGER, INTENT(IN) :: natoms
    ! Output Variables
    INTEGER, DIMENSION(natoms), INTENT(OUT) :: contributes(:)
    INTEGER, INTENT(OUT) :: count
    ! Local Variables
    INTEGER :: i


    ! Initialize to zero charge
    charges = 0.0

    count = 0

    OPEN(22, FILE=TRIM(osc_file), STATUS='old')
    DO i=1, natoms
        read(22,*) contributes(i)
        IF (contributes(i) == 1) THEN
            count = count + 1
    END DO
    CLOSE(22)

END SUBROUTINE read_charges