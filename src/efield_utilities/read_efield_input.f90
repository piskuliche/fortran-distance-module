SUBROUTINE read_efield_input(inputfile, natoms, nframes, n_osc, rc &
                        , charges, oscs, bonds, osc_grps, grp_count &
                        , traj_file, traj_format)
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
    !   bonds: Mapping of bonds in the system (2D Array of size n_osc, 2)
    !   traj_file : Trajectory filename (Character)
    !   traj_format : Integer (1 = xtc, 2 = xyz)
    ! 
    ! Allocates charges, oscs

    IMPLICIT NONE

    ! Input Variables
    CHARACTER(LEN=*), intent(in) :: inputfile

    ! Output Variables

    INTEGER, INTENT(OUT) :: natoms, nframes, n_osc
    REAL, INTENT(OUT) :: rc
    REAL, ALLOCATABLE, INTENT(OUT) :: charges(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: oscs(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: bonds(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: osc_grps(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: grp_count(:)
    ! Local Variables
    INTEGER :: traj_format
    CHARACTER(LEN=40) :: charge_file, osc_file, osc_grp_file
    CHARACTER(LEN=40) :: traj_file, bond_file

    write(*,*) TRIM(inputfile)
    ! Open the input file
    OPEN(21, FILE=TRIM(inputfile), STATUS="OLD")
    
    
    READ(21,*)  ! Trajectory file name, format of trajectory (1:xtc, 2:xyz)
    READ(21,*) traj_file, traj_format
    READ(21,*)  ! charge_file name, oscillator_name
    READ(21,*) charge_file, bond_file
    READ(21,*) 
    READ(21,*) osc_file, osc_grp_file
    READ(21,*) 
    READ(21,*) natoms, nframes
    READ(21,*) 
    READ(21,*) rc

    ALLOCATE(charges(natoms))
    ALLOCATE(oscs(natoms))
    ALLOCATE(osc_grps(natoms))

    charges = 0.0
    oscs = 0
    osc_grps = 0.0


    CALL read_charges(charge_file, natoms, charges)
    CALL read_osc(osc_file, osc_grp_file, natoms, oscs, n_osc, osc_grps, grp_count) ! Allocates grp_count
    CALL read_osc_bonds(bond_file, bonds) ! Allocates bonds


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

SUBROUTINE read_osc(osc_file, osc_grp_file, natoms, contributes, count, osc_grp, grp_count)
    ! This subroutine reads the oscillator file and stores information
    ! about which atoms to calculate the location of the electric field
    ! and the strength and direction on them. 

    IMPLICIT NONE

    ! Input Variables
    CHARACTER(LEN=*), INTENT(IN) :: osc_file, osc_grp_file
    INTEGER, INTENT(IN) :: natoms
    ! Output Variables
    INTEGER, DIMENSION(natoms), INTENT(OUT) :: contributes(:)
    INTEGER, INTENT(OUT) :: count
    INTEGER, DIMENSION(:), INTENT(OUT) :: osc_grp
    INTEGER, ALLOCATABLE, INTENT(OUT) :: grp_count(:)
    ! Local Variables
    INTEGER :: i
    INTEGER :: n_osc_grp
    INTEGER :: max_osc


    ! Initialize to zero charge
    contributes = 0.0

    count = 0

    OPEN(22, FILE=TRIM(osc_file), STATUS='old')
    DO i=1, natoms
        read(22,*) contributes(i)
        IF (contributes(i) == 1) THEN
            count = count + 1
        END IF
    END DO
    CLOSE(22)

    OPEN(22, FILE=trim(osc_grp_file), STATUS='old')

    DO i=1, natoms
        READ(22,*) osc_grp(i)
    END DO

    max_osc = MAXVAL(osc_grp)

    ALLOCATE(grp_count(max_osc))

    grp_count = 0

    DO i=1, natoms
        grp_count(osc_grp(i)) = grp_count(osc_grp(i)) + 1
    END DO

END SUBROUTINE read_osc

SUBROUTINE Molecular_Linked_List(osc_grps, group_list, group_head)
    ! This code takes a list of oscillators and creates a linked list
    ! This is so that each of these oscillators can then be used to
    ! calculate the electric field at each of the atoms in the system.
    ! While maintaining oscillator groups as whole.

    ! This is currently unused.

    IMPlICIT NONE
    INTEGER, PARAMETER :: max_atoms=100
    INTEGER, DIMENSION(:), INTENT(IN) :: osc_grps
    INTEGER, ALLOCATABLE, INTENT(OUT) :: group_list(:), group_head(:)
    
    INTEGER :: n_groups, n_values

    INTEGER :: i


    n_groups = MAXVAL(osc_grps)
    n_values = size(osc_grps)

    IF (.NOT. ALLOCATED(group_list)) THEN
        ALLOCATE(group_list(n_values))
        ALLOCATE(group_head(n_groups))
    END IF

    group_head = 0
    group_list = 0

    DO i=1, n_values
        group_list(i) = group_head(osc_grps(i))
        group_head(osc_grps(i)) = i
    END DO


END SUBROUTINE