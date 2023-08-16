PROGRAM efield_main

    USE mpi_f08
    USE coordinates_module
    USE efield_module
    use distance_module
    use linked_lists

    ! This program calculates the electric field experienced by a subset
    ! of selected atoms, exerted on it by all the charges in the system.
    ! This code is mpi enabled, and it uses the cell-linked-list approach for 
    ! distances within the cutoff. There is furthermore an option to avoid this cutoff
    ! by using a traditional double loop.

    ! General Procedure
    ! (1) Read in input file (x)
    ! (2) Read in the charges (x) 
    ! -> start loop
    ! (3) Read in Coordinates ()
    ! (4) Calculate Electric field ()
    ! (5) Write data to files ()
    ! -> end loop

    IMPLICIT NONE

    ! Command Line would eventually be nice for this one
    CHARACTER(LEN=40) :: inputfile

    ! Loop indices
    INTEGER :: i, j, k
    INTEGER :: frame

    INTEGER :: boxunit, trjunit

    ! Efield Input Subroutine Arguments
    INTEGER :: natoms, nbonds
    INTEGER :: nframes, n_osc, traj_format
    REAL, ALLOCATABLE :: charges(:) 
    INTEGER, ALLOCATABLE :: oscs(:), osc_bnd_indices(:)
    INTEGER, ALLOCATABLE :: bonds(:,:)
    INTEGER, ALLOCATABLE :: osc_grps(:)
    INTEGER, ALLOCATABLE :: grp_count(:)

    REAL :: rc, rc_sq

    CHARACTER(len=40) :: traj_fname

    ! Trajectory variables
    REAL, ALLOCATABLE :: r(:,:), r_osc(:,:)
    REAL, DIMENSION(3) :: coord, L
    REAL, DIMENSION(3,3) :: box
    REAL, DIMENSION(3,3) :: cell

    ! Distance variables
    REAL, ALLOCATABLE :: dr(:)
    INTEGER, ALLOCATABLE :: id1(:), id2(:)
    INTEGER :: drcount, ntmp, ocount
    LOGICAL :: load_balance
    LOGICAL :: include_vector
    REAL, ALLOCATABLE :: drx(:), dry(:), drz(:)

    ! Electric field variables
    REAL, ALLOCATABLE :: efield(:)
    REAL, ALLOCATABLE :: dipole_vec(:,:)
    
    ! Writing Variables
    LOGICAL :: init_field_files


    ! MPI Variables
    INTEGER :: rank, nranks, ierror

    ! Initialize MPI
    CALL MPI_INIT(ierror)

    ! Get the number of processors & the rank of this one
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    ! Set unit numbers
    boxunit = 10
    trjunit = 11
    include_vector=.true.

    inputfile = "test_input.inp"
    ! (1 & 2) Reads in the input file, the charges, and the oscilators
    IF (rank == 0) THEN
        ! Read the input file for the run
        ! This also reads the charges, oscillators, and bonds
        CALL read_efield_input(inputfile & ! Input
       , natoms,  nframes, n_osc, rc, charges, oscs, bonds, osc_grps, grp_count &
       , traj_fname, traj_format) ! Output

        nbonds = size(bonds,1)
        IF (traj_format == 1) THEN
            WRITE(*,*) "Trajectory filetype is XTC"
            WRITE(*,*) "This hasn't been implemented yet - stopping!"
            STOP 1
            !WRITE(*,*) trim(traj_fname)
            !call trj%open(traj_fname)
            !WRITE(*,*) "Opened!"
        ELSE IF (traj_format == 2) THEN
            WRITE(*,*) "Trajectory filetype is XYZ"
            open(trjunit, file=traj_fname, status='old', action='read')
            open(boxunit, file="L.dat", status='old', action='read')
        ELSE
            WRITE(*,*) "Error: Only supported options are [1] or [2] for traj_format"
            WRITE(*,*) "Stopping program."
            STOP
        END IF
        rc_sq = rc**2.
    END IF
    WRITE(*,*) "rank", rank, "has read the input file"

    ! Broadcast input information to all ranks
    CALL MPI_BCAST(natoms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BCAST(nbonds, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

    IF (rank > 0) THEN
        ALLOCATE(r(natoms,3))
        ALLOCATE(bonds(nbonds,2))
    END IF

    ! These are allocated already on rank 0 in the read_efield_input subroutine
    IF (rank > 0) THEN
        ALLOCATE(charges(natoms))
        ALLOCATE(oscs(natoms))
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
    CALL MPI_BCAST(rc_sq, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BCAST(charges, size(charges), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BCAST(oscs, size(oscs), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BCAST(nframes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BCAST(bonds, size(bonds), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

    if (rank == 0) THEN
        WRITE(*,*) "rc_sq: ", rc_sq
        WRITE(*,*) "nframes: ", nframes
        WRITE(*,*) "size(charges): ", size(charges)
        WRITE(*,*) "size(oscs): ", size(oscs)
    END IF

    ! End broadcast
    load_balance  = .true.

    ntmp = 0
    WRITE(*,*) "rank", rank, "has allocated the charges and oscs arrays"
    if (rank == 0) THEN
        WRITE(*,*) "Starting frame loop"
    END IF

    ! Loop over frames
    DO frame=1, nframes
        ! (3) Read Coordinates
        IF (rank == 0) THEN
            write(*,*) "Reading frame", frame
            CALL read_next_xyz(trjunit, r)
            CALL read_next_L(boxunit, L)
            write(*,*) "Read", ntmp, "atoms"
            write(*,*) "L:", L
        END IF

        ! Broadcast the coordinates to all ranks
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        CALL MPI_BCAST(r, size(r), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_BCAST(L, size(L), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

        ! (4a) Calculate the distances
        !CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)
        CALL double_loop_distance(r, r, L, rc_sq, dr, id1, id2 &
                    , drcount, same_array=.true., cell_length=0.0 &
                    , load_balance=load_balance, include_vector=include_vector &
                    , drx=drx, dry=dry, drz=drz, verbosity=0)



        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        ! (4b) Turn distances into the electric field
        IF (rank == 0) THEN
            write(*,*) "Found ", drcount, "Distances Total"

            CALL calculate_field(bonds, drx, dry, drz, dr, id1, id2, charges, osc_grps, grp_count, efield, dipole_vec)
            
        END IF
        
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        ! (5) Write out to a file
        ! Lets rank 0 move onto reading in the next coords
        ! while the other ranks write out.
        if (rank == 0) THEN
            write(*,*) "Finished frame", frame
            
            ! Don't overwrite files on higher frames
            IF (frame == 1) THEN
                init_field_files = .true.
            ELSE
                init_field_files = .false.
            END IF

            CALL Write_Field_Files(efield, dipole_vec, initialize=init_field_files)
        END IF

    END DO 

    deallocate(r, charges, oscs)
    IF (ALLOCATED(osc_grps)) DEALLOCATE(osc_grps)

    IF (ALLOCATED(efield)) THEN
        deallocate(efield)
    END IF

    CALL MPI_FINALIZE(ierror)


END PROGRAM efield_main