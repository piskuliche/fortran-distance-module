PROGRAM efield_main

    USE MPI_F08
    use gmxfort_trajectory
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
    
    type (Trajectory) :: trj

    ! Command Line would eventually be nice for this one
    CHARACTER(LEN=40) :: inputfile

    ! Loop indices
    INTEGER :: i, j, k
    INTEGER :: frame

    ! Efield Input Subroutine Arguments
    INTEGER :: natoms, nframes, n_osc, traj_format
    REAL, ALLOCATABLE :: charges
    INTEGER, ALLOCATABLE :: oscs
    CHARACTER(len=40) :: traj_fname

    ! Trajectory variables
    REAL, ALLOCATABLE :: r_osc(:,:), r(:,:)
    REAL, DIMENSION(3) :: coord, box
    REAL, DIMENSION(3,3) :: trj_cell

    ! MPI Variables
    INTEGER :: rank, nranks, ierror

    ! Initialize MPI
    CALL MPI_INIT(ierror)

    ! Get the number of processors & the rank of this one
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    ! (1 & 2) Reads in the input file, the charges, and the oscilators
    CALL read_efield_input(inputfile & ! Input
       , natoms, nframes, n_osc, charges, oscs, traj_fname, traj_format) ! Output

    IF (traj_format == 1) THEN
        WRITE(*,*) "Trajectory filetype is XTC"
        call trj%open(trim(traj_fname),"index.ndx")
    ELSE IF (traj_format == 2) THEN
        WRITE(*,*) "Support for XYZ filetypes has not been implemented yet"
        WRITE(*,*) "Stopping program."
        STOP
    ELSE
        WRITE(*,*) "Error: Only supported options are [1] or [2] for traj_format"
        WRITE(*,*) "Stopping program."
        STOP
    END IF

    DO frame=1, nframes
        ! (3) Read Coordinates
        IF (rank == 0) THEN
            trj%read_next(1)
            ! Read the 
            DO i=1, natoms 
                coord(:) = trj%x(1,i)
                r(i,:)= coord(:)
            END DO 
        END IF

        ! Broadcast the coordinates to all ranks
        CALL MPI_BCAST(r, size(r), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

        ! (4a) Calculate the distances
        CALL cell_list_distance(r, r, box, cell_length, rc_sq, ll_dr, ll_id1, ll_id2, ll_count, same_array=.true.)

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        ! (4b) Turn distances into the electric field
        
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
        ! (5) Write out to a file
        ! Lets rank 0 move onto reading in the next coords
        ! while the other ranks write out.
        if (rank /= 0) THEN

        END IF

    END DO 



    
    


    

    


END PROGRAM efield_main