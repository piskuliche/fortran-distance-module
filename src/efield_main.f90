PROGRAM efield_main

    USE MPI_F08
    USE efield_module
    use distance_module
    use linked_lists

    ! This program calculates the electric field experienced by a subset
    ! of selected atoms, exerted on it by all the charges in the system.
    ! This code is mpi enabled, and it uses the cell-linked-list approach for 
    ! distances within the cutoff. There is furthermore an option to avoid this cutoff
    ! by using a traditional double loop.

    ! General Procedure
    ! (1) Read in input file
    ! (2) Read in Coordinates
    ! (3) Read in Charges
    ! (4) Calculate Electric field
    ! (5) Write data

    IMPLICIT NONE

    

    


END PROGRAM efield_main