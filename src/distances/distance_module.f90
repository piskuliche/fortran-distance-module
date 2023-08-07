module distance_module
    ! This module contains the functions and subroutines for calculating the distance between
    ! all pairs of atoms in a system (within a cutoff) using a cell linked list.
    !


    implicit none

contains

    include "functions.f90"
    include "double_loop_distance.f90"

end module distance_module