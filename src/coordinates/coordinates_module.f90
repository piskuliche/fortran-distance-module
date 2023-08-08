MODULE coordinates_module
    ! This module contains subroutines for reading coordinates from molecular simulations
    ! For use within these codes.

    IMPLICIT NONE

CONTAINS 

    INCLUDE "read_xyz.f90"
    INCLUDE "coordinate_generator.f90"

END MODULE coordinates_module