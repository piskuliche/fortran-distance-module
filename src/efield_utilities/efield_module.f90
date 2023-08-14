MODULE efield_module

    IMPLICIT NONE

    REAL, PARAMETER :: angperau = 0.52917721092d0

CONTAINS

    INCLUDE "read_efield_input.f90"
    INCLUDE "pick_subset.f90"
    INCLUDE "read_osc_bonds.f90"
    INCLUDE "calculate_field.f90"
    INCLUDE "write_field.f90"
    INCLUDE "read_field.f90"

END MODULE efield_module