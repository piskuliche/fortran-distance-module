module linked_lists
    ! this module builds a linked list and implements its basic operations
    implicit none

    ! Parameters **************************************************************
    integer, parameter :: max_cells = 1000000
    integer, parameter :: max_bins = 500

    ! Variables ****************************************************************

contains

    include "assign_to_grid.f90"
    include "setup_cell_grid.f90"
    include "build_linked_list.f90"
    include "cell_internal_distance.f90"
    include "cell_list_distance.f90"

end module linked_lists
