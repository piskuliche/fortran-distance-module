MODULE spectra_module

    IMPLICIT NONE

    DOUBLE PRECISION :: fs_per_au=2.4188843265857d-2
    DOUBLE PRECISION :: cmi_per_au=2.1947463d5
    DOUBLE PRECISION :: ang_per_au=0.5291772083d0
    DOUBLE PRECISION :: pi=4.0d0*datan(1.0d0)

CONTAINS
    include "ir_tcf.f90"
    include "read_spectra_input.f90"
    include "ir_spectra.f90"
    include "frequency_map.f90"

END MODULE spectra_module