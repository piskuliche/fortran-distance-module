PROGRAM Spectroscopy

    !use MPI_F08

    IMPLICIT NONE

    ! Input file Parameters
    CHARACTER(LEN=100), PARAMETER :: inputfile = "spectra_input.dat"
    INTEGER :: n_osc, ndigits, ntimes, ncorr, nksip
    REAL :: dt, w_resol
    REAL, DIMENSION(3) :: map_w01, map_w12, map_mu
    REAL, DIMENSION(2) :: map_x01, map_x12
    REAL :: T1_relax, T1_by_dt
    INTEGER :: num_Tw
    REAL, ALLOCATABLE :: Tw(:)
    INTEGER :: which_spectra
    REAL :: w1_min, w1_max, w3_min, w3_max

    ! 1) Read Calculation Input File
    ! 2) Read Field Files
    ! 3) Calculate Chosen Spectra
    ! 4) Write Results to file
    

    WRITE(*,*) "PROGRAM START: SPECTROSCOPY"

    ! 1) READ CALCULATION INPUT FILE
    Read_Spectra_Input(inputfile &
                , n_osc, ndigits, ntimes, dt, ncorr, nksip, w_resol &
                , map_w01, map_w12, map_mu, map_x01, map_x12, & 
                , T1_relax, T1_by_dt, num_Tw, Tw, which_spectra & 
                , w1_min, w1_max, w3_min, w3_max)

    ! 2) READ FIELD FILES

    ! 3) CALCULATE CHOSEN SPECTRA

    ! 4) WRITE RESULTS TO FILE

    WRITE(*,*) "END PROGRAM: SPECTROSCOPY"

END PROGRAM Spectroscopy