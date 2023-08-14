PROGRAM Spectroscopy

    !use MPI_F08
    use efield_module

    IMPLICIT NONE

    ! Input file Parameters
    CHARACTER(LEN=100), PARAMETER :: inputfile = "spectra_input.dat"
    INTEGER :: n_osc, ndigits, ntimes, ncorr, nksip
    REAL :: dt, w_resol
    REAL, DIMENSION(3) :: map_w01, map_w12, map_mu
    REAL, DIMENSION(2) :: map_x01, map_x12, map_alpha
    REAL :: T1_relax, T1_by_dt
    INTEGER :: num_Tw
    REAL, ALLOCATABLE :: Tw(:)
    INTEGER :: which_spectra
    REAL :: w1_min, w1_max, w3_min, w3_max

    ! Field Parameters
    INTEGER :: field_unit
    REAL, ALLOCATABLE :: field(:)
    REAL, ALLOCATABLE :: osc_vec(:,:)

    ! TCF parameters
    DOUBLE COMPLEX, ALLOCATABLE :: tcf(:)
    DOUBLE COMPLEX, ALLOCATABLE :: ir_total_tcf(:)


    ! 1) Read Calculation Input File
    ! 2) Read Field Files
    ! 3) Calculate Chosen Spectra
    ! 4) Write Results to file
    

    WRITE(*,*) "PROGRAM START: SPECTROSCOPY"

    ! 1) READ CALCULATION INPUT FILE
    CALL Read_Spectra_Input(inputfile &
                , n_osc, ndigits, ntimes, dt, ncorr, nksip, w_resol &
                , map_w01, map_w12, map_mu, map_x01, map_x12, map_alpha & 
                , T1_relax, T1_by_dt, num_Tw, Tw, which_spectra & 
                , w1_min, w1_max, w3_min, w3_max)

    ! Allocate memory for the field and TCFs
    ALLOCATE(field(ntimes))
    ALLOCATE(osc_vec(ntimes, 3))
    ALLOCATE(ir_total_tcf(ncorr), tcf(ncorr))
    ALLOCATE(w(ntimes), mu(ntimes,3))
    
    ! Calculate TCF based on options set in the input file
    ! Note - Currently there are five planned options, 
    ! but only 1 is implemented thus far:
    ! [X] IR
    ! [ ] Raman
    ! [ ] SFG
    ! [ ] 2D IR
    ! [ ] 2D SFG
    ! *****************************************************
    IF (which_spectra == 1) THEN 
        WRITE(*,*) "OPTION: IR Spectrum Calculation"

    ELSE IF (which_spectra == 2) THEN ! TODO
        WRITE(*,*) "OPTION: RAMAN TCF NOT YET IMPLEMENTED"
        STOP
    
    ELSE IF (which_spectra == 3) THEN ! TODO
        WRITE(*,*) "OPTION: SFG TCF NOT YET IMPLEMENTED"
        STOP

    ELSE IF (which_spectra == 4) THEN ! TODO
        WRITE(*,*) "OPTION: 2D IR TCF NOT YET IMPLEMENTED"
        STOP

    ELSE IF (which_spectra == 5) THEN ! TODO
        WRITE(*,*) "OPTION: 2D SFG TCF NOT YET IMPLEMENTED"
        STOP
    
    END IF

    ! 2) READ FIELD FILES AND CALCULATE OSC TCFS ********************
    DO i=1, n_osc
        ! Read FIELD File
        field_unit = i + 100
        CALL Read_Individual_Field(field_unit, field, osc_vec)

        ! Convert field to frequencies (gives w, mu as output)
        w=0.0; mu = 0.0 ! Zero before the subroutine call
        CALL Apply_Empirical_Map(field, osc_vec, map_w01, map_mu, map_x01, w, mu)

        IF (which_spectra == 1) THEN
            ! Calculate IR TCF
            CALL IR_TCF(ntimes, ncorr, nskip, w, mu, tcf)
            ir_total_tcf(:) = ir_total_tcf(:) + tcf(:)
        ! TODO ADD MORE OPTIONS
        END IF
    END DO
    IF (which_spectra == 1) THEN
        ir_total_tcf(:) = ir_total_tcf(:) / dcmplx(n_osc,0d0)
    END IF

    ! 3) Calculate the final TCFS and Spectra
    IF (which_spectra == 1) THEN
        CALL Calculate_IR_Spectra()
    ! TODO: Add more options
    END IF

    WRITE(*,*) "END PROGRAM: SPECTROSCOPY"

END PROGRAM Spectroscopy