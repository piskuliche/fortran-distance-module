SUBROUTINE Read_Spectra_Input(inputfile &
                , n_osc, ndigits, ntimes, dt, ncorr, nskip, w_resol &
                , map_w01, map_w12, map_mu, map_x01, map_x12, map_alpha & 
                , T1_relax, T1_by_dt, num_Tw, Tw, which_spectra & 
                , w1_min, w1_max, w3_min, w3_max)
    ! *********************************************************************
    ! 
    !
    ! Input:
    !   inputfile   -   Name of the input file (CHARACTER)
    !  
    ! Output:
    ! 
    ! *********************************************************************

    IMPLICIT NONE

    ! Input
    CHARACTER(len=*), INTENT(IN) :: inputfile

    ! Output
    INTEGER, INTENT(OUT) :: n_osc, ndigits
    INTEGER, INTENT(OUT) :: ntimes, ncorr, nskip
    REAL, INTENT(OUT) :: dt, w_resol
    REAL, DIMENSION(3), INTENT(OUT) :: map_w01, map_w12, map_mu
    REAL, DIMENSION(2), INTENT(OUT) :: map_x01, map_x12, map_alpha
    REAL, INTENT(OUT) :: T1_relax, T1_by_dt
    INTEGER, INTENT(OUT) :: num_Tw
    REAL, ALLOCATABLE, INTENT(OUT) :: Tw(:)
    INTEGER, INTENT(OUT) :: which_spectra
    REAL, INTENT(OUT) :: w1_min, w1_max, w3_min, w3_max

    ! Local
    INTEGER :: j

    ! MPI

    OPEN(11, FILE=inputfile, STATUS='OLD', ACTION='READ')

    READ(11,*)
    READ(11,*) n_osc, ndigits                           ! number of oscillators, number of digits in the field file
    READ(11,*)
    READ(11,*) ntimes, dt, ncorr, nskip                 ! number of origins, timestep, length of tcf, spacing between origins
    READ(11,*) 
    READ(11,*) w_resol                                  ! freq resolution (pad TCF w/ zeros)
    READ(11,*)
    READ(11,*) map_w01(1), map_w01(2), map_w01(3)       ! w01 freq map (units of cm^-1 for freq, au for field)
    READ(11,*) 
    READ(11,*) map_w12(1), map_w12(2), map_w12(3)       ! w12 freq map (units of cm^-1 for freq, au for field)
    read(11,*) 
    read(11,*) map_mu(1), map_mu(2), map_mu(3)          ! mu' map (units of au for field)
    read(11,*)
    read(11,*) map_x01(1), map_x01(2)                   ! x01 map (units of cm^-1 for freq, au for x)
    read(11,*)
    read(11,*) map_x12(1), map_x12(2)                   ! x12 map (units of cm^-1 for freq, au for x)
    read(11,*)  
    read(11,*) map_alpha(1), map_alpha(2)               ! Raman alpha map (units of au for field) 
    read(11,*) 
    read(11,*) T1_relax                                 ! T1 relaxation time (in fs)
    read(11,*)
    read(11,*) num_Tw                                   ! Number of waiting times
    read(11,*)
    read(11,*) which_spectra                            ! Which spectra to calculate
    read(11,*)
    ALLOCATE(Tw(num_Tw))
    read(11,*) (Tw(j), j=1, num_Tw)                     ! List of waiting times
    read(11,*)
    read(11,*) w1_min, w1_max, w3_min, w3_max           ! Ranges for the spectrum (cm^-1)

    CLOSE(11)

    ! Convert timestep to atomic units
    dt = dt/fs_per_au

    ! Convert freq. resol to atomic units
    w_resol = w_resol/cmi_per_au

    ! Convert map parameters to atomic units
    map_w01(:) = map_w01(:)/cmi_per_au
    map_w12(:) = map_w12(:)/cmi_per_au

    ! Convert vib. relax time to au
    T1_relax = T1_relax/fs_per_au
    T1_by_dt = T1_relax/dt

    ! Convert waiting times to atomic units
    Tw(:) = Tw(:)/fs_per_au

    ! Convert x01 and x12 parameters
    map_x01(1) = map_x01(1)/ang_per_au
    map_x01(2) = map_x01(2)*cmi_per_au/ang_per_au
    map_x12(1) = map_x12(1)/ang_per_au
    map_x12(2) = map_x12(2)*cmi_per_au/ang_per_au

    ! Convert frequency ranges
    w1_min = w1_min/cmi_per_au
    w1_max = w1_max/cmi_per_au
    w3_min = w3_min/cmi_per_au
    w3_max = w3_max/cmi_per_au
    
END SUBROUTINE Read_Spectra_Input