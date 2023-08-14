SUBROUTINE IR_TCF(ntimes, ncorr, nskip, w, mu, tcf)
    ! *********************************************************************
    ! 
    !
    ! Input:
    !  
    ! Output:
    ! 
    ! *********************************************************************

    IMPLICIT NONE

    ! Input
    INTEGER, INTENT(IN) :: ntimes, ncorr, nskip
    DOUBLE PRECISION, DIMENSION(ntimes) :: w
    DOUBLE PRECISION, DIMENSION(ntimes,3) :: mu

    ! Output
    DOUBLE COMPLEX, DIMENSION(ncorr), INTENT(OUT) :: tcf

    ! Local
    INTEGER :: i, j, ilag
    DOUBLE PRECISION :: origin_count, phase, mu0_dot_mu
    DOUBLE PRECISION, DIMENSION(3) :: mu0
    DOUBLE PRECISION :: phase

    ! Initialize the complex TCF
    tcf = dcmplx(0.0d0,0.0d0)

    origin_count = 0.0d0
    ! Loop over time origins
    DO i = 1, ntimes - ncorr, nskip
        ! Set the original dipole vector
        mu0(:) = mu(i,:)
        origin_count = origin_count + 1d0

        ! Calculate the phase of the TCF
        phase = -w(i)

        ! Loop over the lag times
        DO j = i, i + ncorr
            ! Set lag time index
            ilag = j - i + 1 
            phase = phase + w(j)
            ! Calculate the dot_product of the t0 dipole with the time
            ! t dipole moment.
            mu0_dot_mu = dot_product(mu0, mu(j,:))
            ! Calculate the time correlation function
            tcf(ilag) = tcf(ilag) + mu0_dot_mu * dcmplx(dcos(phase*dt), dsin(phase*dt))
        END DO
    END DO

    ! Normalize the TCF by the number of time origins
    tcf = tcf/dcmplx(origin_count,0d0)

END SUBROUTINE