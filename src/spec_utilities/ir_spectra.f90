SUBROUTINE Calculate_IR_Spectra(total_tcf, T1, dt, ncorr, w_resol,  wavg)
    ! *********************************************************************
    ! This subroutine applies a Fourier Transform to the total TCF
    ! to obtain the IR spectrum.
    !
    ! Input:
    !  
    ! Output:
    ! 
    ! *********************************************************************
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
    include "fftw3.f03"

    ! Input
    DOUBLE COMPLEX, DIMENSION(:), INTENT(INOUT) :: total_tcf
    REAL, INTENT(IN) :: T1, dt, wavg, w_resol
    INTEGER, INTENT(IN) :: ncorr

    ! Output
    !       None!

    ! Local
    INTEGER :: i, p, nt
    TYPE(C_PTR) :: plan
    DOUBLE PRECISION :: w_spec, ti, pi
    complex(C_DOUBLE_COMPLEX), ALLOCATABLE :: input(:), ir_spec(:)

    ! Open the Output File
    OPEN(22, file='IR_Spectrum.dat', status='unknown')

    ! (1) Apply vibrational relaxation factor
    CALL Apply_T1_and_Exp(T1, dt, wavg, ncorr, total_tcf)

    ! (2) Setup the Fourier Transform
    ! Calculate the time domain grid needed to get the frequency resolution
    p = int(dlog(2d0*pi/(dt*w_resol))/dlog(2d0))
    nt = 2**p

    ALLOCATE(input(nt), ir_spec(nt))
    
    ! Plan the FFT
    plan = fftw_plan_dft_1d(nt, input, ir_spec, FFTW_FORWARD, FFTW_ESTIMATE)

    ! (3) Pad the calculated TCF with zeros to get freq resolution before FT
    input = dcmplx(0.0d0, 0.0d0)
    DO i=1, min(nt, ncorr)
        input(i) = total_tcf(i)
    END DO
    
    CALL fftw_execute_dft(plan, input, ir_spec)

    ! (4) Write the IR spectrum to file
    ! When the average frequency is subtracted from the phase, the 
    ! spectrum must be split in them middle.
    ! (4a) Write the frequencies below the average
    DO i = nt/2 + 1, nt
        w_spec = (2d0*pi/dt)*dfloat(i-1-nt)/dfloat(nt) + wavg
        WRITE(22,*) w_spec*cmi_per_au, dreal(ir_spec(i))
    END DO

    ! (4b) Write the frequencies above the average
    DO i = 1, nt/2
        w_spec = (2d0*pi/dt)*dfloat(i-1)/dfloat(nt) + wavg
        WRITE(22,*) w_spec*cmi_per_au, dreal(ir_spec(i))
    END DO

    ! Destroy the FFTW plan
    CALL fftw_destroy_plan(plan)

    ! Deallocate the arrays
    DEALLOCATE(input, ir_spec)
    CLOSE(22)

END SUBROUTINE Calculate_IR_Spectra

SUBROUTINE Apply_T1_and_Exp(T1, dt, wavg, ncorr, tcf)
    IMPLICIT NONE

    ! Input
    REAL, INTENT(IN) :: T1, dt, wavg
    INTEGER, INTENT(IN) :: ncorr
    DOUBLE COMPLEX, INTENT(INOUT) :: tcf(:)


    ! Local
    INTEGER :: i
    DOUBLE PRECISION :: ti

    DO i=1, ncorr
        ti = float(i-1)*dt
        tcf(i) = tcf(i)*dcmplx(dexp(-0.5d0)*ti/T1, 0d0)
    END DO

    DO i=1, size(tcf)
        ti = float(i-1)*dt
        tcf(i) = tcf(i)*dcmplx(dcos(ti*wavg),-dsin(ti*wavg))
    END DO

END SUBROUTINE Apply_T1_and_Exp