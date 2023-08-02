SUBROUTINE CSR_INTEGER_WRITER(io_unit, csr, csr_filename)
    implicit none
    INTEGER, INTENT(IN) :: io_unit
    INTEGER, DIMENSION(:), INTENT(IN) :: csr
    CHARACTER(len=*), INTENT(IN) :: csr_filename 

    INTEGER :: i, status

    ! Open file for CSR Writer
    OPEN(io_unit, file=trim(csr_filename))
    
    ! Write CSR to file
    DO i = 1, size(csr)
        WRITE(io_unit,*) csr(i)
    END DO

    ! Close file for CSR Writer
    CLOSE(io_unit)


END SUBROUTINE CSR_INTEGER_WRITER


SUBROUTINE CSR_REAL_WRITER(io_unit, csr, csr_filename)
    implicit none
    INTEGER, INTENT(IN) :: io_unit
    REAL, DIMENSION(:), INTENT(IN) :: csr
    CHARACTER(len=*), INTENT(IN) :: csr_filename 


    INTEGER :: i, status

    ! Open file for CSR Writer
    OPEN(io_unit, file=trim(csr_filename))
    
    ! Write CSR to file
    DO i = 1, size(csr)
        WRITE(io_unit,*) csr(i)
    END DO

    ! Close file for CSR Writer
    CLOSE(io_unit)

END SUBROUTINE CSR_REAL_WRITER