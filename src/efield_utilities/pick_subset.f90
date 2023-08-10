SUBROUTINE pick_subset(r, osc, r_osc)

    REAL, DIMENSION(:,:), INTENT(IN) :: r
    INTEGER, DIMENSION(:), INTENT(IN) :: osc

    REAL, ALLOCATABLE, INTENT(OUT) :: r_osc(:,:)

    INTEGER :: i, count, n_osc

    IF (size(osc) /= size(r,1)) THEN
        PRINT *, "ERROR: size(osc) /= size(r,1)"
        STOP
    END IF

    n_osc = sum(osc)
    write(*,*) "n_osc = ", n_osc
    IF (.NOT. ALLOCATED(r_osc)) THEN
        ALLOCATE(r_osc(n_osc, 3))
    END IF

    count = 1
    DO i=1, n_osc
        IF (osc(i) == 1) THEN
            r_osc(count,:) = r(i,:)
            count = count + 1
        END IF
    END DO



END SUBROUTINE pick_subset