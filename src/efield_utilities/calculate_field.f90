SUBROUTINE calculate_field(bonds, drx, dry, drz, dr, id1, id2, charges, osc_grps, field, dipole_vec)
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: bonds
    REAL, DIMENSION(:), INTENT(IN) :: drx, dry, drz
    REAL, DIMENSION(:), INTENT(IN) :: dr, charges
    INTEGER, DIMENSION(:), INTENT(IN) :: id1, id2
    INTEGER, DIMENSION(:), INTENT(IN) :: osc_grps


    REAL, ALLOCATABLE, INTENT(OUT) :: field(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: dipole_vec(:,:)

    INTEGER :: i, j, k
    INTEGER :: n_osc

    REAL, ALLOCATABLE :: dr_field(:,:)
    INTEGER, ALLOCATABLE :: osc_bnd_indices(:)
    REAL, DIMENSION(SIZE(drx),3) :: dr_vec
    
    dr_vec(:,1) = drx
    dr_vec(:,2) = dry
    dr_vec(:,3) = drz

    ! Number of oscillators
    n_osc = SIZE(bonds, 1)
    IF (.NOT. ALLOCATED(field)) THEN
        ALLOCATE(dr_field(n_osc,3))
        ALLOCATE(field(n_osc))
        ALLOCATE(dipole_vec(n_osc,3))
    END IF  

    CALL pick_subset(id1, id2, bonds, osc_bnd_indices)

    dr_field = 0.0
    dipole_vec = 0.0
    field = 0.0

    DO i=1, n_osc
        ! Calculate the field at the ith oscillator
        DO j=1, SIZE(id1)
            IF (id1(j) == bonds(i,2) .and. (osc_grps(id1(j)) /= osc_grps(id2(j)))) THEN
                WRITE(*,*) "id1(j) = ", id1(j), "id2(j) = ", id2(j), bonds(i,2)
                dr_field(i,:) = dr_field(i,:) +             charges(id2(j)) * dr_vec(j,:) / sqrt(dr(j))**3
            ELSE IF (id2(j) == bonds(i,2) .and. (osc_grps(id1(j)) /= osc_grps(id2(j)))) THEN
                WRITE(*,*) "id1(j) = ", id1(j), "id2(j) = ", id2(j), bonds(i,2)
                dr_field(i,:) = dr_field(i,:) +  (-1) *     charges(id1(j)) * dr_vec(j,:) / sqrt(dr(j))**3
            END IF
        END DO

        ! Convert units
        dr_field(i,:) = dr_field(i,:) * angperau**2.0
    END DO
    

    DO i=1, n_osc
        dipole_vec(i,:) = -dr_vec(osc_bnd_indices(i),:)
        field(i) = dot_product(dr_field(i,:), dipole_vec(i,:))
    END DO
    
END SUBROUTINE calculate_field