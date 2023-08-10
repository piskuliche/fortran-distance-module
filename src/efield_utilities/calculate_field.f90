SUBROUTINE calculate_field(bonds, drx, dry, drz, dr, id1, id2, charges, field)
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: bonds
    REAL, DIMENSION(:), INTENT(IN) :: drx, dry, drz
    REAL, DIMENSION(:), INTENT(IN) :: dr, charges
    INTEGER, DIMENSION(:), INTENT(IN) :: id1, id2


    REAL, ALLOCATABLE, INTENT(OUT) :: field(:)

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
    END IF  

    CALL pick_subset(id1, id2, bonds, osc_bnd_indices)

    DO i=1, n_osc
        DO j=1, SIZE(id1)
            IF (id1(j) == i) THEN
                dr_field(i,:) = dr_field(i,:) + (-1) * charges(id2(j)) * dr_vec(j,:) / dr(j)**3
            ELSE IF (id2(j) == i) THEN
                dr_field(i,:) = dr_field(i,:) +        charges(id1(j)) * dr_vec(j,:) / dr(j)**3
            END IF
        END DO
    END DO

    DO i=1, n_osc
        field(i) = dot_product(dr_field(i,:), dr_vec(osc_bnd_indices(i),:))
        write(*,*) field(i)
    END DO
    
END SUBROUTINE calculate_field