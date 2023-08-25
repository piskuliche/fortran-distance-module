SUBROUTINE calculate_field(bonds, drx, dry, drz, dr, id1, id2, charges, osc_grps, grp_count, field, dipole_vec)
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: bonds
    REAL, DIMENSION(:), INTENT(IN) :: drx, dry, drz
    REAL, DIMENSION(:), INTENT(IN) :: dr, charges
    INTEGER, DIMENSION(:), INTENT(IN) :: id1, id2
    INTEGER, DIMENSION(:), INTENT(IN) :: osc_grps
    INTEGER, DIMENSION(:), INTENT(IN) :: grp_count ! Array that stores counts of atoms in each group


    REAL, ALLOCATABLE, INTENT(OUT) :: field(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: dipole_vec(:,:)

    INTEGER :: i, j, k
    INTEGER :: n_osc, newcount
    INTEGER :: max_osc

    REAL, ALLOCATABLE :: dr_field(:,:)
    REAL, ALLOCATABLE :: field_contribution(:,:)
    INTEGER, ALLOCATABLE :: osc_bnd_indices(:)
    INTEGER, ALLOCATABLE :: jgroup1(:), jgroup2(:), osc_sum(:)
    REAL, DIMENSION(SIZE(drx),3) :: dr_vec

    ! MPI Variables
    INTEGER :: ierror, nranks, rank
    
    dr_vec(:,1) = drx
    dr_vec(:,2) = dry
    dr_vec(:,3) = drz

    ! Number of oscillators
    n_osc = SIZE(bonds, 1)

    ! Maximum oscillator group number - sets indexing for
    max_osc = MAXVAL(osc_grps)

    IF (.NOT. ALLOCATED(field)) THEN

        ALLOCATE(field(n_osc))
        ALLOCATE(dipole_vec(n_osc,3))

    END IF  

    ALLOCATE(jgroup1(size(id1)), jgroup2(size(id2)))
    ALLOCATE(osc_sum(max_osc))
    ALLOCATE(dr_field(n_osc,3))
    ALLOCATE(field_contribution(max_osc,3))

    ! First key thing to do is to pick out 
    ! the distances that are relevant to the oscillators
    ! pick_subset(id1, id2, bonds, osc_bnd_indices)

    !CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    dr_field = 0.0
    dipole_vec = 0.0
    field = 0.0

    DO i=1, n_osc
        ! Calculate the field at the ith oscillator
        ! 1) Loop over all distances calculated
        newcount = 0
        osc_sum = 0
        field_contribution = 0.0
        DO j=1, SIZE(id1)
            ! 2) For every distance calculated, check if ith oscillator hatom is involved
            !    Also check that the two atoms are not in the same oscillator group  
            ! Then calculate the field contributions
            jgroup1(j) = osc_grps(id1(j))
            jgroup2(j) = osc_grps(id2(j))
            ! Check that the values come from different oscillator groups
            IF (jgroup1(j) /= jgroup2(j)) THEN
                IF (id1(j) == bonds(i,2)) THEN
                    newcount = newcount + 1
                    field_contribution(jgroup2(j),:) = field_contribution(jgroup2(j),:)  & 
                        +    charges(id2(j)) * dr_vec(j,:) / sqrt(dr(j))**3
                    osc_sum(jgroup2(j)) = osc_sum(jgroup2(j)) + 1
                ELSE IF (id2(j) == bonds(i,2)) THEN
                    newcount = newcount + 1
                    field_contribution(jgroup1(j),:) = field_contribution(jgroup1(j),:)  &
                        +    (-1) *     charges(id1(j)) * dr_vec(j,:) / sqrt(dr(j))**3
                    osc_sum(jgroup1(j)) = osc_sum(jgroup1(j))+ 1
                END IF
                ! Add something to check whether all the atoms in a molecule are present
                ! in the list of atoms.
            ELSE
                ! This section grabs the dipole vector for the oscillator
                ! It always points towards the hydrogen atom
                IF (id1(j) == bonds(i,1) .and. id2(j) == bonds(i,2)) THEN
                    dipole_vec(i,:) = -dr_vec(j,:)/sqrt(dr(j))
                ELSE IF (id1(j) == bonds(i,2) .and. id2(j) == bonds(i,1)) THEN
                    dipole_vec(i,:) = dr_vec(j,:)/sqrt(dr(j))
                END IF
            END IF
        END DO

        DO j=1, max_osc
            !write(*,*) osc_sum(j), grp_count(j)
            ! only add the contribution from j if the total number within the osc group
            ! matches what is expected - this is to maintain molecules whole.

            IF (osc_sum(j) == grp_count(j)) THEN
                dr_field(i,:) = dr_field(i,:) + field_contribution(j,:)
            END IF
        END DO 

        
        IF (i == 1 .or. i == 2) THEN
            WRITE(*,*) newcount
            WRITE(*,*) dr_field(i,:)
        ENDIF

        ! Convert units
        dr_field(i,:) = dr_field(i,:) * angperau**2.0
    END DO
    

    DO i=1, n_osc
        !dipole_vec(i,:) = -dr_vec(osc_bnd_indices(i),:)/sqrt(dr(osc_bnd_indices(i)))
        field(i) = dot_product(dr_field(i,:), dipole_vec(i,:))
    END DO

    DEALLOCATE(jgroup1, jgroup2)
    DeALLOCATE(osc_sum)
    DEALLOCATE(dr_field)
    DEALLOCATE(field_contribution)


    
END SUBROUTINE calculate_field