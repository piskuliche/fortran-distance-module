SUBROUTINE calculate_field(bonds, drx, dry, drz, dr, id1, id2, charges, osc_grps, grp_count, field, dipole_vec)
    USE MPI_F08
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
    INTEGER :: n_osc_per_rank, rank_count
    INTEGER :: iter, istart, istop ! Iteration variables for each rank
    REAL, ALLOCATABLE :: field_proc(:) ! MPI field array before gatherv
    REAL, DIMENSION(3) :: dipole_proc
    REAL, ALLOCATABLE :: dpx_proc(:), dpy_proc(:), dpz_proc(:)
    REAL, ALLOCATABLE :: dpx(:), dpy(:), dpz(:)
    INTEGER, ALLOCATABLE :: counts(:), displs(:)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierror)
    
    dr_vec(:,1) = drx
    dr_vec(:,2) = dry
    dr_vec(:,3) = drz

    ! Number of oscillators
    n_osc = SIZE(bonds, 1)
    n_osc_per_rank = CEILING(REAL(n_osc)/REAL(nranks))
    istart = max(1, n_osc_per_rank*rank + 1)
    istop = min(n_osc, n_osc_per_rank*(rank+1))
    rank_count = istop-istart+1
    WRITE(*,*) "test-mpi-field", istart, istop, rank_count
    WRITE(*,*) size(id1)


    ! Maximum oscillator group number - sets indexing for
    max_osc = MAXVAL(osc_grps)
    IF ( rank == 0 ) THEN
        IF (.NOT. ALLOCATED(field)) THEN

            ALLOCATE(field(n_osc))
            ALLOCATE(dipole_vec(n_osc,3))

        END IF  

        field = 0.0
        dipole_vec = 0.0
        
        ALLOCATE(counts(nranks), displs(nranks))
        ALLOCATE(dpx(n_osc), dpy(n_osc), dpz(n_osc))
    END IF

    ALLOCATE(jgroup1(size(id1)), jgroup2(size(id2)))
    ALLOCATE(osc_sum(max_osc))
    ALLOCATE(dr_field(rank_count,3))
    ALLOCATE(field_contribution(max_osc,3))
    ALLOCATE(field_proc(rank_count))
    ALLOCATE(dpx_proc(rank_count), dpy_proc(rank_count), dpz_proc(rank_count))


    ! First key thing to do is to pick out 
    ! the distances that are relevant to the oscillators
    ! pick_subset(id1, id2, bonds, osc_bnd_indices)

    !CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    dr_field = 0.0
    field_proc = 0.0
    dpx_proc = 0.0; dpy_proc = 0.0; dpz_proc = 0.0
    dipole_proc = 0.0

    DO i=istart, istop
        iter = i - istart + 1
        ! Calculate the field at the ith oscillator
        ! 1) Loop over all distances calculated
        osc_sum = 0
        field_contribution = 0.0
        ! Get the field contributions
        CALL Calculate_Field_Contribution(id1, id2, osc_grps, dr_vec, dr & 
                , bonds(i,:), charges, osc_sum, field_contribution, dipole_proc(:))
        ! Set these variables based on dipole_proc
        dpx_proc(iter) = dipole_proc(1)
        dpy_proc(iter) = dipole_proc(2)
        dpz_proc(iter) = dipole_proc(3)
        DO j=1, max_osc
            !write(*,*) osc_sum(j), grp_count(j)
            ! only add the contribution from j if the total number within the osc group
            ! matches what is expected - this is to maintain molecules whole.
            IF (osc_sum(j) == grp_count(j)) THEN
                dr_field(iter,:) = dr_field(iter,:) + field_contribution(j,:)
            END IF
        END DO 

        ! Convert units to Atomic Units
        dr_field(iter,:) = dr_field(iter,:) * angperau**2.0

        ! Calculate the final field value, for the processor
        field_proc(iter) = dot_product(dr_field(iter,:), dipole_proc(:))
    END DO

    ! TODO: Generate counts and displs arrays for gatherv
    IF (rank == 0) THEN
        counts(1) = rank_count
        DO i = 2, nranks
            counts(i) = min(n_osc_per_rank*i, n_osc) - counts(i-1)
        END DO

        displs(1) = 0
        DO i=2, nranks
            displs(i) = sum(counts(1:i-1))
        END DO
        WRITE(*,*) counts
        WRITE(*,*) displs
    END IF

    ! TODO: Gather the field values
    CALL MPI_GATHERV(field_proc, rank_count, MPI_REAL &
        , field, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD)
    ! TODO: Gatehr the dipole values for x, y and z components
    CALL MPI_GATHERV(dpx_proc, rank_count, MPI_REAL &
        , dpx, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD)
    CALL MPI_GATHERV(dpy_proc, rank_count, MPI_REAL &
        , dpy, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD)
    CALL MPI_GATHERV(dpz_proc, rank_count, MPI_REAL &
        , dpz, counts, displs, MPI_REAL, 0, MPI_COMM_WORLD)
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    ! TODO: Set rank 0 values for the dipoole
    IF (rank == 0) THEN
        dipole_vec(:,1) = dpx
        dipole_vec(:,2) = dpy
        dipole_vec(:,3) = dpz
    END IF 

    DEALLOCATE(jgroup1)
    DEALLOCATE(jgroup2)
    DEALLOCATE(osc_sum)
    DEALLOCATE(dr_field)
    DEALLOCATE(field_contribution)
    DEALLOCATE(field_proc)
    DEALLOCATE(dpx_proc, dpy_proc, dpz_proc)
    IF (rank == 0) THEN
        DEALLOCATE(counts)
        DEALLOCATE(displs)
    ENDIF


    
END SUBROUTINE calculate_field


SUBROUTINE Calculate_Field_Contribution(id1, id2, osc_grps, dr_vec, dr, bond_atoms, charges, osc_sum, field_contribution, dipole)
! 2) For every distance calculated, check if ith oscillator hatom is involved
            !    Also check that the two atoms are not in the same oscillator group  
            ! Then calculate the field contributions

    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) :: id1 ! Atom ID of first atom in distance pair
    INTEGER, DIMENSION(:), INTENT(IN) :: id2 ! Atom ID of second atom in distance pair
    INTEGER, DIMENSION(:), INTENT(IN) :: osc_grps ! Oscillator group of each atom
    REAL, DIMENSION(:,:), INTENT(IN) :: dr_vec ! Distance vector between atoms id1 and id2
    REAL, DIMENSION(:), INTENT(IN) :: dr ! Distance between atoms id1 and id2
    INTEGER, DIMENSION(2), INTENT(IN) :: bond_atoms ! Atom IDs of the two atoms in the ith bond
    REAL, DIMENSION(:), INTENT(IN) :: charges ! Charges of each atom

    INTEGER, DIMENSION(:), INTENT(INOUT) :: osc_sum ! Sum variable for field contributions based on osc_grps
    REAL, DIMENSION(:,:), INTENT(INOUT) :: field_contribution ! Field contribution from each oscillator group
    REAL, DIMENSION(3), INTENT(OUT) :: dipole ! Dipole vector for the ith oscillator


    INTEGER :: jgroup1, jgroup2
    INTEGER :: j

    DO j = 1, size(id1)
        dipole = 0.0

        jgroup1 = osc_grps(id1(j))
        jgroup2 = osc_grps(id2(j))

        ! Check that the values come from different oscillator groups
        IF (jgroup1 /= jgroup2 ) THEN
            IF (id1(j) == bond_atoms(2)) THEN
                field_contribution(jgroup2,:) = field_contribution(jgroup2,:)  & 
                    + Add_Field(charges(id2(j)), dr_vec(j,:), dr(j))
                osc_sum(jgroup2) = osc_sum(jgroup2) + 1
            ELSE IF (id2(j) == bond_atoms(2)) THEN
                field_contribution(jgroup1,:) = field_contribution(jgroup1,:)  &
                    - Add_Field(charges(id1(j)), dr_vec(j,:), dr(j))
                osc_sum(jgroup1) = osc_sum(jgroup1)+ 1
            END IF
        ELSE ! Same oscillator group
            ! This section grabs the dipole vector for the oscillator
            ! It always points towards the hydrogen atom
            IF (id1(j) == bond_atoms(1) .and. id2(j) == bond_atoms(2)) THEN
                dipole(:) = -dr_vec(j,:)/sqrt(dr(j))
            ELSE IF (id1(j) == bond_atoms(2) .and. id2(j) == bond_atoms(1)) THEN
                dipole(:) = dr_vec(j,:)/sqrt(dr(j))
            END IF
        END IF
    END DO

END SUBROUTINE Calculate_Field_Contribution

FUNCTION Add_Field(q, dr_vec, dr) RESULT(comp)
    IMPLICIT NONE

    REAL, INTENT(IN) :: q ! Charge of atom
    REAL, DIMENSION(:), INTENT(IN) :: dr_vec ! Vector between atoms
    REAL, INTENT(IN) :: dr ! Distance between atoms

    REAL :: comp(3)

    comp = 0.0

    comp = q * dr_vec / (sqrt(dr)**3)

END FUNCTION Add_Field