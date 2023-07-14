FUNCTION mesh_function(r, q, mesh_points, box) RESULT(rho)

    ! *********************************************************************
    ! Calculates the charge density on a mesh
    ! Adapted from Allen and Tildesley, Computer Simulation of Liquids
    ! *********************************************************************
    ! Inputs:
    !  r:  Positions of charges (n,3)
    !  q:  Charges (n)
    !  mesh_points: Number of mesh points in each direction (3)
    !  box: The box size (3)
    ! Outputs:
    !   rho: Charge density on mesh (sc,sc,sc)
    ! *********************************************************************

    IMPLICIT NONE
    ! Input Variables *****************************************************
    REAL,   DIMENSION(:,:), INTENT(IN) :: r     ! Positions (n,3)
    REAL,   DIMENSION(:),   INTENT(IN) :: q     ! Charges (n)
    INTEGER, DIMENSION(3)   INTENT(IN) :: mesh_points    ! Dimension of mesh
    ! Output Variables ****************************************************
    REAL, DIMENSION(0:mesh_points(1)-1,0:mesh_points(2)-1,0:mesh_points(3)-1) :: rho ! Charge density

    ! Local Variables *****************************************************
    INTEGER :: i, j, k
    INTEGER :: natoms
    INTEGER :: n1, n2, n3 ! Mesh point indices
    REAL :: q1, q2, q3 ! Charge weights

    INTEGER, DIMENSION(3) :: nindex ! Mesh point index
    REAL, DIMENSION(3) :: dr ! Charge positions relative to mesh
    REAL, DIMENSION(3) :: spacing ! Mesh spacing
    REAL, DIMENSION(3, -1:1) :: weights ! Weights in each coordinate direction

    ! *********************************************************************

    natoms = SIZE(q)

    ! Check input variables
    IF ( ANY( SHAPE(r) /= [natoms,3] ) ) THEN
        WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), natoms, 3
        STOP 'Error in mesh_function'
    END IF

    spacing(:) = box(:) / REAL( mesh_points(:) )

    rho = 0.0 ! Zero density

    DO i=1, natoms
        ! Assign points to grid, and use periodic boundary conditions
        nindex(:) = NINT( r(i,:) * mesh_points(:) )  ! Mesh point index
        nindex(:) = MOD( nindex(:), mesh_points(:) )   ! Periodic boundary conditions

        ! Get the vector from the charge to the mesh center
        dr(:) = r(i,:) - REAL( nindex(:) ) * spacing(:)
        dr(:) = dr(:) - box(:) * ANINT( dr(:)/box(:) ) 
        dr(:) = dr(:) / spacing(:)         ! Normalize by cell size   

        ! Weights for three-point assignment scheme
        v(:, -1) = 0.5 * ( 0.5 - dr(:) ) ** 2.0
        v(:,  0) = 0.75 - dr(:) ** 2.0
        v(:, 1) = 0.5 * ( 0.5 + dr(:) ) ** 2.0

        DO i = -1, 1
            q1 = q(i) * v(i, 1)
            n1 = MODULO( nindex(1) + i, mesh_points(:) )

            DO j = -1, 1
                q2 = q1 * v(j, 2)
                n2 = MODULO( nindex(2) + j, mesh_points(:) )

                DO k=-1, 1
                    q3 = q2 * v(k, 3)
                    n3 = MODULO( nindex(3) + k, mesh_points(:) )
                    rho(n1, n2, n3) = rho(n1, n2, n3) + q3
                END DO

            END DO

        END DO

        rho = rho / ( spacing(:) ** 3.0 )

    END DO

    
END FUNCTION mesh_function