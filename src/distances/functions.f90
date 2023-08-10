
    function periodic_distance2(r1, r2, box) result(dr_sq)
    ! This function calculates the square of the distance between two atoms
    ! in a periodic system.
    !
    ! Inputs:
    !   r1: 1D array of atom positions
    !   r2: 1D array of atom positions
    !   box: 1D array of box lengths
    !
    ! Returns:
    !   dr_sq: Square of the distance between the two atoms
    !
        implicit none

        ! Inputs **************************************************************
        REAL, DIMENSION(3), INTENT(in) :: r1, r2, box
        ! Outputs *************************************************************
        REAL :: dr_sq
        ! Local Variables *****************************************************
        INTEGER :: i
        REAL, DIMENSION(3) :: dr_tmp
        ! *********************************************************************

        dr_sq = 0.0; dr_tmp = 0.0
        Do i=1,3
            ! Calculate the distance between the two atoms
            dr_tmp(i) = r1(i) - r2(i)
            ! Apply periodic boundary conditions
            dr_tmp(i) = dr_tmp(i) - box(i)*nint(dr_tmp(i)/box(i))
        EndDo
        dr_sq = sum(dr_tmp**2)

    end function periodic_distance2

    

    function periodic_distance_and_vector(r1, r2, box) result(dr_arr)
    ! This function calculates the square of the distance between two atoms
    ! in a periodic system.
    !
    ! Inputs:
    !   r1: 1D array of atom positions
    !   r2: 1D array of atom positions
    !   box: 1D array of box lengths
    !
    ! Returns:
    !   dr_arr: 1D array that has first three elements as the distance vector
    !           and the last element is the distance.
    !
        implicit none

        ! Inputs **************************************************************
        REAL, DIMENSION(3), INTENT(in) :: r1, r2, box
        ! Outputs *************************************************************
        REAL, DIMENSION(4) :: dr_arr
        ! Local Variables *****************************************************
        INTEGER :: i
        REAL :: dr
        
        REAL, DIMENSION(3) :: dr_tmp
        ! *********************************************************************

        dr = 0.0; dr_arr = 0.0; dr_tmp = 0.0
        Do i=1,3
            ! Calculate the distance between the two atoms
            dr_tmp(i) = r1(i) - r2(i)
            ! Apply periodic boundary conditions
            dr_tmp(i) = dr_tmp(i) - box(i)*nint(dr_tmp(i)/box(i))
        EndDo
        ! Calculate distance
        dr = sqrt(sum(dr_tmp**2))
        ! Assign values to dr_arr and get unit vector
        dr_arr(1:3) = dr_tmp(1:3)/dr
        ! Assign distance to dr_arr
        dr_arr(4) = dr

    end function periodic_distance_and_vector

    function angle_between_points(r1, r2, r3, box) result(theta)
    ! This function calculates the angle between three points in a periodic system.
    !
    implicit none

    ! Inputs **************************************************************
    REAL, DIMENSION(3), INTENT(in) :: r1, r2, r3, box
    ! Outputs *************************************************************
    REAL :: theta
    ! Local Variables *****************************************************
    REAL, DIMENSION(4) :: dr12, dr23
    ! *********************************************************************

    dr12 = 0.0; dr23 = 0.0
    dr12 = periodic_distance_and_vector(r2, r1, box)
    dr23 = periodic_distance_and_vector(r2, r3, box)
    theta = acos(dot_product(dr12(1:3), dr23(1:3)))

    end function angle_between_points