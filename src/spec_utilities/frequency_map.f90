SUBROUTINE Apply_Empirical_Map(efield, osc_vec, map_w01, map_mu, map_x01, w, mu)

    IMPLICIT NONE

    ! Input
    REAL, DIMENSION(:), INTENT(IN) :: efield
    REAL, DIMENSION(:,:), INTENT(IN) :: osc_vec
    REAL, DIMENSION(3), INTENT(IN) :: map_w01, map_mu
    REAL, DIMENSION(2), INTENT(IN) :: map_x01

    ! Output
    DOUBLE PRECISION, DIMENSION(size(efield)), INTENT(OUT) :: w
    DOUBLE PRECISION, DIMENSION(size(efield),3),INTENT(OUT) :: mu

    ! Local
    INTEGER :: i

    w = 0.0; mu = 0.0
    DO i=1, size(efield)
        w(i) = Get_Frequency(map_w01, efield(i))
        mu(i,:) = Get_Transition_Dipole(w(i), map_x01, map_mu, efield(i), osc_vec(i,:))
    END DO

END SUBROUTINE Apply_Empirical_Map

FUNCTION Get_Frequency(map_w01, efield) RESULT(frequency)
    ! This function applies the traditional quadratic frequency
    ! map used for the empiricall mapping approach to get the 
    ! frequency.
    ! **********
    IMPLICIT NONE

    ! Input
    REAL, DIMENSION(3), INTENT(IN) :: map_w01
    REAL, INTENT(IN) :: efield

    ! Output
    DOUBLE PRECISION :: frequency

    frequency = map_w01(1) + map_w01(2)*efield + map_w01(3)*efield**2
    
END FUNCTION Get_Frequency

FUNCTION Get_Transition_Dipole(freq, map_x01, map_mu, efield, osc_vec) RESULT(tr_mu)
    IMPLICIT NONE

    ! Input
    DOUBLE PRECISION, INTENT(IN) :: freq
    REAL, DIMENSION(3), INTENT(IN) :: map_mu
    REAL, DIMENSION(2), INTENT(IN) :: map_x01
    REAL, INTENT(IN) :: efield
    REAL, DIMENSION(3) :: osc_vec

    ! Output
    REAL, DIMENSION(3) :: tr_mu

    ! Local
    REAL :: xtmp, mutmp

    xtmp = map_x01(1) + map_x01(2)*freq

    mutmp = map_mu(1) + map_mu(2)*efield + map_mu(3)*efield**2
    
    tr_mu(:) = osc_vec(:)*mutmp*xtmp
    
END FUNCTION Get_Transition_Dipole

FUNCTION Get_Transition_Polarizability(freq, map_x01, map_alpha, efield) RESULT(tr_pol)
    IMPLICIT NONE

    ! Input
    DOUBLE PRECISION, INTENT(IN) :: freq
    REAL, DIMENSION(2), INTENT(IN) :: map_x01
    REAL, DIMENSION(2), INTENT(IN) :: map_alpha
    REAL, INTENT(IN) :: efield

    ! Output
    REAL :: tr_pol

    ! Local
    REAL :: xtmp

    xtmp = map_x01(1) + map_x01(2)*freq

    tr_pol = xtmp * (map_alpha(1) + map_alpha(2)*efield)
    
END FUNCTION Get_Transition_Polarizability