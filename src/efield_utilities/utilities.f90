FUNCTION calculate_field_at_point(dr_vec, dr, q) RESULT(field)
    IMPLICIT NONE
    REAL, DIMENSION(3) :: dr_vec
    REAL :: dr, q, field
    
    field = q * dr_vec * dr / (dr**3)

END FUNCTION