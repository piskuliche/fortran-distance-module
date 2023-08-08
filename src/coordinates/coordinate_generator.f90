    subroutine coordinate_generator(natoms, box, r)
        implicit none
        ! Input Arguments
        integer, intent(in) :: natoms
        real, dimension(3), intent(in) :: box
        real, dimension(natoms,3), intent(out) :: r
        ! Local Variables
        integer :: i, j
        real, dimension(3) :: rtmp

    
        do i=1, natoms
            rtmp = 0.0
            do j=1, 3
                call random_number(rtmp(j))
            enddo
            r(i,:) = rtmp(:)*box(:)
        end do

    end subroutine coordinate_generator