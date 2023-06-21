module shared

    implicit none

    integer, parameter :: max_cells = 1000000

    integer, dimension(3) :: nbins
    integer, dimension(max_cells) :: list_r1, list_r2
    integer, dimension(:,:,:), allocatable :: head_r1, head_r2

contains
    
    function assign_bins(r, nbins, box) result(arr)
        implicit none
        integer, dimension(3), intent(in) :: nbins
        real, dimension(3), intent(in) :: r, box

        integer :: i
        real, dimension(3) :: arr
        arr = 0
        Do i=1,3
            arr(i) = min(nbins(i), max(1, int(r(i)/box(i)*nbins(i))))
        EndDo

    end function assign_bins

    subroutine build_linked_list(r, nbins, box, head, list)
        implicit none

        integer, dimension(3), intent(in) :: nbins
        real, dimension(3), intent(in) :: box
        real, dimension(:,:), intent(in) :: r
        integer, dimension(:,:,:), intent(out) :: head
        integer, dimension(:), intent(out) :: list

        integer :: i, j
        integer, dimension(3) :: ir


        head=0
        list=0
        Do i=1, size(r,1)
            ! Find cell indices
            ir(:) = assign_bins(r(i,:), nbins, box)
            ! Build linked list
            list(i) = head(ir(1), ir(2), ir(3))
            ! Update Head List
            head(ir(1), ir(2), ir(3)) = i
        EndDo
    end subroutine build_linked_list

    subroutine cell_list_distance(r1, r2, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2)
    ! This subroutine calculates the distance between all pairs of atoms in a system
    ! using a cell linked list. 
    !
    ! Inputs:
    !   r1: 3D array of atom positions
    !   r2: 3D array of atom positions
    !   box: 3D array of box lengths
    !   cell_length: Length of the cell
    !   rc_sq: Square of the distance cutoff

    !  
    ! Outputs:
    !   dr_values: Array of distances
    !   dr_atom1: Array of atom indices
    !   dr_atom2: Array of atom indices

        implicit none
        

        real, intent(in) :: rc_sq, cell_length
        real, dimension(:,:), intent(in) :: r1, r2
        real, dimension(:), intent(in) :: box
        real, dimension(:), intent(out) :: dr_values
        integer, dimension(:), intent(out) :: dr_atom1, dr_atom2

        integer :: i, j, k, l, m, n
        integer :: ii, jj, kk
        integer :: di
        integer :: ihead, jhead
        integer :: count

        integer, dimension(3) :: ir
        real, dimension(3) :: dr_tmp
        real :: rsq

        integer, dimension(3,500) :: map
        

        
        ! Calculate the number of bins in each direction
        map = 0
        do i=1,3
            nbins(i) = max(1, int(box(i)/cell_length))
            if ( nbins(i) > 500 ) then
                write(*,*) "Error: Number of bins in direction ", i, " is greater than 500"
                stop
            endif
            map(i,1) = nbins(i)-1
            map(i,2) = nbins(i)
            do j=1, nbins(i)
                map(i,j+2) = j
            enddo
            map(i,nbins(i)+3) = 1
            map(i,nbins(i)+4) = 2
        enddo
        write(*,*) nbins(:)


        allocate(head_r1(nbins(1), nbins(2), nbins(3)))
        allocate(head_r2(nbins(1), nbins(2), nbins(3)))

        ! Initialize the Arrays
        head_r1 = 0
        head_r2 = 0
        list_r1 = 0
        list_r2 = 0
        dr_values = 0.0
        dr_atom1 = 0
        dr_atom2 = 0
        
        ! Build cell linked list
        call build_linked_list(r1, nbins, box, head_r1, list_r1)
        call build_linked_list(r2, nbins, box, head_r2, list_r2)


        count = 0
        ! Loop over Central Bins
        Do i=1, nbins(1)
        Do j=1, nbins(2)
        Do k=1, nbins(3)
            ! Loop over 125 nearest neighbor cells
            Do l=-2, 2
            Do m=-2, 2
            Do n=-2, 2
                ! Ensure that there is no double counting of cells
                    ii = map(1,i+l+2)
                    jj = map(2,j+m+2)
                    kk = map(3,k+n+2)
                    ! Head of current cell
                    ihead = head_r1(i,j,k) 
                    Do While (ihead /= 0)
                        ! Loop over atoms in the neighboring cell
                        jhead = head_r2(ii, jj, kk)
                        Do While (jhead /= 0 .and. jhead /= ihead)
                            dr_tmp = 0.0
                            Do di=1,3
                                dr_tmp(di) = r1(ihead,di) - r2(jhead,di)
                                dr_tmp(di) = dr_tmp(di) - box(di)*nint(dr_tmp(di)/box(di))
                            EndDo
                            rsq = sum(dr_tmp(:)**2)
                            ! Distance cutoff
                            If (rsq < rc_sq .and. ihead < jhead) Then
                                count = count + 1
                                dr_values(count) = rsq
                                dr_atom1(count) = ihead
                                dr_atom2(count) = jhead
                            EndIf
                            ! Move to the next atom in neighboring cell
                            jhead = list_r2(jhead)
                        EndDo
                        ! Move to nxt atom in current cell
                        ihead = list_r1(ihead)
                    EndDo
                !EndIf
            EndDo !n
            EndDo !m
            EndDo !l
        EndDo !k
        EndDo !j
        EndDo !i
        write(*,*) count
        deallocate(head_r1)
        deallocate(head_r2)
        
    end subroutine cell_list_distance

    subroutine naive_distance(r1, r2, box, rc_sq, dr_values, dr_atom1, dr_atom2, cell_assign_1, cell_assign_2)
    ! This subroutine calculates the distance between all pairs of atoms in a system
        real, intent(in) :: rc_sq
        real, dimension(:,:), intent(in) :: r1, r2
        real, dimension(:), intent(in) :: box
        real, dimension(:), intent(out) :: dr_values
        integer, dimension(:), intent(out) :: dr_atom1, dr_atom2
        integer, dimension(:,:) :: cell_assign_1, cell_assign_2

        
        integer :: i, j, di
        integer :: count
        integer, dimension(3) :: nbins

        real, dimension(3) :: dr_tmp
        real :: rsq
        
        nbins(1) = 58
        nbins(2) = 58
        nbins(3) = 57

        dr_values = 0.0
        dr_atom1 = 0
        dr_atom2 = 0
        count = 0
        Do i=1, size(r1,1)
            Do j=i+1, size(r2,1)
                dr_tmp = 0.0
                Do di=1,3
                    dr_tmp(di) = r1(i,di) - r2(j,di)
                    dr_tmp(di) = dr_tmp(di) - box(di)*nint(dr_tmp(di)/box(di))
                EndDo
                rsq = sum(dr_tmp(:)**2)
                If (rsq < rc_sq .and. i /= j) Then
                    count = count + 1
                    dr_values(count) = rsq
                    dr_atom1(count) = i
                    dr_atom2(count) = j
                    cell_assign_1(count,:) = assign_bins(r1(i,:), nbins, box)
                    cell_assign_2(count,:) = assign_bins(r2(j,:), nbins, box)
                EndIf
            EndDo
        EndDo
        write(*,*) count
    
    end subroutine naive_distance
    
    subroutine test_and_time_distance(natoms, rc)
        implicit none
        integer, intent(in) :: natoms
        real, intent (in) :: rc

        integer :: i, j, ntimes

        real, dimension(3) :: box, rtmp
        real, dimension(natoms,3) :: r

        real :: cell_length, rc_sq

        real, dimension(natoms*500) :: dr_values_naive, dr_values
        integer, dimension(natoms*500) :: dr_atom1, dr_atom2, dr_atom1_naive, dr_atom2_naive
        integer, dimension(natoms*500,3) :: cell_assign_1, cell_assign_2

        real :: start_time, end_time, elapsed_time
        

        call random_seed()

        box(1) = 103
        box(2) = 102
        box(3) = 101

        do i=1, natoms
            rtmp = 0.0
            do j=1, 3
                call random_number(rtmp(j))
            enddo
            r(i,:) = rtmp(:)*box(:)
        end do

        cell_length = rc/2.0
        rc_sq = rc**.2

        ! Time loop
        ntimes = 5
        
        ! Cell-List Approach
        write(*,*) "Cell-List"
        call cpu_time(start_time)
        do i=1, ntimes
            call cell_list_distance(r, r, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2)
        enddo
        call cpu_time(end_time)
        elapsed_time = end_time - start_time

        write(*,*) "Elapsed time: ", elapsed_time/ntimes, " seconds"
        write(*,*) "For 50,000 frames: ", elapsed_time/ntimes*50000/3600, " hours"

        ! Naive Approach
        write(*,*) "Naive-Approach"
        call cpu_time(start_time)
        do i=1, ntimes
            call naive_distance(r, r, box, rc_sq, dr_values_naive, dr_atom1_naive, dr_atom2_naive, cell_assign_1, cell_assign_2)
        enddo
        call cpu_time(end_time)
        elapsed_time = end_time - start_time

        write(*,*) "Elapsed time: ", elapsed_time/ntimes, " seconds"
        write(*,*) "For 50,000 frames: ", elapsed_time/ntimes*50000/3600, " hours"

        open(10, file="map_reg.dat")
        open(11, file="map_naive.dat")
        do i=1, 1000
            write(10,*) dr_values(i), dr_atom1(i), dr_atom2(i), dr_values_naive(i), dr_atom1_naive(i), dr_atom2_naive(i)
            write(11,*) dr_atom1_naive(i), cell_assign_1(i,1), cell_assign_1(i,2), cell_assign_1(i,3)
            write(11,*) dr_atom2_naive(i), cell_assign_2(i,1), cell_assign_2(i,2), cell_assign_2(i,3)
        enddo
        close(10)
        close(11)


    end subroutine test_and_time_distance

    subroutine test_boundary_handling()
        implicit none
        

        real :: cell_length
        real, dimension(3) :: box
        integer, dimension(3,500) :: map

        integer :: i, j, k, ii, jj, kk, l, m, n
        integer :: count


        box(1) = 103
        box(2) = 102
        box(3) = 101
        cell_length = 16.0/2.0

        map = 0
        do i=1,3
            nbins(i) = max(1, int(box(i)/cell_length))
            map(i,1) = nbins(i)-1
            map(i,2) = nbins(i)
            do j=1, nbins(i)
                map(i,j+2) = j
            enddo
            map(i,nbins(i)+3) = 1
            map(i,nbins(i)+4) = 2
        enddo

        

        write(*,*) nbins(:)
        count = 0
         ! Loop over Central Bins
        do i=1, nbins(1)
        do j=1, nbins(2)
        do k=1, nbins(3)
            ! Loop over 24 nearest neighbor cells
            Do l=0, 2
            Do m=-2, 2
            Do n=-2, 2
                ! Ensure that there is no double counting of cells
                    ii = map(1,i+l+2)
                    jj = map(2,j+m+2)
                    kk = map(3,k+n+2)
                    if ((ii == 4 .and. jj == 4 .and. kk == 4) .or. (i == 4 .and. j == 4 .and. k==4)) then
                        count = count + 1
                        write(*,*) i, j, k, ii, jj, kk
                        !write(*,*) i, j, ii, jj
                    endif
            EndDo !n
            EndDo !m
            EndDo !l
        enddo
        enddo
        enddo
            write(*,*) "count", count
            write(*,*) "nbins", nbins(1)*nbins(2)*nbins(3)
        j=12
        m=2
        write(*,*) j+m-1, (j+m-1)/nbins(2), mod(j+m-1, nbins(2)) + 1 - (nbins(2)*(int((j+m-1)/nbins(2))-1))
        write(*,*) map(2,j+m+2)
        j=1
        m=-2
        write(*,*) j+m-1, (j+m-1)/nbins(2), mod(j+m-1, nbins(2)) + 1 - (nbins(2)*(int((j+m-1)/nbins(2))-1))
        write(*,*) map(2,j+m+2)
        
    end subroutine test_boundary_handling

end module shared

program new_distance
    use shared
    implicit none
    integer :: i
    integer, dimension(6) :: natoms
    real :: rc

    rc = 3.5
    
    natoms(1) = 10
    Do i=2, 5
        natoms(i) = natoms(i-1)*10
    EndDo

    !call test_boundary_handling()

    Do i=1,6
        write(*,*) natoms(i), " atoms"
        call test_and_time_distance(natoms(i), rc)
    EndDo
    

end program new_distance
