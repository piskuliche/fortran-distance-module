program test_mpi
    use mpi
    use distance_module
    use mpi_distance
    
    implicit none

    integer, parameter :: natoms = 10000
    integer :: i, j, k
    integer :: ierr, rank, mpisize
    real :: start_time, end_time
    real, dimension(2) :: elapsed_time


    real, dimension(3) :: box 
    real, dimension(3) :: rtmp
    real :: cell_length
    real :: rc_sq
    real, dimension(natoms*5) :: dr_values
    integer, dimension(natoms*5) :: dr_atom1
    integer, dimension(natoms*5) :: dr_atom2

    real, dimension(natoms,3) :: r

    integer, dimension(1,3) :: cell_assign_1, cell_assign_2

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mpisize, ierr)

    box=100
    rc_sq = 16
    cell_length = 2.0
    dr_values = 0.0
    dr_atom1 = 0; dr_atom2 = 0

    if ( rank == 0 ) then
        write(*,*) "Number of MPI processes: ", mpisize
        write(*,*) "Number of atoms: ", natoms
        
        ! Generate random coordinates
        do j=1, natoms
            rtmp = 0.0
            do k=1, 3
                call random_number(rtmp(k))
            enddo
            r(j,:) = rtmp(:)*box(:)
        end do

        ! Cell-List Approach ******************************************************
        call cpu_time(start_time)
        
        ! Send r to all processes

    endif
    call MPI_Bcast(r, natoms*3, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    ! Wait for all processes to reach this point
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Call the MPI distance function
    call mpi_dist(MPI_COMM_WORLD, r, r, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2)

    if ( rank == 0 ) then
        call cpu_time(end_time)
        elapsed_time(1) = end_time - start_time
        write(*,*) "Cell-List Approach: ", elapsed_time(1)
        write(*,*) "Finished Program"
    EndIf
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        dr_values = 0.0
        dr_atom1 = 0; dr_atom2 = 0
        call cell_list_distance(r, r, box, cell_length, rc_sq, dr_values, dr_atom1, dr_atom2, same_array=1)
        write(*,*) count(dr_values>0)
        dr_values = 0.0
        dr_atom1 = 0; dr_atom2 = 0
        call double_loop_distance(r, r, box, rc_sq &
        , dr_values, dr_atom1, dr_atom2, cell_assign_1, cell_assign_2, compare_cell=0, same_array=1, cell_length=cell_length)
        write(*,*) count(dr_values>0)
    end if

    call MPI_Finalize(ierr)

end program test_mpi