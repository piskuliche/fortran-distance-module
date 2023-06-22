module mpi_distance

    use mpi

contains

    
    subroutine mpi_dist(comm, r1, r2, box, cell_length, rc_sq &
        , sparse_mat_values, sparse_mat_row_idx, sparse_mat_col_idx)
        ! **********************************************************************
        ! Calculates the distance between all atoms in r1 and r2, and returns the distances that are less than rc_sq
        ! in sparse matrix format. This subroutine is parallelized using MPI.
        ! Inputs:
        !   comm: MPI communicator
        !   r1: Coordinates of atoms in the first system
        !   r2: Coordinates of atoms in the second system
        !   box: Dimensions of the simulation box
        !   cell_length: Length of the cell in the cell list
        !   rc_sq: Square of the cutoff distance
        !
        ! Outputs:
        !   sparse_mat_values: Values of the sparse matrix
        !   sparse_mat_row_idx: Row indices of the sparse matrix
        !   sparse_mat_col_idx: Column indices of the sparse matrix
        ! **********************************************************************

        use distance_module
        use mpi

        implicit none

        ! Input variables ******************************************************
        integer, intent(in) :: comm
        real, dimension(:,:), intent(in) :: r1
        real, dimension(:,:), intent(in) :: r2
        real, dimension(3), intent(in) :: box
        real, intent(in) :: cell_length
        real, intent(in) :: rc_sq
        ! Output variables *****************************************************
        real, dimension(:), intent(out) :: sparse_mat_values
        integer, dimension(:), intent(out) :: sparse_mat_row_idx
        integer, dimension(:), intent(out) :: sparse_mat_col_idx
        ! Local variables ******************************************************
        integer :: i, j, k
        integer :: mpisize, rank, ierr
        integer :: start_index, end_index
        integer :: nnz_local, nnz_total
        integer :: num_atoms_per_proc
        integer :: size_sparse_mat



        real, allocatable :: sparse_mat_values_local(:)
        integer, allocatable :: sparse_mat_row_idx_local(:)
        integer, allocatable :: sparse_mat_col_idx_local(:)
        integer, allocatable :: nnz(:), displ(:)
        real, allocatable :: r1_local(:,:)
        real, allocatable :: dr_values_all(:)
        integer, allocatable :: dr_atom1_all(:)
        integer, allocatable :: dr_atom2_all(:)
        ! **********************************************************************

        ! Initialize MPI Variables 
        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, mpisize, ierr)
        
        ! Calculate the number of atoms per processor
        num_atoms_per_proc = size(r1, 1) / mpisize
        start_index = rank * num_atoms_per_proc + 1
        end_index = (rank + 1) * num_atoms_per_proc

        ! Set the size_sparse_mat to the number of elements in the sparse matrix
        size_sparse_mat = size(sparse_mat_values, 1)

        write(*,*) "Rank: ", rank, "Num atoms per proc: ", num_atoms_per_proc

        ! Allocate memory for the local variables ******************************
        allocate(r1_local(num_atoms_per_proc, 3))
        allocate(nnz(mpisize))
        allocate(displ(mpisize))
        if (rank == 0) then
            ! dr_values_all allocation
            allocate(dr_values_all(size_sparse_mat))
            allocate(dr_atom1_all(size_sparse_mat))
            allocate(dr_atom2_all(size_sparse_mat))
        endif

        ! sparse_mat allocation
        allocate(sparse_mat_values_local(size_sparse_mat/mpisize))
        allocate(sparse_mat_row_idx_local(size_sparse_mat/mpisize))
        allocate(sparse_mat_col_idx_local(size_sparse_mat/mpisize))
        ! **********************************************************************

        nnz = 0; r1_local = 0
        nnz_total = 0

        ! Calculate the non-zero elements in the local sparse matrix
        nnz = 0
        r1_local = r1(start_index:end_index, :)
        sparse_mat_values_local = 0.0
        sparse_mat_row_idx_local = 0
        sparse_mat_col_idx_local = 0

        call cell_list_distance(r1_local, r2, box, cell_length, rc_sq, &
            sparse_mat_values_local, sparse_mat_row_idx_local, sparse_mat_col_idx_local, same_array=1, offset=start_index)

        nnz_local = count(sparse_mat_values_local /= 0.0)
        write(*,*) nnz_local
        !nnz_local = 1250

        call MPI_Barrier(comm, ierr)

        ! Share the number of non-zero elements on each process
        call MPI_Allgather(nnz_local, 1, MPI_INTEGER, nnz, 1, MPI_INTEGER, comm, ierr)

        ! Calculate the displacement of this process's non-zero elements in the gathered array
        if (rank == 0) then
            displ(1) = 0
            do i = 2, mpisize
                displ(i) = displ(i-1) + nnz(i-1)
            enddo
        endif

        call MPI_Barrier(comm, ierr)

        ! Share the values of displ across all processes
        call MPI_Bcast(displ, mpisize, MPI_INTEGER, 0, comm, ierr)
        call MPI_Allgather(displ(rank+1), 1, MPI_INTEGER, displ, 1, MPI_INTEGER, comm, ierr)

        ! Synchronize all processes
        call MPI_Barrier(comm, ierr)

        ! Print out the values of displ and nnz on process 0
        if (rank == 0) then
            write(*,*) "displ = ", displ
            write(*,*) "nnz = ", nnz
        endif

        write(*,*) rank, size(sparse_mat_values_local,1)

        call MPI_Barrier(comm, ierr)


        call MPI_Gatherv(sparse_mat_values_local, nnz_local, MPI_REAL, &
            dr_values_all, nnz, displ, MPI_REAL, 0, comm, ierr)

        call MPI_Gatherv(sparse_mat_row_idx_local, nnz_local, MPI_INTEGER, &
            dr_atom1_all, nnz, displ, MPI_INTEGER, 0, comm, ierr)

        call MPI_Gatherv(sparse_mat_col_idx_local, nnz_local, MPI_INTEGER, &
            dr_atom2_all, nnz, displ, MPI_INTEGER, 0, comm, ierr)

        if (rank == 0) then

            sparse_mat_values = dr_values_all
            sparse_mat_row_idx = dr_atom1_all
            sparse_mat_col_idx = dr_atom2_all

            write(*,*) "Sparse matrix values: ", count( dr_values_all /= 0)

            deallocate(dr_values_all)
            deallocate(dr_atom1_all)
            deallocate(dr_atom2_all)
        endif

        deallocate(r1_local)


    end subroutine mpi_dist

end module mpi_distance

