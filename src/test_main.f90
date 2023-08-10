Program test
    use gmxfort_trajectory
    implicit none
    type(Trajectory) :: trj
    integer :: ntmp

    call trj%open('../../test_files/temp.xtc', "../..//index.ndx")
    ntmp = trj%read_next(1)
    write(*,*) ntmp
    write(*,*) trj%box(1)
    write(*,*) trj%x(1,1)
    call trj%close()

End Program test