program trial
        use bec_dmc
        implicit none
        real(kind=8), allocatable :: trial_array(:,:,:)
        integer :: N_walk, N_at, N_dim
        integer :: i,j,k
        integer :: io_unit=11

        N_walk=1
        N_at = 127
        N_dim = 3

        allocate(trial_array(N_walk,N_at,N_dim))

        print *, shape(trial_array)
        trial_array = gaussian_rng(N_walk,N_at,N_dim)

        open(unit=io_unit,file='rng.txt',action='Write')
        do i =1,N_walk
            do j =1,N_at
                do k =1,N_dim
                   write(io_unit,*) trial_array(i,j,k)
                end do
            end do
        end do
        close(io_unit)

        deallocate(trial_array)

end program trial
