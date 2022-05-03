program DMC_BEC
      ! use ulimit -s unlimited in the bash to avoid segmentation faults
      use mpi
      use io_handler 
      use bec_dmc
      use bec_vmc
      implicit none

      !mpi_variables
      integer :: ierr
      integer :: nprocs, my_rank

      !date and time 
      character(len=8) :: date
      character(len=10) :: time

      !DMC variables
      integer(kind=8) :: N_at, N_walk, N_max
      integer(kind=8) :: eq_it, samples, dt_sam
      real(kind=8)    :: dt

      !INITIALIZING MPI
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

      !INITIALIZING RNG
      call init_rng(my_rank)

      !PRINTING DATE AND TIME
      if (my_rank .eq. 0) then
         call date_and_time(date,time)
         write(*,*)
         write(*,*) 
         write(*,*) 'DMC run started.'
         write(*,*) 'Date : ', date
         write(*,*) 'Time : ', time
         write(*,*) 
         write(*,*) 
         write(*,*) '---------------------------------------------------'
      end if

      !READING NAME OF INPUT FILE
      if (my_rank .eq. 0) then
              call read_input_file()
              call read_data_dmc(N_at, N_walk, N_max, eq_it, samples, dt_sam, dt)
      end if

      !INFORMATIVE PRINTOUT
      if (my_rank .eq. 0) then
             write(*,*)
             write(*,'(a27,4x,i10)')    'Number of processes      : ', nprocs
             write(*,'(a27,4x,i10)')    'Number of atoms          : ', N_at
             write(*,'(a27,4x,i10)')    'Number of walkers        : ', N_walk
             write(*,'(a27,4x,i10)')    'Max num of walkers       : ', N_max      
             write(*,'(a27,4x,i10)')    'Number of eq. iterations : ', eq_it
             write(*,'(a27,4x,i10)')    'Number of samples        : ', samples
             write(*,'(a27,4x,i10)')    'Steps between two samples: ', dt_sam
             write(*,'(a27,4x,f15.10)') 'Timestep                 : ', dt
      end if
     
      call mpi_finalize(ierr)
end program DMC_BEC

