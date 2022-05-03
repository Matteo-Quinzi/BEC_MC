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

      !DMC variables (global)
      integer(kind=8) :: N_at, N_walk, N_max     !N_atoms, walkers and maximum number of walkers
      integer(kind=8) :: eq_it, samples, dt_sam  !N of iterations fir equilibration, samples and delta_time between samples
      real(kind=8)    :: dt                      !timestep width
      real(kind=8)    :: a, b0, b1 !Scattering length and params of the guiding function
      character(len=50) :: coords_input_file
      real(kind=8), allocatable :: acc_prob(:)   !Overall acc prob for VMC set up

      !DMC variables (local -- for one process)
      integer(kind=8) :: my_walk, my_max !initial number of walkers and max num of walks per process
      real(kind=8), allocatable :: my_configurations(:,:,:)   !max_walkers:atoms:coords
      logical, allocatable :: my_flag(:) !Existance flag for walkers
      real(kind=8), allocatable :: temp_coords(:,:)
      integer(kind=8) :: metro_step
      real(kind=8), allocatable :: my_acc_prob, temp_prob !acceptance probability for VMC in initial condition

      !Loop indexes
      integer(kind=8) :: i

      !CPU-Time variables
      real(kind=8) :: tic, toc  


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
             call read_data_dmc(N_at, N_walk, N_max, eq_it, samples, &
                     dt_sam, dt, a, b0, b1, coords_input_file)
             !Sanity check
                 if (N_walk > N_max) then
                     write(*,*) 'ERROR : N_walk can not be greater than N_max !!!'
                     stop
                 end if
             write(*,*)
             write(*,'(a27,4x,i10)')    'Number of processes      : ', nprocs
             write(*,'(a27,4x,i10)')    'Number of atoms          : ', N_at
             write(*,'(a27,4x,i10)')    'Number of walkers        : ', N_walk
             write(*,'(a27,4x,i10)')    'Max num of walkers       : ', N_max      
             write(*,'(a27,4x,i10)')    'Number of eq. iterations : ', eq_it
             write(*,'(a27,4x,i10)')    'Number of samples        : ', samples
             write(*,'(a27,4x,i10)')    'Steps between two samples: ', dt_sam
             write(*,'(a27,4x,f15.10)') 'Timestep                 : ', dt
             write(*,'(a27,4x,f15.10)') 'Scattering length        : ', a
             write(*,'(a27,4x,f15.10)') 'Guiding func. param b0   : ', b0
             write(*,'(a27,4x,f15.10)') 'Guiding func. param b1   : ', b1
             write(*,'(a27,a50)')       'Coordinates input file   : ', coords_input_file
      end if

      !BROADCASTING INPUT DATA TO ALL PROCS
      call mpi_bcast(N_at, 1, mpi_integer8, 0, mpi_comm_world, ierr)
      call mpi_bcast(N_walk, 1, mpi_integer8, 0, mpi_comm_world, ierr)
      call mpi_bcast(N_max, 1, mpi_integer8, 0, mpi_comm_world, ierr)
      call mpi_bcast(eq_it, 1, mpi_integer8, 0, mpi_comm_world, ierr)
      call mpi_bcast(samples, 1, mpi_integer8, 0, mpi_comm_world, ierr)
      call mpi_bcast(dt_sam, 1, mpi_integer8, 0, mpi_comm_world, ierr)
      call mpi_bcast(dt, 1, mpi_real8, 0, mpi_comm_world, ierr)
      call mpi_bcast(a, 1, mpi_real8, 0, mpi_comm_world, ierr)
      call mpi_bcast(b0, 1, mpi_real8, 0, mpi_comm_world, ierr)
      call mpi_bcast(b1, 1, mpi_real8, 0, mpi_comm_world, ierr)

      !ALLOCATING NUMBER OF WALKER TO EACH PROCESS
      my_walk = N_walk / nprocs
      if (my_rank < mod(N_walk,nprocs)) my_walk = my_walk + 1
      my_max = N_max / nprocs
      if (my_rank < mod(N_max,nprocs)) my_max = my_max + 1

      !ALLOCATING CONFIGURATIONS ARRAY
      allocate(my_configurations(my_max,N_at,3))
      allocate(my_flag(my_max))
      allocate(temp_coords(N_at,3))

      !MASTER RANK READS EQ. COORDS AND BROADCAST THEM TO THE OTHERS
      if (my_rank .eq. 0) then
             allocate(acc_prob(nprocs))
             temp_coords = set_up_from_file(N_at, coords_input_file) 
      end if 
      do i=1,3
          call mpi_bcast(temp_coords(:,i), N_at, mpi_real8, 0, mpi_comm_world, ierr)
      end do


      tic = mpi_wtime()
      !SETTING THE INITIAL CONDITION
      metro_step = 1000*N_at    !Each atom is displaced (on average) at least once
      temp_prob = 0.d0
      do i = 1, my_max
              if (i .le. my_walk) then
                  my_flag(i) = .True.    !Walker is set to alive
                  !This is a trick because in metropolis I need coords to be (1,N_at,3)
                      if (i .eq. 1) then 
                          my_configurations(i,:,:) = temp_coords
                      else
                          my_configurations(i,:,:) = my_configurations(i-1,:,:)
                      end if
                  !Metropolis run is used to change coordinates
                  ! I need just one walker in the run (I am not doing energy averages)!!
                  call no_energy_metropolis(a, b0, b1, N_at, 1, metro_step, my_configurations(i,:,:), temp_prob)
                  my_acc_prob = my_acc_prob + temp_prob/my_walk
              else 
              my_flag(i) = .False.  ! Walker is set to dead
          end if
      end do

      call mpi_gather(my_acc_prob, 1, mpi_real8, &
                      acc_prob, 1, mpi_real8, &
                      0, mpi_comm_world, ierr)
      deallocate(temp_coords) ! I need this no more

      toc = mpi_wtime()

      !DIAGNOSTICS ABOUT INITIAL CONFIGURATION
      if (my_rank .eq. 0) then
          write(*,*)
          write(*,*) '-----------------------------------------------'
          write(*,*) 'Initial Distribution setted '
          write(*,*) 'CPU time : ', toc - tic
          write(*,*) 'Typical acceptance ratio : ', sum(acc_prob)/nprocs
          write(*,*)
      end if 
      
      !DEALLOCATING BEFORE EXITING 
      deallocate(my_configurations)
      deallocate(my_flag)
      call mpi_finalize(ierr)
end program DMC_BEC

