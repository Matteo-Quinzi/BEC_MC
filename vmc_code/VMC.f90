program main
     use mpi 
     use io_handler
     use bec_vmc
 
     implicit none
     
     !mpi variables
     integer :: ierr
     integer :: nprocs, my_rank 
     
     !Input data
     integer(kind=8) :: N_at_row, N_at, N_walk  ! #atoms and #walkers
     integer(kind=8) :: M, N_cyc                ! #step per cycle and #cycles
     real(kind=8) :: a                          ! scattering length
     real(kind=8) :: b0, b1                     ! variational params
     real(kind=8) :: l                          ! gas radius
     integer(kind=8) :: from_file               ! tells if input is from file
     character(len=50) :: coords_input_file     

     !local variables (for each process)
     integer(kind=8) :: N_loc_walk    ! # of local walkers
     integer(kind=8) :: remainder     ! remainder of N_walk/nprocs
     real(kind=8), allocatable :: loc_E(:,:)      ! Local energy (rows are walkers, columns are timesteps)
     real(kind=8), allocatable :: start_coords(:,:)  !initial coordinates 
     real(kind=8), allocatable :: coords(:,:,:)
     real(kind=8) :: acc_prob
     integer(kind=8) :: N_mark
     real(kind=8), allocatable :: energies(:,:), mean_energies(:)
     real(kind=8) :: acc_ratio
     real(kind=8) :: temp_energy, temp_err_energy
     real(kind=8), allocatable :: master_temp_energy(:), master_acc_ratio(:), master_temp_err_energy(:)
     real(kind=8), allocatable :: out_energy(:), out_acc_ratio(:), out_err_energy(:)
     real(kind=8) :: tic, toc, start_tic

     !date and time values
     character(len=8) :: date
     character(len=10) :: time

     !loop indexes
     integer(kind=8) :: i, l1, l2, l3 

     !INITIALIZING MPI 
     call mpi_init(ierr)
     call mpi_comm_size(mpi_comm_world, nprocs, ierr)
     call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

     !Initialize random number generator
     call init_rng(my_rank)

     allocate(master_temp_energy(nprocs)) 
     allocate(master_acc_ratio(nprocs))
     allocate(master_temp_err_energy(nprocs))

     start_tic = mpi_wtime()

     call date_and_time(date,time)

     !Printing date and time 
     if (my_rank .eq. 0) then
         write(*,*)
         write(*,*) 
         write(*,*) 'VMC run started.'
         write(*,*) 'Date : ', date
         write(*,*) 'Time : ', time
         write(*,*) 
         write(*,*) 
         write(*,*) '---------------------------------------------------'
     end if

     ! Rank 0 works as the Master :
     ! It reads input file and broadcast data to all other ranks 
     if (my_rank .eq. 0) then 
         call read_input_file()
         call read_data(N_at, N_walk, M, N_cyc, a,&
                 b0, b1, l, from_file, coords_input_file)
         allocate(out_energy(N_cyc), out_acc_ratio(N_cyc), out_err_energy(N_cyc))
     end if 


     !All procs must make the call to bcast
     call mpi_bcast(N_at, 1, mpi_integer8, 0, mpi_comm_world, ierr)
     call mpi_bcast(N_walk, 1, mpi_integer8, 0, mpi_comm_world, ierr)
     call mpi_bcast(M, 1, mpi_integer8, 0, mpi_comm_world, ierr)
     call mpi_bcast(N_cyc, 1, mpi_integer8, 0, mpi_comm_world, ierr)
     call mpi_bcast(a, 1, mpi_real8, 0, mpi_comm_world, ierr)
     call mpi_bcast(b0, 1, mpi_real8, 0, mpi_comm_world, ierr)
     call mpi_bcast(b1, 1, mpi_real8, 0, mpi_comm_world, ierr)
     call mpi_bcast(l, 1, mpi_real8, 0, mpi_comm_world, ierr)

     !Printing info about the run
     if (my_rank .eq. 0) then
             write(*,*) 
             write(*,*)
             write(*,'(a27,4x,i10)')   'Number of processes      : ', nprocs
             write(*,'(a27,4x,i10)' )  'Number of atoms          : ', N_at
             write(*,'(a27,4x,i10)' )  'Number of walkers        : ', N_walk
             write(*,'(a27,4x,i10)' )  'Number of step per cycle : ', M
             write(*,'(a27,4x,i10)' )  'Number of cycles         : ', N_cyc
             write(*,'(a27,4x,f6.5)')  'Scattering length        : ', a
             write(*,'(a27,4x,f10.6)') 'Variational parameter b0 : ', b0
             write(*,'(a27,4x,f10.6)') 'Variational parameter b1 : ', b1
             write(*,'(a27,4x,f7.5)')  'Max initialization radius: ', l
             write(*,'(a27,4x,i1)')    'Reading from file        : ',from_file
             write(*,'(a27,4x,a50)')   'Coords Input file        : ',coords_input_file
             write(*,*)
             write(*,*)
             write(*,*) '------------------------------------------------------'
     end if

     ! Evaluating the # of walkers for each process
     N_loc_walk = N_walk/nprocs
     remainder = mod(N_walk,nprocs)
     if (my_rank < remainder) N_loc_walk = N_loc_walk + 1

     !Allocating coords 
     allocate(start_coords(N_at,3))
     allocate(coords(N_loc_walk, N_at, 3))

     !Setting up the condensate 
     !Master rank determines the initial coordinates and broadcast them to the 
     !other ranks
     if (my_rank .eq. 0) then
             if (from_file .eq. 0) then
                 start_coords = set_up_boson_gas(N_at, a, l)
                 call save_coords(N_at, start_coords, 'start_coords.txt')
             else if (from_file .eq. 1) then
                 start_coords = set_up_from_file(N_at, coords_input_file)
                 call save_coords(N_at, start_coords, 'start_coords.txt')
             else
                 write(*,*) 'Somethin went wrong with the initial coordinates input ...'
                 stop
             end if    
     endif

     call mpi_bcast(start_coords(:,1), N_at, mpi_real8, 0, mpi_comm_world, ierr)
     call mpi_bcast(start_coords(:,2), N_at, mpi_real8, 0, mpi_comm_world, ierr)
     call mpi_bcast(start_coords(:,3), N_at, mpi_real8, 0, mpi_comm_world, ierr)

     !Initial coordinates are copied for all walkers 
     do i = 1, N_loc_walk
        coords(i,:,:) = start_coords 
     end do

     deallocate(start_coords)
     allocate(mean_energies(M))
      
     do i = 1, N_cyc

         tic = mpi_wtime()

         ! Metropolis run for all walkers 
         ! The outcoming energy is averaged over all walkers
         call metropolis(a, b0, b1, N_at, N_loc_walk, M, coords, mean_energies, acc_ratio)

         !Each process average the energy over all steps
         temp_energy = sum(mean_energies/M)
         temp_err_energy = err_energy(M, mean_energies, temp_energy) 

         !Each process sends the energy to the Master 
         call mpi_gather(temp_energy, 1, mpi_real8, &
                         master_temp_energy, 1, mpi_real8, &
                         0, mpi_comm_world, ierr)
         call mpi_gather(acc_ratio, 1, mpi_real8, &
                         master_acc_ratio, 1, mpi_real8, &
                         0, mpi_comm_world, ierr)
         call mpi_gather(temp_err_energy, 1, mpi_real8, &
                         master_temp_err_energy, 1, mpi_real8, &
                         0, mpi_comm_world, ierr)

         toc = mpi_wtime()

         ! Master average the temp_energy and returns the average value of the local energy
         if (my_rank .eq. 0) then
                 out_energy(i) = sum(master_temp_energy/nprocs)
                 out_acc_ratio(i) = sum(master_acc_ratio/nprocs)
                 out_err_energy(i) = sum(master_temp_err_energy/nprocs)
                 write(*,*) 
                 write(*,'(a30,4x,i10)')   'Terminated cycle             :', i
                 write(*,'(a30,4x,f15.8)') 'Average value of the energy  :', out_energy(i)
                 write(*,'(a30,4x,f15.8)') 'Average error on the energy  :', out_err_energy(i)
                 write(*,'(a30,4x,f15.8)') 'Average acceptance ratio     :', out_acc_ratio(i)
                 write(*,'(a30,4x,f15.8)') 'Time spent within this cycle :', toc-tic
                 write(*,*)
                 write(*,*)
                 write(*,*) '---------------------------------------------------------------'
         end if

     end do

     ! Save final coordinates of one walker
     ! Can be used to reinitialize another run
     if (my_rank == 0) then 
             call save_coords(N_at, coords(1,:,:), 'final_coords.txt')
     end if

     ! Save energy evolution of the last cycle
     call save_rank_energies(my_rank, M, mean_energies)

     ! Printout the final value of the energy 
     if (my_rank .eq. 0) then
              write(*,*) 
              write(*,*) 'Program terminated'
              write(*,'(a6,4x,f15.8)') '!E   :', out_energy(N_cyc)
              write(*,'(a6,4x,f15.8)') '!!err:', out_err_energy(N_cyc)
              write(*,'(a6,4x,f15.8)') '!!a_r:', out_acc_ratio(N_cyc)
              write(*,'(a37,4x,f15.5)') 'Time spent within the whole program :', toc - start_tic 
              write(*,*)
              write(*,*)
              write(*,*) '---------------------------------------------------------------'
              call save_output(N_cyc, out_energy, out_err_energy, out_acc_ratio)
     end if
     
     call mpi_finalize(ierr)

end program main 
