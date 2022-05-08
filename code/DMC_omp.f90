program DMC_omp
      !This is a OpenMP version of the Diffusion Monte Carlo
      ! program for the study of a trapped boson gas.
      !As the number of walkers is varying and we don't want the
      ! processes to remain stucked with different numbers of walkers
      ! an OMP implementation is carried over.
      use io_handler
      use bec_vmc
      use bec_dmc
      use omp_lib
      implicit none

      
      !Date and Time
      character(len=8) :: date
      character(len=10) :: time

      !DMC variables 
      integer(kind=8) :: N_at, N_walk, N_max      !N_atoms, walkers and maximum number of walkers
      integer(kind=8) :: eq_it, samples, dt_sam   !N of iterations fir equilibration, samples and delta_time between samples
      real(kind=8)    :: dt                       !timestep width
      real(kind=8)    :: Er                       !Trial energy
      real(kind=8)    :: a, b0, b1                !Scattering length and params of the guiding function
      character(len=50) :: coords_input_file
      real(kind=8), allocatable :: acc_prob(:)            !Overall acc prob for VMC set up
      real(kind=8), allocatable :: configurations(:,:,:)  !Atoms configurations
      logical, allocatable :: flag(:)

      !Initial configuration variables
      real(kind=8) :: acc_prob_vmc
      real(kind=8), allocatable :: temp_prob(:)
      real(kind=8), allocatable :: temp_coords(:,:)
      integer(kind=8) :: metro_step

      !Equilibration variables
      integer(kind=8), allocatable :: Nt(:)       !Number of walkers at each step
      real(kind=8), allocatable :: walk_en(:), walk_en_old(:)     !Energy for each walker at a given step
      real(kind=8), allocatable :: Et(:)          !Average energy at each step
      real(kind=8), allocatable :: F_driv(:,:,:)  !Driving force for all walkers (N_max:N_at:3)
      real(kind=8) :: acc_prob_eq                 !Acceptance probability in equilibration procedure
      real(kind=8), allocatable :: temp_prob_eq(:)
      logical :: eq_accepted
      integer(kind=8) :: N_walk_start
      real(kind=8) :: N_ratio, alpha

      !CPU time variables
      real(kind=8) :: tic, toc

      !Loops indexes
      integer(kind=8) :: i, it

      !INITIALIZING RNG
      call fix_rng(0)

      !PRINTING DATE AND TIME !
      call date_and_time(date, time)
      write(*,*)
      write(*,*)
      write(*,*) 'DMC run started :'
      write(*,*) 'Date : ', date(:4),'/',date(5:6),'/',date(7:8)
      write(*,*) 'Time : ', time(:2),':',time(3:4),':',time(5:6)
      write(*,*)
      write(*,*) 
      write(*,*) '------------------------------------'

      !READING INPUT DATA
      !Reading the name of the input file
      call read_input_file()
      !Reading input data
      call read_data_dmc(N_at, N_walk, N_max, eq_it, samples, &
                     dt_sam, dt, Er, a, b0, b1, coords_input_file)
      !Sanity check
      if (N_walk > N_max) then
          write(*,*) 'ERROR : N_walk can not be greater than N_max !!!'
          stop
      end if

      !Printing out info about the run
      write(*,*)
      write(*,'(a27,4x,i10)')    'Number of atoms          : ', N_at
      write(*,'(a27,4x,i10)')    'Number of walkers        : ', N_walk
      write(*,'(a27,4x,i10)')    'Max num of walkers       : ', N_max
      write(*,'(a27,4x,i10)')    'Number of eq. iterations : ', eq_it
      write(*,'(a27,4x,i10)')    'Number of samples        : ', samples
      write(*,'(a27,4x,i10)')    'Steps between two samples: ', dt_sam
      write(*,'(a27,4x,f15.10)') 'Timestep                 : ', dt
      write(*,'(a27,4x,f15.10)') 'Energy scale             : ', Er
      write(*,'(a27,4x,f15.10)') 'Scattering length        : ', a
      write(*,'(a27,4x,f15.10)') 'Guiding func. param b0   : ', b0
      write(*,'(a27,4x,f15.10)') 'Guiding func. param b1   : ', b1
      write(*,'(a27,a50)')       'Coordinates input file   : ', coords_input_file

      !ALLOCATING CONFIGURATIONS ARRAYS
      allocate(configurations(N_max, N_at, 3))
      allocate(flag(N_max))
      allocate(temp_coords(N_at,3))

      !READING INITIAL COORDINATES FILE
      temp_coords = set_up_from_file(N_at, coords_input_file)


      tic = omp_get_wtime()
      !CREATING INITIAL COORDINATES DISTRIBUTION WITH A VMC PROCEDURE
      metro_step = 1000*N_at
      acc_prob_vmc = 0.d0
      allocate(temp_prob(N_walk))

      !$OMP PARALLEL DO  
      do i = 1, N_walk
            flag(i) = .True.     !Walker is set to alive
            if (i .eq. 1) then 
                    configurations(1,:,:) = temp_coords
            else
                    configurations(i,:,:) = configurations(i-1,:,:)
            end if
            !Metropolis run is used to change coordinates
            call no_energy_metropolis(a, b0, b1, N_at, 1, metro_step, configurations(i,:,:), temp_prob(i))
      end do
      !$OMP END PARALLEL DO
      acc_prob_vmc = sum(temp_prob) / N_walk

      do i = N_walk + 1, N_max
                flag(i) = .False.    !Walker is set to dead
      end do
      toc = omp_get_wtime()
      deallocate(temp_coords)

      !DIAGNOSTICS ABOUT INITIAL CONFIGURATION
      write(*,*)
      write(*,*) '-----------------------------------------------'
      write(*,*) 'Initial Distribution setted '
      write(*,*) 'CPU time : ', toc - tic
      write(*,*) 'Typical acceptance ratio : ', acc_prob_vmc
      write(*,*)
      write(*,*)

      !EQUILIBRATION
      !All major loops are contained within the main program
      !The most expensive subroutine is the diffusion one, which
      ! evaluates both driving forces and local energies. It is this one
      ! that I want to parallelize.
      allocate(Nt(eq_it), Et(eq_it))
      allocate(walk_en(N_max), walk_en_old(N_max))
      allocate(F_driv(N_max,N_at,3))

      !Evaluating initial driving forces and initial local energies
      !$OMP PARALLEL DO
      do i = 1, N_walk
          F_driv(i,:,:) = driving_force(a, b0, b1, N_at, configurations(i, :, :))
          walk_en(i)    = local_energy(a, b0, b1, N_at, configurations(i,:,:)) 
      end do
      !$OMP END PARALLEL DO

      acc_prob_eq = 0.d0
      allocate(temp_prob_eq(N_max))

      !Saving initial number of walker
      N_walk_start = N_walk

      !Fixing population control parameter
      !alpha = N_at * N_walk_start / dt
      alpha = 1.0d0 / dt

      tic = omp_get_wtime()
      !Equilibration loop
      do it = 1, eq_it

          !Diffusing alive walkers
          !$OMP PARALLEL DO reduction(+:acc_prob_eq)
          do i = 1, N_walk
              !Saving starting walker energy
              walk_en_old(i) = walk_en(i)
              !Each walker diffuse with accept/reject step
              call one_walker_diffusion(a, b0, b1, N_at, dt, configurations(i, :, :), F_driv(i, :, :), walk_en(i), eq_accepted)
              if (eq_accepted .eq. .True.) acc_prob_eq = acc_prob_eq + 1.d0/(N_walk*it)
          end do
          !$OMP END PARALLEL DO

          !Branching
          !N_walk, walk_en, F_driv are going to be updated with the number of replicas
          call branching(N_walk, N_max, N_at, configurations, walk_en, walk_en_old, F_driv, flag, dt, Er)

          !Saving number of walkers at step it 
          Nt(it) = N_walk
          !Evalauting average energy at step it
          Et(it) = sum(walk_en(1:N_walk)) / N_walk

          !Adjusting the energy scale
          if ( mod(it,100) .eq. 0) then
              N_ratio = (N_walk_start*1.d0) / (Nt(it)*1.d0)
              Er = Et(it) + alpha * log(N_ratio)
          end if

          print *, Nt(it), Et(it)

      end do
      toc = omp_get_wtime()
      print *, toc-tic

      deallocate(temp_prob_eq)






      deallocate(Nt, Et)
      deallocate(walk_en, walk_en_old)
      deallocate(F_driv)


      deallocate(configurations)
      deallocate(flag)
end program DMC_omp
