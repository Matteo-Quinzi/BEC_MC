program DMC_omp
      !This is a OpenMP version of the Diffusion Monte Carlo
      ! program for the study of a trapped boson gas.
      !As the number of walkers is varying and we don't want the
      ! processes to remain stucked with different numbers of walkers
      ! an OMP implementation is carried over.
      use io_handler
      use bec_vmc
      use bec_dmc
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
      real(kind=8) :: acc_prob_vmc, temp_prob
      real(kind=8), allocatable :: temp_coords(:,:)
      integer(kind=8) :: metro_step

      !CPU time variables
      real(kind=8) :: tic, toc

      !Loops indexes
      integer(kind=8) :: i

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


      call cpu_time(tic)
      !CREATING INITIAL COORDINATES DISTRIBUTION WITH A VMC PROCEDURE
      metro_step = 1000*N_at
      acc_prob_vmc = 0.d0

      !$OMP PARALLEL DO
      do i = 1, N_walk
            flag(i) = .True.     !Walker is set to alive
            if (i .eq. 1) then 
                    configurations(1,:,:) = temp_coords
            else
                    configurations(i,:,:) = configurations(i-1,:,:)
            end if
            !Metropolis run is used to change coordinates
            call no_energy_metropolis(a, b0, b1, N_at, 1, metro_step, configurations(i,:,:), temp_prob)
            !acc_prob_vmc  = acc_prob_vmc + temp_prob/N_walk 
      end do
      !$OMP END PARALLEL DO

      do i = N_walk + 1, N_max
                flag(i) = .False.    !Walker is set to dead
      end do
      call cpu_time(toc)
      deallocate(temp_coords)

      !DIAGNOSTICS ABOUT INITIAL CONFIGURATION
      write(*,*)
      write(*,*) '-----------------------------------------------'
      write(*,*) 'Initial Distribution setted '
      write(*,*) 'CPU time : ', toc - tic
      write(*,*) 'Typical acceptance ratio : ', acc_prob_vmc
      write(*,*)
      write(*,*)



      deallocate(configurations)
      deallocate(flag)
end program DMC_omp
