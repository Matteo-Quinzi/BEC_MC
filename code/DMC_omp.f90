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

     
      !omp variables
      !It is useful to get the thread id in order to seed the rng
      integer :: thread_id


      !Date and Time
      character(len=8) :: date
      character(len=10) :: time

      !DMC variables 
      integer(kind=8) :: N_at, N_walk, N_max      !N_atoms, walkers and maximum number of walkers
      integer(kind=8) :: eq_it, samples, dt_sam   !N of iterations fir equilibration, samples and delta_time between samples
      real(kind=8)    :: dt                       !timestep width
      real(kind=8)    :: Er                       !Trial energy
      real(kind=8)    :: a, b0, b1                !Scattering length and params of the guiding function
      character(len=50) :: coords_input_file, eq_out_file
      integer(kind=8) :: Nl ! # of points where to evaluate the radial density and density matrix
      real(kind=8) :: r_max !Maximum radius where to look for the OBDM
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
      real(kind=8),allocatable :: acc_prob_eq(:)                 !Acceptance probability in equilibration procedure for each step
      real(kind=8), allocatable :: temp_prob_eq(:)
      logical :: eq_accepted
      integer(kind=8) :: N_walk_start
      real(kind=8) :: N_ratio, alpha

      !Sampling variables
      integer(kind=8) :: sample_it
      real(kind=8), allocatable :: rho_rad(:), r_at(:), r_mesh(:)
      integer(kind=8),allocatable :: rho_rad_walk(:,:)
      real(kind=8) :: mean_Et
      integer(kind=8), allocatable :: obdm_walk(:,:,:)    !OBDM  for all walkers (N_max, Nl, Nl)
      real(kind=8), allocatable :: obdm_lzero(:,:)

      !Variables in OBDM diagonalization
      integer :: M      !number of evaluated eigenvalues
      real(kind=8), allocatable :: occupation_numbers(:)
      integer :: ldz 
      real(kind=8), allocatable :: z(:,:)
      real(kind=8), allocatable :: work(:)
      integer :: lwork
      integer, allocatable :: iwork(:)
      integer, allocatable :: ifail(:)
      integer :: info
      character(len=50) :: eigen_out_file
      real(kind=8) :: check_sum
      integer(kind=8) :: eigenvec_idx
      real(kind=8), allocatable :: identity_mat(:,:)
      real(kind=8) :: r_ghost

      !Pi 
      real(kind=8) :: pi_greek 



      !CPU time variables
      real(kind=8) :: tic, toc

      !Loops indexes
      integer(kind=8) :: i, j, it

      pi_greek = acos(-1.d0)
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
                     dt_sam, dt, Er, a, b0, b1, coords_input_file, eq_out_file,&
                     Nl, r_max)
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
      write(*,'(a27,a50)')       'Equilibration output file: ', eq_out_file
      write(*,'(a27,i10)')       'Meshgrid-side points     : ', Nl
      write(*,'(a27,f15.10)')    'Maximum radius for OBDM  : ', r_max

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
      allocate(Nt(eq_it), Et(eq_it), acc_prob_eq(eq_it))
      allocate(walk_en(N_max), walk_en_old(N_max))
      allocate(F_driv(N_max,N_at,3))

      !Evaluating initial driving forces and initial local energies
      !$OMP PARALLEL DO 
      do i = 1, N_walk
          F_driv(i,:,:) = driving_force(a, b0, b1, N_at, configurations(i, :, :))
          walk_en(i)    = local_energy(a, b0, b1, N_at, configurations(i,:,:)) 
      end do
      !$OMP END PARALLEL DO

      allocate(temp_prob_eq(N_max))

      !Saving initial number of walker
      N_walk_start = N_walk

      !Fixing population control parameter
      !alpha = N_at * N_walk_start / dt
      alpha = 1.0d0 / (N_walk_start * dt)

      tic = omp_get_wtime()
      !Equilibration loop
      do it = 1, eq_it
          acc_prob_eq(it) = 0.d0

          !Diffusing alive walkers
          !$OMP PARALLEL DO reduction(+:acc_prob_eq)
          do i = 1, N_walk
              !Saving starting walker energy
              walk_en_old(i) = walk_en(i)
              !Each walker diffuse with accept/reject step
              call one_walker_diffusion(a, b0, b1, N_at, dt, configurations(i, :, :), F_driv(i, :, :), walk_en(i), eq_accepted)
              if (eq_accepted .eq. .True.) acc_prob_eq(it) = acc_prob_eq(it) + 1.d0/N_walk
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

      end do
      toc = omp_get_wtime()

      !DIAGNOSTIC PRINTOUT AFTER EQUILIBRATION
      write(*,*) '--------------------------------------'
      write(*,*) 'Equilibration concluded'
      write(*,*) 'Equilibration time : ', toc - tic
      write(*,*) 'Average equilibration energy   : ', sum(Et) / eq_it
      write(*,*) 'Average acceptance probability : ', sum(acc_prob_eq) / eq_it
      write(*,*) 'Average number of walkers      : ', 1.d0*(sum(Nt)) / eq_it
      write(*,*)
      write(*,*) 
      
      call save_rank_dmc(0, eq_it, Nt, Et, acc_prob_eq, eq_out_file) 

      deallocate(temp_prob_eq)
      deallocate(Nt, Et, acc_prob_eq)

      !SAMPLING PHASE
      ! After equilibration we continue with the diffusion algorithm
      ! to intermittently sample the configuration space. Each sample 
      ! is used to build up the radial density distribution function
      ! and the (radial part of the) One Body Density Matrix 

      !definig number of total sampling iterations 
      sample_it = samples * dt_sam

      !We evaluate energy only at intermediate steps now

      !Allocating radial quantities
      allocate(rho_rad(Nl), r_at(N_at), r_mesh(Nl))
      allocate(rho_rad_walk(N_max,Nl))
      allocate(obdm_walk(N_max,Nl,Nl))
      allocate(obdm_lzero(Nl, Nl))

      !Initializing mesh and rho
      !r_mesh = define_mesh(Nl, r_max)
      rho_rad = 0.d0  
      obdm_lzero = 0.d0

      tic = omp_get_wtime()
      !Diffusion cycle works as before
      do it = 1, sample_it

          !Diffusing alive walkers
          !$OMP PARALLEL DO 
          do i = 1, N_walk
              !Saving starting walker energy
              walk_en_old(i) = walk_en(i)
              !Each walker diffuse with accept/reject step
              call one_walker_diffusion(a, b0, b1, N_at, dt, configurations(i, :, :), F_driv(i, :, :), walk_en(i), eq_accepted)
          end do
          !$OMP END PARALLEL DO

          !Branching
          !N_walk, walk_en, F_driv are going to be updated with the number of replicas
          call branching(N_walk, N_max, N_at, configurations, walk_en, walk_en_old, F_driv, flag, dt, Er)

          !Adjusting the energy scale
          if ( mod(it,100) .eq. 0) then
              mean_Et = sum(walk_en(1:N_walk)) / N_walk
              N_ratio = (N_walk_start*1.d0) / (N_walk*1.d0)
              Er = mean_Et + alpha * log(N_ratio)
          end if

          !HERE WE ADD THE SAMPLING STEP
          if ( mod(it,dt_sam) .eq. 0) then
                  !$OMP PARALLEL DO private(r_mesh, r_at)
                  do i = 1,N_walk
                      r_mesh = define_mesh(Nl, r_max)
                      r_at = evaluate_atomic_distances(N_at, configurations(i,:,:))
                      rho_rad_walk(i,:) =  one_walk_radial_distribution(N_at, r_at, Nl, r_mesh)
                      obdm_walk(i,:,:) = one_walk_radial_obdm_zero(N_at, r_at, Nl, r_mesh)
                  end do
                  !$OMP END PARALLEL DO

                  !Putting it all together
                  rho_rad(:) = (/( rho_rad(i) + (1.d0*sum(rho_rad_walk(1:N_walk,i)))/(1.d0*N_walk) , i=1,Nl )/)
                  
                  do i = 1,Nl
                      do j = i,Nl
                          obdm_lzero(i,j) = obdm_lzero(i,j) + ( sum(obdm_walk(1:N_walk,i,j))*1.d0 ) / (1.d0*N_walk)
                      end do
                  end do
          end if
      end do

      !Evaluate public mesh_grid
      r_mesh = define_mesh(Nl, r_max)

      toc = omp_get_wtime()

      rho_rad = rho_rad/(N_at*samples*1.d0)
      obdm_lzero = obdm_lzero/(N_at*samples*1.d0) 

      !Regularize the OBDM 
      do i = 1, Nl-1
          do j = i+1, Nl
              obdm_lzero(i,j) =  obdm_lzero(i,j) / (4.d0*pi_greek*r_mesh(i)*r_mesh(j))
          end do
      end do

      check_sum = 0.d0
      do i=1,Nl
       check_sum = check_sum + obdm_lzero(i,i)
      end do


      !DIAGNOSTIC PRINTOUT
      write(*,*) '--------------------------------------'
      write(*,*) 'Sampling procedure completed'
      write(*,*) 'Time taken : ', toc-tic
      write(*,*) 'Saving radial distribution function in the out file'
      write(*,*) 'Radial distribution sums to : ',  sum(rho_rad)
      write(*,*) 'Check sum                   : ', check_sum
      !call save_rho_rad(Nl, r_mesh, rho_rad)
      write(*,*)
      write(*,*)

      !DIAGONALIZING THE OBDM 
      ! The obdm (l=0) is diagonalized to obtain the natural orbitals
      ! and their occupation values. The macroscopically occupated orbital will
      ! give the fraction of atoms within the condensate

      write(*,*) '-----------------------------------------'
      write(*,*) 'Diagonalizing the OBDM (l=0)'


      !I need LDZ=Nl, Z(Nl,Nl)
      !I need Lwork=8*Nl, work(lwork)
      !iwork(5*Nl)
      !ifail(N)
      !info
      allocate(identity_mat(Nl, Nl))
      allocate(occupation_numbers(Nl))
      ldz = Nl
      allocate(z(ldz,Nl))
      lwork = 8*Nl
      allocate(work(lwork))
      allocate(iwork(5*Nl))
      allocate(ifail(Nl))

      identity_mat = 0.d0
      do i = 1,Nl
          identity_mat(i,i) = 1.d0
      end do

      !Redefining obdm to have larger eigenvalues becoming smaller eigenvalues
      !obdm_lzero = identity_mat - obdm_lzero
      write(*,*) (obdm_lzero(i,i), i=1,Nl)

      tic = omp_get_wtime()
      call dsyevx('V', 'A', 'U', Nl, obdm_lzero, Nl, &
                  0.d0, 1.d0, 1, 10, 1.d-10, M, &
                  occupation_numbers, z, ldz, work, lwork, iwork, ifail, info)
      toc = omp_get_wtime()

      !occupation_numbers = 1.d0 - occupation_numbers
      occupation_numbers(M:Nl) = 0.d0

      !DIAGNOSTIC PRINTOUT
      write(*,*)
      write(*,*) 'Diagonalization completed'
      write(*,*) 'Time in dsyevx : ', toc-tic
      write(*,*) 'Number of eigenvalues found : ', M
      write(*,*) 'Largest eigenvalue found : ', maxval(occupation_numbers(:))
      write(*,*) 'Saving eigenvalues and eigenvectors of the OBDM in eigenvec.out'
      write(*,*) 'Summing eigenvalues :', sum(occupation_numbers(1:M))

      write(*,*) 'Index of maximum eigenvalue : ' , maxloc(occupation_numbers(:))
      eigen_out_file = 'eigenvec.out'
      !call save_eigenvectors(Nl, 1, r_mesh, obdm_lzero(:,maxloc(occupation_numbers(:))), eigen_out_file)
      call save_rho_rad(Nl, r_mesh, rho_rad, obdm_lzero(:,maxloc(occupation_numbers))) 

      deallocate(obdm_walk, obdm_lzero)
      deallocate(walk_en, walk_en_old)
      deallocate(F_driv)
      deallocate(configurations)
      deallocate(flag)

      write(*,*) 'PROGRAM COMPLETED CORRECTLY'
end program DMC_omp
