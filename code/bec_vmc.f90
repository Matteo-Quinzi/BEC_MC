module bec_vmc
       use mpi
       ! This module contains the functions/ routines used in 
       ! the VMC study of a boson gas.

       implicit none
       
       real(kind=8), private, parameter :: h=1.d-4 !Finite differences step 

       contains

!--------------------------------------------------------------------------------------------------

       subroutine init_rng(my_rank)
               !Initialize the rng to avoid having same random numbers
               !for all processes and for all runs
               integer, intent(in) :: my_rank
               character(len=8) :: date
               character(len=10) :: time
               integer :: seed_size, time_int, i
               integer, allocatable :: seed_array(:)

               !Convert time in integer
               call date_and_time(date,time)
               read(time(:6),*) time_int

               !Get the seed of the rng
               call random_seed(size=seed_size)
               allocate(seed_array(seed_size))
               call random_seed(get=seed_array)

               !Modify the seed
               seed_array(:) = (/ (seed_array(i) + 100*my_rank - time_int, i=1,seed_size) /)

               !Insert back the new seed
               call random_seed(put=seed_array) 

       end subroutine

! -------------------------------------------------------------------------------------------------

       subroutine fix_rng(my_rank)
               !Initialize the rng to avoid having same random numbers
               !for all procs and for all runs
               ! THIS FIXES THE RNG FOR THE IFORT COMPILER ONLY
               ! NOT PORTABLE
               integer, intent(in) :: my_rank
               integer :: seed_size, time_int, i
               integer :: seed_array(2)
               
               seed_array(1) = 2147483562
               seed_array(2) = 2147483398
               !Put in the fixed seed
               call random_seed(put=seed_array)

       end subroutine fix_rng

! --------------------------------------------------------------------------------------------------
       function set_up_boson_gas(N_at, a, l) result(coords)
               ! This create a configuration of N_at bosons
               ! as the trap we are considering will be spherical 
               ! the extraction of coords is made in such a way to enforce
               ! the spherical symmetry
               integer(kind=8), intent(in) :: N_at    ! # atoms in the gas
               real(kind=8), intent(in) :: a          ! radius of the hard-sphere potential
               real(kind=8), intent(in) :: l          ! maximum radius for coords extraction
               real(kind=8) :: coords(N_at,3)         ! coords of the N_at bosons
               integer(kind=8) :: i,j
               real(kind=8) :: rho, phi, theta        !spherical coords
               real(kind=8) :: x, y, z                !cartesian coords
               real(kind=8) :: rij
               real(kind=8) :: pi = acos(-1.d0)
               integer :: check
               real(kind=8) :: tic, toc
               ! A fixed amount of iterations is fixed to avoid infinite looping
               ! due to bad initializang conditions
               ! (e.g. small radius compared to number and size of atoms)
               integer(kind=8) :: maxit=1000 , counter 
               
               maxit = 1000 * N_at

               call cpu_time(tic)

               ! first move is always good
               call random_number(rho)
               call random_number(phi)
               call random_number(theta)

               rho = rho * l
               phi = phi * 2.d0 * pi
               theta = theta * pi

               x = rho * sin(theta) * cos(phi)
               y = rho * sin(theta) * sin(phi)
               z = rho * cos(theta)

               coords(1,:) = (/x, y, z/)

               do i = 2, N_at
                    
                   counter = 0           
                   check = 1

                   do while (check == 1)
                   counter = counter + 1
                   if (counter >= maxit ) then 
                      write(*,*)
                      write(*,*)
                      write(*,*) 'Could not find a starting configuration'
                      write(*,*) 'Please check if values of initialization radius, number of atoms'
                      write(*,*) ' or scattering lengths are compatible.'
                      write(*,*) 
                      write(*,*)
                      write(*,*) '-------------------------------------------------------------------'
                      stop 
                   end if     
                        !extract three random numbers
                        call random_number(rho)
                        call random_number(phi)
                        call random_number(theta)

                        !transform in polar coords
                        rho = rho * l
                        phi = phi * 2.d0 * pi
                        theta = theta * pi

                        !transform in cartesian coords
                        x = rho * sin(theta) * cos(phi)
                        y = rho * sin(theta) * sin(phi)
                        z = rho * cos(theta)

                        do j = 1, (i-1)
                            rij = sqrt( (x - coords(j,1))**2.d0 + &
                                        (y - coords(j,2))**2.d0 + &
                                        (z - coords(j,3)) **2.d0 )
                            if (rij <= a) check = 0
                        end do

                        if (check == 0) then
                                check = 1 
                                continue
                        else
                                exit   ! exit from while loop if all rij > a
                        end if

                    end do     ! while loop

                    coords(i,:) = (/x, y, z/) ! assign to i-th boson the last value of x,y,z
                end do     !end i loop

                call cpu_time(toc)
               write(*,*)
               write(*,*)
               write(*,*) 'Boson gas has been successfully setted up'
               write(*,*) 'CPU_time in gas creation : ', toc - tic
               write(*,*)
               write(*,*)
               write(*,*) '-------------------------------------------'
       end function set_up_boson_gas

!-----------------------------------------------------------------------------------------------------------------

       function set_up_from_file(N_at, coords_input_file) result (coords)
               ! Initial coords are taken from an input text file
               ! can be useful to avoid re-equilibration
               ! (it seems that equilibration can take a long time with N_at > 10)
               integer(kind=8), intent(in) :: N_at
               character(len=*), intent(in) :: coords_input_file
               real(kind=8) :: coords(N_at,3)
               integer(kind=8) :: i, dummy

               open(unit=10, file=coords_input_file, action='read')
               read(10,*) !Reading the header
               do i =1, N_at
                   read(10,'(i6,4x,f16.13,4x,f16.13,4x,f16.13)') dummy, coords(i,1), coords(i,2), coords(i,3) 
               end do
               close(10)

       end function set_up_from_file

!-----------------------------------------------------------------------------------------------------------------

       function g(b0, b1, r)
              real(kind=8), intent(in) :: b0, b1
              real(kind=8), intent(in) :: r(3)
              real(kind=8) :: g
              real(kind=8) :: r2

              r2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
              g = exp(- ( b0*r2 + b1*r2*r2 ) )    
 
       end function g

!-----------------------------------------------------------------------------------------------------------------

       function first_derg_over_g(b0, b1, r_idx, r2)
               real(kind=8), intent(in) :: b0, b1
               real(kind=8), intent(in) :: r_idx
               real(kind=8), intent(in) :: r2
               real(kind=8) :: first_derg_over_g
               first_derg_over_g = -2.d0*r_idx * (b0 - 2.d0*b1*r2)                  
       end function first_derg_over_g

!------------------------------------------------------------------------------------------------------------------
      
       function f(a, ri, rj)
               real(kind=8), intent(in) :: a
               real(kind=8), intent(in) :: ri(3), rj(3)
               real(kind=8) :: f, rij

               rij = sqrt( (ri(1)-rj(1))**2.d0 + &
                           (ri(2)-rj(2))**2.d0 + &
                           (ri(3)-rj(3))**2.d0 )
               if (rij <= a) then 
                       f = 0.d0
               else 
                       f = 1.d0 - a/rij
               end if
        end function f

!------------------------------------------------------------------------------------------------------------------

        function fast_f(a, rij)
                !You need to check outside if the hard sphere condition is reached
                real(kind=8), intent(in) :: a
                real(kind=8), intent(in) :: rij
                real(kind=8) :: fast_f

                fast_f = 1.d0 - a/rij

        end function fast_f

!-------------------------------------------------------------------------------------------------------------------

        function first_der_f(a, ri, rj, idx)
                !derivative of f wrt ri
                real(kind=8), intent(in) :: a
                real(kind=8), intent(in) :: ri(3), rj(3)
                integer(kind=8), intent(in) :: idx !derivative is made wrt ri-x, j or z
                real(kind=8) :: rij
                real(kind=8) :: first_der_f

                rij = sqrt( (ri(1) - rj(1))**2.d0 + &
                            (ri(2) - rj(2))**2.d0 + &
                            (ri(3) - rj(3))**2.d0 )
                if (rij <= a) then 
                    first_der_f = 0.d0
                else
                    first_der_f =  a * ( ri(idx) - rj(idx) ) / (rij**3.d0 )
                end if
        
        end function first_der_f

!-------------------------------------------------------------------------------------------------------------------

        function fast_first_der_f(a,xij,rij)
                ! You need to check outside if the hard sphere condition is reached
                real(kind=8), intent(in) :: a
                real(kind=8), intent(in) :: xij, rij
                real(kind=8) :: fast_first_der_f

                fast_first_der_f = a * xij / (rij**3.d0)

        end function fast_first_der_f

!--------------------------------------------------------------------------------------------------------------------

        function V_ext(r)
               real(kind=8), intent(in) :: r(3)
               real(kind=8) :: V_ext 
               
               V_ext = 0.5d0 * (r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
        end function

!------------------------------------------------------------------------------------------------------------------
        function lap_single_atom(a, b0, b1, N_at, coords, which) result (nabla)
                ! This evaluates the laplacian associated to the atom 'which'
                ! The laplacian is then normalized with the wavefunction, otherwise
                ! we can get values below the machine precision treshold
                real(kind=8), intent(in) :: a, b0, b1
                integer(kind=8), intent(in) :: N_at 
                real(kind=8), intent(in) :: coords(N_at, 3)
                integer(kind=8), intent(in) :: which
                real(kind=8) :: nabla
                integer(kind=8) :: i,j
                real(kind=8) :: disp(3,3)
                real(kind=8) :: g_still, g_low, g_up, up_contribution, low_contribution
                real(kind=8), dimension(N_at) :: f_still_vec, f_up_vec, f_low_vec 

                disp = reshape((/h, 0.d0, 0.d0, 0.d0, h, 0.d0, 0.d0, 0.d0, h/), shape(disp))
                g_still = g(b0, b1, coords(which,:))

                f_still_vec(:) = (/ ( f(a, coords(which,:), coords(j,:)), j=1,N_at)/)
                f_still_vec(which) = 1.d0 ! so that f_up/f_still = 1 if f_up is also 1

                nabla = 0.d0
                do i =1,3
                    g_up = g(b0, b1, coords(which,:) + disp(i,:))
                    f_up_vec = (/ (f(a, coords(which,:) + disp(i,:), coords(j,:)), j =1,N_at) /)
                    f_up_vec(which) = 1.d0

                    g_low = g(b0, b1, coords(which,:) - disp(i,:))
                    f_low_vec = (/ (f(a, coords(which,:) - disp(i,:), coords(j,:)), j =1,N_at) /)
                    f_low_vec(which) = 1.d0
                    
                    up_contribution = g_up/g_still
                    low_contribution = g_low/g_still

                    do j = 1, N_at
                         up_contribution = up_contribution * (f_up_vec(j)/f_still_vec(j))
                         low_contribution = low_contribution * (f_low_vec(j)/f_still_vec(j))
                    end do

                    nabla = nabla + (up_contribution/N_at -2.d0/N_at + low_contribution/N_at)/(h**2.d0)

                end do

        end function lap_single_atom

!--------------------------------------------------------------------------------------------------------------------

        function local_energy(a, b0, b1, N_at, coords) result (E_loc)
                ! This evaluates the local energy of the boson gas using 
                ! the easier wavefunction
                ! E_loc is given by the sum of N_at kinetic terms + N_at
                ! harmonic terms
                real(kind=8), intent(in) :: a, b0, b1
                integer(kind=8), intent(in) :: N_at
                real(kind=8), intent(in) :: coords(N_at,3)
                real(kind=8) :: E_loc
                integer(kind=8) :: i

                E_loc = 0.d0 
                do i = 1, N_at
                   !THIS IS WITH NUMERIC LAPLACIAN
                   E_loc = E_loc - 0.5d0*lap_single_atom(a, b0, b1, N_at, coords, i) +&
                           V_ext(coords(i,:))/N_at

                   !THIS IS WITH ANALITYC LAPLACIAN
                   !E_loc = E_loc - 0.5d0 * lap_psi_over_psi(a, b0, b1, N_at, coords, i) +&
                   !        V_ext(coords(i,:))/N_at
                end do

        end function local_energy

!--------------------------------------------------------------------------------------------------------------------

        function acceptance_prob(N_at, N_mark, a, b0, b1, coords, trial_r) result (acc_prob)
               integer(kind=8), intent(in) :: N_at, N_mark 
               real(kind=8), intent(in) :: a, b0, b1 
               real(kind=8), intent(in) :: coords(N_at,3), trial_r(3)
               real(kind=8) :: acc_prob
               real(kind=8) :: g_old, g_trial, f_vec(N_at), f_vec_trial(N_at)
               integer :: i,j

               g_old = g(b0, b1, coords(N_mark,:))
               g_trial = g(b0, b1, trial_r)

               f_vec(:) = (/ (f(a, coords(N_mark,:), coords(j,:)), j=1,N_at) /)
               f_vec(N_mark) = 1.d0 ! I want f_vec/f_vec_trial = 1 for N_mark
               f_vec_trial(:) = (/ (f(a, trial_r, coords(j,:)), j=1,N_at )/)
               f_vec_trial(N_mark) = 1.d0

               acc_prob = g_trial / g_old

               do i =1, N_at
                   acc_prob = acc_prob * ( f_vec_trial(i) / f_vec(i))
               end do

               acc_prob = acc_prob**2.d0

        end function acceptance_prob

!--------------------------------------------------------------------------------------------------------------------

        subroutine int_random_number(M, N_walk, N_at, array)
                integer(kind=8), intent(in) :: M, N_walk, N_at
                integer(kind=8), intent(inout) :: array(M,N_walk)
                real(kind=8) :: float_array(M,N_walk)
                integer :: i,j

                call random_number(float_array)

                do i = 1,N_walk
                    array(:,i) = (/  ( 1 + int(float(N_at)*(float_array(j,i))) , j = 1,M ) /) 
                end do  

        end subroutine int_random_number

!--------------------------------------------------------------------------------------------------------------------

        subroutine metropolis(a, b0, b1, N_at, N_walk, M, coords, mean_energies, acc_ratio)
                real(kind=8), intent(in) :: a, b0, b1
                integer(kind=8), intent(in) :: N_at, N_walk, M 
                real(kind=8), intent(inout) :: coords(N_walk,N_at,3)
                real(kind=8), intent(inout) :: mean_energies(M)
                real(kind=8) :: energies(N_walk, M)
                integer(kind=8) :: rd_which(M,N_walk)    !Tells at step M, for walker N which atom is trying to move
                real(kind=8) :: rd_where(M,N_walk,3)     !Gives random coordinates of the displacement
                real(kind=8) :: rd_acc(M,N_walk)
                real(kind=8) :: d                        !half displacement range
                real(kind=8) :: acc_prob(M,N_walk)
                integer(kind=8) :: acc_moves(N_walk)
                real(kind=8) :: acc_ratio
                integer(kind=8) :: i,j
                real(kind=8) :: avg_en(M)
                 

                acc_moves= 0
                d = 1.0d0 

                call int_random_number(M, N_walk, N_at, rd_which)
                call random_number(rd_where)
                rd_where = d * (2.d0*rd_where - 1.d0) 
                call random_number(rd_acc)

                ! Metropolis run of N_walk walkers

                energies(:,1) = (/ (local_energy(a, b0, b1, N_at, coords(j,:,:)), j=1,N_walk )/)
                mean_energies(1) = sum(energies(:,1)/N_walk)

                do i = 2, M
                ! At every step the trial move is evaluated for every walker
                    acc_prob(i,:) = (/ (acceptance_prob(N_at, rd_which(i,j), a, b0, b1, coords(j,:,:),&
                                     coords(j, rd_which(i,j),:) + rd_where(i,j,:)), j=1,N_walk) /)
                    do j = 1, N_walk
                        if (acc_prob(i,j) > rd_acc(i,j)) then
                                coords(j,rd_which(i,j),:) =  coords(j,rd_which(i,j),:) + rd_where(i,j,:)
                                acc_moves(j) = acc_moves(j) + 1

                                energies(j,i) = local_energy(a, b0, b1, N_at, coords(j,:,:))
                        else
                                energies(j,i) = energies(j,i-1)
                        end if
                    end do
                 
                    mean_energies(i) = sum(energies(:,i)/N_walk)

                    !DEBUG PRINTOUT         
                    !print *, 'Step : ',i, 'Average acc_prob', sum(acc_prob(i,:)/N_walk)
                    !print *, 'Step : ',i, 'Average Local Energy', sum(energies(:,i)/N_walk)
                    !print *, 'Step :',i, 'Average Local Energy', mean_energies(i)
                end do 

                acc_ratio = 0.
                do i = 1, N_walk
                    acc_ratio = acc_ratio + acc_moves(i)
                end do

                acc_ratio = acc_ratio/(M*N_walk)

        end subroutine metropolis 

!--------------------------------------------------------------------------------------------------------------------

        subroutine no_energy_metropolis(a, b0, b1, N_at, N_walk, M, coords, acc_ratio)
                !This is a Metropolis that doeasn't evaluate the local energy
                ! It is used to build up the initial configuration for the DMC run
                ! The initial wavefunction is sampled by evolving the coordinates
                real(kind=8), intent(in) :: a, b0, b1
                integer(kind=8), intent(in) :: N_at, N_walk, M 
                real(kind=8), intent(inout) :: coords(N_walk,N_at,3)
                integer(kind=8) :: rd_which(M,N_walk)    !Tells at step M, for walker N which atom is trying to move
                real(kind=8) :: rd_where(M,N_walk,3)     !Gives random coordinates of the displacement
                real(kind=8) :: rd_acc(M,N_walk)
                real(kind=8) :: d                        !half displacement range
                real(kind=8) :: acc_prob(M,N_walk)
                integer(kind=8) :: acc_moves(N_walk)
                real(kind=8) :: acc_ratio
                integer(kind=8) :: i,j

                acc_moves= 0
                d = 1.0d0 

                call int_random_number(M, N_walk, N_at, rd_which)
                call random_number(rd_where)
                rd_where = d * (2.d0*rd_where - 1.d0) 
                call random_number(rd_acc)

                ! Metropolis run of N_walk walkers

                do i = 1, M
                ! At every step the trial move is evaluated for every walker
                    acc_prob(i,:) = (/ (acceptance_prob(N_at, rd_which(i,j), a, b0, b1, coords(j,:,:),&
                                     coords(j, rd_which(i,j),:) + rd_where(i,j,:)), j=1,N_walk) /)
                    do j = 1, N_walk
                        if (acc_prob(i,j) > rd_acc(i,j)) then
                                coords(j,rd_which(i,j),:) =  coords(j,rd_which(i,j),:) + rd_where(i,j,:)
                                acc_moves(j) = acc_moves(j) + 1
                        end if
                    end do
                 
                end do 

                acc_ratio = 0.
                do i = 1, N_walk
                    acc_ratio = acc_ratio + acc_moves(i)
                end do

                acc_ratio = acc_ratio/(M*N_walk)

        end subroutine no_energy_metropolis 
!---------------------------------------------------------------------------------------------------------------------

        function err_energy(M, energies, mean_energy)
                ! Evaluate the error to the energy as the 
                ! std deviation on the number of steps
                integer(kind=8) :: M
                real(kind=8), intent(in) :: energies(M) !here energies should be already 
                                                        !averaged over walkers
                real(kind=8), intent(in) :: mean_energy  
                real(kind=8) :: err_energy
                integer(kind=8) :: i

                err_energy = 0.d0
                do i = 1, M 
                    err_energy = (energies(i) - mean_energy)**2.d0 / (M-1)
                end do
                err_energy = sqrt(err_energy)

        end function err_energy

end module bec_vmc
