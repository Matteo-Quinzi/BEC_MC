module bec_dmc
        !
        use bec_vmc
        implicit none 
        real(kind=8) :: pi=acos(-1.d0)
        !
        contains 
!---------------------------------------------------------------------------------------------------------------------------------
        function gaussian_rng(N_walk, N_at, N_dim) result (grn_matrix)
                ! Returns an array of gaussian random numbers obtained through 
                ! Box-Muller transform
                integer(kind=8), intent(in) :: N_walk, N_at, N_dim
                real(kind=8) :: grn_matrix(N_walk, N_at, N_dim)
                real(kind=8),allocatable :: temp_urn(:), temp_grn(:)
                integer(kind=8) :: L,i
                real(kind=8) :: prefactor 

                !check that you have an even number of numbers
                if (mod(N_walk*N_at*N_dim,2) .eq. 0) then 
                    allocate(temp_grn(N_walk*N_at*N_dim))
                else         
                    allocate(temp_grn(N_walk*N_at*N_dim + 1))
                end if
                L = size(temp_grn)/2

                !extracts N_at*N_walk*N_dim(+1) random numbers
                !sampling the uniform distribution
                call random_number(temp_grn)

                !transform them into gaussian distributed random numbers through 
                !Box-Muller transform
                do i = 0,L-1
                    prefactor = sqrt(-2.d0*log(temp_grn(2*i+1)))
                    temp_grn(2*i+1:2*i+2) = (/ prefactor*cos(2.d0*pi*temp_grn(2*i+2)), prefactor*sin(2.d0*pi*temp_grn(2*i+2))  /) 
                end do

                !Reshaping the array in the correct form
                if (mod(N_walk*N_at*N_dim,2) .eq. 0) then 
                    grn_matrix = reshape(temp_grn, (/N_walk,N_at,N_dim/))
                else
                    grn_matrix = reshape(temp_grn(1:2*L-1), (/N_walk,N_at,N_dim/))
                end if 
                
        end function gaussian_rng
!-----------------------------------------------------------------------------------------------------------------------------------

        function driving_force(a, b0, b1, N_at, coords) result (F)
                ! This is actually half of the driving force
                real(kind=8), intent(in)  :: a, b0, b1
                integer(kind=8), intent(in) :: N_at
                real(kind=8), intent(in)  :: coords(N_at,3)
                real(kind=8) :: F(N_at,3)
                real(kind=8) :: der_f_over_f(N_at,N_at,3)
                real(kind=8) :: r(3), rij, xij(3), temp_f
                real(kind=8) :: r2
                integer :: i,j,k

                !Building up 
                do i = 1, N_at
                    do j = i, N_at
                        ! Evaluate relative coordinates between atom i and atom j
                        xij(:) = (/ (coords(i,k) - coords(j,k) , k=1,3) /)
                        rij = sqrt( xij(1)**2.d0 + xij(2)**2.d0 + xij(3)**2.d0 )

                        !Evaluate f'/f for each spatial dimension
                        if (rij < a ) then 
                                der_f_over_f(i,j,:) = (/ (0.d0 , k=1,3) /)
                        else
                                temp_f = fast_f(a,rij)
                                der_f_over_f(i,j,:) = (/( fast_first_der_f(a,xij(k),rij)/temp_f, k=1,3 )/)
                        end if
                        !Impose symmetry condition
                        der_f_over_f(j,i,:) = der_f_over_f(i,j,:)
                    end do 
                end do

                do i = 1, N_at
                    r = coords(i,:)
                    r2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
                    do k = 1,3
                        F(i,k) = first_derg_over_g(b0, b1, r(k), r2)
                        F(i,k) = F(i,k) + sum(der_f_over_f(i,:,k))
                    end do
                end do

                F = 2.d0 * F

        end function driving_force

!----------------------------------------------------------------------------------------------------------------------------------

        function dmc_acc_prob(a, b0, b1, N_at, coords, new_coords, F_driv, F_driv_new, &
                        E_loc, E_loc_new, dt) result (acc_prob)
                !Evaluate the acceptance probability in the dmc accept/reject step
                !Following the derivation of Ceperley we need to build the ratio as
                !        G(R' -> R) psi_T(R')^2 / G(R -> R') psi_T(R)^2
                ! this way if psi_T(R') = 0 (forbidden configuration) the probability is 
                ! automatically null and the move is surely rejected. This is a nice way
                ! to control the correctness of the hard sphere prescription
                ! after all atoms have been displaced 
                real(kind=8), intent(in) :: a, b0, b1
                integer(kind=8), intent(in) :: N_at
                real(kind=8), intent(in) :: coords(N_at, 3)
                real(kind=8), intent(in) :: new_coords(N_at, 3)
                real(kind=8), intent(in) :: F_driv(N_at,3)
                real(kind=8), intent(in) :: F_driv_new(N_at,3)
                real(kind=8), intent(in) :: E_loc, E_loc_new
                real(kind=8), intent(in) :: dt
                real(kind=8) :: acc_prob
                real(kind=8) :: psi_ratio
                real(kind=8) :: fnew_over_fold(N_at)
                real(kind=8) :: G_old_new_mat(N_at,3), G_new_old_mat(N_at,3)
                real(kind=8) :: G_old_new, G_new_old
                integer(kind=8) :: i, j

                psi_ratio = 1.d0
                do i = 1, N_at
                    psi_ratio = psi_ratio * g(b0, b1, new_coords(i,:)) / g(b0, b1, coords(i,:))
                    do j = 1,N_at
                        if (j .eq. i) then
                            fnew_over_fold(j) = 1.d0
                        else
                            fnew_over_fold(j) = f(a, new_coords(i,:), new_coords(j,:)) / f(a, coords(i,:), coords(j,:))
                        end if
                    psi_ratio = psi_ratio * fnew_over_fold(j)
                    end do
                end do 

                ! This happens when two atoms compenetrate in the new position
                if (psi_ratio .eq. 0.d0) then
                        acc_prob = 0.d0
                else
                    ! If hard sphere condition is respected then go on with the other evaluation
                    G_old_new_mat = new_coords - coords + 0.5d0 * F_driv * dt
                    G_new_old_mat = coords - new_coords + 0.5d0 * F_driv_new * dt
                    G_old_new = 0.d0
                    G_new_old = 0.d0
                    do i = 1,3
                        G_old_new = G_old_new + dot_product(G_old_new_mat(:,i), G_old_new_mat(:,i))
                        G_new_old = G_new_old + dot_product(G_new_old_mat(:,i), G_new_old_mat(:,i))
                    end do
                    acc_prob = - (G_new_old - G_old_new) / (2.d0*dt) - dt * (E_loc_new - E_loc)
                    acc_prob =  exp(acc_prob)
                    acc_prob = acc_prob * psi_ratio
                end if


        end function dmc_acc_prob

!-----------------------------------------------------------------------------------------------------------------------------------

        subroutine one_walker_diffusion(a, b0, b1, N_at, dt, coords, F, E_loc, accepted)
              !Subroutine that changes the coordinates according to the diffusion step
              ! in the DMC algorithm. 
              !The diffusion step takes place as 
              !       R_new = R_old + 1/2 * F(R_old) * dt + e * (dt)**1/2
              ! with F being the guiding force and e a gaussian random number with mean value 0
              ! and variance 1.
              !An accept/reject step is introduced (Ceperley-Alder implementation) to avoid drifts
              ! due to the integration error O(3/2) in dt
              !To optimize the calculation the driving force is evaluated outside and can be changed
              real(kind=8), intent(in) :: a, b0, b1
              integer(kind=8), intent(in) :: N_at
              real(kind=8), intent(in) :: dt
              real(kind=8), intent(inout) :: coords(N_at,3)
              real(kind=8), intent(inout) :: F(N_at,3)
              real(kind=8), intent(inout) :: E_loc
              logical, intent(inout) :: accepted
              real(kind=8) :: e(1,N_at,3)
              real(kind=8) :: rand_acc
              real(kind=8) :: F_new(N_at,3)
              real(kind=8) :: E_loc_new
              real(kind=8) :: new_coords(N_at,3)
              real(kind=8) :: acc_prob
              integer(kind=8) :: i,j

              
              !Extracting a sample of gaussian random numbers
              e = gaussian_rng(1,N_at,3)

              !Extracting a sample of N_at uniform random numbers
              call random_number(rand_acc) ! these will be used in the accept/ reject step

              !Diffusing all atoms
              do i = 1, N_at
                  new_coords(i,:) = (/ (coords(i,j) + 0.5d0 * F(i,j) * dt + e(1,i,j)*sqrt(dt), j=1,3) /)
              end do

              !I need to re-evaluate all driving forces 
              ! and the local energy to evaluate the acceptance prob
              F_new = driving_force(a, b0, b1, N_at, new_coords)
              E_loc_new = local_energy(a, b0, b1, N_at, new_coords) 

              !Evaluating acceptance prob
              acc_prob = dmc_acc_prob(a, b0, b1, N_at, coords, new_coords, &
                                      F, F_new, E_loc, E_loc_new, dt)
              if (acc_prob .gt. rand_acc) then
                      coords = new_coords
                      F = F_new
                      E_loc = E_loc_new
                      accepted = .True.
              else 
                      accepted = .False.
              end if

        end subroutine one_walker_diffusion
!----------------------------------------------------------------------------------------------------------------------------------

        subroutine branching(N_walk, N_max, N_at, configurations, walk_en, F_driv, flag, dt, Er)
              !This subroutine evaluate the branching factor for each alive walker
              ! then add new replicas to the configurations array
              integer(kind=8), intent(inout) :: N_walk
              integer(kind=8), intent(in) :: N_max, N_at
              real(kind=8), intent(inout) :: configurations(N_max, N_at, 3)
              real(kind=8), intent(inout) :: walk_en(N_max)
              real(kind=8), intent(inout) :: F_driv(N_max, N_at, 3)
              logical, intent(inout) :: flag(N_max)
              real(kind=8), intent(in) :: dt, Er
              integer(kind=8) :: W(N_walk)
              real(kind=8) :: rd_shift(N_walk)
              integer(kind=8) :: i, N_walk_new
              logical, dimension(N_walk) :: j0


              !At the beginning N_walk is the number of alive walkers
              N_walk_new = N_walk

              ! Extract N_walk randomd shifts
              call random_number(rd_shift)

              !Evaluate the branching value taking the nearest integer
              W(:) = (/ ( nint(exp( - dt * (walk_en(i) - Er)) + rd_shift) , i=1,N_walk ) /)

              !I need to sign which walkers are dead and which ones are duplicating
              !I am allowing a single walker to generate up to 4 replicas
              do i = 1, N_walk
                  if ( W(i) .eq. 0 ) then
                          !Kills walker at position i
                          flag(i) = .False.
  
                  else if ( W(i) .eq. 2) then
                          ! Adding one walker to the tail
                          N_walk_new = N_walk_new + 1
                          flag(N_walk_new) = .True.             
                          configurations(N_walk_new,:,:) = configurations(i,:,:)
                          walk_en(N_walk_new) = walk_en(i)
                          F_driv(N_walk_new,:,:) = F_driv(i,:,:)

                  else if ( W(i) .eq. 3) then
                          !Adding two walkers to the tail
                          N_walk_new = N_walk_new + 1
                          flag(N_walk_new) = .True.
                          configurations(N_walk_new,:,:) = configurations(i,:,:)
                          walk_en(N_walk_new) = walk_en(i)
                          F_driv(N_walk_new,:,:) = F_driv(i,:,:)

                          N_walk_new = N_walk_new + 1
                          flag(N_walk_new) = .True.
                          configurations(N_walk_new,:,:) = configurations(i,:,:)
                          walk_en(N_walk_new) = walk_en(i)
                          F_driv(N_walk_new,:,:) = F_driv(i,:,:)

                   else if ( W(i) .ge. 4) then
                          !Adding three walkers to the tail
                          N_walk_new = N_walk_new + 1
                          flag(N_walk_new) = .True.
                          configurations(N_walk_new,:,:) = configurations(i,:,:)
                          walk_en(N_walk_new) = walk_en(i)
                          F_driv(N_walk_new,:,:) = F_driv(i,:,:)

                          N_walk_new = N_walk_new + 1
                          flag(N_walk_new) = .True.
                          configurations(N_walk_new,:,:) = configurations(i,:,:)
                          walk_en(N_walk_new) = walk_en(i)
                          F_driv(N_walk_new,:,:) = F_driv(i,:,:)

                          N_walk_new = N_walk_new + 1
                          flag(N_walk_new) = .True.
                          configurations(N_walk_new,:,:) = configurations(i,:,:)
                          walk_en(N_walk_new) = walk_en(i)
                          F_driv(N_walk_new,:,:) = F_driv(i,:,:)

                   end if

              end do

              !Adjust walkers positions
              do i = 1,N_walk
                  if (flag(i) .eq. .False.) then
                      !Takes one walker from the tail and put it here
                      configurations(i,:,:) = configurations(N_walk_new,:,:)
                      walk_en(i) = walk_en(N_walk_new)
                      F_driv(i,:,:) = F_driv(N_walk_new,:,:)
                      flag(N_walk_new) = .False.     !Last walker is set to dead
                      flag(i) = .True. !New walker is set to alive
                      N_walk_new = N_walk_new - 1
                  end if
              end do

              if (N_walk_new .eq. N_max) then
                           write(*,*) '-----------------------'
                           write(*,*) 'ERROR : Reached maximum number of walkers !'
                           write(*,*)
                           stop
              else if (N_walk_new .eq. 0) then
                           write(*,*) '-----------------------'
                           write(*,*) 'All walkers are dead :('
                           write(*,*)
                           stop
              end if

              N_walk = N_walk_new

        end subroutine branching

!-----------------------------------------------------------------------------------------------------------------------------------
end module bec_dmc
