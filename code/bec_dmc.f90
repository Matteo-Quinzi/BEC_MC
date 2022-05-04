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
                integer, intent(in) :: N_walk, N_at, N_dim
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


        end function driving_force


!-----------------------------------------------------------------------------------------------------------------------------------
end module bec_dmc
