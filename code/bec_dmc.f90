module bec_dmc
        !
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
        
!-----------------------------------------------------------------------------------------------------------------------------------
end module bec_dmc
