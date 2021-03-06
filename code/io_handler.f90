      module io_handler 
      implicit none 
      character(len=80) :: arg, input_file, output_file
      integer :: io_unit=10

      contains 

!-------------------------------------------------------------------------------------------------------------

      subroutine read_input_file()
           character(len=128) :: command
           integer(kind=8):: i
    
          call get_command_argument(0,command)
          if (command_argument_count() .ne. 2) then
              write(*,*) 'Usage : ', trim(command), ' Input_file, &
              Output_file' 
              stop 
          else
              call get_command_argument(1,arg)
              read(arg,*) input_file
              call get_command_argument(2,arg)
              read(arg,*) output_file
          end if 
          
          write(*,*) 'Reading : ', input_file
          write(*,*) 'Output will be saved in : ', output_file
      end subroutine 

!--------------------------------------------------------------------------------------------------------------

      subroutine read_data(N_at, N_walk, M, N_cyc, &
                           a, b0, b1, l, f,coords_input_file)
          integer(kind=8) :: N_at, N_walk, M, N_cyc
          real(kind=8) :: a, b0, b1, l
          integer(kind=8) :: f
          character(len=50) :: coords_input_file

          open(unit=10, file=input_file)
          !reading N atoms
          read(10,*)
          read(10, '(4x, i10)') N_at
          read(10, '(4x, i10)') N_walk
          read(10,*)

          !reading step per loop and number of loops
          read(10,*)
          read(10, '(4x, i10)') M
          read(10, '(4x, i10)') N_cyc
          read(10,*)
          
          !reading scattering length
          read(10,*) 
          read(10,'(4x, f15.7)') a 
          read(10,*) 
          
          !reading variational params 
          read(10,*) 
          read(10,'(4x,f15.7)') b0
          read(10,'(4x,f15.7)') b1 
          read(10,*)

          !reading gas maximum radius
          read(10,*)
          read(10,'(4x,f15.7)') l
          read(10,*)

          !reading input from file
          read(10,*)
          read(10,'(14x,i1)') f
          read(10,'(14x,a50)') coords_input_file


          close(unit=10)
      end subroutine read_data

!-------------------------------------------------------------------------------------------------------------------

      subroutine read_data_dmc(N_at, N_walk, N_max, eq_it, samples, &
                               dt_sam, dt, Er, a, b0, b1, &
                               coords_input_file, eq_out_file, &
                               Nl, r_min, r_max)
          integer(kind=8) :: N_at, N_walk, N_max
          integer(kind=8) :: eq_it, samples, dt_sam
          real(kind=8) :: dt, Er
          real(kind=8) :: a, b0, b1
          character(len=50) :: coords_input_file, eq_out_file
          integer(kind=8) :: Nl !number of points on the meshgrid
          real(kind=8) :: r_min
          real(kind=8) :: r_max !Maximum radius where to look for OBDM 

          open(unit=io_unit, file=input_file, action='Read')
              !Reading Atoms and walkers
              read(io_unit,*)
              read(io_unit,'(10x,i10)') N_at
              read(io_unit,'(10x,i10)') N_walk
              read(io_unit,'(10x,i10)') N_max
              read(io_unit,*)

              !Reading iterations, samples, dt ...
              read(io_unit,*)
              read(io_unit,'(10x,i10)') eq_it
              read(io_unit,'(10x,i10)') samples
              read(io_unit,'(10x,i10)') dt_sam
              read(io_unit,'(10x,f15.10)') dt
              read(io_unit,'(10x,f15.10)') Er
              read(io_unit,*)

              !Reading info about the guiding function
              read(io_unit,*)
              read(io_unit,'(10x,f15.10)') a
              read(io_unit,'(10x,f15.10)') b0 
              read(io_unit,'(10x,f15.10)') b1
              read(io_unit,*)

              !Reading coords input file
              read(io_unit,*)
              read(io_unit,'(10x,a50)') coords_input_file
              read(io_unit,'(10x,a50)') eq_out_file
              read(io_unit,*)

              !Reading data for the OBDM and radial distribution
              !function
              read(io_unit,*)
              read(io_unit,'(10x,i10)') Nl
              read(io_unit,'(10x,f15.10)') r_min
              read(io_unit,'(10x,f15.10)') r_max
              read(io_unit,*)

          close(io_unit)
      end subroutine read_data_dmc

!-------------------------------------------------------------------------------------------------------------------

      subroutine print_coords(N,coords)
              integer(kind=8) :: N
              real(kind=8) :: coords(N,3)
              integer(kind=8) :: i
              do i = 1, N
                  print *, 'Atom : ',i, '  coords :  ', coords(i,:) 
              end do
      end subroutine print_coords

!--------------------------------------------------------------------------------------------------------------------

      subroutine save_coords(N_at, coords, output_file_coords)
             ! Save coords of one walker on the output file
             integer(kind=8), intent(in) :: N_at
             real(kind=8), intent(in) :: coords(N_at,3)
             character(len=*), intent(in) :: output_file_coords
             integer(kind=8) :: i

             !n = len(output_file)
             !allocate(character (n) :: output_file_used)
             !output_file_used = output_file

             open(unit=10, file=output_file_coords, action='Write')
             write(10,*) 'Atom          coords(x,y,z)'
             do i=1,N_at
                 write(10, '(i6,4x,f16.13,4x,f16.13,4x,f16.13)') i, &
                         coords(i,1), coords(i,2), coords(i,3)
             end do
             close(unit=10)

      end subroutine save_coords

!---------------------------------------------------------------------------------------------------------------------

      subroutine save_rank_energies(my_rank, M, energies)
              ! Save the energy averged over all walkers 
              ! for each time-step
              integer, intent(in) :: my_rank 
              integer(kind=8), intent(in) :: M
              real(kind=8), intent(in) :: energies(M)
              character(len=20) :: my_file
              integer(kind=8) :: i, out_unit

              out_unit = 5 + my_rank

              write(my_file,'(a10,i1,a4)') 'my_energy_',my_rank,'.txt'
              open(unit=10, file=my_file, action='write')
              write(10,*) 'Step         Average_energy'
              do i =1,M
                  write(10,'(i10,4x,f15.8)') i, energies(i) 
              end do
              close(10)

      end subroutine

!---------------------------------------------------------------------------------------------------------------------

      subroutine save_rank_dmc(my_rank, it, Nt, ene, acc_prob, out_file)
              integer, intent(in) :: my_rank 
              integer(kind=8), intent(in) :: it
              integer(kind=8), intent(in) :: Nt(it)
              real(kind=8), intent(in) :: ene(it)
              real(kind=8), intent(in) :: acc_prob(it)
              character(len=50) :: out_file
              integer(kind=8) :: i
              integer :: io_unit=10
              character(len=25) :: my_file


        open(unit = io_unit, file=out_file, action='write')
        write(io_unit, *) 'Step       Num walks        Ene     acc_prob'
        do i = 1, it
        write(io_unit,'(i10,4x,i10,4x,f15.10,4x,f15.10)') &
              i,Nt(i),ene(i),acc_prob(i) 
        end do
        close(io_unit)

      end subroutine save_rank_dmc 

!---------------------------------------------------------------------------------------------------------------------

      subroutine save_output(N_cyc, energies, err_energies, acc_ratios)
              ! save the output on the given output file
              integer(kind=8), intent(in) :: N_cyc
              real(kind=8), intent(in) :: energies(N_cyc)
              real(kind=8), intent(in) :: err_energies(N_cyc)
              real(kind=8), intent(in) :: acc_ratios(N_cyc)
              integer(kind=8) :: i

              open(unit=10, file=output_file, action='write')
              write(10,*) 'Cycle        E        err_E       a_r'
              do i = 1,N_cyc
                  write(10,'(i10,4x,f15.8,4x,f15.8,4x,f15.8)') i, &
                          energies(i), err_energies(i), acc_ratios(i)
              end do
              close(10)
      end subroutine save_output

!--------------------------------------------------------------------------------------------------------------------------

      subroutine save_rho_rad(Nl, r_mesh, rho_rad, eigenvec)
              integer(kind=8), intent(in) :: Nl
              real(kind=8), intent(in) :: r_mesh(Nl), rho_rad(Nl)
              real(kind=8), intent(in) :: eigenvec(Nl)
              integer:: i

              open(unit=io_unit, file=output_file, action='Write')
              write(io_unit,*) 'r                      rho(r)'
              do i = 1,Nl
              write(io_unit,'(f15.10,4x,f15.10,4x,f15.10)') r_mesh(i),&
                            rho_rad(i), eigenvec(i)
              end do
              close(io_unit)
      end subroutine

!------------------------------------------------------------------------------------------------------------------------------

      subroutine save_eigenvectors(Nl, M, r_mesh, eigenvecs, eigenvec_out_file)
              integer(kind=8), intent(in) :: Nl
              integer, intent(in) ::  M
              real(kind=8) :: r_mesh(Nl)
              real(kind=8), intent(in) :: eigenvecs(Nl,M)
              character(len=50) :: eigenvec_out_file
              integer(kind=8) :: i

              open(unit=io_unit, file=eigenvec_out_file, action='Write')
              do i = 1,Nl
                  write(io_unit,*) r_mesh(i), eigenvecs(i,:)
              end do
              close(io_unit)

      end subroutine save_eigenvectors


      end module io_handler
