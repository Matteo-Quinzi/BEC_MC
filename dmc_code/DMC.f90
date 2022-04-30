program DMC_BEC
      ! use ulimit -s unlimited in the bash to avoid segmentation faults
      use mpi 
      use bec_dmc
      use bec_vmc
      implicit none

      !mpi_variables
      integer :: ierr
      integer :: nprocs, my_rank

      !initialize mpi
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

      print *, 'Total number of procs ',nprocs
      print *, 'My rank ',my_rank
      
      call mpi_finalize(ierr)


end program DMC_BEC

