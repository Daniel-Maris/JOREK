subroutine initialise_pastix
!***********************************************************************
! initialise the pastix library                                        *
!***********************************************************************
use pastix_module
use mumps_module

if (.not. use_pastix) return

write(*,*) '***********************************'
write(*,*) '* initialise PastiX               *'
write(*,*) '***********************************'

pastix_iparm(1)     = 0                                     ! insert default values
pastix_iparm(2)     = 0                                     ! initializse
pastix_iparm(3)     = 0

pastix_iparm(31) = 2
pastix_iparm(35) = 1                  ! thread/mpi
pastix_iparm(39) = 0
pastix_iparm(41) = 1

!call my_pastix(MPI_COMM_N,mumps_par%n,mumps_par%irn,mumps_par%jcn,mumps_par%A,   &
!                   pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

pastix_iparm(6)  = 0
pastix_dparm(6)  = 1.e-12
pastix_dparm(11) = 1.e-32             ! pivot threshold?

return
end