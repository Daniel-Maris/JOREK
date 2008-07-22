subroutine solve_pastix_all_dist(my_id,my_index_min,my_index_max)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using pastix with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use mumps_module
use pastix_module
use global_distributed_matrix
implicit none
include 'mpif.h'

integer                  :: my_index_min, my_index_max       ! global index_min, index_max for this cpu
real*8,allocatable       :: column_local(:)
integer, allocatable     :: pastix_loc2glb(:)
real*8                   :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1
integer                  :: i, k, j, ierr, my_id

write(*,*) my_id,'*********************************'
write(*,*) my_id,'*  solve global matrix (PastiX) *'
write(*,*) my_id,'*********************************'

mumps_par%A_loc   => A_glob(1:nz_glob)
mumps_par%irn_loc => irn_glob(1:nz_glob)
mumps_par%jcn_loc => jcn_glob(1:nz_glob)
mumps_par%rhs     => rhs_glob(1:nz_glob)

mumps_par%n      = ndof_glob
mumps_par%nz_loc = nz_glob

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(mumps_par%N + 1))

call coicsr(mumps_par%N,mumps_par%nz_loc,1,mumps_par%A_loc,mumps_par%IRN_loc,mumps_par%JCN_loc,sparskit_work)

deallocate(sparskit_work)

allocate(pastix_loc2glb(my_index_max-my_index_min+1))

do i=1,my_index_max-my_index_min+1
  pastix_loc2glb(i) = my_index_min + i - 1
enddo

if (allocated(column_scaling))  deallocate(column_scaling)
if (allocated(column_local))    deallocate(column_local)
allocate(column_scaling(mumps_par%N),column_local(mumps_par%N))

column_local = 1.d-20;   column_scaling = 1.d-20
do k=1,mumps_par%nz_loc
  j = mumps_par%jcn_loc(k)
  column_local(j) = max(column_local(j),abs(mumps_par%A_loc(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
do k=1,mumps_par%nz_loc
  j = mumps_par%jcn_loc(k)
  mumps_par%A_loc(k) = mumps_par%A_loc(k) / column_scaling(j)
enddo

call cpu_time(t_analysis_0)

pastix_iparm(2) = 1
pastix_iparm(3) = 3
pastix_iparm(6) = 0          ! refinement : max number of iterations

pastix_iparm(31) = 2
pastix_iparm(35) = 1         !   numthreads : number of threads
pastix_iparm(39) = 0         ! right hand side (0 : use RHS)
pastix_iparm(37) = 1
pastix_iparm(41) = 1

pastix_dparm(6)  = 1.d-20    ! error level refinement
pastix_dparm(11) = 1.d-32    ! pivot threshold?

!call dpastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn_loc, mumps_par%irn_loc, mumps_par%A_loc, &
!                     pastix_loc2glb, pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

call cpu_time(t_analysis_1)
if (my_id .eq. 0)  write(*,'(A,f8.3)') ' PASTIX, analysis  : ',t_analysis_1-t_analysis_0

call cpu_time(t_fact_0)

pastix_iparm(2) = 4
pastix_iparm(3) = 6
pastix_iparm(6) = 1          ! refinement : max number of iterations

!call dpastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn_loc, mumps_par%irn_loc, mumps_par%A_loc, &
!                     pastix_loc2glb, pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

call cpu_time(t_fact_1)

if (my_id .eq. 0) write(*,'(A,f8.3)')  ' PASTIX, fact/solv : ',t_fact_1-t_fact_0

do k=1,mumps_par%n
  deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
!  write(*,*) k,deltas(k)
enddo

return
end
