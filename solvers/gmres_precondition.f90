subroutine gmres_precondition(x,y,i_tor,my_id,my_id_n,MPI_COMM_MASTER,MPI_COMM_N)
!-----------------------------------------------------------------------
! solve step of the local matrices for each toroidal harmonic
! (preconditioner for gmres)
!-----------------------------------------------------------------------
use parameters
use mumps_module
use pastix_module
use global_distributed_matrix
implicit none
include 'mpif.h'
integer             :: my_id, my_id_n, MPI_COMM_MASTER, MPI_COMM_N, ierr, i, i_tor(*), n_dof
real*8              :: x(*), y(*)
real*8, allocatable :: y_tmp(:), Rsnd_buffer(:)
integer             :: index_snd, n_loc_n, n_cpu, n_cpu_n, M_cpu, ifactor, in, j, idisp, index_rcv
integer, allocatable :: send_counts(:), send_disp(:), recv_counts(:), recv_disp(:)
real*8              :: t1, t2, t3, t4, t5, t6

write(*,*) my_id,my_id_n,' GMRES preconditioning ',MPI_COMM_WORLD

call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)     ! the number of cpus

call cpu_time(t1)

n_dof    = ndof_glob

n_loc_n  = n_dof / n_tor
M_cpu    = n_cpu / ((n_tor+1)/2)

ifactor = 2
if (my_id   .eq. 0) ifactor = 1
if (my_id_n .ne. 0) ifactor = 0

if (my_id .eq. 0) then

  allocate(Rsnd_buffer(n_dof),send_counts(n_cpu/M_cpu),send_disp(n_cpu/M_cpu))
  Rsnd_buffer(1:n_loc_n) = x(1:n_dof:n_tor)

  do in=2, (n_tor+1)/2

    index_snd = n_loc_n + (in-2)*2*n_loc_n

    Rsnd_buffer(index_snd+1:index_snd+2*n_loc_n:2) = x(2*in-2:n_dof:n_tor)
    Rsnd_buffer(index_snd+2:index_snd+2*n_loc_n:2) = x(2*in-1:n_dof:n_tor)

  enddo

  send_counts(1) = n_loc_n
  send_disp(1)   = 0
  idisp          = send_counts(1)

  do j=2,n_cpu/M_cpu

    send_counts(j) = 2*n_loc_n
    send_disp(j)   = idisp
    idisp          = idisp + send_counts(j)

  enddo

endif

if (associated(mumps_par%rhs)) deallocate(mumps_par%rhs)

if (my_id_n .eq. 0) then

  allocate(mumps_par%rhs(ifactor*n_loc_n))

  call mpi_scatterv(Rsnd_buffer,send_counts,send_disp,MPI_DOUBLE_PRECISION, &
                    mumps_par%rhs,ifactor*n_loc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_MASTER,ierr)

endif


if (my_id .eq. 0) deallocate(Rsnd_buffer)


if (use_mumps) then

  mumps_par%JOB = 3                                   ! Solve

  call DMUMPS(mumps_par)

else

  pastix_iparm(2) = 5
  pastix_iparm(3) = 6
  pastix_iparm(6) = 1          ! refinement : max number of iterations

  pastix_iparm(31) = 2
  pastix_iparm(35) = 1         ! numthreads   ! number of threads
  pastix_iparm(37) = 1
  pastix_iparm(39) = 0         ! right hand side (0 : use RHS)
  pastix_iparm(41) = 1

  pastix_dparm(6)  = 1.d-20    ! error level refinement
  pastix_dparm(11) = 1.d-32    ! pivot threshold?

  write(*,*) my_id, my_id_n,' PRECONDITIONING using PASTIX '

  if (.not. associated(mumps_par%rhs)) then
    write(*,*) ' gmres: RHS not allocated!',my_id, my_id_n
    allocate(mumps_par%rhs(mumps_par%n))
  endif

  call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)

  write(*,'(2i3,A,2e16.8)') my_id,my_id_n,' precond : rhs before : ',maxval(mumps_par%rhs),minval(mumps_par%rhs)

  call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)


endif

write(*,'(2i3,A,2e16.8)') my_id,my_id_n,' precond : rhs after : ',maxval(mumps_par%rhs),minval(mumps_par%rhs)


if (my_id_n .eq. 0) then

  allocate(y_tmp(n_dof))
  allocate(recv_counts(n_cpu/M_cpu),recv_disp(n_cpu/M_cpu))

  y_tmp(1:n_dof) = 0.d0

  recv_counts(1) = n_loc_n
  do i=2,(n_tor+1)/2
    recv_counts(i) = 2*n_loc_n
  enddo

  recv_disp(1) = 0
  do i=2,(n_tor+1)/2
    recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
  enddo

  call mpi_gatherv(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION, &
                   y_tmp,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,0,MPI_COMM_MASTER,ierr)

  if (my_id .eq. 0) then

    y(1:n_dof:n_tor) = y_tmp(1:n_loc_n)

    do in=2, (n_tor+1)/2

      index_rcv = n_loc_n + (in-2)*2*n_loc_n

      y(2*in-2:n_dof:n_tor) = y_tmp(index_rcv+1:index_rcv+2*n_loc_n:2)
      y(2*in-1:n_dof:n_tor) = y_tmp(index_rcv+2:index_rcv+2*n_loc_n:2)

    enddo

  endif

  deallocate(y_tmp, recv_counts, recv_disp)
endif

call cpu_time(t2)
!if (my_id .eq. 0) write(*,'(i3,A,f14.6)') my_id,' PRECON TOTAL : ',t2-t1

return
end
