subroutine gmres_matrix_vector(x,y,my_id,my_id_n)
!-----------------------------------------------------------------------
! sparse matrix vector product using coordinate scheme
! to be called only on all of MPI_COMM_WORLD
!
! result (y) is only known on id=0
!-----------------------------------------------------------------------
use parameters
use global_distributed_matrix
implicit none
include 'mpif.h'

real*8              :: x(*), y(*), t1, t2, t3, t4, t5
real*8, allocatable :: y_tmp(:)
integer             :: n, i, ir, jc, ierr, my_id, my_id_n
integer             :: n_blocksize, n_blocks, iA_start, ix_start, iy_start

!write(*,*) my_id,my_id_n,' GMRES matrix_vector ',ndof_glob
call cpu_time(t1)

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call cpu_time(t2)

call MPI_BCAST(x,ndof_glob,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

call cpu_time(t3)

allocate(y_tmp(ndof_glob))

y_tmp(1:ndof_glob) = 0.d0

n_blocksize = n_tor * n_var
n_blocks    = nz_glob/n_blocksize**2

do i=1, n_blocks

  iA_start = (i-1) * n_blocksize**2

  ix_start = jcn_glob(iA_start+1)
  iy_start = irn_glob(iA_start+1)

  call dgemv('T',n_blocksize,n_blocksize,1.d0,A_glob(iA_start+1),n_blocksize,x(ix_start),1,1.d0,y_tmp(iy_start),1)

enddo
call cpu_time(t4)

!do i=1,nz_glob
!  ir = irn_glob(i)
!  jc = jcn_glob(i)
!  y_tmp(ir) = y_tmp(ir) + A_glob(i) * x(jc)
!enddo

y(1:ndof_glob) = 0.d0

call MPI_Reduce(y_tmp,y,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

call cpu_time(t5)

deallocate(y_tmp)

!write(*,'(A,i3,3f14.6)') ' M-V timing  barrier: ',my_id,t2-t1
!write(*,'(A,i3,3f14.6)') ' M-V timing  bcast  : ',my_id,t3-t2
!write(*,'(A,i3,3f14.6)') ' M-V timing  dgemv  : ',my_id,t4-t3
!write(*,'(A,i3,3f14.6)') ' M-V timing  reduce : ',my_id,t5-t4
!write(*,'(A,i3,3f14.6)') ' M-V timing  TOTAL  : ',my_id,t5-t1

return
end