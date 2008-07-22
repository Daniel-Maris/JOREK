subroutine scale_global_matrix
!***********************************************************************
!* column scaling of the global matrix                                 *
!***********************************************************************
use global_distributed_matrix
implicit none
include 'mpif.h'
real*8, allocatable :: column_local(:)
integer :: k,j,ierr

if (allocated(column_scaling))  deallocate(column_scaling)
if (allocated(column_local))    deallocate(column_local)
allocate(column_scaling(ndof_glob),column_local(ndof_glob))

column_local = 1.d-20;   column_scaling = 1.d-20

do k=1,nz_glob
  j = jcn_glob(k)
  column_local(j) = max(column_local(j),abs(A_glob(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,ndof_glob,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

do k=1,nz_glob
  j = jcn_glob(k)
  A_glob(k) = A_glob(k) / column_scaling(j)
enddo

deallocate(column_local)

return
end
