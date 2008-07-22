subroutine Broadcast_elements(my_id,element_list)
!----------------------------------------------------------
! subroutine to broadcast all the nodes in the point_list
!----------------------------------------------------------
use data_structure
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
include 'mpif.h'               ! MPI fortran include file
type (type_element)      :: anelement
integer                  :: my_id, ife, ierr, position, bufsize, IDBL_EXT, INT_EXT, ILOG_EXT
integer, allocatable     :: buffer(:)

! type type_element                                    ! type definition for one elements
!    integer :: vertex(n_vertex_max)                   ! the nodes of the corners
!    integer :: neighbours(n_vertex_max)               ! the neighbouring elements
!    real*8  :: size(n_vertex_max,n_order+1)           ! the size of the vectors at each vertex of the element
!  endtype type_element

!  type type_element_list                              ! type definition for a list of elements
!    integer :: n_elements                             ! the number of elements in the list
!    type (type_element)  :: element(n_elements_max)   ! the list of elements
!  endtype type_element_list

!----------------------------------- one line would be enough if only MPI_TYPE_STRUCT would work on IXIA
!call MPI_BCAST(fe_list%fe(1:fe_list%nfe),fe_list%nfe,MPI_element,0,MPI_COMM_WORLD,ierr)



call MPI_TYPE_EXTENT(MPI_INTEGER,INT_EXT,ierr)
call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,IDBL_EXT,ierr)

call MPI_BCAST(element_list%n_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

bufsize = element_list%n_elements * (2*n_vertex_max*INT_EXT + n_vertex_max*(n_order+1)*IDBL_EXT)

allocate(buffer(bufsize/ INT_EXT))

if (my_id .eq. 0) then
  position = 0

  do ife=1,element_list%n_elements

    anelement = element_list%element(ife)

    call MPI_PACK(anelement%vertex,n_vertex_max,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anelement%neighbours,n_vertex_max,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anelement%size,n_vertex_max*(n_order+1),MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  enddo

endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

  position = 0
  do ife=1,element_list%n_elements

    call MPI_UNPACK(buffer,bufsize,position,anelement%vertex,n_vertex_max,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anelement%neighbours,n_vertex_max,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anelement%size,n_vertex_max*(n_order+1),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    element_list%element(ife) = anelement

  enddo

endif

return
end
