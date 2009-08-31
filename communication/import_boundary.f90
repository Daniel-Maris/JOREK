subroutine import_boundary(node_List,boundary_list) 

use data_structure 
implicit none

type(type_boundary_list) :: boundary_list
type(type_node_list)     :: node_list

integer :: i,iv1,iv2,idir1,idir2
 
open(22,file='boundary.txt')
 
read(*,*) boundary_list%n_boundary

do i=1,boundary_list%n_boundary
    iv1   = boundary_list%boundary(i)%vertex(1)
    iv2   = boundary_list%boundary(i)%vertex(2)
    idir1 = boundary_list%boundary(i)%direction(1,2)
    idir2 = boundary_list%boundary(i)%direction(2,2)
    read(22,'(3i6,12e16.8)') i,iv1,iv2, node_list%node(iv1)%x(1,1:2),          &
                                         node_list%node(iv1)%x(idir1,1:2),      &
					 boundary_list%boundary(i)%size(1:2,2), &
					
					 node_list%node(iv2)%x(1,1:2),          &
                                         node_list%node(iv2)%x(idir2,1:2),      &
					 boundary_list%boundary(i)%size(1:2,2) 
					
enddo

close(22)
return
end
