subroutine boundary_check(node_list,rhs)
use data_structure

implicit none

integer :: i,j,in, index1, index2, index3, index_node1, index_node2, index_node3
real*8  :: rhs(*)

type (type_node_list)    :: node_list

write(*,*) '************************************'
write(*,*) '* boundary check                   *'
write(*,*) '************************************'

do i=1,node_list%n_nodes

  if (node_list%node(i)%boundary .ne. 0) then
    write(*,'(A,2i5,20e14.6)') ' bnd node : ',i,node_list%node(i)%boundary, &
              node_list%node(i)%x(1,:),node_list%node(i)%values(1,1:4,1),node_list%node(i)%values(1,1:4,3)
  endif

  index_node1 = node_list%node(i)%index(1)
  index_node2 = node_list%node(i)%index(2)
  index_node3 = node_list%node(i)%index(3)

  if (node_list%node(i)%boundary .eq. 1) then

    do j=1,n_var

      if (j .ne. 4) then

        do in=1,n_tor

          index1 = n_tor*n_var * (index_node1 - 1) + n_tor*(j-1) + in
          index2 = n_tor*n_var * (index_node2 - 1) + n_tor*(j-1) + in

          if ( (abs(rhs(index1)) .gt. 1d-8) .or. (abs(rhs(index2)) .gt. 1d-8) ) then
            write(*,*) ' PROBLEM BOUNDARY CONDITION (1): ',i,j,in,rhs(index1),rhs(index2)
          endif
        enddo

      endif

    enddo

  endif

  if (node_list%node(i)%boundary .eq. 2) then

    do j=1,n_var

      if (j .ne. 4) then

        do in=1,n_tor

          index1 = n_tor*n_var * (index_node1 - 1) + n_tor*(j-1) + in
          index3 = n_tor*n_var * (index_node3 - 1) + n_tor*(j-1) + in

          if ( (abs(rhs(index1)) .gt. 1d-8) .or. (abs(rhs(index3)) .gt. 1d-8) ) then
            write(*,*) ' PROBLEM BOUNDARY CONDITION (2): ',i,j,in,node_list%node(i)%values(in,1,j),node_list%node(i)%values(in,3,j)
          endif
        enddo

      endif

    enddo

  endif

  if (node_list%node(i)%boundary .eq. 3) then

    do j=1,n_var

      if (j .ne. 4) then

        do in=1,n_tor

          index1 = n_tor*n_var * (index_node1 - 1) + n_tor*(j-1) + in
          index2 = n_tor*n_var * (index_node2 - 1) + n_tor*(j-1) + in
          index3 = n_tor*n_var * (index_node3 - 1) + n_tor*(j-1) + in

          if ( (abs(rhs(index1)) .gt. 1d-8) .or. (abs(rhs(index2)) .gt. 1d-8) .or. (abs(rhs(index3)) .gt. 1d-8) ) then
            write(*,*) ' PROBLEM BOUNDARY CONDITION (3): ',i,j,in,node_list%node(i)%values(in,1,j),node_list%node(i)%values(in,2,j)
          endif
        enddo
      endif
    enddo
  endif
enddo
return
end
