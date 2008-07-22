subroutine global_matrix_structure(my_id,node_List,element_list,local_elms,n_local_elms,index_min,index_max)
!***********************************************************************
!* subroutine determines the position of the indices in the global     *
!* matrix                                                              *
!***********************************************************************

use data_structure
use global_distributed_matrix

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: flux_list

integer             :: local_elms(*), index_min, index_max, my_id, n_local_elms
integer             :: i, iv, ik, jv, jk, ielm, inode1, inode2, index1, index2, index1_local, index2_local
integer             :: j_larger, j, ibase, n_max

ndof_glob = -1
do inode1=1,node_list%n_nodes
  ndof_glob = max(ndof_glob,maxval(node_list%node(inode1)%index))
enddo
ndof_glob = ndof_glob * n_tor*n_var

n_max = 800

allocate(ijA_size(index_max-index_min+1))
allocate(irn_jcn(index_max-index_min+1,n_max))

ijA_size    = 0
irn_jcn = 0

do i=1,n_local_elms

  ielm = local_elms(i)

  do iv = 1, n_vertex_max

    inode1 = element_list%element(ielm)%vertex(iv)

    do ik = 1, n_degrees

      index1       = node_list%node(inode1)%index(ik)
      index1_local = index1 - index_min + 1

      if ((index1 .ge. index_min) .and. (index1 .le. index_max)) then

        do jv = 1,n_vertex_max

          inode2 = element_list%element(ielm)%vertex(jv)

          do jk = 1, n_degrees

            index2       = node_list%node(inode2)%index(jk)
            index2_local = index2 - index_min + 1

            if (ijA_size(index1_local) .eq. 0) then

              ijA_size(index1_local) = 1
              irn_jcn(index1_local,1) = index2

            elseif (index2 .gt. irn_jcn(index1_local,ijA_size(index1_local))) then

              irn_jcn(index1_local,ijA_size(index1_local)+1) = index2
              ijA_size(index1_local) = ijA_size(index1_local) + 1

            else

             do j = 1, ijA_size(index1_local)

                if (index2 .le. irn_jcn(index1_local,j) ) then

                  j_larger = j
                  exit

                endif

              enddo

              if (index2 .ne. irn_jcn(index1_local,j_larger) ) then

                do j=ijA_size(index1_local), j_larger, -1

                  irn_jcn(index1_local,j+1) = irn_jcn(index1_local,j)

                enddo

                irn_jcn(index1_local,j_larger) = index2
                ijA_size(index1_local) = ijA_size(index1_local) + 1

                if (ijA_size(index1_local) .gt. n_max) then
                  write(*,*) ' FATAL error : irn_jcn too small ',ijA_size(index1_local)
                endif

              endif

            endif

          enddo
        enddo

      endif

    enddo
  enddo

enddo

allocate(ijA_index(index_max-index_min+1,n_max))

ibase = 0
do i=1,index_max-index_min+1

  do j=1,ijA_size(i)

    ijA_index(i,j) = ibase + 1

    ibase = ibase + (n_tor*n_var)**2

  enddo

enddo

n_glob  = (index_max-index_min+1) * n_tor * n_var

nz_glob = ijA_index(index_max-index_min+1,ijA_size(index_max-index_min+1)) + (n_tor*n_var)**2 - 1

write(*,*) my_id,' size matrices : n, nz = ',n_glob, nz_glob

return
end
