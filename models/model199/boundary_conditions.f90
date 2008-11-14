subroutine boundary_conditions(my_id,node_list,element_list,local_elms,n_local_elms,index_min,index_max, &
                               xpoint2,psi_axis,psi_bnd,Z_xpoint)			    
!---------------------------------------------------------------
! add the boundary condition to the global matrix
!---------------------------------------------------------------
use data_structure
use global_distributed_matrix
use phys_module

implicit none
include 'mpif.h'

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node_list)    :: nodes

integer :: my_id, local_elms(*), n_local_elms, index_min, index_max, index_min_loc, index_max_loc
real*8  :: zbig, psi_axis, psi_bnd, Z_xpoint, T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
real*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
real*8  :: Vpar0_pol_R, Vpar0_pol_Z, Vpol_R, Vpol_Z, znormal_R, znormal_Z
real*8  :: Vpar0_perp, Vpol_perp, Btot, cs_fraction, ratio
real*8  :: grad_s, grad_psi, u0_s, u0_t, u0_x, u0_y
integer :: i_bnd, i, in, ife, iv, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
integer :: index_i, index_large_i, index_large_k, index_node, index_node1, index_node2, i_order, k_order, ic, ielm, ierr
integer :: ijA_position,ijA_position2, nz_AA2, n_AA2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
logical :: xpoint2

zbig = 1.d12

do i=1, n_local_elms

  ielm = local_elms(i)

  do iv=1, n_vertex_max

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .ne. 0) then

      do in=1, n_tor

        do k=1, n_var

!------------------------------------ the open field lines (in case of x-point grid)
          if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then

              index_node = node_list%node(inode)%index(1)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)   = zbig

              endif

              index_node = node_list%node(inode)%index(2)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)    = zbig

              endif

            endif

          endif

!------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
          if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then

              index_node = node_list%node(inode)%index(1)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)   = zbig

              endif

              index_node = node_list%node(inode)%index(3)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)    = zbig

              endif

            endif

          endif

        enddo

      enddo
    endif
  enddo
enddo

return
end
