subroutine construct_matrix(my_id,node_list,element_list,local_elms,n_local_elms,index_min,index_max, &
                            xpoint2,psi_axis,psi_bnd,Z_xpoint)
!---------------------------------------------------------------
! collect the element matrices into one large sparse matrix
! in coordinate format
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

real*8, allocatable :: rhs_loc(:)
integer :: my_id, local_elms(*), n_local_elms, index_min, index_max,index_min_loc, index_max_loc
real*8  :: ELM(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: RHS(n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: ELM2(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: RHS2(n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: zbig, psi_axis, psi_bnd, Z_xpoint, T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
real*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
real*8  :: Vpar0_pol_R, Vpar0_pol_Z, Vpol_R, Vpol_Z, znormal_R, znormal_Z
real*8  :: Vpar0_perp, Vpol_perp, Btot, cs_fraction, ratio
real*8  :: grad_s, grad_psi, u0_s, u0_t, u0_x, u0_y
integer :: i_bnd, i, in, ife, iv, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
integer :: index_i, index_large_i, index_large_k, index_node, index_node1, index_node2, i_order, k_order, ic, ielm, ierr
integer :: ijA_position,ijA_position2, nz_AA2, n_AA2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
logical :: xpoint2

!write(*,*) '****************************************'
!write(*,*) '*  construct matrix                    *'
!write(*,*) '****************************************'
!write(*,*) ' n_elements (local)       : ',my_id,n_local_elms
!write(*,*) ' index_min,index_max      : ',my_id,index_min,index_max

i_bnd = 0

do i=1, n_local_elms

  ielm = local_elms(i)

  do iv=1,n_vertex_max

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .eq. 1) i_bnd = i_bnd + 1
    if (node_list%node(inode)%boundary .eq. 2) i_bnd = i_bnd + 1
    if (node_list%node(inode)%boundary .eq. 3) i_bnd = i_bnd + 2

    index_min_loc = min(index_min_loc, minval(node_list%node(iv)%index))
    index_max_loc = max(index_max_loc, maxval(node_list%node(iv)%index))

  enddo

enddo

if (.not. allocated(A_glob))    allocate(A_glob(nz_glob))
if (.not. allocated(irn_glob))  allocate(irn_glob(nz_glob))
if (.not. allocated(jcn_glob))  allocate(jcn_glob(nz_glob))

if (allocated(rhs_glob))        deallocate(rhs_glob)

allocate(rhs_glob(ndof_glob))

if (.not. allocated(rhs_loc))   allocate(rhs_loc(ndof_glob))

irn_glob = 0
jcn_glob = 0
A_glob   = 0.d0
RHS_glob = 0.d0
RHS_loc  = 0.d0

do ife =1, n_local_elms

  ielm = local_elms(ife)

  element = element_list%element(ielm)

  do iv = 1, n_vertex_max

    inode          = element%vertex(iv)
    nodes%node(iv) = node_list%node(inode)

  enddo

  if (n_tor .gt. 3) then
    call element_matrix_fft(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)      ! use fft for toroidal integration
  else
    call element_matrix(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)           ! use direct integration
  endif

!------------------------------------------------------- comparing two versions of element_matrix
!  if (ife .eq. n_local_elms/2) then
!    call element_matrix_fft(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM2, RHS2)
!    call element_matrix(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)
!    do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
!      if (abs(RHS(i)-RHS2(i))/(abs(RHS(i))+abs(RHS2(i))+1.d0) .gt. 1.d-12) then
!        write(*,'(i3,A,i6,3e16.8)') my_id,' RHS : ',i,RHS(i),RHS2(i),RHS(i)-RHS2(i)
!      endif
!    enddo
!    do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
!      do j=1,n_tor*n_vertex_max*(n_order+1)*n_var
!        if (abs(ELM(i,j)-ELM2(i,j))/(abs(ELM(i,j))+abs(ELM2(i,j))+1.d0) .gt. 1.d-10) then
!          write(*,'(i3,A,2i6,3e16.8)') my_id,' ELM : ',i,j,ELM(i,j),ELM2(i,j),ELM(i,j)-ELM2(i,j)
!        endif
!     enddo
!    enddo
!  endif

!    if (ife .eq. n_local_elms/2) then
!    do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
!      write(*,'(i3,A,i6,4e16.8)') my_id,' RHS : ',i,RHS(i)
!    enddo
!    do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
!      do j=1,n_tor*n_vertex_max*(n_order+1)*n_var
!        write(*,'(i3,A,2i6,8e16.8)') my_id,' ELM : ',i,j,ELM(i,j)
!      enddo
!    enddo
!    endif

  do i=1,n_vertex_max

    inode1         = element%vertex(i)

    do i_order = 1, n_order+1

      index_node1 = node_list%node(inode1)%index(i_order)

      index_large_i = n_tor * n_var * (index_node1 - 1)

      if ((index_node1 .ge. index_min) .and. (index_node1 .le. index_max)) then

        do j = 1, n_var * n_tor

          index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix

          rhs_loc(index_large_i+j) = rhs_loc(index_large_i+j) + RHS(index_ij)

          do k=1,n_vertex_max

            knode         = element%vertex(k)

            do k_order = 1, n_order+1

              index_node2 = node_list%node(knode)%index(k_order)

              index_large_k = n_tor * n_var * (index_node2 - 1)

              call locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position)

              do l = 1, n_var * n_tor

                index_kl = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l   ! index in the ELM matrix

                ilarge2 = ijA_position - 1 + (j-1) * n_var*n_tor + l

                irn_glob(ilarge2) = index_large_i   + j
                jcn_glob(ilarge2) = index_large_k   + l
                A_glob(ilarge2)   = A_glob(ilarge2) + ELM(index_ij,index_kl)

              enddo

            enddo

          enddo
        enddo

      endif

    enddo
  enddo

enddo

!write(*,*) ' nz_aa : ',my_id, nz_aa

!----------------------- boundary conditions

zbig = 1.d10

do i=1, n_local_elms

  ielm = local_elms(i)

  do iv=1, n_vertex_max

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .ne. 0) then

      do in=1, n_tor

        do k=1, n_var

!------------------------------------ the open field lines (in case of x-point grid)
          if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

!            if ((k .eq.   1) .or. (k .eq. 92) .or. (k .eq. 3)  .or. &
!                (k .eq.  4)  .or. (k .eq. 95)  .or. (k .eq. 96) .or. (k .eq.97) ) then

            if ((k .eq.   1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                (k .eq.  4)  .or. (k .eq. 5) .or. (k .eq. 6) .or. (k .eq.97) ) then

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

            if (k .eq. 7) then

              index_node  = node_list%node(inode)%index(1)             ! position of value
              index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

              T0        = node_list%node(inode)%values(1,1,6)
              Vpar0     = node_list%node(inode)%values(1,1,7)
              BigR      = node_list%node(inode)%x(1,1)
              dT0_ds    = node_list%node(inode)%values(1,2,6)
              dVpar0_ds = node_list%node(inode)%values(1,2,7)
              dBigR_ds  = node_list%node(inode)%x(2,1)

              ps0_s     = node_list%node(inode)%values(1,2,1)
              ps0_t     = node_list%node(inode)%values(1,3,1)

              U0_s      = node_list%node(inode)%values(1,2,2)
              U0_t      = node_list%node(inode)%values(1,3,2)

              R_s       = node_list%node(inode)%x(2,1)
              R_t       = node_list%node(inode)%x(3,1)
              Z_s       = node_list%node(inode)%x(2,2)
              Z_t       = node_list%node(inode)%x(3,2)

              xjac  =  R_s*Z_t - R_t*Z_s
              ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
              ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

              u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
              u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

              direction = + ps0_x / abs(ps0_x)             ! temporary solution for lower x-point only

              grad_psi = sqrt(ps0_x**2 + ps0_y**2)

              Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

              if (in .eq. 1) then

!                write(*,'(A,3e14.6,A,e14.6)') ' Boundary : ',Vpar0, -BigR**2 * u0_s/ps0_s, direction*sqrt(GAMMA*T0)/Btot,&
!                                              ' error : ',Vpar0 - BigR**2 * u0_s/ps0_s - direction*sqrt(GAMMA*T0)/Btot

              endif

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)

                index_large_i = n_tor * n_var * (index_node - 1)

                ku = 2
                kv = 7
                kT = 6

                ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                ilarge_vT  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                ilarge_vus = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (ku-1)*n_tor + in

                irn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                jcn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                A_glob(ilarge_vv)   =  zbig

                irn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                jcn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vT)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                irn_glob(ilarge_vus) =  n_tor * n_var * (index_node -1) + (kv-1)*n_tor + in
                jcn_glob(ilarge_vus) =  n_tor * n_var * (index_node2-1) + (ku-1)*n_tor + in
                A_glob(ilarge_vus)   = - zbig * BigR**2 / ps0_s

                RHS_glob(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                         Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*sqrt(GAMMA*T0))

              endif

              index_node = node_list%node(inode)%index(2)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                kv = 7
                kT = 6

                ilarge_vv = ijA_position - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                ilarge_vT = ijA_position - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in

                irn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                jcn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
 !               A_glob(ilarge_vv)   = zbig

                irn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                jcn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kT-1)*n_tor + in
 !               A_glob(ilarge_vT)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

 !               RHS_glob(n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in) = &
 !                       Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction)

              endif

            endif

          endif

!------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
          if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

!            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
!                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 96) .or. (k .eq. 7) ) then
            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) .or. (k .eq. 7) ) then

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

call MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!write(*,*) '******** end construct matrix **********'

deallocate(RHS_loc)

return
end
