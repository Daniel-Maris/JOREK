subroutine Poisson(my_id,itype,node_list,element_list,ivar_in,ivar_out,i_harm,xpoint)
!---------------------------------------------------------------
! collect the element matrices into one large sparse matrix
! in coordinate format
!---------------------------------------------------------------
use data_structure
use mumps_module
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8   :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1))
real*8   :: zbig, Z_xpoint, psi_axis, psi_bnd, psi_xpoint, R_xpoint, s_xpoint, t_xpoint
real*8   :: R_axis, Z_axis, s_axis, t_axis, amix
integer  :: ierr, my_id, itype, ivar_in, ivar_out, i_harm, i_elm_axis, i_elm_xpoint
integer  :: n_AA, nz_AA, nz_AA_old, n_border, ilarge, ife, iv, i,j,k,l
integer  :: n_elements, inode, index_large_i, knode, index_large_k, index_ij, index_kl, index, index_i
logical  :: xpoint

if (my_id .eq. 0) then

  write(*,*) '**************************************'
  write(*,*) '*            Poisson                 *'
  write(*,*) '**************************************'
  write(*,*) ' i_type     : ',itype
  write(*,*) ' n_elements : ',element_list%n_elements
  write(*,*) ' n_nodes    : ',node_list%n_nodes

  nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2

  n_border = 0
  do i=1,node_list%n_nodes
   if (node_list%node(i)%boundary .eq. 1) n_border = n_border+2
   if (node_list%node(i)%boundary .eq. 2) n_border = n_border+2
   if (node_list%node(i)%boundary .eq. 3) n_border = n_border+3
  enddo

  nz_AA = nz_AA + n_border

  n_AA  = node_list%n_nodes * (n_order+1)

  n_AA = 0
  do inode = 1, node_list%n_nodes
    n_AA = max(n_AA,node_list%node(inode)%index(4))
  enddo

  write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)
  write(*,*) ' number of boundary nodes: ',n_border

  if (.not. associated(mumps_par%A))     allocate(mumps_par%A(nz_AA))
  if (.not. associated(mumps_par%rhs))   allocate(mumps_par%rhs(n_AA))
  if (.not. associated(mumps_par%irn))   allocate(mumps_par%irn(nz_AA))
  if (.not. associated(mumps_par%jcn))   allocate(mumps_par%jcn(nz_AA))

  mumps_par%irn = 0
  mumps_par%jcn = 0
  mumps_par%A   = 0.d0
  mumps_par%RHS = 0.d0

  n_elements = element_list%n_elements

  ilarge=0

  if (itype .lt. 0) then
    call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis)

    psi_bnd = 0.d0
    if (xpoint) then
      call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)
      psi_bnd = psi_xpoint
    endif
  endif

  amix = 0.d0
!  if (itype .eq. -1) amix = 0.d0

  do ife =1, n_elements

    element = element_list%element(ife)

    do iv = 1, n_vertex_max

      inode     = element%vertex(iv)
      nodes(iv) = node_list%node(inode)

    enddo

    if (itype .eq. -1) then

      call element_matrix_GS(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

    elseif (itype .eq. -2) then

      call element_matrix_GS_inverse(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

    else

      call element_matrix_Poisson(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

    endif

    do i=1,n_vertex_max

      inode         = element%vertex(i)

      do j=1,n_order+1

        index_ij = (i-1)*(n_order+1) + j     ! index in the ELM matrix

        index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix

        mumps_par%rhs(index_large_i) = mumps_par%rhs(index_large_i) + (1.d0-amix) * RHS(index_ij)

        do k=1,n_vertex_max

          knode         = element%vertex(k)

          do l=1,n_order+1

            index_kl = (k-1)*(n_order+1) + l

            index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix

            ilarge = ilarge +1

            mumps_par%irn(ilarge) = index_large_i
            mumps_par%jcn(ilarge) = index_large_k
            mumps_par%A(ilarge)   = ELM(index_ij,index_kl)

          enddo
        enddo
      enddo
    enddo

  enddo

  nz_AA_old = nz_AA
  nz_AA = ilarge

  zbig = 1.d10

!----------------------- boundary conditions

  do i=1,node_list%n_nodes

    if (node_list%node(i)%boundary .ne. 0) then

      index_i = node_list%node(i)%index(1)  ! base index in the main matrix

      mumps_par%irn(ilarge+1) = index_i
      mumps_par%jcn(ilarge+1) = index_i
      mumps_par%A(ilarge+1)   = zbig
      ilarge = ilarge + 1

      if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then

        index_i = node_list%node(i)%index(2)  ! base index in the main matrix

        mumps_par%irn(ilarge+1) = index_i
        mumps_par%jcn(ilarge+1) = index_i
        mumps_par%A(ilarge+1)   = zbig
        ilarge = ilarge + 1
      endif

      if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3)) then

        index_i = node_list%node(i)%index(3)  ! base index in the main matrix

        mumps_par%irn(ilarge+1) = index_i
        mumps_par%jcn(ilarge+1) = index_i
        mumps_par%A(ilarge+1)   = zbig
        ilarge = ilarge + 1
      endif

    endif
  enddo

  nz_AA_old = nz_AA
  nz_AA     = ilarge

  mumps_par%n  = n_AA
  mumps_par%nz = nz_AA

endif

mumps_par%JOB = 6
mumps_par%SYM = 0
mumps_par%icntl(7) = 4

write(*,*) ' mumps : ',mumps_par%n, mumps_par%nz

call DMUMPS(mumps_par)

if (my_id .eq. 0) then

  do i=1,node_list%n_nodes

    index = node_list%node(i)%index(1)

!    write(*,'(A,2i5,8e16.8)') 'poisson : ',i,index,mumps_par%RHS(index),mumps_par%RHS(index+1),mumps_par%RHS(index+2),mumps_par%RHS(index+3)

    do k=1,n_order+1

      index = node_list%node(i)%index(k)

      node_list%node(i)%values(i_harm,k,ivar_out) = node_list%node(i)%values(i_harm,k,ivar_out) + mumps_par%RHS(index)

    enddo
  enddo

  deallocate(mumps_par%irn,mumps_par%jcn,mumps_par%A,mumps_par%rhs)

endif


return
end
