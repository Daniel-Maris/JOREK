!*******************************************************************************
!* Subroutine: boundary_condition                                              *
!*******************************************************************************
!*                                                                             *
!* Add boundary condition on the matrix.                                       *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   index_min    - Minimal index of local elements                            *
!*   index_max    - Maximal index of local elements                            *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*   gmres        - boolean indicating if we are using GMRES method            *
!*   solve_only   - Indicate if we want to perform only solve                  *
!*                                                                             *
!*******************************************************************************
!* Subroutine: solve_Psi_boundary_eqn                                          *
!*******************************************************************************
!*                                                                             *
!* Solve the differential equation for Psi on boundary, ensuring that n.B = 0  *
!* The solution is then used as an imhomogeneous Dirichlet b.c. for Psi        *
!*                                                                             *
!*******************************************************************************
!* Subroutine: setup_boundary_condition                                        *
!*******************************************************************************
!*                                                                             *
!* Projects the solution of the boundary Psi equation onto the JOREK finite    *
!*  element basis and writes the calculated dofs into the appropriate nodes.   *
!*                                                                             *
!*******************************************************************************
module mod_boundary_conditions
implicit none

real*8, dimension(:), allocatable :: Psi_b, zj_b

private :: Zn, Zn_chi, Zm, Zm_tht, Zm_tht_2, HZn

contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,          &
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2,       & 
                                  R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,  &  
                                  gmres, solve_only, ijA_index, ijA_size, irn_jcn, irn, jcn,          & 
                                  A_mat, i_tor_min, i_tor_max )

    use data_structure
    use phys_module, only: F0, GAMMA, keep_n0_const
    use vacuum, only: is_freebound
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_integer_types

    implicit none

    ! --- Routine parameters
    integer,                            intent(in)    :: my_id
    type (type_node_list),              intent(in)    :: node_list
    type (type_element_list),           intent(in)    :: element_list
    type (type_bnd_node_list),          intent(in)    :: bnd_node_list
    integer,                            intent(in)    :: local_elms(*)
    integer,                            intent(in)    :: n_local_elms
    integer,                            intent(in)    :: index_min
    integer,                            intent(in)    :: index_max
    real*8,                             intent(inout) :: rhs_loc(*)
    logical,                            intent(in)    :: xpoint2
    integer,                            intent(in)    :: xcase2
    real*8,                             intent(in)    :: R_axis
    real*8,                             intent(in)    :: Z_axis
    real*8,                             intent(in)    :: psi_axis
    real*8,                             intent(in)    :: psi_bnd
    real*8,                             intent(in)    :: R_xpoint(2)
    real*8,                             intent(in)    :: Z_xpoint(2)
    real*8,                             intent(in)    :: psi_xpoint(2)
    logical,                            intent(in)    :: gmres
    logical,                            intent(in)    :: solve_only
    integer,                            intent(in)    :: i_tor_min, i_tor_max 
    integer(kind=int_all), allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer(kind=int_all), allocatable, intent(inout) :: irn(:), jcn(:) 
    real*8,                allocatable, intent(inout) :: A_mat(:) 

    ! Internal parameters
    real*8                :: zbig, zbig_backup
    integer               :: i, in, iv, inode, k
    integer               :: ielm
    integer               :: index_node
    integer(kind=int_all) :: ijA_position
    integer               :: ilarge2
    integer               :: ierr, n_tor_local

    n_tor_local = i_tor_max - i_tor_min + 1
    zbig = 1.d12
    zbig_backup = zbig
       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)

             if (node_list%node(inode)%boundary .ne. 0) then

                do in=i_tor_min, i_tor_max 
                  if (keep_n0_const  .and.  in .eq. 1 ) then
                    zbig = 1.d15
                  else
                    zbig = zbig_backup
                  endif

                   do k=1, n_var

                      !------------------------------------ the open field lines (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                              (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then
 
                          if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)
                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local  & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)   = zbig

                            endif

                            index_node = node_list%node(inode)%index(2)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local   & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)    = zbig

                            endif

                          endif
                        endif
                      endif

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local   & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)   = zbig

                            endif

                            index_node = node_list%node(inode)%index(3)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local   & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)   = zbig
                            end if

                         endif

                      endif

                   enddo

                enddo
             endif
          enddo
       enddo

    return
  end subroutine boundary_conditions
  
  subroutine solve_Psi_boundary_eqn(node_list, element_list, boundary_list)
    use constants, only: pi
    use mod_parameters, only: n_period, n_plane, n_order, n_tor, n_coord_tor
    use phys_module, only: F0, m_pol_bc
    use gauss
    use basis_at_gaussian
    use data_structure
    use mod_chi
    use mod_SolveMN
    implicit none

    type(type_node_list),        intent(in) :: node_list
    type(type_element_list),     intent(in) :: element_list
    type(type_bnd_element_list), intent(in) :: boundary_list

    real*8  :: phi, theta, element_size_ij, element_size_perp, BigR, v, grad_chi(3), Jgrad_ps(3), xjac
    real*8  :: Psi_tht, Psi_tht_2, Bv2, grad_theta(3), grad_Bv2(3), grad_xjac(3), Lap_theta, alt_Lap_theta, Bv_dll_gt_gc, zj, chi_p
    real*8  :: gst, gtt, gtp, J2gst_s, J2gtt_t, J2gtp_p, J_s, J_t, J_p, alt_grad_theta(3)
    real*8  :: x_s_x, x_s_y, x_s_p, x_t_x, x_t_y, x_t_p, x_p_x, x_p_y, x_p_p
    real*8  :: y_s_x, y_s_y, y_s_p, y_t_x, y_t_y, y_t_p, y_p_x, y_p_y, y_p_p
    integer :: N_tht, mp, ms, ielm, i, i2, i3, j, j2, j3, n, m, nn, mm, in, in2, iv, iv2, ind1, ind2, info, direction_perp(2)

    type(type_bnd_element) :: bnd_element
    type(type_element)     :: element
    type(type_node)        :: node, node2

    real*8, dimension(:,:,:), allocatable :: x_g, x_s, x_t, x_p, x_ss, x_tt, x_st, x_sp, x_tp, x_pp
    real*8, dimension(:,:,:), allocatable :: y_g, y_s, y_t, y_p, y_ss, y_tt, y_st, y_sp, y_tp, y_pp
    real*8, dimension(:,:),   allocatable :: Amat
    real*8, dimension(:),     allocatable :: RHS

    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

    write(*,*) "***************************************"
    write(*,*) "*    Solving boundary Psi equation    *"
    write(*,*) "***************************************"
    write(*,'(A,I4,A,I4)') "Toroidal modes: ", n_tor
    write(*,'(A,I4,A,I4)') "Poloidal modes: ", m_pol_bc

    N_tht = boundary_list%n_bnd_elements

    allocate(x_g(n_plane,n_gauss,N_tht)); allocate(x_s, x_t, x_p, x_ss, x_tt, x_st, x_sp, x_tp, x_pp, mold=x_g)
    allocate(y_g(n_plane,n_gauss,N_tht)); allocate(y_s, y_t, y_p, y_ss, y_tt, y_st, y_sp, y_tp, y_pp, mold=y_g)
    allocate(Amat(n_tor*m_pol_bc,n_tor*m_pol_bc))
    allocate(RHS(n_tor*m_pol_bc))

    x_g  = 0.d0; x_s = 0.d0; x_t = 0.d0; x_p = 0.d0; x_ss = 0.d0; x_tt = 0.d0; x_st = 0.d0; x_sp = 0.d0; x_tp = 0.d0; x_pp = 0.d0
    y_g  = 0.d0; y_s = 0.d0; y_t = 0.d0; y_p = 0.d0; y_ss = 0.d0; y_tt = 0.d0; y_st = 0.d0; y_sp = 0.d0; y_tp = 0.d0; y_pp = 0.d0

    do ielm=1,N_tht
      bnd_element = boundary_list%bnd_element(ielm)
      element = element_list%element(bnd_element%element)

      do i=1,2
        node = node_list%node(bnd_element%vertex(i))
        direction_perp = (/ 6/bnd_element%direction(i,2), 4 /)
        i2 = findloc(element%vertex, bnd_element%vertex(i), 1)
        if (i .eq. 1) then
          if (bnd_element%side .le. 2) then
            i3 = mod(mod(mod(i2,4)+1,4)+1,4)+1
          else
            i3 = mod(i2,4)+1
          end if
        else
          if (bnd_element%side .le. 2) then
            i3 = mod(i2,4)+1
          else
            i3 = mod(mod(mod(i2,4)+1,4)+1,4)+1
          end if
        end if
        node2 = node_list%node(element%vertex(i3))

        do j=1,2
          j2 = bnd_element%direction(i,j)
          element_size_ij = bnd_element%size(i,j)
          j3 = direction_perp(j)
          if (bnd_element%side .eq. 1) then
            element_size_perp =  3.d0*element%size(i2,direction_perp(1))
          else
            element_size_perp = -3.d0*element%size(i2,direction_perp(1))
          end if

          do ms=1,n_gauss
            do mp=1,n_plane
              do in=1,n_coord_tor
                x_g(mp,ms,ielm)  = x_g(mp,ms,ielm)  + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)
                x_s(mp,ms,ielm)  = x_s(mp,ms,ielm)  + node%x(in,j2,1)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)
                x_t(mp,ms,ielm)  = x_t(mp,ms,ielm)  + node%x(in,j3,1)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)*element_size_perp
                x_p(mp,ms,ielm)  = x_p(mp,ms,ielm)  + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)
                x_ss(mp,ms,ielm) = x_ss(mp,ms,ielm) + node%x(in,j2,1)*element_size_ij*H1_ss(i,j,ms)*HZ_coord(in,mp)
                x_sp(mp,ms,ielm) = x_sp(mp,ms,ielm) + node%x(in,j2,1)*element_size_ij*H1_s(i,j,ms) *HZ_coord_p(in,mp)
                x_st(mp,ms,ielm) = x_st(mp,ms,ielm) + node%x(in,j3,1)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)*element_size_perp
                x_tp(mp,ms,ielm) = x_tp(mp,ms,ielm) + node%x(in,j3,1)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)*element_size_perp
                x_pp(mp,ms,ielm) = x_pp(mp,ms,ielm) + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord_pp(in,mp)
                x_tt(mp,ms,ielm) = x_tt(mp,ms,ielm) - 6.d0*(node%x(in,j2,1)*element_size_ij + 2.d0*node%x(in,j3,1)*element%size(i2,j3) &
                                 - node2%x(in,j2,1)*element%size(i3,j2) - node2%x(in,j3,1)*element%size(i3,j3))*H1(i,j,ms)*HZ_coord(in,mp)

                y_g(mp,ms,ielm)  = y_g(mp,ms,ielm)  + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)
                y_s(mp,ms,ielm)  = y_s(mp,ms,ielm)  + node%x(in,j2,2)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)
                y_t(mp,ms,ielm)  = y_t(mp,ms,ielm)  + node%x(in,j3,2)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)*element_size_perp
                y_p(mp,ms,ielm)  = y_p(mp,ms,ielm)  + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)
                y_ss(mp,ms,ielm) = y_ss(mp,ms,ielm) + node%x(in,j2,2)*element_size_ij*H1_ss(i,j,ms)*HZ_coord(in,mp)
                y_sp(mp,ms,ielm) = y_sp(mp,ms,ielm) + node%x(in,j2,2)*element_size_ij*H1_s(i,j,ms) *HZ_coord_p(in,mp)
                y_st(mp,ms,ielm) = y_st(mp,ms,ielm) + node%x(in,j3,2)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)*element_size_perp
                y_tp(mp,ms,ielm) = y_tp(mp,ms,ielm) + node%x(in,j3,2)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)*element_size_perp
                y_pp(mp,ms,ielm) = y_pp(mp,ms,ielm) + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord_pp(in,mp)
                y_tt(mp,ms,ielm) = y_tt(mp,ms,ielm) - 6.d0*(node%x(in,j2,2)*element_size_ij + 2.d0*node%x(in,j3,2)*element%size(i2,j3) &
                                 - node2%x(in,j2,2)*element%size(i3,j2) - node2%x(in,j3,2)*element%size(i3,j3))*H1(i,j,ms)*HZ_coord(in,mp)
              end do
            end do
          end do
        end do
      end do
    end do

    Amat = 0.d0; RHS = 0.d0

    do ielm=1,N_tht
      do ms=1,n_gauss
        theta = 2.d0*pi*(float(ielm-1) + xgauss(ms))/float(N_tht)

        do mp=1,n_plane
          BigR = x_g(mp,ms,ielm)
          phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
          chi = get_chi(x_g(mp,ms,ielm),y_g(mp,ms,ielm),phi)
          grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/BigR /)
          Jgrad_ps = (/ -BigR*y_s(mp,ms,ielm), BigR*x_s(mp,ms,ielm), x_p(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_p(mp,ms,ielm) /)
          Jgrad_ps = Jgrad_ps*N_tht/(2.d0*pi)

          ind1 = 1
          do n=1,n_tor
            do m=1,m_pol_bc
              RHS(ind1) = RHS(ind1) + dot_product(Jgrad_ps,grad_chi)*Zn(n,chi(0,0,0))*Zm(m,theta)*wgauss(ms)
              ind1 = ind1 + 1
            end do
          end do
        end do
      end do
    end do
    RHS = RHS*(2.d0*pi)**2/float(n_period*n_plane*N_tht)

    ind1 = 1
    do n=1,n_tor
      do m=1,m_pol_bc
        ind2 = 1
        do nn=1,n_tor
          do mm=1,m_pol_bc
            if (n .eq. nn) then
              if (mod(m,2) .eq. 1 .and. mm .eq. m + 1) then
                Amat(ind1,ind2) =  pi*(mm/2)
              else if (mod(m,2) .eq. 0 .and. mm .eq. m - 1) then
                Amat(ind1,ind2) = -pi*(m/2)
              end if

              if (n .eq. 1) then
                Amat(ind1,ind2) = Amat(ind1,ind2)*2.d0*pi*F0/n_period
              else
                Amat(ind1,ind2) = Amat(ind1,ind2)*pi*F0/n_period
              end if
            end if

            ind2 = ind2 + 1
          end do
        end do

        ind1 = ind1 + 1
      end do
    end do

    call SolveMN(Amat, info)
    if (info .ne. 0) then
      write(*,'(A,I3)') "ERROR: solve_Psi_boundary_eqn, matrix inversion failed: ", info
      stop
    end if

    if (allocated(Psi_b)) deallocate(Psi_b)
    allocate(Psi_b(n_tor*m_pol_bc))

    Psi_b = matmul(Amat,RHS)

    RHS = 0.d0

    open(25,file="zj.txt",action="write",status="replace")
    do ielm=1,N_tht
      do ms=1,n_gauss
        theta = 2.d0*pi*(float(ielm-1) + xgauss(ms))/float(N_tht)

        do mp=1,n_plane
          BigR = x_g(mp,ms,ielm)
          phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
          chi = get_chi(x_g(mp,ms,ielm),y_g(mp,ms,ielm),phi)
          chi_p = chi(0,0,1) + x_p(mp,ms,ielm)*chi(1,0,0) + y_p(mp,ms,ielm)*chi(0,1,0)
          grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/BigR /)
          Bv2 = dot_product(grad_chi,grad_chi)
          xjac = x_t(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_t(mp,ms,ielm)
          grad_theta = (/ -y_t(mp,ms,ielm), x_t(mp,ms,ielm), (x_p(mp,ms,ielm)*y_t(mp,ms,ielm) - x_t(mp,ms,ielm)*y_p(mp,ms,ielm))/BigR /)
          grad_theta = grad_theta*2.d0*pi/(xjac*N_tht)
          grad_Bv2 = 2.d0*(/ chi(1,0,0)*chi(2,0,0) + chi(0,1,0)*chi(1,1,0) + chi(0,0,1)*chi(1,0,1)/BigR**2 - chi(0,0,1)**2/BigR**3, &
                             chi(1,0,0)*chi(1,1,0) + chi(0,1,0)*chi(0,2,0) + chi(0,0,1)*chi(0,1,1)/BigR**2, &
                            (chi(1,0,0)*chi(1,0,1) + chi(0,1,0)*chi(0,1,1) + chi(0,0,1)*chi(0,0,2)/BigR**2)/BigR /)

          x_s_x = (x_st(mp,ms,ielm)*y_s(mp,ms,ielm) - x_ss(mp,ms,ielm)*y_t(mp,ms,ielm))/xjac
          x_s_y = (x_ss(mp,ms,ielm)*x_t(mp,ms,ielm) - x_st(mp,ms,ielm)*x_s(mp,ms,ielm))/xjac
          x_s_p = x_sp(mp,ms,ielm) - x_p(mp,ms,ielm)*x_s_x - y_p(mp,ms,ielm)*x_s_y
          x_t_x = (x_tt(mp,ms,ielm)*y_s(mp,ms,ielm) - x_st(mp,ms,ielm)*y_t(mp,ms,ielm))/xjac
          x_t_y = (x_st(mp,ms,ielm)*x_t(mp,ms,ielm) - x_tt(mp,ms,ielm)*x_s(mp,ms,ielm))/xjac
          x_t_p = x_tp(mp,ms,ielm) - x_p(mp,ms,ielm)*x_t_x - y_p(mp,ms,ielm)*x_t_y
          x_p_x = (x_tp(mp,ms,ielm)*y_s(mp,ms,ielm) - x_sp(mp,ms,ielm)*y_t(mp,ms,ielm))/xjac
          x_p_y = (x_sp(mp,ms,ielm)*x_t(mp,ms,ielm) - x_tp(mp,ms,ielm)*x_s(mp,ms,ielm))/xjac
          x_p_p = x_pp(mp,ms,ielm) - x_p(mp,ms,ielm)*x_p_x - y_p(mp,ms,ielm)*x_p_y
          y_s_x = (y_st(mp,ms,ielm)*y_s(mp,ms,ielm) - y_ss(mp,ms,ielm)*y_t(mp,ms,ielm))/xjac
          y_s_y = (y_ss(mp,ms,ielm)*x_t(mp,ms,ielm) - y_st(mp,ms,ielm)*x_s(mp,ms,ielm))/xjac
          y_s_p = y_sp(mp,ms,ielm) - x_p(mp,ms,ielm)*y_s_x - y_p(mp,ms,ielm)*y_s_y
          y_t_x = (y_tt(mp,ms,ielm)*y_s(mp,ms,ielm) - y_st(mp,ms,ielm)*y_t(mp,ms,ielm))/xjac
          y_t_y = (y_st(mp,ms,ielm)*x_t(mp,ms,ielm) - y_tt(mp,ms,ielm)*x_s(mp,ms,ielm))/xjac
          y_t_p = y_tp(mp,ms,ielm) - x_p(mp,ms,ielm)*y_t_x - y_p(mp,ms,ielm)*y_t_y
          y_p_x = (y_tp(mp,ms,ielm)*y_s(mp,ms,ielm) - y_sp(mp,ms,ielm)*y_t(mp,ms,ielm))/xjac
          y_p_y = (y_sp(mp,ms,ielm)*x_t(mp,ms,ielm) - y_tp(mp,ms,ielm)*x_s(mp,ms,ielm))/xjac
          y_p_p = y_pp(mp,ms,ielm) - x_p(mp,ms,ielm)*y_p_x - y_p(mp,ms,ielm)*y_p_y

          grad_xjac = (/ x_t_x*y_s(mp,ms,ielm) + x_t(mp,ms,ielm)*y_s_x - x_s_x*y_t(mp,ms,ielm) - x_s(mp,ms,ielm)*y_t_x, &
                         x_t_y*y_s(mp,ms,ielm) + x_t(mp,ms,ielm)*y_s_y - x_s_y*y_t(mp,ms,ielm) - x_s(mp,ms,ielm)*y_t_y, &
                        (y_s(mp,ms,ielm)*x_t_p + x_t(mp,ms,ielm)*y_s_p - y_t(mp,ms,ielm)*x_s_p - x_s(mp,ms,ielm)*y_t_p)/BigR /)
          Lap_theta = -dot_product(grad_theta,grad_xjac)/xjac - 2.d0*pi*(y_t_x + y_t(mp,ms,ielm)/BigR - x_t_y &
                    + (x_t_p*y_p(mp,ms,ielm) + x_t(mp,ms,ielm)*y_p_p - x_p_p*y_t(mp,ms,ielm) - x_p(mp,ms,ielm)*y_t_p)/BigR**2)/(xjac*N_tht)
          Bv_dll_gt_gc = dot_product(grad_Bv2,grad_theta)/2.d0 + grad_theta(1)*grad_chi(3)**2/BigR &
              - 2.d0*pi*(chi(1,0,0)*(chi(1,0,0)*y_t_x + chi(0,1,0)*y_t_y + chi(0,0,1)*y_t_p/BigR**2) &
                       - chi(0,1,0)*(chi(1,0,0)*x_t_x + chi(0,1,0)*x_t_y + chi(0,0,1)*x_t_p/BigR**2) &
                       + chi(0,0,1)*(chi(1,0,0)*x_t_x + chi(0,1,0)*x_t_y + chi(0,0,1)*x_t_p/BigR**2)*y_p(mp,ms,ielm)/BigR**2 &
                       + chi(0,0,1)*(chi(1,0,0)*y_p_x + chi(0,1,0)*y_p_y + chi(0,0,1)*y_p_p/BigR**2)*x_t(mp,ms,ielm)/BigR**2 &
                       - chi(0,0,1)*(chi(1,0,0)*x_p_x + chi(0,1,0)*x_p_y + chi(0,0,1)*x_p_p/BigR**2)*y_t(mp,ms,ielm)/BigR**2 &
                       - chi(0,0,1)*(chi(1,0,0)*y_t_x + chi(0,1,0)*y_t_y + chi(0,0,1)*y_t_p/BigR**2)*x_p(mp,ms,ielm)/BigR**2)/(xjac*N_tht) &
                       - dot_product(grad_theta,grad_chi)*dot_product(grad_chi,grad_xjac)/xjac - 2.d0*grad_chi(1)*grad_chi(3)*grad_theta(3)/BigR
          J2gst_s = ((-BigR**2 - x_p(mp,ms,ielm)**2)*y_t(mp,ms,ielm) + x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_p(mp,ms,ielm))*y_st(mp,ms,ielm) &
                  + ((-2*BigR*x_t(mp,ms,ielm) - 2*x_p(mp,ms,ielm)*x_tp(mp,ms,ielm))*y_t(mp,ms,ielm) + (-BigR**2 - x_p(mp,ms,ielm)**2)*y_tt(mp,ms,ielm) &
                  + x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_tp(mp,ms,ielm) + x_p(mp,ms,ielm)*x_tt(mp,ms,ielm)*y_p(mp,ms,ielm) &
                  + x_t(mp,ms,ielm)*y_p(mp,ms,ielm)*x_tp(mp,ms,ielm))*y_s(mp,ms,ielm) - BigR**2*x_t(mp,ms,ielm)*x_st(mp,ms,ielm) &
                  - BigR**2*x_tt(mp,ms,ielm)*x_s(mp,ms,ielm) - 2*BigR*x_t(mp,ms,ielm)**2*x_s(mp,ms,ielm) &
                  + x_p(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_tt(mp,ms,ielm) + x_p(mp,ms,ielm)*x_s(mp,ms,ielm)*y_t(mp,ms,ielm)*y_tp(mp,ms,ielm) &
                  + x_p(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_st(mp,ms,ielm) - 2*x_t(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_tp(mp,ms,ielm) &
                  - x_t(mp,ms,ielm)*y_p(mp,ms,ielm)**2*x_st(mp,ms,ielm) - x_tt(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)**2 &
                  + x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_tp(mp,ms,ielm)
          J2gtt_t = (2*BigR*x_s(mp,ms,ielm) + 2*x_p(mp,ms,ielm)*x_sp(mp,ms,ielm))*y_t(mp,ms,ielm)**2 &
                  + 2*(BigR**2 + x_p(mp,ms,ielm)**2)*y_t(mp,ms,ielm)*y_st(mp,ms,ielm) + 2*BigR**2*x_t(mp,ms,ielm)*x_st(mp,ms,ielm) &
                  + 2*BigR*x_t(mp,ms,ielm)**2*x_s(mp,ms,ielm) - 2*x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_p(mp,ms,ielm)*y_st(mp,ms,ielm) &
                  - 2*x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_t(mp,ms,ielm)*y_sp(mp,ms,ielm) - 2*x_p(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_st(mp,ms,ielm) &
                  + 2*x_t(mp,ms,ielm)**2*y_p(mp,ms,ielm)*y_sp(mp,ms,ielm) + 2*x_t(mp,ms,ielm)*y_p(mp,ms,ielm)**2*x_st(mp,ms,ielm) &
                  - 2*x_t(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_sp(mp,ms,ielm)
          J2gtp_p = (x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_t(mp,ms,ielm) - x_t(mp,ms,ielm)**2*y_p(mp,ms,ielm))*y_sp(mp,ms,ielm) &
                  + (x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_tp(mp,ms,ielm) + x_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_tp(mp,ms,ielm) &
                  + x_pp(mp,ms,ielm)*x_t(mp,ms,ielm)*y_t(mp,ms,ielm) - x_t(mp,ms,ielm)**2*y_pp(mp,ms,ielm) &
                  - 2*x_t(mp,ms,ielm)*y_p(mp,ms,ielm)*x_tp(mp,ms,ielm))*y_s(mp,ms,ielm) - 2*x_p(mp,ms,ielm)*x_s(mp,ms,ielm)*y_t(mp,ms,ielm)*y_tp(mp,ms,ielm) &
                  - x_p(mp,ms,ielm)*y_t(mp,ms,ielm)**2*x_sp(mp,ms,ielm) - x_pp(mp,ms,ielm)*x_s(mp,ms,ielm)*y_t(mp,ms,ielm)**2 &
                  + x_t(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_tp(mp,ms,ielm) + x_t(mp,ms,ielm)*x_s(mp,ms,ielm)*y_pp(mp,ms,ielm)*y_t(mp,ms,ielm) &
                  + x_t(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_sp(mp,ms,ielm) + x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm)*x_tp(mp,ms,ielm)
          gst = (((-BigR**2 - x_p(mp,ms,ielm)**2)*y_t(mp,ms,ielm) + x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_p(mp,ms,ielm))*y_s(mp,ms,ielm) &
              - BigR**2*x_t(mp,ms,ielm)*x_s(mp,ms,ielm) + x_p(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm) &
              - x_t(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)**2)/(BigR**2*xjac**2)
          gtt = ((BigR**2 + x_p(mp,ms,ielm)**2)*y_t(mp,ms,ielm)**2 + BigR**2*x_t(mp,ms,ielm)**2 &
              - 2*x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm) + x_t(mp,ms,ielm)**2*y_p(mp,ms,ielm)**2)/(BigR**2*xjac**2)
          gtp = ((x_p(mp,ms,ielm)*x_t(mp,ms,ielm)*y_t(mp,ms,ielm) - x_t(mp,ms,ielm)**2*y_p(mp,ms,ielm))*y_s(mp,ms,ielm) &
              - x_p(mp,ms,ielm)*x_s(mp,ms,ielm)*y_t(mp,ms,ielm)**2 + x_t(mp,ms,ielm)*x_s(mp,ms,ielm)*y_p(mp,ms,ielm)*y_t(mp,ms,ielm))/(BigR**2*xjac**2)
          J_s = (x_t(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_t(mp,ms,ielm))*x_t(mp,ms,ielm) + (x_t(mp,ms,ielm)*y_st(mp,ms,ielm) &
              + x_tt(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_tt(mp,ms,ielm) - y_t(mp,ms,ielm)*x_st(mp,ms,ielm))*BigR
          J_t = (x_t(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_t(mp,ms,ielm))*x_s(mp,ms,ielm) + (x_t(mp,ms,ielm)*y_ss(mp,ms,ielm) &
              - x_s(mp,ms,ielm)*y_st(mp,ms,ielm) - x_ss(mp,ms,ielm)*y_t(mp,ms,ielm) + y_s(mp,ms,ielm)*x_st(mp,ms,ielm))*BigR
          J_p = (x_t(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_t(mp,ms,ielm))*x_p(mp,ms,ielm) + (x_t(mp,ms,ielm)*y_sp(mp,ms,ielm) &
              - x_s(mp,ms,ielm)*y_tp(mp,ms,ielm) - y_t(mp,ms,ielm)*x_sp(mp,ms,ielm) + y_s(mp,ms,ielm)*x_tp(mp,ms,ielm))*BigR
          alt_Lap_theta = 2.d0*pi*((J2gst_s + J2gtt_t + J2gtp_p)/(BigR**2*xjac**2) - (gst*J_s + gtt*J_t + gtp*J_p)/(BigR*xjac))/N_tht
          alt_grad_theta = 2.d0*pi*(gst*(/ x_t(mp,ms,ielm), y_t(mp,ms,ielm), 0.d0 /) + gtt*(/ x_s(mp,ms,ielm), y_s(mp,ms,ielm), 0.d0 /) &
                                  + gtp*(/ x_p(mp,ms,ielm), y_p(mp,ms,ielm), BigR /))/N_tht

          Psi_tht = 0.d0
          Psi_tht_2 = 0.d0

          ind1 = 1
          do n=1,n_tor
            do m=1,m_pol_bc
              Psi_tht = Psi_tht + Psi_b(ind1)*Zn(n,chi(0,0,0))*Zm_tht(m,theta)
              Psi_tht_2 = Psi_tht_2 + Psi_b(ind1)*Zn(n,chi(0,0,0))*Zm_tht_2(m,theta)
              ind1 = ind1 + 1
            end do
          end do
          
          zj = Psi_tht_2*(dot_product(grad_theta,grad_theta) - dot_product(grad_theta,grad_chi)**2/Bv2) &
             + Psi_tht*(dot_product(grad_Bv2,grad_theta)/Bv2 + Lap_theta - Bv_dll_gt_gc/Bv2)
          write(25,'(4E14.6)') theta, phi, zj, (Lap_theta - alt_Lap_theta)
          ind1 = 1
          do n=1,n_tor
            do m=1,m_pol_bc
              RHS(ind1) = RHS(ind1) + chi_p*zj*Zn(n,chi(0,0,0))*Zm(m,theta)*wgauss(ms)
              ind1 = ind1 + 1
            end do
          end do
        end do
      end do
    end do
    close(25)
    RHS = RHS*(2.d0*pi)**2/float(n_period*n_plane*N_tht)

    if (allocated(zj_b)) deallocate(zj_b)
    allocate(zj_b(n_tor*m_pol_bc))

    ind1 = 1
    do n=1,n_tor
      do m=1,m_pol_bc
        zj_b(ind1) = RHS(ind1)*n_period/(F0*pi**2)
        if (n .eq. 1) zj_b(ind1) = zj_b(ind1)/2.d0
        ind1 = ind1 + 1
      end do
    end do
  end subroutine solve_Psi_boundary_eqn

  subroutine setup_boundary_condition(node_list, bnd_node_list)
    use constants, only: pi
    use mod_parameters, only: n_tor, n_order, n_plane, n_period, n_coord_tor, var_Psi
    use phys_module, only: F0, m_pol_bc
    use basis_at_gaussian
    use data_structure
    use mod_chi
    implicit none
    
    type(type_node_list),     intent(inout) :: node_list
    type(type_bnd_node_list), intent(in)    :: bnd_node_list
    
    integer :: ielm, i, j, mp, im, in, n, m, ind
    real*8  :: R, z, R_s, z_s, phi, theta, Psi, Psi_tht, Psi_chi, zj, zj_tht, zj_chi, chi_tht
    
    real*8, dimension(n_tor) :: Psi_dof, Psi_dof2, zj_dof, zj_dof2
    
    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi
    
    do i=1,bnd_node_list%n_bnd_nodes
      in = bnd_node_list%bnd_node(i)%index_jorek
      j  = bnd_node_list%bnd_node(i)%direction(2)
      theta = 2.d0*pi*float(i-1)/float(bnd_node_list%n_bnd_nodes)
      Psi_dof = 0.d0; Psi_dof2 = 0.d0
      zj_dof = 0.d0; zj_dof2 = 0.d0
      
      do mp=1,n_plane
        phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)

        R = 0.d0; R_s = 0.d0
        z = 0.d0; z_s = 0.d0
        Psi = 0.d0; Psi_tht = 0.d0; Psi_chi = 0.d0
        zj = 0.d0; zj_tht = 0.d0; zj_chi = 0.d0
          
        do im=1,n_coord_tor
          R   = R   + node_list%node(in)%x(im,1,1)*HZ_coord(im,mp)
          R_s = R_s + node_list%node(in)%x(im,j,1)*HZ_coord(im,mp)*3.d0
          z   = z   + node_list%node(in)%x(im,1,2)*HZ_coord(im,mp)
          z_s = z_s + node_list%node(in)%x(im,j,2)*HZ_coord(im,mp)*3.d0
        end do
        
        chi = get_chi(R,z,phi)
        chi_tht = (R_s*chi(1,0,0) + z_s*chi(0,1,0))*bnd_node_list%n_bnd_nodes/(2.d0*pi)
        
        ind = 1
        do n=1,n_tor
          do m=1,m_pol_bc
            Psi = Psi + Psi_b(ind)*Zn(n,chi(0,0,0))*Zm(m,theta)
            Psi_tht = Psi_tht + Psi_b(ind)*Zn(n,chi(0,0,0))*Zm_tht(m,theta)
            Psi_chi = Psi_chi + Psi_b(ind)*Zn_chi(n,chi(0,0,0))*Zm(m,theta)
            zj = zj + zj_b(ind)*Zn(n,chi(0,0,0))*Zm(m,theta)
            zj_tht = zj_tht + zj_b(ind)*Zn(n,chi(0,0,0))*Zm_tht(m,theta)
            zj_chi = zj_chi + zj_b(ind)*Zn_chi(n,chi(0,0,0))*Zm(m,theta)
            ind = ind + 1
          end do
        end do
        
        do im=1,n_tor
          Psi_dof(im) = Psi_dof(im) + Psi*HZn(im,phi)
          Psi_dof2(im) = Psi_dof2(im) + (Psi_tht + Psi_chi*chi_tht)*HZn(im,phi)
          zj_dof(im) = zj_dof(im) + zj*HZn(im,phi)
          zj_dof2(im) = zj_dof2(im) + (zj_tht + zj_chi*chi_tht)*HZn(im,phi)
        end do
      end do
      
      Psi_dof2 = Psi_dof2*2.d0*pi/bnd_node_list%n_bnd_nodes
      zj_dof2 = zj_dof2*2.d0*pi/bnd_node_list%n_bnd_nodes
      
      Psi_dof(1) = Psi_dof(1)/n_plane
      Psi_dof2(1) = Psi_dof2(1)/(3.d0*n_plane)
      zj_dof(1) = zj_dof(1)/n_plane
      zj_dof2(1) = zj_dof2(1)/(3.d0*n_plane)
      if (n_tor .gt. 1) then
        Psi_dof(2:) = Psi_dof(2:)*2.d0/n_plane
        Psi_dof2(2:) = Psi_dof2(2:)*2.d0/(3.d0*n_plane)
        zj_dof(2:) = zj_dof(2:)*2.d0/n_plane
        zj_dof2(2:) = zj_dof2(2:)*2.d0/(3.d0*n_plane)
      end if
      
      node_list%node(in)%values(:,1,var_Psi) = F0*Psi_dof
      node_list%node(in)%values(:,j,var_Psi) = F0*Psi_dof2
      node_list%node(in)%values(:,1,var_zj) = F0*zj_dof
      node_list%node(in)%values(:,j,var_zj) = F0*zj_dof2
    end do
  end subroutine setup_boundary_condition

  pure real*8 function Zn(n, chi)
    use mod_parameters, only: n_period
    use phys_module, only: F0
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: chi

    if (n .eq. 1) then
      Zn = 1.d0
      return
    end if

    if (mod(n,2) .eq. 0) then
      Zn = cos(n_period*int(n/2)*chi/F0)
    else
      Zn = sin(n_period*int(n/2)*chi/F0)
    end if
  end function Zn

  pure real*8 function Zn_chi(n, chi)
    use mod_parameters, only: n_period
    use phys_module, only: F0
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: chi

    if (n .eq. 1) then
      Zn_chi = 0.d0
      return
    end if

    if (mod(n,2) .eq. 0) then
      Zn_chi = -n_period*int(n/2)*sin(n_period*int(n/2)*chi/F0)/F0
    else
      Zn_chi =  n_period*int(n/2)*cos(n_period*int(n/2)*chi/F0)/F0
    end if
  end function Zn_chi

  pure real*8 function Zm(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta

    if (mod(m,2) .eq. 1) then
      Zm = cos(int((m+1)/2)*theta)
    else
      Zm = sin(int((m+1)/2)*theta)
    end if
  end function Zm

  pure real*8 function Zm_tht(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta

    if (mod(m,2) .eq. 1) then
      Zm_tht = -int((m+1)/2)*sin(int((m+1)/2)*theta)
    else
      Zm_tht =  int((m+1)/2)*cos(int((m+1)/2)*theta)
    end if
  end function Zm_tht

  pure real*8 function Zm_tht_2(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta

    if (mod(m,2) .eq. 1) then
      Zm_tht_2 = -int((m+1)/2)**2*cos(int((m+1)/2)*theta)
    else
      Zm_tht_2 = -int((m+1)/2)**2*sin(int((m+1)/2)*theta)
    end if
  end function Zm_tht_2

  pure real*8 function HZn(n, phi)
    use mod_parameters, only: n_period
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: phi

    if (n .eq. 1) then
      HZn = 1.d0
      return
    end if

    if (mod(n,2) .eq. 0) then
      HZn = cos(n_period*int(n/2)*phi)
    else
      HZn = sin(n_period*int(n/2)*phi)
    end if
  end function HZn
end module mod_boundary_conditions
