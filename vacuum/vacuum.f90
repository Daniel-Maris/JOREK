subroutine vacuum(my_id,node_list,element_list,boundary_list,index_min,index_max, &
                  xpoint2,psi_axis,psi_bnd,Z_xpoint)	
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundary integral of the 
! induction equation using the vacuum response from STARWALL
!---------------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use global_distributed_matrix

implicit none

type (type_node_list)     :: node_list
type (type_element_list)  :: element_list
type (type_boundary_list) :: boundary_list

real*8     :: MK(2,2,2,2,boundary_list%n_boundary)

real*8     :: x_g(n_gauss), x_s(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

integer    :: index_min, index_max
integer    :: my_id, ibnd,i, j, ms, mp, kp, kbnd, k, l, jdir, ldir, index_node, index_node2, ilarge_pp, ijA_position
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: ws, xjac,  BigR, PI, phi, eps_cyl
real*8     :: psi_axis, psi_bnd, Z_xpoint
real*8     :: rhs_glob_1, A_glob_11
real*8     :: psi_norm, theta, zeta, gamma_sheeth

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, r0, T0  
real*8     :: psi, psi_s, rho,  T, eta_T  
real*8     :: amat_61, amat_65, amat_66, amat_67
logical    :: xpoint2


PI    = 2.d0*asin(1.d0)

theta = 0.5d0; zeta = 0.d0          ! Crank-Nicholson parameter
!theta = 1.0d0  ; zeta = 0.0d0       ! Euler scheme 
!theta = 1.0d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0;  x_ss  = 0.d0; 
y_g  = 0.d0; y_s  = 0.d0;  y_ss  = 0.d0; 
eq_g = 0.d0; eq_s = 0.d0;  eq_ss = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; 

do ibnd = 1, boundary_list%n_boundary                     ! loop over all boundary elements

  do i=1,2                                                ! loop over two corners of each boundary element
  
    do j=1,2                                              ! loop over the two basis functions (H, H_s)

      do ms=1, n_gauss                                    ! loop over Gaussian points, construct coordinates and values

        x_g(ms)  = x_g(ms)  + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%x(j,1) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1(i,j,ms)

        x_s(ms)  = x_s(ms)  + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%x(j,1) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1_s(i,j,ms)

        y_g(ms)  = y_g(ms)  + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%x(j,2) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1(i,j,ms)

        y_s(ms)  = y_s(ms)  + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%x(j,2) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1_s(i,j,ms)

        do mp=1,n_plane                                   ! loop over toroidal planes

          do k=1,n_var                                    ! loop over all variables

            do in=1,n_tor                                 ! loop over toroidal harmonics

              eq_g(mp,k,ms)  = eq_g(mp,k,ms) + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%values(in,j,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1(i,j,ms)   * HZ(in,mp)

              eq_s(mp,k,ms)  = eq_s(mp,k,ms) + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%values(in,j,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1_s(i,j,ms) * HZ(in,mp)

              eq_p(mp,k,ms)  = eq_p(mp,k,ms) + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%values(in,j,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1(i,j,ms)   * HZ_p(in,mp)

              delta_g(mp,k,ms) = delta_g(mp,k,ms)  + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%deltas(in,j,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1(i,j,ms)   * HZ(in,mp)
         
	      delta_s(mp,k,ms) = delta_s(mp,k,ms)  + node_list%node(boundary_list%boundary(ibnd)%vertex(i))%deltas(in,j,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1_s(i,j,ms)   * HZ(in,mp)

            enddo                                         ! end loop over toroidal harmonics
          enddo                                           ! end loop over all variables
        enddo                                             ! end loop over toroidal planes

      enddo                                               ! end loop over Gaussian points
    enddo                                                 ! end loop over the two basis function
  enddo                                                   ! end loop over two corners of each boundary element
 

!--------------------------------------------------- sum over the Gaussian integration points
  do ms=1, n_gauss                                        ! loop over Gaussian points

    ws = wgauss(ms)
!    xjac  =  x_s(ms)*y_t(ms) - x_t(ms)*y_s(ms)
    
    do mp = 1, n_plane                                    ! sum over toroidal planes (direct integration in toroidal direction)

      ps0   = eq_g(mp,1,ms)
      ps0_s = eq_s(mp,1,ms)

      r0    = eq_g(mp,5,ms)
      T0    = eq_g(mp,6,ms)
     
      BigR = x_g(ms)
      
!      grad_s2 = (x_t(ms)**2 + y_t(ms)**2) / xjac**2
     
      eta_T   = eta   * (abs(T0)/T_0)**(-1.5d0)                                 ! temperature dependent resistivity

      do i=1,2                                                                  ! loop over nodes of this piece of boundary
     
        do j=1,2                                                                ! loop over basis functions

          jdir = boundary_list%boundary(ibnd)%direction(i,j)
          
          index_node = node_list%node(boundary_list%boundary(ibnd)%vertex(i))%index(jdir)

          if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then  
  
            do im=1,n_tor                                                       ! loop over toroidal harmonics
 
              v   =  H1(i,j,ms) * boundary_list%boundary(ibnd)%size(i,j) * HZ(im,mp)         ! test function



              do kbnd = 1, boundary_list%n_boundary                     ! loop over all boundary elements


                do k=1,2                                                          ! loop over nodes

                  do l=1,2                                                        ! loop over basis functions

                    ldir = boundary_list%boundary(kbnd)%direction(k,l)
          
                    index_node2 = node_list%node(boundary_list%boundary(kbnd)%vertex(k))%index(ldir)
		
                    call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position)
	
                    kp = 1

                    do in = 1, n_tor                                                       ! loop over toroidal harmonics
  
                      ilarge_pp  = ijA_position  - 1 + ((kp-1)*n_tor + im-1) * n_var*n_tor + (kp-1)*n_tor + in  ! to be verified (in<->im)

                      psi   =  H1(k,l,ms)   * boundary_list%boundary(kbnd)%size(k,l) * HZ(in,mp)         ! test function
                      psi_s =  H1_s(k,l,ms) * boundary_list%boundary(kbnd)%size(k,l) * HZ(in,mp)         ! test function derivative

                      rho = psi
		      T   = psi   
		    
		      A_glob_11 = - v * eta_T * MK(k,l,i,j,kbnd) * psi_s / BigR * xjac * tstep ! (how to avoid xjac)
 
                      irn_glob(ilarge_pp) =  n_tor * n_var * (index_node-1)  + (kp-1)*n_tor + im
                      jcn_glob(ilarge_pp) =  n_tor * n_var * (index_node2-1) + (kp-1)*n_tor + in
                    
		      A_glob(ilarge_pp)   =  A_glob(ilarge_pp) + ws * A_glob_11

                    enddo   ! end loop over toroidal harmonics
                  enddo     ! end loop over basis functions
                enddo       ! end loop over nodes

              enddo         ! end of loop over all boundary elements

            enddo         ! end loop over toroidal harmonics  
	  
	  endif           ! endif selection of local indices
	  
        enddo             ! end of loop over basis functions        
	
      enddo               ! end of loop over nodes     
    enddo                 ! end of loop over toroidal planes
  enddo                   ! end of loop over Gaussian points

enddo                     ! end of loop over all boundary elements

return
end
