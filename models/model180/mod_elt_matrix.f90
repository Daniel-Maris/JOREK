module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
  !---------------------------------------------------------------
  ! calculates the matrix contribution of one element
  !---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use tr_module
use diffusivities, only: get_dperp, get_zkperp    
use corr_neg
use vacuum, only: freeb_fact

implicit none

type (type_element), intent(in)   :: element
type (type_node)   , intent(in)   :: nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in) :: tid

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, kl1, kl2, kl3, kl4, kl5, kl6
real*8     :: wst,  xjac, xjac_x, xjac_y, xjac_s, xjac_t, BigR, r2, phi, eps_cyl
real*8     :: Bv, Bv_x, Bv_y, Bv_p
real*8     :: chi_x, chi_y, chi_p
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: D_prof, ZK_prof, psi_norm
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6
real*8     :: theta, zeta, delta_u_x, delta_u_y

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xs, v_ys, v_xt, v_yt, v_xy
real*8     :: ps0, ps0_x, ps0_y, ps0_p, ps0_s, ps0_t, ps0_ss, ps0_st, ps0_tt, ps0_xx, ps0_xy, ps0_yy, ps0_pp, ps0_xp, ps0_yp
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t,  r0_hat, r0_x_hat, r0_y_hat, T0, T0_x, T0_y, T0_p, T0_s, T0_t
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_yy, psi_pp, psi_xy, psi_xp, psi_yp
real*8     :: u, u_x, u_y, u_p, u_s, u_t
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, T, T_x, T_y, T_s, T_t, T_p
real*8     :: P0, P0_x, P0_y, P0_s, P0_t
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_32, amat_44, amat_42
real*8     :: amat_51, amat_52, amat_55, amat_61, amat_62, amat_66
real*8     :: amat_stab_11, amat_stab_12, amat_stab_13, amat_stab_14 ,amat_stab_21,amat_stab_22, amat_stab_23, amat_stab_24
real*8     :: amat_stab_31, amat_stab_32, amat_stab_33, amat_stab_34 ,amat_stab_41,amat_stab_42, amat_stab_43, amat_stab_44

logical    :: xpoint2

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_st, eq_ss, eq_tt
real*8, dimension(:,:,:,:) , pointer :: eq_p, eq_pp, eq_sp, eq_tp
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t

eq_g    => thread_struct(tid)%eq_g   
eq_s    => thread_struct(tid)%eq_s   
eq_t    => thread_struct(tid)%eq_t   
eq_p    => thread_struct(tid)%eq_p
eq_pp   => thread_struct(tid)%eq_pp   
eq_sp   => thread_struct(tid)%eq_sp
eq_tp   => thread_struct(tid)%eq_tp
eq_ss   => thread_struct(tid)%eq_ss   
eq_st   => thread_struct(tid)%eq_st   
eq_tt   => thread_struct(tid)%eq_tt   
delta_g => thread_struct(tid)%delta_g
delta_s => thread_struct(tid)%delta_s
delta_t => thread_struct(tid)%delta_t


ELM = 0.d0
RHS = 0.d0

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
zeta  = time_evol_zeta

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s   = 0.d0; x_t   = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s   = 0.d0; y_t   = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s  = 0.d0; eq_t  = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0;
eq_p = 0.d0; eq_pp = 0.d0; eq_sp = 0.d0; eq_tp = 0.d0

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source     = 0.d0

do i=1,n_vertex_max
 do j=1,n_order+1

   do ms=1, n_gauss
     do mt=1, n_gauss

       x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
       x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
       x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)
       x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
       x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_st(i,j,ms,mt)
       x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

       y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
       y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
       y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)
       y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
       y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_st(i,j,ms,mt)
       y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

       do mp=1,n_plane

         do k=1,n_var

           do in=1,n_tor

             eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)

             eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)

             eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)

             eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
             
             eq_pp(mp,k,ms,mt) = eq_pp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ_pp(in,mp)
             eq_sp(mp,k,ms,mt) = eq_sp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ_p(in,mp)
             eq_tp(mp,k,ms,mt) = eq_tp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ_p(in,mp)
             
             eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
             eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
             eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)

             delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
             delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
             delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)

           enddo

         enddo

       enddo

       if (keep_current_prof) &
         call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
       call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))

     enddo
   enddo
 enddo
enddo

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac    = x_s(ms,mt)*y_t(ms,mt)  - x_t(ms,mt)*y_s(ms,mt)
   
   xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
	         + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
	         + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

   xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
	         + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
	         + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

   BigR    = x_g(ms,mt)
   BigR_x  = 1.d0
   
   Bv      = F0/BigR
   Bv_x    = - F0**2 / BigR**3
   Bv_y    = 0.
   Bv_p    = 0.
   chi_x   = 0.
   chi_y   = 0.
   chi_p   = 0.

   eps_cyl = 1.d0          ! for cylinder geometry : epscyl = eps

   do mp = 1, n_plane

     ps0    = eq_g(mp,1,ms,mt)
     ps0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_p  = eq_p(mp,1,ms,mt)
     ps0_s  = eq_s(mp,1,ms,mt)
     ps0_t  = eq_t(mp,1,ms,mt)
     ps0_pp = eq_pp(mp,1,ms,mt)
     ps0_ss = eq_ss(mp,1,ms,mt)
     ps0_tt = eq_tt(mp,1,ms,mt)
     ps0_st = eq_st(mp,1,ms,mt)
     ps0_xp = (   y_t(ms,mt) * eq_sp(mp,1,ms,mt) - y_s(ms,mt) * eq_tp(mp,1,ms,mt) ) / xjac
     ps0_yp = ( - x_t(ms,mt) * eq_sp(mp,1,ms,mt) + x_s(ms,mt) * eq_tp(mp,1,ms,mt) ) / xjac

     u0    = eq_g(mp,2,ms,mt)
     u0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_p  = eq_p(mp,2,ms,mt)
     u0_s  = eq_s(mp,2,ms,mt)
     u0_t  = eq_t(mp,2,ms,mt)

     r0    = abs(eq_g(mp,3,ms,mt))
     r0_x  = (   y_t(ms,mt) * eq_s(mp,3,ms,mt) - y_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     r0_y  = ( - x_t(ms,mt) * eq_s(mp,3,ms,mt) + x_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     r0_p  = eq_p(mp,3,ms,mt)
     r0_s  = eq_s(mp,3,ms,mt)
     r0_t  = eq_t(mp,3,ms,mt)

     T0    = abs(eq_g(mp,4,ms,mt))
     T0_x  = (   y_t(ms,mt) * eq_s(mp,4,ms,mt) - y_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     T0_y  = ( - x_t(ms,mt) * eq_s(mp,4,ms,mt) + x_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     T0_p  = eq_p(mp,4,ms,mt)
     T0_s  = eq_s(mp,4,ms,mt)
     T0_t  = eq_t(mp,4,ms,mt)

     P0    = r0 * T0
     P0_x  = r0_x * T0 + r0 * T0_x
     P0_y  = r0_y * T0 + r0 * T0_y
     P0_s  = r0_s * T0 + r0 * T0_s
     P0_t  = r0_t * T0 + r0 * T0_t
     
     ps0_xx = (ps0_ss * y_t(ms,mt)**2 - 2.d0*ps0_st * y_s(ms,mt)*y_t(ms,mt) + ps0_tt * y_s(ms,mt)**2 &
             + ps0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
             + ps0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &
             - xjac_x * (ps0_s* y_t(ms,mt) - ps0_t * y_s(ms,mt))  / xjac**2

     ps0_yy = (ps0_ss * x_t(ms,mt)**2 - 2.d0*ps0_st * x_s(ms,mt)*x_t(ms,mt) + ps0_tt * x_s(ms,mt)**2 &
             + ps0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
             + ps0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             &
             - xjac_y * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )  / xjac**2

     ps0_xy = (- ps0_ss * y_t(ms,mt)*x_t(ms,mt) - ps0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
              + ps0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                          &
              - ps0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                          &
              - ps0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
              - xjac_x * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )   / xjac**2
     
     ! --- Temperature dependent resistivity
     if ( eta_T_dependent ) then
       eta_T     =   eta * (corr_neg_temp(T0)/T_0)**(-1.5d0)
       deta_dT   = - eta * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
       d2eta_d2T =   eta * (3.75d0) * corr_neg_temp(T0)**(-3.5d0) * T_0**(1.5d0)
     else
       eta_T     = eta
       deta_dT   = 0.d0
       d2eta_d2T = 0.d0
     end if
     
     ! --- Temperature dependent viscosity
     if ( visco_T_dependent ) then
       visco_T   =   visco * (corr_neg_temp(T0)/T_0)**(-1.5d0)
       dvisco_dT = - visco * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
     else
       visco_T   = visco
       dvisco_dT = 0.d0
     end if
     
     psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
     if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
         psi_norm = 2.d0 - psi_norm
       endif
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
         psi_norm = 2.d0 - psi_norm
       endif
     endif

     D_prof  = get_dperp (psi_norm)
     ZK_prof = get_zkperp(psi_norm)
     
     do i=1,n_vertex_max

       do j=1,n_order+1

         do im=1,n_tor

           index_ij = n_tor*n_var*(n_order+1)*(i-1) + n_tor * n_var * (j-1) + im   ! index in the ELM matrix

           v   =  H(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HZ(im,mp)
           v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HZ(im,mp)

           v_s = h_s(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_t = h_t(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_p = H(i,j,ms,mt)   * element%size(i,j) * HZ_p(im,mp)
           
           v_ss = h_ss(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_tt = h_tt(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_st = h_st(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)

	       v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &
	           	+ v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
	            + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &	
		        - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

	       v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &
	           	+ v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
	            + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &	
	           	- xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

           v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                    &
       	        + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &
                - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &
	            - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2           &
                - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2

           rhs_ij_1 = v * BigR * (chi_x*u0_x + chi_y*u0_y + (F0+chi_p)*u0_p/BigR**2) / Bv**2                     * xjac * tstep &
                    + v * BigR &
                    * pbrac(ps0_x, ps0_y, ps0_p/BigR, u0_x, u0_y, u0_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv * xjac * tstep &
                    + v * BigR * eta_T * (ps0_xx + ps0_x/BigR + ps0_pp/BigR**2 + ps0_yy &
                    - (chi_x      *      (chi_x*ps0_xx           +            (chi_p+F0)*ps0_xp/BigR**2 + chi_y*ps0_xy - (chi_p+F0)*ps0_p/BigR**3) &
                    +  ((chi_p+F0)/BigR)*(chi_x*(ps0_xp/BigR-ps0_p/BigR**2) + (chi_p+F0)*ps0_pp/BigR**3 + chi_y*ps0_yp/BigR + (chi_p+F0)*ps0_x/BigR**2) &
                    +  chi_y      *      (chi_x*ps0_xy           +            (chi_p+F0)*ps0_yp/BigR**2 + chi_y*ps0_yy)) / Bv**2 &
                    + (ps0_x*Bv_x + ps0_y*Bv_y + ps0_p*Bv_p/BigR**2) / Bv)                                       * xjac * tstep &
                    + zeta * v * BigR * delta_g(mp,1,ms,mt)                                                      * xjac

           rhs_ij_2 = 0
           
           rhs_ij_3 = v * BigR * particle_source(ms,mt)                                                          * xjac * tstep &
                    - v * BigR &
                    * pbrac(r0_x, r0_y, r0_p/BigR, u0_x, u0_y, u0_p/BigR, 0., 0., F0/BigR) / Bv                  * xjac * tstep &
                    + 2 * v * BigR * r0 &
                    * pbrac(Bv_x, Bv_y, Bv_p/BigR, u0_x, u0_y, u0_p/BigR, 0., 0., F0/BigR) / Bv**2               * xjac * tstep &
                    - D_prof * BigR  * (v_x*r0_x + v_y*r0_y + v_p*r0_p * eps_cyl**2 /BigR**2 )                   * xjac * tstep &
                    + zeta * v * BigR * delta_g(mp,3,ms,mt)                                                      * xjac
           
           rhs_ij_4 = 0

           ij1 = index_ij
           ij2 = index_ij + 1*n_tor
           ij3 = index_ij + 2*n_tor
           ij4 = index_ij + 3*n_tor

           RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
           RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst
           RHS(ij3) = RHS(ij3) + rhs_ij_3 * wst
           RHS(ij4) = RHS(ij4) + rhs_ij_4 * wst

           do k=1,n_vertex_max

             do l=1,n_order+1

               do in = 1, n_tor

                 psi   = H(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 psi_x =  (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac    &
                              * element%size(k,l) * HZ(in,mp)

                 psi_y =  ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) )  / xjac   &
                              * element%size(k,l) * HZ(in,mp)

                 psi_xp = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac    &
                              * element%size(k,l) * HZ_p(in,mp)

                 psi_yp = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) )  / xjac   &
                              * element%size(k,l) * HZ_p(in,mp)

                 psi_p = H(k,l,ms,mt)   * element%size(k,l) * HZ_p(in,mp)
                 psi_pp = H(k,l,ms,mt)  * element%size(k,l) * HZ_pp(in,mp)
                 psi_s = h_s(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_t = h_t(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 
                 psi_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2  &
		                    + psi_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
	                        + psi_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &	
		                    - xjac_x * (psi_s * y_t(ms,mt) - psi_t * y_s(ms,mt)) / xjac**2

	             psi_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2  &
		                    + psi_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
	                        + psi_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2               &
		                    - xjac_y * (- psi_s * x_t(ms,mt) + psi_t * x_s(ms,mt) ) / xjac**2
		                    
		         psi_xy = (- psi_ss * y_t(ms,mt)*x_t(ms,mt) - psi_tt * x_s(ms,mt)*y_s(ms,mt)                      &
     	                + psi_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                             &
                        - psi_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                             &
	                    - psi_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2               &
                        - xjac_x * (- psi_s * x_t(ms,mt) + psi_t * x_s(ms,mt) )   / xjac**2

                 u   = psi   ;    rho   = psi    ;    T   = psi
                 u_x = psi_x ;    rho_x = psi_x  ;    T_x = psi_x
                 u_y = psi_y ;    rho_y = psi_y  ;    T_y = psi_y
                 u_p = psi_p ;    rho_p = psi_p  ;    T_p = psi_p
                 u_s = psi_s ;    rho_s = psi_s  ;    T_s = psi_s
                 u_t = psi_t ;    rho_t = psi_t  ;    T_t = psi_t
                 
                 index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix

!---------------------------------------------------------------- equation 1
                 amat_11 = v * psi * BigR * xjac * (1.d0 + zeta) &
                         - v * BigR * pbrac(psi_x, psi_y, psi_p/BigR, u0_x, u0_y, u0_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv * xjac * theta * tstep &
                         - v * BigR * eta_T * (psi_xx + psi_x/BigR + psi_pp/BigR**2 + psi_yy &
                         - (chi_x      *      (chi_x*psi_xx           +            (chi_p+F0)*psi_xp/BigR**2 + chi_y*psi_xy - (chi_p+F0)*psi_p/BigR**3) &
                         +  ((chi_p+F0)/BigR)*(chi_x*(psi_xp/BigR-psi_p/BigR**2) + (chi_p+F0)*psi_pp/BigR**3 + chi_y*psi_yp/BigR + (chi_p+F0)*psi_x/BigR**2) &
                         +  chi_y      *      (chi_x*psi_xy           +            (chi_p+F0)*psi_yp/BigR**2 + chi_y*psi_yy)) / Bv**2 &
                         + (psi_x*Bv_x + psi_y*Bv_y + psi_p*Bv_p/BigR**2) / Bv)                                                  * xjac * theta * tstep

                 amat_12 = -v * BigR * (chi_x*u_x + chi_y*u_y + (F0+chi_p)*u_p/BigR**2) / Bv**2                                  * xjac * theta * tstep &
                         -  v * BigR * pbrac(ps0_x, ps0_y, ps0_p/BigR, u_x, u_y, u_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv   * xjac * theta * tstep

!---------------------------------------------------------------- equation 2
                 amat_22 = 1.

                 amat_21 = 0.
                 amat_23 = 0.

                 amat_24 = 0.

!---------------------------------------------------------------- equation 3
                 amat_32 = v * BigR &
                         * pbrac(r0_x, r0_y, r0_p/BigR, u_x, u_y, u_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv       * xjac * theta * tstep &
                         - 2 * v * BigR * r0 &
                         * pbrac(Bv_x, Bv_y, Bv_p/BigR, u_x, u_y, u_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv**2    * xjac * theta * tstep
                 
                 amat_33 = v * rho * BigR * xjac * (1.d0 + zeta) &
                         + v * BigR &
                         * pbrac(rho_x, rho_y, rho_p/BigR, u0_x, u0_y, u0_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv * xjac * theta * tstep &
                         - 2 * v * BigR * rho &
                         * pbrac(Bv_x, Bv_y, Bv_p/BigR, u0_x, u0_y, u0_p/BigR, chi_x, chi_y, (F0+chi_p)/BigR) / Bv**2 * xjac * theta * tstep &
                         + D_prof * BigR * (v_x*rho_x + v_y*rho_y + v_p*rho_p * eps_cyl**2 /BigR**2 )                 * xjac * theta * tstep

!---------------------------------------------------------------- equation 4
                 amat_44 = 1.
                 amat_42 = 0.

                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor
                 kl3 = index_kl + 2*n_tor
                 kl4 = index_kl + 3*n_tor

                 ELM(ij1,kl1) =  ELM(ij1,kl1) + wst * amat_11
                 ELM(ij1,kl2) =  ELM(ij1,kl2) + wst * amat_12

                 ELM(ij2,kl1) =  ELM(ij2,kl1) + wst * amat_21
                 ELM(ij2,kl2) =  ELM(ij2,kl2) + wst * amat_22
                 ELM(ij2,kl3) =  ELM(ij2,kl3) + wst * amat_23
                 ELM(ij2,kl4) =  ELM(ij2,kl4) + wst * amat_24

                 ELM(ij3,kl2) =  ELM(ij3,kl2) + wst * amat_32
                 ELM(ij3,kl3) =  ELM(ij3,kl3) + wst * amat_33

                 ELM(ij4,kl2) =  ELM(ij4,kl2) + wst * amat_42
                 ELM(ij4,kl4) =  ELM(ij4,kl4) + wst * amat_44

               enddo
             enddo
           enddo

         enddo
       enddo

     enddo
   enddo

 enddo
enddo

return

end subroutine element_matrix

!> Calculates the Poisson bracket given the gradients of scalar functions a and b, and vacuum magnetic field components
!> a_p and b_p are projections of the gradient on the phi direction ((1/R)*(da/dphi)), rather than just da/dphi
pure real*8 function pbrac(a_x, a_y, a_p, b_x, b_y, b_p, Bvx, Bvy, Bvp)
  implicit none
  real*8, intent(in) :: a_x, a_y, a_p, b_x, b_y, b_p, Bvx, Bvy, Bvp
  real*8             :: Bv
  
  Bv = sqrt(Bvx**2 + Bvy**2 + Bvp**2)
  pbrac = (Bvx/Bv) * (a_y*b_p - a_p*b_y) &
         + (Bvy/Bv) * (a_x*b_p - a_p*b_x) &
         + (Bvp/Bv) * (a_x*b_y - a_y*b_x)
end function pbrac

!> Calculates the inner product of two gradient perpendicular to the vacuum field,
!> given the gradients of scalar functions a and b, and vacuum magnetic field components
!> a_p and b_p are projections of the gradient on the phi direction ((1/R)*(da/dphi)), rather than just da/dphi
pure real*8 function inprod(a_x, a_y, a_p, b_x, b_y, b_p, Bvx, Bvy, Bvp)
  implicit none
  real*8, intent(in) :: a_x, a_y, a_p, b_x, b_y, b_p, Bvx, Bvy, Bvp
  real*8             :: Bv, ab, bb, ap_x, ap_y, ap_p, bp_x, bp_y, bp_p
  
  Bv = sqrt(Bvx**2 + Bvy**2 + Bvp**2)
  ab = (Bvx*a_x + Bvy*a_y + Bvp*a_p)/Bv
  bb = (Bvx*b_x + Bvy*b_y + Bvp*b_p)/Bv
  
  ap_x = a_x - (Bvx/Bv)*ab
  ap_y = a_y - (Bvy/Bv)*ab
  ap_p = a_p - (Bvp/Bv)*ab
  bp_x = b_x - (Bvx/Bv)*bb
  bp_y = b_y - (Bvy/Bv)*bb
  bp_p = b_p - (Bvp/Bv)*bb
  
  inprod = ap_x*bp_x + ap_y*bp_y + ap_p*bp_p
end function inprod
end module mod_elt_matrix
