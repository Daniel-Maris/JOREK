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
use mod_semianalytical
use mod_equations

implicit none

type (type_element), intent(in)   :: element
type (type_node)   , intent(in)   :: nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in) :: tid

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, kl1, kl2, kl3, kl4, kl5, kl6
real*8     :: wst,  xjac, xjac_x, xjac_y, xjac_s, xjac_t, BigR, r2, phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: psi_norm
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6
real*8     :: theta, zeta, delta_u_x, delta_u_y

real*8     :: amat_11, amat_12, amat_13, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31
real*8     :: amat_41, amat_42, amat_43, amat_44, amat_52, amat_55, amat_66

real*8     :: F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz, FFprime_prof, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz

logical    :: xpoint2

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt
#ifdef altcs
real*8, dimension(n_gauss,n_gauss)    :: psieq, psieq_s, psieq_t
real*8                                :: psieq_x, psieq_y
#endif

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_st, eq_ss, eq_tt
real*8, dimension(:,:,:,:) , pointer :: eq_p, eq_pp, eq_sp, eq_tp
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t, delta_p

real*8, dimension(:,:,:,:), pointer :: eq

real*8 :: rhsdt0, rhsdt1, amatdt0, amatdt1

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
delta_p => thread_struct(tid)%delta_p

eq => thread_eq(tid)%eq

ELM = 0.d0
RHS = 0.d0

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
zeta  = time_evol_zeta

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
#ifdef altcs
psieq = 0.d0; psieq_s = 0.d0; psieq_t = 0.d0
#endif
x_g  = 0.d0; x_s   = 0.d0; x_t   = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s   = 0.d0; y_t   = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s  = 0.d0; eq_t  = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0;
eq_p = 0.d0; eq_pp = 0.d0; eq_sp = 0.d0; eq_tp = 0.d0

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0; delta_p = 0.d0

eq = 0.d0

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
       
#ifdef altcs
       psieq(ms,mt)   = psieq(ms,mt)   + nodes(i)%psi_eq(j)*element%size(i,j)*H(i,j,ms,mt)
       psieq_s(ms,mt) = psieq_s(ms,mt) + nodes(i)%psi_eq(j)*element%size(i,j)*H_s(i,j,ms,mt)
       psieq_t(ms,mt) = psieq_t(ms,mt) + nodes(i)%psi_eq(j)*element%size(i,j)*H_t(i,j,ms,mt)
#endif

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

   BigR = x_g(ms,mt)
   
#ifdef altcs
   psieq_x = ( y_t(ms,mt)*psieq_s(ms,mt) - y_s(ms,mt)*psieq_t(ms,mt))/xjac
   psieq_y = (-x_t(ms,mt)*psieq_s(ms,mt) + x_s(ms,mt)*psieq_t(ms,mt))/xjac
#endif

   do mp = 1, n_plane
     phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)

     ! delta_u^n
     eq(1:n_var,0,0,0) = eq_g(mp,:,ms,mt)
     eq(1:n_var,1,0,0) = (y_t(ms,mt)*eq_s(mp,:,ms,mt) - y_s(ms,mt)*eq_t(mp,:,ms,mt))/xjac
     eq(1:n_var,0,1,0) = (-x_t(ms,mt)*eq_s(mp,:,ms,mt) + x_s(ms,mt)*eq_t(mp,:,ms,mt))/xjac
     eq(1:n_var,0,0,1) = eq_p(mp,:,ms,mt)
     eq(1:n_var,2,0,0) = (eq_ss(mp,:,ms,mt)*y_t(ms,mt)**2 - 2.d0*eq_st(mp,:,ms,mt)*y_s(ms,mt)*y_t(ms,mt)                      &
                       + eq_tt(mp,:,ms,mt)*y_s(ms,mt)**2 + eq_s(mp,:,ms,mt)*(y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt)) &
                       + eq_t(mp,:,ms,mt)*(y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt)))/xjac**2                          &
                       - xjac_x*(eq_s(mp,:,ms,mt)*y_t(ms,mt) - eq_t(mp,:,ms,mt)*y_s(ms,mt))/xjac**2
     eq(1:n_var,0,2,0) = (eq_ss(mp,:,ms,mt)*x_t(ms,mt)**2 - 2.d0*eq_st(mp,:,ms,mt)*x_s(ms,mt)*x_t(ms,mt)                      &
                       + eq_tt(mp,:,ms,mt)*x_s(ms,mt)**2 + eq_s(mp,:,ms,mt)*(x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt)) &
                       + eq_t(mp,:,ms,mt)*(x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt)))/xjac**2                          &
                       - xjac_y*(-eq_s(mp,:,ms,mt)*x_t(ms,mt) + eq_t(mp,:,ms,mt)*x_s(ms,mt))/xjac**2
     eq(1:n_var,1,1,0) = (-eq_ss(mp,:,ms,mt)*y_t(ms,mt)*x_t(ms,mt) - eq_tt(mp,:,ms,mt)*x_s(ms,mt)*y_s(ms,mt) &
                       + eq_st(mp,:,ms,mt)*(y_s(ms,mt)*x_t(ms,mt) + y_t(ms,mt)*x_s(ms,mt))                   &
                       - eq_s(mp,:,ms,mt)*(x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt))                  &
                       - eq_t(mp,:,ms,mt)*(x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt)))/xjac**2         &
                       - xjac_x*(-eq_s(mp,:,ms,mt)*x_t(ms,mt) + eq_t(mp,:,ms,mt)*x_s(ms,mt))/xjac**2
     eq(1:n_var,0,0,2) = eq_pp(mp,:,ms,mt)
     eq(1:n_var,1,0,1) = (y_t(ms,mt)*eq_sp(mp,:,ms,mt) - y_s(ms,mt)*eq_tp(mp,:,ms,mt))/xjac
     eq(1:n_var,0,1,1) = (-x_t(ms,mt)*eq_sp(mp,:,ms,mt) + x_s(ms,mt)*eq_tp(mp,:,ms,mt))/xjac
     
     ! delta_u^(n-1)
     eq(n_var+1:2*n_var,0,0,0) = delta_g(mp,:,ms,mt)
     eq(n_var+1:2*n_var,1,0,0) = (y_t(ms,mt)*delta_s(mp,:,ms,mt) - y_s(ms,mt)*delta_t(mp,:,ms,mt))/xjac
     eq(n_var+1:2*n_var,0,1,0) = (-x_t(ms,mt)*delta_s(mp,:,ms,mt) + x_s(ms,mt)*delta_t(mp,:,ms,mt))/xjac
     eq(n_var+1:2*n_var,0,0,1) = delta_p(mp,:,ms,mt)
     
     eq(2*n_var+3,0,0,0) = F0*phi                     ! chi
     eq(2*n_var+3,0,0,1) = F0
     eq(2*n_var+4,0,0,0) = x_g(ms,mt)                 ! psi_v
     eq(2*n_var+4,1,0,0) = 1.d0
     eq(2*n_var+5,0,0,0) = x_g(ms,mt)                 ! R
     eq(2*n_var+5,1,0,0) = 1.d0
     
     psi_norm = (eq(1,0,0,0) - psi_axis)/(psi_bnd - psi_axis)
     if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
         psi_norm = 2.d0 - psi_norm
       endif
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
         psi_norm = 2.d0 - psi_norm
       endif
     endif
     
     ! The Psi in the equations differs by a factor of F0 from the normal JOREK Psi
     eq(1,:,:,:) = eq(1,:,:,:)/F0
     eq(n_var+1,:,:,:) = eq(n_var+1,:,:,:)/F0
     eq(3,:,:,:) = eq(3,:,:,:)/F0
     eq(n_var+3,:,:,:) = eq(n_var+3,:,:,:)/F0

     eq(2*n_var+6,0,0,0) = get_dperp(psi_norm)       ! D_perp
     eq(2*n_var+7,0,0,0) = get_zkperp(psi_norm)      ! k_perp
     if (mp .eq. 1) then ! these quantities are the same on all poloidal planes
       eq(2*n_var+8,0,0,0) = particle_source(ms,mt)    ! S_rho
       eq(2*n_var+9,0,0,0) = heat_source(ms,mt)       ! S_e
       eq(2*n_var+10,0,0,0) = current_source(ms,mt)/F0 ! S_j
#ifdef altcs
       call F_profile(xpoint2,xcase2,y_g(ms,mt),Z_xpoint,psieq(ms,mt),psi_axis,psi_bnd,F_prof,dF_dpsi,dF_dz,dF_dpsi2,dF_dz2, &
                      dF_dpsi_dz,FFprime_prof,dFF_dpsi,dFF_dz,dFF_dpsi2,dFF_dz2,dFF_dpsi_dz)
       eq(2*n_var+11,0,0,0) = F_prof
       eq(2*n_var+11,1,0,0) = dF_dpsi*psieq_x
       eq(2*n_var+11,0,1,0) = dF_dpsi*psieq_y
#else
       call F_profile(xpoint2,xcase2,y_g(ms,mt),Z_xpoint,F0*eq(1,0,0,0),psi_axis,psi_bnd,F_prof,dF_dpsi,dF_dz,dF_dpsi2,dF_dz2, &
                      dF_dpsi_dz,FFprime_prof,dFF_dpsi,dFF_dz,dFF_dpsi2,dFF_dz2,dFF_dpsi_dz)
       eq(2*n_var+11,0,0,0) = F_prof
       eq(2*n_var+11,1,0,0) = dF_dpsi*F0*eq(1,1,0,0)
       eq(2*n_var+11,0,1,0) = dF_dpsi*F0*eq(1,0,1,0)
#endif
     end if
     ! Resistivity
     if (eta_T_dependent) then
       eq(2*n_var+12,0,0,0) = eta*(corr_neg_temp(eq(6,0,0,0))/T_0)**(-1.5d0)               ! eta
       eq(2*n_var+13,0,0,0) = -1.5d0*eta*corr_neg_temp(eq(6,0,0,0))**(-2.5d0)*T_0**(1.5d0) ! deta/dT
     else
       eq(2*n_var+12,0,0,0) = eta
       eq(2*n_var+13,0,0,0) = 0.d0
     end if
     ! Viscosity
     if (visco_T_dependent) then
       eq(2*n_var+14,0,0,0) = visco*(corr_neg_temp(eq(6,0,0,0))/T_0)**(-1.5d0)               ! visco
       eq(2*n_var+15,0,0,0) = -1.5d0*visco*corr_neg_temp(eq(6,0,0,0))**(-2.5d0)*T_0**(1.5d0) ! dvisco/dT
     else
       eq(2*n_var+14,0,0,0) = visco
       eq(2*n_var+15,0,0,0) = 0.d0
     end if
     
     ! Auxiliary variables (aux)
#ifdef DEBUG
     eq(2*n_var+16,0,0,0) = eval(thread_eq(tid)%aBv2seq); eq(2*n_var+16,1,0,0) = eval(thread_eq(tid)%aBv2xseq)
     eq(2*n_var+16,0,1,0) = eval(thread_eq(tid)%aBv2yseq); eq(2*n_var+16,0,0,1) = eval(thread_eq(tid)%aBv2pseq)
     eq(2*n_var+17,0,0,0) = eval(thread_eq(tid)%aj0xseq)
     eq(2*n_var+18,0,0,0) = eval(thread_eq(tid)%aj0yseq)
     eq(2*n_var+19,0,0,0) = eval(thread_eq(tid)%aj0pseq)
     eq(2*n_var+20,0,0,0) = eval(thread_eq(tid)%aj0chiseq)
#else
#include "aux_unreadable.h"
#endif
     
     do i=1,n_vertex_max

       do j=1,n_order+1

         do im=1,n_tor

           index_ij = n_tor*n_var*(n_order+1)*(i-1) + n_tor*n_var*(j-1) + im   ! index in the ELM matrix

           eq(2*n_var+1,0,0,0) =  H(i,j,ms,mt)*element%size(i,j)*HZ(im,mp)
           eq(2*n_var+1,1,0,0) = (y_t(ms,mt)*h_s(i,j,ms,mt) - y_s(ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac
           eq(2*n_var+1,0,1,0) = (-x_t(ms,mt)*h_s(i,j,ms,mt) + x_s(ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac
           eq(2*n_var+1,0,0,1) = H(i,j,ms,mt)*element%size(i,j)*HZ_p(im,mp)
           eq(2*n_var+1,2,0,0) = (h_ss(i,j,ms,mt)*y_t(ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*y_s(ms,mt)*y_t(ms,mt)                           &
	                             + h_tt(i,j,ms,mt)*y_s(ms,mt)**2 + h_s(i,j,ms,mt)*(y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt))      &
                               + h_t(i,j,ms,mt)*(y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt)))*element%size(i,j)*HZ(im,mp)/xjac**2 &
                               - xjac_x*(h_s(i,j,ms,mt)*y_t(ms,mt) - h_t(i,j,ms,mt)*y_s(ms,mt))*element%size(i,j)*HZ(im,mp)/xjac**2
           eq(2*n_var+1,0,2,0) = (h_ss(i,j,ms,mt)*x_t(ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*x_s(ms,mt)*x_t(ms,mt)                           &
                               + h_tt(i,j,ms,mt)*x_s(ms,mt)**2 + h_s(i,j,ms,mt)*(x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt))      &
                               + h_t(i,j,ms,mt)*(x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt)))*element%size(i,j)*HZ(im,mp)/xjac**2 &
           	                   - xjac_y*(-h_s(i,j,ms,mt)*x_t(ms,mt) + h_t(i,j,ms,mt)*x_s(ms,mt))*element%size(i,j)*HZ(im,mp)/xjac**2
           eq(2*n_var+1,0,0,2) = H(i,j,ms,mt)*element%size(i,j)*HZ_pp(im,mp)
           eq(2*n_var+1,1,1,0) = (-h_ss(i,j,ms,mt)*y_t(ms,mt)*x_t(ms,mt) - h_tt(i,j,ms,mt)*x_s(ms,mt)*y_s(ms,mt)                       &
       	                       + h_st(i,j,ms,mt)*(y_s(ms,mt)*x_t(ms,mt) + y_t(ms,mt)*x_s(ms,mt))                                       &
                               - h_s(i,j,ms,mt)*(x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt))                                      &
                               - h_t(i,j,ms,mt)*(x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt)))*element%size(i,j)*HZ(im,mp)/xjac**2 &
                               - xjac_x*(-h_s(i,j,ms,mt)*x_t(ms,mt) + h_t(i,j,ms,mt)*x_s(ms,mt))*element%size(i,j)*HZ(im,mp)/xjac**2
           eq(2*n_var+1,1,0,1) = (y_t(ms,mt)*h_s(i,j,ms,mt) - y_s(ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ_p(im,mp)/xjac
           eq(2*n_var+1,0,1,1) = (-x_t(ms,mt)*h_s(i,j,ms,mt) + x_s(ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ_p(im,mp)/xjac

#ifdef DEBUG
           rhs_ij_1 = eval(thread_eq(tid)%rhs1seq)*BigR*xjac
           rhs_ij_2 = eval(thread_eq(tid)%rhs2seq)*BigR*xjac
           rhs_ij_3 = 0.d0 ! eval(thread_eq(tid)%rhs3seq)*xjac*freeb_fact/BigR
           rhs_ij_4 = 0.d0 ! (eval(thread_eq(tid)%rhs4dt0seq) + eval(thread_eq(tid)%rhs4dt1seq)*tstep)*xjac
           rhs_ij_5 = eval(thread_eq(tid)%rhs5seq)*BigR*xjac
           rhs_ij_6 = 0.d0
#else
#include "rhs_unreadable.h"

           rhs_ij_1 = rhs_ij_1*BigR*xjac
           rhs_ij_2 = rhs_ij_2*BigR*xjac
           rhs_ij_3 = 0.d0
           rhs_ij_4 = 0.d0
           rhs_ij_5 = rhs_ij_5*BigR*xjac
           rhs_ij_6 = 0.d0
#endif

           ij1 = index_ij
           ij2 = index_ij + 1*n_tor
           ij3 = index_ij + 2*n_tor
           ij4 = index_ij + 3*n_tor
           ij5 = index_ij + 4*n_tor
           ij6 = index_ij + 5*n_tor

           RHS(ij1) = RHS(ij1) + rhs_ij_1*wst
           RHS(ij2) = RHS(ij2) + rhs_ij_2*wst
           RHS(ij3) = RHS(ij3) + rhs_ij_3*wst
           RHS(ij4) = RHS(ij4) + rhs_ij_4*wst
           RHS(ij5) = RHS(ij5) + rhs_ij_5*wst
           RHS(ij6) = RHS(ij6) + rhs_ij_6*wst

           do k=1,n_vertex_max

             do l=1,n_order+1

               do in = 1, n_tor

                 eq(2*n_var+2,0,0,0) = H(k,l,ms,mt)*element%size(k,l)*HZ(in,mp)
                 eq(2*n_var+2,1,0,0) = (y_t(ms,mt)*h_s(k,l,ms,mt) - y_s(ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac
                 eq(2*n_var+2,0,1,0) = (-x_t(ms,mt)*h_s(k,l,ms,mt) + x_s(ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac
                 eq(2*n_var+2,0,0,1) = H(k,l,ms,mt)*element%size(k,l)*HZ_p(in,mp)
                 eq(2*n_var+2,2,0,0) = (h_ss(k,l,ms,mt)*y_t(ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*y_s(ms,mt)*y_t(ms,mt) + h_tt(k,l,ms,mt)*y_s(ms,mt)**2 &
                                     + h_s(k,l,ms,mt)*(y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt))                                            &
                                     + h_t(k,l,ms,mt)*(y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt)))*element%size(k,l)*HZ(in,mp)/xjac**2       &	
                                     - xjac_x*(h_s(k,l,ms,mt)*y_t(ms,mt) - h_t(k,l,ms,mt)*y_s(ms,mt))*element%size(k,l)*HZ(in,mp)/xjac**2
                 eq(2*n_var+2,0,2,0) = (h_ss(k,l,ms,mt)*x_t(ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*x_s(ms,mt)*x_t(ms,mt) + h_tt(k,l,ms,mt)*x_s(ms,mt)**2 &
                                     + h_s(k,l,ms,mt)*(x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt))                                            &
                                     + h_t(k,l,ms,mt)*(x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt)))*element%size(k,l)*HZ(in,mp)/xjac**2       &
                                     - xjac_y*(-h_s(k,l,ms,mt)*x_t(ms,mt) + h_t(k,l,ms,mt)*x_s(ms,mt))*element%size(k,l)*HZ(in,mp)/xjac**2
                 eq(2*n_var+2,0,0,2) = H(k,l,ms,mt)*element%size(k,l)*HZ_pp(in,mp)
                 eq(2*n_var+2,1,1,0) = (-h_ss(k,l,ms,mt)*y_t(ms,mt)*x_t(ms,mt) - h_tt(k,l,ms,mt)*x_s(ms,mt)*y_s(ms,mt)                       &
     	                               + h_st(k,l,ms,mt)*(y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt))                                      &
                                     - h_s(k,l,ms,mt)*(x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt))                                      &
                                     - h_t(k,l,ms,mt)*(x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt)))*element%size(k,l)*HZ(in,mp)/xjac**2 &
                                     - xjac_x*(-h_s(k,l,ms,mt)*x_t(ms,mt) + h_t(k,l,ms,mt)*x_s(ms,mt))*element%size(k,l)*HZ(in,mp)/xjac**2
                 eq(2*n_var+2,1,0,1) = (y_t(ms,mt)*h_s(k,l,ms,mt) - y_s(ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ_p(in,mp)/xjac
                 eq(2*n_var+2,0,1,1) = (-x_t(ms,mt)*h_s(k,l,ms,mt) + x_s(ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ_p(in,mp)/xjac
                 
                 index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix
                 
#ifdef DEBUG
!---------------------------------------------------------------- Auxiliary variables involving unknowns (aux2)
                 eq(2*n_var+21,0,0,0) = eval(thread_eq(tid)%atjxseq)
                 eq(2*n_var+22,0,0,0) = eval(thread_eq(tid)%atjyseq)
                 eq(2*n_var+23,0,0,0) = eval(thread_eq(tid)%atjpseq)
                 eq(2*n_var+24,0,0,0) = eval(thread_eq(tid)%atjchiseq)
                 
!---------------------------------------------------------------- equation 1
                 amat_11 = eval(thread_eq(tid)%amat11seq)*BigR*xjac/F0
                 amat_12 = eval(thread_eq(tid)%amat12seq)*BigR*xjac
                 amat_13 = eval(thread_eq(tid)%amat13seq)*BigR*xjac/F0

!---------------------------------------------------------------- equation 2
                 amat_21 = eval(thread_eq(tid)%amat21seq)*BigR*xjac/F0
                 amat_22 = eval(thread_eq(tid)%amat22seq)*BigR*xjac
                 amat_23 = eval(thread_eq(tid)%amat23seq)*BigR*xjac/F0
                 amat_24 = eval(thread_eq(tid)%amat24seq)*BigR*xjac
                 amat_25 = eval(thread_eq(tid)%amat25seq)*BigR*xjac

!---------------------------------------------------------------- equation 3
                 amat_31 = eval(thread_eq(tid)%amat31seq)*xjac/(BigR*F0)
                 amat_33 = eval(thread_eq(tid)%amat33seq)*xjac/(BigR*F0)

!---------------------------------------------------------------- equation 4
                 amat_42 = eval(thread_eq(tid)%amat42seq)*BigR*xjac
                 amat_44 = eval(thread_eq(tid)%amat44seq)*BigR*xjac
                 
                 amat_52 = eval(thread_eq(tid)%amat52seq)*BigR*xjac
                 amat_55 = eval(thread_eq(tid)%amat55seq)*BigR*xjac
                 
                 amat_66 = 1.d0
#else
#include "aux2_unreadable.h"
#include "amat_unreadable.h"

                 amat_11 = amat_11*BigR*xjac/F0; amat_12 = amat_12*BigR*xjac; amat_13 = amat_13*BigR*xjac/F0
                 amat_21 = amat_21*BigR*xjac/F0; amat_22 = amat_22*BigR*xjac; amat_23 = amat_23*BigR*xjac/F0; amat_24 = amat_24*BigR*xjac
                 amat_25 = amat_25*BigR*xjac
                 amat_31 = amat_31*xjac/(BigR*F0); amat_33 = amat_33*xjac/(BigR*F0)
                 amat_42 = amat_42*BigR*xjac; amat_44 = amat_44*BigR*xjac
                 amat_52 = amat_52*BigR*xjac; amat_55 = amat_55*BigR*xjac
                 amat_66 = 1.d0
#endif

                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor
                 kl3 = index_kl + 2*n_tor
                 kl4 = index_kl + 3*n_tor
                 kl5 = index_kl + 4*n_tor
                 kl6 = index_kl + 5*n_tor

                 ELM(ij1,kl1) =  ELM(ij1,kl1) + wst*amat_11
                 ELM(ij1,kl2) =  ELM(ij1,kl2) + wst*amat_12
                 ELM(ij1,kl3) =  ELM(ij1,kl3) + wst*amat_13

                 ELM(ij2,kl1) =  ELM(ij2,kl1) + wst*amat_21
                 ELM(ij2,kl2) =  ELM(ij2,kl2) + wst*amat_22
                 ELM(ij2,kl3) =  ELM(ij2,kl3) + wst*amat_23
                 ELM(ij2,kl4) =  ELM(ij2,kl4) + wst*amat_24
                 ELM(ij2,kl5) =  ELM(ij2,kl5) + wst*amat_25

                 ELM(ij3,kl1) =  ELM(ij3,kl1) + wst*amat_31
                 ELM(ij3,kl3) =  ELM(ij3,kl3) + wst*amat_33

                 ELM(ij4,kl2) =  ELM(ij4,kl2) + wst*amat_42
                 ELM(ij4,kl4) =  ELM(ij4,kl4) + wst*amat_44
                 
                 ELM(ij5,kl2) =  ELM(ij5,kl2) + wst*amat_52
                 ELM(ij5,kl5) =  ELM(ij5,kl5) + wst*amat_55
                 
                 ELM(ij6,kl6) =  ELM(ij6,kl6) + wst*amat_66

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
end module mod_elt_matrix
