module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element,nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use mod_semianalytical
use mod_equations

implicit none

type(algexpr), parameter :: u0       = algexpr(basic=.true.,var=1)
type(algexpr), parameter :: w0       = algexpr(basic=.true.,var=2)
type(algexpr), parameter :: delta_u  = algexpr(basic=.true.,var=3)
type(algexpr), parameter :: delta_w  = algexpr(basic=.true.,var=4)
type(algexpr), parameter :: v        = algexpr(basic=.true.,var=5)
type(algexpr), parameter :: u        = algexpr(basic=.true.,var=6)
type(algexpr), parameter :: w        = algexpr(basic=.true.,var=6)

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in) :: tid

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, kl1, kl2
integer    :: last
real*8     :: wst, xjac, xjac_x, xjac_y
real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), psi_norm
real*8     :: rhs_ij_1,   rhs_ij_2
real*8     :: theta, zeta

real*8     :: amat_11, amat_12, amat_21, amat_22

real*8     :: R,Z,WW, WW_R, WW_Z, WW_RR, WW_ZZ, WW_RZ
logical    :: xpoint2

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_p
real*8, dimension(:,:,:,:) , pointer :: eq_ss, eq_st, eq_tt   
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t

real*8, dimension(:,:,:,:), pointer :: eq

type(algexpr) :: rhs1
type(algexpr) :: amat11, amat12, amat21, amat22

type(action), dimension(:), allocatable, target :: rhs1seq, amat11seq, amat12seq, amat21seq, amat22seq

eq_g    => thread_struct(tid)%eq_g   
eq_s    => thread_struct(tid)%eq_s   
eq_t    => thread_struct(tid)%eq_t   
eq_p    => thread_struct(tid)%eq_p   
eq_ss   => thread_struct(tid)%eq_ss  
eq_st   => thread_struct(tid)%eq_st  
eq_tt   => thread_struct(tid)%eq_tt  
delta_g => thread_struct(tid)%delta_g
delta_s => thread_struct(tid)%delta_s
delta_t => thread_struct(tid)%delta_t

eq => thread_eq(tid)%eq

!allocate(eq(2*n_var+2,0:n_order-1,0:n_order-1,0:n_order-1))

ELM = 0.d0
RHS = 0.d0

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
zeta  = time_evol_zeta

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_ss  = 0.d0; x_st  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_ss  = 0.d0; y_st  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_ss = 0.d0; eq_st = 0.d0; eq_tt = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

mp = 1
in = 1

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

        do k=1,n_var

          eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  
          eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
          eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)

          eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)
          eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)
          eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)

          delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  
          delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
          delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) 

        enddo

      enddo
    enddo 
      
  enddo
enddo

!rhs1 = -tstep*w0*pbrack(v,u0) - visco*tstep*inprod(v,w0) - (0.25d0*tstep**2)*pbrack(w0,u0)*pbrack(v,u0) - zeta*inprod(v,delta_u)
!allocate(rhs1seq(countsubexprs(rhs1)))
!last = 0
!call buildsequence(rhs1, rhs1seq, eq, last)

!amat11 = -(1.d0 + zeta)*inprod(v,u) + theta*tstep*w0*pbrack(v,u) &
!       + (0.25d0*theta*tstep**2)*(pbrack(w0,u)*pbrack(v,u0) + pbrack(w0,u0)*pbrack(v,u))
!allocate(amat11seq(countsubexprs(amat11)))
!last = 0
!call buildsequence(amat11, amat11seq, eq, last)

!amat12 = theta*tstep*w*pbrack(v,u0) + tstep*theta*visco*inprod(v,w) + (0.25d0*theta*tstep**2)*pbrack(w,u0)*pbrack(v,u0)
!allocate(amat12seq(countsubexprs(amat12)))
!last = 0
!call buildsequence(amat12, amat12seq, eq, last)

!amat21 = inprod(v,u)
!allocate(amat21seq(countsubexprs(amat21)))
!last = 0
!call buildsequence(amat21, amat21seq, eq, last)

!amat22 = v*w
!allocate(amat22seq(countsubexprs(amat22)))
!last = 0
!call buildsequence(amat22, amat22seq, eq, last)

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  do mt=1, n_gauss

    wst = wgauss(ms)*wgauss(mt)

    xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
    
    xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)  &           
            + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                              &
            + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac
	   
    xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)  &           
            + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                              &
            + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

    eq(1:n_var,0,0,0) = eq_g(mp,:,ms,mt)
    eq(1:n_var,1,0,0) = (y_t(ms,mt)*eq_s(mp,:,ms,mt) - y_s(ms,mt)*eq_t(mp,:,ms,mt))/xjac
    eq(1:n_var,0,1,0) = (-x_t(ms,mt)*eq_s(mp,:,ms,mt) + x_s(ms,mt)*eq_t(mp,:,ms,mt))/xjac
    eq(1:n_var,2,0,0) = (eq_ss(mp,:,ms,mt)*y_t(ms,mt)**2 - 2.d0*eq_st(mp,:,ms,mt)*y_s(ms,mt)*y_t(ms,mt) &
                      + eq_tt(mp,:,ms,mt)*y_s(ms,mt)**2                                                 &
                      + eq_s(mp,:,ms,mt)*(y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt))              &
                      + eq_t(mp,:,ms,mt)*(y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt)))/xjac**2     &
                      - xjac_x*(eq_s(mp,:,ms,mt)*y_t(ms,mt) - eq_t(mp,:,ms,mt)*y_s(ms,mt))/xjac**2
    eq(1:n_var,0,2,0) = (eq_ss(mp,:,ms,mt)*x_t(ms,mt)**2 - 2.d0*eq_st(mp,:,ms,mt)*x_s(ms,mt)*x_t(ms,mt) &
                      + eq_tt(mp,:,ms,mt)*x_s(ms,mt)**2                                                 &
                      + eq_s(mp,:,ms,mt)*(x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt))              &
                      + eq_t(mp,:,ms,mt)*(x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt)))/xjac**2     &
                      - xjac_y*(-eq_s(mp,:,ms,mt)*x_t(ms,mt) + eq_t(mp,:,ms,mt)*x_s(ms,mt))/xjac**2
    eq(1:n_var,1,1,0) = (-eq_ss(mp,:,ms,mt)*y_t(ms,mt)*x_t(ms,mt) - eq_tt(mp,:,ms,mt)*x_s(ms,mt)*y_s(ms,mt) &
                      + eq_st(mp,:,ms,mt)*(y_s(ms,mt)*x_t(ms,mt) + y_t(ms,mt)*x_s(ms,mt))                   &
                      - eq_s(mp,:,ms,mt)*(x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt))                  &
                      - eq_t(mp,:,ms,mt)*(x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt)))/xjac**2         &
                      - xjac_x*(-eq_s(mp,:,ms,mt)*x_t(ms,mt) + eq_t(mp,:,ms,mt)*x_s(ms,mt))/xjac**2

    eq(n_var+1:2*n_var,0,0,0) = delta_g(mp,:,ms,mt)
    eq(n_var+1:2*n_var,1,0,0) = (y_t(ms,mt)*delta_s(mp,:,ms,mt) - y_s(ms,mt)*delta_t(mp,:,ms,mt))/xjac
    eq(n_var+1:2*n_var,0,1,0) = (-x_t(ms,mt)*delta_s(mp,:,ms,mt) + x_s(ms,mt)*delta_t(mp,:,ms,mt))/xjac

    do i=1,n_vertex_max

      do j=1,n_order+1

        index_ij = n_var*(n_order+1)*(i-1) + n_var * (j-1) + 1   ! index in the ELM matrix

        eq(2*n_var+1,0,0,0) = h(i,j,ms,mt)*element%size(i,j)
        eq(2*n_var+1,1,0,0) = (y_t(ms,mt)*h_s(i,j,ms,mt) - y_s(ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)/xjac
        eq(2*n_var+1,0,1,0) = (-x_t(ms,mt)*h_s(i,j,ms,mt) + x_s(ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)/xjac
        eq(2*n_var+1,2,0,0) = (h_ss(i,j,ms,mt)*y_t(ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*y_s(ms,mt)*y_t(ms,mt)                 &
                            + h_tt(i,j,ms,mt)*y_s(ms,mt)**2                                                               &
                            + h_s(i,j,ms,mt)*(y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt))                            &
                            + h_t(i,j,ms,mt)*(y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt)))*element%size(i,j)/xjac**2 &
                            - xjac_x*(h_s(i,j,ms,mt)*y_t(ms,mt) - h_t(i,j,ms,mt)*y_s(ms,mt))*element%size(i,j)/xjac**2
        eq(2*n_var+1,0,2,0) = (h_ss(i,j,ms,mt)*x_t(ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*x_s(ms,mt)*x_t(ms,mt)                 &
                            + h_tt(i,j,ms,mt)*x_s(ms,mt)**2                                                               &
                            + h_s(i,j,ms,mt)*(x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt))                            &
                            + h_t(i,j,ms,mt)*(x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt)))*element%size(i,j)/xjac**2 &
                            - xjac_y*(-h_s(i,j,ms,mt)*x_t(ms,mt) + h_t(i,j,ms,mt)*x_s(ms,mt))*element%size(i,j)/xjac**2


        rhs_ij_1 = (eval(thread_eq(tid)%rhs1dt0seq) + eval(thread_eq(tid)%rhs1dt1seq)*tstep + eval(thread_eq(tid)%rhs1dt2seq)*tstep**2)*xjac
                 ! -eq(2,0,0,0)*(h_s(i,j,ms,mt)*eq_t(mp,1,ms,mt) - h_t(i,j,ms,mt)*eq_s(mp,1,ms,mt))*element%size(i,j)*tstep &
                 ! - visco*(eq(2*n_var+1,1,0,0)*eq(2,1,0,0) + eq(2*n_var+1,0,1,0)*eq(2,0,1,0))*xjac*tstep &
                 ! - 0.25d0*(eq(2,1,0,0)*eq(1,0,1,0) - eq(2,0,1,0)*eq(1,1,0,0))*(eq(2*n_var+1,1,0,0)*eq(1,0,1,0) - eq(2*n_var+1,0,1,0)*eq(1,1,0,0))*xjac*tstep**2 ! &
                 ! - zeta*(eq(2*n_var+1,1,0,0)*eq(3,1,0,0) + eq(2*n_var+1,0,1,0)*eq(3,0,1,0))*xjac
           
        rhs_ij_2 = 0.d0 

        ij1 = index_ij
        ij2 = index_ij + 1

        RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
        RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst

        do k=1,n_vertex_max

          do l=1,n_order+1
            
            eq(2*n_var+2,0,0,0) = h(k,l,ms,mt)*element%size(k,l)
            eq(2*n_var+2,1,0,0) = (y_t(ms,mt)*h_s(k,l,ms,mt) - y_s(ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)/xjac
            eq(2*n_var+2,0,1,0) = (-x_t(ms,mt)*h_s(k,l,ms,mt) + x_s(ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)/xjac
            eq(2*n_var+2,2,0,0) = (h_ss(k,l,ms,mt)*y_t(ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*y_s(ms,mt)*y_t(ms,mt)                 &
                                + h_tt(k,l,ms,mt)*y_s(ms,mt)**2                                                               &
                                + h_s(k,l,ms,mt)*(y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt))                            &
                                + h_t(k,l,ms,mt)*(y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt)))*element%size(k,l)/xjac**2 &
                                - xjac_x*(h_s(k,l,ms,mt)*y_t(ms,mt) - h_t(k,l,ms,mt)*y_s(ms,mt))*element%size(k,l)/xjac**2
            eq(2*n_var+2,0,2,0) = (h_ss(k,l,ms,mt)*x_t(ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*x_s(ms,mt)*x_t(ms,mt)                 &
                                + h_tt(k,l,ms,mt)*x_s(ms,mt)**2                                                               &
                                + h_s(k,l,ms,mt)*(x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt))                            &
                                + h_t(k,l,ms,mt)*(x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt)))*element%size(k,l)/xjac**2 &
                                - xjac_y*(-h_s(k,l,ms,mt)*x_t(ms,mt) + h_t(k,l,ms,mt)*x_s(ms,mt))*element%size(k,l)/xjac**2
            
		
            index_kl = n_var*(n_order+1)*(k-1) + n_var * (l-1) + 1   ! index in the ELM matrix

!---------------------------------------------------------------- equation 1
 		      
            amat_11 = (eval(thread_eq(tid)%amat11dt0seq) + eval(thread_eq(tid)%amat11dt1seq)*tstep + eval(thread_eq(tid)%amat11dt2seq)*tstep**2)*xjac
                  ! - (eq(2*n_var+1,1,0,0)*eq(2*n_var+2,1,0,0) + eq(2*n_var+1,0,1,0)*eq(2*n_var+2,0,1,0))*xjac*(1.d0 + zeta)                 &  
                  ! + eq(2,0,0,0)*(h_s(i,j,ms,mt)*h_t(k,l,ms,mt) - h_t(i,j,ms,mt)*h_s(k,l,ms,mt))*element%size(i,j)*element%size(k,l)*theta*tstep &
                  ! + 0.25d0*(eq(2,1,0,0)*eq(2*n_var+2,0,1,0) - eq(2,0,1,0)*eq(2*n_var+2,1,0,0))*(eq(2*n_var+1,1,0,0)*eq(1,0,1,0) - eq(2*n_var+1,0,1,0)*eq(1,1,0,0))*xjac*theta*tstep**2 &
		              ! + 0.25d0*(eq(2,1,0,0)*eq(1,0,1,0) - eq(2,0,1,0)*eq(1,1,0,0))*(eq(2*n_var+1,1,0,0)*eq(2*n_var+2,0,1,0) - eq(2*n_var+1,0,1,0)*eq(2*n_var+2,1,0,0))*xjac*theta*tstep**2
     
            amat_12 = (eval(thread_eq(tid)%amat12dt1seq)*tstep + eval(thread_eq(tid)%amat12dt2seq)*tstep**2)*xjac
                    ! eq(2*n_var+2,0,0,0)*(h_s(i,j,ms,mt)*eq_t(mp,1,ms,mt) - h_t(i,j,ms,mt)*eq_s(mp,1,ms,mt))*element%size(i,j)*theta*tstep &
                    ! + visco*(eq(2*n_var+1,1,0,0)*eq(2*n_var+2,1,0,0) + eq(2*n_var+1,0,1,0)*eq(2*n_var+2,0,1,0))*xjac*theta*tstep &
                    ! + 0.25d0*(eq(2*n_var+2,1,0,0)*eq(1,0,1,0) - eq(2*n_var+2,0,1,0)*eq(1,1,0,0))*(eq(2*n_var+1,1,0,0)*eq(1,0,1,0) - eq(2*n_var+1,0,1,0)*eq(1,1,0,0))*xjac*theta*tstep**2
		      
!---------------------------------------------------------------- equation 2
            amat_22 = eval(thread_eq(tid)%amat22dt0seq)*xjac
            amat_21 = eval(thread_eq(tid)%amat21dt0seq)*xjac

            kl1 = index_kl
            kl2 = index_kl + 1

            ELM(ij1,kl1) =  ELM(ij1,kl1) + wst * amat_11
            ELM(ij1,kl2) =  ELM(ij1,kl2) + wst * amat_12

            ELM(ij2,kl1) =  ELM(ij2,kl1) + wst * amat_21
            ELM(ij2,kl2) =  ELM(ij2,kl2) + wst * amat_22

          enddo
        enddo

      enddo
    enddo
  enddo
enddo

return

!contains

!type(algexpr) function pbrack(a,b)
!  implicit none
!  type(algexpr), intent(in) :: a, b
  
!  pbrack = dx(a)*dy(b) - dy(a)*dx(b)
!end function pbrack

!type(algexpr) function inprod(a,b)
!  implicit none
!  type(algexpr), intent(in) :: a, b
  
!  inprod = dx(a)*dx(b) + dy(a)*dy(b)
!end function inprod
end subroutine element_matrix
end module mod_elt_matrix
