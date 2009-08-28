module basis_at_gaussian
!------------------------------------------------------------
! module containing the basis function on the Gaussian points
! index 1  : vertex
! index 2  : basis function (i.e. p0,u,v,w)
! index 3  : s_gauss
! index 4  : t_gauss
!-------------------------------------------------------------
use parameters
use Gauss
  real*8 :: H(n_vertex_max,n_order+1,n_gauss,n_gauss)       ! basis functions in poloidal plane
  real*8 :: H_s(n_vertex_max,n_order+1,n_gauss,n_gauss)     ! derivative with respect to first coordinate
  real*8 :: H_t(n_vertex_max,n_order+1,n_gauss,n_gauss)     ! derivative with respect to second coordinate
  real*8 :: H_st(n_vertex_max,n_order+1,n_gauss,n_gauss)    ! cross-derivative 
  real*8 :: H_ss(n_vertex_max,n_order+1,n_gauss,n_gauss)    ! second derivative, first coordinate
  real*8 :: H_tt(n_vertex_max,n_order+1,n_gauss,n_gauss)    ! second derivative, second coordinate

  real*8 :: HZ(n_tor,n_plane), HZ_p(n_tor,n_plane)          ! basis functions in toroidal direction

  real*8 :: H1(2,n_order+1,n_gauss)                         ! one dimensional basis function
  real*8 :: H1_s(2,n_order+1,n_gauss)                       ! first derivative
  real*8 :: H1_ss(2,n_order+1,n_gauss)                      ! second derivative
end module
