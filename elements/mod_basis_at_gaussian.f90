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
  real*8 :: H(n_vertex_max,n_order+1,n_gauss,n_gauss)
  real*8 :: H_s(n_vertex_max,n_order+1,n_gauss,n_gauss)
  real*8 :: H_t(n_vertex_max,n_order+1,n_gauss,n_gauss)
  real*8 :: H_st(n_vertex_max,n_order+1,n_gauss,n_gauss)
  real*8 :: H_ss(n_vertex_max,n_order+1,n_gauss,n_gauss)
  real*8 :: H_tt(n_vertex_max,n_order+1,n_gauss,n_gauss)
  real*8 :: HZ(n_tor,n_plane), HZ_p(n_tor,n_plane)
end module