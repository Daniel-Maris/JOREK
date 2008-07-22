module phys_module
use parameters
  real*8  :: tstep, tstep_in, eta, visco, visco_par
  real*8  :: R_begin, R_end, Z_begin, Z_end, R_geo, Z_geo, amin, fbnd(1026), fpsi(1026)
  real*8  :: ellip,tria_u,tria_l,quad_u,quad_l,xampl,xwidth,xsig,xtheta,xshift,xleft
  real*8  :: xr1, sig1, xr2, sig2
  real*8  :: F0, GAMMA
  real*8  :: zjz_0, zjz_1,  zj_coef(10)
  real*8  :: T_0,   T_1,    T_coef(10)
  real*8  :: rho_0, rho_1,  rho_coef(10)
  real*8  :: FF_0,  FF_1,   FF_coef(10)
  real*8  :: particlesource, heatsource, ZK_perp(10), ZK_par, D_perp(10), D_par, eta_num, visco_num
  real*8  :: t_start, t_now
  integer :: nstep, n_R, n_Z, n_radial, n_pol, n_tht, n_flux
  integer :: n_open, n_private, n_leg
  integer :: mf, index_start, index_now, mode(n_tor), nout
  logical :: restart, regrid, import_equil, xpoint
  real*8, allocatable   :: xtime(:), energies(:,:,:)
end module