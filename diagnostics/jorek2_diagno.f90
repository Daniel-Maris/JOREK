!**********************************************************************
!* program to extract data from a JOREK2 restart file                 *
!**********************************************************************

program jorek2_diagno
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

namelist /in1/  tstep, nstep, eta, visco, visco_par,                &
                restart,  regrid,                                   &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open,n_private,n_leg,  nout,                      &
                xr1, sig1, xr2, sig2,                               &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                F0,                                                 &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_perp, D_par, D_perp,                     &
                particlesource, heatsource,                         &
                eta_num, visco_num,                                 &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint


call initialise_basis                              ! define the basis functions at the Gaussian points

call export_helena(node_list,element_list)

!----------------------------------------- plot profiles
call begplt('profiles.ps')

call plot_velocity_profile(node_list,element_list, 3.d0, 0.d0, 3.d0, 2.d0)

call finplt

end


