!> Calculates n.B on the plasma boundary

subroutine determine_boundary_flux(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint,q,rad)

use constants
use tr_module 
use data_structure
use phys_module
use equil_info, only : get_psi_n
use mod_interp


implicit none

! --- Gaussian points between (-1.,1.) for Gauss-integration
real*8, parameter :: xgs(4) = (/-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /)
real*8, parameter :: wgs(4) = (/ 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /)

! --- Input parameters.
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
type (type_surface_list), intent(in)    :: surface_list
real*8,                   intent(in)    :: psi_axis
real*8,                   intent(in)    :: psi_xpoint(2)
real*8,                   intent(in)    :: Z_xpoint(2)
real*8,                   intent(inout) :: q(surface_list%n_psi)
real*8,                   intent(inout) :: rad(surface_list%n_psi)

! --- Local variables
integer :: n_int
integer :: i_elm, j, k, n1, n2, n3
real*8  :: p, t,rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi, dA
real*8  :: RRgi, dRRgi_dr, dRRgi_ds, dRRgi_dp, ZZgi, dZZgi_dr, dZZgi_ds, dZZgi_dp, dRRgi_dt, dZZgi_dt
real*8  :: PSgi, dPSgi_dr, dPSgi_ds, PSI_R, PSI_Z, RZJAC, grad_psi, psi_n
real*8  :: Fgi,dFgi_dr,dFgi_ds,dFgi_drs,dFgi_drr,dFgi_dss
real*8  :: sum_dA, B_tot2, delta_phi
real*8  :: dRRgi_drs,dRRgi_drr,dRRgi_dss, dRRgi_drp,dRRgi_dsp,dRRgi_dpp
real*8  :: dZZgi_drs,dZZgi_drr,dZZgi_dss, dZZgi_drp,dZZgi_dsp,dZZgi_dpp, dPSgi_drs,dPSgi_drr,dPSgi_dss
integer :: i,m, ig1, ig2, i_plane
real*8  :: s_phi, c_phi, ndotB, cross_deriv(3), n_perp(3), B_boundary(3)
real*8  :: ndotB_max=0.0


write(*,*) "*********************************"
write(*,*) "*    Determine Boundary Flux    *"
write(*,*) "*********************************"

delta_phi = 2 * PI / float(n_plane+1) / float(n_period)
rad(:) = 0.d0
q(:)   = 0.d0

! Loop through surfaces
i=surface_list%n_psi
rad(i) = 0.d0
q(i)   = 0.d0
sum_dA = 0.d0
! Loop through pieces in poloidal plane
open(21, file='normal_points.dat')
do i_elm=(n_flux-2)*n_tht+1, (n_flux-1)*n_tht
  do i_plane=1,n_plane+1
    do ig1 = 1, 4
      si = 0.5 * (xgs(ig1) + 1.0)
      ri = 1.0

      do ig2 = 1, 4
        p = (i_plane - 1 + 0.5 * (xgs(ig2) + 1.0)) * delta_phi

        call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)

        call interp_RZP(node_list,element_list,i_elm,ri,si,p,   &
                        RRgi,dRRgi_dr,dRRgi_ds,dRRgi_dp,dRRgi_drs,dRRgi_drr,dRRgi_dss,dRRgi_drp, dRRgi_dsp, dRRgi_dpp, &
                        ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_dp,dZZgi_drs,dZZgi_drr,dZZgi_dss,dZZgi_drp, dZZgi_dsp, dZZgi_dpp)

        ! Radial component discarded because grid is already flux surface aligned
        dRRgi_dt = dRRgi_ds  ! + dRRgi_dr * dri
        dZZgi_dt = dZZgi_ds  ! + dZZgi_dr * dri
        
        ! Calculate normal to boundary
        n_perp = (/-dZZgi_dt, dRRgi_dt*dZZgi_dp-dRRgi_dp*dZZgi_dt, dRRgi_dt/)
        n_perp = n_perp * 1 / RRgi
        n_perp = n_perp / sqrt(sum(n_perp*n_perp))
        if ((i_plane .eq. 1) .and. (ig2 .eq. 1)) write(21,'(5ES16.8)') RRgi, ZZgi, n_perp

        ! Calculate surface area contribution from covariant components
        c_phi = cos(p)
        s_phi = sin(p)
        cross_deriv = (/(dRRgi_dt*s_phi*dZZgi_dp-(dRRgi_dp*s_phi+RRgi*c_phi)*dZZgi_dt),   &
                      (-dRRgi_dt*dZZgi_dp*c_phi+(dRRgi_dp*c_phi-RRgi*s_phi)*dZZgi_dt),  &
                      (dRRgi_dt * RRgi)  /)
        dA = sqrt(cross_deriv(1)*cross_deriv(1) + cross_deriv(2)*cross_deriv(2) + cross_deriv(3)*cross_deriv(3)) 
        !dA = sqrt(dRRgi_dt**2 + dZZgi_dt**2 + (dRRgi_dt * dZZgi_dp - dZZgi_dt * dRRgi_dp) ** 2)

        ! Calculate n.B
        RZjac  = DRRgi_dr * dZZgi_ds - dRRgi_ds * dZZgi_dr
        PSI_R = (   dPSgi_dr * dZZgi_ds - dPSgi_ds * dZZgi_dr ) / RZjac
        PSI_Z = ( - dPSgi_dr * dRRgi_ds + dPSgi_ds * dRRgi_dr ) / RZjac
        B_boundary = (/ -1/RRgi * Psi_Z, F0/RRgi, 1.0/RRgi * Psi_R /)
        ndotB = sum(n_perp*B_boundary)      
        ndotB_max = max(ndotB, ndotB_max)

        ! Factors of 0.5 for conversion of weights from -1 to 1 to weights between 0 and 1
        sum_dA = sum_dA +  wgs(ig1) * 0.5 * wgs(ig2) * 0.5 * delta_phi * dA * ndotB 
      end do
    end do
  end do
end do
close(21)

write(*,*)  "Max n.B: ", ndotB_max
write(*, *) "Total Boundary Flux: ", n_period * sum_dA, "Tm^2"

end subroutine determine_boundary_flux
