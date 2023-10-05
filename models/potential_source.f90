subroutine potential_source(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,potential_profile, &
                       dPhi_dpsi,dPhi_dz,dPhi_dpsi2,dPhi_dz2,dPhi_dpsi_dz,dPhi_dpsi3,dPhi_dpsi_dz2, dPhi_dpsi2_dz)
use phys_module
use vacuum, only: current_FB_fact

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: potential_profile, dPhi_dpsi, dPhi_dz, dPhi_dpsi2, dPhi_dz2, &
                        dPhi_dpsi_dz, dPhi_dpsi3, dPhi_dpsi_dz2, dPhi_dpsi2_dz

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, sig_T, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, factor
real*8  :: cosh1, cosh2, cosh3, cosh3_u
real*8  :: tanh1, tanh2, tanh2_u
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi

psi_n = max( min(psi_n, 2.), 0. )


!!!!!
if (.not. num_Phi) then
        write(*,*) "Potential source must be numeric"
        stop 

else
        left = 1
        right = num_Phi_len
        do
                if(right == left + 1) exit
                mid = left+right/2
                if (num_Phi_x(mid) >= psi_n)
                        right = mid
                else
                        left = mid
                end if                    
        end do        
        !weights
        aux1 = (psi_n - num_Phi_x(left))/(num_Phi_x(right) - num_Phi_x(left)) !d from desired point to left point, normalized2 step size
        aux2 = (1. - aux1)
        prof1 = num_Phi_y0(left)*aux2 + num_Phi_y0(right)*aux1 !Phi at the desired point
        dprof1_dpsi = (num_Phi_y1(left)*aux2 - num_Phi_y1(right)*aux1) / delta_psi !first derivative at desired point
        dprof2_dpsi = (num_Phi_y2(left)*aux2 - num_Phi_y2(right)*aux1) / delta_psi**2
        dprof3_dpsi = (num_Phi_y3(left)*aux2 - num_Phi_y3(right)*aux1) / delta_psi**3 

end if
!!!!!!!!!!!!!!!!!
write(*,*)   "Potential source is currently not implemented!"
stop

return
end subroutine potential_source
