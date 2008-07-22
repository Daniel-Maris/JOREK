subroutine current(xpoint2,R,Z,Z_xpoint,psi,psi_axis,psi_bnd,zjz)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: R, Z, Z_xpoint, psi, psi_axis, psi_bnd, zj0, zjz, dj_dpsi, psi_n, sigz
real*8  :: zjz_p, zjz_q, dj_dz, ss, dj_ds
real*8  :: zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8  :: zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
real*8  :: zFFprime, dFFprime_dpsi, dFFprime_dz, dFFprime_dpsi_dz,dFFprime_dpsi2,dFFprime_dz2
real*8  :: dc(4), ddc(4), bar_start, bar_mid, bar_end

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

call density(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,&
             zn,dn_dpsi,dn_dz,dn_dpsi2, dn_dz2, dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)

call temperature(xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                 zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)

call FFprime(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
             zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)

zjz   = zFFprime - R*R * (zn * dT_dpsi + dn_dpsi * zT)

return
end
