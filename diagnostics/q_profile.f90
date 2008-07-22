subroutine q_profile(node_list,element_list,surface_list)
!----------------------------------------------------------------------
! subroutine calculats the q-profile from the flux surface representation
! (adapted from helena20)
!---------------------------------------------------------------------
use data_structure
use phys_module

implicit none

!--------------------------------------- gaussian points between (-1.,1.)
real*8 :: xgs(4), wgs(4)
data xgs /-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /
data wgs / 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer :: n_int
integer :: i_elm, j, k, n1, n2, n3
real*8  :: t,rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi, dl
real*8  :: RRgi, dRRgi_dr, dRRgi_ds, ZZgi, dZZgi_dr, dZZgi_ds, dRRgi_dt, dZZgi_dt
real*8  :: PSgi, dPSgi_dr, dPSgi_ds, PSI_R, PSI_Z, RZJAC, grad_psi, psi_n
real*8  :: sum_dl, B_tot2, PI
real*8  :: dRRgi_drs,dRRgi_drr,dRRgi_dss, dZZgi_drs,dZZgi_drr,dZZgi_dss, dPSgi_drs,dPSgi_drr,dPSgi_dss
real*8, allocatable :: q(:)
integer :: i,m, ig, ip

write(*,*) '**********************************'
write(*,*) '*        q-profile               *'
write(*,*) '**********************************'

allocate(q(surface_list%n_psi))

PI = 2.d0 *asin(1.d0)

do i=2, surface_list%n_psi

  q(i)   = 0.d0
  sum_dl = 0.d0

  do k=1, surface_list%flux_surfaces(i)%n_pieces

    do ig = 1, 4

      t = xgs(ig)

      rr1  = surface_list%flux_surfaces(i)%s(1,k)
      drr1 = surface_list%flux_surfaces(i)%s(2,k)
      rr2  = surface_list%flux_surfaces(i)%s(3,k)
      drr2 = surface_list%flux_surfaces(i)%s(4,k)

      ss1  = surface_list%flux_surfaces(i)%t(1,k)
      dss1 = surface_list%flux_surfaces(i)%t(2,k)
      ss2  = surface_list%flux_surfaces(i)%t(3,k)
      dss2 = surface_list%flux_surfaces(i)%t(4,k)

      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

      i_elm = surface_list%flux_surfaces(i)%elm(k)

      call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)

      call interp_RZ(node_list,element_list,i_elm,ri,si,RRgi,dRRgi_dr,dRRgi_ds,dRRgi_drs,dRRgi_drr,dRRgi_dss, &
                                                        ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_drs,dZZgi_drr,dZZgi_dss)

      dRRgi_dt = dRRgi_dr * dri + dRRgi_ds * dsi
      dZZgi_dt = dZZgi_dr * dri + dZZgi_ds * dsi

      dl = sqrt(dRRgi_dt**2 + dZZgi_dt**2)

      RZjac  = DRRgi_dr * dZZgi_ds - dRRgi_ds * dZZgi_dr

      PSI_R = (   dPSgi_dr * dZZgi_ds - dPSgi_ds * dZZgi_dr ) / RZjac
      PSI_Z = ( - dPSgi_dr * dRRgi_ds + dPSgi_ds * dRRgi_dr ) / RZjac

      grad_psi = sqrt(PSI_R * PSI_R + PSI_Z * PSI_Z)

      B_tot2 =  (F0 / RRgi)**2 + (grad_psi/ RRgi)**2

      sum_dl = sum_dl +  wgs(ig) * dl

      q(i) = q(i) +  wgs(ig) / (RRgi * grad_psi) * dl

    enddo

  enddo

  q(i) = F0 * q(i) / (2.d0 * PI)

  write(*,'(A,i5,3e16.8)') ' psi, q : ',i,surface_list%psi_values(i),q(i),sum_dl

enddo


!----------------------------------- values on axis
!q(1) =  PI / sqrt(CRR_axis*CZZ_Axis)

call lplot(2,1,1,surface_list%psi_values(2),q(2),surface_list%n_psi-1,1,'q-profile',9,'flux',4,'q',1)

deallocate(q)

return
end