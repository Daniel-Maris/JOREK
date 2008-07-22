subroutine Define_Boundary
!---------------------------------------------------------------------
! subroutine to define the spline coefficients for the shape of the
! plasma boundary
!   ELLIP : ellipticity
!   TRIA  : triangularity
! output is contained in the module boundary
!---------------------------------------------------------------------
use phys_module

implicit none

integer             :: n_bnd, i, j, m
real*8, allocatable :: r_bnd(:),psi_bnd(:),dr_bnd(:),dpsi_bnd(:),tht_bnd(:)
real*8, allocatable :: Work(:)
real*8              :: PI, Vr(4), Vpsi(4), RP, ZP, theta, tht_i

write(*,*) '*******************************************'
write(*,*) '*    Defining boundary                    *'
write(*,*) '*******************************************'

PI = 2.d0*asin(1.d0)

!------------------------------- boundary given by ellip, tria etc as splined in r_bnd, psi_bnd (module boundary)
if (mf .le. 0) then

  write(*,'(A18,f8.4)')  '  ellipticity     : ',ellip
  write(*,'(A18,2f8.4)') '  triangularity   : ',tria_u,tria_l
  write(*,'(A18,2f8.4)') '  quadrangularity : ',quad_u,quad_l
  write(*,'(A18,f8.4)')  '  Xpoint(ampl)    : ',xampl
  write(*,'(A18,f8.4)')  '  Xpoint(width)   : ',xwidth
  write(*,'(A18,f8.4)')  '  Xpoint(sigma)   : ',xsig
  write(*,'(A18,f8.4)')  '  Xpoint(angle)   : ',xtheta

  n_bnd = 256
  allocate(tht_bnd(n_bnd),r_bnd(n_bnd),dr_bnd(n_bnd),psi_bnd(n_bnd),dpsi_bnd(n_bnd))

  do i=1,n_bnd

    theta = 2.d0*PI * real(i-1)/real(n_bnd-1)

    if (theta .lt. pi) then

      RP = cos(theta + tria_u*sin(theta) + quad_u*sin(2.d0*theta))

    else

      RP = cos(theta + tria_l*sin(theta) + quad_l*sin(2.d0*theta))

    endif

    ZP = ellip * sin(theta)

    tht_bnd(i) = atan2(ZP,RP)
    if (tht_bnd(i) .lt. 0.d0*pi)    tht_bnd(i) = tht_bnd(i) + 2.d0*pi

    tht_i = tht_bnd(i)
    if (theta      .lt. 0.5d0*pi)   tht_i = tht_i + 2.d0*pi

    r_bnd(i)   = sqrt(rp**2+zp**2)
    psi_bnd(i) =  - xshift * sin(tht_i) + xleft * cos(tht_i) &
               + xampl*(-1.d0 + (xwidth*(tht_i-xtheta)/xsig)**2)* exp( - ((tht_i-xtheta)/xsig)**2)

  enddo

  r_bnd(n_bnd)   = r_bnd(1)
  psi_bnd(n_bnd) = psi_bnd(1)

  allocate(work(3*n_bnd))

  call TB15A(n_bnd,tht_bnd,r_bnd,dr_bnd,work,6)           ! periodic spline of the radius
  call TB15A(n_bnd,tht_bnd,psi_bnd,dpsi_bnd,work,6)       ! periodic spline of flux

  call lplot6(1,1,tht_bnd,psi_bnd,n_bnd,'psi at boundary')

  deallocate(work)

  mf = 256

  do j=1, mf

    theta = 2.d0 * PI * float(j-1)/float(mf)

    call  TG02A(-1,n_bnd,tht_bnd,r_bnd,dr_bnd,theta,Vr)
    call  TG02A(-1,n_bnd,tht_bnd,psi_bnd,dpsi_bnd,theta,Vpsi)

    fbnd(j) = Vr(1)
    fpsi(j) = Vpsi(1)

  enddo

  call rft2(fbnd,mf,1)
  call rft2(fpsi,mf,1)

  do m=1,mf
    fbnd(m) = 2.d0 * fbnd(m) / float(mf)
    fpsi(m) = 2.d0 * fpsi(m) / float(mf)
  enddo
  do m=2,mf,2
    fbnd(m) = - fbnd(m)
    fpsi(m) = - fpsi(m)
  enddo

else
  write(*,*) ' boundary defined by Fourier series : ',mf
endif

return

end