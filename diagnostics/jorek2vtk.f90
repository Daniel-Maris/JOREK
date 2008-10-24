!**********************************************************************
!* program to convert a JOREK2 restart file into binary VTK format    *
!**********************************************************************

program jorek2vtk
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
integer,allocatable   :: ien (:,:)
integer               :: i, j, k, m, etype, ivtk, irst, int, i_var, i_tor, index, index_node
character             :: buffer*80, lf*1, str1*12, str2*12
character*8, allocatable :: scalar_names(:), vector_names(:)
real*4                :: float
real*8                :: s, t
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt, ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt, W,W_s,W_t,W_st,W_ss,W_tt
real*8                :: U,U_s,U_t,U_st,U_ss,U_tt, RHO,RH_s,RH_t,RH_st,RH_ss,RH_tt, TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: u0_x, u0_y, xjac, v_perp, Psi_J, R_p, error, zj_x, zj_y, ps_x, ps_y
real*8                :: V, V_s, V_t, V_st, V_ss, V_tt

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

read(5,in1)               ! read the namelist input file

ivtk = 22                 ! an arbitrary unit number for the VTK output file

nsub  = 11                 ! the number of subdivisions of the cubic finite elements into linear pieces
i_tor = -1

n_scalars = n_var         ! number of scalars to write to the VTK output file
n_vectors = 3

allocate(scalar_names(n_scalars), vector_names(n_vectors))

scalar_names = (/ 'flux    ','U       ','current ','W       ','density ','T       ','Vpar     '/)
vector_names = (/ 'v_perp  ','v_par   ','V_tot   '/)

call import_restart(node_list,element_list)

nnos = nsub*nsub*node_list%n_nodes
allocate(xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))

nnoel = 4
nel   = (nsub-1)*(nsub-1)*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0

do i=1,element_list%n_elements

  do j=1,nsub
    s = float(j-1)/float(nsub-1)
    do k=1,nsub
      t = float(k-1)/float(nsub-1)

      call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

      xjac  = R_s * Z_t - R_t * Z_s

      inode = inode+1

      xyz(1:3,inode) = (/ R, Z, 0.d0/)

      if ((i_tor .ge. 1) .and. (i_tor .le. n_tor)) then

        do m=1,n_var
          call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
          scalars(inode,m) = P
        enddo

        if ((xjac .gt. 1.d-6)) then

          call interp(node_list,element_list,i,1,i_tor,s,t,Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)
          call interp(node_list,element_list,i,2,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
          call interp(node_list,element_list,i,3,i_tor,s,t,ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt)
          call interp(node_list,element_list,i,4,i_tor,s,t,W,W_s,W_t,W_st,W_ss,W_tt)
          call interp(node_list,element_list,i,5,i_tor,s,t,RHO,RH_s,RH_t,RH_st,RH_ss,RH_tt)
          call interp(node_list,element_list,i,6,i_tor,s,t,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
          call interp(node_list,element_list,i,7,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)

	  u0_x  = (   Z_t * U_s - Z_s * U_t ) / xjac
          u0_y  = ( - R_t * U_s + R_s * U_t ) / xjac

          ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac
          ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac

          zj_x  = (   Z_t * ZJ_s - Z_s * ZJ_t ) / xjac
          zj_y  = ( - R_t * ZJ_s + R_s * ZJ_t ) / xjac

	  v_perp = R * sqrt(u0_x*u0_x + u0_y * u0_y)

	  vectors(inode,:,1) = (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
          vectors(inode,:,2) = (/ + ps_y /R * V, - ps_x /R * V, F0/R * V /)
          vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * V, + R * u0_x - ps_x /R * V, F0/R * V /)

          psi_J = (Ps_s * ZJ_t - PS_t * ZJ_s ) / xjac
          R_p   = (2.d0 * R * (R_s * (RH_t * TT + RHO * TT_t) - R_t * (RH_s * TT + RHO * TT_s) )) / xjac
          error = psi_J - R_p

        endif

      else

        do i_tor = 1,n_tor

	  do m=1,n_var
            call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,m) = scalars(inode,m) + P * HZ(i_tor,1)
          enddo

          if ((xjac .gt. 1.d-6)) then

            call interp(node_list,element_list,i,1,i_tor,s,t,Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)
            call interp(node_list,element_list,i,2,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
            call interp(node_list,element_list,i,3,i_tor,s,t,ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt)
            call interp(node_list,element_list,i,4,i_tor,s,t,W,W_s,W_t,W_st,W_ss,W_tt)
            call interp(node_list,element_list,i,5,i_tor,s,t,RHO,RH_s,RH_t,RH_st,RH_ss,RH_tt)
            call interp(node_list,element_list,i,6,i_tor,s,t,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
            call interp(node_list,element_list,i,7,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)

            u0_x  = (   Z_t * U_s - Z_s * U_t ) / xjac * HZ(i_tor,1)
            u0_y  = ( - R_t * U_s + R_s * U_t ) / xjac * HZ(i_tor,1)

            ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac * HZ(i_tor,1)
            ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac * HZ(i_tor,1)

            zj_x  = (   Z_t * ZJ_s - Z_s * ZJ_t ) / xjac * HZ(i_tor,1)
            zj_y  = ( - R_t * ZJ_s + R_s * ZJ_t ) / xjac * HZ(i_tor,1)

	    v_perp = R * sqrt(u0_x*u0_x + u0_y * u0_y)

!            vectors(inode,:,1) = vectors(inode,:,1) + (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
!            vectors(inode,:,2) = vectors(inode,:,2) + (/ + ps_y /R * V, - ps_x /R * V, F0/R * V /)
!            vectors(inode,:,3) = vectors(inode,:,3) + (/ - R * u0_y + ps_y /R * V, + R * u0_x - ps_x /R * V, F0/R * V /)

            vectors(inode,:,1) = vectors(inode,:,1) + (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
            vectors(inode,:,2) = vectors(inode,:,2) + (/ + ps_y /R * V, - ps_x /R * V, 0.d0 /)
            vectors(inode,:,3) = vectors(inode,:,3) + (/ - R * u0_y + ps_y /R * V, + R * u0_x - ps_x /R * V, 0.d0 /)

	  endif

	enddo

      endif

    enddo
  enddo

  do j=1,nsub-1
    do k=1,nsub-1
      ielm        = ielm+1
      ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
      ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
      ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
      ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
    enddo
  enddo

enddo


!--------------------------------------------------- write the binary VTK file
etype = 9  ! for vtk_quad

lf = char(10) ! line feed character

open(unit=ivtk,file='jorek_tmp.vtk',form='binary',convert='BIG_ENDIAN')

buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:12),'(i12)') nnos
buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) ((xyz(i,j),i=1,3),j=1,nnos)

! CELLS SECTION
write(str1(1:12),'(i12)') nel            ! number of elements (cells)
write(str2(1:12),'(i12)') nel*(1+nnoel)  ! size of the following element list (nel*(nnoel+1))
buffer = lf//lf//'CELLS '//str1//' '//str2//lf                                        ; write(ivtk) trim(buffer)
write(ivtk) (nnoel,(ien(i,j),i=1,nnoel),j=1,nel)

! CELL_TYPES SECTION
write(str1(1:12),'(i12)') nel   ! number of elements (cells)
buffer = lf//lf//'CELL_TYPES'//str1//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) (etype,i=1,nel)

! POINT_DATA SECTION
write(str1(1:12),'(i12)') nnos
buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
  write(ivtk) (scalars(i,i_var),i=1,nnos)
enddo

do i_var =1, n_vectors
  buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf                                          ; write(ivtk) trim(buffer)
  write(ivtk) ((vectors(j,i,i_var),i=1,3),j=1,nnos)
enddo

close(ivtk)

end

