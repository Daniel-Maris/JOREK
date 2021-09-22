program recalc_egr
  use data_structure
  use nodes_elements
  use basis_at_gaussian
  use mod_import_restart
  use mod_expression
!  use mod_integrals3D_nompi
  use phys_module
  use mod_chi
  implicit none
  
  integer, parameter :: ts1 = 52
  integer, parameter :: nts = 20
  
  character*14 :: filein
  integer      :: ts, ierr, i_mag, i_kin, imf
  
  real*8, dimension(:), allocatable              :: res
  real*8, dimension(nts)                         :: time
  real*8, dimension(1+int(n_coord_period/2),nts) :: Wmag, Wkin
  
  call det_modes
  call initialise_basis
  call init_chi_basis
  call initialise_parameters(0,  "__NO_FILENAME__")
  call init_expr
  
  allocate(res(exprs_all_int%n_expr+1))
  res = 0.d0
  
  index_now = ts1
  
  i_mag = get_expr_num_int(exprs_all_int,'Wmag_tot')
  i_kin = get_expr_num_int(exprs_all_int,'Kin_perp_tot')
  
  do ts=1,nts
    write(filein,'(A,i5.5)') 'jorek', index_now
    call import_restart(node_list,element_list, filein, rst_format, ierr, .true.)
    call energy3D(node_list,element_list,Wmag(:,ts),Wkin(:,ts))
!    call int3d_new(0, node_list, element_list, bnd_node_list, bnd_elm_list, exprs_all_int, res, 1)
    time(ts) = xtime(index_now)
!    Wmag(ts) = res(i_mag+1)
!    Wkin(ts) = res(i_kin+1)
    index_now = index_now + 1
  end do
  
  open(21,file="energies_growth_rates.dat",action="write",status="replace")
  do ts=1,nts
    write(21,'(12E18.8)',advance='no') time(ts), Wmag(:,ts), Wkin(:,ts)
    if (ts .gt. 1) then
      do imf=1,1+int(n_coord_period/2)
        if (Wmag(imf,ts) .gt. 0. .and. Wmag(imf,ts-1) .gt. 0.) then
          write(21,'(E18.8)',advance='no') 0.5d0*(log(Wmag(imf,ts)) - log(Wmag(imf,ts-1)))/(time(ts) - time(ts-1))
        else
          write(21,'(E18.8)',advance='no') 0.d0
        end if
        if (Wkin(imf,ts) .gt. 0. .and. Wkin(imf,ts-1) .gt. 0.) then
          write(21,'(E18.8)',advance='no') 0.5d0*(log(Wkin(imf,ts)) - log(Wkin(imf,ts-1)))/(time(ts) - time(ts-1))
        else
          write(21,'(E18.8)',advance='no') 0.d0
        end if
      end do
      write(21,*)
    else
      write(21,'(12E18.8)') (0.d0,imf=1,2+2*int(n_coord_period/2))
    end if
  end do
  close(21)
end program recalc_egr