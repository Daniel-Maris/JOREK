program algexpr2fort
#ifdef SEMIANALYTICAL
  use phys_module
  use mod_semianalytical
  use mod_equations
  use mod_parameters
  implicit none
  
  integer :: i
  integer, parameter :: line_length = 130
  
  character(3)  :: model_num
  
  character(:),  allocatable       :: varname, full
  character(8),  allocatable       :: varname_rhs(:)
  character(7),  allocatable       :: varname_amat(:)
  character(14), dimension(n_aux)  :: varname_aux
  
  type(algexpr), allocatable       :: rhs(:)
  type(algexpr), allocatable       :: amat(:)
  
  type(algexpr), dimension(n_aux)  :: aux
  
  varname = "eq"
  
  call preset_parameters()
  time_evol_zeta = 0.
  time_evol_theta = 1.
  call init_equations()
  !allocate(rhs(n_rhs), varname_rhs(n_rhs))
  !allocate(amat(n_amat), varname_amat(n_amat))
!  if ( with_TiTe ) then
!    allocate(rhs(n_rhs), varname_rhs(7))
!    allocate(amat(30), varname_amat(30))
!  else
!    allocate(rhs(6), varname_rhs(6))
!    allocate(amat(22), varname_amat(22))
!  end if

  call get_rhs(rhs, varname_rhs)
  call get_amat(amat, varname_amat)
  if (n_aux .ne. 0) call get_aux(aux, varname_aux)
  
  write(model_num,'(I3.3)') jorek_model
  
  open(10, file="models/model"//model_num//"/rhs_unreadable.h", action="write", status="replace")
  do i=1,n_rhs
    full = varname_rhs(i) // "=" // gencode(rhs(i), varname)
    call write_long_string(10,full)
  end do
  close(10)
  
  open(20, file="models/model"//model_num//"/amat_unreadable.h", action="write", status="replace")
  do i=1,n_amat
    full = varname_amat(i) // " = " // gencode(amat(i), varname)
    call write_long_string(20,full)
  end do
  close(20)
  
  if (n_aux .ne. 0) then
    open(30, file="models/model"//model_num//"/aux_unreadable.h", action="write", status="replace")
    do i=1,n_aux
      full = varname_aux(i) // "=" // gencode(aux(i), varname)
      call write_long_string(30,full)
    end do
    close(30)
  end if
  
  contains
  
  subroutine write_long_string(id, str)
    implicit none
    integer, intent(in) :: id
    character(:), allocatable, intent(in) :: str
    integer :: j, length, n, e
    
    length = len(str)
    n = length/line_length
    if (n .gt. 1) then
      e = merge(n-1,n,n*line_length.eq.length)
      write(id,'(A,A)') str(1:line_length), "&"
      do j=2,e
        write(id,'(A,A,A)') "&", str((j-1)*line_length+1:j*line_length), "&"
      end do
      write(id,'(A,A)') "&", str(e*line_length+1:length)
    else
      write(id,'(A)') str
    end if
  end subroutine write_long_string
#else
  write(*,*) ">> Code generation is only for semianalytical models. <<"
#endif
end program algexpr2fort
