program algexpr2fort
#ifdef SEMIANALYTICAL
  use phys_module
  use mod_semianalytical
  use mod_equations
  use mod_parameters
  use mod_equations
  implicit none
  
  integer :: i_var, j_var, i_aux
  integer, parameter :: line_length = 130
  
  character(3)  :: model_num
  
  character(:),  allocatable       :: varname, full
  character(8),  dimension(n_var)  :: index_names(n_var)
  character(19), dimension(n_aux)  :: varname_aux
  
  type(algexpr), dimension(n_aux)  :: aux
  
  varname = "eq"
  
  call preset_parameters()
  time_evol_zeta = 0.
  time_evol_theta = 1.
  call init_equations()
  call get_varnames(index_names)
  if (n_aux .ne. 0) call get_aux(aux, varname_aux)
  
  write(model_num,'(I3.3)') jorek_model
  
  ! Write out RHS for included variables
  open(10, file="models/model"//model_num//"/rhs_unreadable.h", action="write", status="replace")
  do i_var=1,n_var
    if ((.not. associated(rhs_semianalytic(i_var)%operand1)) .and. (.not. associated(rhs_semianalytic(i_var)%operand2))) then
      cycle
    else
      full = "rhs_ij("// index_names(i_var) // ") = " // gencode(rhs_semianalytic(i_var), varname)
      call write_long_string(10,full)
    endif
  end do
  close(10)
  
  ! Write out AMAT for included variables
  open(20, file="models/model"//model_num//"/amat_unreadable.h", action="write", status="replace")
  do i_var = 1, n_var
    do j_var = 1, n_var
      if ((.not. associated(amat_semianalytic(i_var, j_var)%operand1)) .and. (.not. associated(amat_semianalytic(i_var, j_var)%operand2))) then
        cycle
      else
        full = "amat_ij("// index_names(i_var) // "," // index_names(j_var) // ") = " // gencode(amat_semianalytic(i_var, j_var), varname) 
        call write_long_string(20,full)
      endif 
    enddo
  enddo
  close(20)
  
  if (n_aux .ne. 0) then
    open(30, file="models/model"//model_num//"/aux_unreadable.h", action="write", status="replace")
    do i_aux=1,n_aux
      full = varname_aux(i_aux) // "=" // gencode(aux(i_aux), varname)
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
