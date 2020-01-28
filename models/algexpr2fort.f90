program algexpr2fort
  use phys_module
  use mod_semianalytical
  use mod_equations
  use mod_parameters
  implicit none
  
  integer :: i, j, length, n, e
  integer, parameter :: line_length = 130
  
  character(3) :: model_num
  
  character(:),  allocatable       :: varname, full
  character(8),  dimension(n_rhs)  :: varname_rhs
  character(7),  dimension(n_amat) :: varname_amat
  character(12), dimension(n_aux)  :: varname_aux
  
  type(algexpr), dimension(n_rhs)  :: rhs
  type(algexpr), dimension(n_amat) :: amat
  type(algexpr), dimension(n_aux)  :: aux
  
  varname = "eq"
  
  call preset_parameters()
  time_evol_zeta = 0.
  time_evol_theta = 1.
  call init_equations()
  call get_rhs(rhs, varname_rhs)
  call get_amat(amat, varname_amat)
  if (n_aux .ne. 0) call get_aux(aux, varname_aux)
  
  write(model_num,'(I3.3)') jorek_model
  
  open(10, file="models/model"//model_num//"/rhs_unreadable.h", action="write", status="replace")
  do i=1,n_rhs
    full = varname_rhs(i) // "=" // gencode(rhs(i), varname)
    length = len(full)
    n = length/line_length
    if (n .gt. 1) then
      e = merge(n-1,n,n*line_length.eq.length)
      write(10,'(A,A)') full(1:line_length), "&"
      do j=2,e
        write(10,'(A,A,A)') "&", full((j-1)*line_length+1:j*line_length), "&"
      end do
      write(10,'(A,A)') "&", full(e*line_length+1:length)
    else
      write(10,'(A)') full
    end if
  end do
  close(10)
  
  open(20, file="models/model"//model_num//"/amat_unreadable.h", action="write", status="replace")
  do i=1,n_amat
    full = varname_amat(i) // " = " // gencode(amat(i), varname)
    length = len(full)
    n = length/line_length
    if (n .gt. 1) then
      e = merge(n-1,n,n*line_length.eq.length)
      write(20,'(A,A)') full(1:line_length), "&"
      do j=2,e
        write(20,'(A,A,A)') "&", full((j-1)*line_length+1:j*line_length), "&"
      end do
      write(20,'(A,A)') "&", full(e*line_length+1:length)
    else
      write(20,'(A)') full
    end if
  end do
  close(20)
  
  if (n_aux .ne. 0) then
    open(30, file="models/model"//model_num//"/aux_unreadable.h", action="write", status="replace")
    do i=1,n_aux
      full = varname_aux(i) // "=" // gencode(aux(i), varname)
      length = len(full)
      n = length/line_length
      if (n .gt. 1) then
        e = merge(n-1,n,n*line_length.eq.length)
        write(30,'(A,A)') full(1:line_length), "&"
        do j=2,e
          write(30,'(A,A,A)') "&", full((j-1)*line_length+1:j*line_length), "&"
        end do
        write(30,'(A,A)') "&", full(e*line_length+1:length)
      else
        write(30,'(A)') full
      end if
    end do
    close(30)
  end if
end program algexpr2fort
