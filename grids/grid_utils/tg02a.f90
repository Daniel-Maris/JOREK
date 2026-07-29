subroutine tg02a(ix,n,x0,s,d,x,v)
  
  implicit none
  !------------------------------------------------------------------
  ! Routine for interpolation of splines and their derivative
  !------------------------------------------------------------------
  !    n  : number of points
  !    ix : negative 0 -> no initial guess for where xi is
  !        positive -> gues for index close to value x
  !   x0  : the coordinates of the spline points
  !   f   : the function values of the spline points
  !   d   : the derivatives on the spline points
  !   x   : the coordinate where the output is wanted
  !   v(1-4) : value and derivatives of the spline interpolation
  !------------------------------------------------------------------

  real*8, intent(in)   :: x
  integer, intent(in)  :: ix,n
  real*8, intent(inout)   :: d(*),s(*),x0(*),v(4)

  ! Local variables
  integer :: j
  logical :: found
  real*8  :: h00, h10, h01, h11, theta, theta2, theta3, t, eps, dx

  found = .false.
  eps = 1.d-30

  if (x .le. x0(1)) then
     if (dabs(x-x0(1)) .lt. eps*dabs(x0(n)-x0(1))) then
        j = 1
     else
        v(:) = 0.d0
        return
     end if
  else if (x .ge. x0(n)) then
     if (dabs(x-x0(n)) .lt. eps * dabs(x0(n)-x0(1))) then
        j = n-1
     else
        v(:) = 0.d0
        return
     end if
  else 
     if (ix .gt. 0) j = max(min(ix,n-1),1)
     if (ix .le. 0) j = (x-x0(1))/(x0(n)-x0(1))*(n-1)+1 ! initial guess of value
     
     do while  (.not. found)
        if (x .gt. x0(j+1)) then
           j = j + 1
        else if (x .ge. x0(j)) then
           found = .true.
        else
           j = j - 1
        end if
     end do
  end if

  dx = x0(j+1) - x0(j)
  theta = (x-x0(j))/dx
  theta2 = theta*theta
  theta3 = theta2*theta


  h00 = 2.d0*theta3 - 3*theta2 + 1
  h10 = theta3 - 2.d0*theta2 + theta
  h01 = -2.d0*theta3 + 3*theta2
  h11 = theta3 - theta2

  ! Calculate spline values and their derivative
  v(1) = h00*s(j) + h10*dx*d(j) + h01*s(j+1) + h11*dx*d(j+1)
  v(2) = 6.d0*(theta2 - theta)*s(j) + (3.d0*theta2 - 4.d0*theta + 1)*dx*d(j) + 6.d0*(-theta2 + theta)*s(j+1) + (3.d0*theta2 - 2.d0*theta)*dx*d(j+1)
  v(2) = v(2)/dx
  v(3) = (12.d0*theta - 6.d0)*s(j) + (6.d0*theta - 4.d0)*dx*d(j) + (-12.d0*theta + 6.d0)*s(j+1) + (6.d0*theta - 2.d0)*dx*d(j+1)
  v(3) = v(3)/(dx**2)
  v(4) = 6.d0*(2*s(j) - 2.d0*s(j+1) + dx*d(j) + dx*d(j+1))/(dx**3)


end subroutine tg02a
