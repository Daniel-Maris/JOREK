subroutine tg02a(ix,n,u,s,d,x,v)
  implicit none
  !------------------------------------------------------------------
  !    n  : number of points
  !    ix : negative 0 -> no initial guess for where xi is
  !        positive -> gues for index close to value x
  !   u   : the coordinates of the spline points
  !   s   : the function values of the spline points
  !   d   : the derivatives on the spline points
  !   x   : the coordinate where the output is wanted
  !   v(1-4) : value and derivatives of the spline interpolation
  !------------------------------------------------------------------
  real*8, intent(in)   :: x
  integer, intent(in)  :: ix,n
  real*8, intent(inout)   :: d(*),s(*),u(*),v(4)

  integer :: j
  logical :: found
  real*8  :: h00, h10, h01, h11, theta, theta2, theta3, t, eps,h
  found = .false.
  eps = 1.d-30
  if (x .le. u(1)) then
     if (dabs(x-u(1)) .lt. eps*dabs(u(n)-u(1))) then
        j = 1
     else
        v(:) = 0.d0
        return
     end if

  else if (x .ge. u(n)) then
     if (dabs(x-u(n)) .lt. eps * dabs(u(n)-u(1))) then
        j = n-1
     else
        v(:) = 0.d0
        return
     end if
  else 
     if (ix .gt. 0.d0) j = max(min(ix,n-1),1)
     if (ix .lt. 0.d0) j = (x-u(1))/(u(n)-u(1))*(n-1)+1 ! initial guess of value

     do while  (.not. found)
        if (x .gt. u(j+1)) then
           j = j + 1
        else if (x .ge. u(j)) then
           found = .true.
        else
           j = j - 1
        end if
     end do
  end if

  h = u(j+1) - u(j)
  theta = (x-u(j))/h
  theta2 = theta*theta
  theta3 = theta2*theta


  h00 = 2.d0*theta3 - 3*theta2 + 1
  h10 = theta3 - 2.d0*theta2 + theta
  h01 = -2.d0*theta3 + 3*theta2
  h11 = theta3 - theta2

  ! Calculate spline values and their derivative
  v(1) = h00*s(j) + h10*h*d(j) + h01*s(j+1) + h11*h*d(j+1)
  v(2) = 6.d0*(theta2 - theta)*s(j) + (3.d0*theta2 - 4.d0*theta + 1)*h*d(j) + 6.d0*(-theta2 + theta)*s(j+1) + (3.d0*theta2 - 2.d0*theta)*h*d(j+1)
  v(2) = v(2)/h
  v(3) = (12.d0*theta - 6.d0)*s(j) + (6.d0*theta - 4.d0)*h*d(j) + (-12.d0*theta + 6.d0)*s(j+1) + (6.d0*theta - 2.d0)*h*d(j+1)
  v(3) = v(3)/(h**2)
  v(4) = 6.d0*(2*s(j) - 2.d0*s(j+1) + h*d(j) + h*d(j+1))/(h**3)


end subroutine tg02a
