subroutine tb15c(n,x,f,d,w,lp)
  implicit none
  !------------------------------------------------------------------
  ! hsl routine for cubic spline with periodic boundary conditions
  ! first point must be the same as last : f(1)=f(n)
  !    n : number of points
  !    x : coordinate (input)
  !    f : the function values to be splined (input)
  !    d : the derivatives at the points (output)
  !    w : workspace (dimension 3n)
  !   lp : unit number for output
  !------------------------------------------------------------------
  real*8     :: zero,one,two,three
  parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
  integer    :: lp,n
  real*8     :: d(n),f(n),w(*),x(n)
  real*8     :: a3n1,f1,f2,h1,h2,p
  integer    :: i,j,k,n2

  ! Local arrays for tridiagonal matrix coefficients and right-hand side
  real*8 :: a(n),b(n),c(n),rhs(n)

  ! Compute coefficients for the tridiagonal system
  do i = 1, n
    h1 = x(i) - x(i-1)  ! Step size between current and previous point
    if (i == n) then
      h2 = x(2) - x(1)  ! Step size for periodic boundary condition
    else
      h2 =  x(i+1) - x(i)  ! Step size between current and next point
    endif 
    a(i) = 1/h1  ! Lower diagonal coefficient
    b(i) = 2*(1/h1 + 1/h2)  ! Main diagonal coefficient
    c(i) = 1/h2  ! Upper diagonal coefficient
    if (i == 1) then
      rhs(i) = 3*((f(2) - f(1))/(h2*h2) + (f(n) - f(n-1))/(h1*h1))  ! Right-hand side for first point
    elseif (i == n) then
      rhs(i) = 3*((f(2) - f(1))/(h1*h1) + (f(n) - f(n-1))/(h2*h2))  ! Right-hand side for last point
    else
      rhs(i) = 3*((f(i+1) - f(i))/(h2*h2) + (f(i) - f(i-1))/(h1*h1))  ! Right-hand side for other points
    end if
  end do

  ! Solve the cyclic tridiagonal system
  call cyclic(a,b,c,a(1),c(n),rhs,x,n)
  
end subroutine


subroutine tridag(a, b, c, r, u, n)
  implicit none
  !------------------------------------------------------------------
  ! Solves a tridiagonal system of equations using the Thomas algorithm
  !    a : lower diagonal (input)
  !    b : main diagonal (input)
  !    c : upper diagonal (input)
  !    r : right-hand side (input)
  !    u : solution vector (output)
  !    n : number of equations
  !------------------------------------------------------------------
  integer, intent(in) :: n
  real*8, dimension(n), intent(in) :: a, b, c, r
  real*8, dimension(n), intent(out) :: u
  integer :: j
  real*8 :: bet
  real*8, dimension(:), allocatable :: gam

  ! Allocate workspace vector
  allocate(gam(n))

  ! Check for invalid input
  if (b(1) == 0.0) then
    print *, "tridag: rewrite equations"
    stop
  end if

  ! Forward substitution
  bet = b(1)
  u(1) = r(1) / bet
  do j = 2, n
    gam(j) = c(j - 1) / bet  ! Compute gamma
    bet = b(j) - a(j) * gam(j)  ! Update beta
    if (bet == 0.0) then
      print *, "tridag failed"
      deallocate(gam)
      stop
    end if
    u(j) = (r(j) - a(j) * u(j - 1)) / bet  ! Update solution
  end do

  ! Back substitution
  do j = n - 1, 1, -1
    u(j) = u(j) - gam(j + 1) * u(j + 1)  ! Update solution
  end do

  ! Deallocate workspace vector
  deallocate(gam)

  return
end subroutine tridag
 

subroutine cyclic(a,b,c,alpha,beta,r,x,n)
  !------------------------------------------------------------------
  ! Solves a cyclic tridiagonal system of equations
  !    a : lower diagonal (input)
  !    b : main diagonal (input)
  !    c : upper diagonal (input)
  ! alpha : coefficient for first equation (input)
  !  beta : coefficient for last equation (input)
  !    r : right-hand side (input)
  !    x : solution vector (output)
  !    n : number of equations
  !------------------------------------------------------------------
  integer n,nmax
  real*8:: alpha,beta,a(n),b(n),c(n),r(n),x(n)
  parameter (nmax=500)

  integer i
  real*8:: fact,gamma,bb(nmax),u(nmax),z(nmax)

  ! Modify the main diagonal for the cyclic system
  gamma=-b(1) 
  bb(1)=b(1)-gamma 
  bb(n)=b(n)-alpha*beta/gamma
  do  i=2,n-1
    bb(i)=b(i)
  enddo 

  ! Solve the modified tridiagonal system
  call tridag(a,bb,c,r,x,n) 

  ! Solve for the correction vector
  u(1)=gamma 
  u(n)=alpha
  do  i=2,n-1
    u(i)=0.
  enddo 
  call tridag(a,bb,c,u,z,n) 

  ! Apply the correction
  fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
  do  i=1,n 
    x(i)=x(i)-fact*z(i)
  enddo 

  return
end subroutine cyclic

subroutine tb15c(n,x,f,d,w,lp)

  implicit none
  
  integer,intent(in) :: n,lp
  real(8),intent(in) :: x(n),f(n)
  real(8),intent(out):: d(n)
  real(8),intent(inout):: w(*)
  
  real(8) :: a(n-1),b(n-1),c(n-1),rhs(n-1)
  real(8) :: sol(n-1)
  
  real(8) :: h1,h2,f1,f2
  integer :: i,m
  
  
  if (n < 4) then
     write(lp,*) "n too small"
     return
  endif
  
  
  m=n-1
  
  
  !------------------------------------------------------
  ! Build equations for unknowns:
  !
  ! sol(1)=d(2)
  ! ...
  ! sol(n-1)=d(n)
  !
  !------------------------------------------------------
  
 do i=2,n

    if(i==n) then
        h2 = 1.d0/(x(2)-x(1))
        f2 = f(2)
    else
        h2 = 1.d0/(x(i+1)-x(i))
        f2 = f(i+1)
    endif


    h1 = 1.d0/(x(i)-x(i-1))

    f1=f(i-1)


    a(i-1)=h1

    b(i-1)=2.d0*(h1+h2)

    c(i-1)=h2


    rhs(i-1)=3.d0*( &
          f2*h2*h2 &
        + f(i)*(h1*h1-h2*h2) &
        - f1*h1*h1 )

enddo
  
  
  !------------------------------------------------------
  ! Cyclic corner terms
  !
  ! first row:
  ! b1*x1+c1*x2+beta*xm
  !
  ! last row:
  ! alpha*x1+a_m*x_{m-1}+b_m*xm
  !
  !------------------------------------------------------
  
  call cyclic_solve(m,a,b,c,a(1),c(m),rhs,sol)
  
  
  do i=2,n
      d(i)=sol(i-1)
  enddo
  
  d(1)=d(n)
  
  
  end subroutine tb15c

  subroutine cyclic_solve(n,a,b,c,alpha,beta,r,x)

    implicit none
    
    integer,intent(in)::n
    
    real(8),intent(in)::a(n),b(n),c(n)
    real(8),intent(in)::alpha,beta
    real(8),intent(in)::r(n)
    
    real(8),intent(out)::x(n)
    
    real(8)::bb(n),u(n),z(n)
    real(8)::gamma,factor
    integer::i
    
    
    gamma=-b(1)
    
    
    bb=b
    bb(1)=b(1)-gamma
    bb(n)=b(n)-alpha*beta/gamma
    
    
    call tridag(n,a,bb,c,r,x)
    
    
    u=0.d0
    u(1)=gamma
    u(n)=alpha
    
    
    call tridag(n,a,bb,c,u,z)
    
    
    factor=(x(1)+beta*x(n)/gamma) / &
           (1.d0+z(1)+beta*z(n)/gamma)
    
    
    do i=1,n
        x(i)=x(i)-factor*z(i)
    enddo
    
    
    end subroutine cyclic_solve
    subroutine tridag(n,a,b,c,r,u)

implicit none

integer,intent(in)::n

real(8),intent(in)::a(n),b(n),c(n),r(n)
real(8),intent(out)::u(n)

real(8)::gam(n)
real(8)::bet
integer::i


bet=b(1)

if(bet==0.d0) stop "zero pivot"


u(1)=r(1)/bet


do i=2,n

    gam(i)=c(i-1)/bet

    bet=b(i)-a(i)*gam(i)

    if(bet==0.d0) stop "zero pivot"

    u(i)=(r(i)-a(i)*u(i-1))/bet

enddo


do i=n-1,1,-1

    u(i)=u(i)-gam(i+1)*u(i+1)

enddo


end subroutine tridag