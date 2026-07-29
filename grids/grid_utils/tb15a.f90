! Todo: make into module
subroutine tb15a(n,x,f,d,w,lp)
  !-------------------------------
  ! Calculate periodic spline parameters
  !-------------------------------
  implicit none

  integer,intent(in) :: n,lp      ! Number of points, unit number 
  real(8),intent(in) :: x(n),f(n) ! evaluation points and function values
  real(8),intent(out):: d(n)      ! Result: Derivatives at knots
  real(8),intent(inout):: w(*)    ! Workspace for comptability reasons

  real(8) :: a(n-1),b(n-1),c(n-1),rhs(n-1)
  real(8) :: sol(n-1)

  real(8) :: h1,h2,f1,f2
  integer :: i,m



  if (f(n) .ne. f(1)) then
     write(*,*) "Function values must be periodic", f(n), f(1)
     stop
  end if

  if (n < 4) then
     write(*,*) "More than 4 points needed for spline calculation"
     stop
  endif

  do i = 2, n
     if (x(i) .le. x(i-1)) then
        write(*,*) "Coordinate points must be increasing"
        stop
     end if
  end do

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

     a(i-1)=h1             ! Subdiagonal
     b(i-1)=2.d0*(h1+h2)   ! Diagonal
     c(i-1)=h2             ! Supdiagonal


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


end subroutine tb15a

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
