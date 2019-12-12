module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
    type(action), dimension(:), allocatable :: rhs1seq
    type(action), dimension(:), allocatable :: amat11seq, amat12seq
    type(action), dimension(:), allocatable :: amat21seq, amat22seq
    
    real*8, dimension(:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  type(algexpr), parameter, private :: u0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: w0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: delta_u  = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: delta_w  = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: v        = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: u        = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: w        = algexpr(basic=.true.,var=6)
  
  type(algexpr), private :: rhs1
  type(algexpr), private :: amat11, amat12, amat21, amat22
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, visco, tstep
    implicit none
    integer :: i, last
    real*8  :: zeta
    real*8  :: theta
    
    zeta  = time_evol_zeta
    theta = time_evol_theta
    
    rhs1 = -tstep*(w0*pbrack(v,u0) + visco*inprod(v,w0)) - (0.25d0*tstep**2)*pbrack(w0,u0)*pbrack(v,u0) - zeta*inprod(v,delta_u)
    
    amat11 = -(1.d0 + zeta)*inprod(v,u) + theta*tstep*w0*pbrack(v,u) &
           + (0.25d0*theta*tstep**2)*(pbrack(w0,u)*pbrack(v,u0) + pbrack(w0,u0)*pbrack(v,u))
    
    amat12 = theta*tstep*(w*pbrack(v,u0) + visco*inprod(v,w)) + (0.25d0*theta*tstep**2)*pbrack(w,u0)*pbrack(v,u0)
    
    amat21 = inprod(v,u)
    
    amat22 = v*w
    
    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        allocate(thread_eq(i)%eq(2*n_var+2,0:n_order-1,0:n_order-1,0:n_order-1))
        allocate(thread_eq(i)%rhs1seq(countsubexprs(rhs1)))
        allocate(thread_eq(i)%amat11seq(countsubexprs(amat11)))
        allocate(thread_eq(i)%amat12seq(countsubexprs(amat12)))
        allocate(thread_eq(i)%amat21seq(countsubexprs(amat21)))
        allocate(thread_eq(i)%amat22seq(countsubexprs(amat22)))
      end do
    end if
    
    do i=1,nbthreads
      last = 0; call buildsequence(rhs1, thread_eq(i)%rhs1seq, thread_eq(i)%eq, last)
      last = 0; call buildsequence(amat11, thread_eq(i)%amat11seq, thread_eq(i)%eq, last)
      last = 0; call buildsequence(amat12, thread_eq(i)%amat12seq, thread_eq(i)%eq, last)
      last = 0; call buildsequence(amat21, thread_eq(i)%amat21seq, thread_eq(i)%eq, last)
      last = 0; call buildsequence(amat22, thread_eq(i)%amat22seq, thread_eq(i)%eq, last)
    end do
  end subroutine init_equations
  
  type(algexpr) function pbrack(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
  
    pbrack = dx(a)*dy(b) - dy(a)*dx(b)
  end function pbrack

  type(algexpr) function inprod(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
  
    inprod = dx(a)*dx(b) + dy(a)*dy(b)
  end function inprod
end module mod_equations