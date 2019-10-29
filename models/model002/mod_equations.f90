module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
    type(action), dimension(:), allocatable :: rhs1dt0seq, rhs1dt1seq, rhs1dt2seq
    type(action), dimension(:), allocatable :: amat11dt0seq, amat11dt1seq, amat11dt2seq, amat12dt1seq, amat12dt2seq
    type(action), dimension(:), allocatable :: amat21dt0seq, amat22dt0seq
    
    real*8, dimension(:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  type(algexpr), parameter, private :: u0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: w0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: delta_u  = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: delta_w  = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: v        = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: u        = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: w        = algexpr(basic=.true.,var=6)
  
  type(algexpr), private :: rhs1dt0, rhs1dt1, rhs1dt2
  type(algexpr), private :: amat11dt0, amat11dt1, amat11dt2, amat12dt1, amat12dt2, amat21dt0, amat22dt0
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, visco
    implicit none
    integer :: i, last
    real*8  :: zeta
    real*8  :: theta
    
    zeta  = time_evol_zeta
    theta = time_evol_theta
    
    rhs1dt0 = -zeta*inprod(v,delta_u); rhs1dt1 = -w0*pbrack(v,u0) - visco*inprod(v,w0); rhs1dt2 = -0.25d0*pbrack(w0,u0)*pbrack(v,u0)
    
    amat11dt0 = -(1.d0 + zeta)*inprod(v,u); amat11dt1 = theta*w0*pbrack(v,u)
    amat11dt2 = 0.25d0*theta*(pbrack(w0,u)*pbrack(v,u0) + pbrack(w0,u0)*pbrack(v,u))
    
    amat12dt1 = theta*w*pbrack(v,u0) + theta*visco*inprod(v,w); amat12dt2 = 0.25d0*theta*pbrack(w,u0)*pbrack(v,u0)
    
    amat21dt0 = inprod(v,u)
    
    amat22dt0 = v*w
    
    allocate(thread_eq(nbthreads))
    do i=1,nbthreads
      allocate(thread_eq(i)%eq(2*n_var+2,0:n_order-1,0:n_order-1,0:n_order-1))
    
      allocate(thread_eq(i)%rhs1dt0seq(countsubexprs(rhs1dt0)))
      last = 0; call buildsequence(rhs1dt0, thread_eq(i)%rhs1dt0seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%rhs1dt1seq(countsubexprs(rhs1dt1)))
      last = 0; call buildsequence(rhs1dt1, thread_eq(i)%rhs1dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%rhs1dt2seq(countsubexprs(rhs1dt2)))
      last = 0; call buildsequence(rhs1dt2, thread_eq(i)%rhs1dt2seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat11dt0seq(countsubexprs(amat11dt0)))
      last = 0; call buildsequence(amat11dt0, thread_eq(i)%amat11dt0seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat11dt1seq(countsubexprs(amat11dt1)))
      last = 0; call buildsequence(amat11dt1, thread_eq(i)%amat11dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat11dt2seq(countsubexprs(amat11dt2)))
      last = 0; call buildsequence(amat11dt2, thread_eq(i)%amat11dt2seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat12dt1seq(countsubexprs(amat12dt1)))
      last = 0; call buildsequence(amat12dt1, thread_eq(i)%amat12dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat12dt2seq(countsubexprs(amat12dt2)))
      last = 0; call buildsequence(amat12dt2, thread_eq(i)%amat12dt2seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat21dt0seq(countsubexprs(amat21dt0)))
      last = 0; call buildsequence(amat21dt0, thread_eq(i)%amat21dt0seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat22dt0seq(countsubexprs(amat22dt0)))
      last = 0; call buildsequence(amat22dt0, thread_eq(i)%amat22dt0seq, thread_eq(i)%eq, last)
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