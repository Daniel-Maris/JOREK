module mod_semianalytical
  implicit none

  ! A simple algebraic expression with two operands and an arithmetical operator (+, -, * or \)
  ! The operands can either be basics variables or algebraic expressions themselves
  type algexpr
    type(algexpr), pointer :: operand1 => NULL(), operand2 => NULL()
    character              :: oprtr = ''
    integer                :: dx = 0, dy = 0, dp = 0    ! Orders of differential operators acting on this expression
    logical                :: basic = .false.           ! .true. if there are no more sub-operands in the expression
    integer                :: var = 0                   ! The index of the basic variable. Only initialized if basic .eq. .true.
    real*8                 :: factor = 1.d0, add = 0.d0 ! A numerical multiplicative factor and an additive constant
    
    ! These pointers are only for debugging purposes and should be removed in production
    ! up points to the expression for which this expression is an operand
    ! origin points to the un-expanded expression from which this expression was obtained
    type(algexpr), pointer :: up => NULL(), origin => NULL()
  end type algexpr

  ! An instruction to add, subtract, multiply or divide two numbers
  ! Algebraic expressions are compiled to sequences (arrays) of actions
  type action
    real*8, pointer :: v1, v2      ! Point to two values, which are to be added, subtracted, multiplied or divided
    integer*1       :: c1, c2, c3  ! These three constants collectively determine the type of arithmetic operation to be performed
    real*8          :: reslt       ! Stores the result
    real*8          :: f1, f2      ! Numerical factors individually multiplying v1 and v2
    real*8          :: factor, add ! A numerical multiplicative factor (for entire expression) and an additive constant
    
    ! Only for debugging purposes: points to the exact expression from which this action was obtained
    type(algexpr), pointer :: origin => NULL()
  end type action
  
  type(algexpr), parameter :: one = algexpr(basic = .true.,var=1,factor=0,add=1)

  ! Operators for making algebraic expressions
  interface operator (+)
    procedure addexpr
    procedure addexprn
    procedure addnexpr
  end interface operator (+)
  
  interface operator (-)
    procedure subexpr
    procedure subexprn
    procedure subnexpr
    procedure negate
  end interface operator (-)
  
  interface operator (*)
    procedure multexpr
    procedure multexprn
    procedure multnexpr
  end interface operator (*)
  
  interface operator (/)
    procedure divexpr
    procedure divexprn
    procedure divnexpr
  end interface operator (/)
  
contains

  ! Functions implementing the operators
  type(algexpr) function addexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(addexpr%operand1)
    allocate(addexpr%operand2)
    
    addexpr%operand1 = e1
    addexpr%operand2 = e2
    addexpr%oprtr    = '+'
  end function addexpr
  
  type(algexpr) function addexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    
    addexprn = e1
    addexprn%add = addexprn%add + n2
  end function addexprn
  
  type(algexpr) function addnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    
    addnexpr = e2
    addnexpr%add = addnexpr%add + n1
  end function addnexpr
  
  type(algexpr) function subexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(subexpr%operand1)
    allocate(subexpr%operand2)
    
    subexpr%operand1 = e1
    subexpr%operand2 = e2
    subexpr%oprtr    = '-'
  end function subexpr
  
  type(algexpr) function subexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    
    subexprn = e1
    subexprn%add = subexprn%add - n2
  end function subexprn
  
  type(algexpr) function subnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    
    subnexpr = e2
    subnexpr%factor = -subnexpr%factor
    subnexpr%add = n1 - subnexpr%add
  end function subnexpr
  
  type(algexpr) function negate(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    
    negate = expr
    negate%factor = -negate%factor
    negate%add = -negate%add
  end function negate
  
  type(algexpr) function multexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(multexpr%operand1)
    allocate(multexpr%operand2)
    
    multexpr%operand1 = e1
    multexpr%operand2 = e2
    multexpr%oprtr    = '*'
  end function multexpr
  
  type(algexpr) function multexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    
    multexprn = e1
    multexprn%factor = n2*multexprn%factor
    multexprn%add = n2*multexprn%add
  end function multexprn
  
  type(algexpr) function multnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    
    multnexpr = e2
    multnexpr%factor = n1*multnexpr%factor
    multnexpr%add = n1*multnexpr%add
  end function multnexpr
  
  type(algexpr) function divexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(divexpr%operand1)
    allocate(divexpr%operand2)
    
    divexpr%operand1 = e1
    divexpr%operand2 = e2
    divexpr%oprtr    = '/'
  end function divexpr
  
  type(algexpr) function divexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    
    divexprn = e1
    divexprn%factor = divexprn%factor/n2
    divexprn%add = divexprn%add/n2
  end function divexprn
  
  type(algexpr) function divnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    
    allocate(divnexpr%operand1)
    allocate(divnexpr%operand2)
    
    divnexpr%operand1        = one
    divnexpr%operand1%factor = n1
    divnexpr%operand2        = e2
    divnexpr%oprtr           = '/'
  end function divnexpr
  
  ! Functions applying differential operators to the expression
  type(algexpr) function dx(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    dx    = expr
    dx%dx = dx%dx + 1
  end function dx
  
  type(algexpr) function dy(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    dy    = expr
    dy%dy = dy%dy + 1
  end function dy
  
  type(algexpr) function dp(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    dp    = expr
    dp%dp = dp%dp + 1
  end function dp
  
  ! This function makes a deep copy of the entire expression tree
  type(algexpr) recursive function deepcopy(expr) result(res)
    implicit none
    type(algexpr), intent(in) :: expr
    
    if (expr%basic) then
      res = expr
    else
      res%oprtr  = expr%oprtr
      res%dx     = expr%dx
      res%dy     = expr%dy
      res%dp     = expr%dp
      res%var    = expr%var
      res%factor = expr%factor
      res%add    = expr%add
      
      res%origin => expr%origin
      res%up     => expr%up
      
      allocate(res%operand1)
      allocate(res%operand2)
      res%operand1 = deepcopy(expr%operand1)
      res%operand2 = deepcopy(expr%operand2)
    end if
  end function deepcopy
  
  ! This subroutine initializes the up pointers in an expression tree (only useful for debugging)
  recursive subroutine inituptree(expr)
    implicit none
    type(algexpr), target, intent(inout) :: expr
    
    if (.not. expr%basic) then
      expr%operand1%up => expr
      expr%operand2%up => expr
      call inituptree(expr%operand1)
      call inituptree(expr%operand2)
    end if
  end subroutine inituptree

  type(algexpr) recursive function Dexpand(expr) result(res)
    implicit none
    type(algexpr), target, intent(in) :: expr
    type(algexpr), pointer            :: oldop1, oldop2, rescp
    
    if (.not. expr%basic) then
      res = deepcopy(expr)
      if (expr%dx .ne. 0) then
        select case (expr%oprtr)
          case ('+')
            res = dx(res%operand1) + dx(res%operand2)
          case ('-')
            res = dx(res%operand1) - dx(res%operand2)
          case ('*')
            res = dx(res%operand1)*res%operand2 + res%operand1*dx(res%operand2)
          case ('/')
            res = (res%operand2*dx(res%operand1) - res%operand1*dx(res%operand2))/(res%operand2*res%operand2)
        end select
        res%dx = expr%dx - 1
        res%dy = expr%dy
        res%dp = expr%dp
        res%add = 0
        res%origin => expr
      else if (expr%dy .ne. 0) then
        select case (expr%oprtr)
          case ('+')
            res = dy(res%operand1) + dy(res%operand2)
          case ('-')
            res = dy(res%operand1) - dy(res%operand2)
          case ('*')
            res = dy(res%operand1)*res%operand2 + res%operand1*dy(res%operand2)
          case ('/')
            res = (res%operand2*dy(res%operand1) - res%operand1*dy(res%operand2))/(res%operand2*res%operand2)
        end select
        res%dx = expr%dx
        res%dy = expr%dy - 1
        res%dp = expr%dp
        res%add = 0
        res%origin => expr
      else if (expr%dp .ne. 0) then
        select case (expr%oprtr)
          case ('+')
            res = dp(res%operand1) + dp(res%operand2)
          case ('-')
            res = dp(res%operand1) - dp(res%operand2)
          case ('*')
            res = dp(res%operand1)*res%operand2 + res%operand1*dp(res%operand2)
          case ('/')
            res = (res%operand2*dp(res%operand1) - res%operand1*dp(res%operand2))/(res%operand2*res%operand2)
        end select
        res%dx = expr%dx
        res%dy = expr%dy
        res%dp = expr%dp - 1
        res%add = 0
        res%origin => expr
      end if
      res%factor = expr%factor
      
      if (res%dx .ne. 0 .or. res%dy .ne. 0 .or. res%dp .ne. 0) then
        allocate(rescp)
        rescp = res
        res = Dexpand(rescp)
      else
        oldop1 => res%operand1
        oldop2 => res%operand2
        allocate(res%operand1)
        allocate(res%operand2)
        res%operand1 = Dexpand(oldop1)
        res%operand2 = Dexpand(oldop2)
      end if
    else
      res = expr
      if (res%dx .ne. 0 .or. res%dy .ne. 0 .or. res%dp .ne. 0) res%add = 0
    end if
  end function Dexpand
  
  ! Print an expression using basic text
  recursive subroutine printexpr(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    
    if (expr%factor .ne. 1) write(*,"(F10.3,A)",advance='no') expr%factor,"*"
    
    if (expr%dx .ne. 0) write(*,"(A,I1)",advance='no') "dx",expr%dx
    if (expr%dy .ne. 0) write(*,"(A,I1)",advance='no') "dy",expr%dy
    if (expr%dp .ne. 0) write(*,"(A,I1)",advance='no') "dp",expr%dp
    
    if (expr%basic) then
      write(*,"(A,I1,A)",advance='no') "[var",expr%var,"]"
    else
      write(*,"(A)",advance='no') "("
      call printexpr(expr%operand1)
      write(*,"(A)",advance='no') expr%oprtr
      call printexpr(expr%operand2)
      write(*,"(A)",advance='no') ")"
    end if
    
    if (expr%add .ne. 0) write(*,"(A,F10.3)",advance='no') " + ",expr%add
  end subroutine printexpr
  
  ! Print LaTeX code for an expression
  ! varsymb: a string containing one-character symbols (in proper order) for each variable
  recursive subroutine printlatex(expr,varsymb)
    implicit none
    type(algexpr),    intent(in) :: expr
    character(len=*), intent(in) :: varsymb
    
    if (expr%factor .ne. 1) write(*,"(F10.3,A)",advance='no') expr%factor,"*"
    
    if (expr%dx .eq. 1) then
      write(*,"(A)",advance='no') "\frac{\partial}{\partial x}"
    else if (expr%dx .gt. 1) then
      write(*,"(A,I1,A,I1,A)",advance='no') "\frac{\partial^",expr%dx,"}{\partial x^",expr%dx,"}"
    end if
    
    if (expr%dy .eq. 1) then
      write(*,"(A)",advance='no') "\frac{\partial}{\partial y}"
    else if (expr%dy .gt. 1) then
      write(*,"(A,I1,A,I1,A)",advance='no') "\frac{\partial^",expr%dy,"}{\partial y^",expr%dy,"}"
    end if
    
    if (expr%dp .eq. 1) then
      write(*,"(A)",advance='no') "\frac{\partial}{\partial p}"
    else if (expr%dp .gt. 1) then
      write(*,"(A,I1,A,I1,A)",advance='no') "\frac{\partial^",expr%dp,"}{\partial p^",expr%dp,"}"
    end if
    
    if (expr%basic) then
      write(*,"(A)",advance='no') varsymb(expr%var:expr%var)
    else
      if (expr%oprtr .eq. '/') then
        write(*,"(A)",advance='no') "\frac{"
        call printlatex(expr%operand1,varsymb)
        write(*,"(A)",advance='no') "}{"
        call printlatex(expr%operand2,varsymb)
        write(*,"(A)",advance='no') "}"
      else if (expr%oprtr .eq. '*') then
        write(*,"(A)",advance='no') "\left("
        call printlatex(expr%operand1,varsymb)
        call printlatex(expr%operand2,varsymb)
        write(*,"(A)",advance='no') "\right)"
      else
        write(*,"(A)",advance='no') "\left("
        call printlatex(expr%operand1,varsymb)
        write(*,"(A)",advance='no') expr%oprtr
        call printlatex(expr%operand2,varsymb)
        write(*,"(A)",advance='no') "\right)"
      end if
    end if
    
    if (expr%add .ne. 0) write(*,"(A,F10.3)",advance='no') " + ",expr%add
  end subroutine printlatex

  integer recursive function countsubexprs(expr) result(res)
    implicit none
    type(algexpr), intent(in) :: expr
    
    if (expr%basic) then
      res = 0
    else
      res = 1 + countsubexprs(expr%operand1) + countsubexprs(expr%operand2)
    end if
  end function countsubexprs

  ! Use an algebraic expression to build a sequence of instructions
  recursive subroutine buildsequence(expr, actseq, eq, last)
    implicit none
    type(algexpr),                       target, intent(in)    :: expr
    type(action), dimension(:),          target, intent(inout) :: actseq
    real*8,       dimension(:,0:,0:,0:), target, intent(in)    :: eq
    integer,                                     intent(inout) :: last
    type(action),                        target                :: act
    
    act%factor =  expr%factor
    act%add    =  expr%add
    act%origin => expr
    
    if (expr%operand1%basic) then
      act%v1  => eq(expr%operand1%var,expr%operand1%dx,expr%operand1%dy,expr%operand1%dp)
      act%add =  act%add + expr%operand1%add
      act%f1  =  expr%operand1%factor
    else
      call buildsequence(expr%operand1,actseq,eq,last)
      act%v1 => actseq(last)%reslt
      act%f1 =  1.d0
    end if
    
    if (expr%operand2%basic) then
      act%v2  => eq(expr%operand2%var,expr%operand2%dx,expr%operand2%dy,expr%operand2%dp)
      act%add =  act%add + expr%operand2%add
      act%f2  =  expr%operand2%factor
    else
      call buildsequence(expr%operand2,actseq,eq,last)
      act%v2 => actseq(last)%reslt
      act%f2 =  1.d0
    end if
    
    select case (expr%oprtr)
      case ('+')
        act%c1 = 1
        act%c2 = 1
        act%c3 = 1
      case ('-')
        act%c1 =  1
        act%c2 =  1
        act%c3 = -1
      case ('*')
        act%c1 = 0
        act%c2 = 1
        act%c3 = 0
      case ('/')
        act%c1 = 1
        act%c2 = 0
        act%c3 = 0
    end select
    last = last + 1
    actseq(last) = act
  end subroutine buildsequence

  ! Execute the instruction sequence
  ! Only this function is called in the loops
  real*8 function eval(actseq)
    implicit none
    type(action), dimension(:), intent(inout) :: actseq
    integer                                   :: i
    
    !dir$ noparallel
    !dir$ novector
    do i=1,size(actseq)
      actseq(i)%reslt = actseq(i)%factor*((actseq(i)%f1*actseq(i)%v1)*((1 - actseq(i)%c1)*(actseq(i)%f2*actseq(i)%v2) + actseq(i)%c1)/((1 - actseq(i)%c2)*(actseq(i)%f2*actseq(i)%v2) + actseq(i)%c2) + actseq(i)%c3*(actseq(i)%f2*actseq(i)%v2)) + actseq(i)%add
    end do
    
    eval = actseq(size(actseq))%reslt
  end function eval
end module mod_semianalytical
