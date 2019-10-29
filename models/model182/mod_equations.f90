module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
    type(action), dimension(:), allocatable :: rhs1dt0seq, rhs1dt1seq, rhs2dt0seq, rhs2dt1seq, rhs3dt0seq, rhs5dt0seq, rhs5dt1seq
    type(action), dimension(:), allocatable :: amat11dt0seq, amat11dt1seq, amat12dt1seq, amat13dt1seq, amat14dt1seq
    type(action), dimension(:), allocatable :: amat21dt1seq, amat22dt0seq, amat22dt1seq, amat23dt1seq, amat24dt1seq
    type(action), dimension(:), allocatable :: amat31dt0seq, amat33dt0seq
    type(action), dimension(:), allocatable :: amat52dt1seq, amat55dt0seq, amat55dt1seq
    
    real*8, dimension(:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: zj0        = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: w0         = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: delta_Psi  = algexpr(basic=.true.,var=7)
  type(algexpr), parameter, private :: delta_Phi  = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: delta_zj   = algexpr(basic=.true.,var=9)
  type(algexpr), parameter, private :: delta_w    = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: delta_rho  = algexpr(basic=.true.,var=11)
  type(algexpr), parameter, private :: delta_T    = algexpr(basic=.true.,var=12)
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=13)
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: zj         = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: w          = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: sqrtT0c    = algexpr(basic=.true.,var=15)
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=16)
  type(algexpr), parameter, private :: psi_v      = algexpr(basic=.true.,var=17)
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=18)
  type(algexpr), parameter, private :: D_perp     = algexpr(basic=.true.,var=19)
  type(algexpr), parameter, private :: k_perp     = algexpr(basic=.true.,var=20)
  type(algexpr), parameter, private :: S_rho      = algexpr(basic=.true.,var=21)
  type(algexpr), parameter, private :: S_e        = algexpr(basic=.true.,var=22)
  
  type(algexpr), private :: Bv2
  
  type(algexpr), private :: rhs1dt0, rhs1dt1, rhs2dt0, rhs2dt1, rhs3dt0, rhs3dt1, rhs4dt0, rhs4dt1, rhs5dt0, rhs5dt1
  type(algexpr), private :: amat11dt0, amat11dt1, amat12dt1, amat13dt1, amat14dt1
  type(algexpr), private :: amat21dt1, amat22dt0, amat22dt1, amat23dt1, amat24dt1
  type(algexpr), private :: amat31dt0, amat33dt0
  type(algexpr), private :: amat41dt0, amat41dt1, amat42dt0, amat42dt1, amat43dt0, amat43dt1, amat44dt0, amat44dt1
  type(algexpr), private :: amat52dt1, amat55dt0, amat55dt1
  
  type(algexpr), private :: rhs1dt1e, rhs2dt1e, rhs3dt1e, rhs4dt1e, rhs5dt1e
  type(algexpr), private :: amat11dt1e, amat12dt1e
  type(algexpr), private :: amat21dt1e, amat22dt1e, amat23dt1e
  type(algexpr), private :: amat32dt1e, amat33dt1e
  type(algexpr), private :: amat41dt1e, amat42dt1e, amat43dt1e, amat44dt1e
  type(algexpr), private :: amat52dt1e, amat55dt1e
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, eta, visco, zk_par, eta_T_dependent, visco_T_dependent, T_0, gamma
    implicit none
    integer       :: i, last
    real*8        :: zeta
    real*8        :: theta
    type(algexpr) :: j0x, j0y, j0p
    type(algexpr) :: tjx, tjy, tjp
    
    zeta  = time_evol_zeta
    theta = time_evol_theta
    
    Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)
    
    j0x =  -(zj0 + 2.d0*dx(Psi0)/R)*dx(chi) + Bv_parderiv(dx(Psi0)) - gradprod(Psi0,dx(chi))
    j0y =  -(zj0 + 2.d0*dx(Psi0)/R)*dy(chi) + Bv_parderiv(dy(Psi0)) - gradprod(Psi0,dy(chi))
    j0p = (-(zj0 + 2.d0*dx(Psi0)/R)*dp(chi) + Bv_parderiv(dp(Psi0)) - gradprod(Psi0,dp(chi)))/R + 2.d0*(dx(Psi0)*dp(chi) - dp(Psi0)*dx(chi))/(R*R)
    
    tjx =  -2.d0*dx(Psi)*dx(chi)/R + Bv_parderiv(dx(Psi)) - gradprod(Psi,dx(chi))
    tjy =  -2.d0*dx(Psi)*dy(chi)/R + Bv_parderiv(dy(Psi)) - gradprod(Psi,dy(chi))
    tjp = (-2.d0*dx(Psi)*dp(chi)/R + Bv_parderiv(dp(Psi)) - gradprod(Psi,dp(chi)))/R + 2.d0*(dx(Psi)*dp(chi) - dp(Psi)*dx(chi))/(R*R)
    
    rhs1dt0 = zeta*v*Bv_pbrack(psi_v,delta_Psi)
!    if (eta_T_dependent) then
!      rhs1dt1 = v*Bv_pbrack((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))/Bv2,psi_v) &
!              - (eta*T_0**(1.5d0))*(dx(v)*(dy(psi_v)*j0p - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*j0p) &
!              + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)/(sqrtT0c*sqrtT0c*sqrtT0c)
!    else
      rhs1dt1 = & ! v*Bv_pbrack((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))/Bv2,psi_v) &
              - eta*(dx(v)*(dy(psi_v)*j0p - dp(psi_v)*j0y/R) &
              + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*j0p) + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)
!    end if
    
    rhs2dt0 = -zeta*inprod(v,delta_Phi)
    if (visco_T_dependent) then
      rhs2dt1 = pLap(Phi0)*Bv_pbrack(v,Phi0)/Bv2 - Bv2*(j0x*dx(v) + j0y*dy(v) + j0p*dp(v)/R)/rho0 &
              + (dx(chi)*j0x + dy(chi)*j0y + dp(chi)*j0p/R)*Bv_parderiv(v)/rho0 &
              - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*Bv_pbrack(v,Psi0)/rho0 - D_perp*gradprod(rho0,inprod(v,Phi0)/rho0) &
              + S_rho*inprod(v,Phi0)/rho0 - v*Bv_pbrack(rho0,T0)/rho0 + (visco*T_0**(1.5d0))*Lap(v)*pLap(Phi0)/(sqrtT0c*sqrtT0c*sqrtT0c)
    else
      rhs2dt1 = pLap(Phi0)*Bv_pbrack(v,Phi0)/Bv2 - Bv2*(j0x*dx(v) + j0y*dy(v) + j0p*dp(v)/R)/rho0 &
              + (dx(chi)*j0x + dy(chi)*j0y + dp(chi)*j0p/R)*Bv_parderiv(v)/rho0 &
              - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*Bv_pbrack(v,Psi0)/rho0 - D_perp*gradprod(rho0,inprod(v,Phi0)/rho0) &
              + S_rho*inprod(v,Phi0)/rho0 - v*Bv_pbrack(rho0,T0)/rho0 + visco*Lap(v)*pLap(Phi0)
    end if
    
    rhs3dt0 = v*zj0 + gradprod(v,Psi0)
    
    rhs5dt0 = zeta*v*delta_rho; rhs5dt1 = -v*Bv_pbrack(rho0/Bv2,Phi0) - D_perp*gradprod(v,rho0) + S_rho*v
    
    
    amat11dt0 = (1.d0 + zeta)*v*Bv_pbrack(psi_v,Psi)
!    if (eta_T_dependent) then
!      amat11dt1 = -theta*v*Bv_pbrack(Bv_pbrack(Psi,Phi0)/Bv2,psi_v) + (theta*eta*T_0**(1.5d0))*(dx(v)*(dy(psi_v)*tjp - dp(psi_v)*tjy/R) &
!                + dy(v)*(dp(psi_v)*tjx/R - dx(psi_v)*tjp) + dp(v)*(dx(psi_v)*tjy - dy(psi_v)*tjx)/R)/(sqrtT0c*sqrtT0c*sqrtT0c)
!    else
      amat11dt1 = & ! -theta*v*Bv_pbrack(Bv_pbrack(Psi,Phi0)/Bv2,psi_v) + &
                  theta*eta*(dx(v)*(dy(psi_v)*tjp - dp(psi_v)*tjy/R) &
                + dy(v)*(dp(psi_v)*tjx/R - dx(psi_v)*tjp) + dp(v)*(dx(psi_v)*tjy - dy(psi_v)*tjx)/R)
!    end if
    amat12dt1 = -theta*v*Bv_pbrack((Bv_pbrack(Psi0,Phi) - Bv_parderiv(Phi))/Bv2,psi_v)
    amat13dt1 = -theta*eta*zj*Bv_pbrack(v,psi_v)
    
    amat31dt0 = -gradprod(v,Psi)
    amat33dt0 = -v*zj
    
    amat52dt1 = theta*v*Bv_pbrack(rho0/Bv2,Phi)
    amat55dt0 = (1.d0 + zeta)*v*rho; amat55dt1 = theta*(v*Bv_pbrack(rho/Bv2,Phi0) + D_perp*gradprod(v,rho))
    
    rhs1dt1e = Dexpand(deepcopy(rhs1dt1))
    rhs2dt1e = Dexpand(deepcopy(rhs2dt1))
    rhs5dt1e = Dexpand(deepcopy(rhs5dt1))
    
    amat11dt1e = Dexpand(deepcopy(amat11dt1)); amat12dt1e = Dexpand(deepcopy(amat12dt1))
    amat52dt1e = Dexpand(deepcopy(amat52dt1)); amat55dt1e = Dexpand(deepcopy(amat55dt1))
    
    
    allocate(thread_eq(nbthreads))
    do i=1,nbthreads
      allocate(thread_eq(i)%eq(2*n_var+11,0:n_order-1,0:n_order-1,0:n_order-1))
    
      allocate(thread_eq(i)%rhs1dt0seq(countsubexprs(rhs1dt0)))
      last = 0; call buildsequence(rhs1dt0, thread_eq(i)%rhs1dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%rhs1dt1seq(countsubexprs(rhs1dt1e)))
      last = 0; call buildsequence(rhs1dt1e, thread_eq(i)%rhs1dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%rhs2dt0seq(countsubexprs(rhs2dt0)))
      last = 0; call buildsequence(rhs2dt0, thread_eq(i)%rhs2dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%rhs2dt1seq(countsubexprs(rhs2dt1e)))
      last = 0; call buildsequence(rhs2dt1e, thread_eq(i)%rhs2dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%rhs3dt0seq(countsubexprs(rhs3dt0)))
      last = 0; call buildsequence(rhs3dt0, thread_eq(i)%rhs3dt0seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%rhs5dt0seq(countsubexprs(rhs5dt0)))
      last = 0; call buildsequence(rhs5dt0, thread_eq(i)%rhs5dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%rhs5dt1seq(countsubexprs(rhs5dt1e)))
      last = 0; call buildsequence(rhs5dt1e, thread_eq(i)%rhs5dt1seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat11dt0seq(countsubexprs(amat11dt0)))
      last = 0; call buildsequence(amat11dt0, thread_eq(i)%amat11dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat11dt1seq(countsubexprs(amat11dt1e)))
      last = 0; call buildsequence(amat11dt1e, thread_eq(i)%amat11dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat12dt1seq(countsubexprs(amat12dt1e)))
      last = 0; call buildsequence(amat12dt1e, thread_eq(i)%amat12dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat13dt1seq(countsubexprs(amat13dt1)))
      last = 0; call buildsequence(amat13dt1, thread_eq(i)%amat13dt1seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat31dt0seq(countsubexprs(amat31dt0)))
      last = 0; call buildsequence(amat31dt0, thread_eq(i)%amat31dt0seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat33dt0seq(countsubexprs(amat33dt0)))
      last = 0; call buildsequence(amat33dt0, thread_eq(i)%amat33dt0seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat52dt1seq(countsubexprs(amat52dt1e)))
      last = 0; call buildsequence(amat52dt1e, thread_eq(i)%amat52dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat55dt0seq(countsubexprs(amat55dt0)))
      last = 0; call buildsequence(amat55dt0, thread_eq(i)%amat55dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat55dt1seq(countsubexprs(amat55dt1e)))
      last = 0; call buildsequence(amat55dt1e, thread_eq(i)%amat55dt1seq, thread_eq(i)%eq, last)
    end do
  end subroutine init_equations
  
  type(algexpr) function Bv_pbrack(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
  
    Bv_pbrack = ((dy(a)*dp(b) - dp(a)*dy(b))*dx(chi) + (dp(a)*dx(b) - dx(a)*dp(b))*dy(chi) + (dx(a)*dy(b) - dy(a)*dx(b))*dp(chi))/R
  end function Bv_pbrack
  
  type(algexpr) function Bv_parderiv(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    Bv_parderiv = dx(a)*dx(chi) + dy(a)*dy(chi) + dp(a)*dp(chi)/(R*R)
  end function Bv_parderiv
  
  type(algexpr) function gradprod(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
    
    gradprod = dx(a)*dx(b) + dy(a)*dy(b) + dp(a)*dp(b)/(R*R)
  end function gradprod

  type(algexpr) function inprod(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
  
    inprod = gradprod(a,b) - Bv_parderiv(a)*Bv_parderiv(b)/Bv2
  end function inprod
  
  type(algexpr) function Lap(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    Lap = dx(R*dx(a))/R + dy(dy(a)) + dp(dp(a))/(R*R)
  end function Lap
  
  type(algexpr) function pLap(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    pLap = Lap(a) - Bv_parderiv(Bv_parderiv(a)/Bv2)
  end function pLap
end module mod_equations