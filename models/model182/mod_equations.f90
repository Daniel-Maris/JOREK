module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
    type(action), dimension(:), allocatable :: rhs1seq, rhs2seq, rhs3seq, rhs5seq
    type(action), dimension(:), allocatable :: amat11seq, amat12seq, amat13seq, amat14seq
    type(action), dimension(:), allocatable :: amat21seq, amat22seq, amat23seq, amat24seq, amat25seq
    type(action), dimension(:), allocatable :: amat31seq, amat33seq
    type(action), dimension(:), allocatable :: amat42seq, amat44seq
    type(action), dimension(:), allocatable :: amat52seq, amat55seq
    
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
  type(algexpr), parameter, private :: S_j        = algexpr(basic=.true.,var=23)
  
  type(algexpr), private :: Bv2
  
  type(algexpr), private :: rhs1, rhs2, rhs3, rhs4, rhs5
  type(algexpr), private :: amat11, amat12, amat13, amat14
  type(algexpr), private :: amat21, amat22, amat23, amat24, amat25
  type(algexpr), private :: amat31, amat33
  type(algexpr), private :: amat42, amat44
  type(algexpr), private :: amat52, amat55
  
  type(algexpr), private :: rhs1e, rhs2e, rhs3e, rhs4e, rhs5e
  type(algexpr), private :: amat11e, amat12e, amat13e
  type(algexpr), private :: amat21e, amat22e, amat24e, amat25e
  type(algexpr), private :: amat52e, amat55e
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, eta, visco, zk_par, eta_T_dependent, visco_T_dependent, T_0, gamma, tstep
    implicit none
    integer       :: i, last
    real*8        :: zeta
    real*8        :: theta
    type(algexpr) :: j0x, j0y, j0p, j0chi
    type(algexpr) :: tjx, tjy, tjp, tjchi
    
    zeta  = time_evol_zeta
    theta = time_evol_theta
    
    Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)
    
    j0x =  -(zj0 + 2.d0*dx(Psi0)/R)*dx(chi) + Bv_parderiv(dx(Psi0)) - gradprod(Psi0,dx(chi))
    j0y =  -(zj0 + 2.d0*dx(Psi0)/R)*dy(chi) + Bv_parderiv(dy(Psi0)) - gradprod(Psi0,dy(chi))
    j0p = (-(zj0 + 2.d0*dx(Psi0)/R)*dp(chi) + Bv_parderiv(dp(Psi0)) - gradprod(Psi0,dp(chi)))/R + 2.d0*(dx(Psi0)*dp(chi) - dp(Psi0)*dx(chi))/(R*R)
    j0chi = -(Bv2*zj0 + 2.d0*Bv2*dx(Psi0)/R - Bv_parderiv(Bv_parderiv(Psi0)) + Bv_parderiv(Bv2)*Bv_parderiv(Psi0)/Bv2 + inprod(Bv2,Psi0))
    
    tjx =  -2.d0*dx(Psi)*dx(chi)/R + Bv_parderiv(dx(Psi)) - gradprod(Psi,dx(chi))
    tjy =  -2.d0*dx(Psi)*dy(chi)/R + Bv_parderiv(dy(Psi)) - gradprod(Psi,dy(chi))
    tjp = (-2.d0*dx(Psi)*dp(chi)/R + Bv_parderiv(dp(Psi)) - gradprod(Psi,dp(chi)))/R + 2.d0*(dx(Psi)*dp(chi) - dp(Psi)*dx(chi))/(R*R)
    tjchi  = -(2.d0*Bv2*dx(Psi)/R - Bv_parderiv(Bv_parderiv(Psi)) + Bv_parderiv(Bv2)*Bv_parderiv(Psi)/Bv2 + inprod(Bv2,Psi))
    
    if (eta_T_dependent) then
      rhs1 = -tstep*((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))*Bv_pbrack(v,psi_v)/Bv2 &
           + (eta*T_0**(1.5d0))*(dx(v)*(dy(psi_v)*j0p - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*j0p) &
           + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)/(sqrtT0c*sqrtT0c*sqrtT0c)) &
           + zeta*v*Bv_pbrack(psi_v,delta_Psi)
    else
      rhs1 = -tstep*((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))*Bv_pbrack(v,psi_v)/Bv2 + &
             eta*(dx(v)*(dy(psi_v)*(j0p+S_j*dp(chi)/R) - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*(j0p+S_j*dp(chi)/R)) &
           + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)) &
           + zeta*v*Bv_pbrack(psi_v,delta_Psi)
    end if
    
    if (visco_T_dependent) then
      rhs2 = tstep*(w0*Bv_pbrack(v,Phi0)/Bv2 - Bv2*(j0x*dx(v) + j0y*dy(v) + j0p*dp(v)/R)/rho0 &
           + j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0))/rho0 - D_perp*gradprod(rho0,inprod(v,Phi0)/rho0) &
           + S_rho*inprod(v,Phi0)/rho0 &! - v*Bv_pbrack(rho0,T0)/rho0 
           - (visco*T_0**(1.5d0))*gradprod(v,w0)/(sqrtT0c*sqrtT0c*sqrtT0c)) &
           - zeta*inprod(v,delta_Phi)
    else
      rhs2 = tstep*(w0*Bv_pbrack(v,Phi0)/Bv2 - Bv2*(j0x*dx(v) + j0y*dy(v) + j0p*dp(v)/R)/rho0 &
           + j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0))/rho0 - D_perp*gradprod(rho0,inprod(v,Phi0)/rho0) &
           + S_rho*inprod(v,Phi0)/rho0 &! - v*Bv_pbrack(rho0,T0)/rho0 
           - visco*gradprod(v,w0)) - zeta*inprod(v,delta_Phi) !+ (0.25d0*tstep**2)*(Bv_pbrack(w0/Bv2,Phi0)*Bv_pbrack(v,Phi0) - w0*w0*inprod(v,Phi0))/Bv2
    end if
    
    rhs3 = v*zj0 + gradprod(v,Psi0)
    
    rhs4 = v*w0 + inprod(v,Phi0)
    
    rhs5 = -tstep*(v*Bv_pbrack(rho0/Bv2,Phi0) + D_perp*gradprod(v,rho0) - S_rho*v) + zeta*v*delta_rho
    
    
    if (eta_T_dependent) then
      amat11 = (1.d0 + zeta)*v*Bv_pbrack(psi_v,Psi) &
             + tstep*theta*(Bv_pbrack(Psi,Phi0)*Bv_pbrack(v,psi_v)/Bv2 + (eta*T_0**(1.5d0))*(dx(v)*(dy(psi_v)*tjp - dp(psi_v)*tjy/R) &
             + dy(v)*(dp(psi_v)*tjx/R - dx(psi_v)*tjp) + dp(v)*(dx(psi_v)*tjy - dy(psi_v)*tjx)/R)/(sqrtT0c*sqrtT0c*sqrtT0c))
    else
      amat11 = (1.d0 + zeta)*v*Bv_pbrack(psi_v,Psi) + tstep*theta*(Bv_pbrack(Psi,Phi0)*Bv_pbrack(v,psi_v)/Bv2 + &
               eta*(dx(v)*(dy(psi_v)*tjp - dp(psi_v)*tjy/R) &
             + dy(v)*(dp(psi_v)*tjx/R - dx(psi_v)*tjp) + dp(v)*(dx(psi_v)*tjy - dy(psi_v)*tjx)/R))
    end if
    amat12 = tstep*theta*(Bv_pbrack(Psi0,Phi) - Bv_parderiv(Phi))*Bv_pbrack(v,psi_v)/Bv2
    amat13 = (-tstep)*theta*eta*zj*Bv_pbrack(v,psi_v)
    
    amat21 = tstep*theta*(Bv2*(tjx*dx(v) + tjy*dy(v) + tjp*dp(v)/R) - tjchi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0)) - j0chi*Bv_pbrack(v,Psi))/rho0
    amat22 = -(1.d0 + zeta)*inprod(v,Phi) - tstep*theta*(w0*Bv_pbrack(v,Phi)/Bv2 - D_perp*gradprod(inprod(v,Phi)/rho0,rho0) + S_rho*inprod(v,Phi)/rho0) !&
           !- (0.25d0*tstep**2)*(Bv_pbrack(w0/Bv2,Phi)*Bv_pbrack(v,Phi0) + Bv_pbrack(w0/Bv2,Phi0)*Bv_pbrack(v,Phi) - w0*w0*inprod(v,Phi))/Bv2
    amat23 = tstep*theta*Bv2*zj*Bv_pbrack(v,Psi0)/rho0
    if (visco_T_dependent) then
      amat24 = -tstep*theta*w*Bv_pbrack(v,Phi0)/Bv2 + (tstep*theta*visco*T_0**(1.5d0))*gradprod(v,w)/(sqrtT0c*sqrtT0c*sqrtT0c)
    else
      amat24 = -tstep*theta*w*Bv_pbrack(v,Phi0)/Bv2 + tstep*theta*visco*gradprod(v,w) !- (0.25d0*tstep**2)*(Bv_pbrack(w/Bv2,Phi0)*Bv_pbrack(v,Phi0) - 2.d0*w0*w*inprod(v,Phi0))/Bv2
    end if
    amat25 = -tstep*theta*(Bv2*(j0x*dx(v) + j0y*dy(v) + j0p*dp(v)/R) - j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0)) &
           - S_rho*inprod(v,Phi0) &! + v*Bv_pbrack(rho0,T0)
             )*rho/(rho0*rho0) + tstep*theta*(D_perp*gradprod(inprod(v,Phi0)/rho0,rho) &
           - D_perp*gradprod(inprod(v,Phi0)*rho/(rho0*rho0),rho0) )! + v*Bv_pbrack(rho,T0)/rho0)
    
    amat31 = gradprod(v,Psi)
    amat33 = v*zj
    
    amat42 = inprod(v,Phi)
    amat44 = v*w
    
    amat52 = tstep*theta*v*Bv_pbrack(rho0/Bv2,Phi)
    amat55 = (1.d0 + zeta)*v*rho + tstep*theta*(v*Bv_pbrack(rho/Bv2,Phi0) + D_perp*gradprod(v,rho))
    
    rhs1e = Dexpand(deepcopy(rhs1))
    rhs2e = Dexpand(deepcopy(rhs2))
    rhs5e = Dexpand(deepcopy(rhs5))
    
    amat11e = Dexpand(deepcopy(amat11)); amat12e = Dexpand(deepcopy(amat12))
    amat21e = Dexpand(deepcopy(amat21)); amat22e = Dexpand(deepcopy(amat22)); amat24e = Dexpand(deepcopy(amat24)); amat25e = Dexpand(deepcopy(amat25))
    amat52e = Dexpand(deepcopy(amat52)); amat55e = Dexpand(deepcopy(amat55))
    
    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        allocate(thread_eq(i)%eq(2*n_var+11,0:n_order-1,0:n_order-1,0:n_order-1))
        allocate(thread_eq(i)%rhs1seq(countsubexprs(rhs1e)))
        allocate(thread_eq(i)%rhs2seq(countsubexprs(rhs2e)))
        allocate(thread_eq(i)%rhs3seq(countsubexprs(rhs3)))
        allocate(thread_eq(i)%rhs5seq(countsubexprs(rhs5e)))
        allocate(thread_eq(i)%amat11seq(countsubexprs(amat11e)))
        allocate(thread_eq(i)%amat12seq(countsubexprs(amat12e)))
        allocate(thread_eq(i)%amat13seq(countsubexprs(amat13)))
        allocate(thread_eq(i)%amat21seq(countsubexprs(amat21e)))
        allocate(thread_eq(i)%amat22seq(countsubexprs(amat22e)))
        allocate(thread_eq(i)%amat23seq(countsubexprs(amat23)))
        allocate(thread_eq(i)%amat24seq(countsubexprs(amat24e)))
        allocate(thread_eq(i)%amat25seq(countsubexprs(amat25e)))
        allocate(thread_eq(i)%amat31seq(countsubexprs(amat31)))
        allocate(thread_eq(i)%amat33seq(countsubexprs(amat33)))
        allocate(thread_eq(i)%amat42seq(countsubexprs(amat42)))
        allocate(thread_eq(i)%amat44seq(countsubexprs(amat44)))
        allocate(thread_eq(i)%amat52seq(countsubexprs(amat52e)))
        allocate(thread_eq(i)%amat55seq(countsubexprs(amat55e)))
      end do
    end if
    
    do i=1,nbthreads
      call buildsequence(rhs1e, thread_eq(i)%rhs1seq, thread_eq(i)%eq)
      call buildsequence(rhs2e, thread_eq(i)%rhs2seq, thread_eq(i)%eq)
      call buildsequence(rhs3, thread_eq(i)%rhs3seq, thread_eq(i)%eq)
      call buildsequence(rhs5e, thread_eq(i)%rhs5seq, thread_eq(i)%eq)
      
      call buildsequence(amat11e, thread_eq(i)%amat11seq, thread_eq(i)%eq)
      call buildsequence(amat12e, thread_eq(i)%amat12seq, thread_eq(i)%eq)
      call buildsequence(amat13, thread_eq(i)%amat13seq, thread_eq(i)%eq)
      
      call buildsequence(amat21e, thread_eq(i)%amat21seq, thread_eq(i)%eq)
      call buildsequence(amat22e, thread_eq(i)%amat22seq, thread_eq(i)%eq)
      call buildsequence(amat23, thread_eq(i)%amat23seq, thread_eq(i)%eq)
      call buildsequence(amat24e, thread_eq(i)%amat24seq, thread_eq(i)%eq)
      call buildsequence(amat25e, thread_eq(i)%amat25seq, thread_eq(i)%eq)
      
      call buildsequence(amat31, thread_eq(i)%amat31seq, thread_eq(i)%eq)
      call buildsequence(amat33, thread_eq(i)%amat33seq, thread_eq(i)%eq)
      
      call buildsequence(amat42, thread_eq(i)%amat42seq, thread_eq(i)%eq)
      call buildsequence(amat44, thread_eq(i)%amat44seq, thread_eq(i)%eq)
      
      call buildsequence(amat52e, thread_eq(i)%amat52seq, thread_eq(i)%eq)
      call buildsequence(amat55e, thread_eq(i)%amat55seq, thread_eq(i)%eq)
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