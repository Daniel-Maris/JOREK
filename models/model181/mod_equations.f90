module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
    type(action), dimension(:), allocatable :: rhs1dt0seq, rhs1dt1seq, rhs2dt0seq, rhs2dt1seq, rhs3dt0seq, rhs3dt1seq, rhs4dt0seq
    type(action), dimension(:), allocatable :: rhs4dt1seq
    type(action), dimension(:), allocatable :: amat11dt0seq, amat11dt1seq, amat12dt1seq, amat14dt1seq
    type(action), dimension(:), allocatable :: amat21dt1seq, amat22dt0seq, amat22dt1seq, amat23dt1seq, amat24dt1seq
    type(action), dimension(:), allocatable :: amat32dt1seq, amat33dt0seq, amat33dt1seq
    type(action), dimension(:), allocatable :: amat41dt0seq, amat41dt1seq, amat42dt0seq, amat42dt1seq, amat43dt0seq, amat43dt1seq
    type(action), dimension(:), allocatable :: amat44dt0seq, amat44dt1seq
    
    real*8, dimension(:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: delta_Psi  = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: delta_Phi  = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: delta_rho  = algexpr(basic=.true.,var=7)
  type(algexpr), parameter, private :: delta_T    = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=9)
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: sqrtT0c    = algexpr(basic=.true.,var=11)
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=12)
  type(algexpr), parameter, private :: psi_v      = algexpr(basic=.true.,var=13)
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: D_perp     = algexpr(basic=.true.,var=15)
  type(algexpr), parameter, private :: k_perp     = algexpr(basic=.true.,var=16)
  type(algexpr), parameter, private :: S_rho      = algexpr(basic=.true.,var=17)
  type(algexpr), parameter, private :: S_e        = algexpr(basic=.true.,var=18)
  type(algexpr), parameter, private :: S_j        = algexpr(basic=.true.,var=19)
  
  type(algexpr), private :: Bv2
  
  type(algexpr), private :: rhs1dt0, rhs1dt1, rhs2dt0, rhs2dt1, rhs3dt0, rhs3dt1, rhs4dt0, rhs4dt1
  type(algexpr), private :: amat11dt0, amat11dt1, amat12dt1, amat14dt1
  type(algexpr), private :: amat21dt1, amat22dt0, amat22dt1, amat23dt1, amat24dt1
  type(algexpr), private :: amat32dt1, amat33dt0, amat33dt1
  type(algexpr), private :: amat41dt0, amat41dt1, amat42dt0, amat42dt1, amat43dt0, amat43dt1, amat44dt0, amat44dt1
  
  type(algexpr), private :: rhs1dt1e, rhs2dt1e, rhs3dt1e, rhs4dt1e
  type(algexpr), private :: amat11dt1e, amat12dt1e
  type(algexpr), private :: amat21dt1e, amat22dt1e, amat23dt1e
  type(algexpr), private :: amat32dt1e, amat33dt1e
  type(algexpr), private :: amat41dt1e, amat42dt1e, amat43dt1e, amat44dt1e
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, eta, visco, zk_par, eta_T_dependent, visco_T_dependent, T_0, gamma
    implicit none
    integer       :: i, last
    real*8        :: zeta
    real*8        :: theta
    type(algexpr) :: j0x, j0y, j0p
    type(algexpr) :: jx, jy, jp
    
    zeta  = time_evol_zeta
    theta = time_evol_theta
    
    Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)
    
    j0x =  -Lap(Psi0)*dx(chi) + Bv_parderiv(dx(Psi0)) - gradprod(Psi0,dx(chi))
    j0y =  -Lap(Psi0)*dy(chi) + Bv_parderiv(dy(Psi0)) - gradprod(Psi0,dy(chi))
    j0p = (-Lap(Psi0)*dp(chi) + Bv_parderiv(dp(Psi0)) - gradprod(Psi0,dp(chi)))/R + 2.d0*(dx(Psi0)*dp(chi) - dp(Psi0)*dx(chi))/(R*R)
    
    jx =  -Lap(Psi)*dx(chi) + Bv_parderiv(dx(Psi)) - gradprod(Psi,dx(chi))
    jy =  -Lap(Psi)*dy(chi) + Bv_parderiv(dy(Psi)) - gradprod(Psi,dy(chi))
    jp = (-Lap(Psi)*dp(chi) + Bv_parderiv(dp(Psi)) - gradprod(Psi,dp(chi)))/R + 2.d0*(dx(Psi)*dp(chi) - dp(Psi)*dx(chi))/(R*R)
    
    rhs1dt0 = zeta*v*Bv_pbrack(psi_v,delta_Psi)
!    if (eta_T_dependent) then
!      rhs1dt1 = v*Bv_pbrack((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))/Bv2,psi_v) &
!              - (eta*T_0**(1.5d0))*(dx(v)*(dy(psi_v)*j0p - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*j0p) &
!              + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)/(sqrtT0c*sqrtT0c*sqrtT0c)
!    else
      rhs1dt1 = & ! v*Bv_pbrack((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))/Bv2,psi_v) 
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
    
    rhs3dt0 = zeta*v*delta_rho; rhs3dt1 = -v*Bv_pbrack(rho0/Bv2,Phi0) - D_perp*gradprod(v,rho0) + S_rho*v
    
    rhs4dt0 = zeta*v*(delta_rho*inprod(Phi0,Phi0)/(2.d0*Bv2) + delta_rho*T0/(gamma-1.d0) + rho0*delta_T/(gamma-1.d0) &
            + rho0*inprod(Phi0,delta_Phi)/Bv2 + Bv2*inprod(Psi0,delta_Psi))
    if (eta_T_dependent) then
      rhs4dt1 = -v*Bv_pbrack(rho0*inprod(Phi0,Phi0)/(2.d0*Bv2*Bv2) + gamma*rho0*T0/((gamma-1.d0)*Bv2) + inprod(Psi0,Psi0),Phi0) &
              +  v*Bv_parderiv(inprod(Phi0,Psi0)) + v*Bv_pbrack(inprod(Phi0,Psi0),Psi0) &
              + (eta*T_0**(1.5d0))*((dx(v)*(j0y*dp(chi)/R - j0p*dy(chi)) + dy(v)*(j0p*dx(chi) - j0x*dp(chi)/R) &
              + dp(v)*(j0x*dy(chi) - j0y*dx(chi))/R) - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*gradprod(v,Psi0) &
              - (j0x*dx(Psi0) + j0y*dy(Psi0) + j0p*dp(Psi0)/R)*Bv_parderiv(v))/(sqrtT0c*sqrtT0c*sqrtT0c) - k_perp*gradprod(T0,v) &
              - T0*D_perp*gradprod(rho0,v)/(gamma-1.d0) - zk_par*Bv_parderiv(v)*(Bv_parderiv(T0) + Bv_pbrack(T0,Psi0))/Bv2 &
              + D_perp*gradprod(v*inprod(Phi0,Phi0)/Bv2,rho0)/2.d0 - v*inprod(Phi0,Phi0)*S_rho/(2.d0*Bv2) + v*S_e
    else
      rhs4dt1 = -v*Bv_pbrack(rho0*inprod(Phi0,Phi0)/(2.d0*Bv2*Bv2) + gamma*rho0*T0/((gamma-1.d0)*Bv2) + inprod(Psi0,Psi0),Phi0) &
              +  v*Bv_parderiv(inprod(Phi0,Psi0)) + v*Bv_pbrack(inprod(Phi0,Psi0),Psi0) &
              + eta*((dx(v)*(j0y*dp(chi)/R - j0p*dy(chi)) + dy(v)*(j0p*dx(chi) - j0x*dp(chi)/R) + dp(v)*(j0x*dy(chi) - j0y*dx(chi))/R) &
              - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*gradprod(v,Psi0) - (j0x*dx(Psi0) + j0y*dy(Psi0) + j0p*dp(Psi0)/R)*Bv_parderiv(v)) &
              - k_perp*gradprod(T0,v) - T0*D_perp*gradprod(rho0,v)/(gamma-1.d0) &
              - zk_par*Bv_parderiv(v)*(Bv_parderiv(T0) + Bv_pbrack(T0,Psi0))/Bv2 + D_perp*gradprod(v*inprod(Phi0,Phi0)/Bv2,rho0)/2.d0 &
              - v*inprod(Phi0,Phi0)*S_rho/(2.d0*Bv2) + v*S_e
    end if
    
    
    amat11dt0 = (1.d0 + zeta)*v*Bv_pbrack(psi_v,Psi)
!    if (eta_T_dependent) then
!      amat11dt1 = -theta*v*Bv_pbrack(Bv_pbrack(Psi,Phi0)/Bv2,psi_v) + (theta*eta*T_0**(1.5d0))*(dx(v)*(dy(psi_v)*jp - dp(psi_v)*jy/R) &
!                + dy(v)*(dp(psi_v)*jx/R - dx(psi_v)*jp) + dp(v)*(dx(psi_v)*jy - dy(psi_v)*jx)/R)/(sqrtT0c*sqrtT0c*sqrtT0c)
!    else
      amat11dt1 = & ! -theta*v*Bv_pbrack(Bv_pbrack(Psi,Phi0)/Bv2,psi_v) + 
                  theta*eta*(dx(v)*(dy(psi_v)*jp - dp(psi_v)*jy/R) &
                + dy(v)*(dp(psi_v)*jx/R - dx(psi_v)*jp) + dp(v)*(dx(psi_v)*jy - dy(psi_v)*jx)/R)
!    end if
    amat12dt1 = -theta*v*Bv_pbrack((Bv_pbrack(Psi0,Phi) - Bv_parderiv(Phi))/Bv2,psi_v)
    if (eta_T_dependent) then
      amat14dt1 = (-1.5d0*theta*eta*T_0**(1.5d0))*T*(dx(v)*(dy(psi_v)*j0p - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*j0p) &
                + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)/(sqrtT0c*sqrtT0c*sqrtT0c*sqrtT0c*sqrtT0c)
    else
      amat14dt1 = 0.d0*(one + one)
    end if
    
    amat21dt1 = theta*(Bv2*(jx*dx(v) + jy*dy(v) + jp*dp(v)/R)/rho0 - (dx(chi)*jx + dy(chi)*jy + dp(chi)*jp/R)*Bv_parderiv(v)/rho0 &
              + (Bv2*pLap(Psi) + inprod(Psi,Bv2))*Bv_pbrack(v,Psi0)/rho0 + (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*Bv_pbrack(v,Psi)/rho0)
    amat22dt0 = (1.d0 + zeta)*v*pLap(Phi)
    if (visco_T_dependent) then
      amat22dt1 = -theta*(pLap(Phi)*Bv_pbrack(v,Phi0)/Bv2 - pLap(Phi0)*Bv_pbrack(v,Phi) + D_perp*gradprod(inprod(v,Phi)/rho0,rho0) &
                - S_rho*inprod(v,Phi)/rho0 - (visco*T_0**(1.5d0))*Lap(v)*pLap(Phi0)/(sqrtT0c*sqrtT0c*sqrtT0c))
    else
      amat22dt1 = -theta*(pLap(Phi)*Bv_pbrack(v,Phi0)/Bv2 - pLap(Phi0)*Bv_pbrack(v,Phi) + D_perp*gradprod(inprod(v,Phi)/rho0,rho0) &
                - S_rho*inprod(v,Phi)/rho0 - visco*Lap(v)*pLap(Phi0))
    end if
    amat23dt1 = theta*(-Bv2*(j0x*dx(v) + j0y*dy(v) + j0p*dp(v)/R)*rho/(rho0*rho0) &
              + (dx(chi)*j0x + dy(chi)*j0y + dp(chi)*j0p/R)*Bv_parderiv(v)*rho/(rho0*rho0) &
              - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*Bv_pbrack(v,Psi0)*rho/(rho0*rho0) &
              - D_perp*gradprod(inprod(v,Phi0)*rho/(rho0*rho0),rho0) + S_rho*inprod(v,Phi0)*rho/(rho0*rho0) &
              + D_perp*gradprod(inprod(v,Phi0)/rho0,rho) - v*Bv_pbrack(rho0,T0)*rho/(rho0*rho0) + v*Bv_pbrack(rho,T0)/rho0)
    amat24dt1 = theta*v*Bv_pbrack(rho0,T)/rho0
    
    amat32dt1 = theta*v*Bv_pbrack(rho0/Bv2,Phi)
    amat33dt0 = (1.d0 + zeta)*v*rho; amat33dt1 = theta*(v*Bv_pbrack(rho/Bv2,Phi0) + D_perp*gradprod(v,rho))
    
    amat41dt0 = (1.d0 + zeta)*v*Bv2*inprod(Psi0,Psi)
    if (eta_T_dependent) then
      amat41dt1 = theta*(v*Bv_pbrack(inprod(Psi0,Psi),Phi0) - v*Bv_parderiv(inprod(Phi0,Psi)) - v*Bv_pbrack(inprod(Phi0,Psi),Psi0) &
                - v*Bv_pbrack(inprod(Phi0,Psi0),Psi) - (eta*T_0**(1.5d0))*((dx(v)*(jy*dp(chi)/R - jp*dy(chi)) &
                + dy(v)*(jp*dx(chi) - jx*dp(chi)/R) + dp(v)*(jx*dy(chi) - jy*dx(chi))/R) &
                - (Bv2*pLap(Psi)  + inprod(Psi,Bv2))*gradprod(v,Psi0) - (jx*dx(Psi0) + jy*dy(Psi0) + jp*dp(Psi0)/R)*Bv_parderiv(v) &
                - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*gradprod(v,Psi) - (j0x*dx(Psi) + j0y*dy(Psi) + j0p*dp(Psi)/R)*Bv_parderiv(v))/(sqrtT0c*sqrtT0c*sqrtT0c) &
                + zk_par*Bv_parderiv(v)*Bv_pbrack(T0,Psi)/Bv2)
    else
      amat41dt1 = theta*(v*Bv_pbrack(inprod(Psi0,Psi),Phi0) - v*Bv_parderiv(inprod(Phi0,Psi)) - v*Bv_pbrack(inprod(Phi0,Psi),Psi0) &
                - v*Bv_pbrack(inprod(Phi0,Psi0),Psi) - eta*((dx(v)*(jy*dp(chi)/R - jp*dy(chi)) &
                + dy(v)*(jp*dx(chi) - jx*dp(chi)/R) + dp(v)*(jx*dy(chi) - jy*dx(chi))/R) &
                - (Bv2*pLap(Psi)  + inprod(Psi,Bv2))*gradprod(v,Psi0) - (jx*dx(Psi0) + jy*dy(Psi0) + jp*dp(Psi0)/R)*Bv_parderiv(v) &
                - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*gradprod(v,Psi) - (j0x*dx(Psi) + j0y*dy(Psi) + j0p*dp(Psi)/R)*Bv_parderiv(v)) &
                + zk_par*Bv_parderiv(v)*Bv_pbrack(T0,Psi)/Bv2)
    end if
    amat42dt0 = (1.d0 + zeta)*v*rho0*inprod(Phi0,Phi)/Bv2
    amat42dt1 = theta*v*(Bv_pbrack(rho0*inprod(Phi0,Phi)/(Bv2*Bv2),Phi0) + Bv_pbrack(rho0*inprod(Phi0,Phi0)/(2.d0*Bv2*Bv2),Phi) &
              - Bv_parderiv(inprod(Phi,Psi0)) - Bv_pbrack(inprod(Phi,Psi0),Psi0) - inprod(Phi0,Phi)*S_rho/Bv2) &
              - theta*D_perp*gradprod(v*inprod(Phi0,Phi)/Bv2,rho0)
    amat43dt0 = (1.d0 + zeta)*v*(rho*inprod(Phi0,Phi0)/(2.d0*Bv2) + rho*T0/(gamma-1.d0))
    amat43dt1 = theta*(v*Bv_pbrack(rho*inprod(Phi0,Phi0)/(2.d0*Bv2*Bv2) + gamma*rho*T0/((gamma-1.d0)*Bv2),Phi0) &
              + D_perp*T0*gradprod(v,rho)/(gamma-1.d0) - D_perp*gradprod(v*inprod(Phi0,Phi0)/(2.d0*Bv2),rho))
    amat44dt0 = (1.d0 + zeta)*v*rho0*T/(gamma-1.d0)
    if (eta_T_dependent) then
      amat44dt1 = theta*(v*Bv_pbrack(gamma*rho0*T/((gamma-1.d0)*Bv2),Phi0) &
                + (1.5d0*eta*T_0**(1.5d0))*T*((dx(v)*(j0y*dp(chi)/R - j0p*dy(chi)) + dy(v)*(j0p*dx(chi) - j0x*dp(chi)/R) &
                + dp(v)*(j0x*dy(chi) - j0y*dx(chi))/R) - (Bv2*pLap(Psi0) + inprod(Psi0,Bv2))*gradprod(v,Psi0) &
                - (j0x*dx(Psi0) + j0y*dy(Psi0) + j0p*dp(Psi0)/R)*Bv_parderiv(v))/(sqrtT0c*sqrtT0c*sqrtT0c*sqrtT0c*sqrtT0c) &
                - k_perp*gradprod(v,T) - D_perp*T*gradprod(v,rho0)/(gamma-1.d0) &
                - zk_par*Bv_parderiv(v)*(Bv_parderiv(T) + Bv_pbrack(T,Psi0))/Bv2)
    else
      amat44dt1 = theta*(v*Bv_pbrack(gamma*rho0*T/((gamma-1.d0)*Bv2),Phi0) &
                - k_perp*gradprod(v,T) - D_perp*T*gradprod(v,rho0)/(gamma-1.d0) &
                - zk_par*Bv_parderiv(v)*(Bv_parderiv(T) + Bv_pbrack(T,Psi0))/Bv2)
    end if
    
    rhs1dt1e = Dexpand(deepcopy(rhs1dt1))
    rhs2dt1e = Dexpand(deepcopy(rhs2dt1))
    rhs3dt1e = Dexpand(deepcopy(rhs3dt1))
    rhs4dt1e = Dexpand(deepcopy(rhs4dt1))
    
    amat11dt1e = Dexpand(deepcopy(amat11dt1)); amat12dt1e = Dexpand(deepcopy(amat12dt1))
    amat21dt1e = Dexpand(deepcopy(amat21dt1)); amat22dt1e = Dexpand(deepcopy(amat22dt1)); amat23dt1e = Dexpand(deepcopy(amat23dt1))
    amat32dt1e = Dexpand(deepcopy(amat32dt1)); amat33dt1e = Dexpand(deepcopy(amat33dt1))
    amat41dt1e = Dexpand(deepcopy(amat41dt1)); amat42dt1e = Dexpand(deepcopy(amat42dt1)); amat43dt1e = Dexpand(deepcopy(amat43dt1))
    amat44dt1e = Dexpand(deepcopy(amat44dt1))
    
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
      allocate(thread_eq(i)%rhs3dt1seq(countsubexprs(rhs3dt1e)))
      last = 0; call buildsequence(rhs3dt1e, thread_eq(i)%rhs3dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%rhs4dt0seq(countsubexprs(rhs4dt0)))
      last = 0; call buildsequence(rhs4dt0, thread_eq(i)%rhs4dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%rhs4dt1seq(countsubexprs(rhs4dt1e)))
      last = 0; call buildsequence(rhs4dt1e, thread_eq(i)%rhs4dt1seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat11dt0seq(countsubexprs(amat11dt0)))
      last = 0; call buildsequence(amat11dt0, thread_eq(i)%amat11dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat11dt1seq(countsubexprs(amat11dt1e)))
      last = 0; call buildsequence(amat11dt1e, thread_eq(i)%amat11dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat12dt1seq(countsubexprs(amat12dt1e)))
      last = 0; call buildsequence(amat12dt1e, thread_eq(i)%amat12dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat14dt1seq(countsubexprs(amat14dt1)))
      last = 0; call buildsequence(amat14dt1, thread_eq(i)%amat14dt1seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat21dt1seq(countsubexprs(amat21dt1e)))
      last = 0; call buildsequence(amat21dt1e, thread_eq(i)%amat21dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat22dt0seq(countsubexprs(amat22dt0)))
      last = 0; call buildsequence(amat22dt0, thread_eq(i)%amat22dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat22dt1seq(countsubexprs(amat22dt1e)))
      last = 0; call buildsequence(amat22dt1e, thread_eq(i)%amat22dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat23dt1seq(countsubexprs(amat23dt1e)))
      last = 0; call buildsequence(amat23dt1e, thread_eq(i)%amat23dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat24dt1seq(countsubexprs(amat24dt1)))
      last = 0; call buildsequence(amat24dt1, thread_eq(i)%amat24dt1seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat32dt1seq(countsubexprs(amat32dt1e)))
      last = 0; call buildsequence(amat32dt1e, thread_eq(i)%amat32dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat33dt0seq(countsubexprs(amat33dt0)))
      last = 0; call buildsequence(amat33dt0, thread_eq(i)%amat33dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat33dt1seq(countsubexprs(amat33dt1e)))
      last = 0; call buildsequence(amat33dt1e, thread_eq(i)%amat33dt1seq, thread_eq(i)%eq, last)
      
      
      allocate(thread_eq(i)%amat41dt0seq(countsubexprs(amat41dt0)))
      last = 0; call buildsequence(amat41dt0, thread_eq(i)%amat41dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat41dt1seq(countsubexprs(amat41dt1e)))
      last = 0; call buildsequence(amat41dt1e, thread_eq(i)%amat41dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat42dt0seq(countsubexprs(amat42dt0)))
      last = 0; call buildsequence(amat42dt0, thread_eq(i)%amat42dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat42dt1seq(countsubexprs(amat42dt1e)))
      last = 0; call buildsequence(amat42dt1e, thread_eq(i)%amat42dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat43dt0seq(countsubexprs(amat43dt0)))
      last = 0; call buildsequence(amat43dt0, thread_eq(i)%amat43dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat43dt1seq(countsubexprs(amat43dt1e)))
      last = 0; call buildsequence(amat43dt1e, thread_eq(i)%amat43dt1seq, thread_eq(i)%eq, last)
      
      allocate(thread_eq(i)%amat44dt0seq(countsubexprs(amat44dt0)))
      last = 0; call buildsequence(amat44dt0, thread_eq(i)%amat44dt0seq, thread_eq(i)%eq, last)
      allocate(thread_eq(i)%amat44dt1seq(countsubexprs(amat44dt1e)))
      last = 0; call buildsequence(amat44dt1e, thread_eq(i)%amat44dt1seq, thread_eq(i)%eq, last)
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