module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads


  implicit none
  
  type type_thread_eq
#ifdef DEBUG
    type(action), dimension(:), allocatable :: rhs1seq, rhs2seq, rhs3seq, rhs4seq, rhs5seq, rhs6seq
    type(action), dimension(:), allocatable :: amat11seq, amat12seq, amat13seq, amat16seq
    type(action), dimension(:), allocatable :: amat21seq, amat22seq, amat23seq, amat24seq, amat25seq, amat26seq
    type(action), dimension(:), allocatable :: amat31seq, amat33seq
    type(action), dimension(:), allocatable :: amat42seq, amat44seq
    type(action), dimension(:), allocatable :: amat51seq, amat52seq, amat55seq
    type(action), dimension(:), allocatable :: amat61seq, amat62seq, amat63seq, amat65seq, amat66seq, amat67seq
    type(action), dimension(:), allocatable :: amat71seq, amat72seq, amat73seq, amat75seq, amat76seq, amat77seq
    type(action), dimension(:), allocatable :: aBv2seq, aBv2xseq, aBv2yseq, aBv2pseq, aB2seq
#endif

    ! Indices in eq array: variable index (see algexpr's below), R derivative order, z derivative order, phi derivative order,
    !   separation of terms with covariant phi derivatives in test function and unknown (FFT)
    real*8, dimension(:,:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  ! Variables at current time step
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: zj0        = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: w0         = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: T0_i       = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: T0_e       = algexpr(basic=.true.,var=7)
  ! Changes since previous time step
  type(algexpr), parameter, private :: delta_Psi  = algexpr(basic=.true.,var=n_var+1)
  type(algexpr), parameter, private :: delta_Phi  = algexpr(basic=.true.,var=n_var+2)
  type(algexpr), parameter, private :: delta_zj   = algexpr(basic=.true.,var=n_var+3)
  type(algexpr), parameter, private :: delta_w    = algexpr(basic=.true.,var=n_var+4)
  type(algexpr), parameter, private :: delta_rho  = algexpr(basic=.true.,var=n_var+5)
  type(algexpr), parameter, private :: delta_T    = algexpr(basic=.true.,var=n_var+6)
  type(algexpr), parameter, private :: delta_T_i  = algexpr(basic=.true.,var=n_var+6)
  type(algexpr), parameter, private :: delta_T_e  = algexpr(basic=.true.,var=n_var+7)  
  ! Test function
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=2*n_var+1)
  ! Unknowns
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: zj         = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: w          = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: T_i        = algexpr(basic=.true.,var=2*n_var+2)
  type(algexpr), parameter, private :: T_e        = algexpr(basic=.true.,var=2*n_var+2)
  ! Other quantities
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=15)
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=16)
  type(algexpr), parameter, private :: D_perp     = algexpr(basic=.true.,var=17)
  type(algexpr), parameter, private :: S_rho      = algexpr(basic=.true.,var=18)
  type(algexpr), parameter, private :: S_j        = algexpr(basic=.true.,var=19)
  type(algexpr), parameter, private :: eta        = algexpr(basic=.true.,var=20)
  type(algexpr), parameter, private :: deta_dT    = algexpr(basic=.true.,var=21)
  type(algexpr), parameter, private :: visco      = algexpr(basic=.true.,var=22)
  type(algexpr), parameter, private :: dvisco_dT  = algexpr(basic=.true.,var=23)
  type(algexpr), parameter, private :: k_perp     = algexpr(basic=.true.,var=24)
  type(algexpr), parameter, private :: k_perp_i   = algexpr(basic=.true.,var=25)
  type(algexpr), parameter, private :: k_perp_e   = algexpr(basic=.true.,var=26)
  type(algexpr), parameter, private :: S_e        = algexpr(basic=.true.,var=27)
  type(algexpr), parameter, private :: S_e_i      = algexpr(basic=.true.,var=28)
  type(algexpr), parameter, private :: S_e_e      = algexpr(basic=.true.,var=29)
  type(algexpr), parameter, private :: k_par      = algexpr(basic=.true.,var=30) 
  type(algexpr), parameter, private :: k_par_i    = algexpr(basic=.true.,var=31)
  type(algexpr), parameter, private :: k_par_e    = algexpr(basic=.true.,var=32)
  type(algexpr), parameter, private :: dk_par_dT  = algexpr(basic=.true.,var=33)
  type(algexpr), parameter, private :: dk_par_dT_i= algexpr(basic=.true.,var=34)
  type(algexpr), parameter, private :: dk_par_dT_e= algexpr(basic=.true.,var=35)
  type(algexpr), parameter, private :: dTe_i      = algexpr(basic=.true.,var=36)
  type(algexpr), parameter, private :: ddTe_i_dT_i= algexpr(basic=.true.,var=37)
  type(algexpr), parameter, private :: ddTe_i_dT_e= algexpr(basic=.true.,var=38)
  type(algexpr), parameter, private :: ddTe_i_drho= algexpr(basic=.true.,var=39)

  ! Auxiliary variables (aux)
  type(algexpr), parameter, private :: Bv2        = algexpr(basic=.true.,var=40)
  type(algexpr), parameter, private :: B2         = algexpr(basic=.true.,var=41)
  type(const), private :: tstep, zeta, theta, visco_num, eta_num, D_perp_num, k_perp_num, gamma, reta
  
  type(algexpr), private :: rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7
  type(algexpr), private :: amat11, amat12, amat13, amat16,  amat17
  type(algexpr), private :: amat21, amat22, amat23, amat24, amat25, amat26, amat27
  type(algexpr), private :: amat31, amat33
  type(algexpr), private :: amat42, amat44
  type(algexpr), private :: amat51, amat52, amat55
  type(algexpr), private :: amat61, amat62, amat63, amat65, amat66, amat67
  type(algexpr), private :: amat71, amat72, amat73, amat75, amat76, amat77
  type(algexpr), private :: a_Bv2, a_B2
  
  integer            :: n_rhs, n_amat
  integer, parameter :: n_aux  = 5
  
  type(algexpr), private :: rhs2e, rhs3e, rhs4e, rhs5e, rhs6e, rhs7e
  type(algexpr), private :: amat22e, amat24e, amat25e, amat26e, amat27e
  type(algexpr), private :: amat31e
  type(algexpr), private :: amat52e, amat55e
  type(algexpr), private :: amat62e, amat65e, amat66e, amat67e
  type(algexpr), private :: amat72e, amat75e, amat76e, amat77e
  type(algexpr), private :: ea_Bv2x, ea_Bv2y, ea_Bv2p

  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains

  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, Igamma => gamma, Itstep => tstep, Ivisco_num => visco_num, Ieta_num => eta_num, &
                           ID_perp_num => D_perp_num, zk_perp_num, Ieta => eta, eta_ohmic
    implicit none

    tstep      = const(value = Itstep,          token = "tstep")
    zeta       = const(value = time_evol_zeta,  token = "zeta")
    theta      = const(value = time_evol_theta, token = "theta")
    visco_num  = const(value = Ivisco_num,      token = "visco_num")
    eta_num    = const(value = Ieta_num,        token = "eta_num")
    D_perp_num = const(value = ID_perp_num,     token = "D_perp_num")
    k_perp_num = const(value = zk_perp_num,     token = "zk_perp_num")
    gamma      = const(value = Igamma,          token = "gamma")
    if (Ieta .ne. 0.d0) then
      reta     = const(value = eta_ohmic/Ieta,  token = "reta")
    else
      reta     = const(value = 0.d0,            token = "reta")
    end if
    
    a_Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)
    a_B2 = Bv2 + Bv2*inprod(Psi0,Psi0)
    
    if (with_TiTe) then
      
      ! --- RHS
      rhs1 = tstep*v*((Bv_parderiv(Phi0) - Bv_pbrack(Psi0,Phi0))/Bv2 + eta*(zj0 - S_j)) + tstep*eta_num*inprod(v,zj0) + zeta*v*delta_Psi
      
      rhs2 = -tstep*((Bv_pbrack(rho0/Bv2,v)*inprod(Phi0,Phi0)/2.d0 - Bv_pbrack(v,Phi0)*rho0*w0/Bv2 - Bv_pbrack(rho0/Bv2,Phi0)*inprod(v,Phi0) &
           + Bv_pbrack(v,rho0*(T0_i+T0_e)))/Bv2 - v*Bv_parderiv(zj0) - v*Bv_pbrack(zj0,Psi0) + visco*inprod(v,w0) + visco_num*Lap(v)*Lap(w0)) &
           - zeta*(rho0*inprod(v,delta_Phi) + delta_rho*inprod(v,Phi0))/Bv2

      rhs3 = -Bv2*inprod(v,Psi0) - v*Bv2*zj0

      rhs4 = -inprod(v,Phi0) - v*w0

      rhs5 = -tstep*(v*Bv_pbrack(rho0/Bv2,Phi0) + D_perp*gradgrad_perp(v,rho0) - S_rho*v) + zeta*v*delta_rho

      rhs6 = -tstep*(v*Bv_pbrack(rho0*T0_i, Phi0)/Bv2 - gamma*v*rho0*T0_i*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) + k_perp_i*gradprod(v,T0_i) &
           + (k_par_i-k_perp_i)*B0_parderiv(v)*B0_parderiv(T0_i)/B2 + k_perp_num*Lap(v)*Lap(T0_i) + D_perp*T0_i*gradgrad_perp(v, rho0) &
           - (gamma - 1.d0)*reta*eta*v*Bv2*zj0*zj0 - v*S_e) + zeta*v*(rho0*delta_T_i + T0_i*delta_rho) &
           + v*dTe_i
      rhs7 = -tstep*(v*Bv_pbrack(rho0*T0_e, Phi0)/Bv2 - gamma*v*rho0*T0_e*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) + k_perp_i*gradprod(v,T0_e) &
           + (k_par_e-k_perp_e)*B0_parderiv(v)*B0_parderiv(T0_e)/B2 + k_perp_num*Lap(v)*Lap(T0_e) + D_perp*T0_e*gradgrad_perp(v, rho0) &
           - (gamma - 1.d0)*reta*eta*v*Bv2*zj0*zj0 - v*S_e) + zeta*v*(rho0*delta_T_e + T0_e*delta_rho) &
           - v*dTe_i
     
      ! --- LHS
      amat11 = (1.d0 + zeta)*v*Psi + tstep*theta*v*Bv_pbrack(Psi,Phi0)/Bv2
      amat12 = (-tstep*theta)*v*(Bv_parderiv(Phi) - Bv_pbrack(Psi0,Phi))/Bv2
      amat13 = (-tstep*theta)*(eta*v*zj + eta_num*inprod(v,zj))
      amat17 = (-tstep*theta)*v*deta_dT*T_e*zj0

      amat21 = (-tstep*theta)*v*Bv_pbrack(zj0,Psi)
      amat22 = -(1.d0 + zeta)*rho0*inprod(v,Phi)/Bv2 + tstep*theta*(Bv_pbrack(rho0/Bv2,v)*inprod(Phi0,Phi) - rho0*w0*Bv_pbrack(v,Phi)/Bv2 &
             - Bv_pbrack(rho0/Bv2,Phi)*inprod(v,Phi0) - Bv_pbrack(rho0/Bv2,Phi0)*inprod(v,Phi))/Bv2
      amat23 = (-tstep*theta)*v*(Bv_parderiv(zj) + Bv_pbrack(zj,Psi0))
      amat24 = -tstep*theta*rho0*w*Bv_pbrack(v,Phi0)/(Bv2*Bv2) + tstep*theta*(visco*inprod(v,w) + visco_num*Lap(v)*Lap(w))
      amat25 = -(1.d0 + zeta)*rho*inprod(v,Phi0)/Bv2 + tstep*theta*(Bv_pbrack(rho/Bv2,v)*inprod(Phi0,Phi0)/2.d0 - rho*w0*Bv_pbrack(v,Phi0)/Bv2 &
             - Bv_pbrack(rho/Bv2,Phi0)*inprod(v,Phi0) + Bv_pbrack(v,rho*(T0_i+T0_e)))/Bv2
      amat26 = tstep*theta*Bv_pbrack(v,rho0*T_i)/Bv2 
      amat27 = tstep*theta*Bv_pbrack(v,rho0*T_e)/Bv2 + tstep*theta*dvisco_dT*T_e*inprod(v,w0)

      amat31 = theta*Bv2*inprod(v,Psi)
      amat33 = theta*v*Bv2*zj

#ifdef DEBUG
      amat42 = theta*inprod(v,Phi) + 0.d0*one
#else
      amat42 = theta*inprod(v,Phi)
#endif
      amat44 = theta*v*w

      amat51 = (-tstep*theta)*D_perp*gradDgrad_par(v,rho0)
      amat52 = tstep*theta*v*Bv_pbrack(rho0/Bv2,Phi)
      amat55 = (1.d0 + zeta)*v*rho + tstep*theta*(v*Bv_pbrack(rho/Bv2,Phi0) + D_perp*gradgrad_perp(v,rho))
   
      amat61 = tstep*theta*((k_par_i - k_perp_i)*gradDgrad_par(v,T0_i) - D_perp*T0_i*gradDgrad_par(v,rho0))
      amat62 = tstep*theta*v*(Bv_pbrack(rho0*T0_i,Phi) - gamma*rho0*T0_i*Bv_pbrack(Bv2,Phi)/Bv2)/Bv2
      amat63 = -2.d0*tstep*theta*(gamma - 1.d0)*v*reta*eta*Bv2*zj0*zj
      amat65 = (1.d0 + zeta)*v*rho*T0_i + tstep*theta*(v*Bv_pbrack(rho*T0_i,Phi0)/Bv2 - gamma*v*rho*T0_i*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
             + D_perp*T0_i*gradgrad_perp(v,rho)) + ddTe_i_drho
      amat66 = (1.d0 + zeta)*v*rho0*T_i + tstep*theta*(v*Bv_pbrack(rho0*T_i,Phi0)/Bv2 - gamma*v*rho0*T_i*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
             + k_perp_i*gradprod(v,T_i) + (k_par_i - k_perp_i)*B0_parderiv(v)*B0_parderiv(T)/B2 &
             + dk_par_dT_i*T_i*B0_parderiv(v)*B0_parderiv(T0_i)/B2 + k_perp_num*Lap(v)*Lap(T_i) &
             + D_perp*T_i*gradgrad_perp(v,rho0) + ddTe_i_dT_i)
      amat67 = tstep*theta*ddTe_i_dT_e

      amat71 = tstep*theta*((k_par_e - k_perp_e)*gradDgrad_par(v,T0_e) - D_perp*T0_e*gradDgrad_par(v,rho0))
      amat72 = tstep*theta*v*(Bv_pbrack(rho0*T0_e,Phi) - gamma*rho0*T0_e*Bv_pbrack(Bv2,Phi)/Bv2)/Bv2
      amat73 = -2.d0*tstep*theta*(gamma - 1.d0)*v*reta*eta*Bv2*zj0*zj
      amat75 = (1.d0 + zeta)*v*rho*T0_e + tstep*theta*(v*Bv_pbrack(rho*T0_e,Phi0)/Bv2 - gamma*v*rho*T0_e*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
             + D_perp*T0_e*gradgrad_perp(v,rho))
      amat76 = - tstep * theta * ddTe_i_dT_i
      amat77 = (1.d0 + zeta)*v*rho0*T_e + tstep*theta*(v*Bv_pbrack(rho0*T_e,Phi0)/Bv2 - gamma*v*rho0*T_e*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
             + k_perp_e*gradprod(v,T_e) + (k_par_e - k_perp_e)*B0_parderiv(v)*B0_parderiv(T_e)/B2 &
             + dk_par_dT_e*T_e*B0_parderiv(v)*B0_parderiv(T0_e)/B2 + k_perp_num*Lap(v)*Lap(T_e) &
             + D_perp*T_e*gradgrad_perp(v,rho0) - v*reta*deta_dT*T_e*Bv2*zj0*zj0 - ddTe_i_dT_e)

      rhs2e = Dexpand(deepcopy(rhs2))
      rhs5e = Dexpand(deepcopy(rhs5))
      rhs6e = Dexpand(deepcopy(rhs6))
      rhs7e = Dexpand(deepcopy(rhs7))

      amat22e = Dexpand(deepcopy(amat22)); amat24e = Dexpand(deepcopy(amat24)); amat25e = Dexpand(deepcopy(amat25)); 
      amat26e = Dexpand(deepcopy(amat26)); amat27e = Dexpand(deepcopy(amat27))
      amat31e = Dexpand(deepcopy(amat31))
      amat52e = Dexpand(deepcopy(amat52)); amat55e = Dexpand(deepcopy(amat55))
      amat62e = Dexpand(deepcopy(amat62)); amat65e = Dexpand(deepcopy(amat65)); amat66e = Dexpand(deepcopy(amat66));  
      amat67e = Dexpand(deepcopy(amat67))
      amat72e = Dexpand(deepcopy(amat72)); amat75e = Dexpand(deepcopy(amat75)); amat76e = Dexpand(deepcopy(amat76));
      amat77e = Dexpand(deepcopy(amat77))
    else     
      ! --- RHS
      rhs1 = tstep*v*((Bv_parderiv(Phi0) - Bv_pbrack(Psi0,Phi0))/Bv2 + eta*(zj0 - S_j)) + tstep*eta_num*inprod(v,zj0) + zeta*v*delta_Psi

      rhs2 = -tstep*((Bv_pbrack(rho0/Bv2,v)*inprod(Phi0,Phi0)/2.d0 - Bv_pbrack(v,Phi0)*rho0*w0/Bv2 - Bv_pbrack(rho0/Bv2,Phi0)*inprod(v,Phi0) &
           + Bv_pbrack(v,rho0*T0))/Bv2 - v*Bv_parderiv(zj0) - v*Bv_pbrack(zj0,Psi0) + visco*inprod(v,w0) + visco_num*Lap(v)*Lap(w0)) &
           - zeta*(rho0*inprod(v,delta_Phi) + delta_rho*inprod(v,Phi0))/Bv2

      rhs3 = -Bv2*inprod(v,Psi0) - v*Bv2*zj0

      rhs4 = -inprod(v,Phi0) - v*w0

      rhs5 = -tstep*(v*Bv_pbrack(rho0/Bv2,Phi0) + D_perp*gradgrad_perp(v,rho0) - S_rho*v) + zeta*v*delta_rho

      rhs6 = -tstep*(v*Bv_pbrack(rho0*T0,Phi0)/Bv2 - gamma*v*rho0*T0*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) + k_perp*gradprod(v,T0) &
           + (k_par - k_perp)*B0_parderiv(v)*B0_parderiv(T0)/B2 + k_perp_num*Lap(v)*Lap(T0) + D_perp*T0*gradgrad_perp(v,rho0) &
           - (gamma - 1.d0)*reta*eta*v*Bv2*zj0*zj0 - v*S_e) + zeta*v*(rho0*delta_T + T0*delta_rho)

      ! --- RHS
      amat11 = (1.d0 + zeta)*v*Psi + tstep*theta*v*Bv_pbrack(Psi,Phi0)/Bv2
      amat12 = (-tstep*theta)*v*(Bv_parderiv(Phi) - Bv_pbrack(Psi0,Phi))/Bv2
      amat13 = (-tstep*theta)*(eta*v*zj + eta_num*inprod(v,zj))
      amat16 = (-tstep*theta)*v*deta_dT*T*zj0

      amat21 = (-tstep*theta)*v*Bv_pbrack(zj0,Psi)
      amat22 = -(1.d0 + zeta)*rho0*inprod(v,Phi)/Bv2 + tstep*theta*(Bv_pbrack(rho0/Bv2,v)*inprod(Phi0,Phi) - rho0*w0*Bv_pbrack(v,Phi)/Bv2 &
             - Bv_pbrack(rho0/Bv2,Phi)*inprod(v,Phi0) - Bv_pbrack(rho0/Bv2,Phi0)*inprod(v,Phi))/Bv2
      amat23 = (-tstep*theta)*v*(Bv_parderiv(zj) + Bv_pbrack(zj,Psi0))
      amat24 = -tstep*theta*rho0*w*Bv_pbrack(v,Phi0)/(Bv2*Bv2) + tstep*theta*(visco*inprod(v,w) + visco_num*Lap(v)*Lap(w))
      amat25 = -(1.d0 + zeta)*rho*inprod(v,Phi0)/Bv2 + tstep*theta*(Bv_pbrack(rho/Bv2,v)*inprod(Phi0,Phi0)/2.d0 - rho*w0*Bv_pbrack(v,Phi0)/Bv2 &
             - Bv_pbrack(rho/Bv2,Phi0)*inprod(v,Phi0) + Bv_pbrack(v,rho*T0))/Bv2
      amat26 = tstep*theta*Bv_pbrack(v,rho0*T)/Bv2 + tstep*theta*dvisco_dT*T*inprod(v,w0)

      amat31 = theta*Bv2*inprod(v,Psi)
      amat33 = theta*v*Bv2*zj

#ifdef DEBUG
      amat42 = theta*inprod(v,Phi) + 0.d0*one
#else
      amat42 = theta*inprod(v,Phi)
#endif
      amat44 = theta*v*w

      amat51 = (-tstep*theta)*D_perp*gradDgrad_par(v,rho0)
      amat52 = tstep*theta*v*Bv_pbrack(rho0/Bv2,Phi)
      amat55 = (1.d0 + zeta)*v*rho + tstep*theta*(v*Bv_pbrack(rho/Bv2,Phi0) + D_perp*gradgrad_perp(v,rho))

      amat61 = tstep*theta*((k_par - k_perp)*gradDgrad_par(v,T0) - D_perp*T0*gradDgrad_par(v,rho0))
      amat62 = tstep*theta*v*(Bv_pbrack(rho0*T0,Phi) - gamma*rho0*T0*Bv_pbrack(Bv2,Phi)/Bv2)/Bv2
      amat63 = -2.d0*tstep*theta*(gamma - 1.d0)*v*reta*eta*Bv2*zj0*zj
      amat65 = (1.d0 + zeta)*v*rho*T0 + tstep*theta*(v*Bv_pbrack(rho*T0,Phi0)/Bv2 - gamma*v*rho*T0*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
             + D_perp*T0*gradgrad_perp(v,rho))
      amat66 = (1.d0 + zeta)*v*rho0*T + tstep*theta*(v*Bv_pbrack(rho0*T,Phi0)/Bv2 - gamma*v*rho0*T*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
             + k_perp*gradprod(v,T) + (k_par - k_perp)*B0_parderiv(v)*B0_parderiv(T)/B2 + dk_par_dT*T*B0_parderiv(v)*B0_parderiv(T0)/B2 &
             + k_perp_num*Lap(v)*Lap(T) + D_perp*T*gradgrad_perp(v,rho0) - v*reta*deta_dT*T*Bv2*zj0*zj0)

      rhs2e = Dexpand(deepcopy(rhs2))
      rhs5e = Dexpand(deepcopy(rhs5))
      rhs6e = Dexpand(deepcopy(rhs6))

      amat22e = Dexpand(deepcopy(amat22)); amat24e = Dexpand(deepcopy(amat24)); amat25e = Dexpand(deepcopy(amat25)); amat26e = Dexpand(deepcopy(amat26))
      amat31e = Dexpand(deepcopy(amat31))
      amat52e = Dexpand(deepcopy(amat52)); amat55e = Dexpand(deepcopy(amat55))
      amat62e = Dexpand(deepcopy(amat62)); amat65e = Dexpand(deepcopy(amat65)); amat66e = Dexpand(deepcopy(amat66))
    end if
  
    ea_Bv2x = Dexpand(deepcopy(dx(a_Bv2))); ea_Bv2y = Dexpand(deepcopy(dy(a_Bv2))); ea_Bv2p = Dexpand(deepcopy(dp(a_Bv2)))
  end subroutine init_equations
  
  subroutine init_eq_struct()
    use data_structure, only: nbthreads
    implicit none
    integer :: i
    
    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        if (with_TiTe) then !> better to reduce for no TiTe?
          allocate(thread_eq(i)%eq(41,0:n_order-1,0:n_order-1,0:n_order-1,4))
        else
          allocate(thread_eq(i)%eq(41,0:n_order-1,0:n_order-1,0:n_order-1,4))
        end if
#ifdef DEBUG
        allocate(thread_eq(i)%rhs1seq(countsubexprs(rhs1)))
        allocate(thread_eq(i)%rhs2seq(countsubexprs(rhs2e)))
        allocate(thread_eq(i)%rhs3seq(countsubexprs(rhs3)))
        allocate(thread_eq(i)%rhs4seq(countsubexprs(rhs4)))
        allocate(thread_eq(i)%rhs5seq(countsubexprs(rhs5e)))
        allocate(thread_eq(i)%rhs6seq(countsubexprs(rhs6e)))
        allocate(thread_eq(i)%amat11seq(countsubexprs(amat11)))
        allocate(thread_eq(i)%amat12seq(countsubexprs(amat12)))
        allocate(thread_eq(i)%amat13seq(countsubexprs(amat13)))
        allocate(thread_eq(i)%amat21seq(countsubexprs(amat21)))
        allocate(thread_eq(i)%amat22seq(countsubexprs(amat22e)))
        allocate(thread_eq(i)%amat23seq(countsubexprs(amat23)))
        allocate(thread_eq(i)%amat24seq(countsubexprs(amat24e)))
        allocate(thread_eq(i)%amat25seq(countsubexprs(amat25e)))
        allocate(thread_eq(i)%amat26seq(countsubexprs(amat26e)))
        allocate(thread_eq(i)%amat31seq(countsubexprs(amat31e)))
        allocate(thread_eq(i)%amat33seq(countsubexprs(amat33)))
        allocate(thread_eq(i)%amat42seq(countsubexprs(amat42)))
        allocate(thread_eq(i)%amat44seq(countsubexprs(amat44)))
        allocate(thread_eq(i)%amat51seq(countsubexprs(amat51)))
        allocate(thread_eq(i)%amat52seq(countsubexprs(amat52e)))
        allocate(thread_eq(i)%amat55seq(countsubexprs(amat55e)))
        allocate(thread_eq(i)%amat61seq(countsubexprs(amat61)))
        allocate(thread_eq(i)%amat62seq(countsubexprs(amat62e)))
        allocate(thread_eq(i)%amat63seq(countsubexprs(amat63)))
        allocate(thread_eq(i)%amat65seq(countsubexprs(amat65e)))
        allocate(thread_eq(i)%amat66seq(countsubexprs(amat66e)))
        allocate(thread_eq(i)%aBv2seq(countsubexprs(a_Bv2)))
        allocate(thread_eq(i)%aBv2xseq(countsubexprs(ea_Bv2x)))
        allocate(thread_eq(i)%aBv2yseq(countsubexprs(ea_Bv2y)))
        allocate(thread_eq(i)%aBv2pseq(countsubexprs(ea_Bv2p)))
        allocate(thread_eq(i)%aB2seqtype(algexpr), allocatable, intent(out) :: rhs(:)
    character(8),  allocatable, intent(out) :: varnames(:)(countsubexprs(a_B2)))
        if (with_TiTe) then
          allocate(thread_eq(i)%rhs7seq(countsubexprs(rhs7e)))
          allocate(thread_eq(i)%amat17seq(countsubexprs(amat17)))
          allocate(thread_eq(i)%amat27seq(countsubexprs(amat27e)))
          allocate(thread_eq(i)%amat67seq(countsubexprs(amat67e)))
          allocate(thread_eq(i)%amat71seq(countsubexprs(amat71)))
          allocate(thread_eq(i)%amat72seq(countsubexprs(amat72e)))
          allocate(thread_eq(i)%amat73seq(countsubexprs(amat73)))
          allocate(thread_eq(i)%amat75seq(countsubexprs(amat75e)))
          allocate(thread_eq(i)%amat76seq(countsubexprs(amat76e)))
          allocate(thread_eq(i)%amat77seq(countsubexprs(amat77e)))
        end if
#endif
      end do
    end if
  end subroutine init_eq_struct
  
#ifdef DEBUG
  subroutine build_all_seq()
    use data_structure, only: nbthreads
    implicit none
    integer :: i
    
    do i=1,nbthreads
      call buildsequence(rhs1, thread_eq(i)%rhs1seq, thread_eq(i)%eq)
      call buildsequence(rhs2e, thread_eq(i)%rhs2seq, thread_eq(i)%eq)
      call buildsequence(rhs3,  thread_eq(i)%rhs3seq, thread_eq(i)%eq)
      call buildsequence(rhs4,  thread_eq(i)%rhs4seq, thread_eq(i)%eq)
      call buildsequence(rhs5e, thread_eq(i)%rhs5seq, thread_eq(i)%eq)
      call buildsequence(rhs6e, thread_eq(i)%rhs6seq, thread_eq(i)%eq)

      call buildsequence(amat11, thread_eq(i)%amat11seq, thread_eq(i)%eq)
      call buildsequence(amat12, thread_eq(i)%amat12seq, thread_eq(i)%eq)
      call buildsequence(amat13, thread_eq(i)%amat13seq, thread_eq(i)%eq)
      if (with_TiTe) then
        call buildsequence(rhs7e, thread_eq(i)%rhs7seq, thread_eq(i)%eq)
        call buildsequence(amat17, thread_eq(i)%amat17seq, thread_eq(i)%eq)
      else
        call buildsequence(amat16, thread_eq(i)%amat16seq, thread_eq(i)%eq)
      end if

      call buildsequence(amat21, thread_eq(i)%amat21seq, thread_eq(i)%eq)
      call buildsequence(amat22e, thread_eq(i)%amat22seq, thread_eq(i)%eq)
      call buildsequence(amat23, thread_eq(i)%amat23seq, thread_eq(i)%eq)
      call buildsequence(amat24e, thread_eq(i)%amat24seq, thread_eq(i)%eq)
      call buildsequence(amat25e, thread_eq(i)%amat25seq, thread_eq(i)%eq)
      call buildsequence(amat26e, thread_eq(i)%amat26seq, thread_eq(i)%eq)
      call buildsequence(amat33, thread_eq(i)%amat33seq, thread_eq(i)%eq)
      
      call buildsequence(amat42, thread_eq(i)%amat42seq, thread_eq(i)%eq)
      call buildsequence(amat44, thread_eq(i)%amat44seq, thread_eq(i)%eq)
      
      call buildsequence(amat51,  thread_eq(i)%amat51seq, thread_eq(i)%eq)
      call buildsequence(amat52e, thread_eq(i)%amat52seq, thread_eq(i)%eq)
      call buildsequence(amat55e, thread_eq(i)%amat55seq, thread_eq(i)%eq)
      
      call buildsequence(amat61, thread_eq(i)%amat61seq, thread_eq(i)%eq)
      call buildsequence(amat62e, thread_eq(i)%amat62seq, thread_eq(i)%eq)
      call buildsequence(amat63, thread_eq(i)%amat63seq, thread_eq(i)%eq)
      call buildsequence(amat65e, thread_eq(i)%amat65seq, thread_eq(i)%eq)
      call buildsequence(amat66e, thread_eq(i)%amat66seq, thread_eq(i)%eq)

      if (with_TiTe) then
        call buildsequence(amat27e, thread_eq(i)%amat27seq, thread_eq(i)%eq)
        call buildsequence(amat67e, thread_eq(i)%amat67seq, thread_eq(i)%eq)
        call buildsequence(amat71, thread_eq(i)%amat71seq, thread_eq(i)%eq)
        call buildsequence(amat72e, thread_eq(i)%amat72seq, thread_eq(i)%eq)
        call buildsequence(amat73, thread_eq(i)%amat73seq, thread_eq(i)%eq)
        call buildsequence(amat75e, thread_eq(i)%amat75seq, thread_eq(i)%eq)
        call buildsequence(amat76e, thread_eq(i)%amat76seq, thread_eq(i)%eq)
        call buildsequence(amat77e, thread_eq(i)%amat77seq, thread_eq(i)%eq)
      end if

      call buildsequence(a_Bv2, thread_eq(i)%aBv2seq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2x, thread_eq(i)%aBv2xseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2y, thread_eq(i)%aBv2yseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2p, thread_eq(i)%aBv2pseq, thread_eq(i)%eq)
      
      call buildsequence(a_B2, thread_eq(i)%aB2seq, thread_eq(i)%eq)
    end do
  end subroutine build_all_seq
#endif
  
  subroutine get_rhs(rhs,varnames)
    implicit none
    
!    type(algexpr), dimension(n_rhs), intent(out) :: rhs
!    character(8),  dimension(n_rhs), intent(out) :: varnames

    type(algexpr), allocatable, intent(out) :: rhs(:)
    character(8),  allocatable, intent(out) :: varnames(:)
    if (with_TiTe) then
      n_rhs = 7
      allocate(rhs(n_rhs), varnames(n_rhs)) 
      rhs = (/ rhs1, rhs2e, rhs3, rhs4, rhs5e, rhs6e, rhs7e /)
      varnames = (/ "rhs_ij_1", "rhs_ij_2", "rhs_ij_3", "rhs_ij_4", "rhs_ij_5", "rhs_ij_6", "rhs_ij_7" /)
    else    
      n_rhs = 6
      allocate(rhs(n_rhs), varnames(n_rhs))
      rhs = (/ rhs1, rhs2e, rhs3, rhs4, rhs5e, rhs6e /)
      varnames = (/ "rhs_ij_1", "rhs_ij_2", "rhs_ij_3", "rhs_ij_4", "rhs_ij_5", "rhs_ij_6" /)
    end if
  end subroutine get_rhs
  
  subroutine get_amat(amat,varnames)
    implicit none
    
!    type(algexpr), dimension(n_amat), intent(out) :: amat
!    character(7),  dimension(n_amat), intent(out) :: varnames

    type(algexpr), allocatable, intent(out) :: amat(:)
    character(7),  allocatable, intent(out) :: varnames(:)
    
    if ( with_TiTe) then
      n_amat = 30
      allocate(amat(n_amat), varnames(n_amat))
      amat = (/ amat11,  amat12,  amat13,                            amat17,   &
                amat21,  amat22e, amat23, amat24e, amat25e, amat26e, amat27e,  &
                amat31e,          amat33,                                      &
                         amat42,          amat44,                              &
                amat51,  amat52e,                  amat55e,                    &
                amat61,  amat62e, amat63,          amat65e,  amat66e, amat67e, &
                amat71,  amat72e, amat73,          amat75e,  amat76e, amat77e  /)
      varnames = (/ "amat_11", "amat_12", "amat_13",                                  "amat_17", &
                    "amat_21", "amat_22", "amat_23", "amat_24", "amat_25", "amat_26", "amat_27", &
                    "amat_31",            "amat_33", &
                               "amat_42",            "amat_44", &
                    "amat_51", "amat_52",                       "amat_55", &
                    "amat_61", "amat_62", "amat_63",            "amat_65", "amat_66", "amat_67", &
                    "amat_71", "amat_72", "amat_73",            "amat_75", "amat_76", "amat_77"  /)
    else
      n_amat = 22
      allocate(amat(n_amat), varnames(n_amat))
      amat = (/ amat11, amat12,  amat13,                   amat16,  &
                amat21, amat22e, amat23, amat24e, amat25e, amat26e, &
                amat31e,           amat33,                          &
                        amat42,           amat44,                   &
                amat51, amat52e,                   amat55e,         &
                amat61, amat62e, amat63,           amat65e, amat66e /)
      varnames = (/ "amat_11", "amat_12", "amat_13",                       "amat_16", &
                    "amat_21", "amat_22", "amat_23", "amat_24", "amat_25", "amat_26", &
                    "amat_31",            "amat_33",                                  &
                               "amat_42",            "amat_44",                       &
                    "amat_51", "amat_52",                       "amat_55",            &
                    "amat_61", "amat_62", "amat_63",            "amat_65", "amat_66" /)
    end if 
 end subroutine get_amat
  
  subroutine get_aux(aux,varnames)
    implicit none
    type(algexpr), dimension(n_aux), intent(out) :: aux
    character(14), dimension(n_aux), intent(out) :: varnames
    integer      :: i
    character(2) :: num
    
    aux = (/ a_Bv2, ea_Bv2x, ea_Bv2y, ea_Bv2p, a_B2 /)
    varnames = (/ "eq(38,0,0,0,:)", "eq(38,1,0,0,:)", "eq(38,0,1,0,:)", "eq(38,0,0,1,:)", "eq(39,0,0,0,:)" /)
  end subroutine get_aux
  
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
  
  type(algexpr) function B0_parderiv(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    B0_parderiv = Bv_parderiv(a) + Bv_pbrack(a,Psi0)
  end function B0_parderiv
  
  type(algexpr) function B_parderiv(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    B_parderiv = Bv_pbrack(a,Psi)
  end function B_parderiv
  
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
  
  type(algexpr) function gradgrad_perp(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
  
    gradgrad_perp = gradprod(a,b) - B0_parderiv(a)*B0_parderiv(b)/B2
  end function gradgrad_perp
  
  type(algexpr) function gradDgrad_par(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b
    
    gradDgrad_par = (B_parderiv(a)*B0_parderiv(b) + B0_parderiv(a)*B_parderiv(b) - 2.d0*Bv2*inprod(Psi0,Psi)*B0_parderiv(a)*B0_parderiv(b)/B2)/B2
  end function gradDgrad_par
  
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
