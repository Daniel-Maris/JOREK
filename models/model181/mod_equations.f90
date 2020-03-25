module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
#ifdef DEBUG
    type(action), dimension(:), allocatable :: rhs1seq, rhs2seq, rhs3seq, rhs4seq
    type(action), dimension(:), allocatable :: amat11seq, amat12seq, amat14seq
    type(action), dimension(:), allocatable :: amat21seq, amat22seq, amat23seq, amat24seq
    type(action), dimension(:), allocatable :: amat32seq, amat33seq
    type(action), dimension(:), allocatable :: amat41seq, amat42seq, amat43seq, amat44seq
    type(action), dimension(:), allocatable :: aBv2seq, aBv2xseq, aBv2yseq, aBv2pseq
    type(action), dimension(:), allocatable :: aj0xseq, aj0yseq, aj0pseq, aj0chiseq, ajxseq, ajyseq, ajpseq, ajchiseq
#endif
    
    real*8, dimension(:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  ! Variables at the current time step
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=4)
  ! Changes since previous time step
  type(algexpr), parameter, private :: delta_Psi  = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: delta_Phi  = algexpr(basic=.true.,var=6)
  type(algexpr), parameter, private :: delta_rho  = algexpr(basic=.true.,var=7)
  type(algexpr), parameter, private :: delta_T    = algexpr(basic=.true.,var=8)
  ! Test function
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=9)
  ! Unknowns
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=10)
  ! Other quantities
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=11)
  type(algexpr), parameter, private :: psi_v      = algexpr(basic=.true.,var=12)
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=13)
  type(algexpr), parameter, private :: D_perp     = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: k_perp     = algexpr(basic=.true.,var=15)
  type(algexpr), parameter, private :: S_rho      = algexpr(basic=.true.,var=16)
  type(algexpr), parameter, private :: S_e        = algexpr(basic=.true.,var=17)
  type(algexpr), parameter, private :: S_j        = algexpr(basic=.true.,var=18)
  type(algexpr), parameter, private :: F          = algexpr(basic=.true.,var=19)
  type(algexpr), parameter, private :: eta        = algexpr(basic=.true.,var=20)
  type(algexpr), parameter, private :: deta_dT    = algexpr(basic=.true.,var=21)
  type(algexpr), parameter, private :: visco      = algexpr(basic=.true.,var=22)
  type(algexpr), parameter, private :: dvisco_dT  = algexpr(basic=.true.,var=23)
  type(algexpr), parameter, private :: k_par      = algexpr(basic=.true.,var=24)
  type(algexpr), parameter, private :: dk_par_dT  = algexpr(basic=.true.,var=25)
  ! Auxiliary variables (aux)
  type(algexpr), parameter, private :: Bv2        = algexpr(basic=.true.,var=26)
  type(algexpr), parameter, private :: j0x        = algexpr(basic=.true.,var=27)
  type(algexpr), parameter, private :: j0y        = algexpr(basic=.true.,var=28)
  type(algexpr), parameter, private :: j0p        = algexpr(basic=.true.,var=29)
  type(algexpr), parameter, private :: j0chi      = algexpr(basic=.true.,var=30)
  ! Auxiliary variables involving unknowns (aux2)
  type(algexpr), parameter, private :: jx         = algexpr(basic=.true.,var=31)
  type(algexpr), parameter, private :: jy         = algexpr(basic=.true.,var=32)
  type(algexpr), parameter, private :: jp         = algexpr(basic=.true.,var=33)
  type(algexpr), parameter, private :: jchi       = algexpr(basic=.true.,var=34)
  
  type(const), private :: tstep, zeta, theta, visco_num, eta_num, gamma, reta
  
  type(algexpr), private :: rhs1, rhs2, rhs3, rhs4
  type(algexpr), private :: amat11, amat12, amat14
  type(algexpr), private :: amat21, amat22, amat23, amat24
  type(algexpr), private :: amat32, amat33
  type(algexpr), private :: amat41, amat42, amat43, amat44
  type(algexpr), private :: a_Bv2, a_j0x, a_j0y, a_j0p, a_j0chi
  type(algexpr), private :: a_jx, a_jy, a_jp, a_jchi
  
  integer, parameter :: n_rhs = 4, n_amat = 13, n_aux = 8, n_aux2 = 4
  
  type(algexpr), private :: rhs2e, rhs3e, rhs4e
  type(algexpr), private :: amat22e, amat23e, amat24e
  type(algexpr), private :: amat32e, amat33e
  type(algexpr), private :: amat42e, amat43e, amat44e
  type(algexpr), private :: ea_Bv2x, ea_Bv2y, ea_Bv2p
  type(algexpr), private :: ea_j0chi, ea_jchi
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, Igamma => gamma, Itstep => tstep, Ivisco_num => visco_num, Ieta_num => eta_num, &
                           Ieta => eta, eta_ohmic
    implicit none
    
    tstep     = const(value = Itstep,          token = "tstep")
    zeta      = const(value = time_evol_zeta,  token = "time_evol_zeta")
    theta     = const(value = time_evol_theta, token = "time_evol_theta")
    visco_num = const(value = Ivisco_num,      token = "visco_num")
    eta_num   = const(value = Ieta_num,        token = "eta_num")
    gamma     = const(value = Igamma,          token = "gamma")
    if (Ieta .ne. 0.d0) then
      reta    = const(value = eta_ohmic/Ieta,  token = "reta")
    else
      reta    = const(value = 0.d0,            token = "reta")
    end if
    
    a_Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)
    
    a_j0x =  -Lap(Psi0)*dx(chi) + Bv_parderiv(dx(Psi0)) - gradprod(Psi0,dx(chi))
    a_j0y =  -Lap(Psi0)*dy(chi) + Bv_parderiv(dy(Psi0)) - gradprod(Psi0,dy(chi))
    a_j0p = (-Lap(Psi0)*dp(chi) + Bv_parderiv(dp(Psi0)) - gradprod(Psi0,dp(chi)))/R + 2.d0*(dx(Psi0)*dp(chi) - dp(Psi0)*dx(chi))/(R*R)
!    a_j0chi = -(Bv2*pLap(Psi0) + inprod(Bv2,Psi0))
    a_j0chi = -Bv2*Lap(Psi0) + Bv_parderiv(Bv_parderiv(Psi0)) - Bv_parderiv(Bv2)*Bv_parderiv(Psi0)/Bv2 - inprod(Bv2,Psi0)
    
    a_jx =  -Lap(Psi)*dx(chi) + Bv_parderiv(dx(Psi)) - gradprod(Psi,dx(chi))
    a_jy =  -Lap(Psi)*dy(chi) + Bv_parderiv(dy(Psi)) - gradprod(Psi,dy(chi))
    a_jp = (-Lap(Psi)*dp(chi) + Bv_parderiv(dp(Psi)) - gradprod(Psi,dp(chi)))/R + 2.d0*(dx(Psi)*dp(chi) - dp(Psi)*dx(chi))/(R*R)
!    a_jchi  = -(Bv2*pLap(Psi) + inprod(Bv2,Psi))
    a_jchi  = -Bv2*Lap(Psi) + Bv_parderiv(Bv_parderiv(Psi)) - Bv_parderiv(Bv2)*Bv_parderiv(Psi)/Bv2 - inprod(Bv2,Psi)
    
    rhs1 = -tstep*((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))*Bv_pbrack(v,psi_v)/Bv2 &
         + eta*(dx(v)*(dy(psi_v)*(j0p+S_j*dp(chi)/R) - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*(j0p+S_j*dp(chi)/R)) &
         + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)) + zeta*v*Bv_pbrack(psi_v,delta_Psi)
    
    rhs2 = tstep*(pLap(Phi0)*Bv_pbrack(v,Phi0)/Bv2 - Bv2*((j0x+dy(F)/R)*dx(v) + (j0y-dx(F)/R)*dy(v) + j0p*dp(v)/R)/rho0 &
         + j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0))/rho0 - D_perp*gradprod(rho0,inprod(v,Phi0)/rho0) + S_rho*inprod(v,Phi0)/rho0 &
         !- v*Bv_pbrack(rho0,T0)/rho0 
         + visco*Lap(v)*pLap(Phi0)) - zeta*inprod(v,delta_Phi)
    
    rhs3 = -tstep*(v*Bv_pbrack(rho0/Bv2,Phi0) + D_perp*gradprod(v,rho0) - S_rho*v) + zeta*v*delta_rho
    
    rhs4 = -tstep*(v*Bv_pbrack(rho0*T0,Phi0)/Bv2 - gamma*v*rho0*T0*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) + k_perp*gradprod(v,T0) &
         + k_par*Bv_parderiv(v)*(Bv_parderiv(T0) + Bv_pbrack(T0,Psi0))/Bv2 + D_perp*T0*gradprod(v,rho0) &
         - reta*eta*v*(j0x*j0x + j0y*j0y + j0p*j0p) - v*S_e) + zeta*v*(rho0*delta_T + T0*delta_rho)
    
    
    amat11 = (1.d0 + zeta)*v*Bv_pbrack(psi_v,Psi) &
           + tstep*theta*(Bv_pbrack(Psi,Phi0)*Bv_pbrack(v,psi_v)/Bv2 + eta*(dx(v)*(dy(psi_v)*jp - dp(psi_v)*jy/R) &
           + dy(v)*(dp(psi_v)*jx/R - dx(psi_v)*jp) + dp(v)*(dx(psi_v)*jy - dy(psi_v)*jx)/R))
    amat12 = tstep*theta*(Bv_pbrack(Psi0,Phi) - Bv_parderiv(Phi))*Bv_pbrack(v,psi_v)/Bv2
    amat14 = tstep*theta*deta_dT*T*(dx(v)*(dy(psi_v)*(j0p+S_j*dp(chi)/R) - dp(psi_v)*j0y/R) &
           + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*(j0p+S_j*dp(chi)/R)) + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)
    
    amat21 = tstep*theta*(Bv2*(jx*dx(v) + jy*dy(v) + jp*dp(v)/R) - jchi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0)) - j0chi*Bv_pbrack(v,Psi))/rho0
    amat22 = -(1.d0 + zeta)*inprod(v,Phi) - tstep*theta*(pLap(Phi)*Bv_pbrack(v,Phi0)/Bv2 + pLap(Phi0)*Bv_pbrack(v,Phi)/Bv2 &
           - D_perp*gradprod(inprod(v,Phi)/rho0,rho0) + S_rho*inprod(v,Phi)/rho0 + visco*Lap(v)*pLap(Phi))
    amat23 = -tstep*theta*(Bv2*((j0x+dy(F)/R)*dx(v) + (j0y-dx(F)/R)*dy(v) + j0p*dp(v)/R) - j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0)) &
           - S_rho*inprod(v,Phi0) &! + v*Bv_pbrack(rho0,T0)
           )*rho/(rho0*rho0) + tstep*theta*(D_perp*gradprod(inprod(v,Phi0)/rho0,rho) &
           - D_perp*gradprod(inprod(v,Phi0)*rho/(rho0*rho0),rho0) )! + v*Bv_pbrack(rho,T0)/rho0)
    amat24 = tstep*theta*v*Bv_pbrack(rho0,T)/rho0 - tstep*theta*dvisco_dT*T*pLap(Phi0)*Lap(v)
    
    amat32 = tstep*theta*v*Bv_pbrack(rho0/Bv2,Phi)
    amat33 = (1.d0 + zeta)*v*rho + tstep*theta*(v*Bv_pbrack(rho/Bv2,Phi0) + D_perp*gradprod(v,rho))
    
    amat41 = tstep*theta*k_par*Bv_parderiv(v)*Bv_pbrack(T0,Psi)/Bv2 - 2.d0*tstep*theta*reta*eta*v*(j0x*jx + j0y*jy + j0p*jp)
    amat42 = tstep*theta*v*(Bv_pbrack(rho0*T0,Phi) - gamma*rho0*T0*Bv_pbrack(Bv2,Phi)/(Bv2*Bv2))
    amat43 = (1.d0 + zeta)*v*T0*rho + tstep*theta*(v*Bv_pbrack(rho*T0,Phi0)/Bv2 - gamma*v*rho*T0*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) &
           + D_perp*T0*gradprod(v,rho))
    amat44 = (1.d0 + zeta)*v*rho0*T + tstep*theta*(v*Bv_pbrack(rho0*T,Phi0)/Bv2 - gamma*v*rho0*T*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2) + k_perp*gradprod(v,T) &
           + k_par*Bv_parderiv(v)*(Bv_parderiv(T) + Bv_pbrack(T,Psi0))/Bv2 + dk_par_dT*T*Bv_parderiv(v)*(Bv_parderiv(T0) + Bv_pbrack(T0,Psi0))/Bv2 &
           + D_perp*T*gradprod(v,rho0) + v*reta*deta_dT*T*(j0x*j0x + j0y*j0y + j0p*j0p))
    
    rhs2e = Dexpand(deepcopy(rhs2))
    rhs3e = Dexpand(deepcopy(rhs3))
    rhs4e = Dexpand(deepcopy(rhs4))
    
    amat22e = Dexpand(deepcopy(amat22)); amat23e = Dexpand(deepcopy(amat23)); amat24e = Dexpand(deepcopy(amat24))
    amat32e = Dexpand(deepcopy(amat32)); amat33e = Dexpand(deepcopy(amat33))
    amat42e = Dexpand(deepcopy(amat42)); amat43e = Dexpand(deepcopy(amat43)); amat44e = Dexpand(deepcopy(amat44))
    
    ea_Bv2x = Dexpand(deepcopy(dx(a_Bv2))); ea_Bv2y = Dexpand(deepcopy(dy(a_Bv2))); ea_Bv2p = Dexpand(deepcopy(dp(a_Bv2)))
    ea_j0chi = Dexpand(deepcopy(a_j0chi)); ea_jchi = Dexpand(deepcopy(a_jchi))
  end subroutine init_equations
  
  subroutine init_eq_struct()
    use data_structure, only: nbthreads
    implicit none
    integer :: i
    
    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        allocate(thread_eq(i)%eq(2*n_var+26,0:n_order-1,0:n_order-1,0:n_order-1))
#ifdef DEBUG
        allocate(thread_eq(i)%rhs1seq(countsubexprs(rhs1)))
        allocate(thread_eq(i)%rhs2seq(countsubexprs(rhs2e)))
        allocate(thread_eq(i)%rhs3seq(countsubexprs(rhs3e)))
        allocate(thread_eq(i)%rhs4seq(countsubexprs(rhs4e)))
        allocate(thread_eq(i)%amat11seq(countsubexprs(amat11)))
        allocate(thread_eq(i)%amat12seq(countsubexprs(amat12)))
        allocate(thread_eq(i)%amat14seq(countsubexprs(amat14)))
        allocate(thread_eq(i)%amat21seq(countsubexprs(amat21)))
        allocate(thread_eq(i)%amat22seq(countsubexprs(amat22e)))
        allocate(thread_eq(i)%amat23seq(countsubexprs(amat23e)))
        allocate(thread_eq(i)%amat24seq(countsubexprs(amat24e)))
        allocate(thread_eq(i)%amat32seq(countsubexprs(amat32e)))
        allocate(thread_eq(i)%amat33seq(countsubexprs(amat33e)))
        allocate(thread_eq(i)%amat41seq(countsubexprs(amat41)))
        allocate(thread_eq(i)%amat42seq(countsubexprs(amat42e)))
        allocate(thread_eq(i)%amat43seq(countsubexprs(amat43e)))
        allocate(thread_eq(i)%amat44seq(countsubexprs(amat44e)))
        allocate(thread_eq(i)%aBv2seq(countsubexprs(a_Bv2)))
        allocate(thread_eq(i)%aBv2xseq(countsubexprs(ea_Bv2x)))
        allocate(thread_eq(i)%aBv2yseq(countsubexprs(ea_Bv2y)))
        allocate(thread_eq(i)%aBv2pseq(countsubexprs(ea_Bv2p)))
        allocate(thread_eq(i)%aj0xseq(countsubexprs(a_j0x)))
        allocate(thread_eq(i)%aj0yseq(countsubexprs(a_j0y)))
        allocate(thread_eq(i)%aj0pseq(countsubexprs(a_j0p)))
        allocate(thread_eq(i)%aj0chiseq(countsubexprs(ea_j0chi)))
        allocate(thread_eq(i)%ajxseq(countsubexprs(a_jx)))
        allocate(thread_eq(i)%ajyseq(countsubexprs(a_jy)))
        allocate(thread_eq(i)%ajpseq(countsubexprs(a_jp)))
        allocate(thread_eq(i)%ajchiseq(countsubexprs(ea_jchi)))
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
      call buildsequence(rhs3e, thread_eq(i)%rhs3seq, thread_eq(i)%eq)
      call buildsequence(rhs4e, thread_eq(i)%rhs4seq, thread_eq(i)%eq)
      
      call buildsequence(amat11, thread_eq(i)%amat11seq, thread_eq(i)%eq)
      call buildsequence(amat12, thread_eq(i)%amat12seq, thread_eq(i)%eq)
      call buildsequence(amat14, thread_eq(i)%amat14seq, thread_eq(i)%eq)

      call buildsequence(amat21, thread_eq(i)%amat21seq, thread_eq(i)%eq)
      call buildsequence(amat22e, thread_eq(i)%amat22seq, thread_eq(i)%eq)
      call buildsequence(amat23e, thread_eq(i)%amat23seq, thread_eq(i)%eq)
      call buildsequence(amat24e, thread_eq(i)%amat24seq, thread_eq(i)%eq)
      
      call buildsequence(amat32e, thread_eq(i)%amat32seq, thread_eq(i)%eq)
      call buildsequence(amat33e, thread_eq(i)%amat33seq, thread_eq(i)%eq)
      
      call buildsequence(amat41, thread_eq(i)%amat41seq, thread_eq(i)%eq)
      call buildsequence(amat42e, thread_eq(i)%amat42seq, thread_eq(i)%eq)
      call buildsequence(amat43e, thread_eq(i)%amat43seq, thread_eq(i)%eq)
      call buildsequence(amat44e, thread_eq(i)%amat44seq, thread_eq(i)%eq)
      
      call buildsequence(a_Bv2, thread_eq(i)%aBv2seq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2x, thread_eq(i)%aBv2xseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2y, thread_eq(i)%aBv2yseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2p, thread_eq(i)%aBv2pseq, thread_eq(i)%eq)
      
      call buildsequence(a_j0x, thread_eq(i)%aj0xseq, thread_eq(i)%eq)
      call buildsequence(a_j0y, thread_eq(i)%aj0yseq, thread_eq(i)%eq)
      call buildsequence(a_j0p, thread_eq(i)%aj0pseq, thread_eq(i)%eq)
      call buildsequence(ea_j0chi, thread_eq(i)%aj0chiseq, thread_eq(i)%eq)
      call buildsequence(a_jx, thread_eq(i)%ajxseq, thread_eq(i)%eq)
      call buildsequence(a_jy, thread_eq(i)%ajyseq, thread_eq(i)%eq)
      call buildsequence(a_jp, thread_eq(i)%ajpseq, thread_eq(i)%eq)
      call buildsequence(ea_jchi, thread_eq(i)%ajchiseq, thread_eq(i)%eq)
    end do
  end subroutine build_all_seq
#endif
  
  subroutine get_rhs(rhs,varnames)
    implicit none
    type(algexpr), dimension(n_rhs), intent(out) :: rhs
    character(8),  dimension(n_rhs), intent(out) :: varnames
    
    rhs = (/ rhs1, rhs2e, rhs3e, rhs4e /)
    varnames = (/ "rhs_ij_1", "rhs_ij_2", "rhs_ij_3", "rhs_ij_4" /)
  end subroutine get_rhs
  
  subroutine get_amat(amat,varnames)
    implicit none
    type(algexpr), dimension(n_amat), intent(out) :: amat
    character(7),  dimension(n_amat), intent(out) :: varnames
    
    amat = (/ amat11, amat12,           amat14,  &
              amat21, amat22e, amat23e, amat24e, &
                      amat32e, amat33e,          &
              amat41, amat42e, amat43e, amat44e   /)
    varnames = (/ "amat_11", "amat_12",            "amat_14", &
                  "amat_21", "amat_22", "amat_23", "amat_24", &
                             "amat_32", "amat_33",            &
                  "amat_41", "amat_42", "amat_43", "amat_44"   /)
  end subroutine get_amat
  
  subroutine get_aux(aux,varnames)
    implicit none
    type(algexpr), dimension(n_aux), intent(out) :: aux
    character(12), dimension(n_aux), intent(out) :: varnames
    integer      :: i
    character(2) :: num
    
    aux = (/ a_Bv2, ea_Bv2x, ea_Bv2y, ea_Bv2p, a_j0x, a_j0y, a_j0p, ea_j0chi /)
    varnames(1:4) = (/ "eq(26,0,0,0)", "eq(26,1,0,0)", "eq(26,0,1,0)", "eq(26,0,0,1)" /)
    do i=5,n_aux
      write(num,'(I2)') ((i-3)+2*n_var+17)
      varnames(i) = "eq(" // num // ",0,0,0)"
    end do
  end subroutine get_aux
  
  subroutine get_aux2(aux2,varnames)
    implicit none
    type(algexpr), dimension(n_aux2), intent(out) :: aux2
    character(12), dimension(n_aux2), intent(out) :: varnames
    integer      :: i
    character(2) :: num
    
    aux2 = (/ a_jx, a_jy, a_jp, ea_jchi /)
    do i=1,n_aux2
      write(num,'(I2)') (i+2*n_var+22)
      varnames(i) = "eq(" // num // ",0,0,0)"
    end do
  end subroutine get_aux2
  
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
    
    Lap = dx(dx(a)) + dx(a)/R + dy(dy(a)) + dp(dp(a))/(R*R)
  end function Lap
  
  type(algexpr) function pLap(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    pLap = Lap(a) - Bv_parderiv(Bv_parderiv(a)/Bv2)
  end function pLap
end module mod_equations
