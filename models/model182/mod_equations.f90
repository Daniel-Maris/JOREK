module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq
#ifdef DEBUG
    type(action), dimension(:), allocatable :: rhs1seq, rhs2seq, rhs3seq, rhs5seq
    type(action), dimension(:), allocatable :: amat11seq, amat12seq, amat13seq, amat14seq
    type(action), dimension(:), allocatable :: amat21seq, amat22seq, amat23seq, amat24seq, amat25seq
    type(action), dimension(:), allocatable :: amat31seq, amat33seq
    type(action), dimension(:), allocatable :: amat42seq, amat44seq
    type(action), dimension(:), allocatable :: amat52seq, amat55seq
    type(action), dimension(:), allocatable :: aBv2seq, aBv2xseq, aBv2yseq, aBv2pseq
    type(action), dimension(:), allocatable :: aj0xseq, aj0yseq, aj0pseq, aj0chiseq, atjxseq, atjyseq, atjpseq, atjchiseq
#endif
    
    real*8, dimension(:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  ! Variables at current time step
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: zj0        = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: w0         = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=6)
  ! Changes since previous time step
  type(algexpr), parameter, private :: delta_Psi  = algexpr(basic=.true.,var=7)
  type(algexpr), parameter, private :: delta_Phi  = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: delta_zj   = algexpr(basic=.true.,var=9)
  type(algexpr), parameter, private :: delta_w    = algexpr(basic=.true.,var=10)
  type(algexpr), parameter, private :: delta_rho  = algexpr(basic=.true.,var=11)
  type(algexpr), parameter, private :: delta_T    = algexpr(basic=.true.,var=12)
  ! Test function
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=13)
  ! Unknowns
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: zj         = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: w          = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=14)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=14)
  ! Other quantities
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=15)
  type(algexpr), parameter, private :: psi_v      = algexpr(basic=.true.,var=16)
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=17)
  type(algexpr), parameter, private :: D_perp     = algexpr(basic=.true.,var=18)
  type(algexpr), parameter, private :: k_perp     = algexpr(basic=.true.,var=19)
  type(algexpr), parameter, private :: S_rho      = algexpr(basic=.true.,var=20)
  type(algexpr), parameter, private :: S_e        = algexpr(basic=.true.,var=21)
  type(algexpr), parameter, private :: S_j        = algexpr(basic=.true.,var=22)
  type(algexpr), parameter, private :: F          = algexpr(basic=.true.,var=23)
  type(algexpr), parameter, private :: eta        = algexpr(basic=.true.,var=24)
  type(algexpr), parameter, private :: deta_dT    = algexpr(basic=.true.,var=25)
  type(algexpr), parameter, private :: visco      = algexpr(basic=.true.,var=26)
  type(algexpr), parameter, private :: dvisco_dT  = algexpr(basic=.true.,var=27)
  ! Auxiliary variables (aux)
  type(algexpr), parameter, private :: Bv2        = algexpr(basic=.true.,var=28)
  type(algexpr), parameter, private :: j0x        = algexpr(basic=.true.,var=29)
  type(algexpr), parameter, private :: j0y        = algexpr(basic=.true.,var=30)
  type(algexpr), parameter, private :: j0p        = algexpr(basic=.true.,var=31)
  type(algexpr), parameter, private :: j0chi      = algexpr(basic=.true.,var=32)
  ! Auxiliary variables involving unknowns (aux2)
  type(algexpr), parameter, private :: tjx        = algexpr(basic=.true.,var=33)
  type(algexpr), parameter, private :: tjy        = algexpr(basic=.true.,var=34)
  type(algexpr), parameter, private :: tjp        = algexpr(basic=.true.,var=35)
  type(algexpr), parameter, private :: tjchi      = algexpr(basic=.true.,var=36)
  
  type(const), private :: tstep, zeta, theta, visco_num, gamma, k_par
  
  type(algexpr), private :: rhs1, rhs2, rhs3, rhs4, rhs5
  type(algexpr), private :: amat11, amat12, amat13, amat14
  type(algexpr), private :: amat21, amat22, amat23, amat24, amat25
  type(algexpr), private :: amat31, amat33
  type(algexpr), private :: amat42, amat44
  type(algexpr), private :: amat52, amat55
  type(algexpr), private :: a_Bv2, a_j0x, a_j0y, a_j0p, a_j0chi
  type(algexpr), private :: a_tjx, a_tjy, a_tjp, a_tjchi
  
  integer, parameter :: n_rhs = 3, n_amat = 14, n_aux = 8, n_aux2 = 4
  
  type(algexpr), private :: rhs1e, rhs2e, rhs3e, rhs4e, rhs5e
  type(algexpr), private :: amat11e, amat12e, amat13e
  type(algexpr), private :: amat21e, amat22e, amat24e, amat25e
  type(algexpr), private :: amat52e, amat55e
  type(algexpr), private :: ea_Bv2x, ea_Bv2y, ea_Bv2p, ea_j0chi
  type(algexpr), private :: ea_tjchi
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, zk_par, Igamma => gamma, Itstep => tstep, Ivisco_num => visco_num
    implicit none
    
    tstep     = const(value = Itstep,          token = "tstep")
    zeta      = const(value = time_evol_zeta,  token = "time_evol_zeta")
    theta     = const(value = time_evol_theta, token = "time_evol_theta")
    visco_num = const(value = Ivisco_num,      token = "visco_num")
    gamma     = const(value = Igamma,          token = "gamma")
    k_par     = const(value = zk_par,          token = "zk_par")
    
    a_Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)
    
    a_j0x =  -(zj0 + 2.d0*dx(Psi0)/R)*dx(chi) + Bv_parderiv(dx(Psi0)) - gradprod(Psi0,dx(chi))
    a_j0y =  -(zj0 + 2.d0*dx(Psi0)/R)*dy(chi) + Bv_parderiv(dy(Psi0)) - gradprod(Psi0,dy(chi))
    a_j0p = (-(zj0 + 2.d0*dx(Psi0)/R)*dp(chi) + Bv_parderiv(dp(Psi0)) - gradprod(Psi0,dp(chi)))/R + 2.d0*(dx(Psi0)*dp(chi) - dp(Psi0)*dx(chi))/(R*R)
    a_j0chi = -(Bv2*zj0 + 2.d0*Bv2*dx(Psi0)/R - Bv_parderiv(Bv_parderiv(Psi0)) + Bv_parderiv(Bv2)*Bv_parderiv(Psi0)/Bv2 + inprod(Bv2,Psi0))
    
    a_tjx =  -2.d0*dx(Psi)*dx(chi)/R + Bv_parderiv(dx(Psi)) - gradprod(Psi,dx(chi))
    a_tjy =  -2.d0*dx(Psi)*dy(chi)/R + Bv_parderiv(dy(Psi)) - gradprod(Psi,dy(chi))
    a_tjp = (-2.d0*dx(Psi)*dp(chi)/R + Bv_parderiv(dp(Psi)) - gradprod(Psi,dp(chi)))/R + 2.d0*(dx(Psi)*dp(chi) - dp(Psi)*dx(chi))/(R*R)
    a_tjchi  = -(2.d0*Bv2*dx(Psi)/R - Bv_parderiv(Bv_parderiv(Psi)) + Bv_parderiv(Bv2)*Bv_parderiv(Psi)/Bv2 + inprod(Bv2,Psi))
    
    rhs1 = -tstep*((Bv_pbrack(Psi0,Phi0) - Bv_parderiv(Phi0))*Bv_pbrack(v,psi_v)/Bv2 &
         + eta*(dx(v)*(dy(psi_v)*(j0p+S_j*dp(chi)/R) - dp(psi_v)*j0y/R) + dy(v)*(dp(psi_v)*j0x/R - dx(psi_v)*(j0p+S_j*dp(chi)/R)) &
         + dp(v)*(dx(psi_v)*j0y - dy(psi_v)*j0x)/R)) &
         + zeta*v*Bv_pbrack(psi_v,delta_Psi)
    
    rhs2 = tstep*(w0*Bv_pbrack(v,Phi0)/Bv2 - Bv2*((j0x+dy(F)/R)*dx(v) + (j0y-dx(F)/R)*dy(v) + j0p*dp(v)/R)/rho0 &
         + j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0))/rho0 - D_perp*gradprod(rho0,inprod(v,Phi0)/rho0) &
         + S_rho*inprod(v,Phi0)/rho0 &! - v*Bv_pbrack(rho0,T0)/rho0 
         - visco*gradprod(v,w0) - visco_num*Lap(v)*Lap(w0)/R) - zeta*inprod(v,delta_Phi)
    
    rhs3 = v*zj0 + gradprod(v,Psi0)
    
    rhs4 = v*w0 + inprod(v,Phi0)
    
    rhs5 = -tstep*(v*Bv_pbrack(rho0/Bv2,Phi0) + D_perp*gradprod(v,rho0) - S_rho*v) + zeta*v*delta_rho
    
    
    amat11 = (1.d0 + zeta)*v*Bv_pbrack(psi_v,Psi) + tstep*theta*(Bv_pbrack(Psi,Phi0)*Bv_pbrack(v,psi_v)/Bv2 + &
             eta*(dx(v)*(dy(psi_v)*tjp - dp(psi_v)*tjy/R) &
           + dy(v)*(dp(psi_v)*tjx/R - dx(psi_v)*tjp) + dp(v)*(dx(psi_v)*tjy - dy(psi_v)*tjx)/R))
    amat12 = tstep*theta*(Bv_pbrack(Psi0,Phi) - Bv_parderiv(Phi))*Bv_pbrack(v,psi_v)/Bv2
    amat13 = (-tstep)*theta*eta*zj*Bv_pbrack(v,psi_v)
    
    amat21 = tstep*theta*(Bv2*(tjx*dx(v) + tjy*dy(v) + tjp*dp(v)/R) - tjchi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0)) - j0chi*Bv_pbrack(v,Psi))/rho0
    amat22 = -(1.d0 + zeta)*inprod(v,Phi) - tstep*theta*(w0*Bv_pbrack(v,Phi)/Bv2 - D_perp*gradprod(inprod(v,Phi)/rho0,rho0) + S_rho*inprod(v,Phi)/rho0)
    amat23 = tstep*theta*Bv2*zj*Bv_pbrack(v,Psi0)/rho0
    amat24 = -tstep*theta*w*Bv_pbrack(v,Phi0)/Bv2 + tstep*theta*(visco*gradprod(v,w) + visco_num*Lap(v)*Lap(w)/R)
    amat25 = -tstep*theta*(Bv2*((j0x+dy(F)/R)*dx(v) + (j0y-dx(F)/R)*dy(v) + j0p*dp(v)/R) - j0chi*(Bv_parderiv(v) + Bv_pbrack(v,Psi0)) &
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
    
    ea_Bv2x = Dexpand(deepcopy(dx(a_Bv2))); ea_Bv2y = Dexpand(deepcopy(dy(a_Bv2))); ea_Bv2p = Dexpand(deepcopy(dp(a_Bv2)))
    ea_j0chi = Dexpand(deepcopy(a_j0chi)); ea_tjchi = Dexpand(deepcopy(a_tjchi))
  end subroutine init_equations
  
  subroutine init_eq_struct()
    use data_structure, only: nbthreads
    implicit none
    integer :: i
    
    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        allocate(thread_eq(i)%eq(2*n_var+24,0:n_order-1,0:n_order-1,0:n_order-1))
#ifdef DEBUG
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
        allocate(thread_eq(i)%aBv2seq(countsubexprs(a_Bv2)))
        allocate(thread_eq(i)%aBv2xseq(countsubexprs(ea_Bv2x)))
        allocate(thread_eq(i)%aBv2yseq(countsubexprs(ea_Bv2y)))
        allocate(thread_eq(i)%aBv2pseq(countsubexprs(ea_Bv2p)))
        allocate(thread_eq(i)%aj0xseq(countsubexprs(a_j0x)))
        allocate(thread_eq(i)%aj0yseq(countsubexprs(a_j0y)))
        allocate(thread_eq(i)%aj0pseq(countsubexprs(a_j0p)))
        allocate(thread_eq(i)%aj0chiseq(countsubexprs(ea_j0chi)))
        allocate(thread_eq(i)%atjxseq(countsubexprs(a_tjx)))
        allocate(thread_eq(i)%atjyseq(countsubexprs(a_tjy)))
        allocate(thread_eq(i)%atjpseq(countsubexprs(a_tjp)))
        allocate(thread_eq(i)%atjchiseq(countsubexprs(ea_tjchi)))
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
      
      call buildsequence(a_Bv2, thread_eq(i)%aBv2seq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2x, thread_eq(i)%aBv2xseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2y, thread_eq(i)%aBv2yseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2p, thread_eq(i)%aBv2pseq, thread_eq(i)%eq)
      
      call buildsequence(a_j0x, thread_eq(i)%aj0xseq, thread_eq(i)%eq)
      call buildsequence(a_j0y, thread_eq(i)%aj0yseq, thread_eq(i)%eq)
      call buildsequence(a_j0p, thread_eq(i)%aj0pseq, thread_eq(i)%eq)
      call buildsequence(ea_j0chi, thread_eq(i)%aj0chiseq, thread_eq(i)%eq)
      call buildsequence(a_tjx, thread_eq(i)%atjxseq, thread_eq(i)%eq)
      call buildsequence(a_tjy, thread_eq(i)%atjyseq, thread_eq(i)%eq)
      call buildsequence(a_tjp, thread_eq(i)%atjpseq, thread_eq(i)%eq)
      call buildsequence(ea_tjchi, thread_eq(i)%atjchiseq, thread_eq(i)%eq)
    end do
  end subroutine build_all_seq
#endif
  
  subroutine get_rhs(rhs,varnames)
    implicit none
    type(algexpr), dimension(n_rhs), intent(out) :: rhs
    character(8),  dimension(n_rhs), intent(out) :: varnames
    
    rhs = (/ rhs1e, rhs2e, rhs5e /)
    varnames = (/ "rhs_ij_1", "rhs_ij_2", "rhs_ij_5" /)
  end subroutine get_rhs
  
  subroutine get_amat(amat,varnames)
    implicit none
    type(algexpr), dimension(n_amat), intent(out) :: amat
    character(7),  dimension(n_amat), intent(out) :: varnames
    
    amat = (/ amat11e, amat12e, amat13,                   &
              amat21e, amat22e, amat23, amat24e, amat25e, &
              amat31,           amat33, &
                      amat42,           amat44, &
                      amat52e,                   amat55e  /)
    varnames = (/ "amat_11", "amat_12", "amat_13",                       &
                  "amat_21", "amat_22", "amat_23", "amat_24", "amat_25", &
                  "amat_31",            "amat_33", &
                             "amat_42",            "amat_44", &
                             "amat_52",                       "amat_55"  /)
  end subroutine get_amat
  
  subroutine get_aux(aux,varnames)
    implicit none
    type(algexpr), dimension(n_aux), intent(out) :: aux
    character(12), dimension(n_aux), intent(out) :: varnames
    integer      :: i
    character(2) :: num
    
    aux = (/ a_Bv2, ea_Bv2x, ea_Bv2y, ea_Bv2p, a_j0x, a_j0y, a_j0p, ea_j0chi /)
    varnames(1:4) = (/ "eq(28,0,0,0)", "eq(28,1,0,0)", "eq(28,0,1,0)", "eq(28,0,0,1)" /)
    do i=5,n_aux
      write(num,'(I2)') ((i-3)+2*n_var+15)
      varnames(i) = "eq(" // num // ",0,0,0)"
    end do
  end subroutine get_aux
  
  subroutine get_aux2(aux2,varnames)
    implicit none
    type(algexpr), dimension(n_aux2), intent(out) :: aux2
    character(12), dimension(n_aux2), intent(out) :: varnames
    integer      :: i
    character(2) :: num
    
    aux2 = (/ a_tjx, a_tjy, a_tjp, ea_tjchi /)
    do i=1,n_aux2
      write(num,'(I2)') (i+2*n_var+20)
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
    
    Lap = dx(R*dx(a))/R + dy(dy(a)) + dp(dp(a))/(R*R)
  end function Lap
  
  type(algexpr) function pLap(a)
    implicit none
    type(algexpr), intent(in) :: a
    
    pLap = Lap(a) - Bv_parderiv(Bv_parderiv(a)/Bv2)
  end function pLap
end module mod_equations
