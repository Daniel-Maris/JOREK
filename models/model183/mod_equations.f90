module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none
  
  type type_thread_eq

    ! Indices in eq array: variable index (see algexpr's below), R derivative order, z derivative order, phi derivative order,
    !   separation of terms with covariant phi derivatives in test function and unknown (FFT)
    real*8, dimension(:,:,:,:,:), allocatable :: eq
  end type type_thread_eq
  
  ! Indices
  integer, parameter  :: var_dPsi        = n_var+var_Psi
  integer, parameter  :: var_dPhi        = n_var+var_Phi
  integer, parameter  :: var_dzj         = n_var+var_zj
  integer, parameter  :: var_dw          = n_var+var_w
  integer, parameter  :: var_drho        = n_var+var_rho
  integer, parameter  :: var_dT          = n_var+var_T
  integer, parameter  :: var_dvpar       = n_var+var_vpar
  integer, parameter  :: var_dTi         = n_var+var_Ti
  integer, parameter  :: var_dTe         = n_var+var_Te
  integer, parameter  :: var_v           = 2*n_var+1
  integer, parameter  :: var_varStar     = 2*n_var+2
  integer, parameter  :: var_chi         = 2*n_var+3
  integer, parameter  :: var_R           = 2*n_var+4 
  integer, parameter  :: var_D_perp      = 2*n_var+5 
  integer, parameter  :: var_S_rho       = 2*n_var+6 
  integer, parameter  :: var_D_par       = 2*n_var+7 
  integer, parameter  :: var_S_j         = 2*n_var+8 
  integer, parameter  :: var_eta         = 2*n_var+9 
  integer, parameter  :: var_deta_dT     = 2*n_var+10
  integer, parameter  :: var_visco       = 2*n_var+11
  integer, parameter  :: var_dvisco_dT   = 2*n_var+12
  integer, parameter  :: var_k_perp      = 2*n_var+13
  integer, parameter  :: var_k_perp_i    = 2*n_var+14
  integer, parameter  :: var_k_perp_e    = 2*n_var+15
  integer, parameter  :: var_S_e         = 2*n_var+16
  integer, parameter  :: var_S_e_i       = 2*n_var+17
  integer, parameter  :: var_S_e_e       = 2*n_var+18
  integer, parameter  :: var_k_par       = 2*n_var+19
  integer, parameter  :: var_k_par_i     = 2*n_var+20
  integer, parameter  :: var_k_par_e     = 2*n_var+21
  integer, parameter  :: var_dk_par_dT   = 2*n_var+22
  integer, parameter  :: var_dk_par_dT_i = 2*n_var+23
  integer, parameter  :: var_dk_par_dT_e = 2*n_var+24
  integer, parameter  :: var_dTe_i       = 2*n_var+25
  integer, parameter  :: var_ddTe_i_dT_i = 2*n_var+26
  integer, parameter  :: var_ddTe_i_dT_e = 2*n_var+27
  integer, parameter  :: var_ddTe_i_drho = 2*n_var+28
  integer, parameter  :: var_Bv2         = 2*n_var+29
  integer, parameter  :: var_B2          = 2*n_var+30

  ! Variables at current time step
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=var_Psi)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=var_Phi)
  type(algexpr), parameter, private :: zj0        = algexpr(basic=.true.,var=var_zj)
  type(algexpr), parameter, private :: w0         = algexpr(basic=.true.,var=var_w)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=var_rho)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=var_T)
  type(algexpr), parameter, private :: vpar0      = algexpr(basic=.true.,var=var_vpar)
  type(algexpr), parameter, private :: T0_i       = algexpr(basic=.true.,var=var_Ti)
  type(algexpr), parameter, private :: T0_e       = algexpr(basic=.true.,var=var_Te)
  ! Changes since previous time step
  type(algexpr), parameter, private :: delta_Psi  = algexpr(basic=.true.,var=var_dPsi )
  type(algexpr), parameter, private :: delta_Phi  = algexpr(basic=.true.,var=var_dPhi )
  type(algexpr), parameter, private :: delta_zj   = algexpr(basic=.true.,var=var_dzj  )
  type(algexpr), parameter, private :: delta_w    = algexpr(basic=.true.,var=var_dw   )
  type(algexpr), parameter, private :: delta_rho  = algexpr(basic=.true.,var=var_drho )
  type(algexpr), parameter, private :: delta_T    = algexpr(basic=.true.,var=var_dT   )
  type(algexpr), parameter, private :: delta_vpar = algexpr(basic=.true.,var=var_dvpar)
  type(algexpr), parameter, private :: delta_T_i  = algexpr(basic=.true.,var=var_dTi  )
  type(algexpr), parameter, private :: delta_T_e  = algexpr(basic=.true.,var=var_dTe  )  
  ! Test function
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=var_v)
  ! Unknowns
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: zj         = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: w          = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: vpar       = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: T_i        = algexpr(basic=.true.,var=var_varStar)
  type(algexpr), parameter, private :: T_e        = algexpr(basic=.true.,var=var_varStar)
  ! Other quantities
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=var_chi        )
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=var_R          )
  type(algexpr), parameter, private :: D_perp     = algexpr(basic=.true.,var=var_D_perp     )
  type(algexpr), parameter, private :: S_rho      = algexpr(basic=.true.,var=var_S_rho      )
  type(algexpr), parameter, private :: D_par      = algexpr(basic=.true.,var=var_D_par      )
  type(algexpr), parameter, private :: S_j        = algexpr(basic=.true.,var=var_S_j        )
  type(algexpr), parameter, private :: eta        = algexpr(basic=.true.,var=var_eta        )
  type(algexpr), parameter, private :: deta_dT    = algexpr(basic=.true.,var=var_deta_dT    )
  type(algexpr), parameter, private :: visco      = algexpr(basic=.true.,var=var_visco      )
  type(algexpr), parameter, private :: dvisco_dT  = algexpr(basic=.true.,var=var_dvisco_dT  )
  type(algexpr), parameter, private :: k_perp     = algexpr(basic=.true.,var=var_k_perp     )
  type(algexpr), parameter, private :: k_perp_i   = algexpr(basic=.true.,var=var_k_perp_i   )
  type(algexpr), parameter, private :: k_perp_e   = algexpr(basic=.true.,var=var_k_perp_e   )
  type(algexpr), parameter, private :: S_e        = algexpr(basic=.true.,var=var_S_e        )
  type(algexpr), parameter, private :: S_e_i      = algexpr(basic=.true.,var=var_S_e_i      )
  type(algexpr), parameter, private :: S_e_e      = algexpr(basic=.true.,var=var_S_e_e      )
  type(algexpr), parameter, private :: k_par      = algexpr(basic=.true.,var=var_k_par      ) 
  type(algexpr), parameter, private :: k_par_i    = algexpr(basic=.true.,var=var_k_par_i    )
  type(algexpr), parameter, private :: k_par_e    = algexpr(basic=.true.,var=var_k_par_e    )
  type(algexpr), parameter, private :: dk_par_dT  = algexpr(basic=.true.,var=var_dk_par_dT  )
  type(algexpr), parameter, private :: dk_par_dT_i= algexpr(basic=.true.,var=var_dk_par_dT_i)
  type(algexpr), parameter, private :: dk_par_dT_e= algexpr(basic=.true.,var=var_dk_par_dT_e)
  type(algexpr), parameter, private :: dTe_i      = algexpr(basic=.true.,var=var_dTe_i      )
  type(algexpr), parameter, private :: ddTe_i_dT_i= algexpr(basic=.true.,var=var_ddTe_i_dT_i)
  type(algexpr), parameter, private :: ddTe_i_dT_e= algexpr(basic=.true.,var=var_ddTe_i_dT_e)
  type(algexpr), parameter, private :: ddTe_i_drho= algexpr(basic=.true.,var=var_ddTe_i_drho)

  ! Auxiliary variables (aux)
  type(algexpr), parameter, private :: Bv2        = algexpr(basic=.true.,var=var_Bv2)
  type(algexpr), parameter, private :: B2         = algexpr(basic=.true.,var=var_B2)

  type(const), private :: tstep, zeta, theta, visco_num, eta_num, D_perp_num, k_perp_num, gamma, reta
  
  type(algexpr), public  :: rhs_semianalytic(n_var)
  type(algexpr), public  :: amat_semianalytic(n_var, n_var)
  type(algexpr), private :: a_Bv2, a_B2
  
  integer, parameter :: n_aux  = 5
  
  type(algexpr), private :: ea_Bv2x, ea_Bv2y, ea_Bv2p
  
  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq
  
  contains
  
  subroutine init_equations()
    use phys_module, only: time_evol_zeta, time_evol_theta, Igamma => gamma, Itstep => tstep, Ivisco_num => visco_num, Ieta_num => eta_num, &
                           ID_perp_num => D_perp_num, zk_perp_num, Ieta => eta, eta_ohmic
    implicit none
    
    integer  :: i, i_var, j_var

#ifdef WITH_TiTe
#define DIMT 2
#else
#define DIMT 1
#endif
    integer                           :: num_T, T_indices(DIMT)
    type(algexpr), dimension(DIMT, 7) :: T_expr
    type(algexpr)                     :: i_T0, i_delta_T, i_T, i_k_perp, i_k_par, i_dk_par_dT, i_S_e

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
    
    !###################################################################################################
    !#  Auxiliary vacuum and total magnetic field                                                      #
    !###################################################################################################
    a_Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)                               ! Magnitude of vacuum field squared
    a_B2 = Bv2 + Bv2*inprod(Psi0,Psi0)                                                              ! Magnitude of total field squared
    
    !###################################################################################################
    !#  Induction Equation                                                                             #
    !#     Questions:                                                                                  #
    !#        - Should the deta_dT term include S_j?                                                   #
    !###################################################################################################
    rhs_semianalytic(var_Psi) = tstep*v*((Bv_parderiv(Phi0) - Bv_pbrack(Psi0,Phi0))/Bv2     &       ! v x B ideal component
                              + eta*(zj0 - S_j))                                            &       ! Resistivity and current source
                              + tstep*eta_num*inprod(v,zj0)                                 &       ! Hyper resistivity
                              + zeta*v*delta_Psi                                                    ! dPsi_dt
                                                                                                    
    amat_semianalytic(var_Psi, var_Psi) = (1.d0 + zeta)*v*Psi                               &       ! dPsi_dt
                                        + tstep*theta*v*Bv_pbrack(Psi,Phi0)/Bv2                     ! v x B ideal component
    amat_semianalytic(var_Psi, var_Phi) = (-tstep*theta)*v*(Bv_parderiv(Phi)                &       ! v x B ideal component
                                        - Bv_pbrack(Psi0,Phi))/Bv2                                  ! v x B ideal component
    amat_semianalytic(var_Psi, var_zj ) = (-tstep*theta)*(eta*v*zj                          &       ! Resistivity and current source
                                        + eta_num*inprod(v,zj))                                     ! Hyper resistivity
    if (with_TiTe) then                                                                             
      amat_semianalytic(var_Psi, var_Te) = (-tstep*theta)*v*deta_dT*T_e*(zj0 - S_j)                 ! Resistivity and current source
    else                                                                                            
      amat_semianalytic(var_Psi, var_T) = (-tstep*theta)*v*deta_dT*T*(zj0 - S_j)                    ! Resistivity and current source
    end if

    !###################################################################################################
    !#  Perpendicular Momentum Equation                                                                #
    !#                                                                                                 #
    !#  Missing terms:                                                                                 #
    !#     - d(v_par B)_dt:       change in parallel momentum                                          #
    !###################################################################################################
    rhs_semianalytic(var_Phi) = -tstep*((Bv_pbrack(rho0/Bv2,v)*inprod(Phi0,Phi0)/2.d0   &            ! 1/2 rho grad(v^2)  
                              - Bv_pbrack(v,Phi0)*rho0*w0/Bv2                           &            ! rho omega x v 
                              - Bv_pbrack(rho0/Bv2,Phi0)*inprod(v,Phi0))/Bv2            &            ! density source
                              - v*Bv_parderiv(zj0)                                      &            ! j x B component
                              - v*Bv_pbrack(zj0,Psi0)                                   &            ! j x B component
                              + visco*inprod(v,w0)                                      &            ! Ad-hoc viscous tensor
                              + visco_num*Lap(v)*Lap(w0))                               &            ! Hyper viscosity
                              - zeta*(rho0*inprod(v,delta_Phi)                          &            ! d(rho v)_dt component
                              + delta_rho*inprod(v,Phi0))/Bv2                                        ! d(rho v)_dt component
                                                                                                     ! density source
    if (with_TiTe) then                                                                              
      rhs_semianalytic(var_Phi) = rhs_semianalytic(var_Phi)                             &            
                                -tstep * Bv_pbrack(v,rho0*(T0_i+T0_e))/Bv2                           ! grad(p) component
    else                                                                                             
      rhs_semianalytic(var_Phi) = rhs_semianalytic(var_Phi)                             &            
                                  -tstep * Bv_pbrack(v,rho0*T0)/Bv2                                  ! grad(p) component
    end if                                                                                           
                                                                                                     
    amat_semianalytic(var_Phi, var_Psi) = (-tstep*theta)*v*Bv_pbrack(zj0,Psi)                        ! j x B component
                                                                                                     
    amat_semianalytic(var_Phi, var_Phi) = -(1.d0 + zeta)*rho0*inprod(v,Phi)/Bv2               &      ! d(rho v)_dt component
                                        + tstep*theta*(Bv_pbrack(rho0/Bv2,v)*inprod(Phi0,Phi) &      ! 1/2 rho grad(v^2)
                                        - rho0*w0*Bv_pbrack(v,Phi)/Bv2                        &      ! rho omega x v
                                        - Bv_pbrack(rho0/Bv2,Phi)*inprod(v,Phi0)              &      ! density source
                                        - Bv_pbrack(rho0/Bv2,Phi0)*inprod(v,Phi))/Bv2                ! density source
                                                                                                     
    amat_semianalytic(var_Phi,  var_zj) = (-tstep*theta)*v*(Bv_parderiv(zj)                   &      ! j x B component
                                        + Bv_pbrack(zj,Psi0))                                        ! j x B component
                                                                                                     
    amat_semianalytic(var_Phi,   var_w) = -tstep*theta*rho0*w*Bv_pbrack(v,Phi0)/(Bv2*Bv2)     &      ! rho omega x v
                                        + tstep*theta*(visco*inprod(v,w)                      &      ! Ad-hoc viscous tensor
                                        + visco_num*Lap(v)*Lap(w))                                   ! Hyper viscosity
    
    amat_semianalytic(var_Phi, var_rho) = -(1.d0 + zeta)*rho*inprod(v,Phi0)/Bv2                    & ! d(rho v)_dt component
                                        + tstep*theta*(Bv_pbrack(rho/Bv2,v)*inprod(Phi0,Phi0)/2.d0 & ! 1/2 rho grad(v^2)
                                        - rho*w0*Bv_pbrack(v,Phi0)/Bv2                             & ! rho omega x v
                                        - Bv_pbrack(rho/Bv2,Phi0)*inprod(v,Phi0))/Bv2                ! density source
    if (with_TiTe) then 
      amat_semianalytic(var_Phi, var_rho) = amat_semianalytic(var_Phi, var_rho)                    &
                                          + tstep*theta*Bv_pbrack(v,rho*(T0_i+T0_e))/Bv2             ! grad(p) component

      amat_semianalytic(var_Phi, var_Ti) = tstep*theta*Bv_pbrack(v,rho0*T_i)/Bv2                     ! grad(p) component

      amat_semianalytic(var_Phi, var_Te) = tstep*theta*Bv_pbrack(v,rho0*T_e)/Bv2                   & ! grad(p) component
                                         + tstep*theta*dvisco_dT*T_e*inprod(v,w0)                    ! dvisco_dT_e
    else
      amat_semianalytic(var_Phi, var_rho) = amat_semianalytic(var_Phi, var_rho)                    &
                                          + tstep*theta*Bv_pbrack(v,rho*(T0))/Bv2                    ! grad(p) component

      amat_semianalytic(var_Phi,   var_T) = tstep*theta*Bv_pbrack(v,rho0*T)/Bv2                    & ! grad(p) component
                                          + tstep*theta*dvisco_dT*T*inprod(v,w0)                     ! dvisco_dT
    end if

    !###################################################################################################
    !#  Current Definition Equation                                                                    #
    !###################################################################################################
    rhs_semianalytic(var_zj) = -Bv2*inprod(v,Psi0)                                                  & ! Lap(Psi) 
                             - v*Bv2*zj0                                                              ! current zj
    
    amat_semianalytic(var_zj, var_Psi) = theta*Bv2*inprod(v,Psi)                                      ! change in Lap(Psi)
    amat_semianalytic(var_zj,  var_zj) = theta*v*Bv2*zj                                               ! change in zj

    !###################################################################################################
    !#  Vorticity Definition Equation                                                                  #
    !###################################################################################################
    rhs_semianalytic(var_w) = -inprod(v,Phi0)                                                       & ! Lap(Phi) 
                            - v*w0                                                                    ! current w
    
    amat_semianalytic(var_w, var_Phi) = theta*inprod(v,Phi)                                           ! change in Lap(Phi)
    amat_semianalytic(var_w,   var_w) = theta*v*w                                                     ! change in w

    !###################################################################################################
    !#  Density Equation                                                                               #
    !###################################################################################################
    rhs_semianalytic(var_rho) = -tstep*(v*Bv_pbrack(rho0/Bv2,Phi0)                               &    ! div(rho v) component
                              + D_perp*gradprod(v,rho0)                                          &    ! D_perp grad(rho)
                              + (D_par - D_perp)*B0_parderiv(v)*B0_parderiv(rho0)/B2             &    ! (D_par - D_perp) * grad_par(rho)
                              - S_rho*v) + zeta*v*delta_rho                                           ! Density source
    
    amat_semianalytic(var_rho, var_Psi) = tstep*theta*(D_par - D_perp)*gradDgrad_par(v,rho0)          ! (D_par - D_perp) * grad_par(rho)
    amat_semianalytic(var_rho, var_Phi) = tstep*theta*v*Bv_pbrack(rho0/Bv2,Phi)                       ! div(rho v) component
    amat_semianalytic(var_rho, var_rho) = (1.d0 + zeta)*v*rho                                    &    ! drho_dt 
                                        + tstep*theta*(v*Bv_pbrack(rho/Bv2,Phi0)                 &    ! div(rho v) component
                                        + D_perp*gradprod(v,rho)                                 &    ! D_perp grad(rho)
                                        + (D_par - D_perp)*B0_parderiv(v)*B0_parderiv(rho)/B2)        ! (D_par - D_perp) * grad_par(rho)

    !###################################################################################################
    !#  Pressure Equation                                                                              #
    !###################################################################################################
    ! index, T0, T, k_perp, k_par, S_e
#ifdef WITH_TiTe
      num_T = 2
      T_expr(1,:)   = (/ T0_i, delta_T_i, T_i, k_perp_i, k_par_i, dk_par_dT_i, S_e_i /); 
      T_expr(2,:)   = (/ T0_e, delta_T_e, T_e, k_perp_e, k_par_e, dk_par_dT_e, S_e_e /)
      T_indices(:) = (/ var_Ti, var_Te /)
#else
      num_T = 1
      T_expr(1,:)   = (/ T0, delta_T, T, k_perp, k_par, dk_par_dT, S_e /)
      T_indices(1) = var_T
#endif

    do i = 1, num_T
      i_var = T_indices(i)
      i_T0     = T_expr(i,1); i_delta_T    = T_expr(i,2); i_T = T_expr(i,3); 
      i_k_perp = T_expr(i,4); i_k_par      = T_expr(i,5); i_dk_par_dT = T_expr(i,6)
      i_S_e    = T_expr(i,7); 

      rhs_semianalytic(i_var) = -tstep*(v*Bv_pbrack(rho0*i_T0, Phi0)/Bv2                             & ! v.grad(p) component
                               - gamma*v*rho0*i_T0*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2)                     & ! gamma p div(v) component
                               + i_k_perp*gradprod(v,i_T0)                                           & ! K_perp grad(T)
                               + (i_k_par-i_k_perp)*B0_parderiv(v)*B0_parderiv(i_T0)/B2              & ! (K_par - K_perp) grad_par(T)
                               + k_perp_num*Lap(v)*Lap(i_T0)                                       & ! ad-hoc hyper-conductivity
                               + D_perp*i_T0*gradprod(v, rho0)                                       & ! D_perp T grad(rho)
                               + (D_par - D_perp)*i_T0*B0_parderiv(v)*B0_parderiv(rho0)/B2           & ! (D_par - D_perp) T grad_par(rho) 
                               - v*i_S_e)                                                            & ! heat source
                               + zeta*v*(rho0*i_delta_T + i_T0*delta_rho)                              ! dp_dt

      amat_semianalytic(i_var, var_Psi) = tstep*theta*((i_k_par - i_k_perp)*gradDgrad_par(v,i_T0)    & ! (K_par - K_perp) grad_par(T)
                                        + (D_par - D_perp)*i_T0*gradDgrad_par(v,rho0))                 ! (D_par - D_perp) T grad_par(rho)

      amat_semianalytic(i_var, var_Phi) = tstep*theta*v*(Bv_pbrack(rho0*i_T0,Phi)                    & ! v.grad(p) component
                                        - gamma*rho0*i_T0*Bv_pbrack(Bv2,Phi)/Bv2)/Bv2                  ! gamma p div(v) component

      amat_semianalytic(i_var, var_rho) = (1.d0 + zeta)*v*rho*i_T0                                   & ! dp_dt
                                         + tstep*theta*(v*Bv_pbrack(rho*i_T0,Phi0)/Bv2               & ! v.grad(p) component
                                         - gamma*v*rho*i_T0*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2)            & ! gamma p div(v) component
                                         + D_perp*i_T0*gradprod(v,rho)                               & ! D_perp T grad(rho)
                                         + (D_par - D_perp)*i_T0*B0_parderiv(v)*B0_parderiv(rho)/B2)   ! (D_par - D_perp) T grad_par(rho)

      amat_semianalytic(i_var, i_var)   = (1.d0 + zeta)*v*rho0*i_T                                   &  ! dp_dt
                                        + tstep*theta*(v*Bv_pbrack(rho0*i_T,Phi0)/Bv2                &  ! v.grad(p) component
                                        - gamma*v*rho0*i_T*Bv_pbrack(Bv2,Phi0)/(Bv2*Bv2)             &  ! gamma p div(v) component
                                        + i_k_perp*gradprod(v,i_T)                                   &  ! K_perp grad(T)
                                        + (i_k_par - i_k_perp)*B0_parderiv(v)*B0_parderiv(i_T)/B2    &  ! (K_par - K_perp) grad_par(T)
                                        + i_dk_par_dT*i_T*B0_parderiv(v)*B0_parderiv(i_T0)/B2        &  ! (K_par - K_perp) grad_par(T)
                                        + k_perp_num*Lap(v)*Lap(i_T)                               &  ! ad-hoc hyper-conductivity
                                        + D_perp*i_T*gradprod(v,rho0)                                &  ! D_perp T grad(rho)
                                        + (D_par - D_perp)*i_T*B0_parderiv(v)*B0_parderiv(rho0)/B2)     ! (D_par - D_perp) T grad_par(rho)
    enddo

    if (with_TiTe) then
      rhs_semianalytic(var_Ti) = rhs_semianalytic(var_Ti) + tstep*v*dTe_i                                     ! thermalization
    
      amat_semianalytic(var_Ti, var_rho) = amat_semianalytic(var_Ti, var_rho) - tstep*theta*v*ddTe_i_drho*rho ! thermalization
      amat_semianalytic(var_Ti, var_Ti)  = amat_semianalytic(var_Ti, var_Ti)  - tstep*theta*v*ddTe_i_dT_i*T_i ! thermalization
      amat_semianalytic(var_Ti, var_Te)  = - tstep*theta*v*ddTe_i_dT_e*T_e                                    ! thermalization
      
      rhs_semianalytic(var_Te) = rhs_semianalytic(var_Te) - tstep*v*dTe_i                           &         ! thermalization
                               + tstep * (gamma - 1.d0)*reta*eta*v*Bv2*zj0*zj0

      amat_semianalytic(var_Te,  var_zj) = -2.d0*tstep*theta*(gamma - 1.d0)*v*reta*eta*Bv2*zj0*zj             ! ohmic heating

      amat_semianalytic(var_Te, var_rho) = amat_semianalytic(var_Te, var_rho) + tstep*theta*v*ddTe_i_drho*rho ! thermalization
      amat_semianalytic(var_Te,  var_Ti) = tstep * theta * v * ddTe_i_dT_i * T_i                              ! thermalization
      amat_semianalytic(var_Te,  var_Te) = amat_semianalytic(var_Te,  var_Te)                       &         
                                         + tstep*theta*(v*ddTe_i_dT_e*T_e                           &         ! Thermalization
                                         - v*reta*deta_dT*T_e*Bv2*zj0*zj0)                                    ! ohmic heating
    else            
      rhs_semianalytic(var_T) = rhs_semianalytic(var_T) + tstep*(gamma - 1.d0)*reta*eta*v*Bv2*zj0*zj0         ! ohmic heating

      amat_semianalytic(var_T,  var_zj) = -2.d0*tstep*theta*(gamma - 1.d0)*v*reta*eta*Bv2*zj0*zj              ! ohmic heating
      
      amat_semianalytic(var_T,   var_T) = amat_semianalytic(var_T,   var_T)                         &
                                        - tstep*theta*v*reta*deta_dT*T*Bv2*zj0*zj0                            ! ohmic heating
    end if

    !###################################################################################################
    !#  Parallel Momentum Equation                                                                     #
    !###################################################################################################
    if (with_vpar) then
      amat_semianalytic(var_vpar, var_vpar) = v*vpar
    endif

    ! Expansion of differential operators
    do i_var = 1, n_var
      if ((associated(rhs_semianalytic(i_var)%operand1)) .and. (associated(rhs_semianalytic(i_var)%operand2))) then
        rhs_semianalytic(i_var) = Dexpand(deepcopy(rhs_semianalytic(i_var)))
      endif
      do j_var = 1, n_var  
        if ((associated(amat_semianalytic(i_var, j_var)%operand1)) .and. (associated(amat_semianalytic(i_var, j_var)%operand2))) then
          amat_semianalytic(i_var, j_var) = Dexpand(deepcopy(amat_semianalytic(i_var, j_var)))
        endif
      enddo
    enddo
    ea_Bv2x = Dexpand(deepcopy(dx(a_Bv2))); ea_Bv2y = Dexpand(deepcopy(dy(a_Bv2))); ea_Bv2p = Dexpand(deepcopy(dp(a_Bv2)))
  end subroutine init_equations
  
  subroutine init_eq_struct()
    use data_structure, only: nbthreads
    implicit none
    integer :: i
    
    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        allocate(thread_eq(i)%eq(2*n_var+30,0:n_order-1,0:n_order-1,0:n_order-1,4))
      end do
    end if
  end subroutine init_eq_struct
  
  subroutine get_varnames(varnames)
    implicit none
    character(8), dimension(n_var), intent(out) :: varnames
      
    varnames(var_Psi)    = " var_Psi"
    varnames(var_Phi)    = " var_Phi"
    varnames( var_zj)    = "  var_zj"
    varnames(  var_w)    = "   var_w"
    varnames(var_rho)    = " var_rho"
    if (with_vpar) then
      varnames(var_vpar) = "var_vpar"
    endif
    if (with_TiTe) then
      varnames( var_Ti)    = "  var_Ti"
      varnames( var_Te)    = "  var_Te"
    else
      varnames( var_T)     = "   var_T"
    endif
  end subroutine get_varnames

  subroutine get_aux(aux,varnames)
    implicit none
    type(algexpr), dimension(n_aux), intent(out) :: aux
    character(19), dimension(n_aux), intent(out) :: varnames
    integer      :: i
    character(2) :: num
    
    aux = (/ a_Bv2, ea_Bv2x, ea_Bv2y, ea_Bv2p, a_B2 /)
    varnames = (/ "eq(var_Bv2,0,0,0,:)", "eq(var_Bv2,1,0,0,:)", "eq(var_Bv2,0,1,0,:)", "eq(var_Bv2,0,0,1,:)", "eq( var_B2,0,0,0,:)" /)
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
