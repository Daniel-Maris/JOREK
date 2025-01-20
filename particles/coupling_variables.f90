!> storage of the variables required in each of the coupling schemes
!> used to construct the indices of the feedback from the particle evolution
!> and the indices of the aux_node_list
!> PLEASE CAREFULLY OBSERVE EXISTING PARAMETERS TO AVOID OVERLAP WHEN ADDING A NEW COUPLING SCHEME
module coupling_variables
  implicit none

  integer, parameter :: n_aux_var_max = 200

  ! ====== variables names and number of coupling schemes ====== !

  integer, parameter :: var_name_len  = 15
  character(len=var_name_len), dimension(n_aux_var_max) :: coupling_vars

  ! NCS
  character(len=var_name_len), dimension(3) :: ncs_var_names = [&
    "rho",    & !> density
    "Vpar",   & !> parallel velocity
    "T"     & !> temperature
  ]

  ! CCS
  character(len=var_name_len), dimension(4) :: ccs_var_names = [&
    "q",     & !> charge density
    "j_R",    & !> R   component of current
    "j_Z",    & !> Z   component of current
    "j_Phi"   & !> Phi compoment of current  
  ]


  ! =========== Storage variables for kinetic coupling indices ======== !

  !> kn denotes kinetics
  ! variables indices
  integer :: rho_idx_kn   = 0
  integer :: Vpar_idx_kn  = 0
  integer :: T_idx_kn     = 0
  integer :: q_idx_kn     = 0
  integer :: j_R_idx_kn   = 0
  integer :: j_Z_idx_kn   = 0
  integer :: j_Phi_idx_kn = 0

end module coupling_variables