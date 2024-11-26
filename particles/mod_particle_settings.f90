module mod_particle_settings
    use mod_model_settings

    implicit none

    ! determined by scanning and combining group settings, maybe through a config file
    logical :: use_ncs = .true. ! use neutral coupling
    logical :: use_ics = .true.  ! use impurity coupling
    logical :: use_ccs = .false. ! use current coupling scheme for fast particles
    logical :: use_pcs = .false. ! use pressure coupling scheme for fast particles
    integer :: n_aux_var = n_var
    integer :: n_diag_var = n_var
    
    ! aux_node_list indices
    integer :: var_S_rho = 1
    integer :: var_S_mom = 2
    integer :: var_S_E = 3

    ! diag_node_list indices
    integer :: var_sigma = 1



end module mod_particle_settings