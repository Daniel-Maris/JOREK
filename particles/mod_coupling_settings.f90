!> variables and functions related to settings for the coupling between kinetic particles and the fluid
module mod_coupling_settings
use mod_model_settings, only: n_var
use phys_module, only: n_part_groups, n_part_groups_max, particle_group_configs
use phys_module, only: n_aux_var, n_diag_var
use coupling_variables

implicit none
logical :: use_ncs         = .false. !< use neutral kinetic particles
logical :: use_ccs         = .false. !< use current coupling scheme for fast particles
logical :: use_pcs         = .false. !< use pressure coupling scheme for fast particles
logical :: use_pcs_full    = .false. !< use full tensor pressure coupling scheme for fast particles
contains

    

!> Scans over all the particle groups and determines how the coupling scheme booleans should be initialized
subroutine determine_coupling_schemes()
    implicit none
    integer        :: group_num

    do group_num=1, n_part_groups
    select case (particle_group_configs(group_num)%coupling_scheme)
        case ('ncs')
        use_ncs = .true.
        case ('ccs')
        use_ccs = .true.
        case ('pcs')
        use_pcs = .true.
        case ('pcf')
        use_pcs_full = .true.
    end select
    enddo 
    ! maybe some write out here to provide info?
end subroutine determine_coupling_schemes

!> used for accumulating variables for coupling
subroutine assess_and_accumulate_variable(assessed_var, coupling_var_idx, coupling_vars)
    implicit none
    character(len=var_name_len), intent(in) :: assessed_var
    integer, intent(inout)                  :: coupling_var_idx
    character(len=var_name_len), dimension(n_aux_var_max), intent(inout) :: coupling_vars


    if (.not. (any(coupling_vars == assessed_var))) then
        coupling_var_idx = coupling_var_idx + 1
        coupling_vars(coupling_var_idx) = assessed_var
    endif

end subroutine assess_and_accumulate_variable


subroutine determine_coupling_variables()
    implicit none
    integer :: i
    integer :: coupling_var_idx, final_var_idx

    coupling_vars = ""
    coupling_var_idx = 0
    final_var_idx    = 0

    !> construct a list of unique coupling variables required
    if (use_ncs) then 
    do i=1, size(ncs_var_names)
        call assess_and_accumulate_variable(ncs_var_names(i), coupling_var_idx, coupling_vars)
    enddo
    endif
    
    if (use_ccs) then
    do i=1, size(ccs_var_names)
        call assess_and_accumulate_variable(ccs_var_names(i), coupling_var_idx, coupling_vars)
    enddo
    endif

    !use_pcs

    !use_pcs_full
    
    !> assign indices to the coupling variables and determine n_aux_var
    do i=1, coupling_var_idx
    final_var_idx = final_var_idx + 1

    select case (trim(coupling_vars(i)))
        case ("rho")
        rho_idx_kn = final_var_idx
        case ("Vpar")
        Vpar_idx_kn = final_var_idx
        case ("T")
        T_idx_kn = final_var_idx
        case ("q")
        q_idx_kn = final_var_idx
        case ("j_R")
        j_R_idx_kn = final_var_idx
        case ("j_Z")
        j_Z_idx_kn = final_var_idx
        case ("j_Phi")
        j_Phi_idx_kn = final_var_idx
        case default
        write(*,*) "Error: no match found for coupling variable, please check coupling_variables.f90 and recompile"
        stop 1
    end select
    enddo

    n_aux_var = final_var_idx
    n_aux_var = n_aux_var + 1 ! temporary as diag projections not yet created
    ! maybe some write out here to provide info?

end subroutine determine_coupling_variables



end module mod_coupling_settings