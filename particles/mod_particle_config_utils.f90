
!> functions related to extracting model settings from particle_config
module mod_particle_config_utils
use mod_model_settings, only: n_var
use phys_module, only: use_ncs, use_ccs, use_pcs, use_pcs_full
use phys_module, only: n_part_groups, n_part_groups_max, particle_configs
use phys_module, only: part_groups_in_use
use phys_module, only: n_aux_var, n_diag_var
use coupling_variables
implicit none

integer :: id_counter = 0
character(len=3), dimension(:), allocatable :: part_group_ids

contains


!> Scans over all the particle groups and determines how the coupling scheme booleans should be initialized
subroutine determine_coupling_schemes()
    implicit none
    integer        :: group_num

    do group_num=1, n_part_groups
    select case (particle_configs(group_num)%coupling_scheme)
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

subroutine assign_part_group_ids()
    implicit none
    integer :: i

    do i=1, n_part_groups
        if (particle_configs(i)%id(1:1) == 'P') then
            write(*,*) "Error: Self assigned particle ids cannot start with 'P'  " // &
                       "as it is reserved for system assigned ids."
            stop
        endif

        if (particle_configs(i)%id == 'non') then
            if (part_groups_in_use(i) /= 'non') then 
                write(*,*) "Error: part_group_in_use is defined, which requires id to be explicitly " // &
                           "defined for all members of part_configs."
                stop
            endif
            call generate_part_group_id(particle_configs(i)%id)
        endif
    enddo

    if (part_groups_in_use(1) == 'non') then ! part_group_in_use is not assigned
        do i=1, n_part_groups
            part_groups_in_use(i) = particle_configs(i)%id
        enddo
    endif

end subroutine assign_part_group_ids

subroutine generate_part_group_id(id)
    implicit none
    character(len=3), intent(inout) :: id
    character(len=2)                :: temp

    
    id_counter = id_counter + 1
    write(temp, '(I2.2)') id_counter
    id = 'P' // temp
end subroutine generate_part_group_id

end module mod_particle_config_utils