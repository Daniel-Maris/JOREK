
!> functions related to extracting model settings from particle_config
module mod_particle_config_utils
use phys_module, only: part_groups_in_use, n_part_groups, particle_configs
implicit none

integer :: id_counter = 0

contains

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