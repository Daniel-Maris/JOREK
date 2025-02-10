
!> functions related to assigning ids for particle groups
module mod_particle_group_id
use phys_module, only: part_groups_in_use, n_part_groups, particle_group_configs, n_part_groups_max
implicit none

integer :: id_counter = 0  !< how many ids have been automatically generated. Used to keep IDs unique.

contains

  !> if part_groups_in_use is manually defined, checks whether the number and order of
  !> groups specified matches with the particle_group_configs
  subroutine check_part_groups_in_use_matches_configs()
    implicit none
    integer :: i, defined_groups

    defined_groups = 0
    if (trim(part_groups_in_use(1)) /= 'non') then !< checking if part_groups_in_use is manually defined

      do i=1, n_part_groups_max
        if (part_groups_in_use(i) /= particle_group_configs(i)%id) then
          write(*,*) "ERROR: if manually defining particle groups using part_groups_in_use, the" 
          write(*,*) " groups in particle_group_configs() must match in number and id and order"
          stop
        endif

        if (part_groups_in_use(i) /= 'non') defined_groups = defined_groups + 1
      enddo

      if (n_part_groups /= defined_groups) then
        write(*,*) "ERROR: if manually defining particle groups using part_groups_in_use, the" 
        write(*,*) " groups in particle_group_configs() must match in number and id and order"
        stop
      endif

    endif
    
  end subroutine check_part_groups_in_use_matches_configs

  subroutine assign_part_group_ids()
    implicit none
    integer :: i


    ! check if part_groups_in_use is manually defined
    if (part_groups_in_use(1) == 'non') then ! not manually defined
       
      do i=1, n_part_groups
        if (particle_group_configs(i)%id == 'non') call generate_part_group_id(particle_group_configs(i)%id)
        part_groups_in_use(i) = particle_group_configs(i)%id
      enddo
    else ! is manually defined 

    endif

    do i=1, n_part_groups

      if (part_groups_in_use(i) /= particle_group_configs(i)%id) then

        if (particle_group_configs(i)%id(1:1) == 'P') then
          write(*,*) "Error: Self assigned particle ids cannot start with 'P'  " // &
                "as it is reserved for system assigned ids."
          stop
        endif

        if (particle_group_configs(i)%id == 'non') then
            write(*,*) "Error: part_groups_in_use is defined, which requires id to be explicitly " // &
                  "defined for all members of part_configs."
            stop
          
        endif

      endif
    enddo



  end subroutine assign_part_group_ids

  subroutine generate_part_group_id(id)
    implicit none
    character(len=3), intent(inout) :: id
    character(len=2)        :: temp

    id_counter = id_counter + 1
    write(temp, '(I2.2)') id_counter
    id = 'P' // temp
  end subroutine generate_part_group_id

end module mod_particle_group_id