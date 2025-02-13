
!> This module handles the checking of the validity of, assigning and matching 
!> of particle group ids provided by the input namelist arrays 'part_groups_in_use' 
!> and 'part_group_configs'
module mod_particle_group_id
use phys_module, only: part_groups_in_use, n_part_groups, part_group_configs, n_part_groups_max
implicit none

integer :: id_counter = 0  !< how many ids have been automatically generated. Used to keep IDs unique.

!> the ith element of this array stores the matching part_group_configs index
!> for ith group in part_groups_array (i.e converts part_groups_in_use index to 
!> corresponding part_group_configs index)
integer, dimension(n_part_groups_max) :: matching_part_config_indices = 1 ! only the first n_part_groups value of this should be accessed

public matching_part_config_indices, match_part_groups_and_configs, generate_part_groups_in_use
private generate_part_group_id

contains

  !> matches the groups requested in part_groups_in_use with a corresponding 
  !> group in part_group_configs. Also checks that the ids in part_groups_in_use
  !> and part_groups_configs are unique
  subroutine match_part_groups_and_configs()
    implicit none
    integer :: i, j, matched_num, matched_idx

    do i=1, n_part_groups ! loop over defined groups in part_groups_in_use

      !> check that the id is unique
      if ( (count(part_groups_in_use == part_groups_in_use(i)) > 1) .and. part_groups_in_use(i) /= 'non') then
        write(*,*) "ERROR: The entries in part_groups_in_use are not unique!"
        stop
      endif
      
      !> check that a matching group exists in part_group_configs
      matched_num = 0
      matched_idx = -1
      do j=1, n_part_groups_max ! loop over part_group_configs
        if (trim(part_group_configs(j)%id) == trim(part_groups_in_use(i))) then
          matched_num = matched_num + 1
          matched_idx = j
        endif
      enddo

      if (matched_num == 0) then
        write(*,*) "ERROR: No matching part_group_configs entry found for group id: "
        write(*,*) " '", part_groups_in_use(i),"' defined in 'part_groups_in_use'. "
        stop
      else if (matched_num > 1) then
        write(*,*) "ERROR: More than one group in part_group_configs has the id: '", part_groups_in_use(i), "'"
        stop
      endif

      matching_part_config_indices(i) = matched_idx
    enddo
  end subroutine match_part_groups_and_configs

  !> generates part_groups_in_use based on part_groups_configs if it is
  !> not manually defined. A unique system generated ID is assigned if no
  !> input for part_groups_config%id is found.
  subroutine generate_part_groups_in_use()
    implicit none

    integer  ::  i

    do i=1, n_part_groups_max
      !> loop over actually defined configs
      if ( (part_group_configs(i)%type /= 'none') .or. (part_group_configs(i)%n_particles > 0)) then

        !> check if the ID of the group has been manually assigned
        if (trim(part_group_configs(i)%id) == 'non') then
          part_group_configs(i)%id = generate_part_group_id()
          write(*,*) "WARNING: No ID defined for particle group in slot: ", i
          write(*,*) " Assigning it the system generated ID: '", part_group_configs(i)%id, "'."

        else
        !> if manually assigned, ensure that it does not start with 'P'
        !> (Which are reserved for system generated groups IDs)
          if (part_group_configs(i)%id(1:1) == 'P') then
            write(*,*) "Error: Self assigned particle ids cannot start with 'P'  " // &
                  "as it is reserved for system assigned ids."
            stop
          endif
        endif 

        !> check if the assigned ID has been used before
        if (any(part_groups_in_use == trim(part_group_configs(i)%id))) then
          write(*,*) "ERROR: The id: '", trim(part_group_configs(i)%id), "' is being assigned to multiple particle groups."
          stop
        endif

        part_groups_in_use(i) = trim(part_group_configs(i)%id)
      endif
    enddo
  end subroutine generate_part_groups_in_use

  !> generates a unique particle group ID starting with 'P', 
  !> followed by two numerical digits
  function generate_part_group_id() result(id)
    implicit none
    character(len=3)        :: id
    character(len=2)        :: temp

    id_counter = id_counter + 1
    write(temp, '(I2.2)') id_counter
    id = 'P' // temp
  end function generate_part_group_id

end module mod_particle_group_id