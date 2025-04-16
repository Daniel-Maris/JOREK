
!> This module handles the checking of the validity of, assigning and matching 
!> of particle group ids provided by the input namelist arrays 'part_groups_in_use' 
!> and 'part_group_configs'
module mod_particle_group_id
use phys_module, only: part_groups_in_use, n_part_groups, part_group_configs, n_part_groups_max
implicit none

integer :: id_counter = 0  !< how many ids have been automatically generated. Used to keep IDs unique.

! Index translation arrays:
! the indices in part_group_configs do not have to correspond to those in sim%groups because part_group_configs can be dropped based on part_groups_in_use
! therefore we need the following two translation arrays to get the one from the other:

!> used to translate from sim%groups(group_num) to part_group_configs%(config_num)
!> as config_num = matching_part_config_indices(group_num)
integer, dimension(n_part_groups_max) :: matching_part_config_indices = n_part_groups_max ! only the first n_part_groups value of this should be accessed

!> used to translate from part_group_configs%(config_num) to sim%groups(group_num)
!> as group_num = matching_sim_groups_indices(config_num)
integer, dimension(n_part_groups_max) :: matching_sim_groups_indices  = n_part_groups_max

private
public match_part_groups_and_configs, generate_part_groups_in_use
public matching_part_config_indices, matching_sim_groups_indices
public group_num_from_id, config_num_from_id

contains

  !> matches the groups requested in part_groups_in_use with a corresponding 
  !> group in part_group_configs. Also checks that the ids in part_groups_in_use
  !> and part_groups_configs are unique
  subroutine match_part_groups_and_configs()
    
    implicit none
    integer :: group_num, j, matched_num, matched_config_num

    do group_num=1, n_part_groups ! loop over defined groups in part_groups_in_use

      !> check that the id is unique
      if ( (count(part_groups_in_use == part_groups_in_use(group_num)) > 1) .and. part_groups_in_use(group_num) /= 'non') then
        write(*,*) "ERROR: The entries in part_groups_in_use are not unique!"
        stop
      endif
      
      !> check that a matching group exists in part_group_configs
      matched_num = 0
      matched_config_num = -1
      do j=1, n_part_groups_max ! loop over part_group_configs
        if (trim(part_group_configs(j)%id) == trim(part_groups_in_use(group_num))) then
          matched_num = matched_num + 1
          matched_config_num = j
        endif
      enddo

      if (matched_num == 0) then
        write(*,*) "ERROR: No matching part_group_configs entry found for group id: "
        write(*,*) " '", part_groups_in_use(group_num),"' defined in 'part_groups_in_use'. "
        stop
      else if (matched_num > 1) then
        write(*,*) "ERROR: More than one group in part_group_configs has the id: '", part_groups_in_use(group_num), "'"
        stop
      endif

      matching_part_config_indices(group_num) = matched_config_num
      matching_sim_groups_indices(matched_config_num) = group_num
    enddo

    !logging the matching indices
    write(*,"(A,20I4)") 'matching_part_config_indices = ', matching_part_config_indices
    write(*,"(A,20I4)") 'matching_sim_groups_indices  = ', matching_sim_groups_indices  
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

  !> returns the group_num which satisfies sim%groups(group_num)%id = id
  function group_num_from_id(sim,id) result(group_num)
    use mod_particle_sim, only: particle_sim

    type(particle_sim), intent(in) :: sim
    character(len=3),   intent(in) :: id !< particle group %id
    integer :: group_num !< the number sim%groups(group_num)
    integer :: i

    group_num = -1

    do i=1,size(sim%groups,1)
      if(sim%groups(i)%id == id) then ! matching id is found
        group_num = i
        return
      end if
    end do

    ! if at the end the matching id is not found, then the input id is not actually a valid used id in the sim
    if(sim%my_id == 0) write(*,"(3A)") "ERROR: id ",id," not found among sim%groups(:)%id (group_num_from_id)"

  end function group_num_from_id

  !> returns the group_num which satisfies part_group_configs(config_num)%id = id
  !> note that this particle group does not necessarily have to be in use!
  function config_num_from_id(id) result(config_num)
    use phys_module, only: part_group_configs

    character(len=3),   intent(in) :: id !< particle group %id
    integer :: config_num !< the number part_group_configs(config_num)%id
    integer :: i

    config_num = -1

    do i=1,size(part_group_configs,1)
      if(part_group_configs(i)%id == id) then ! matching id is found
        config_num = i
        return
      end if
    end do

    ! if at the end the matching id is not found, then the input id is not actually a valid used id in the sim
    write(*,"(3A)") "ERROR: id ",id," not found among part_group_configs(:)%id (config_num_from_id)"

  end function config_num_from_id

end module mod_particle_group_id