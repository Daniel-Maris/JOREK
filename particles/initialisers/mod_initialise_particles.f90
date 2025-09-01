!> Master module which handles the creation of the initial distribution of
!> superparticles for specific particle groups the for kinetic simulations.
!> This module mainly handles the interface between the initialisers and the simulation
!> The function of the specific initialiser subroutines are contained in the 
!> "initialisers_*.f90" files
module mod_initialise_particles
  use initialisers_RE
  use initialisers_base
  use phys_module, only: part_group_configs, type_part_group_config, n_part_groups
  use mod_particle_group_id, only: matching_part_config_indices
  
  implicit none
  contains

  subroutine initialise_particles_for_sim(sim)
    implicit none
    class(particle_sim), intent(inout)       :: sim
    integer                                  :: i, j

    do i=1, n_part_groups ! loop over part_groups_in_use
      j = matching_part_config_indices(i) ! get the matching part_group_config index
      
      if (sim%my_id == 0) write(*,*) "----- Initialising particles for group '", part_group_configs(j)%id, "' with coupling scheme '", part_group_configs(j)%coupling_scheme, "' -----"

      select case(trim(part_group_configs(j)%coupling_scheme))

      !> initialisation of runaway electrons
      case('rep')
        call initialise_group_RE(sim, i)
      !> initialisation of energetic particles
      case ('epc', 'epp', 'epf')
        call initialise_group_EP(sim, i)
      !> default case - give error
      case default
        if (sim%my_id == 0) write(*,*) "ERROR : No particle coupling scheme selected for group '", part_group_configs(j)%id, "'"
        stop 1
      end select
    enddo

  end subroutine initialise_particles_for_sim

  subroutine initialise_group_RE(sim, group_num)
    use mod_pcg32_rng

    implicit none
    class(particle_sim), intent(inout) :: sim
    integer,             intent(in)    :: group_num
    character(len=50)                  :: init_function_name, init_pdf_name
    type(type_part_group_config)       :: config

    config             = part_group_configs(matching_part_config_indices(group_num))
    init_function_name = config%init_function
    init_pdf_name      = config%init_pdf

    select case (trim(init_function_name))
      case ("basic")
        !> Call initialiser subroutine
        if (sim%my_id == 0) write(*,*) "  Using the 'basic_initialization' function, with PDF: ", trim(init_pdf_name)
        call basic_initialization(sim, group_num, pcg32_rng(), init_pdf_name, config%re_energy, config%re_pitch, config%re_std_energy)

        !> Set particle charge and weight
        select type (particles => sim%groups(group_num)%particles)
        type is (particle_kinetic_relativistic)
          particles(:)%q = -1  !> default electron charge
          particles(:)%weight = config%num_re / config%n_particles
        end select

        if (sim%my_id == 0) then 
          write(*,*) "----- Finished initialisation for group '", config%id, "' with coupling scheme '", config%coupling_scheme, "' -----"
          write(*,*) ""
        endif
      case default
        if (sim%my_id == 0) then
          write(*,*) "ERROR: ", trim(init_function_name), " is not a valid initialisation function "
          write(*,*) "  for group '", config%id, "' with coupling scheme: '", config%coupling_scheme, "'"
        endif
        stop 1
    end select

  end subroutine initialise_group_RE

  subroutine initialise_group_EP(sim, group_num)
    use mod_pcg32_rng

    implicit none
    class(particle_sim), intent(inout) :: sim
    integer,             intent(in)    :: group_num
    character(len=50)                  :: init_function_name
    type(type_part_group_config)       :: config
    real*8                             :: T_maxwell
    integer                            :: n_phi_planes_in
    real*8                            :: n_particles_total

    config = part_group_configs(matching_part_config_indices(group_num))
    init_function_name = config%init_function
    T_maxwell = config%T_maxwell
    n_phi_planes_in = config%n_phi_planes
    n_particles_total = config%n_particles_total

    select case (trim(init_function_name))
    case ("maxwell")

      !> Initialise particles via a maxwellian distribution
      if (sim%my_id == 0) write(*,*) " Using the initialise_particles_h_mu_psi_phiplanes subroutine"
      call initialise_particles_H_mu_psi_phiplanes(sim%groups(group_num)%particles, sim%fields, pcg32_rng(),sim%groups(group_num)%mass, T_maxwell=T_maxwell, n_phi_planes_in=n_phi_planes_in)

      !> Adjust particle weights so sum(weights) = n_particles_total
      call adjust_particle_weights(sim%groups(group_num)%particles, n_particles_total)
      if (sim%my_id == 0) then
        write(*,*) "----- Finished initialisation for group '", config%id, "' with coupling scheme '", config%coupling_scheme, "' -----"
        write(*,'(A,ES12.2)')   " Number of super particles    : ", config%n_particles
        write(*,'(A,ES12.2)')   " Total number of particles    : ", n_particles_total
        write(*,'(A,ES12.2)')   " Particle weights adjusted to : ", sim%groups(group_num)%particles(1)%weight
        write(*,'(A,ES12.2,A)') " Temperature of Maxwellian    : ", T_maxwell, "[eV]"
        write(*,*) ""
      endif
    case default
      if (sim%my_id == 0) then
        write(*,*) "ERROR: ", trim(init_function_name), " is not a valid initialisation function "
        write(*,*) "  for group '", config%id, "' with coupling scheme: '", config%coupling_scheme, "'"
      endif
      stop 1
    end select
  end subroutine initialise_group_EP

end module mod_initialise_particles