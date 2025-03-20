!> Module for particle to particle and fluid to particle interactions at the wall.
!> 
!> The main object is the wall action, in which the wall_action%type is a string
!> determining the type of interaction. Examples of wall actions are a plasma fluid 
!> species sputtering one particle group (e.g. D plasma sputtering W impurities),
!> or particles from a particular group (e.g. N) reflecting against the wall
!> Every such interaction needs it's own object, and internally the right routine's
!> are then called when wall_action%do(sim) is called.
!> 
!> The currently implemented interaction types are: 
!> "self sputter" (e.g. W -> W), "fluid sputter" (e.g. fluid D+ -> W), "reflection" 
!> (e.g. kinetic D -> D) and "wall recomb" (e.g. kinetic D+ -> D)
!> 
!> Eckstein coefficients are used to determine the yield of the interaction and the 
!> resulting energy of the resulting particles. These yields are automatically loaded
!> from the simulation folder based on the original and target species symbols.
!>
!> Actions have global diagnostics (e.g. total particles of group i reflected off the 
!> wall) printed out in the logfile, and sputter diagnostics also have local vtk
!> diagnostics using mod_edge_elements
!> 
!> Limitations:
!> - The incoming angle of the particle/fluid is not taken into account (it is hardcoded 
!>   to 0). Implementing this correctly would require some estimation of surface roughness.
!> - The sampling from and integrating the fluid integrals is done using mod_edge_elements 
!>   rather than using a bezier FE description like the fluid (same is true for wall projections)
!> - A simplified model is used for the energy of sampled particles from the fluid in 
!>   fluid2part actions
module mod_particle_wall_interaction
  use mod_edge_elements
  use mod_io_actions, only: io_action
  use mod_sampling
  use mod_particle_types
  use mod_eckstein_y_ye
  use constants
  use mod_rng, only: type_rng, setup_shared_rngs
  use mod_boundary, only: wall_normal_vector
  use mod_interp
  use mod_atomic_elements !< chemical elements
  use mod_particle_sim
  use mod_event
  use mod_particle_allocation, only: calc_n_particles_per_mpi
  use equil_info, only:find_xpoint
  use mod_particle_create, only: part_create_scheme, type_part_create_scheme
  !$ use omp_lib 
  
  implicit none
   
  private
  public :: wall_action, wall_actions_from_config

  type, extends(io_action) :: wall_action
    integer           :: origin_group  !< index specifying which group is undergoing this wall interaction. Either particle group number or fluid group number (if it is a fluid to particle interaction type)
    integer           :: target_group  !< which particle group this wall interaction affects
    character(len=20) :: type = "none" !< type of the wall interaction, namely "self sputter" (e.g. W -> W), "fluid sputter" (e.g. fluid D+ -> W), "other sputter" (e.g. kinetic N -> W), "reflection" (e.g. kinetic D -> D) or "wall recomb" (e.g. kinetic D+ -> D)

    ! internal variables to determine which kind of backend function needs to be called, depending on whether the origin group is fluid or not and the target group is the origin group or not
    logical :: part2self = .false., part2other = .false., fluid2part = .false.

    type(eckstein_sputter_yield)          :: yield  !< eckstein coefficients for the wall interaction yield
    type(eckstein_sputtered_energy_coeff) :: energy !< eckstein coefficients for determining energy of the resulting particle
    type(thompson_dist)                   :: E_dist = thompson_dist(E_b = 8.7d0, n=2) !< produces energies in eV (value for W default)
    logical :: use_thompson = .false. !< Use a thompson distribution for the energy of sputtered particles
    logical :: use_Yn_func  = .false. !< Use Ecksteins interpolating functions instead of interpolating manually
    
    class(type_rng), dimension(:), allocatable :: rng !< one RNG per openmp thread
    
    !> when the origin group is a fluid species
    integer             :: fluid_Z         = -999     !< Z of this fluid species (e.g. -2 for D)
    type(edge_elements) :: fluid_yield_integral       !< the yield (of the specified interaction type) integrated over f(v) for this fluid species

    type(type_part_create_scheme) :: create_scheme    !< super particles create scheme
    real*8  :: weight_factor = 1.d0 !< additional weight factor of the yield (e.g. useful to split a single plasma fluid into D and T neutrals upon wall recombination)

    ! diagnostics
    logical             :: do_wall_projection=.true. !< whether to do wall projections for this interaction
    type(edge_elements) :: wall_projection           !< diagnostic to keep track of particle- and heat fluxes, sputtering yields, etc. resolved on the wall (1D for 2D simulations, 2D for 3D simulations)
    integer             :: i_step_diag = 0           !< how many steps have been taken between the previous diagnostic output and now
    integer             :: n_step_diag               !< after how many timesteps the wall projection should be saved (as vtk file in the simulation folder)
    real*8              :: last_diag_time = -9.d99   !< Last time of output of diagnostics
    integer             :: n_project_extra           !< Number of extra projection diagnostics for this interaction on top of n_project_general
    integer             :: n_project_tot             !< Total number of wall projections for this diagnostic (i.e. shorthand for n_project_general + n_project_extra)
    integer             :: n_project_part = -1       !< Number of particle projections (as those need to be MPI reduced)
    real*8              :: delta_t                   !< [s] tstep in SI
 
    logical             :: constructed =.false.      !< whether the constructor has been called (this is used as assert in the do action) 
  contains
    procedure :: do => do_wall_action
    procedure :: load_eckstein_data
  end type wall_action
  
  !> indices of different diagnostics in the global diagnostics array which is used for the output file
  !> number of super particles is intentionally stored in a real, to easily handle all diagnostics simultaneously (in omp reductions and in MPI_reduce)
  integer, parameter :: n_global_diagnostics=8, i_wall_part_in=1, i_wall_flux_in=2, i_wall_heat_in=3, i_wall_part_out=4, i_wall_flux_out=5, i_wall_heat_out=6, i_wall_flux_refl=7, i_wall_heat_refl=8
  
  integer, parameter :: n_project_general=4 !< number of general projections (on top of the number of interaction type specific interactions)

  real*8  :: supers_ratio_wall_default = 5.d-4                 !< if none of the wall_act_configs(i)%supers_..._wall options are set, supers_to_create will be calculated
                                                               !< as supers_ratio_wall_default * part_group_config(this%target_group)%n_particles
                                                               !< In this case this default value overrides the value from preset_parameters.f90
contains

!> Constructor for the particle_sputter type, setting the io_action parameters and sputtering parameters.
subroutine construct_wall_action(this, sim, origin_group, target_group, type, edge_element_template, fluid_Z, wall_config_num, filename, basename, decimal_digits, fractional_digits, rng, input_identifier)
  use mod_pcg32_rng, only: pcg32_rng
  use mod_random_seed, only: random_seed
  use phys_module, only: nout_projection, n_fluid_groups, n_part_groups, fluid_configs, type_wall_act_config
  use mod_particle_group_id, only: matching_part_config_indices

  implicit none
  type(wall_action),   intent(inout) :: this      !< the new wall_action object. Inout because it may need some settings already
  type(particle_sim),  intent(in) :: sim
  integer,             intent(in) :: origin_group !< index specifying which group is undergoing this wall interaction. Either particle group number or fluid group number (if it is a fluid to particle interaction type)
  integer,             intent(in) :: target_group !< which particle group this wall interaction affects
  character(len=*),    intent(in) :: type         !< type of the wall interaction, namely "self sputter" (e.g. W -> W), "fluid sputter" (e.g. fluid D+ -> W), "other sputter" (e.g. kinetic N -> W), "reflection" (e.g. kinetic D -> D) or "wall recomb" (e.g. kinetic D+ -> D)
  type(edge_elements), intent(in) :: edge_element_template !< a prepared set of edge elements
  integer,             intent(in), optional :: fluid_Z          !< Z of this fluid species (e.g. -2 for D)
  integer,             intent(in), optional :: wall_config_num  !< which wall_act_config this action is based on, [fluid]/[part_groups]_configs(i)%wall_act_configs(wall_config_num)
  character(len=*),    intent(in), optional :: filename         !< where to save the diagnostics
  character(len=*),    intent(in), optional :: basename
  integer,             intent(in), optional :: decimal_digits
  integer,             intent(in), optional :: fractional_digits
  class(type_rng),     intent(in), optional :: rng              !< random-number generator to use (default PCG32)
  character(len=*),    intent(in), optional :: input_identifier !< extra message on stops, to determine which construct_wall_action call had wrong input
 
  character(len=100) :: name
  integer :: my_seed, i, j, target_group_loc
  character(len=14), dimension(:), allocatable :: extra_proj_scalar_names !< additional scalar names on top of the normal ones
  character(len=1000) :: msg !< error message
  character(len=1000) :: identifier
  type(type_wall_act_config) :: config

  ! setting the identifier (it is optional for the sake of using construct_wall_action directly from inside a program rather than through the namelist)
  if(present(input_identifier)) then
    identifier = input_identifier
  else
    identifier = ""
  end if

  ! --- determining the interaction type
  this%type = trim(type)

  select case(trim(type))
  case("self sputter")
    this%part2self = .true.
  case("fluid sputter")
    this%fluid2part = .true.
  case("other sputter")
    this%part2other = .true.
    write(msg,*) 'type "other sputter" is still to be implemented'
    call wrong_input(msg, sim%my_id, identifier)
  case("reflection")
    this%part2self = .true.
  case("wall recomb")
    this%fluid2part = .true.
  case default
    call wrong_interaction_type(type, identifier)
  end select
  
  ! --- general checks on input
  if (this%fluid2part) then
    if(origin_group < 1 .or. origin_group > n_fluid_groups) then
      write(msg,*) "fluid origin group is not valid (origin_group/max)",origin_group,n_fluid_groups
      call wrong_input(msg, sim%my_id, identifier)
    end if
    this%origin_group = origin_group

    if(present(fluid_Z)) then
      ! check whether Z is sensical
      if(fluid_Z < lbound(element_symbols,1) .or. fluid_Z > ubound(element_symbols,1)) then
        write(msg,"(A,3I5)") "fluid Z not in bound (Z/min/max)",fluid_Z,lbound(element_symbols,1),ubound(element_symbols,1)
        call wrong_input(msg, sim%my_id, identifier)
      end if
      this%fluid_Z = fluid_Z
    else ! fluid Z must be present
      write(msg,"(A,3I5)") "fluid Z must be specified for fluid type interaction"
      call wrong_input(msg, sim%my_id, identifier)
    end if
  else ! origin_group is a particle group
    if(origin_group < 1 .or. origin_group > n_part_groups) then
      write(msg,*) "particle origin group is not valid (origin_group/max)",origin_group,n_part_groups
      call wrong_input(msg, sim%my_id, identifier)
    end if

    this%origin_group = matching_part_config_indices(origin_group)
  end if
  
  !self interactions default to have the same target as origin if the wall_act_config%target_group is the unchanged namelist input value -1
  if(this%part2self .and. target_group == -1) then
    target_group_loc = origin_group
  else
    target_group_loc = target_group
  end if

  !checking and setting target group
  if(target_group_loc < 1 .or. target_group_loc > n_part_groups) then
    write(msg,*) "target group is not valid (target_group/max)",target_group_loc,n_part_groups
    call wrong_input(msg, sim%my_id, identifier)
  end if
  this%target_group = matching_part_config_indices(target_group_loc)

  !checking whether the user set origin group and target group differently while it is a self interaction
  if (this%part2self) then
    call check_self_type(this, sim%my_id, identifier)
  end if

  !setting the creation scheme
  if (this%fluid2part) then
    if(.not. present(wall_config_num)) call wrong_input("wall_config_num should be present for particle creating routines in construct_wall_action", sim%my_id, identifier)
    config = fluid_configs(this%origin_group)%wall_act_configs(wall_config_num)
    this%create_scheme = part_create_scheme(config%supers_num_wall,config%supers_weight_wall,config%supers_ratio_wall, &
        sim%groups(target_group)%n_particles,default=supers_ratio_wall_default,my_id=sim%my_id,identifier=identifier)
  end if

  ! --- diagnostics
  if(this%fluid2part) extra_proj_scalar_names = ["n_e           ","T_e           ","cos_alpha     ","Psi_n         ","fluid_flux    ","fluid_heatflux","fluid_yield   "]
  
  ! if there are no extra projections, set the allocatable to 0
  if(.not. allocated(extra_proj_scalar_names)) allocate(extra_proj_scalar_names(0))
  if(this%n_project_part < 0) this%n_project_part = n_project_general

  call this%load_eckstein_data(sim)

  ! initialising the edge_element objects from the template
  if (.not. allocated(edge_element_template%patch(1)%xyz)) then
    write(msg,*) 'Edge element template needs to be prepared, exiting'
    call wrong_input(msg, sim%my_id, identifier)
  end if
  
  this%wall_projection = edge_element_template
  this%fluid_yield_integral = edge_element_template
  ! Clean up the passed edge elements
  do i=1,size(edge_element_template%patch,1)
    if (allocated(edge_element_template%patch(i)%scalars)) then
      deallocate(this%wall_projection%patch(i)%scalars, &
                 this%fluid_yield_integral%patch(i)%scalars)
    end if
    if (allocated(edge_element_template%patch(i)%scalar_names)) then
      deallocate(this%wall_projection%patch(i)%scalar_names, &
                 this%fluid_yield_integral%patch(i)%scalar_names)
    end if
  end do

  ! settings for the diagnostics
  write(name, "(A,I2.2,A,I2.2)") spaces2underscore(type)//"_", origin_group, "_to_", this%target_group
  this%n_step_diag = nout_projection
  this%basename = trim(name)//"_"
  if (present(filename)) this%filename = filename
  if (present(basename)) this%basename = basename
  if (present(decimal_digits)) this%decimal_digits = decimal_digits
  if (present(fractional_digits)) this%fractional_digits = fractional_digits
  this%extension = '.vtk'
  this%name = trim(name)
  this%log = .true.

  ! Set up scalars and scalar names for the diagnostic projections
  this%n_project_extra = size(extra_proj_scalar_names, dim=1)
  this%n_project_tot = n_project_general + this%n_project_extra
  do i=1,size(this%wall_projection%patch,1)
    allocate(this%wall_projection%patch(i)%scalars(size(this%wall_projection%patch(i)%st,2), this%n_project_tot))
    this%wall_projection%patch(i)%scalars(:,:) = 0.d0 ! initialising scalars
    
    allocate(this%wall_projection%patch(i)%scalar_names(this%n_project_tot))
    
    ! defining the scalar names
    associate (sn => this%wall_projection%patch(i)%scalar_names)
      sn(1) = "part_flux"
      sn(2) = "part_heatflux"
      sn(3) = "part_promptflux"
      sn(4) = "part_yield"
      do j=1,this%n_project_extra
        sn(n_project_general+j) = trim(extra_proj_scalar_names(j))
      enddo
    end associate
  end do

  ! Allocate scalars for the fluid_yield_integral
  do i=1,size(this%fluid_yield_integral%patch,1)
    allocate(this%fluid_yield_integral%patch(i)%scalars( &
      size(this%fluid_yield_integral%patch(i)%st,2), 1))
    
    this%fluid_yield_integral%patch(i)%scalars = -1
  end do

  ! --- allocate random seed for sampling
  my_seed = random_seed()
  if (present(rng)) then
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=rng, rngs=this%rng)
  else
    ! default to pcg32_rng
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=pcg32_rng(), rngs=this%rng)
  end if

  !constructor finished
  this%constructed = .true.
end subroutine construct_wall_action


!> set up the wall actions array from the configs of the namelist
function wall_actions_from_config(sim, edge_element_template) result(wall_actions)
  use phys_module, only: part_group_configs, n_part_groups, n_part_groups_max, type_wall_act_config
  use phys_module, only: fluid_configs, n_fluid_groups, n_fluid_groups_max
  use mod_particle_group_id, only: matching_part_config_indices

  implicit none

  type(particle_sim),  intent(in) :: sim
  type(edge_elements), intent(in) :: edge_element_template !< a prepared set of edge elements
  type(wall_action), dimension(:), allocatable :: wall_actions

  type(type_wall_act_config) :: config

  character(len=1000) :: identifier
  integer :: i, j, Z, config_num_i, config_num_j
  integer :: n_wall_acts !< number of wall_action objects to make
  integer :: i_wall_acts !< ith wall_action to be put in the wall_actions

  !determining number of wall actions necessary
  n_wall_acts = 0
  
  do i=1,n_part_groups !from particles
    do j=1,n_part_groups_max
      config = part_group_configs(i)%wall_act_configs(j)
      if(trim(config%type) == "none") cycle

      n_wall_acts = n_wall_acts + 1
    end do
  end do
  
  do i=1,n_fluid_groups_max !from the fluid
    do j=1,n_part_groups_max
      config = fluid_configs(i)%wall_act_configs(j)
      if(trim(config%type) == "none") cycle

      n_wall_acts = n_wall_acts + 1
    end do
  end do

  if(n_wall_acts > n_part_groups_max**2 + n_part_groups_max*n_fluid_groups_max) then
    write(*,*) "ERROR: number of wall actions is bigger than should be possible?"
    stop
  end if

  !setting up the wall actions
  allocate(wall_actions(n_wall_acts))

  i_wall_acts = 0
  
  !from particles
  do i=1,n_part_groups !loop over particle groups
    do j=1,n_part_groups_max !loop over wall action configs
      config_num_i = matching_part_config_indices(i)
      config_num_j = matching_part_config_indices(j)
      config = part_group_configs(config_num_i)%wall_act_configs(config_num_j)
      if(trim(config%type) == "none") cycle

      ! being here means it is a wall action that should be used
      i_wall_acts = i_wall_acts + 1
      write(identifier,"(A,I2,A,I2,A,I3,A)") " for input namelist: particle_group_config(",config_num_i,")%wall_act_configs(",config_num_j,"). (This corresponds to wall_action: ",i_wall_acts,")"
      call construct_wall_action(wall_actions(i_wall_acts),sim,i,config%target_group,config%type,edge_element_template,input_identifier=identifier)      
    end do
  end do

  !from the fluid
  do i=1,n_fluid_groups_max !loop over fluid groups
    Z = fluid_configs(i)%Z
    
    do j=1,n_part_groups_max !loop over wall action configs
      config_num_j = matching_part_config_indices(j)
      config = fluid_configs(i)%wall_act_configs(config_num_j)
      if(trim(config%type) == "none") cycle

      ! being here means it is a wall action that should be used
      i_wall_acts = i_wall_acts + 1
      write(identifier,"(A,I2,A,I2,A,I3,A)") " for input namelist: fluid_config(",i,")%wall_act_configs(",config_num_j,"). (This corresponds to wall_action: ",i_wall_acts,")"
      call construct_wall_action(wall_actions(i_wall_acts),sim,i,config%target_group,config%type,edge_element_template,fluid_Z=Z,wall_config_num=j,input_identifier=identifier)
    end do
  end do

end function wall_actions_from_config


!> Load eckstein sputtering yields and energy coefficients for this interaction
subroutine load_eckstein_data(this, sim)
  implicit none

  class(wall_action), intent(inout) :: this
  type(particle_sim), intent(in)    :: sim

  integer :: Z_origin, Z_target

  ! determining Z of origin and target (needed to read the correct data file)
  if (this%fluid2part) then
    Z_origin = this%fluid_Z
  else !< origin is particle group
    Z_origin = sim%groups(this%origin_group)%Z
  end if

  Z_target = sim%groups(this%target_group)%Z

  ! setting the yield object
  this%yield%Z_ion    = Z_origin
  this%yield%Z_target = Z_target
  this%yield%use_Yn_func = this%use_Yn_func

  ! reading the yield data
  call this%yield%read()

  if (.not. this%use_thompson) then ! use eckstein coefficients
    ! setting the energy object
    this%energy%Z_ion    = Z_origin
    this%energy%Z_target = Z_target
    this%energy%use_Yn_func = this%use_Yn_func

    !reading the energy data
    call this%energy%read()
  end if
end subroutine load_eckstein_data


!> Perform the wall interaction according to the setting in the wall_action object
!> How and when to run depends on the details of the interaction
subroutine do_wall_action(this, sim, ev)
  use mod_atomic_elements, only: element_symbols
  use mod_parameters, only: n_plane, n_period
  use phys_module, only: use_manual_random_seed, tstep, central_mass, central_density
  
  class(wall_action), intent(inout)    :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev   !< this is here so that it is compatible with the event structure, but we can remove it if we want to move away from it

  integer :: i

  if(sim%my_id == 0) write(*,"(A)") "--- wall action: "//trim(this%name)//" --- "

  ! --- setup
  ! updating the timestep in SI (as tstep can change)
  this%delta_t = (tstep*sqrt((MU_ZERO * CENTRAL_MASS * MASS_PROTON * CENTRAL_DENSITY * 1.d20)))

  ! check whether the constructor was used (so that all other sanity checks can be done once in the constructor)
  if (.not. this%constructed) then
    write(*,*)'=======================ERROR!!=================================='
    write(*,*)"particle wall action object of type "//trim(this%type)
    write(*,"(A,I2,A,I2, A)") "with origin ",this%origin_group," and target ", this%target_group, "has not finished it's construction"
    write(*,*)'please use the constructor'
    stop
  end if

  ! --- underlying function calls
  if(this%fluid2part) call fluid2part_action(this, sim)

  if(this%part2self) call part2self_action(this, sim)
  
  ! in the future add part2other

  ! --- area for writing the projected diagnostic
  this%i_step_diag = this%i_step_diag + 1
  if (this%i_step_diag .ge. this%n_step_diag) then
    call write_wall_project_vtk(this, sim)
  end if
end subroutine do_wall_action


!> Routine to do the fluid to particle wall interactions.
!> Models fluid sputtering ("fluid sputter") and wall 
!> recombination ("wall recomb")
!>
!> The fluid-particle interaction is done by calculating 
!> the yield by integrating over the velocity distribution.
!> Then particles are sampled using these local yields to represent
!> the incoming flux (with the weight already adjusted for particles 
!> resulting from the interaction). 
!> The particles then undergo a particle-particle interaction to 
!> determine their new energy, and the new particles are stored in
!> free slots in the group's particle array
subroutine fluid2part_action(this, sim)
  use mpi_mod
  use mod_atomic_elements, only: element_symbols
  use mod_interp, only: interp_RZ
  use mod_particle_create, only: free_particle_indices
  use mod_particle_types, only: initialize_particle_to_zero
  use mod_edge_elements, only: sample_edge_elements, integrate_edge_elements, type_cdf_data
  
  class(wall_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  
  type(type_cdf_data) :: res !< data on the cumulative distribution function calculated at the integral which is needed when sampling
  type(particle_kinetic_leapfrog) :: particle
  integer :: j, i_p, n_supers, n_supers_loc, i_rng
  integer :: q, Z
  real*8 :: E !< [eV] particle energy  (eV because of eckstein coeffs).
  real*8 :: n_e, T_e, T_eV
  real*8 :: integral !< total weight of all created particles in this wall action
  real*8,  allocatable :: xyz_sampled(:,:), st_sampled(:,:), rng_sample(:,:) !< (3,n_supers_loc), (2,n_supers_loc), (3,n_supers_loc)
  integer, allocatable :: i_elm_sampled(:) !< (n_supers_loc)

  !> For check free particles
  integer, allocatable, dimension(:) :: i_free
  integer :: n_free

  ! diagnostics
  real*8, dimension(n_global_diagnostics) :: diagnostics         !< diagnostics for the global wall loads
  real*8, dimension(n_global_diagnostics) :: diagnostics_all_mpi !< MPI reduced version of diagnostics
  real*8  :: mol_binding_E=3.526d-19, ion_binding_E=2.18d-18  !< ! (J) default values are only true for hydrogen, should be the sum of ionisation energies from 0 to q for impurities.
  
  !> this subroutine will calculate the incident ion flux over every fluid species on edge domain
  !> And the resulting yield of created particles (in atoms/m^2 during delta_t)
  call project_sputter_vars_on_edge(this, sim)

  ! determine integral over the domain
  call integrate_edge_elements(this%fluid_yield_integral, 1, integral, res)

  if (integral .le. 1d-12) then
    if(sim%my_id == 0) write(*,"(A,I2,A)") "fluid2part wall action "//trim(this%type)//" with Z=",this%fluid_Z," has 0 yield. returning"
    return ! Move along, nothing to do
  end if

  ! determine how many particles to initialise on this MPI proces
  n_supers = this%create_scheme%supers_to_create(sim%my_id,integral*this%weight_factor)
  n_supers_loc = calc_n_particles_per_mpi(n_supers, sim%n_mpi, sim%my_id)

  ! determine indices of free particles
  call free_particle_indices(sim%groups(this%target_group)%particles, n_free, i_free, n_needed=n_supers_loc)  

  allocate(rng_sample(3,size(i_free)))
  allocate(xyz_sampled(3,size(i_free)))
  allocate(st_sampled(2,size(i_free)))
  allocate(i_elm_sampled(size(i_free)))

  diagnostics = 0.d0

  ! We need to properly use all RNGS here to avoid missing numbers
  ! needs default(shared) for gfortran
#ifdef __GFORTRAN__
  !$omp parallel default(shared) &
#else
  !$omp parallel default(none) &
  !$omp shared(this, rng_sample, n_supers_loc) &
#endif 
  !$omp private(i_rng, j)
  i_rng = 1
  !$ i_rng = omp_get_thread_num()+1
  !$omp do schedule(static,1)
  do j=1,n_supers_loc
    call this%rng(i_rng)%next(rng_sample(:,j))
  end do
  !$omp end do
  !$omp end parallel

  q = this%fluid_Z
  if (q .le. 0) q = 1 ! deuterium, tritium special case
  q = min(q, 4) ! limit to 4 for divertor conditions
  Z = this%fluid_Z

  call sample_edge_elements(this%fluid_yield_integral, res, 1, n_supers_loc, rng_sample(1:2,:), xyz_sampled, st_sampled, i_elm_sampled)

  if (sim%my_id .eq. 0) then
    write(*,"(A,i8,A,A,A,A,A,i2,A,i2,A,es16.6,A,es16.6)") "fluid2wall will create ", n_supers," ", element_symbols(sim%groups(this%target_group)%Z),&
      " from ", element_symbols(Z), " in group ", this%target_group, &
    " (Z=", sim%groups(this%target_group)%Z, ") with total weight ", integral*this%weight_factor, "  particles flux #/s : ", integral*this%weight_factor/this%delta_t
  end if

  select type (pa => sim%groups(this%target_group)%particles)
  type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
  !$omp parallel default(shared) &
#else
  !$omp parallel default(none) &
  !$omp shared(this, sim, rng_sample, xyz_sampled, st_sampled, i_elm_sampled, i_free, &
  !$omp integral, q, Z, n_supers_loc) &
#endif
  !$omp private(i_rng, j, E, T_e, T_eV, n_e, particle, i_p) &
  !$omp reduction(+:diagnostics)
  i_rng = 1
  !$ i_rng = omp_get_thread_num()+1
  !$omp do schedule(static,1)
  do j=1,n_supers_loc
    !> make a new particle which at the end of the do loop will be written into a free particle in the array
    call initialize_particle_to_zero(particle)

    particle%q = int(q,1)
    particle%i_elm = i_elm_sampled(j)
    if (i_elm_sampled(j) .le. 0) cycle
    particle%st = st_sampled(:,j)
    call interp_RZ(sim%fields%node_list, sim%fields%element_list, i_elm_sampled(j), &
      st_sampled(1,j), st_sampled(2,j), &
      particle%x(1), &
      particle%x(2))
    particle%x(3) = xyz_sampled(3,j) ! phi coordinate from sampling

    !> weight of fluid particle is equally distributed as a fraction of the incoming flux. Such that the sum of all incoming fluid particles ,
    !> is the total amount of incoming particles over the edge domain area * delta_t
    !> multiplication with this%weight will be done in single_self_interaction
    particle%weight = integral/(n_supers_loc*sim%n_mpi)

    ! Calculate temperature at this position to determine particle energy
    call sim%fields%calc_NeTe(sim%time, i_elm_sampled(j), st_sampled(:,j), xyz_sampled(3,j), n_e, T_e)
    T_eV = T_e * K_BOLTZ / EL_CHG

    select case(trim(this%type))
    case("wall recomb")
      ! determine E
      call sample_fluid_particle_energy(T_eV, rng_sample(1:3,j), Z, E)

      ! determine outcoming particle
      call single_self_interaction(this, sim, particle, this%rng(i_rng), diagnostics, E, "reflection")
    case("fluid sputter")
      ! The yield at a specific position is given by
      ! \[
      !   \int_v Y(E) f(v) dv
      ! \]
      ! where $f(v)$ is a maxwellian and $E$ includes the sheath potential and the Bohm outflow
      ! condition additionally.
      !
      ! To now calculate the energy of the sputtered particle we multiply the sputtered energy
      ! coefficient with E of a particle sampled from f(v).
      ! Taking the sputtered energy coefficient * Y as a weight factor and sampling from f(v) will do the trick.
      ! we need to normalize with the sputter yield at that position, which we have calculated before.
      ! Basically this is a weighted average of Y_E(E) * E, weighted with Y(E) f(E)
      
      ! If sampling from the incoming energy distribution function, the
      ! sputtered energy coefficient needs to be reweighed with the sputtering
      ! yield at this energy (since the tail contributes more)
      ! This is commented below since we have simplified the model for now to
      ! work at a fixed energy of 3 q T_e + 2 T_i, so we don't need to do this
      ! anymore. The extension to realistic IEDFs should be done later, so 
      ! I've kept some of the code around.

      E = 2 * T_eV !< from the bohm criterion, E = E_sheath_entrance + E_sheath_acceleration = 2 T_i + 3 q T_e, but for now T_i = T_e
      ! so E_sheath_entrance  = 2 T_i = 2 T, and E_sheath_acceleration will be added later
      call single_self_interaction(this, sim, particle, this%rng(i_rng), diagnostics, E, "self sputter", .true.)

      ! Non-implemented alternative to the above method:
      ! We could sample directly from Y(E) Y_E(E) f(E), but I don't know how to do this generally.
      ! That would have the advantage of better distribution of statistics (more uniform weights).

      !sputtering_yield = this%yield%interp(E, theta)
      ! Workaround if sputtered energy coeff threshold is lower than sputtering
      ! threshold: use sputtered energy coeff just above threshold instead
      ! (note: all this doesn't take into account theta properly)

      !av_yield = fluid_sputtering_yield(this%yield, T_eV, Z, theta)
      ! we could probably avoid the calculation of fluid_sputtering_yield by
      ! using the discretisation we just sampled from (if theta is constant)
      !if (av_yield .le. 1d-18) av_yield = 1d-6 ! does not matter since then sputtering_yield must be 0, just to avoid a NaN below
      
      ! now we weigh the particles with the prevalence of this energy in sputtered particles, i.e. sputtering_yield
      ! over the integral of sputtering_yield, which we calculate (again)
      !particle%weight = &
      !particle%weight * &
        !sputtering_yield / av_yield
    case default
      call wrong_interaction_type(this%type)
    end select

    ! write the created particle to a free slot in the array
    i_p = i_free(j)
    pa(i_p) = particle ! assignment(=+ operator is defined for particle_base as copy, so this works as you would intuitively think
  end do
  !$omp end do
  !$omp end parallel
  class default
    write(*,*) 'Target particle type not implemented for fluid2part actions (origin/target)', this%origin_group, this%target_group
    stop
  end select

  call write_global_diag(this, sim, diagnostics)

  ! on top of the standard global diagnostics, a few extra ones are printed
  if (sim%my_id .eq. 0) then
    write(*,'(A50,1es16.8)') "atom wall-assisted recombination power [W] = ",      diagnostics(i_wall_flux_in)                                  * ion_binding_E / this%delta_t
    write(*,'(A50,1es16.8)') "molecule wall-assisted recombination power [W] = ", (diagnostics(i_wall_flux_in) - diagnostics(i_wall_flux_refl)) * mol_binding_E / this%delta_t
    write(*,'(A50,1es16.8)') "Power to (fast) reflected atoms [W] = ",             diagnostics(i_wall_heat_refl)                                                / this%delta_t
  endif
  !< atom wall-assisted recombination power = recycled flux * 13.6 eV. All ions are neutralized on the wall. This increaes the heat load on the wall ~stangeby2000 p.653
  !< molecule wall-assisted recombination power = thermal desorption flux*2.2 eV. When neutrals on the wall form neutrals, the wall heat load is increased by 2.2 eV per molecule. ~stangeby2000 p.653
  !< power to (fast) reflected atoms = energy retained by reflected neutrals. This energy is not deposited on the wall, thus decreases the plasma heat load.
  !  From ITER PFPO-1 test in 2D : energy_reflected_all > enery_wall_recombi_all >> energy_mol_recombi_all

  deallocate(rng_sample, xyz_sampled, st_sampled, i_elm_sampled)
  deallocate(i_free)
end subroutine fluid2part_action


!> calls single_self_interaction() for all particles in the specified this%origin_group
!> also prints the global diagnostics to the output file
subroutine part2self_action(this, sim)
  use phys_module, only: use_manual_random_seed
  
  class(wall_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  
  real*8, dimension(n_global_diagnostics) :: diagnostics         !< diagnostics for the global wall loads
  real*8, dimension(n_global_diagnostics) :: diagnostics_all_mpi !< MPI reduced version of diagnostics

  integer :: j, i_rng
  
  diagnostics = 0.d0

  select type (pa => sim%groups(this%origin_group)%particles)
  type is (particle_kinetic_leapfrog)
  if(use_manual_random_seed) then
    !$ call omp_set_schedule(omp_sched_static,10)
  else
    !$ call omp_set_schedule(omp_sched_dynamic,10)
  end if
  
#ifdef __GFORTRAN__
  !$omp parallel default(shared) & ! workaround for Error: ‘__vtab_mod_pcg32_rng_Pcg32_rng’ not specified in enclosing ‘parallel’
#else
  !$omp parallel default(none) &
  !$omp shared(this, sim)      & 
#endif
  !$omp private(i_rng)         &
  !$omp reduction(+:diagnostics)

  i_rng = 1
  !$ i_rng = omp_get_thread_num()+1
  !$omp do schedule(runtime)
  do j = 1,size(sim%groups(this%origin_group)%particles,1)
    ! Skip if this particle is not lost in a specific location (i_elm .eq. 0 means lost 'somewhere')
    if (pa(j)%i_elm .ge. 0) cycle
    
    !> Place particle back into domain
    pa(j)%i_elm = -pa(j)%i_elm
    
    !> do single particle wall interaction
    call single_self_interaction(this, sim, pa(j), this%rng(i_rng), diagnostics)
  end do
  !$omp end do
  !$omp end parallel
  class default
    write(*,*) "part2self_action not implemented for this kinetic type, group=",this%origin_group
    call exit(13)
  end select
  
  call write_global_diag(this, sim, diagnostics)
end subroutine part2self_action


!> The interaction of a single particle with the wall, only affecting that super particle (self sputter or reflect)
subroutine single_self_interaction(this, sim, particle, rng, diagnostics, E_in, type_in, weight_preadjusted)
  implicit none

  class(wall_action),                      intent(inout) :: this
  type(particle_sim),                      intent(in)    :: sim
  type(particle_kinetic_leapfrog),         intent(inout) :: particle           !< particle to undergo interaction
  class(type_rng),                         intent(inout) :: rng                !< RNG object of the current openmp thread
  real*8, dimension(n_global_diagnostics), intent(inout) :: diagnostics        !< diagnostics for the global wall loads
  real*8,            optional,             intent(in)    :: E_in               !< [eV] energy of the incoming particle (if not specified will be determined from particle%v)
  character(len=*),  optional,             intent(in)    :: type_in            !< type of single interaction (either "self sputter" or "reflection") (if not specified will be set to this%type) 
  logical,           optional,             intent(in)    :: weight_preadjusted !< whether the weight was already adjusted beforehand to take the yield into account (true) or not (false, default)

  real*8 :: n_e, T_e, theta
  real*8 :: E !<[eV] particle energy. E is in [eV] in this subroutine, because of eckstein coeffs.
  real*8 :: vector_normal(3)
  logical :: fast_reflection !< whether the reflection is a fast reflection or a thermal desorption (not that release is instant, but the energy of the reflected particle is different)
  real*8 :: yield, energy_coeff, T_eV, fast_reflect_chance, v_new
  real*8 :: u(2)
  character(len=20) :: local_type !< which single particle interaction to do, used to call self interaction from within fluid2part_action (=type_in if present, else =this%type)
  logical :: skip_yield !< if weight_preadjusted = true, then the yield calculation should be skipped
  
  ! set the incoming particle energy
  if(present(E_in)) then
    E = E_in
  else
    ! calculate the energy associated with the velocity of the particle (in eV)
    E = 0.5d0*sim%groups(this%origin_group)%mass*ATOMIC_MASS_UNIT*dot_product(particle%v, particle%v)/EL_CHG !< must be in eV
  end if

  ! determine the type, this can be different from this%type if single_self_interaction is called from within fluid2part
  if(present(type_in)) then
    local_type = type_in
  else
    local_type = this%type
  end if

  ! determine whether to calculate the yield
  if(present(weight_preadjusted)) then
    skip_yield = weight_preadjusted
  else
    skip_yield = .false.
  end if

  ! use normal vector and velocity of particle to determine incoming angle
  ! cos(theta) = (n . v)/ (||n||.||v||)
  vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, particle%i_elm, particle%st(1), particle%st(2))

  ! Hard-code theta to 0 to fix issues with sputtering module at strange angles
  ! the angle calculation should be revisited. Before using theta != 0 the
  ! surface roughness should be estimated, as this gives a distribution of
  ! impact angles as well
  theta = 0.d0 
  !> old theta 
  ! theta = acos(dot_product(-vector_normal,particle%v)/norm2(particle%v))*180.d0/PI !< acos gives results in radians
  ! ! theta must be in degrees as the theta_star is also in degrees
  ! if (abs(theta) .gt. 91) then
  !   ! This is like an assert, it cannot really happen... but it does
  !   !!$omp critical
  !   !write(*,*) 'incoming angle warning', theta, vector_normal, particle%v
  !   !!$omp end critical
  ! end if

  ! Update the particle energy from the potential drop in the sheath
  call sim%fields%calc_NeTe(sim%time, particle%i_elm, particle%st, particle%x(3), n_e, T_e)
  T_eV = T_e*K_BOLTZ/EL_CHG
  E = E + simple_potential_drop(int(particle%q,4),T_eV)

  ! store this particle's contribution to incoming particle, heatflux and flux onto the wall
  diagnostics(i_wall_part_in) = diagnostics(i_wall_part_in) + 1
  diagnostics(i_wall_flux_in) = diagnostics(i_wall_flux_in) + particle%weight
  diagnostics(i_wall_heat_in) = diagnostics(i_wall_heat_in) + particle%weight * E *EL_CHG

  !> determining interaction yield and new energy depending on interaction type
  select case (trim(local_type))
  case ("reflection")
    !> a particle can either bounce of the wall (fast_reflection=.true.) or be thermally released
    !> whether a particle reflects directly is determined through eckstein coefficients set for this goal
    fast_reflect_chance = this%yield%interp(E,theta)
    
    call rng%next(u)
    if (u(1) .le. fast_reflect_chance) then
      fast_reflection = .true.
    else
      fast_reflection = .false. 
    end if
    
    !> assume wall saturation (pumping implementation is done separately)
    yield = 1.d0

    !> storing this particle's contribution on a 2D edge element patch grid as diagnostic
    call particle_projection_diagnostic(this, sim, particle, E, yield)

    !> determine new energy
    if (fast_reflection) then
      ! still some energy and momentum can be lost at the reflection against the wall, this is modelled using another set of eckstein coefficients
      energy_coeff = this%energy%interp(E,theta)
      E = energy_coeff * E

      ! since we have wall_flux_in, and wall_flux_in = wall_flux_refl + wall_flux_therm, we also know wall_flux_thermal. Similarly we know wall_heat_thermal
      diagnostics(i_wall_flux_refl)   = diagnostics(i_wall_flux_refl) + particle%weight
      diagnostics(i_wall_heat_refl)   = diagnostics(i_wall_heat_refl) + particle%weight * E * EL_CHG
    else ! thermal release
      E = (800.d0 + 273.d0) *K_BOLTZ/EL_CHG! must be in eV (800 degrees celsius)
    endif  
  case ("self sputter")
    ! use eckstein sputtering coefficients to determine both the sputter yield and resulting energy
    
    if(skip_yield) then
      yield = 1.d0
    else
      yield = this%yield%interp(E,theta)
    end if

    !> exponential self sputtering for yield > 1
    if (yield .gt. 1.d0 + 1.d-12) then
      !$omp critical
      write(*,"(A,f5.0,A,f8.3)") "> 1 self-sputtering detected, E=", E, "yield=", yield
      !$omp end critical
    end if

    !> storing this particle's contribution on a 2D edge element patch grid as diagnostic
    call particle_projection_diagnostic(this, sim, particle, E, yield)

    !> determining the energy of the particle post sputtering
    if (this%use_thompson) then
      call rng%next(u)
      ! Option below to remove the highest 2% of the distribution by clipping u (hacky)
      ! u = min(u, 0.98d0)
      E = sample_dist(this%E_dist, u(1))
    else
      !> avoiding numerical issues with E being too small to calculate energy_coeff
      if (E < this%energy%E_threshold + 1d0) then
        !$omp critical
        write(*,*) "WARNING: E too small for yields",E,this%energy%E_threshold,"setting E to just above threshold, please expand coefficients range"
        !$omp end critical
        E = this%energy%E_threshold + 1d0
      end if

      energy_coeff = this%energy%interp(E,theta)
      E = energy_coeff * E
    end if

  case default
    write(*,*) "ERROR: unknown single_self_interaction type",local_type
  end select
  
  ! update weight of simulated particle after the wall interaction
  particle%weight = this%weight_factor * yield * particle%weight 

  ! use E from previous section to calculate velocity in one 
  v_new = sqrt(2.d0* E *EL_CHG/(sim%groups(this%origin_group)%mass * ATOMIC_MASS_UNIT)) !< TODO: is this wrong, shouldn't it be sqrt(3 kb T / m)? (leave for now for regtest)
  
  ! give particle a new direction:
  ! Calculate vector normal and select a random vector with a cosine distribution in angle between the normal and itself
  call rng%next(u)
  particle%v =  v_new * sample_cosine(u(1:2),vector_normal) 
  ! [[not sure what this comment is about]] Since it is a neutral the half-step for boris method does not matter at all

  ! wall interactions typically neutralise the particles if they used to have charge
  particle%q = 0_1

  ! after the wall interaction, the particle is now considered a new particle, so update i_life and t_birth
  particle%i_life = particle%i_life + 1
  particle%t_birth = sim%time
  ! For particle-particle sputtering we might want them to have the same identifiers
  ! if so comment the line above

  !> nan check (in fortran, for x=nan, x == x will return false)
  if (any(particle%x .ne. particle%x) .or. E .ne. E .or. particle%weight .ne. particle%weight) then
    !$omp critical
    write(*,*) 'ERROR: removing particle with nans in function single_self_interaction() (x,E,w,i_elm):', particle%x, E, particle%weight, particle%i_elm
    !$omp end critical
    particle%i_elm = 0 ! skip this one since sputtering went wrong
  end if

  ! store this particle's contribution to outgoing particle, heatflux and flux onto the wall
  diagnostics(i_wall_part_out) = diagnostics(i_wall_part_out) + 1
  diagnostics(i_wall_flux_out) = diagnostics(i_wall_flux_out) + particle%weight
  diagnostics(i_wall_heat_out) = diagnostics(i_wall_heat_out) + particle%weight * E *EL_CHG
  
end subroutine single_self_interaction


!> The potential drop from a debye sheath. Could support two-temperature model later
pure function debye_potential_drop(q, T_eV) result(U_drop)
  use constants, only: TWOPI, ATOMIC_MASS_UNIT, MASS_ELECTRON
  use phys_module, only: central_mass
  integer, intent(in) :: q
  real*8, intent(in) :: T_eV !< Local temperature in eV
  real*8 :: T_i, T_e, U_drop
  ! Equal temperatures
  T_i = T_eV
  T_e = T_eV

  !> Potential drop in eV
  U_drop = 0.5d0 * log((TWOPI * MASS_ELECTRON/(central_mass * ATOMIC_MASS_UNIT))*(1.d0+T_i/T_e))
end function debye_potential_drop


!> Calculate the energy gain of a potential drop from a sheath in the simplest model possible
pure function simple_potential_drop(q, T_eV) result(ion_energy)
  integer, intent(in) :: q
  real*8, intent(in) :: T_eV !< Local temperature in eV
  real*8 :: ion_energy !< The energy of the outgoing ion in eV
  
  ion_energy = 3.d0*real(q,8)*T_eV !< for sputtering from fluid perspective. Add the original energy E to this
end function simple_potential_drop


!> Integrate the sputtering yield over the distribution of incoming velocities.
pure function fluid_sputtering_yield(coeff, T_eV, Z, theta) result(yield)
  use gauss
  class(eckstein_coeff_set), intent(in) :: coeff
  real*8,                    intent(in) :: T_eV  !< Plasma temperature in eV
  integer,                   intent(in) :: Z     !< Atomic number of the incoming particles
  real*8,                    intent(in) :: theta !< angle of impact (usually assumed 0) in degrees
  real*8                                :: yield !< The sputter yield in atoms/ion

  real*8 :: U_drop
  integer :: q

  if (Z .le. 0) then
    q = 1
  else
    q = min(Z, 4) ! cap to 4 for divertor conditions
  end if 

  ! We use a simplified model for now! 2 T_i + 3 q T_e
  U_drop = simple_potential_drop(q, T_eV) ! assume particle has full charge
  yield = coeff%interp(2*T_eV + U_drop, theta)

  ! Alternative but unused version not assuming the simplified model:
  
  ! Since we use inverse transform sampling on u to calculate the energy we can
  ! just integrate over u from 0 to 1 to cover the whole distribution.
  ! Do this with n subelements, using gaussian quadrature in each element
  ! the subintervals go from 1/2 to n-1/2 to avoid using 0 and 1, since 1 should
  ! lead to infinity for sampling from a gaussian. Skipping the first part is reasonable
  ! since the sputtering yield will be very low there. For the high energies a maxwellian
  ! is perhaps not even the best approximation so that is probably not so bad either.
  
  ! The returned yield is averaged over the maxwellian at T_eV + the potential drop
  
  ! real*8 :: E
  ! integer :: i, j, k
  ! integer, parameter :: n_interval = 4 !< number of intervals to calculate. (using 1 is already pretty good)
  ! real*8, parameter :: idu = 1.d0/real(n_interval,8) !< interval size
  ! real*8 :: u(3) !< the integration point

  ! yield = 0.d0
  ! if (T_eV .le. 1d-1) return
  ! do i=0,n_interval*n_gauss-1
  !   u(1) = (real(i/n_gauss,8) + xgauss(mod(i,n_gauss)+1))*idu
  !   do j=0,n_interval*n_gauss-1
  !     u(2) = (real(j/n_gauss,8) + xgauss(mod(j,n_gauss)+1))*idu
  !     do k=0,1 ! we only use the sign of this one
  !       u(3) = real(k,8)

  !       ! add to this energy the plasma sheath potential
  !       call sample_fluid_particle_energy(T_eV, u, Z, E)
  !       U_drop = simple_potential_drop(q, T_eV) ! assume particle has full charge
  !       yield = yield + coeff%interp(E + U_drop, theta)
  !     end do
  !   end do
  ! end do
  ! yield = yield / (2*n_interval**2)
end function fluid_sputtering_yield


!> Calculate the flux to and some diagnostics for fluid flux in
!> a period delta_t.
!> 
!> Fluid_yield_integral contains a single scalar, the incoming 
!> fluid flux. Assume all particles are moving at the same
!> velocity, so multiplying with the relative density is enough 
!> to get the flux of a specific species.
!>
!> Assume that the impact angle of all particles is 0
subroutine project_sputter_vars_on_edge(this, sim)
  use mod_atomic_elements, only: atomic_weights
  use phys_module, only: central_mass, xpoint, xcase, min_sheath_angle, gamma
  
  type(wall_action),  intent(inout) :: this
  type(particle_sim), intent(in)    :: sim
  
  integer :: q, i, i_patch, Z
  real*8 :: vector_normal(3), cos_alpha, mass_ion, c_s, Gamma_d
  real*8 :: T_i, T_e, n_e, yield, vpar
  real*8, dimension(3) :: E, B, B_hat
  real*8 :: m, psi, U
  real*8 :: c_angle !< min_sheath_angle but then in radians, same as in mod_boundary_matrix_open

  real*8 :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_xpoint(2), psi_limit, R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
  integer :: i_elm_axis, ifail, i_elm_xpoint(2)

  c_angle = min_sheath_angle * PI/180.d0

  ! projection diagnostic
  ! Preparation (force my_id to 1 to suppress message)
  ! Note that this does not do proper time interpolation! We should probably
  ! have a proper function on the simulation to obtain those parameters
  ! for a rough estimate it will work however
  !t_xpoint = 0.d0
  !s_xpoint= 0.d0
  call find_axis(1,sim%fields%node_list,sim%fields%element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (xpoint) then
    call find_xpoint(1,sim%fields%node_list,sim%fields%element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_limit  = psi_xpoint(1)
    if((xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)))) then
      psi_limit = psi_xpoint(2)
    end if
  else
    if (sim%my_id .eq. 0) then
      write(*,*) "WARNING: limiter config for sputtering unsupported, use at your own risk"
    end if
    psi_limit = 0.d0 ! not really supported
  end if

  ! resetting fluid yield integral scalars
  do i=1,size(this%fluid_yield_integral%patch,1)
    this%fluid_yield_integral%patch(i)%scalars = -1
  end do

  do i_patch = 1, size(this%fluid_yield_integral%patch,1) !< different parts of edge domain
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp shared(this, sim, gamma, &
    !$omp i_patch, central_mass, psi_axis, psi_limit, c_angle) &
#endif
    !$omp private(i, n_e, T_e, vpar, E, B, psi, U, vector_normal, B_hat, cos_alpha, q, T_i, mass_ion, c_s, m, Gamma_d, &
    !$omp         yield, Z) schedule(static)
    do i = 1, size(this%fluid_yield_integral%patch(i_patch)%xyz, 2) !< over all nodes
      call sim%fields%calc_NeTevpar(sim%time, this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), this%fluid_yield_integral%patch(i_patch)%st(:,i), &
        real(this%fluid_yield_integral%patch(i_patch)%xyz(3,i), 8), n_e, T_e, vpar)
      
      call sim%fields%calc_EBpsiU(sim%time, this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), &
           this%fluid_yield_integral%patch(i_patch)%st(:,i), &
           real(this%fluid_yield_integral%patch(i_patch)%xyz(3,i), 8), &
           E, B, psi, U)
      
      !> normal vector calculation
      vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, &
          this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), &
          this%fluid_yield_integral%patch(i_patch)%st(1,i), &
          this%fluid_yield_integral%patch(i_patch)%st(2,i))
      
      !alpha = acos( dot_product(vector_normal,NORM2(B,dim=1))) !< acos is in radians
      ! the flux is given by the velocity along B dot n
      B_hat = B/norm2(B)
      cos_alpha = abs(dot_product(vector_normal,B_hat))
        
      q = 1 ! for calculation of sound speed
      T_i = T_e !< not made for model 400 [K]
      mass_ion = central_mass* ATOMIC_MASS_UNIT !< now we use only the deuterium soundspeed
      ! c_s = sqrt((k_boltz/mass_ion)*(T_e + gamma * T_i)) ! m/s !< gamma *(Te+Ti) in model303 and 307
      c_s = sqrt((k_boltz/mass_ion)*(gamma * (T_i+T_e))) !< IF model =303 / 307
      !<TODO: test c_s is vpar0, as this should account for all models
      
      Z = this%fluid_Z
      m = atomic_weights(Z) * ATOMIC_MASS_UNIT
      
      Gamma_d = n_e * abs(vpar) * norm2(B) * cos_alpha + n_e * c_s * c_angle

      ! Assume an impact angle of 0!
      ! need the abs here because we cheat using negative numbers to indicate D, T
      ! cap ionisation level to 4
      q = min(abs(this%fluid_Z), 4)
      select case(trim(this%type))
      case("wall recomb")
        yield = 1.d0 !<assuming complete wall saturation
      case("fluid sputter")
        yield = fluid_sputtering_yield(this%yield, T_e * K_BOLTZ/EL_CHG, q, 0.d0)
      case default
        call wrong_interaction_type(trim(this%type))
      end select

      this%fluid_yield_integral%patch(i_patch)%scalars(i,1) = Gamma_d * this%delta_t * yield !< particles / m^2 in this timestep

      if (this%do_wall_projection) then
        !associate (sc => this%wall_projection%patch(i_patch)%scalars) ! associate is nice to make more readable but cannot be used in OMP before version 4.5 (so not in OneAPI's OMP)
        ! These are all also multiplied by delta_t so we can make an average
        ! over the diagnostics period. Disregard the time in the units.
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+1) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+1) + n_e * this%delta_t ! n_e [m^-3]

        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+2) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+2) + T_e * K_BOLTZ / EL_CHG * this%delta_t ! T_e [eV]

        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+3) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+3) + cos_alpha * this%delta_t ! dimensionless, cosine of angle between wall normal and fieldline B

        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+4) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+4) + (psi - psi_axis)/(psi_limit - psi_axis) * this%delta_t ! normalized psi, dimensionless
        
        ! incoming fluid projections, should be similar to incoming particle projections, so this is like a sanity check
        ! number of particles incoming
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+5) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+5) + Gamma_d * this%delta_t !< incident particle flux (particles/m^2)

        ! incident energy integrated over delta_t
        ! where we assume the ion energy to be 2 k T_i + 3 q k T_e as in the ! sputtering calculation above
        ! J/m^2
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+6) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+6) + &
          Gamma_d * this%delta_t * (2.d0 * k_boltz * T_i + 3.d0 * k_boltz * q * T_e)

        ! sputtering yield in this time interval at this location [particles/m^2]
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+7) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+7) + &
          Gamma_d * this%delta_t * yield
        !end associate
      end if
    end do
    !$omp end parallel do
  end do
end subroutine project_sputter_vars_on_edge


!> Sample the energy of a particle with charge Z_ion in the plasma (before the sheath)
!> from the local temperature
!>
!> For the ions treated as a fluid we make the assumption that they travel at the
!> background plasma sound speed.
!>
!> The energy is determined by the criterion $<v_par> > c_s$ with $c_s$ the background
!> plasma sound speed. That leads to the factor sqrt(m_ion/central_mass)
!>
!> The calculation here proceeds as follows:
!> 1. Calculate total energy from chi_squared(3) distribution, equal to chi_squared(2) + chi_squared(1)
!> 2. Calculate ratio of perpendicular and total energies muB = perp/(perp + par) (see https://en.wikipedia.org/wiki/Chi-squared_distribution#Relation_to_other_distributions)
!> 3. Calculate new parallel energy from the square of E +- cs, with + or - 50/50
!> 4. Add all energies together
!> 5. Correct for atomic weight, assuming all velocities are central_mass velocities
pure subroutine sample_fluid_particle_energy(T_eV, u, Z_ion, E, E_threshold)
  use phys_module, only: central_mass
  use mod_sampling, only: sample_chi_squared_3
  use mod_atomic_elements, only: atomic_weights

  real*8, intent(in)             :: T_eV !< Temperature in eV
  real*8, intent(in)             :: u(3) !< random numbers for sampling
  integer, intent(in)            :: Z_ion
  real*8, intent(out)            :: E !< Energy in eV
  real*8, intent(in), optional   :: E_threshold !< Theshold energy in eV, not to sample particles below this energy

  real*8                         :: beta, v

  ! Sample an energy at the local temperature
  E = T_eV*0.5d0*sample_chi_squared_3(u(1)) ! in eV
  ! Solve now for u = 1-sqrt(1-x) (CDF of beta(1,1/2) distribution)
  beta = 2.d0*u(2)-u(2)**2
  ! this is also the ratio between perpendicular and total energies
  ! the parallel energy is then given by
  ! E*(1-beta)
  ! and we take the square root of that to get a parallel velocity
  ! the direction of this is either + or - with 50/50 probability.
  ! Add the soundspeed (positive) to this and calculate the new energy
  ! v = sqrt(2E/m) (+ or - with 50/50 prob)
  v = sign(sqrt(2.d0*E*EL_CHG*(1.d0-beta)/(central_mass*ATOMIC_MASS_UNIT)), u(3)-0.5d0) ! m/s
  ! the sound speed is sqrt(k (1+gamma) T/m) = sqrt(T_eV*EL_CHG/m)
  v = v + sqrt(T_eV*EL_CHG/(central_mass*ATOMIC_MASS_UNIT)) ! m/s
  E = E * beta + 0.5d0 * central_mass*ATOMIC_MASS_UNIT * v**2 / EL_CHG

  E = E*sqrt(atomic_weights(Z_ion)/central_mass) ! correct for atomic weight
end subroutine sample_fluid_particle_energy


!> find in which edge element patch index the particle is
function elm_in_patch(i_elm, edge_element_obj) result(i_patch)
  implicit none
  integer, intent(in) :: i_elm
  type(edge_elements), intent(in) :: edge_element_obj
  integer :: i_patch

  logical :: found

  found = .false.
  i_patch = -1

  do i_patch = 1,size(edge_element_obj%patch,1)
    ! if i_elm in the i_elm list of this edge domain exit the loop
    ! Note that this has issues at sharp corners, where particles may be
    ! lost in a different patch but at the same element number!
    if (any(i_elm .eq. edge_element_obj%patch(i_patch)%i_elm_jorek_edge(:))) then
      found = .true.
      exit
    endif
  end do
  ! i_patch should now be the first patch with correct element number, unless it wasn't found
  
  if (.not. found) then
    i_patch = -1 ! impossible number
    return
  end if

end function elm_in_patch


!> adds this particle's contribution to the wall_projection diagnostic tool
subroutine particle_projection_diagnostic(this, sim, particle, E, sputtering_yield)
  use phys_module, only: n_period, n_plane

  implicit none

  class(wall_action),              intent(inout) :: this
  class(particle_sim),             intent(in)    :: sim
  type(particle_kinetic_leapfrog), intent(in)    :: particle !< particle to undergo interaction
  real*8,                          intent(in)    :: E !< old energy of particle in eV
  real*8,                          intent(in)    :: sputtering_yield
  
  integer :: k, i_patch
 
  !> Prompt loss calculation
  integer :: is_prompt_loss
  real*8 :: Efield(3), B(3), pot, psi
  
  integer :: i_edge_elm, i_edge_nodes(4)
  real*8 :: area(4), dphi
  !> for mpi_reduce of particle contributions
  integer :: toroidal_offset !< Number of elements in the toroidal direction
  
  if(.not. this%do_wall_projection) return

  !> find in which patch the particle is lost
  i_patch = elm_in_patch(particle%i_elm, this%fluid_yield_integral)
  if (i_patch < 0) then
    write(*,*) "ERR in particle_self_reflection elm_in_patch, particle lost to somewhere unknown"
    return
  end if

  !> Write several diagnostics for the particle-particle sputtering
    ! the projection of a variable into the edge elements is simply a weighted addition to four points around an element
    ! Calculate the weight factors first and then store the relevant diagnostics
    ! find the corner point of the edge element we'll add the diagnostics to
  i_edge_elm = find_edge_element(this%wall_projection%patch(i_patch), particle%i_elm, particle%st(1), particle%st(2), particle%x(3))
  if (i_edge_elm .le. 0) then
    !$omp critical
    write(*,*) "ERROR: cannot find edge element for particle lost in this patch", particle%i_elm, particle%x(1), particle%x(2), i_edge_elm
    ! call flush(6)
    !$omp end critical
    return
  end if
  ! the weighting is done by inverse area
  ! 3-------|-----------4
  ! |   2   |   k=1     |
  ! |       |           |
  ! --------X------------
  ! |   4   |     3     |
  ! 1-------|-----------2
  ! in real space. i.e. calculate for each of the four quadrants above the surface area of the element
  ! and give them a fraction opposite area / total each.
  !
  ! The integrals are simple, since the elements are linear. It is given by
  ! \[
  !   \int_{l_0}^{l_1} \int_{\phi_0}^{\phi_1} R dl dphi
  ! \]
  ! The phi-integral drops out since it does not depend on l (they are orthogonal)
  ! and the other integral can be simplified since dl is along a straight line.
  ! this has as answer: 
  ! \[
  !   \left(r_0 l + \frac{1}{2} l^2 \frac{dr}{dl}\right) * (\phi_1 - \phi_0)
  ! \]
  ! with dr/dl = delta r / delta l (i.e. bounded between 0 and 1), 1 for purely outwards.
  !this%wall_projection%patch(i_patch)%scalars(index_node,5) = E * particle%weight
  toroidal_offset = this%wall_projection%patch(i_patch)%nsub_toroidal*n_plane
  if (toroidal_offset .eq. 1) toroidal_offset = 0 ! special case for fully axisymmetric
  i_edge_nodes = [i_edge_elm, i_edge_elm+1, &
      i_edge_elm + toroidal_offset,  &
      i_edge_elm + toroidal_offset + 1]

  ! area = r_0 l + (r_1-r_0) l / 2 = (r_1 + r_0) l / 2
  ! The indices k are as above shown, i.e. of the area opposite the node
  ! this is related to the edge nodes as
  ! 1 <-> 4 and 2 <-> 3, so 5-i
  do k=1,4
    if (i_edge_nodes(5-k) .gt. size(this%wall_projection%patch(i_patch)%xyz(1,:))) then
      write(*,*) "ERROR indexing problem in mod_wall_actioning",k,i_edge_elm,toroidal_offset, i_edge_nodes(5-k), size(this%wall_projection%patch(i_patch)%xyz(1,:))
      write(*,*) "ERROR temporary fix: set i_edge_nodes(5-k) = 1"
      i_edge_nodes(5-k) = 1
    end if

    area(k) = (this%wall_projection%patch(i_patch)%xyz(1,i_edge_nodes(5-k)) + particle%x(1)) &
       * norm2(this%wall_projection%patch(i_patch)%xyz(1:2,i_edge_nodes(5-k))-particle%x(1:2), dim=1) * 0.5d0
  end do
  ! multiply with delta-phi part
  ! we assume below that the particle is in this element (as it came from find_edge_element)
  dphi = TWOPI / (n_period * n_plane)
  area(1:2) = area(1:2) * modulo(dphi - particle%x(3), dphi) ! distance from X to top row
  area(3:4) = area(3:4) * modulo(particle%x(3) - dphi, dphi) ! distance from X to bottom row

  ! Multiply by this below (I might be guilty of some premature optimization here)
  is_prompt_loss = 0
  call sim%fields%calc_EBpsiU(sim%time, particle%i_elm, particle%st, particle%x(3), Efield, B, psi, pot)
  ! If the age of this particle is less than an a gyroperiod at the local magnetic field strength
  ! this particle is considered a prompt loss and will be written down below
  if ((sim%time - particle%t_birth) .lt. TWOPI * sim%groups(this%origin_group)%mass*ATOMIC_MASS_UNIT/(EL_CHG * norm2(B))) is_prompt_loss = 1

  ! we need to loop here since omp atomic cannot set an array at once
  !associate (sc => this%wall_projection%patch(i_patch)%scalars) ! associate is nice to make more readable but cannot be used in OMP before version 4.5 (so not in OneAPI's OMP)
  do k=1,4
    ! particle flux
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),1) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),1) + particle%weight * area(k)/sum(area)**2
    
    ! particle heat flux on edge elements (including sheath potential)
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),2) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),2) + particle%weight * E * EL_CHG * area(k)/sum(area)**2
    
    ! particle flux from prompt redeposition (i.e. from particles younger than 2 pi / omega_c)
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),3) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),3) + particle%weight * is_prompt_loss * area(k)/sum(area)**2
    
    ! sputtering yield
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),4) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),4) + particle%weight * sputtering_yield * area(k)/sum(area)**2
    
  end do
  !end associate
end subroutine particle_projection_diagnostic


subroutine write_wall_project_vtk(this, sim)
  use mpi_mod

  implicit none

  type(wall_action),  intent(inout) :: this
  type(particle_sim), intent(in)    :: sim
  
  integer :: nnos, i, ierr
  real*4, allocatable :: scalars(:,:) !< for mpi_reduce of particle contributions
  character(len=120)  :: filename

  ! if not initialised, setting initial value of this%last_diag_time
  if (this%last_diag_time < 0) then
    this%last_diag_time = sim%time - this%delta_t
  end if

  ! determine filename
  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if
  
  ! time normalising and MPI reducing the quantities
  do i = 1,size(this%wall_projection%patch,1)
    ! Turn all quantities from fluences into fluxes by dividing by the time since the last diagnostics output
    ! some of these (like T_e and n_e) were actually not fluences, but
    ! multiply those by this%delta_t anyway so this normalisation works and we get
    ! a decent time average
    this%wall_projection%patch(i)%scalars(:,:) = &
        this%wall_projection%patch(i)%scalars(:,:) / real(sim%time - this%last_diag_time,4)

    nnos = size(this%wall_projection%patch(i)%scalars,1)
    if (sim%my_id .eq. 0) then
      allocate(scalars(nnos,this%n_project_part))
    else
      allocate(scalars(0,0))
    end if
    
    ! Calculate the sum across mpi procs
    ! this needs to be done for all particle-quantities only (i.e. 1:this%n_project_part)
    call MPI_Reduce(this%wall_projection%patch(i)%scalars(:,1:this%n_project_part), &
        scalars, &
        nnos*this%n_project_part, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if (sim%my_id .eq. 0) then
      this%wall_projection%patch(i)%scalars(:,1:this%n_project_part) = scalars
    end if
    
    deallocate(scalars)
  end do 

  ! writing the projection
  if (sim%my_id .eq. 0) write(*,*) 'Writing wall projection diagnostics to ', trim(filename)
  call this%wall_projection%write_vtk_projection(filename)

  ! Reset diagnostic
  do i = 1,size(this%wall_projection%patch,1)
    this%wall_projection%patch(i)%scalars(:,:) = 0.d0
  end do
  this%i_step_diag = 0
  this%last_diag_time = sim%time
end subroutine write_wall_project_vtk


!> centralised routine to write out what input was wrong, and then stop the program
subroutine wrong_input(message, my_id, identifier)
  implicit none

  character(len=*), intent(in) :: message
  integer,          intent(in) :: my_id
  character(len=*), intent(in) :: identifier

  if(my_id == 0) then
    write(*,"(3A)") trim(message), " ", trim(identifier)
  end if

  stop
end subroutine wrong_input


!> centralised routine to throw the error of unsupported type and exit (please change this when you add a new type)
subroutine wrong_interaction_type(type, identifier)
  implicit none
  character(len=*), intent(in) :: type !< type which is not supported (will be trimmed in this subroutine)
  character(len=*), intent(in), optional :: identifier

  write(*,"(A)") "Wall interaction type "//trim(type)//" not supported (mod_wall_interaction.f90)"
  write(*,"(A)") 'Available types: "self sputter", "fluid sputter", "other sputter", "reflection" or "wall recomb" '
  if(present(identifier)) write(*,"(2A)") "Error detected ",trim(identifier)
  stop
end subroutine wrong_interaction_type


!> exit with message if origin_group is not target_group
subroutine check_self_type(this, my_id, identifier)
  implicit none
  type(wall_action), intent(in) :: this
  integer,           intent(in) :: my_id
  character(len=*),  intent(in) :: identifier
  
  character(len=1000) :: msg

  if(this%origin_group /= this%target_group) then
    write(msg,*) "type "//trim(this%type)//" is a self interaction type so origin_group should be target_group"
    call wrong_input(msg, my_id, identifier)
  end if
end subroutine


!> MPI reduces and writes the normal global diagnostics, 
!> returns the MPI reduced diagnostics back into diagnostics
subroutine write_global_diag(this,sim,diagnostics)
  use mpi_mod

  implicit none

  type(wall_action),   intent(in) :: this
  type(particle_sim),  intent(in) :: sim
  real*8, dimension(n_global_diagnostics), intent(inout) :: diagnostics !< diagnostics for the global wall loads
  
  real*8, dimension(n_global_diagnostics)                :: diagnostics_all_mpi !< MPI reduced diagnostics for the global wall loads
  integer :: ierr

  ! MPI reduce can be done at once for all diagnostics  
  call MPI_REDUCE(diagnostics, diagnostics_all_mpi, n_global_diagnostics, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  ! write standard diagnostics to logfile
  if (sim%my_id .eq. 0) then
    write(*,'(A,2f14.0)' ) "superparticles going (in/out) = ", diagnostics_all_mpi(i_wall_part_in),             diagnostics_all_mpi(i_wall_part_out) 
    write(*,'(A,2es16.6)') "particle flux (in/out) [#/s]  = ", diagnostics_all_mpi(i_wall_flux_in)/this%delta_t,diagnostics_all_mpi(i_wall_flux_out)/this%delta_t 
    write(*,'(A,2es16.6)') "heatflux (in/out) [W]         = ", diagnostics_all_mpi(i_wall_heat_in)/this%delta_t,diagnostics_all_mpi(i_wall_heat_out)/this%delta_t 
  endif

  ! write MPI reduced diagnostics back to diagnostics so other things can be printed if the user wants it
  diagnostics = diagnostics_all_mpi
end subroutine


!> first trims and then replaces spaces by underscores (_) in the string
function spaces2underscore(string_in) result(string_out)
  implicit none

  character(len=*), intent(in)  :: string_in !< string in which to replace spaces for _
  character(len=:), allocatable :: string_out

  character(len=:), allocatable :: string
  integer :: i

  string = trim(string_in)

  do i=1,len(string)
    if(string(i:i) == " ") string(i:i) = "_"
  end do

  string_out = string
end function

end module mod_particle_wall_interaction
