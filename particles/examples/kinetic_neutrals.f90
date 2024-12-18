!> Standard neutral atomic particle example in 2D
!>
!> The physics model includes puffing, recombination, ionisation, recycling and charge exchange for atomic Deuterium (no molecules).
!> OpenADAS is used for atomic physics
!> Plasma and neutral wall interaction are based on SDTRIM coefficients. 
!> External files y_DD.dat and ye_DD.dat are used to determine wall recombination of plasma into atomic neutral deuterium. 
!> These are based on interaction with a W wall
!>
!> To adjust the puff to your scenario, see  "! --- Setting up puffing" below
!>
!> To use a particle restart file: use restart_particles=.t. in the input file.
!>
!> To change to 3D recombination:
!> Change subroutine do1particlerecombination: use mod_integrate_recomb3D, only : integrate_recombination
!> Particle puffing is axisymmetric by default.

program kinetic_neutral_loop

use particle_tracer
use mod_particle_diagnostics
use mpi
use mod_interp
use mod_atomic_elements
use mod_particle_evolution
use mod_particle_recomb
use mod_particle_conservation
use mod_particle_io
use mod_event
use mod_project_particles
use mod_jorek_timestepping
use mod_random_seed
use mod_basisfunctions
use nodes_elements
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
use mod_projection_functions, only: proj_f_combined_density, proj_f_combined_energy, proj_f_combined_par_momentum
use mod_particle_puffing
use mod_edge_domain
use mod_edge_elements
use mod_atomic_coeff_deuterium, only: ad_deuterium 
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use mod_boundary,   only: boundary_from_grid
use equil_info

use phys_module, only: tstep,tstep_n,restart_particles, restart, t_start, nout
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
use phys_module, only: n_part_groups, n_aux_var
use phys_module, only: nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: deuterium_adas,sqrt_mu0_over_rho0
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: puff_rate, n_puff, valves
use phys_module, only: use_manual_random_seed, manual_seed

!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader
type(event), dimension(:), allocatable            :: sputter_events, puff_events ! can also be not allocatable and have size n_part_groups_max
type(particle_sputter)                            :: sputter_source
type(event)                                       :: gas_puff_event, gas_puff2_event
type(event), target                               :: project_jorek_feedback, jorek_stepper_event
type(pcg32_rng), dimension(:), allocatable        :: rng
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback
type(jorek_timestep_action), target               :: jorek_stepper
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: sputter_edge
type(particle_puffing)                            :: gas_puff, gas_puff2

real*8    :: rho_norm, t_norm, n_norm, tstep_fluid_si 
real*8    :: tstep_part_adj !< tstep_particles adjusted so that an integer amount of steps (nstep_particles) fit into a fluid step (tstep)
!$ real*8 :: w0, w1, mmm(3)

integer   :: n_reflect
integer   :: i, j, istep, group_num
integer   :: seed, i_rng, n_stream
integer   :: sputter_counter = 0
integer   :: recomb_counter  = 0
integer, dimension(:), allocatable :: recomb_groups

! Puffing parameters
real*8  :: t_puff_start          !< [s] time to start ramping the puff rate if puff_t_dependent=.true.
real*8  :: t_puff_slope          !< [s] time over which the puff rated is ramped to puff_rate (input) from t_puff_start where the puff rate was still puffing_rate_start
real*8  :: puffing_rate_start    !< [atoms/s] initial puff rate before the ramp if puff_t_dependent=.true.
real*8  :: poly_R(4)             !< [m] R coordinates of the quadrangular puffing valve if boxpuff=.true.
real*8  :: poly_Z(4)             !< [m] Z coordinates of the quadrangular puffing valve if boxpuff=.true.
real*8  :: poly_R2(4),poly_Z2(4) !  [m] second puffing valve location
logical :: puff_t_dependent      !< puff time dependent using a flat - ramp - flat pattern (=.true.) or no time dependence at all (.false.) 
logical :: boxpuff               !< whether to puff in a simple (=.false., uses r_valve etc) or quadrangular (=.true., uses pol_R, poly_Z) puff valve

!***********************************************************************
!*                            intialisation                            *
!***********************************************************************

! Start up MPI, jorek
call sim%initialize()

! Set up the field reader < can this be moved to sim%initialize
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)
tstep = tstep_n(1) !< the field reader overwrites tstep for some reason, this resets that

! setting up the particles
if (restart_particles) then
  ! reading the particles from a file
  if (sim%my_id == 0) write(*,*) 'INFO: READING PARTICLES RESTART FILE'
  partreader = event(read_action(filename='part_restart.h5'))
  call with(sim, partreader) !<defines sim%groups and the corresponding particles

  call configure_particle_group(sim)

  !TODO? Sven: We should make an option to use partreader but increase n_particles; may be similar to phi_zero_whrite to a sim_in and sim_out but with different allocation size.
else
  if (sim%my_id == 0) write(*,*) 'INFO: INITIALIZING PARTICLES', sim%n_cpu, " cpus "

  !> is this needed for neutrals?
  if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

  call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase )

  ! Setting up particle characteristics and allocation
  call configure_particle_group(sim)
  call allocate_particles(sim)
endif ! (restart_particles)

! Read Open ADAS data for plasma fluid
if (deuterium_adas .and. sim%groups(1)%use_kn_recombination) ad_deuterium =  read_adf11(sim%my_id,'96_h') !< move to core (jorek2_main for particles)

! --- Setting up random numbers for ionisation probability
seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
write(*,*) "id, n_cpu, n_stream",sim%my_id, sim%n_cpu, n_stream
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(1, seed, n_stream, i)
end do

! --- Check if the user tried to use nstep_particles rather than tstep_particles to define the particle timestepping
if (nstep_particles .ne. 0 .and. sim%my_id .eq. 0) then
  write(*,*) "ERROR: nstep_particles is defined in the input file, while for this example the combination tstep_particles, nstep and tstep define nstep_particles. Please remove nstep_particles from your input file to avoid ambiguity. Stopping now."
  stop
end if

!***********************************************************************
!*                         setting up physics                          *
!***********************************************************************

! --- Calculating normalisation constants
n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek 

! --- Setting up sputtering and recombination
allocate(sputter_events(n_part_groups), recomb_groups(n_part_groups)) 
do group_num=1, n_part_groups
  if (sim%groups(group_num)%use_kn_sputtering) then
    sputter_counter = sputter_counter + 1
    n_reflect = ceiling(sim%groups(group_num)%n_particles * sim%groups(group_num)%n_reflect_ratio)
    sputter_source = initialise_sputtering(sim%fields%node_list, sim%fields%element_list, group_num, n_reflect)
    sputter_events(sputter_counter) = event(sputter_source)
  endif

  if (sim%groups(group_num)%use_kn_recombination) then

    ! add group to the list of groups requiring recombination
    recomb_counter = recomb_counter + 1
    recomb_groups(recomb_counter) = group_num
  endif
enddo 

! --- Setting up recombination



!> Adapt the following to customize the time dependent puff rate:
!> puffing_rate_start = initial puffing rate [atoms/s]
!> puff_rate = final puffing rate [atoms/s] <input parameter>
!> t_puff_start = At what time the puffing rate starts to increase [s]
!> t_puff_slope = How much time it takes to increase linearly from puffing_rate_start to puff_rate [s]
!> n_puff = number of super particles puffed per valve per jorek timestep (should be small fraction of total number of super particles) <input parameter>
puff_t_dependent = .true. 
puffing_rate_start = puff_rate/1.5d0 !< initial puffing rate [atoms/s]

!> puff location can be determined for a circular valve by setting input parameters: 
!> r_valve (valve radius), R_valve_loc, Z_valve (R,Z, coordinates of simple valve)
!> if boxpuff=.true., poly_R(4),poly_Z(4) are the vertices of the quadrangular puffing valve
boxpuff = .true. !< whether to puff using 

if (sim%groups(1)%use_kn_puffing) then
  t_puff_start = 5000*t_norm !< start puffing after this amount of seconds, t_SI = t_jorek*t_norm jorek time units
  t_puff_slope = 4.d-3       !< [s] linearly ramps up the puffing during this time

  gas_puff = particle_puffing(n_puff, puff_rate/2.d0, valves(1), puff_t_dependent=puff_t_dependent,t_puff_start=t_puff_start,t_puff_slope=t_puff_slope, & 
      puffing_rate_start=puffing_rate_start/2.d0)
  gas_puff2 = particle_puffing(n_puff, puff_rate/2.d0, valves(1), puff_t_dependent=puff_t_dependent,t_puff_start=t_puff_start,t_puff_slope=t_puff_slope, &
      puffing_rate_start=puffing_rate_start/2.d0)
  
  gas_puff_event  = event(gas_puff)
  gas_puff2_event = event(gas_puff2)
  
  if (sim%my_id .eq.0) then
    write(*,*) "Gas puffing rate [#/s] : ", puff_rate
    write(*,*) "puff_t_dependent : ",puff_t_dependent, "with puff slope",t_puff_slope,"starting at", t_puff_start, "s"
  endif
else 
  gas_puff = particle_puffing(0, 5d20, valves(1))
  gas_puff2 = particle_puffing(0, 5d20, valves(1))
endif

! --- Set up feedback to the plasma (does not currently include recombination)
jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0=filter_par_n0,      &
                                filter = filter_perp, filter_hyper = filter_hyper, filter_parallel=filter_par, fractional_digits = 9, &
                                do_zonal = .false., calc_integrals=.false., to_vtk=.TRUE., to_h5 = .false., basename='projections', nsub=2)
aux_node_list => jorek_feedback%node_list

!> define feedback size dependent on the number of variables required for coupling
allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, n_aux_var))

jorek_feedback%rhs = 0.d0

! --- Setting up jorek timestepper
! For proper timestepping, the projections need to be defined before the jorek timestepper
jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

project_jorek_feedback = new_event_ptr(jorek_feedback,   start = sim%time)
jorek_stepper_event    = new_event_ptr(jorek_stepper,    start = sim%time)


!***********************************************************************
!*                           main loop                                 *
!***********************************************************************

istep = 0
do while (.not. sim%stop_now)
  istep = istep + 1
  if(sim%my_id .eq. 0) write(*,'(A100)'  ) "===================================================================================================="
  if(sim%my_id .eq. 0) write(*,'(A37,I6)') "Starting main loop iteration istep = ",istep
  if(sim%my_id .eq. 0) write(*,'(A100)'  ) "===================================================================================================="

  ! --- Determining the time stepping for this fluid step
  tstep = get_tstep_n(istep) ! tstep is also set in stepper, but tstep is already used in the calls before the stepper
  tstep_fluid_si = tstep*t_norm
  sim%time = sim%time + tstep_fluid_si ! carries the time at the end of the current step

  nstep_particles = ceiling(tstep_fluid_si / tstep_particles) ! ceiling makes sure tstep_part_adj is never bigger than tstep_particles
  tstep_part_adj = tstep_fluid_si / nstep_particles ! slightly smaller tstep_particles to fit an exact integer amount in one fluid timestep
  
  if (sim%my_id .eq. 0) then
     write(*,*) "PARTICLE : tstep_particles : ",tstep_particles
     write(*,*) "PARTICLE : tstep_part_adj  : ",tstep_part_adj
     write(*,*) "PARTICLE : sim%time        : ",sim%time
     write(*,*) "PARTICLE : nstep_particles : ",nstep_particles
     write(*,*) "PARTICLE : tstep_fluid_si  : ",tstep_fluid_si
     write(*,*) "PARTICLE : n*dt_part - dt  : ",nstep_particles*tstep_part_adj - tstep_fluid_si
  endif


  ! --- Interactions that happen on the fluid timestep (creating kinetic particles)
  
  !> The sputtering modules actually contains 3 different effects: sputtering (plasma to W, which is not used in this example), 
  !> kinetic particle reflection off the wall, and wall recombination of the plasma into kinetic neutrals (i.e. recycling)
  !> Do this call before recombination and puffing. Otherwise to-be-reflected particles can be overwritten.
  do i=1, sputter_counter    
    call write_to_outputfile(sim%my_id, "Sputtering")
    call with(sim, sputter_events(i))
  enddo
  
  do i=1, recomb_counter
    call write_to_outputfile(sim%my_id, "Recombination")
    call do_1particle_recombination(element_list,node_list, recomb_groups(i), jorek_stepper,rng, tstep_fluid_si) 
  enddo
    
  if (sim%groups(1)%use_kn_puffing) then
    call write_to_outputfile(sim%my_id, "Puffing")
    call with(sim, gas_puff_event) 
    call with(sim, gas_puff2_event)
  endif ! use_kn_puffing  


  ! --- Interactions that happen on the particle timesteps
  
  !> ionisation + CX + pushing the particles + calculating the feedback
  call write_to_outputfile(sim%my_id, "Particle loop")

  do group_num=1, n_part_groups
    call evolve_particle_group(sim, group_num, jorek_feedback, rng, tstep_part_adj)
  enddo  

  ! --- Update the fluid
  
  !> Project the collected feedback from the particles onto the finite element grid so that the MHD solver can use it
  !> Also writes the projection.vtk file which contains the interaction terms (particle, energy and momentum exchange to the fluid) and neutral density
  call write_to_outputfile(sim%my_id, "Projecting feedback from particles to fluid FE grid")
  call with(sim, project_jorek_feedback)
  
  !> Calls the MHD solver which timesteps the MHD fluid based on the fluid itself using the projected
  !> feedback of the particles as sources and sinks in the MHD equations 
  !> Also writes .h5 file, updates tstep and sets sim%stop_now = .true. if all fluid steps are done
  call write_to_outputfile(sim%my_id, "Fluid stepper")
  call with(sim, jorek_stepper_event) 
  

  ! -- Finalising the fluid timestep
  
  !Writing interim particle restart files every 500 fluid steps done. Overwrites previous restart file to save space
  if ( mod(istep,500) .eq. 0 ) then
    call write_to_outputfile(sim%my_id, "Writing interim_part_restart.h5")
    call write_simulation_hdf5(sim, 'interim_part_restart.h5')
  endif

  ! Writing some conservation checks to the ouput file
  call write_to_outputfile(sim%my_id, "Conservation checks")
  call conservation_checks(sim)

end do ! while

!***********************************************************************
!*                          end of simulation                          *
!***********************************************************************
call write_to_outputfile(sim%my_id, "End of simulation")
  
call write_simulation_hdf5(sim, 'part_restart.h5')

deallocate(sputter_events, recomb_groups)

call sim%finalize

!***********************************************************************
!*                          end of main program                        *
!***********************************************************************

contains

subroutine write_to_outputfile(id,what)
  implicit none
  
  integer, intent(in) :: id
  character(len=*),intent(in) :: what

  if(id .ne. 0) return

  write(*,'(A100)') "===================================================================================================="
  write(*,*) what
  write(*,'(A100)') "===================================================================================================="

end subroutine

function initialise_sputtering(node_list, element_list, target_group, n_reflect) result(sputter_source)

  use mod_edge_domain
  use mod_edge_elements

  type(type_node_list), intent(in)    :: node_list
  type(type_element_list)             :: element_list
  integer, intent(in)                 :: target_group
  integer, intent(in)                 :: n_reflect
  type(particle_sputter)              :: sputter_source
  !real*8, allocatable, dimension(:)   :: wall_albedo
  type(type_edge_domain), allocatable, dimension(:) :: edge_domains

  ! number of particles to sputter per species (should be renormalized to yield)

  call find_edge_domains(node_list,element_list, edge_domains)!, discont_corner=.true.)
  if (sim%my_id .eq. 0) write(*,*) "n_domains = ", size(edge_domains,1)
  
  call sputter_edge%prepare(node_list, element_list, edge_domains, nsub=6, nsub_toroidal=1)!,wall_albedo=wall_albedo)

  ! target group, number of particles per mpi task, densities, Zs, basename
  sputter_source = particle_sputter(sputter_edge, target_group, n_reflect)
  sputter_source%use_Yn_func = .false.
  sputter_source%n_save = nout !10 ! or nout
  sputter_source%albedo_for_neutrals = 1.d0
  sputter_source%sputtered_particle_weight_threshold = 1.d0

end function

pure function f_adapted(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:4)
  real*4 :: f

  coeff(0)=0.53
  coeff(1)=0.3
  coeff(2)=0.2
  coeff(3)=0.52
  coeff(4)=0.26

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  f = (f - coeff(4)) / (1.d0 - coeff(4))

end function f_adapted

pure function f_original(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_original

pure function f_toroidal_flux(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 :: s, psi_norm, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  s = 0.957 * psi_norm + 0.043 * psi_norm**2 

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_toroidal_flux


end program kinetic_neutral_loop