!> Testing the coupling of the projections of particles to JOREK

program recombination_loop

!use mod_integrate_recombination
use mod_H2_AMJUEL_coeff !compilen
use particle_tracer
use mod_particle_diagnostics
use mpi
use mod_interp
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_project_particles
!use mod_particle_loop
use mod_jorek_timestepping
use mod_random_seed
use mod_basisfunctions
use nodes_elements
use phys_module, only: tstep,restart_particles, restart, t_start, nout
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: use_ncs, use_pcs, use_ccs, deuterium_adas,sqrt_mu0_over_rho0
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
! use phys_module, only: use_sputtering , use_cx, use_ionisation, use_sputtering

use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG

use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
use mod_projection_functions, only: proj_f_combined_density, &
                                    proj_f_combined_energy, proj_f_combined_par_momentum
! use mod_radiation, only : proj_PLT
use mod_particle_puffing
use mod_edge_domain
use mod_edge_elements

 use mod_atomic_coeff_deuterium, only: ad_deuterium 

use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use mod_boundary,   only: boundary_from_grid
use equil_info

!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader
type(event)                                       :: D_sputter_event,gas_puff_event ,gas_puff2_event, gas_puff3_event!, partwriter
type(adf11_all)                                   :: adas
type(pcg32_rng), dimension(:), allocatable        :: rng
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback, project_atom_density,project_mol_density, project_current
type(jorek_timestep_action), target               :: jorek_stepper
type(particle_sputter)                            :: D_sputter_source
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge
type(particle_puffing)                            :: gas_puff
type(particle_puffing)                            :: gas_puff2,gas_puff3
type(write_particle_diagnostics)                  :: diag

real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
real*8    :: target_time, projection_time
real*8    :: physical_particles, weight
real*8    :: tstep_keep,oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm(3)
real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, rho_part, v2, tstart_jorek
real*8    :: momentum_conserv(3)
!$ real*8 :: w0, w1, mmm(3)

integer   :: n_particles_local,n_particles_local_mol, n_particles_local_mol_ion ,n_reflect,ifail
integer   :: i, j, k, l, m, n_steps, i_elm_old,ierr
integer   :: seed, i_rng, n_stream
integer   :: atoms, molecules

! Puffing parameters
real*8  :: r_valve, R_valve_loc, Z_valve,  R_valve_loc2, Z_valve2, puff_rate,t_puff_start,t_puff_slope, fueling_rate_start
real*8   ::r_valve3, R_valve_loc3, Z_valve3,puff_rate3,poly_R(4),poly_Z(4),poly_R2(4),poly_Z2(4),poly_R3(4),poly_Z3(4)
integer :: n_puff
logical :: puff_t_dependent,boxpuff


!use physics
logical :: use_recombination, use_puffing, use_cx, use_ionisation , use_sputtering,use_line_radiation, USE_H2plus_DISSOCIATION
logical :: use_molecules, use_dissociation, use_dissionisation, use_nondissionisation, use_mol_cx
logical  :: run_stepper, run_rec !, one_rec_only !< when recombination is used

! diagnostics
real*8    :: density_tot, density_in, density_out,  pressure, pressure_in, pressure_out
real*8    :: mom_par_tot, mom_par_in, mom_par_out, kin_par_tot, kin_par_out, kin_par_in
real*8    :: particles_remaining, momentum_remaining, energy_remaining, all_particles, all_momentum, all_energy, lost_particle_weights
real*8    :: totallostparticles, lost_energy, total_lost_energy, momentuminplasma(3), allmomentuminplasma(3), parallelmomentuminplasma, allparallelmomentuminplasma
real*8    :: all_momentum_conserv(3), all_totallostparticles, all_total_lost_energy
integer   :: superparticles_remaining,all_superparticles,closest_iteration 
!integer   :: particles_per_element


! Start up MPI, jorek
call sim%initialize(num_groups=2)

!> make sure tstep from namelist doesn't get overwritten
tstep_keep        = tstep
timesteps         = tstep_particles

!> saving part_restart every part_n_save steps
!part_i_save = 1
!part_n_save = 500

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)
tstep = tstep_keep !< fieldreader overwrites tstep, do this to counter that

! if (restart_particles) then 

   ! if (sim%my_id == 0) write(*,*) 'INFO: READING PARTICLES RESTART FILE'
   ! partreader = event(read_action(filename='part_restart.h5'))
   ! call with(sim, partreader) !<defines sim%groups en particles, if more than one group, change num_groups

   ! n_particles_local = size(sim%groups(atoms)%particles(:)) !< sputtering and puffing amount is function 
   ! write(*,*) "n_particles_local = ", n_particles_local
   ! !We should make an option to use partreader but increase n_particles
   ! !may be similar to phi_zero_whrite to a sim_in and sim_out but with different allocation size.
      
! else
	if (sim%my_id == 0) write(*,*) 'INFO: INITIALIZING PARTICLES', sim%n_cpu, " cpus "

	! Read Open ADAS data
	write(*,*) "deuterium_adas (12)",  deuterium_adas
	adas = read_adf11('12_h')

	!> is this needed for neutrals?
	if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
	call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

	call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase )
	!<
	
	!> may be write it as atom_group = 1 and molecule_group = 2
	!> sim%groups(atom_group)%Z , makes it more readable??
  ! Set up particle group characteristics
  atoms=1
	sim%groups(atoms)%Z    = -2 !< for deuterium 1
	sim%groups(atoms)%mass = atomic_weights(-2) !< atomic mass units
	sim%groups(atoms)%ad   = adas
  if (sim%my_id == 0) write(*,*) 'Group 1 = D'
  
  !Set diagnostic for particles exiting the plasma
  
  write(*,*) "test125 size(sim%groups(atoms)%particles)", size(sim%groups(atoms)%particles)

  ! Set up hydrogenic molecule group characterisitcs
  molecules=2
	sim%groups(molecules)%Z    = -2 !< for deuterium 1
	sim%groups(molecules)%mass = 2.d0*atomic_weights(-2) !< atomic mass units
	! sim%groups(molecules)%AMJUEL   = adas
  if (sim%my_id == 0) write(*,*) 'Group 2 = D2'
  

	
	! setting up particles per MPI node
	n_particles_local = int(n_particles/sim%n_cpu) 
  allocate(particle_kinetic_leapfrog::sim%groups(atoms)%particles(n_particles_local))
  n_particles_local_mol = int(100/sim%n_cpu) !< TODO: make new n_particles for molecules or make n_particles(array)
  !n_particles_local_mol = int(10000/sim%n_cpu) !< TODO: make new n_particles for molecules or make n_particles(array)
  allocate(particle_kinetic_leapfrog::sim%groups(molecules)%particles(n_particles_local_mol))


	!>initialise particles here if needed
	! call initialise_particles_H_mu_psi 
	! call adjust weights

  ! Set up particles
  ! sim%groups(molecules)%Z    = -2
  ! sim%groups(molecules)%mass = atomic_weights(-2) !< atomic mass units
  ! sim%groups(molecules)%ad   = adas

  !allocate(particle_kinetic_leapfrog::sim%groups(molecules)%particles(n_particles_local))

  call initialise_particles_H_mu_psi(sim%groups(molecules)%particles, sim%fields, pcg32_rng(), sim%groups(molecules)%mass, &
            uniform_space=.true., uniform_space_rej_f=f_psi_inside, uniform_space_rej_vars=[1], charge = 0)

  physical_particles = 1.d17
  !physical_particles = 0.d0
  weight = physical_particles/(n_particles_local_mol * sim%n_cpu) !n_cpu is nu 1
  !weight = 0.d0
  v_kin_temp = 0.d0! sqrt( (2.d0 * 1d5) / (sim%groups(molecules)%mass* ATOMIC_MASS_UNIT) / physical_particles) !kinetische temp. Kan met 0 geinitalisserd worden



  select type (p => sim%groups(molecules)%particles)
  type is (particle_kinetic_leapfrog)

    p(:)%q      = 0
    p(:)%weight = weight
    write(*,*) 'test123 184 p-> simgroupsmoleculesparticles'
    do j=1,size(p,1)
      call sim%fields%calc_EBpsiU(sim%time , p(j)%i_elm, p(j)%st, p(j)%x(3), E, B, psi, U)
      B_norm = B/norm2(B)
      p(j)%v(1)  = v_kin_temp * B_norm(1)
      p(j)%v(2)  = v_kin_temp * B_norm(2)
      p(j)%v(3)  = v_kin_temp * B_norm(3)
!todo: 0 maken
  !    p(j)%weight = p(j)%weight* (1.d0 + 0.8d0*cos(    p(j)%x(3))  +  0.d0*sin(     p(j)%x(3)) &
  !                                     + 11.d0*cos(2.d0*p(j)%x(3)) + 13.d0*sin(2.d0*p(j)%x(3)) &
  !                                     + 17.d0*cos(3.d0*p(j)%x(3)) + 19.d0*sin(3.d0*p(j)%x(3)) &
  !                                     + 23.d0*cos(4.d0*p(j)%x(3)) + 27.d0*sin(4.d0*p(j)%x(3))  )
                                    
    end do

    call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(molecules)%mass, sim%fields, sim%time, timesteps)

  end select



	select type (p => sim%groups(atoms)%particles)
	type is (particle_kinetic_leapfrog)  
	!> only set everything to zero when you do not initialise particles
	p(:)%q      = 0 !< for neutrals
	p(:)%weight = 0.0!weight
	p(:)%i_elm  = 0
	p(:)%v(1)   = 0.d0 
	p(:)%v(2)   = 0.d0
	p(:)%v(3)   = 0.d0
	! call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(atoms)%mass, sim%fields, sim%time, timesteps)
  end select

  ! 	select type (p => sim%groups(molecules)%particles)
	! type is (particle_kinetic_leapfrog)  
	! !> only set everything to zero when you do not initialise particles
	! p(:)%q      = 0 !< for neutrals
	! p(:)%weight = 0.0!weight
	! p(:)%i_elm  = 0
	! p(:)%v(1)   = 0.d0 
	! p(:)%v(2)   = 0.d0
	! p(:)%v(3)   = 0.d0
	! ! call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(atoms)%mass, sim%fields, sim%time, timesteps)
  ! end select
  
! endif ! (restart_particles)

! setting up particles per MPI node and timestep
rho_part    = 1.195d19 !(corrected value to obtain density=1.441e17 (as in benchmark, for original profile with toroidal flux) 
! n_particles_local = int(n_particles/sim%n_cpu) 
! timesteps         = tstep_particles
! tstep_keep        = tstep

! selecting physics (should be done in input file) !working: cx, ionisation, line_rad
use_puffing           = .false. !.false. 
use_cx                = .false. !.true.
use_ionisation        = .false. !.false.!.false.
use_sputtering        = .false. !.false. !false
use_recombination     = .false.  !
use_line_radiation    = .false.
use_molecules         = .true.
use_dissociation      = .false.
use_nondissionisation = .false.
use_dissionisation    = .false. 
use_mol_cx            = .false.
USE_H2plus_DISSOCIATION   = .true.


! Read Open ADAS data for plasma fluid
 if (deuterium_adas .and. use_recombination) ad_deuterium =  read_adf11('96_h') !< move to core (jorek2_main for particles)
 
n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek 
 
! Setting up edge_elements and amount of sputtered super particles per event
!TODO: set sputtering up for molecules as well.
if (use_sputtering) then  
  n_reflect = int(n_particles_local* sim%n_cpu * 5.d-4) !1.d-3 int(n_particles_local * 2.d-3)
  D_sputter_source = initialise_sputtering(sim%fields%node_list, sim%fields%element_list, n_reflect)
  D_sputter_event = event(D_sputter_source)
endif

! setting up particle puffing Top puff
! puff_t_dependent = .true. !.true. !< select if you want time dependent puffing
! puff_rate = 40.d21 !40.d21!100.d21 !8.85d21 !4.d21 !8.d22 !4.d22 !4.d21
! fueling_rate_start = 10.d21 !10 40 worked, 20 before
! r_valve     = 0.02d0!              0.01d0 !0.04d0 !0.02d0 !0.04d0 !.005d0
! R_valve_loc = 4.27d0!               4.4d0 !4.42787 !4.42787!2.33!2.6!2.1 !< for JET test !1.98991!2.58888  or 1.98991
! Z_valve     = -3.74d0!             -3.8d0 !-3.7 !-3.77948! -1.86 !-1.0!-1.75 !-0.550736!1.86579   or -0.550736

! R_valve_loc2 = 5.55d0!                  5.4d0 !5.46d0
! Z_valve2     = -4.35d0!                  -4.19d0 !-4.2d0

! puff_rate3 = 160.d21 !136.d21 ! 109.d21 !72.d21 !160.d21 !160.d21!85.d21
! R_valve_loc3 = 6.05d0!                  5.4d0 !5.46d0
! Z_valve3     = 4.15d0! 
! r_valve3    = 0.10d0!  .12

!Bot puff
puff_t_dependent = .true. !.true. !< select if you want time dependent puffing
puff_rate = 70.d21 !70.d21 !280.d21 !160.d21 !40.d21!100.d21 !8.85d21 !4.d21 !8.d22 !4.d22 !4.d21
fueling_rate_start = 40.d21 !40.d21 !40 40 worked, 20 before
r_valve     = 0.05d0 !0.02d0!              0.01d0 !0.04d0 !0.02d0 !0.04d0 !.005d0
R_valve_loc = 4.3d0 !4.27d0!               4.4d0 !4.42787 !4.42787!2.33!2.6!2.1 !< for JET test !1.98991!2.58888  or 1.98991
Z_valve     = -3.8d0 !-3.74d0!             -3.8d0 !-3.7 !-3.77948! -1.86 !-1.0!-1.75 !-0.550736!1.86579   or -0.550736
poly_R = (/4.2566d0 ,4.474d0 ,4.237d0 ,4.4917d0 /)
poly_Z= (/-3.727d0 ,-3.629d0 ,-3.7738d0 ,-3.6587d0 /)


R_valve_loc2 = 5.5d0 !5.55d0!                  5.4d0 !5.46d0
Z_valve2     = -4.35d0!                  -4.19d0 !-4.2d0
poly_R2 = (/5.426d0 ,5.559d0 ,5.455d0 ,5.5586d0 /)
poly_Z2= (/-4.155d0 ,-4.4005d0 ,-4.0803d0 ,-4.330d0 /)


puff_rate3 = 0.d21 !136.d21 ! 109.d21 !72.d21 !160.d21 !160.d21!85.d21
R_valve_loc3 = 6.05d0!                  5.4d0 !5.46d0
Z_valve3     = 4.15d0! 
r_valve3    = 0.10d0!  .12
poly_R3 = (/5.77d0 ,6.735d0 ,5.72d0 ,6.68d0 /)
poly_Z3= (/4.51d0 ,3.760d0 ,4.46d0 ,3.71d0 /)
boxpuff = .true.


!R_valve_loc = 4.307! touching leg
!Z_valve     = -3.7898!
if (use_puffing) then  
  n_puff      = int(5.d-5*n_particles_local* sim%n_cpu) !0.25 0.5d-4 !< now total n_puff
  if (puff_t_dependent) then
	t_puff_start = 5000*t_norm !25000*t_norm !34995*t_norm !5000*t_norm !< start puffing after this amount of seconds, t_SI = t_jorek*t_norm jorek time units
	t_puff_slope = 4.d-3 !4.d-3 !< linearly ramps up the puffing during this time
	!fueling_rate_start = 5.d21 !40 worked, 20 before
	!gas_puff = particle_puffing(n_puff, puff_rate/2.d0, r_valve, R_valve_loc, Z_valve, puff_t_dependent, t_puff_start, t_puff_slope)
	!gas_puff2 = particle_puffing(n_puff, puff_rate/2.d0, r_valve, R_valve_loc2, Z_valve2, puff_t_dependent, t_puff_start, t_puff_slope)
    gas_puff = particle_puffing(n_puff, puff_rate/2.d0, r_valve, R_valve_loc, Z_valve,target_group=molecules, puff_t_dependent=puff_t_dependent,t_puff_start=t_puff_start,t_puff_slope=t_puff_slope, & 
	           fueling_rate_start=fueling_rate_start/2.d0,poly_R=poly_R,poly_Z=poly_Z,boxpuff=boxpuff)
	gas_puff2 = particle_puffing(n_puff, puff_rate/2.d0, r_valve, R_valve_loc2, Z_valve2,target_group=molecules, puff_t_dependent=puff_t_dependent,t_puff_start=t_puff_start,t_puff_slope=t_puff_slope, &
	            fueling_rate_start=fueling_rate_start/2.d0,poly_R=poly_R2,poly_Z=poly_Z2,boxpuff=boxpuff)
	! gas_puff3 = particle_puffing(n_puff/2, puff_rate3   , r_valve3, R_valve_loc3, Z_valve3,target_group=2, puff_t_dependent=puff_t_dependent,t_puff_start=t_puff_start,t_puff_slope=t_puff_slope, &
				! fueling_rate_start=0.d21,poly_R=poly_R3,poly_Z=poly_Z3,boxpuff=boxpuff) !20.d21
  else

	gas_puff = particle_puffing(n_puff, puff_rate, r_valve, R_valve_loc, Z_valve,target_group=molecules) ! was 1
	!gas_puff2 = particle_puffing(n_puff, puff_rate/2.d0, r_valve, R_valve_loc2, Z_valve2,target_group=2)!-0.0) !-1.77 ! jet 2.8d0, -1.77
	! gas_puff3 = particle_puffing(n_puff, 20.d21, 0.12d0, 6.05, 4.15)
  end if
  gas_puff_event = new_event_ptr(gas_puff) !< new_event_ptr  allows for changing the gas puff within the even structure without redefining the event
  !gas_puff2_event = new_event_ptr(gas_puff2)
  ! gas_puff3_event = new_event_ptr(gas_puff3)
	!gas_puff = particle_puffing(n_puff, 5d22, r_valve, R_valve_loc, Z_valve)
	
	if (sim%my_id .eq.0) then
	write(*,*) "Gas puffing rate [#/s] : ", puff_rate
	write(*,*) "puff_t_dependent : ",puff_t_dependent, "with puff slope",t_puff_slope,"starting at", t_puff_start, "s"
	endif
else 
	n_puff = 0.d0
	gas_puff = particle_puffing(n_puff, 5d20, r_valve, R_valve_loc, Z_valve)
	gas_puff2 = particle_puffing(n_puff, 5d20, r_valve, R_valve_loc, Z_valve)
	! gas_puff3 = particle_puffing(n_puff, 5d20, r_valve, R_valve_loc, Z_valve)
endif

! write(*,*) 'main : t_start = ',t_start

! n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
! rho_norm  = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
! t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek

tstep_si  = tstep * t_norm
n_steps   = floor(tstep_si / timesteps)
timesteps = tstep_si / n_steps
n_steps   = tstep_si / timesteps

if (sim%my_id .eq.0) then
  write(*,*) 'main : t_start = ',t_start
  write(*,*) ' adapt time step to be multiple of jorek time step'
  write(*,*) "tstep = ", tstep_si, n_steps, timesteps
  write(*,*) "check :", n_steps, tstep_si - n_steps*timesteps
endif
!< sim%time = t_start ? 

!partwriter = event(write_action()) !< event writing particle restart files
! Set up feedback
jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0=filter_par_n0,      &
                                filter = filter_perp, filter_hyper = filter_hyper, filter_parallel=filter_par, fractional_digits = 9, &
                                do_zonal = .false., calc_integrals=.false., to_vtk=.TRUE., to_h5 = .false., basename='projections', nsub=2)

! jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     ! filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
                     ! filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                     ! fractional_digits = 9,  to_vtk=.TRUE., to_h5 = .FALSE., basename='projections')
aux_node_list => jorek_feedback%node_list

!> define feedback size as function of the coupling scheme
if (use_ncs) then
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 4)) !< stacksize should be big enough
elseif (use_pcs) then  ! not implemented yet!
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))
elseif (use_ccs) then  ! not implemented yet!
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 4))
else
  stop 'define use_ncs, use_pcs or use_ccs'
endif
jorek_feedback%rhs = 0.d0

!Setting up projections 
!project_density = new_projection(sim%fields%node_list, sim%fields%element_list,   &
!                      filter = 0d-3, filter_hyper = 1d-5, filter_parallel = 0.d0, &
!                      f=[proj_f(proj_one, group = 1)], fractional_digits = 9,     &
!                      calc_integrals=.true., to_vtk=.true., to_h5=.false., basename='density', nsub=5)

project_atom_density = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
                     filter_n0 = 5.d-6, filter_hyper_n0 = 2.d-11, filter_parallel_n0 = filter_par_n0, &
                     f=[proj_f(proj_one, group = 1)], &
                     fractional_digits = 9,  to_vtk=.TRUE., to_h5=.FALSE., basename='atom_density', nsub=2)

call with(sim, project_atom_density)

project_mol_density = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
                     filter_n0 = 5.d-6, filter_hyper_n0 = 2.d-11, filter_parallel_n0 = filter_par_n0, &
                     f=[proj_f(proj_one, group = 2)], &
                     fractional_digits = 9,  to_vtk=.TRUE., to_h5=.FALSE., basename='mol_density', nsub=2)
					 
call with(sim, project_mol_density)					 

! if (use_line_radiation) then
	! project_PLT = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     ! filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
                     ! filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                     ! f=[proj_f(proj_PLT, group = 1)], &
                     ! fractional_digits = 9,calc_integrals=.true.,  to_vtk=.TRUE., to_h5=.FALSE., basename='linerad', nsub=5)
! endif	 !use_line_radiation				 
!project_current = new_projection(sim%fields%node_list, sim%fields%element_list,   &
!                      filter = 0d-3, filter_hyper = 1d-5, filter_parallel = 0.d0, &
!                      f=[proj_f(proj_jPhi, group = 1)], fractional_digits = 9,    &
!                      calc_integrals=.true., to_vtk=.false., to_h5=.false., basename='current', nsub=5)

! For proper timestepping, the projections need to be defined before the jorek timestepper
jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

diag = write_particle_diagnostics(filename='diag.h5', append=.true.)

if (restart_particles) then
   tstart_jorek = sim%time + tstep_si
else
   tstart_jorek = sim%time
endif

if (sim%my_id .eq. 0) write(*,*) 'tstart_jorek : ',tstart_jorek

diag_time = 1.1*n_steps*timesteps
events = [ new_event_ptr(jorek_feedback,   start = tstart_jorek),            &
           new_event_ptr(jorek_stepper,    start = tstart_jorek),            & 
!		   event(D_sputter_source      ,start = tstart_jorek+tstep_si/2.d0, step=tstep_si),&
!		   event(gas_puff, step = 5.d-6),                                &
!		   event(gas_puff2, step = 5.d-6),                                & !
!          new_event_ptr(D_sputter_source, start = tstart_jorek, step=diag_time), & !20.d-7 1 sputter per jorek timestep step=1.d-6), &
!          event(count_action(),           start = tstart_jorek, step=1d-5), &
!          event(write_particle_diagnostics(filename='diag.h5'), step=diag_time), &
!         event(write_action(), step=diag_time),                        &
!           new_event_ptr(project_density, step=0.11d-6),                  &
!		   new_event_ptr(project_density, step=2.5d-6),                  &
!          event(diag, step=diag_time)									 &
           event(stop_action(), start=1d12)                              &
        ]


!< physics and projection like gas_puff, sputtering and projection density can be used in the events[] list as well.
!< We've decided to call them separately as this then doesn't interrupt the timestepping and may be altering the result.  ~Sven
!< look at the "if (run_stepper)then " part to see how to call events or projections as function of the jorek timestepper 
!< without interupting the simulation


! if(.not. restart_particles) then
jorek_stepper%extra_event => events(1) !< is used as first event before enetering particle loop (skipped if particles_restart)
! endif !restart_particles
!================================================================================================
!                                      MAIN PARTICLE LOOP
!================================================================================================

totallostparticles = 0.d0
total_lost_energy = 0.d0


! Set up random numbers for ionisation probability
seed = 1 !random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(1, seed, n_stream, i)
end do
!write(*,*) "test123 471, stream"
! Call events at sim%time once to help event scheduler, before entering particle loop
step_rest_time = 0.d0
call with(sim, events, at=sim%time)

! do recombination to initialise particles. Maintains conservation with first jorek timestep
if (.not. restart_particles) then 
	if (sim%my_id .eq. 0) write(*,*) "Do 1 particle recombination"
	particle_step_time = 0.d0 !< results in first step Nans for recombination diagnostics
	if (use_recombination) call do_1particle_recombination(element_list,node_list,jorek_stepper,rng,particle_step_time) 

	call with(sim, project_atom_density) !< directly project the first recombination at t=0
endif !.not. restart_particles

do while (.not. sim%stop_now) !begin loop

  target_time = next_event_at(sim, events) 
  particle_start_time = (sim%time - step_rest_time)
  particle_step_time  = target_time - particle_start_time
  n_steps             = particle_step_time/timesteps
  step_rest_time      = particle_step_time - real(n_steps,8) * timesteps
  write(*,*) 'n_steps = ', n_steps
  if (sim%my_id .eq. 0) then
     if (n_steps < 10) write(*,*) 'low n_steps,', n_steps
     write(*,*) 'Time difference between particles and jorek: ', step_rest_time
     write(*,*) "PARTICLE : target time         : ",target_time
     write(*,*) "PARTICLE : timesteps           : ",timesteps
     write(*,*) "PARTICLE : sim%time            : ",sim%time
     write(*,*) "PARTICLE : particle_start_time : ",particle_start_time
     write(*,*) "PARTICLE : particle_step_time  : ",particle_step_time
     write(*,*) "PARTICLE : n_steps             : ",n_steps
     write(*,*) "PARTICLE : step_rest_time      : ",step_rest_time
  endif
   
  !------------ Call Kinetic loops --------------------
  
  call with(sim, counter) !count particles in all groups and write diagnostic

  ! in this area we should call all the kinetic loops.
  ! we should be carefull to not reset the rhs to zero in every kinetic loop.
  !> may be only set it to 0.d0 here??
  jorek_feedback%rhs = 0.d0
  
  
  !> ionisation + CX + pushing the particles + calculating the feedback
  call loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)
  
  all_totallostparticles = 0.d0
  all_total_lost_energy = 0.d0
  all_momentum_conserv = 0.d0
  lost_particle_weights = 0.d0 
  lost_energy = 0.d0
  momentum_conserv = 0.d0
  select type (particles => sim%groups(atoms)%particles)
  type is (particle_kinetic_leapfrog)
    do j=1,size(particles,1) 
      if (particles(j)%i_elm .le. 0) then !sum weights of lost particles
        !write(*,*) 'test21 i_elm', particles(j)%i_elm, 'weight =', particles(j)%weight
        lost_particle_weights = lost_particle_weights + particles(j)%weight

        lost_energy    = lost_energy    + particles(j)%weight * dot_product(particles(j)%v,particles(j)%v) * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT /2.d0
        
        
        momentum_conserv(1) = momentum_conserv(1) + particles(j)%weight * (cos(particles(j)%x(3))*particles(j)%v(1)-sin(particles(j)%x(3)) * particles(j)%v(3) )* sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
        momentum_conserv(2) = momentum_conserv(2) + particles(j)%weight * (-sin(particles(j)%x(3))*particles(j)%v(1)-cos(particles(j)%x(3))  * particles(j)%v(3) ) *sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
        momentum_conserv(3) = momentum_conserv(3) + particles(j)%weight * 1.d0 *sim%groups(atoms)%mass* particles(j)%v(2) * ATOMIC_MASS_UNIT
  
        !if (particles(j)%v(1) .ne. 0.d0) then
        !  write(*,*) 'test1256', sqrt(abs(dot_product(particles(j)%v,particles(j)%v)))
        !endif
        !lost_energy    = lost_energy    + particles(j)%weight * 4.5d0 * EL_CHG !dot_product(particles(j)%v,particles(j)%v) * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT /2.d0

        particles(j)%weight = 0.d0 !test to see if particles get overwritten
        particles(j)%v = 0.d0

        

          
      
        
      
      

      endif
      !if (particles(j)%i_elm .gt. 0) write(*,*) 'j - weight - ielm:', j, particles(j)%weight, particles(j)%i_elm
  
    enddo
  end select
  totallostparticles = totallostparticles + lost_particle_weights
  total_lost_energy = total_lost_energy + lost_energy
  call MPI_REDUCE(totallostparticles, all_totallostparticles, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(total_lost_energy, all_total_lost_energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_conserv, all_momentum_conserv, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  !write(*,*) 'time',sim%time,'test213 Lost particle weights', lost_particle_weights, 'Sum lost particles: ', totallostparticles
  if (sim%my_id .eq. 0) write(*,'(a,1e16.8,a,1e16.8,a,1e16.8,a,3e16.8)') 'Time total diag ', sim%time, ' Total lost particles ', all_totallostparticles, ' Total lost energy ', all_total_lost_energy, ' Momentum conserv ', all_momentum_conserv
  !write(*,'(a,1e16.8,a,1e16.8)') 'Time', sim%time, 'Total lost particles', totallostparticles
  !write(*,'(a,1e16.8,a,1e16.8)') 'Time', sim%time, 'Total lost energy', total_lost_energy
  !write(*,'(a,1e16.8,a,3e16.8)') 'Time', sim%time, 'Momentum conserv', momentum_conserv
  if (use_molecules) then
    !write(*,*) 'test123 we are in the if molecules statement'
    call loop_particle_kinetic_molecule_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)
  endif !use_molecules



  !? add molecule loop here. !
  !call loop_particle_kinetic_molecule_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)
  !< Does the order matter here? First atoms than molecules? or first molecules than atoms??
  
  !----------------------------------------------------
  
  !> do particle sources if next event is jorek_stepper (event(2))
  !> run_stepper tells if particle sources will run. run_stepper must be before call with(sim..
  run_stepper = events(2)%run_at(target_time) !
  projection_time = target_time !< = target_time, but for the purpose of projections. needed for call to project diagnostics such as project_density
  sim%time = target_time !< set time to exactly target_time for calling next event
  
  call with(sim, events, at=sim%time) !< gives new target time
  !doet feedback en nieuwe tijdstap

  !> run particle source routines directly after the jorek_stepper
  !> Density projection added which now run every nout steps
  !> You can put anything in here that you want to solely depend on the jorek timestep.
  if (run_stepper)then

	!> call projection only every nout jorek steps. useful for longer runs
	closest_iteration = nint((projection_time - tstart_jorek)/(tstep_si*nout)) !< very similar to run_at function. May be put this in a function?
	if ( (abs((tstart_jorek +closest_iteration*tstep_si*nout) -projection_time) .le. 1.d-13) .or. sim%stop_now) then !< == true every tstep * nout steps
		call with(sim, project_atom_density)
		call with(sim, project_mol_density)
		
		 ! if (use_line_radiation) call with(sim, project_PLT)
	endif !< write projection or diagnostics
	
	if ( (abs((tstart_jorek +nint((projection_time - tstart_jorek)/(tstep_si*500))*tstep_si*500) -projection_time) .le. 1.d-13)) then !< == true every tstep * 100 steps
		call write_simulation_hdf5(sim, 'interum_part_restart.h5') !trim(this%get_filename(sim%time)))
		

	endif !< write interim particle restart file every 100 tsteps
	!if (part_i_save .ge. part_n_save) then
	!	call with(sim, partwriter)
	!	part_i_save = 0
	!endif
	!part_i_save = part_i_save + 1
	
	if (use_recombination) then
	  !call recombination
	  call do_1particle_recombination(element_list,node_list,jorek_stepper,rng,particle_step_time) 
    endif !use_recombination
  
    


	if (use_sputtering) then
	  ! call sputtering
	   call with(sim, D_sputter_event) !event(D_sputter_source))
	endif   !use_sputtering
	  
	if (use_puffing) then
      call with(sim, gas_puff_event) 
	  !call with(sim, gas_puff2_event)
	  ! call with(sim, gas_puff3_event)
    endif ! use_puffing	
  endif !run_stepper
  


!>  separate subroutine?
!=======================================================================
!                     Run diagnostics for conservation check for atoms
!======================================================================== 
  call Integrals_3D(sim%my_id, sim%fields%node_list, sim%fields%element_list, density_tot, density_in, density_out, &
                    pressure, pressure_in, pressure_out, kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out)
!======================================================================== DIT VOOR MOLECULEN GAAN DOEN
  particles_remaining = 0.d0
  momentum_remaining  = 0.d0
  energy_remaining    = 0.d0
  momentuminplasma = 0.d0
  parallelmomentuminplasma = 0.d0

  superparticles_remaining = 0 !< just for debugging. Use count_action in actual simulation.!.d0!.d0

  select type (particles => sim%groups(atoms)%particles)
  type is (particle_kinetic_leapfrog)
  
    !$omp parallel do default(none) &
    !$omp reduction(+:particles_remaining, momentum_remaining, energy_remaining,superparticles_remaining,lost_particle_weights,momentuminplasma, parallelmomentuminplasma) &
    !$omp shared(sim, particles, atoms) &
    !$omp private(j, E, B, psi, U, B_norm)
    !write(*,*) 'test123 we are in the if atom diagnostics statement'
    do j=1,size(particles,1) 

      !is nu na de molecule loop waar deeltjes metene weer gebruikt en overschreven worden
      ! if (particles(j)%i_elm .le. 0) then !sum weights of lost particles 
      !   !write(*,*) 'test21 i_elm', particles(j)%i_elm, 'weight =', particles(j)%weight
      !   lost_particle_weights = lost_particle_weights + particles(j)%weight
      ! endif

      if (particles(j)%i_elm .le. 0) cycle

      call sim%fields%calc_EBpsiU(sim%time , particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
      B_norm = B/norm2(B)

      particles_remaining = particles_remaining + particles(j)%weight
      momentum_remaining  = momentum_remaining  + particles(j)%weight * dot_product(B_norm,particles(j)%v) *sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
      !vx = cos(particles(j)%x(3))*particles(j)%v(1)-particles(j)%x(1)*sin(particles(j)%x(3))
      !vy = -sin(particles(j)%x(3))*particles(j)%v(1)-particles(j)%x(1)*cos(particles(j)%x(3))
      !vz = 1.d0
      momentuminplasma(1) = momentuminplasma(1) + particles(j)%weight * (cos(particles(j)%x(3))*particles(j)%v(1)-sin(particles(j)%x(3)) * particles(j)%v(3)) *sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
      momentuminplasma(2) = momentuminplasma(2) + particles(j)%weight * (-sin(particles(j)%x(3))*particles(j)%v(1)-cos(particles(j)%x(3)) * particles(j)%v(3)) *sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
      momentuminplasma(3) = momentuminplasma(3) + particles(j)%weight * particles(j)%v(2) * 1.d0 * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
      
      parallelmomentuminplasma = parallelmomentuminplasma + particles(j)%weight * dot_product(particles(j)%v, B_norm)*sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
      
      

      energy_remaining    = energy_remaining    + particles(j)%weight * dot_product(particles(j)%v,particles(j)%v) *sim%groups(atoms)%mass * ATOMIC_MASS_UNIT /2.d0
!      energy_remaining    = energy_remaining    + particles(j)%weight * 2.18d-15
      
      superparticles_remaining = superparticles_remaining + 1

      !ONLY USE FOR MOMENTUM CHECK
      !particles(j)%i_elm = 0 !parallel momentum is only conserved at the moment of creation
    enddo !j
    ! if (lost_particle_weights .gt. 0.d0) then !write if particles are lost
    !   write(*,*) 'test21 lost_particle_weights atoms', lost_particle_weights
    ! endif
	!omp end parallel do'


  end select
!========================================================================

  call MPI_REDUCE(particles_remaining, all_particles, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_remaining,  all_momentum,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_remaining,    all_energy,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(superparticles_remaining,all_superparticles,1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentuminplasma,allmomentuminplasma,3 , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(parallelmomentuminplasma, allparallelmomentuminplasma,1 , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

   
  !write(33,'(A,126e16.8)') ' TOTAL P% info : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
  !                                      mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
  !                                      pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot 
  
  if (sim%my_id .eq. 0) then
    write(*,'(A,3e16.8)') 'REMAINING (START) ATOMS : ',all_particles, all_momentum, all_energy

    write(*,'(A,126e16.8)') ' TOTAL1 : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
                                        mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
                                        pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot

    write (*,'(A,3e16.8)') 'Atom momentum in plasma: ', allmomentuminplasma
	!write(33,'(A,126e16.8)') ' TOTAL P% info : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
    !                                    mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
    !                                    pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot 


    write (*,'(A,3e16.8)') 'Parallel momentum given to atoms at this timestep: ', allparallelmomentuminplasma
	write(*,'(A,I13,A,E8.2,A,F13.10,A)') 'Superparticles in use :',all_superparticles,' of ', n_particles, '| in use :', &
				real(all_superparticles)/n_particles*100.d0,'%'
				
	if ( all_superparticles .gt. 0 )	then
		write(*,'(A,2E16.8)') 'Average weight of particles',(all_particles)/all_superparticles
	endif	!real(count(
	
  endif !(sim%my_id .eq. 0)  
!===================================================  End diagnostics for conservation  
!======================================================= Start diagnostics for molecules
  particles_remaining = 0.d0
  momentum_remaining  = 0.d0
  energy_remaining    = 0.d0
  superparticles_remaining = 0 !< just for debugging. Use count_action in actual simulation.!.d0!.d0
  select type (particles => sim%groups(molecules)%particles)
  type is (particle_kinetic_leapfrog)
  
    !$omp parallel do default(none) &
    !$omp reduction(+:particles_remaining, momentum_remaining, energy_remaining,superparticles_remaining) &
    !$omp shared(sim, particles, atoms,molecules) &
    !$omp private(j, E, B, psi, U, B_norm)
    do j=1,size(particles,1)

      if (particles(j)%i_elm .le. 0) cycle

      call sim%fields%calc_EBpsiU(sim%time , particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
      B_norm = B/norm2(B)

      particles_remaining = particles_remaining + particles(j)%weight
      momentum_remaining  = momentum_remaining  + particles(j)%weight * dot_product(B_norm,particles(j)%v) *sim%groups(molecules)%mass * ATOMIC_MASS_UNIT
      energy_remaining    = energy_remaining    + particles(j)%weight * dot_product(particles(j)%v,particles(j)%v) *sim%groups(molecules)%mass * ATOMIC_MASS_UNIT /2.d0
!      energy_remaining    = energy_remaining    + particles(j)%weight * 2.18d-15
      
      superparticles_remaining = superparticles_remaining + 1

    enddo !j
	!omp end parallel do


  end select
!========================================================================

  call MPI_REDUCE(particles_remaining, all_particles, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_remaining,  all_momentum,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_remaining,    all_energy,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(superparticles_remaining,all_superparticles,1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)   
  !write(33,'(A,126e16.8)') ' TOTAL P% info : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
  !                                      mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
  !                                      pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot 
  
  if (sim%my_id .eq. 0) then
    write(*,'(A,3e16.8)') 'REMAINING (START) MOLECULES : ',all_particles, all_momentum, all_energy

    ! write(*,'(A,126e16.8)') ' TOTAL1 : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
    !                                     mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
    !                                     pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot
										
	!write(33,'(A,126e16.8)') ' TOTAL P% info : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
    !                                    mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
    !                                    pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot 

	write(*,'(A,I13,A,E8.2,A,F13.10,A)') 'Molecular superparticles in use :',all_superparticles,' of ', n_particles, '| in use :', &
				real(all_superparticles)/n_particles*100.d0,'%'
				
	if ( all_superparticles .gt. 0 )	then
		write(*,'(A,2E16.8)') 'Average weight of molecular particles',(all_particles)/all_superparticles
	endif	!real(count(
	
  endif !(sim%my_id .eq. 0)  
!===================================================  End diagnostics for conservation    
end do ! while

call write_simulation_hdf5(sim, 'part_restart.h5')
! call write_simulation_hdf5(sim, 'restart_part.h5') !< make sure part_restart won't be overwritten (guido solution)

call sim%finalize

contains

!==================================================================================
!                                ATOM LOOP
!==================================================================================
subroutine loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)
use mod_project_particles
use mod_random_seed
use mod_interp, only: mode_moivre
use mod_basisfunctions
use mod_particle_types, only: copy_particle_kinetic_leapfrog
use mod_sampling, only: boxmueller_transform,sample_chi_squared_3

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(projection), target, intent(inout)                   :: jorek_feedback
type(count_action)                                        :: counter
type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
type(particle_kinetic_leapfrog)                           :: particle_tmp

real*8, intent(in)     :: timesteps, particle_start_time 
real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2)
! real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3)!, v
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: ion_rate, ion_source, ion_prob, ion_ran(1), cx_ran(8),st_ran(2), cx_source, cx_energy ,PLT
real*8    :: cx_prob, CX_rate
real*8    :: kinetic_energy, ion_energy,line_rad_energy
real*8    :: n_lost_ion, n_lost_ion_all, p_plt_lost,p_plt_lost_all,p_cx_lost,p_cx_lost_all,p_lost_ion,p_lost_ion_all
integer   :: n_super_ionized, n_super_ionized_all
real*8    :: particle_source, velocity_par_source, energy_source
real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v, v_v, v_E,extra_proj
real*8    :: vvector(3),sum_ran(3), E_th, v_th,ran_norm(4)

!$ real*8 :: w0, w1, mmm(3)

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_elm_old, i_elm 
integer   :: seed, i_rng, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: n_particles, ifail
logical   :: limits
real*8,allocatable :: feedback_rhs(:,:,:,:,:)

!$ w0 = omp_get_wtime()

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

  n_lost_ion = 0.d0
  n_lost_ion_all = 0.d0
  p_lost_ion   = 0.d0
  p_lost_ion_all   = 0.d0
  p_plt_lost  = 0.d0
  p_plt_lost_all  = 0.d0
  p_cx_lost   = 0.d0
  p_cx_lost_all   = 0.d0
  
  n_super_ionized = 0
  n_super_ionized_all = 0
  
  
jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps

allocate(feedback_rhs,source=jorek_feedback%rhs)

! jorek_feedback%rhs = 0.d0
feedback_rhs       = 0.d0

!call with(sim, counter) 

select type (particles => sim%groups(atoms)%particles)
type is (particle_kinetic_leapfrog)!momentum
#ifdef __GFORTRAN__
 !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else 
!$omp parallel do default(none) &
#endif
 !$omp schedule(dynamic,10) &
 !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
 !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm,                           &
 !$omp use_cx, use_ionisation,use_line_radiation, atoms,                                                   &
 !$omp CENTRAL_DENSITY, CENTRAL_MASS)                                              &
 !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old,    &
 !$omp i_elm_old, i_elm, n_e, T_e,                                                 &
 !$omp PLT,ion_rate, ion_prob, ion_ran, ion_source, ion_energy, kinetic_energy, line_rad_energy,       &  
 !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail,limits,    &
 !$omp CX_rate, CX_prob, CX_source, CX_energy, v, v_E, v_v,extra_proj,                        &
 !$omp particle_source, velocity_par_source, energy_source, v_temp, K_eV, T_eV, cx_ran,&
 !$omp E_th, v_th,sum_ran,vvector,ran_norm)                                                                 &
 !$omp reduction(+:feedback_rhs,n_lost_ion,p_plt_lost,p_cx_lost,p_lost_ion,n_super_ionized)
 
 ! shared jorek_feedback
 !private 
!velocity_par_source = 0.d0
 do j=1,size(particles,1)

    call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)
     i_rng = 1
     !$ i_rng = omp_get_thread_num()+1
	if (particle_tmp%i_elm .le. 0) cycle
    do k=1,n_steps

      if (particle_tmp%i_elm .le. 0) exit

      t = particle_start_time + (k-1)*timesteps

      call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
      rz_old    = particle_tmp%x(1:2)
      st_old    = particle_tmp%st
      i_elm_old = particle_tmp%i_elm
	  
	  !> in use ionisation as well?
	  call sim%fields%calc_NeTe(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), n_e, T_e)
	  ! call calc_ U ne Te vpar
	  ion_source = 0.d0
	  ion_energy = 0.d0
	  
	  call sim%fields%calc_vvector(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), vvector)
	  !vvector is fluid flow velocity [v_R, v_Z, v_phi] m/s
	  !TODO: add upper limits if necessary
	  limits = n_e .le. 1e14 .or. T_e * K_BOLTZ / EL_CHG .le. 1.d0 !ADAS limits
	  if (particle_tmp%weight .lt. 0.0d0) write(*,*) "Negative particle weight p(j)%w=", particle_tmp%weight
	  
	  !>for impurities, bremsstrahlung and CX radiation can be added here as well. (see W_rad_example)
	  line_rad_energy = 0.d0
	  if (use_line_radiation .and. .not. limits) then !< before or after Ionisation and CX ??
			call sim%groups(atoms)%ad%PLT%interp(int(particle_tmp%q), log10(n_e), log10(T_e), PLT) ! [J m^3/s]
			! call ad_deuterium%plt%interp( 1, ne_si_log10, Te_si_log10, LradDrays_T, dLradDrays_dT)
			line_rad_energy = n_e * particle_tmp%weight * PLT * timesteps
	  endif ! use_line_radiation
	  
	  if (use_ionisation .and. .not. limits) then
       
          call sim%groups(atoms)%ad%SCD%interp(int(particle_tmp%q), log10(n_e), log10(T_e), ion_rate) ! [m^3/s]
          ion_prob = 1.d0 - exp(-ion_rate * n_e * timesteps) ! [0] poisson point process, exponential 

          ! If the weight is to small throw away the particle with the probability, else reduce weight with ionising probability
          ion_source = 0.d0

          if (particle_tmp%weight .le. 1.0d9) then !1.0d9 !1.0d10 1.0d7
            call rng(i_rng)%next(ion_ran)
            if (ion_ran(1) .le. ion_prob) then
              particle_tmp%i_elm  = 0
              ion_source = particle_tmp%weight
			  !superparticles ionized
			  n_super_ionized = n_super_ionized +1
            else
              ion_source = 0.d0
            endif

          else 
            ion_source = particle_tmp%weight * ion_prob
            particle_tmp%weight = particle_tmp%weight * (1.d0 - ion_prob)
          endif 

          kinetic_energy = dot_product(particle_tmp%v,particle_tmp%v) *sim%groups(atoms)%mass * ATOMIC_MASS_UNIT /2.d0

          ion_energy     = kinetic_energy - binding_energy !<binding energy should be here
		  !<including binding energy will make ion_energy negative, so it becomes a sink for the plasma

	  endif ! use_ionisation

	  
	  ! Charge Exchange
	  ! It is assumed that we will have a exchange between hydrogen isotopes
	  v_temp    = particle_tmp%v
	  cx_source = 0.d0
	  cx_energy = 0.d0
	  
	  if (use_cx  .and. .not. limits) then !< CX uses adas as well. Te limit could be lower.
	  
          call sim%groups(atoms)%ad%CCD%interp(int(particle_tmp%q+1), log10(n_e), log10(T_e), CX_rate) ! [m^3/s]
          CX_prob = 1.d0 - exp(-CX_rate * n_e * timesteps)

          call rng(i_rng)%next(cx_ran)
           if (cx_ran(1) .le. CX_prob) then
            ! sample boltzman, randomize velocity
            T_eV = T_e * K_BOLTZ / EL_CHG !< T_eV = electron T in [eV]

			!============== NEW CX PARTICLE
			  !Box-Mueller sample velocities with st.dev=1
			  ran_norm = boxmueller_transform(cx_ran(2:5))
			  !>v_temp = sqrt(kT/m) * ran_norm
			  v_temp = sqrt(T_e * K_BOLTZ/(sim%groups(atoms)%mass * ATOMIC_MASS_UNIT))*ran_norm(2:4)
			  !write(*,*) "vtemp", v_temp
			  !>add bulk fluid flow
			  v_temp = v_temp + vvector 

              CX_source = particle_tmp%weight
              CX_energy   = 0.5d0 * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT *  (dot_product(particle_tmp%v,particle_tmp%v) - dot_product(v_temp,v_temp))
			
              !write(*,*) "neTe",n_e,T_e			
			  !write(*,*) "CX", vvector
          endif ! cx_ran
	  endif ! use_cx
	  
	  if (isnan(ion_source * ion_energy + cx_source * cx_energy - line_rad_energy)) then
		write(*,*) "ion_energy", ion_energy
		write(*,*) "cx_energy", cx_energy
		write(*,*) "line_rad_energy", line_rad_energy
		particle_tmp%i_elm  = 0
		CYCLE !< don't feed this particle into the feedback
		
	  endif
	  
	  
	  ! feedback from each particle at each timestep
	  energy_source       = ion_source * ion_energy + cx_source * cx_energy - line_rad_energy
	  particle_source     = ion_source * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT !< mass source in SI
	  velocity_par_source = ion_source * dot_product(B, particle_tmp%v) * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT &	
			+ CX_source  * dot_product(B, particle_tmp%v - v_temp) * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT 
			   
	  particle_tmp%v = v_temp 
	  n_lost_ion = n_lost_ion + ion_source	!< local sum #particles lost due to ionisation
	  p_lost_ion = p_lost_ion + ion_source * ion_energy
	  p_plt_lost = p_plt_lost + line_rad_energy
	  p_cx_lost  = p_cx_lost + cx_source * cx_energy
	  !Calculate the projection of the ion source in real-time
		call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
		call mode_moivre(particle_tmp%x(3), HZ)
			  
		do l=1,n_vertex_max
		  do m=1,n_order+1

			index_lm = (l-1)*(n_order+1) + m

			v   = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_source     * t_norm / rho_norm
			v_E = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source       * t_norm / E_norm
			v_v = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * velocity_par_source * t_norm / m_norm
			extra_proj = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) *particle_tmp%weight * 1.d0/real(n_steps,8) !<average density over jorek timesteps!real(floor(k/n_steps))!1.d0 !<density proj

			do i_tor=1,n_tor
			  feedback_rhs(m,l,i_elm_old,i_tor,1) = feedback_rhs(m,l,i_elm_old,i_tor,1) + HZ(i_tor) * v
			  feedback_rhs(m,l,i_elm_old,i_tor,2) = feedback_rhs(m,l,i_elm_old,i_tor,2) + HZ(i_tor) * v_E
			  feedback_rhs(m,l,i_elm_old,i_tor,3) = feedback_rhs(m,l,i_elm_old,i_tor,3) + HZ(i_tor) * v_v
			  feedback_rhs(m,l,i_elm_old,i_tor,4) = feedback_rhs(m,l,i_elm_old,i_tor,4) + HZ(i_tor) * extra_proj !< buiten de steps loop
			enddo

		  enddo
		enddo
	  
	  
      if (particle_tmp%i_elm .gt. 0) then
        ! Push the particle and determine it's new location.
        call boris_push_cylindrical(particle_tmp, sim%groups(atoms)%mass, E, B, timesteps)
        

        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                            particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
        if (ifail .lt. 0) then !if outside grid in new position, set i_elm =0 
          particle_tmp%i_elm=0
          !write(*,*) 'test125 particle at boundary, ifail=', ifail, 'weight', particle_tmp%weight
        endif

      endif
   
    enddo ! steps 

    call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))

  
  enddo   ! particles
  !$omp end parallel do
  
end select

if (use_ncs) then
    write(*,*) 'GATHER TIME : ',jorek_feedback%rhs_gather_time
    !jorek_feedback%rhs = feedback_rhs / jorek_feedback%rhs_gather_time !* TWOPI
	jorek_feedback%rhs(:,:,:,:,1:3) = feedback_rhs(:,:,:,:,1:3) / jorek_feedback%rhs_gather_time !* TWOPI
	jorek_feedback%rhs(:,:,:,:,4) = feedback_rhs(:,:,:,:,4)
    jorek_feedback%rhs_gather_time = 0.d0
else
    jorek_feedback%rhs = feedback_rhs 
  endif
  
deallocate(feedback_rhs)

call MPI_REDUCE(n_lost_ion, n_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(p_lost_ion, p_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(p_plt_lost, p_plt_lost_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(p_cx_lost, p_cx_lost_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(n_super_ionized, n_super_ionized_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

if (sim%my_id .eq. 0) write(*,'(A46,E14.6,I6)') "Lost superparticles at t due to ionisation: ", sim%time, n_super_ionized_all
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') " Lost particles at t due to ionisation: ", sim%time, n_lost_ion_all
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') " Ionization rate at time t [#/s]: ", sim%time, n_lost_ion_all / (timesteps * n_steps)
p_lost_ion_all = p_lost_ion_all / (timesteps * n_steps)
p_plt_lost_all = p_plt_lost_all / (timesteps * n_steps)
p_cx_lost_all = p_cx_lost_all / (timesteps * n_steps)
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Energy exchange to plasma [W] at t due to ionisation: ", sim%time, p_lost_ion_all ! energy gain
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost energy [W] at t due to line radiation: ", sim%time, p_plt_lost_all
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Energy exchange to plasma [W] at t due to CX radiation: ", sim%time, p_cx_lost_all
if (sim%my_id .eq. 0) write(*,'(A17,5E14.6)') 'TOTAL Exchange , delta t: ' ,sim%time,p_lost_ion_all, -p_plt_lost_all, p_cx_lost_all, timesteps * n_steps
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Total energy exchange to plasma [W]: ", sim%time, p_lost_ion_all -p_plt_lost_all+ p_cx_lost_all

! if (sim%my_id .eq. 0) write(*,*) " Lost energy [J] at t due to line radiation: ", sim%time, p_plt_lost_all
!$ w1 = omp_get_wtime()
!$ mmm = mpi_minmeanmax(w1-w0)
!$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"


!  write(*,*) 'CAREFUL: averaging over n_steps : ',n_steps
!  jorek_feedback%rhs = jorek_feedback%rhs / real(n_steps,8)
if (sim%my_id .eq. 0) write(*,*) 'done loop_particle_kinetic_local'

end subroutine

!==================================================================================
!                                MOLECULE LOOP
!==================================================================================
!feedback: op elke node, verschillende waardes voor terugkoppeling energy, momentum en dichtheid. RNG: random number generator. 
subroutine loop_particle_kinetic_molecule_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)
use mod_project_particles
use mod_random_seed
use mod_interp, only: mode_moivre
use mod_basisfunctions
use mod_particle_types, only: copy_particle_kinetic_leapfrog
use mod_sampling, only: boxmueller_transform,sample_chi_squared_3
use mod_find_free_particle, only: find_free_particles

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(projection), target, intent(inout)                   :: jorek_feedback
type(count_action)                                        :: counter
type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
type(particle_kinetic_leapfrog)                           :: particle_tmp
type(particle_kinetic_leapfrog)                           :: H2plus_tmp(2) 

!deze zijn atomair specifiek. Moet nog worden aangepast!. 
real*8, intent(in)     :: timesteps, particle_start_time 
real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2)
! real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3)!, v
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: ion_rate, ion_source, ion_prob, ion_ran(1), cx_ran(8),mol_cx_ran(5),st_ran(2), ion_diss_ran(4), cx_source, cx_energy ,PLT
real*8    :: cx_prob, CX_rate, diss_ran(4), diss_ion_ran(4), non_diss_ion_ran(1)
real*8    :: kinetic_energy, ion_energy,line_rad_energy
real*8    :: n_lost_ion, n_lost_ion_all, n_lost_diss_ion, p_plt_lost,p_plt_lost_all,p_cx_lost,p_cx_lost_all,p_lost_ion,p_lost_ion_all
integer   :: n_super_ionized, n_super_ionized_all
real*8    :: particle_source, velocity_par_source, energy_source
real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v, v_v, v_E,extra_proj
real*8    :: vvector(3),sum_ran(3), E_th, v_th,ran_norm(4)

! moleculair specifiek
real*8    :: electron_cooling_rate,electron_radiation_rate, combined_rate, kinetic_energy_H2_initial, potential_energy_H2_eV
real*8    :: diss_prob, diss_source, potential_energy_H2_J, atom_final_energy, atom_diss_final_speed, atom_diss_final_energy
real*8    :: electron_cooling_rate_diss_ion,electron_radiation_rate_diss_ion,  diss_ion_prob, diss_ion_source, atom_diss_ion_final_energy, atom_diss_ion_final_speed
real*8    :: non_diss_ion_prob, non_diss_ion_source, electron_cooling_rate_non_diss_ion, electron_radiation_rate_non_diss_ion, H2ionisationenergy
real*8    :: mol_cx_prob, mol_cx_source, mol_CX_energy
real*8    :: ion_diss_exc_prob, ion_diss_ion_prob, ion_diss_recomb_prob, sum_probabilities, ion_diss_exc_source(2), electron_cooling_rate_ion_diss_exc
real*8    :: ion_diss_exc_final_energy, ion_diss_exc_final_speed, v_ion_diss_exc(2,3), ion_diss_recomb_source(2), electron_cooling_rate_ion_diss_recomb(2)
real*8    :: ion_diss_recomb_final_energy(2), ion_diss_recomb_final_speed(2), ion_diss_ion_source(2), electron_cooling_rate_ion_diss_ion, ion_diss_ion_final_energy
real*8    :: ion_diss_ion_final_speed, ion1_diss_ion_velocity(2,3), ion2_diss_ion_velocity(2,3), ion_diss_exc_final_energy_coupling(2), ion2_diss_ion_energycoupling(2), ion1_diss_ion_energycoupling(2)


real*8, dimension(3)    :: rand_direc_vec,  rand_direc_vec_diss_ion, v_final_atom, v_final_molecular_ion, v_Hplus_tmp, v_mol_tmp
real*8, dimension(3)    :: random(3), v_diss_ion, Hvtemp(3)


!$ real*8 :: w0, w1, mmm(3)

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_p, i_f, i_elm_old, i_elm, counter2, ii
integer   :: seed, i_rng, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: n_particles, ifail
logical   :: limits, limitsAMJUEL
integer, dimension(:), allocatable ::  i_free
real*8,allocatable :: feedback_rhs(:,:,:,:,:)

!$ w0 = omp_get_wtime()
!write(*,*) 'test123 passed the init'
n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

Hvtemp = 0.d0
! OD diagnostics. Sommaties om mass en energie conservation bij te houden. 
  n_lost_ion = 0.d0 !hoeveel weight deeltjes is geioniseerd op 1 mpi process. 
  n_lost_ion_all = 0.d0 !over alle mpi processes 
  n_lost_diss_ion = 0.d0
  p_lost_ion   = 0.d0
  p_lost_ion_all   = 0.d0
  p_plt_lost  = 0.d0
  p_plt_lost_all  = 0.d0
  p_cx_lost   = 0.d0
  p_cx_lost_all   = 0.d0
  
  n_super_ionized = 0
  n_super_ionized_all = 0
  
    !============== Finding free particles 
  i_free = find_free_particles(sim%groups(atoms)%particles)
  !write(*,*) 'test124 = find_free_particles(sim%groups(atoms)%particles) i_free = ', i_free(1),i_free(2)
  !write(*,*) 'test124 = find_free_particles(sim%groups(atoms)%particles)', find_free_particles(sim%groups(atoms)%particles)
  !write(*,*) 'test126 = size(sim%groups(atoms)%particles, 1)', size(sim%groups(atoms)%particles, 1)
! ==================
  ! i_f = 1 !initialisatie voor de loop

#ifdef _OPENMP
i_f = -1 !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
!write(*,*) 'omp i_f =', i_f
#else
i_f = 1
!write(*,*) 'serial i_f =', i_f
#endif
!:q  write(*,*) 'first i_f =', i_f

  
jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps !kinetische tijdstappen voor normale tijdstap

allocate(feedback_rhs,source=jorek_feedback%rhs)

! jorek_feedback%rhs = 0.d0
feedback_rhs       = 0.d0

!call with(sim, counter)
!omp: miltithreading. Shared: elke cpu zelfde info. Private: anders. !$omp aanroep. 
select type (Hatom => sim%groups(atoms)%particles) !atoom maken/
type is (particle_kinetic_leapfrog)

select type (particles => sim%groups(molecules)%particles) !< group molecule (in this example it is group 2)
type is (particle_kinetic_leapfrog)



#ifdef __GFORTRAN__
 !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�

#else
 !$omp parallel do default(none) &
#endif
 !$omp schedule(dynamic,10) &
 !$omp firstprivate(i_f) &
 !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,i_free,        &
 !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm,                           &
 !$omp use_cx, use_ionisation,use_line_radiation,atoms,molecules, H2ionisationenergy,use_mol_cx, USE_H2plus_DISSOCIATION,use_molecules,                      &
 !$omp CENTRAL_DENSITY, CENTRAL_MASS,use_dissociation, potential_energy_H2_eV,potential_energy_H2_J,H2_NON_DISS_ION, use_dissionisation, use_nondissionisation,&
 !$omp H2plus_ELEC_COOL, H2plus_DISS_REC, H2plus_DISS_ION, H2plus_DISS_EXC, H2_ELEC_COOL,H2_DISS, H2_DISS_ION, H2_ion_con, H2_ELEC_RAD, Hatom  )                                              &
 !$omp private(particle_tmp, H2plus_tmp, i_rng, i,ii,j,k,l,m,i_p, t, E, B, psi, U, rz_old, st_old, counter2, v_diss_ion,  &
 !$omp i_elm_old, i_elm, n_e, T_e, Hvtemp,                                                &
 !$omp PLT,ion_rate, ion_prob, ion_ran, ion_source, ion_energy, kinetic_energy, line_rad_energy, diss_prob, diss_source, diss_ion_prob, diss_ion_source,      &
 !$omp non_diss_ion_prob, non_diss_ion_source, electron_cooling_rate_non_diss_ion, electron_radiation_rate_non_diss_ion, MOL_CX_prob, mol_CX_source,        &
 !$omp ion_diss_exc_prob, ion_diss_ion_prob, ion_diss_recomb_prob, sum_probabilities, ion_diss_exc_source, electron_cooling_rate_ion_diss_exc, &
 !$omp ion_diss_exc_final_energy, ion_diss_exc_final_speed, v_ion_diss_exc, ion_diss_recomb_source, electron_cooling_rate_ion_diss_recomb, &
 !$omp ion_diss_recomb_final_energy, ion_diss_recomb_final_speed, ion_diss_ion_source, electron_cooling_rate_ion_diss_ion, ion_diss_ion_final_energy, &
 !$omp ion_diss_ion_final_speed, ion1_diss_ion_velocity, ion2_diss_ion_velocity, &
 !$omp v_final_atom, v_final_molecular_ion, v_Hplus_tmp,         &
 !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail,limits,limitsAMJUEL,    &
 !$omp CX_rate, CX_prob, CX_source, CX_energy, v, v_E, v_v,extra_proj, v_mol_tmp,                 &
 !$omp particle_source, velocity_par_source, energy_source, v_temp, K_eV, T_eV, cx_ran,mol_cx_ran ,ion_diss_ran,diss_ran,diss_ion_ran,non_diss_ion_ran, &
 !$omp E_th, v_th,sum_ran,vvector,ran_norm, ion_diss_exc_final_energy_coupling, ion2_diss_ion_energycoupling, ion1_diss_ion_energycoupling,      &             
 !$omp electron_cooling_rate, electron_radiation_rate, electron_cooling_rate_diss_ion, electron_radiation_rate_diss_ion,   &
 !$omp combined_rate,atom_diss_final_energy,atom_diss_final_speed,random,rand_direc_vec,rand_direc_vec_diss_ion, atom_diss_ion_final_energy, atom_diss_ion_final_speed, mol_CX_energy)  &                                                        
 !$omp reduction(+:feedback_rhs,n_lost_ion,n_lost_diss_ion,p_plt_lost,p_cx_lost,p_lost_ion,n_super_ionized) 
  
! shared jorek_feedback
 !private reduction
!write(*,*) 'test123 1139: size(particles,1)', size(particles,1)

 

 do j=1,size(particles,1) ! loop over alle particles,.
  !$ if (i_f .eq. -1) i_f = omp_get_thread_num()+1
 !
  !write(*,*) 'i_f = ', i_f, 'omp_get_thread_num', omp_get_thread_num(), 'j', j
!write(*,*) "test123 in loop over all particles now: particle j =", j



 !write(*,*) 'test123 1140 i_f', i_f
    call copy_particle_kinetic_leapfrog(particles(j),particle_tmp) !particle temp/tmp: voorkomt telkens doorlopen van grote array
    call copy_particle_kinetic_leapfrog(particle_tmp, H2plus_tmp(1))
    call copy_particle_kinetic_leapfrog(particle_tmp, H2plus_tmp(2))
    
  !write(*,*) 'test123 n_steps for kinetic timesteps = ', n_steps
  !write(*,*) 'test123 particle_tmp%i_elm', particle_tmp%i_elm
  !write(*,*) 'test124 n_particles', n_particles
      i_rng = 1
 !$ i_rng = omp_get_thread_num()+1
  if (particle_tmp%i_elm .le. 0) cycle
  if (particle_tmp%weight .eq. 0.d0) cycle
    do k=1,n_steps !aantal kin steps per fluid timestep. 
      !initialisation of rates when not used
      electron_cooling_rate = 0.d0
      electron_radiation_rate = 0.d0
      electron_cooling_rate_diss_ion = 0.d0
      electron_radiation_rate_diss_ion = 0.d0
      electron_cooling_rate_non_diss_ion = 0.d0
      combined_rate = 1.d0 !not 0 because we divide by it!
      line_rad_energy = 0.d0
      cx_source = 0.d0
      cx_energy = 0.d0
      diss_source = 0.d0
      diss_ion_source = 0.d0
      non_diss_ion_source = 0.d0
      v_diss_ion = 0.d0
      atom_diss_ion_final_energy = 0.d0
      non_diss_ion_source = 0.d0
      electron_cooling_rate_non_diss_ion = 0.d0
      electron_radiation_rate_non_diss_ion = 0.d0
      mol_cx_source = 0.d0
      mol_cx_energy = 0.d0
      ion_diss_ion_source = 0.d0
      electron_cooling_rate_non_diss_ion = 0.d0
      ion_diss_ion_final_energy= 0.d0
      ion1_diss_ion_velocity= 0.d0
      ion2_diss_ion_velocity= 0.d0
      ion_diss_exc_source= 0.d0
      v_ion_diss_exc= 0.d0
      ion_diss_recomb_source = 0.d0
      ion_diss_recomb_final_energy = 0.d0
      electron_cooling_rate_ion_diss_ion = 0.d0
      ion_diss_exc_final_energy = 0.d0
      electron_cooling_rate_ion_diss_exc = 0.d0
      
      

      !write(*,*) 'test123 in kinetic loop now'
      !write(*,*) 'test123 k =', k
      !write(*,*) 'test123 inloop i_f =', i_f
      !write(*,*) 'test123 particle_tmp%i_elm', particle_tmp%i_elm
      if (particle_tmp%weight .le. 0.d0) then
        particle_tmp%i_elm = 0

      endif
      if (i_f .ge. size(sim%groups(atoms)%particles, 1)) write(*,*) 'test123 exiting because i_f > n_part'
      if (i_f .ge. size(sim%groups(atoms)%particles, 1)) exit !TODO adapt to size of atoms list:
      if (i_f .le. 0) write(*,*) 'i_f below zero!', i_f
      if (particle_tmp%i_elm .le. 0) write (*,*) ' exiting because particle tmp outside grid'
      if (particle_tmp%i_elm .le. 0) exit !zit niet in grid, dus doe er niks mee. !i_elm. Element nummer in fin element grid. 


      !write(*,*) 'before find free particles function: i_f = ', i_f
      !i_f is index van deeltjes die we willen creëeren. (1 tm 1000 bijv)


      !write(*,*) 'test123 after find free particles function: i_p = ', i_p
      t = particle_start_time + (k-1)*timesteps
	!x hieronde ris locatie in realspace. st is locatie in local space. x(3) derde richting (phi)
	  call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U) !uitlezen background field
      rz_old    = particle_tmp%x(1:2)
      st_old    = particle_tmp%st
      i_elm_old = particle_tmp%i_elm
	  !write(*,*) 'test123 line 1192'
	  !> in use ionisation as well?
	  call sim%fields%calc_NeTe(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), n_e, T_e)
	  ! call calc_ U ne Te vpar
	  ion_source = 0.d0
	  ion_energy = 0.d0
	  !write(*,*) 'test123 line 1198'
	  call sim%fields%calc_vvector(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), vvector)
	  !vvector is fluid flow velocity [v_R, v_Z, v_phi] m/s
	  !TODO: add upper limits if necessary
    limits = n_e .le. 1e14 .or. T_e * K_BOLTZ / EL_CHG .le. 1.d0 !ADAS limits
    limitsAMJUEL = .false. !todo: add limits
	  if (particle_tmp%weight .lt. 0.0d0) write(*,*) "Negative particle weight p(j)%w=", particle_tmp%weight
	  !write(*,*) 'test123 line 1205'

	  ! ------------------------------- ADD PHYSICS REACTION HERE ---------	  

    if (use_molecules .and. T_e*K_BOLTZ/EL_CHG .ge. 1.d0 ) then! .and. .not. limitsAMJUEL) then
    if (use_dissociation) then  

    call rng(i_rng)%next(diss_ran)
    
      !If the weight is to small throw away the particle with the probability, else reduce weight with ionising probability

      if (particle_tmp%weight .gt. 5.d16) then 
        
          diss_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_DISS,n_e,T_e) * n_e * timesteps) ! [0] poisson point process, exponential !amjuel rate toegevoegd
          diss_source = particle_tmp%weight * diss_prob !zoveel ioniseren we 
          particle_tmp%weight = particle_tmp%weight * (1.d0 - diss_prob) !nieuwe particle wieght
          !write(*,*) 'i_p test', i_free(1)

          i_p = i_free(i_f)
          !write(*,*) 'i_p,', i_p
          Hatom(i_p)%weight = diss_source  

          electron_cooling_rate = AMJUEL_rate_coeff_neTe(H2_ELEC_COOL,n_e,T_e)*EL_CHG ! [Js-1m-3]

          ! write(*,*) 'electron_cooling_rate', electron_cooling_rate


          electron_radiation_rate = AMJUEL_rate_coeff_neTe(H2_ELEC_RAD,n_e,T_e)*EL_CHG ! [Js-1m-3]
          combined_rate = AMJUEL_rate_coeff_neTe(H2_DISS,n_e,T_e)+AMJUEL_rate_coeff_neTe(H2_NON_DISS_ION,n_e,T_e)
          ! write(*,*) 'after combined rate'
          potential_energy_H2_eV = 4.48 ![eV]
          potential_energy_H2_J = potential_energy_H2_eV*EL_CHG
          atom_diss_final_energy      = 0.5d0*(electron_cooling_rate/combined_rate - electron_radiation_rate/combined_rate - potential_energy_H2_J) !<binding energy should be here !plasma meer energie: positief!
          
          atom_diss_final_speed       = sqrt(2.d0*atom_diss_final_energy/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))


          ! write(*,*) 'voor rand vec'
          !call random_number(random) !array met 3 random waarden (x,y,z) richting
          !random = random - 0.5d0
          rand_direc_vec = diss_ran(2:4)/sqrt((dot_product(diss_ran(2:4),diss_ran(2:4)))) !normaliseer naar unit vector
          write(*,*) 'randvec', rand_direc_vec
          ! write(*,*) 'na randvec'


          !write(*,*) 'rand_direc_vec length = ', sqrt(rand_direc_vec(1)*rand_direc_vec(1)+rand_direc_vec(2)*rand_direc_vec(2)+rand_direc_vec(3)*rand_direc_vec(3))
          Hatom(i_p)%v = particle_tmp%v+rand_direc_vec*atom_diss_final_speed    !geef Hatom snelheid van molecuul + snelheid van botsing in random richting
          Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
          Hatom(i_p)%st    = st_old
          Hatom(i_p)%i_elm = i_elm_old
          Hatom(i_p)%q = 0
          Hatom(i_p)%v(3) = Hatom(i_p)%v(3)!+1000.d0
          !Hvtemp(1) = 

          #ifdef _OPENMP
          !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
          !write(*,*) 'omp i_f =', i_f
          #else
          i_f =i_f +1 !if no OMP
          #endif

          
          

          i_p = i_free(i_f) !switch naar nieuw vrij atoom
          Hatom(i_p)%weight = diss_source  
          Hatom(i_p)%v =  particle_tmp%v-rand_direc_vec*atom_diss_final_speed    !geef Hatom snelheid van molecuul - snelheid van botsing in random richting (momentum cons)
          !write(*,*) 'test12345', Hatom(i_p)%v, Hatom(i_p)%v+Hvtemp

      
          Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
          Hatom(i_p)%st    = st_old
          Hatom(i_p)%i_elm = i_elm_old
          Hatom(i_p)%q = 0

          #ifdef _OPENMP
          !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
          !write(*,*) 'omp i_f =', i_f
          #else
          i_f =i_f +1 !if no OMP
          #endif

          if (i_f .ge. size(sim%groups(atoms)%particles, 1)) write(*,*) 'exited because of indice outside of array'
          if (i_f .ge. size(sim%groups(atoms)%particles, 1)) exit  


          !write(*,'(a,2e16.8,a,3e16.8)') 'atom_diss_final_energy energie van 1 deeltje = ', atom_diss_final_energy/EL_CHG,1.d0*ATOMIC_MASS_UNIT*dot_product(Hatom(i_p)%v,Hatom(i_p)%v)/EL_CHG, 'Hatom%v',Hatom(i_p)%v
          !i_p = i_free(i_f) !switch naar nieuw vrij atoom 
           !particle_tmp%weight .gt. 1.0d15) then 

          else!if (particle_tmp%weight .le. 1.0d15) then 
            !write(*,*) 'test12345 enter else'
          !call random_number(random) 
          !random(1) = 0.d0 !to always enter loop
          diss_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_DISS,n_e,T_e) * n_e * timesteps)
          
              if (diss_ran(1) .le. diss_prob) then !
                  diss_prob = 1.d0
                  diss_source = particle_tmp%weight * diss_prob 
                  particle_tmp%weight = 0.d0
                  

                  i_p = i_free(i_f)
                  Hatom(i_p)%weight = diss_source  
        
                  electron_cooling_rate = AMJUEL_rate_coeff_neTe(H2_ELEC_COOL,n_e,T_e)*EL_CHG ! [Js-1m-3]
                  electron_radiation_rate = AMJUEL_rate_coeff_neTe(H2_ELEC_RAD,n_e,T_e)*EL_CHG ! [Js-1m-3]
                  combined_rate = AMJUEL_rate_coeff_neTe(H2_DISS,n_e,T_e)+AMJUEL_rate_coeff_neTe(H2_NON_DISS_ION,n_e,T_e)
        
                  potential_energy_H2_eV = 4.48 ![eV]
                  potential_energy_H2_J = potential_energy_H2_eV*EL_CHG
                  atom_diss_final_energy      = 0.5d0*(electron_cooling_rate/combined_rate - electron_radiation_rate/combined_rate - potential_energy_H2_J) !<binding energy should be here !plasma meer energie: positief!
                  
                  atom_diss_final_speed       = sqrt(2.d0*atom_diss_final_energy/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))
           
                  
                  !random = random - 0.5d0
                  rand_direc_vec = diss_ran(2:4)/sqrt((dot_product(diss_ran(2:4),diss_ran(2:4)))) !normaliseer naar unit vector
                  
        
        
                  !write(*,*) 'rand_direc_vec length = ', sqrt(rand_direc_vec(1)*rand_direc_vec(1)+rand_direc_vec(2)*rand_direc_vec(2)+rand_direc_vec(3)*rand_direc_vec(3))
                  Hatom(i_p)%v = particle_tmp%v+rand_direc_vec*atom_diss_final_speed    !geef Hatom snelheid van molecuul + snelheid van botsing in random richting
                  Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
                  Hatom(i_p)%st    = st_old
                  Hatom(i_p)%i_elm = i_elm_old
                  Hatom(i_p)%q = 0

        
                  #ifdef _OPENMP
                  !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
                  !write(*,*) 'omp i_f =', i_f
                  #else
                  i_f =i_f +1 !if no OMP
                  #endif

                  i_p = i_free(i_f) !switch naar nieuw vrij atoom
                  Hatom(i_p)%weight = diss_source  
                  Hatom(i_p)%v =  particle_tmp%v-rand_direc_vec*atom_diss_final_speed    !geef Hatom snelheid van molecuul - snelheid van botsing in random richting (momentum cons)
                  !write(*,*) 'test12345', Hatom(i_p)%v, Hatom(i_p)%v+Hvtemp
        
              
                  Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
                  Hatom(i_p)%st    = st_old
                  Hatom(i_p)%i_elm = i_elm_old
                  Hatom(i_p)%q = 0
                  if (i_f .ge. size(sim%groups(atoms)%particles, 1)) write(*,*) 'exited because of indice outside of array'
                  if (i_f .ge. size(sim%groups(atoms)%particles, 1)) exit  
                  
                  #ifdef _OPENMP
                  !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
                  !write(*,*) 'omp i_f =', i_f
                  #else
                  i_f =i_f +1 !if no OMP
                  #endif

        
                  !write(*,'(a,2e16.8,a,3e16.8)') 'atom_diss_final_energy energie van 1 deeltje = ', atom_diss_final_energy/EL_CHG,1.d0*ATOMIC_MASS_UNIT*dot_product(Hatom(i_p)%v,Hatom(i_p)%v)/EL_CHG, 'Hatom%v',Hatom(i_p)%v
                  !i_p = i_free(i_f) !switch naar nieuw vrij atoom 
                  particle_tmp%i_elm  = 0 !bring particle out of domain
               !(ion_ran(1) .le. diss_ion_prob)
                
              

              endif !(ion_ran(1) .le. diss_prob) then ! als probability hoger dan random waarde
      endif !(particle_tmp%weight .le. 1.0d9) then
      
      endif !  (use_dissociation) then  



      if (use_dissionisation) then
      
      call rng(i_rng)%next(diss_ion_ran)


      if (particle_tmp%weight .gt. 5.d16) then !5.d13
        diss_ion_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_DISS_ION,n_e,T_e) * n_e * timesteps) ! [0] poisson point process, exponential !amjuel rate toegevoegd
        diss_ion_source = particle_tmp%weight * diss_ion_prob !zoveel ioniseren we 
        particle_tmp%weight = particle_tmp%weight * (1.d0 - diss_ion_prob) !nieuwe particle wieght

     

        i_p = i_free(i_f)

        Hatom(i_p)%weight = diss_ion_source  

        electron_cooling_rate_diss_ion = 2.808d1*EL_CHG ! [Js-1m-3]
        electron_radiation_rate_diss_ion = 0.d0!  ! [Js-1m-3]
      
        atom_diss_ion_final_energy      =  5.d0*EL_CHG ![J] 5eV
        atom_diss_ion_final_speed       = sqrt(2.d0*atom_diss_ion_final_energy/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))
 
        
        diss_ion_ran(2:4) = diss_ion_ran(2:4) - 0.5d0 !change random values to between -0.5 and 0.5
        rand_direc_vec_diss_ion = diss_ion_ran(2:4)/sqrt((dot_product(diss_ion_ran(2:4),diss_ion_ran(2:4)))) !normaliseer naar unit vector     

        Hatom(i_p)%v = particle_tmp%v+rand_direc_vec_diss_ion*atom_diss_ion_final_speed    !geef Hatom snelheid van molecuul + snelheid van botsing in random richting
        v_diss_ion = particle_tmp%v-rand_direc_vec_diss_ion*atom_diss_ion_final_speed
        Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
        Hatom(i_p)%st    = st_old
        Hatom(i_p)%i_elm = i_elm_old
        Hatom(i_p)%q = 0


        if (i_f .ge. size(sim%groups(atoms)%particles, 1)) write(*,*) 'exited because of indice outside of array'
        if (i_f .ge. size(sim%groups(atoms)%particles, 1)) exit  
        #ifdef _OPENMP
        !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
        !write(*,*) 'omp i_f =', i_f
        #else
        i_f =i_f +1 !if no OMP
        #endif
        


        else!if (particle_tmp%weight .le. 1.0d15) then 
          
        
        diss_ion_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_DISS_ION,n_e,T_e) * n_e * timesteps) ! 
            if (diss_ion_ran(1) .le. diss_ion_prob) then !
                diss_ion_prob = 1.d0
                diss_ion_source = particle_tmp%weight * diss_ion_prob 
                particle_tmp%weight = 0.d0
                
                i_p = i_free(i_f)
                Hatom(i_p)%weight = diss_ion_source  
      
                electron_cooling_rate_diss_ion = 2.808d1*EL_CHG ! [Js-1m-3]
                electron_radiation_rate_diss_ion = 0.d0!  ! [Js-1m-3]
              
                atom_diss_ion_final_energy      =  5.d0*EL_CHG ![J]
                atom_diss_ion_final_speed       = sqrt(2.d0*atom_diss_ion_final_energy/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))
         
                
                diss_ion_ran(2:4) = diss_ion_ran(2:4) - 0.5d0 !change random values to between -0.5 and 0.5
                rand_direc_vec_diss_ion = diss_ion_ran(2:4)/sqrt((dot_product(diss_ion_ran(2:4),diss_ion_ran(2:4)))) !normaliseer naar unit vector

                Hatom(i_p)%v = particle_tmp%v+rand_direc_vec_diss_ion*atom_diss_ion_final_speed    !geef Hatom snelheid van molecuul + snelheid van botsing in random richting
                v_diss_ion = particle_tmp%v-rand_direc_vec_diss_ion*atom_diss_ion_final_speed
                Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
                Hatom(i_p)%st    = st_old
                Hatom(i_p)%i_elm = i_elm_old
                Hatom(i_p)%q = 0
      
 
      
      
                !write(*,'(a,2e16.8,a,3e16.8)') 'atom_diss_final_energy energie van 1 deeltje = ', atom_diss_final_energy/EL_CHG,1.d0*ATOMIC_MASS_UNIT*dot_product(Hatom(i_p)%v,Hatom(i_p)%v)/EL_CHG, 'Hatom%v',Hatom(i_p)%v
                !i_p = i_free(i_f) !switch naar nieuw vrij atoom 
                particle_tmp%i_elm  = 0 !bring particle out of domain
                #ifdef _OPENMP
                !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
                !write(*,*) 'omp i_f =', i_f
                #else
                i_f =i_f +1 !if no OMP
                #endif

             !(ion_ran(1) .le. diss_ion_prob)v_temp
              
            

            endif !(ion_ran(1) .le. diss_prob) then ! als probability hoger dan random waarde
    endif !(particle_tmp%weight .le. 1.0d9) then

      endif !use_dissionisation then
        !write(*,*) 'test987, before diss'


    if (use_nondissionisation) then      
          if (particle_tmp%weight .gt. 5.d16) then !1.d6 
            !write(*,*) exp(-AMJUEL_rate_coeff_neTe(H2_NON_DISS_ION,n_e,T_e))
            !write(*,*) 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_NON_DISS_ION,n_e,T_e) * n_e * timesteps)
        
              non_diss_ion_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_NON_DISS_ION,n_e,T_e) * n_e * timesteps) ! [0] poisson point process, exponential !amjuel rate toegevoegd
              non_diss_ion_source = particle_tmp%weight * non_diss_ion_prob !zoveel ioniseren we 
              particle_tmp%weight = particle_tmp%weight * (1.d0 - non_diss_ion_prob) !nieuwe particle wieght
                
              electron_cooling_rate_non_diss_ion = 1.5386d1*EL_CHG ! [Js-1m-3]
              electron_radiation_rate_non_diss_ion = 0.d0*EL_CHG ! [Js-1m-3]
            
              !H2ionisationenergy = 1.5386d1*EL_CHG
              
              H2plus_tmp(1)%weight = non_diss_ion_source 
              H2plus_tmp(1)%v = particle_tmp%v !geef Hatom snelheid van molecuul + snelheid van botsing in random richting
              H2plus_tmp(1)%q = 1
            
              if (i_f .ge. size(sim%groups(atoms)%particles, 1)) write(*,*) 'exited because of indice outside of array'
              if (i_f .ge. size(sim%groups(atoms)%particles, 1)) exit  
    
              else!if (particle_tmp%weight .le. 1.0d15) then 

                  call rng(i_rng)%next(non_diss_ion_ran)
                  non_diss_ion_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2_NON_DISS_ION,n_e,T_e) * n_e * timesteps)

                  if ((non_diss_ion_ran(1)) .le. non_diss_ion_prob) then !
                      non_diss_ion_prob = 1.d0
                      non_diss_ion_source = particle_tmp%weight * non_diss_ion_prob 
                      particle_tmp%weight = 0.d0
                      
                      electron_cooling_rate_non_diss_ion = 1.5386d1*EL_CHG ! [Js-1m-3]
                      electron_radiation_rate_non_diss_ion = 0.d0*EL_CHG ! [Js-1m-3]
                    
                      !H2ionisationenergy = 1.5386d1*EL_CHG
                      
                      H2plus_tmp(1)%weight = non_diss_ion_source 
                      H2plus_tmp(1)%v = particle_tmp%v !geef Hatom snelheid van molecuul + snelheid van botsing in random richting
                      H2plus_tmp(1)%q = 1
                      
                      particle_tmp%i_elm  = 0 !bring particle out of domain
                    
                  
    
                  endif !(ion_ran(1) .le. non_diss_ion_prob) then ! als probability hoger dan random waarde
          endif !(particle_tmp%weight .le. 1.0d9) then
          
          endif !  (use_nondissionisation) then  



                      ! Charge Exchange
          ! It is assumed that we will have a exchange between hydrogen isotopes

          mol_cx_source = 0.d0
          mol_cx_energy = 0.d0
          v_Hplus_tmp = 0.d0
          
          if (use_mol_cx) then !< 
            
            mol_CX_prob = 1.d0 - exp(-AMJUEL_rate_coeff_Te(H2_ION_CON,T_e) * n_e * timesteps)
            call rng(i_rng)%next(mol_cx_ran)
  
            if (mol_cx_ran(1) .le. mol_CX_prob) then
            ! sample boltzman, randomize velocity
            !============== NEW CX PARTICLE
              !Box-Mueller sample velocities with st.dev=1
              !ran_norm = boxmueller_transform(cx_ran(2:5))
              !>v_temp = sqrt(kT/m) * ran_norm

            mol_CX_source = particle_tmp%weight
            v_mol_tmp    = particle_tmp%v
            
            ran_norm = boxmueller_transform(mol_cx_ran(2:5))

            v_Hplus_tmp = sqrt(T_e * K_BOLTZ/(sim%groups(atoms)%mass * ATOMIC_MASS_UNIT))*ran_norm(2:4)+vvector

            v_final_atom = v_mol_tmp * sim%groups(molecules)%mass/sim%groups(atoms)%mass
            v_final_molecular_ion = (v_Hplus_tmp) * sim%groups(atoms)%mass/sim%groups(molecules)%mass
            
            
            i_p = i_free(i_f)

            Hatom(i_p)%weight = mol_CX_source 
            Hatom(i_p)%v = v_final_atom
            Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
            Hatom(i_p)%st    = st_old
            Hatom(i_p)%i_elm = i_elm_old
            Hatom(i_p)%q = 0


            H2plus_tmp(2)%v = v_final_molecular_ion

            
            mol_CX_energy   =  0.5d0 * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT * dot_product(v_Hplus_tmp,v_Hplus_tmp)


            particle_tmp%weight = 0.d0
            particle_tmp%i_elm  = 0 
            #ifdef _OPENMP
            !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
            !write(*,*) 'omp i_f =', i_f
            #else
            i_f =i_f +1 !if no OMP
            #endif

            endif ! (random(1) .le. mol_CX_prob) 
          endif ! (use_mol_cx)

          if (use_H2plus_dissociation) then 

            !REMOVE WHEN RUNNING DIFFERENT THAN DISSOCIATION CHECKS

            !H2plus_tmp(1)%weight = 1.d10
            !H2plus_tmp(2)%weight = 1.d10
            !H2plus_tmp(1)%v = 0.d0
            !H2plus_tmp(2)%v = 0.d0

            !H2plus_tmp(1)%v(3) = -1.d2
            !H2plus_tmp(2)%v(3) = -1.d2




            !END REMOVE


            do ii = 1, 2
            
            ion_diss_exc_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2plus_DISS_EXC,n_e,T_e) * n_e * timesteps)
            ion_diss_ion_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2plus_DISS_ION,n_e,T_e) * n_e * timesteps)
            ion_diss_recomb_prob = 1.d0 - exp(-AMJUEL_rate_coeff_neTe(H2plus_DISS_REC,n_e,T_e) * n_e * timesteps) 

            sum_probabilities = ion_diss_exc_prob + ion_diss_ion_prob + ion_diss_recomb_prob
            !write(*,*) 'Sum, dissexc, dission, dissrecomb', sum_probabilities, ion_diss_exc_prob, ion_diss_ion_prob, ion_diss_recomb_prob
            call rng(i_rng)%next(ion_diss_ran)
            ion_diss_ran(2:4) = ion_diss_ran(2:4) - 0.5d0
            rand_direc_vec = ion_diss_ran(2:4)/sqrt((dot_product(ion_diss_ran(2:4),ion_diss_ran(2:4)))) !normaliseer naar unit vector

            !write(*,*) 'testverdel1' ,ion_diss_ran(1)*sum_probabilities
            !write(*,*) 'exc, ion, recomb', ion_diss_exc_prob, ion_diss_ion_prob, ion_diss_recomb_probe
!            if (ion_diss_ran(1)*sum_probabilities .le. ion_diss_exc_prob) then
            if (ion_diss_ran(1)*sum_probabilities .le. ion_diss_exc_prob) then
            


              ion_diss_exc_source(ii) = H2plus_tmp(ii)%weight * 1.d0! ion_diss_exc_prob !zoveel ioniseren we 
              electron_cooling_rate_ion_diss_exc = 1.05d1*EL_CHG ! [Js-1m-3]


              
 
            
              ion_diss_exc_final_energy      =  4.3d0*EL_CHG ![J]
              ion_diss_exc_final_speed       = sqrt(2.d0*ion_diss_exc_final_energy/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))

              i_p = i_free(i_f)
              Hatom(i_p)%weight = ion_diss_exc_source(ii) 
              Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
              Hatom(i_p)%st    = st_old
              Hatom(i_p)%i_elm = i_elm_old
              Hatom(i_p)%q = 0
              Hatom(i_p)%v = H2plus_tmp(ii)%v+rand_direc_vec*ion_diss_exc_final_speed    !geef Hatom snelheid van molecuul + snelheid van botsing in random richting


              

              v_ion_diss_exc(ii,1:3) = particle_tmp%v-rand_direc_vec*ion_diss_exc_final_speed
              ion_diss_exc_final_energy_coupling(ii) = 0.5d0 * dot_product(v_ion_diss_exc(ii,1:3),v_ion_diss_exc(ii,1:3)) * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT

              #ifdef _OPENMP
              !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
             ! write(*,*) 'omp i_f =', i_f
              #else
              i_f =i_f +1 !if no OMP
              #endif

              


              H2plus_tmp(ii)%weight = particle_tmp%weight * (1.d0 - ion_diss_exc_prob) !nieuwe particle wieght (=0)
  


            
              elseif ((ion_diss_ran(1)*sum_probabilities .gt. ion_diss_exc_prob) .and. ion_diss_ran(1)*sum_probabilities .le. (ion_diss_exc_prob + ion_diss_recomb_prob)) then
             
            

              ion_diss_recomb_source(ii) = H2plus_tmp(ii)%weight * 1.d0 !ion_diss_recomb_prob !zoveel ioniseren we 
              electron_cooling_rate_ion_diss_recomb(ii) =  (AMJUEL_rate_coeff_Te(H2plus_ELEC_COOL,T_e)/AMJUEL_rate_coeff_neTe(H2plus_DISS_REC,n_e,T_e))*EL_CHG ! [Js-1m-3]
              ion_diss_recomb_final_energy(ii) = 0.5d0*max((electron_cooling_rate_ion_diss_recomb(ii)-(1.35d0+1.36d1/(1.5d0**2))*EL_CHG),0.d0) ![J]
              ion_diss_recomb_final_speed(ii)       = sqrt(2.d0*ion_diss_recomb_final_energy(ii)/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))


              
              i_p = i_free(i_f)
              Hatom(i_p)%weight = ion_diss_recomb_source(ii)
              Hatom(i_p)%v = H2plus_tmp(ii)%v + ion_diss_recomb_final_speed(ii) * rand_direc_vec
              Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
              Hatom(i_p)%st    = st_old
              Hatom(i_p)%i_elm = i_elm_old
              Hatom(i_p)%q = 0


              #ifdef _OPENMP
              !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
              !write(*,*) 'omp i_f =', i_f
              #else
              i_f =i_f +1 !if no OMP
              #endif


              i_p = i_free(i_f)

              Hatom(i_p)%weight = ion_diss_recomb_source(ii)
              Hatom(i_p)%v = H2plus_tmp(ii)%v - ion_diss_recomb_final_speed(ii) * rand_direc_vec
              Hatom(i_p)%x(1:3)   = particle_tmp%x(1:3)
              Hatom(i_p)%st    = st_old
              Hatom(i_p)%i_elm = i_elm_old
              Hatom(i_p)%q = 0

              H2plus_tmp(ii)%weight = particle_tmp%weight * (1.d0 - ion_diss_recomb_prob)
              


              #ifdef _OPENMP
              !$ i_f = i_f + omp_get_num_threads() !< if we run omp, we need a special first i_f value to then share it over the num_threads with unique values
              !write(*,*) 'omp i_f =', i_f
              #else
              i_f =i_f +1 !if no OMP
              #endif



            else


              ion_diss_ion_source(ii) = H2plus_tmp(ii)%weight * 1.d0 !ion_diss_ion_prob 
              electron_cooling_rate_ion_diss_ion = 1.55d1*EL_CHG ! [Js-1m-3]
              ion_diss_ion_final_energy = 2.5d-1 * EL_CHG ![J]
              ion_diss_ion_final_speed       = sqrt(2.d0*ion_diss_ion_final_energy/(sim%groups(atoms)%mass*ATOMIC_MASS_UNIT))

              
              ion1_diss_ion_velocity(ii,1:3) = H2plus_tmp(ii)%v + ion_diss_ion_final_speed * rand_direc_vec
              ion2_diss_ion_velocity(ii,1:3) = H2plus_tmp(ii)%v - ion_diss_ion_final_speed * rand_direc_vec

              ion1_diss_ion_energycoupling(ii) = 0.5d0*dot_product(ion1_diss_ion_velocity(ii,1:3),ion1_diss_ion_velocity(ii,1:3))* sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
              ion2_diss_ion_energycoupling(ii) = 0.5d0*dot_product(ion2_diss_ion_velocity(ii,1:3),ion2_diss_ion_velocity(ii,1:3))* sim%groups(atoms)%mass * ATOMIC_MASS_UNIT
              H2plus_tmp(ii)%weight = particle_tmp%weight * (1.d0 - ion_diss_ion_prob) !nieuwe particle wieght (=0)



              
            endif !(random(1)*sum_probabilities .le. ion_diss_exc_prob) 

            enddo !ii = 1, 2

          endif !if (use_H2_dissociation) 



    
      endif ! use_molecules
      


	  
	  ! ------------------------------- END OF PHYSICS REACTION PART ---------	  
	  
	  !> check for NaNs !voeg 
	  ! if (isnan(diss_source * atom_diss_final_energy + cx_source * cx_energy - line_rad_energy)) then
		! write(*,*) "ion_energy", ion_energy
		! write(*,*) "cx_energy", cx_energy
		! write(*,*) "line_rad_energy", line_rad_energy
		! particle_tmp%i_elm  = 0
		! CYCLE !< don't feed this particle into the feedback
	  
	  ! endif !isnan(diss_source * atom_diss_final_energy + cx_source * cx_energy - line_rad_energy)) 
	  !write(*,*) 'after NAN Check'
    !write(*,*) "diss_source", diss_source
    !write(*,*) "2nd Eterm ", diss_ion_source * (-electron_cooling_rate_diss_ion + atom_diss_ion_final_energy)
    !write(*,*) "electron_cooling_rate/combined_rate ", electron_cooling_rate/combined_rate 
    !write(*,*) 'before energy source'
    ! feedback from each particle at each timestep
    !write(*,*)           -diss_source * electron_cooling_rate/combined_rate 
    !write(*,*)   + diss_ion_source * -electron_cooling_rate_diss_ion  
    !write(*,*)  + atom_diss_ion_final_energy 
    !write(*,*)  + ion_diss_recomb_source * (-ion_diss_recomb_final_energy)
    !write(*,*)   - mol_cx_source * mol_CX_energy 
    !write(*,*)   - non_diss_ion_source*electron_cooling_rate_non_diss_ion 

    !write(*,*) 'before rror?'
    !write(*,*) ion_diss_ion_source
    !write(*,*) ion_diss_ion_final_energy
    !write(*,*) electron_cooling_rate_ion_diss_ion
    !write(*,*)    + ion_diss_ion_source*(2.d0*ion_diss_ion_final_energy -electron_cooling_rate_ion_diss_ion) 
    !write(*,*)   + ion_diss_exc_source * (ion_diss_exc_final_energy - electron_cooling_rate_ion_diss_exc) 


    energy_source       = -diss_source * electron_cooling_rate/combined_rate &
                        + diss_ion_source * (-electron_cooling_rate_diss_ion + atom_diss_ion_final_energy) &
                        + ion_diss_recomb_source(1) * (-electron_cooling_rate_ion_diss_recomb(1)) &
                        + ion_diss_recomb_source(2) * (-electron_cooling_rate_ion_diss_recomb(2)) &
                        - mol_cx_source * mol_CX_energy &
                        - non_diss_ion_source*electron_cooling_rate_non_diss_ion &
                        + ion_diss_ion_source(1)*(ion1_diss_ion_energycoupling(1) + ion2_diss_ion_energycoupling(1) - electron_cooling_rate_ion_diss_ion) &
                        + ion_diss_ion_source(2)*(ion1_diss_ion_energycoupling(2) + ion2_diss_ion_energycoupling(2) - electron_cooling_rate_ion_diss_ion) &
                        + ion_diss_exc_source(1) * (ion_diss_exc_final_energy_coupling(1) - electron_cooling_rate_ion_diss_exc) &
                        + ion_diss_exc_source(2) * (ion_diss_exc_final_energy_coupling(2) - electron_cooling_rate_ion_diss_exc)        
    particle_source     = (diss_ion_source  - mol_cx_source + ion_diss_ion_source(1)*2.d0 + ion_diss_ion_source(2)*2.d0 + ion_diss_exc_source(1)+ion_diss_exc_source(2)) &
                        * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT  !
    velocity_par_source = (diss_ion_source * dot_product(B, v_diss_ion ) &
                        - mol_cx_source * dot_product(B,v_Hplus_tmp)  &
                        + ion_diss_ion_source(1) * (dot_product(B, ion1_diss_ion_velocity(1,1:3)) + dot_product(B, ion2_diss_ion_velocity(1,1:3) )) &
                        + ion_diss_ion_source(2) * (dot_product(B, ion1_diss_ion_velocity(2,1:3)) + dot_product(B, ion2_diss_ion_velocity(2,1:3) )) &
                        + ion_diss_exc_source(1) * dot_product(B,v_ion_diss_exc(1,1:3))+ ion_diss_exc_source(2) * dot_product(B,v_ion_diss_exc(2,1:3)))&
                        * sim%groups(atoms)%mass * ATOMIC_MASS_UNIT !&	!=0, electronen geen mom



!    write(*,*) 'n_lost_ion'
    n_lost_ion = n_lost_ion + ion_source	!< local sum #particles lost due to ionisation
    
    n_lost_diss_ion = n_lost_diss_ion + diss_ion_source

 !   write(*,*) 'n_lost_ion'
    p_lost_ion = p_lost_ion + ion_source * ion_energy
 !   write(*,*) 'n_lost_ion'
    p_plt_lost = p_plt_lost + line_rad_energy
!    write(*,*) 'n_lost_ion'
    p_cx_lost  = p_cx_lost + cx_source * cx_energy
!    write(*,*) 'n_lost_ion'
    !Calculate the projection of the ion source in real-time
    !write(*,*) 'particletmp%st', particle_tmp%st(1), particle_tmp%st(2)
 !   write(*,*) 'test987', particle_tmp%st(1), particle_tmp%st(2)
		call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
		call mode_moivre(particle_tmp%x(3), HZ)
			  
		do l=1,n_vertex_max
		  do m=1,n_order+1

			index_lm = (l-1)*(n_order+1) + m

			v   = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_source     * t_norm / rho_norm
			v_E = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source       * t_norm / E_norm
			v_v = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * velocity_par_source * t_norm / m_norm
			extra_proj = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) *particle_tmp%weight * 1.d0/real(n_steps,8) !<average density over jorek timesteps!real(floor(k/n_steps))!1.d0 !<density proj

			do i_tor=1,n_tor
			  feedback_rhs(m,l,i_elm_old,i_tor,1) = feedback_rhs(m,l,i_elm_old,i_tor,1) + HZ(i_tor) * v
			  feedback_rhs(m,l,i_elm_old,i_tor,2) = feedback_rhs(m,l,i_elm_old,i_tor,2) + HZ(i_tor) * v_E
			  feedback_rhs(m,l,i_elm_old,i_tor,3) = feedback_rhs(m,l,i_elm_old,i_tor,3) + HZ(i_tor) * v_v
			  feedback_rhs(m,l,i_elm_old,i_tor,4) = feedback_rhs(m,l,i_elm_old,i_tor,4) + HZ(i_tor) * extra_proj !< buiten de steps loop
			enddo

		  enddo
		enddo
	  
	  
      if ((particle_tmp%i_elm .gt. 0) .and. (particle_tmp%weight .gt. 0.d0)) then
        ! Push the particle and determine it's new location.
        call boris_push_cylindrical(particle_tmp, sim%groups(molecules)%mass, E, B, timesteps)

        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                            particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
      end if
   
    end do ! steps 

    call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))

  
  end do   ! particles
  !$omp end parallel do
  
end select
end select 

if (use_ncs) then
    write(*,*) 'GATHER TIME : ',jorek_feedback%rhs_gather_time
    !jorek_feedback%rhs = feedback_rhs / jorek_feedback%rhs_gather_time !* TWOPI
	jorek_feedback%rhs(:,:,:,:,1:3) = feedback_rhs(:,:,:,:,1:3) / jorek_feedback%rhs_gather_time !* TWOPI
	jorek_feedback%rhs(:,:,:,:,4) = feedback_rhs(:,:,:,:,4)
    jorek_feedback%rhs_gather_time = 0.d0
else
    jorek_feedback%rhs = feedback_rhs 
  endif
  
deallocate(feedback_rhs)
!MPI_REDUCE diagnostic. 
call MPI_REDUCE(n_lost_ion, n_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(p_lost_ion, p_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(p_plt_lost, p_plt_lost_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(p_cx_lost, p_cx_lost_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(n_super_ionized, n_super_ionized_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

if (sim%my_id .eq. 0) write(*,'(A46,E14.6,I6)') "Lost molecular superparticles at t due to dissociation: ", sim%time, n_super_ionized_all
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') " Lost molecular particles at t due to diss: ", sim%time, n_lost_ion_all
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "  molecular dissocoation rate at time t [#/s]: ", sim%time, n_lost_ion_all / (timesteps * n_steps)
p_lost_ion_all = p_lost_ion_all / (timesteps * n_steps)
p_plt_lost_all = p_plt_lost_all / (timesteps * n_steps)
p_cx_lost_all = p_cx_lost_all / (timesteps * n_steps)
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Energy exchange to plasma [W] at t due to molecular diss: ", sim%time, p_lost_ion_all ! energy gain
!if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost energy [W] at t due to line radiation: ", sim%time, p_plt_lost_all
!if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Energy exchange to plasma [W] at t due to CX radiation: ", sim%time, p_cx_lost_all
if (sim%my_id .eq. 0) write(*,'(A17,5E14.6)') 'TOTAL Exchange , delta t: ' ,sim%time,p_lost_ion_all, -p_plt_lost_all, p_cx_lost_all, timesteps * n_steps
if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Total energy exchange to plasma [W]: ", sim%time, p_lost_ion_all -p_plt_lost_all+ p_cx_lost_all

! if (sim%my_id .eq. 0) write(*,*) " Lost energy [J] at t due to line radiation: ", sim%time, p_plt_lost_all
!$ w1 = omp_get_wtime()
!$ mmm = mpi_minmeanmax(w1-w0)
!$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"
!if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"

!  write(*,*) 'CAREFUL: averaging over n_steps : ',n_steps
!  jorek_feedback%rhs = jorek_feedback%rhs / real(n_steps,8)
if (sim%my_id .eq. 0) write(*,*) 'done loop_particle_kinetic_local_molecule'

end subroutine ! molecule loop




!AFBLIJVEN

!================================================================================
!                                 RECOMBINATION
!================================================================================
subroutine do_1particle_recombination(element_list,node_list,jorek_stepper,rng,particle_step_time)
use mod_jorek_timestepping !< gives us access to sim?
! use mod_ionisation_recombination, only : rec_rate_local, rec_rate_global, rec_mom_local,rec_energy_local, rec_v_R, rec_v_Z, rec_v_phi
use particle_tracer
!use mod_particle_diagnostics
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_integrate_recomb, only : integrate_recombination
use mod_find_free_particle, only: find_free_particles
!mod_integrate_recomb.f90
implicit none

!class(particle_sim), target, intent(inout)                :: sim
!type(pcg32_rng), dimension(:), allocatable      :: rng
type(pcg32_rng), dimension(:), intent(inout)    :: rng
type(jorek_timestep_action),target           :: jorek_stepper !target
TYPE (type_node_list),         intent(in)     :: node_list
TYPE (type_element_list),      intent(in)     :: element_list
real*8, intent(in)                            :: particle_step_time ! in seconds
	
!internal variables
type (type_element)               :: element
logical, allocatable, dimension(:) :: is_free
integer, allocatable, dimension(:) :: i_free
integer             :: Nrec_part, particles_per_element
real*8              :: total_rec,total_rec_all ,total_volume,total_volume_all
real*8              :: total_Erec_neutral,total_Erec_neutral_all, total_Erec_rad,total_Erec_rad_all
integer             :: n_free,i, k,ielm,ife, i_rng!, element_loc
real*8              :: s, t,R, Z, st_ran(2)

!debug rec
real*8                              :: sanity_rec_local,total_sanity_rec
!rec variables
real*8, dimension(:), allocatable  :: rec_rate_local , rec_v_R, rec_v_Z, rec_v_phi 
real*8, dimension(:), allocatable  :: volume_check, energy_neutrals, energy_radiation  

!Call mod_integrate_recombination
call integrate_recombination(sim%my_id,sim%n_cpu, rec_rate_local, rec_v_R, rec_v_Z, rec_v_phi,volume_check, energy_neutrals, energy_radiation)

sanity_rec_local = 0.d0
!calculate total recombination per mpi proces
total_volume = sum( volume_check(:) )
total_Erec_neutral = sum( energy_neutrals(:) )
total_Erec_rad = sum( energy_radiation(:) )
total_rec = sum( rec_rate_local(:) )
! total recombination
call MPI_REDUCE(total_rec, total_rec_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(total_volume, total_volume_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(total_Erec_neutral, total_Erec_neutral_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(total_Erec_rad, total_Erec_rad_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
if (sim%my_id .eq. 0 .and. sim%time .gt. 0.d0) then
    write(*,'(A30,2E16.8)') 'total recombination weight : ' , sim%time,total_rec_all* central_density* 1.d20 
	write(*,'(A30,2E16.8)') 'Recombination rate  [#/s] : ', sim%time, total_rec_all* central_density* 1.d20 /particle_step_time
	write(*,*) 'total energy to recombined neutrals [J] : ' , total_Erec_neutral_all *1.5d0 / MU_ZERO
	write(*,*) 'total power to recombined neutrals [MW]: ' , total_Erec_neutral_all *1.5d0 / MU_ZERO/particle_step_time /1.d6
	write(*,*) 'total energy lost to Prb [J]: ' , total_Erec_rad_all *1.5d0 / MU_ZERO
	write(*,*) 'total power lost to Prb [MW]: ' , total_Erec_rad_all *1.5d0 / MU_ZERO/particle_step_time /1.d6
	write(*,*) 'total volume : ' , total_volume_all
	write(*,'(A15,6E14.6)') 'TOTAL RECOMB: ',sim%time, total_rec_all* central_density* 1.d20 , total_Erec_neutral_all *1.5d0 / MU_ZERO, total_Erec_neutral_all *1.5d0 / MU_ZERO/particle_step_time /1.d6, &
	                           total_Erec_rad_all *1.5d0 / MU_ZERO, total_Erec_rad_all *1.5d0 / MU_ZERO/particle_step_time /1.d6
endif
!Nrec_part amount of particles needed for this amount of recombination
Nrec_part = int( max(n_particles * 1.d-2 ,total_rec/1.d14 ) )!< assumed average weight per particle (not necesarily the actual weight, as that depends on Srec)
!< limited to 1% of the total initialized particles


!============== Finding free particles 
!> # is_free > n_elements * particles_per_element 
i_free = find_free_particles(sim%groups(atoms)%particles)
!==============================

! loop over all elements
k = 0 !< first free particle
particles_per_element = 1	
!write(*,*) "Doing 1 particle recombination over total n_elemnts", element_list%n_elements
!write(*,*) "Doing 1 particle recombination over n_local_elms", jorek_stepper%n_local_elms
!write(*,*) "Size rec_rate_local", SHAPE(rec_rate_local)
select type (particles => sim%groups(atoms)%particles)
type is (particle_kinetic_leapfrog)
!omp
!$omp parallel do default(shared) &
!$omp schedule(dynamic,10)      &
!$omp shared(sim, particles,jorek_stepper, element_list, node_list, rec_v_R,rec_v_Z,rec_v_phi, &
!$omp i_free,rng,rec_rate_local, &
!$omp CENTRAL_DENSITY, CENTRAL_MASS,sqrt_mu0_over_rho0,particles_per_element ) &
!$omp private(ife,ielm,k,i,element,s,t,R, Z , &
!$omp st_ran, i_rng ) &
!$omp reduction(+:sanity_rec_local)
do ife = 1, size(rec_rate_local) ! loop over all local elements

	if (isnan(rec_v_R(ife)) .or. isnan(rec_v_Z(ife)) .or. isnan(rec_v_phi(ife))) CYCLE !NaN check
	if (rec_rate_local(ife)* central_density* 1.d20 .le. 1.d3) CYCLE
  
  i_rng = 1
	!$ i_rng = omp_get_thread_num()+1
	!if (rec_rate_local(ife) / real(particles_per_element)* central_density* 1.d20 .le. 1.d7) cycle
   
	k = ife !< every OMP thread gets different values
	!< every MPI process has it's own list of i_free.
   
		! --- Get element
	!ielm = jorek_stepper%local_elms(ife) !< actual element number
	ielm    = (sim%my_id+1) + sim%n_cpu*(ife - 1)
	element = element_list%element(ielm)
	
	! initialise particle in the element with Position, Weight, Energy, Momentum			
	do i = 1, particles_per_element
	    k = k *i !< update free particle index ! at begin of loop as k is initialized at k =0
		particles(i_free(k))%weight = rec_rate_local(ife) / real(particles_per_element,8)* central_density* 1.d20 !< rec_rate = in jorek units?
		particles(i_free(k))%i_elm  = ielm  !x, i_elm, st
		particles(i_free(k))%q      = 0
		
		sanity_rec_local = sanity_rec_local + particles(i_free(k))%weight
        !write(31,*) "ielm,",ielm, "k,",k,"i_free(k)",i_free(k) , "particles(i_free(k))%weight,",particles(i_free(k))%weight
		
		!call rng(1)%next(st_ran) !< i_rng should be thread dependent
		call rng(i_rng)%next(st_ran)
		! sample random st combination
		!particles(i_free(k))%st(1:2) = st_ran(2)! [s, t] !< dummi for later
		particles(i_free(k))%st(1) = 0.5d0
		particles(i_free(k))%st(2) = 0.5d0
		
		s = particles(i_free(k))%st(1)
		t = particles(i_free(k))%st(2)
		
		!> uses i_elm and s,t to give us R,Z
		call interp_RZ(node_list,element_list,ielm,s,t,R,Z)
		particles(i_free(k))%x(1:2)  = [R, Z]!  = [R, Z, phi] no phi for axisymmetrix particles
		
		!> distribute directly fluid velocity?
		particles(i_free(k))%v(1)  = rec_v_R(ife)   / (particles(i_free(k))%weight * CENTRAL_MASS * ATOMIC_MASS_UNIT )/ sqrt_mu0_over_rho0 !m/s
		particles(i_free(k))%v(2)  = rec_v_Z(ife)   / (particles(i_free(k))%weight * CENTRAL_MASS * ATOMIC_MASS_UNIT )/ sqrt_mu0_over_rho0
		particles(i_free(k))%v(3)  = rec_v_phi(ife) / (particles(i_free(k))%weight * CENTRAL_MASS * ATOMIC_MASS_UNIT )/ sqrt_mu0_over_rho0
        !< v = momentum fluid lost to recombination / (mass of superparticle)
	end do ! parts_per_element
 
enddo   !ife 
!$omp end parallel do
end select
!end omp

call MPI_REDUCE(sanity_rec_local, total_sanity_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
if (sim%my_id .eq. 0 .and. sim%time .gt. 0.d0) then
  write(*,*) 'SANITY recombination weight : ' , total_sanity_rec 
  write(*,*) 'SANITY Recombination rate  [#/s] : ' , total_sanity_rec /particle_step_time
endif		
!!!!!--------------------------------------------------------------------------
! end subroutine !do_1particle_recombination
end subroutine !do_1particle_recombination

function initialise_sputtering(node_list, element_list, n_reflect) result(D_sputter_source)

  use mod_edge_domain
  use mod_edge_elements

  type(type_node_list), intent(in)    :: node_list
  type(type_element_list)             :: element_list
  type(particle_sputter)              :: D_sputter_source
  integer                             :: n_reflect
  !real*8, allocatable, dimension(:)   :: wall_albedo
  type(type_edge_domain), allocatable, dimension(:) :: edge_domains

  ! number of particles to sputter per species (should be renormalized to yield)

  call find_edge_domains(node_list,element_list, edge_domains)!, discont_corner=.true.)
  if (sim%my_id .eq. 0) write(*,*) "n_domains = ", size(edge_domains,1)
  ! allocate(wall_albedo(size(edge_domains,1)))
  ! wall_albedo(:) = 1.d0 !0.9d0
  ! wall_albedo(size(edge_domains,1)) = 0.1d0
  
  call D_edge%prepare(node_list, element_list, edge_domains, nsub=6, nsub_toroidal=1)!,wall_albedo=wall_albedo)

  ! target group, number of particles per mpi task, densities, Zs, basename
  D_sputter_source = particle_sputter(D_edge, 1, n_reflect, basename='D_reflect')
  D_sputter_source%use_Yn_func = .false.
  D_sputter_source%n_save = nout !10 ! or nout
  D_sputter_source%albedo_for_neutrals = 1.d0
  D_sputter_source%sputtered_particle_weight_threshold = 1.d0
  
  !

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

pure function f_psi_inside(n, P, grad_P) result(f) 
!deeltjes te sampelen met een probability
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*4 :: f, psi_norm
  f = 5e-2 !sampling probability 
  
  psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)
  ! if (P(1) .lt. -0.26 .and. P(1) .gt. -0.35) f = 1e0
  if (psi_norm .lt. 0.9d0 .and. psi_norm .gt. 0.d0) f = 1e0
end function f_psi_inside

end program recombination_loop
