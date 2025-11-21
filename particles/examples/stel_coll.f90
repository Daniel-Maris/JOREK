!>#Example 1: single particle in JOREK field
!> This example follows a particle in a static JOREK equilibrium field
program ex1_stel
use particle_tracer
!use mod_interp, only: interp_gvec
use mod_particle_io
use mod_interp, only: mode_moivre, interp, interp_RZ, interp_RZP
use mod_basisfunctions
use mod_jorek_timestepping
use omp_lib
use mpi
use mod_pcg32_rng, only: pcg32_rng
use mod_ionisation_recombination, only: new_charge
use mod_collisions
!use mod_project_particles
use constants
use phys_module, only: restart, restart_particles
use phys_module, only: tstep_particles, nstep_particles, nstep_n
use phys_module, only: use_manual_random_seed, nout, index_start
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0

use mod_random_seed
use mod_particle_allocation
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff


implicit none

interface
  subroutine find_RZP(node_list,element_list,R_find,Z_find,phi_find,R_out,Z_out,ielm_out,s_out,t_out,ifail,checked_elms)
    use data_structure
    type (type_node_list),    intent(in)     :: node_list
    type (type_element_list), intent(in)     :: element_list
    real*8,                   intent(in)     :: R_find, Z_find, phi_find
    real*8,                   intent(out)    :: R_out,Z_out,s_out,t_out
    integer,                  intent(inout)  :: ielm_out
    integer,                  intent(out)    :: ifail,checked_elms
  end subroutine find_RZP
end interface

! --- Set up the simulation variables containing
!     sim: particles, time, and io.
!     events: halting points for the pushers and actions to run.
class(*), pointer :: p
integer :: i, j, k, n_steps, n_lost, istep, ierr
real*8  :: target_time, time
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U
type(event) :: fieldreader, partreader
character(len=1024) :: filename, output_filename
integer :: i_elm, ifail, i_elm_old, q_old, checked_elms
real*8  :: dummy, DUMMY_R, DUMMY_Z, s, t
real*8 :: current_time, start_time, time_per_step, eta
real*8, parameter  :: R=2.044, Z=0.0, phi=0.0, sigma_r=.0004d0, sigma_z=.004d0, A=1.d0, sigma_phi=PI/32.d0
real*8, parameter :: vR=0.d0, vZ=0.d0, vphi=0.d0
real*8, dimension(:,:), allocatable :: positions ! For gaussian ditributed markers

! Coronal equilibrium
real*8 :: n_e, T_e
real*8    :: ionize_ran_imp(2)
type(pcg32_rng), dimension(:), allocatable :: rng
integer   :: n_stream, l, seed_size, i_rng, seed

! To restart particles
real*8 :: n_norm, rho_norm, t_norm
real*8 :: tstep, tstep_fluid_si, tstep_part_adj
character (len=29) :: hdf5_file_name
! w MPI logic
integer, dimension(:), allocatable :: n_particles_per_mpi, global_start_index
integer :: current_offset

! Collisions
real*8, dimension(1) :: P_, P_s, P_t, P_phi, P_time
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, R_, Z_, R_phi, Z_phi
logical :: use_coll
integer, parameter :: n_coll = 5     ! 5
integer*1 :: q_b
real*8 :: n_b, m_b, kTb, coulomb_log
real*8 :: grad_T_e(3), q(3), ran2(6, n_coll), v_b(3, n_coll), ran(6)
real*8 :: Du,v(3),nv

use_coll = .true.

!***********************************************************************
!*                           initialisation                            *
!***********************************************************************

! Start up MPI, jorek
call sim%initialize()

! Loading the jorek fields
if (restart) then
    fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
    call with(sim, fieldreader)
  else
    if (sim%my_id == 0) write(*,*) 'ERROR: using this program without restarting from a jorek field is not possible. Please set restart=.t. in the namelist and provide a jorek_restart.h5 file'
    stop
end if

! --------- Set up the particles ---------
if (restart_particles) then
  write(*,*) "--------------------------- restarting particles ---------------------------"

  ! Read particles from a file
  if (sim%my_id == 0) write(*,*) '========== INFO: READING PARTICLES RESTART FILE =========='
  partreader = event(read_action(filename='part_restart.h5'))
  call with(sim, partreader) !<defines sim%groups and the corresponding particles
else
  call allocate_particles_for_sim(sim) ! populate the particle arrays in the particle groups

  ! Generate Gaussian distributed blob of particles
  allocate(positions(3, int(sim%groups(1)%n_particles)))
  if (sim%my_id == 0) positions = generate_3d_gaussian(int(sim%groups(1)%n_particles), R, Z, phi, sigma_r, sigma_z, sigma_phi, A)

  !write(*,*) "sim%my_id n_particles", sim%my_id, size(sim%groups(1)%particles)

  ! Broadcast (share) positions to all MPI threads
  call MPI_BCAST(positions, 3 * int(sim%groups(1)%n_particles), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
    write(*,*) 'Error in MPI_BCAST on rank', sim%my_id
  end if

  ! Retrieving ADAS data
  call read_adf11_test(sim%groups(1)%ad, 0,'50_w')

  n_particles_per_mpi = calc_n_particles_per_mpi_array(int(sim%groups(1)%n_particles), sim%n_mpi)
  !if(sim%my_id .eq. 0) write(*,*) "n_particles_per_mpi", n_particles_per_mpi

  allocate(global_start_index(sim%n_mpi))
  current_offset = 1
  do j=1, sim%n_mpi
    global_start_index(j) = current_offset
    !write(*,*) "j, current_offset, global_start_index", j, current_offset, global_start_index(j)
    current_offset = current_offset + n_particles_per_mpi(j)
  end do

  ! Get particle location in element, s, t, space (parallelized with MPI and OpenMP)
  !$omp parallel do default(none) shared(sim, positions, global_start_index, n_particles_per_mpi) private(p, DUMMY_R, DUMMY_Z, i_elm, s, t, ifail, j, E, B, checked_elms)
  do j = 1, n_particles_per_mpi(sim%my_id + 1)  ! each MPI thread allocates a position for its local number of particles
    select type (p => sim%groups(1)%particles(j))
    type is (particle_kinetic_leapfrog)
      call find_RZP(sim%fields%node_list, sim%fields%element_list,  &
      !              positions(1,j), positions(2,j), positions(3,j), &
                    positions(1, global_start_index(sim%my_id + 1) + (j - 1)), & 
                    positions(2, global_start_index(sim%my_id + 1) + (j - 1)), &
                    positions(3, global_start_index(sim%my_id + 1) + (j - 1)), &
                    DUMMY_R, DUMMY_Z, i_elm, s, t, ifail, checked_elms)
      !p%x      = [R,Z,phi]
      !p%x      = positions(:, j)
      p%x      = positions(:, global_start_index(sim%my_id + 1) + (j - 1))
      p%i_elm  = i_elm
      p%st     = [s,t]
      p%v      = [vR, vZ, vphi]
      p%q      = 0
      p%weight = 1.0d15
      p%T_e    = 0.
      p%n_e    = 0.
      if (j==1) write(*,*) "sim%my_id, p%x", sim%my_id, p%x ! if they are the same the rng has same seed for different mpis
    type is (particle_fieldline)
      p%x = [R,Z,phi]
      p%i_elm = i_elm
      p%st = [s, t]
      E = 0.0
      B = 0.0
    end select
  end do
  !$omp end parallel do

  deallocate(positions)
  deallocate(global_start_index)
  deallocate(n_particles_per_mpi)

  if(sim%my_id .eq. 0) write(*,*) "--- particles initialized"
end if


! --- Set dpsi_dt to be zero as treating JOREK as an equilibrium field
call sim%fields%set_flag_dpsidt(.true.)


! --- Setting up random numbers for ionisation probability
seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
write(*,*) "id, n_mpi, n_stream",sim%my_id, sim%n_mpi, n_stream
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(1, seed, n_stream, i)
end do



!***********************************************************************
!*                         setting up physics                          *
!***********************************************************************

! --- Calculating normalisation constants
n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek 


!***********************************************************************
!*                           main loop                                 *
!***********************************************************************

! --- Loop until the simulation requests a stop
if(sim%my_id .eq. 0) write(*,*) "--------------------------- Start tracing ---------------------------"
istep = 0
n_lost = 0
do while (.not. sim%stop_now)
  istep = istep + 1
  if(sim%my_id .eq. 0) write(*,'(A80)'  ) "================================================================================"
  if(sim%my_id .eq. 0) write(*,'(A37,I6)') "Starting main loop iteration istep = ",istep
  if(sim%my_id .eq. 0) write(*,'(A80)'  ) "================================================================================"

  ! --- Determining the time stepping for this fluid step
  tstep = get_tstep_n(istep)            ! fluid dt in JOREK units
  tstep_fluid_si = tstep*t_norm         ! fuild dt in SI units
  sim%time = sim%time + tstep_fluid_si  ! carries the time at the end of the current step

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

  ! --- Loop over all particle groups
  do i=1,1
    if(use_manual_random_seed) then
      !$ call omp_set_schedule(omp_sched_static,10)
    else
      !$ call omp_set_schedule(omp_sched_dynamic,10)
    end if
    !$omp parallel do schedule(runtime) default(none) shared(sim, n_steps, istep, i, rng, tstep_particles, nstep_particles, tstep_part_adj, use_coll, v, nv, Du) &
    !$omp private(j, k, t, p, i_elm_old, q_old, E, B, psi, U, rz_old, st_old, n_e, T_e, ionize_ran_imp, ifail, filename, i_rng, l, &
    !$omp grad_T_e, P_, P_s, P_t, P_phi, P_time, R_, Z_, R_s, R_t, Z_s, Z_t, R_phi, Z_phi, kTb, n_b, q_b, m_b, q, coulomb_log, ran2, v_b, ran) &
    !$omp reduction(+:n_lost)
    do j=1,size(sim%groups(i)%particles)
      !$ i_rng = omp_get_thread_num()+1
      do k=1, nstep_particles

        t = sim%time + (k-1)*tstep_part_adj

        ! --- Integrating the particle with corresponding time evolution scheme
        select type (p => sim%groups(i)%particles(j))
          type is (particle_kinetic_leapfrog)
            if (p%i_elm .gt. 0) then
              call sim%fields%calc_EBpsiU(t, p%i_elm, p%st, p%x(3), E, B, psi, U)
              rz_old    = p%x(1:2)
              st_old    = p%st
              i_elm_old = p%i_elm
              q_old     = p%q

              !> exit evolution loop if particle is outside domain
              if (p%i_elm .le. 0) exit

              !> check that particle weight is non negative
              if (p%weight .lt. 0.0d0) write(*,*) "Negative particle weight p(j)%w=", p%weight

              !> calculate n_i and T_e (T_e in K, n_e in m-3) (jorek model assumption: n_e = n_i)
              call sim%fields%calc_NeTe(t, p%i_elm, p%st, p%x(3), n_e=n_e, T_e=T_e, grad_T_e=grad_T_e)  ! &
              ! always use grad_T_e, otherwise T2_1 is NaN matrix when sampling shifted Maxwellian
              p%n_e = n_e
              p%T_e = T_e * K_BOLTZ / EL_CHG !< T_e saved in eV

              !> IONISATION AND RECOMBINATION
              call rng(i_rng)%next(ionize_ran_imp)
              p%q = int(new_charge(int(q_old,4), sim%groups(i)%ad, log10(n_e), log10(T_e), tstep_part_adj, ionize_ran_imp(1:2)),1)
              !< IONISATION AND RECOMBINATION
              
              !> COLLISIONS
              if (use_coll .eq. .true.) then
                if (p%q .gt. 0) then
                  ! variables
                  kTb = T_e * K_BOLTZ
                  n_b = n_e
                  q_b = 1
                  m_b = 2.
                
                  ! double-check this!
                  q = q_homma2013(kTb, grad_T_e*K_BOLTZ, B, n_b, m_b, q_b)
                  
                  ! double-check this!
                  coulomb_log = coulomb_logarithm(kTb, n_b, p%q, q_b, sim%groups(1)%mass, m_b)
                  coulomb_log = max(10.d0, coulomb_log)
                  coulomb_log = min(20.d0, coulomb_log)
                
                  ! double-check this!
                  call sim%fields%interp_PRZP_1(t, p%i_elm, [var_Vpar], 1, p%st(1), p%st(2), p%x(3), &
                                                P_, P_s, P_t, P_phi, P_time, R_, R_s, R_t, R_phi, Z_, Z_s, Z_t, Z_phi)
                
                  do l=1,n_coll
                    call rng(i_rng)%next(ran2(:,l))
                  end do
                  
                  ! double-check this!
                  !write(*,*) "P_(1)*B/norm2(B)/sim%t_norm", P_(1)*B/norm2(B)/sim%t_norm
                  !write(*,*) "P_(1)*B/sim%t_norm", P_(1)*B/sim%t_norm
                  call sample_velocity_dist_magnetized(n_coll, ran2(1:6,:), kTb, q, n_b, m_b, q_b, P_(1)*B/norm2(B)/sim%t_norm, v_b) ! <- why different in mod_particle_evolution?

                  do l=1,n_coll
                    call rng(i_rng)%next(ran)

                    !> Check particle timestepping
                    !v = p%v - v_b(:,l)
                    !nv = norm2(v)
                    !Du = real(int(p%q,4)**2 * int(q_b,4)**2,8) * EL_CHG**4 * n_b * coulomb_log / &
                    !    (8.d0 * PI * EPS_ZERO**2 * (sim%groups(1)%mass*m_b/(sim%groups(1)%mass+m_b))**2 * ATOMIC_MASS_UNIT**2 * (nv**3))

                    !if (sim%my_id == 0) write(*,*) "Du, dt, Du*dt", Du, tstep_part_adj/real(n_coll,8), Du*tstep_part_adj/real(n_coll,8)
                    !< Check particle timestepping

                    ! double-check this!
                    call collide_particles(ran(1:3), p%q, sim%groups(1)%mass, p%v, &
                            q_b, m_b, v_b(:,l), n_b, coulomb_log, tstep_part_adj/real(n_coll,8))
                  end do
                end if
              end if !< COLLISIONS

              !> Push particles and find out where they are next
              call boris_push_cylindrical(p, m=sim%groups(i)%mass, E=E, B=B, dt=tstep_part_adj) ! timesteps(i)
              call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, p%x(1), p%x(2), p%st(1), p%st(2), p%i_elm, ifail, p%x(3))
              !call interp_gvec(sim%fields%node_list, sim%fields%element_list, p%i_elm, 4, 1, 1, p%st(1), p%st(2), psi_norm, dummy, dummy, dummy, dummy, dummy) !s_norm, BRg_s, BRg_t, BRg_st, BRg_ss, BRg_tt
              !psi_norm = psi_norm * psi_norm
            end if
          type is (particle_fieldline)
            call field_line_runge_kutta_fixed_dt_push_jorek(sim%fields, p, sim%time, t)
            !call interp_gvec(sim%fields%node_list, sim%fields%element_list, p%i_elm, 4, 1, 1, p%st(1), p%st(2), s_norm, dummy, dummy, dummy, dummy, dummy)
        end select

        if (sim%groups(i)%particles(j)%i_elm .le. 0) then
          n_lost = n_lost + 1
          write(*,*) "Particle ", j, " lost at step ", k
          exit
        endif
      end do ! steps
    end do ! particles
    !$omp end parallel do
    if (sim%my_id == 0) write(*,*) "number/% of lost particles: ", n_lost, n_lost/sim%groups(1)%n_particles*100.0
  end do ! groups
  
  ! --- Write restart files
  if (mod(istep, nout) .eq. 0) then
    write(*,*) 'Writing particle restart file'
    write(hdf5_file_name, '(A, I0.7, A)') 'coll_1e5_jorek_part', index_start + istep, '.h5'
    call write_simulation_hdf5(sim, hdf5_file_name)
  end if

  ! --- Finish loop according to nstep_n variable, which can be assigned in namelist
  if (istep .ge. nstep_n(1)) then
    sim%stop_now = .true.
    call write_simulation_hdf5(sim, 'part_restart.h5') ! Write last particle restart file as part_restart.h5
  end if

end do ! while

deallocate(rng)

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

  write(*,'(A80)') "================================================================================"
  write(*,*) what
  write(*,'(A80)') "================================================================================"

end subroutine

!> Generates positions following a 3D Gaussian distribution
!! Creates an array of points with x, y, and z coordinates sampled from
!! independent Gaussian distributions.
!!
!! @param n_points Number of positions to generate
!! @param x0       Mean of the x-coordinate Gaussian distribution
!! @param y0       Mean of the y-coordinate Gaussian distribution
!! @param z0       Mean of the z-coordinate Gaussian distribution
!! @param sigma_x  Standard deviation for the x-coordinate
!! @param param sigma_y Standard deviation for the y-coordinate
!! @param sigma_z  Standard deviation for the z-coordinate
!! @param A        Amplitude parameter (not used in current implementation)
!! @return positions Array containing [x, y, z] coordinates for each point
function generate_3d_gaussian(n_points, x0, y0, z0, sigma_x, sigma_y, sigma_z, A) result(positions)
  use constants, only: pi
  implicit none

  ! Argument declarations
  integer, intent(in) :: n_points
  real(8), intent(in) :: x0, y0, z0, sigma_x, sigma_y, sigma_z, A
  real(8), dimension(3,n_points) :: positions

  ! Local variable declarations
  integer :: i
  real(8) :: u1, u2, r_gauss, theta_gauss, norm1, norm2
  real(8) :: u3, u4, r2_gauss, theta2_gauss, norm3

  do i = 1, n_points
  ! Generate two independent standard normals for x and y using Box–Muller transform
    call random_number(u1)
    call random_number(u2)
    if(u1 == 0.0d0) u1 = 1.0d-10  ! Avoid log(0)
    r_gauss = sqrt(-2.0d0 * log(u1))
    theta_gauss = 2.0d0 * pi * u2
    norm1 = r_gauss * cos(theta_gauss)
    norm2 = r_gauss * sin(theta_gauss)
      
    positions(1,i) = x0 + sigma_x * norm1
    positions(2,i) = y0 + sigma_y * norm2

    ! Generate one standard normal for z using another Box–Muller transform
    call random_number(u3)
    call random_number(u4)
    if(u3 == 0.0d0) u3 = 1.0d-10  ! Avoid log(0)
    r2_gauss = sqrt(-2.0d0 * log(u3))
    theta2_gauss = 2.0d0 * pi * u4
    norm3 = r2_gauss * cos(theta2_gauss)
      
    positions(3,i) = z0 + sigma_z * norm3
  end do
end function generate_3d_gaussian

end program ex1_stel
              
