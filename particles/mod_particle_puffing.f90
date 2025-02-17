!> Module for puffing gas into the plasma
!> This module will create new particles at the locations where gas valves will be.

module mod_particle_puffing
  use mod_edge_elements
  use mod_io_actions, only: io_action
  use mod_sampling
  use mod_particle_types
  use constants, only: TWOPI, K_BOLTZ, ATOMIC_MASS_UNIT
  use mod_rng, only: type_rng, setup_shared_rngs
  use mod_boundary, only: wall_normal_vector
  use mod_atomic_elements 
  use mod_particle_sim
  use mod_event
  use mod_find_rz_nearby, only: find_rz_nearby
  use phys_module, only: type_valve, valves, part_group_configs, type_puff_ctrl, n_puff_segment_max 

  implicit none

  private
  public  :: particle_puffing

  ! Extend type
  type, extends(io_action) :: particle_puffing
    
    class(type_rng), dimension(:), allocatable :: rng    !< one RNG per openmp thread
    !> Valve
    type(type_valve) :: puff_valve      !< determines the location of the puffing
    integer          :: valve_num       !< the index of the valve used for this puffing event
    !> Puff ctrl 
    type(type_puff_ctrl) :: puff_ctrl   !< the controller for this puffing action
    !> Variables required for piecewise linear time dependent puffing control 
    integer          :: current_puff_seg = 0   
    integer          :: last_puff_seg = n_puff_segment_max !< the segment number after the last defined puffing keyframe
    !> Target particle group
    integer          :: target_group
    
  contains
    procedure :: do => do_particle_puffing
  end type particle_puffing

  interface particle_puffing
    module procedure new_particle_puffing
  end interface particle_puffing
contains

function new_particle_puffing(sim, target_group, valve_num, rng, seed) result(new)

  use mod_pcg32_rng,   only: pcg32_rng
  use mod_random_seed, only: random_seed
  
  type(particle_puffing)    :: new

  type(particle_sim), intent(in)  :: sim           
  integer, intent(in)             :: target_group
  integer, intent(in)             :: valve_num           !< the valve number to use for this puffing
    
  class(type_rng), intent(in), optional :: rng !< random number generator to use (deafult PCG32)
  integer, intent(in), optional         :: seed !< Seed for the RNG (default random_seed() on my_id + bcast)
  integer                               :: my_seed, i
  
  new%target_group     = target_group
  new%puff_valve       = valves(valve_num)
  new%valve_num        = valve_num
  new%puff_ctrl        = part_group_configs(target_group)%puff_ctrl(valve_num)

  ! determine the current puffing segment and the last puffing segment
  do i=1, n_puff_segment_max-1
    !> find last puffing segment
    if ((new%puff_ctrl%times(i) /= -1) .and. (new%puff_ctrl%times(i+1) == -1)) new%last_puff_seg = i
    !> find current puffing segment
    if (sim%time > new%puff_ctrl%times(i) .and. (new%puff_ctrl%times(i) /= -1)) new%current_puff_seg = i
  enddo

  !> allocate random seed for sampling
  if (present(seed)) then
    my_seed = seed
  else
    my_seed = random_seed()
  end if
  if (present(rng)) then
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=rng, rngs=new%rng)
  else
    ! default to pcg32_rng for reflection
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=pcg32_rng(), rngs=new%rng)
  end if

end function new_particle_puffing

!> Actually puff gass
subroutine do_particle_puffing(this,sim, ev)
  use mpi_mod
  use phys_module, only: tstep, central_mass, central_density
  use constants, only: MASS_PROTON, MU_ZERO
  ! !$ use omp_lib

  class(particle_puffing) , intent(inout) :: this
  type(particle_sim), intent(inout)       :: sim
  type(event), intent(inout), optional    :: ev 

  integer :: ierr,i_scalar, n_free, j, k, n_group, i_elm, i_elm_new, ifail, i_p, to_puff, supers_per_puff_local, i_rng
  logical, allocatable, dimension(:) :: is_free
  integer, allocatable, dimension(:) :: i_free
  real*8  :: tstep_fluid_si, c, R, Z, phi, s, t
  real*8  :: R_new, Z_new, s_new, t_new, r_valve, theta
  real*8  :: vector_normal(3), u(5)
  real*8  :: puff_rate !< possibly time dependent fueling rate

  ! variables for piecewise linear time dependent puffing
  real*8  :: puff_rate_0, puff_rate_1, puff_time_0, puff_time_1

  integer ::    puffed_this_step_local, all_puffed_this_step
  real*8  ::    puff_weight_local, all_puff_weight

  tstep_fluid_si = tstep*sqrt((MU_ZERO * CENTRAL_MASS * MASS_PROTON * CENTRAL_DENSITY * 1.d20))

  if (sim%my_id .eq. 0) write(*,'(A,A,A,I1,A)') "--- Started puffing for Group: ", sim%groups(this%target_group)%id, ", Valve: ", this%valve_num, " ---"
  
  if (this%puff_ctrl%supers_num_puff .le. -1.d-6) then ! 0.d0
    if (sim%my_id .eq. 0) write(*,*) "ERROR: No puff quota set, i.e. puff_ctrl%supers_num_puff = 0 for this group"  
    stop
  end if
  
  supers_per_puff_local = this%puff_ctrl%supers_num_puff / sim%n_mpi !supers_per_puff_local is the amount of superparticles that will be puffed per MPI process.

  !> check if the simulation has advanced to the next puffing segment
  if (this%current_puff_seg /= this%last_puff_seg) then
    if (sim%time > this%puff_ctrl%times(this%current_puff_seg + 1)) this%current_puff_seg = this%current_puff_seg + 1 
  endif

  !> set the bounding values of the puffing segment
  if (this%current_puff_seg == 0) then 
    !> for t < puff_ctrl%times(1), we keep puff_rate constant at puff_ctrl%rates(1)
    puff_rate_0  = this%puff_ctrl%rates(1)
    puff_rate_1  = this%puff_ctrl%rates(1) 

    !> the puff_ctrl%times no longer matter in this case but (puff_ctrl%times_1 - puff_ctrl%times_0) has to be non zero 
    puff_time_0 = this%puff_ctrl%times(1) 
    puff_time_1 = this%puff_ctrl%times(1) + 1 
  else if (this%current_puff_seg == this%last_puff_seg) then
    !> for t > puff_ctrl%times(last_puff_seg), we keep puff_rate constant at puff_ctrl%rates(last_puff_seg)
    puff_rate_0  = this%puff_ctrl%rates(this%last_puff_seg)
    puff_rate_1  = this%puff_ctrl%rates(this%last_puff_seg) 

    !> the puff_ctrl%times no longer matter in this case but (puff_ctrl%times_1 - puff_ctrl%times_0) has to be non zero 
    puff_time_0 = this%puff_ctrl%times(this%last_puff_seg) 
    puff_time_1 = this%puff_ctrl%times(this%last_puff_seg) + 1 
  else
    !> in general, puff_ctrl%times and puff_ctrl%rates are the defined values bounding the segment
    puff_rate_0  = this%puff_ctrl%rates(this%current_puff_seg)
    puff_rate_1  = this%puff_ctrl%rates(this%current_puff_seg + 1)
    puff_time_0 = this%puff_ctrl%times(this%current_puff_seg)
    puff_time_1 = this%puff_ctrl%times(this%current_puff_seg + 1)
  endif

!============== Finding free particles !< make into a function?
allocate(is_free(size(sim%groups(this%target_group)%particles,1))) 
!$omp parallel do default(none) shared(sim, this, n_free, i_free, is_free) &
!$omp private(j) schedule(dynamic, 100)
do j=1,size(sim%groups(this%target_group)%particles,1) !sim%groups(this%target_group)%particles
  is_free(j) = sim%groups(this%target_group)%particles(j)%i_elm .le. 0  !< array T/F is particle is free
end do
!$omp end parallel do
!$omp barrier
n_free = count(is_free)
allocate(i_free(n_free))
k = 1
do j=1,size(is_free,1)
  if (is_free(j)) then
    i_free(k) = j !< i_free(k) has index of free particle in  sim%groups(this%target_group)%particles(j)
    k = k+1
    !if (sim%my_id .eq. 0) write(*,*) "Adding to the list number: ", j
  end if
end do
! ==================
  
  
  ! the current set up may only work for puffing Hydrogen
  ! Assuming the incoming gas at T=300C and a diatomic gas
  c = sqrt((7.d0/5.d0)*(300.d0+273.d0)*K_BOLTZ/(2.d0*sim%groups(this%target_group)%mass*ATOMIC_MASS_UNIT))
  if (this%puff_valve%type == "circ") then
    call find_RZ(sim%fields%node_list, sim%fields%element_list, this%puff_valve%R_valve_loc, this%puff_valve%Z_valve_loc, R, Z, &
           i_elm, s, t ,ifail)
           
  else
    call find_RZ(sim%fields%node_list, sim%fields%element_list, sum(this%puff_valve%poly_R(1:2))/2.d0, sum(this%puff_valve%poly_Z(1:2))/2.d0, R, Z, &
           i_elm, s, t ,ifail)
  endif
  if (ifail .ne. 0) then
    if (sim%my_id .eq. 0) write(*,*) "Warning: The valve location for puffing could not be found, maybe it was placed outside of the grid?"
    stop
  end if

  vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, &
          i_elm, s, t)

!------------- Decide how many superparticles to initiate 

  to_puff = supers_per_puff_local

  !> calculate time dependent puffing rate
  puff_rate = calc_puff_rate_linear(sim%time, puff_rate_0, puff_rate_1, puff_time_0, puff_time_1)

  !> output time dependent puffing details
  if (sim%my_id .eq. 0) then
    write(*,"(A,G12.6,A)")           "Puffing details for time t=", sim%time, ":"
    write(*,"(2X,A12)") "Set-up:     "
    write(*,"(4X,A15, ' = ', I12)")        "puff segment   ", this%current_puff_seg
    write(*,"(4X,A15, ' = ', G12.6)")      "puff_rate_0    " , puff_rate_0
    write(*,"(4X,A15, ' = ', G12.6)")      "puff_rate_1    " , puff_rate_1
    write(*,"(4X,A15, ' = ', G12.6,A)")    "puff_rate      "   , puff_rate, " atoms/s"
    write(*,*) ""
  endif
      
  if (to_puff .ge. n_free) then
    write(*,*) "Warning could not puff the requested amount."
    to_puff = n_free
  end if

!-------------  
  
  
  puffed_this_step_local = 0
  puff_weight_local      = 0.d0
  select type (pa => sim%groups(this%target_group)%particles)
  type is (particle_kinetic_leapfrog)
  !> To do: Parallellize loop 
  !> This loop initializes the to be puffed particles and places then in the computational domain.
  !> It counts the amount of marker particles and weight that was initialized.
  !> This loop can be OMP parallel. there was a small bug in the OMP implementation.
  !> The almost done loop is written below.
  !> reduction of puffed_this_step_local,puff_weight_local for diagnostics
  !>
  ! #ifdef __GFORTRAN__
  !  !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
  ! #else
  ! !$omp parallel do default(shared) &
  ! !$omp schedule(dynamic,10) &
  ! !$omp shared(sim, pa, this,i_free,c, vector_normal,                       &
  ! !$omp   to_puff,supers_per_puff_local, tstep_fluid_si,puff_rate )                        &
  ! !$omp private(i_p, i_rng, j,k,u , R,Z,s,t,R_new,Z_new,s_new,t_new,     &
  ! !$omp  i_elm,i_elm_new,r_valve, theta,                                         &
  ! !$omp ifail)                                                                    &
  ! !$omp reduction(+:puffed_this_step_local,puff_weight_local)
    do j = 1, to_puff
      i_p = i_free(j)
      do 
    
        !      !$ i_rng = omp_get_thread_num()+1
        call this%rng(1)%next(u) !rng(1)
        if (this%puff_valve%type == "circ") then
          r_valve = this%puff_valve%r_valve*sample_piecewise_linear(2, [0.d0, 1.d0], [1.d0, 0.d0], u(1))
          theta = TWOPI * u(2)
          R_new = this%puff_valve%R_valve_loc + r_valve * cos(theta)
          Z_new = this%puff_valve%Z_valve_loc + r_valve * sin(theta)
        else
          s = u(1)
          t = u(2)
          R_new = this%puff_valve%poly_R(1)*(1.d0-s)*(1.d0-t) + this%puff_valve%poly_R(2)*s*(1.d0-t) + this%puff_valve%poly_R(3)*(1.d0-s)*t + this%puff_valve%poly_R(4)*s*t
          Z_new = this%puff_valve%poly_Z(1)*(1.d0-s)*(1.d0-t) + this%puff_valve%poly_Z(2)*s*(1.d0-t) + this%puff_valve%poly_Z(3)*(1.d0-s)*t + this%puff_valve%poly_Z(4)*s*t
        endif
    
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, R, Z, s, t, i_elm, &
        R_new, Z_new, s_new, t_new, i_elm_new, ifail)
        if (ifail .ge. 0) exit
      end do
      R     = R_new
      Z     = Z_new
      s     = s_new
      t     = t_new
      i_elm = i_elm_new
      if (this%puff_valve%phi .lt. 0.d0) then
        pa(i_p)%x(3) = TWOPI*u(3)
      else 
        pa(i_p)%x(3) = this%puff_valve%phi 
      end if
      pa(i_p)%x(1:2)  = [R, Z]
      pa(i_p)%st(1:2) = [s, t]
      pa(i_p)%i_elm   = i_elm
      pa(i_p)%weight  = real(1.d0/supers_per_puff_local) * tstep_fluid_si * puff_rate
   
      pa(i_p)%v       = c * sample_cosine(u(4:5), vector_normal)   
      pa(i_p)%q       = 0_1
      if (sim%groups(this%target_group)%particles(i_p)%weight  .le. 1.d-2) then ! if the weight is too low. 
        sim%groups(this%target_group)%particles(i_p)%i_elm = 0
        cycle       
      end if
    puffed_this_step_local = puffed_this_step_local+1
    puff_weight_local      = puff_weight_local + pa(i_p)%weight 
    end do
  ! !$omp end parallel do  
  class default
    write(*,*) 'Particle type not implemented for gas fueling.'
    stop
  end select

! puffed_this_step_local
call MPI_REDUCE(puffed_this_step_local,all_puffed_this_step,1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)   
call MPI_REDUCE(puff_weight_local,all_puff_weight,1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)  
if (sim%my_id .eq. 0) then
  write(*,"(2X,A12)") "Puffed:     "
  write(*,"(4X,A15, ' = ', I12)")    "Superparticles ", all_puffed_this_step
  write(*,"(4X,A15, ' = ', E12.6)")  "Total Weight   ", all_puff_weight
  write(*,"(A)") "---------------------------------------------"
endif
  
  
end subroutine do_particle_puffing

!< linearly interpolates the puff rate between two determined points (t0, y0), (t1, y1)
!< y = y0 + [(y1 - y0)/(t1 - t0)] * (t - t0)
function calc_puff_rate_linear(time, puff_rate_0, puff_rate_1, puff_time_0, puff_time_1) result(puff_rate)

  implicit none
  real*8, intent(in)    :: time          !< current simulation time (t)
  real*8, intent(in)    :: puff_rate_0   !< the left set value of the linear function (y0)
  real*8, intent(in)    :: puff_rate_1   !< the right set value of the linear function (y1)
  real*8, intent(in)    :: puff_time_0   !< t0
  real*8, intent(in)    :: puff_time_1   !< t1
  real*8                :: puff_rate     !< the puff rate at time t (y)

  puff_rate = puff_rate_0 + (puff_rate_1 - puff_rate_0) / (puff_time_1 - puff_time_0) * (time - puff_time_0)

end function calc_puff_rate_linear

end module