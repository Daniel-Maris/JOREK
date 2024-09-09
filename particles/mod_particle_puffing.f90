!> Module for puffing gas into the plasma
!> This module will create new particles at the locations where gas valves will be.

module mod_particle_puffing
  use mod_edge_elements
  use mod_io_actions, only: io_action
  use mod_sampling
  use mod_particle_types
  use constants, only: TWOPI, K_BOLTZ, ATOMIC_MASS_UNIT!, MU_ZERO
  ! use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
  use mod_rng, only: type_rng, setup_shared_rngs
  use mod_boundary, only: wall_normal_vector
  use mod_atomic_elements !mod_elements !< Chemical elements
  use mod_particle_sim
  use mod_event
  use mod_find_rz_nearby, only: find_rz_nearby


  

  implicit none

  private
  public  :: particle_puffing

  ! Extend type
  type, extends(io_action) :: particle_puffing
   
    class(type_rng), dimension(:), allocatable :: rng  !< one RNG per openmp thread
   
    ! number of simulation particles/s to puff across all processes
    integer :: n_puff = -1 
    ! Average fueling rate: 9.7d22; max fueling rate 18d22
    real*8  :: fueling_rate = -1.d0
    real*8  :: R = -1.d0, Z = -1.d0, phi = -1.d0
    real*8  :: valve_r = -1.d0  !< radius of gas valve
    real*8  :: last_time = 0.d0 !< When did we puff last 
    real*8  :: last_diag_time = 0.d0 !< Last time of output of diagnostics
	
    !> Should maybe go into a shape function?
    !box volume puff, define 4 RZ points to determine volume
    real*8  :: poly_R(4) = -1.d0
    real*8	:: poly_Z(4) = -1.d0
    logical :: boxpuff = .false.
    
    
    !Time dependent puffing
    logical :: puff_t_dependent = .false. !< flat ramp flat
    real*8  :: fueling_rate_start = 0.d0
    real*8  :: t_puff_start = 0.d0 !< defined in JOREK time units
    real*8  :: t_puff_slope = 0.d0 !<defined in SI
    !Getting the puffingrate from a file
    logical                          :: rate_from_file = .false.
    character(len=64)                :: puff_rate_filename = 'none'
    real*8,dimension(:), allocatable :: t_pts         !< ascending list with time in (s)
    real*8,dimension(:), allocatable :: puff_rate_pts !< list with corresponding puffrates at those times
    integer                          :: it_last = 2   !< index of t_pts at last puffing event
    
  contains
    procedure :: do => do_particle_puffing
  end type particle_puffing

  interface particle_puffing
    module procedure new_particle_puffing
  end interface particle_puffing
contains

function new_particle_puffing(n_puff, fueling_rate, valve_r, R, Z, phi, rng, seed, puff_t_dependent, &
  t_puff_start,t_puff_slope,fueling_rate_start,poly_R,poly_Z,boxpuff,rate_from_file,puff_rate_filename,sim_id) result(new)
  use mod_pcg32_rng,   only: pcg32_rng
  use mod_random_seed, only: random_seed
  use profiles,        only: readProf

  type(particle_puffing)    :: new

  integer, intent(in)           :: n_puff
  real*8, intent(in), optional  :: fueling_rate
  real*8, intent(in), optional  :: valve_r
  real*8, intent(in), optional  :: R, Z
  real*8, intent(in), optional  :: phi ! If no phi is given axisymmetric puffing will be excecuted.
  
  real*8, intent(in), optional  :: poly_R(4)
  real*8, intent(in), optional	:: poly_Z(4)
  logical, intent(in), optional :: boxpuff
  
  logical, intent(in), optional :: puff_t_dependent
  real*8, intent(in), optional  :: t_puff_start,t_puff_slope
  real*8, intent(in), optional  :: fueling_rate_start 
  
  logical, intent(in), optional :: rate_from_file
  character(len=64), intent(in), optional :: puff_rate_filename
  
  class(type_rng), intent(in), optional :: rng !< random number generator to use (deafult PCG32)
  integer, intent(in), optional         :: seed !< Seed for the RNG (default random_seed() on my_id + bcast)
  integer                               :: my_seed
  
  integer, intent(in), optional         :: sim_id

  !local variables
  integer :: len, i
  logical :: print_info
  real*8,dimension(:), allocatable :: t_pts
  real*8,dimension(:), allocatable :: puff_rate_pts

  new%n_puff = n_puff
  if (present(fueling_rate))       new%fueling_rate = fueling_rate
  if (present(R))                  new%R = R
  if (present(Z))                  new%Z = Z
  if (present(valve_r))            new%valve_r = valve_r
  if (present(phi))                new%phi = phi
  if (present(puff_t_dependent))   new%puff_t_dependent  = puff_t_dependent
  if (present(t_puff_start))       new%t_puff_start = t_puff_start
  if (present(t_puff_slope))       new%t_puff_slope = t_puff_slope
  if (present(fueling_rate_start)) new%fueling_rate_start = fueling_rate_start
  if (present(poly_R))             new%poly_R = poly_R
  if (present(poly_Z))             new%poly_Z = poly_Z
  if (present(boxpuff))            new%boxpuff = boxpuff
  if (present(rate_from_file))     new%rate_from_file = rate_from_file
  if (present(puff_rate_filename)) new%puff_rate_filename = puff_rate_filename
  
  if(present(sim_id)) then
    if(sim_id .eq. 0) then
      print_info = .true.
    else
      print_info = .false.
    end if
  else
    print_info = .true.
  end if
  
  if(print_info) write(*,*) "---------puff event info-----------"
  
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

  if(new%rate_from_file .and. present(puff_rate_filename)) then
    call readProf(t_pts, puff_rate_pts, len, puff_rate_filename)
    !some sanity checks, based on checks in read_num_profiles
    if ( len < 2 ) then 
      if (len == 1) then
        if(print_info) write(*,'(A80, ES14.4)') 'WARNING: puff rate input file only has 1 point. Assuming constant puff rate of ',puff_rate_pts
        t_pts = [0.d0,1.d0]
        puff_rate_pts = [puff_rate_pts, puff_rate_pts]
        !write(*,*) new%t_pts, new%puff_rate_pts
      else
        write(*,*) '  ERROR: Could not read the numerical profile for '//trim(puff_rate_filename)
        stop
      end if
    end if
    if ( t_pts(1)>1.d-6 ) then
      if(print_info) write(*,*) 'WARNING: Numerical '//trim(puff_rate_filename)//' input does not start at t=0'
    end if
    if ( len > 500 ) then
      if(print_info) write(*,*) 'WARNING: Numerical '//trim(puff_rate_filename)//' input has a very large number of points: ',len
    end if
    do i = 1, len-1
      if ( t_pts(i+1)<=t_pts(i) ) then
        write(*,*) 'ERROR: Numerical '//trim(puff_rate_filename)//' input not correct'
        write(*,*) '  t values do not increase in a strictly monotonic way'
        stop
      end if
    end do
    if ( (minval(t_pts) /= minval(t_pts)) .or. (minval(puff_rate_pts) /= minval(puff_rate_pts)) ) then
      write(*,*) 'ERROR: Numerical '//trim(puff_rate_filename)//' input contains NaNs'
      stop
    end if
    if ( minval(puff_rate_pts) < 0.d0 ) then
      write(*,*) 'ERROR: Numerical '//trim(puff_rate_filename)//' input has non-positive values'
      stop
    end if
    
    new%t_pts         = t_pts
    new%puff_rate_pts = puff_rate_pts
  end if

  !print where the puffing will be done from
  if (present(R) .and. present(Z) .and. present(valve_r) .and. (.not. new%boxpuff)) then 
    if(print_info) write(*,'(A30,3E12.2)') "puffing at valve (R/Z/radius)",R,Z,valve_r 
  else if (present(poly_R) .and. present(poly_Z) .and. new%boxpuff) then
    if(print_info) write(*,'(A23,4F7.2,A8, 4F7.2)') "puffing at box: poly_R",poly_R," poly_Z",poly_Z 
  else
    write(*,*) "Invalid puffing location set. Please fill (R & Z & valve_r) or (poly_R & poly_Z)"
    stop
  end if

  !print what the rate will be
  if(new%rate_from_file .and. present(puff_rate_filename)) then
    if(print_info) write(*,*) "puffing rate from file named ",trim(puff_rate_filename)
  else if((.not. new%rate_from_file) .and. new%puff_t_dependent .and. present(t_puff_start) .and. &
      present(t_puff_slope) .and. present(fueling_rate) ) then
    if(print_info) write(*,*) "puffing flat-ramp-flat with puff slope",t_puff_slope,"starting at", t_puff_start, " and ending at ",fueling_rate
  else if((.not. new%rate_from_file) .and. (.not. new%puff_t_dependent) .and. present(fueling_rate)) then
    if(print_info) write(*,*) "constant puffing rate ",fueling_rate
  else
    write(*,*) "Invalid combination of puffing rate variables set up"
    stop
  end if

end function new_particle_puffing

  ! Actually puff gass
subroutine do_particle_puffing(this,sim, ev)
	use mpi_mod
    ! !$ use omp_lib
	
	
  class(particle_puffing) , intent(inout) :: this
  type(particle_sim), intent(inout)       :: sim
  type(event), intent(inout), optional    :: ev !! <STIJN> is this nececary?

  integer :: ierr,i_scalar, n_free, j, k, n_group, i_elm, i_elm_new, ifail, i_p, to_puff, n_puff_local,i_rng
  logical, allocatable, dimension(:) :: is_free
  integer, allocatable, dimension(:) :: i_free
  real*8  :: delta_t, c, R, Z, phi, s, t
  real*8  :: R_new, Z_new, s_new, t_new, r_valve, theta
  real*8  :: vector_normal(3), u(5)
  real*8  :: fueling_rate_t !< possibly time dependent fueling rate

  integer ::	  puffed_this_step_local, all_puffed_this_step
  real*8  ::	  puff_weight_local, all_puff_weight

  integer :: it  !< index for loop through t_pts
  real*8  :: time    !< sim%time (s)

  if (sim%my_id .eq. 0) write(*,*) "Started puffing!"
  
  if (this%last_time .eq. 0.d0) then
    this%last_time = sim%time
    this%last_diag_time = sim%time
    if (abs(sim%time) .le. 1d-10) then
      this%last_time = 1d-10 ! so we actually do most of the sputtering should we start at 0
    end if
    return
  end if
  delta_t = sim%time - this%last_time
  this%last_time = sim%time
  
  if (this%n_puff .le. -1.d-6) then ! 0.d0
    if (sim%my_id .eq. 0) write(*,*) 'No puff quota set, exiting. --- n_puff == 0 this will now stop the program'
    stop
  end if
  
  n_puff_local = this%n_puff / sim%n_cpu !n_puff local is the amount of superparticles that will be puffed per MPI process.
 
  ! allocate(is_free(size(sim%groups(1)%particles,1))) 

  ! might be replaced with omp workshare, or just the array expression.
  ! there is an issue with derived type arrays in gfortran though, and this works
  ! !$omp parallel do default(none) shared(sim, this, n_free, i_free, is_free) &
  ! !$omp private(j)
  ! do j=1,size(sim%groups(1)%particles,1)
    ! is_free(j) = sim%groups(1)%particles(j)%i_elm .le. 0 
  ! end do
  ! !$omp end parallel do
  ! !$omp barrier
  ! n_free = count(is_free)
  ! allocate(i_free(n_free))
  ! k = 1
  ! do j=1,size(is_free,1)
    ! if (is_free(j)) then
      ! i_free(k) = j
      ! k = k+1
      ! if (sim%my_id .eq. 0) write(*,*) "Adding to the list number: ", j
    ! end if
  ! end do

  !============== Finding free particles !< make into a function?
  allocate(is_free(size(sim%groups(1)%particles,1))) 
  !$omp parallel do default(none) shared(sim, n_free, i_free, is_free) &
  !$omp private(j) schedule(dynamic, 100)
  do j=1,size(sim%groups(1)%particles,1) !sim%groups(1)%particles
    is_free(j) = sim%groups(1)%particles(j)%i_elm .le. 0  !< array T/F is particle is free
  end do
  !$omp end parallel do
  !$omp barrier
  n_free = count(is_free)
  allocate(i_free(n_free))
  k = 1
  do j=1,size(is_free,1)
    if (is_free(j)) then
      i_free(k) = j !< i_free(k) has index of free particle in  sim%groups(1)%particles(j)
      k = k+1
      !if (sim%my_id .eq. 0) write(*,*) "Adding to the list number: ", j
    end if
  end do
  ! ==================
  
  
  n_group = 1   ! Puffing Hydrogen (or actually the element at groups(1)) only
  ! Assuming the incoming gas at T=300K and a diatomic gas
  c = sqrt((7.d0/5.d0)*(300.d0+273.d0)*K_BOLTZ/(2.d0*sim%groups(n_group)%mass*ATOMIC_MASS_UNIT))
	if (.not. this%boxpuff) then
		call find_RZ(sim%fields%node_list, sim%fields%element_list, this%R, this%Z, R, Z, &
				   i_elm, s, t ,ifail)
				   
	else
		call find_RZ(sim%fields%node_list, sim%fields%element_list, sum(this%poly_R(1:2))/2.d0, sum(this%poly_Z(1:2))/2.d0, R, Z, &
				   i_elm, s, t ,ifail)
	endif
  if (ifail .ne. 0) then
    if (sim%my_id .eq. 0) write(*,*) "Warning: The valve location for puffing could not be found, maybe it was placed outside of the grid?"
    stop
  end if

  vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, &
          i_elm, s, t)

  !------------- Decide how many superparticles to initiate		 
  !Adjust amount of superparticles + fueling rate if we use time dependent puffing 
  ! same way as in model500
  ! write(*,*) "puff_t_dependent", this%puff_t_dependent
  if (this%rate_from_file) then
    to_puff = n_puff_local
    
    !find the linearly interpolated puffing rate at time sim%time from the input profile
    time = sim%time
    fueling_rate_t = -1.d0
    do it=this%it_last,size(this%t_pts) !< starting value of this%it_last = 2, thus this loops for all linear interpolations between t_pts(1) and t_pts(end)
      if(time .le. this%t_pts(it)) then
        if(it .eq. 2) then
          if(time .le. this%t_pts(1)) then !< if t<t_pts(1), take puff_rate = puff_rate_pts(1)
            fueling_rate_t = this%puff_rate_pts(1)
            exit
          end if
        end if
        !linear interpolation between two (t,puffrate) points:
        fueling_rate_t = this%puff_rate_pts(it-1) + (time - this%t_pts(it-1))*(this%puff_rate_pts(it)-this%puff_rate_pts(it-1))/(this%t_pts(it)-this%t_pts(it-1))
        this%it_last = it !< reduces computation time of this do loop, making rate from file limited by data, not computation, thus many points can be used
        exit
      end if
    end do
    if(fueling_rate_t < 0.d0) then !< if t not in between t_pts(1) and t_pts(end)
      if(time .gt. this%t_pts(size(this%t_pts))) then !< default to last point
        fueling_rate_t = this%puff_rate_pts(size(this%t_pts))
      else
        write(*,*) 'problem in loading rate from file',this%puff_rate_filename,time,this%it_last,size(this%t_pts),this%puff_rate_pts(size(this%t_pts))
        fueling_rate_t = 0.d0 
      end if
  end if
  else if (this%puff_t_dependent) then
    to_puff        = n_puff_local !int( maxval((/ time_dependent_puff(real(n_puff_local,8)       ,sim%time, this%t_puff_start,this%t_puff_slope) ,10.d0 /)))
	  fueling_rate_t = time_dependent_puff(this%fueling_rate ,sim%time, this%t_puff_start,this%t_puff_slope, this%fueling_rate_start)
    !write(*,*) "n_puff", this%n_puff, "to_puff", to_puff, "fueling_rate_t", fueling_rate_t
    !write(*,*) "to_puff_real" , maxval((/ time_dependent_puff(real(n_puff_local,8)       ,sim%time, this%t_puff_start,this%t_puff_slope) ,10.d0 /))
  else
    fueling_rate_t = this%fueling_rate
	  to_puff = n_puff_local
  end if !< time dependent puffing
  
  if (to_puff .gt. n_free) then
    write(*,*) "problem: not enough superparticles free for puffing (wanted/free)",to_puff,n_free
    to_puff = n_free
  end if

  if (sim%my_id .eq. 0) write(*,"(A,g12.4,A,g16.8, A)") "Planned puffing rate at time t:", sim%time, " is fueling_rate_t:",fueling_rate_t, "atoms/s"
  !-------------  
  
  
  puffed_this_step_local = 0
  puff_weight_local      = 0.d0
  select type (pa => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)
! #ifdef __GFORTRAN__
!  !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
! #else
 ! !$omp parallel do default(shared) &
 ! !$omp schedule(dynamic,10) &
 ! !$omp shared(sim, pa, this,i_free,c, vector_normal,                       &
 ! !$omp   to_puff,n_puff_local, delta_t,fueling_rate_t )                        &
 ! !$omp private(i_p, i_rng, j,k,u , R,Z,s,t,R_new,Z_new,s_new,t_new,     &
 ! !$omp  i_elm,i_elm_new,r_valve, theta,                                         &
 ! !$omp ifail)                                                                    &
 ! !$omp reduction(+:puffed_this_step_local,puff_weight_local)
    do j = 1, to_puff
      i_p = i_free(j)
      do !basically a do while .true. loop
	  
        !	    !$ i_rng = omp_get_thread_num()+1
        call this%rng(1)%next(u) !rng(1)
        if (.not. this%boxpuff) then
          r_valve = this%valve_r*sample_piecewise_linear(2, [0.d0, 1.d0], [1.d0, 0.d0], u(1))
          theta = TWOPI * u(2)
          R_new = this%R + r_valve * cos(theta) !< this puffs only on the edge of the circle? I thought a valve puff means in a circular domain (so also it's centre)
          Z_new = this%Z + r_valve * sin(theta)
        else
          s = u(1)
          t = u(2)
          R_new = this%poly_R(1)*(1.d0-s)*(1.d0-t) + this%poly_R(2)*s*(1.d0-t) + this%poly_R(3)*(1.d0-s)*t + this%poly_R(4)*s*t
          Z_new = this%poly_Z(1)*(1.d0-s)*(1.d0-t) + this%poly_Z(2)*s*(1.d0-t) + this%poly_Z(3)*(1.d0-s)*t + this%poly_Z(4)*s*t
        endif
        
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, R, Z, s, t, i_elm, &
        R_new, Z_new, s_new, t_new, i_elm_new, ifail)
        !call find_RZ(sim%fields%node_list, sim%fields%element_list, R_new, Z_new, R, Z, &
        !       i_elm_new, s_new, t_new, ifail)
        if (ifail .ge. 0) exit
      end do
      R     = R_new
      Z     = Z_new
      s     = s_new
      t     = t_new
      i_elm = i_elm_new
      if (this%phi .lt. 0.d0) then
        pa(i_p)%x(3) = TWOPI*u(3)
      else 
        pa(i_p)%x(3) = this%phi 
      end if
      pa(i_p)%x(1:2)  = [R, Z]
      pa(i_p)%st(1:2) = [s, t]
      pa(i_p)%i_elm   = i_elm
      pa(i_p)%weight  = 1.d0/real(n_puff_local,8) * delta_t * fueling_rate_t / real(sim%n_cpu,8)
      pa(i_p)%v       = c * sample_cosine(u(4:5), vector_normal)   ! <STIJN> Maybe this needs to be isotropic, or 1+cos like?
      pa(i_p)%q       = 0_1
      if (sim%groups(1)%particles(i_p)%weight  .le. 1.d-2) then ! if the weight is too low. 
        if (sim%groups(1)%particles(i_p)%weight .ne. 0.d0) write(*,*) "problem: particle dumped in mod particle puffing due to extremely low weight",sim%groups(1)%particles(i_p)%weight
        sim%groups(1)%particles(i_p)%i_elm = 0
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
write(*,'(A60,I7,E14.6)') "Superparticles, weight puffed this puffing action action = ", all_puffed_this_step, all_puff_weight
endif
  
  
end subroutine do_particle_puffing

pure function time_dependent_puff(max_puff,time, t_puff_start,t_puff_slope, min_puff) result(to_puff)
real*8,intent(in)   :: max_puff, min_puff
real*8              :: to_puff
real*8,intent(in)    :: t_puff_start,t_puff_slope
real*8,intent(in)    :: time

! real*8               :: t_norm

! n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
! rho_norm  = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
! t_norm    = sqrt((MU_ZERO * CENTRAL_MASS * ATOMIC_MASS_UNIT * CENTRAL_DENSITY * 1.d20    ))                           ! t_SI   = t_norm * t_jorek


if (time-(t_puff_start+t_puff_slope) .ge. 0.d0) then
	to_puff = max_puff
elseif (time-t_puff_start .ge. 0.d0) then
	to_puff = min_puff+ (max_puff -min_puff) * (time-t_puff_start)/(t_puff_slope)  
else
    to_puff = min_puff !default = 0.d0
endif
end function time_dependent_puff

end module