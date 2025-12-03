!> Module for routines doing general tasks in the output file such as general diagnostics (non-interaction specific) and writing new header lines
module mod_output_file_routines
  use mod_particle_sim
  use mod_particle_types
  use constants, only: ATOMIC_MASS_UNIT
  use mod_event, only: mpi_minmeanmax
  use mpi
  use phys_module, only: tracked_group_id
  !$ use omp_lib

  implicit none
  private
  public write_to_outputfile

  integer,parameter :: num_tracked=7      !< number of tracked variables (see conserv_obj)

  real*8  :: conserv_obj(num_tracked)=0.d0 !< superparticles, real particles, momentum (3 directions), energy
  real*8  :: last_time=0                   !< omp_get_wtime
  logical :: next_call_only_header=.false. !< whether the coming block should only have a header (.true.) or also have conservation and/or timing information (.false.)

contains

!> Write the new header line for the coming block, with possibly extra general output of the previous block such as
!> how much computational time the block previous block took, and the changes in particles/momentum/energy of a tracked species
subroutine write_to_outputfile(sim,what,next_block_only_header)
  implicit none

  type(particle_sim), intent(in) :: sim
  character(len=*),   intent(in) :: what
  logical, optional,  intent(in) :: next_block_only_header !< whether the coming block should only have a header (.true.) or also have conservation and/or timing information (.false.). Default .false.

  if(.not. next_call_only_header) then
    if (tracked_group_id /= "non") call conservation_block(sim)

    call cpu_time_block(sim)
  endif

  if(present(next_block_only_header)) then
    next_call_only_header=next_block_only_header
  else
    next_call_only_header=.false.
  endif
  
  if(sim%my_id .ne. 0) return

  write(*,'(A100)') "===================================================================================================="
  write(*,"(2X,A)") what
  write(*,'(A100)') "===================================================================================================="

end subroutine


!> keeps track of global superparticle, particle, momentum and energy balances for 1 group,
!> and writes the change and new values to the logfile, and times the previous block
subroutine conservation_block(sim)
  implicit none
  
  class(particle_sim), target, intent(in) :: sim

  real*8  :: old(num_tracked), diff(num_tracked)
  integer :: j, ierr, group_num_tracked
  real*8    :: particles_remaining, particles_elm_lt0, momentum_remaining(3), energy_remaining, all_particles, all_momentum(3), all_energy, all_elm_lt0
  integer   :: superparticles_remaining,all_superparticles,closest_iteration

  group_num_tracked = group_num_from_id(sim, tracked_group_id)

  particles_remaining = 0.d0
  particles_elm_lt0   = 0.d0
  momentum_remaining  = 0.d0
  energy_remaining    = 0.d0
  superparticles_remaining = 0

  select type (particles => sim%groups(group_num_tracked)%particles)
  type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
    !$omp parallel do default(none) &
    !$omp shared(sim,group_num_tracked) &
#endif
    !$omp reduction(+:particles_remaining, particles_elm_lt0, momentum_remaining, energy_remaining,superparticles_remaining)
      do j=1,size(particles,1)

        if (particles(j)%i_elm .lt. 0) particles_elm_lt0 = particles_elm_lt0 + particles(j)%weight


        if (particles(j)%i_elm .le. 0) cycle

        particles_remaining = particles_remaining + particles(j)%weight
        momentum_remaining  = momentum_remaining  + particles(j)%weight * particles(j)%v *sim%groups(group_num_tracked)%mass * ATOMIC_MASS_UNIT
        energy_remaining    = energy_remaining    + particles(j)%weight * 0.5d0 * sim%groups(group_num_tracked)%mass * ATOMIC_MASS_UNIT * dot_product(particles(j)%v,particles(j)%v)
        superparticles_remaining = superparticles_remaining + 1

      enddo !j
    !omp end parallel do
  class default
      if(sim%my_id == 0) write(*,*) "conservation_block() only implemented for particle_kinetic_leapfrog, not for particle type of group ",sim%groups(group_num_tracked)%id
      return
  end select

  call MPI_REDUCE(particles_remaining, all_particles,         1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(particles_elm_lt0,   all_elm_lt0,           1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_remaining,  all_momentum,          3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_remaining,    all_energy,            1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(superparticles_remaining,all_superparticles,1, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  
  if (sim%my_id .eq. 0) then
    old = conserv_obj

    conserv_obj(1:3) = [real(all_superparticles),all_particles, all_elm_lt0]
    conserv_obj(4:6) = all_momentum
    conserv_obj(7)   = all_energy

    diff = 0.d0
    diff = conserv_obj - old
    
    write(*,"(A)") "conservation checks --------------------------------------------------------------"
    write(*,"(A,7A15)") "qty: ", "superparticles", "particles", "w ielm < 0", "momentum R", "momentum Z", "momentum phi", "energy"
    write(*,"(A,7es15.5)") "diff ",diff
    write(*,"(A,7es15.5)") "new  ",conserv_obj
  endif !(sim%my_id .eq. 0)

end subroutine conservation_block

!> determines the cpu time spent since the last time this function was called (so it times the previous block)
!> and writes this to the output file
subroutine cpu_time_block(sim)
  implicit none
  
  class(particle_sim), target, intent(in) :: sim

  real*8 :: now, mmm(3)

  ! the following are executed only if omp is used, they are not actual comments

  !$ if(last_time < 0.d0) last_time = 0
  !$ now = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(now - last_time)
  !$ last_time = now

  !$ if(sim%my_id == 0) write(*,"(A,3f17.5,A)") "block done in (min/mean/max) ", mmm, " s"
end subroutine cpu_time_block

end module mod_output_file_routines
