!> functionality to do with the allocation of the particle arrays for each particle group
!> including load balancing. 
!> In the future this module will be a user of the different initialization schemes defined
!> in mod_initialise_particles.f90 (or similar) once that system has been reworked.

module mod_particle_allocation
use mod_particle_sim, only: particle_sim
use mpi

implicit none
contains

  !> allocate an empty array of a specified particle type for a specific particle group for each mpi process
  subroutine allocate_particle_array(sim, group_num, particle_type_str, n_particles_per_proc, mpi_comm_loc)
    use mod_particle_types
  
    implicit none
  
    type(particle_sim),            intent(inout) :: sim
    integer,                       intent(in)    :: group_num
    character(len=:), allocatable, intent(in)    :: particle_type_str
    integer,                       intent(in)    :: n_particles_per_proc
    integer, optional,             intent(in)    :: mpi_comm_loc
    integer                                      :: errorcode, ierr
  
    select case (trim(particle_type_str))
    case ("particle_kinetic")
      allocate(particle_kinetic::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_kinetic_leapfrog")
      allocate(particle_kinetic_leapfrog::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_gc")
      allocate(particle_gc::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_gc_vpar")
      allocate(particle_gc_vpar::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_gc_Qin")
      allocate(particle_gc_Qin::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_fieldline")
      allocate(particle_fieldline::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_kinetic_relativistic")
      allocate(particle_kinetic_relativistic::sim%groups(group_num)%particles(n_particles_per_proc))
    case ("particle_gc_relativistic")
      allocate(particle_gc_relativistic::sim%groups(group_num)%particles(n_particles_per_proc))
    case default
      write(*,*) "Error: missing type name declaration ",trim(particle_type_str)," for reading: ABORT!"
      if (present(mpi_comm_loc)) then
        call MPI_Abort(mpi_comm_loc,errorcode,ierr)
      else
        stop
      endif
    end select
  
  end subroutine allocate_particle_array

  !> calculate the distribution of particles for a group across mpi processers (load balancing)
  subroutine calc_n_particles_per_mpi(n_particles_tot, n_mpi, n_particles_per_mpi, master_task)
    implicit none

    integer,  intent(in)                                  :: n_particles_tot       !< total number of particles for a particle group
    integer,  intent(in)                                  :: n_mpi                 !< number of mpi processes being used
    integer,  dimension(:),  allocatable,  intent(inout)  :: n_particles_per_mpi   !< array with the number of particles per mpi process
    integer,  optional                                    :: master_task           !< which mpi process is the "master" and will hence hold
                                                                                   !< the remainder of the particles after division (default 0)

    if (.not. present(master_task)) master_task = 0

    n_particles_per_mpi = int(n_particles_tot)/n_mpi
    n_particles_per_mpi(master_task+1) = int(n_particles_tot) - (n_mpi-1)*n_particles_per_mpi(master_task+1)

  end subroutine calc_n_particles_per_mpi




end module mod_particle_allocation