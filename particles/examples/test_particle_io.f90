! Program use for testing the particle IO
program test_particle_io
  use mpi
  use fruit
  use fruit_mpi
  use hdf5,                           only: H5F_ACC_TRUNC_F
  use mod_gnu_rng,                    only: set_seed_sys_time
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_particle_types,             only: particle_fieldline_id,particle_fieldline
  use mod_particle_types,             only: particle_gc_id,particle_gc_vpar_id
  use mod_particle_types,             only: particle_kinetic_id,particle_kinetic_leapfrog_id
  use mod_particle_types,             only: particle_kinetic_relativistic_id
  use mod_particle_types,             only: particle_gc_relativistic_id
  use mod_particle_sim,               only: particle_sim
  use mod_event,                      only: event,with
  use mod_io_actions,                 only: read_action,write_action
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: fill_particles
  use mod_particle_assert_equal,      only: assert_equal_particle
  implicit none
  !> Variable declaration ----------------------------------------
  integer,parameter                :: master_task=0
  type(particle_sim)               :: sim_write,sim_read
  type(event),dimension(1),target  :: event_write,event_read
  type(read_action)                :: read_particle
  type(write_action)               :: write_particle
  integer                          :: ii,jj,ifail
  integer                          :: n_groups,fill_particle_type
  integer                          :: mpi_comm_io,mpi_info_io
  integer                          :: file_access_write
  integer,dimension(2)             :: Z_interval,seed_interval
  integer,dimension(:),allocatable :: p_types,n_particles
  real*8,dimension(2)              :: mass_interval
  logical                          :: do_write,do_read
  logical                          :: use_hdf5_access_properties
  logical                          :: use_native_hdf5_mpio
  logical                          :: mpio_collective
  logical                          :: remove_file,file_exist
  character(len=24)                :: filename
  !> Define test parameters --------------------------------------
  filename           = "particle_restart_test.h5"
  do_write           = .true.
  do_read            = .true.
  ifail              = 0 
  n_groups           = 7  !< number of groups = number of particle types
  allocate(n_particles(n_groups)); allocate(p_types(n_groups))
  n_particles        = [1000,2000,1000,3000,4000,7000,1000] !< number of particles
  seed_interval      = [12,31244553]
  mass_interval      = [1d-5,1d0] 
  Z_interval         = [-1,100] 
  p_types = [particle_fieldline_id,particle_gc_id,particle_gc_vpar_id,&
            particle_kinetic_id,particle_kinetic_leapfrog_id,&
            particle_kinetic_relativistic_id,particle_gc_relativistic_id]
  mpi_comm_io                = MPI_COMM_WORLD
  mpi_info_io                = MPI_INFO_NULL
  use_hdf5_access_properties = .true.
  use_native_hdf5_mpio       = .true.
  file_access_write          = H5F_ACC_TRUNC_F
  mpio_collective            = .true.
  remove_file                = .false.
  !> Initialise --------------------------------------------------
  call sim_write%initialize(n_groups,.true.,do_jorek_init_in=.false.)
  write_particle = write_action(filename,file_access_in=file_access_write,&
  use_native_hdf5_mpio_in=use_native_hdf5_mpio,&
  use_hdf5_access_properties_in=use_hdf5_access_properties,&
  mpi_comm_in=mpi_comm_io,mpi_info_in=mpi_info_io,mpio_collective_in=mpio_collective)
  if(do_write) then
    call allocate_one_particle_list_type(n_groups,n_particles,p_types,sim_write%groups,ifail)
    call fill_particles(n_groups,sim_write%groups,fill_particle_in=fill_particle_type,&
    rank_in=sim_write%my_id)
    call set_seed_sys_time(seed_interval)
    do ii=1,n_groups
      call gnu_rng_interval(mass_interval,sim_write%groups(ii)%mass)
      call gnu_rng_interval(Z_interval,sim_write%groups(ii)%Z)
    enddo
  endif
  !> Read
  read_particle  = read_action(filename=filename,mpi_comm_in=mpi_comm_io,&
  mpi_info_in=mpi_info_io,use_hdf5_access_properties_in=use_hdf5_access_properties)
  if(do_read) then
    event_write = [event(write_particle)]; event_read = [event(read_particle)];
    call sim_read%initialize(n_groups,.true.,my_id=sim_write%my_id,&
    n_cpu=sim_write%n_cpu,do_jorek_init_in=.false.)
  endif
  !> Read - write operations -------------------------------------
  if(do_write) call with(sim_write,event_write)
  if(do_read)  call with(sim_read,event_read)
  !> Check that files are equals ---------------------------------
  if(do_write.and.do_read) then
    do ii=2,n_groups
      select type (p=>sim_write%groups(ii)%particles)
      type is (particle_fieldline)
        select type(p2=>sim_read%groups(ii)%particles)
          type is (particle_fieldline)
          do jj=1,size(p)
            p(jj)%B_hat_prev = 0d0;
          enddo
        endselect
      end select
      write(*,*) "Size particle old: ",size(sim_write%groups(ii)%particles),&
      " size new: ",size(sim_read%groups(ii)%particles)
      call init_fruit
      call fruit_init_mpi_xml(sim_write%my_id)
      call fruit_hide_dots
      call assert_equals(sim_write%groups(ii)%Z,sim_read%groups(ii)%Z,"List Z mismatch!")
      call assert_equals(sim_write%groups(ii)%mass,&
      sim_read%groups(ii)%mass,"List mass mismatch!")
      call assert_equal_particle(size(sim_write%groups(ii)%particles),&
      sim_write%groups(ii)%particles,sim_read%groups(ii)%particles)
      call fruit_summary_mpi(sim_write%n_cpu,sim_write%my_id)
      call fruit_summary_mpi_xml(sim_write%n_cpu,sim_write%my_id)
      call fruit_finalize_mpi(sim_write%n_cpu,sim_write%my_id)
    enddo
  endif
  !> Clean-up ---------------------------------------------------- 
  if(allocated(p_types))     deallocate(p_types); 
  if(allocated(n_particles)) deallocate(n_particles);
  if((sim_write%my_id.eq.master_task).and.remove_file) then
    inquire(file=filename,exist=file_exist)
    if(file_exist) call system("rm -rf "//trim(filename))
  endif
  call sim_write%finalize
end program test_particle_io

