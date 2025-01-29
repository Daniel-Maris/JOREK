!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
!> The documentation of the module is at link: https://www.jorek.eu/wiki/doku.php?id=jorek-particles_i_o
module mod_particle_io
implicit none
private
public write_simulation_hdf5,read_simulation_hdf5,get_simulation_hdf5_time

!> Module wide variables
integer,parameter :: master_task=0
integer,parameter :: n_cpu_1=1
integer,parameter :: group_name_len=12
contains

!> Export all particles using HDF5 IO
!> Parallel support to writing operations is provided
!> by MPI enabled HDF5 procedures.
!> inputs:
!>   sim:                        (particle_sim) particle simulation object
!>   filename:                   (character)(N) name of the output file
!>   file_access_in:             (integer)(optional) define how the HDF5 should behave:
!>                               when opening or creating a new file.
!>                               default: H5F_ACC_TRUNC_F (truncate fille if already exists
!>                               create otherwise)
!>   use_native_hdf5_mpio:       (logical)(optional) if true, the native hdf5-mpio is used
!>                               for parallel writing, otherwise data are first gathered 
!>                               in the master task node and then written in serial,
!>                               default: .false.
!>   use_hdf5_access_properties: (logical)(optional) HDF5 file access property
!>                               if must be set to .false. for parallel I/O
!>                               default: .true. 
!>   collective_mpio_in:         (logical)(optional) define whether MPIO collective calls are
!>                               performed default: true
!>   mpi_comm_in:                (integer)(optional) MPI communicator identifier
!>   mpi_info_in:                (integer)(optional) MPI info structre for parallel IO
subroutine write_simulation_hdf5(sim,filename,file_access_in,use_native_hdf5_mpio_in,&
use_hdf5_access_properties,collective_mpio_in,mpi_comm_in,mpi_info_in)
  use phys_module
  use mpi
  use hdf5,               only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F
  use hdf5,               only: H5Gcreate_f,H5Gclose_f
  use hdf5_io_module,     only: HDF5_open_or_create,HDF5_close
  use hdf5_io_module,     only: HDF5_integer_saving 
  use hdf5_io_module,     only: HDF5_real_saving,HDF5_char_saving
  use hdf5_io_module,     only: HDF5_array1D_saving_int_native_or_gatherv
  use hdf5_io_module,     only: HDF5_array1D_saving_r4_native_or_gatherv
  use hdf5_io_module,     only: HDF5_array1D_saving_native_or_gatherv
  use hdf5_io_module,     only: HDF5_array2D_saving_native_or_gatherv
  use hdf5_io_module,     only: HDF5_array3D_saving_native_or_gatherv
  use mod_particle_types, only: particle_arrays_from_list
  use mod_particle_types, only: deallocate_particle_arrays
  use mod_particle_sim,   only: particle_sim
  implicit none
  !> parameters
  integer,parameter             :: master_rank=0
  integer(HSIZE_T),parameter    :: i0_HSIZE_T=int(0,kind=HSIZE_T)
  !> input variables
  type(particle_sim),intent(in) :: sim
  character(len=*),  intent(in) :: filename 
  integer, intent(in), optional :: file_access_in,mpi_comm_in,mpi_info_in
  logical, intent(in), optional :: use_native_hdf5_mpio_in
  logical, intent(in), optional :: use_hdf5_access_properties,collective_mpio_in
  !> variables
  integer                       :: file_access_loc
  integer                       :: mpi_comm_loc,mpi_info_loc 
  integer                       :: ii,jj,ierr,h5err,n_groups,n_particles
  integer                       :: n_particles_per_group
  integer(HID_T)                :: file_id,group_id
  integer(HSIZE_T)              :: n_particles_offset
  integer,  dimension(:),    allocatable :: n_particles_loc,particle_displacement
  integer,  dimension(:),    allocatable :: i_elm_arr,i_life_arr
  integer,  dimension(:),    allocatable :: q_arr
  integer,  dimension(:,:),  allocatable :: n_particles_glob
  real*4,   dimension(:),    allocatable :: t_birth_arr
  real*8,   dimension(:),    allocatable :: weight_arr,v_1d_arr
  real*8,   dimension(:),    allocatable :: E_arr,mu_arr,vpar_arr
  real*8,   dimension(:),    allocatable :: B_norm_arr,vpar_m_arr,Bn_k_arr
  real*8,   dimension(:,:),  allocatable :: st_arr,x_arr,B_hat_prev_arr,v_2d_arr
  real*8,   dimension(:,:),  allocatable :: x_m_arr,Astar_m_arr,Astar_k_arr
  real*8,   dimension(:,:),  allocatable :: dBn_k_arr,Bnorm_k_arr,E_k_arr
  real*8,   dimension(:,:,:),allocatable :: dAstar_k_arr
  logical                                :: use_hdf5_parallel,use_gatherv_mpio
  logical                                :: create_access_plist,collective_mpio_loc
  character(len=group_name_len)          :: group_name
  character(len=:),          allocatable :: particle_type_str

  !> preparation
  h5err = 0; use_hdf5_parallel = .false.     !< do not use parallel HDF5 by default
  use_gatherv_mpio  = .not.use_hdf5_parallel !< use MPI gatherv for collecting all data for writing
  if(present(use_native_hdf5_mpio_in)) then 
    use_hdf5_parallel = use_native_hdf5_mpio_in
    use_gatherv_mpio  = .not.use_native_hdf5_mpio_in
  endif
<<<<<<< HEAD

  ! Create group to write particle groups in
  call h5gcreate_f(file, "/groups", group_id, hdferr)
  call h5gclose_f(group_id, hdferr)

  ! Write the time
  call HDF5_real_saving(file,sim%time,'/time')

  ! Write n_part_groups and part_groups_in_use
  call HDF5_integer_saving(file,n_part_groups,'/n_part_groups')
end if



! saving for each group
if (allocated(sim%groups)) then 
  do i=1,size(sim%groups,1)
    if (.not. allocated(sim%groups(i)%particles)) then
      write(*,*) "WARNING: group ", i, " not allocated, exiting"
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    ! Find the number of particles on each node
    n_here = size(sim%groups(i)%particles,1)
    call MPI_Gather(n_here,1,MPI_INTEGER,particles_per_proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_total = sum(particles_per_proc,1) ! set it on other processes as well (to
    ! 0) so they can allocate the useless arrays

    if (my_id .eq. 0) then
      ! Create group to write in
      write(group_name,"(A,i0.3,A)") "/groups/", i, "/"
      call h5gcreate_f(file, group_name, group_id, hdferr)
      call h5gclose_f(group_id, hdferr)
    end if

    ! particle_base properties
    ! x
    allocate(x(3,n_here), x_all(3,n_total))
    do j=1,n_here
      x(:,j) = sim%groups(i)%particles(j)%x
    end do
!    call MPI_Gatherv(x(:,:), 3*n_here, MPI_REAL8, &
!      x_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    call MPI_Gatherv(x(1,:), n_here, MPI_REAL8, &
                     x_all(1,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(x(2,:), n_here, MPI_REAL8, &
                     x_all(2,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(x(3,:), n_here, MPI_REAL8, &
                     x_all(3,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! st
    allocate(st(2,n_here), st_all(2,n_total))
    do j=1,n_here
      st(:,j) = sim%groups(i)%particles(j)%st
    end do
!    call MPI_Gatherv(st(:,:), 2*n_here, MPI_REAL8, &
!      st_all(:,:), particles_per_proc*2, [(sum(particles_per_proc(0:i-1),1)*2, i=0,n_cpu-1)], &
!      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(st(1,:), n_here, MPI_REAL8, &
                     st_all(1,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(st(2,:), n_here, MPI_REAL8, &
                     st_all(2,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! weight
    allocate(weight(n_here), weight_all(n_total))
    do j=1,n_here
      weight(j) = sim%groups(i)%particles(j)%weight
    end do
    call MPI_Gatherv(weight(:), n_here, MPI_REAL8, &
      weight_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! i_elm
    allocate(i_elm(n_here), i_elm_all(n_total))
    do j=1,n_here
      i_elm(j) = sim%groups(i)%particles(j)%i_elm
    end do
    call MPI_Gatherv(i_elm(:), n_here, MPI_INTEGER, &
      i_elm_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! i_life
    allocate(i_life(n_here), i_life_all(n_total))
    do j=1,n_here
      i_life(j) = sim%groups(i)%particles(j)%i_life
    end do
    call MPI_Gatherv(i_life(:), n_here, MPI_INTEGER, &
      i_life_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! t_birth
    allocate(t_birth(n_here), t_birth_all(n_total))
    do j=1,n_here
      t_birth(j) = sim%groups(i)%particles(j)%t_birth
    end do
    call MPI_Gatherv(t_birth(:), n_here, MPI_REAL4, &
      t_birth_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

    ! Write out stuff depending on particle type
    select type (p => sim%groups(i)%particles)

    type is (particle_kinetic)
      particle_type_name = 'particle_kinetic'

      ! v
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%v
      end do
!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      call MPI_Gatherv(v(1,:), n_here, MPI_REAL8, &
                       v_all(1,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(2,:), n_here, MPI_REAL8, &
                       v_all(2,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(3,:), n_here, MPI_REAL8, &
                       v_all(3,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_kinetic_leapfrog)
      particle_type_name = 'particle_kinetic_leapfrog'

      ! v
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%v
      end do
!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      call MPI_Gatherv(v(1,:), n_here, MPI_REAL8, &
                       v_all(1,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(2,:), n_here, MPI_REAL8, &
                       v_all(2,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(3,:), n_here, MPI_REAL8, &
                       v_all(3,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_gc)
      particle_type_name = 'particle_gc'

      ! E
      allocate(E(n_here), E_all(n_total))
      do j=1,n_here
        E(j) = p(j)%E
      end do
      call MPI_Gatherv(E(:), n_here, MPI_REAL8, &
        E_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,E_all,n_total,group_name//"E")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(E,mu,q,E_all,mu_all,q_all)

    type is (particle_gc_vpar)
      particle_type_name = 'particle_gc_vpar'

      ! Vpar
      allocate(Vpar(n_here), Vpar_all(n_total))
      do j=1,n_here
        Vpar(j) = p(j)%Vpar
      end do
      call MPI_Gatherv(Vpar(:), n_here, MPI_REAL8, &
        Vpar_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! B_norm
      allocate(B_norm(n_here), B_norm_all(n_total))
      do j=1,n_here
        B_norm(j) = p(j)%B_norm
      end do
      call MPI_Gatherv(B_norm(:), n_here, MPI_REAL8, &
        B_norm_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,Vpar_all,n_total,group_name//"Vpar")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving(file,B_norm_all,n_total,group_name//"B_norm")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(Vpar,mu,B_norm,q,Vpar_all,mu_all,B_norm_all,q_all)

    type is (particle_fieldline)
      particle_type_name = 'particle_fieldline'
      ! v
      allocate(v1(n_here), v1_all(n_total))
      do j=1,n_here
        v1(j) = p(j)%v
      end do
      call MPI_Gatherv(v1(:), n_here, MPI_REAL8, &
        v1_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,v1_all,n_total,group_name//"v")
      end if
      deallocate(v1, v1_all)

    type is (particle_kinetic_relativistic)
      particle_type_name = 'particle_kinetic_relativistic'
 
      ! momenta
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%p
      end do

!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(1,:), n_here, MPI_REAL8, &
                       v_all(1,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(2,:), n_here, MPI_REAL8, &
                       v_all(2,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(3,:), n_here, MPI_REAL8, &
                       v_all(3,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
 
      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v") 
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q") 
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_gc_relativistic)
      particle_type_name = 'particle_gc_relativistic'

      ! momenta (parallel momentum and magnetic moment)
      allocate(v(2,n_here), v_all(2,n_total))
      do j=1,n_here
        v(:,j) = p(j)%p
      end do

!      call MPI_Gatherv(v(:,:), 2*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*2, [(sum(particles_per_proc(0:i-1),1)*2, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(1,:), n_here, MPI_REAL8, &
                       v_all(1,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(2,:), n_here, MPI_REAL8, &
                       v_all(2,:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,2,n_total,group_name//"v") 
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q") 
      end if
      deallocate(v,q,v_all,q_all)

    class default
      write(*,*) "error: missing type name declaration for write"
      call exit(1)
    end select


    ! It is important to do the gathering first, because that is the collective part
    if (my_id .eq. 0) then
      call HDF5_array2D_saving(file,x_all,3,n_total,group_name//"x")
      call HDF5_array2D_saving(file,st_all,2,n_total,group_name//"st")
      call HDF5_array1D_saving(file,weight_all,n_total,group_name//"weight")
      call HDF5_array1D_saving_int(file,i_elm_all,n_total,group_name//"i_elm")
      call HDF5_array1D_saving_int(file,i_life_all,n_total,group_name//"i_life")
      call HDF5_array1D_saving_r4(file,t_birth_all,n_total,group_name//"t_birth")
=======
  file_access_loc = H5F_ACC_TRUNC_F !< truncate the file by default
  if(present(file_access_in)) file_access_loc = file_access_in;
  create_access_plist = .false. !< serial access by default
  if(present(use_hdf5_access_properties)) create_access_plist = .not.use_hdf5_access_properties
  if(use_hdf5_parallel) create_access_plist = .true.
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_info_loc = mpi_info_in
  collective_mpio_loc = .true. !< enable collective MPIO applications by default
  if(present(collective_mpio_in)) collective_mpio_loc = collective_mpio_in
  !> allocate the gatherv displacement array if required
  if(use_gatherv_mpio) allocate(particle_displacement(sim%n_cpu),source=0) 
  !> create the hdf5 file and the groups fields
  if(use_gatherv_mpio.and.(sim%my_id.eq.master_rank)) call HDF5_open_or_create(&
  trim(filename),file_id,ierr=h5err,file_access=file_access_loc)
  if(use_hdf5_parallel) call HDF5_open_or_create(trim(filename),file_id,&
  ierr=h5err,file_access=file_access_loc,create_access_plist_in=create_access_plist,& 
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  if(h5err.gt.0) then
    write(*,*) "Failed to create or open the ",filename," file: ",h5err,&
    " MPI task: ",sim%my_id,", ABORT!"
    call MPI_Abort(mpi_comm_loc,-1,ierr)
  endif
  if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.use_hdf5_parallel) then
    call H5Gcreate_f(file_id,"/groups",group_id,h5err) !< create particle groups
    call H5Gclose_f(group_id,h5err)
    !> write the time in HDF5 file, we assume that each MPI task reached the same physical time
    !> if HDF5-MPIO is used, the routine must be executed by all tasks for avoiding deadlocks
    call HDF5_real_saving(file_id,sim%time,"/time") 
  endif
  !> check if loops are allocated and loop on them
  if(allocated(sim%groups)) then
    !> it is assumed that all processors has the same number of groups but
    !> different size of particle lists per each group so we gather the
    !> particle list size for all groups firstly
    n_groups = size(sim%groups,1); allocate(n_particles_loc(n_groups)); n_particles_loc=0;
    do ii=1,n_groups
      n_particles_loc(ii) = size(sim%groups(ii)%particles,1)
    enddo
    allocate(n_particles_glob(sim%n_cpu*n_groups,1)); n_particles_glob=0;
    call MPI_Allgather(n_particles_loc,n_groups,MPI_INTEGER,n_particles_glob(:,1),n_groups,&
    MPI_INTEGER,mpi_comm_loc,ierr)
    n_particles_glob = transpose(reshape(n_particles_glob,(/n_groups,sim%n_cpu/)))
    !> write particle list data into HDF5 file
    do ii=1,n_groups
      if(.not.allocated(sim%groups(ii)%particles)) then
        if(sim%my_id.eq.master_task) write(*,*) "WARNING: particle list N# ",&
        ii,"  not allocated: skip!"; cycle;
      endif
      !> number of total particles for the group
      n_particles_per_group = sum(n_particles_glob(:,ii))
      if(use_gatherv_mpio) then
        n_particles_offset    = i0_HSIZE_T
        particle_displacement(2:sim%n_cpu) = [(sum(n_particles_glob(1:jj,ii)),jj=1,sim%n_cpu-1)]
      else
        n_particles_offset  = int(sum(n_particles_glob(1:sim%my_id,ii)),kind=HSIZE_T)
      endif
      !> create the HDF5 group for the particle list
      write(group_name,"(A,i0.3,A)") "/groups/",ii,"/"
      if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.use_hdf5_parallel) then
        call H5Gcreate_f(file_id,trim(group_name),group_id,h5err)
        call H5Gclose_f(group_id,h5err)
      endif
      !> reorganize and store the particle data in congruent arrays
      call particle_arrays_from_list(sim%groups(ii)%particles,n_particles,&
      i_elm_arr,i_life_arr,q_arr,t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,&
      vpar_arr,B_norm_arr,vpar_m_arr,st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,&
      Astar_m_arr,Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr,&
      particle_type_str)
      !> write data in HDF5 file
      if(allocated(i_elm_arr)) call HDF5_array1D_saving_int_native_or_gatherv(&
      file_id,i_elm_arr,n_particles_per_group,trim(trim(group_name)//"i_elm"),&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)
>>>>>>> 46d8f4abddc8a950131af8d646356ad3e65b16d2
      
      if(allocated(i_life_arr)) call HDF5_array1D_saving_int_native_or_gatherv(&
      file_id,i_life_arr,n_particles_per_group,trim(group_name)//"i_life",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)
      
      if(allocated(q_arr)) call HDF5_array1D_saving_int_native_or_gatherv(&
      file_id,q_arr,n_particles_per_group,trim(group_name)//"q",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)
      
      if(allocated(t_birth_arr)) call HDF5_array1D_saving_r4_native_or_gatherv(&
      file_id,t_birth_arr,n_particles_per_group,trim(group_name)//"t_birth",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(weight_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,weight_arr,n_particles_per_group,trim(group_name)//"weight",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(v_1d_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,v_1d_arr,n_particles_per_group,trim(group_name)//"v",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(E_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,E_arr,n_particles_per_group,trim(group_name)//"E",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(mu_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,mu_arr,n_particles_per_group,trim(group_name)//"mu",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(vpar_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,vpar_arr,n_particles_per_group,trim(group_name)//"Vpar",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(B_norm_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,B_norm_arr,n_particles_per_group,trim(group_name)//"B_norm",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc) 

      if(allocated(vpar_m_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,vpar_m_arr,n_particles_per_group,trim(group_name)//"Vpar_m",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc) 

      if(allocated(st_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,st_arr,size(st_arr,1),n_particles_per_group,trim(group_name)//"st",&
      use_gatherv_mpio,dim2_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[i0_HSIZE_T,n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)  

      if(allocated(x_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,x_arr,size(x_arr,1),n_particles_per_group,trim(group_name)//"x",&
      use_gatherv_mpio,dim2_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[i0_HSIZE_T,n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc) 

      if(allocated(B_hat_prev_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,B_hat_prev_arr,size(B_hat_prev_arr,1),n_particles_per_group,&
      trim(group_name)//"B_hat_prev",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc) 

      if(allocated(v_2d_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,v_2d_arr,size(v_2d_arr,1),n_particles_per_group,&
      trim(group_name)//"v",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc) 

      if(allocated(x_m_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,x_m_arr,size(x_m_arr,1),n_particles_per_group,&
      trim(group_name)//"x_m",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc) 

      if(allocated(Astar_m_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,Astar_m_arr,size(Astar_m_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_m",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc) 

      if(allocated(Astar_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,Astar_k_arr,size(Astar_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_k",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc) 

      if(allocated(Bn_k_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,Bn_k_arr,n_particles_per_group,trim(group_name)//"Bn_k",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      if(allocated(dBn_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,dBn_k_arr,size(dBn_k_arr,1),n_particles_per_group,&
      trim(group_name)//"dBn_k",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc) 

      if(allocated(Bnorm_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,Bnorm_k_arr,size(Bnorm_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Bnorm_k",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],use_hdf5_parallel_in=use_hdf5_parallel,&
      mpio_collective_in=collective_mpio_loc)

      if(allocated(E_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,E_k_arr,size(E_k_arr,1),n_particles_per_group,trim(group_name)//"E_k",&
      use_gatherv_mpio,dim2_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[i0_HSIZE_T,n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)
      
      if(allocated(dAstar_k_arr)) call HDF5_array3D_saving_native_or_gatherv(&
      file_id,dAstar_k_arr,size(dAstar_k_arr,1),size(dAstar_k_arr,2),&
      n_particles_per_group,trim(group_name)//"dAstar_k",use_gatherv_mpio,&
      dim3_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,i0_HSIZE_T,n_particles_offset],&
      use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=collective_mpio_loc)

      !> Write particle group attributes in HDF5 file, we assume that the attributes
      !> of the same group index for all tasks are equals. Therefore, we write the
      !> group attributes of only the master task. For HDF5-MPIO, the routines must
      !> be executed by all tasks for avoiding deadlocks
      if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.use_hdf5_parallel) then
        if(allocated(particle_type_str)) call HDF5_char_saving(file_id,&
        particle_type_str,trim(group_name)//"type")
        call HDF5_char_saving(file_id,sim%groups(ii)%ad%suffix,trim(group_name)//"adas_suffix")
        call HDF5_integer_saving(file_id,sim%groups(ii)%Z,trim(group_name)//"Z")
        call HDF5_real_saving(file_id,sim%groups(ii)%mass,trim(group_name)//"mass")
      endif
      !> deallocate structures
      call deallocate_particle_arrays(n_particles,i_elm_arr,i_life_arr,q_arr,&
      t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,vpar_arr,B_norm_arr,vpar_m_arr,&
      st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,Astar_m_arr,Astar_k_arr,&
      Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr)
    enddo
  else
    if(sim%my_id.eq.master_task) write(*,*) "WARNING: sim particle groups is not allocated!"
  endif
  !> cleanups
  if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.&
  use_hdf5_parallel) call HDF5_close(file_id)
  if(allocated(n_particles_loc))       deallocate(n_particles_loc)
  if(allocated(n_particles_glob))      deallocate(n_particles_glob)
  if(allocated(particle_displacement)) deallocate(particle_displacement)
  if(allocated(particle_type_str))     deallocate(particle_type_str)
end subroutine write_simulation_hdf5

!> Parallel read of particle HDF5 restart file
!> inputs:
!>   filename:       (character)(N) name of the output file
!>   sim:            (particle_sim) particle simulation object
!>   use_hdf5_access_properties: (logical)(optional) HDF5 file access property
!>                   if must be set to .false. for parallel I/O
!>                   default: .true. 
!>   mpi_comm_in:    (integer)(optional) MPI communicator identifier
!>   mpi_info_in:    (integer)(optional) MPI info structre for parallel IO
!> outputs:
!>   sim: (particle_sim) particle simulation object
subroutine read_simulation_hdf5(sim,filename,use_hdf5_access_properties,&
mpi_comm_in,mpi_info_in,test_in)
  use phys_module
  use mpi
  use hdf5,               only: HSIZE_T,HID_T
  use hdf5,               only: H5Gopen_f,H5Gget_info_f
  use hdf5,               only: H5Gclose_f
  use hdf5_io_module,     only: HDF5_open,HDF5_close
  use hdf5_io_module,     only: HDF5_char_reading
  use hdf5_io_module,     only: HDF5_allocatable_char_reading
  use hdf5_io_module,     only: HDF5_real_reading,HDF5_integer_reading
  use hdf5_io_module,     only: HDF5_get_dataset_rank_dims
  use hdf5_io_module,     only: HDF5_allocatable_array1D_reading_int
  use hdf5_io_module,     only: HDF5_allocatable_array1D_reading_r4
  use hdf5_io_module,     only: HDF5_allocatable_array1D_reading
  use hdf5_io_module,     only: HDF5_allocatable_array2D_reading
  use hdf5_io_module,     only: HDF5_allocatable_array3D_reading
  use mod_particle_types, only: particle_list_from_arrays
  use mod_particle_types, only: initialize_particle_list_to_zero
  use mod_particle_types, only: deallocate_particle_arrays
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_gc_vpar
  use mod_particle_types, only: particle_gc_Qin
  use mod_particle_types, only: particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_sim,   only: particle_sim, configure_particle_groups
  use mod_coronal,        only: coronal
  use mod_openadas,       only: read_adf11
  implicit none
  !> parameters:
  integer(HSIZE_T),parameter  :: i0_HSIZE_T=int(0,kind=HSIZE_T)
  integer(HSIZE_T),parameter  :: n1_HSIZE_T=int(-1,kind=HSIZE_T)
  !> inputs:
  character(len=*),intent(in) :: filename
  logical,intent(in),optional :: use_hdf5_access_properties
  integer,intent(in),optional :: mpi_comm_in,mpi_info_in
  logical,intent(in),optional :: test_in
  !> inputs-outputs:
  class(particle_sim),intent(inout) :: sim  
  !> variables:
  integer                                        :: ii,ierr,h5err,errorcode,n_groups_old
  integer                                        :: mpi_comm_loc,mpi_info_loc
  integer                                        :: storage_type,max_corder,rank
  integer(HID_T)                                 :: file_id,group_id
  integer(HSIZE_T)                               :: offset,n_particles_hsizet
  integer,          dimension(:),    allocatable :: n_particles_per_proc
  integer*4,        dimension(:),    allocatable :: i_elm_arr,i_life_arr
  integer*4,        dimension(:),    allocatable :: q_arr
  integer(HSIZE_T), dimension(:),    allocatable :: n_particles_tot,n_particles_max
  real*4,           dimension(:),    allocatable :: t_birth_arr
  real*8,           dimension(:),    allocatable :: weight_arr,v_1d_arr
  real*8,           dimension(:),    allocatable :: E_arr,mu_arr,vpar_arr
  real*8,           dimension(:),    allocatable :: B_norm_arr,vpar_m_arr,Bn_k_arr
  real*8,           dimension(:,:),  allocatable :: st_arr,x_arr,B_hat_prev_arr,v_2d_arr
  real*8,           dimension(:,:),  allocatable :: x_m_arr,Astar_m_arr,Astar_k_arr
  real*8,           dimension(:,:),  allocatable :: dBn_k_arr,Bnorm_k_arr,E_k_arr
  real*8,           dimension(:,:,:),allocatable :: dAstar_k_arr
  character(len=group_name_len)                  :: group_name
  character(len=:),                  allocatable :: particle_type_str
  logical                                        :: create_access_plist,test
  !> initialisation
  allocate(n_particles_per_proc(sim%n_cpu))
  !> set optional parameters
  create_access_plist = .false.
  if(present(use_hdf5_access_properties)) create_access_plist = .not.use_hdf5_access_properties
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_info_loc = mpi_info_in
  test = .false.; if(present(test_in)) test = test_in;
  !> open HDF5 file 
  call HDF5_open(filename,file_id,ierr,create_access_plist_in=create_access_plist,&
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  !> read the simulation time
  call HDF5_real_reading(file_id,sim%time,"/time",mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu)
  !> get number of groups
  call H5Gopen_f(file_id,"/groups/",group_id,h5err)   
  call H5Gget_info_f(group_id,storage_type,n_groups_old,max_corder,h5err) ! n_groups_old as n_groups can change in the future based on user input
  call H5Gclose_f(group_id,h5err)  
  !> check if n_part_groups from restart fits in n_part_group_max
  if (n_groups_old > n_part_groups_max) then
    write(*,*) "Error: n_part_groups being imported from restart exceeds n_part_groups_max"
    stop
  endif
  !> check if the number of groups requested differs from the number of groups in part_restart.h5
  !> will be removed after adding feature to support changing num groups
  if (n_groups_old /= n_part_groups) then
    write(*,*) "Error: mismatch between n_part_groups requested and that found in part_restart.h5"
    stop
  endif
  !> the line below should only be the case in certain unit tests, otherwise the
  !> allocation of groups and their configuration should be done in sim%initialize()
  if (.not. (allocated(sim%groups))) call sim%allocate_groups(n_part_groups)

<<<<<<< HEAD


!> Import all particles.
!> Reads the number of particles from a file, determines a
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine read_simulation_hdf5(sim, filename)
use mod_openadas, only: read_adf11
use mod_coronal
use phys_module
implicit none
type(particle_sim) , intent(inout) :: sim
character*(*)      , intent(in)  :: filename

integer                            :: my_id, n_cpu, ierr, rank
integer, allocatable, dimension(:) :: particles_per_proc ! particles per mpi proc for each group (group, processor_id)

! For HDF5 reading
integer(HID_T)    :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)    :: data_type
integer(HID_T)    :: group_id
integer(HID_T)    :: time_set_id
integer(HSIZE_T)  :: i_here
integer           :: n_here, n_part_groups_old
integer           :: storage_type, max_corder
character(len=12) :: group_name
character(len=3), allocatable :: part_groups_in_use_old(:)
character(len=3)  :: dropped_groups(n_part_groups_max) 
character(len=3)  :: group_id_tmp
integer           :: i, j, k, hdferr, n_alive, id
integer, allocatable :: n_alive_all(:)
logical           :: exists
! temp group attributes
integer            :: tmp_Z, tmp_n_particles, dropped_groups_counter
real*8             :: tmp_mass
character(len=3)   :: tmp_cs
character(len=8) :: tmp_ad_suffix 
character(len=particle_type_name_length)  :: tmp_part_type

type(c_ptr) :: p_ptr
integer*8, dimension(1:2) :: tmp, maxdims
real*8, dimension(:,:), allocatable :: real8_2D
integer*4, dimension(:), allocatable :: int4_1D
real*4, dimension(:), allocatable :: real4_1D
real*8, dimension(:), allocatable :: real8_1D

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! read the time
call HDF5_real_reading(file,sim%time,'/time')

! get old number of groups
call HDF5_integer_reading(file, n_part_groups_old, '/n_part_groups') 

! check if n_part_groups from restart fits in n_part_group_max
if (n_part_groups_old > n_part_groups_max) then
  write(*,*) "Error: n_part_groups being imported from restart exceeds n_part_groups_max"
  stop
endif

if (n_part_groups < n_part_groups_old) then
  write(*,*) "Error: n_part_groups specified is smaller than the number of groups from part_restart.h5. " // &
  "If you would like to remove particle groups, please specifiy the group ids to keep using part_group_in_use. "
else if (n_part_groups > n_part_groups_old) then
  write(*,*) "Warning: n_part_groups greater than number of groups from part_restart.h5, will be initializing new particle groups."
endif

allocate(part_groups_in_use_old(n_part_groups_old))
call HDF5_array1D_reading_char_len(file, part_groups_in_use_old, 3, '/part_groups_in_use')
if (part_groups_in_use(1) == 'non') part_groups_in_use(1:n_part_groups_old) = part_groups_in_use_old 

dropped_groups_counter = 0
do i=1, n_part_groups ! loop over groups in part_groups_in_use
  if (sim%my_id == 0) write(*,*) "Group ID to match: ", part_groups_in_use(i)
  
  do k=1, n_part_groups_old ! loop over the saved particle groups
    write(group_name,'(A,i0.3,A)') '/groups/', k, '/'
    call HDF5_char_reading(file,group_id_tmp,group_name//"id")

    if (group_id_tmp == part_groups_in_use(i)) then
      
      if (sim%my_id == 0) write(*,*) "Matching restart group found. Loading data."

      ! Open the dataset for x -----------------------------
      call h5dopen_f(file, group_name//"x", dset, hdferr)
    
      ! Open the file dataspace ----------------------------
      call h5dget_space_f(dset, file_space, hdferr)

      ! Reading particle group attributes ---------------------
      call HDF5_integer_reading(file,tmp_Z,group_name//"Z")
      call HDF5_real_reading(file,tmp_mass,group_name//"mass")
      call HDF5_char_reading(file,tmp_ad_suffix,group_name//"adas_suffix")
      call HDF5_char_reading(file,tmp_cs,group_name//"coupling_scheme")
      call HDF5_char_reading(file,tmp_part_type,group_name//"type")
      
      ! n_particles
      call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
      call h5sget_simple_extent_dims_f(file_space, tmp, maxdims, hdferr)
      tmp_n_particles = int(tmp(2),4)
      
      ! check if read attributes matches with specified attributes 

      if (particle_group_configs(i)%Z /= tmp_Z) then
        write(*,*) "IMPORT ERROR: Attribute 'Z' mismatch between namelist and imported data for group: ", group_id_tmp
        stop
      endif

      if (particle_group_configs(i)%mass /= tmp_mass) then
        write(*,*) "IMPORT ERROR: Attribute 'mass' mismatch between namelist and imported data for group: ", group_id_tmp
        stop
      endif

      if (trim(particle_group_configs(i)%coupling_scheme) /= trim(tmp_cs)) then
        write(*,*) "IMPORT ERROR: Attribute 'coupling_scheme' mismatch between namelist and imported data for group: ", group_id_tmp
        stop
      endif

      if (trim(particle_group_configs(i)%type) /= trim(tmp_part_type)) then
        write(*,*) "IMPORT ERROR: Attribute 'type' mismatch between namelist and imported data for group: ", group_id_tmp
        write(*,*) "config type: ", trim(particle_group_configs(i)%type)
        write(*,*) "load type: ", trim(tmp_part_type)
        stop
      endif

      if (particle_group_configs(i)%n_particles /= tmp_n_particles) then ! could allow to be changed in the future?
        write(*,*) "IMPORT ERROR: Attribute 'n_particles' mismatch between namelist and imported data for group: ", group_id_tmp
        stop
      endif

      ! attributes that are allowed to change
      if (trim(particle_group_configs(i)%atom_data_suffix) /= trim(tmp_ad_suffix)) then
        write(*,*) "IMPORT WARNING: Attribute 'atom_data_suffix' mismatch between namelist "
        write(*,*) " and imported data for group: ", group_id_tmp
        write(*,*) "Namelist: ", trim(particle_group_configs(k)%atom_data_suffix), " | Imported: ", trim(tmp_ad_suffix)
        write(*,*) "Will proceed with suffix given by namelist."
      endif
    
      ! Set attributes --------------------------------------
      sim%groups(i)%Z = tmp_Z
      sim%groups(i)%mass = tmp_mass
      sim%groups(i)%coupling_scheme = tmp_cs

      ! adas data
      sim%groups(i)%ad%suffix = tmp_ad_suffix

      if (len_trim(sim%groups(k)%ad%suffix) .gt. 0) then
        sim%groups(i)%ad = read_adf11(my_id,sim%groups(i)%ad%suffix)
        sim%groups(i)%cor = coronal(sim%groups(i)%ad)
      end if

      ! n_particles and load balancing
      sim%groups(i)%n_particles = tmp_n_particles
      call h5sclose_f(file_space, hdferr)
      call h5dclose_f(dset, hdferr)
      particles_per_proc(0)         = tmp_n_particles/n_cpu + modulo(tmp_n_particles, n_cpu)
      particles_per_proc(1:n_cpu-1) = tmp_n_particles/n_cpu
    
      i_here = sum(particles_per_proc(0:my_id-1))
      n_here = particles_per_proc(my_id)

      ! particle type
      ierr = 0
      select case (trim(tmp_part_type))
      case ('particle_kinetic')
        allocate(particle_kinetic::sim%groups(i)%particles(n_here), stat=ierr)
      case ('particle_kinetic_leapfrog')
        allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(n_here), stat=ierr)
      case ('particle_gc')
        allocate(particle_gc::sim%groups(i)%particles(n_here), stat=ierr)
      case ('particle_gc_vpar')
        allocate(particle_gc_vpar::sim%groups(i)%particles(n_here), stat=ierr)
      case ('particle_fieldline')
        allocate(particle_fieldline::sim%groups(i)%particles(n_here), stat=ierr)
      case ('particle_kinetic_relativistic')
        allocate(particle_kinetic_relativistic::sim%groups(i)%particles(n_here), stat=ierr)
      case ('particle_gc_relativistic')
        allocate(particle_gc_relativistic::sim%groups(i)%particles(n_here), stat=ierr)
      case default
        write(*,*) "error: missing type name declaration ", trim(tmp_part_type), " for read"
        call exit(1)
      end select
      if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
          "unable to allocate particles(", particles_per_proc(my_id), ")"
  
    
      ! Read particle specific properties -----------------
      ! x
      allocate(real8_2D(3,n_here))
      call HDF5_array2D_reading(file, real8_2D, group_name//"x",start=[0_HSIZE_T,i_here])
      do j=1,n_here
        sim%groups(i)%particles(j)%x = real8_2D(:,j)
      end do
      deallocate(real8_2D)
    
      ! st
      allocate(real8_2D(2,n_here))
      call HDF5_array2D_reading(file, real8_2D, group_name//"st",start=[0_HSIZE_T,i_here])
      do j=1,n_here
        sim%groups(i)%particles(j)%st = real8_2D(:,j)
      end do
      deallocate(real8_2D)
    
      ! weight
      allocate(real8_1D(n_here))
      call HDF5_array1D_reading(file, real8_1D, group_name//"weight",start=[i_here])
      do j=1,n_here
        sim%groups(i)%particles(j)%weight = real8_1D(j)
      end do
      deallocate(real8_1D)
      ! i_elm
      allocate(int4_1D(n_here))
      call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_elm", start=[i_here])
      do j=1,n_here
        sim%groups(i)%particles(j)%i_elm = int4_1D(j)
      end do
      deallocate(int4_1D)
    
      ! The following two are relatively new, so might not always be present
      ! i_life
      call h5lexists_f(file, group_name//"i_life", exists, ierr)
      if (exists) then
        allocate(int4_1D(n_here))
        int4_1D = 0 ! preset to 0 in case not present (i.e. old restart files)
        ! will still give a nasty error message, but should work,
        ! though who knows what hdf5 does with our array if reading fails...
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_life", start=[i_here])
        do j=1,n_here
          sim%groups(i)%particles(j)%i_life = int4_1D(j)
        end do
        deallocate(int4_1D)
      end if
    
      ! t_birth
      call h5lexists_f(file, group_name//"t_birth", exists, ierr)
      if (exists) then
        allocate(real4_1D(n_here))
        real4_1D = 0.0 ! preset to 0 in case not present (i.e. old restart files)
        call HDF5_array1D_reading_r4(file, real4_1D, group_name//"t_birth",start=[i_here])
        do j=1,n_here
          sim%groups(i)%particles(j)%t_birth = real4_1D(j)
        end do
        deallocate(real4_1D)
      end if
    
      select type (p => sim%groups(i)%particles)
      type is (particle_kinetic)
        ! v
        allocate(real8_2D(3,n_here))
        call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
        do j=1,n_here
          p(j)%v = real8_2D(:,j)
        end do
        deallocate(real8_2D)
    
        ! q
        allocate(int4_1D(n_here))
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
        do j=1,n_here
          p(j)%q = int4_1D(j)
        end do
        deallocate(int4_1D)
    
      type is (particle_kinetic_leapfrog)
        ! v
        allocate(real8_2D(3,n_here))
        call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
        do j=1,n_here
          p(j)%v = real8_2D(:,j)
        end do
        deallocate(real8_2D)
    
        ! q
        allocate(int4_1D(n_here))
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
        do j=1,n_here
          p(j)%q = int4_1D(j)
        end do
        deallocate(int4_1D)
    
      type is (particle_gc)
        ! E
        allocate(real8_1D(n_here))
        call HDF5_array1D_reading(file, real8_1D, group_name//"E",start=[i_here])
        do j=1,n_here
          p(j)%E = real8_1D(j)
        end do
        deallocate(real8_1D)
        ! mu
        allocate(real8_1D(n_here))
        call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
        do j=1,n_here
          p(j)%mu = real8_1D(j)
        end do
        deallocate(real8_1D)
        ! q
        allocate(int4_1D(n_here))
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
        do j=1,n_here
          p(j)%q = int4_1D(j)
        end do
        deallocate(int4_1D)
      
      type is (particle_gc_vpar)
        ! 
        allocate(real8_1D(n_here))
        call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar",start=[i_here])
        do j=1,n_here
          p(j)%Vpar = real8_1D(j)
        end do
        deallocate(real8_1D)
        ! mu
        allocate(real8_1D(n_here))
        call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
        do j=1,n_here
          p(j)%mu = real8_1D(j)
        end do
        deallocate(real8_1D)
        ! Bnorm
        allocate(real8_1D(n_here))
        call HDF5_array1D_reading(file, real8_1D, group_name//"B_norm",start=[i_here])
        do j=1,n_here
          p(j)%B_norm = real8_1D(j)
        end do
        deallocate(real8_1D)
        ! q
        allocate(int4_1D(n_here))
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
        do j=1,n_here
          p(j)%q = int4_1D(j)
        end do
        deallocate(int4_1D)
      
      type is (particle_fieldline)
        ! v
        allocate(real8_1D(n_here))
        call HDF5_array1D_reading(file, real8_1D, group_name//"v",start=[i_here])
        do j=1,n_here
          p(j)%v = real8_1D(j)
        end do
        deallocate(real8_1D)
    
      type is (particle_kinetic_relativistic)
    
        ! momenta [AMU*m/s]
        allocate(real8_2D(3,n_here))
        call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
        do j=1,n_here
          p(j)%p = real8_2D(:,j)
        end do
        deallocate(real8_2D)
    
        ! electric charge q
        allocate(int4_1D(n_here))
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
        do j=1,n_here
          p(j)%q = int4_1D(j)
        end do
        deallocate(int4_1D)    
    
      type is (particle_gc_relativistic)
    
        ! momenta (parallel momentum and magnetic moment)
        allocate(real8_2D(2,n_here))
        call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
        do j=1,n_here
          p(j)%p = real8_2D(:,j)
        end do
        deallocate(real8_2D)
    
        ! electric charge q
        allocate(int4_1D(n_here))
        call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
        do j=1,n_here
          p(j)%q = int4_1D(j)
        end do
        deallocate(int4_1D)    
    
      end select
    
      ! Check if the balance between processors is okay, by comparing the number of
      ! alive particles
      n_alive = count(sim%groups(i)%particles(:)%i_elm .gt. 0)
      allocate(n_alive_all(sim%n_cpu))
      call MPI_Gather(n_alive, 1, MPI_INTEGER, n_alive_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
      ! Check if the imbalance is not too great
      if (sim%my_id .eq. 0) then
        if (maxval(n_alive_all) .gt. minval(n_alive_all) * 1.5) then
          write(*,*) "WARNING: ", (maxval(n_alive_all)*100)/minval(n_alive_all), '% imbalance between CPUs, counts:'
          write(*,*) n_alive_all
        end if
      end if
      deallocate(n_alive_all)
    
    else 
      dropped_groups_counter = dropped_groups_counter + 1
      dropped_groups(dropped_groups_counter) = group_id_tmp
    endif ! group id match

  enddo ! n_part_groups_old
enddo ! groups in part_groups_in_use 

if (dropped_groups_counter > 0) then
  if (sim%my_id == 0) then
    write(*, "(1X,A, ' = ')", advance="no") "Imported particle groups now dropped: "
    do i = 1, dropped_groups_counter
      write(*, "(A, A)", advance="no") "'", trim(dropped_groups(i)) // "' "
    end do
    write(*,*)
  endif
endif


! dealloc temp arrays
if (allocated(part_groups_in_use_old)) deallocate(part_groups_in_use_old)

! Close everything else
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

=======
  do ii=1,n_groups_old
    !> read and load group datasets. We assume that the ith-particle group of 
    !> all tasks is defined by the same unique value stored in the hdf5 file
    ierr = 0; write(group_name,'(A,i0.3,A)') "/groups/",ii,"/";
    call HDF5_allocatable_char_reading(file_id,particle_type_str,trim(group_name)//"type")
    call HDF5_integer_reading(file_id,sim%groups(ii)%Z,trim(group_name)//"Z")
    call HDF5_real_reading(file_id,sim%groups(ii)%mass,trim(group_name)//"mass")
    call HDF5_char_reading(file_id,sim%groups(ii)%coupling_scheme,trim(group_name)//"coupling_scheme")
    call HDF5_char_reading(file_id,sim%groups(ii)%ad%suffix,trim(group_name)//"adas_suffix")
    if((len_trim(sim%groups(ii)%ad%suffix).gt.0).and.(.not.test)) then
      sim%groups(ii)%ad  = read_adf11(sim%my_id,sim%groups(ii)%ad%suffix)
      sim%groups(ii)%cor = coronal(sim%groups(ii)%ad) 
    endif
    !> compute the number of particles per processor and allocate particle array
    call HDF5_get_dataset_rank_dims(file_id,trim(group_name)//"i_elm",&
    rank,n_particles_tot,n_particles_max)
    n_particles_per_proc = int(n_particles_tot(1))/sim%n_cpu
    n_particles_per_proc(master_task+1) = int(n_particles_tot(1)) - &
    (sim%n_cpu-1)*n_particles_per_proc(master_task+1)
    offset = int(sum(n_particles_per_proc(1:sim%my_id)),kind=HSIZE_T)
    n_particles_hsizet = int(n_particles_per_proc(sim%my_id+1),kind=HSIZE_T)
    !> set n_particles for the group
    sim%groups(ii)%n_particles = int(n_particles_tot(1))
    !> allocate particle list and initialise to 0
    select case (trim(particle_type_str))
    case ("particle_kinetic")
      allocate(particle_kinetic::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_kinetic_leapfrog")
      allocate(particle_kinetic_leapfrog::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc")
      allocate(particle_gc::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc_vpar")
      allocate(particle_gc_vpar::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc_Qin")
      allocate(particle_gc_Qin::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_fieldline")
      allocate(particle_fieldline::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_kinetic_relativistic")
      allocate(particle_kinetic_relativistic::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc_relativistic")
      allocate(particle_gc_relativistic::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case default
      write(*,*) "Error: missing type name declaration ",trim(particle_type_str)," for reading: ABORT!"
      call MPI_Abort(mpi_comm_loc,errorcode,ierr)
    end select
    call initialize_particle_list_to_zero(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr)
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 1D array
    call HDF5_allocatable_array1D_reading_int(file_id,i_elm_arr,trim(group_name)//"i_elm",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading_int(file_id,i_life_arr,trim(group_name)//"i_life",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading_int(file_id,q_arr,trim(group_name)//"q",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: float 1D array
    call HDF5_allocatable_array1D_reading_r4(file_id,t_birth_arr,trim(group_name)//"t_birth",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: double 1D array
    call HDF5_allocatable_array1D_reading(file_id,weight_arr,trim(group_name)//"weight",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,v_1d_arr,trim(group_name)//"v",&
    reqdims_in=[n_particles_hsizet],start=[offset])   
    call HDF5_allocatable_array1D_reading(file_id,E_arr,trim(group_name)//"E",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,mu_arr,trim(group_name)//"mu",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,vpar_arr,trim(group_name)//"Vpar",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,B_norm_arr,trim(group_name)//"B_norm",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,vpar_m_arr,trim(group_name)//"Vpar_m",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,Bn_k_arr,trim(group_name)//"Bn_k",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 2D array
    call HDF5_allocatable_array2D_reading(file_id,st_arr,trim(group_name)//"st",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,x_arr,trim(group_name)//"x",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,B_hat_prev_arr,trim(group_name)//"B_hat_prev",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,x_m_arr,trim(group_name)//"x_m",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,Astar_m_arr,trim(group_name)//"Astar_m",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,Astar_k_arr,trim(group_name)//"Astar_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,dBn_k_arr,trim(group_name)//"dBn_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,Bnorm_k_arr,trim(group_name)//"Bnorm_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,E_k_arr,trim(group_name)//"E_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,v_2d_arr,trim(group_name)//"v",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 3D array
    call HDF5_allocatable_array3D_reading(file_id,dAstar_k_arr,trim(group_name)//"dAstar_k",&
    reqdims_in=[n1_HSIZE_T,n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,i0_HSIZE_T,offset])
    !> fill particle list from arrays
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    i_elm_arr=i_elm_arr,i_life_arr=i_life_arr,t_birth_arr=t_birth_arr,weight_arr=weight_arr,&
    x_arr=x_arr,st_arr=st_arr,q_arr=q_arr,v_1d_arr=v_1d_arr,E_arr=E_arr,mu_arr=mu_arr,&
    vpar_arr=vpar_arr,B_norm_arr=B_norm_arr,vpar_m_arr=vpar_m_arr,B_hat_prev_arr=B_hat_prev_arr,&
    v_2d_arr=v_2d_arr,x_m_arr=x_m_arr,Astar_m_arr=Astar_m_arr,Astar_k_arr=Astar_k_arr,&
    Bn_k_arr=Bn_k_arr,dBn_k_arr=dBn_k_arr,Bnorm_k_arr=Bnorm_k_arr,&
    E_k_arr=E_k_arr,dAstar_k_arr=dAstar_k_arr)
    !> deallocate structures
    call deallocate_particle_arrays(n_particles_per_proc(sim%my_id+1),i_elm_arr,&
    i_life_arr,q_arr,t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,vpar_arr,&
    B_norm_arr,vpar_m_arr,st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,&
    Astar_m_arr,Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr)   
  enddo
  !> clean-up
  call HDF5_close(file_id)
  if(allocated(n_particles_tot))      deallocate(n_particles_tot)
  if(allocated(n_particles_max))      deallocate(n_particles_max)
  if(allocated(n_particles_per_proc)) deallocate(n_particles_per_proc)
  if(allocated(particle_type_str))    deallocate(particle_type_str)
>>>>>>> 46d8f4abddc8a950131af8d646356ad3e65b16d2
end subroutine read_simulation_hdf5

!> !> Get '/time' from a file. Does not alter the units in any way
!> code works for jorek and particle restart files, and returns values in
!> different units for both
!> inputs:
!>   filename:            (character) name of the HDF5 file
!>   create_access_plist: (logical) create HDF parameter list (for MPIO)
!>                        default: false (H5P_DEFAULT_F), if true create
!>                        H5P_FILE_ACCESS_F
!>   mpi_comm:            (integer) MPI communicator identifier (for MPIO)
!>   mpi_info:            (integer) MPI parameter object (for MPIO)
!>   my_id:               (integer) MPI task identifier (for MPIO)
!>   n_cpu:               (integer) number of MPI task (for MPIO)
!> outputs:
!>   time:     (real8) restart simulation time
function get_simulation_hdf5_time(filename,use_hdf5_access_properties,&
mpi_comm_loc,mpi_info_loc,my_id,n_cpu) result(time)
  use mpi
  use hdf5,           only: HID_T
  use hdf5_io_module, only: HDF5_open,HDF5_close,HDF5_real_reading
  implicit none
  character(len=*),intent(in) :: filename
  integer,intent(in),optional :: mpi_comm_loc,mpi_info_loc
  integer,intent(in),optional :: my_id,n_cpu
  logical,intent(in),optional :: use_hdf5_access_properties
  real*8                      :: time
  integer                     :: h5err
  integer(HID_T)              :: file_id
  if(present(use_hdf5_access_properties).and.present(mpi_comm_loc).and.&
  present(mpi_info_loc).and.present(my_id).and.present(n_cpu)) then
    call HDF5_open(trim(filename),file_id,h5err,&
    create_access_plist_in=.not.use_hdf5_access_properties,&
    mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
    call HDF5_real_reading(file_id,time,"/time",mpi_rank=my_id,n_mpi_tasks=n_cpu)
  else
    call HDF5_open(trim(filename),file_id,h5err)
    call HDF5_real_reading(file_id,time,"/time")
  endif
  call HDF5_close(file_id)
end function get_simulation_hdf5_time

end module mod_particle_io
