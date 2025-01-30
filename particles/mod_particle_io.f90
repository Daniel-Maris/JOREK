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


! allocate an empty array of a specified particle type for a specific particle group for each mpi process
subroutine allocate_particle_array(sim, group_num, particle_type_str, n_particles_per_proc, mpi_comm_loc)
  use mod_particle_sim, only: particle_sim
  use mod_particle_types
  use mpi

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
  use hdf5_io_module,     only: HDF5_array1D_saving_char
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

    ! Write n_part_groups and part_groups_in_use
    call HDF5_integer_saving(file_id,n_part_groups,'/n_part_groups')
    call HDF5_array1D_saving_char(file_id, part_groups_in_use, n_part_groups, '/part_groups_in_use')
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
      write(group_name,"(A,A,A)") "/groups/",trim(sim%groups(ii)%id),"/"
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
        call HDF5_char_saving(file_id,sim%groups(ii)%coupling_scheme,trim(group_name)//"coupling_scheme")
        call HDF5_char_saving(file_id,sim%groups(ii)%id,trim(group_name)//"id")
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
  use hdf5_io_module,     only: HDF5_array1D_reading_char
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
  integer                                        :: i,ii,ierr,h5err,errorcode
  integer                                        :: mpi_comm_loc,mpi_info_loc
  integer                                        :: storage_type,max_corder,rank
  integer                                        :: n_groups_old
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
  character(len=3)                               :: dropped_groups(n_part_groups_max)
  character(len=3)                               :: part_group_id
  character(len=3),                  allocatable :: part_groups_in_use_old(:) 

  ! temp group attributes
  integer            :: tmp_Z, dropped_groups_counter
  real*8             :: tmp_mass
  character(len=3)   :: tmp_cs
  character(len=8)   :: tmp_ad_suffix 

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

  ! checks and warnings about changes to n_part_groups

  !> reading in part_groups_in_use from part_restart.h5
  allocate(part_groups_in_use_old(n_groups_old))
  call HDF5_array1D_reading_char(file_id, part_groups_in_use_old, '/part_groups_in_use')
  !> overwrite part_groups_in_use with part_groups_in_use_old read in from 
  !> restart if it has not been manually defined
  if (part_groups_in_use(1) == 'non') part_groups_in_use(1:n_groups_old) = part_groups_in_use_old 

  !> the line below should only be the case in certain unit tests, otherwise the
  !> allocation of groups and their configuration should be done in sim%initialize()
  if (.not. (allocated(sim%groups))) call sim%allocate_groups(n_part_groups)

  if (sim%my_id == 0) then
    write(*,*) "" 
    write(*,*) "====== LOADING PARTICLE GROUPS ======"
  endif 

  dropped_groups_counter = 0
  do i=1, n_part_groups ! loop over groups in part_groups_in_use

    if (sim%my_id == 0) then
      write(*,*) ""
      write(*,*) " --- Group ID to match: ", part_groups_in_use(i), " --- "
    endif

    do ii=1,n_groups_old ! loop over the saved particle groups
      !> check for id match between id from part_groups_in_use and part_restart.h5
      ierr = 0; write(group_name,'(A,A,A)') "/groups/",trim(part_groups_in_use_old(ii)),"/";
      call HDF5_char_reading(file_id, part_group_id, trim(group_name)//"id")

      if (part_group_id == part_groups_in_use(i)) then
        if (sim%my_id == 0) write(*,*) "Matching restart group found. Loading data..."
        !> read and load group datasets. We assume that the ith-particle group of 
        !> all tasks is defined by the same unique value stored in the hdf5 file

        !> Reading particle group attributes ----------------
        call HDF5_allocatable_char_reading(file_id,particle_type_str,trim(group_name)//"type")
        call HDF5_integer_reading(file_id,tmp_Z,trim(group_name)//"Z")
        call HDF5_real_reading(file_id,tmp_mass,trim(group_name)//"mass")
        call HDF5_char_reading(file_id,tmp_cs,trim(group_name)//"coupling_scheme")
        call HDF5_char_reading(file_id,tmp_ad_suffix,trim(group_name)//"adas_suffix")
        call HDF5_get_dataset_rank_dims(file_id,trim(group_name)//"i_elm",&
        rank,n_particles_tot,n_particles_max) ! n_particles_tot

        !> Check if read attributes matches with specified attributes -----------
        if (particle_group_configs(i)%Z /= tmp_Z) then
          write(*,*) "IMPORT ERROR: Attribute 'Z' mismatch between namelist and imported data for group: ", part_group_id
          stop
        endif

        if (particle_group_configs(i)%mass /= tmp_mass) then
          write(*,*) "IMPORT ERROR: Attribute 'mass' mismatch between namelist and imported data for group: ", part_group_id
          stop
        endif

        if (trim(particle_group_configs(i)%coupling_scheme) /= trim(tmp_cs)) then
          write(*,*) "IMPORT ERROR: Attribute 'coupling_scheme' mismatch between namelist and imported data for group: ", part_group_id
          stop
        endif

        if (trim(particle_group_configs(i)%type) /= trim(particle_type_str)) then
          write(*,*) "IMPORT ERROR: Attribute 'type' mismatch between namelist and imported data for group: ", part_group_id
          write(*,*) "Namelist type: ", trim(particle_group_configs(i)%type)
          write(*,*) "Imported type: ", trim(particle_type_str)
          stop
        endif

        !> attributes that are allowed to change
        if (trim(particle_group_configs(i)%atom_data_suffix) /= trim(tmp_ad_suffix)) then
          write(*,*) "IMPORT WARNING: Attribute 'atom_data_suffix' mismatch between namelist and imported data for group: ", part_group_id
          write(*,*) " Namelist value: ", trim(particle_group_configs(i)%atom_data_suffix) 
          write(*,*) " Imported value: ", trim(tmp_ad_suffix)
          write(*,*) " Will proceed with suffix given by namelist."
        endif

        if (particle_group_configs(i)%n_particles /= n_particles_tot(1)) then
          write(*,*) "IMPORT WARNING: Attribute 'n_particles' mismatch between namelist and imported data for group: ", part_group_id
          write(*,*) " Namelist value: ", particle_group_configs(i)%n_particles 
          write(*,*) " Imported value: ", n_particles_tot(1)
          write(*,*) " Will proceed with namelist value. (Resizing particles array for group)"
        endif

        !> will probably need some more checks here when particle_group_configs(i)%n_particles < n_particles_tot(1)
    
        !> Setting attributes ------------------------
        sim%groups(i)%Z = tmp_Z
        sim%groups(i)%mass = tmp_mass
        sim%groups(i)%coupling_scheme = tmp_cs
        sim%groups(i)%n_particles = n_particles_tot(1)
  
        ! adas data
        sim%groups(i)%ad%suffix = tmp_ad_suffix
        if((len_trim(sim%groups(ii)%ad%suffix).gt.0).and.(.not.test)) then
          sim%groups(ii)%ad  = read_adf11(sim%my_id,sim%groups(ii)%ad%suffix)
          sim%groups(ii)%cor = coronal(sim%groups(ii)%ad) 
        endif

        !> compute the number of particles per processor and allocate particle array ---------
        n_particles_per_proc = int(n_particles_tot(1))/sim%n_cpu
        n_particles_per_proc(master_task+1) = int(n_particles_tot(1)) - &
        (sim%n_cpu-1)*n_particles_per_proc(master_task+1)
        offset = int(sum(n_particles_per_proc(1:sim%my_id)),kind=HSIZE_T)
        n_particles_hsizet = int(n_particles_per_proc(sim%my_id+1),kind=HSIZE_T)
        !> set n_particles for the group
        sim%groups(ii)%n_particles = int(n_particles_tot(1))

        !> allocate particle list and initialise to 0 -----------------
        call allocate_particle_array(sim, ii, particle_type_str, n_particles_per_proc(sim%my_id+1), mpi_comm_loc)

        call initialize_particle_list_to_zero(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr)


        !> Read particle base datasets from HDF5 and fill the particle lists: integer 1D array -------------
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
      else
        write(*,*) "WARNING: No matching group found in restart file for '", part_groups_in_use(i), "'"
        write(*,*) "A new group will be created and initialized. "
        
        ! calculate load balancing ------
        ! turn below into a function
        n_particles_per_proc = int(particle_group_configs(i)%n_particles)/sim%n_cpu
        n_particles_per_proc(master_task+1) = int(particle_group_configs(i)%n_particles) - &
        (sim%n_cpu-1)*n_particles_per_proc(master_task+1)

        ! getting particle type --------- 
        ! should clean up a bit here
        if(allocated(particle_type_str)) deallocate(particle_type_str)
        allocate(character(len=len(trim(particle_group_configs(i)%type)))::particle_type_str)
        particle_type_str = trim(particle_group_configs(i)%type)

        call allocate_particle_array(sim, i, particle_type_str, n_particles_per_proc(sim%my_id+1), mpi_comm_loc)
        
      endif ! (part_group_id == part_groups_in_use(i))
    enddo ! n_groups_old (particle groups in restart)
  enddo ! n_part_groups (particle groups requested in part_groups_in_use)

  if (sim%my_id == 0) then 
    write(*,*) ""
    write(*,*) "====== PARTICLE GROUPS LOADED ====== "
    write(*,*) ""
  endif
  ! inform about the groups that have been dropped


  !> clean-up
  call HDF5_close(file_id)
  if(allocated(n_particles_tot))      deallocate(n_particles_tot)
  if(allocated(n_particles_max))      deallocate(n_particles_max)
  if(allocated(n_particles_per_proc)) deallocate(n_particles_per_proc)
  if(allocated(particle_type_str))    deallocate(particle_type_str)
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
