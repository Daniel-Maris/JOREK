!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
implicit none
private
public write_simulation_hdf5,read_simulation_hdf5,get_simulation_hdf5_time
public read_simulation_hdf5_original,get_simulation_hdf5_time_original

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
  use mpi
  use hdf5,               only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F
  use hdf5,               only: H5P_DEFAULT_F,H5Gcreate_f,H5Gclose_f
  use hdf5_io_module,     only: HDF5_set_parallel_io_properties
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
  integer(HID_T)                :: file_id,group_id,transfer_property
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
  logical                                :: use_gatherv_mpio
  logical                                :: create_access_plist,collective_mpio_loc
  character(len=group_name_len)          :: group_name
  character(len=:),          allocatable :: particle_type_str

  !> preparation
  use_gatherv_mpio  = .true. !< use MPI gatherv for collecting all data for writing
  if(present(use_native_hdf5_mpio_in)) use_gatherv_mpio=.not.use_native_hdf5_mpio_in
  file_access_loc = H5F_ACC_TRUNC_F !< truncate the file by default
  if(present(file_access_in)) file_access_loc = file_access_in;
  create_access_plist = .false. !< serial access by default
  if(present(use_hdf5_access_properties)) create_access_plist = .not.use_hdf5_access_properties
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_info_loc = mpi_info_in
  collective_mpio_loc = .true. !< enable collective MPIO applications by default
  if(present(collective_mpio_in)) collective_mpio_loc = collective_mpio_in
  transfer_property = H5P_DEFAULT_F;
  if(.not.use_gatherv_mpio) call HDF5_set_parallel_io_properties(transfer_property,collective_mpio_loc)
  !> allocate the gatherv displacement array if required
  if(use_gatherv_mpio) allocate(particle_displacement(sim%n_cpu),source=0) 
  !> create the hdf5 file and the groups fields
  if(use_gatherv_mpio.and.(sim%my_id.eq.master_rank)) call HDF5_open_or_create(&
  trim(filename),file_id,ierr=h5err,file_access=file_access_loc)
  if(.not.use_gatherv_mpio) call HDF5_open_or_create(trim(filename),file_id,&
  ierr=h5err,file_access=file_access_loc,create_access_plist_in=create_access_plist,& 
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  if(h5err.gt.0) then
    if(sim%my_id.eq.master_task) write(*,*) "Failed to create or open the ",&
    filename," file: ",h5err,", ABORT!"
    call MPI_Abort(mpi_comm_loc,-1,ierr)
  endif
  if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.(.not.use_gatherv_mpio)) then
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
      if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.(.not.use_gatherv_mpio)) then
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
      transfer_property_in=transfer_property)
      
      if(allocated(i_life_arr)) call HDF5_array1D_saving_int_native_or_gatherv(&
      file_id,i_life_arr,n_particles_per_group,trim(group_name)//"i_life",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)
      
      if(allocated(q_arr)) call HDF5_array1D_saving_int_native_or_gatherv(&
      file_id,q_arr,n_particles_per_group,trim(group_name)//"q",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)
      
      if(allocated(t_birth_arr)) call HDF5_array1D_saving_r4_native_or_gatherv(&
      file_id,t_birth_arr,n_particles_per_group,trim(group_name)//"t_birth",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(weight_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,weight_arr,n_particles_per_group,trim(group_name)//"weight",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(v_1d_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,v_1d_arr,n_particles_per_group,trim(group_name)//"v",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(E_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,E_arr,n_particles_per_group,trim(group_name)//"E",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(mu_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,mu_arr,n_particles_per_group,trim(group_name)//"mu",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(vpar_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,vpar_arr,n_particles_per_group,trim(group_name)//"Vpar",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(B_norm_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,B_norm_arr,n_particles_per_group,trim(group_name)//"B_norm",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property) 

      if(allocated(vpar_m_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,vpar_m_arr,n_particles_per_group,trim(group_name)//"Vpar_m",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property) 

      if(allocated(st_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,st_arr,size(st_arr,1),n_particles_per_group,trim(group_name)//"st",&
      use_gatherv_mpio,dim2_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[i0_HSIZE_T,n_particles_offset],&
      transfer_property_in=transfer_property)  

      if(allocated(x_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,x_arr,size(x_arr,1),n_particles_per_group,trim(group_name)//"x",&
      use_gatherv_mpio,dim2_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[i0_HSIZE_T,n_particles_offset],&
      transfer_property_in=transfer_property) 

      if(allocated(B_hat_prev_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,B_hat_prev_arr,size(B_hat_prev_arr,1),n_particles_per_group,&
      trim(group_name)//"B_hat_prev",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property) 

      if(allocated(v_2d_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,v_2d_arr,size(v_2d_arr,1),n_particles_per_group,&
      trim(group_name)//"v",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property) 

      if(allocated(x_m_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,x_m_arr,size(x_m_arr,1),n_particles_per_group,&
      trim(group_name)//"x_m",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property) 

      if(allocated(Astar_m_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,Astar_m_arr,size(Astar_m_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_m",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property) 

      if(allocated(Astar_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,Astar_k_arr,size(Astar_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_k",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property) 

      if(allocated(Bn_k_arr)) call HDF5_array1D_saving_native_or_gatherv(&
      file_id,Bn_k_arr,n_particles_per_group,trim(group_name)//"Bn_k",&
      use_gatherv_mpio,dim1_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[n_particles_offset],&
      transfer_property_in=transfer_property)

      if(allocated(dBn_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,dBn_k_arr,size(dBn_k_arr,1),n_particles_per_group,&
      trim(group_name)//"dBn_k",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property) 

      if(allocated(Bnorm_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,Bnorm_k_arr,size(Bnorm_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Bnorm_k",use_gatherv_mpio,&
      dim2_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,n_particles_offset],transfer_property_in=transfer_property)

      if(allocated(E_k_arr)) call HDF5_array2D_saving_native_or_gatherv(&
      file_id,E_k_arr,size(E_k_arr,1),n_particles_per_group,trim(group_name)//"E_k",&
      use_gatherv_mpio,dim2_all_tasks=n_particles_glob(:,ii),&
      displs=particle_displacement,mpi_rank=sim%my_id,n_cpu=sim%n_cpu,&
      mpi_comm_loc=mpi_comm_loc,start=[i0_HSIZE_T,n_particles_offset],&
      transfer_property_in=transfer_property)
      
      if(allocated(dAstar_k_arr)) call HDF5_array3D_saving_native_or_gatherv(&
      file_id,dAstar_k_arr,size(dAstar_k_arr,1),size(dAstar_k_arr,2),&
      n_particles_per_group,trim(group_name)//"dAstar_k",use_gatherv_mpio,&
      dim3_all_tasks=n_particles_glob(:,ii),displs=particle_displacement,&
      mpi_rank=sim%my_id,n_cpu=sim%n_cpu,mpi_comm_loc=mpi_comm_loc,&
      start=[i0_HSIZE_T,i0_HSIZE_T,n_particles_offset],&
      transfer_property_in=transfer_property)

      !> Write particle group attributes in HDF5 file, we assume that the attributes
      !> of the same group index for all tasks are equals. Therefore, we write the
      !> group attributes of only the master task. For HDF5-MPIO, the routines must
      !> be executed by all tasks for avoiding deadlocks
      if((use_gatherv_mpio.and.(sim%my_id.eq.master_rank)).or.(.not.use_gatherv_mpio)) then
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
  (.not.use_gatherv_mpio)) call HDF5_close(file_id)
  if(allocated(n_particles_loc))       deallocate(n_particles_loc)
  if(allocated(n_particles_glob))      deallocate(n_particles_glob)
  if(allocated(particle_displacement)) deallocate(particle_displacement)
  if(allocated(particle_type_str))     deallocate(particle_type_str)
end subroutine write_simulation_hdf5

!> Original procedure for importing all particles.
!> Reads the number of particles from a file, determines a
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine read_simulation_hdf5_original(sim, filename, test_in)
use iso_c_binding, only: c_ptr
use mpi
use hdf5, only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F,H5P_FILE_ACCESS_F,H5F_ACC_RDONLY_F
use hdf5, only: h5open_f,h5pcreate_f,h5pset_fapl_mpio_f,h5fopen_f,h5pclose_f
use hdf5, only: h5gopen_f,h5gget_info_f,h5gclose_f,h5dopen_f,h5dget_space_f
use hdf5, only: h5sget_simple_extent_ndims_f,h5sget_simple_extent_dims_f
use hdf5, only: h5sclose_f,h5dclose_f,h5lexists_f,h5fclose_f,h5close_f
use hdf5_io_module, only: HDF5_real_reading,HDF5_allocatable_char_reading,HDF5_integer_reading
use hdf5_io_module, only: HDF5_array1D_reading,HDF5_array2D_reading,HDF5_array3D_reading
use hdf5_io_module, only: HDF5_array1D_reading_int,HDF5_array1D_reading_r4
use hdf5_io_module, only: HDF5_allocatable_array1D_reading
use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
use mod_particle_types, only: particle_gc,particle_gc_vpar
use mod_particle_types, only: particle_gc_Qin
use mod_particle_types, only: particle_fieldline
use mod_particle_types, only: particle_kinetic_relativistic
use mod_particle_types, only: particle_gc_relativistic
use mod_openadas,   only: read_adf11
use mod_coronal
use mod_particle_sim, only: particle_sim
type(particle_sim) , intent(inout)         :: sim
character*(*)      , intent(in)            :: filename
logical            , intent(in) , optional :: test_in

integer                            :: my_id, n_cpu, ierr, rank, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)               :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)               :: data_type
integer(HID_T)               :: group_id
integer(HID_T)               :: time_set_id
integer(HSIZE_T)             :: i_here
integer                      :: n_here
integer                      :: storage_type, max_corder
character(len=12)            :: group_name
character(len=:),allocatable :: particle_type_name,adas_suffix
integer                      :: i, j, n, hdferr, n_alive
integer, allocatable         :: n_alive_all(:)
logical                      :: exists,test

type(c_ptr) :: p_ptr
integer*8, dimension(1:2) :: tmp, maxdims
real*8,    dimension(:,:,:), allocatable :: real8_3D
real*8,    dimension(:,:),   allocatable :: real8_2D
integer*4, dimension(:),     allocatable :: int4_1D
real*4,    dimension(:),     allocatable :: real4_1D
real*8,    dimension(:),     allocatable :: real8_1D

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))
test = .false.; if(present(test_in)) test = test_in;

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! read the time
call HDF5_allocatable_array1D_reading(file,real8_1D,'/time')
sim%time = real8_1D(1); if(allocated(real8_1D)) deallocate(real8_1D);

! Get the number of groups
call h5gopen_f(file, '/groups/', group_id, hdferr)
call h5gget_info_f(group_id, storage_type, n, max_corder, hdferr)
call h5gclose_f(group_id, hdferr)

! Reallocate groups if necessary
if (allocated(sim%groups)) deallocate(sim%groups)
allocate(sim%groups(n))
do i=1,n
  ! Open the dataset for x
  write(group_name,'(A,i0.3,A)') '/groups/', i, '/'
  call h5dopen_f(file, group_name//"x", dset, hdferr)

  ! Open the file dataspace
  call h5dget_space_f(dset, file_space, hdferr)

  ! Get the number of particles
  call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
  call h5sget_simple_extent_dims_f(file_space, tmp, maxdims, hdferr)
  n_particles = int(tmp(2),4)
  call h5sclose_f(file_space, hdferr)
  call h5dclose_f(dset, hdferr)

  ! Divide particles over processors
  particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
  particles_per_proc(1:n_cpu-1) = n_particles/n_cpu
  i_here = sum(particles_per_proc(0:my_id-1))
  n_here = particles_per_proc(my_id)

  ! Get the particle type from the attribute
  call HDF5_allocatable_char_reading(file,particle_type_name,group_name//"type",&
  mpi_rank=0,n_mpi_tasks=1)

  ierr = 0
  select case (trim(particle_type_name))
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
  case ('particle_gc_Qin')
    allocate(particle_gc_Qin::sim%groups(i)%particles(n_here), stat=ierr)
  case default
    write(*,*) "error: missing type name declaration ", trim(particle_type_name), " for read"
    call exit(1)
  end select

  if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
      "unable to allocate particles(", particles_per_proc(my_id), ")"
  call HDF5_integer_reading(file,sim%groups(i)%Z,group_name//"Z")
  call HDF5_real_reading(file,sim%groups(i)%mass,group_name//"mass")
  call HDF5_allocatable_char_reading(file,adas_suffix,group_name//"adas_suffix",&
  mpi_rank=0,n_mpi_tasks=1)
  sim%groups(i)%ad%suffix = adas_suffix; 
  if ((len_trim(sim%groups(i)%ad%suffix) .gt. 0) .and. (.not.test)) then
    sim%groups(i)%ad = read_adf11(my_id,sim%groups(i)%ad%suffix)
    sim%groups(i)%cor = coronal(sim%groups(i)%ad)
  end if

  ! Read base particle attributes
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
    ! vpar 
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

  type is (particle_gc_Qin)
    ! vpar 
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar",start=[i_here])
    do j=1,n_here
      p(j)%vpar = real8_1D(j)
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
    ! x_m
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"x_m",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%x_m = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! Vpar_m
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar_m",start=[i_here])
    do j=1,n_here
      p(j)%vpar_m = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! Astar_m
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"Astar_m",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%Astar_m = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! Astar_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"Astar_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%Astar_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! dAstar_k
    allocate(real8_3D(3,3,n_here))
    call HDF5_array3D_reading(file, real8_3D, group_name//"dAstar_k",start=[0_HSIZE_T,0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%dAstar_k = real8_3D(:,:,j)
    end do
    deallocate(real8_3D)
    ! Bn_k 
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Bn_k",start=[i_here])
    do j=1,n_here
      p(j)%Bn_k = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! dBn_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"dBn_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%dBn_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! Bnorm_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"Bnorm_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%Bnorm_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! E_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"E_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%E_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
   
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
  if(allocated(particle_type_name)) deallocate(particle_type_name)
  if(allocated(adas_suffix))        deallocate(adas_suffix);
end do

! Close everything else
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine read_simulation_hdf5_original


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
  use mod_particle_sim,   only: particle_sim
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
  integer                                        :: ii,ierr,h5err,errorcode,n_groups
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
  call H5Gget_info_f(group_id,storage_type,n_groups,max_corder,h5err)
  call H5Gclose_f(group_id,h5err)  
  !> allocate simulation groups
  if(allocated(sim%groups)) deallocate(sim%groups); allocate(sim%groups(n_groups));
  do ii=1,n_groups
    !> read and load group datasets. We assume that the ith-particle group of 
    !> all tasks is defined by the same unique value stored in the hdf5 file
    ierr = 0; write(group_name,'(A,i0.3,A)') "/groups/",ii,"/";
    call HDF5_allocatable_char_reading(file_id,particle_type_str,trim(group_name)//"type")
    call HDF5_integer_reading(file_id,sim%groups(ii)%Z,trim(group_name)//"Z")
    call HDF5_real_reading(file_id,sim%groups(ii)%mass,trim(group_name)//"mass")
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
end subroutine read_simulation_hdf5

!> Original procedure for getting '/time' from a file. 
!> Does not alter the units in any way code works for jorek 
!> and particle restart files, and returns values in different units for both
function get_simulation_hdf5_time_original(filename) result(time)
  use hdf5, only: HID_T,H5F_ACC_RDONLY_F
  use hdf5, only: h5open_f,h5fopen_f,h5fclose_f,h5close_f
  use hdf5_io_module, only: HDF5_real_reading
  character*(*), intent(in)  :: filename
  real*8 :: time
  integer(HID_T) :: file
  integer :: hdferr
  call h5open_f(hdferr)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  call HDF5_real_reading(file,time,'/time')
  call h5fclose_f(file,hdferr)
  call h5close_f(hdferr)
end function get_simulation_hdf5_time_original

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
    write(*,*) "filename: ",trim(filename) !< DEBG DEBUG DEBUG
    call HDF5_real_reading(file_id,time,"/time",mpi_rank=my_id,n_mpi_tasks=n_cpu)
  else
    call HDF5_open(trim(filename),file_id,h5err)
    call HDF5_real_reading(file_id,time,"/time")
  endif
  call HDF5_close(file_id)
end function get_simulation_hdf5_time

end module mod_particle_io
