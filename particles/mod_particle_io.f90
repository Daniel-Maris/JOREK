!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
implicit none
private
public write_simulation_hdf5, read_simulation_hdf5, get_simulation_hdf5_time

!> Module wide variables
integer,parameter :: master_task=0
integer,parameter :: n_cpu_1=1
integer,parameter :: group_name_len=12

contains

!> Export all particles using HDF5 IO
!> Parallel support to writing operations is provided
!> by MPI enabled HDF5 procedures.
!> inputs:
!>   sim:                       (particle_sim) particle simulation object
!>   filename:                  (character)(N) name of the output file
!>   file_access_in:            (integer)(optional) define the file access type:
!>                              default: H5F_ACC_TRUNC_F
!>   access_type_in:            (integer)(optional) HDF5 file access property
!>                              default: H5P_FILE_ACCESS_F
!>   type_dataset_transfert_in: (integer)(optional) method for HDF5 parallel IO
!>                              default: 1 (parallel collective IO)
!>   mpi_comm_in:               (integer)(optional) MPI communicator identifier
!>   mpi_info_in:               (integer)(optional) MPI info structre for parallel IO
subroutine write_simulation_hdf5(sim,filename,file_access_in,&
access_type_in,type_dataset_transfert_in,mpi_comm_in,mpi_info_in)
  use mpi
  use hdf5,               only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F
  use hdf5,               only: H5Gcreate_f,H5Gclose_f
  use hdf5_io_module,     only: HDF5_open_or_create,HDF5_close
  use hdf5_io_module,     only: HDF5_integer_saving 
  use hdf5_io_module,     only: HDF5_real_saving,HDF5_char_saving
  use hdf5_io_module,     only: HDF5_array1D_saving_int
  use hdf5_io_module,     only: HDF5_array1D_saving_r4
  use hdf5_io_module,     only: HDF5_array1D_saving
  use hdf5_io_module,     only: HDF5_array2D_saving
  use hdf5_io_module,     only: HDF5_array3D_saving
  use mod_particle_types, only: particle_arrays_from_list
  use mod_particle_types, only: deallocate_particle_arrays
  use mod_particle_sim,   only: particle_sim
  implicit none
  !> parameters
  integer(HSIZE_T),parameter    :: i0_HSIZE_T=int(0,kind=HSIZE_T)
  !> input variables
  type(particle_sim),intent(in) :: sim
  character(len=*),  intent(in) :: filename 
  integer, intent(in), optional :: file_access_in,access_type_in
  integer, intent(in), optional :: mpi_comm_in,mpi_info_in
  integer, intent(in), optional :: type_dataset_transfert_in
  !> variables
  integer                       :: file_access_loc,access_type_loc
  integer                       :: mpi_comm_loc,mpi_info_loc
  integer                       :: type_dataset_transfert_loc
  integer                       :: ii,ierr,h5err,n_groups,n_particles
  integer                       :: n_particles_per_group
  integer(HID_T)                :: file_id,group_id
  integer(HSIZE_T)              :: n_particles_offset
  integer,  dimension(:),    allocatable :: n_particles_loc
  integer*4,dimension(:),    allocatable :: i_elm_arr,i_life_arr
  integer*4,dimension(:),    allocatable :: q_arr
  integer,  dimension(:,:),  allocatable :: n_particles_glob
  real*4,   dimension(:),    allocatable :: t_birth_arr
  real*8,   dimension(:),    allocatable :: weight_arr,v_1d_arr
  real*8,   dimension(:),    allocatable :: E_arr,mu_arr,vpar_arr
  real*8,   dimension(:),    allocatable :: B_norm_arr,vpar_m_arr,Bn_k_arr
  real*8,   dimension(:,:),  allocatable :: st_arr,x_arr,B_hat_prev_arr,v_2d_arr
  real*8,   dimension(:,:),  allocatable :: x_m_arr,Astar_m_arr,Astar_k_arr
  real*8,   dimension(:,:),  allocatable :: dBn_k_arr,Bnorm_k_arr,E_k_arr
  real*8,   dimension(:,:,:),allocatable :: dAstar_k_arr
  character(len=group_name_len)          :: group_name
  character(len=:),          allocatable :: particle_type_str

  !> preparation
  file_access_loc = H5F_ACC_TRUNC_F !< truncate the file by default
  if(present(file_access_in)) file_access_loc = file_access_in;
  access_type_loc = 1 !< parallel access by default
  if(present(access_type_in)) access_type_loc = access_type_in
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_info_loc = mpi_info_in
  type_dataset_transfert_loc = 1 !< enable parallel dataset transfert by default
  if(present(type_dataset_transfert_in)) type_dataset_transfert_loc = type_dataset_transfert_in
  !> create the hdf5 file and the groups fields
  call HDF5_open_or_create(filename,file_id,h5err,&
  file_access=file_access_loc,access_type_in=access_type_loc,& 
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  if(h5err.gt.0) then
    if(sim%my_id.eq.master_task) write(*,*) "Failed to create or open the ",&
    filename," file: ",h5err,", ABORT!"
    call MPI_Abort(mpi_comm_loc,-1,ierr)
  endif
  call H5Gcreate_f(file_id,"/groups",group_id,h5err) !< create particle groups
  call H5Gclose_f(group_id,h5err) 
  call HDF5_real_saving(file_id,sim%time,"/time",sim%my_id,sim%n_cpu,&
  type_dataset_transfert_loc) !< write the time in HDF5 file
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
      n_particles_offset    = int(sum(n_particles_glob(1:sim%my_id,ii)),kind=HSIZE_T)
      !> create the HDF5 group for the particle list
      write(group_name,"(A,i0.3,A)") "/groups/",ii,"/"
      call H5Gcreate_f(file_id,trim(group_name),group_id,h5err)
      call H5Gclose_f(group_id,h5err)
      !> reorganize and store the particle data in congruent arrays
      call particle_arrays_from_list(sim%groups(ii)%particles,n_particles,&
      i_elm_arr,i_life_arr,q_arr,t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,&
      vpar_arr,B_norm_arr,vpar_m_arr,st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,&
      Astar_m_arr,Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr,&
      particle_type_str)
      !> write data in HDF5 file
      if(allocated(i_elm_arr)) call HDF5_array1D_saving_int(file_id,i_elm_arr,&
      n_particles_per_group,trim(group_name)//"i_elm",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(i_life_arr)) call HDF5_array1D_saving_int(file_id,i_life_arr,&
      n_particles_per_group,trim(group_name)//"i_life",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(q_arr)) call HDF5_array1D_saving_int(file_id,q_arr,&
      n_particles_per_group,trim(group_name)//"q",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(t_birth_arr)) call HDF5_array1D_saving_r4(file_id,t_birth_arr,&
      n_particles_per_group,trim(group_name)//"t_birth",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(weight_arr)) call HDF5_array1D_saving(file_id,weight_arr,&
      n_particles_per_group,trim(group_name)//"weight",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(v_1d_arr)) call HDF5_array1D_saving(file_id,v_1d_arr,&
      n_particles_per_group,trim(group_name)//"v",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)     
      if(allocated(E_arr)) call HDF5_array1D_saving(file_id,E_arr,&
      n_particles_per_group,trim(group_name)//"E",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(mu_arr)) call HDF5_array1D_saving(file_id,mu_arr,&
      n_particles_per_group,trim(group_name)//"mu",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(vpar_arr)) call HDF5_array1D_saving(file_id,vpar_arr,&
      n_particles_per_group,trim(group_name)//"Vpar",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(B_norm_arr)) call HDF5_array1D_saving(file_id,B_norm_arr,&
      n_particles_per_group,trim(group_name)//"B_norm",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(vpar_m_arr)) call HDF5_array1D_saving(file_id,vpar_m_arr,&
      n_particles_per_group,trim(group_name)//"Vpar_m",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(st_arr)) call HDF5_array2D_saving(file_id,st_arr,size(st_arr,1),&
      n_particles_per_group,trim(group_name)//"st",start=[i0_HSIZE_T,&
      n_particles_offset],type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(x_arr)) call HDF5_array2D_saving(file_id,x_arr,size(x_arr,1),&
      n_particles_per_group,trim(group_name)//"x",start=[i0_HSIZE_T,&
      n_particles_offset],type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(B_hat_prev_arr)) call HDF5_array2D_saving(&
      file_id,B_hat_prev_arr,size(B_hat_prev_arr,1),n_particles_per_group,&
      trim(group_name)//"B_hat_prev",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(v_2d_arr)) call HDF5_array2D_saving(&
      file_id,v_2d_arr,size(v_2d_arr,1),n_particles_per_group,&
      trim(group_name)//"v",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(x_m_arr)) call HDF5_array2D_saving(&
      file_id,x_m_arr,size(x_m_arr,1),n_particles_per_group,&
      trim(group_name)//"x_m",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(Astar_m_arr)) call HDF5_array2D_saving(&
      file_id,Astar_m_arr,size(Astar_m_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_m",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(Astar_k_arr)) call HDF5_array2D_saving(&
      file_id,Astar_k_arr,size(Astar_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_k",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(Bn_k_arr)) call HDF5_array1D_saving(file_id,Bn_k_arr,&
      n_particles_per_group,trim(group_name)//"Bn_k",start=[n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(dBn_k_arr)) call HDF5_array2D_saving(file_id,&
      dBn_k_arr,size(dBn_k_arr,1),n_particles_per_group,&
      trim(group_name)//"dBn_k",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(Bnorm_k_arr)) call HDF5_array2D_saving(file_id,&
      Bnorm_k_arr,size(Bnorm_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Bnorm_k",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(E_k_arr)) call HDF5_array2D_saving(&
      file_id,E_k_arr,size(E_k_arr,1),n_particles_per_group,&
      trim(group_name)//"E_k",start=[i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(dAstar_k_arr)) call HDF5_array3D_saving(file_id,&
      dAstar_k_arr,size(dAstar_k_arr,1),size(dAstar_k_arr,2),n_particles_per_group,&
      trim(group_name)//"dAstar_k",start=[i0_HSIZE_T,i0_HSIZE_T,n_particles_offset],&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      if(allocated(particle_type_str)) call HDF5_char_saving(file_id,particle_type_str,&
      trim(group_name)//"type",mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu,&
      type_dataset_transfert_in=type_dataset_transfert_loc)
      call HDF5_char_saving(file_id,sim%groups(ii)%ad%suffix,trim(group_name)//"adas_suffix",&
      mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu,type_dataset_transfert_in=type_dataset_transfert_loc)
      call HDF5_integer_saving(file_id,sim%groups(ii)%Z,trim(group_name)//"Z",&
      mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu,type_dataset_transfert_in=type_dataset_transfert_loc)
      call HDF5_real_saving(file_id,sim%groups(ii)%mass,trim(group_name)//"mass",&
      mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu,type_dataset_transfert_in=type_dataset_transfert_loc)
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
  call HDF5_close(file_id)
  if(allocated(n_particles_loc))   deallocate(n_particles_loc)
  if(allocated(n_particles_glob))  deallocate(n_particles_glob)
  if(allocated(particle_type_str)) deallocate(particle_type_str)
end subroutine write_simulation_hdf5

!> Parallel read of particle HDF5 restart file
!> inputs:
!>   filename:       (character)(N) name of the output file
!>   sim:            (particle_sim) particle simulation object
!>   access_type_in: (integer)(optional) HDF5 file access property
!>                   default: H5P_FILE_ACCESS_F
!>   mpi_comm_in:    (integer)(optional) MPI communicator identifier
!>   mpi_info_in:    (integer)(optional) MPI info structre for parallel IO
!> outputs:
!>   sim: (particle_sim) particle simulation object
subroutine read_simulation_hdf5(sim,filename,access_type_in,mpi_comm,mpi_info)
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
  integer,intent(in),optional :: access_type_in
  integer,intent(in),optional :: mpi_comm,mpi_info
  !> inputs-outputs:
  class(particle_sim),intent(inout) :: sim  
  !> variables:
  integer                                        :: ii,ierr,h5err,errorcode,n_groups
  integer                                        :: access_type_loc
  integer                                        :: mpi_comm_loc,mpi_info_loc
  integer                                        :: storage_type,max_corder,rank
  integer(HID_T)                                 :: file_id,group_id
  integer(HSIZE_T)                               :: offset,n_particles_hsizet
  integer,          dimension(:),    allocatable :: n_particles_per_proc
  integer*4,        dimension(:),    allocatable :: array1D_int
  integer(HSIZE_T), dimension(:),    allocatable :: n_particles_tot,n_particles_max
  real*4,           dimension(:),    allocatable :: array1D_r4
  real*8,           dimension(:),    allocatable :: array1D_r8
  real*8,           dimension(:,:),  allocatable :: array2D_r8
  real*8,           dimension(:,:,:),allocatable :: array3D_r8
  character(len=group_name_len)                  :: group_name
  character(len=:),                  allocatable :: particle_type_str
  !> initialisation
  allocate(n_particles_per_proc(sim%n_cpu))
  !> set optional parameters
  access_type_loc = 1
  if(present(access_type_in)) access_type_loc = access_type_in
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm)) mpi_comm_loc = mpi_comm
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info)) mpi_info_loc = mpi_info
  !> open HDF5 file 
  call HDF5_open(filename,file_id,ierr,access_type_in=access_type_loc,&
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
    !> read and load group datasets
    ierr = 0; write(group_name,'(A,i0.3,A)') "/groups/",ii,"/";
    call HDF5_allocatable_char_reading(file_id,particle_type_str,&
    trim(group_name)//"type",mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu)
    call HDF5_integer_reading(file_id,sim%groups(ii)%Z,trim(group_name)//"Z",&
    mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu)
    call HDF5_real_reading(file_id,sim%groups(ii)%mass,trim(group_name)//"mass",&
    mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu)
    call HDF5_char_reading(file_id,sim%groups(ii)%ad%suffix,trim(group_name)//&
    "adas_suffix",mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu)
    if(len_trim(sim%groups(ii)%ad%suffix).gt.0) then
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
    !> allocate particle list
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
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 1D array
    call HDF5_allocatable_array1D_reading_int(file_id,array1D_int,trim(group_name)//"i_elm",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    i_elm_arr=array1D_int)
    call HDF5_allocatable_array1D_reading_int(file_id,array1D_int,trim(group_name)//"i_life",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    i_life_arr=array1D_int)
    call HDF5_allocatable_array1D_reading_int(file_id,array1D_int,trim(group_name)//"q",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    q_arr=array1D_int)
    !> Read particle base datasets from HDF5 and fill the particle lists: float 1D array
    call HDF5_allocatable_array1D_reading_r4(file_id,array1D_r4,trim(group_name)//"t_birth",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: double 1D array
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    t_birth_arr=array1D_r4)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"weight",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    weight_arr=array1D_r8)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"v",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    v_1d_arr=array1D_r8)   
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"E",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    E_arr=array1D_r8)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"mu",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    mu_arr=array1D_r8)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"Vpar",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    vpar_arr=array1D_r8)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"B_norm",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    B_norm_arr=array1D_r8)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"Vpar_m",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    vpar_m_arr=array1D_r8)
    call HDF5_allocatable_array1D_reading(file_id,array1D_r8,trim(group_name)//"Bn_k",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    Bn_k_arr=array1D_r8)
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 2D array
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"st",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    st_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"x",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    x_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"B_hat_prev",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    B_hat_prev_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"x_m",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    x_m_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"Astar_m",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    Astar_m_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"Astar_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    Astar_k_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"dBn_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    dBn_k_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"Bnorm_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    Bnorm_k_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"E_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    E_k_arr=array2D_r8)
    call HDF5_allocatable_array2D_reading(file_id,array2D_r8,trim(group_name)//"v",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    v_2d_arr=array2D_r8)
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 3D array
    call HDF5_allocatable_array3D_reading(file_id,array3D_r8,trim(group_name)//"dAstar_k",&
    reqdims_in=[n1_HSIZE_T,n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,i0_HSIZE_T,offset])
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),&
    sim%groups(ii)%particles,ierr,dAstar_k_arr=array3D_r8)
  enddo
  !> clean-up
  call HDF5_close(file_id)
  if(allocated(n_particles_tot))      deallocate(n_particles_tot)
  if(allocated(n_particles_max))      deallocate(n_particles_max)
  if(allocated(n_particles_per_proc)) deallocate(n_particles_per_proc)
  if(allocated(array1D_int))          deallocate(array1D_int)
  if(allocated(array1D_r4))           deallocate(array1D_r4)
  if(allocated(array1D_r8))           deallocate(array1D_r8)
  if(allocated(array2D_r8))           deallocate(array2D_r8)
  if(allocated(array3D_r8))           deallocate(array3D_r8)
  if(allocated(particle_type_str))    deallocate(particle_type_str)
end subroutine read_simulation_hdf5

!> !> Get '/time' from a file. Does not alter the units in any way
!> code works for jorek and particle restart files, and returns values in
!> different units for both
!> inputs:
!>   filename: (character) name of the HDF5 file
!> outputs:
!>   time:     (real8) restart simulation time
function get_simulation_hdf5_time(filename,access_type_in,mpi_comm_in,&
mpi_info_in,my_id_in,n_cpu_in) result(time)
  use mpi
  use hdf5, only: HID_T
  use hdf5_io_module, only: HDF5_open,HDF5_close,HDF5_real_reading
  implicit none
  character(len=*),intent(in) :: filename
  integer,intent(in),optional :: access_type_in,mpi_comm_in,mpi_info_in 
  integer,intent(in),optional :: my_id_in,n_cpu_in
  real*8                      :: time
  integer                     :: h5err,access_type_loc,mpi_comm_loc
  integer                     :: mpi_info_loc,my_id,n_cpu
  integer(HID_T)              :: file_id
  !> set optional parameters
  access_type_loc = 1
  if(present(access_type_in)) access_type_loc = access_type_in
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_comm_loc = mpi_comm_in
  my_id = master_task; if(present(my_id_in)) my_id = my_id_in;
  n_cpu = n_cpu_1;     if(present(n_cpu_in)) n_cpu = n_cpu_in;
  call HDF5_open(trim(filename),file_id,h5err,access_type_in=access_type_loc,&
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_real_reading(file_id,time,"/time",mpi_rank=my_id,n_mpi_tasks=n_cpu)
  call HDF5_close(file_id)
end function get_simulation_hdf5_time

end module mod_particle_io
