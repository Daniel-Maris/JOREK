!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
implicit none
private
public write_simulation_hdf5, read_simulation_hdf5, get_simulation_hdf5_time

!> Module wide variables
integer,parameter :: master_task=0
integer,parameter :: group_name_len=12

contains

!> Export all particles using HDF5 IO
!> Parallel support to writing operations is provided
!> by MPI enabled HDF5 procedures.
!> inputs:
!>   sim:      (particle_sim) particle simulation object
!>   filename: (character)(N) name of the output file
subroutine write_simulation_hdf5(sim,filename,file_access_in,&
access_type_in,mpi_comm_in,mpi_info_in)
  use hdf5_io_module, only:
  use hdf5, only:
  use mpi
  use mod_particle_types, only: particle_arrays_from_list
  use mod_particle_types, only: deallocate_particle_arrays
  use mod_particle_sim,   only: particle_sim
  implicit none
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
  integer                       :: ii,jj,ierr,h5err,n_particles
  integer(HID_T)                :: file_id,group,group_id
  integer*4,dimension(:),    allocatable :: i_elm_arr,i_life_arr
  integer*4,dimension(:),    allocatable :: q_arr
  real*4,   dimension(:),    allocatable :: t_birth_arr
  real*8,   dimension(:),    allocatable :: weight_arr,v_1d_arr
  real*8,   dimension(:),    allocatable :: E_arr,mu_arr,vpar_arr
  real*8,   dimension(:),    allocatable :: B_norm_arr,vpar_m_arr,Bn_k_arr
  real*8,   dimension(:,:),  allocatable :: st_arr,x_arr,B_hat_prev_arr,v_2d_arr
  real*8,   dimension(:,:),  allocatable :: x_m_arr,Astar_m_arr,Astar_k_arr
  real*8,   dimension(:,:),  allocatable :: Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr
  real*8,   dimension(:,:,:),allocatable :: dAstar_k_arr
  character(len=group_name_len) :: group_name

  !> preparation
  file_access_loc = H5F_TRUNC_F !< truncate the file by default
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
  if(h5eff.gt.0) then
    if(sim%my_id.eq.master_task) write(*,*) "Failed to create or open the ",&
    filename," file: ",h5err,", ABORT!"
    call MPI_Abort(mpi_comm_loc,-1,ierr)
  endif
  call H5Gcreate_f(file_id,"/groups",group_id,h5err) !< create particle groups
  call H5Gclose_f(group_id,h5err) 
  call HDF5_real_saving(file_id,sim%time,"/time",sim%my_id,sim%n_cpus,&
  type_dataset_transfert_loc) !< write the time in HDF5 file
  !> check if loops are allocated and loop on them
  if(allocated(sim%groups)) then
    do ii=1,size(sim%groups,1)
      if(.not.allocated(sim%groups(ii)%particles)) then
        if(sim%my_id.eq.master_task) write(*,*) "WARNING: particle list N# ",&
        ii,"  not allocated: skip!"; cycle;
      endif
      !> create the HDF5 group for the particle list
      write(group_name,"(A,i0.3,A)") "/groups/",ii,"/"
      call H5Gcreate_f(file_id,group_name,group_id,h5err)
      call H5Gclose_f(group_id,h5err)
      !> reorganize and store the particle data in congruent arrays
      call particle_arrays_from_list(sim%groups(ii)%particles,n_particles,&
      i_elm_arr,i_life_arr,q_arr,t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,&
      vpar_arr,B_norm_arr,vpar_m_arr,st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,&
      Astar_m_arr,Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr)
      !> write data in HDF5 file
      !> deallocate structures
      call deallocate_particle_arrays(n_particles,i_elm_arr,i_life_arr,q_arr,&
      t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,vpar_arr,B_norm_arr,vpar_m_arr,&
      st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,Astar_m_arr,&
      Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr)
    enddo
  else
    if(sim%my_id.eq.master_rank) write(*,*) "WARNING: sim particle groups is not allocated!"
  endif
end subroutine write_simulation_hdf5

end module mod_particle_io
