module mod_hdf5_io_module_mpi_test
use mpi, only: MPI_INFO_NULL,MPI_COMM_WORLD
use mpi, only: MPI_Info_create,MPI_Info_free
use fruit
use fruit_mpi
use hdf5_io_module
implicit none
private
public :: run_fruit_hdf5_io_module_mpi
!> Variables --------------------------------------------------
integer,parameter :: master_rank=0
integer,parameter :: mpi_comm=MPI_COMM_WORLD
integer,parameter :: mpi_info=MPI_INFO_NULL
integer,parameter :: access_hdf5_parallel=1
character(len=13) :: filename_base="test_hdf5_file"
character(len=3)  :: extension=".h5"
character(len=1)  :: rank_format
integer :: rank_loc,n_tasks_loc,ifail_loc
integer :: mpi_comm_loc,mpi_info_loc
!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
subroutine run_fruit_hdf5_io_module_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  if(rank.eq.master_rank) write(*,'(/A)') " ... setting-up: hdf5 IO module mpi tests"
  call setup(rank,n_tasks,ifail) 
  if(rank.eq.master_rank) write(*,'(/A)') " ... running: hdf5 IO module mpi tests" 
  call run_test_case(test_create_hdf5_file,'test_create_hdf5_file')
  call run_test_case(test_open_hdf5_file,'test_open_hdf5_file')
  if(rank.eq.master_rank) write(*,'(/A)') " ... tearing-down: hdf5 IO module mpi tests" 
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_hdf5_io_module_mpi

!> Set-up and tear-down ---------------------------------------
subroutine setup(rank,n_tasks,ifail)
  implicit none
  integer,intent(in)    :: rank,n_tasks
  integer,intent(inout) :: ifail
  rank_loc = rank; n_tasks_loc = n_tasks; ifail_loc = ifail;
  mpi_info_loc = mpi_info; mpi_comm_loc = mpi_comm;
  if(mpi_info_loc.ne.MPI_INFO_NULL) call MPI_Info_create(mpi_info_loc,ifail_loc)
  rank_format = '1' 
  if(rank.gt.0) write(rank_format,'(I1)') int(log10(real(rank_loc,kind=8)))+1
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  implicit none
  integer,intent(in)    :: rank,n_tasks
  integer,intent(inout) :: ifail
  if(mpi_info_loc.ne.MPI_INFO_NULL) call MPI_Info_free(mpi_info_loc,ifail_loc)
  ifail = ifail_loc; rank_loc = -1; n_tasks_loc = -1;
end subroutine teardown
!> Tests ------------------------------------------------------
!> Test procedure for exclusively creating HDF5 files
subroutine test_create_hdf5_file()
  implicit none
  integer(HID_T)     :: file_id
  character(len=100) :: filename
  logical            :: file_exists
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  call HDF5_create(trim(filename),file_id,ifail_loc)
  call HDF5_close(file_id)
  inquire(file=trim(filename),exist=file_exists)
  call assert_true(file_exists,"Error test create HDF5 file posix: file "//&
  trim(filename)//" not created!")
  call remove_file(filename,file_exists_in=file_exists)
  filename = trim(filename_base)//trim(extension)
  call HDF5_create(filename,file_id,access_type_in=access_hdf5_parallel,&
  mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_close(file_id)
  inquire(file=trim(filename),exist=file_exists)
  call assert_true(file_exists,"Error test create HDF5 file access FILE_ACCESS: file "//&
  trim(filename)//" not created!")
  call remove_file(filename,file_exists_in=file_exists)
end subroutine test_create_hdf5_file

!> Test procedure for exclusively opening HDF5 file
subroutine test_open_hdf5_file()
  implicit none
  integer(HID_T)     :: file_id
  character(len=100) :: filename
  logical            :: file_exists
  write(filename,'(A,A,I'//trim(rank_format)//',A)') & 
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_create(filename,file_id,ifail_loc); call HDF5_close(file_id); 
  call HDF5_open(filename,file_id,ifail_loc); call HDF5_close(file_id);
  call assert_equals(ifail_loc,0,"Error test open HDF5 file posix: file "//&
  trim(filename)//" not opened!"); call remove_file(filename);
  filename = trim(filename_base)//trim(extension); ifail_loc=0;
  call HDF5_create(filename,file_id,access_type_in=access_hdf5_parallel,&
  mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc); call HDF5_close(file_id); 
  call HDF5_create(filename,file_id,ierr=ifail_loc,access_type_in=access_hdf5_parallel,&
  mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc); call HDF5_close(file_id);
  call assert_equals(ifail_loc,0,"Error test open HDF5 file access FILE_ACCESS: file "//&
  trim(filename)//" not opened!")
  call remove_file(filename)
end subroutine test_open_hdf5_file

!> Tools ------------------------------------------------------
subroutine remove_file(filename,file_exists_in)
  implicit none
  character(len=*),intent(in) :: filename
  logical,optional,intent(in) :: file_exists_in
  logical                     :: file_exists
  if(present(file_exists_in)) then
    file_exists = file_exists_in
  else
    inquire(file=filename,exist=file_exists)
  endif
  if(file_exists) call system("rm -rf "//trim(filename))
end subroutine remove_file

!> ------------------------------------------------------------
end module mod_hdf5_io_module_mpi_test
