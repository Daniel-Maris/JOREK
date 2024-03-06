module mod_hdf5_io_module_mpi_test
use mpi,           only: MPI_INFO_NULL,MPI_COMM_WORLD
use mpi,           only: MPI_Info_create,MPI_Info_free
use mod_pcg32_rng, only: pcg32_rng
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
integer,parameter :: type_dataset_transfert_mpi=1
integer,parameter :: ndims_tot=5
integer,dimension(ndims_tot),parameter :: n_elements=[13,20,18,7,4]
real*8,parameter  :: tol_r8=1.d-16
character(len=14) :: filename_base="test_hdf5_file"
character(len=3)  :: extension=".h5"
character(len=1)  :: rank_format
integer           :: rank_loc,n_tasks_loc,ifail_loc
integer           :: mpi_comm_loc,mpi_info_loc
real*8,dimension(n_elements(1),n_elements(2),n_elements(3),&
n_elements(4),n_elements(5)) :: array_sol
type(pcg32_rng)   :: rng
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
  call run_test_case(test_create_open_hdf5_file,'test_create_open_hdf5_file')
  call run_test_case(test_HDF5_array1D_saving_int,"test_HDF5_array1D_saving_int")
  call run_test_case(test_HDF5_array2D_saving_int,"test_HDF5_array2D_saving_int")
  call run_test_case(test_HDF5_array1D_saving_r4,"test_HDF5_array1D_saving_r4")
  call run_test_case(test_HDF5_real_saving,"test_HDF5_real_saving")
  call run_test_case(test_HDF5_array1D_saving_r8,"test_HDF5_array1D_saving_r8")
  call run_test_case(test_HDF5_array2D_saving_r8,"test_HDF5_array2D_saving_r8")
  call run_test_case(test_HDF5_array3D_saving_r8,"test_HDF5_array3D_saving_r8")
  call run_test_case(test_HDF5_array4D_saving_r8,"test_HDF5_array4D_saving_r8")
  call run_test_case(test_HDF5_array5D_saving_r8,"test_HDF5_array5D_saving_r8")
  if(rank.eq.master_rank) write(*,'(/A)') " ... tearing-down: hdf5 IO module mpi tests" 
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_hdf5_io_module_mpi

!> Set-up and tear-down ---------------------------------------
subroutine setup(rank,n_tasks,ifail)
  use mod_random_seed, only: random_seed
  implicit none
  integer,intent(in)    :: rank,n_tasks
  integer,intent(inout) :: ifail
  integer               :: ii,jj,kk,pp
  rank_loc = rank; n_tasks_loc = n_tasks; ifail_loc = ifail;
  mpi_info_loc = mpi_info; mpi_comm_loc = mpi_comm;
  if(mpi_info_loc.ne.MPI_INFO_NULL) call MPI_Info_create(mpi_info_loc,ifail_loc)
  rank_format = '1' 
  if(rank.gt.0) write(rank_format,'(I1)') int(log10(real(rank_loc,kind=8)))+1
  !> setup the rng
  call rng%initialize(n_elements(1),random_seed(),n_tasks_loc,rank_loc,ifail_loc)
  do pp=1,n_elements(5)
    do kk=1,n_elements(4)
      do jj=1,n_elements(3)
        do ii=1,n_elements(2)
          call rng%next(array_sol(:,ii,jj,kk,pp))
        enddo
      enddo
    enddo
  enddo
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
  file_exists=.false.; inquire(file=trim(filename),exist=file_exists)
  call assert_true(file_exists,"Error test create HDF5 file posix: file "//&
  trim(filename)//" not created!")
  call remove_file(filename,file_exists_in=file_exists)
  filename = ''; filename = trim(filename_base)//trim(extension)
  call HDF5_create(filename,file_id,access_type_in=access_hdf5_parallel,&
  mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_close(file_id)
  file_exists=.false.; inquire(file=trim(filename),exist=file_exists);
  call assert_true(file_exists,"Error test create HDF5 file access FILE_ACCESS: file "//&
  trim(filename)//" not created!")
  call remove_file(filename,file_exists_in=file_exists)
end subroutine test_create_hdf5_file

!> Test procedure for exclusively opening HDF5 file
subroutine test_open_hdf5_file()
  implicit none
  integer(HID_T)     :: file_id
  character(len=100) :: filename
  write(filename,'(A,A,I'//trim(rank_format)//',A)') & 
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_create(filename,file_id,ifail_loc); call HDF5_close(file_id); 
  call HDF5_open(filename,file_id,ifail_loc); call HDF5_close(file_id);
  call assert_equals(ifail_loc,0,"Error test open HDF5 file posix: file "//&
  trim(filename)//" not opened!"); call remove_file(filename);
  filename = ''; filename = trim(filename_base)//trim(extension); ifail_loc=0;
  call HDF5_create(filename,file_id,access_type_in=access_hdf5_parallel,&
  mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc); call HDF5_close(file_id); 
  call HDF5_create(filename,file_id,ierr=ifail_loc,access_type_in=access_hdf5_parallel,&
  mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc); call HDF5_close(file_id);
  call assert_equals(ifail_loc,0,"Error test open HDF5 file access FILE_ACCESS: file "//&
  trim(filename)//" not opened!"); call remove_file(filename);
end subroutine test_open_hdf5_file

!> test combined procedure for creating and opening HDF5 files
subroutine test_create_open_hdf5_file()
  implicit none
  integer(HID_T)     :: file_id
  character(len=100) :: filename
  logical            :: file_exists
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_close(file_id); file_exists=.false.; 
  inquire(file=filename,exist=file_exists);
  call assert_true(file_exists,"Error test create-open HDF5 file posix: file "//&
  trim(filename)//" not created!"); ifail_loc=0;
  call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_close(file_id);
  call assert_equals(ifail_loc,0,"Error test create-open HDF5 file posix: file "//&
  trim(filename)//" not opened!"); call remove_file(trim(filename)); 
  filename=''; filename = trim(filename_base)//trim(extension); ifail_loc=0;
  call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_close(file_id); file_exists=.false.; 
  inquire(file=filename,exist=file_exists);
  call assert_true(file_exists,"Error test create-open HDF5 file access FILE_ACCESS: file "//&
  trim(filename)//" not created!");
  ifail_loc=0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_close(file_id);
  call assert_equals(ifail_loc,0,"Error test create-open HDF5 file access FILE_ACCESS: file "//&
  trim(filename)//" not opened!"); call remove_file(filename);
end subroutine test_create_open_hdf5_file

!> the the posix and collective writing / reading HDF5 file
!> of integer 1D array
subroutine test_HDF5_array1D_saving_int()
  implicit none
  character(len=11),parameter      :: datasetname='array1D_int'
  integer,dimension(n_elements(1)) :: test_array,result_array
  integer(HID_T)                   :: file_id
  integer(HID_T),dimension(1)      :: offset
  character(len=100)               :: filename 
  !> initialise posix test
  test_array = 0; result_array = int(1d3*array_sol(:,1,1,1,1));
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array1D_saving_int(file_id,result_array,n_elements(1),datasetname)
  call HDF5_array1D_reading_int(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),&
  "Error test HDF5 I/O 1D integer posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[rank_loc*n_elements(1)];
  test_array = 0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array1D_saving_int(file_id,result_array,n_tasks_loc*n_elements(1),&
  datasetname,start=offset,type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array1D_reading_int(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),&
  "Error test HDF5 I/O 1D integer MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array1D_saving_int

!> the the posix and collective writing / reading HDF5 file
!> of integer 2D array
subroutine test_HDF5_array2D_saving_int()
  implicit none
  character(len=11),parameter                    :: datasetname='array2D_int'
  integer,dimension(n_elements(1),n_elements(2)) :: test_array,result_array
  integer(HID_T)                                 :: file_id
  integer(HID_T),dimension(2)                    :: offset
  character(len=100)                             :: filename 
  !> initialise posix test
  test_array = 0; result_array = int(1d3*array_sol(:,:,1,1,1));
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array2D_saving_int(file_id,result_array,n_elements(1),n_elements(2),&
  datasetname); call HDF5_array2D_reading_int(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),n_elements(2),&
  "Error test HDF5 I/O 2D integer posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[0,rank_loc*n_elements(2)];
  test_array = 0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array2D_saving_int(file_id,result_array,n_elements(1),n_tasks_loc*n_elements(2),&
  datasetname,start=offset,type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array2D_reading_int(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),n_elements(2),&
  "Error test HDF5 I/O 2D integer MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array2D_saving_int

!> the the posix and collective writing / reading HDF5 file
!> of floats 1D array
subroutine test_HDF5_array1D_saving_r4()
  implicit none
  character(len=11),parameter      :: datasetname='array1D_r4'
  real*4,dimension(n_elements(1))  :: test_array,result_array
  integer(HID_T)                   :: file_id
  integer(HID_T),dimension(1)      :: offset
  character(len=100)               :: filename 
  !> initialise posix test
  test_array = 0.0; result_array = real(array_sol(:,1,1,1,1),kind=4);
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array1D_saving_r4(file_id,result_array,n_elements(1),datasetname)
  call HDF5_array1D_reading_r4(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),&
  "Error test HDF5 I/O 1D float posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[rank_loc*n_elements(1)];
  test_array = 0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array1D_saving_r4(file_id,result_array,n_tasks_loc*n_elements(1),&
  datasetname,start=offset,type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array1D_reading_r4(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),&
  "Error test HDF5 I/O 1D float MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array1D_saving_r4

!> the the posix and collective writing / reading HDF5 file of a single double
subroutine test_HDF5_real_saving()
  implicit none
  character(len=12),parameter :: datasetname='double_value'
  real*8                      :: test_value,result_value
  integer(HID_T)              :: file_id
  character(len=100)          :: filename 
  !> initialise posix test
  test_value = 0d0; result_value = array_sol(1,1,1,1,1);
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_real_saving(file_id,result_value,datasetname)
  call HDF5_real_reading(file_id,test_value,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_value,result_value,&
  "Error test HDF5 I/O 0D double posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension);
  test_value = 0d0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_real_saving(file_id,result_value,datasetname,mpi_rank=rank_loc,&
  n_mpi_tasks=n_tasks_loc,type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_real_reading(file_id,test_value,datasetname,&
  mpi_rank=rank_loc,n_mpi_tasks=n_tasks_loc)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_value,result_value,&
  "Error test HDF5 I/O 1D double MPI collective: test and result array mismatch!")
end subroutine test_HDF5_real_saving

!> the the posix and collective writing / reading HDF5 file
!> of double 1D array
subroutine test_HDF5_array1D_saving_r8()
  implicit none
  character(len=11),parameter      :: datasetname='array1D_r8'
  real*8,dimension(n_elements(1))  :: test_array,result_array
  integer(HID_T)                   :: file_id
  integer(HID_T),dimension(1)      :: offset
  character(len=100)               :: filename 
  !> initialise posix test
  test_array = 0d0; result_array = array_sol(:,1,1,1,1);
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array1D_saving(file_id,result_array,n_elements(1),datasetname)
  call HDF5_array1D_reading(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),&
  "Error test HDF5 I/O 1D double posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[rank_loc*n_elements(1)];
  test_array = 0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array1D_saving(file_id,result_array,n_tasks_loc*n_elements(1),&
  datasetname,start=offset,type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array1D_reading(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),&
  "Error test HDF5 I/O 1D double MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array1D_saving_r8

!> the the posix and collective writing / reading HDF5 file
!> of double 2D array
subroutine test_HDF5_array2D_saving_r8()
  implicit none
  character(len=11),parameter                   :: datasetname='array2D_r8'
  real*8,dimension(n_elements(1),n_elements(2)) :: test_array,result_array
  integer(HID_T)                                :: file_id
  integer(HID_T),dimension(2)                   :: offset
  character(len=100)                            :: filename 
  !> initialise posix test
  test_array = 0d0; result_array = array_sol(:,:,1,1,1);
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array2D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  datasetname); call HDF5_array2D_reading(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),n_elements(2),&
  "Error test HDF5 I/O 2D double posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[0,rank_loc*n_elements(2)];
  test_array = 0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array2D_saving(file_id,result_array,n_elements(1),n_tasks_loc*n_elements(2),&
  datasetname,start=offset,type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array2D_reading(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals(test_array,result_array,n_elements(1),n_elements(2),&
  "Error test HDF5 I/O 2D double MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array2D_saving_r8

!> the the posix and collective writing / reading HDF5 file
!> of double 3D array
subroutine test_HDF5_array3D_saving_r8()
  use mod_assert_equals_tools, only: assert_equals_extended
  implicit none
  character(len=11),parameter                                 :: datasetname='array3D_r8'
  real*8,dimension(n_elements(1),n_elements(2),n_elements(3)) :: test_array,result_array
  integer(HID_T)                                              :: file_id
  integer(HID_T),dimension(3)                                 :: offset
  character(len=100)                                          :: filename 
  !> initialise posix test
  test_array = 0d0; result_array = array_sol(:,:,:,1,1);
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array3D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  n_elements(3),datasetname); call HDF5_array3D_reading(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals_extended(n_elements(1),n_elements(2),n_elements(3),test_array,&
  result_array,tol_r8,"Error test HDF5 I/O 3D double posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[0,0,rank_loc*n_elements(3)];
  test_array = 0d0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array3D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  n_tasks_loc*n_elements(3),datasetname,start=offset,&
  type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array3D_reading(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals_extended(n_elements(1),n_elements(2),n_elements(3),test_array,&
  result_array,tol_r8,"Error test HDF5 I/O 3D double MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array3D_saving_r8

!> the the posix and collective writing / reading HDF5 file
!> of double 4D array
subroutine test_HDF5_array4D_saving_r8()
  use mod_assert_equals_tools, only: assert_equals_extended
  implicit none
  character(len=11),parameter                     :: datasetname='array4D_r8'
  real*8,dimension(n_elements(1),n_elements(2),&
  n_elements(3),n_elements(4))                    :: test_array,result_array
  integer(HID_T)                                  :: file_id
  integer(HID_T),dimension(4)                     :: offset
  character(len=100)                              :: filename 
  !> initialise posix test
  test_array = 0d0; result_array = array_sol(:,:,:,:,1);
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array4D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  n_elements(3),n_elements(4),datasetname); 
  call HDF5_array4D_reading(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals_extended(n_elements(1),n_elements(2),n_elements(3),n_elements(4),test_array,&
  result_array,tol_r8,"Error test HDF5 I/O 4D double posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[0,0,0,rank_loc*n_elements(4)];
  test_array = 0d0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array4D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  n_elements(3),n_tasks_loc*n_elements(4),datasetname,start=offset,&
  type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array4D_reading(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals_extended(n_elements(1),n_elements(2),n_elements(3),n_elements(4),test_array,&
  result_array,tol_r8,"Error test HDF5 I/O 4D double MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array4D_saving_r8

!> the the posix and collective writing / reading HDF5 file
!> of double 5D array
subroutine test_HDF5_array5D_saving_r8()
  use mod_assert_equals_tools, only: assert_equals_extended
  implicit none
  character(len=11),parameter                     :: datasetname='array5D_r8'
  real*8,dimension(n_elements(1),n_elements(2),&
  n_elements(3),n_elements(4),n_elements(5))      :: test_array,result_array
  integer(HID_T)                                  :: file_id
  integer(HID_T),dimension(5)                     :: offset
  character(len=100)                              :: filename 
  !> initialise posix test
  test_array = 0d0; result_array = array_sol;
  write(filename,'(A,A,I'//trim(rank_format)//',A)') &
  trim(filename_base),'_rank',rank_loc,trim(extension)
  ifail_loc=0; call HDF5_open_or_create(trim(filename),file_id,ierr=ifail_loc);
  call HDF5_array5D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  n_elements(3),n_elements(4),n_elements(5),datasetname); 
  call HDF5_array5D_reading(file_id,test_array,datasetname)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals_extended(n_elements(1),n_elements(2),n_elements(3),&
  n_elements(4),n_elements(5),test_array,result_array,tol_r8,&
  "Error test HDF5 I/O 5D double posix: test and result array mismatch!")
  filename = trim(filename_base)//trim(extension); offset=[0,0,0,0,rank_loc*n_elements(5)];
  test_array = 0d0; call HDF5_open_or_create(filename,file_id,ierr=ifail_loc,&
  access_type_in=access_hdf5_parallel,mpi_comm=mpi_comm_loc,mpi_info=mpi_info_loc)
  call HDF5_array5D_saving(file_id,result_array,n_elements(1),n_elements(2),&
  n_elements(3),n_elements(4),n_tasks_loc*n_elements(5),datasetname,start=offset,&
  type_dataset_transfert_in=type_dataset_transfert_mpi)
  call HDF5_array5D_reading(file_id,test_array,datasetname,start=offset)
  call HDF5_close(file_id); call remove_file(filename);
  call assert_equals_extended(n_elements(1),n_elements(2),n_elements(3),&
  n_elements(4),n_elements(5),test_array,result_array,tol_r8,&
  "Error test HDF5 I/O 5D double MPI collective: test and result array mismatch!")
end subroutine test_HDF5_array5D_saving_r8

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
