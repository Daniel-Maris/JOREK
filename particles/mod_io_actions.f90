!> Particle input-output event module
module mod_io_actions
use mod_particle_io
use mod_event
use mod_particle_sim
implicit none
private
public io_action, read_action, write_action, &
    get_filename

type, extends(action), abstract :: io_action
  character(len=120), public :: filename = ''         !< Filename to use
  character(len=80), public  :: basename = 'part'     !< If no filename, use basename + digits
  integer, public            :: decimal_digits = 3    !< Number of decimals before the point in timestamp
  integer, public            :: fractional_digits = 8 !< Number of decimals after the point
  character(len=5), public   :: extension = '.h5'     !< I/O file type extension
  integer, public            :: access_type = 1       !< type of access to the hdf5 file 1: mpi collective
  integer, public            :: mpi_comm_io           !< mpi communicator
  integer, public            :: mpi_info_io           !< mpi information 
  logical, public            :: original = .false.    !< if true, use original I/O procedures (needed only for legacy)
  contains
    procedure :: get_filename
end type io_action

type, extends(io_action) :: read_action
  real*8  :: time             !< used with the formats from io_action if filename is unset
  logical :: legacy = .false. !< if true read old io files
  logical :: test   = .false. !< if true avoid to read adas for unit testing
contains
  procedure :: do => do_read_action
end type read_action
interface read_action
  module procedure new_read_action
end interface read_action

type, extends(io_action) :: write_action
  integer :: file_access
  integer :: type_dataset_transfert = 1 !< mpi collective dataset transfert
contains
  procedure :: do => do_write_action
end type write_action
interface write_action
  module procedure new_write_action
end interface write_action

contains
!> Write a filename consisting of this%basename and time as floating-point
!> number, without spaces (tricky), with this%decimal_digits before the `.` and
!> this%fractional_digits behind the `.`. Set these two to 0 to just set the name.
function get_filename(this, time) result(filename)
  class(io_action), intent(in) :: this
  real*8, intent(in)           :: time
  character(len=120)           :: filename
  character(len=20)            :: format
  if (this%decimal_digits .eq. 0 .and. this%fractional_digits .eq. 0) then
    write(filename,'(A,A)') trim(this%basename), this%extension
  else if (this%fractional_digits .eq. 0) then
    write(format,'(A,I0,A)') '(A,I', this%decimal_digits, ',A)'
    write(filename,trim(format)) trim(this%basename), int(time), this%extension
  else
    write(format,'(A,I0,A,I0,A,I0,A)') '(A,I', this%decimal_digits, '.', this%decimal_digits, &
        ',f0.', this%fractional_digits, ',A)'
    write(filename,trim(format)) trim(this%basename), int(time), time-real(int(time)), this%extension
  end if
end function get_filename

!> Constructor for read_action.
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_read_action(filename, basename, decimal_digits, fractional_digits, extension, &
access_type_in, mpi_comm_in, mpi_info_in, legacy_in,original_in, test_in)
  use mpi
  implicit none
  type(read_action) :: new_read_action
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  character(len=*), intent(in), optional :: extension
  integer,          intent(in), optional :: access_type_in,mpi_comm_in,mpi_info_in
  logical,          intent(in), optional :: legacy_in,original_in,test_in
  new_read_action%mpi_comm_io = MPI_COMM_WORLD
  new_read_action%mpi_info_io = MPI_INFO_NULL
  if (present(filename)) new_read_action%filename = filename
  if (present(basename)) new_read_action%basename = basename
  if (present(decimal_digits)) new_read_action%decimal_digits = decimal_digits
  if (present(fractional_digits)) new_read_action%fractional_digits = fractional_digits
  if (present(extension)) new_read_action%extension = extension
  if (present(access_type_in)) new_read_action%access_type = access_type_in
  if (present(mpi_comm_in)) new_read_action%mpi_comm_io = mpi_comm_in
  if (present(mpi_info_in)) new_read_action%mpi_info_io = mpi_info_in
  if (present(legacy_in)) new_read_action%legacy = legacy_in
  if (present(original_in)) new_read_action%original = original_in
  if (present(test_in)) new_read_action%test = test_in
  new_read_action%name = "ReadAction"
  new_read_action%log = .true.
end function new_read_action

!> Action for reading the simulation
subroutine do_read_action(this, sim, ev)
  class(read_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  if (len_trim(this%filename) .eq. 0) then
    if (this%original) then
      call read_simulation_hdf5_original(sim, trim(this%get_filename(this%time)),&
      test_in=this%test)
    else
      call read_simulation_hdf5(sim, trim(this%get_filename(this%time)), &
      access_type_in=this%access_type, mpi_comm_in=this%mpi_comm_io, &
      mpi_info_in=this%mpi_info_io, legacy_in=this%legacy, test_in=this%test)
    end if
  else
    if (this%original) then
      call read_simulation_hdf5_original(sim, trim(this%filename), &
      test_in=this%test)
    else
      call read_simulation_hdf5(sim, trim(this%filename), &
      access_type_in=this%access_type, mpi_comm_in=this%mpi_comm_io, &
      mpi_info_in=this%mpi_info_io, legacy_in=this%legacy, test_in=this%test)
    end if
  end if
end subroutine do_read_action


!> Constructor for write_action
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_write_action(filename, basename, decimal_digits, fractional_digits, extension, &
file_access_in, access_type_in, mpi_comm_in, mpi_info_in, type_dataset_transfert_in,original_in)
  use mpi
  use hdf5, only: H5F_ACC_TRUNC_F
  implicit none
  type(write_action) :: new_write_action
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  character(len=*), intent(in), optional :: extension
  integer,          intent(in), optional :: file_access_in,mpi_comm_in,mpi_info_in
  integer,          intent(in), optional :: access_type_in,type_dataset_transfert_in
  logical,          intent(in), optional :: original_in
  new_write_action%name = "WriteAction"
  new_write_action%log = .true.
  new_write_action%file_access = H5F_ACC_TRUNC_F
  new_write_action%mpi_comm_io = MPI_COMM_WORLD
  new_write_action%mpi_info_io = MPI_INFO_NULL
  if (present(filename)) new_write_action%filename = filename
  if (present(basename)) new_write_action%basename = basename
  if (present(decimal_digits)) new_write_action%decimal_digits = decimal_digits
  if (present(fractional_digits)) new_write_action%fractional_digits = fractional_digits
  if (present(extension)) new_write_action%extension = extension
  if (present(file_access_in)) new_write_action%file_access = file_access_in
  if (present(access_type_in)) new_write_action%access_type = access_type_in
  if (present(mpi_comm_in)) new_write_action%mpi_comm_io = mpi_comm_in
  if (present(mpi_info_in)) new_write_action%mpi_info_io = mpi_info_in
  if (present(type_dataset_transfert_in)) new_write_action%type_dataset_transfert = type_dataset_transfert_in
  if (present(original_in)) new_write_action%original = original_in
end function new_write_action

!> Action for writing the simulation
subroutine do_write_action(this, sim, ev)
  class(write_action), intent(inout) :: this
  type(particle_sim), intent(inout)  :: sim
  type(event), intent(inout), optional :: ev
  if (len_trim(this%filename) .eq. 0) then
    if (this%original) then
      call write_simulation_hdf5_original(sim, trim(this%get_filename(sim%time))) 
    else
      call write_simulation_hdf5(sim, trim(this%get_filename(sim%time)), &
      file_access_in=this%file_access, access_type_in=this%access_type, &
      type_dataset_transfert_in=this%type_dataset_transfert, &
      mpi_comm_in=this%mpi_comm_io, mpi_info_in=this%mpi_info_io)
    end if
  else
    if (this%original) then
      call write_simulation_hdf5_original(sim, trim(this%filename)) 
    else
      call write_simulation_hdf5(sim, trim(this%filename), &
      file_access_in=this%file_access, access_type_in=this%access_type, &
      type_dataset_transfert_in=this%type_dataset_transfert, &
      mpi_comm_in=this%mpi_comm_io, mpi_info_in=this%mpi_info_io)
    end if
  end if
end subroutine do_write_action
end module mod_io_actions
