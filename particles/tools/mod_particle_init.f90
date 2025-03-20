!> module to contain common tools for initialisation and destruction of particles
module mod_particle_init
  
  implicit none
   
  private
  public :: free_particle_indices, determine_supers_create_scheme

contains


!> checks the set parameters for the create scheme to determine the scheme, or stop for wrong input 
subroutine determine_supers_create_scheme(supers_num, supers_weight, supers_ratio, scheme, default, my_id, identifier)
  implicit none

  integer,           intent(in)           :: supers_num
  real*8,            intent(in)           :: supers_weight 
  real*8,            intent(inout)        :: supers_ratio  !< inout because it will take the default value if nothing is specified
  character(len=10), intent(out)          :: scheme        !< chosen scheme (num, weight or ratio)
  real*8,            intent(in), optional :: default       !< default for supers_ratio for this scheme
  integer,           intent(in), optional :: my_id         !< for printing errors
  character(len=*),  intent(in), optional :: identifier    !< for tracing back the origin of the error
  

  integer :: id, n_create_schemes
  real*8  :: supers_ratio_default
  character(len=1000) :: identifier_local

  id = 0
  if(present(my_id)) id = my_id
  
  identifier_local = ""
  if(present(identifier)) identifier_local = identifier

  supers_ratio_default = 1.d-4
  if(present(default)) supers_ratio_default = default
  
  !> determine supers_create_scheme (which scheme is used to calculate supers_to_create) -----
  n_create_schemes = 0

  ! check that the supers_... settings are non-negative
  if (supers_num /= -1.d0) then 
    if (supers_num < 0) then
      if (id == 0) write(*,"(A,A,A)") "ERROR: ",trim(identifier_local),"%supers_num is negative"
      stop
    endif 
    scheme = "num"
    n_create_schemes = n_create_schemes + 1
  endif
  if (supers_weight /= -1.d0) then
    if (supers_weight < 0) then
      if (id == 0) write(*,"(A,A,A)") "ERROR: ",trim(identifier_local),"%supers_weight is negative"
      stop
    endif 
    scheme = "weight"
    n_create_schemes = n_create_schemes + 1
  endif
  if (supers_ratio /= -1.d0) then
    if (supers_ratio < 0) then
      if (id == 0) write(*,"(A,A,A)") "ERROR: ",trim(identifier_local),"%supers_ratio is negative"
      stop
    endif
    scheme = "ratio"
    n_create_schemes = n_create_schemes + 1
  endif

  ! check that only 1 types of supers_... settings are set, or else use default setting
  if (n_create_schemes > 1) then
    if (id == 0) then 
      write(*,"(A,A)") "ERROR: in create scheme ",trim(identifier_local)
      write(*,*) "  Only one type of supers_... can be used per creation event"
    endif
    stop
  else if (n_create_schemes == 0) then
    scheme = "ratio"
    supers_ratio = supers_ratio_default
    if (id == 0) then
      write(*,"(A,I2,A,I2,A)") "WARNING: in create scheme ",trim(identifier_local)
      write(*,*) "  no scheme for determining the number of superparticles (supers_...)"
      write(*,*) "  has been assigned. Using the default"
      write(*,*) "  setting of supers_ratio = ", supers_ratio_default
    endif
  endif

end subroutine determine_supers_create_scheme


!> determines the indices of all free particles for a particles array (typically sim%particle_group(group_num)%particles)
subroutine free_particle_indices(part_arr, n_free, i_free, n_needed)
  use mod_particle_types, only: particle_base
  use phys_module, only: use_manual_random_seed
  !$ use omp_lib 

  implicit none
  
  class(particle_base), dimension(:), allocatable, intent(in)  :: part_arr          !< typically sim%particle_group(group_num)%particles
  integer,                                         intent(out) :: n_free            !< number of free particles
  integer,             dimension(:), allocatable,  intent(out) :: i_free            !< indices of free particles (size: n_free or n_needed if specified)
  integer,                           optional,     intent(in)  :: n_needed          !< number of free particles needed. If this is specified, routine returns first n_needed free particles instead 
  
  logical, dimension(:), allocatable :: is_free
  integer, dimension(:), allocatable :: i_free_tmp
  integer :: i, j, k, k_thread, n_part
  
  n_part = size(part_arr)

  ! Step 1: find which particles are free using a mask is_free
  allocate(is_free(n_part)) 

  ! might be replaced with omp workshare, or just the array expression.
  ! there is an issue with derived type arrays in gfortran though, and this works
  if(use_manual_random_seed) then
    !$ call omp_set_schedule(omp_sched_static,100)
  else
    !$ call omp_set_schedule(omp_sched_dynamic,100)
  end if
  !$omp parallel do default(none) shared(part_arr, n_part, is_free) &
  !$omp private(j) schedule(runtime)
  do j=1,n_part
    ! i_elm < 0 technically means the particle left the domain rather than that it is free, thus to easier catch bugs, only i_elm = 0 particles are considered free by default 
    !however for reg test purpose we need .le.
    !is_free(j) = part_arr(j)%i_elm .eq. 0 
    is_free(j) = part_arr(j)%i_elm .le. 0 
  end do
  !$omp end parallel do
  
  ! what is the function of this barrier?
  !$omp barrier
  
  n_free = count(is_free)
  allocate(i_free(n_free))
  
  ! Step 2: write their indices in an array
  k = 1

  ! omp parallelisation breaks reg test because of rng
  ! if(use_manual_random_seed) then
  !   !$ call omp_set_schedule(omp_sched_static,1)
  ! else
  !   !$ call omp_set_schedule(omp_sched_dynamic,100)
  ! end if
  ! !$omp parallel do default(none) shared(is_free, i_free, n_part, k) &
  ! !$omp private(j, k_thread) schedule(runtime)
  do j=1,n_part
    if (is_free(j)) then
      ! !$omp atomic capture
      k_thread = k
      k = k+1
      ! !$omp end atomic
      i_free(k_thread) = j
    end if
  end do
  ! !$omp end parallel do

  ! Step 3: give first n_needed elements back if n_needed was specified, else return array with all free particle indices
  ! TODO: could copying be sped up using omp?
  if(present(n_needed)) then
    if(n_free < n_needed) then
      write(*,*) "ERROR: More free particles needed than available, returning only available free particles (avail/needed)", n_free, n_needed
      ! then just send the full i_free
      return
    end if
    ! make i_free smaller
    i_free_tmp = i_free(1:n_needed)
    deallocate(i_free)
    i_free = i_free_tmp
  end if

end subroutine free_particle_indices

end module mod_particle_init
