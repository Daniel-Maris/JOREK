!> module to contain common tools for initialisation and destruction of particles
module mod_particle_init
  
  implicit none
   
  private
  public :: free_particle_indices

contains

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
