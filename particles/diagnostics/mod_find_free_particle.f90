module mod_find_free_particle
    use mod_particle_types
    !$ use omp_lib
    implicit none
    private
    public find_free_particles

contains


function find_free_particles(particles) result(i_free)

    class (particle_base), intent(in), dimension(:)    :: particles
    integer, dimension(:), allocatable ::  i_free
  
    logical, dimension(:), allocatable :: is_free
    integer                            :: j,k, n_free
  
    allocate(is_free(size(particles,1))) 
  !$omp parallel do default(none) shared(particles, n_free, i_free, is_free) &
  !$omp private(j) &
  !$omp schedule(dynamic, 100)
  do j=1,size(particles,1) !sim%groups(1)%particles
      is_free(j) = particles(j)%i_elm .le. 0  !< array T/F is particle is free
  end do
  !$omp end parallel do
  !$omp barrier
  n_free = count(is_free)
  allocate(i_free(n_free))
  k = 1
  do j=1,size(is_free,1)
      if (is_free(j)) then
        i_free(k) = j !< i_free(k) has index of free particle in  sim%groups(atoms)%particles(j)
        k = k+1
        !if (sim%my_id .eq. 0) write(*,*) "Adding to the list number: ", j
      end if
  end do
  
  end function find_free_particles

end module mod_find_free_particle