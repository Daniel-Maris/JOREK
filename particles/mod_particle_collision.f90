!> Module for neutral-neutral collisions
!> Nearest neighbour approach
module mod_particle_collision
  use mod_edge_elements
  use mod_particle_types
  use constants
  use mod_pcg32_rng,   only: pcg32_rng
  use mod_particle_sim
  use phys_module, only: tstep, n_period, use_manual_random_seed !< we need the intrinsic fortran gamma function so have to use only
  use mod_interp
  use mod_event, only: mpi_minmeanmax
  use mpi
  
  !$ use omp_lib

  implicit none
   
  private
  public :: neutral_self_collision

  !> an object containing the array of particle indices of particles in this element 
  !> s.t. an array of these objects can store all particle indices as array per element
  type :: indices_in_elm
    integer, dimension(:), allocatable :: pa_ind
  end type

contains

!> neutral-neutral collision within a species 
!> nearest neighbour approach where each element has a number of subcells in the toroidal direction
!> the number of subcells in the poloidal direction are done by dividing the element up in an s,t grid
!> which is determined automatically by the amount of particles and the length of the element along s and t 
subroutine neutral_self_collision(sim, rng, dt)
  implicit none

  type(particle_sim),            intent(inout) :: sim
  type(pcg32_rng), dimension(:), intent(inout) :: rng

  real*8,                        intent(in)    :: dt !< timestep over which the collisions must be calculated

  type(particle_kinetic_leapfrog), dimension(:),     allocatable :: pa_subbin !< array of particles in this collisional subbin
  type(indices_in_elm),            dimension(:),     allocatable :: sorted_ind_arr !< object containing all particle indices as arrays per element number
  type(indices_in_elm),            dimension(:,:,:), allocatable :: i_pa_bin !< object containing all this element's particle indices as arrays per bin (s bin,t bin, phi bin)
  
  integer :: n_phi !< number of toroidal bins to do collisions in
  integer :: n_elm
  integer :: i_pa_group, i_phi, is, it, ns, nt, nst, n_pa_bin, n_pa_subbin, aim_pa_per_subbin, i_pair
  integer :: i_subbin !< the bin index of the first particle in this subbin (so if this is the second subbin for aim_per_pa=50, i_subbin=51, etc)
  real*8, parameter :: subbin_factor = 2 !< factor by which a collisional bin is allowed to be larger than aim_pa_per_bin before splitting into subbings
  integer :: i,j, i_elm, i_loc, i_global
  
  real*8 :: P_coll  !< [chance] P = 2 w_g n \sigma_T v_r dt
  real*8 :: w2, w1  !< [physical particles] weight of super particle 1
  real*8 :: w_s     !< [physical particles] weight of smaller super particle in binary collision w_s = min(w1,w2)
  real*8 :: w_g     !< [physical particles] weight of greater super particle in binary collision w_g = max(w1, w2)
  real*8 :: v_r     !< [m/s] scalar relative velocity of collisional pair
  real*8 :: sigma_T !< [m^2] total collisional cross section of the collision pair

  integer :: pa_in_elm !< number of particles in the element under consideration
  integer, parameter :: aim_pa_per_bin = 50 !< wanted amount of super particles in one collisional bin (element is split up to satisfy this)
  integer, dimension(:),       allocatable :: i_pa_elm     !< global particle indices of the particles in this element
  !integer, dimension(:,:,:,:), allocatable :: i_pa_bin     !< global particle indices of the particles in this collisional bin (bin particle, s bin,t bin, phi bin)
  integer, dimension(:,:,:),   allocatable :: n_pa_bin_arr !< number of particles in each collisional bin of this element (s bin, t bin, phi bin)
  integer, dimension(:,:,:),   allocatable :: i_loc_bin    !< local index of particles in each collisional bin of this element (s bin, t bin, phi bin)
  logical, dimension(:),       allocatable :: paired       !< array to keep track of which particle is not paired yet
  integer, dimension(:),       allocatable :: pa_elm_arr   !< precalculated pa(:)%i_elm for faster masking over i_elm 
  integer, dimension(:),       allocatable :: pa_in_elm_arr !< number of particles in the element (n_elm)
  integer, dimension(:),       allocatable :: i_loc_arr !< local index in the element array of particle indices (n_elm)
  
  
  real :: R, R_s, R_t, Z, Z_s, Z_t
  real :: lt,ls !< [m] physical size of element / bin in s direction
  real :: V_c !< [m^3] grid cell (collisional bin) volume 
  real*8 :: RN(3) !< random numbers for chance to collide, impact parameter and scattering plane angle
  integer :: i_rng

  !nearest neighbour
  real*8  :: RZPhi(3), RZPhi_try(3), d2_try, d2, d, w2_phi !< squared weight to scale a distance in phi to a distance in R, Z
  integer :: i_pa1, i_pa2
  real*8  :: n_bin !< [physical particles/m^3] number density in the bin

  !elastic collision parameters
  real*8 :: alpha  !< scattering plane angle to Z-axis around *v1*
  real*8 :: Theta  !< scattering angle in normal coordinates
  real*8 :: m2, m1 !< mass of particle 1
  real*8 :: v1f(3) !< final velocity of particle 1
  real*8 :: v2f(3) !< final velocity of particle 2
  
  !dbg
  real*8 :: t(6)  !< cpu times (begin, end, diff)
  real*8 :: t_mask, t_rest, t_elm, t_priv(25), t_priv_tot(25)
  integer :: max_n_pa(2)
  real*8 :: P_max, P_max_global
  integer :: ierr
  !$ real*8 :: w(2), mmm(3)

  ! --- start of code

  !$ w(1) = omp_get_wtime()
  
  t = 0.d0
  call cpu_time(t(1))

  n_elm = sim%fields%element_list%n_elements
  
  i_rng = 1 !default if not using OMP

  do i_pa_group = 1,size(sim%groups,1)
  
    select type (pa => sim%groups(i_pa_group)%particles)
    type is (particle_kinetic_leapfrog)
      
      !allocate here so that different groups can have different resolutions in the future, this should probably just be based on n_plane
      n_phi = 1 !4 !n_plane
      
      ! --- sorting the particle indices into array by their element number
      call cpu_time(t(4))
      
      ! allocate(pa_ind_arr(size(pa)))
      if(allocated(pa_elm_arr)) deallocate(pa_elm_arr)
      allocate(pa_elm_arr(size(pa)))
      if(.not. allocated(pa_in_elm_arr)) allocate(pa_in_elm_arr(n_elm))
      pa_elm_arr = 0
      pa_in_elm_arr = 0
      
      !$omp parallel do default(none)  &
      !$omp shared(pa,pa_elm_arr)      &
      !$omp private(i_elm)             &
      !$omp reduction(+:pa_in_elm_arr)
      do i=1,size(pa)
        i_elm = pa(i)%i_elm
        if(i_elm < 1) cycle
        pa_elm_arr(i) = i_elm
        pa_in_elm_arr(i_elm) = pa_in_elm_arr(i_elm) + 1
      end do
      !$omp end parallel do

      !allocating sorted_ind_arr object
      if(.not. allocated(sorted_ind_arr)) allocate(sorted_ind_arr(n_elm))
      !$omp parallel do default(none) &
      !$omp shared(pa_in_elm_arr, sorted_ind_arr)
      do i=1,size(sorted_ind_arr)
        if(allocated(sorted_ind_arr(i)%pa_ind)) deallocate(sorted_ind_arr(i)%pa_ind) ! this can be done here as we know size(sorted_ind_arr)=n_elm is fixed
        allocate(sorted_ind_arr(i)%pa_ind(pa_in_elm_arr(i)))
      end do
      !$omp end parallel do

      !filling sorted_ind_arr object
      allocate(i_loc_arr(n_elm))
      i_loc_arr = 1
      !$omp parallel do default(none)                     &
      !$omp shared(pa_elm_arr, i_loc_arr, sorted_ind_arr) &
      !$omp private(i_elm, i_loc)
      do i=1,size(pa_elm_arr)
        i_elm = pa_elm_arr(i)
        if(i_elm < 1) cycle
        !$omp atomic capture
        i_loc = i_loc_arr(i_elm)
        i_loc_arr(i_elm) = i_loc_arr(i_elm) + 1
        !$omp end atomic
        sorted_ind_arr(i_elm)%pa_ind(i_loc) = i
      end do
      !$omp end parallel do
      deallocate(i_loc_arr)

      call cpu_time(t(5))
      t(6) = t(5) - t(1) !t(6) + t(5) - t(4)
      
      t_priv = 0.d0
      t_priv_tot = 0.d0
      t_mask = 0.d0
      t_rest = 0.d0
      t_elm  = 0.d0
      max_n_pa = 0
      P_max = 0
      ! Start loop over each finite element number
      if(use_manual_random_seed) then
        !$ call omp_set_schedule(omp_sched_static,1)
      else
        !$ call omp_set_schedule(omp_sched_dynamic,1)
      end if
#ifdef __GFORTRAN__
      !$omp parallel do default(shared)                                                            &
#else
      !$omp parallel do default(none)                                                              &
      !$omp shared(sim, rng, dt,                                                                   &
      !$omp n_elm, n_phi, i_pa_group, sorted_ind_arr, pa_in_elm_arr)                               &
#endif
      !$omp schedule(runtime)                                                                      &
      !$omp private(i_pa_elm, pa_in_elm, i, i_global, i_phi, it, is, nst, ns, nt, R, R_s, R_t, Z, Z_s, Z_t,  &
      !$omp ls, lt, i_pa_bin, n_pa_bin_arr, n_pa_bin, i_subbin, n_pa_subbin, pa_subbin, n_bin, w2_phi, V_c,  &
      !$omp paired, i_pair, i_pa1, i_pa2, RZPhi, d2, RZPhi_try, d2_try, w1, w2, w_s, w_g, v_r, m1, m2,  &
      !$omp sigma_T, P_coll, i_rng, RN, Theta, alpha, v1f, v2f, t_priv, i_loc, i_loc_bin, aim_pa_per_subbin) &
      !$omp reduction(+:t_mask, t_rest, t_elm, t_priv_tot) reduction(max:max_n_pa, P_max)
      do i_elm=1,n_elm
        t_priv=0
        call cpu_time(t_priv(1))
        !$ i_rng = omp_get_thread_num()+1
        t_priv_tot=0

        max_n_pa = 0
        P_max=0
        
        if(allocated(i_pa_elm)) deallocate(i_pa_elm)
        allocate(i_pa_elm(pa_in_elm_arr(i_elm)))
        i_pa_elm = sorted_ind_arr(i_elm)%pa_ind(:)
        
        pa_in_elm = size(i_pa_elm,1)        
        
        max_n_pa(1) = pa_in_elm

        call cpu_time(t_priv(2))
        t_mask = t_priv(2) - t_priv(1)
        t_elm  = t_mask ! overwrite later if rest of loop is done
        t_rest = 0.d0

        !> skip collisions if there are nearly no particles
        if(pa_in_elm .lt. 2*n_phi) cycle
        
        call cpu_time(t_priv(4))
        t_priv_tot(4) = t_priv(4)-t_priv(2)

        !determine how many poloidal bins for this element. floor to make sure the average number of particles does not get much smaller than aim_pa_per_bin
        nst = floor((pa_in_elm/real(aim_pa_per_bin * n_phi))) !< number of bins in poloidal plane

        if(nst .gt. 1) then        
          !determine how to distribute the poloidal bins based on what size the element is along the s and t coordinates
          call interp_RZ(sim%fields%node_list,sim%fields%element_list,i_elm,0.5d0,0.5d0,R,R_s,R_t,Z,Z_s,Z_t)
          ls = sqrt(R_s**2 + Z_s**2)
          lt = sqrt(R_t**2 + Z_t**2)
          
          !aim for ns*nt = nst together with (ls/lt)nt = ns and (ls/lt)ns = nt
          ns = nint(sqrt(nst*ls/lt))
          nt = nint(sqrt(nst*lt/ls))

          !if the element is very elongated, make sure there's always at least one bin in both directions
          if(ns .le. 1) then
            ns = 1
            nt = nst
          end if
          if(nt .le. 1) then
            nt = 1
            ns = nst
          end if
        else 
          ns = 1
          nt = 1
        end if

        ! !$omp critical
        !   write(*,"(A,7I14)") "DBG ",i_elm,ns,nt,n_phi,ns*nt*n_phi,ns*nt*n_phi*aim_pa_per_bin,pa_in_elm
        ! !$omp end critical

        call cpu_time(t_priv(18))
        t_priv_tot(18) = t_priv(18)-t_priv(4)

        allocate(i_pa_bin(ns,nt,n_phi))
        allocate(n_pa_bin_arr(ns,nt,n_phi))
        
        call cpu_time(t_priv(20))
        t_priv_tot(20) = t_priv(20)-t_priv(18)

        ! i_pa_bin = 0
        n_pa_bin_arr = 0

        call cpu_time(t_priv(5))
        t_priv_tot(5) = t_priv(5)-t_priv(4)

        t_priv_tot(19) = t_priv(5)-t_priv(18)
        t_priv_tot(1)  = t_priv(5)-t_priv(20)

        !loop to get amount of particles per bin
        do i=1,pa_in_elm
          i_global = i_pa_elm(i)
          is = ith_bin(pa(i_global)%st(1),ns)
          it = ith_bin(pa(i_global)%st(2),nt)
          
          i_phi = ceiling(pa(i_global)%x(3)/ (2.d0 * PI / float(n_period*n_phi)))
          if (i_phi .gt. n_phi) i_phi = mod(i_phi,n_phi)
          if (i_phi .lt. 1)     i_phi = mod(i_phi,n_phi) + n_phi
          
          n_pa_bin_arr(is,it,i_phi) = n_pa_bin_arr(is,it,i_phi) + 1
        end do !i_pa_elm

        call cpu_time(t_priv(21))
        t_priv_tot(21) = t_priv(21)-t_priv(5)

        !loop to allocate bin sizes in i_pa_bin
        do i_phi=1,n_phi
          do it=1,nt
            do is=1,ns
              allocate(i_pa_bin(is,it,i_phi)%pa_ind(n_pa_bin_arr(is,it,i_phi)))
            end do
          end do
        end do

        call cpu_time(t_priv(22))
        t_priv_tot(22) = t_priv(22)-t_priv(21)

        !fill i_pa_bin
        allocate(i_loc_bin(ns,nt,n_phi))
        i_loc_bin = 1
        do i=1,pa_in_elm
          i_global = i_pa_elm(i)
          is = ith_bin(pa(i_global)%st(1),ns)
          it = ith_bin(pa(i_global)%st(2),nt)
          
          i_phi = ceiling(pa(i_global)%x(3)/ (2.d0 * PI / float(n_period*n_phi)))
          if (i_phi .gt. n_phi) i_phi = mod(i_phi,n_phi)
          if (i_phi .lt. 1)     i_phi = mod(i_phi,n_phi) + n_phi
          
          i_loc = i_loc_bin(is,it,i_phi)
          i_loc_bin(is,it,i_phi) = i_loc_bin(is,it,i_phi) + 1
          i_pa_bin(is,it,i_phi)%pa_ind(i_loc) = i_global
        end do !i_pa_elm
        
        call cpu_time(t_priv(6))
        t_priv_tot(6) = t_priv(6)-t_priv(5)

        t_priv_tot(2) = t_priv(6)-t_priv(22)

        !loop through each collisional bin
        do i_phi=1,n_phi 
          do it=1,nt
            do is=1,ns 
              call cpu_time(t_priv(7))
              
              n_pa_bin = n_pa_bin_arr(is,it,i_phi) ! shorthand notation
              
              if(n_pa_bin > max_n_pa(2)) max_n_pa(2) = n_pa_bin
              if(n_pa_bin .le. 1) cycle !< can't collide 0 or 1 particles
              
              ! detemine quantities relevant to all subbins: n_bin and w2_phi
              call interp_RZ(sim%fields%node_list,sim%fields%element_list,i_elm,(real(is)+0.5d0)/real(ns),(real(it)+0.5d0)/real(nt),R,R_s,R_t,Z,Z_s,Z_t)
              ls = sqrt(R_s**2 + Z_s**2)/real(ns)
              lt = sqrt(R_t**2 + Z_t**2)/real(nt)
              V_c = ls*lt*(2.d0*PI/real(n_period*n_phi))
              !n_bin = sum(pa(i_pa_bin(is,it,i_phi)%pa_ind(:))%weight)/V_c !n = N/V_c, V_c cylindrical approximation temporary until we can read from projections?
              
              !determine weight of phi direction necessary for nearest neighbour calculation based on bin size
              if(n_period .eq. 1 .and. n_phi .eq. 1) then !< axi symmetry
                w2_phi = 0.d0
              else
                w2_phi = ((ls+lt)/(2.d0*PI/real(n_period*n_phi)))**2
              end if

              ! If many particles are in exactly the same location, they will all be in the same bin, meaning the nearest neighbour approach is very slow, 
              ! thus subbins are used to spread these highly localised particles into several sub bins for computation sake
              if(n_pa_bin > ceiling(subbin_factor*aim_pa_per_bin)) then !much more particles than expected in this bin, so use subbinning
                aim_pa_per_subbin = aim_pa_per_bin
              else !only one subbin, effectively skipping subbin implementation
                aim_pa_per_subbin = n_pa_bin
              end if

              call cpu_time(t_priv(24))
              t_priv_tot(24) = t_priv_tot(24) + t_priv(24) - t_priv(7)

              do i_subbin=1,n_pa_bin,aim_pa_per_subbin
                call cpu_time(t_priv(23))

                !last subbin of size <= aim_pa_per_subbin
                if(i_subbin + aim_pa_per_subbin < n_pa_bin) then
                  n_pa_subbin = aim_pa_per_subbin
                else
                  n_pa_subbin = n_pa_bin - i_subbin + 1
                end if

                if(n_pa_subbin .le. 1) cycle !< can't collide 1 particle

                call cpu_time(t_priv(8))
                t_priv_tot(8) = t_priv_tot(8) + t_priv(8)-t_priv(23)

                allocate(paired(n_pa_subbin))
                paired(:) = .false.

                !make a local copy of the particles in this subbin
                allocate(pa_subbin(n_pa_subbin))
                do i=1,n_pa_subbin
                  call copy_particle_kinetic_leapfrog(pa(i_pa_bin(is,it,i_phi)%pa_ind(i_subbin+i-1)),pa_subbin(i)) !copy from MPI pa array to subbin arr
                end do
                
                call cpu_time(t_priv(9))
                t_priv_tot(9) = t_priv_tot(9) + t_priv(9)-t_priv(8)

                do i_pair = 1,n_pa_subbin/2 !< dummy variable for loop, for odd n_pa_bin, there are (n_pa_bin-1)/2 pairs possible, so integer division here gives the intended result
                  call cpu_time(t_priv(12))
                
                  ! Find first unpaired particle
                  do i=i_pair,n_pa_subbin !start looking from i_pair, as all i<i_pair is definitely already used
                    if(paired(i)) cycle
                    i_pa1 = i
                    paired(i) = .true.
                    exit
                  end do
                  
                  call cpu_time(t_priv(13))
                  t_priv_tot(13) = t_priv_tot(13) + t_priv(13)-t_priv(12)

                  !> Find nearest neighbour pa_subbin(i_pa2) to pa_subbin(i_pa1)
                  RZPhi=pa_subbin(i_pa1)%x
                  d2 = 1.d99
                  i_pa2 = -1
                  do i=1,n_pa_subbin !< finds weighted nearest neighbour by looping through all particles in bin
                    if(paired(i)) cycle !< skip particles that have been paired already
                    RZPhi_try = pa_subbin(i)%x 
                    !> find lowest squared distance rather than real distance
                    d2_try = (RZPHI(1) - RZPhi_try(1))**2 + (RZPHI(2) - RZPhi_try(2))**2 + w2_phi*(RZPHI(3) - RZPhi_try(3))**2
                    if(d2_try .lt. d2) then
                      d2=d2_try
                      i_pa2=i
                    end if
                  end do !loop for nearest neighbour
                  paired(i_pa2) = .true.
                  
                  call cpu_time(t_priv(14))
                  t_priv_tot(14) = t_priv_tot(14) + t_priv(14)-t_priv(13)

                  ! calculate P for this collision pair
                  w1 = pa_subbin(i_pa1)%weight
                  w2 = pa_subbin(i_pa2)%weight
                  w_g = max(w1, w2)
                  w_s = min(w1, w2)
                  v_r = norm2(pa_subbin(i_pa1)%v - pa_subbin(i_pa2)%v)
                  m1 = sim%groups(i_pa_group)%mass
                  m2 = sim%groups(i_pa_group)%mass
                  sigma_T = calc_sigma_T(v_r, m1, m2)

                  call cpu_time(t_priv(15))
                  t_priv_tot(15) = t_priv_tot(15) + t_priv(15)-t_priv(14)
                  
                  P_coll = w_g * n_pa_bin * sigma_T * v_r * dt / V_c
                  if (P_coll > P_max) P_max = P_coll
                  if(P_coll .gt. 1.d0) &
                    write(*,"(A,7es15.5)") "ERROR in NNC: P_coll > 1 (P,n,sigma,v_r,dt,w1,w2)",P_coll, n_bin, sigma_T, v_r, dt, w1, w2

                  !generate random number
                  call rng(i_rng)%next(RN)

                  call cpu_time(t_priv(16))
                  t_priv_tot(16) = t_priv_tot(16) + t_priv(16)-t_priv(15)

                  if (RN(1) .le. P_coll) then ! do binary collision

                    Theta = PI - 2*asin(RN(2))
                    alpha = TWOPI*RN(3)

                    ! updates velocities of colliding part of super particles
                    call binary_elastic_collision(m1,m2,pa_subbin(i_pa1)%v,pa_subbin(i_pa2)%v,Theta,alpha,v1f,v2f)

                    !correct new velocities for non-colliding component in heaviest particle
                    if(w1 .gt. w2) then
                      pa_subbin(i_pa1)%v = (w_s/w1) * v1f + (1-w_s/w1) * pa_subbin(i_pa1)%v
                      pa_subbin(i_pa2)%v = v2f
                    else
                      pa_subbin(i_pa1)%v = v1f
                      pa_subbin(i_pa2)%v = (w_s/w2) * v2f + (1-w_s/w2) * pa_subbin(i_pa2)%v
                    end if 

                  end if !collision happens

                  call cpu_time(t_priv(17))
                  t_priv_tot(17) = t_priv_tot(17) + t_priv(17)-t_priv(16)
                
                end do !loop over i_pair

                call cpu_time(t_priv(10))
                t_priv_tot(10) = t_priv_tot(10) + t_priv(10)-t_priv(9)
  
                ! copy back into MPI pa array
                do i=1,n_pa_subbin
                  call copy_particle_kinetic_leapfrog(pa_subbin(i), pa(i_pa_bin(is,it,i_phi)%pa_ind(i_subbin+i-1))) 
                end do
  
                deallocate(pa_subbin)
                
                deallocate(paired)
  
                call cpu_time(t_priv(11))
                t_priv_tot(11) = t_priv_tot(11) + t_priv(11)-t_priv(10)
                t_priv_tot(23) = t_priv_tot(23) + t_priv(11)-t_priv(23) 

              end do !loop over i_subbin
            
            call cpu_time(t_priv(25))
            t_priv_tot(25) = t_priv_tot(25) + t_priv(25) - t_priv(24)

            end do !is
          end do !it
        end do !i_phi toroidal bins
        
        if(allocated(n_pa_bin_arr)) deallocate(n_pa_bin_arr)
        deallocate(i_loc_bin)
        do i_phi=1,n_phi
          do it=1,nt
            do is=1,ns
              deallocate(i_pa_bin(is,it,i_phi)%pa_ind)
            end do
          end do
        end do
        deallocate(i_pa_bin)
        
        call cpu_time(t_priv(3))
        t_rest = t_priv(3) - t_priv(2)
        t_elm  = t_priv(3) - t_priv(1)
        t_priv_tot(3) = t_priv(3) - t_priv(6)
      end do !i_elm
      !$omp end parallel do
        
    class default
      write(*,*) "neutral-neutral self collisions not implemented for this type, group=", i
      call exit(13)
    end select

  end do !i_particle_group
  call cpu_time(t(2))
  t(3) = t(2) - t(1) ! cpu time spent by program
  if(sim%my_id .eq. 0) write(*,"(A,30es15.5)") "NNC cpu times (s): (tot, t_mask, t_rest, t_elm, sort, t_priv_tot)", t(3), t_mask, t_rest, t_elm, t(6), t_priv_tot
  call MPI_REDUCE(P_max, P_max_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  if(sim%my_id .eq. 0) write(*,"(A,2I8,es15.5)") "max (pa in elm/pa in bin/P_coll) =",max_n_pa,P_max_global
  !$ w(2) = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w(2)-w(1))
  !$ if (sim%my_id .eq. 0) write(*,"(A,3f9.4,A)") "Neutral self collision complete in (min/mean/max) ", mmm, " s"

end subroutine neutral_self_collision

!> calculates the elastic collisional cross section for D + D elastic collisions
!> variable hard sphere based on [1] eq 4.63
!> [1]: Bird, G. A.. (1994). Molecular Gas Dynamics and the Direct Simulation of Gas Flows.
function calc_sigma_T(v_r, m1, m2) result(sigma_T)
  
  use constants, only: PI, K_BOLTZ
  
  implicit none
  
  real*8, intent(in)  :: v_r     !< relative velocity
  real*8, intent(in)  :: m1      !< mass of first particle
  real*8, intent(in)  :: m2      !< mass of second particle

  real*8              :: sigma_T !< total collisional cross section sigma_T

  ! local variables
  real*8 :: d     !< diameter
  real*8 :: T_ref !< reference temperature for d_ref and omega
  real*8 :: d_ref !< diameter at reference temperature 
  real*8 :: omega !< viscosity index
  real*8 :: m_r   !< reduced mass m=m1*m2/(m1+m2)

  if (abs(v_r) .lt. 1.d-12) then
    sigma_T = 0.d0
    ! write(*,*) "NNC sigma_T(v_r=0) is set to 0"
    return
  end if

  ! setting local variables
  d_ref = 2*1.2d-10*0.8 !1.2Å vdWaals radius for H, but that should be nearly equal to D. For H2 d_ref=2.92d-10 m according to [1] appendix A, and He d_ref=2.33d-10 while r_vdWaals = 1.40d-10m for He. So factor 0.8 is to handwavingly convert from r_vdWaals to d_ref
  T_ref = 273
  omega = 0.68 !Varoutis 2017 paper 
  m_r   = m1*m2/(m1+m2)

  d = d_ref*sqrt(((2*K_BOLTZ*T_ref/(m_r*v_r**2))**(omega-0.5d0) )/ gamma(2.5d0-omega)) !< gamma here is the gamma function, not phys_module gamma
  sigma_T = PI*d**2

end function calc_sigma_T

!> Calculates the resulting velocity vectors after a binary elastic collision 
!> Uses CM scattering angle Theta and angle of the scattering plane alpha
subroutine binary_elastic_collision(m1,m2,v1i,v2i,Theta,alpha,v1f,v2f)
  use constants, only: PI

  implicit none
  
  real*8, intent(in)  :: m1      !< mass of particle 1
  real*8, intent(in)  :: m2      !< mass of particle 2
  real*8, intent(in)  :: Theta   !< scattering angle of centre-of-mass system
  real*8, intent(in)  :: alpha   !< scattering plane angle to Z-axis around v1i
  real*8, intent(in)  :: v1i(3)  !< initial velocity of particle 1
  real*8, intent(in)  :: v2i(3)  !< initial velocity of particle 2
  real*8, intent(out) :: v1f(3)  !< final velocity of particle 1
  real*8, intent(out) :: v2f(3)  !< final velocity of particle 2

  !local parameters
  ! collision velocities in scattering plane [parallel to v1, perp to v1] where v1=v_r, v2=0, and p stands for prime, the final velocity
  real*8 :: v_r(3)     !< relative velocity in R,Z,phi
  real*8 :: v_r_hat(3) !< direction of v_r
  real*8 :: scalar_v1, scalar_v1p, scalar_v2p, v1p(2), v2p(2) 
  real*8 :: theta_1, theta_2, perp1(3), perp2(3)

  !sanity
  real*8 :: dp(3), dE, tol
  real*8 :: db(20) !< debug checks, should all be 0
  
  db(:)=0.d0

  !this calculation follows the logic in Liebermann (2005), Principles of Plasma Discharges and Materials Processing, section 3.2, applied to 3D

  !go to 3D rest frame of particle 2
  v_r = v1i - v2i
  !go to scattering plane, 2D rest frame of particle 2 (i.e. v2i = 0 in this frame, while v1i = [v_r,0])
  scalar_v1 = norm2(v_r)

  !calculate scattering angles of particles in scattering plane from CM frame scattering angle Theta
  theta_2 = -(PI-Theta)/2
  theta_1 = atan(sin(Theta)/(m1/m2 + cos(Theta)))

  !calculate the velocities in scattering plane
  scalar_v1p = scalar_v1 * sqrt(1.d0/(1 + (m1/m2)*(sin(theta_1)/sin(theta_2))**2)) !< from energy balance, substituting scalar_v2p (see next line)
  scalar_v2p = -(m1/m2) * (sin(theta_1) / sin(theta_2)) * scalar_v1p               !< from perpendicular momentum balance

  v2p = scalar_v2p*[cos(theta_2), sin(theta_2)]
  v1p = scalar_v1p*[cos(theta_1), sin(theta_1)]

  !transform back to 3D rest frame of particle 2
  !we need 2 orthonormal vectors to v_r to determine the scattering plane using alpha
  perp1 = cross_product(v_r,[1.d0,0.d0,0.d0]) !< cross product with R axis to get orthogonal vector
  perp2 = cross_product(v_r,perp1)            !< cross product with first perpendicular axis
  
  v_r_hat = v_r  /norm2(v_r)   !< normalising to get orthonormal basis vectors
  perp1   = perp1/norm2(perp1)
  perp2   = perp2/norm2(perp2)

  !calculate the resulting velocity from parallel velocity component and perpendicular split up into the two perpendicular directions using alpha
  v1f = v1p(1)*v_r_hat + v1p(2)*cos(alpha)*perp1 + v1p(2)*sin(alpha)*perp2
  v2f = v2p(1)*v_r_hat + v2p(2)*cos(alpha)*perp1 + v2p(2)*sin(alpha)*perp2

  db(1)  = m1*scalar_v1**2 - m1*norm2(v1f)**2 - m2*norm2(v2f)**2
  db(2) = norm2(v1p)**2 - norm2(v1f)**2
  db(3) = norm2(v2p)**2 - norm2(v2f)**2
  db(4) = m1*v1p(1) + m2*v2p(1) - m1*scalar_v1
  db(5) = m1*v1p(2) + m2*v2p(2)
  
  !transform back from rest frame of particle 2 to R,Z,phi frame
  v1f = v1f + v2i
  v2f = v2f + v2i

  ! check conservation quantities
  dp = (m1*(v1f-v1i) + m2*(v2f-v2i))/max(1.d0,norm2(m1*v1i + m2*v2i)) !< relative change in momentum due to collision, should be 0, m is in amu
  dE = (m1*(norm2(v1f)**2-norm2(v1i)**2) + m2*(norm2(v2f)**2-norm2(v2i)**2)) / max(1.d0,m1*norm2(v1i)**2 + m2*norm2(v2i)**2) !< relative change in energy [J/(kg/amu)]

  db(6)  = m1*scalar_v1**2 - m1*scalar_v1p**2 - m2*scalar_v2p**2
  db(7)  = m1*scalar_v1**2 - m1*norm2(v1p)**2 - m2*norm2(v2p)**2
  db(8)  = norm2(v_r_hat) - 1.d0
  db(9)  = norm2(perp1)   - 1.d0
  db(10) = norm2(perp2)   - 1.d0
  db(11) = dot_product(v_r_hat,perp1)
  db(12) = dot_product(v_r_hat,perp2)
  db(13) = dot_product(perp1,perp2)
  db(14) = dp(1)
  db(15) = dp(2)
  db(16) = dp(3)
  db(17) = dE

  tol = 1.d-10
  if(abs(dp(1)) .gt. tol .or. abs(dp(1)) .gt. tol .or. abs(dp(1)) .gt. tol .or. abs(dE) .gt. tol) then
    write(*,*) "ERROR, momentum or energy not conserved in binary elastic collision (m1,m2,v1i,v2i,Theta,alpha,v1f,v2f,dp,dE)",m1,m2,v1i,v2i,Theta,alpha,v1f,v2f,dp,dE
    write(*,*) "debug: (should all be 0)"
    write(*,*) db
  end if
  
end subroutine binary_elastic_collision

!> orthonormal basis cross product
pure function cross_product(v1,v2) result(v_perp)
  implicit none
  real*8, intent(in) :: v1(3), v2(3)
  real*8 :: v_perp(3)

  v_perp(1) = v1(2)*v2(3) - v1(3)*v2(2)
  v_perp(2) = v1(3)*v2(1) - v1(1)*v2(3)
  v_perp(3) = v1(1)*v2(2) - v1(2)*v2(1)

end function cross_product

!> calculates in which bin (out of n bins) a value x=[0,1] should fall
pure function ith_bin(x_in,n) result(i)
  implicit none
  real*8,  intent(in) :: x_in !< value between 0 and 1, of which the corresponding bin should be determined
  integer, intent(in) :: n !< number of bins
  integer :: i !< ith element in the bin

  real*8 :: x !< copy of x_in
  real*8,parameter :: tol=1.d-10 !< numerical tolerance
  
  x = x_in

  ! avoid x exactly 0 or 1
  if (x < tol)        x = tol
  if (x > 1.d0 - tol) x = 1.d0 - tol

  i = ceiling(x*n)
end function ith_bin

end module mod_particle_collision
