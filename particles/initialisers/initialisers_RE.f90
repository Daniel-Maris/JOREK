module initialisers_RE
  use mod_particle_types
  use mod_particle_sim
  use mod_rng
  use mod_initialise_particles
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT, MASS_ELECTRON, TWOPI
  use mod_pusher_tools, only: get_orthonormals
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  implicit none

  contains

  ! Quick and rough function to sample markers based on RZ-coordinates
  pure function RZ_pdf(var) result(p)
  real*8, intent(in)  :: var(2) ! var(1)=j
  real*8              :: p
  real*8              :: minor_r
  real*8              :: R_ax, Z_ax

  R_ax = 10.d0
  Z_ax = 0.d0

  minor_r = sqrt((var(1)-Z_ax)**2 + (var(2)-R_ax)**2)

  p = 1 / (1 + minor_r)**2

  end function RZ_pdf
    
  subroutine basic_initialization(sim, rng, energy, pitch, std_energy)
    use phys_module, only: tstep_particles
    use mod_kinetic_relativistic
    use mod_sampling, only: boxmueller_transform

    type(particle_sim),                   intent(inout) :: sim
    class(type_rng),                      intent(in)    :: rng
    real*8,                               intent(in)    :: energy, pitch ! Kinetic energy in units of eV and pitch
    real*8,               optional,       intent(in)    :: std_energy
    real*8,               allocatable                   :: p_tot(:), p_par(:), p_perp(:)
    integer                                             :: j
    real(kind=8)                                        :: psi, U, gyro_angle
    real(kind=8),         dimension(3)                  :: E, B, B_cart, B_norm
    real*8                                              :: e1(3), e2(3) 
    real*8                                              :: R_bound(2), Z_bound(2)
    integer                                             :: num_part
    real*8,               allocatable                   :: ran_uniform(:), ran_gaussian(:)

    R_bound = [1.9d0, 3.3d0]
    Z_bound = [-1.d0, 0.8d0]
    ! Initialise every particle
    !call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[var_psi], transform=psi_pdf)
    call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-2,-1], transform=RZ_pdf)
    !call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-1,var_zj], transform=current_pdf) 
    !call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-2,-1], transform=analytical_pdf)
    !call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-1, 9], transform=normalized_nre_pdf, normalize=.false.) ! i_var = 7 for JET case, 9 for ITER case
    ! call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-1, 7], transform=new_nre_pdf, normalize=.true., Rbound=R_bound, Zbound=Z_bound) ! i_var = 7 for JET case, 9 for ITER case

    !call initialise_particles(sim%groups(1)%particles, sim%fields%node_list, sim%fields%element_list, rng, normalize=.false.) ! i_var = 7 for JET case, 9 for ITER case
    !call weigh_with_interp_f(sim%fields%node_list, sim%fields%element_list, sim%groups(1)%particles, [-2, -1], RZ_pdf)

    num_part = size(sim%groups(1)%particles,1) 

    allocate(p_tot(num_part))
    allocate(p_par(num_part))
    allocate(p_perp(num_part))

    ! Generate gassian distributed energy
    if (present(std_energy)) then

      allocate(ran_uniform(num_part + mod(num_part,2)))
      allocate(ran_gaussian(num_part + mod(num_part,2)))

      call random_number(ran_uniform)
      ran_gaussian = boxmueller_transform(ran_uniform)

      p_tot               = sqrt(((energy + std_energy*ran_gaussian(:num_part))*EL_CHG/SPEED_OF_LIGHT + MASS_ELECTRON*SPEED_OF_LIGHT)**2 - (MASS_ELECTRON*SPEED_OF_LIGHT)**2)/ATOMIC_MASS_UNIT ! [AMU*m/s]

      deallocate(ran_uniform)
      deallocate(ran_gaussian)

    else

      p_tot               = sqrt((energy*EL_CHG/SPEED_OF_LIGHT + MASS_ELECTRON*SPEED_OF_LIGHT)**2 - (MASS_ELECTRON*SPEED_OF_LIGHT)**2)/ATOMIC_MASS_UNIT ! [AMU*m/s]

    end if

    p_par               = pitch * p_tot
    p_perp              = (1-pitch) * p_tot

    ! Set particle momentum
    select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_relativistic)
      !$omp parallel do default(none) &
      !$omp private(E, B, psi, U, B_cart, B_norm, e1, e2, gyro_angle, j) &
      !$omp shared (sim, tstep_particles, p_par, p_perp, num_part)
      do j=1,num_part

        ! Extract magnetic field and convert to cartesian coordinates
        call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
        B_cart = vector_cylindrical_to_cartesian(particles(j)%x(3),B)

        B_norm = B_cart/norm2(B_cart)

        ! Generate perpendicular component based on sampled gyro angle
        call get_orthonormals(B_norm, e1, e2)
        call random_number(gyro_angle)
        gyro_angle = gyro_angle * TWOPI

        !particles(j)%p = p_par * B_norm + p_perp*(e1*cos(gyro_angle) + e2*sin(gyro_angle))
        particles(j)%p = -p_par(j) * B_norm + p_perp(j)*(e1*cos(gyro_angle) + e2*sin(gyro_angle))

      end do
      !$omp end parallel do 
    end select

    deallocate(p_tot)
    deallocate(p_par)
    deallocate(p_perp)

  end subroutine basic_initialization

end module initialisers_RE