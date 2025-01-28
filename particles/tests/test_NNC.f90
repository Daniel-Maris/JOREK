!> Benchmark program to test the neutral neutral collisions
program test_NNC

  use data_structure
  use mpi
  use mod_pcg32_rng
  use mod_random_seed
  use mod_interp
  use mod_parameters, only: n_degrees
  use basis_at_gaussian, only: initialise_basis
  use mod_particle_collision
  use mod_particle_sim
  use phys_module
  use mod_atomic_elements
  use mod_particle_types
  use nodes_elements
  use mod_boundary,   only: boundary_from_grid
  use equil_info
  use mod_fields_linear, only: jorek_fields_interp_linear
  use constants
  use mod_sampling, only: boxmueller_transform
  use mod_event, only: mpi_minmeanmax, count_action, with, event, new_event_ptr
  use mod_log_params
  use mod_project_particles

  !$ use omp_lib

  implicit none

  type(particle_sim)                         :: sim
  type(pcg32_rng), dimension(:), allocatable :: rng
  type(count_action)                         :: counter
  type(projection), target                   :: projections
  type(event), target                        :: project_diagnostics
  
  integer, parameter :: n_diag=8 !< how many diagnostic fields to project
  !< rho_1, rho_2, T, p_ii (p_RR, p_ZZ, p_\phi\phi), 1 empty

  !parameters of the grid
  real*8, parameter :: R_0 = 100.d0, Z_0 = 0.d0, length = 0.1d0
  integer, parameter :: n = 10 !100 ! number of nodes in r, z directions

  !variables
  integer :: i, i_elm, i_step
  integer :: seed, i_rng, n_stream
  real*8 :: s, t, phi, R, Z
  real*8 :: RN(8), weight, T_av !< [K] average intial temperature of gas
  character(len=100) :: i_step_char
  !$ real*8 :: w(2), mmm(3)

  ! --- start of code

  !$ w(1) = omp_get_wtime()

  ! Start up MPI, jorek
  call sim%initialize(num_groups=1)

  ! --- Write out all parameters defined in parameters and the namelist input file.
  call log_parameters(sim%my_id)

  ! Setup the grid
  call init_node_list(node_list, n_nodes_max, 0, n_var)
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n, n, R_0-length,R_0+length, Z_0-length, Z_0+length, .true., node_list, element_list)

  ! setup sim
  allocate(jorek_fields_interp_linear::sim%fields)
  if (.not. associated(sim%fields%node_list))    sim%fields%node_list    => node_list
  if (.not. associated(sim%fields%element_list)) sim%fields%element_list => element_list
  
  if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

  ! dummy fields
  do i=1,size(sim%fields%node_list%node)
    sim%fields%node_list%node(i)%values(:,:,:) = 0
    sim%fields%node_list%node(i)%deltas(:,:,:) = 0
  end do
  
  call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase )

  ! --- Plot the grid
  if (sim%my_id == 0) then
    call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.true.,.false.,'initial')
  end if

  ! Set up particle group characteristics
  sim%groups(1)%Z    = -2 !< for deuterium 1
  sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units
  !sim%groups(1)%ad   = read_adf11(sim%my_id,'12_h')

  ! setting up particles per MPI node
  allocate(particle_kinetic_leapfrog::sim%groups(1)%particles( ceiling(n_particles/sim%n_cpu) ))

  ! --- Setting up random numbers
  seed = random_seed()
  n_stream = 1
  !$ n_stream = omp_get_max_threads()
  write(*,*) "id, n_cpu, n_stream",sim%my_id, sim%n_cpu, n_stream
  allocate(rng(n_stream))
  do i=1,n_stream
    call rng(i)%initialize(1, seed, n_stream, i)
  end do

  ! --- setting up initial particle array
  weight = 1.d25/(n_particles)
  T_av=1*EL_CHG/K_BOLTZ
  select type (p => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)  
    !$omp parallel do default(none)  &
    !$omp shared(sim, p, rng, T_Av, weight)       &
    !$omp private(i_rng, RN, R, Z, s, t, i_elm)
    do i=1,size(p)
      !$ i_rng = omp_get_thread_num()+1
      call rng(i_rng)%next(RN)

      p(i)%q      = 0 !< for neutrals
      p(i)%weight = weight
      p(i)%t_birth = 0.d0
      
      s = RN(1)
      t = RN(2)
      i_elm = ceiling(sim%fields%element_list%n_elements * RN(3))
      p(i)%i_elm  = i_elm      
      p(i)%st     = [s,t]
      call interp_RZ(sim%fields%node_list,sim%fields%element_list,i_elm,s,t,R,Z)
      p(i)%x      = [R,Z,RN(4)*TWOPI/n_period]
      
      RN(5:8) = boxmueller_transform(RN(5:8))
      p(i)%v      = 0.d0 + sqrt(T_av*K_BOLTZ/(sim%groups(1)%mass * ATOMIC_MASS_UNIT))*RN(5:7)

      if(R < R_0) then
        p(i)%i_life = 1
      else
        p(i)%i_life = 2
      end if
    end do
    !$omp end parallel do
  end select

  !projections for feedback purposes
  projections = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                do_zonal = .false., calc_integrals=.false., to_vtk=.TRUE., to_h5 = .false., basename='projections',   &
                                index_h5=.false., nsub=10)
  
  allocate(projections%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, n_diag))
  projections%rhs = 0.d0

  project_diagnostics = new_event_ptr(projections,   start = sim%time)


  ! -------------------------------------- main loop
  do i_step=1,nstep
    write(i_step_char,*) "i_step=",i_step
    call write_to_outputfile(sim%my_id, trim(i_step_char))
    
    call write_to_outputfile(sim%my_id, "Kinetic pusher")
    call push_particle(sim, rng, projections)
    call with(sim, counter)
    call with(sim, project_diagnostics)

    call write_to_outputfile(sim%my_id, "Neutral self collision")
    call neutral_self_collision(sim,rng,tstep_particles*nstep_particles)

    sim%time = sim%time + tstep_particles*nstep_particles
  end do

  ! --- end
  call write_to_outputfile(sim%my_id, "End of sim")
  !$ w(2) = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w(2)-w(1))
  !$ if (sim%my_id .eq. 0) write(*,"(A,3f9.4,A)") "test_NNC done in (min/mean/max) ", mmm, " s"
  call sim%finalize()
contains

  subroutine push_particle(sim,rng,projections)
    use mod_sampling, only: sample_cosine
    use mod_find_rz_nearby
    use mod_boris, only: boris_push_cylindrical
    use mod_boundary, only: wall_normal_vector
    use mod_basisfunctions
    
    implicit none
    
    type(particle_sim), intent(inout) :: sim
    type(pcg32_rng), dimension(:), intent(inout) :: rng
    type(projection), target, intent(inout)   :: projections
    

    type(particle_kinetic_leapfrog) :: particle_tmp
    real*8 :: rz_old(2), st_old(2), vector_normal(3), RN(2)
    integer :: i_elm_old, ifail, i_rng, i, j

    !diagnostics
    integer :: m, l, k, i_tor, i_elm
    real*8 :: HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
    real*8 :: qty(n_diag), w, mass
    real*8 :: factor !< the contribution of this particle's weight to the number density as projected onto this order of the basis function
    real*8, allocatable :: proj(:,:,:,:,:)
    
    !dbg:
    logical :: to_be_reflected
    integer :: is,it,ns,nt
    real*8  :: a,b

    allocate(proj,source=projections%rhs)

    projections%rhs = 0.d0
    proj           = 0.d0
    
    mass = sim%groups(1)%mass * ATOMIC_MASS_UNIT

    select type (p => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)  
      ! !$omp parallel do default(private)  &
      ! !$omp shared(sim, p, tstep_particles, mass) &
      ! !$omp reduction(+:proj)
      do i=1,size(p)
        i_rng=1
        ! !$ i_rng = omp_get_thread_num()+1
        call copy_particle_kinetic_leapfrog(p(i),particle_tmp)

        do j=1,nstep_particles
          to_be_reflected = particle_tmp%i_elm < 0 

          !reflections of domain boundary
          if(particle_tmp%i_elm < 0) then
            particle_tmp%i_elm = -particle_tmp%i_elm
            vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, particle_tmp%i_elm, particle_tmp%st(1), particle_tmp%st(2))
            call rng(i_rng)%next(RN)
            particle_tmp%v = norm2(particle_tmp%v) * sample_cosine(RN(1:2), vector_normal)
          end if

          if (particle_tmp%i_elm .gt. 0) then
            rz_old    = particle_tmp%x(1:2)
            st_old    = particle_tmp%st
            i_elm_old = particle_tmp%i_elm
            
            !contribution to diagnostics
            
            !Calculate the projection of the ion source in real-time
            call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
            call mode_moivre(particle_tmp%x(3), HZ)
            
            do l=1,n_vertex_max ! corners of the elements
              do m=1,n_order+1  ! order of the polynomial
                qty = 0.d0
                factor = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_tmp%weight / nstep_particles
                qty(particle_tmp%i_life) = factor ! species density (qty(1) = n_1, qty(2) = n_2)
                qty(3)   = factor * mass * norm2(particle_tmp%v)**2 / (3*EL_CHG) ! T (E = 1/2 m v^2 = 3/2 k_B T, k_B T / e = m v^2 / (3 e))
                qty(4:6) = factor * mass * particle_tmp%v * particle_tmp%v ! directional pressure P_ii = dp_i/dt /A = mv_i^2 dt A n /dt /A = mv_i^2 n
                qty(7)   = factor ! total species density
                qty(8)   = particle_tmp%weight !< to normalise later on
                
                do i_tor=1,n_tor ! toroidal harmonic
                  do k=1,n_diag  ! diagnostic quantity
                    proj(m,l,i_elm_old,i_tor,k) = proj(m,l,i_elm_old,i_tor,k) + HZ(i_tor) * qty(k)
                  enddo
                enddo
              enddo
            enddo
            
            ! Push the particle and determine it's new location. Since q=0, the Boris method simplifies to v=constant, so this acts as a straight leapfrog pusher
            call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, [1.d0,0.d0,0.d0], [1.d0,0.d0,0.d0], tstep_particles)
            
            call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                                particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
            if(to_be_reflected .and. particle_tmp%i_elm < 0) then
              ! !$omp critical
              !  write(*,*) "ERROR: particle lost (i,x,st,i_elm,ifail)",i,particle_tmp%x,particle_tmp%st,particle_tmp%i_elm,ifail
              ! !$omp end critical  
            end if
          else
            ! !$omp critical
              write(*,*) "ERROR: particle lost?",i,particle_tmp%x,particle_tmp%st,particle_tmp%i_elm
            ! !$omp end critical  
          end if
        end do ! loop over particle steps

        call copy_particle_kinetic_leapfrog(particle_tmp, p(i))

      end do ! loop over particles
      ! !$omp end parallel do

      ! rescaling the projections using the total weight at the location
      do l=1,n_vertex_max
        do m=1,n_order+1
          do i_elm=1,sim%fields%element_list%n_elements
            do i_tor=1,n_tor
              w = proj(m,l,i_elm,i_tor,1)+proj(m,l,i_elm,i_tor,2)
              proj(m,l,i_elm,i_tor,3) = proj(m,l,i_elm,i_tor,3)/max(proj(m,l,i_elm,i_tor,8),1.d0) !< divide by total weight
            enddo
          enddo
        enddo
      enddo

      projections%rhs = proj

    end select
  end subroutine push_particle

  subroutine write_to_outputfile(id,what)
    implicit none
    
    integer, intent(in) :: id
    character(len=*),intent(in) :: what
  
    if(id .ne. 0) return
  
    write(*,'(A100)') "===================================================================================================="
    write(*,*) what
    write(*,'(A100)') "===================================================================================================="
  
  end subroutine
end program test_NNC
