!> Benchmark program to test the neutral neutral collisions
program test_NNC

  use data_structure
  use mpi
  use mod_pcg32_rng
  use mod_random_seed
  use mod_interp
  use mod_parameters, only: n_degrees
  use basis_at_gaussian, only: initialise_basis
  use mod_particle_allocation, only: allocate_particles_for_sim
  use mod_particle_collision
  use mod_particle_sim
  use mod_particle_io, only: write_simulation_hdf5
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
  use particle_tracer, only: sim

  !$ use omp_lib

  implicit none

  type(pcg32_rng), dimension(:), allocatable :: rng
  type(count_action)                         :: counter
  type(projection), target                   :: projections
  type(event), target                        :: project_diagnostics
  type(type_neutral_collision)               :: neutral_collision

  integer, parameter :: n_diag=8 !< how many diagnostic fields to project
  !< rho_1, rho_2, T, p_ii (p_RR, p_ZZ, p_\phi\phi), 1 empty

  !parameters of the grid
  real*8, parameter :: R_0 = 200.d0, Z_0 = 0.d0, length = 1.d0
  integer, parameter :: n = 6 !10!2!100 ! number of nodes in r, z directions

  !variables
  integer :: i, i_elm, i_step, nstep_proj
  integer :: seed, i_rng, n_stream
  real*8 :: s, t, phi, R, Z
  real*8 :: RN(8), weight, T_av !< [K] average intial temperature of gas
  character(len=100) :: i_step_char, filename, MSD_file
  real*8 :: last_time
  !$ real*8 :: w(2), mmm(3)
  
  !distribution w_i(R,t) diagnostic
  real*8, allocatable :: w_iRt(:,:,:)  !< weights per species as a function of R and t, dim(species,R_bins,0:nstep) 
  real*8, allocatable :: t_arr(:)      !< timesteps, dim(0:nstep)
  integer, parameter  :: R_bins = 1000 !< number of bins along R for the w_i(R,t) diagnostic

  !conservation
  real*8 :: conserv_obj(6)

  ! --- start of code
  last_time = -1.d99
  !$ last_time = omp_get_wtime()
  !$ w(1) = last_time

  ! Start up MPI, jorek
  call sim%initialize()

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

  write(MSD_file,"(A,I1,A)") 'p',sim%my_id,'.csv'
  open(unit=13+sim%my_id, file=MSD_file, status='replace')

  ! Set up particle group characteristics
  sim%groups(1)%Z    = -2 !< for deuterium 1
  sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units
  !sim%groups(1)%ad   = read_adf11(sim%my_id,'12_h')

  ! setting up particles per MPI node
  allocate(particle_kinetic_leapfrog::sim%groups(1)%particles( ceiling(part_group_configs(1)%n_particles/sim%n_mpi) ))

  ! --- Setting up random numbers
  seed = random_seed()
  n_stream = 1
  !$ n_stream = omp_get_max_threads()
  write(*,*) "id, n_mpi, n_stream",sim%my_id, sim%n_mpi, n_stream
  allocate(rng(n_stream))
  do i=1,n_stream
    call rng(i)%initialize(1, seed, n_stream, i)
  end do

  ! --- initialize the neutral collision action
  call neutral_collision%initialize(sim)

  ! --- setting up initial particle array
  weight = particlesource/(sim%groups(1)%n_particles) !< abusing particlesource to mean total weight of all particles in the sim for easy input iterations
  T_av=1*EL_CHG/K_BOLTZ
  select type (p => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)  
    !$omp parallel do default(none)  &
    !$omp shared(sim, p, rng, T_av, weight)       &
    !$omp private(i_rng, RN, R, Z, s, t, i_elm)
    do i=1,size(p)
      !$ i_rng = omp_get_thread_num()+1
      call rng(i_rng)%next(RN)

      p(i)%q      = 0 !< for neutrals
      p(i)%weight = weight
      
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

    !setting the MSD particle
    p(1)%x = [R_0,Z_0,0.d0]
  end select

  !projections for feedback purposes
  projections = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                do_zonal = .false., calc_integrals=.false., to_vtk=.TRUE., to_h5 = .false., basename='projections',   &
                                index_h5=.false., nsub=10)
  
  allocate(projections%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, n_diag))
  projections%rhs = 0.d0

  project_diagnostics = new_event_ptr(projections,   start = sim%time)

  ! setup up w_i(R,t) diagnostic
  if(nout_projection < 1) then
    if(sim%my_id == 0) write(*,*) "nout_projection not set, producing output every i_step"
    nout_projection = 1
  end if
  nstep_proj = nstep/nout_projection !< intended integer division
  allocate(t_arr(0:nstep_proj))
  allocate(w_iRt(2,R_bins,0:nstep_proj))
  w_iRt = 0.d0
  t_arr = 0.d0
  i_step = 0
  call with(sim, counter)
  call fill_w_iRt(sim,w_iRt,t_arr,i_step)

  conserv_obj(:) = 0.d0
  call conservation_checks(sim,conserv_obj)

  i_step=0
  ! write(filename,"(A,I2.2)") "i_step_",i_step
  ! call write_simulation_hdf5(sim,filename)
  ! -------------------------------------- main loop
  do i_step=1,nstep
    write(i_step_char,*) "i_step=",i_step
    call log_block(sim%my_id, trim(i_step_char), last_time)

    sim%time = sim%time + tstep_particles*nstep_particles

    call conservation_checks(sim,conserv_obj)

    call log_block(sim%my_id, "Kinetic pusher")
    call push_particle(sim, rng, projections)

    call conservation_checks(sim,conserv_obj)

    call global_av_T(sim)

    if(mod(i_step,nout_projection)==0) then
      call log_block(sim%my_id, "Diagnostics", last_time)
      call with(sim, counter)
      call fill_w_iRt(sim,w_iRt,t_arr,i_step)
      call with(sim, project_diagnostics)
      ! call conservation_checks(sim,conserv_obj)
    end if

    call log_block(sim%my_id, "Neutral self collision", last_time)
    call neutral_collision%do(sim,tstep_particles*nstep_particles)

    ! ! printout sim%particles
    ! select type (pa => sim%groups(1)%particles)
    ! type is (particle_kinetic_leapfrog)
    !   do i=1,size(pa,1)
    !     write(*,"(A,I10,7es20.10)") "DBG: pa(i) i,%x,%v",i,pa(i)%x,pa(i)%v
    !   end do
    ! end select

    ! if(mod(i_step,5000)==0) then
    !   call log_block(sim%my_id, "")
    !   call reset_part_labels(sim)
    ! end if
    ! call conservation_checks(sim,conserv_obj)

    ! write(filename,"(A,I2.2)") "i_step_",i_step
    ! call write_simulation_hdf5(sim,filename)  

    select type (pa => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)  
    do i=1,1
      write(13+sim%my_id,"(2I10, 8es20.10)") i, i_step, sim%time, pa(i)%x, pa(i)%v
    end do
    end select
  end do

  ! --- end
  call log_block(sim%my_id, "End of sim", last_time)
  call write_diag_output(sim%my_id,w_iRt,t_arr)
  close(13+sim%my_id)
  !$ w(2) = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w(2)-w(1))
  !$ if (sim%my_id .eq. 0) write(*,"(A,3es13.3,A)") "test_NNC done in (min/mean/max) ", mmm, " s"
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
      !$omp parallel do default(none)  &
      !$omp shared(sim, p, tstep_particles, mass, rng, nstep_particles) &
      !$omp private(i_rng, particle_tmp, to_be_reflected, vector_normal, &
      !$omp rz_old, st_old, i_elm_old, HH, HH_s, HH_t, HZ, l, m, qty,    &
      !$omp factor, i_tor, k, ifail, RN, j) &
      !$omp reduction(+:proj)
      do i=1,size(p)
        i_rng=1
        !$ i_rng = omp_get_thread_num()+1
        call copy_particle_kinetic_leapfrog(p(i),particle_tmp)

        do j=1,nstep_particles
          to_be_reflected = particle_tmp%i_elm < 0 

          if (particle_tmp%i_elm .gt. 0) then
            rz_old    = particle_tmp%x(1:2)
            st_old    = particle_tmp%st
            i_elm_old = particle_tmp%i_elm
            
            !contribution to diagnostics
            call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
            call mode_moivre(particle_tmp%x(3), HZ)
            
            do l=1,n_vertex_max ! corners of the elements
              do m=1,n_order+1  ! order of the polynomial
                qty = 0.d0
                factor = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_tmp%weight / nstep_particles
                qty(particle_tmp%i_life) = factor ! species density [m^-3] (qty(1) = n_1, qty(2) = n_2)
                qty(3)   = factor * mass * dot_product(particle_tmp%v,particle_tmp%v) / (3*EL_CHG) ! T [eV] (E = 1/2 m v^2 = 3/2 k_B T, k_B T / e = m v^2 / (3 e))
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
            
            !reflections of domain boundary
            if(particle_tmp%i_elm < 0) then
              particle_tmp%i_elm = -particle_tmp%i_elm
              vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, particle_tmp%i_elm, particle_tmp%st(1), particle_tmp%st(2))
              call rng(i_rng)%next(RN)
              particle_tmp%v = norm2(particle_tmp%v) * sample_cosine(RN(1:2), vector_normal)
            end if
            
            if(to_be_reflected .and. particle_tmp%i_elm < 0) then
              !$omp critical
               write(*,*) "ERROR: particle lost a (i,x,st,i_elm,ifail)",i,particle_tmp%x,particle_tmp%st,particle_tmp%i_elm,ifail
              !$omp end critical  
            end if
          else
            !$omp critical
              write(*,*) "ERROR: particle lost b",i,particle_tmp%x,particle_tmp%st,particle_tmp%i_elm
            !$omp end critical  
          end if
        end do ! loop over particle steps

        call copy_particle_kinetic_leapfrog(particle_tmp, p(i))

      end do ! loop over particles
      !$omp end parallel do

      ! rescaling the projections using the total weight at the location
      do l=1,n_vertex_max
        do m=1,n_order+1
          do i_elm=1,sim%fields%element_list%n_elements
            do i_tor=1,n_tor
              w = proj(m,l,i_elm,i_tor,1)+proj(m,l,i_elm,i_tor,2)
              proj(m,l,i_elm,i_tor,3) = proj(m,l,i_elm,i_tor,3)/max(proj(m,l,i_elm,i_tor,8),1.d0)/sim%n_mpi !< divide by total weight
            enddo
          enddo
        enddo
      enddo

      projections%rhs = proj

    end select
  end subroutine push_particle

  !> writes a new header in the log file, and if last_time is supplied, it times the previous block
  subroutine log_block(id,what,last_time)
    implicit none
    
    integer, intent(in) :: id
    character(len=*),intent(in) :: what !< what to call the new block
    real*8, intent(inout), optional :: last_time !< omp_get_wtime

    real*8 :: now, mmm(3)
    logical :: timing !< whether to write and update timing info, default .false., set to .true. if last_time is specified

    timing = .false.
    if(present(last_time)) timing = .true.
    
    if(timing) then
      !$ now = omp_get_wtime()
      !$ mmm = mpi_minmeanmax(now - last_time)
      !$ last_time = now
    endif

    if(id .ne. 0) return

    !$ if(timing) write(*,"(A,3f9.4,A)") "block done in (min/mean/max) ", mmm, " s"

    write(*,'(A100)') "===================================================================================================="
    write(*,*) what
    write(*,'(A100)') "===================================================================================================="
  
  end subroutine

  subroutine fill_w_iRt(sim,w_iRt,t_arr,i_step)
    implicit none

    type(particle_sim), intent(inout)  :: sim
    real*8, allocatable, intent(inout) :: w_iRt(:,:,:), t_arr(:)
    integer, intent(in) :: i_step

    integer :: iR, it
    real*8  :: R
    real*8, allocatable :: w_iR(:,:) !< w_iRt(:,:,it) for this diagnostic timestep

    it = i_step/nout_projection !< intended integer division
    if(sim%my_id == 0) write(*,*) "appending t_arr,w_iRt at it=",it
    t_arr(it) = sim%time
    
    w_iR=w_iRt(:,:,it)

    select type (p => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)  
      !$omp parallel do default(none)  &
      !$omp shared(sim, p) &
      !$omp private(iR, R) &
      !$omp reduction(+:w_iR)
      do i=1,size(p)
        if(p(i)%i_elm .le. 0) then
          ! !$omp critical
          ! write(*,*) "ERROR w_iRt",sim%my_id,p(i)%i_elm,p(i)%x
          ! !$omp end critical  
          cycle
        end if
        R = p(i)%x(1)
        iR = nint(((R - (R_0 - length)) / (2*length)) * R_bins + 0.5d0)
        if(iR < 1) then
          !$omp critical
          write(*,*) "error in fill_w_iRt, iR < 1, iR/i/R=",iR,i,R
          !$omp end critical  
          iR = 1
        end if
        if(iR > R_bins) then
          !$omp critical
          write(*,*) "error in fill_w_iRt, iR > R_bins, iR/i/R=",iR,i,R
          !$omp end critical  
          iR = R_bins
        end if
        

        w_iR(p(i)%i_life,iR) = w_iR(p(i)%i_life,iR) + p(i)%weight
      end do
      !$omp end parallel do
    end select

    w_iRt(:,:,it) = w_iR

  end subroutine fill_w_iRt

  subroutine write_diag_output(id, w_iRt_in,t_arr)
    implicit none

    integer, intent(in) :: id
    real*8, allocatable, intent(in) :: w_iRt_in(:,:,:), t_arr(:)
    
    real*8, allocatable :: w_iRt_tot(:,:,:)
    integer :: i,iR,it,ierr,size,nstep_proj

    allocate(w_iRt_tot, mold=w_iRt_in)
    nstep_proj = nstep/nout_projection
    size = 2*(nstep_proj + 1)*R_bins
    call MPI_REDUCE(w_iRt_in, w_iRt_tot, size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if(id == 0) then

      open(unit=10, file='w_iRt.csv', status='replace')
      do i=1,2
        do iR=1,R_bins
          do it=0,nstep_proj
            write(10,"(3I10, es20.10)") i, iR, it, w_iRt_tot(i,iR,it)
          end do
        end do
      end do
      close(10)

      open(unit=11, file='t_arr.csv', status='replace')
      do it=0,nstep_proj
        write(11,*) t_arr(it)
      end do
      close(11)
    
    end if
  end subroutine

  subroutine conservation_checks(sim, conserv_obj)
    implicit none
    
    integer,parameter :: num=6
    class(particle_sim), target, intent(in) :: sim
    real*8, intent(inout) :: conserv_obj(num) ! superparticles, real particles, momentum (3 directions), energy
    
    real*8  :: old(num), diff(num)
    integer :: j, ierr
    real*8    :: particles_remaining, momentum_remaining(3), energy_remaining, all_particles, all_momentum(3), all_energy
    integer   :: superparticles_remaining,all_superparticles,closest_iteration!, part_i_save,part_n_save
    
    
    particles_remaining = 0.d0
    momentum_remaining  = 0.d0
    energy_remaining    = 0.d0
    superparticles_remaining = 0 !< just for debugging. Use count_action in actual simulation.!.d0!.d0

  
    select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
      !$omp parallel do default(none) &
      !$omp shared(sim, particles) &
#endif
      !$omp reduction(+:particles_remaining, momentum_remaining, energy_remaining,superparticles_remaining)
        do j=1,size(particles,1)
  
          if (particles(j)%i_elm .le. 0) cycle
  
          particles_remaining = particles_remaining + particles(j)%weight
          momentum_remaining  = momentum_remaining  + particles(j)%weight * particles(j)%v *sim%groups(1)%mass * ATOMIC_MASS_UNIT
          energy_remaining    = energy_remaining    + particles(j)%weight * 0.5d0 * sim%groups(1)%mass * ATOMIC_MASS_UNIT * dot_product(particles(j)%v,particles(j)%v)
          superparticles_remaining = superparticles_remaining + 1
  
        enddo !j
      !omp end parallel do
    end select
  
    call MPI_REDUCE(particles_remaining, all_particles,         1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(momentum_remaining,  all_momentum,          3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(energy_remaining,    all_energy,            1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(superparticles_remaining,all_superparticles,1, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
    
    if (sim%my_id .eq. 0) then
      old = conserv_obj

      conserv_obj(1:2) = [real(all_superparticles),all_particles]
      conserv_obj(3:5) = all_momentum
      conserv_obj(6)   = all_energy

      diff = 0.d0
      diff = conserv_obj - old
      
      write(*,"(A)") "conservation checks --------------------------------------------------------------"
      write(*,"(A,6A15)") "qty: ", "superparticles", "particles", "momentum R", "momentum Z", "momentum phi", "energy"
      write(*,"(A,6es15.5)") "diff ",diff
      write(*,"(A,6es15.5)") "new  ",conserv_obj
      
  
    endif !(sim%my_id .eq. 0)
  end subroutine conservation_checks

  !> calculates and prints the global average temperature (in eV)
  subroutine global_av_T(sim)
    implicit none
    type(particle_sim), intent(inout)  :: sim

    real*8 :: T_av_measured, T_av_measured_red
    integer :: ierr
    
    T_av_measured = 0
    select type (p => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)  
      !$omp parallel do default(none)  &
      !$omp shared(p,sim) reduction(+:T_av_measured)
      do i=1,size(p)
        T_av_measured = T_av_measured + (sim%groups(1)%mass * ATOMIC_MASS_UNIT) * dot_product(p(i)%v,p(i)%v) / (3*EL_CHG)
      end do
      !$omp end parallel do
      T_av_measured = T_av_measured/size(p)
    end select
    call MPI_REDUCE(T_av_measured, T_av_measured_red,       1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(sim%my_id==0) then
      T_av_measured_red = T_av_measured_red/sim%n_mpi
      write(*,"(A,es15.5)") "Average measured temperature (eV)",T_av_measured_red
    end if
  end subroutine global_av_T

  subroutine reset_part_labels(sim)
    implicit none

    class(particle_sim), target, intent(inout) :: sim

    select type (p => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)  
      !$omp parallel do default(none)  &
      !$omp shared(sim, p, rng, T_Av, weight)       &
      !$omp private(i_rng, RN, R, Z, s, t, i_elm)
      do i=1,size(p)
        if(p(i)%x(1) < R_0) then
          p(i)%i_life = 1
        else
          p(i)%i_life = 2
        end if
      end do
      !$omp end parallel do
    end select
  end subroutine reset_part_labels

end program test_NNC
