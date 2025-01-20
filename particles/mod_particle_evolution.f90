!> Handles the evolution of the particle groups per time step
module mod_particle_evolution
    use particle_tracer
    use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
    use phys_module, only: nstep_particles, use_manual_random_seed
    use mod_coupling_settings
    use mod_project_particles
    use mod_random_seed
    use mod_interp, only: mode_moivre
    use mod_basisfunctions
    use mod_particle_types, only: copy_particle_kinetic_leapfrog
    use mod_sampling, only: boxmueller_transform,sample_chi_squared_3
    !$ use omp_lib

    implicit none
    real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)


contains


  !> Evolves the particles with every interaction that only depends on the particles except for wall reflections,
  !> 1. calculates interaction of the group with its environment (e.g. CX, ionisation... for neutrals)
  !> 2. pushes the particle
  !> 3. calculates the feedback for coupling with the fluid plasma
  subroutine evolve_particle_group(sim, group_num, jorek_feedback, rng, tstep_part_adj)

    use coupling_variables

    implicit none

    class(particle_sim), target, intent(inout)                :: sim
    integer, intent(in)                                       :: group_num
    type(projection), target, intent(inout)                   :: jorek_feedback
    type(count_action)                                        :: counter
    type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
    type(particle_kinetic_leapfrog)                           :: particle_tmp
  
    real*8, intent(in)     :: tstep_part_adj
  
    real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
    real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
    real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2)
    real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
    real*8    :: ion_rate, ion_source, ion_prob, ion_ran(1), cx_ran(8),st_ran(2), cx_source, cx_energy ,PLT
    real*8    :: cx_prob, CX_rate
    real*8    :: kinetic_energy, ion_energy,line_rad_energy
    real*8    :: n_lost_ion, n_lost_ion_all, p_plt_lost,p_plt_lost_all,p_cx_lost,p_cx_lost_all,p_lost_ion,p_lost_ion_all
    integer   :: n_super_ionized, n_super_ionized_all
    real*8    :: particle_source, velocity_par_source, energy_source
    real*8    :: v_temp(3), K_eV, v_kin_temp, B_norm(3), v, v_v, v_E,extra_proj
    real*8    :: vvector(3),sum_ran(3), E_th, v_th,ran_norm(4)
    !$ real*8 :: w0, w1, mmm(3)
  
    integer   :: i, j, k, l, m, i_elm_old, i_elm 
    integer   :: seed, i_rng, n_stream, ierr, nthreads
    integer   :: i_tor, index_lm, i_elm_temp
    integer   :: ifail
    logical   :: limits
    real*8,allocatable :: feedback_rhs(:,:,:,:,:)
  
    !$ w0 = omp_get_wtime()
  
    n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
    rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
    t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
    v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
    E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
    M_norm   = rho_norm * v_norm                                    ! momentum normalisation
  
    n_lost_ion = 0.d0
    n_lost_ion_all = 0.d0
    p_lost_ion   = 0.d0
    p_lost_ion_all   = 0.d0
    p_plt_lost  = 0.d0
    p_plt_lost_all  = 0.d0
    p_cx_lost   = 0.d0
    p_cx_lost_all   = 0.d0
    
    n_super_ionized = 0
    n_super_ionized_all = 0
      
      
    jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + nstep_particles * tstep_part_adj
  
    allocate(feedback_rhs,source=jorek_feedback%rhs)
  
    jorek_feedback%rhs = 0.d0
    feedback_rhs       = 0.d0
  
    call with(sim, counter)
  
    select type (particles => sim%groups(group_num)%particles)
    type is (particle_kinetic_leapfrog)
    if(use_manual_random_seed) then
      !$ call omp_set_schedule(omp_sched_static,10)
    else
      !$ call omp_set_schedule(omp_sched_dynamic,10)
    end if
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
    !$omp parallel do default(none) &
    !$omp shared(sim, group_num, particles, nstep_particles, tstep_part_adj, rng,                              &
    !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm,                                         &
    !$omp rho_idx_kn, Vpar_idx_kn, T_idx_kn,                                                        &
  
    !$omp CENTRAL_DENSITY, CENTRAL_MASS)                                                            &
#endif
    !$omp schedule(runtime)                                                                         &
    !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old,                  &
    !$omp i_elm_old, i_elm, n_e, T_e,                                                               &
    !$omp PLT,ion_rate, ion_prob, ion_ran, ion_source, ion_energy, kinetic_energy, line_rad_energy, &
    !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail,limits,           &
    !$omp CX_rate, CX_prob, CX_source, CX_energy, v, v_E, v_v,extra_proj,                           &
    !$omp particle_source, velocity_par_source, energy_source, v_temp, K_eV, cx_ran,                &
    !$omp E_th, v_th,sum_ran,vvector,ran_norm)                                                      &
    !$omp reduction(+:feedback_rhs,n_lost_ion,p_plt_lost,p_cx_lost,p_lost_ion,n_super_ionized)
    
    do j=1,size(particles,1)
  
      call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)
      !$ i_rng = omp_get_thread_num()+1
      
      do k=1,nstep_particles
  
        if (particle_tmp%i_elm .le. 0) exit
  
        t = sim%time + (k-1)*tstep_part_adj
  
        call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
        rz_old    = particle_tmp%x(1:2)
        st_old    = particle_tmp%st
        i_elm_old = particle_tmp%i_elm
        
        !> in use ionisation as well?
        call sim%fields%calc_NeTe(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), n_e, T_e)
        ! call calc_ U ne Te vpar
        ion_source = 0.d0
        ion_energy = 0.d0
        
        call sim%fields%calc_vvector(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), vvector)
        !vvector is fluid flow velocity [v_R, v_Z, v_phi] m/s
        !TODO: add upper limits if necessary
        limits = n_e .le. 1e14 .or. T_e * K_BOLTZ / EL_CHG .le. 1.d0 !ADAS limits
        if (particle_tmp%weight .lt. 0.0d0) write(*,*) "Negative particle weight p(j)%w=", particle_tmp%weight
        
        !>for impurities, bremsstrahlung and CX radiation can be added here as well. (see W_rad_example)
        line_rad_energy = 0.d0
        if (sim%groups(group_num)%use_kn_line_radiation .and. .not. limits) then !< before or after Ionisation and CX ??
          call sim%groups(group_num)%ad%PLT%interp(int(particle_tmp%q), log10(n_e), log10(T_e), PLT) ! [J m^3/s]
          ! call ad_deuterium%plt%interp( 1, ne_si_log10, Te_si_log10, LradDrays_T, dLradDrays_dT)
          line_rad_energy = n_e * particle_tmp%weight * PLT * tstep_part_adj
        endif ! use_kn_line_radiation
        
        if (sim%groups(group_num)%use_kn_ionisation .and. .not. limits) then
          
          call sim%groups(group_num)%ad%SCD%interp(int(particle_tmp%q), log10(n_e), log10(T_e), ion_rate) ! [m^3/s]
          ion_prob = 1.d0 - exp(-ion_rate * n_e * tstep_part_adj) ! [0] poisson point process, exponential 
  
          ! If the weight is too small throw away the particle with the probability, else reduce weight with ionising probability
          ion_source = 0.d0
  
          if (particle_tmp%weight .le. 1.0d9) then
            call rng(i_rng)%next(ion_ran)
            if (ion_ran(1) .le. ion_prob) then
              particle_tmp%i_elm  = 0
              ion_source = particle_tmp%weight
              !superparticles ionized
              n_super_ionized = n_super_ionized +1
            else
              ion_source = 0.d0
            endif
  
          else 
            ion_source = particle_tmp%weight * ion_prob
            particle_tmp%weight = particle_tmp%weight * (1.d0 - ion_prob)
          endif 
  
          kinetic_energy = dot_product(particle_tmp%v,particle_tmp%v) *sim%groups(group_num)%mass * ATOMIC_MASS_UNIT /2.d0
  
          ion_energy     = kinetic_energy - binding_energy !<binding energy should be here
          !<including binding energy will make ion_energy negative, so it becomes a sink for the plasma
  
        endif ! use_kn_ionisation
  
        
        ! Charge Exchange
        ! It is assumed that we will have a exchange between hydrogen isotopes
        v_temp    = particle_tmp%v
        cx_source = 0.d0
        cx_energy = 0.d0
        
        if (sim%groups(group_num)%use_kn_cx  .and. .not. limits) then !< CX uses adas as well. Te limit could be lower.
    
          call sim%groups(group_num)%ad%CCD%interp(int(particle_tmp%q+1), log10(n_e), log10(T_e), CX_rate) ! [m^3/s]
          CX_prob = 1.d0 - exp(-CX_rate * n_e * tstep_part_adj)
  
          call rng(i_rng)%next(cx_ran)
          if (cx_ran(1) .le. CX_prob) then
            ! sample boltzman, randomize velocity
  
            !============== NEW CX PARTICLE
            !Box-Mueller sample velocities with st.dev=1
            ran_norm = boxmueller_transform(cx_ran(2:5))
            v_temp = sqrt(T_e * K_BOLTZ/(sim%groups(group_num)%mass * ATOMIC_MASS_UNIT))*ran_norm(2:4)
  
            !>add bulk fluid flow
            v_temp = v_temp + vvector 
  
            CX_source = particle_tmp%weight
            CX_energy   = 0.5d0 * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT *  (dot_product(particle_tmp%v,particle_tmp%v) - dot_product(v_temp,v_temp))
          
  
            endif ! cx_ran
        endif ! use_kn_cx
        
        if (isnan(ion_source * ion_energy + cx_source * cx_energy - line_rad_energy)) then
          write(*,*) "ion_energy", ion_energy
          write(*,*) "cx_energy", cx_energy
          write(*,*) "line_rad_energy", line_rad_energy
          particle_tmp%i_elm  = 0
          cycle !< don't feed this particle into the feedback
        endif
        
        
        ! feedback from each particle at each timestep
        energy_source       = ion_source * ion_energy + cx_source * cx_energy - line_rad_energy
        particle_source     = ion_source * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT !< mass source in SI
        velocity_par_source = ion_source * dot_product(B, particle_tmp%v) * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT &
              + CX_source  * dot_product(B, particle_tmp%v - v_temp) * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT 
            
        particle_tmp%v = v_temp 
        n_lost_ion = n_lost_ion + ion_source  !< local sum #particles lost due to ionisation
        p_lost_ion = p_lost_ion + ion_source * ion_energy
        p_plt_lost = p_plt_lost + line_rad_energy
        p_cx_lost  = p_cx_lost + cx_source * cx_energy
        !Calculate the projection of the ion source in real-time
        call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
        call mode_moivre(particle_tmp%x(3), HZ)
            
        do l=1,n_vertex_max
          do m=1,n_order+1
  
            index_lm = (l-1)*(n_order+1) + m
  
            v   = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_source     * t_norm / rho_norm
            v_E = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source       * t_norm / E_norm
            v_v = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * velocity_par_source * t_norm / m_norm
            extra_proj = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) *particle_tmp%weight * 1.d0/real(nstep_particles,8) !<average density over jorek tstep_part_adj!real(floor(k/nstep_particles))!1.d0 !<density proj
  
            do i_tor=1,n_tor
              feedback_rhs(m,l,i_elm_old,i_tor,rho_idx_kn) = feedback_rhs(m,l,i_elm_old,i_tor,rho_idx_kn) + HZ(i_tor) * v
              feedback_rhs(m,l,i_elm_old,i_tor,Vpar_idx_kn) = feedback_rhs(m,l,i_elm_old,i_tor,Vpar_idx_kn) + HZ(i_tor) * v_v
              feedback_rhs(m,l,i_elm_old,i_tor,T_idx_kn) = feedback_rhs(m,l,i_elm_old,i_tor,T_idx_kn) + HZ(i_tor) * v_E
              feedback_rhs(m,l,i_elm_old,i_tor,4) = feedback_rhs(m,l,i_elm_old,i_tor,4) + HZ(i_tor) * extra_proj !< buiten de steps loop
            enddo
          enddo
        enddo
        
        
        if (particle_tmp%i_elm .gt. 0) then
          ! Push the particle and determine it's new location.
          call boris_push_cylindrical(particle_tmp, sim%groups(group_num)%mass, E, B, tstep_part_adj)
  
          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                              particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
        end if
      
      end do ! steps 
  
      call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))
  
      
    end do   ! particles
    !$omp end parallel do
      
    end select
  
    if (use_ncs) then
      write(*,*) 'GATHER TIME : ',jorek_feedback%rhs_gather_time
      jorek_feedback%rhs(:,:,:,:,rho_idx_kn) = feedback_rhs(:,:,:,:,rho_idx_kn) / jorek_feedback%rhs_gather_time !* TWOPI
      jorek_feedback%rhs(:,:,:,:,Vpar_idx_kn) = feedback_rhs(:,:,:,:,Vpar_idx_kn) / jorek_feedback%rhs_gather_time !* TWOPI
      jorek_feedback%rhs(:,:,:,:,T_idx_kn) = feedback_rhs(:,:,:,:,T_idx_kn) / jorek_feedback%rhs_gather_time !* TWOPI
      jorek_feedback%rhs(:,:,:,:,4) = feedback_rhs(:,:,:,:,4) ! additional project for density, will be moved to diagnostics in the future
      jorek_feedback%rhs_gather_time = 0.d0
    else
      jorek_feedback%rhs = feedback_rhs 
    endif
      
    deallocate(feedback_rhs)
  
    call MPI_REDUCE(n_lost_ion, n_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(p_lost_ion, p_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(p_plt_lost, p_plt_lost_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(p_cx_lost, p_cx_lost_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(n_super_ionized, n_super_ionized_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
    if (sim%my_id .eq. 0) write(*,'(A46,E14.6,I6)') "Lost superparticles at t due to ionisation: ", sim%time, n_super_ionized_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') " Lost particles at t due to ionisation: ", sim%time, n_lost_ion_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') " Ionization rate at time t [#/s]: ", sim%time, n_lost_ion_all / (tstep_part_adj * nstep_particles)
    p_lost_ion_all = p_lost_ion_all / (tstep_part_adj * nstep_particles)
    p_plt_lost_all = p_plt_lost_all / (tstep_part_adj * nstep_particles)
    p_cx_lost_all = p_cx_lost_all / (tstep_part_adj * nstep_particles)
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Energy exchange to plasma [W] at t due to ionisation: ", sim%time, p_lost_ion_all ! energy gain
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost energy [W] at t due to line radiation: ", sim%time, p_plt_lost_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Energy exchange to plasma [W] at t due to CX radiation: ", sim%time, p_cx_lost_all
    if (sim%my_id .eq. 0) write(*,'(A17,5E14.6)') 'TOTAL Exchange , delta t: ' ,sim%time,p_lost_ion_all, -p_plt_lost_all, p_cx_lost_all, tstep_part_adj * nstep_particles
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Total energy exchange to plasma [W]: ", sim%time, p_lost_ion_all -p_plt_lost_all+ p_cx_lost_all
  
    ! if (sim%my_id .eq. 0) write(*,*) " Lost energy [J] at t due to line radiation: ", sim%time, p_plt_lost_all
    !$ w1 = omp_get_wtime()
    !$ mmm = mpi_minmeanmax(w1-w0)
    !$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"
  
  
    if (sim%my_id .eq. 0) write(*,*) 'done evolving particle group: ', group_num
  
  end subroutine 

end module mod_particle_evolution