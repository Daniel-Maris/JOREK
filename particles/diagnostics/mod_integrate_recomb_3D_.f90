module mod_integrate_recomb_3D

  implicit none
  contains

#include "corr_neg_include.f90"

  !-------------------------------------------------------------------------------------------------
  ! Subroutine: integrate_recombination
  ! Purpose:    Calculates the recombination rate, momentum, and energy source terms for each
  !             element in a 3D domain.
  !-------------------------------------------------------------------------------------------------
  subroutine integrate_recombination(my_id, n_mpi, rec_rate_local, rec_v_R, rec_v_Z, rec_v_phi, volume_check, energy_neutrals, energy_radiation)

    use data_structure
    use gauss
    use basis_at_gaussian
    use phys_module
    use nodes_elements
    use mpi
    use constants, only: TWOPI
    use mod_atomic_coeff_deuterium, only: rec_rate_to_kinetic
    !$ use omp_lib

    implicit none

    ! --- Argument Declarations ---
    integer :: n_mpi, my_id
    real*8, dimension(:, :), allocatable, intent(out) :: rec_rate_local              !< [out] Local recombination rate per element, per plane
    real*8, dimension(:, :), allocatable, intent(out) :: rec_v_R, rec_v_Z, rec_v_phi   !< [out] Recombination momentum source terms (R, Z, phi)
    real*8, dimension(:, :), allocatable, optional, intent(out) :: volume_check      !< [out, optional] Integrated volume of each element for verification
    real*8, dimension(:, :), allocatable, optional, intent(out) :: energy_neutrals   !< [out, optional] Neutral kinetic energy from recombination
    real*8, dimension(:, :), allocatable, optional, intent(out) :: energy_radiation  !< [out, optional] Radiated energy from recombination

    ! --- Local Type Declarations ---
    type (type_element) :: element
    type (type_node)    :: nodes(n_vertex_max)

    ! --- Local Variables ---
    ! Geometry and field values at Gaussian points
    real*8 :: x_g(n_gauss, n_gauss), x_s(n_gauss, n_gauss), x_t(n_gauss, n_gauss)
    real*8 :: y_g(n_gauss, n_gauss), y_s(n_gauss, n_gauss), y_t(n_gauss, n_gauss)
    real*8, dimension(n_plane, n_var, n_gauss, n_gauss) :: eq_g, eq_s, eq_t

    ! Loop counters and indices
    integer :: i, j, k, in, ms, mt, mp, iv, inode, ife, ielm, n_elements
    integer :: missing, loc_rec_elms

    ! Physical quantities at an integration point
    real*8 :: xjac, BigR, wst, delta_phi
    real*8 :: ps0_x, ps0_y, u0_x, u0_y, r0, T0, vpar0
    real*8 :: T0_corr, r0_corr

    ! Atomic rate coefficients
    real*8 :: Sion_T, dSion_dT
    real*8 :: Srec_T, dSrec_dT
    real*8 :: LradDcont_T, dLradDcont_T
    real*8 :: LradDcont_corr, dLradDcont_dT_corr

    ! MPI-related variables
    integer, dimension(:), allocatable :: local_rec_elements

    !================================================================================
    ! Subroutine Body
    !================================================================================

    ! --- Distribute elements across MPI processes for load balancing ---
    allocate(local_rec_elements(n_mpi))
    n_elements = element_list%n_elements
    missing = mod(n_elements, n_mpi) !< Number of elements left after even distribution

    ! Distribute elements as evenly as possible. The first 'missing' processes get one extra.
    do i = 1, n_mpi
      local_rec_elements(i) = floor(n_elements / real(n_mpi, 8))
      if (missing > 0) then
        missing = missing - 1
        local_rec_elements(i) = local_rec_elements(i) + 1
      endif
    enddo

    ! --- Allocate local result arrays if they are not already allocated ---
    if (.not. allocated(rec_rate_local)) then
      allocate(rec_rate_local(local_rec_elements(my_id + 1), n_plane))
      allocate(rec_v_R(local_rec_elements(my_id + 1), n_plane))
      allocate(rec_v_Z(local_rec_elements(my_id + 1), n_plane))
      allocate(rec_v_phi(local_rec_elements(my_id + 1), n_plane))
      allocate(volume_check(local_rec_elements(my_id + 1), n_plane))
      allocate(energy_neutrals(local_rec_elements(my_id + 1), n_plane))
      allocate(energy_radiation(local_rec_elements(my_id + 1), n_plane))
    endif

    ! --- Initialize local arrays to zero ---
    rec_rate_local(:, :)   = 0.d0
    rec_v_R(:, :)          = 0.d0
    rec_v_Z(:, :)          = 0.d0
    rec_v_phi(:, :)        = 0.d0
    volume_check(:, :)     = 0.d0
    energy_neutrals(:, :)  = 0.d0
    energy_radiation(:, :) = 0.d0

    ! --- Define the toroidal angle step for 3D integration ---
    delta_phi = 2.d0 * PI / real(n_plane, 8) / real(n_period, 8)

    !$omp parallel do default(none)                                                              &
    !$omp   schedule(static, 100)                                                                &
    !$omp   shared(local_rec_elements, my_id, n_mpi, element_list, node_list,                     &
    !$omp          H, H_s, H_t, HZ, tstep, F0, delta_phi, gamma,                                  &
    !$omp          rec_rate_local, rec_v_R, rec_v_Z, rec_v_phi,                                   &
    !$omp          volume_check, energy_neutrals, energy_radiation                                &
    !$omp          )                                                                              &
    !$omp   private(ife, ielm, iv, i, j, k, ms, mt, mp, in, inode,                                &
    !$omp           nodes, element, x_g, y_g, x_s, y_s, x_t, y_t,                                 &
    !$omp           xjac, eq_g, eq_s, eq_t, wst, BigR, r0, T0, vpar0,                             &
    !$omp           ps0_x, ps0_y, u0_x, u0_y, r0_corr, T0_corr,                                   &
    !$omp           Sion_T, dSion_dT, Srec_T, dSrec_dT, LradDcont_T,                              &
    !$omp           dLradDcont_dT, LradDcont_corr, dLradDcont_dT_corr                             &
    !$omp           )
    ! Loop over elements assigned to this MPI process
    do ife = 1, local_rec_elements(my_id + 1)

      ! Map local loop index 'ife' to global element index 'ielm'.
      ! This scheme distributes adjacent elements to different MPI ranks for better load balance.
      ielm = (my_id + 1) + n_mpi * (ife - 1)
      element = element_list%element(ielm)

      ! Create a thread-private deep copy of node data to prevent race conditions in OpenMP.
      do iv = 1, n_vertex_max
        inode = element%vertex(iv)
        call make_deep_copy_node(node_list%node(inode), nodes(iv))
      enddo

      ! --- Interpolate geometry and physical variables from nodes to Gaussian points ---
      x_g = 0.d0; x_s = 0.d0; x_t = 0.d0
      y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
      eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0

      do i = 1, n_vertex_max
        do j = 1, n_order + 1
          do ms = 1, n_gauss
            do mt = 1, n_gauss
              x_g(ms, mt) = x_g(ms, mt) + nodes(i)%x(1, j, 1) * element%size(i, j) * H(i, j, ms, mt)
              x_s(ms, mt) = x_s(ms, mt) + nodes(i)%x(1, j, 1) * element%size(i, j) * H_s(i, j, ms, mt)
              x_t(ms, mt) = x_t(ms, mt) + nodes(i)%x(1, j, 1) * element%size(i, j) * H_t(i, j, ms, mt)
              y_g(ms, mt) = y_g(ms, mt) + nodes(i)%x(1, j, 2) * element%size(i, j) * H(i, j, ms, mt)
              y_s(ms, mt) = y_s(ms, mt) + nodes(i)%x(1, j, 2) * element%size(i, j) * H_s(i, j, ms, mt)
              y_t(ms, mt) = y_t(ms, mt) + nodes(i)%x(1, j, 2) * element%size(i, j) * H_t(i, j, ms, mt)
            end do
          end do
          do ms = 1, n_gauss
            do mt = 1, n_gauss
              do k = 1, n_var
                do in = 1, n_tor
                  do mp = 1, n_plane
                    eq_g(mp, k, ms, mt) = eq_g(mp, k, ms, mt) + nodes(i)%values(in, j, k) * element%size(i, j) * H(i, j, ms, mt) * HZ(in, mp)
                    eq_s(mp, k, ms, mt) = eq_s(mp, k, ms, mt) + nodes(i)%values(in, j, k) * element%size(i, j) * H_s(i, j, ms, mt) * HZ(in, mp)
                    eq_t(mp, k, ms, mt) = eq_t(mp, k, ms, mt) + nodes(i)%values(in, j, k) * element%size(i, j) * H_t(i, j, ms, mt) * HZ(in, mp)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      ! --- Integrate over Gaussian points and toroidal planes ---
      do mp = 1, n_plane
        do ms = 1, n_gauss
          do mt = 1, n_gauss
            wst = wgauss(ms) * wgauss(mt)                               !< Gaussian weight
            xjac = x_s(ms, mt) * y_t(ms, mt) - x_t(ms, mt) * y_s(ms, mt) !< Jacobian of the coordinate transformation
            BigR = x_g(ms, mt)                                          !< Major radius at the integration point

            ! Gradients of poloidal flux (psi) and potential (u)
            ps0_x = (y_t(ms, mt) * eq_s(mp, 1, ms, mt) - y_s(ms, mt) * eq_t(mp, 1, ms, mt)) / xjac
            ps0_y = (-x_t(ms, mt) * eq_s(mp, 1, ms, mt) + x_s(ms, mt) * eq_t(mp, 1, ms, mt)) / xjac
            u0_x  = (y_t(ms, mt) * eq_s(mp, 2, ms, mt) - y_s(ms, mt) * eq_t(mp, 2, ms, mt)) / xjac
            u0_y  = (-x_t(ms, mt) * eq_s(mp, 2, ms, mt) + x_s(ms, mt) * eq_t(mp, 2, ms, mt)) / xjac

            ! Plasma density, temperature, and parallel velocity
            r0 = eq_g(mp, 5, ms, mt)
            r0_corr = corr_neg_dens1(r0)
            T0 = eq_g(mp, 6, ms, mt)
            T0_corr = corr_neg_temp1(T0)
            vpar0 = eq_g(mp, 7, ms, mt)

            ! Calculate atomic recombination and radiation rates
            call rec_rate_to_kinetic(r0, 0.5d0 * T0, Sion_T, dSion_dT, Srec_T, dSrec_dT, &
                                     LradDcont_T, dLradDcont_dT, LradDcont_corr, dLradDcont_dT_corr)
            
            ! Transform derivatives from electron temperature (Te) to total temperature (T)
            dSrec_dT = dSrec_dT / 2.d0
            dLradDcont_dT = dLradDcont_dT / 2.d0
            dLradDcont_dT_corr = dLradDcont_dT_corr / 2.d0

            ! --- Perform the integration using the correct 3D volume element (BigR * xjac * delta_phi) ---
            rec_rate_local(ife, mp)   = rec_rate_local(ife, mp) + (Srec_T * r0_corr * r0_corr) * BigR * xjac * tstep * delta_phi * wst !< Neutral density gain
            rec_v_R(ife, mp)          = rec_v_R(ife, mp) + (Srec_T * r0_corr * r0_corr) * (-BigR * u0_y + vpar0 / BigR * ps0_y) * BigR * xjac * tstep * delta_phi * wst !< Momentum source in R
            rec_v_Z(ife, mp)          = rec_v_Z(ife, mp) + (Srec_T * r0_corr * r0_corr) * (+BigR * u0_x - vpar0 / BigR * ps0_x) * BigR * xjac * tstep * delta_phi * wst !< Momentum source in Z
            rec_v_phi(ife, mp)        = rec_v_phi(ife, mp) + (Srec_T * r0_corr * r0_corr) * (+F0 * vpar0 / BigR) * BigR * xjac * tstep * delta_phi * wst !< Momentum source in phi
            volume_check(ife, mp)     = volume_check(ife, mp) + (1.d0) * BigR * xjac * delta_phi * wst !< Element volume contribution
            energy_neutrals(ife, mp)  = energy_neutrals(ife, mp) + (gamma - 1.d0) * 0.5d0 * T0 * r0_corr * r0_corr * Srec_T * BigR * xjac * tstep * delta_phi * wst !< Kinetic energy gain
            energy_radiation(ife, mp) = energy_radiation(ife, mp) + r0_corr * r0_corr * LradDcont_corr * BigR * xjac * tstep * delta_phi * wst !< Radiated energy loss
          enddo
        enddo
      enddo

      ! Free the memory of the private node copy to prevent memory leaks
      do iv = 1, n_vertex_max
        call dealloc_node(nodes(iv))
      enddo

    enddo ! End of loop over local elements
    !$omp end parallel do

  end subroutine

end module mod_integrate_recomb_3D