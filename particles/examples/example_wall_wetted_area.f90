!> Assess wall loads using relativistic guiding centers.
!>
!> This script requires that markers are given as an input in a part_restart.h5 file. You will
!> also need the wall input. See the wiki page "How to assess wall loads with the particle tracker?"
!> for details.
!>
Program  example_wallload

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff
use mod_gc_relativistic
use mod_kinetic_relativistic
use mod_wall_collision
use hdf5_io_module
use constants, only: PI
use phys_module, only : sqrt_mu0_rho0, F0
                  
implicit none

! Set up the simulation variables
real(kind=8)                      :: timesteps(1) = [0.1d-11] 
real(kind=8)                      :: target_time, t, R, phi, v_part
integer(kind=4)                   :: n_part, i, j, k, l, n_steps, ifail, max_depth, wall_id
integer                           :: n_tri    
integer(HID_T)                    :: file_id                       
type(write_particle_diagnostics)  :: diag
real(kind=8),dimension(3)         :: pos_prev, wall_pos, xyz_tria, norm_tria, v21, v31, Bphi_v
real*8, allocatable :: iangle(:,:), l_part(:)

type(particle_kinetic_relativistic), allocatable :: prtkin(:)
real*8 :: rnd(1), psi, U, B(3), E(3)

type(octree_node) :: wall
type(octree_triangle), allocatable :: triangles(:)
real*8, allocatable                :: wtmp(:)

! --- Read wall with octree
max_depth = 6
call mod_wall_collision_init('wall.h5',max_depth,wall)

! Read wall to initialize particles
file_id = 1
call HDF5_open('wall.h5',file_id)
call HDF5_integer_reading(file_id,n_tri,'ntriangle')

allocate(wtmp(9*n_tri))
call HDF5_array1D_reading(file_id, wtmp, 'nodes')
call HDF5_close(file_id)

! Convert wall array into triangles
allocate(triangles(n_tri))
do i=1,n_tri
   triangles(i)%triangle_id = i
   triangles(i)%v0 = (/ wtmp( (i-1)*9 + 1 ), wtmp( (i-1)*9 + 2 ), wtmp( (i-1)*9 + 3 ) /)
   triangles(i)%v1 = (/ wtmp( (i-1)*9 + 4 ), wtmp( (i-1)*9 + 5 ), wtmp( (i-1)*9 + 6 ) /)
   triangles(i)%v2 = (/ wtmp( (i-1)*9 + 7 ), wtmp( (i-1)*9 + 8 ), wtmp( (i-1)*9 + 9 ) /)
end do

call sim%initialize(num_groups=1)

! Set events to write output data and stop the simulation.
! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time))), & 
     !event(diag,start=sim%time,step=1d-8),         &
     event(stop_action(),start=sim%time+7.d-4)]

! Run first event to read the JOREK fields
call with(sim, events, at=0.d0)

! ------------- Initialize particles -----------------------------------------------------------------------------
sim%groups(:)%Z    = 1
sim%groups(:)%mass = 2.d0 !< atomic mass units

n_part = n_tri

allocate(particle_gc_relativistic::sim%groups(1)%particles(n_tri))

do i=1, n_tri
  
  select type (p=>sim%groups(1)%particles(i))
  
  type is (particle_gc_relativistic)
  
   !--- Calculate wall normals
   v21 = triangles(i)%v1 - triangles(i)%v0 
   v31 = triangles(i)%v2 - triangles(i)%v0 

   norm_tria = [v21(2)*v31(3)- v21(3)*v31(2), v21(3)*v31(1)- v21(1)*v31(3), v21(1)*v31(2)- v21(2)*v31(1)]

   norm_tria = norm_tria / norm2(norm_tria)
   
   ! Assign particle position to triangle center
   xyz_tria = (triangles(i)%v0 + triangles(i)%v1 + triangles(i)%v2) / 3.d0

   !--- Move the point 1 mm away from the triangle
   xyz_tria  = xyz_tria + norm_tria * 0.001

   R   = sqrt(xyz_tria(1)**2.0 + xyz_tria(2)**2.0)
   phi = atan2(xyz_tria(2),xyz_tria(1)) 
   if (phi < 0) phi = phi  + 2.0*PI

   p%x = [R,xyz_tria(3),phi]

   call find_RZ(sim%fields%node_list, sim%fields%element_list, &
            p%x(1), p%x(2), &
            p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)
   
   if (ifail /= 0) write(*,*) 'Initial position of triangle id = ', i, ' not found in JOREK grid, move the wall inside'

   Bphi_v = [-sin(phi), -cos(phi), 0.d0] * F0 / R

   ! --- Particle momentum in the field direction (must be corrected to point towards the domain)
   p%weight = 1.0
   p%q      = 1
   p%p(1)   = sign(1.d0,dot_product(Bphi_v,norm_tria))   ! parallel momentum aligns with normal (towards plasma)
   p%p(1)   = p%p(1) * 2.014 * 6.9d4    ! thermal velocity for 100 eV deuterium
   p%p(2)   = 0.d0    ! No perp momentum, fully parallel

   end select
  
end do

allocate(iangle(1,n_part))
allocate(l_part(n_part))
iangle = 0
! ----------- end paticle initialization -----------------------------------------------------------------------


! Set up the diagnostics output
!diag = write_particle_diagnostics(filename='diag.h5',only=[1,2,6,12,13,14,15]) ! store total and kinetic energies, p_phi, ielm, phi, R, Z



call check_and_fix_timesteps(timesteps, events)

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)

  do i=1,1
    n_steps = nint((target_time - sim%time)/timesteps(i))

    select type (particles => sim%groups(i)%particles)
    type is (particle_gc_relativistic)	
       l_part = 0.d0
      !$omp parallel do default(private) &
      !$omp shared (i, n_steps, timesteps, sim, wall, iangle, l_part)
       do j=1,size(particles,1)

          if (l_part(j) > 6.0) cycle  !--- If the particle is able to move away X m from the wall, it is a wetted area
          do k=1,n_steps
             if(particles(j)%i_elm .le. 0) exit

             pos_prev = particles(j)%x

             call runge_kutta_fixed_dt_gc_push_jorek(sim%fields,sim%time,timesteps(i), &
                  sim%groups(i)%mass,particles(j))

             l_part(j) = l_part(j) + norm2(pos_prev-particles(j)%x)

             if (particles(j)%i_elm .le. 0) exit

             if (k < 3) cycle  ! don't check collisions during the initial steps to allow the particle to move away from the wall
             call mod_wall_collision_check(pos_prev, particles(j)%x, wall, wall_id, wall_pos, iangle(i,j))
             if(wall_id .gt. 0) then
                particles(j)%x      = wall_pos
                particles(j)%i_elm  = -wall_id
             end if
          end do
       end do
       !$omp end parallel do

    type is (particle_kinetic_relativistic)
       !$omp parallel do default(private) &
       !$omp shared (i, n_steps, timesteps, sim, wall, iangle)
       do j=1,size(particles,1)
          do k=1,n_steps
             if (particles(j)%i_elm .le. 0) exit
             pos_prev = particles(j)%x

             call volume_preserving_push_jorek(particles(j),sim%fields,sim%groups(i)%mass,sim%time,timesteps(i),ifail)
             if (particles(j)%i_elm .le. 0) exit
             call mod_wall_collision_check(pos_prev, particles(j)%x, wall, wall_id, wall_pos, iangle(i,j))
             if(wall_id .gt. 0) then
                particles(j)%x      = wall_pos
                particles(j)%i_elm  = -wall_id
             end if
          end do
       end do
       !$omp end parallel do
    end select

  enddo

  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo

call mod_wall_collision_free(wall)

call write_simulation_hdf5(sim, 'part_out.h5')

call mod_wall_collision_export(sim, 'wallload.h5', iangle)

deallocate(iangle)

! Finalize the simulation
call sim%finalize

end program example_wallload
