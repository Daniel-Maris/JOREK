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
use mod_coordinate_transforms, only: cylindrical_to_cartesian, cartesian_to_cylindrical
use mod_wall_collision
use hdf5_io_module
use constants, only: PI
use phys_module, only : sqrt_mu0_rho0, F0
                  
implicit none

! Set up the simulation variables
real(kind=8)                      :: timesteps(1) = [1d-7] 
real(kind=8)                      :: target_time, t, R, phi, v_part, norm_R, dir_sign, tstep_field=1.d-2
integer(kind=4)                   :: n_part, i, j, k, l, ifail, max_depth, wall_id, n_steps = 1000
integer(kind=4),allocatable       :: indices(:,:)
integer                           :: n_tri, n_nodes
integer             :: filehandle = 60
integer(HID_T)                    :: file_id                       
type(write_particle_diagnostics)  :: diag
real(kind=8),dimension(3)         :: pos_prev, wall_pos, xyz_tria, norm_tria, v21, v31, Bphi_v, xyz_prev, xyz
real*8, allocatable :: iangle(:,:), l_part(:)

type(particle_kinetic_relativistic), allocatable :: prtkin(:)
real*8 :: rnd(1), psi, U, B(3), E(3)

type(octree_node) :: wall
type(octree_triangle), allocatable :: triangles(:)
real*8, allocatable                :: wtmp(:),nodes_xyz(:,:), normals_all(:,:)

! --- Read wall with octree
max_depth = 6
call mod_wall_collision_init('wall.h5',max_depth,wall)

! Read wall to initialize particles
file_id = 1
call HDF5_open('wall_to_load.h5',file_id)
call HDF5_integer_reading(file_id,n_tri,'ntriangle')

allocate(wtmp(9*n_tri))
call HDF5_array1D_reading(file_id, wtmp, 'nodes')

allocate(indices(n_tri,3))
call HDF5_array2D_reading_int(file_id, indices, 'indices')
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
! events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time))), & 
!      !event(diag,start=sim%time,step=1d-8),         &
!      event(stop_action(),start=sim%time+7.d-4)]

events = [event(read_jorek_fields_interp_linear(i=-1)), & 
event(stop_action(),start=sim%time+7.d-4)]

! Run first event to read the JOREK fields
call with(sim, events, at=0.d0)

! ------------- Initialize particles -----------------------------------------------------------------------------
sim%groups(:)%Z    = 1
sim%groups(:)%mass = 2.d0 !< atomic mass units

n_part = n_tri

allocate(particle_fieldline::sim%groups(1)%particles(n_tri))
allocate(normals_all(n_tri,3))

do i=1, n_tri
  
  select type (p=>sim%groups(1)%particles(i))
  
  type is (particle_fieldline)
  
   !--- Calculate wall normals (we assume they point towards the domain)
   v21 = triangles(i)%v1 - triangles(i)%v0 
   v31 = triangles(i)%v2 - triangles(i)%v0 

   norm_tria = [v21(2)*v31(3)- v21(3)*v31(2), v21(3)*v31(1)- v21(1)*v31(3), v21(1)*v31(2)- v21(2)*v31(1)]

   norm_tria        = norm_tria / norm2(norm_tria)
   normals_all(i,:) = norm_tria
   
   ! Assign particle position to triangle center
   xyz_tria = (triangles(i)%v0 + triangles(i)%v1 + triangles(i)%v2) / 3.d0

   !--- Move the point 1 mm away from the triangle
   !xyz_tria  = xyz_tria + norm_tria * 0.0001

   p%x = cartesian_to_cylindrical(xyz_tria)

   norm_R = (xyz_tria(1)*norm_tria(1) + xyz_tria(2)*norm_tria(2))/p%x(1)
   write(121,'(6ES16.6)') p%x(1), p%x(2), p%x(3), norm_R, norm_tria(3)

   call find_RZ(sim%fields%node_list, sim%fields%element_list, &
            p%x(1), p%x(2), &
            p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)
   
   if (ifail /= 0) write(*,*) 'Initial position of triangle id = ', i, ' not found in JOREK grid, move the wall inside'

   p%weight = 1.0
   p%v      = 1.d0 ! --- not really used

   end select
  
end do

allocate(iangle(1,n_part))
allocate(l_part(n_part))
iangle = 0
! ----------- end paticle initialization -----------------------------------------------------------------------

call check_and_fix_timesteps(timesteps, events)

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)

  do i=1,1

    select type (particles => sim%groups(i)%particles)
    type is (particle_fieldline)	
       l_part = 0.d0
      !$omp parallel do default(private) &
      !$omp shared (i, n_steps, timesteps, sim, wall, iangle, l_part, normals_all, tstep_field)
       do j=1,size(particles,1)

          do k=1,n_steps

             if(particles(j)%i_elm .le. 0) exit
             if (l_part(j) > 5.0) exit  !--- If the particle is able to move away X m from the wall, it is a wetted area

             pos_prev = particles(j)%x

             ! --- Do a step and check if the particle moves in the normal direction (away from the wall), otherwise correct direction
             if (k==1) then
               call field_line_runge_kutta_fixed_dt_push_jorek(sim%fields, particles(j), sim%time, 1.d-3)!*particles(j)%v)
               xyz_prev = cylindrical_to_cartesian(pos_prev)
               xyz      = cylindrical_to_cartesian(particles(j)%x)

               if (dot_product(xyz-xyz_prev,normals_all(j,:)) > 0.d0) then
                  dir_sign =  1.d0 
               else
                  dir_sign = -1.d0
               endif
             endif

             call field_line_runge_kutta_fixed_dt_push_jorek(sim%fields, particles(j), sim%time, tstep_field*dir_sign)

             xyz_prev = cylindrical_to_cartesian(pos_prev)
             xyz      = cylindrical_to_cartesian(particles(j)%x)

             l_part(j) = l_part(j) + norm2(xyz-xyz_prev)

             if (particles(j)%i_elm .le. 0) exit

             if (k < 3) cycle  ! don't check collisions during the initial steps to allow the particle to move away from the wall
             call mod_wall_collision_check(pos_prev, particles(j)%x, wall, wall_id, wall_pos, iangle(i,j))
             if(wall_id .gt. 0) then
                particles(j)%x      = wall_pos
                particles(j)%i_elm  = -wall_id
             end if
          end do
       end do

    end select

  enddo

  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo

call mod_wall_collision_free(wall)

call write_simulation_hdf5(sim, 'part_out.h5')

call mod_wall_collision_export(sim, 'wallload.h5', iangle)

! --- Write results to vtk
open(filehandle, file='wetted_area.vtk', status='replace', action='write')
140 format(a)
141 format(a,i8,a)
142 format(3es16.8)
143 format(a,2i8)
144 format(4i8)
write(filehandle,140) '# vtk DataFile Version 2.0'
write(filehandle,140) 'testdata'
write(filehandle,140) 'ASCII'
write(filehandle,140) 'DATASET POLYDATA'

n_nodes = maxval(indices) + 1
write(*,*) n_nodes
allocate(nodes_xyz(n_nodes,3))
do i = 1, n_tri
   nodes_xyz(indices(i,1)+1,:) = triangles(i)%v0
   nodes_xyz(indices(i,2)+1,:) = triangles(i)%v1
   nodes_xyz(indices(i,3)+1,:) = triangles(i)%v2
end do

! --- Triangle node positions
write(filehandle,141) 'POINTS', n_nodes, ' float'
do i = 1, n_nodes
  write(filehandle,142) nodes_xyz(i,:)*1000
end do

! --- Node indices corresponding to triangles
write(filehandle,143) 'POLYGONS', n_tri, n_tri * 4
do i = 1, n_tri
  write(filehandle,144) 3, indices(i,:) 
end do

! --- Write cell data
write(filehandle,141) 'CELL_DATA', n_tri
write(filehandle,140) 'SCALARS L_pre_collision float'
write(filehandle,140) 'LOOKUP_TABLE default'

do i = 1, n_tri
   if (l_part(i)>4.9) then
      write(filehandle,142) 1.d0
   else 
      write(filehandle,142) 0.d0
   endif
 end do

! --- Close file, clean up
close(filehandle)

deallocate(iangle)

! Finalize the simulation
call sim%finalize

end program example_wallload
