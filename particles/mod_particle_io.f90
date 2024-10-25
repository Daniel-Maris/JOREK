!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
implicit none
private
public write_simulation_hdf5,read_simulation_hdf5,get_simulation_hdf5_time
public write_simulation_hdf5_original,read_simulation_hdf5_original,get_simulation_hdf5_time_original

!> Module wide variables
integer,parameter :: master_task=0
integer,parameter :: n_cpu_1=1
integer,parameter :: group_name_len=12
contains

!> Export all particles using HDF5 IO
!> This communicates the particles internally before writing since this gave
!> better performance in our tests. A parallel version is in the git history.
!> Note that reading is still performed in parallel as this is not considered
!> time-critical and simpler to write.
subroutine write_simulation_hdf5_original(sim, filename)
use iso_c_binding, only: c_ptr
use mpi
use hdf5, only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F
use hdf5, only: h5open_f,h5fcreate_f,h5gcreate_f,h5gclose_f,h5fclose_f,h5close_f
use hdf5_io_module,     only: HDF5_real_saving,HDF5_integer_saving,HDF5_char_saving
use hdf5_io_module,     only: HDF5_array1D_saving,HDF5_array2D_saving,HDF5_array3D_saving
use hdf5_io_module,     only: HDF5_array1D_saving_r4,HDF5_array1D_saving_int
use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
use mod_particle_types, only: particle_gc,particle_gc_vpar
use mod_particle_types, only: particle_gc_Qin
use mod_particle_types, only: particle_fieldline
use mod_particle_types, only: particle_kinetic_relativistic
use mod_particle_types, only: particle_gc_relativistic
use mod_particle_sim,   only: particle_sim
integer(HSIZE_T), parameter :: particle_type_name_length = 40 !< length of the string used to identify a specific type of particle
character(len=*), parameter :: particle_type_name_field_name = 'particle_type' !< Name of the field containing the particle_type_name
type(particle_sim)   , intent(in) :: sim
character*(*)        , intent(in) :: filename

integer :: my_id, n_cpu, ierr
integer :: n_here, n_total
integer(HSIZE_T) :: i_here
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, create_file_space, write_file_space, dset, plist ! handles
integer(HID_T)                :: group_id
integer(HID_T)                :: data_type
integer(HID_T)                :: time_space_id, time_set_id
character(len=12)             :: group_name
character(len=particle_type_name_length) :: particle_type_name
integer                       :: i, j, hdferr, n_total_subarray
integer                       :: subarray2dn2type, doublesize
integer                       :: subarray2dn3type, subarray3d3x3type
integer                       :: resized2dn2type, resized2dn3type, resized3d3x3type
type(c_ptr) :: p_ptr
real*8,  dimension(:,:),   allocatable :: x, v, x_all, v_all, st, st_all
real*8,  dimension(:,:),   allocatable :: x_m, x_m_all, Astar_m, Astar_m_all
real*8,  dimension(:,:),   allocatable :: Astar_k, Astar_k_all, dBn_k, dBn_k_all 
real*8,  dimension(:,:),   allocatable :: Bnorm_k, Bnorm_k_all, E_k, E_k_all
real*8,  dimension(:,:,:), allocatable :: dAstar_k, dAstar_k_all
real*8,  dimension(:),     allocatable :: weight, weight_all, Vpar, E, mu, v1, B_norm
real*8,  dimension(:),     allocatable :: E_all, mu_all, v1_all, Vpar_all, B_norm_all
real*8,  dimension(:),     allocatable :: Vpar_m, Vpar_m_all, Bn_k, Bn_k_all
real*4,  dimension(:),     allocatable :: t_birth, t_birth_all
integer, dimension(:),     allocatable :: i_elm, i_elm_all, i_life, i_life_all
integer, dimension(:),     allocatable :: q, q_all, lost, lost_all

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))
particles_per_proc = 0


if (my_id .eq. 0) then
  ! Create file
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)
  if (hdferr .gt. 0) then
    write(*,*) "file open failed:", hdferr
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  endif

  ! Create group to write particle groups in
  call h5gcreate_f(file, "/groups", group_id, hdferr)
  call h5gclose_f(group_id, hdferr)

  ! Write the time
  call HDF5_real_saving(file,sim%time,'/time')
end if


if (allocated(sim%groups)) then
  do i=1,size(sim%groups,1)
    if (.not. allocated(sim%groups(i)%particles)) then
      write(*,*) "WARNING: group ", i, " not allocated, exiting"
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    ! Find the number of particles on each node
    n_here = size(sim%groups(i)%particles,1)
    call MPI_Gather(n_here,1,MPI_INTEGER,particles_per_proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_total = sum(particles_per_proc,1) ! set it on other processes as well (to
    ! 0) so they can allocate the useless arrays
    n_total_subarray = n_total
    call MPI_Bcast(n_total_subarray,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if (my_id .eq. 0) then
      ! Create group to write in
      write(group_name,"(A,i0.3,A)") "/groups/", i, "/"
      call h5gcreate_f(file, group_name, group_id, hdferr)
      call h5gclose_f(group_id, hdferr)
    end if

    ! create MPI types for transfering 2D arrays [3,n_particles] and 3d arrays [3,3,n_particles]
    call MPI_Type_size(MPI_DOUBLE_PRECISION,doublesize,ierr);
    ! 2D n=3 array type
    call MPI_Type_create_subarray(2,[2,n_total_subarray],[2,n_here],[0,0],&
    MPI_ORDER_FORTRAN,MPI_REAL8,subarray2dn2type,ierr)
    call MPI_Type_commit(subarray2dn2type,ierr);
    call MPI_Type_create_resized(subarray2dn2type,int(0,kind=MPI_Address_kind),&
    int(2*doublesize,kind=MPI_Address_kind),resized2dn2type,ierr)
    call MPI_Type_commit(resized2dn2type,ierr)
    ! 2D n=3 array type
    call MPI_Type_create_subarray(2,[3,n_total_subarray],[3,n_here],[0,0],&
    MPI_ORDER_FORTRAN,MPI_REAL8,subarray2dn3type,ierr)
    call MPI_Type_commit(subarray2dn3type,ierr);
    call MPI_Type_create_resized(subarray2dn3type,int(0,kind=MPI_Address_kind),&
    int(3*doublesize,kind=MPI_Address_kind),resized2dn3type,ierr)
    call MPI_Type_commit(resized2dn3type,ierr)
    ! 3D 3x3 array type
    call MPI_Type_create_subarray(3,[3,3,n_total_subarray],[3,3,n_here],[0,0,0],&
    MPI_ORDER_FORTRAN,MPI_REAL8,subarray3d3x3type,ierr) 
    call MPI_Type_commit(subarray3d3x3type,ierr);
    call MPI_Type_create_resized(subarray3d3x3type,int(0,kind=MPI_Address_kind),&
    int(9*doublesize,kind=MPI_Address_kind),resized3d3x3type,ierr)
    call MPI_type_commit(resized3d3x3type,ierr)  

    ! particle_base properties
    ! x
    allocate(x(3,n_here), x_all(3,n_total))
    do j=1,n_here
      x(:,j) = sim%groups(i)%particles(j)%x
    end do
!    call MPI_Gatherv(x(:,:), 3*n_here, MPI_REAL8, &
!      x_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_gatherv(x,3*n_here,MPI_REAL8,x_all,[(1,i=0,n_cpu-1)],&
    [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

    ! st
    allocate(st(2,n_here), st_all(2,n_total))
    do j=1,n_here
      st(:,j) = sim%groups(i)%particles(j)%st
    end do
!    call MPI_Gatherv(st(:,:), 2*n_here, MPI_REAL8, &
!      st_all(:,:), particles_per_proc*2, [(sum(particles_per_proc(0:i-1),1)*2, i=0,n_cpu-1)], &
!      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_gatherv(st,2*n_here,MPI_REAL8,st_all,[(1,i=0,n_cpu-1)],&
    [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn2type,0,MPI_COMM_WORLD,ierr)

    ! weight
    allocate(weight(n_here), weight_all(n_total))
    do j=1,n_here
      weight(j) = sim%groups(i)%particles(j)%weight
    end do
    call MPI_Gatherv(weight(:), n_here, MPI_REAL8, &
      weight_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! i_elm
    allocate(i_elm(n_here), i_elm_all(n_total))
    do j=1,n_here
      i_elm(j) = sim%groups(i)%particles(j)%i_elm
    end do
    call MPI_Gatherv(i_elm(:), n_here, MPI_INTEGER, &
      i_elm_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! i_life
    allocate(i_life(n_here), i_life_all(n_total))
    do j=1,n_here
      i_life(j) = sim%groups(i)%particles(j)%i_life
    end do
    call MPI_Gatherv(i_life(:), n_here, MPI_INTEGER, &
      i_life_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! t_birth
    allocate(t_birth(n_here), t_birth_all(n_total))
    do j=1,n_here
      t_birth(j) = sim%groups(i)%particles(j)%t_birth
    end do
    call MPI_Gatherv(t_birth(:), n_here, MPI_REAL4, &
      t_birth_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
      MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

    ! Write out stuff depending on particle type
    select type (p => sim%groups(i)%particles)

    type is (particle_kinetic)
      particle_type_name = 'particle_kinetic'

      ! v
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%v
      end do
!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_gatherv(v,3*n_here,MPI_REAL8,v_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

    ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_kinetic_leapfrog)
      particle_type_name = 'particle_kinetic_leapfrog'

      ! v
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%v
      end do
!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_gatherv(v,3*n_here,MPI_REAL8,v_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_gc)
      particle_type_name = 'particle_gc'

      ! E
      allocate(E(n_here), E_all(n_total))
      do j=1,n_here
        E(j) = p(j)%E
      end do
      call MPI_Gatherv(E(:), n_here, MPI_REAL8, &
        E_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,E_all,n_total,group_name//"E")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(E,mu,q,E_all,mu_all,q_all)

    type is (particle_gc_vpar)
      particle_type_name = 'particle_gc_vpar'

      ! Vpar
      allocate(Vpar(n_here), Vpar_all(n_total))
      do j=1,n_here
        Vpar(j) = p(j)%Vpar
      end do
      call MPI_Gatherv(Vpar(:), n_here, MPI_REAL8, &
        Vpar_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! B_norm
      allocate(B_norm(n_here), B_norm_all(n_total))
      do j=1,n_here
        B_norm(j) = p(j)%B_norm
      end do
      call MPI_Gatherv(B_norm(:), n_here, MPI_REAL8, &
        B_norm_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,Vpar_all,n_total,group_name//"Vpar")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving(file,B_norm_all,n_total,group_name//"B_norm")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(Vpar,mu,B_norm,q,Vpar_all,mu_all,B_norm_all,q_all)

    type is (particle_gc_Qin)
      particle_type_name = 'particle_gc_Qin'

      ! Vpar
      allocate(Vpar(n_here), Vpar_all(n_total))
      do j=1,n_here
        Vpar(j) = p(j)%Vpar
      end do
      call MPI_Gatherv(Vpar(:), n_here, MPI_REAL8, &
        Vpar_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! B_norm
      allocate(B_norm(n_here), B_norm_all(n_total))
      do j=1,n_here
        B_norm(j) = p(j)%B_norm
      end do
      call MPI_Gatherv(B_norm(:), n_here, MPI_REAL8, &
        B_norm_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ! x_m
      allocate(x_m(3,n_here), x_m_all(3,n_total))
      do j=1,n_here
        x_m(:,j) = p(j)%x_m
      enddo      
      call MPI_gatherv(x_m,3*n_here,MPI_REAL8,x_m_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr) 

      ! Vpar_m
      allocate(Vpar_m(n_here), Vpar_m_all(n_total))
      do j=1,n_here
        Vpar_m(j) = p(j)%Vpar_m
      end do
      call MPI_Gatherv(Vpar_m(:), n_here, MPI_REAL8, &
        Vpar_m_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! Astar_m
      allocate(Astar_m(3,n_here), Astar_m_all(3,n_total))
      do j=1,n_here
        Astar_m(:,j) = p(j)%Astar_m
      enddo      
      call MPI_gatherv(Astar_m,3*n_here,MPI_REAL8,Astar_m_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      ! Astar_k
      allocate(Astar_k(3,n_here), Astar_k_all(3,n_total))
      do j=1,n_here
        Astar_k(:,j) = p(j)%Astar_k
      enddo      
      call MPI_gatherv(Astar_k,3*n_here,MPI_REAL8,Astar_k_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      ! dAstar_k
      allocate(dAstar_k(3,3,n_here), dAstar_k_all(3,3,n_total))
      do j=1,n_here
        dAstar_k(:,:,j) = p(j)%dAstar_k
      enddo
      call MPI_gatherv(dAstar_k,9*n_here,MPI_REAL8,dAstar_k_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized3d3x3type,0,MPI_COMM_WORLD,ierr)

      ! Bn_k
      allocate(Bn_k(n_here), Bn_k_all(n_total))
      do j=1,n_here
        Bn_k(j) = p(j)%Bn_k
      end do
      call MPI_Gatherv(Bn_k(:), n_here, MPI_REAL8, &
        Bn_k_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! dBn_k
      allocate(dBn_k(3,n_here), dBn_k_all(3,n_total))
      do j=1,n_here
        dBn_k(:,j) = p(j)%dBn_k
      enddo      
      call MPI_gatherv(dBn_k,3*n_here,MPI_REAL8,dBn_k_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      ! Bnorm_k
      allocate(Bnorm_k(3,n_here), Bnorm_k_all(3,n_total))
      do j=1,n_here
        Bnorm_k(:,j) = p(j)%Bnorm_k
      enddo      
      call MPI_gatherv(Bnorm_k,3*n_here,MPI_REAL8,Bnorm_k_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      ! E_k
      allocate(E_k(3,n_here), E_k_all(3,n_total))
      do j=1,n_here
        E_k(:,j) = p(j)%E_k
      enddo      
      call MPI_gatherv(E_k,3*n_here,MPI_REAL8,E_k_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,Vpar_all,n_total,group_name//"Vpar")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving(file,B_norm_all,n_total,group_name//"B_norm")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
        call HDF5_array2D_saving(file,x_m_all,3,n_total,group_name//"x_m")
        call HDF5_array1D_saving(file,Vpar_m_all,n_total,group_name//"Vpar_m")
        call HDF5_array2D_saving(file,Astar_m_all,3,n_total,group_name//"Astar_m")
        call HDF5_array2D_saving(file,Astar_k_all,3,n_total,group_name//"Astar_k")
        call HDF5_array3D_saving(file,dAstar_k_all,3,3,n_total,group_name//"dAstar_k")
        call HDF5_array1D_saving(file,Bn_k_all,n_total,group_name//"Bn_k")
        call HDF5_array2D_saving(file,dBn_k_all,3,n_total,group_name//"dBn_k")
        call HDF5_array2D_saving(file,Bnorm_k_all,3,n_total,group_name//"Bnorm_k")
        call HDF5_array2D_saving(file,E_k_all,3,n_total,group_name//"E_k")
        
      end if
      deallocate(Vpar,mu,B_norm,q,Vpar_all,mu_all,B_norm_all,q_all,x_m,x_m_all,&
      Vpar_m,Vpar_m_all,Astar_m,Astar_m_all,Astar_k,Astar_k_all,dAstar_k,dAstar_k_all,&
      dBn_k,dBn_k_all,Bnorm_k,Bnorm_k_all,E_k,E_k_all)

    type is (particle_fieldline)
      particle_type_name = 'particle_fieldline'
      ! v
      allocate(v1(n_here), v1_all(n_total))
      do j=1,n_here
        v1(j) = p(j)%v
      end do
      call MPI_Gatherv(v1(:), n_here, MPI_REAL8, &
        v1_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,v1_all,n_total,group_name//"v")
      end if
      deallocate(v1, v1_all)

    type is (particle_kinetic_relativistic)
      particle_type_name = 'particle_kinetic_relativistic'
 
      ! momenta
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%p
      end do

!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(0:i-1),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
       call MPI_gatherv(v,3*n_here,MPI_REAL8,v_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn3type,0,MPI_COMM_WORLD,ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v") 
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q") 
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_gc_relativistic)
      particle_type_name = 'particle_gc_relativistic'

      ! momenta (parallel momentum and magnetic moment)
      allocate(v(2,n_here), v_all(2,n_total))
      do j=1,n_here
        v(:,j) = p(j)%p
      end do

!      call MPI_Gatherv(v(:,:), 2*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*2, [(sum(particles_per_proc(0:i-1),1)*2, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_gatherv(v,2*n_here,MPI_REAL8,v_all,[(1,i=0,n_cpu-1)],&
      [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)],resized2dn2type,0,MPI_COMM_WORLD,ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,2,n_total,group_name//"v") 
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q") 
      end if
      deallocate(v,q,v_all,q_all)

    class default
      write(*,*) "error: missing type name declaration for write"
      call exit(1)
    end select


    ! It is important to do the gathering first, because that is the collective part
    if (my_id .eq. 0) then
      call HDF5_array2D_saving(file,x_all,3,n_total,group_name//"x")
      call HDF5_array2D_saving(file,st_all,2,n_total,group_name//"st")
      call HDF5_array1D_saving(file,weight_all,n_total,group_name//"weight")
      call HDF5_array1D_saving_int(file,i_elm_all,n_total,group_name//"i_elm")
      call HDF5_array1D_saving_int(file,i_life_all,n_total,group_name//"i_life")
      call HDF5_array1D_saving_r4(file,t_birth_all,n_total,group_name//"t_birth")

      call HDF5_char_saving(file,particle_type_name,group_name//"type")
      call HDF5_integer_saving(file,sim%groups(i)%Z,group_name//"Z")
      call HDF5_real_saving(file,sim%groups(i)%mass,group_name//"mass")

      call HDF5_char_saving(file,sim%groups(i)%ad%suffix,group_name//"adas_suffix")
    end if
    deallocate(x,x_all,st,st_all,weight,weight_all,t_birth,t_birth_all,i_elm,i_elm_all,i_life,i_life_all)
    call MPI_Type_free(subarray2dn2type,ierr);  call MPI_Type_free(resized2dn2type,ierr);
    call MPI_Type_free(subarray2dn3type,ierr);  call MPI_Type_free(resized2dn3type,ierr);
    call MPI_Type_free(subarray3d3x3type,ierr); call MPI_Type_free(resized3d3x3type,ierr);
  end do
end if

! Close everything
if (my_id .eq. 0) then
  call h5fclose_f(file, hdferr)
  call h5close_f(hdferr)

  write(*,*) "Writing particle output file to ", filename, " succeeded"
end if

end subroutine write_simulation_hdf5_original

!> Export all particles using HDF5 IO
!> Parallel support to writing operations is provided
!> by MPI enabled HDF5 procedures.
!> inputs:
!>   sim:                        (particle_sim) particle simulation object
!>   filename:                   (character)(N) name of the output file
!>   file_access_in:             (integer)(optional) define how the HDF5 should behave:
!>                               when opening or creating a new file.
!>                               default: H5F_ACC_TRUNC_F (truncate fille if already exists
!>                               create otherwise)
!>   use_hdf5_access_properties: (logical)(optional) HDF5 file access property
!>                               if must be set to .false. for parallel I/O
!>                               default: .true. 
!>   collective_mpio_in:         (logical)(optional) define whether MPIO collective calls are
!>                               performed default: true
!>   mpi_comm_in:                (integer)(optional) MPI communicator identifier
!>   mpi_info_in:                (integer)(optional) MPI info structre for parallel IO
subroutine write_simulation_hdf5(sim,filename,file_access_in,&
use_hdf5_access_properties,collective_mpio_in,mpi_comm_in,mpi_info_in)
  use mpi
  use hdf5,               only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F
  use hdf5,               only: H5Gcreate_f,H5Gclose_f
  use hdf5_io_module,     only: HDF5_open_or_create,HDF5_close
  use hdf5_io_module,     only: HDF5_integer_saving 
  use hdf5_io_module,     only: HDF5_real_saving,HDF5_char_saving
  use hdf5_io_module,     only: HDF5_array1D_saving_int
  use hdf5_io_module,     only: HDF5_array1D_saving_r4
  use hdf5_io_module,     only: HDF5_array1D_saving
  use hdf5_io_module,     only: HDF5_array2D_saving
  use hdf5_io_module,     only: HDF5_array3D_saving
  use mod_particle_types, only: particle_arrays_from_list
  use mod_particle_types, only: deallocate_particle_arrays
  use mod_particle_sim,   only: particle_sim
  implicit none
  !> parameters
  integer,parameter             :: master_rank=0
  integer(HSIZE_T),parameter    :: i0_HSIZE_T=int(0,kind=HSIZE_T)
  !> input variables
  type(particle_sim),intent(in) :: sim
  character(len=*),  intent(in) :: filename 
  integer, intent(in), optional :: file_access_in,mpi_comm_in,mpi_info_in
  logical, intent(in), optional :: use_hdf5_access_properties,collective_mpio_in
  !> variables
  integer                       :: file_access_loc
  integer                       :: mpi_comm_loc,mpi_info_loc 
  integer                       :: ii,ierr,h5err,n_groups,n_particles
  integer                       :: n_particles_per_group
  integer(HID_T)                :: file_id,group_id
  integer(HSIZE_T)              :: n_particles_offset
  integer,  dimension(:),    allocatable :: n_particles_loc
  integer*4,dimension(:),    allocatable :: i_elm_arr,i_life_arr
  integer*4,dimension(:),    allocatable :: q_arr
  integer,  dimension(:,:),  allocatable :: n_particles_glob
  real*4,   dimension(:),    allocatable :: t_birth_arr
  real*8,   dimension(:),    allocatable :: weight_arr,v_1d_arr
  real*8,   dimension(:),    allocatable :: E_arr,mu_arr,vpar_arr
  real*8,   dimension(:),    allocatable :: B_norm_arr,vpar_m_arr,Bn_k_arr
  real*8,   dimension(:,:),  allocatable :: st_arr,x_arr,B_hat_prev_arr,v_2d_arr
  real*8,   dimension(:,:),  allocatable :: x_m_arr,Astar_m_arr,Astar_k_arr
  real*8,   dimension(:,:),  allocatable :: dBn_k_arr,Bnorm_k_arr,E_k_arr
  real*8,   dimension(:,:,:),allocatable :: dAstar_k_arr
  logical                                :: create_access_plist,collective_mpio_loc
  character(len=group_name_len)          :: group_name
  character(len=:),          allocatable :: particle_type_str

  !> preparation
  file_access_loc = H5F_ACC_TRUNC_F !< truncate the file by default
  if(present(file_access_in)) file_access_loc = file_access_in;
  create_access_plist = .false. !< serial access by default
  if(present(use_hdf5_access_properties)) create_access_plist = .not.use_hdf5_access_properties
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_info_loc = mpi_info_in
  collective_mpio_loc = .true. !< enable collective MPIO applications by default
  if(present(collective_mpio_in)) collective_mpio_loc = collective_mpio_in
  !> create the hdf5 file and the groups fields
  call HDF5_open_or_create(filename,file_id,h5err,&
  file_access=file_access_loc,create_access_plist_in=create_access_plist,& 
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  if(h5err.gt.0) then
    if(sim%my_id.eq.master_task) write(*,*) "Failed to create or open the ",&
    filename," file: ",h5err,", ABORT!"
    call MPI_Abort(mpi_comm_loc,-1,ierr)
  endif
  call H5Gcreate_f(file_id,"/groups",group_id,h5err) !< create particle groups
  call H5Gclose_f(group_id,h5err)
  !> write the time in HDF5 file, we assume that each MPI task reached the same physical time
  call HDF5_real_saving(file_id,sim%time,"/time") !< TODO find a way that only master write the data standard if condition deadlocks
  !> check if loops are allocated and loop on them
  if(allocated(sim%groups)) then
    !> it is assumed that all processors has the same number of groups but
    !> different size of particle lists per each group so we gather the
    !> particle list size for all groups firstly
    n_groups = size(sim%groups,1); allocate(n_particles_loc(n_groups)); n_particles_loc=0;
    do ii=1,n_groups
      n_particles_loc(ii) = size(sim%groups(ii)%particles,1)
    enddo
    allocate(n_particles_glob(sim%n_cpu*n_groups,1)); n_particles_glob=0;
    call MPI_Allgather(n_particles_loc,n_groups,MPI_INTEGER,n_particles_glob(:,1),n_groups,&
    MPI_INTEGER,mpi_comm_loc,ierr)
    n_particles_glob = transpose(reshape(n_particles_glob,(/n_groups,sim%n_cpu/)))
    !> write particle list data into HDF5 file
    do ii=1,n_groups
      if(.not.allocated(sim%groups(ii)%particles)) then
        if(sim%my_id.eq.master_task) write(*,*) "WARNING: particle list N# ",&
        ii,"  not allocated: skip!"; cycle;
      endif
      !> number of total particles for the group
      n_particles_per_group = sum(n_particles_glob(:,ii))
      n_particles_offset    = int(sum(n_particles_glob(1:sim%my_id,ii)),kind=HSIZE_T)
      !> create the HDF5 group for the particle list
      write(group_name,"(A,i0.3,A)") "/groups/",ii,"/"
      call H5Gcreate_f(file_id,trim(group_name),group_id,h5err)
      call H5Gclose_f(group_id,h5err)
      !> reorganize and store the particle data in congruent arrays
      call particle_arrays_from_list(sim%groups(ii)%particles,n_particles,&
      i_elm_arr,i_life_arr,q_arr,t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,&
      vpar_arr,B_norm_arr,vpar_m_arr,st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,&
      Astar_m_arr,Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr,&
      particle_type_str)
      !> write data in HDF5 file
      !> TODO for allowing parallel applications on implementation having
      !> having only serial HDF5 installed create a wrapped for each of 
      !> this function in which one can choose to use the parallel implementation
      !> of HDF5 or the previously implemented MPI gather + serial HDF5 writing
      if(allocated(i_elm_arr)) call HDF5_array1D_saving_int(file_id,i_elm_arr,&
      n_particles_per_group,trim(group_name)//"i_elm",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(i_life_arr)) call HDF5_array1D_saving_int(file_id,i_life_arr,&
      n_particles_per_group,trim(group_name)//"i_life",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(q_arr)) call HDF5_array1D_saving_int(file_id,q_arr,&
      n_particles_per_group,trim(group_name)//"q",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(t_birth_arr)) call HDF5_array1D_saving_r4(file_id,t_birth_arr,&
      n_particles_per_group,trim(group_name)//"t_birth",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(weight_arr)) call HDF5_array1D_saving(file_id,weight_arr,&
      n_particles_per_group,trim(group_name)//"weight",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(v_1d_arr)) call HDF5_array1D_saving(file_id,v_1d_arr,&
      n_particles_per_group,trim(group_name)//"v",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)     
      if(allocated(E_arr)) call HDF5_array1D_saving(file_id,E_arr,&
      n_particles_per_group,trim(group_name)//"E",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(mu_arr)) call HDF5_array1D_saving(file_id,mu_arr,&
      n_particles_per_group,trim(group_name)//"mu",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(vpar_arr)) call HDF5_array1D_saving(file_id,vpar_arr,&
      n_particles_per_group,trim(group_name)//"Vpar",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(B_norm_arr)) call HDF5_array1D_saving(file_id,B_norm_arr,&
      n_particles_per_group,trim(group_name)//"B_norm",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(vpar_m_arr)) call HDF5_array1D_saving(file_id,vpar_m_arr,&
      n_particles_per_group,trim(group_name)//"Vpar_m",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(st_arr)) call HDF5_array2D_saving(file_id,st_arr,size(st_arr,1),&
      n_particles_per_group,trim(group_name)//"st",start=[i0_HSIZE_T,&
      n_particles_offset],mpio_collective_in=collective_mpio_loc)
      if(allocated(x_arr)) call HDF5_array2D_saving(file_id,x_arr,size(x_arr,1),&
      n_particles_per_group,trim(group_name)//"x",start=[i0_HSIZE_T,&
      n_particles_offset],mpio_collective_in=collective_mpio_loc)
      if(allocated(B_hat_prev_arr)) call HDF5_array2D_saving(&
      file_id,B_hat_prev_arr,size(B_hat_prev_arr,1),n_particles_per_group,&
      trim(group_name)//"B_hat_prev",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(v_2d_arr)) call HDF5_array2D_saving(&
      file_id,v_2d_arr,size(v_2d_arr,1),n_particles_per_group,&
      trim(group_name)//"v",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(x_m_arr)) call HDF5_array2D_saving(&
      file_id,x_m_arr,size(x_m_arr,1),n_particles_per_group,&
      trim(group_name)//"x_m",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(Astar_m_arr)) call HDF5_array2D_saving(&
      file_id,Astar_m_arr,size(Astar_m_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_m",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(Astar_k_arr)) call HDF5_array2D_saving(&
      file_id,Astar_k_arr,size(Astar_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Astar_k",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(Bn_k_arr)) call HDF5_array1D_saving(file_id,Bn_k_arr,&
      n_particles_per_group,trim(group_name)//"Bn_k",start=[n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(dBn_k_arr)) call HDF5_array2D_saving(file_id,&
      dBn_k_arr,size(dBn_k_arr,1),n_particles_per_group,&
      trim(group_name)//"dBn_k",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(Bnorm_k_arr)) call HDF5_array2D_saving(file_id,&
      Bnorm_k_arr,size(Bnorm_k_arr,1),n_particles_per_group,&
      trim(group_name)//"Bnorm_k",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(E_k_arr)) call HDF5_array2D_saving(&
      file_id,E_k_arr,size(E_k_arr,1),n_particles_per_group,&
      trim(group_name)//"E_k",start=[i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      if(allocated(dAstar_k_arr)) call HDF5_array3D_saving(file_id,&
      dAstar_k_arr,size(dAstar_k_arr,1),size(dAstar_k_arr,2),n_particles_per_group,&
      trim(group_name)//"dAstar_k",start=[i0_HSIZE_T,i0_HSIZE_T,n_particles_offset],&
      mpio_collective_in=collective_mpio_loc)
      !> Write particle group attributes in HDF5 file, we assume that the attributes
      !> of the same group index for all tasks are equals. Therefore, we write the
      !> group attributes of only the master task
      !> TODO find a way that only master write the data standard if condition deadlocks
      if(allocated(particle_type_str)) call HDF5_char_saving(file_id,&
      particle_type_str,trim(group_name)//"type")
      call HDF5_char_saving(file_id,sim%groups(ii)%ad%suffix,trim(group_name)//"adas_suffix")
      call HDF5_integer_saving(file_id,sim%groups(ii)%Z,trim(group_name)//"Z")
      call HDF5_real_saving(file_id,sim%groups(ii)%mass,trim(group_name)//"mass")
      !> deallocate structures
      call deallocate_particle_arrays(n_particles,i_elm_arr,i_life_arr,q_arr,&
      t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,vpar_arr,B_norm_arr,vpar_m_arr,&
      st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,Astar_m_arr,Astar_k_arr,&
      Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr)
    enddo
  else
    if(sim%my_id.eq.master_task) write(*,*) "WARNING: sim particle groups is not allocated!"
  endif
  !> cleanups
  call HDF5_close(file_id)
  if(allocated(n_particles_loc))   deallocate(n_particles_loc)
  if(allocated(n_particles_glob))  deallocate(n_particles_glob)
  if(allocated(particle_type_str)) deallocate(particle_type_str)
end subroutine write_simulation_hdf5

!> Original procedure for importing all particles.
!> Reads the number of particles from a file, determines a
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine read_simulation_hdf5_original(sim, filename, test_in)
use iso_c_binding, only: c_ptr
use mpi
use hdf5, only: HSIZE_T,HID_T,H5F_ACC_TRUNC_F,H5P_FILE_ACCESS_F,H5F_ACC_RDONLY_F
use hdf5, only: h5open_f,h5pcreate_f,h5pset_fapl_mpio_f,h5fopen_f,h5pclose_f
use hdf5, only: h5gopen_f,h5gget_info_f,h5gclose_f,h5dopen_f,h5dget_space_f
use hdf5, only: h5sget_simple_extent_ndims_f,h5sget_simple_extent_dims_f
use hdf5, only: h5sclose_f,h5dclose_f,h5lexists_f,h5fclose_f,h5close_f
use hdf5_io_module, only: HDF5_real_reading,HDF5_allocatable_char_reading,HDF5_integer_reading
use hdf5_io_module, only: HDF5_array1D_reading,HDF5_array2D_reading,HDF5_array3D_reading
use hdf5_io_module, only: HDF5_array1D_reading_int,HDF5_array1D_reading_r4
use hdf5_io_module, only: HDF5_allocatable_array1D_reading
use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
use mod_particle_types, only: particle_gc,particle_gc_vpar
use mod_particle_types, only: particle_gc_Qin
use mod_particle_types, only: particle_fieldline
use mod_particle_types, only: particle_kinetic_relativistic
use mod_particle_types, only: particle_gc_relativistic
use mod_openadas,   only: read_adf11
use mod_coronal
use mod_particle_sim, only: particle_sim
type(particle_sim) , intent(inout)         :: sim
character*(*)      , intent(in)            :: filename
logical            , intent(in) , optional :: test_in

integer                            :: my_id, n_cpu, ierr, rank, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)               :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)               :: data_type
integer(HID_T)               :: group_id
integer(HID_T)               :: time_set_id
integer(HSIZE_T)             :: i_here
integer                      :: n_here
integer                      :: storage_type, max_corder
character(len=12)            :: group_name
character(len=:),allocatable :: particle_type_name,adas_suffix
integer                      :: i, j, n, hdferr, n_alive
integer, allocatable         :: n_alive_all(:)
logical                      :: exists,test

type(c_ptr) :: p_ptr
integer*8, dimension(1:2) :: tmp, maxdims
real*8,    dimension(:,:,:), allocatable :: real8_3D
real*8,    dimension(:,:),   allocatable :: real8_2D
integer*4, dimension(:),     allocatable :: int4_1D
real*4,    dimension(:),     allocatable :: real4_1D
real*8,    dimension(:),     allocatable :: real8_1D

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))
test = .false.; if(present(test_in)) test = test_in;

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! read the time
call HDF5_allocatable_array1D_reading(file,real8_1D,'/time')
sim%time = real8_1D(1); if(allocated(real8_1D)) deallocate(real8_1D);

! Get the number of groups
call h5gopen_f(file, '/groups/', group_id, hdferr)
call h5gget_info_f(group_id, storage_type, n, max_corder, hdferr)
call h5gclose_f(group_id, hdferr)

! Reallocate groups if necessary
if (allocated(sim%groups)) deallocate(sim%groups)
allocate(sim%groups(n))
do i=1,n
  ! Open the dataset for x
  write(group_name,'(A,i0.3,A)') '/groups/', i, '/'
  call h5dopen_f(file, group_name//"x", dset, hdferr)

  ! Open the file dataspace
  call h5dget_space_f(dset, file_space, hdferr)

  ! Get the number of particles
  call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
  call h5sget_simple_extent_dims_f(file_space, tmp, maxdims, hdferr)
  n_particles = int(tmp(2),4)
  call h5sclose_f(file_space, hdferr)
  call h5dclose_f(dset, hdferr)

  ! Divide particles over processors
  particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
  particles_per_proc(1:n_cpu-1) = n_particles/n_cpu
  i_here = sum(particles_per_proc(0:my_id-1))
  n_here = particles_per_proc(my_id)

  ! Get the particle type from the attribute
  call HDF5_allocatable_char_reading(file,particle_type_name,group_name//"type",&
  mpi_rank=0,n_mpi_tasks=1)

  ierr = 0
  select case (trim(particle_type_name))
  case ('particle_kinetic')
    allocate(particle_kinetic::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_kinetic_leapfrog')
    allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc')
    allocate(particle_gc::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc_vpar')
    allocate(particle_gc_vpar::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_fieldline')
    allocate(particle_fieldline::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_kinetic_relativistic')
    allocate(particle_kinetic_relativistic::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc_relativistic')
    allocate(particle_gc_relativistic::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc_Qin')
    allocate(particle_gc_Qin::sim%groups(i)%particles(n_here), stat=ierr)
  case default
    write(*,*) "error: missing type name declaration ", trim(particle_type_name), " for read"
    call exit(1)
  end select

  if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
      "unable to allocate particles(", particles_per_proc(my_id), ")"
  call HDF5_integer_reading(file,sim%groups(i)%Z,group_name//"Z")
  call HDF5_real_reading(file,sim%groups(i)%mass,group_name//"mass")
  call HDF5_allocatable_char_reading(file,adas_suffix,group_name//"adas_suffix",&
  mpi_rank=0,n_mpi_tasks=1)
  sim%groups(i)%ad%suffix = adas_suffix; 
  if ((len_trim(sim%groups(i)%ad%suffix) .gt. 0) .and. (.not.test)) then
    sim%groups(i)%ad = read_adf11(my_id,sim%groups(i)%ad%suffix)
    sim%groups(i)%cor = coronal(sim%groups(i)%ad)
  end if

  ! Read base particle attributes
  ! x
  allocate(real8_2D(3,n_here))
  call HDF5_array2D_reading(file, real8_2D, group_name//"x",start=[0_HSIZE_T,i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%x = real8_2D(:,j)
  end do
  deallocate(real8_2D)

  ! st
  allocate(real8_2D(2,n_here))
  call HDF5_array2D_reading(file, real8_2D, group_name//"st",start=[0_HSIZE_T,i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%st = real8_2D(:,j)
  end do
  deallocate(real8_2D)

  ! weight
  allocate(real8_1D(n_here))
  call HDF5_array1D_reading(file, real8_1D, group_name//"weight",start=[i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%weight = real8_1D(j)
  end do
  deallocate(real8_1D)

 ! i_elm
  allocate(int4_1D(n_here))
  call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_elm", start=[i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%i_elm = int4_1D(j)
  end do
  deallocate(int4_1D)

  ! The following two are relatively new, so might not always be present
  ! i_life
  call h5lexists_f(file, group_name//"i_life", exists, ierr)
  if (exists) then
    allocate(int4_1D(n_here))
    int4_1D = 0 ! preset to 0 in case not present (i.e. old restart files)
    ! will still give a nasty error message, but should work,
    ! though who knows what hdf5 does with our array if reading fails...
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_life", start=[i_here])
    do j=1,n_here
      sim%groups(i)%particles(j)%i_life = int4_1D(j)
    end do
    deallocate(int4_1D)
  end if

  ! t_birth
  call h5lexists_f(file, group_name//"t_birth", exists, ierr)
  if (exists) then
    allocate(real4_1D(n_here))
    real4_1D = 0.0 ! preset to 0 in case not present (i.e. old restart files)
    call HDF5_array1D_reading_r4(file, real4_1D, group_name//"t_birth",start=[i_here])
    do j=1,n_here
      sim%groups(i)%particles(j)%t_birth = real4_1D(j)
    end do
    deallocate(real4_1D)
  end if

  select type (p => sim%groups(i)%particles)
  type is (particle_kinetic)
    ! v
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%v = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)

  type is (particle_kinetic_leapfrog)
    ! v
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%v = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)

  type is (particle_gc)
    ! E
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"E",start=[i_here])
    do j=1,n_here
      p(j)%E = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! mu
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
    do j=1,n_here
      p(j)%mu = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)

  type is (particle_gc_vpar)
    ! vpar 
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar",start=[i_here])
    do j=1,n_here
      p(j)%Vpar = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! mu
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
    do j=1,n_here
      p(j)%mu = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! Bnorm
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"B_norm",start=[i_here])
    do j=1,n_here
      p(j)%B_norm = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)

  type is (particle_gc_Qin)
    ! vpar 
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar",start=[i_here])
    do j=1,n_here
      p(j)%vpar = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! mu
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
    do j=1,n_here
      p(j)%mu = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! Bnorm
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"B_norm",start=[i_here])
    do j=1,n_here
      p(j)%B_norm = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)
    ! x_m
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"x_m",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%x_m = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! Vpar_m
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar_m",start=[i_here])
    do j=1,n_here
      p(j)%vpar_m = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! Astar_m
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"Astar_m",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%Astar_m = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! Astar_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"Astar_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%Astar_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! dAstar_k
    allocate(real8_3D(3,3,n_here))
    call HDF5_array3D_reading(file, real8_3D, group_name//"dAstar_k",start=[0_HSIZE_T,0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%dAstar_k = real8_3D(:,:,j)
    end do
    deallocate(real8_3D)
    ! Bn_k 
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Bn_k",start=[i_here])
    do j=1,n_here
      p(j)%Bn_k = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! dBn_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"dBn_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%dBn_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! Bnorm_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"Bnorm_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%Bnorm_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
    ! E_k
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"E_k",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%E_k = real8_2D(:,j)
    end do
    deallocate(real8_2D)
   
  type is (particle_fieldline)
    ! v
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"v",start=[i_here])
    do j=1,n_here
      p(j)%v = real8_1D(j)
    end do
    deallocate(real8_1D)

  type is (particle_kinetic_relativistic)

    ! momenta [AMU*m/s]
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%p = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! electric charge q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)    

  type is (particle_gc_relativistic)
    ! momenta (parallel momentum and magnetic moment)
    allocate(real8_2D(2,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%p = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! electric charge q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)    

  end select
  ! Check if the balance between processors is okay, by comparing the number of
  ! alive particles
  n_alive = count(sim%groups(i)%particles(:)%i_elm .gt. 0)
  allocate(n_alive_all(sim%n_cpu))
  call MPI_Gather(n_alive, 1, MPI_INTEGER, n_alive_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  ! Check if the imbalance is not too great
  if (sim%my_id .eq. 0) then
    if (maxval(n_alive_all) .gt. minval(n_alive_all) * 1.5) then
      write(*,*) "WARNING: ", (maxval(n_alive_all)*100)/minval(n_alive_all), '% imbalance between CPUs, counts:'
      write(*,*) n_alive_all
    end if
  end if
  deallocate(n_alive_all)
  if(allocated(particle_type_name)) deallocate(particle_type_name)
  if(allocated(adas_suffix))        deallocate(adas_suffix);
end do

! Close everything else
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine read_simulation_hdf5_original


!> Parallel read of particle HDF5 restart file
!> inputs:
!>   filename:       (character)(N) name of the output file
!>   sim:            (particle_sim) particle simulation object
!>   use_hdf5_access_properties: (logical)(optional) HDF5 file access property
!>                   if must be set to .false. for parallel I/O
!>                   default: .true. 
!>   mpi_comm_in:    (integer)(optional) MPI communicator identifier
!>   mpi_info_in:    (integer)(optional) MPI info structre for parallel IO
!> outputs:
!>   sim: (particle_sim) particle simulation object
subroutine read_simulation_hdf5(sim,filename,use_hdf5_access_properties,&
mpi_comm_in,mpi_info_in,test_in)
  use mpi
  use hdf5,               only: HSIZE_T,HID_T
  use hdf5,               only: H5Gopen_f,H5Gget_info_f
  use hdf5,               only: H5Gclose_f
  use hdf5_io_module,     only: HDF5_open,HDF5_close
  use hdf5_io_module,     only: HDF5_char_reading
  use hdf5_io_module,     only: HDF5_allocatable_char_reading
  use hdf5_io_module,     only: HDF5_real_reading,HDF5_integer_reading
  use hdf5_io_module,     only: HDF5_get_dataset_rank_dims
  use hdf5_io_module,     only: HDF5_allocatable_array1D_reading_int
  use hdf5_io_module,     only: HDF5_allocatable_array1D_reading_r4
  use hdf5_io_module,     only: HDF5_allocatable_array1D_reading
  use hdf5_io_module,     only: HDF5_allocatable_array2D_reading
  use hdf5_io_module,     only: HDF5_allocatable_array3D_reading
  use mod_particle_types, only: particle_list_from_arrays
  use mod_particle_types, only: initialize_particle_list_to_zero
  use mod_particle_types, only: deallocate_particle_arrays
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_gc_vpar
  use mod_particle_types, only: particle_gc_Qin
  use mod_particle_types, only: particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_sim,   only: particle_sim
  use mod_coronal,        only: coronal
  use mod_openadas,       only: read_adf11
  implicit none
  !> parameters:
  integer(HSIZE_T),parameter  :: i0_HSIZE_T=int(0,kind=HSIZE_T)
  integer(HSIZE_T),parameter  :: n1_HSIZE_T=int(-1,kind=HSIZE_T)
  !> inputs:
  character(len=*),intent(in) :: filename
  logical,intent(in),optional :: use_hdf5_access_properties
  integer,intent(in),optional :: mpi_comm_in,mpi_info_in
  logical,intent(in),optional :: test_in
  !> inputs-outputs:
  class(particle_sim),intent(inout) :: sim  
  !> variables:
  integer                                        :: ii,ierr,h5err,errorcode,n_groups
  integer                                        :: mpi_comm_loc,mpi_info_loc
  integer                                        :: storage_type,max_corder,rank
  integer(HID_T)                                 :: file_id,group_id
  integer(HSIZE_T)                               :: offset,n_particles_hsizet
  integer,          dimension(:),    allocatable :: n_particles_per_proc
  integer*4,        dimension(:),    allocatable :: i_elm_arr,i_life_arr
  integer*4,        dimension(:),    allocatable :: q_arr
  integer(HSIZE_T), dimension(:),    allocatable :: n_particles_tot,n_particles_max
  real*4,           dimension(:),    allocatable :: t_birth_arr
  real*8,           dimension(:),    allocatable :: weight_arr,v_1d_arr
  real*8,           dimension(:),    allocatable :: E_arr,mu_arr,vpar_arr
  real*8,           dimension(:),    allocatable :: B_norm_arr,vpar_m_arr,Bn_k_arr
  real*8,           dimension(:,:),  allocatable :: st_arr,x_arr,B_hat_prev_arr,v_2d_arr
  real*8,           dimension(:,:),  allocatable :: x_m_arr,Astar_m_arr,Astar_k_arr
  real*8,           dimension(:,:),  allocatable :: dBn_k_arr,Bnorm_k_arr,E_k_arr
  real*8,           dimension(:,:,:),allocatable :: dAstar_k_arr
  character(len=group_name_len)                  :: group_name
  character(len=:),                  allocatable :: particle_type_str
  logical                                        :: create_access_plist,test
  !> initialisation
  allocate(n_particles_per_proc(sim%n_cpu))
  !> set optional parameters
  create_access_plist = .false.
  if(present(use_hdf5_access_properties)) create_access_plist = .not.use_hdf5_access_properties
  mpi_comm_loc = MPI_COMM_WORLD
  if(present(mpi_comm_in)) mpi_comm_loc = mpi_comm_in
  mpi_info_loc = MPI_INFO_NULL
  if(present(mpi_info_in)) mpi_info_loc = mpi_info_in
  test = .false.; if(present(test_in)) test = test_in;
  !> open HDF5 file 
  call HDF5_open(filename,file_id,ierr,create_access_plist_in=create_access_plist,&
  mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
  !> read the simulation time
  call HDF5_real_reading(file_id,sim%time,"/time",mpi_rank=sim%my_id,n_mpi_tasks=sim%n_cpu)
  !> get number of groups
  call H5Gopen_f(file_id,"/groups/",group_id,h5err)   
  call H5Gget_info_f(group_id,storage_type,n_groups,max_corder,h5err)
  call H5Gclose_f(group_id,h5err)  
  !> allocate simulation groups
  if(allocated(sim%groups)) deallocate(sim%groups); allocate(sim%groups(n_groups));
  do ii=1,n_groups
    !> read and load group datasets. We assume that the ith-particle group of 
    !> all tasks is defined by the same unique value stored in the hdf5 file
    ierr = 0; write(group_name,'(A,i0.3,A)') "/groups/",ii,"/";
    call HDF5_allocatable_char_reading(file_id,particle_type_str,trim(group_name)//"type")
    call HDF5_integer_reading(file_id,sim%groups(ii)%Z,trim(group_name)//"Z")
    call HDF5_real_reading(file_id,sim%groups(ii)%mass,trim(group_name)//"mass")
    call HDF5_char_reading(file_id,sim%groups(ii)%ad%suffix,trim(group_name)//"adas_suffix")
    if((len_trim(sim%groups(ii)%ad%suffix).gt.0).and.(.not.test)) then
      sim%groups(ii)%ad  = read_adf11(sim%my_id,sim%groups(ii)%ad%suffix)
      sim%groups(ii)%cor = coronal(sim%groups(ii)%ad) 
    endif
    !> compute the number of particles per processor and allocate particle array
    call HDF5_get_dataset_rank_dims(file_id,trim(group_name)//"i_elm",&
    rank,n_particles_tot,n_particles_max)
    n_particles_per_proc = int(n_particles_tot(1))/sim%n_cpu
    n_particles_per_proc(master_task+1) = int(n_particles_tot(1)) - &
    (sim%n_cpu-1)*n_particles_per_proc(master_task+1)
    offset = int(sum(n_particles_per_proc(1:sim%my_id)),kind=HSIZE_T)
    n_particles_hsizet = int(n_particles_per_proc(sim%my_id+1),kind=HSIZE_T)
    !> allocate particle list and initialise to 0
    select case (trim(particle_type_str))
    case ("particle_kinetic")
      allocate(particle_kinetic::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_kinetic_leapfrog")
      allocate(particle_kinetic_leapfrog::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc")
      allocate(particle_gc::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc_vpar")
      allocate(particle_gc_vpar::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc_Qin")
      allocate(particle_gc_Qin::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_fieldline")
      allocate(particle_fieldline::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_kinetic_relativistic")
      allocate(particle_kinetic_relativistic::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case ("particle_gc_relativistic")
      allocate(particle_gc_relativistic::sim%groups(ii)%particles(n_particles_per_proc(sim%my_id+1)))
    case default
      write(*,*) "Error: missing type name declaration ",trim(particle_type_str)," for reading: ABORT!"
      call MPI_Abort(mpi_comm_loc,errorcode,ierr)
    end select
    call initialize_particle_list_to_zero(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr)
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 1D array
    call HDF5_allocatable_array1D_reading_int(file_id,i_elm_arr,trim(group_name)//"i_elm",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading_int(file_id,i_life_arr,trim(group_name)//"i_life",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading_int(file_id,q_arr,trim(group_name)//"q",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: float 1D array
    call HDF5_allocatable_array1D_reading_r4(file_id,t_birth_arr,trim(group_name)//"t_birth",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: double 1D array
    call HDF5_allocatable_array1D_reading(file_id,weight_arr,trim(group_name)//"weight",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,v_1d_arr,trim(group_name)//"v",&
    reqdims_in=[n_particles_hsizet],start=[offset])   
    call HDF5_allocatable_array1D_reading(file_id,E_arr,trim(group_name)//"E",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,mu_arr,trim(group_name)//"mu",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,vpar_arr,trim(group_name)//"Vpar",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,B_norm_arr,trim(group_name)//"B_norm",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,vpar_m_arr,trim(group_name)//"Vpar_m",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    call HDF5_allocatable_array1D_reading(file_id,Bn_k_arr,trim(group_name)//"Bn_k",&
    reqdims_in=[n_particles_hsizet],start=[offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 2D array
    call HDF5_allocatable_array2D_reading(file_id,st_arr,trim(group_name)//"st",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,x_arr,trim(group_name)//"x",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,B_hat_prev_arr,trim(group_name)//"B_hat_prev",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,x_m_arr,trim(group_name)//"x_m",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,Astar_m_arr,trim(group_name)//"Astar_m",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,Astar_k_arr,trim(group_name)//"Astar_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,dBn_k_arr,trim(group_name)//"dBn_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,Bnorm_k_arr,trim(group_name)//"Bnorm_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,E_k_arr,trim(group_name)//"E_k",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    call HDF5_allocatable_array2D_reading(file_id,v_2d_arr,trim(group_name)//"v",&
    reqdims_in=[n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,offset])
    !> Read particle base datasets from HDF5 and fill the particle lists: integer 3D array
    call HDF5_allocatable_array3D_reading(file_id,dAstar_k_arr,trim(group_name)//"dAstar_k",&
    reqdims_in=[n1_HSIZE_T,n1_HSIZE_T,n_particles_hsizet],start=[i0_HSIZE_T,i0_HSIZE_T,offset])
    !> fill particle list from arrays
    call particle_list_from_arrays(n_particles_per_proc(sim%my_id+1),sim%groups(ii)%particles,ierr,&
    i_elm_arr=i_elm_arr,i_life_arr=i_life_arr,t_birth_arr=t_birth_arr,weight_arr=weight_arr,&
    x_arr=x_arr,st_arr=st_arr,q_arr=q_arr,v_1d_arr=v_1d_arr,E_arr=E_arr,mu_arr=mu_arr,&
    vpar_arr=vpar_arr,B_norm_arr=B_norm_arr,vpar_m_arr=vpar_m_arr,B_hat_prev_arr=B_hat_prev_arr,&
    v_2d_arr=v_2d_arr,x_m_arr=x_m_arr,Astar_m_arr=Astar_m_arr,Astar_k_arr=Astar_k_arr,&
    Bn_k_arr=Bn_k_arr,dBn_k_arr=dBn_k_arr,Bnorm_k_arr=Bnorm_k_arr,&
    E_k_arr=E_k_arr,dAstar_k_arr=dAstar_k_arr)
    !> deallocate structures
    call deallocate_particle_arrays(n_particles_per_proc(sim%my_id+1),i_elm_arr,&
    i_life_arr,q_arr,t_birth_arr,weight_arr,v_1d_arr,E_arr,mu_arr,vpar_arr,&
    B_norm_arr,vpar_m_arr,st_arr,x_arr,B_hat_prev_arr,v_2d_arr,x_m_arr,&
    Astar_m_arr,Astar_k_arr,Bn_k_arr,dBn_k_arr,Bnorm_k_arr,E_k_arr,dAstar_k_arr)   
  enddo
  !> clean-up
  call HDF5_close(file_id)
  if(allocated(n_particles_tot))      deallocate(n_particles_tot)
  if(allocated(n_particles_max))      deallocate(n_particles_max)
  if(allocated(n_particles_per_proc)) deallocate(n_particles_per_proc)
  if(allocated(particle_type_str))    deallocate(particle_type_str)
end subroutine read_simulation_hdf5

!> Original procedure for getting '/time' from a file. 
!> Does not alter the units in any way code works for jorek 
!> and particle restart files, and returns values in different units for both
function get_simulation_hdf5_time_original(filename) result(time)
  use hdf5, only: HID_T,H5F_ACC_RDONLY_F
  use hdf5, only: h5open_f,h5fopen_f,h5fclose_f,h5close_f
  use hdf5_io_module, only: HDF5_real_reading
  character*(*), intent(in)  :: filename
  real*8 :: time
  integer(HID_T) :: file
  integer :: hdferr
  call h5open_f(hdferr)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  call HDF5_real_reading(file,time,'/time')
  call h5fclose_f(file,hdferr)
  call h5close_f(hdferr)
end function get_simulation_hdf5_time_original

!> !> Get '/time' from a file. Does not alter the units in any way
!> code works for jorek and particle restart files, and returns values in
!> different units for both
!> inputs:
!>   filename:            (character) name of the HDF5 file
!>   create_access_plist: (logical) create HDF parameter list (for MPIO)
!>                        default: false (H5P_DEFAULT_F), if true create
!>                        H5P_FILE_ACCESS_F
!>   mpi_comm:            (integer) MPI communicator identifier (for MPIO)
!>   mpi_info:            (integer) MPI parameter object (for MPIO)
!>   my_id:               (integer) MPI task identifier (for MPIO)
!>   n_cpu:               (integer) number of MPI task (for MPIO)
!> outputs:
!>   time:     (real8) restart simulation time
function get_simulation_hdf5_time(filename,use_hdf5_access_properties,&
mpi_comm_loc,mpi_info_loc,my_id,n_cpu) result(time)
  use mpi
  use hdf5,           only: HID_T
  use hdf5_io_module, only: HDF5_open,HDF5_close,HDF5_real_reading
  implicit none
  character(len=*),intent(in) :: filename
  integer,intent(in),optional :: mpi_comm_loc,mpi_info_loc
  integer,intent(in),optional :: my_id,n_cpu
  logical,intent(in),optional :: use_hdf5_access_properties
  real*8                      :: time
  integer                     :: h5err
  integer(HID_T)              :: file_id
  if(present(use_hdf5_access_properties).and.present(mpi_comm_loc).and.&
  present(mpi_info_loc).and.present(my_id).and.present(n_cpu)) then
    call HDF5_open(trim(filename),file_id,h5err,&
    create_access_plist_in=.not.use_hdf5_access_properties,&
    mpi_comm_in=mpi_comm_loc,mpi_info=mpi_info_loc)
    call HDF5_real_reading(file_id,time,"/time",mpi_rank=my_id,n_mpi_tasks=n_cpu)
  else
    call HDF5_open(trim(filename),file_id,h5err)
    call HDF5_real_reading(file_id,time,"/time")
  endif
  call HDF5_close(file_id)
end function get_simulation_hdf5_time

end module mod_particle_io
