!> get SOL currents from a JOREK restart file
program jorek2_SOLcurrent

use data_structure
use phys_module
use basis_at_gaussian
use elements_nodes_neighbours
use mod_neighbours
use mod_import_restart
use mod_log_params
use equil_info, only : get_psi_n, ES
use mod_interp
use mod_boundary
use mod_position

implicit none

!type (type_node_list)   ,     pointer :: node_list
!type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list    

type(t_pol_pos_list), target :: initial_bnd_pos

character(len=512) :: s
real*8, allocatable :: Phi_start_list(:)
real*8 :: R_1, R_2, Z_1, Z_2, Phi_1, Phi_2 !<start and end points of the traced fieldline
real*8 :: TT_1, TT_2, ne_1, ne_2 !<quantities at the beginning and end
real*8 :: L                                       !<integrals along the traced line: connection length L
integer :: my_id
!real*8  :: rr, zz, psi
integer :: i, j, iside_i, iside_j, ip, np, i_tor, i_harm, i_var_psi = 1,i_var_n=5,i_var_T=6
integer :: i_elm, ifail, n_phi_steps, i_elm_out, i_elm_prev, i_elm_tmp,i_steps
integer :: n_phi0, i_phi0, n_elm_pts
real*8  :: R_line, Z_line, s_line, t_line, p_line, R_mid, Z_mid, s_mid, t_mid, p_mid, s_out, t_out
!real*8, allocatable :: R_start(:), Z_start(:), P_start(:)
real*8  :: R, R_s, R_t, Z, Z_s, Z_t 
real*8,dimension(3) :: P, P_s, P_t, P_phi
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local
real*8  :: atmp, cur_pert
real*8  :: psi_out
integer :: ierr
real*8  :: dr,R_old,Z_old !< integration variables
logical :: in_domain !< keeps track of whether tracer position is inside the domain or on its boundary
real*8  :: B(2)!B(3) !< Magnetic field [T]
real*8  :: Bdotn !< B dot outward pointing normal to the boundary
real*8  :: inv_st_jac, psi_R, psi_Z
logical :: debug = .false.!.true.


write(*,*) '***************************************'
write(*,*) '* JOREK2_SOLcurrent                   *'
write(*,*) '***************************************'

!allocate(node_list)
!allocate(element_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

my_id=0

! --- Initialize mode and mode_type arrays
call det_modes()

call initialise_parameters(my_id,  "__NO_FILENAME__")
call log_parameters(my_id)

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis

call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

allocate(element_neighbours(4,element_list%n_elements))

element_neighbours = 0

do i=1,element_list%n_elements
  
  do j=i+1,element_list%n_elements
    
    if (neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j)) then
      element_neighbours(iside_i,i) = j
      element_neighbours(iside_j,j) = i
    endif
    
  enddo
enddo

! NUMERICAL PARAMETERS
n_phi_steps   = 1500 !number of steps in the toroidal direction
n_phi0        = 2    !number of toroidal angles to start from for each poloidal bnd point
n_elm_pts     = 1!10   !number of starting points along each bnd element

! step at maximum at delta_phi (smaller step is element boundary is crossed)
delta_phi = 2.d0 * PI / float(n_period*n_phi_steps)

allocate(Phi_start_list(n_phi0))
do i_phi0 = 1, n_phi0
  Phi_start_list(i_phi0) = float(i_phi0-1)/float(n_phi0) * 2.d0 * PI / float(n_period)
end do

initial_bnd_pos = bnd_pos(node_list, element_list, bnd_node_list, bnd_elm_list, n_elm_pts)

np = size(initial_bnd_pos%pos(1,:))

write(*,*) 'toroidal angles per starting position: ',n_phi0
write(*,*) 'number of fieldlines to trace:         ',n_phi0*np

i_var_psi = 1

!Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
!do i=1,node_list%n_nodes
!  Rmin = min(Rmin,node_list%node(i)%x(1,1))
!  Rmax = max(Rmax,node_list%node(i)%x(1,1))
!  Zmin = min(Zmin,node_list%node(i)%x(1,2))
!  Zmax = max(Zmax,node_list%node(i)%x(1,2))
!enddo

mode(1) = 0
do i=1,(n_tor-1)/2
  mode(2*i)   = i * n_period
  mode(2*i+1) = i * n_period
enddo
write(*,*) ' modes   : ',mode
write(*,*) ' nperiod : ',n_period

  
!call begplt('poincare.ps')

! --- Open the output file to which the SOL current data will be written in ascii format
open(21,file='SOLcurrents.dat')
write(21,*) '#  R_1  Z_1 Phi_1 R_2 Z_2 Phi_2 ne_1 ne_2 T_1 T_2 L' !add other terms and total SOL current here
!call nframe(21,11,1,Rmin,Rmax,Zmin,Zmax,'Poincare',8,'R [m]',4,'Z [m]',4)

! --- Trace the fieldlines
L_p: do ip=1,np !loop over the starting positions
  
  L_phi0: do i_phi0 = 1, n_phi0
    if(debug) write(*,*) 'phi loop start, angle ',i_phi0,'/',n_phi0
  
    s_line = initial_bnd_pos%pos(1,ip)%s
    t_line = initial_bnd_pos%pos(1,ip)%t
    i_elm  = initial_bnd_pos%pos(1,ip)%ielm

    Phi_1 = Phi_start_list(i_phi0)
    p_line = Phi_1
    
    !call interp_RZ(node_list,element_list,i_elm,s_line,t_line, R, R_s, R_t, Z, Z_s, Z_t) !gets R and Z at the starting location
    ! Calculate the derivatives to R and Z
    call interp_PRZ(node_list,element_list,i_elm,[i_var_psi,i_var_n,i_var_T],3,s_line,t_line,p_line,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
    ne_1 = P(2)
    TT_1 = P(3)
    
    R_1 = R
    Z_1 = Z
    R_old = R
    Z_old = Z

    if(debug) write(*,*) 'R,Z',R,Z
    
    !determine sign of delta_phi
    
    inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
    psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
    psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac

    B     = [+psi_Z, -psi_R] * 1.d0/R
    Bdotn = dot_product(B,initial_bnd_pos%pos(i,ip)%bnd_normal)
    delta_phi = sign(delta_phi, Bdotn)
    if(debug) write(*,*) 'B:',B, 'Bdotn:',Bdotn, 'delta_phi:',delta_phi,'F0:', F0, 'Bphi',F0/R
    
    !resetting running integrals
    L = 0.d0
  
    i_steps = 0

    in_domain = .true.
    do while (in_domain)
    
      i_steps = i_steps + 1
      if(i_steps .ge. 1000*n_phi_steps) then
        write(*,*) 'ERROR: could not find connecting wall'
        cycle
      end if
    
      
      call step(i_elm,s_line,t_line,p_line,delta_phi,delta_s,delta_t)

      s_mid = s_line + 0.5d0 * delta_s
      t_mid = t_line + 0.5d0 * delta_t
      p_mid = p_line + 0.5d0 * delta_phi
          
      call step(i_elm,s_mid,t_mid,p_mid,delta_phi,delta_s,delta_t)
      
      small_delta_s = 1.d0
    
      if  (s_line + delta_s .gt. 1.d0) then
    
        small_delta_s = (1.d0 - s_line)/delta_s

      elseif  (s_line + delta_s .lt. 0.d0) then

        small_delta_s = abs(s_line/delta_s)

      endif
    
      small_delta_t = 1.d0

      if  (t_line + delta_t .gt. 1.d0)  then
    
        small_delta_t = (1.d0 - t_line)/delta_t

      elseif  (t_line + delta_t .lt. 0.d0)  then

        small_delta_t = abs(t_line/delta_t)
    
      endif      

      small_delta = min(small_delta_s, small_delta_t)

      !	write(*,'(A,5e16.8)') ' small delta : ',small_delta,delta_s,delta_t,s_line,t_line

      if (small_delta .lt. 1.d0)  then ! moving the tracked position by delta_t and delta_s moves it out of the element so correct that the tracked position ends up on the element instead

        !s_mid = s_line + 0.5d0 * small_delta * delta_s
        !t_mid = t_line + 0.5d0 * small_delta * delta_t
        !p_mid = p_line + 0.5d0 * small_delta * delta_phi
          
        !call step(i_elm,s_mid,t_mid,p_mid,delta_phi,delta_s,delta_t) ! why is this in here? better accuracy at the cost of getting a possibly getting new position s,t \in [0,1], meaning the position is not updated as it's inside the wrong if statement
    
        if (small_delta_s .lt. small_delta_t) then

          if (s_line + delta_s .gt. 1.d0) then
    
            s_line = 1.d0
            t_line = t_line + small_delta * delta_t
            p_line = p_line + small_delta * delta_phi
      
            call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

            i_elm_prev = i_elm      
            i_elm      = element_neighbours(2,i_elm_prev)
            if ( i_elm == 0 ) then !tracer is at domain boundary
              if(debug) write(*,*) 'end of field line tracing (0)'
              in_domain = .false.
            else
              i_elm_tmp  = element_neighbours(4,i_elm)
        
              if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'

              s_line = 0.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
        
              if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
            endif
          elseif (s_line + delta_s .lt. 0.d0) then

            s_line = 0.d0
            t_line = t_line + small_delta * delta_t
            p_line = p_line + small_delta * delta_phi

            call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)
    
            i_elm_prev = i_elm      
            i_elm      = element_neighbours(4,i_elm_prev)
            if ( i_elm == 0 ) then ! at domain boundary
              if(debug) write(*,*) 'end of field line tracing (1)'
              in_domain = .false.
            else
              i_elm_tmp  = element_neighbours(2,i_elm)
        
              if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

              s_line = 1.d0
              
              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
        
              if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
            end if

          endif

        else

          if (t_line + delta_t .gt. 1.d0) then
    
            s_line = s_line + small_delta * delta_s
            t_line = 1.d0
            p_line = p_line + small_delta * delta_phi

            call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

            i_elm_prev = i_elm      
            i_elm      = element_neighbours(3,i_elm_prev)
            if ( i_elm == 0 ) then
              if(debug) write(*,*) 'end of field line tracing (2)'
              in_domain = .false.
            else
              i_elm_tmp  = element_neighbours(1,i_elm)
        
              if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'

              t_line = 0.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
        
              if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
            end if

          elseif (t_line + delta_t .lt. 0.d0) then

            s_line = s_line + small_delta * delta_s	
            t_line = 0.d0
            p_line = p_line + small_delta * delta_phi

            call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

            i_elm_prev = i_elm      
            i_elm      = element_neighbours(1,i_elm_prev)
            if ( i_elm == 0 ) then
              if(debug) write(*,*) 'end of field line tracing (3)'
              in_domain = .false.
            else
              i_elm_tmp  = element_neighbours(3,i_elm)
        
              if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)',i_elm_prev,i_elm

              t_line = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
        
              if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
            end if

          endif
      
        endif
        
      else      

        s_line = s_line + delta_s
        t_line = t_line + delta_t
        p_line = p_line + delta_phi

        small_delta = 1.d0 !this should not be necessary
    
      endif
      
      delta_phi_local = small_delta * delta_phi
      !if(debug) write(*,*) 'prior to dr: ',i_elm,s_line,t_line
      call interp_RZ(node_list,element_list,i_elm_prev,s_line,t_line,R,Z)
      dr = sqrt((R-R_old)**2 + (Z-Z_old)**2 + (R*delta_phi_local)**2)
      !write(*,*) 'prior to L: ',R, R_old, Z, Z_old, delta_phi_local,i_elm_prev, L, dr

      L = L +  1.d0 * dr !connection length
      ! other running integrals

      R_old = R
      Z_old = Z

      if ( i_elm == 0 ) then
        if(debug) write(*,*) '(4), do nothing'
      end if

      !write(*,'(A,5e16.8)') ' s,t : ',s_line,t_line
    
      
      !if (i_steps .gt. 8) write(*,'(A,5i6)') ' WARNING : isteps ',i_lines,i_turn,i_phi,i_steps,i_elm
            
    enddo ! traced a fieldline for one bnd point for 1 starting angle phi

    call interp_RZ(node_list,element_list,i_elm_prev,s_line,t_line,R_2,Z_2)
    Phi_2 = p_line
    
    call interp(node_list,element_list,i_elm_prev,[i_var_psi,i_var_n,i_var_T],3,s_line,t_line,P,P_s,P_t,P_st,P_ss,P_tt)
    ne_2 = P(2)
    TT_2 = P(3)

    !call var_value(i_elm,1,s_line,t_line,p_line,psi_out)

    !Tp(ip)  = atan2( Z_line - ES%Z_axis, R_line - ES%R_axis)
    !Pp(ip)  = get_psi_n(psi_out, Z_line)

    !if ( i_elm == 0 ) then
    !  write(*,*) 'ERROR, i_elm == 0'
    !  exit
    !end if
    
    !put write statements here
    write(21,'(7e18.8)') R_1, Z_1, Phi_1, R_2, Z_2, Phi_2, ne_1, ne_2, TT_1, TT_2, L

  enddo L_phi0 !traced all starting angles for one bnd point
  
  !if(debug) write(*,*) 'ip/np',ip,np
  if (mod(ip,(np/10+1)) == 0) then
    write(*,*) 'starting point = ',ip,' of ',np,' at (R,Z) = (',R_1,',',Z_1,') and L = ',L
  endif
  
  !write(21,*)
  
  !call lincol(mod(i_lines,8))
   
  !if (iplot_type .eq. 1) then
  !  call pplot(1,1,Rp,Zp,ip,1)
  !else
  !  call pplot(1,1,Pp,Tp,ip,1)
  !endif
  
end do L_p !traced all starting points and angles

close(21)
!close(22)
!call finplt

end program jorek2_SOLcurrent



subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t)
use mod_parameters
use elements_nodes_neighbours
use phys_module
use mod_interp

implicit none

integer :: i_var_psi, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
real*8 :: R,R_s,R_t,Z,Z_s,Z_t
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, Zjac
real*8 :: AR0_Z, AR0_p, AR0_s, AR0_t, AZ0_R, AZ0_p, AZ0_s, AZ0_t, A30_R, A30_Z, BR0, BZ0, Bp0, Fprof

i_var_psi = 1

call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,Z,Z_s,Z_t)

Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

psi_s = P0_s 
psi_t = P0_t 

#ifdef fullmhd
  call interp(node_list,element_list,i_elm,var_AR,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  AR0_s = P0_s 
  AR0_t = P0_t 

  call interp(node_list,element_list,i_elm,var_AZ,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  AZ0_s = P0_s 
  AZ0_t = P0_t 

  AR0_p = 0.d0
  AZ0_p = 0.d0

  call interp(node_list,element_list,i_elm,710,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  Fprof = P0
#endif

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
  psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
  psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)

#ifdef fullmhd
  call interp(node_list,element_list,i_elm,var_AR,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
  AR0_s = AR0_s + Pcos_s * cos(mode(i_harm)*p_in)
  AR0_t = AR0_t + Pcos_t * cos(mode(i_harm)*p_in)
  AR0_p = AR0_p - Pcos   * sin(mode(i_harm)*p_in) * mode(i_harm)
  call interp(node_list,element_list,i_elm,var_AR,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
  AR0_s = AR0_s + Psin_s * sin(mode(i_harm+1)*p_in)
  AR0_t = AR0_t + Psin_t * sin(mode(i_harm+1)*p_in)
  AR0_p = AR0_p + Psin   * cos(mode(i_harm+1)*p_in) * mode(i_harm+1)

  call interp(node_list,element_list,i_elm,var_AZ,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
  AZ0_s = AZ0_s + Pcos_s * cos(mode(i_harm)*p_in)
  AZ0_t = AZ0_t + Pcos_t * cos(mode(i_harm)*p_in)
  AZ0_p = AZ0_p - Pcos   * sin(mode(i_harm)*p_in) * mode(i_harm)
  call interp(node_list,element_list,i_elm,var_AZ,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
  AZ0_s = AZ0_s + Psin_s * sin(mode(i_harm+1)*p_in)
  AZ0_t = AZ0_t + Psin_t * sin(mode(i_harm+1)*p_in)
  AZ0_p = AZ0_p + Psin   * cos(mode(i_harm+1)*p_in) * mode(i_harm+1)
#endif

enddo

#ifdef fullmhd
AR0_Z = ( - R_t * AR0_s  + R_s * AR0_t ) / Zjac
AZ0_R = (   Z_t * AZ0_s  - Z_s * AZ0_t ) / Zjac
A30_R = (   Z_t * psi_s  - Z_s * psi_t ) / Zjac
A30_Z = ( - R_t * psi_s  + R_s * psi_t ) / Zjac

BR0 = ( A30_Z - AZ0_p )/ R
BZ0 = ( AR0_p - A30_R )/ R
Bp0 = ( AZ0_R - AR0_Z )       +   Fprof / R

! dR/Rdphi = B_R / B_phi ; dZ/Rdphi = B_Z / B_phi
! ds = (Z_t dR - R_t dZ) / Zjac ; dt = ( -Z_s dR + R_s dZ) / Zjac
delta_s =  ( Z_t*BR0 - R_t*BZ0) / ( Bp0 * Zjac ) * R * delta_p
delta_t =  (-Z_s*BR0 + R_s*BZ0) / ( Bp0 * Zjac ) * R * delta_p     

#else
delta_s =   psi_t * R / (Zjac * F0) * delta_p
delta_t = - psi_s * R / (Zjac * F0) * delta_p
#endif

return
end subroutine step



subroutine var_value(i_elm,i_var,s_in,t_in,p_in,value_out)
use mod_parameters
use elements_nodes_neighbours
use phys_module
use mod_interp, only: interp

implicit none

integer :: i_var, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt
real*8 :: value_out

!call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

value_out = P0

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  value_out = value_out + Pcos * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  value_out = value_out + Psin * sin(mode(i_harm+1)*p_in)

enddo

return
end subroutine var_value
