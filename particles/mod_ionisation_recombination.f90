!> This module contains routines to calculate ionisation
!> and recombination probabilities of particles in a specific time
module mod_ionisation_recombination
use mod_openadas
!implicit none
!private
!public new_charge
!public fields_interp_ne_Te

!contains

  ! use mod_edge_elements
  ! use mod_io_actions, only: io_action
  ! use mod_sampling
  ! use mod_particle_types
  ! use constants, only: TWOPI, K_BOLTZ, ATOMIC_MASS_UNIT
  ! use mod_rng, only: type_rng, setup_shared_rngs
  ! use mod_boundary, only: wall_normal_vector
  ! use mod_atomic_elements !mod_elements !< Chemical elements
  ! use mod_particle_sim
  ! use mod_event
  ! use mod_find_rz_nearby, only: find_rz_nearby
  use mod_parameters, only : n_vertex_max, n_elements_max, n_order,n_nodes_max
  use nodes_elements
  use data_structure, only : type_node, type_element
  


  implicit none

    type (type_element)      :: element
 type (type_node)         :: nodes(n_vertex_max)
  
  private
  public rec_rate_global
  public rec_rate_local !(ife or i_elm) ! size n_elements_max
  public rec_mom_local
  public rec_energy_local
  public rec_v_R
  public rec_v_Z  
  public rec_v_phi
  
  public new_charge
  public fields_interp_ne_Te
  !public  :: particle_recombination

  ! Extend type
  ! type, extends(io_action) :: particle_recombination
   
    ! class(type_rng), dimension(:), allocatable :: rng  !< one RNG per openmp thread
   
    ! ! number of simulation particles/s to puff across all processes
    ! integer :: n_puff = -1 
    ! ! Average fueling rate: 9.7d22; max fueling rate 18d22
    ! real*8  :: fueling_rate = -1.d0
    ! real*8  :: R = -1.d0, Z = -1.d0, phi = -1.d0
    ! real*8  :: valve_r = -1.d0  !< radius of gas valve
    ! real*8  :: last_time = 0.d0 !< When did we puff last 
    ! real*8 :: last_diag_time = 0.d0 !< Last time of output of diagnostics

  ! contains
    ! procedure :: do => do_particle_recombination
  ! end type particle_recombination

  ! interface particle_recombination
    ! module procedure new_particle_recombination
  ! end interface particle_recombination
  
  real*8, dimension(n_vertex_max, n_order+1, n_nodes_max) :: rec_rate_global
  real*8, dimension(n_elements_max) :: rec_rate_local, rec_mom_local, rec_energy_local
  real*8, dimension(n_elements_max) :: rec_v_R, rec_v_Z, rec_v_phi  


  
contains

  	! rec_rate_local = 0.d0
	! rec_mom_local  = 0.d0
	! rec_energy_local= 0.d0
	! rec_v_R= 0.d0
	! rec_v_Z= 0.d0
	! rec_v_phi= 0.d0









! already existing --------------->>>>

!> Calculate new charge state at a specific density, temperature and timestep
function new_charge(z, ad, electron_density, electron_temperature, timestep, ran2) result(z_new)
implicit none

integer, intent(in)              :: z !< Old charge state
type(ADF11_all), intent(in)      :: ad !< ADF11 data
real*8, intent(in)               :: electron_density !< log10 Electron density in m^-3
real*8, intent(in)               :: electron_temperature !< log10 Electron temperature in K
real*8, intent(in)               :: timestep !< Timestep in s
real*8, intent(in), dimension(2) :: ran2 !< Random numbers to use to select a new charge or not
integer :: z_new

real*8 :: prob(2), acd, scd

z_new = z
! probabilities of recombination and ionisation events
call ad%ACD%interp(z,   electron_density, electron_temperature, acd)
call ad%SCD%interp(z,   electron_density, electron_temperature, scd)
prob = 1.d0 - exp(-[acd, scd] * 10.d0**electron_density * timestep)
z_new = z_new + sum([-1, 1], mask=(prob .gt. ran2))

end function new_charge

pure subroutine fields_interp_ne_Te(fields, time, s, t, phi, i_elm, n_e, T_e)
use mod_fields
use phys_module, only: central_density
use constants
class(fields_base), intent(in)                    :: fields
real*8, intent(in)                                :: time, s, t, phi
integer, intent(in)                               :: i_elm
real*8, intent(out)                               :: n_e !< electron density [m^-3]
real*8, intent(out)                               :: T_e !< electron temperature [K]

real*8, dimension(2) :: P, P_s, P_t, P_phi, P_time
real*8               :: R, R_s, R_t, Z, Z_s, Z_t
call fields%interp_PRZ(time,i_elm,&
#if (JOREK_MODEL == 400)
      [5,8],& ! electron temperature
#else
      [5,6],& ! electron temperature + ion temperature (assumed equal)
#endif
          2,s,t,phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)

n_e = max(P(1) * central_density * 1d20,1d16)                           ! plasma density [1/m^3], capped against negative
T_e = max(P(2)/(2.d0*MU_ZERO*central_density*1.d20)/K_BOLTZ, 1.d0) ! temperature capped against going negative
#if (JOREK_MODEL == 400)
T_e = T_e*2d0 ! P(1) contains the electron temperature, reverse previous correction
#endif
end subroutine fields_interp_ne_Te
end module mod_ionisation_recombination
