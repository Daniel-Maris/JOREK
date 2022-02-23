!> Module for controlling a paramtere in the plasma
!> As an example, this module will control the time-varying gas puff rate at the locations where gas valves will be.
!> But it is aimed at being usable for various other applications

!note: all variables that you allocate here cannot be allocated at the same time in the module where you use mod_controller

module mod_controller
    use mod_particle_puffing
    use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
    use constants,   only: MU_ZERO, MASS_PROTON

    implicit none
        
    ! controller parameters
    logical                :: use_controller
    
    ! usage parameters
    type(particle_puffing) :: gas_puff
    real*8                 :: rho_norm, t_norm, n_norm

contains

subroutine controller_function(use_controller)
    implicit none 
  
    logical,intent(in) :: use_controller
    
    n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
    rho_norm  = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
    t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek 

    ! determine what the controller does
    ! note: if you want to change parameters that are part of an action within an event, the event must be a pointer (use new_event_ptr() instead of event())
    if (use_controller) then
        write(*,*) "The controller_function works. The controller is on"
        gas_puff%fueling_rate = 60.d21    ! we can adjust this%fuelling_rate here depending on the required input signal
        gas_puff%t_puff_slope = 1.d-3
        gas_puff%t_puff_start = 10*t_norm
        
        ! call the setpoint on this timestep 
        ! measure value
        ! determine error
        ! calculate output controller
        ! save error and integral
        ! calculate new output = setpoint + controller_output 
    else
        write(*,*) "The controller_function works. The controller is off"
    endif !(use_controller)
end subroutine controller_function

end module