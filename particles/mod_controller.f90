!> Module for controlling a parameter in the plasma
!> As an example, this module will control the time-varying gas puff rate at the locations where a gas valve is present
!> But it is aimed at being usable for various other applications

!note: all variables that you allocate here cannot be allocated at the same time in the module where you use mod_controller

module mod_controller
    
    use mod_particle_puffing
    use mod_particle_sim
    use phys_module, only: index_now

    implicit none
        
    ! controller parameters
    logical                          :: use_controller
    character(len=256)               :: controller_timedependentsignal_file

    ! use-case specific parameters
    type(particle_puffing)           :: gas_puff
    class(particle_puffing), pointer :: this
    real*8                           :: t_norm

    !defining a type for a general time dependent signal
    type :: time_dependent_signal
    integer            :: len      !< Number of points in numerical time trace.
    real*8, allocatable:: time(:)  !< time-values of numerical time trace
    real*8, allocatable:: signal(:)  !< current-values of numerical time trace
    end type time_dependent_signal

contains

subroutine controller_function(use_controller,this,sim, time_dependent_signal_controller)
    use profiles, only: readProf
    
    implicit none 
    
    class(particle_puffing) , intent(inout) :: this
    type(particle_sim), intent(inout)       :: sim
    type(time_dependent_signal), intent(inout) ::time_dependent_signal_controller
    logical,intent(in) :: use_controller

    controller_timedependentsignal_file = '/home/ITER/vanhooe/Documents/Datafiles/datafile.dat' ! Note: check the directory to the datafile
    
    !> The following if-statement determines what the controller does
    !> note: if you want to change parameters that are part of an action within an event, the event must be a pointer (use new_event_ptr() instead of event())
    if (use_controller) then
        write(*,*) "The controller_function works. The controller is on"
        !> Example how to change the defined !time-dependent! signal                       the lines below were used for run 16 and part of commit 10
        !gas_puff%fueling_rate = 60.d21    ! we can adjust this%fuelling_rate here depending on the required input signal
        !gas_puff%t_puff_slope = 1.d-3
        !gas_puff%t_puff_start = 10*t_norm
        
        !test to make sure sim%time and index_now work properly and are usable to make the controller function time dependent
        write(*,"(A,g12.4)") "test voor controller, this is the time now", sim%time 
        write(*,"(A,g12.4)") "test voor controller, this is index_now", index_now 
       
        !> Example of how to define a certain time dependent fuelling rate signal from within the controller function. 
        !> Currently it is the same function as the time-dependent_puff function within mod_particle_puffing (but copied and renamed below)
        !> To use it, puff_t_dependent must be set .false. in my_example.
        !gas_puff%fueling_rate = time_dependent_puff_controller(25.d21,sim%time, 10*t_norm,500*t_norm, 20.d21)
         
        !> Example of how to import a time dependent signal using a datafile
        call readProf(time_dependent_signal_controller%time, time_dependent_signal_controller%signal, &
        time_dependent_signal_controller%len, controller_timedependentsignal_file)
        write(*,"(A,g12.4,A,g12.4)") "current time signal", time_dependent_signal_controller%time(index_now), "current signal", time_dependent_signal_controller%signal(index_now) 
        gas_puff%fueling_rate = time_dependent_signal_controller%signal(index_now)

        ! Next steps to implement in the controller function:
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

pure function time_dependent_puff_controller(max_puff,time, t_puff_start,t_puff_slope, min_puff) result(to_puff)
real*8,intent(in)   :: max_puff, min_puff
real*8              :: to_puff
real*8,intent(in)    :: t_puff_start,t_puff_slope
real*8,intent(in)    :: time

if (time-(t_puff_start+t_puff_slope) .ge. 0.d0) then
	to_puff = max_puff
elseif (time-t_puff_start .ge. 0.d0) then
	to_puff = min_puff+ (max_puff -min_puff) * (time-t_puff_start)/(t_puff_slope)  
else
    to_puff = min_puff !default = 0.d0
endif
end function time_dependent_puff_controller

end module