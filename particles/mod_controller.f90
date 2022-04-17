!> Module for controlling a parameter in the plasma
!> As an example, this module will control the time-varying gas puff fueling rate at the locations where a gas valve is present
!> But it is aimed at being usable for various other applications

module mod_controller
    
    use mod_particle_puffing
    use mod_particle_sim
    use phys_module, only: index_now, central_density
    use data_structure
    use mod_fields
    
    implicit none

    !note: all variables that you allocate right here cannot be allocated at the same time in the module where you use mod_controller
    !> so it is preferred to allocate all used variables in the subroutine itself

    !defining a type for a general time dependent signal
    type :: time_dependent_signal
    integer                          :: len        !< Number of points in numerical time trace.
    real*8, allocatable              :: time(:)    !< time-values of numerical time trace
    real*8, allocatable              :: signal(:)  !< functionvalues of numerical time trace
    end type time_dependent_signal

contains

subroutine controller_function(use_controller,  sim, t_dep_signal_controller, contr_selfdefined, contr_usedatafile, &
                                contr_analytical, control_t_dep_signal_file, analytical_expression, analytical_len, &
                                analytical_tmax, controllerhasbeencalledbefore, previous_time_controller, controller_K_p, &
                                controller_K_i, controller_K_d, node_list, element_list, puff_t_dependent, t_norm, &
                                setpoint, max_value, min_value, controller_type, actuator_signal)
    
    implicit none 
    
    ! Allocating the function arguments
    type(particle_sim), intent(inout)          :: sim
    type(time_dependent_signal), intent(inout) :: t_dep_signal_controller
    type(type_node_list), intent(in)           :: node_list
    type (type_element_list), intent(in)       :: element_list
    logical, intent(in)                        :: use_controller
    logical, intent(in)                        :: contr_selfdefined 
    logical, intent(in)                        :: contr_usedatafile 
    logical, intent(in)                        :: contr_analytical 
    character(len=512), intent(in)             :: control_t_dep_signal_file
    character(len=512), intent(in)             :: analytical_expression
    integer, intent(in)                        :: analytical_len
    real*8, intent(in)                         :: analytical_tmax
    logical, intent(inout)                     :: controllerhasbeencalledbefore
    real*8, intent(inout)                      :: previous_time_controller
    real*8, intent(in)                         :: controller_K_p
    real*8, intent(in)                         :: controller_K_i
    real*8, intent(in)                         :: controller_K_d
    real*8, intent(in)                         :: setpoint
    logical, intent(in)                        :: puff_t_dependent
    real*8, intent(in)                         :: t_norm
    real*8, intent(in)                         :: max_value
    real*8, intent(in)                         :: min_value
    character, intent(in)                      :: controller_type
    real*8, intent(out)                        :: actuator_signal

    real*8                                     :: controller_error_prev
    real*8                                     :: input_signal
    real*8                                     :: controller_output
    
    !> The following if-statement determines what the controller does when it is called
    !> note: if you want to change parameters that are part of an action within an event, the event must be a pointer (use new_event_ptr() instead of event())
    if (use_controller) then
        if (sim%my_id .eq. 0) write(*,*) "The controller_function works. The controller is on"
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is the time now:", sim%time 
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is index_now:", index_now 
        
        if (controller_type == 'openloop') then

            call controller_input(contr_selfdefined, contr_usedatafile, contr_analytical, sim, t_norm, &
                                            controllerhasbeencalledbefore, t_dep_signal_controller, control_t_dep_signal_file, &
                                            analytical_expression, analytical_len, analytical_tmax, input_signal)

            actuator_signal = max(input_signal, min_value)
            actuator_signal = min(input_signal, max_value)

        else if (controller_type == 'closedloop') then

            if (.not. controllerhasbeencalledbefore) then
                controller_error_prev = 0.d0
            endif

            call closed_loop(sim, element_list, node_list, previous_time_controller, setpoint, controller_K_p, controller_K_i, controller_K_d, &
                            controller_output, controller_error_prev) 

            actuator_signal = max(controller_output, min_value)
            actuator_signal = min(controller_output, max_value)

        else if (controller_type == 'feedforward') then

            call controller_input(contr_selfdefined, contr_usedatafile, contr_analytical, sim, t_norm, &
                                            controllerhasbeencalledbefore, t_dep_signal_controller, control_t_dep_signal_file, &
                                            analytical_expression, analytical_len, analytical_tmax, input_signal)

            if (.not. controllerhasbeencalledbefore) then
                controller_error_prev = 0.d0
            endif

            call closed_loop(sim, element_list, node_list, previous_time_controller, setpoint, controller_K_p, controller_K_i, controller_K_d, &
                            controller_output, controller_error_prev) 

            actuator_signal = max(input_signal + controller_output, min_value)
            actuator_signal = min(input_signal + controller_output, max_value)
        else
            if (sim%my_id .eq. 0) write(*,*) "ERROR: The controller_type must be specified when use_controller is true."
            stop
        endif

        if (actuator_signal .eq. min_value) then
            if (sim%my_id .eq. 0) write(*,*) "WARNING: The actuator_signal is set to the minimal value because it cannot be smaller than that"
        endif
        if (actuator_signal .eq. max_value) then
            if (sim%my_id .eq. 0) write(*,*) "WARNING: The actuator_signal is set to the maximal value because it cannot be larger than that"
        endif

        if (sim%my_id .eq. 0) write(*,"(E16.8,E16.8,A52)") sim%time, actuator_signal, "#simulation time and actuator_signal"

        if (sim%my_id .eq. 0) write(*,*) "test for controller, end of the controller loop"
    else
        if (sim%my_id .eq. 0) write(*,*) "The controller_function works. The controller is called but it is not used, as the controller is turned off"
    endif !(use_controller)

    !> This makes sure that certain actions are only performed once and not during every timestep in which the controller is called
    controllerhasbeencalledbefore = .true.

end subroutine controller_function

subroutine controller_input(contr_selfdefined, contr_usedatafile, contr_analytical, sim, t_norm, controllerhasbeencalledbefore, &
                            t_dep_signal_controller, control_t_dep_signal_file, analytical_expression, analytical_len, &
                            analytical_tmax, input_signal)

use profiles, only: readProf, interpolProf
    
implicit none 
    
! Allocating the function arguments
type(particle_sim), intent(inout)          :: sim
type(time_dependent_signal), intent(inout) :: t_dep_signal_controller
logical, intent(in)                        :: contr_selfdefined 
logical, intent(in)                        :: contr_usedatafile 
logical, intent(in)                        :: contr_analytical 
character(len=512), intent(in)             :: control_t_dep_signal_file
character(len=512), intent(in)             :: analytical_expression
integer, intent(in)                        :: analytical_len
real*8, intent(in)                         :: analytical_tmax
logical, intent(in)                        :: controllerhasbeencalledbefore
real*8, intent(in)                         :: t_norm
real*8, intent(out)                        :: input_signal

! Allocating parameters necessary for calculating the analytical expression using python
integer                                    :: el, err
character(len=60)                          :: es, filename
real*8                                     :: er

if (contr_selfdefined .and. (contr_usedatafile .or. contr_analytical)) then
    if (sim%my_id .eq. 0) write(*,*) "ERROR: You cannot define the input signal in (more than) two ways at once - 1."
    stop
endif
  
if (contr_usedatafile .and. contr_analytical) then
    if (sim%my_id .eq. 0) write(*,*) "ERROR: You cannot define the input signal in (more than) two ways at once - 2."
    stop
endif

!> Example of how to define a certain time dependent fueling rate signal from within the controller function. 
!> Currently it is the same function as the time-dependent_puff function within mod_particle_puffing (but copied and renamed below)
if (contr_selfdefined) then
    input_signal = time_dependent_puff_controller(40.d21, sim%time, 10*t_norm, 500*t_norm, 20.d21) ! < change the parameters here
    if (sim%my_id .eq. 0) write(*,*) sim%time, input_signal, "#sim time & input signal selfdefined function"

    !> Example of how to import a time dependent signal using a datafile
else if (contr_usedatafile) then
    if (.not. controllerhasbeencalledbefore) then ! first time controller is called --> read datafile 
        call readProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, control_t_dep_signal_file)
        if (sim%my_id .eq. 0) write(*,*) "During the first timestep in the controller the datafile is imported."
    endif 
    
    ! First make sure simulation will run with a constant final value when the data file time is shorter than the simulation time
    if (sim%time .ge.t_dep_signal_controller%time(t_dep_signal_controller%len)) then 
        input_signal = t_dep_signal_controller%signal(t_dep_signal_controller%len)
        if (sim%my_id .eq. 0) write(*,*) "The simulation time is larger than the final datafile time. The final value in the datafile is kept as constant."
    else !Interpolate the data when necessary
        input_signal = interpolProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, sim%time)
        if (sim%my_id .eq. 0) write(*,*) "controller input signal data interpolation done."
    endif
        
    if (sim%my_id .eq. 0) write(*,*) sim%time, input_signal, "#sim time & input_signal using datafile"

!> Example of how to create and use a datafile for a time dependent signal using Python
else if (contr_analytical) then
    if (.not. controllerhasbeencalledbefore) then ! first time controller is called --> produce datafile using analytical expression.
        ! --- Python script - copied from vacuum.f90 and adjusted 
        call random_seed()
        err = 1
        do while ( err /= 0 )
            call random_number(er)
            el = er * 99999999
            write(es,*) el
            filename='./jorek_controller_expr_'//trim(adjustl(es))//'.py'
            open(42, file=trim(filename), status='new', iostat=err)
        end do
        111 format(2a)
        112 format(a,i16)
        113 format(a,es25.16) ! check of es hier es moet zijn of s
        !The following lines are the Python code
        write(42,111) 'from math import *'
        write(42,111) 'import numpy as np'
        write(42,111) 'def f(t):'
        write(42,111) '  return ', trim(analytical_expression)
        write(42,112) 'len=', analytical_len
        write(42,113) 'tmin=', 0.d0
        write(42,113) 'tmax=', analytical_tmax
        write(42,111) 'for x in range(1,len):'
        write(42,111) '  t=tmin+(x-1)/float(len-1)*(tmax-tmin)'
        write(42,111) '  s = "%25.16e"%t'
        write(42,111) '  s += "%25.16e"%f(t)'
        write(42,111) '  print(s)'
        close(42)

        ! --- Call Python
        call system('python ./jorek_controller_expr_'//trim(adjustl(es))//'.py > ./jorek_controller_expr_'//trim(adjustl(es))//'.dat')

        ! --- Read the result
        call readProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, &
                        t_dep_signal_controller%len, './jorek_controller_expr_'//trim(adjustl(es))//'.dat')
        if (sim%my_id .eq. 0) write(*,*) "During the first timestep in the controller the datafile is made using Python and imported."
        ! --- Delete temporary files ! not used here because it is useful to check the generated datafile
        !call system('rm ./jorek_controller_expr_'//trim(adjustl(es))//'.py ./jorek_controller_expr_'//trim(adjustl(es))//'.dat')
    endif ! first time controller is called --> produce datafile using analytical expression.

    ! First make sure simulation will run with a constant final value when the data file time is shorter than the simulation time
    if (sim%time .ge.t_dep_signal_controller%time(t_dep_signal_controller%len)) then
        input_signal = t_dep_signal_controller%signal(t_dep_signal_controller%len)
        if (sim%my_id .eq. 0) write(*,*) "The simulation time is larger than the final datafile time. The final value in the datafile is kept as constant."
    else !Interpolate the data when necessary
        input_signal = interpolProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, sim%time)
        if (sim%my_id .eq. 0) write(*,*) "controller input signal data interpolation done."
    endif      
    
    if (sim%my_id .eq. 0) write(*,*) sim%time, input_signal, "#sim time & input_signal using analytical expression"  

else
    if (sim%my_id .eq. 0) write(*,*) "ERROR: when you use the controller you need to specify what the input signal is. One of the logicals must be true."
    stop
endif ! this if statement determines the input signal: selfdefined / usedatafile / analytical
end subroutine controller_input

subroutine closed_loop(sim, element_list, node_list, previous_time_controller, setpoint, controller_K_p, controller_K_i, controller_K_d, controller_output, &
                    controller_error_prev) 

implicit none 
    
! Allocating the function arguments
type(particle_sim), intent(inout)          :: sim
type(type_node_list), intent(in)           :: node_list
type (type_element_list), intent(in)       :: element_list
real*8, intent(inout)                      :: previous_time_controller
real*8, intent(in)                         :: controller_K_p
real*8, intent(in)                         :: controller_K_i
real*8, intent(in)                         :: controller_K_d
real*8, intent(in)                         :: setpoint
real*8, intent(out)                        :: controller_output
real*8, intent(inout)                      :: controller_error_prev

! Allocating parameters necessary for using function find_RZ()
integer                                    :: i
real*8                                     :: R_find, Z_find
real*8                                     :: R_out, Z_out, s_out, t_out
integer                                    :: ifail

! Allocating parameters necessary for using function interp_PRZ()
integer                                    :: i_elm_out, i_elm, n_v
integer                                    :: i_v(1) !interp_PRZ is used for 1 variable here, that is why the length of i_v(1) is 1
real*8                                     :: time, phi, s, t
real*8                                     :: density_controller(1), P(1), P_s(1), P_t(1), P_time(1) !interp_PRZ is used for 1 variable here, that is why the length is 1
real*8                                     :: R, R_s, R_t, Z, Z_s, Z_t
real*8                                     :: P_phi(1) !interp_PRZ is used for 1 variable here, that is why the length of P_phi(1) is 1

! closed loop controller parameters
real*8                                     :: measured_value
real*8                                     :: controller_error
real*8                                     :: controller_tstep
real*8                                     :: controller_P
real*8                                     :: controller_I = 0.d0
real*8                                     :: controller_D

!> The part below calculates the controller time step
controller_tstep = sim%time - previous_time_controller ! calculate the timestep in which the controller is used
if (sim%my_id .eq. 0) write(*,*) "this is controller_tstep:", controller_tstep

!> The part below is the measurement of the controlled parameter
real*8 :: measurement_array_R(:)
real*8 :: measurement_array_Z(:)
measurement_array_R = 4.173d0 !make sure input for R_find and Z_find is given as a real!
measurement_array_Z = -3.0d0, -3.5d0,-4.0d0 !make sure input for R_find and Z_find is given as a real!
if (sim%my_id .eq. 0) write(*,*) "this is measurement_array_R", measurement_array_R
if (sim%my_id .eq. 0) write(*,*) "this is measurement_array_Z", measurement_array_Z
real*8 :: measurements(size(measurement_array_R),size(measurement_array_Z))
do i=1,size(measurement_array_R)
    do j=1,size(measurement_array_Z)
        if (sim%my_id .eq. 0) write(*,*) "this is i", i, "this is j", j
        if (sim%my_id .eq. 0) write(*,*) "this is R_in", measurement_array_R(i), "this is Z_in" measurement_array_Z(j)
        call find_RZ(sim%fields%node_list,sim%fields%element_list,measurement_array_R(i),measurement_array_Z(j),R_out,Z_out,i_elm_out,s_out,t_out,ifail))
        if (sim%my_id .eq. 0) write(*,*) "this is R and Z after call findRZ:", R_out, Z_out
        if (sim%my_id .eq. 0) write(*,*) "this is ielm_out, t_out, s_out:", i_elm_out, s_out, t_out
        if ((i_elm_out .le. 0) .and. (sim%my_id .eq. 0)) write(*,*) "WARNING: i_elm_out < 0"
        call sim%fields%interp_PRZ(sim%time, i_elm_out, [5],1, s_out, t_out, phi, density_controller, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t) !old version measurement: node_list%node(541)%values(1,1,var_rho) ! note the node number
        if (sim%my_id .eq. 0) write(*,*) "this is R and Z after call interp_PRZ:", R, Z
        if (sim%my_id .eq. 0) write(*,*) "this is the measured value after interp_PRZ:", density_controller
        if (sim%my_id .eq. 0) write(*,*) sim%time, density_controller, "#time & measured value"
        measurements(i,j) = density_controller
        if (sim%my_id .eq. 0) write(*,*) "this is i, j, measurement", i, j, measurements(i,j)
    enddo
enddo
if (sim%my_id .eq. 0) write(*,*) "this is the measurement array", measurements
heatflux_measurement = max(measurements) !change this when no or another condition on the measurement must hold
if (sim%my_id .eq. 0) write(*,*) "this is the max heatflux", heatflux_measurement

stop !for testing
!> The part below specifies the closed loop PID controller
if (sim%my_id .eq. 0) write(*,*) "start closed loop controller"
measured_value = density_controller(1) 
if (sim%my_id .eq. 0) write(*,"(E16.8,A25)") setpoint, "this is the setpoint" 
controller_error = (setpoint - measured_value)*central_density*1.d20 
if (sim%my_id .eq. 0) write(*,"(E16.8,E16.8,A52)") sim%time, controller_error, "this is the error between setpoint and measurement"
controller_P = controller_error  
controller_I = controller_I + controller_error*controller_tstep 
if (sim%my_id .eq. 0) write(*,"(E16.8,E16.8,A40)") sim%time, controller_I, "this is the integral over the error"
if (controller_error_prev .eq. 0.d0) then
    controller_D = 0
else
    controller_D = (controller_error - controller_error_prev)/controller_tstep 
endif
if (sim%my_id .eq. 0) write(*,"(E16.8,E16.8,A40)") sim%time, controller_D, "this is the derivative of the error"
controller_output = controller_K_p*controller_P + controller_K_i*controller_I + controller_K_d*controller_D 
if (sim%my_id .eq. 0) write(*,"(E16.8,E16.8,A35)") sim%time, controller_output, "this is the controller output"
    
!> The part below saves the values that are necessary in the next timestep
controller_error_prev = controller_error  
if (sim%my_id .eq. 0) write(*,*) "this is the saved controller_error_prev:", controller_error_prev
previous_time_controller = sim%time
if (sim%my_id .eq. 0) write(*,*) "this is the saved previous time:", previous_time_controller
end subroutine closed_loop

!This is the same function as defined in mod_particle_puffing. it can be adjusted to contain any other time dependent function
pure function time_dependent_puff_controller(max_puff,time, t_puff_start,t_puff_slope, min_puff) result(to_puff)
real*8,intent(in)   :: max_puff, min_puff
real*8              :: to_puff
real*8,intent(in)   :: t_puff_start,t_puff_slope
real*8,intent(in)   :: time

if (time-(t_puff_start+t_puff_slope) .ge. 0.d0) then
	to_puff = max_puff
elseif (time-t_puff_start .ge. 0.d0) then
	to_puff = min_puff+ (max_puff -min_puff) * (time-t_puff_start)/(t_puff_slope)  
else
    to_puff = min_puff 
endif
end function time_dependent_puff_controller

end module