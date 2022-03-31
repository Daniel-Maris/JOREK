!> Module for controlling a parameter in the plasma
!> As an example, this module will control the time-varying gas puff fueling rate at the locations where a gas valve is present
!> But it is aimed at being usable for various other applications

!note: all variables that you allocate here cannot be allocated at the same time in the module where you use mod_controller

module mod_controller
    
    use mod_particle_puffing
    use mod_particle_sim
    use phys_module, only: index_now, central_density
    use data_structure
    use mod_fields
    
    implicit none
    
    !defining a type for a general time dependent signal
    type :: time_dependent_signal
    integer                          :: len        !< Number of points in numerical time trace.
    real*8, allocatable              :: time(:)    !< time-values of numerical time trace
    real*8, allocatable              :: signal(:)  !< functionvalues of numerical time trace
    end type time_dependent_signal

contains

subroutine controller_function(use_controller,  sim, t_dep_signal_controller, contr_change_t_indep, &
                                contr_change_t_dep, contr_selfdefined, contr_usedatafile, contr_analytical, &
                                control_t_dep_signal_file, analytical_expression, analytical_len, analytical_tmax, &
                                controllerhasbeencalledbefore, previous_time_controller, controller_K_p, controller_K_i, &
                                controller_K_d,node_list,element_list,puff_t_dependent,t_norm,gas_puff,setpoint)
    
    use profiles, only: readProf, interpolProf
    
    implicit none 
    
    ! Allocate the function arguments
    type(particle_sim), intent(inout)          :: sim
    type(time_dependent_signal), intent(inout) :: t_dep_signal_controller
    type(type_node_list), intent(in)           :: node_list
    type (type_element_list), intent(in)       :: element_list
    type(particle_puffing), intent(inout)      :: gas_puff
    logical,intent(in)                         :: use_controller
    logical, intent(in)                        :: contr_change_t_indep
    logical, intent(in)                        :: contr_change_t_dep 
    logical, intent(in)                        :: contr_selfdefined 
    logical, intent(in)                        :: contr_usedatafile 
    logical, intent(in)                        :: contr_analytical 
    character(len=512),intent(in)              :: control_t_dep_signal_file
    character(len=512),intent(in)              :: analytical_expression
    integer, intent(in)                        :: analytical_len
    real*8, intent(in)                         :: analytical_tmax
    logical, intent(inout)                     :: controllerhasbeencalledbefore
    real*8,intent(inout)                       :: previous_time_controller
    real*8, intent(in)                         :: controller_K_p
    real*8, intent(in)                         :: controller_K_i
    real*8, intent(in)                         :: controller_K_d
    real*8, intent(in)                         :: setpoint
    logical,intent(in)                         :: puff_t_dependent
    real*8, intent(in)                         :: t_norm
    
    ! parameters necessary for using function find_RZ()
    integer                                    :: i
    real*8                                     :: R_find, Z_find
    real*8                                     :: R_out,Z_out,s_out,t_out
    integer                                    :: ifail

    ! parameters necessary for using function interp_PRZ()
    integer                                    :: i_elm_out, i_elm, n_v
    integer                                    :: i_v(1)
    real*8                                     :: time, phi, s, t
    real*8                                     :: density_controller(1), P(1), P_s(1), P_t(1), P_time(1)
    real*8                                     :: R, R_s, R_t, Z, Z_s, Z_t
    real*8                                     :: P_phi(1)

    ! parameters necessary for calculating the analytical expression using python
    integer                                    :: el, err
    character(len=60)                          :: es, filename
    real*8                                     :: er

    ! closed loop controller parameters
    real*8                                     :: measured_value
    real*8                                     :: controller_error
    real*8                                     :: controller_tstep
    real*8                                     :: controller_P
    real*8                                     :: controller_I = 0.d0
    real*8                                     :: controller_D
    real*8                                     :: controller_error_prev = 0.d0
    real*8                                     :: controller_output 

    !> The following if-statement determines what the controller does when it is called
    !> note: if you want to change parameters that are part of an action within an event, the event must be a pointer (use new_event_ptr() instead of event())
    if (use_controller) then
        if (sim%my_id .eq. 0) write(*,*) "The controller_function works. The controller is on"
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is the time now:", sim%time 
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is index_now:", index_now 
        
        if (puff_t_dependent .and. .not. contr_change_t_dep) then
            if (sim%my_id .eq. 0) write(*,*) "ERROR: puff_t_dependent cannot be true if you do not want to use the controller for changing the input values of the time_dependent_puff function."
            stop
        endif 

        if (contr_change_t_indep .and. (contr_change_t_dep .or. contr_selfdefined .or. contr_usedatafile .or. contr_analytical)) then
            if (sim%my_id .eq. 0) write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 1."
            stop
        endif
          
        if (contr_change_t_dep .and. (contr_selfdefined .or. contr_usedatafile .or. contr_analytical)) then
            if (sim%my_id .eq. 0) write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 2."
            stop
        endif
          
        if (contr_selfdefined .and. (contr_usedatafile .or. contr_analytical)) then
            if (sim%my_id .eq. 0) write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 3."
            stop
        endif
          
        if (contr_usedatafile .and. contr_analytical) then
            if (sim%my_id .eq. 0) write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 4."
            stop
        endif

        if (contr_change_t_dep .and. .not. puff_t_dependent) then
            if (sim%my_id .eq. 0) write(*,*) "ERROR: If you want to change the predefined time dependent signal values using the controller, puff_t_dependent must be true."
            stop
        endif 

        !> The first part below is the open loop controller, that determines the setpoint values.

        !> Example on how to change the fueling rate of a time independent puff
        if (contr_change_t_indep) then
            gas_puff%fueling_rate = 10.d21 
        
        !> Example how to change the defined !!time-dependent!! signal  
        else if (contr_change_t_dep .and. puff_t_dependent) then                     
            gas_puff%fueling_rate = 60.d21    
            gas_puff%t_puff_slope = 500*t_norm
            gas_puff%t_puff_start = 10*t_norm
        
        !> Example of how to define a certain time dependent fuelling rate signal from within the controller function. 
        !> Currently it is the same function as the time-dependent_puff function within mod_particle_puffing (but copied and renamed below)
        else if (contr_selfdefined) then
            gas_puff%fueling_rate = time_dependent_puff_controller(40.d21, sim%time, 10*t_norm, 500*t_norm, 20.d21) ! < change the parameters here

        !> Example of how to import a time dependent signal using a datafile
        else if (contr_usedatafile) then
            if (.not. controllerhasbeencalledbefore) then ! first time controller is called --> read datafile 
                call readProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, control_t_dep_signal_file)
                if (sim%my_id .eq. 0) write(*,*) "During the first timestep in the controller the datafile is imported."
            endif 
            
            ! First make sure simulation will run with a constant final value when the data file time is shorter than the simulation time
            if (sim%time .ge.t_dep_signal_controller%time(t_dep_signal_controller%len)) then 
                gas_puff%fueling_rate = t_dep_signal_controller%signal(t_dep_signal_controller%len)
                if (sim%my_id .eq. 0) write(*,*) "The simulation time is larger than the final datafile time. The final value in the datafile is kept as constant."
            else !Interpolate the data when necessary
                gas_puff%fueling_rate = interpolProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, sim%time)
                if (sim%my_id .eq. 0) write(*,*) "controller interpolation done. current time", sim%time, "current signal", gas_puff%fueling_rate
            endif
                
            if (sim%my_id .eq. 0) write(*,*) sim%time, gas_puff%fueling_rate, "#time & fuelingratesignal"

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
                if (sim%my_id .eq. 0) write(*,*) "During the first timestep in the controller the datafile is made using Python."
                ! --- Delete temporary files ! not used here because it is useful to check the generated datafile
                !call system('rm ./jorek_controller_expr_'//trim(adjustl(es))//'.py ./jorek_controller_expr_'//trim(adjustl(es))//'.dat')
            endif ! first time controller is called --> produce datafile using analytical expression.

            ! First make sure simulation will run with a constant final value when the data file time is shorter than the simulation time
            if (sim%time .ge.t_dep_signal_controller%time(t_dep_signal_controller%len)) then
                gas_puff%fueling_rate = t_dep_signal_controller%signal(t_dep_signal_controller%len)
                if (sim%my_id .eq. 0) write(*,*) "The simulation time is larger than the final datafile time. The final value in the datafile is kept as constant."
            else !Interpolate the data when necessary
                gas_puff%fueling_rate = interpolProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, sim%time)
                if (sim%my_id .eq. 0) write(*,*) "controller interpolation done. current time", sim%time, "current signal", gas_puff%fueling_rate
            endif        
        else
            if (sim%my_id .eq. 0) write(*,*) "ERROR: when you use the controller you need to specify what happens in the controller function. One of the logicals must be true"
            stop
        endif ! this if statement determines the controllerfunction: change_t_indep / change_t_dep / selfdefined / usedatafile or analytical

        ! test to make a list of node numbers and their respective coordinates (R,Z)
        !!! for simple plasma, not iter plasma!!!
        !if (.not. controllerhasbeencalledbefore) then
        !    open(42, file=trim('./finaldatanodelistITER.dat'), status='new')
        !    do i=1,16490 !3315=65*51 = n_radial*n_pol in my case
        !        if (sim%my_id .eq. 0) write(42,*) i , node_list%node(i)%x(1,:)
        !    end do
        !    close(42)
        !endif

        !> The part below is the closed loop controller
        controller_tstep = sim%time - previous_time_controller ! calculate the timestep in which the controller is active
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is controller_tstep:", controller_tstep

        call find_RZ(sim%fields%node_list,sim%fields%element_list,8.173d0,-0.05d0,R_out,Z_out,i_elm_out,s_out,t_out,ifail) !make sure input for R_find and Z_find is given as a real!
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is R and Z after controller called findRZ:", R_out, Z_out
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is ielm_out t_out s_out:", i_elm_out, s_out, t_out
        
        ! warning bijschrijven als je een i_elm < 0 krijgt
        
        ! USE statement mod_fields aan, dan:
        call sim%fields%interp_PRZ(sim%time, i_elm_out, [5],1, s_out, t_out, phi, density_controller, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
        
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is R and Z after interp_PRZ:", R, Z
        if (sim%my_id .eq. 0) write(*,*) "this is P(in this case the density) after interp_PRZ:", density_controller
        if (sim%my_id .eq. 0) write(*,*) sim%time, density_controller, "#time & density"
        !if (sim%my_id .eq. 0) then
        !    open(42, file=trim('./densitydata.dat'),action='write', status='old')
        !    write(42,*) sim%time, density_controller
        !    close(42)
        !endif
        if (sim%my_id .eq. 0) write(*,*) "start closed loop controller"
        measured_value = density_controller(1) !node_list%node(541)%values(1,1,var_rho) ! let op node number
        !if (sim%my_id .eq. 0) write(*,*) "test for controller, this is the density at node 541 now:", measured_value
        !if (sim%my_id .eq. 0) write(*,*) "test for controller, this is the position (R,Z,phi):", node_list%node(541)%x(1,:)
        if (sim%my_id .eq. 0) write(*,"(E14.8,A23)") setpoint, "this is the setpoint" 
        controller_error = (setpoint - measured_value)*central_density*10.d20 ! determine the error between setpoint and measured value
        if (sim%my_id .eq. 0) write(*,"(E14.8,E15.8,A40)") sim%time, controller_error, "test for controller, this is the error"
        controller_P = controller_error ! calculate proportional term  
        controller_I = controller_I + controller_error*controller_tstep ! calculate integral over error
        if (sim%my_id .eq. 0) write(*,"(E14.8,E15.8,A58)") sim%time, controller_I, "test for controller, this is integral na aanpassing"
        if (sim%my_id .eq. 0) write(*,"(E16.8,A70)") controller_error_prev, "test for controller, this is controller error previous voor D-term"
        if (controller_error_prev .eq. 0.d0) then
            controller_D = 0
        else
            controller_D = (controller_error - controller_error_prev)/controller_tstep ! calculate derivative of error
        endif
        if (sim%my_id .eq. 0) write(*,"(E14.8,E15.8,A45)") sim%time, controller_D, "test for controller, this is the derivative"
        controller_output = controller_K_p*controller_P + controller_K_i*controller_I + controller_K_d*controller_D ! calculate output controller
        if (sim%my_id .eq. 0) write(*,"(E14.8,E15.8,A55)") sim%time, controller_output, "test for controller, this is the controller output"
        gas_puff%fueling_rate = max(gas_puff%fueling_rate + controller_output,0.d0) ! calculate new output = setpoint + controller_output
        if (gas_puff%fueling_rate .eq. 0.d0) then
            if (sim%my_id .eq. 0) write(*,*) "WARNING: The output signal is set to zero because it cannot be negative"
        endif
        if (sim%my_id .eq. 0) write(*,"(E14.8,E15.8,A52)") sim%time, gas_puff%fueling_rate, "test for controller, this is the new fueling rate"
        
        controller_error_prev = controller_error  ! save error
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is the saved controller_error_prev:", controller_error_prev
        previous_time_controller = sim%time       ! save old sim%time to calculate controller_tstep
        if (sim%my_id .eq. 0) write(*,*) "test for controller, this is the saved previous time:", previous_time_controller

        controllerhasbeencalledbefore = .true. ! set this to true so that datafiles are only read once and not during every timestep the controller is called
        if (sim%my_id .eq. 0) write(*,*) "test for controller, hij heeft de controller loop afgemaakt"
    else
        if (sim%my_id .eq. 0) write(*,*) "The controller_function works. The controller is off"
    endif !(use_controller)
end subroutine controller_function

!This is the same function as defined in mod_particle_puffing. it can be adjusted to contain any other time dependent function
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
    to_puff = min_puff 
endif
end function time_dependent_puff_controller

end module