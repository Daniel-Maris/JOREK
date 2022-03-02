!> Module for controlling a parameter in the plasma
!> As an example, this module will control the time-varying gas puff fueling rate at the locations where a gas valve is present
!> But it is aimed at being usable for various other applications

!note: all variables that you allocate here cannot be allocated at the same time in the module where you use mod_controller

module mod_controller
    
    use mod_particle_puffing
    use mod_particle_sim
    use phys_module, only: index_now

    implicit none
        
    ! use-case specific parameters
    logical                          :: puff_t_dependent
    type(particle_puffing)           :: gas_puff
    class(particle_puffing), pointer :: this
    real*8                           :: t_norm

    ! parameters used for interpolation of dataranges
    integer :: left, right, mid
    real*8  :: aux1, aux2, currenttime

    ! parameters necessary for calculating the analytical expression using python
    integer :: el, err
    character(len=60) :: s, filename
    real*8 :: r

    !defining a type for a general time dependent signal
    type :: time_dependent_signal
    integer             :: len        !< Number of points in numerical time trace.
    real*8, allocatable :: time(:)    !< time-values of numerical time trace
    real*8, allocatable :: signal(:)  !< functionvalues of numerical time trace
    end type time_dependent_signal

contains

subroutine controller_function(use_controller, this, sim, t_dep_signal_controller, contr_change_t_indep, &
                                contr_change_t_dep, contr_selfdefined, contr_usedatafile, contr_analytical, &
                                control_t_dep_signal_file, analytical_expression)
    use profiles, only: readProf
    
    implicit none 
    
    class(particle_puffing) , intent(inout)    :: this
    type(particle_sim), intent(inout)          :: sim
    type(time_dependent_signal), intent(inout) :: t_dep_signal_controller
    logical,intent(in)                         :: use_controller
    logical, intent(in)                        :: contr_change_t_indep
    logical, intent(in)                        :: contr_change_t_dep 
    logical, intent(in)                        :: contr_selfdefined 
    logical, intent(in)                        :: contr_usedatafile 
    logical, intent(in)                        :: contr_analytical 
    character(len=512),intent(in)              :: control_t_dep_signal_file
    character(len=512),intent(in)              :: analytical_expression

    !> The following if-statement determines what the controller does
    !> note: if you want to change parameters that are part of an action within an event, the event must be a pointer (use new_event_ptr() instead of event())
    if (use_controller) then
        write(*,*) "The controller_function works. The controller is on"
        write(*,"(A,g12.4)") "test voor controller, this is the time now", sim%time 
        write(*,"(A,g12.4)") "test voor controller, this is index_now", index_now 

        if (contr_change_t_indep .and. (contr_change_t_dep .or. contr_selfdefined .or. contr_usedatafile .or. contr_analytical)) then
            write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 1."
            stop
        endif
          
        if (contr_change_t_dep .and. (contr_selfdefined .or. contr_usedatafile .or. contr_analytical)) then
            write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 2."
            stop
        endif
          
        if (contr_selfdefined .and. (contr_usedatafile .or. contr_analytical)) then
            write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 3."
            stop
        endif
          
        if (contr_usedatafile .and. contr_analytical) then
            write(*,*) "ERROR: You cannot use the controller for changing the signal in (more then) two ways at once - 4."
            stop
        endif

        if (contr_change_t_dep .and. .not. puff_t_dependent) then
            write(*,*) "ERROR: If you want to change the predefined time dependent signal values using the controller, puff_t_dependent must be true."
            stop
        endif 

        !> Example on how to change the fueling rate of a time independent puff
        if (contr_change_t_indep) then
            gas_puff%fueling_rate = 40.d21 

        !> Example how to change the defined !!time-dependent!! signal  
        else if (contr_change_t_dep .and. puff_t_dependent) then                     
            gas_puff%fueling_rate = 60.d21    
            gas_puff%t_puff_slope = 500*t_norm
            gas_puff%t_puff_start = 10*t_norm
        
        !> Example of how to define a certain time dependent fuelling rate signal from within the controller function. 
        !> Currently it is the same function as the time-dependent_puff function within mod_particle_puffing (but copied and renamed below)
        !> To use it, puff_t_dependent must be set .false. in my_example.
        else if (contr_selfdefined) then
            gas_puff%fueling_rate = time_dependent_puff_controller(40.d21, sim%time, 10*t_norm, 500*t_norm, 20.d21) ! < change the parameters here

        !> Example of how to import a time dependent signal using a datafile
        else if (contr_usedatafile) then
            if (index_now .eq. 2) then
                call readProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, t_dep_signal_controller%len, control_t_dep_signal_file)
                write(*,*) "During the first timestep in the controller the datafile is imported."
            endif !index_now = 2 --> read datafile during first timestep in which the controller is called
            
            ! First make sure simulation will run with a constant final value when the data file time is shorter than the simulation time
            if (sim%time .ge.t_dep_signal_controller%time(t_dep_signal_controller%len)) then 
                gas_puff%fueling_rate = t_dep_signal_controller%signal(t_dep_signal_controller%len)
                write(*,*) "The simulation time is larger than the final datafile time. The final value in the datafile is kept as constant."
            else !Interpolate the data when necessary
                left  = 1
                right = t_dep_signal_controller%len
                currenttime = sim%time
                do !Search for the two datapoints in the datafile where in between the current time lies, in order for linear interpolation
                    if ( right == left + 1 ) exit
                    mid = (left + right) / 2
                    if ( t_dep_signal_controller%time(mid) >= currenttime ) then
                      right = mid
                    else
                    left = mid
                    end if
                end do
                aux1 = (currenttime - t_dep_signal_controller%time(left)) / (t_dep_signal_controller%time(right) - t_dep_signal_controller%time(left))
                aux2 = (1. - aux1)
                gas_puff%fueling_rate = t_dep_signal_controller%signal(left) * aux2 + t_dep_signal_controller%signal(right) * aux1
                write(*,"(A,g12.4,A,g12.4)") "previous time", t_dep_signal_controller%time(left), "previous signal", t_dep_signal_controller%signal(left)
                write(*,"(A,g12.4,A,g12.4)") "interpolation done. current time", sim%time, "current signal", gas_puff%fueling_rate
                write(*,"(A,g12.4,A,g12.4)") "next time", t_dep_signal_controller%time(right), "next signal", t_dep_signal_controller%signal(right)
            endif !use final value if simtime> datafiletime, or use interpolation
                
        !> Example of how to create and use a datafile for a time dependent signal using Python
        else if (contr_analytical) then
            if (index_now .eq. 2) then
                write(*,*) "During the first timestep in the controller the datafile is made using Python."
                ! --- Python script - copied from vacuum.f90 and adjusted 
                call random_seed()
                err = 1
                do while ( err /= 0 )
                    call random_number(r)
                    el = r * 99999999
                    write(s,*) el
                    filename='./jorek_controller_expr_'//trim(adjustl(s))//'.py'
                    open(42, file=trim(filename), status='new', iostat=err)
                end do
                111 format(2a)
                112 format(a,i16)
                113 format(a,es25.16)
                !The following lines are the Python code
                write(42,111) 'from math import *'
                write(42,111) 'def f(t):'
                write(42,111) '  return ', trim(analytical_expression) !old expression for tests: '1200*exp(-(t-1000.)**2/(200.)**2)')
                write(42,112) 'len=', 144
                write(42,113) 'tmin=', 0
                write(42,113) 'tmax=', 2.d-4
                write(42,111) 'for x in range(1,len):'
                write(42,111) '  t=tmin+(x-1)/float(len-1)*(tmax-tmin)'
                write(42,111) '  s = "%25.16e"%t'
                write(42,111) '  s += "%25.16e"%f(t)'
                write(42,111) '  print(s)'
                close(42)
        
                ! --- Call Python
                call system('python ./jorek_controller_expr_'//trim(adjustl(s))//'.py > ./jorek_controller_expr_'//trim(adjustl(s))//'.dat')
        
                ! --- Read the result
                call readProf(t_dep_signal_controller%time, t_dep_signal_controller%signal, &
                                t_dep_signal_controller%len, './jorek_controller_expr_'//trim(adjustl(s))//'.dat')
        
                ! --- Delete temporary files
                call system('rm ./jorek_controller_expr_'//trim(adjustl(s))//'.py ./jorek_controller_expr_'//trim(adjustl(s))//'.dat')
            endif ! index_now = 2 --> produce datafile using analytical expression.

            ! First make sure simulation will run with a constant final value when the data file time is shorter than the simulation time
            if (sim%time .ge.t_dep_signal_controller%time(t_dep_signal_controller%len)) then
                gas_puff%fueling_rate = t_dep_signal_controller%signal(t_dep_signal_controller%len)
                write(*,*) "The simulation time is larger than the final datafile time. The final value in the datafile is kept as constant."
            else !Interpolate the data when necessary
                left  = 1
                right = t_dep_signal_controller%len
                currenttime = sim%time
                do !Search for the two datapoints in the datafile where in between the current time lies, in order for linear interpolation
                    if ( right == left + 1 ) exit
                    mid = (left + right) / 2
                    if ( t_dep_signal_controller%time(mid) >= currenttime ) then
                      right = mid
                    else
                    left = mid
                    end if
                end do
                aux1 = (currenttime - t_dep_signal_controller%time(left)) / (t_dep_signal_controller%time(right) - t_dep_signal_controller%time(left))
                aux2 = (1. - aux1)
                gas_puff%fueling_rate = t_dep_signal_controller%signal(left) * aux2 + t_dep_signal_controller%signal(right) * aux1
                write(*,"(A,g12.4,A,g12.4)") "previous time", t_dep_signal_controller%time(left), "previous signal", t_dep_signal_controller%signal(left)
                write(*,"(A,g12.4,A,g12.4)") "interpolation done. current time", sim%time, "current signal", gas_puff%fueling_rate
                write(*,"(A,g12.4,A,g12.4)") "next time", t_dep_signal_controller%time(right), "next signal", t_dep_signal_controller%signal(right)
            endif ! use final value if simtime> datafiletime. otherwise: use interpolation        
        else
            write(*,*) "ERROR: when you use the controller you need to specify what happens in the controller function. One of the logicals must be true"
            stop
        endif !this if statement determines the controllerfunction: change_t_indep / change_t_dep / selfdefined / usedatafile or analytical

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