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
    ! the ones below do not have to be defined as parameters in phys-module and stuf... why??
    logical                          :: contr_change_timedepsignal 
    logical                          :: contr_selfdefinedsignal 
    logical                          :: contr_usedatafile 
    logical                          :: contr_analytical 
    logical                          :: puff_t_dependent
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

subroutine controller_function(use_controller,this,sim, time_dependent_signal_controller,contr_change_timedepsignal,contr_selfdefinedsignal,contr_usedatafile,contr_analytical)
    use profiles, only: readProf
    
    implicit none 
    
    class(particle_puffing) , intent(inout)    :: this
    type(particle_sim), intent(inout)          :: sim
    type(time_dependent_signal), intent(inout) ::time_dependent_signal_controller
    logical,intent(in)                         :: use_controller
    logical, intent(in)                        :: contr_change_timedepsignal 
    logical, intent(in)                        :: contr_selfdefinedsignal 
    logical, intent(in)                        :: contr_usedatafile 
    logical, intent(in)                        :: contr_analytical 

    integer :: l, err
    character(len=60) :: s, filename
    real*8 :: r

    controller_timedependentsignal_file = '/home/ITER/vanhooe/Documents/Datafiles/datafile.dat' ! Note: check the directory to the datafile

    !> The following if-statement determines what the controller does
    !> note: if you want to change parameters that are part of an action within an event, the event must be a pointer (use new_event_ptr() instead of event())
    if (use_controller) then
        write(*,*) "The controller_function works. The controller is on"
        !test to make sure sim%time and index_now work properly and are usable to make the controller function time dependent
        write(*,"(A,g12.4)") "test voor controller, this is the time now", sim%time 
        write(*,"(A,g12.4)") "test voor controller, this is index_now", index_now 

        !> Example how to change the defined !!time-dependent!! signal  
        if (contr_change_timedepsignal .and. puff_t_dependent) then                     
            gas_puff%fueling_rate = 60.d21    ! we can adjust this%fuelling_rate here depending on the required input signal
            gas_puff%t_puff_slope = 1.d-3
            gas_puff%t_puff_start = 10*t_norm
        
        !> Example of how to define a certain time dependent fuelling rate signal from within the controller function. 
        !> Currently it is the same function as the time-dependent_puff function within mod_particle_puffing (but copied and renamed below)
        !> To use it, puff_t_dependent must be set .false. in my_example.
        else if (contr_selfdefinedsignal) then
            gas_puff%fueling_rate = time_dependent_puff_controller(25.d21,sim%time, 10*t_norm,500*t_norm, 20.d21)
         
        !> Example of how to import a time dependent signal using a datafile
        else if (contr_usedatafile) then
            call readProf(time_dependent_signal_controller%time, time_dependent_signal_controller%signal, &
            time_dependent_signal_controller%len, controller_timedependentsignal_file)
            write(*,"(A,g12.4,A,g12.4)") "current time signal", time_dependent_signal_controller%time(index_now), "current signal", time_dependent_signal_controller%signal(index_now) 
            if (index_now .le. time_dependent_signal_controller%len) then
                gas_puff%fueling_rate = time_dependent_signal_controller%signal(index_now)
            else
                write(*,*) "The datafile is not long enough for this simulation. There are no more datapoints left. Simulation is stopped."
                stop
            endif

        !> Example of how to create and use a datafile for a time dependent signal using Python
        else if (contr_analytical) then
            if (index_now .eq. 2) then
                write(*,*) "During the first timestep in the controller the datafile is made using Python."
                ! --- Python script - copied from vacuum.f90 and adjusted afterwards
                call random_seed()
                err = 1
                do while ( err /= 0 )
                    call random_number(r)
                    l = r * 99999999
                    write(s,*) l
                    filename='./jorek_controller_expr_'//trim(adjustl(s))//'.py'
                    open(42, file=trim(filename), status='new', iostat=err)
                end do
                111 format(2a)
                112 format(a,i16)
                113 format(a,es25.16)
                !The following lines are the Python code
                write(42,111) 'from math import *'
                write(42,111) 'def f(t):'
                write(42,111) '  return ', trim('1200*exp(-(t-1000.)**2/(200.)**2)')
                write(42,112) 'len=', 144
                write(42,113) 'tmin=', 0
                write(42,113) 'tmax=', 1.d-4
                write(42,111) 'for x in range(1,len):'
                write(42,111) '  t=tmin+(x-1)/float(len-1)*(tmax-tmin)'
                write(42,111) '  s = "%25.16e"%t'
                write(42,111) '  s += "%25.16e"%f(t)'
                write(42,111) '  print(s)'
                close(42)
        
                ! --- Call Python
                call system('python ./jorek_controller_expr_'//trim(adjustl(s))//'.py > ./jorek_controller_expr_'//trim(adjustl(s))//'.dat')
        
                ! --- Read the result
                call readProf(time_dependent_signal_controller%time, time_dependent_signal_controller%signal, &
                time_dependent_signal_controller%len, './jorek_controller_expr_'//trim(adjustl(s))//'.dat')
        
                ! --- Delete temporary files
                call system('rm ./jorek_curr_expr_'//trim(adjustl(s))//'.py ./jorek_curr_expr_'//trim(adjustl(s))//'.dat')
            elseif (index_now .eq. 3) then
                write(*,*) "The first timestep in the controller has ended so the datafile is already generated and does not have to be generated again"
            endif

            !test to check the working of the controller
            write(*,"(A,g12.4,A,g12.4)") "current time signal", time_dependent_signal_controller%time(index_now), "current signal", time_dependent_signal_controller%signal(index_now) 

            if (index_now .le. time_dependent_signal_controller%len) then
                gas_puff%fueling_rate = time_dependent_signal_controller%signal(index_now)
            else
                write(*,*) "The datafile is not long enough for this simulation. There are no more datapoints left. Simulation is stopped."
                stop
            endif
        endif !this if statement determines change_timedepsignal / selfdefinedsignal / usedatafile or analytical

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