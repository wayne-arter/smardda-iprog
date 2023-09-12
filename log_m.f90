module log_m

  use const_kind_m
  use date_time_m
  use log_h
  use smitermpi_h !> Needed for global rank
    

  implicit none

  private
  
  ! public subroutines
  public :: &
       log_init,   & !< create log file and initialise
       log_error,  & !< add error message to log
       log_value,  & !< add name = value to log
       log_alloc_check,  & !< check array allocation
       log_open_check,  & !< check file opening
       log_read_check,  & !< check result of read
       log_write_check,  & !< check result of write
       log_getunit,  & !< get unit number for logging (namelists)
       log_flusherror,  & !< flushes the error log to file
       log_close     !< close log

  interface log_value
     module procedure log_value_ki2
     module procedure log_value_ki4
     module procedure log_value_kr4
     module procedure log_value_kr8
     module procedure log_value_char
     module procedure log_value_kl
  end interface log_value

  ! public variables
  integer, public, parameter :: error_fatal=0 !< public variable
  integer, public, parameter :: error_serious=1 !< public variable
  integer, public, parameter :: error_warning=2 !< public variable
  integer, public, parameter :: error_info=3 !< public variable
  integer, public, parameter :: error_debug=4 !< public variable
  integer, public, parameter :: log_info=5 !< public variable


  ! private variables
  type(date_time_t) :: timedate !< timestamp
  character(*), parameter :: m_name='log_m' !< module name
  character(*), dimension(0:5), parameter :: errorname = & !< type of "error"
       &(/'Fatal  ','Serious','Warning', 'Info   ', 'Debug  ','Info   '/)
  character(128) :: logfile !< name of error logging file
  integer :: nlog    !< error logging unit number
  integer(ki4) :: i       !< loop counter
  integer(ki4) :: errorno       !<  number of error messages
  integer(ki4) :: seriouserrors !< number of serious errors
  integer(ki4),parameter :: maxseriouserrors=10 !< max no of serious errors
  integer(ki4),parameter :: ebufinitsize=32 !< initial size of the error buffer
  integer(ki4),parameter :: ebufmaxsize=100000 !< maximum size of the error buffer (about 7MB)

  type(error_buffer_t) :: error_buffer
  
contains
  !---------------------------------------------------------------------
  !> create log file and initialise
  subroutine log_init(fileroot,timestamp,decomp)

    !! arguments
    character(*), intent(inout) :: fileroot !<root of log file name
    type(date_time_t), intent(in) :: timestamp !< timestamp
    type(decomp_t), intent(in), optional :: decomp !< decomp data structure
    !! local
    logical :: unitused !< flag to test unit is available
    integer(ki4) :: ilen  !< length of string
    integer(ki4) :: ierr !< needed for MPI Barrier if requred

    ! Allocate error buffer with a default size 
    ! this can grow dynamically
    if( .not. allocated(error_buffer%error_line)) then
       allocate(error_buffer%error_line(ebufinitsize))
    endif
    error_buffer%size = ebufinitsize
    
    !! initialise counters
    errorno=0
    seriouserrors=0
    error_buffer%filled = 0
    error_buffer%errorno = 0

    !! create filename
    timedate=timestamp
    ! strip ctl from end of name
    ilen=len_trim(fileroot)
    if (ilen>4) then
       if (fileroot(ilen-3:ilen)=='.ctl') ilen=ilen-4
    end if
    logfile=fileroot(1:ilen)//'.log'
    ! and reset file root
    fileroot=logfile(1:ilen)

    !! get file unit
    do i=99,1,-1
       inquire(i,opened=unitused)
       if(.not.unitused)then
          nlog=i
          exit
       end if
    end do
    !! open file
    if(myrank_log == 0) open(unit=nlog,file=logfile,status='REPLACE')
#ifdef WITH_MPI
    if(present(decomp)) then ! If the decomnp structure is included then this is a call from powcal
       call MPI_BARRIER(decomp%commpowcal,ierr) ! Need to wait until the root process creates the log file       
       if(myrank_log > 0) open(unit=nlog,file=logfile,status='OLD')
    endif
#endif
    

    !! write initial header

    write(nlog,'(/,a,a,/)') 'LOG FILE: ',trim(logfile)
    write(nlog, '(a,/)') timedate%long


  end subroutine log_init
  !---------------------------------------------------------------------
  !>write message to error and info log
  subroutine log_error(modname,subname,point,severity,errormessage,writeout)

    !! arguments
    character(*), intent(in) :: modname  !< module name
    character(*), intent(in) :: subname  !< subprogram name
    integer(ki4), intent(in) :: point    !< calling point
    integer(ki4), intent(in) :: severity !< error severity
    character(*), intent(in) :: errormessage  !< error message
    logical, intent(in), optional :: writeout

    !! local
    logical :: flush = .True. ! errors default to write to screen, change to false to write to a buffer that can be periodically flushed

    !! If there is only 1 process then we don't want to store this in the buffer but flush directly
    if(nproc_log == 1) flush = .True.

    !! Or it can be forced to flush manually
    if(present(writeout)) flush = writeout

    if( flush .eqv. .False. ) then
       call add_error_line(modname,subname,point,severity,errormessage)          
    else if(severity<log_info) then       
       errorno=errorno+1   !<  This errorno is only updated during a flush, errorno is also stored independently in the error_buffer
       write(nlog, '(i7.7,a,i2,a)') &
            &   errorno, ': '//errorname(severity)//': '// trim(modname)//':'//trim(subname)//':', point, ':- '//trim(errormessage)
    else
       write(nlog, '("Log  ",a,i2,a)') &
            &   ': '//errorname(severity)//': '// trim(modname)//':'//trim(subname)//':', point, ':- '//trim(errormessage)
    end if

    !! fatal error
    if(severity==error_fatal) then
       !!close down
       seriouserrors=seriouserrors+1
       write(nlog,'(a)') 'FATAL ERROR - run terminated'
       call log_close
       stop 1

       !!serious error
    else if(severity==error_serious) then
       seriouserrors=seriouserrors+1
       if(seriouserrors>=maxseriouserrors)then
          !!close down
          write(nlog,'(a)') 'TOO MANY SERIOUS ERRORS - run terminated'
          call log_close
          stop 1

       end if


       !!error warning
    else

       !!error info

       !!error debug

       !! log_info

    end if

  end subroutine log_error
  !---------------------------------------------------------------------
  !> check for errors in array allocation
  subroutine log_alloc_check(modname,subname,point,status)

    ! arguments
    character(len=*), intent(in) :: modname  !< module name
    character(len=*), intent(in) :: subname  !< subprogram name
    integer, intent(in) :: point    !< calling point
    integer, intent(in) :: status   !< error status flag

    if(status/=0)then
       call log_error(modname,subname,point,error_fatal,'allocation failed')
    end if


  end subroutine log_alloc_check
  !---------------------------------------------------------------------
  !> check for errors in file opening
  subroutine log_open_check(modname,subname,point,status)

    ! arguments
    character(len=*), intent(in) :: modname  !< module name
    character(len=*), intent(in) :: subname  !< subprogram name
    integer, intent(in) :: point    !< calling point
    integer, intent(in) :: status   !< error status flag

    if(status/=0)then
       call log_error(modname,subname,point,error_fatal,'file opening failed')
    end if


  end subroutine log_open_check
  !---------------------------------------------------------------------
  !> check for errors in read
  subroutine log_read_check(modname,subname,point,status)

    ! arguments
    character(len=*), intent(in) :: modname  !< module name
    character(len=*), intent(in) :: subname  !< subprogram name
    integer, intent(in) :: point    !< calling point
    integer, intent(in) :: status   !< error status flag

    if(status/=0)then
       call log_error(modname,subname,point,error_fatal,'read failed')
    end if


  end subroutine log_read_check
  !---------------------------------------------------------------------
  !> check for errors in write
  subroutine log_write_check(modname,subname,point,status)

    ! arguments
    character(len=*), intent(in) :: modname  !< module name
    character(len=*), intent(in) :: subname  !< subprogram name
    integer, intent(in) :: point    !< calling point
    integer, intent(in) :: status   !< error status flag

    if(status/=0)then
       call log_error(modname,subname,point,error_fatal,'write failed')
    end if


  end subroutine log_write_check
  !---------------------------------------------------------------------
  !> get unit number for logging (namelists)
  subroutine log_getunit(kunit)

    !! arguments
    integer(ki4), intent(out) :: kunit    !< log unit number

    kunit=nlog

  end subroutine log_getunit
  !---------------------------------------------------------------------
  !> close log file
  subroutine log_close

    if(myrank_log .eq. 0) then
       ! Flush the error buffer first
       call log_flusherror
    
       !! arguments
       
       write(nlog, '(//,a,/,a)') ' Error Summary ',' ------------- '
       write(nlog, '(a,i7)') ' total number of errors   = ',errorno
       write(nlog, '(a,i7)') ' number of serious errors = ',seriouserrors
       write(nlog,'(a,/)') date_time_long()
       write(nlog,'(a,/)') 'END OF LOG FILE'
       close(unit=nlog)
    endif
    
    deallocate(error_buffer%error_line)
    error_buffer%size = 0
    error_buffer%filled = 0
    error_buffer%errorno = 0

  end subroutine log_close
  !---------------------------------------------------------------------
  !> log name=value
  subroutine log_value_ki2(varname,value,units)
    !! arguments
    character(*), intent(in) :: varname  !< variable name
    integer(ki2par), intent(in) :: value    !< variable value
    character(*), intent(in),optional :: units  !< units name
    
    !!output
    if(present(units)) then
       write(nlog, '("Log  : ",a," = ",i10,1x,a)') varname,value,units
    else
       write(nlog, '("Log  : ",a," = ",i10)') varname,value
    end if
  end subroutine log_value_ki2
  subroutine log_value_ki4(varname,value,units)
    !! arguments
    character(*), intent(in) :: varname  !< variable name
    integer(ki4), intent(in) :: value    !< variable value
    character(*), intent(in),optional :: units  !< units name

    if(myrank_log .gt. 0) return
    
    !!output
    if(present(units)) then
       write(nlog, '("Log  : ",a," = ",i10,1x,a)') varname,value,units
    else
       write(nlog, '("Log  : ",a," = ",i10)') varname,value
    end if

  end subroutine log_value_ki4
  
  subroutine log_value_kr4(varname,value,units)
    !! arguments
    character(*), intent(in) :: varname  !< variable name
    real(ki4), intent(in) :: value    !< variable value
    character(*), intent(in), optional :: units  !< units name

    if(myrank_log .gt. 0) return
        
    !!output
    if(present(units)) then
       write(nlog, '("Log  : ",a32,"  = ",g12.5,1x,a)') varname,value,units
    else
       write(nlog, '("Log  : ",a," = ",g12.5)') varname,value
    end if
  end subroutine log_value_kr4
  
  subroutine log_value_kr8(varname,value,units)
    !! arguments
    character(*), intent(in) :: varname  !< variable name
    real(ki8), intent(in) :: value    !< variable value
    character(*), intent(in),optional :: units  !< units name

    if(myrank_log .gt. 0) return
    
    !!output

    if(present(units)) then
       write(nlog, '("Log  : ",a32,"  = ",1pe15.8,2x,a)') varname,value,units
    else
       write(nlog, '("Log  : ",a," = ",g20.8)') varname,value
    end if
  end subroutine log_value_kr8
  
  subroutine log_value_kl(varname,value)
    !! arguments
    character(*), intent(in) :: varname  !< variable name
    logical, intent(in) :: value    !< variable value

    if(myrank_log .gt. 0) return
    
    !!output
    write(nlog, '("Log  : ",a," = ",l2)') varname,value
  end subroutine log_value_kl
  
  subroutine log_value_char(varname,value)
    !! arguments
    character(*), intent(in) :: varname  !< variable name
    character(*), intent(in) :: value    !< variable value

    if(myrank_log .gt. 0) return

    !!output
    write(nlog, '("Log  : ",a," = ",a)') varname,value
  end subroutine log_value_char

!-----------------------------------------------------------------------------------------
!> Routines for the error_buffer
!> Adds an error to the error_buffer
  subroutine add_error_line(modname,subname,point,severity,errormessage)

    character(*), intent(in) :: modname  !< module name
    character(*), intent(in) :: subname  !< subprogram name
    integer(ki4), intent(in) :: point    !< calling point
    integer(ki4), intent(in) :: severity !< error severity
    character(*), intent(in) :: errormessage  !< error message

    ! Check if the error buffer is full
    if(error_buffer%filled >= error_buffer%size) then
       call grow_error_buffer
    endif

    error_buffer%filled = error_buffer%filled + 1
    error_buffer%errorno = error_buffer%errorno + 1 ! This is cumulative and not flushed

    error_buffer%error_line(error_buffer%filled)%modname = modname
    error_buffer%error_line(error_buffer%filled)%subname = subname
    error_buffer%error_line(error_buffer%filled)%point = point
    error_buffer%error_line(error_buffer%filled)%severity = severity
    error_buffer%error_line(error_buffer%filled)%message = errormessage
    
  end subroutine add_error_line


  ! Doubles the size of the error buffer if it is full
  subroutine grow_error_buffer

    type(error_line_t),dimension(:), allocatable :: new_error_line
    integer(ki4) :: newsize

    newsize = 2 * error_buffer%size

    ! First check if the buffer is too big and flush and return if it is
    if(newsize .gt. ebufmaxsize) then
       call log_flusherror()
       return
    endif
    
    allocate(new_error_line(newsize))

    new_error_line(1:error_buffer%size) = error_buffer%error_line
    deallocate(error_buffer%error_line)
    call move_alloc(new_error_line,error_buffer%error_line)

    error_buffer%size = newsize

  end subroutine grow_error_buffer

! Flushes the current error buffer to file
  subroutine log_flusherror

    integer(ki4) :: i

    do i=1,error_buffer%filled
       call log_error(error_buffer%error_line(i)%modname,error_buffer%error_line(i)%subname,error_buffer%error_line(i)%point, &
            error_buffer%error_line(i)%severity,error_buffer%error_line(i)%message,.True.)
    end do

    error_buffer%filled = 0
    
  end subroutine log_flusherror  
  
end module log_m
