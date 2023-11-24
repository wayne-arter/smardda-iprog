program p_sumtot

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use log_h
  use log_m
  use clock_m
  use misc_m
  use gfile_m
  use vfile_m
  use icontrol_h
  use icontrol_m
  use halton_m
  use sum_h
  use sum_m
  use sumtot_h
  use sumtot_m
  use smitermpi_h

  implicit none

! Local variables
  character(*), parameter :: m_name='sumtot' !< module name
  type(iparams_t)     :: param      !< run control parameters
  type(ifiles_t)     :: file      !< names of files
  type(iplots_t)     :: plot      !< diagnostic plot selectors
  type(sumtot_t)  :: sumtot      !< sumtot object
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run

  integer(ki4):: nin=0 !< unit for other data
  integer(ki4):: nplot=0 !< unit for vtk files
  integer(ki4):: nprint=0 !< unit for gnuplot files
  integer(ki4):: i !< loop variable
  integer(ki4):: j !< loop variable

  character(len=80) :: ibuf !< character workspace
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(40)
  call clock_start(1,'run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'sumtot: integral estimates '
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! file root

!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run sumtot type at the command line:'
     print *, '   sumtot fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'control initialisation time')
  call icontrol_init(fileroot)
  call icontrol_read(file,param,sumtot%n,sumtot%sum%n,plot)
  call clock_stop(2)
!--------------------------------------------------------------------------
!! other data for sumtot and/or read sumtot controls directly

  call clock_start(15,'initialisation time')
  call sumtot_initfile(file%sumtotdata,nin)
  !call sumtot_readcon(sumtot%n)
  call clock_stop(15)
!--------------------------------------------------------------------------
!! do the main work

  call clock_start(21,'evaluation time')

  select case(param%control)
  case default
     call sumtot_solve(sumtot)
  end select

  call clock_stop(21)
!--------------------------------------------------------------------------
!! data diagnostics to log file - optional

  call clock_start(31,'diagnostics time')

  select case(param%control)
  case default
     call sumtot_dia(sumtot)
  end select

  call clock_stop(31)
!--------------------------------------------------------------------------
!! vtk plot file(s) - optional

  call clock_start(32,'vtk plot time')

  select case(param%control)
  case default
     if(plot%vtk) then
        call vfile_init(file%vtk,'file header',nplot)
        call sumtot_writev(sumtot,'selector',nplot)
        call vfile_close
     end if
  end select

  call clock_stop(32)
!--------------------------------------------------------------------------
!! gnuplot file(s) - optional

  call clock_start(33,'gnuplot time')

  select case(param%control)
  case default
     if(plot%gnu) then
        call gfile_init(trim(file%gnu),'file header',nprint)
        call sumtot_writeg(sumtot,'selector',nprint)
        call gfile_close
     end if
  end select

  call clock_stop(33)
!--------------------------------------------------------------------------
!! output file

  call clock_start(34,'outfile_init time')
  call sumtot_initwrite(fileroot)
  call sumtot_write(sumtot)
  call sumtot_closewrite()
  call clock_stop(34)
!--------------------------------------------------------------------------
!! cleanup and closedown

  call sumtot_delete(sumtot)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program p_sumtot
