module icontrol_m

  use const_kind_m
  use log_m
  use misc_m
  use icontrol_h
  use sumtot_h
  use sumtot_m
  use sum_h
  use sum_m

  implicit none
  private

! public subroutines
  public :: &
 &icontrol_init, & !< open input control data file
 &icontrol_close, & !< close input control data file
 &icontrol_read !< read data for this run

! private variables
  character(*), parameter :: m_name='icontrol_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: nin      !< control file unit number
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine icontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='icontrol_init' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: controlfile !< control file name


! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  call misc_getfileunit(nin)

  !! open file
  controlfile=trim(fileroot)//".ctl"
  root=fileroot
  call log_value("Control data file",trim(controlfile))
  open(unit=nin,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if


end  subroutine icontrol_init
!---------------------------------------------------------------------
!> close input control data file
subroutine icontrol_close

  !! local
  character(*), parameter :: s_name='icontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin)

end  subroutine icontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine icontrol_read(file,param,inumerics,sunumerics,plot)

  !! arguments
  type(ifiles_t), intent(out) :: file !< file names
  type(sunumerics_t), intent(out) :: sunumerics !< controls for box and function
  type(iparams_t), intent(out) :: param !< control parameters
  type(inumerics_t), intent(out) :: inumerics !< controls for sumtot
  type(iplots_t), intent(out) :: plot !< plot selectors

  !!local
  character(*), parameter :: s_name='icontrol_read' !< subroutine name
  logical :: filefound !< true if file exists

  character(len=80) :: prog_input_file  !< sumtot data input file
  character(len=80) :: prog_output_file  !< sumtot data output file

  character(len=80) :: prog_control  !< option control parameter
  real(kr8) :: prog_realpar !< real control parameter
  integer(ki4) :: prog_intpar !< integer control parameter
  logical :: prog_logicpar !< logical control parameter

  logical :: plot_sumtotout !< sumtot output data selector
  logical :: plot_vtk !< vtk plot selector
  logical :: plot_gnu !< gnuplot plot selector

  !! file names
  namelist /progfiles/ &
 &prog_input_file, &
 &prog_output_file

  !! misc parameters
  namelist /miscparameters/ &
 &prog_control,&
 &prog_realpar, prog_intpar, &
 &prog_logicpar

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_sumtotout, &
 &plot_vtk, &
 &plot_gnu

  !---------------------------------------------------------------------
  !! read input file names
  prog_input_file='null'
  prog_output_file='null'
  read(nin,nml=progfiles,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
     print '("Fatal error reading input filenames")'
  end if

  file%sumtotdata = prog_input_file
  file%out = prog_output_file

  !!check file exists

  call log_value("sumtot data file, prog_input_file",trim(file%sumtotdata))
  if(file%sumtotdata.NE.'null') then
     inquire(file=prog_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find sumtot data file, ",a)',prog_input_file
        call log_error(m_name,s_name,2,error_fatal,'sumtot data file not found')
     end if
  end if

  !! create output file names from root

  !!sumtot file root
  file%sumtotout = trim(root)//"_sumtotout"
  !!vtk file
  file%vtk = trim(root)//"_vtk"
  !!gnu file roots
  file%gnu = trim(root)//"_gnu"

  !---------------------------------------------------------------------
  !! set default misc parameters
  prog_control = 'standard'
  prog_realpar = .0001_kr8
  prog_intpar = 1
  prog_logicpar = .false.

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
     print '("Fatal error reading misc parameters")'
  end if

  !! check for valid data
  if(prog_control/='standard') then
     call log_error(m_name,s_name,11,error_fatal,'only standard control allowed')
  end if
  ! positive real parameter
  if(prog_realpar<0._kr8) then
     call log_error(m_name,s_name,20,error_fatal,'realpar must be non-negative')
  end if
  ! positive integer parameter
  if(prog_intpar<=0) then
     call log_error(m_name,s_name,30,error_fatal,'intpar must be positive')
  end if
  if(prog_logicpar) &
 &call log_error(m_name,s_name,40,error_warning,'logicpar is true')

  param%control=prog_control
  param%realpar=prog_realpar
  param%intpar=prog_intpar
  param%logicpar=prog_logicpar

  !---------------------------------------------------------------------
  !! set default plot selections
  plot_sumtotout = .false.
  plot_vtk = .false.
  plot_gnu = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,50,error_fatal,'Error reading plot selections')
     print '("Fatal error reading plot selections")'
  end if

  !! store values
  call sumtot_readcon(inumerics,nin)
  
  call sum_readcon(sunumerics,nin)
  sunumerics%lambda_q = inumerics%lambda_q
  sunumerics%a = inumerics%a
  sunumerics%b = inumerics%b
  sunumerics%beta = inumerics%beta
  sunumerics%epsx = inumerics%epsx
  sunumerics%epsy = inumerics%epsy
  sunumerics%hx = inumerics%hx
  sunumerics%hy = inumerics%hy
  sunumerics%f = inumerics%f
  plot%sumtotout = plot_sumtotout
  plot%vtk = plot_vtk
  plot%gnu = plot_gnu

end  subroutine icontrol_read

end module icontrol_m
