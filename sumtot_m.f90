module sumtot_m

  use sumtot_h
  use log_m
  use const_numphys_h
  use const_kind_m
  use misc_m
  use halton_m
  use sum_h
  use sum_m

  implicit none
  private

! public subroutines
  public :: &
  sumtot_initfile,  & !< open file
  sumtot_readcon,  & !< read data from file
  sumtot_solve,  & !< generic subroutine
  sumtot_userdefined,  & !< user-defined function
  sumtot_fn, &  !< general external function call
  sumtot_dia, &  !< object diagnostics to log file
  sumtot_initwrite, & !< open new file, making up name
  sumtot_write, &  !< write out object
  sumtot_writeg, &  !< write out object as gnuplot
  sumtot_writev, &  !< write out object as vtk
  sumtot_delete, & !< delete object
  sumtot_close, & !< close file
  sumtot_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='sumtot_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: nini=5     !< control file unit number
  integer(ki4), save  :: nouti=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4)  :: ilog      !< for namelist dump after error

  contains
!---------------------------------------------------------------------
!> open file
subroutine sumtot_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='sumtot_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     return
  end if

! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nini=i if (present(channel)) channel=i exit end if end do

  call misc_getfileunit(nini)
  if (present(channel)) channel=nini

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  open(unit=nini,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if

end subroutine sumtot_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine sumtot_readcon(selfn,channel)

  !! arguments
  type(inumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='sumtot_readcon' !< subroutine name
  character(len=80) :: sumtot_formula !< formula to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< variable with meaningful name
  integer(ki4) :: pattern_of_triangles !< self-explanatory local
  real(kr8) :: profile_falloff_lengthscale !< self-explanatory local
  real(kr8) :: box_horizontal_size !< self-explanatory local
  real(kr8) :: box_vertical_size !< self-explanatory local
  real(kr8) :: rotation_from_vertical !< self-explanatory local
  real(kr8) :: horizontal_mesh_displacement !< self-explanatory local
  real(kr8) :: vertical_mesh_displacement !< self-explanatory local
  real(kr8) :: horizontal_mesh_spacing !< self-explanatory local
  real(kr8) :: vertical_mesh_spacing !< self-explanatory local

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< local variable
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable

  !! sumtot parameters
  namelist /sumtotparameters/ &
 & pattern_of_triangles , &
 & profile_falloff_lengthscale , &
 & box_horizontal_size , &
 & box_vertical_size , &
 & rotation_from_vertical , &
 & horizontal_mesh_displacement , &
 & vertical_mesh_displacement , &
 & horizontal_mesh_spacing , &
 & vertical_mesh_spacing , &
 &power_split, sumtot_formula, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  !! set default sumtot parameters
  power_split=1.0_kr8

  sumtot_formula='additional'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  ! number of samples to take
  number_of_integer_parameters=1
  general_integer_parameters(1)=100
  pattern_of_triangles = 0
  profile_falloff_lengthscale = 0.1
  box_horizontal_size = 1
  box_vertical_size = 2
  rotation_from_vertical = 0
  horizontal_mesh_displacement = 0.01
  vertical_mesh_displacement = 0.01
  horizontal_mesh_spacing = 0.01
  vertical_mesh_spacing = 0.01

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     nini=channel
  end if

  !!read sumtot parameters
  read(nini,nml=sumtotparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading sumtot parameters")'
     call log_getunit(ilog)
     write(ilog,nml=sumtotparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading sumtot parameters')
  end if

  call lowor(sumtot_formula,1,len_trim(sumtot_formula))
  !! check for valid data

  formula_chosen: select case (sumtot_formula)
  case('unset','exp')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')

  case('expdouble')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,21,error_fatal,'power_split must be >=0 and <=1')

  case('userdefined','additional')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,44,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,45,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,46,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,47,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,48,error_fatal,'no parameters set')

  end select formula_chosen

  !! store values
  selfn%formula=sumtot_formula

  selfn%f=power_split

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters

  formula_allocate: select case (sumtot_formula)

  case('userdefined','additional')
     if (number_of_real_parameters>0) then
        allocate(selfn%rpar(number_of_real_parameters), stat=status)
        call log_alloc_check(m_name,s_name,65,status)
        selfn%rpar=general_real_parameters(:number_of_real_parameters)
     end if
     if (number_of_integer_parameters>0) then
        allocate(selfn%npar(number_of_integer_parameters), stat=status)
        call log_alloc_check(m_name,s_name,66,status)
        selfn%npar=general_integer_parameters(:number_of_integer_parameters)
     end if
  case default
  end select formula_allocate
  selfn%npat = pattern_of_triangles
  selfn%lambda_q = profile_falloff_lengthscale
  selfn%a = box_horizontal_size
  selfn%b = box_vertical_size
  selfn%beta = rotation_from_vertical
  selfn%epsx = horizontal_mesh_displacement
  selfn%epsy = vertical_mesh_displacement
  selfn%hx = horizontal_mesh_spacing
  selfn%hy = vertical_mesh_spacing

end  subroutine sumtot_readcon
!---------------------------------------------------------------------
!> generic subroutine
subroutine sumtot_solve(self)

  !! arguments
  type(sumtot_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sumtot_solve' !< subroutine name
  integer(ki4) :: isamp
  real(kr8) :: zepsta
  real(kr8) :: delteps

  !! Analytic result
  self%sum%sumexact = 1-self%n%lambda_q*sum_fn(self%sum,self%n%a)
  call log_value("exact value  of integral",self%sum%sumexact)

  isamp=self%n%npar(1)
  zepsta=-self%n%hx/2
  delteps=0
  if (isamp>1) delteps=self%n%hx/(isamp-1)

  ! loop over displacements
  do j=1,isamp
  self%sum%n%epsx=zepsta+(j-1)*delteps
  !! Numerical integral
  mesh_pattern: select case (self%n%npat)
  case(0)
    ! union jack
    call sum_solve(self%sum)
  case(1)
    ! hexagons
    call sum_solvehex(self%sum)
  case(10)
    ! random numbers
    call sum_solverand(self%sum)
  case(11)
    ! quasirandom numbers
    call sum_solveqrand(self%sum)
  end select mesh_pattern
  end do

end subroutine sumtot_solve
!---------------------------------------------------------------------
!> output to log file
subroutine sumtot_dia(self)

  !! arguments
  type(sumtot_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sumtot_dia' !< subroutine name

  call sum_dia(self%sum)

end subroutine sumtot_dia
!---------------------------------------------------------------------
!> userdefined function
function sumtot_userdefined(self,psi)

  !! arguments
  type(sumtot_t), intent(in) :: self !< module object
  real(kr8) :: sumtot_userdefined !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='sumtot_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !> user defines \Tt{pow} here
  !! .....
  !! return sumtot
  sumtot_userdefined=pow

end function sumtot_userdefined
!---------------------------------------------------------------------
!> general external function call
function sumtot_fn(self,psi)

  !! arguments
  type(sumtot_t), intent(in) :: self !< module object
  real(kr8) :: sumtot_fn !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='sumtot_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  pow=0._kr8
  !! select sumtot
  formula_chosen: select case (self%n%formula)
  case('userdefined','additional')
     pow=sumtot_userdefined(self,psi)
  end select formula_chosen

  !! return sumtot
  sumtot_fn=pow

end function sumtot_fn
!---------------------------------------------------------------------
!> open new file, making up name
subroutine sumtot_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='sumtot_initwrite' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then if (present(channel)) channel=i exit end if end do nouti=i

  call misc_getfileunit(nouti)
  if (present(channel)) channel=nouti

  !! open file
  outputfile=trim(fileroot)//"_sumtot.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=nouti,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=nouti,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine sumtot_initwrite
!---------------------------------------------------------------------
!> write sumtot data
subroutine sumtot_write(self,channel)

  !! arguments
  type(sumtot_t), intent(in) :: self   !< sumtot data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for sumtot data structure

  !! local
  character(*), parameter :: s_name='sumtot_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for sumtot data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=nouti
  end if

  write(iout,*,iostat=status) 'sumtot_formula'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%n%formula
  call log_write_check(m_name,s_name,19,status)
  write(iout,*,iostat=status) 'f'
  call log_write_check(m_name,s_name,20,status)
  write(iout,*,iostat=status) self%n%f
  call log_write_check(m_name,s_name,21,status)
  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) self%n%nrpams
  call log_write_check(m_name,s_name,47,status)
  if (self%n%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,48,status)
     write(iout,*,iostat=status) self%n%rpar
     call log_write_check(m_name,s_name,49,status)
  end if
  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,50,status)
  write(iout,*,iostat=status) self%n%nipams
  call log_write_check(m_name,s_name,51,status)
  if (self%n%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,52,status)
     write(iout,*,iostat=status) self%n%npar
     call log_write_check(m_name,s_name,53,status)
  end if

end subroutine sumtot_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine sumtot_writeg(self,select,channel)

  !! arguments
  type(sumtot_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for sumtot data structure

  !! local
  character(*), parameter :: s_name='sumtot_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for sumtot data structure

  call log_error(m_name,s_name,1,log_info,'gnuplot file produced')

  plot_type: select case(select)
  case('cartesian')

  case default
  
  mesh_pattern: select case (self%n%npat)
  case(0)
    ! union jack
    call sum_writeg(self%sum,select,channel)
  case(1)
    ! hexagons
    call sum_writehexg(self%sum,select,channel)
  end select mesh_pattern

  end select plot_type

end subroutine sumtot_writeg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine sumtot_writev(self,select,channel)

  !! arguments
  type(sumtot_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for sumtot data structure

  !! local
  character(*), parameter :: s_name='sumtot_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for sumtot data structure

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine sumtot_writev
!---------------------------------------------------------------------
!> close write file
subroutine sumtot_closewrite

  !! local
  character(*), parameter :: s_name='sumtot_closewrite' !< subroutine name

  !! close file
  close(unit=nouti,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine sumtot_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine sumtot_delete(self)

  !! arguments
  type(sumtot_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sumtot_delete' !< subroutine name

  formula_deallocate: select case (self%n%formula)
  case('userdefined','additional')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate

end subroutine sumtot_delete
!---------------------------------------------------------------------
!> close file
subroutine sumtot_close

  !! local
  character(*), parameter :: s_name='sumtot_close' !< subroutine name

  !! close file
  close(unit=nini,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine sumtot_close

end module sumtot_m
