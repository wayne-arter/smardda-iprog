module sum_m

  use sum_h
  use log_m
  use const_numphys_h
  use const_kind_m
  use misc_m

  implicit none
  private

! public subroutines
  public :: &
  sum_initfile,  & !< open file
  sum_readcon,  & !< read data from file
  sum_solve,  & !<  perform UJ quadrature
  sum_solvehex,  & !<  perform  quadrature on hexagons
  sum_userdefined,  & !< user-defined function
  sum_exp,  & !< exp function
  sum_fn, &  !< general external function call
  sum_dia, &  !< object diagnostics to log file
  sum_initwrite, & !< open new file, making up name
  sum_write, &  !< write out object
  sum_writeg, &  !< write out object as gnuplot
  sum_writehexg, &  !< write out object as gnuplot
  sum_writev, &  !< write out object as vtk
  sum_delete, & !< delete object
  sum_close, & !< close file
  sum_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='sum_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: ninsu=5     !< control file unit number
  integer(ki4), save  :: noutsu=6      !< output file unit number
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
subroutine sum_initfile(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='sum_initfile' !< subroutine name
  logical :: unitused !< flag to test unit is available

  if (trim(file)=='null') then
     call log_error(m_name,s_name,1,log_info,'null filename ignored')
     return
  end if

! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then ninsu=i if (present(channel)) channel=i exit end if end do

  call misc_getfileunit(ninsu)
  if (present(channel)) channel=ninsu

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  open(unit=ninsu,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,2,error_fatal,'Cannot open control data file')
     stop
  end if

end subroutine sum_initfile
!---------------------------------------------------------------------
!> read data from file
subroutine sum_readcon(selfn,channel)

  !! arguments
  type(sunumerics_t), intent(out) :: selfn !< type which data will be assigned to
  integer(ki4), intent(in),optional :: channel   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='sum_readcon' !< subroutine name
  character(len=80) :: sum_formula !< formula to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< variable with meaningful name

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< local variable
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable

  !! sum parameters
  namelist /sumparameters/ &
 &power_split, sum_formula, &
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters

  !! set default sum parameters
  power_split=0.5_kr8

  sum_formula='unset'
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     ninsu=channel
  end if

  !!read sum parameters
  read(ninsu,nml=sumparameters,iostat=status)
  if(status/=0) then
     !!dump namelist contents to logfile to assist error location
     print '("Fatal error reading sum parameters")'
     call log_getunit(ilog)
     write(ilog,nml=sumparameters)
     call log_error(m_name,s_name,1,error_fatal,'Error reading sum parameters')
  end if

  call lowor(sum_formula,1,len_trim(sum_formula))
  !! check for valid data

  formula_chosen: select case (sum_formula)
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
  selfn%formula=sum_formula

  selfn%f=power_split

  !! allocate arrays and assign

  selfn%nrpams=number_of_real_parameters
  selfn%nipams=number_of_integer_parameters

  formula_allocate: select case (sum_formula)

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

end  subroutine sum_readcon
!---------------------------------------------------------------------
!>  perform  UJ quadrature
subroutine sum_solve(self)

  !! arguments
  type(sum_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sum_solve' !< subroutine name
  real(kr8) :: cosb !< local variable
  real(kr8) :: sinb !< local variable
  real(kr8), dimension(8) :: x !< local variable
  real(kr8), dimension(8) :: y !< local variable
  real(kr8) :: xtfm !< local variable
  real(kr8) :: ytfm !< local variable
  real(kr8) :: ah !< local variable
  real(kr8) :: bh !< local variable
  integer(ki4) :: i6 !< local variable
  integer(ki4) :: j6 !< local variable
  integer(ki4) :: ipts !< local variable
  logical :: linx !< local variable
  logical :: liny !< local variable

  cosb=cos(self%n%beta*const_degrad)
  sinb=sin(self%n%beta*const_degrad)
  self%ah=0.5d0*(self%n%a*cosb+self%n%b*sinb)
  self%bh=0.5d0*(self%n%a*sinb+self%n%b*cosb)
  self%nh=1+self%ah/self%n%hx
  self%mh=1+self%bh/self%n%hy
  self%hxsxh=self%n%hx/6
  self%hysxh=self%n%hy/6
  ah=self%n%a/2
  bh=self%n%b/2
  self%sumest=0
  ipts=0
  do j = -self%mh,self%mh
     j6=6*j
     do i = -self%nh,self%nh
        i6=6*i

        x(1)= (i6-2)*self%hxsxh+self%n%epsx
        y(1)= (j6-1)*self%hysxh+self%n%epsy

        x(2)= (i6-2)*self%hxsxh+self%n%epsx
        y(2)= (j6+1)*self%hysxh+self%n%epsy

        x(3)= (i6-1)*self%hxsxh+self%n%epsx
        y(3)= (j6+2)*self%hysxh+self%n%epsy

        x(4)= (i6+1)*self%hxsxh+self%n%epsx
        y(4)= (j6+2)*self%hysxh+self%n%epsy

        x(5)= (i6+2)*self%hxsxh+self%n%epsx
        y(5)= (j6+1)*self%hysxh+self%n%epsy

        x(6)= (i6+2)*self%hxsxh+self%n%epsx
        y(6)= (j6-1)*self%hysxh+self%n%epsy

        x(7)= (i6+1)*self%hxsxh+self%n%epsx
        y(7)= (j6-2)*self%hysxh+self%n%epsy

        x(8)= (i6-1)*self%hxsxh+self%n%epsx
        y(8)= (j6-2)*self%hysxh+self%n%epsy

        do k=1,8
!WA15          write(15,*) x(k),y(k) !WA15 output all points
           xtfm= x(k)*cosb+y(k)*sinb + ah
           linx = ( xtfm*(self%n%a-xtfm) > 0)
           if (linx) then
              ytfm= -x(k)*sinb+y(k)*cosb + bh
              liny = ( ytfm*(self%n%b-ytfm) > 0)
              if (liny) then
                 ipts=ipts+1
                 self%sumest=self%sumest+sum_fn(self,xtfm)
!WA                 write(10,*) xtfm,ytfm  !WA output selected transformed points
!WA                 write(11,*) x(k),y(k)  !WA output selected points
              end if
           end if
        end do
     end do
  end do
  !! normalise
  call log_value("displacement ",self%n%epsx)
  call log_value("number of quadrature points ",ipts)
  call log_value("unnormalised estimate ",self%sumest)
  self%sumest=self%sumest*self%n%hx*self%n%hy/(8*self%n%b)
  call log_value("normalised estimate ",self%sumest)

end subroutine sum_solve
!---------------------------------------------------------------------
!>  perform  quadrature on hexagons
subroutine sum_solvehex(self)

  !! arguments
  type(sum_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sum_solvehex' !< subroutine name
  real(kr8) :: cosb !< local variable
  real(kr8) :: sinb !< local variable
  real(kr8), dimension(8) :: x !< local variable
  real(kr8), dimension(8) :: y !< local variable
  real(kr8) :: xtfm !< local variable
  real(kr8) :: ytfm !< local variable
  real(kr8) :: ah !< local variable
  real(kr8) :: bh !< local variable
  integer(ki4) :: i6 !< local variable
  integer(ki4) :: j6 !< local variable
  integer(ki4) :: ipts !< local variable
  logical :: linx !< local variable
  logical :: liny !< local variable
  real(kr8) :: zxh !< local variable
  real(kr8) :: zyh !< local variable
  real(kr8) :: zxq !< local variable
  real(kr8) :: zyq !< local variable

  cosb=cos(self%n%beta*const_degrad)
  sinb=sin(self%n%beta*const_degrad)
  self%ah=0.5d0*(self%n%a*cosb+self%n%b*sinb)
  self%bh=0.5d0*(self%n%a*sinb+self%n%b*cosb)
  self%nh=1.5+self%ah/self%n%hx
  self%mh=1.5+self%bh/self%n%hy
  self%hxsxh=self%n%hx/6
  self%hysxh=self%n%hy/6
  zxq=self%n%hx/4
  zyq=self%n%hy/4
  zxh=self%n%hx/2
  zyh=self%n%hy/2
  ah=self%n%a/2
  bh=self%n%b/2
  self%sumest=0
  ipts=0
  do j = -self%mh,self%mh
     j6=6*j
     do i = -self%nh,self%nh
        i6=6*i

        x(1)= (i6-3)*self%hxsxh+self%n%epsx+zxh
        y(1)= (j6-1)*self%hysxh+self%n%epsy

        x(2)= (i6-3)*self%hxsxh+self%n%epsx+zxh
        y(2)= (j6+1)*self%hysxh+self%n%epsy

        x(3)= (i6-3)*self%hxsxh+self%n%epsx+zxh+zxq
        y(3)= (j6+2)*self%hysxh+self%n%epsy

        x(4)= i6*self%hxsxh+self%n%epsx+zxh
        y(4)= (j6+1)*self%hysxh+self%n%epsy

        x(5)= i6*self%hxsxh+self%n%epsx+zxh+zxq
        y(5)= (j6+2)*self%hysxh+self%n%epsy

        x(6)= i6*self%hxsxh+self%n%epsx+zxh
        y(6)= (j6-1)*self%hysxh+self%n%epsy

        x(7)= i6*self%hxsxh+self%n%epsx+zxq
        y(7)= (j6-2)*self%hysxh+self%n%epsy

        x(8)= (i6-3)*self%hxsxh+self%n%epsx+zxq
        y(8)= (j6-2)*self%hysxh+self%n%epsy

        do k=1,8
!WA15          write(15,*) x(k),y(k) !WA15
           xtfm= x(k)*cosb+y(k)*sinb + ah
           linx = ( xtfm*(self%n%a-xtfm) > 0)
           if (linx) then
              ytfm= -x(k)*sinb+y(k)*cosb + bh
              liny = ( ytfm*(self%n%b-ytfm) > 0)
              if (liny) then
                 ipts=ipts+1
                 self%sumest=self%sumest+sum_fn(self,xtfm)
!WA                 write(10,*) xtfm,ytfm !WA
!WA                 write(11,*) x(k),y(k) !WA
              end if
           end if
        end do
     end do
  end do
  !! normalise
  call log_value("displacement ",self%n%epsx)
  call log_value("number of quadrature points ",ipts)
  call log_value("unnormalised estimate ",self%sumest)
  self%sumest=self%sumest*self%n%hx*self%n%hy/(8*self%n%b)
  call log_value("normalised estimate ",self%sumest)

end subroutine sum_solvehex
!---------------------------------------------------------------------
!> output to log file
subroutine sum_dia(self)

  !! arguments
  type(sum_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sum_dia' !< subroutine name

  call log_value("estimate of integral",self%sumest)

end subroutine sum_dia
!---------------------------------------------------------------------
!> exp function
function sum_exp(self,psi)

  !! arguments
  type(sum_t), intent(in) :: self !< module object
  real(kr8) :: sum_exp !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='sum_exp' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=exp(-zpos/self%n%lambda_q)/self%n%lambda_q
  sum_exp=pow

end function sum_exp
!---------------------------------------------------------------------
!> userdefined function
function sum_userdefined(self,psi)

  !! arguments
  type(sum_t), intent(in) :: self !< module object
  real(kr8) :: sum_userdefined !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='sum_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  integer(ki4) :: ilocal !< local integer variable

  zpos=psi
  pow=0._kr8
  !! TBD
  sum_userdefined=pow

end function sum_userdefined
!---------------------------------------------------------------------
!> general external function call
function sum_fn(self,psi)

  !! arguments
  type(sum_t), intent(in) :: self !< module object
  real(kr8) :: sum_fn !< local variable
  real(kr8), intent(in) :: psi !< position in \f$ \psi \f$

  !! local variables
  character(*), parameter :: s_name='sum_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  pow=0._kr8
  !! select sum
  formula_chosen: select case (self%n%formula)
  case('exp')
     pow=sum_exp(self,psi)
  case('userdefined')
     pow=sum_userdefined(self,psi)
  end select formula_chosen

  !! return sum
  sum_fn=pow

end function sum_fn
!---------------------------------------------------------------------
!> open new file, making up name
subroutine sum_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='sum_initwrite' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then if (present(channel)) channel=i exit end if end do noutsu=i

  call misc_getfileunit(noutsu)
  if (present(channel)) channel=noutsu

  !! open file
  outputfile=trim(fileroot)//"_sum.out"
  call log_value("Control data file",trim(outputfile))
  open(unit=noutsu,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutsu,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if

end subroutine sum_initwrite
!---------------------------------------------------------------------
!> write sum data
subroutine sum_write(self,channel)

  !! arguments
  type(sum_t), intent(in) :: self   !< sum data structure
  integer(ki4), intent(in), optional :: channel   !< output channel for sum data structure

  !! local
  character(*), parameter :: s_name='sum_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for sum data structure

  !! sort out unit
  if(present(channel)) then
     iout=channel
  else
     iout=noutsu
  end if

  write(iout,*,iostat=status) 'sum_formula'
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

end subroutine sum_write
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine sum_writeg(self,select,channel)

  !! arguments
  type(sum_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for sum data structure

  !! local
  character(*), parameter :: s_name='sum_writeg' !< subroutine name
  integer(ki4) :: iout   !< output channel for sum data structure
  real(kr8) :: cosb !< local variable
  real(kr8) :: sinb !< local variable
  real(kr8), dimension(8) :: x !< local variable
  real(kr8), dimension(8) :: y !< local variable
  real(kr8), dimension(2,0:8) :: xl !< local variable
  real(kr8), dimension(2,0:8) :: yl !< local variable
  real(kr8) :: xtfm !< local variable
  real(kr8) :: ytfm !< local variable
  real(kr8) :: zxh !< local variable
  real(kr8) :: zyh !< local variable
  real(kr8) :: ah !< local variable
  real(kr8) :: bh !< local variable
  integer(ki4) :: ii !< local variable
  integer(ki4) :: i6 !< local variable
  integer(ki4) :: j6 !< local variable
  integer(ki4) :: i2 !< local variable
  integer(ki4) :: j2 !< local variable
  integer(ki4) :: ipts !< local variable
  logical :: linx !< local variable
  logical :: liny !< local variable

  plot_type: select case(select)
  case('cartesian')

  case default
     cosb=cos(self%n%beta*const_degrad)
     sinb=sin(self%n%beta*const_degrad)
     ah=self%n%a/2
     bh=self%n%b/2
     zxh=self%n%hx/2
     zyh=self%n%hy/2
     ipts=0
     do j = -self%mh,self%mh
        j2=2*j
        j6=6*j
        do i = -self%nh,self%nh
           i2=2*i
           i6=6*i

           ii=0
           xl(1,ii)=i2*zxh+self%n%epsx
           yl(1,ii)=j2*zyh+self%n%epsy
           xl(2,ii)=(i2-1)*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           x(1)= (i6-2)*self%hxsxh+self%n%epsx
           y(1)= (j6-1)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2-1)*zxh+self%n%epsx
           yl(2,ii)=j2*zyh+self%n%epsy

           x(2)= (i6-2)*self%hxsxh+self%n%epsx
           y(2)= (j6+1)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2-1)*zxh+self%n%epsx
           yl(2,ii)=(j2+1)*zyh+self%n%epsy

           x(3)= (i6-1)*self%hxsxh+self%n%epsx
           y(3)= (j6+2)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=i2*zxh+self%n%epsx
           yl(2,ii)=(j2+1)*zyh+self%n%epsy

           x(4)= (i6+1)*self%hxsxh+self%n%epsx
           y(4)= (j6+2)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2+1)*zxh+self%n%epsx
           yl(2,ii)=(j2+1)*zyh+self%n%epsy

           x(5)= (i6+2)*self%hxsxh+self%n%epsx
           y(5)= (j6+1)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2+1)*zxh+self%n%epsx
           yl(2,ii)=j2*zyh+self%n%epsy

           x(6)= (i6+2)*self%hxsxh+self%n%epsx
           y(6)= (j6-1)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2+1)*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           x(7)= (i6+1)*self%hxsxh+self%n%epsx
           y(7)= (j6-2)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=i2*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           x(8)= (i6-1)*self%hxsxh+self%n%epsx
           y(8)= (j6-2)*self%hysxh+self%n%epsy
           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2-1)*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           do k=1,8
              !           xtfm= x(k)*cosb+y(k)*sinb + ah
              !           linx = ( xtfm*(self%n%a-xtfm) > 0)
              !           if (linx) then
              !              ytfm= -x(k)*sinb+y(k)*cosb + bh
              !              liny = ( ytfm*(self%n%b-ytfm) > 0)
              !              if (liny) then
              !                 ipts=ipts+1
              write(channel,*) xl(1,k),yl(1,k)
              write(channel,*) xl(2,k),yl(2,k)
              write(channel,*) xl(1,0),yl(1,0)
              !              end if
              !           end if
           end do
        end do
        write(channel,*) ' '
        !write(channel,*) xl(2,4),yl(2,4)
     end do

  end select plot_type

end subroutine sum_writeg
!---------------------------------------------------------------------
!> write object data as gnuplot
subroutine sum_writehexg(self,select,channel)

  !! arguments
  type(sum_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for sum data structure

  !! local
  character(*), parameter :: s_name='sum_writehexg' !< subroutine name
  integer(ki4) :: iout   !< output channel for sum data structure
  real(kr8) :: cosb !< local variable
  real(kr8) :: sinb !< local variable
  real(kr8), dimension(8) :: x !< local variable
  real(kr8), dimension(8) :: y !< local variable
  real(kr8), dimension(2,0:8) :: xl !< local variable
  real(kr8), dimension(2,0:8) :: yl !< local variable
  real(kr8) :: xtfm !< local variable
  real(kr8) :: ytfm !< local variable
  real(kr8) :: zxh !< local variable
  real(kr8) :: zyh !< local variable
  real(kr8) :: zxq !< local variable
  real(kr8) :: zyq !< local variable
  real(kr8) :: ah !< local variable
  real(kr8) :: bh !< local variable
  integer(ki4) :: ii !< local variable
  integer(ki4) :: i6 !< local variable
  integer(ki4) :: j6 !< local variable
  integer(ki4) :: i2 !< local variable
  integer(ki4) :: j2 !< local variable
  integer(ki4) :: ipts !< local variable
  logical :: linx !< local variable
  logical :: liny !< local variable

  plot_type: select case(select)
  case('cartesian')

  case default
     cosb=cos(self%n%beta*const_degrad)
     sinb=sin(self%n%beta*const_degrad)
     ah=self%n%a/2
     bh=self%n%b/2
     zxh=self%n%hx/2
     zyh=self%n%hy/2
     zxq=self%n%hx/4
     zyq=self%n%hy/4
     ipts=0
     do j = -self%mh,self%mh
        j2=2*j
        j6=6*j
        do i = -self%nh,self%nh
           i2=2*i
           i6=6*i

           ii=0
           xl(1,ii)=i2*zxh+self%n%epsx+zxq
           yl(1,ii)=j2*zyh+self%n%epsy
           xl(2,ii)=(i2-1)*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2-1)*zxh+self%n%epsx+zxq
           yl(2,ii)=j2*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2-1)*zxh+self%n%epsx+zxh
           yl(2,ii)=(j2+1)*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=i2*zxh+self%n%epsx+zxh
           yl(2,ii)=(j2+1)*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2+1)*zxh+self%n%epsx+zxh
           yl(2,ii)=(j2+1)*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2+1)*zxh+self%n%epsx+zxq
           yl(2,ii)=j2*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2+1)*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=i2*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           ii=ii+1
           xl(1,ii)=xl(2,ii-1)
           yl(1,ii)=yl(2,ii-1)
           xl(2,ii)=(i2-1)*zxh+self%n%epsx
           yl(2,ii)=(j2-1)*zyh+self%n%epsy

           do k=1,8
              !           xtfm= x(k)*cosb+y(k)*sinb + ah
              !           linx = ( xtfm*(self%n%a-xtfm) > 0)
              !           if (linx) then
              !              ytfm= -x(k)*sinb+y(k)*cosb + bh
              !              liny = ( ytfm*(self%n%b-ytfm) > 0)
              !              if (liny) then
              !                 ipts=ipts+1
              write(channel,*) xl(1,k),yl(1,k)
              write(channel,*) xl(2,k),yl(2,k)
              if (k==4) then
                 write(channel,*) xl(2,3),yl(2,3)
                 write(channel,*) xl(2,5),yl(2,5)
               else if (k==8) then
                 write(channel,*) xl(2,7),yl(2,7)
                 write(channel,*) xl(2,1),yl(2,1)
               else
                 write(channel,*) xl(1,0),yl(1,0)
               end if
              !              end if
              !           end if
           end do
           write(channel,*) ' '
        end do
        write(channel,*) ' '
        !write(channel,*) xl(2,4),yl(2,4)
     end do

  end select plot_type

end subroutine sum_writehexg
!---------------------------------------------------------------------
!> write object data as vtk
subroutine sum_writev(self,select,channel)

  !! arguments
  type(sum_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer(ki4), intent(in), optional :: channel   !< output channel for sum data structure

  !! local
  character(*), parameter :: s_name='sum_writev' !< subroutine name
  integer(ki4) :: iout   !< output channel for sum data structure

  call log_error(m_name,s_name,1,log_info,'vtk file produced')

  plot_type: select case(select)
  case('cartesian')

  case default

  end select plot_type

end subroutine sum_writev
!---------------------------------------------------------------------
!> close write file
subroutine sum_closewrite

  !! local
  character(*), parameter :: s_name='sum_closewrite' !< subroutine name

  !! close file
  close(unit=noutsu,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine sum_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine sum_delete(self)

  !! arguments
  type(sum_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='sum_delete' !< subroutine name

  formula_deallocate: select case (self%n%formula)
  case('userdefined')
     if (self%n%nrpams>0) deallocate(self%n%rpar)
     if (self%n%nipams>0) deallocate(self%n%npar)
  case default
  end select formula_deallocate

end subroutine sum_delete
!---------------------------------------------------------------------
!> close file
subroutine sum_close

  !! local
  character(*), parameter :: s_name='sum_close' !< subroutine name

  !! close file
  close(unit=ninsu,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine sum_close

end module sum_m
