module misc_m

  use const_kind_m
  use const_numphys_h
  use log_m

  implicit none
  private

! public subroutines
  public :: &
 &misc_getfileunit, & !< find new file unit number
 &misc_close, &  !< close file with error handling
 &misc_line2d, &  !< sample 2-D straight line between end-points
 &misc_lines2d, &  !< sample 2-D straight lines between end-points
 &misc_anglevec, & !< return angle in degrees between two vectors
 &misc_fsuffixget, & !< get file suffix as lower case string
 &misc_fileafter, & !< move to point in file after line beginning with string
 &misc_countnos

! public types

! private variables
  character(*), parameter :: m_name='misc_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer :: nread  !< file unit for read
  integer :: status  !< status flag
  integer :: istatus  !<  status flag
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> find new file unit number
subroutine misc_getfileunit(kunit)

  !! arguments
  integer, intent(out) :: kunit !< file unit

  !! local
  integer :: i !< loop counter
  logical :: unitused !< flag to test unit is available

  kunit=199
  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        kunit=i
        exit
     end if
  end do

end subroutine misc_getfileunit
!---------------------------------------------------------------------
!> close file with error handling
subroutine misc_close(kunit)

  !! local
  character(*), parameter :: s_name='misc_close' !< subroutine name
  integer, intent(in) :: kunit    !< unit number

  ibuf1=' '
  inquire(unit=kunit,named=iltest,name=ibuf1, iostat=istatus)

  !! file has a name, so use it
  if (iltest) call log_value('Closing file name',ibuf1)

  close(kunit,iostat=status)

  if(status/=0)then
     !! error closing file
     call log_value('Unit number',kunit)
     call log_value('Iostat',status)
     call log_error(m_name,s_name,1,error_fatal,'Error closing file')
  else
     call log_error(m_name,s_name,2,log_info,'File successfully closed')
  end if

end subroutine misc_close
!---------------------------------------------------------------------
!> sample 2-D straight line between end-points
subroutine misc_line2d(rpos,zpos,stpos,finpos,ldiv)

  !! arguments
  real(kr8), dimension(:), allocatable, intent(out) :: rpos !< positions in 1 coordinate
  real(kr8), dimension(:), allocatable, intent(out) :: zpos !< positions in 2 coordinate
  real(kr8), dimension(3), intent(in) :: stpos !< start position
  real(kr8), dimension(3), intent(in) :: finpos !< finish position
  integer(ki4), intent(in) :: ldiv !< number of divisions in straight line joining start and finish positions

  !! local
  character(*), parameter :: s_name='misc_line2d' !< subroutine name
  real(kr8) :: zr1 !< value of \f$ R \f$
  real(kr8) :: zr2 !< value of \f$ R \f$
  real(kr8) :: zz1 !< value of \f$ Z \f$
  real(kr8) :: zz2 !< value of \f$ Z \f$
  real(kr8) :: zdelr !< value of \f$ \Delta R \f$
  real(kr8) :: zdelz !< value of \f$ \Delta Z \f$

  !! uniformly divided line between stpos and finpos (as polars (R,Z) )
  zr1=stpos(1)
  zr2=finpos(1)
  zz1=stpos(2)
  zz2=finpos(2)
  allocate(rpos(ldiv+1),zpos(ldiv+1),stat=status)
  call log_alloc_check(m_name,s_name,51,status)
  zdelr=(zr2-zr1)/ldiv
  zdelz=(zz2-zz1)/ldiv
  do j=1,ldiv+1
     rpos(j)=zr1+(j-1)*zdelr
     zpos(j)=zz1+(j-1)*zdelz
  end do

end subroutine misc_line2d
!---------------------------------------------------------------------
!> sample 2-D straight lines between end-points
subroutine misc_lines2d(rpos,zpos,npos,ldiv,div)

  !! arguments
  real(kr8), dimension(:), allocatable, intent(inout) :: rpos !< positions in 1 coordinate
  real(kr8), dimension(:), allocatable, intent(inout) :: zpos !< positions in 2 coordinate
  integer(ki4), intent(inout)  :: npos !< number of entries in rpos and zpos arrays
  integer(ki4), intent(in) :: ldiv !< number of subdivisions in straight line joining each position in rpos and zpos
  real(kr8), intent(in), optional :: div !< approx size of division in straight line joining each position in rpos and zpos (INERT)

  !! local
  character(*), parameter :: s_name='misc_lines2d' !< subroutine name
  real(kr8), dimension(:), allocatable :: rposn !< new positions in 1 coordinate
  real(kr8), dimension(:), allocatable :: zposn !< new positions in 2 coordinate
  integer(ki4) :: iposn !< number of positions in output array
  real(kr8) :: zr1 !< value of \f$ R \f$
  real(kr8) :: zr2 !< value of \f$ R \f$
  real(kr8) :: zz1 !< value of \f$ Z \f$
  real(kr8) :: zz2 !< value of \f$ Z \f$
  real(kr8) :: zdelr !< value of \f$ \Delta R \f$
  real(kr8) :: zdelz !< value of \f$ \Delta Z \f$

  !! nothing to do if no subdivisions
  if (ldiv==0) return

  iposn=(npos-1)*ldiv+1
  allocate(rposn(iposn),zposn(iposn),stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  !! uniformly divided line between each input point (as polars (R,Z) )
  i=1
  zr1=rpos(1)
  zz1=zpos(1)
  rposn(1)=zr1
  zposn(1)=zz1
  do l=2,npos
     zr2=rpos(l)
     zz2=zpos(l)
     zdelr=(zr2-zr1)/ldiv
     zdelz=(zz2-zz1)/ldiv
     do j=1,ldiv
        i=i+1
        rposn(i)=zr1+j*zdelr
        zposn(i)=zz1+j*zdelz
     end do
     zr1=zr2
     zz1=zz2
  end do
  ! now move  back to subroutine arguments
  deallocate(rpos,zpos)
  allocate(rpos(iposn),zpos(iposn),stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  rpos=rposn
  zpos=zposn
  npos=iposn
  deallocate(rposn,zposn)

end subroutine misc_lines2d
!---------------------------------------------------------------------
!> return angle in degrees between two vectors
function misc_anglevec(pa,pb)
  !! arguments
  real(kr4), dimension(3), intent(in) :: pa !< local variable
  real(kr4), dimension(3), intent(in) :: pb !< local variable
  !! local
  real(kr4) :: zabsa !< local variable
  real(kr4) :: zabsb !< local variable
  real(kr4) :: zadotb !< local variable
  real(kr4) :: zangle !< local variable
  real(kr4) :: misc_anglevec !< return variable

  zangle = 0
  zabsa = sqrt(dot_product(pa,pa))
  !print *, 'zabsa:',zabsa !print

  zabsb = sqrt(dot_product(pb,pb))
  !print *, 'zabsb:',zabsb !print

  zadotb = dot_product(pa,pb)
  !print *, 'product',zadotb !print
  if (zabsa *zabsb > 0 ) then
     zangle = acos(zadotb/(zabsb*zabsa))
  end if
  misc_anglevec = zangle*180/const_pi

  !print *, misc_anglevec !print
  return
end function misc_anglevec
!---------------------------------------------------------------------
!> get file suffix as lower case string
subroutine misc_fsuffixget(filename,filesuffix,kerr)
  !! arguments
  character(*), intent(in) :: filename !< file name
  character(*), intent(out) :: filesuffix !< file suffix
  integer(ki4), intent(out)  :: kerr !< return status

  !! local
  integer(ki4) :: indot !< position of dot in filename
  integer(ki4) :: ilenf !< length of input string
  integer(ki4) :: ilent !< length of suffix string
  character(len=80) :: icsuf  !< file suffix local variable

  ibuf1=adjustl(filename)
  ilenf=len_trim(ibuf1)
  indot=index(ibuf1,'.',.TRUE.)
  if (indot>=ilenf.OR.indot==0) then
     kerr=error_warning
     filesuffix='nul'
  else
     kerr=0
     icsuf=ibuf1(indot+1:ilenf)
     ilent=len_trim(icsuf)
     call lowor(icsuf,1,ilent)
     filesuffix=icsuf
  end if

end subroutine misc_fsuffixget
!---------------------------------------------------------------------
subroutine misc_countnos(bigbuf,kfmt)
  character(len=*),intent(in) :: bigbuf !< buffer for input
  integer(ki4), intent(out) :: kfmt !< format of buffer - number of 3-vectors
  character(len=132) :: ibuf !< buffer for input/output
  integer(ki4) :: ilen !< length of string
  integer(ki4) :: iblan !< number of blank substrings
  integer(ki4) :: isw !< switch on if last character was not blank
  integer(ki4) :: ji !< loop variable
  iblan=0
  ibuf=adjustl(bigbuf)
  ilen=len_trim(ibuf)
  isw=1
  do ji=1,ilen
     if (ibuf(ji:ji)==' ') then
        if (isw/=0) then
           iblan=iblan+1
           isw=0
        end if
     else
        isw=1
     end if
  end do
  kfmt=(iblan+1)/3
end subroutine misc_countnos
!>---------------------------------------------------------------------
subroutine misc_fileafter(kfind,kin)
  !< move to point in file after line beginning with string
  !! arguments
  character(len=*), intent(in) :: kfind   !< character string to find
  integer, intent(in) :: kin   !< file unit for reading

  !! local
  character(*), parameter :: s_name='misc_fileafter' !< subroutine name
  character(len=80) :: ibuff   !< character string

  do
     read(kin,fmt='(a)',iostat=status) ibuff
     call log_read_check(m_name,s_name,1,status)
     if (adjustl(ibuff)==kfind) exit
  end do

  if (adjustl(ibuff)/=kfind) then
     ! jm string not found
     call log_error(m_name,s_name,10,error_fatal,'Error reading object data')
  end if

end subroutine misc_fileafter

end module misc_m
