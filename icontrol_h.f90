module icontrol_h

  use const_kind_m

!! public types

!> run control parameters
  type, public :: iparams_t
     character(len=80) :: control !< option control parameter
     real(kr8) :: realpar !< real control parameter
     integer(ki4) :: intpar !< integer control parameter
     logical :: logicpar !< logical control parameter
  end type iparams_t

!> file names
  type, public :: ifiles_t
     character(len=80)  :: out       !< output data
     character(len=80)  :: log           !< log file
     character(len=80)  :: sumtotdata         !< sumtot input data file
     character(len=80)  :: sumtotout         !< sumtot output data file
     character(len=80)  :: vtk   !< vtk file
     character(len=80)  :: gnu !< gnuplot file
  end type ifiles_t

!> plot output selectors
  type, public :: iplots_t
     logical  :: sumtotout !< sumtot output data selector
     logical  :: vtk   !< vtk plot selector
     logical  :: gnu !< gnuplot plot selector
  end type iplots_t

end module icontrol_h
