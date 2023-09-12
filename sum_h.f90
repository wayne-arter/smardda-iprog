module sum_h

  use const_kind_m


!> parameters describing how to construct object
  type, public :: sunumerics_t
     integer(ki4) :: npat !< pattern_of_triangles
     real(kr8) :: lambda_q !< profile_falloff_lengthscale
     real(kr8) :: a !< box_horizontal_size
     real(kr8) :: b !< box_vertical_size
     real(kr8) :: beta !< rotation_from_vertical
     real(kr8) :: epsx !< horizontal_mesh_displacement
     real(kr8) :: epsy !< vertical_mesh_displacement
     real(kr8) :: hx !< horizontal_mesh_spacing
     real(kr8) :: hy !< vertical_mesh_spacing
     character(len=80) :: formula !< sum formula
     real(kr8) :: f !< power split (math variable name allowed)
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters
  end type sunumerics_t

  type, public :: sum_t
     real(kr8) :: pow !< power
     real(kr8) :: aext !< box horizontal extent
     real(kr8) :: bext !< box vertical extent
     integer(ki4) :: nh !< number of mesh cells in horizontal
     integer(ki4) :: mh !< number of mesh cells in vertical
     real(kr8) :: ah !< half box width
     real(kr8) :: bh !< half box height
     real(kr8) :: hxsxh !< sixth horizontal mesh spacing
     real(kr8) :: hysxh !< sixth vertical mesh spacing
     real(kr8) :: sumest !< estimate of integral
     real(kr8) :: sumexact !< analytic value of integral
     
     type(sunumerics_t) :: n !< control  parameters
  end type sum_t

end module sum_h
