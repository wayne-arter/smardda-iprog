module sumtot_h

  use const_kind_m
  use sum_h


!> parameters describing how to construct object
  type, public :: inumerics_t
     integer(ki4) :: npat !< pattern_of_triangles
     real(kr8) :: lambda_q !< profile_falloff_lengthscale
     real(kr8) :: a !< box_horizontal_size
     real(kr8) :: b !< box_vertical_size
     real(kr8) :: beta !< rotation_from_vertical
     real(kr8) :: epsx !< horizontal_mesh_displacement
     real(kr8) :: epsy !< vertical_mesh_displacement
     real(kr8) :: hx !< horizontal_mesh_spacing
     real(kr8) :: hy !< vertical_mesh_spacing
     character(len=80) :: formula !< sumtot formula
     real(kr8) :: f !< power split (math variable name allowed)
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters
  end type inumerics_t

  type, public :: sumtot_t
     real(kr8) :: pow !< power
     type(inumerics_t) :: n !< control  parameters
     type(sum_t) :: sum !< box and function
  end type sumtot_t

end module sumtot_h
