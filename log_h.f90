module log_h

  use const_kind_m

  type, public :: error_line_t
     character(32) :: modname
     character(32) :: subname
     integer(ki4)  :: point
     integer(ki4)  :: severity
     character(80) :: message
  end type error_line_t

  type, public :: error_buffer_t
     type(error_line_t),dimension(:),allocatable :: error_line !< contains the individual error information
     integer(ki4) :: size !< The current size of the error buffer
     integer(ki4) :: filled  !< The number of errors currently stored in the error buffer
     integer(ki4) :: errorno !< The cumulative number of errors added to the buffer, this is not flushed
  end type error_buffer_t
  
end module log_h
