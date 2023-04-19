! Sample ODE codes for numerical methods by Aleks Donev
module Precision
    use ISO_C_BINDING
    implicit none
    private
    
    integer, parameter, public :: sp=kind(0.0) ! Single precision floating-point numbers
    integer, parameter, public :: dp=kind(0.0D0) ! Double precision floating-point numbers
    ! Quad precision floating-point numbers (not very portable)
    integer, parameter, public :: qp=C_LONG_DOUBLE !C_FLOAT128 in gfortran ensures IEEE quad precision
    
    integer, parameter, public :: wp = dp ! Choose working precision (usually double)

end module Precision

