module constants
  use iso_fortran_env, only : dp => real64
  implicit none

  public

  complex(dp),   parameter :: cplx_zero = (0.0_dp, 0.0_dp)
  real(kind=dp), parameter :: particle_mass = 206.76_dp
  real(dp),      parameter :: pi = 3.14159265358979_dp 
  real(kind=dp), parameter :: speed_of_light = 137.0_dp
  real(kind=dp), parameter :: fine_structure = 0.007297352_dp
  real(kind=dp), parameter :: amu = 1.66e-27_dp
  real(kind=dp), parameter :: electron_mass = 9.11e-31_dp
  real(kind=dp), parameter :: second =1.0_dp / 2.41888e-17_dp
  real(dp),      parameter :: factorial_table(0:20) = (/&
   1.0_dp, &
   1.0_dp, &
   2.0_dp, &
   6.0_dp, &
   24.0_dp, &
   120.0_dp, &
   720.0_dp, &
   5040.0_dp, &
   40320.0_dp, &
   362880.0_dp, &
   3628800.0_dp, &
   39916800.0_dp, &
   479001600.0_dp, &
   6227020800.0_dp, &
   87178291200.0_dp, &
   1307674368000.0_dp, &
   20922789888000.0_dp, &
   355687428096000.0_dp, &
   6402373705728000.0_dp, &
   121645100408832000.0_dp, &
   2432902008176640000.0_dp/)


  type schrodinger_state
    integer       :: n ! primary
    integer       :: l ! orbital
    integer       :: Z ! nuclear charge
    real(kind=dp), allocatable, dimension(:) :: m ! All m values
    real(kind=dp), allocatable, dimension(:) :: norms
    real(kind=dp)                            :: E
    real(kind=dp)                            :: mass, red_mass
  end type schrodinger_state

  type, extends(schrodinger_state) :: dirac_state 
    real(kind=dp) :: j ! total
    integer       :: k ! dirac
  end type dirac_state

  type grid_type
    integer :: N ! num of points
    real(kind=dp) :: dx
    real(kind=dp) :: x_min ! max value in 1 dimension
    real(kind=dp), dimension(:), allocatable :: x, y, z, weight
  end type grid_type
 end module constants
