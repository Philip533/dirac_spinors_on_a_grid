module hydrogenic
  use utils
  implicit none
  private

  public :: hydrogenic_schro_energy
  public :: hydrogenic_schro_wvfn
  public :: hydrogenic_unbound_wvfn

contains
  function hydrogenic_schro_energy(state) result(energy)
    implicit none

    type(schrodinger_state), intent(in) :: state
    real(kind=dp)                 :: energy

    energy = -state%red_mass * state%Z**2 / (2.0_dp * real(state%n,dp)**2)
  end function hydrogenic_schro_energy
  function hydrogenic_schro_wvfn(r, state) result(wvfn)
    implicit none

    real(kind=dp), intent(in) :: r
    type(schrodinger_state), intent(in) :: state
    real(kind=dp)             :: wvfn, prefactor, arg, laguerre_poly
    real(kind=dp)             :: Z, mu
    integer                   :: n, l

    n = state%n; l = state%l; Z = state%Z; mu = state%red_mass
    
    arg = Z * mu / n

    prefactor = sqrt((2.0_dp * arg)**3 * factorial(n-l-1)/(2.0_dp * real(n,dp) * factorial(n+l)))&
                & * exp(-arg * r) * (2.0_dp * arg * r) ** l

    laguerre_poly = gen_laguerre_poly(n-l-1,real(2*l+1,dp), 2.0_dp * arg * r)

    wvfn = prefactor * laguerre_poly
    
    
  end function

  function hydrogenic_unbound_wvfn(r, E, Z) result(wvfn)
    implicit none

    real(kind=dp),           intent(in) :: r(3)
    real(kind=dp),           intent(in) :: E
    integer,                 intent(in) :: Z
    real(kind=dp)                    :: wvfn(4)

    ! Continuum parameter of the electron. Dependent on Z and delta E
    real(fgsl_double) :: y

    ! Calculate up to L = 3
    real(fgsl_double):: fc_array(4), norm_array(4)

    real(fgsl_double) :: l_min, x(3), f_exponent, k_max
    real(fgsl_double)     :: dE
    integer(fgsl_int) :: status1
    integer :: i
    

    dE = real(E* hartree, fgsl_double)

    ! Continuum parameter
    ! IDENTICAL IN PRACTICE
    y = real(Z* fine_structure / sqrt(2.0_fgsl_double * dE + dE**2), fgsl_double)
    ! y = real(Z / (dE/(speed_of_light*bohr_radius)), fgsl_double)
    ! write(*,*) "Y = ", y

    ! We always start at 0
    l_min = 0.0_fgsl_double

    ! position to the correct type
    x = real(r, fgsl_double)

    ! Calculate the F function which will have the same number of elements
    ! as the length of fc_array
    status1 = fgsl_sf_coulomb_wave_F_array(l_min, y, norm2(x), fc_array, f_exponent)  
    ! status1 = fgsl_sf_coulomb_CL_array(l_min, y, norm_array);


    ! We have obtained the radial part
    wvfn(:) = fc_array !/ norm_array


  end function
end module hydrogenic
