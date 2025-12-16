module hydrogenic
  use utils
  implicit none
  private

  public :: hydrogenic_schro_energy
  public :: hydrogenic_schro_wvfn

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
end module hydrogenic
