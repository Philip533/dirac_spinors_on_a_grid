program dirac
  use iso_fortran_env, only : dp => real64
  implicit none

  ! Contains the quantum numbers
  type dirac_state
    integer       :: n ! primary
    integer       :: l ! orbital
    real(kind=dp) :: j ! total
    integer       :: k ! dirac
    integer       :: Z ! nuclear charge
    real(kind=dp), allocatable, dimension(:) :: m ! All m_j values
    real(kind=dp), allocatable, dimension(:) :: norms
    real(kind=dp)                            :: E
    real(kind=dp)                            :: mass, red_mass
  end type dirac_state

  type grid_type
    integer :: N ! num of points
    real(kind=dp) :: dx
    real(kind=dp) :: x_min ! max value in 1 dimension
    real(kind=dp), dimension(:), allocatable :: x, y, z, weight
  end type grid_type

  complex(dp), parameter :: cplx_zero = (0.0_dp, 0.0_dp)
  real(kind=dp), parameter :: particle_mass = 206.76_dp
  integer,  parameter :: nuclear_Z = 2
  ! real(kind=dp), parameter :: particle_mass = 1.0_dp
  real(dp), parameter :: pi = 3.14159265358979_dp 
  real(kind=dp), parameter :: speed_of_light = 137.0_dp
  real(kind=dp), parameter :: fine_structure = 0.007297352_dp
  real(kind=dp), parameter :: amu = 1.66e-27_dp
  real(kind=dp), parameter :: electron_mass = 9.11e-31_dp
  real(kind=dp), parameter :: second =1.0_dp / 2.41888e-17_dp
  real(dp), parameter :: factorial_table(0:20) = (/&
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

  integer :: i, j, k, nuc_z, N, m1, m2

  ! Define the Pauli \sigma vector
  complex(kind=dp), dimension(3, 2, 2) :: pauli_vec

  ! Define the Dirac \alpha vector
  complex(kind=dp), dimension(3, 4, 4) :: alpha_vec

  ! Dirac spinor
  complex(kind=dp), dimension(4)  :: dirac_spinor1, dirac_spinor2
  complex(kind=dp),dimension(3)   :: wvfn_product
  real(kind=dp) :: dx, transition_energy, red_mass, x_min, final_average
  real(kind=dp) :: bispin_integrand, normalisation1
  complex(kind=dp) :: integrand

  type(dirac_state) :: initial_state, final_state
  type(grid_type) :: grid
  real(kind=dp), dimension(2) :: bispinor1, bispinor2, PQ_vec
  real(kind=dp), dimension(3) :: r
  real(kind=dp)                 :: gamma_k, gauss_norm
  real(kind=dp), dimension(2,2) :: prefactor_matrix

  ! Interface to the routine for calculating the radial part R(r)
  ! For a Dirac spinor, this would be P/r and Q/r
  abstract interface 
    subroutine radial_function(r, ds, f) 
      import :: dp
      import :: dirac_state

      real(kind=dp), dimension(3), intent(in)  :: r
      type(dirac_state), intent(in)            :: ds
      real(kind=dp), intent(out), dimension(2) :: f
    end subroutine
  end interface 

  procedure(radial_function), pointer :: radial 

  ! n, l, j, m, Z
  call set_quantum_numbers(initial_state, 2, 1, 1.5_dp, particle_mass, nuclear_Z)
  call set_quantum_numbers(final_state, 1, 0, 0.5_dp, particle_mass, nuclear_Z)

  ! Initialise Pauli matrices in the vector
  pauli_vec(1,:,:) = reshape((/ 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp /), shape(pauli_vec(1,:,:)))
  pauli_vec(2,:,:) = reshape((/ (0.0_dp,0.0_dp), (0.0_dp,-1.0_dp), (0.0_dp,1.0_dp), (0.0_dp,0.0_dp) /), shape(pauli_vec(1,:,:)))
  pauli_vec(3,:,:) = reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp /), shape(pauli_vec(1,:,:)))

  ! Define the \alpha vector
  do i = 1, 3
    alpha_vec(i,1:2, 1:2) = 0.0_dp
    alpha_vec(i,3:4, 3:4) = 0.0_dp
    alpha_vec(i,3:4, 1:2) = pauli_vec(i,:,:)
    alpha_vec(i,1:2, 3:4) = pauli_vec(i,:,:)
  end do

  ! radial => gaussian_radial  
  ! radial => constant_radial  
  radial => hydrogenic_u_bispinor
  ! call radial((/0.0_dp, 0.0_dp, 0.0_dp /),bispinor)

  call read_in_grid("lebedev_he.dat", grid)
  
  initial_state%norms = check_spinor_normalisation(initial_state, grid)
  write(*,*) "Initial state norms for each m value = ", initial_state%norms
  final_state%norms = check_spinor_normalisation(final_state, grid)
  write(*,*) "Final state norms for each m value = ", final_state%norms

  transition_energy =  -hydrogenic_dirac_energy(final_state)&
  & + hydrogenic_dirac_energy(initial_state)

  final_average = 0.0_dp
  do m1 = 1, size(initial_state%m)
    do m2 = 1, size(final_state%m)
      if (abs(initial_state%m(m1) - final_state%m(m2)) > 1) then
        write(*,*) "magnetic transition m1 = ", initial_state%m(m1) , " to m2 = ", &
        & final_state%m(m2) ," is not allowed for dipole transition"
        cycle
      end if
      integrand = 0.0_dp
      !$omp parallel do reduction(+:integrand) shared(alpha_vec,initial_state,final_state,m1,m2, &
      !$omp transition_energy, grid)&
      !$omp private(dirac_spinor1,dirac_spinor2, wvfn_product, r)  default(none) 
      do i = 1, size(grid%x)
        r = (/ grid%x(i), grid%y(i), grid%z(i) /)
        call hydrogenic_dirac_spinor(r, initial_state, initial_state%m(m1), dirac_spinor1)
        call hydrogenic_dirac_spinor(r, final_state, final_state%m(m2),  dirac_spinor2)

        wvfn_product(1) = dot_product((dirac_spinor2), matmul(alpha_vec(1,:,:), dirac_spinor1)) 
        wvfn_product(2) = dot_product((dirac_spinor2), matmul(alpha_vec(2,:,:), dirac_spinor1))
        wvfn_product(3) = dot_product((dirac_spinor2), matmul(alpha_vec(3,:,:), dirac_spinor1))
        wvfn_product = wvfn_product * grid%weight(i)
        integrand = integrand + dot_product(wvfn_product, wvfn_product)
        ! write(71,*) r(1), r(3), real(dot_product(wvfn_product, wvfn_product))
      end do

      write(*,'(A, F4.1, A, F4.1, ES15.5)') "Total rate for m1 =  ", initial_state%m(m1), " m2 = ",final_state%m(m2),  &
      &integrand  *(transition_energy) * second / speed_of_light * 4.0_dp / 3.0_dp
      final_average = final_average + integrand  *(transition_energy) / speed_of_light *second * 4.0_dp / 3.0_dp
    end do
  end do

  write(*,'(A, ES15.5)') "Final total rate = ", final_average / (2.0_dp * initial_state%j + 1.0_dp)

contains

  ! Returns <psi|psi> on the grid
  function check_spinor_normalisation(state, grid) result(norm)
    implicit none

    type(dirac_state), intent(in) :: state
    type(grid_type), intent(in) :: grid
    real(kind=dp), dimension(:), allocatable   :: norm

    complex(kind=dp), dimension(4) :: spinor 
    integer :: m1, i, j, k
    
    allocate(norm(size(state%m)))
    ! Lets check normalisation of our states
    norm = 0.0_dp
    do m1 = 1, size(state%m)
      !$omp parallel do default(none) shared(state,m1,grid)&
      !$omp private(spinor, r)  reduction(+:norm) 
      do i = 1, size(grid%x)
        ! r = (/-grid%x_min+real(i,dp)*grid%dx, -grid%x_min+real(j,dp)*grid%dx,  -grid%x_min+real(k,dp)*grid%dx/)
        r = (/ grid%x(i), grid%y(i), grid%z(i) /)
        ! r = (/-grid%x_min+real(i,dp)*grid%dx, 0.0_dp,  -grid%x_min+real(k,dp)*grid%dx/)
        call hydrogenic_dirac_spinor(r, state, state%m(m1), spinor)

        norm(m1)= norm(m1) + dot_product((spinor),spinor) * grid%weight(i)
        ! write(28,*) r(1), r(3), real(dot_product(spinor, spinor))
      end do
    end do
  end function
  ! Setter for the quantum numbers, and also calculates k
  subroutine set_quantum_numbers(state, n, l, j, mass, Z)
    implicit none

    type(dirac_state), intent(out) :: state

    integer,           intent(in)  :: n, l, Z
    real(kind=dp),     intent(in)  :: j, mass

    integer :: i
    state%n = n
    state%l = l
    state%j = j
    state%Z = Z
    state%mass = mass

    if(state%j == (state%l - 0.5_dp)) then
      state%k = state%l
    else if (state%j == (state%l + 0.5_dp)) then
      state%k = -(state%l + 1)
    else
      stop "J must be equal to l += 0.5"
    end if

    allocate(state%m(int(2*state%j + 1)))
    allocate(state%norms(int(2*state%j + 1)))
    do i = 1, int(2*state%j) + 1
      state%m(i) = -state%j + real(i-1,dp)
    end do

    state%red_mass = mass * real(Z,dp) * amu / electron_mass / (mass + real(Z, dp) * amu / electron_mass)
    state%E = hydrogenic_dirac_energy(state)
    write(*,*) "E = ", state%E, " for n  = ", state%n, ", l = ", state%l, ", j = ", state%j, ", k = ", state%k, " redmass = ",&
    state%red_mass
  end subroutine set_quantum_numbers

  ! Hydrogenlike Dirac spinor, including both major and
  ! minor components
  subroutine hydrogenic_dirac_spinor(r, state, m, dirac_spinor)
    implicit none

    real(kind=dp), dimension(3), intent(in)     :: r
    type(dirac_state), intent(in)               :: state
    real(kind=dp), intent(in)                   :: m
    complex(kind=dp), dimension(4), intent(out) :: dirac_spinor 
    real(kind=dp), dimension(2)   :: radial_spinor

    real(kind=dp), dimension(2) :: u_bispinor

    real(kind=dp), dimension(2,2) :: prefactor_matrix
    real(kind=dp)                 :: gamma_k, norm_r
    integer                       :: sgn_k

    complex(kind=dp), dimension(2)                 :: weyl_spinor

    call radial(r, state, u_bispinor)

    gamma_k = sqrt(state%k**2 - fine_structure**2)
    prefactor_matrix(1,:) = (/ fine_structure, -(state%k - gamma_k) /)
    prefactor_matrix(2,:) = (/ -(state%k - gamma_k), fine_structure /)

    radial_spinor = 1.0_dp/(sqrt(2.0_dp*state%k*(state%k-gamma_k))) * matmul(prefactor_matrix, u_bispinor)

    call spin_spherical_harmonic(state%l, state%k, real(m,dp), r, weyl_spinor)

    norm_r = norm2(r)
    if(norm_r == 0.0_dp) then
      norm_r = tiny(1.0_dp)
    end if
    ! We now need to multiply these by the appropriate 1/r, and the spin spherical harmonics
    ! dirac_spinor(1:2) = radial_spinor(1) / norm_r * weyl_spinor
    dirac_spinor(1:2) = radial_spinor(1) * weyl_spinor  

    if(state%k > 0) then
      sgn_k = 1
    else
      sgn_k = -1
    end if
    call spin_spherical_harmonic(state%l-sgn_k, -state%k, real(m,dp), r, weyl_spinor)
    ! dirac_spinor(3:4) = -complex(0.0_dp, 1.0_dp) * radial_spinor(2) / norm_r * weyl_spinor
    dirac_spinor(3:4) = -complex(0.0_dp, 1.0_dp) * radial_spinor(2)  * weyl_spinor 

  end subroutine

  subroutine gaussian_radial(r, state, f)
    implicit none

    real(kind=dp), intent(in), dimension(3) :: r
    type(dirac_state), intent(in)           :: state
    real(kind=dp), intent(out), dimension(2) :: f

    f = exp(-norm2(r)**2/2.0_dp)/sqrt(gamma(1.5_dp))
  end subroutine
  subroutine constant_radial(r, state, f)
    implicit none

    real(kind=dp), intent(in), dimension(3) :: r
    type(dirac_state), intent(in)           :: state
    real(kind=dp), intent(out), dimension(2) :: f

    if(norm2(r) < 1.0_dp) then
      f = 1.0_dp
    else
      f = 0.0_dp
    end if
  end subroutine
  ! Hydrogenlike bispinor u
  subroutine hydrogenic_u_bispinor(r, state, u_bispinor)
    implicit none

    real(kind=dp), dimension(3), intent(in)  :: r
    type(dirac_state), intent(in)            :: state
    real(kind=dp), dimension(2), intent(out) :: u_bispinor

    ! sqrt(m^2 - E^2). Equivalent to K in mudirac
    real(kind=dp) :: C_nk

    ! sqrt(k^2 - alpha^2)
    real(kind=dp) :: gamma_k

    ! n - |k|
    integer       :: n_r, k

    real(kind=dp) :: A_plus_nk, A_minus_nk

    real(kind=dp) :: vector_prefactor, norm_r, red_mass

    ! Absolute value of r
    norm_r = norm2(r)

    if(norm_r == 0.0_dp) then
      norm_r = tiny(1.0_dp)
    end if

    n_r = state%n - abs(state%k)
    k = state%k

    C_nk    = sqrt(state%red_mass**2*speed_of_light**2 - state%E**2/speed_of_light**2)
    gamma_k = sqrt(state%k**2 - fine_structure**2)
    
    if(n_r == 0) then
      A_plus_nk = 0.0_dp
    else
      A_plus_nk = sqrt((C_nk * factorial(n_r - 1))/(gamma(n_r + 2.0_dp * gamma_k + 1.0_dp))) *&
                & sqrt((n_r + gamma_k + gamma_k/k * sqrt(fine_structure**2 + (n_r + gamma_k)**2))/(2.0_dp * gamma_k**2/k**2 * &
                & (fine_structure**2 + (n_r + gamma_k)**2))) 
    end if

    A_minus_nk = sqrt((C_nk * factorial(n_r))/(gamma(n_r + 2.0_dp * gamma_k))) *&
              & sqrt((n_r + gamma_k - gamma_k/k * sqrt(fine_structure**2 + (n_r + gamma_k)**2))/(2.0_dp * gamma_k**2/k**2 * &
              & (fine_structure**2 + (n_r + gamma_k)**2))) 

    vector_prefactor = (2.0_dp * C_nk * norm_r)**gamma_k * exp(-C_nk * norm_r)
    ! write(*,*) norm_r, vector_prefactor, A_plus_nk, A_minus_nk, C_nk, gamma_k

    if(A_plus_nk == 0.0_dp) then

      u_bispinor(1) = 0.0_dp
    else
      u_bispinor(1) = vector_prefactor * A_plus_nk * 2.0_dp * C_nk * norm_r * gen_laguerre_poly(n_r - 1, 2.0_dp*gamma_k+1.0_dp,&
                    & 2.0_dp*C_nk*norm_r) /norm_r
    end if
    ! write(*,*) gen_laguerre_poly(n_r, 2.0_dp*gamma_k-1.0_dp,2.0_dp*C_nk*norm_r)
    u_bispinor(2) = vector_prefactor * A_minus_nk * gen_laguerre_poly(n_r, 2.0_dp*gamma_k-1.0_dp,&
                  & 2.0_dp*C_nk*norm_r) /norm_r

  end subroutine

  ! Multiplies correct spherical harmonics by the appropriate
  ! CG coefficients to form the spin spherical harmonics, which are
  ! a 2-vector
  subroutine spin_spherical_harmonic(l, k, m, x, spin_harmonic) 
    implicit none

    integer,                        intent(in)  :: l, k
    real(kind=dp),                  intent(in)  :: m
    real(kind=dp), dimension(3),    intent(in)  :: x
    complex(kind=dp), dimension(2), intent(out) :: spin_harmonic


    ! Up spin
    if(abs(int(m-0.5_dp)) > l ) then
      spin_harmonic(1) = 0.0_dp
    else
      spin_harmonic(1) = clebsch_gordon_spin_spherical_harmonics(k,m,1) * SphericalYCartesian(l, int(m-0.5), x)
    end if

    ! Down spin
    if(abs(int(m+0.5_dp)) > l ) then
      spin_harmonic(2) = 0.0_dp
    else
      spin_harmonic(2) = clebsch_gordon_spin_spherical_harmonics(k,m,0) * SphericalYCartesian(l, int(m+0.5), x)
    end if


  end subroutine

  ! Take in the Dirac quantum number, magnetic number, and particle spin
  ! to calculate the CG coefficient used to couple together the orbital
  ! and spin angular momentum of the particle
  function clebsch_gordon_spin_spherical_harmonics(k, m, s) result(coeff)
    implicit none

    integer,       intent(in) :: k, s
    real(kind=dp), intent(in) :: m

    complex(kind=dp) :: coeff
    integer       :: sgn_k
    complex(kind=dp) :: prefactor, denom

    if (k > 0) then
      sgn_k = 1
    else
      sgn_k = -1
    end if

    denom = complex(2.0_dp * real(k,dp) + 1.0_dp,0.0_dp)
    ! Spin up
    if (s == 1) then
      prefactor = complex(real(k,dp) + 0.5_dp -m, 0.0_dp)
      coeff =  -sgn_k * sqrt(prefactor)/sqrt(denom)

    ! Spin down
    else
      prefactor = complex(real(k,dp) + 0.5_dp +m, 0.0_dp)
      coeff = sqrt(prefactor)/sqrt(denom)
    end if



  end function clebsch_gordon_spin_spherical_harmonics

  function SolidRCartesian(l, m, x)

    complex(dp) :: SolidRCartesian
    integer, intent(in) :: l, m
    real(dp), intent(in) :: x(3)
    integer :: p, q, s

    SolidRCartesian = CPLX_ZERO

    do p = 0, l
      q = p - m
      s = l - p - q

      if ((q >= 0) .and. (s >= 0)) then
         SolidRCartesian = SolidRCartesian + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                           * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                           * (x(3)**s) &
                                           / (factorial(p) * factorial(q) * factorial(s)))
      end if
    end do

    SolidRCartesian = SolidRCartesian * sqrt(factorial(l + m) * factorial(l - m))

  end function SolidRCartesian
  function SphericalYCartesian(l, m, x)

    complex(dp) :: SphericalYCartesian
    integer, intent(in) :: l, m
    real(dp), intent(in) :: x(3)

    if(vector_normsq(x) /= 0.0_dp) then
      SphericalYCartesian = SolidRCartesian(l, m, x) * sqrt(((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                                                    * (vector_normsq(x)**(-0.5_dp * l))
    else
      SphericalYCartesian = 0.0_dp
    end if


  end function SphericalYCartesian
  elemental function factorial(n) result(res)

    ! factorial_real

    integer, intent(in) :: n
    real(dp)            :: res
    integer :: i

    if (n<0) then
      res = huge(1.0_dp)
      res = res + 1.0_dp
    else if(n <= 20) then
      res = factorial_table(n)
    else
      res=1.0_dp
      do i=2,n
        res = res*i
      end do
    end if

  end function factorial
  pure function vector_normsq(vector) result(normsq)

    real(dp), intent(in), dimension(:) :: vector
    real(dp)             :: normsq

    normsq = dot_product(vector,vector)

  end function vector_normsq
  function hydrogenic_dirac_energy(ds) result(energy)
    implicit none

    type(dirac_state), intent(in)  :: ds
    real(kind=dp)                  :: energy

    energy = ds%red_mass * speed_of_light **2 / &
             & (sqrt(1.0_dp + ((ds%Z * fine_structure)/(ds%n - abs(ds%k) + sqrt(ds%k**2-ds%Z**2*fine_structure**2)))**2))
    !- mu * speed_of_light**2
  end function hydrogenic_dirac_energy
  recursive function gen_laguerre_poly(n, alpha, x) result(val)

    integer,       intent(in)     :: n
    real(kind=dp), intent(in)     :: x, alpha
    real(kind=dp)                 :: val

    if (n < 0) then
      stop "n for generalised Laguerre polynomial must be greater than or equal to 0"
    end if

    select case (n)

      case (0)
        val = 1.0_dp
        return

      case (1)
        val =  1.0_dp + alpha  - x
        return

      case (2)
        val = 0.5_dp * (x**2 - 2.0_dp * x * (alpha + 2.0_dp) + (alpha + 1.0_dp) * (alpha + 2.0_dp))
        return

      case (3)
        val = 1.0_dp / 6.0_dp * (- x**3 + 3.0_dp * x**2 * (alpha + 3.0_dp) - 3.0_dp * x * (alpha + 2.0_dp) * (alpha + 3.0_dp)&
        & + (alpha + 1.0_dp) * (alpha + 2.0_dp) * (alpha + 3.0_dp))
        return

      case default
        val = ((2.0_dp * real(n,dp) - 1.0_dp + alpha - x) * gen_laguerre_poly(n - 1, alpha, x) - (real(n,dp) - 1.0_dp + alpha)&
        & * gen_laguerre_poly(n-2, alpha, x))/real(n,dp)
          
    end select
  end function gen_laguerre_poly

  function spherical_bessel_zero(k, r) result(func)
    implicit none

    real(kind=dp), intent(in)  :: k, r
    real(kind=dp) :: func

    if(k*r == 0.0_dp) then
      func = 1.0_dp
    else
      func = sin(k*r) / (k*r)
    end if

  end function

  ! Reads in a Lebedev grid with positions and weights
  subroutine read_in_grid(filename, grid)
    implicit none

    character(len=*), intent(in)   :: filename
    type(grid_type), intent(inout) :: grid

    integer :: i, istat, file_length
    real(kind=dp) :: val

    open(30, file=filename)

    file_length = 0
    do 
      read(30,*, iostat=istat) val
      if(istat /= 0) then
        exit
      else
        file_length = file_length + 1
      end if
    end do

    rewind(30)
    allocate(grid%x(file_length))
    allocate(grid%y(file_length))
    allocate(grid%z(file_length))
    allocate(grid%weight(file_length))

    do i = 1, file_length
      read(30, *) grid%x(i), grid%y(i), grid%z(i), grid%weight(i)
    end do
    grid%x = grid%x
    grid%y = grid%y
    grid%z = grid%z
    grid%weight = grid%weight
  end subroutine
end program dirac
