program dirac
  use constants
  use utils
  use hydrogenic
  implicit none

  integer :: i, j, k, nuc_z, N, m1, m2, n1, l1, n2, l2

  ! Define the Dirac \alpha vector
  complex(kind=dp), dimension(3, 4, 4) :: alpha_vec

  ! Dirac spinor
  complex(kind=dp), dimension(4)  :: dirac_spinor1, dirac_spinor2
  real(kind=dp) :: dx, transition_energy, red_mass, x_min, final_average
  complex(kind=dp) :: integrand

  type(schrodinger_state) :: initial_state_schro, final_state_schro
  type(dirac_state) :: initial_state_dir, final_state_dir
  type(grid_type) :: grid
  real(kind=dp), dimension(3) :: r
  integer       :: nuclear_Z, num_args
  character(len=12), dimension(:), allocatable :: args
  character(len=7)                             :: iupac_transition
  character(len=3)                             :: state1, state2
  logical                                      :: s1, s2

  ! Interface to the routine for calculating the radial part R(r)
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

  ! Transition matrix
  type(t_matrix) :: transition_matrix


  ! First argument is Z.
  ! Second argument is the transition in IUPAC notation
  num_args = command_argument_count()
  allocate(args(num_args))
  do i = 1, num_args
    call get_command_argument(i, args(i))
  end do

  read(args(1), '(I2)') nuclear_Z
  iupac_transition = args(2)
  state1 = split_string(iupac_transition,"-", 2)
  state2 = split_string(iupac_transition,"-", 1)
  call iupac_to_atomic(state1, n1, l1, s1)
  call iupac_to_atomic(state2, n2, l2, s2)

  if(l1 >= n1) then
    stop "l1 cannot be greater than n1"
  else if(l2 >= n2) then
    stop "l1 cannot be greater than n1"
  end if

  ! Initialise the Dirac \alpha vector
  call set_alpha_vector(alpha_vec)

  ! Choose the radial function
  ! radial => hydrogenic_u_bispinor
  radial => mudirac_bispinor

  ! Read in the atomic grid
  call read_in_grid("lebedev_h.dat", grid)
  
  call set_dirac_quantum_numbers(initial_state_dir, n1, l1, s1, particle_mass, nuclear_Z)
  call set_dirac_quantum_numbers(final_state_dir, n2, l2, s2, particle_mass, nuclear_Z)
  initial_state_dir%E = hydrogenic_dirac_energy(initial_state_dir)
  final_state_dir%E = hydrogenic_dirac_energy(final_state_dir)

  ! First dimension is the first state, second dimension is second state
  allocate(transition_matrix%T(size(initial_state_dir%m), size(final_state_dir%m)))
  transition_matrix%k1 = initial_state_dir%k; transition_matrix%k2 = final_state_dir%k
  transition_matrix%m1 = initial_state_dir%m; transition_matrix%m2 = final_state_dir%m

  ! Check if we have dipole or quadrupole
  if (abs(l2 - l1) == 1) then

    call calculate_dirac_dipole(initial_state_dir, final_state_dir, transition_matrix)

  else if (abs(l2 - l1) == 2 .or. abs(l2 - l1) == 0) then

    call calculate_dirac_quadrupole(initial_state_dir, final_state_dir, transition_matrix)

  end if

  ! We have a total rate that needs to be averaged over the number of starting m states
  transition_matrix%total_rate = transition_matrix%total_rate/(size(initial_state_dir%m))

  call write_transition_matrix(transition_matrix)

  deallocate(grid%x, grid%y, grid%z, grid%weight)

contains

  subroutine calculate_dirac_quadrupole(initial_state, final_state, transition_matrix)
    implicit none

    type(dirac_state),                          intent(in)    :: initial_state, final_state
    type(t_matrix),                             intent(inout) :: transition_matrix

    real(kind=dp) :: transition_energy, integrand, rate
    complex(kind=dp),dimension(3)   :: wvfn_product
    complex(kind=dp) :: xax, yax, zax, xay, yay, zay, xaz, yaz, zaz

    transition_energy =  -hydrogenic_dirac_energy(final_state)&
    & + hydrogenic_dirac_energy(initial_state)
    rate = 0.0_dp
    do m1 = 1, size(initial_state%m)
      do m2 = 1, size(final_state%m)
        xax = 0.0_dp; yax = 0.0_dp; zax = 0.0_dp
        xay = 0.0_dp; yay = 0.0_dp; zay = 0.0_dp
        xaz = 0.0_dp; yaz = 0.0_dp; zaz = 0.0_dp
        if (abs(initial_state%m(m1) - final_state%m(m2)) > 2) then
          transition_matrix%T(m1, m2) = 0.0_dp
          cycle
        end if
        integrand = 0.0_dp
        do i = 1, size(grid%x)

          ! Put the position into a 3-vector
          r = (/ grid%x(i), grid%y(i), grid%z(i) /)

          ! Evaluate the wavefunctions for the initial and final states
          call hydrogenic_dirac_spinor(r, initial_state, initial_state%m(m1), dirac_spinor1)
          call hydrogenic_dirac_spinor(r, final_state, final_state%m(m2),  dirac_spinor2)

          xax = xax + dot_product(dirac_spinor2, matmul(alpha_vec(1,:,:), dirac_spinor1)) * r(1) * grid%weight(i)
          yax = yax + dot_product(dirac_spinor2, matmul(alpha_vec(1,:,:), dirac_spinor1)) * r(2) * grid%weight(i)
          zax = zax + dot_product(dirac_spinor2, matmul(alpha_vec(1,:,:), dirac_spinor1)) * r(3) * grid%weight(i)

          xay = xay + dot_product(dirac_spinor2, matmul(alpha_vec(2,:,:), dirac_spinor1)) * r(1) * grid%weight(i)
          yay = yay + dot_product(dirac_spinor2, matmul(alpha_vec(2,:,:), dirac_spinor1)) * r(2) * grid%weight(i)
          zay = zay + dot_product(dirac_spinor2, matmul(alpha_vec(2,:,:), dirac_spinor1)) * r(3) * grid%weight(i)

          xaz = xaz + dot_product(dirac_spinor2, matmul(alpha_vec(3,:,:), dirac_spinor1)) * r(1) * grid%weight(i)
          yaz = yaz + dot_product(dirac_spinor2, matmul(alpha_vec(3,:,:), dirac_spinor1)) * r(2) * grid%weight(i)
          zaz = zaz + dot_product(dirac_spinor2, matmul(alpha_vec(3,:,:), dirac_spinor1)) * r(3) * grid%weight(i)

        end do

        integrand =  (conjg(xax) * xax + 2.0_dp * conjg(xay) * xay + 2.0_dp * conjg(xaz) * xaz  &
                  & - xay * yax + conjg(yay) * yay - xaz * zax &
                  & + 2.0_dp * (conjg(yax) * yax + conjg(yaz) * yaz + conjg(zax) * zax) &
                  & - yaz * zay + 2.0_dp * conjg(zay) * zay - yay * zaz &
                  & + conjg(zaz) * zaz - xax * yay -  xax * zaz)
        integrand = 0.5_dp*(4.0_dp * xay * conjg(xay) + 4.0_dp * xaz * conjg(xaz) - conjg(xay) * yax &
                  & - xay * conjg(yax)+ 4.0_dp*yax*conjg(yax) - conjg(xax) * yay + 2.0_dp * yay * conjg(yay) &
                  & + 4.0_dp * yaz * conjg(yaz) - conjg(xaz) * zax - xaz * conjg(zax) &
                  & + 4.0_dp * zax * conjg(zax) - conjg(yaz) * zay - yaz * conjg(zay) &
                  & + 4.0_dp  * zay * conjg(zay) - conjg(xax) * zaz - conjg(yay) * zaz &
                  & + xax * (2.0_dp * conjg(xax) - conjg(yay) - conjg(zaz)) - yay * conjg(zaz) & 
                  & + 2.0_dp * zaz * conjg(zaz))
        rate = rate + integrand 
        transition_matrix%T(m1,m2) = integrand  *(transition_energy) / speed_of_light / (2.0_dp * pi) * (8.0_dp * pi / 15.0_dp)
        
        
      end do
    end do

    transition_matrix%total_rate = rate * transition_energy/speed_of_light / (2.0_dp * pi) * second * (8.0_dp * pi / 15.0_dp)
  end subroutine
  subroutine calculate_schrodinger_dipole(initial_state, final_state, rate)
    implicit none

    type(schrodinger_state), intent(in)  :: initial_state, final_state
    real(kind=dp),     intent(out) :: rate

    real(kind=dp) :: transition_energy
    complex(kind=dp) :: integrand, wvfn1, wvfn2
    complex(kind=dp),dimension(3)   :: wvfn_product

    integer :: m1, m2

    ! Calculation the transition energy
    transition_energy =  -hydrogenic_schro_energy(final_state)&
    & + hydrogenic_schro_energy(initial_state)
    
    final_average = 0.0_dp
    do m1 = 1, size(initial_state%m)
      do m2 = 1, size(final_state%m)
        if (abs(initial_state%m(m1) - final_state%m(m2)) > 1) then
          write(*,*) "magnetic transition m1 = ", initial_state%m(m1) , " to m2 = ", &
          & final_state%m(m2) ," is not allowed for dipole transition"
          cycle
        end if
        integrand = 0.0_dp
        do i = 1, size(grid%x)

          ! Put the position into a 3-vector
          r = (/ grid%x(i), grid%y(i), grid%z(i) /)

          wvfn1 = hydrogenic_schro_wvfn(norm2(r), initial_state)
          wvfn1 = wvfn1 * SphericalYCartesian(initial_state%l, int(initial_state%m(m1)), r)

          wvfn2 = hydrogenic_schro_wvfn(norm2(r), final_state)
          wvfn2 = wvfn2 * SphericalYCartesian(final_state%l, int(final_state%m(m2)), r)

          wvfn_product = conjg(wvfn1) * wvfn2 *  r

          integrand = integrand + sum(wvfn_product) * grid%weight(i)

        end do

        write(*,'(A, F4.1, A, F4.1, ES15.5)') "Total rate for m1 =  ", initial_state%m(m1), " m2 = ",final_state%m(m2),  &
        &integrand*conjg(integrand)  *(transition_energy) * second / speed_of_light * 4.0_dp / 3.0_dp
        rate = rate + integrand*conjg(integrand)  *(transition_energy)**3 * 4.0_dp * fine_structure**3 / (3.0_dp) * second
      end do
    end do
  end subroutine
  subroutine calculate_dirac_dipole(initial_state, final_state, transition_matrix)
    implicit none

    type(dirac_state),   intent(in)  :: initial_state, final_state
    type(t_matrix),    intent(inout) :: transition_matrix

    real(kind=dp) :: transition_energy, final_average
    complex(kind=dp) :: integrand
    complex(kind=dp),dimension(3)   :: wvfn_product

    integer :: m1, m2

    ! Calculation the transition energy
    transition_energy =  -hydrogenic_dirac_energy(final_state)&
    & + hydrogenic_dirac_energy(initial_state)
    
    final_average = 0.0_dp
    do m1 = 1, size(initial_state%m)
      do m2 = 1, size(final_state%m)
        if (abs(initial_state%m(m1) - final_state%m(m2)) > 1) then
          cycle
        end if
        integrand = 0.0_dp
        do i = 1, size(grid%x)

          ! Put the position into a 3-vector
          r = (/ grid%x(i), grid%y(i), grid%z(i) /)

          ! Evaluate the wavefunctions for the initial and final states
          call hydrogenic_dirac_spinor(r, initial_state, initial_state%m(m1), dirac_spinor1)
          call hydrogenic_dirac_spinor(r, final_state, final_state%m(m2),  dirac_spinor2)

          ! Find the matrix element for each \alpha
          wvfn_product(1) = dot_product((dirac_spinor2), matmul(alpha_vec(1,:,:), dirac_spinor1)) 
          wvfn_product(2) = dot_product((dirac_spinor2), matmul(alpha_vec(2,:,:), dirac_spinor1))
          wvfn_product(3) = dot_product((dirac_spinor2), matmul(alpha_vec(3,:,:), dirac_spinor1))

          ! Sum up each matrix element with the Lebedev weight
          integrand = integrand + sum(wvfn_product) * grid%weight(i) &
          & * spherical_bessel_zero(transition_energy/speed_of_light, norm2(r))

        end do

        final_average = final_average + integrand * conjg(integrand) * transition_energy/speed_of_light*4.0_dp/3.0_dp*second
        transition_matrix%T(m1,m2) = integrand  *conjg(integrand) *transition_energy / speed_of_light * 4.0_dp / 3.0_dp
      end do
    end do
    transition_matrix%total_rate = final_average

  end subroutine
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
      end do
    end do
  end function
  ! Setter for the quantum numbers, and also calculates k


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

    gamma_k = sqrt(state%k**2 - state%Z**2*fine_structure**2)
    ! prefactor_matrix(1,:) = (/ state%Z*fine_structure, -(state%k - gamma_k) /)
    ! prefactor_matrix(2,:) = (/ -(state%k - gamma_k), state%Z*fine_structure /)

    ! radial_spinor = 1.0_dp/(sqrt(2.0_dp*state%k*(state%k-gamma_k))) * matmul(prefactor_matrix, u_bispinor)
    radial_spinor = u_bispinor

    call spin_spherical_harmonic(state%l, state%k, real(m,dp), r, weyl_spinor)

    norm_r = norm2(r)
    if(norm_r == 0.0_dp) then
      norm_r = tiny(1.0_dp)
    end if
    ! We now need to multiply these by the appropriate 1/r, and the spin spherical harmonics
    ! dirac_spinor(1:2) = radial_spinor(1) / norm_r * weyl_spinor
    dirac_spinor(1:2) = radial_spinor(1) * weyl_spinor  /norm_r

    if(state%k > 0) then
      sgn_k = 1
    else
      sgn_k = -1
    end if
    call spin_spherical_harmonic(state%l-sgn_k, -state%k, real(m,dp), r, weyl_spinor)
    ! dirac_spinor(3:4) = -complex(0.0_dp, 1.0_dp) * radial_spinor(2) / norm_r * weyl_spinor
    dirac_spinor(3:4) = complex(0.0_dp, 1.0_dp) * radial_spinor(2)  * weyl_spinor  / norm_r

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
    gamma_k = sqrt(state%k**2 - state%Z**2*fine_structure**2)
    
    if(n_r == 0) then
      A_plus_nk = 0.0_dp
    else
      A_plus_nk = sqrt((C_nk * factorial(n_r - 1))/(gamma(n_r + 2.0_dp * gamma_k + 1.0_dp))) *&
                & sqrt((n_r + gamma_k + gamma_k/k * sqrt(state%Z**2*fine_structure**2 + (n_r + gamma_k)**2))&
                /(2.0_dp * gamma_k**2/k**2 *(state%Z**2*fine_structure**2 + (n_r + gamma_k)**2))) 
      ! write(*,*) "APLUS = ", a_plus_nk 
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
                    & 2.0_dp*C_nk*norm_r) 
    end if
    u_bispinor(2) = vector_prefactor * A_minus_nk * gen_laguerre_poly(n_r, 2.0_dp*gamma_k-1.0_dp,&
                  & 2.0_dp*C_nk*norm_r) 

  end subroutine
  subroutine mudirac_bispinor(r, state, u_bispinor)
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

    real(kind=dp) :: A_plus_nk, A_minus_nk, rho, A, rhodep, E_k

    real(kind=dp) :: vector_prefactor, norm_r, red_mass, lagP, lagM

    ! Absolute value of r
    norm_r = norm2(r)

    if(norm_r == 0.0_dp) then
      norm_r = tiny(1.0_dp)
    end if

    n_r = state%n - abs(state%k)
    k = state%k

    C_nk    = sqrt(state%red_mass**2*speed_of_light**2 - state%E**2/speed_of_light**2)
    gamma_k = sqrt(state%k**2 - state%Z**2*fine_structure**2)
    
    rho = 2.0_dp * C_nk * norm_r
    rhodep = (rho)**gamma_k * exp(-0.5 * rho)


    if(n_r == 0) then
      A = sqrt(C_nk / (2.0_dp * state%n * (state%n + gamma_k) * gamma_k * gamma(2.0_dp * gamma_k)))
      u_bispinor(1) = A * (state%n + gamma_k) * rhodep 
      u_bispinor(2) = -A * state%Z * fine_structure * rhodep 
    else
      E_k = state%E * state%k / (gamma_k * state%red_mass * speed_of_light**2)
      A = sqrt(C_nk * factorial(state%n - abs(state%k) - 1)/(4.0_dp * state%k * (state%k - gamma_k) * &
      & (state%n-abs(state%k)+gamma_k)*gamma(state%n-abs(state%k) + 2.0_dp * gamma_k + 1.0_dp)) * (E_k + E_k**2))
      lagP = rho * gen_laguerre_poly(state%n - abs(state%k) - 1, 2.0_dp * gamma_k + 1.0_dp, rho)
      lagM = (gamma_k * state%red_mass * speed_of_light**2 - state%k*state%E) /(speed_of_light * C_nk) *&
      & gen_laguerre_poly(state%n - abs(state%k), 2.0_dp * gamma_k - 1.0_dp, rho)
      u_bispinor(1) = A*rhodep* (state%Z * fine_structure * lagP +(gamma_k - k) * lagM)
      u_bispinor(2) = -A*rhodep* (state%Z * fine_structure * lagM +(gamma_k - k) * lagP)
    end if
    ! write(*,*) norm_r, rhodep, C_nk, A, gamma_k

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
    if(m == (k+0.5_dp) ) then
      spin_harmonic(1) = 0.0_dp
    else
      spin_harmonic(1) = clebsch_gordon_spin_spherical_harmonics(k,m,1) * SphericalYCartesian(l, int(m-0.5), x)
      ! write(*,*)clebsch_gordon_spin_spherical_harmonics(k,m,1) , " Up", l, k, m
    end if

    ! Down spin
    if(m ==(-k - 0.5_dp)) then
      spin_harmonic(2) = 0.0_dp
    else
      spin_harmonic(2) = clebsch_gordon_spin_spherical_harmonics(k,m,0) * SphericalYCartesian(l, int(m+0.5), x)
      ! write(*,*)clebsch_gordon_spin_spherical_harmonics(k,m,2), " Down", l, k, m
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

  function hydrogenic_dirac_energy(ds) result(energy)
    implicit none

    type(dirac_state), intent(in)  :: ds
    real(kind=dp)                  :: energy

    energy = ds%red_mass * speed_of_light **2 / &
             & (sqrt(1.0_dp + ((ds%Z * fine_structure)/(ds%n - abs(ds%k) + sqrt(ds%k**2-ds%Z**2*fine_structure**2)))**2))
    !- mu * speed_of_light**2
  end function hydrogenic_dirac_energy

end program dirac
