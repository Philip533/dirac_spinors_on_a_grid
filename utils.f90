module utils
  use constants
  implicit none

contains

  subroutine write_transition_matrix(transition_matrix)
    implicit none

    type(t_matrix),    intent(in)              :: transition_matrix


    integer :: i, j

    write(*,'(A)') '#####################################################'
    write(*,'(A,I2, A, I2)') '# TransitionMatrix from state with k = ',transition_matrix%k1," to state with k = &
                             & ", transition_matrix%k2
    write(*,'(A, ES11.5, A)') '# Total rate = ', transition_matrix%total_rate, " s^-1"
    write(*,'(A)') '#####################################################'
    write(*,'(15X)', advance="no")
    do i = 0, size(transition_matrix%m1)
      
      do j = 1, size(transition_matrix%m2)
        if(i == 0) then
          if(j == size(transition_matrix%m2)) then
            write(*,'(F4.1)',advance="yes")  transition_matrix%m2(j)
          else
            write(*,'(F4.1, 11X)',advance="no") transition_matrix%m2(j)
          end if
        end if
      end do

        if(i /= 0) then
          do j = 0, size(transition_matrix%m2)
            if(j==0) then
              write(*,'(F4.1, 4X)', advance="no") transition_matrix%m1(i)
            else if(j/=0 .and. j /= size(transition_matrix%m2)) then

              if(transition_matrix%T(i,j) == 0) then
                write(*,'(10X,I1, 4X)', advance="no") int(transition_matrix%T(i,j))
              else
                write(*,'(ES11.5, 4X)', advance="no") transition_matrix%T(i,j)
              end if
            else if(j== size(transition_matrix%m2)) then
              if(transition_matrix%T(i,j) == 0) then
                write(*,'(10X,I1, 4X)', advance="yes") int(transition_matrix%T(i,j))
              else
                write(*,'(ES11.5, 4X)', advance="yes") transition_matrix%T(i,j)
              end if
            end if

          end do
      end if 

    end do
  end subroutine write_transition_matrix

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
  end subroutine
  subroutine set_schrodinger_quantum_numbers(state, n, l, mass, Z)
    implicit none

    type(schrodinger_state), intent(out) :: state

    integer,           intent(in)  :: n, l, Z
    real(kind=dp),     intent(in)  :: mass

    integer :: i

    state%n = n
    state%l = l
    state%Z = Z
    state%mass = mass

    allocate(state%norms(int(2*state%l + 1)))
    allocate(state%m(int(2*state%l + 1)))
    do i = 1, int(2*state%l) + 1
      state%m(i) = -state%l + real(i-1,dp)
    end do

    state%red_mass = mass * real(Z,dp) * amu / electron_mass / (mass + real(Z, dp) * amu / electron_mass)
  end subroutine set_schrodinger_quantum_numbers

  subroutine set_dirac_quantum_numbers(state, n, l, s, mass, Z)
    implicit none

    type(dirac_state), intent(out) :: state
    integer,           intent(in)  :: n, l, Z
    real(kind=dp),     intent(in)  :: mass
    logical,           intent(in)  :: s

    integer :: i

    state%n = n
    state%l = l
    state%Z = Z
    state%mass = mass

    if(s) then
      state%j = real(state%l,dp) + 0.5_dp
      state%k = -(state%l + 1)
    else
      state%j = real(state%l,dp) - 0.5_dp
      state%k = state%l
    end if

    allocate(state%norms(int(2*state%j + 1)))
    allocate(state%m(int(2*state%j + 1)))
    do i = 1, int(2*state%j) + 1
      state%m(i) = -state%j + real(i-1,dp)
    end do
    state%red_mass = mass * real(Z,dp) * amu / electron_mass / (mass + real(Z, dp) * amu / electron_mass)
  end subroutine
  subroutine set_alpha_vector(alpha_vec)
    implicit none

    complex(kind=dp), dimension(3,4,4), intent(out) :: alpha_vec

    complex(kind=dp), dimension(3,2,2)              :: pauli_vec

    integer :: i
    ! Initialise Pauli matrices in the vector
    pauli_vec(1,:,:) = reshape((/ 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp /), shape(pauli_vec(1,:,:)))
    pauli_vec(2,:,:) = reshape((/ (0.0_dp,0.0_dp), (0.0_dp,1.0_dp), (0.0_dp,-1.0_dp), (0.0_dp,0.0_dp) /), shape(pauli_vec(2,:,:)))
    pauli_vec(3,:,:) = reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp /), shape(pauli_vec(3,:,:)))

    ! Define the \alpha vector
    do i = 1, 3
      alpha_vec(i,1:2, 1:2) = 0.0_dp
      alpha_vec(i,3:4, 3:4) = 0.0_dp
      alpha_vec(i,3:4, 1:2) = pauli_vec(i,:,:)
      alpha_vec(i,1:2, 3:4) = pauli_vec(i,:,:)
    end do
  end subroutine
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
  pure function vector_normsq(vector) result(normsq)

    real(dp), intent(in), dimension(:) :: vector
    real(dp)             :: normsq

    normsq = dot_product(vector,vector)

  end function vector_normsq
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
  function split_string(string,sep,column)
    implicit none

    character(len=*), intent(in)   :: string
    character(len=1), intent(in)   :: sep
    integer, intent(in)            :: column
    character(len=:), allocatable  :: split_string

    character(len=:), allocatable :: temp
    character(len=:), allocatable :: work

    integer :: i, cnt, lastpos
    
    cnt = 0
    lastpos = 0
    split_string = ""
   
    temp = string//sep
    do i = 1, len_trim(temp)
      if(temp(i:i) == sep)then
        cnt = cnt + 1
        if(cnt == column) then
          work = temp(lastpos+1:i-1)
          exit
        end if
        lastpos = i
      end if
    end do
    split_string = work
  end function split_string
    subroutine iupac_to_atomic(istate, n, l, s)
    implicit none

    character(len=*), intent(in)  :: istate
    integer,          intent(out) :: n, l
    logical,          intent(out) :: s
    integer                       :: length
    character(len=1)              :: shell
    character(len=2)              :: trimstate
    integer                       :: orbit

    length = len(istate)

    ! Catch the edge case for K
    if (length < 2) then
      if (istate == "K") then
        n = 1
        l = 0
        s = .false.
      else
        stop "Invalid IUPAC state"
      end if
    end if

    ! Trim the whitespace and split the IUPAC into letter and number
    trimstate = trim(adjustl(istate))
    shell = trimstate(1:1)
    read(trimstate(2:2),*)orbit

    ! Calculate the quantum numbers
    n = ichar(shell) - 74
    l = orbit/2
    s = (mod(orbit,2)) == 1
  end subroutine
end module utils
