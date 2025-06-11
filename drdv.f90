module dvdr_constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nrad = 400
  integer, parameter :: na = nrad * 3
end module dvdr_constants

subroutine dvdr(r, rv, z, nmt, r2dvdr)
  use dvdr_constants
  implicit none
  ! Dummy arguments
  real(dp), intent(in) :: r(nrad), z
  integer, intent(in) :: nmt
  real(dp), intent(inout) :: rv(nrad)
  real(dp), intent(out) :: r2dvdr(nrad)
  ! Local variables
  real(dp), allocatable :: r2v(:), a(:)
  integer :: i, ia_err

  ! Allocate arrays
  allocate(r2v(nmt), a(na), stat=ia_err)
  if (ia_err /= 0) stop 'dvdr: allocation failed'

  ! Compute r2v = r * rv
  do i = 1, nmt
    r2v(i) = r(i) * rv(i)
  end do

  ! Calculate spline derivative
  call derspl(nmt, r, r2v, r2dvdr, a)

  ! Compute r2dvdr = d(r2v)/dr - 2*rv and update rv = rv/r
  do i = 1, nmt
    r2dvdr(i) = r2dvdr(i) - 2.0_dp * rv(i)
    rv(i) = rv(i) / r(i)
  end do

  ! Zero out r2dvdr for i > nmt
  r2dvdr(nmt+1:nrad) = 0.0_dp

  ! Clean up
  deallocate(r2v, a)
end subroutine dvdr

subroutine derspl(n, x, f, d, a)
  use dvdr_constants
  implicit none
  ! Dummy arguments
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n), f(n)
  real(dp), intent(out) :: d(n)
  real(dp), intent(inout) :: a(na)
  ! Local variables
  real(dp) :: h1, h2, p
  integer :: i, j, k

  ! Check for monotonicity
  do i = 2, n
    if (x(i) <= x(i-1)) then
      write(6, '(a,i3,a)') 'return from derspl  ', i, ' out of order'
      a(1) = 1.0_dp
      return
    end if
  end do

  ! Set up tridiagonal system
  do i = 1, n
    if (i == 1 .or. i == n) then
      j = merge(2, n-1, i == 1)
      h1 = 1.0_dp / (x(j) - x(j-1))
      h2 = 1.0_dp / (x(j+1) - x(j))
      a(3*i-2) = h1 * h1
      a(3*i-1) = h1 * h1 - h2 * h2
      a(3*i) = -h2 * h2
      d(i) = 2.0_dp * (f(j) * (h2 * h2 * h2 + h1 * h1 * h1) - &
                       f(j+1) * h2 * h2 * h2 - f(j-1) * h1 * h1 * h1)
    else
      h1 = 1.0_dp / (x(i) - x(i-1))
      h2 = 1.0_dp / (x(i+1) - x(i))
      a(3*i-2) = h1
      a(3*i-1) = 2.0_dp * (h1 + h2)
      a(3*i) = h2
      d(i) = 3.0_dp * (f(i+1) * h2 * h2 + f(i) * (h1 * h1 - h2 * h2) - &
                       f(i-1) * h1 * h1)
    end if
  end do

  ! Forward elimination
  p = a(4) / a(1)
  a(5) = a(5) - p * a(2)
  a(6) = a(6) - p * a(3)
  d(2) = d(2) - p * d(1)
  do i = 3, n
    k = 3 * i - 4
    p = a(k+2) / a(k)
    a(k+3) = a(k+3) - p * a(k+1)
    d(i) = d(i) - p * d(i-1)
    if (i == n-1) then
      p = a(k+5) / a(k)
      a(k+5) = a(k+6) - p * a(k+1)
      a(k+6) = a(k+7)
      d(n) = d(n) - p * d(n-2)
    end if
  end do

  ! Back substitution
  d(n) = d(n) / a(3*n-1)
  do i = 3, n
    j = n + 2 - i
    d(j) = (d(j) - a(3*j) * d(j+1)) / a(3*j-1)
  end do
  d(1) = (d(1) - d(2) * a(2) - d(3) * a(3)) / a(1)
  a(1) = 0.0_dp
end subroutine derspl
