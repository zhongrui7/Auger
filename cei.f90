module cei_constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nrad = 250
  real(dp), parameter :: tiny = 1.0e-20_dp
end module cei_constants

subroutine cei(f2, f4, f1, f3, r, n, lambda, l2, l4, l1, l3, dx, rws, result)
  use cei_constants
  use omp_lib
  implicit none
  ! Dummy arguments
  real(dp), intent(in) :: f1(n), f2(n), f3(n), f4(n), r(n)
  integer, intent(in) :: n, lambda
  real(dp), intent(in) :: l1, l2, l3, l4, dx, rws
  real(dp), intent(out) :: result
  ! Local variables
  real(dp), allocatable :: rl(:), rmlp1(:), arg1(:), arg2(:), arg3(:), arg4(:)
  real(dp) :: ll
  integer :: i, ia_err

  ! Allocate arrays
  allocate(rl(n), rmlp1(n), arg1(n), arg2(n), arg3(n), arg4(n), stat=ia_err)
  if (ia_err /= 0) stop 'cei: allocation failed'

  ! Compute initial arrays
  !$omp parallel do private(i)
  do i = 1, n
    rl(i) = r(i)**lambda
    rmlp1(i) = 1.0_dp / (rl(i) * r(i))
    arg1(i) = f2(i) * f4(i) * rl(i)
  end do
  !$omp end parallel do

  ! First integration
  ll = l2 + l4 + real(lambda, dp)
  call dintg(ll, 1, arg1, arg2, r, dx, n, rws)

  ! Update arg2
  ll = ll + 1.0_dp
  !$omp parallel do private(i)
  do i = 1, n
    arg2(i) = arg2(i) * f1(i) * f3(i) * rmlp1(i)
  end do
  !$omp end parallel do

  ! Second integration preparation
  ll = ll + l1 + l3 - real(lambda, dp) - 1.0_dp
  !$omp parallel do private(i)
  do i = 1, n
    arg3(i) = f2(i) * f4(i) * rmlp1(i)
  end do
  !$omp end parallel do

  ! Second integration
  ll = l2 + l4 - real(lambda, dp) - 1.0_dp
  call dintg(ll, -1, arg3, arg4, r, dx, n, rws)

  ! Update arg4
  ll = ll + 1.0_dp
  !$omp parallel do private(i)
  do i = 1, n
    arg4(i) = arg4(i) * f1(i) * f3(i) * rl(i)
  end do
  !$omp end parallel do

  ! Combine results
  ll = ll + l1 + l3 + real(lambda, dp)
  !$omp parallel do private(i)
  do i = 1, n
    arg2(i) = arg2(i) + arg4(i)
  end do
  !$omp end parallel do

  ! Final integration
  result = sintg(ll, arg2, r, dx, n)

  ! Clean up
  deallocate(rl, rmlp1, arg1, arg2, arg3, arg4)
end subroutine cei

function sintg(ll, fct, r, dx, n) result(res)
  use cei_constants
  use omp_lib
  implicit none
  ! Dummy arguments
  real(dp), intent(in) :: ll, fct(n), r(n), dx
  integer, intent(in) :: n
  real(dp) :: res
  ! Local variables
  real(dp) :: sum, rll, fact, x0, x1, x2
  integer :: i

  ! Initial integration term
  if (ll + 1.0_dp <= tiny) then
    write(9, *) 'sintg: ll+1=', ll + 1
    sum = 0.0_dp
  else
    rll = r(1)**ll
    if (rll <= tiny) then
      sum = 0.0_dp
    else
      fact = fct(1) / rll
      sum = fact * (r(1)**(ll + 1)) / (ll + 1)
    end if
  end if

  ! Simpson's rule integration
  x0 = fct(1) * r(1)
  !$omp parallel do private(i, x1, x2) reduction(+:sum)
  do i = 3, n, 2
    x1 = fct(i - 1) * r(i - 1)
    x2 = fct(i) * r(i)
    sum = sum + dx * (x0 + 4.0_dp * x1 + x2)
    x0 = x2
  end do
  !$omp end parallel do
  res = sum / 3.0_dp

  ! Handle even number of points
  if (mod(n, 2) == 0) then
    res = res + (fct(n) + fct(n - 1)) / 2.0_dp * (r(n) - r(n - 1))
  end if
end function sintg

subroutine dintg(ll, idir, fct, yint, r, dx, n, rws)
  use cei_constants
  use omp_lib
  implicit none
  ! Dummy arguments
  real(dp), intent(in) :: ll, fct(n), r(n), dx, rws
  integer, intent(in) :: idir, n
  real(dp), intent(out) :: yint(n)
  ! Local variables
  real(dp) :: sum, rll, fact, x0, x1, corr
  integer :: i, j

  if (idir > 0) then
    ! Forward integration
    if (ll + 1.0_dp <= tiny) then
      write(9, *) 'dintg: ll+1=', ll + 1
      sum = 0.0_dp
    else
      rll = r(1)**ll
      if (rll <= tiny) then
        sum = 0.0_dp
      else
        fact = fct(1) / rll
        sum = fact * (r(1)**(ll + 1)) / (ll + 1)
      end if
    end if
    yint(1) = sum
    x0 = fct(1) * r
    !$omp parallel do private(i, x1) reduction(+:sum)
    do i = 2, n
      x1 = fct(i) * r(i)
      sum = sum + dx * (x0 + x1) / 2.0_dp
      yint(i) = sum
      x0 = x1
    end do
    !$omp end parallel do
  else
    ! Backward integration
    yint(n) = 0.0_dp
    sum = 0.0_dp
    x1 = fct(n) * r(n)
    !$omp parallel do private(i, x0, corr, j) reduction(+:sum)
    do i = n - 1, 1, -1
      x0 = fct(i) * r(i)
      sum = sum + dx * (x0 + x1) / 2.0_dp
      yint(i) = sum
      if (i == n - 3) then
        call inter1(r(n - 3), yint(n - 3), 4, 1, rws, corr)
        sum = sum - corr
        do j = n - 3, n
          yint(j) = yint(j) - corr
        end do
      end if
      x1 = x0
    end do
    !$omp end parallel do
  end if
end subroutine dintg

subroutine inter1(r, p, n, id, rs, ps)
  use cei_constants
  implicit none
  ! Dummy arguments
  real(dp), intent(in) :: r(n), p(n), rs
  integer, intent(in) :: n, id
  real(dp), intent(out) :: ps
  ! Local variables
  real(dp) :: term, denom
  integer :: i, j

  ps = 0.0_dp
  do j = 1, n, id
    term = 1.0_dp
    denom = 1.0_dp
    do i = 1, n, id
      if (i == j) cycle
      denom = denom * (r(j) - r(i))
      term = term * (rs - r(i))
    end do
    ps = ps + term * p(j) / denom
  end do
end subroutine inter1
