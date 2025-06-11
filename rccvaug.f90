module constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: cfac = 13.606_dp
  real(dp), parameter :: pi = 3.141596_dp
  integer, parameter :: ne0 = 500, ne1 = 301, nem0 = 71
  integer, parameter :: lmax0 = 3
  integer, parameter :: ngc = 100
end module constants

module types
  implicit none
  type :: spectrum_data
    real(8), allocatable :: e0(:), e1(:), em0(:)
    real(8), allocatable :: pdos(:,:), pspect(:,:), stot(:), shelp(:)
    real(8), allocatable :: ptot(:), phelp(:)
    real(8), allocatable :: pmat(:,:), z(:)
  end type spectrum_data
end module types

program rccvaug
  use constants
  use types
  use omp_lib
  implicit none

  character(len=80) :: head
  character(len=1) :: name(10), iddos(10), idcore1(5), idcore3(5)
  type(spectrum_data) :: spec
  real(dp) :: hws, hwv, hwc, ef, decore, a0, conv, norm = 0.0_dp
  integer :: krel, lmax, maxd, maxm, mine1, n1, ios
  integer :: dos = 1, mat = 2, tty = 5, prt = 3, plot = 7
  integer :: i, k, l

  ! Allocate arrays
  allocate(spec%e0(ne0), spec%e1(ne1), spec%em0(nem0))
  allocate(spec%pdos(ne0,-lmax0-1:lmax0), spec%pspect(ne0,-lmax0-1:lmax0))
  allocate(spec%stot(ne0), spec%shelp(ne0))
  allocate(spec%ptot(ne1), spec%phelp(ne1))
  allocate(spec%pmat(nem0,-lmax0-1:lmax0), spec%z(-lmax0-1:lmax0))

  ! Print header
  print '(/,a)', ' *****************************************************'
  print '(a)',   '           this is relativistic CCV Auger  '
  print '(a)',   ' *****************************************************'

  ! Open files
  call ccvopen(dos, mat, tty, prt, plot)

  ! Read input
  read(tty, '(a80)') head
  write(plot, '(1x,a80)') head
  read(tty, *) hws
  read(tty, *) hwv
  read(tty, *) hwc
  read(tty, *) krel

  ! Read matrix elements
  read(mat, *)
  read(mat, '(10a1)') name
  read(mat, *)
  read(mat, '(5a1)') idcore1
  read(mat, *)
  read(mat, '(5a1)') idcore3
  read(mat, *)
  read(mat, *) decore
  read(mat, *)
  read(mat, *) lmax
  read(mat, *)
  read(mat, *) a0
  read(mat, *)

  conv = 2.0_dp * pi / a0
  conv = conv * conv
  if (krel == 0) conv = 1.0_dp

  ! Read matrix elements
  maxm = 0
  do
    maxm = maxm + 1
    read(mat, *, iostat=ios) spec%em0(maxm), &
                           (spec%pmat(maxm,k), k=-lmax-1,-1), &
                           (spec%pmat(maxm,k), k=1,lmax)
    if (ios /= 0) exit
  end do
  maxm = maxm - 1

  ! Read DOS
  read(dos, '(10a1)') iddos
  read(dos, *)
  read(dos, *) ef
  ef = conv * ef
  read(dos, *) lmax
  read(dos, *) maxd
  read(dos, *)

  ! Parallelized DOS reading
  !$omp parallel do private(i, l, bla, fl, fl1, fl2)
  do i = 1, maxd
    if (krel == 0) then
      read(dos, *) spec%e0(i), (spec%z(l), l=0,lmax), bla
      spec%e0(i) = spec%e0(i) * conv
      spec%pdos(i,-1) = spec%z(0)
      do l = 1, lmax
        fl = real(2*l+1, dp)
        fl = fl + fl
        fl1 = real(2*l, dp) / fl
        fl2 = real(2*l+2, dp) / fl
        spec%pdos(i,l) = spec%z(l) * fl1
        spec%pd ш
        spec%pdos(i,-l-1) = spec%z(l) * fl2
      end do
    else
      read(dos, *) spec%e0(i), bla, spec%pdos(i,-1),  &
                 (spec%pdos(i,l), spec%pdos(i,-l-1), l=1,lmax)
      spec%e0(i) = spec%e0(i) * conv
    end if
  end do
  !$omp end parallel do

  ! Printout
  write(prt, '(/,a,10a1,/,a,5a1,/,a,5a1,/,a,f10.3,a,/,a,f10.3,a,/,a,f10.3,a,/,a,f10.3,a,/,a,f10.3,a)') &
    ' Relativistic CCV Auger spectrum of ', iddos, &
    '                   final core state ', idcore1, &
    '                 initial core state ', idcore3, &
    '                  transition energy ', decore, ' ryd', &
    '                       fermi-energy ', ef, ' ryd', &
    '              spectrometer (fwhm) : ', hws, '  Ev', &
    '         valence life-time (fwhm) : ', hwv, '  Ev', &
    '      core state life-time (fwhm) : ', hwc, '  Ev'
  write(prt, '(/,a)') ' Matrixelements'
  do i = 1, maxm
    write(prt, '(f11.5,8e13.5)') spec%em0(i), (spec%pmat(i,k), k=-lmax-1,-1), &
                                (spec%pmat(i,k), k=1,lmax)
  end do
  write(prt, '(/,a)') ' DOS'
  do i = 1, maxd
    write(prt, '(10f11.5)') spec%e0(i), spec%pdos(i,-1), &
                          (spec%pdos(i,l), spec%pdos(i,-l-1), l=1,lmax)
  end do

  ! Synspec
  print *, 'starting synspec'
  call synspec(spec%e0, spec%pdos, ne0, maxd, spec%pspect, spec%stot, &
               spec%em0, spec%pmat, nem0, maxm, lmax0, lmax, prt, 1)

  write(prt, '(/,a)') '  total unbroadened spectrum'
  do i = 1, maxd
    write(prt, '(f11.5,e13.5)') (spec%e0(i) - ef) * cfac, spec%stot(i)
  end do

  ! Broadening
  print *, 'starting broad'
  call broad(spec%e0, spec%stot, ne0, maxd, mine1, prt, hwc, hwv, hws, &
             spec%e1, spec%ptot, spec%phelp, ne1, ef, cfac, norm)

  ! Partial spectra broadening
  !$omp parallel do private(i, k)
  do k = -lmax-1, lmax
    if (k == 0) cycle
    spec%shelp(1:ne0) = spec%pspect(1:ne0,k)
    print *, 'starting broad'
    call broad(spec%e0, spec%shelp, ne0, maxd, mine1, prt, hwc, hwv, hws, &
               spec%e1, spec%pspect(1:ne0,k), spec%phelp, ne1, ef, cfac, norm)
  end do
  !$omp end parallel do

  ! Final output
  n1 = ne1 - mine1 + 1
  write(plot, '(2x,5a1,2x,5a1,/,f15.5)') idcore1, idcore3, decore*cfac
  write(plot, '(5x,a,5x,a,4x,a,6(5x,i2,6x))') 'e', 'Total', '-1', (l,-l-1,l=1,lmax)
  do i = mine1, ne1
    write(plot, '(f11.5,8e13.5)') spec%e1(i), spec%ptot(i)/norm, spec%pspect(i,-1)/norm, &
                                (spec%pspect(i,l)/norm, spec%pspect(i,-l-1)/norm, l=1,lmax)
  end do
  write(plot, '(5x,a,5x,a,5x,a,6(5x,i2,2x))') 'e', 'Total', '-1', (l,-l-1,l=1,lmax)
  do i = mine1, ne1
    write(plot, '(f11.5,8f9.3)') spec%e1(i), spec%ptot(i), spec%pspect(i,-1), &
                               (spec%pspect(i,l), spec%pspect(i,-l-1), l=1,lmax)
  end do

  ! Deallocate arrays
  deallocate(spec%e0, spec%e1, spec%em0)
  deallocate(spec%pdos, spec%pspect, spec%stot, spec%shelp)
  deallocate(spec%ptot, spec%phelp, spec%pmat, spec%z)

  stop 'end rccvaug'
end program rccvaug

subroutine ccvopen(dos, mat, tty, prt, plot)
  implicit none
  integer, intent(out) :: dos, mat, tty, prt, plot
  character(len=40) :: name

  dos = 1; mat = 2; prt = 3; tty = 5; plot = 7
  open(unit=tty, file='rccvaug.in')
  read(tty, '(a40)') name
  open(unit=dos, file=trim(name))
  read(tty, '(a40)') name
  open(unit=mat, file=trim(name))
  read(tty, '(a40)') name
  open(unit=prt, file=trim(name))
  read(tty, '(a40)') name
  open(unit=plot, file=trim(name))
end subroutine ccvopen

subroutine synspec(e0, pdos, ne0, maxd, pspect, stot, em0, pmat, nem0, maxm, ldim, lmax, prt, ipr)
  use constants
  use omp_lib
  implicit none
  integer, intent(in) :: ne0, maxd, nem0, maxm, ldim, lmax, prt, ipr
  real(dp), intent(in) :: e0(:), em0(:), pdos(:,:), pmat(:,:)
  real(dp), intent(out) :: pspect(:,:), stot(:)
  real(dp) :: xmat, sum
  integer :: i, k, iex, j, nc, iw, k2

  !$omp parallel do private(i, k, xmat)
  do k = -lmax-1, lmax
    if (k == 0) cycle
    do i = 1, maxd
      xmat = ylag(e0(i), em0, pmat(1:maxm,k), 0, 3, maxm, iex)
      pspect(i,k) = xmat * pdos(i,k)
    end do
  end do
  !$omp end parallel do

  if (ipr > 0) write(prt, '(/,a,/,a)') &
    ' synspec: total unbroadened spectrum', &
    '   e     intensity      e     intensity      e     intensity      e     intensity'

  !$omp parallel do private(i, sum)
  do i = 1, maxd
    sum = sum(pspect(i,-lmax-1:lmax), mask=(/(k,k=-lmax-1,lmax)/ /= 0))
    stot(i) = sum
  end do
  !$omp end parallel do

  if (ipr > 0) then
    nc = 4
    iw = maxd/nc + 1
    do i = 1, iw
      k = nc * iw + i
      if (k > maxd) k = maxd
      write(prt, '(1x,4(0pf6.3,1pe12.4,3x))') (e0(j), stot(j), j=i,k,iw)
    end do
  end if
end subroutine synspec

function ylag(xi, x, y, ind1, n1, imax, iex) result(res)
  use constants
  implicit none
  real(dp) :: res
  real(dp), intent(in) :: xi, x(:), y(:)
  integer, intent(in) :: ind1, n1, imax
  integer, intent(out) :: iex
  integer :: ind, n, j, i, inl, inu
  real(dp) :: s, p, d, xd

  ind = ind1
  n = n1
  iex = 0
  if (n > imax) then
    n = imax
    iex = n
  end if
  if (ind <= 0) then
    do j = 1, imax
      if (xi == x(j)) then
        res = y(j)
        return
      else if (xi < x(j)) then
        ind = j
        exit
      end if
    end do
    if (ind == 0) iex = 1
  end if
  if (ind <= 1) iex = -1
  inl = ind - (n+1)/2
  if (inl <= 0) inl = 1
  inu = inl + n - 1
  if (inu > imax) then
    inl = imax - n + 1
    inu = imax
  end if
  s = 0.0_dp
  p = 1.0_dp
  do j = inl, inu
    p = p * (xi - x(j))
    d = 1.0_dp
    do i = inl, inu
      if (i == j) then
        xd = xi
      else
        xd = x(j)
      end if
      d = d * (xd - x(i))
    end do
    s = s + y(j) / d
  end do
  res = s * p
end function ylag

subroutine broad(e0, a0, n0, max, mine1, prt, hwc, hwv, hws, e1, a1, ai, n1, ef, cfac, snorm)
  use constants
  use omp_lib
  implicit none
  integer, intent(in) :: n0, max, prt, n1
  integer, intent(out) :: mine1
  real(dp), intent(in) :: e0(:), a0(:), hwc, hwv, hws, ef, cfac
  real(dp), intent(inout) :: snorm
  real(dp), intent(out) :: e1(:), a1(:), ai(:)
  real(dp) :: gc(-ngc:ngc), ah(-ngc:ngc)
  real(dp) :: emin0, xmin = -14.0_dp, estep = 0.05_dp
  real(dp) :: rate1, gam, smax, afac, bfac, delta, delta2
  integer :: ix0 = 281, i, iex

  if (snorm <= 0.0_dp) then
    !$omp parallel do
    do i = 1, max
      e0(i) = (e0(i) - ef) * cfac
    end do
    !$omp end parallel do

    do i = 1, n1
      e1(i) = xmin + (i-1) * estep
    end do
    do i = 1, n1
      if (a0(i) > 0.0_dp) then
        emin0 = e0(i)
        exit
      end if
    end do
    do i = 1, n1
      if (e1(i) > emin0) then
        mine1 = i
        exit
      end if
    end do
  end if

  ! Interpolate
  !$omp parallel do private(iex)
  do i = mine1, n1
    ai(i) = ylag(e1(i), e0, a0, 0, 3, max, iex)
  end do
  !$omp end parallel do
  ai(1:mine1-1) = 0.0_dp
  ai(ix0+1:n1) = 0.0_dp

  ! Lifetime broadening
  !$omp parallel do private(rate1, gam)
  do i = mine1, n1
    rate1 = e1(i) / e1(1)
    gam = (hwv * rate1 * rate1 + hwc + hwc) / 2.0_dp
    if (e1(i) > 0.0_dp) gam = hwc
    a1(i) = convlo(gam, ai, i, estep, mine1, n1, pi)
  end do
  !$omp end parallel do

  ai(mine1:n1) = a1(mine1:n1)

  ! Gaussian convolution
  afac = -log(0.5_dp) / (hws / 2.0_dp)**2
  bfac = sqrt(afac / pi)
  !$omp parallel do private(delta, delta2)
  do i = -ngc, ngc
    delta = i * estep
    delta2 = delta * delta
    gc(i) = exp(-afac * delta2) * bfac
  end do
  !$omp end parallel do

  !$omp parallel do
  do i = mine1, n1
    a1(i) = convgau(ai, i, gc, ah, estep, ngc, n1)
  end do
  !$omp end parallel do

  if (snorm <= 0.0_dp) then
    smax = maxval(a1(mine1:n1))
    snorm = 100.0_dp / smax
  end if

  !$omp parallel do
  do i = mine1, n1
    a1(i) = a1(i) * snorm
  end do
  !$omp end parallel do
end subroutine broad

function convlo(gam, a1, ind, estep, mine1, n1, pi) result(res)
  use constants
  implicit none
  real(dp) :: res
  real(dp), intent(in) :: gam, a1(:), estep, pi
  integer, intent(in) :: ind, mine1, n1
  real(dp) :: sum, atp, atm
  integer :: i

  sum = 0.0_dp
  do i = mine1, n1
    atp = atan(estep * ((i-ind) + 0.5_dp) / gam)
    atm = atan(estep * ((i-ind) - 0.5_dp) / gam)
    sum = sum + a1(i) * (atp - atm)
  end do
  res = sum / pi
end function convlo

function convgau(a1, ind, gc, ah, estep, ngc, n1) result(res)
  use constants
  implicit none
  real(dp) :: res
  real(dp), intent(in) :: a1(:), gc(-ngc:), estep
  real(dp), intent(out) :: ah(-ngc:)
  integer, intent(in) :: ind, ngc, n1
  integer :: i, imx

  ah(-ngc:ngc) = 0.0_dp
  do i = 0, ngc
    imx = ind + i
    if (imx > n1) exit
    ah(i) = a1(imx) * gc(i)
  end do
  do i = -1, -ngc, -1
    imx = ind + i
    if (imx < 1) exit
    ah(i) = a1(imx) * gc(i)
  end do
  res = sum(ah(-ngc:ngc)) * estep
end function convgau
