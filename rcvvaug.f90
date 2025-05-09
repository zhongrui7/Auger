! Program to calculate relativistic core-valence-valence Auger spectra
! Parallelized with OpenMP by Grok, May 2025

program rcvvaug
  use omp_lib
  implicit none

  ! Constants
  integer, parameter :: nemax = 500, nemax2 = 1000, n1m = 300, nemf = 15, nemf2 = 30
  integer, parameter :: km = 3, kmtot = 28
  real(8), parameter :: pi = 3.1415926536, efac = 13.606, tiny = 1.0e-5

  ! Variables
  character(len=50) :: name, filen
  character(len=78) :: title
  character(len=10) :: namorb
  real(8) :: e1(n1m), e(nemax), efile(nemf), efile2(nemf2), e2(nemax2)
  integer :: irange(nemf,nemf2)
  real(8) :: z(-km-1:km), udos(nemax,-km-1:km), dos(nemax,-km-1:km)
  real(8) :: crsfile(nemf,nemf2,-km-1:km,-km-1:km), crs(nemax,-km-1:km,-km-1:km)
  real(8) :: p(nemax2,kmtot), psum(kmtot), pdiag(nemax2), ptot(nemax2)
  integer :: iord(kmtot), ikap(kmtot), ikapp(kmtot)
  real(8) :: y(nemax), phelp(n1m), p1(n1m,4), p1tot(n1m), p1diag(n1m)

  ! I/O units
  integer :: in = 1, jdos = 2, jmat = 3, jprt = 7, jspect = 8

  ! Local variables
  real(8) :: hws, hwlev, hwlr, hwc, hwlf, u, ecore, e0, de
  real(8) :: ef, a0, conv, et0, e2min, e2max, de2, pmax, psumtot
  real(8) :: emax, emin, pdiagg, ptott, smax, snorm, e2maxev
  real(8) :: energy, bla, fl, fl1, fl2, eep, ee, dum
  integer :: nume, nume2, lmax, maxd, krel, necount
  integer :: i, j, ii, jj, k, kp, ie2, ne, nef, ktot, ktest
  integer :: ifst, ilast, nval, idoc, iminb, n1

  ! Open I/O files
  open(unit=in, file='rcvvaug.in')
  read(in,'(a50)') filen
  print '(a50)', filen
  open(unit=jdos, file=filen)

  read(in,'(a50)') filen
  print '(a50)', filen
  open(unit=jmat, file=filen)

  read(in,'(a50)') filen
  print '(a50)', filen
  open(unit=jprt, file=filen)

  read(in,'(a50)') filen
  print '(a50)', filen
  open(unit=jspect, file=filen)

  ! Read input parameters
  read(in,'(a78)') title
  print '(a78)', title

  read(in,*) hws
  print *, hws

  read(in,*) hwlev
  print *, hwlev
  hwlr = hwlev/efac

  read(in,*) hwc
  print *, hwc

  read(in,*) hwlf
  print *, hwlf

  read(in,*) krel
  print *, krel

  read(in,*) u
  print *, u

  ! Read matrix elements
  read(jmat,*)
  read(jmat,'(a50)') name
  read(jmat,*)
  read(jmat,'(a10)') namorb
  read(jmat,*)
  read(jmat,*) ecore
  read(jmat,*)
  read(jmat,*) e0, de, nume
  read(jmat,*)
  read(jmat,*) lmax
  read(jmat,*)
  read(jmat,*) a0
  read(jmat,*)
  read(jmat,*)

  write(jprt,*) ' Matrixelements'
  write(jprt,'(a50)') name
  write(jprt,'(a10)') namorb
  write(jprt,'(/1x,"i2",2x,"i",3x,"e2",4x,"ep",4x,"e",2x,28(i6,i3,"    "))') &
    ((kap,kapp,kapp=-lmax-1,kap),kap=-lmax-1,-1), &
    ((kap,kapp,kapp=-lmax-1,-1),(kap,kapp,kapp=1,kap),kap=1,lmax)

  write(jspect,'(a78)') title
  write(jspect,'(a50)') name
  write(jspect,'("   HWS=",f5.2," eV","    HWL=",f5.2," eV", &
                "   HWC=",f5.2," eV","   HWLF=",f5.2," eV")') hws, hwlev, hwc, hwlf

  conv = 2.0*pi/a0
  conv = conv*conv
  if (krel == 0) conv = 1.0

  nume2 = 2*nume - 1
  necount = nume*nume
  efile(1:nume) = [(e0 + de*(i-1), i=1,nume)]

  et0 = 2.0*e0 - ecore
  efile2(1:nume2) = [(et0 + de*(i-1), i=1,nume2)]

  do i = 1, nume2
    do j = 1, nume
      jp = i + 1 - j
      irange(j,i) = merge(1, 0, jp >= 1 .and. jp <= nume)
      if (irange(j,i) == 0) cycle

      read(jmat,*) ii, jj, dum, eep, ee, &
        ((crsfile(jj,ii,kap,kapp), kapp=-lmax-1,kap), kap=-lmax-1,-1), &
        ((crsfile(jj,ii,kap,kapp), kapp=-lmax-1,-1), &
         (crsfile(jj,ii,kap,kapp), kapp=1,kap), kap=1,lmax)

      write(jprt,'(2i3,3f6.2,28d13.5)') i, j, efile2(i), efile(jp), efile(j), &
        ((crsfile(j,i,kap,kapp), kapp=-lmax-1,kap), kap=-lmax-1,-1), &
        ((crsfile(j,i,kap,kapp), kapp=-lmax-1,-1), &
         (crsfile(j,i,kap,kapp), kapp=1,kap), kap=1,lmax)

      if (ii /= i .or. jj /= j .or. abs(eep-efile(jp)) > tiny .or. &
          abs(ee-efile(j)) > tiny) then
        print *, i, j, efile(jp), efile(j)
        print *, ii, jj, eep, ee
        stop 'Error in matrix element data'
      end if
    end do
  end do

  ! Read density of states
  read(jdos,'(a50)') name
  read(jdos,*)
  read(jdos,*) ef
  ef = conv*ef
  read(jdos,*) lmax
  read(jdos,*) maxd
  read(jdos,*)

  ii = 0
  do i = 1, maxd
    if (krel == 0) then
      read(jdos,*) energy, (z(l), l=0,lmax), bla
      if (bla < tiny .and. ii == 0) cycle
      ii = ii + 1
      e(ii) = energy*conv
      udos(ii,-1) = z(0)
      do l = 1, lmax
        fl = dble(2*l+1)
        fl = fl + fl
        fl1 = dble(2*l)/fl
        fl2 = dble(2*l+2)/fl
        udos(ii,l) = z(l)*fl1
        udos(ii,-l-1) = z(l)*fl2
      end do
    else
      read(jdos,*) energy, bla, z(-1), (z(l), z(-l-1), l=1,lmax)
      if (bla < tiny .and. ii == 0) cycle
      ii = ii + 1
      e(ii) = energy*conv
      udos(ii,-lmax-1:lmax) = z(-lmax-1:lmax)
    end if
  end do
  ne = ii

  write(jprt,*)
  write(jprt,*) ' Densities of states'
  do i = 1, ne
    write(jprt,'(f11.5,6f9.3)') e(i), udos(i,-1), &
      (udos(i,l), udos(i,-l-1), l=1,lmax)
  end do

  nef = iright(e, ne, ef)

  ! Check energy ranges
  if (e(1) < efile(1) .or. e(nef) > efile(nume)) then
    print *, 'Energy range of matrix elements must exceed that of DOS!'
    stop
  end if

  ! Lifetime broadening of DOS
  call dosbroad(ef, nef, e, udos, hwlr, dos, lmax)

  write(jprt,*)
  write(jprt,*) ' Densities of states after broadening'
  do i = 1, nef
    write(jprt,'(f11.5,6f9.3)') e(i), dos(i,-1), &
      (dos(i,l), dos(i,-l-1), l=1,lmax)
  end do

  ! Set up index arrays
  ii = 0
  do k = -lmax-1, lmax
    if (k == 0) cycle
    do kp = -lmax-1, k
      if (kp == 0) cycle
      ii = ii + 1
      ikap(ii) = k
      ikapp(ii) = kp
    end do
  end do
  ktot = ii
  print *, ' ktot=', ktot

  ! Determine boundaries for energy loop
  e1ryd = ecore
  e2min = 2.0*e(1) - e1ryd
  e2max = 2.0*e(nef) - e1ryd

  ! Set up final energy scale
  de = e(2) - e(1)
  de2 = e(3) - e(1)
  ne2 = int((e2max - e2min)/de2) + 1
  e2(1:ne2) = [(e2min + de2*(i-1), i=1,ne2)]

  ktest = ne2/6

  ! Main loop over final energies, parallelized with OpenMP
  !$omp parallel do private(ie2, emax, emin, ilast, ifst, nval, idoc, y, pdiagg, ptott, ii, k, kp, i) &
  !$omp shared(ne2, e2, e1ryd, e, nef, de, ktest, crsfile, efile, nume, efile2, nume2, nemf, nemf2, lmax, &
  !$omp irange, dos, p, pdiag, ptot, ikap, ikapp, ktot, jprt)
  do ie2 = 1, ne2
    ! Each thread prints its progress
    !$omp critical
    print *, 'Thread ', omp_get_thread_num(), ': ie2=', ie2, (e2(ie2) + e1ryd)/2.0
    !$omp end critical

    emax = min(e1ryd + e2(ie2) - e(1), e(nef))
    emin = max(e1ryd + e2(ie2) - e(nef), e(1))
    ilast = iright(e, nef, emax)
    ifst = leftind(e, nef, emin)
    nval = ilast - ifst + 1

    idoc = merge(1, 0, mod(ie2,ktest) == 0 .or. ie2 <= 2 .or. ie2 >= ne2-1)

    ! Allocate thread-private crs array to avoid race conditions
    block
      real(8) :: crs(nemax,-km-1:km,-km-1:km)
      call matinter(crsfile, efile, nume, efile2, nume2, nemf, nemf2, lmax, &
                    irange, e2(ie2), emin, emax, e, nemax, crs, ifst, ilast, idoc, jprt)

      ! Compute spectrum
      ii = 0
      do k = -lmax-1, lmax
        if (k == 0) cycle
        do kp = -lmax-1, k
          if (kp == 0) cycle
          ii = ii + 1
          y(ifst:ilast) = [(crs(i,k,kp)*dos(i,k)*dos(ilast+ifst-i,kp), i=ifst,ilast)]
          p(ie2,ii) = simpson(y(ifst:ilast), nval, de)
        end do
      end do
    end block

    pdiagg = sum(p(ie2,1:ktot), mask=ikap(1:ktot) == ikapp(1:ktot))
    ptott = sum(p(ie2,1:ktot))
    pdiag(ie2) = pdiagg
    ptot(ie2) = ptott
  end do
  !$omp end parallel do

  ! Compute relative magnitudes
  psumtot = sum(ptot(1:ne2))
  psum(1:ktot) = [(sum(p(1:ne2,ii))/psumtot, ii=1,ktot)]

  pmax = maxval(ptot(1:ne2))

  call order(psum, ktot, iord)

  ! Broadening, parallelized with OpenMP
  snorm = 0.0
  print *, ' broad'
  call broad(e2, ptot, nemax2, ne2, iminb, e1, p1tot, phelp, n1m, n1, snorm, smax, hws, hwlf, hwc, u)

  ! Parallelize the loop over the four most important contributions
  !$omp parallel do private(i, dummy) shared(e2, p, nemax2, ne2, iminb, e1, p1, phelp, n1m, n1, snorm, hws, hwlf, hwc, u, iord)
  do i = 1, 4
    !$omp critical
    print *, 'Thread ', omp_get_thread_num(), ': broad i=', i
    !$omp end critical
    call broad(e2, p(:,iord(i)), nemax2, ne2, iminb, e1, p1(:,i), phelp, n1m, n1, snorm, dummy, hws, hwlf, hwc, u)
  end do
  !$omp end parallel do

  print *, ' broad'
  call broad(e2, pdiag, nemax2, ne2, iminb, e1, p1diag, phelp, n1m, n1, snorm, dummy, hws, hwlf, hwc, u)

  ! Output spectra
  e2maxev = e2max*efac

  write(jprt,*) ' unbroadened spectra'
  write(jprt,'(4x,"e2",5x,5x,"Total",3x,4(5x,i2,1x,i2,3x),5x,"Diag.")') &
    (ikap(iord(i)), ikapp(iord(i)), i=1,4)
  do i = 1, ne2
    write(jprt,'(f11.5,6e13.5)') e2(i)-e2max, ptot(i), &
      (p(i,iord(ii)), ii=1,4), pdiag(i)
  end do

  write(jspect,'(a10)') namorb
  write(jspect,'(f11.5)') ecore*efac
  write(jspect,'(4x,"e2",5x,5x,"Total",3x,4(5x,i2,1x,i2,3x),5x,"Diag.")') &
    (ikap(iord(i)), ikapp(iord(i)), i=1,4)
  do i = iminb, n1
    write(jspect,'(f11.5,6e13.5)') e1(i)-e2maxev, p1tot(i)/snorm, &
      (p1(i,ii)/snorm, ii=1,4), p1diag(i)/snorm
  end do

  write(jspect,'(4x,"e2",5x,4x,"Total",4(4x,i2,1x,i2),4x,"Diag.")') &
    (ikap(iord(i)), ikapp(iord(i)), i=1,4)
  do i = iminb, n1
    write(jspect,'(f11.5,6e13.5)') e1(i)-e2maxev, p1tot(i), &
      (p1(i,ii), ii=1,4), p1diag(i)
  end do

contains

  subroutine dosbroad(ef, nef, e, udos, hwlr, dos, lmax)
    implicit none
    real(8), intent(in) :: ef, hwlr
    integer, intent(in) :: nef, lmax
    real(8), intent(in) :: e(nemax), udos(nemax,-km-1:km)
    real(8), intent(out) :: dos(nemax,-km-1:km)
    real(8) :: a(nemax), gam, estep
    integer :: i, j, k
    real(8), parameter :: pi = 3.141596, hwlr0 = 0.1
    integer, parameter :: ilor = 2

    if (hwlr < hwlr0) then
      do k = -lmax-1, lmax
        if (k == 0) cycle
        dos(1:nef,k) = udos(1:nef,k)
      end do
      return
    end if

    estep = e(2) - e(1)
    !$omp parallel do private(i, j, a, gam) shared(dos, udos, e, ef, nef, lmax, estep, pi, ilor, hwlr)
    do k = -lmax-1, lmax
      if (k == 0) cycle
      do i = 1, nef
        a(1:nef) = udos(1:nef,k)
        if (ilor == 1) then
          gam = hwlr*(ef-e(i))/(ef-e(1))/2.0
        else
          gam = hwlr*((ef-e(i))/(ef-e(1)))**2/2.0
        end if
        dos(i,k) = convlo(gam, a, i, estep, 1, nef, pi)
      end do
    end do
    !$omp end parallel do
  end subroutine dosbroad

  subroutine broad(ein, a0, n0, max0, mine1, e1, a1, ai, n1m, n1, snorm, smax, hws, hwl, hwc, u)
    implicit none
    integer, intent(in) :: n0, max0, n1m
    real(8), intent(in) :: ein(n0), a0(n0), hws, hwl, hwc, u
    integer, intent(out) :: mine1, n1
    real(8), intent(out) :: e1(n1m), a1(n1m), ai(n1m), snorm, smax
    real(8) :: e0(nmax), gc(-ngc:ngc), ah(-ngc:ngc)
    real(8) :: de0, e0max, e00, evmin, evmax, estep, emin0, afac, bfac, dele, gam
    integer :: i, j, max, nlst, i0
    integer, parameter :: ngc = 100, nmax = 500
    real(8), parameter :: hwc0 = 5.e-3, cfac = 13.606, estep_def = 0.05, pi = 3.141596
    integer, parameter :: ix0 = 201

    e0(1:max0) = ein(1:max0)*cfac - u
    de0 = e0(2) - e0(1)
    e0max = e0(max0)
    e00 = e0max
    max = max0
    do i = max0+1, nmax
      e00 = e00 + de0
      e0(i) = e00
      a0(i) = 0.0
      if (e00 - e0max > u) exit
      max = i
    end do

    if (snorm <= 0.0) then
      evmin = e0(1) - 1.0
      evmax = e0(max) + 1.0
      estep = estep_def
      do
        n1 = int((evmax-evmin)/estep) + 1
        if (n1 <= n1m) exit
        estep = estep + estep
      end do
      i0 = int(evmin/estep)
      evmin = i0*estep

      nlst = n1
      e1(1:n1) = [(evmin + (i-1)*estep, i=1,n1)]
      do i = 1, n1
        if (e1(i) > e0(max) .and. nlst == n1) nlst = i
      end do

      do i = 1, n0
        if (a0(i) > 0.0) exit
      end do
      emin0 = e0(i)
      do i = 1, n1
        if (e1(i) > emin0) exit
      end do
      mine1 = i
    end if

    ! Interpolate
    !$omp parallel do private(j)
    do i = mine1, nlst
      j = 1
      do while (e0(j) < e1(i) .and. j < max)
        j = j + 1
      end do
      j = max(1, j-2)
      if (j > max-3) j = max-3
      call inter(e0(j:j+3), a0(j:j+3), 4, e1(i), ai(i))
    end do
    !$omp end parallel do
    ai(1:mine1-1) = 0.0
    ai(nlst+1:n1) = 0.0

    ! Convolute
    !$omp parallel do private(dele, gam)
    do i = mine1, n1
      dele = (e1(nlst) - e1(i))/(e1(nlst) - e1(mine1))
      gam = (hwl*dele*dele + hwc)/2.0
      a1(i) = convlo(gam, ai, i, estep, mine1, n1, pi)
    end do
    !$omp end parallel do
    ai(1:n1) = a1(1:n1)

    if (hws /= 0.0) then
      afac = -log(0.5)/(hws/2.0)**2
      bfac = sqrt(afac/pi)
      gc(-ngc:ngc) = [(exp(-afac*(i*estep)**2)*bfac, i=-ngc,ngc)]
      !$omp parallel do
      do i = mine1, n1
        a1(i) = convgau(ai, i, gc, ah, estep, ngc, n1)
      end do
      !$omp end parallel do
    end if

    if (snorm <= 0.0) then
      smax = maxval(a1(mine1:n1))
      snorm = 100.0/smax
    end if

    a1(mine1:n1) = a1(mine1:n1)*snorm
  end subroutine broad

  real(8) function convlo(gam, a1, ind, estep, mine1, n1, pi)
    implicit none
    real(8), intent(in) :: gam, a1(n1), estep, pi
    integer, intent(in) :: ind, mine1, n1
    real(8) :: sum, atp, atm
    integer :: i
    sum = 0.0
    do i = mine1, n1
      atp = atan(estep*((i-ind)+0.5)/gam)
      atm = atan(estep*((i-ind)-0.5)/gam)
      sum = sum + a1(i)*(atp-atm)
    end do
    convlo = sum/pi
  end function convlo

  real(8) function convgau(a1, ind, gc, ah, estep, ngc, n1)
    implicit none
    integer, intent(in) :: ind, ngc, n1
    real(8), intent(in) :: a1(n1), gc(-ngc:ngc), estep
    real(8), intent(out) :: ah(-ngc:ngc)
    real(8) :: sum
    integer :: i, imx
    ah = 0.0
    do i = 0, ngc
      imx = ind + i
      if (imx > n1) exit
      ah(i) = a1(imx)*gc(i)
    end do
    do Burnham,  i = -1, -ngc, -1
      imx = ind + i
      if (imx < 1) exit
      ah(i) = a1(imx)*gc(i)
    end do
    convgau = sum(ah(-ngc:ngc))*estep
  end function convgau

  subroutine matinter(crsfile, efile, nume, efile2, nume2, nemf, nemf2, km, irange, e2, &
                      emin, emax, e, nemax, crs, ifst, ilast, idoc, jprt)
    implicit none
    integer, intent(in) :: nume, nume2, nemf, nemf2, km, nemax, ifst, ilast, idoc, jprt
    real(8), intent(in) :: crsfile(nemf,nemf2,-km-1:km,-km-1:km), efile(nemf), efile2(nemf2)
    integer, intent(in) :: irange(nemf,nemf2)
    real(8), intent(in) :: e2, emin, emax, e(nemax)
    real(8), intent(out) :: crs(nemax,-km-1:km,-km-1:km)
    real(8) :: ehelp(nm), chelp(nm), crs1(nm,-km-1:km,-km-1:km), ebar
    integer :: if2l, if2r, ifl, ifr, insleft, insright, nint, iextra, nf
    integer :: i, if, l, lp
    integer, parameter :: kml = 3, nm = 50

    if2l = leftind(efile2, nume2, e2)
    if2r = iright(efile2, nume2, e2)
    ifl = leftind(efile, nume, emin)
    ifr = iright(efile, nume, emax)

    insleft = merge(1, 0, irange(ifl,if2r) == 0)
    insright = merge(1, 0, irange(ifr,if2l) == 0)

    nint = min(5, nume)
    if (insleft == 1) then
      ebar = (efile2(nume2) + efile2(1))/2.0
      ehelp(ifl) = e2 - ebar + efile(1)
      do l = -km-1, km
        if (l == 0) cycle
        do lp = -km-1, l
          if (lp == 0) cycle
          chelp(1:nume) = [(crsfile(i,nume2-nume+i,l,lp), i=nume,1,-1)]
          call lagrange(efile(1:nume), chelp(1:nume), nume, nint, ehelp(ifl), &
                        crs1(ifl,l,lp), iextra)
        end do
      end do
      if (idoc == 1) write(jprt,'(" non-grid value for lower e boundary")')
    end if

    if (insright == 1) then
      ehelp(ifr) = e2 - efile2(1) + efile(1)
      do l = -km-1, km
        if (l == 0) cycle
        do lp = -km-1, l
          if (lp == 0) cycle
          chelp(1:nume) = crsfile(1:nume,1,l,lp)
          call lagrange(efile(1:nume), chelp(1:nume), nume, nint, ehelp(ifr), &
                        crs1(ifr,l,lp), iextra)
        end do
      end do
      if (idoc == 1) write(jprt,'(" non-grid value for upper e boundary")')
    end if

    if (idoc == 1) write(jprt,'(" ifl,ifr: ",2i4)') ifl, ifr

    do if = ifl, ifr
      if (if == ifl .and. insleft == 1 .or. if == ifr .and. insright == 1) cycle
      ehelp(if) = efile(if)
      do l = -km-1, km
        if (l == 0) cycle
        do lp = -km-1, l
          if (lp == 0) cycle
          chelp(if:if+nume-1) = crsfile(if,if:if+nume-1,l,lp)
          call lagrange(efile2(if:if+nume-1), chelp(if:if+nume-1), nume, nint, &
                        e2, crs1(if,l,lp), iextra)
        end do
      end do
    end do

    if (idoc == 1) then
      write(jprt,'(" crs1:")')
      do i = ifl, ifr
        write(jprt,'(i3,f10.5,28e13.5)') i, ehelp(i), &
          ((crs1(i,l,lp), lp=-km-1,l), l=-km-1,-1), &
          ((crs1(i,l,lp), lp=-km-1,-1), (crs1(i,l,lp), lp=1,l), l=1,km)
      end do
    end if

    nf = ifr - ifl + 1
    nint = min(6, nf)
    !$omp parallel do private(l, lp, ie, chelp, iextra) shared(crs, e, ifst, ilast, ifl, ifr, crs1, efile, nf, nint, idoc, jprt)
    do ie = ifst, ilast
      do l = -km-1, km
        if (l == 0) cycle
        do lp = -km-1, l
          if (lp == 0) cycle
          chelp(ifl:ifr) = crs1(ifl:ifr,l,lp)
          call lagrange(efile(ifl:ifr), chelp(ifl:ifr), nf, nint, e(ie), &
                        crs(ie,l,lp), iextra)
          if (iextra == 1 .and. idoc == 1) then
            !$omp critical
            write(jprt,'(" extrapolation: crs, ie=",i3," e(ie)=",f10.5/28e13.5)') &
              ie, e(ie), ((crs(ie,k,kp), kp=-km-1,k), k=-km-1,-1), &
              ((crs(ie,k,kp), kp=-km-1,-1), (crs(ie,k,kp), kp=1,k), k=1,km)
            !$omp end critical
          end if
        end do
      end do
    end do
    !$omp end parallel do

    if (idoc == 1) then
      !$omp critical
      write(jprt,'(" e2: ",f10.5)') e2
      do ie = ifst, ilast
        write(jprt,'(i3,f10.5,28e13.4)') ie, e(ie), &
          ((crs(ie,l,lp), lp=-km-1,l), l=-km-1,-1), &
          ((crs(ie,l,lp), lp=-km-1,-1), (crs(ie,l,lp), lp=1,l), l=1,km)
      end do
      write(jprt,'(120("-"))')
      !$omp end critical
    end if
  end subroutine matinter

  subroutine lagrange(x, y, ndim, npkt, xval, yinter, iextra)
    implicit none
    integer, intent(in) :: ndim, npkt
    real(8), intent(in) :: x(ndim), y(ndim), xval
    real(8), intent(out) :: yinter
    integer, intent(out) :: iextra
    integer :: igr, ileft, ifst
    iextra = 0
    if (xval < x(1)) iextra = -1
    if (xval > x(ndim)) iextra = 1
    igr = 1
    do while (x(igr) < xval .and. igr < ndim)
      if (x(igr) == xval) then
        yinter = y(igr)
        return
      end if
      igr = igr + 1
    end do
    ileft = npkt/2
    if (mod(npkt,2) /= 0) ileft = ileft + 1
    ifst = igr - ileft
    ifst = max(1, min(ifst, ndim-npkt+1))
    call inter(x(ifst:ifst+npkt-1), y(ifst:ifst+npkt-1), npkt, xval, yinter)
  end subroutine lagrange

  subroutine inter(r, p, n, rs, ps)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: r(n), p(n), rs
    real(8), intent(out) :: ps
    real(8) :: term, denom
    integer :: i, j
    ps = 0.0
    do j = 1, n
      term = 1.0
      denom = 1.0
      do i = 1, n
        if (i == j) cycle
        denom = denom*(r(j) - r(i))
        term = term*(rs - r(i))
      end do
      ps = ps + term*p(j)/denom
    end do
  end subroutine inter

  integer function leftind(field, ndim, value)
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: field(ndim), value
    integer :: i
    do i = ndim, 1, -1
      if (field(i) < value) exit
    end do
    leftind = max(1, i)
  end function leftind

  integer function iright(field, ndim, value)
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: field(ndim), value
    integer :: i
    do i = 1, ndim
      if (field(i) > value) exit
    end do
    iright = min(ndim, i)
  end function iright

  real(8) function simpson(fct, n, dx)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: fct(n), dx
    real(8) :: sum, x0, x1, x2
    integer :: i
    sum = 0.0
    x0 = fct(1)
    do i = 3, n, 2
      x1 = fct(i-1)
      x2 = fct(i)
      sum = sum + dx*(x0 + 4*x1 + x2)
      x0 = x2
    end do
    simpson = sum/3.0
    if (mod(n,2) == 0) simpson = simpson + (fct(n) + fct(n-1))*dx/2.0
  end function simpson

  subroutine order(f, n, ior)
    implicit none
    integer, intent(in) :: n
    real(8), intent(inout) :: f(n)
    integer, intent(out) :: ior(n)
    real(8) :: g
    integer :: i, j, l
    ior = [(i, i=1,n)]
    do i = 1, n-1
      do j = i+1, n
        if (f(j) > f(i)) then
          g = f(j)
          f(j) = f(i)
          f(i) = g
          l = ior(j)
          ior(j) = ior(i)
          ior(i) = l
        end if
      end do
    end do
  end subroutine order

end program rcvvaug
