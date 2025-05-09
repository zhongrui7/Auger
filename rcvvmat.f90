module parameters
  implicit none
  integer, parameter :: nrad = 250
  integer, parameter :: nemax = 15
  integer, parameter :: nemax2 = 2 * nemax
  integer, parameter :: lmax = 10  ! Assumed; adjust based on l.par
  integer, parameter :: lvmax = 5   ! Assumed; adjust based on l.par
  integer, parameter :: lmax2 = 2 * lmax
  integer, parameter :: ndim = 100
  real(8), parameter :: tiny = 1.0d-10
  real(8), parameter :: small = 0.01d0
  real(8), parameter :: half = 0.5d0
  real(8), parameter :: pi = 3.1415926535897932d0
  real(8), parameter :: ev = 13.606d0
  real(8), parameter :: csmall = 274.0746d0
  real(8), parameter :: chuge = 1.0d8
end module parameters

program rcvvmat
  use parameters
  implicit none

#ifdef _OPENMP
  include 'omp_lib.h'
#endif

  ! Variable declarations
  real(8) :: v(nrad), rs(nrad), xs(nrad)
  real(8) :: earr(nemax), earr2(nemax2)
  integer :: irange(nemax, nemax2)
  real(8) :: gc1(nrad), fc1(nrad)
  real(8) :: gvwf(nrad, -lmax-1:lmax), fvwf(nrad, -lmax-1:lmax)
  real(8) :: gvwfp(nrad, -lmax-1:lmax), fvwfp(nrad, -lmax-1:lmax)
  real(8) :: gcwf(nrad, -lmax-1:lmax), fcwf(nrad, -lmax-1:lmax)
  real(8) :: sd(-lvmax-1:lvmax, -lvmax-1:lvmax)
  real(8) :: sde(-lvmax-1:lvmax, -lvmax-1:lvmax)
  real(8) :: ssum(-lvmax-1:lvmax, -lvmax-1:lvmax, nemax)
  real(8) :: crsecs(-lvmax-1:lvmax, -lvmax-1:lvmax, nemax, nemax2)
  real(8) :: idgg(-lmax-1:lmax, 0:lmax), idgf(-lmax-1:lmax, 0:lmax)
  real(8) :: idfg(-lmax-1:lmax, 0:lmax), idff(-lmax-1:lmax, 0:lmax)
  real(8) :: iegg(-lmax-1:lmax, 0:lmax), iegf(-lmax-1:lmax, 0:lmax)
  real(8) :: iefg(-lmax-1:lmax, 0:lmax), ieff(-lmax-1:lmax, 0:lmax)
  real(8) :: cdgg(-lmax-1:lmax, 0:lmax), cdgf(-lmax-1:lmax, 0:lmax)
  real(8) :: cdfg(-lmax-1:lmax, 0:lmax), cdff(-lmax-1:lmax, 0:lmax)
  real(8) :: cegg(-lmax-1:lmax, 0:lmax), cegf(-lmax-1:lmax, 0:lmax)
  real(8) :: cefg(-lmax-1:lmax, 0:lmax), ceff(-lmax-1:lmax, 0:lmax)
  integer :: rl(-lmax-1:lmax), rlb(-lmax-1:lmax), rj(-lmax-1:lmax)
  real(8) :: triad(-lmax-1:lmax, -lmax-1:lmax, 0:lmax2)
  real(8) :: w6j1(0:lvmax, 0:lmax, 0:lmax, 0:lmax)
  real(8) :: w6j3(0:lvmax, 0:lmax, 0:lmax, 0:lmax)
  real(8) :: w6j5(0:lvmax, 0:lmax, 0:lmax, 0:lmax)
  real(8) :: w6j7(0:lvmax, 0:lmax, 0:lmax, 0:lmax)
  integer :: in_unit, wf_unit, po_unit, mat_unit, pun_unit, zatom
  character(1) :: name(40), norb1(5), norb(5)
  character(50) :: filen
  real(8) :: e0, de, ec1, x0, dx, rmt, rws, a0, zeromt
  integer :: lval, ne, nonrel, norm, nmt, nws, kap1, nwf
  real(8) :: clight, coef, signv, x, e, ep, e2, emin2
  integer :: i, j, kap, kapp, lam, kap2, nei, nei2, neip, kout, kout1
  integer :: l, lp, lkap, lkapp, jj, jjp, l1, j1, jj1
  real(8) :: y, yp, y1, y2, xj1, xj3, t1v, tm1mv, t1vp, tm1mvp
  real(8) :: tv2, tmvm2, tvp2, tmvpm2, result, sigd, sigde, contr
  real(8) :: sum1, sum2, llam, llamp, dummy
  integer :: ii, jp, l2, j3

  ! File unit numbers
  in_unit = 1
  wf_unit = 2
  po_unit = 3
  mat_unit = 7
  pun_unit = 8

  ! Open input and error files
  open(unit=in_unit, file='rcvvmat.in', status='old', action='read')
  open(unit=9, file='rcvvmat.err', status='unknown', action='write')

  ! Main input loop
  do
    ! Read file names
    read(in_unit, '(a50)', iostat=i) filen
    if (i /= 0) exit  ! End of file
    open(unit=pun_unit, file=trim(filen), status='unknown', action='write')

    read(in_unit, '(a50)') filen
    open(unit=po_unit, file=trim(filen), status='old', action='read')

    read(in_unit, '(a50)') filen
    open(unit=wf_unit, file=trim(filen), status='old', action='read')

    read(in_unit, '(a50)') filen
    open(unit=mat_unit, file=trim(filen), status='unknown', action='write')

    ! Initialize relativistic quantum numbers
    do kap = -lmax-1, -1
      l = -kap - 1
      rl(kap) = l
      rlb(kap) = l + 1
      rj(kap) = 2 * l + 1
    end do
    do kap = 1, lmax
      l = kap
      rl(kap) = l
      rlb(kap) = l - 1
      rj(kap) = 2 * l - 1
    end do

    ! Read maximum angular momentum for valence states
    read(in_unit, *) lval
    write(pun_unit, '(1x,a,t20,i4,/)') 'LVAL=', lval

    ! Read energy panel parameters
    read(in_unit, *) e0, de, ne

    ! Check for non-relativistic limit
    read(in_unit, *) nonrel
    if (nonrel == 1) then
      clight = chuge
    else
      clight = csmall
    end if

    ! Read normalization flag
    read(in_unit, *) norm
    write(pun_unit, '(1x,a,i4,/)') &
      'NORMALIZATION OF VALENCE WF. (0 FOR MT, 1 FOR WS):', norm

    ! Read core state identifier
    read(in_unit, '(5a1)') norb1

    read(in_unit, *)  ! Skip blank line

    ! Read potential and mesh parameters
    read(po_unit, '(40a1)') name
    read(po_unit, *) zatom, x0, nmt, dx, rws, a0, zeromt
    read(po_unit, *) v(1:nmt)

    coef = 2.0d0 * zatom / clight
    coef = coef * coef

    ! Set up logarithmic and radial mesh
    x = x0
    xs(1) = x0
    rs(1) = exp(x)
    signv = 1.0d0
    if (v(1) > 0.0d0) signv = -1.0d0
    v(1) = signv * v(1) / rs(1) - zeromt

    do j = 2, nmt
      x = x + dx
      xs(j) = x
      rs(j) = exp(x)
      v(j) = signv * v(j) / rs(j) - zeromt
    end do
    rmt = rs(nmt)

    nws = 0
    do j = nmt + 1, nrad
      x = x + dx
      xs(j) = x
      rs(j) = exp(x)
      if (nws == 0 .and. rws < rs(j)) nws = j
      v(j) = 0.0d0
    end do

    write(*, '(40a1)') name
    write(pun_unit, '(40a1)') name
    write(pun_unit, '(/)')
    write(pun_unit, '(1x,a,t20,i4/,1x,a,t20,f5.1/,1x,a,t20,f10.7/, &
      &1x,a,t20,i4/,1x,a,t20,f10.7/,1x,a,t20,i4/,1x,a,t20,f10.7/, &
      &1x,a,t20,f5.1//1x,a)') &
      ' ZATOM =', zatom, &
      '    X0 =', x0, &
      '    DX =', dx, &
      '   NMT =', nmt, &
      '   RMT =', rmt, &
      '   NWS =', nws, &
      '   RWS =', rws, &
      'ZEROMT =', zeromt, &
      'RADIAL MESH AND POTENTIAL*R'
    write(pun_unit, '(4d20.10)') (rs(j), v(j) * rs(j), j=1, nws)

    ! Read core wavefunctions
    do
      read(wf_unit, '(5a1)') norb
      read(wf_unit, *) ec1
      ec1 = 2.0d0 * ec1
      read(wf_unit, *) kap1
      read(wf_unit, *) nwf
      if (nonrel == 1) then
        do i = 1, nwf
          read(wf_unit, *) dummy, gc1(i)
          fc1(i) = 0.0d0
        end do
      else
        do i = 1, nwf
          read(wf_unit, *) dummy, gc1(i), fc1(i)
        end do
      end if
      if (all(norb(1:3) == norb1(1:3))) exit
    end do

    write(pun_unit, '(//1x,a/1x,a,t20,5a1/1x,a,t20,f12.6/1x,a,t20,i4//)') &
      'CORE WAVEFUNCTION', &
      'ORB. =', norb, &
      ' EC1 =', ec1, &
      'KAP1 =', kap1
    write(pun_unit, '(3d20.10)') (rs(j), gc1(j), fc1(j), j=1, nws)

    write(pun_unit, '(//a//)') ' ****************** END INPUT *******************'

    ! Set up energy panels
    ii = 0
    do i = 1, ne
      earr(i) = e0 + (i - 1) * de
    end do
    ne2 = 2 * ne - 1
    emin2 = 2.0d0 * e0 - ec1
    write(pun_unit, '(1x,a//(7x,15f7.2)/)') 'ENERGY PLANE', earr(1:ne)
    do i = 1, ne2
      earr2(i) = emin2 + (i - 1) * de
      do j = 1, ne
        jp = i + 1 - j
        if (jp >= 1 .and. jp <= ne) then
          ii = ii + 1
          irange(j, i) = ii
        else
          irange(j, i) = 0
        end if
      end do
      write(pun_unit, '(f7.2,15(2x,i3,2x))') earr2(i), irange(:, i)
    end do
    write(pun_unit, '(/)')

    ! Initialize angular integration coefficients
    l1 = rl(kap1)
    j1 = rj(kap1)
    jj1 = j1 + 1
    xj1 = j1 * half
    y1 = sqrt(kap1 * kap1 - coef)
    if (nonrel == 1) y1 = real(l1 + 1, 8)

    call triadini(triad, rl, rlb, rj, pun_unit)

    do l = 0, 3
      j3 = 2 * l + 1
      xj3 = 0.5d0 * real(j3, 8)
      if (l == 0) call wigini(w6j1, xj1, xj3, pun_unit)
      if (l == 1) call wigini(w6j3, xj1, xj3, pun_unit)
      if (l == 2) call wigini(w6j5, xj1, xj3, pun_unit)
      if (l == 3) call wigini(w6j7, xj1, xj3, pun_unit)
    end do

    ! Initialize shared arrays
    ssum = 0.0d0
    crsecs = 0.0d0

    ! Continuum energy loop (parallelized)
    !$omp parallel do schedule(dynamic) &
    !$omp private(nei, neip, e, ep, kout, kout1, &
    !$omp         kap, kapp, lam, kap2, y, yp, jj, jjp, lkap, lkapp, &
    !$omp         idgg, idgf, idfg, idff, iegg, iegf, iefg, ieff, &
    !$omp         cdgg, cdgf, cdfg, cdff, cegg, cegf, cefg, ceff, &
    !$omp         t1v, tm1mv, t1vp, tm1mvp, tv2, tmvm2, tvp2, tmvpm2, &
    !$omp         sigd, sigde, sum1, sum2, llam, llamp, contr, result, &
    !$omp         sd, sde, l2, jj2, i, ip) &
    !$omp shared(ne2, earr2, irange, earr, nonrel, ne, lval, &
    !$omp        rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, norm, &
    !$omp        gcwf, fcwf, gvwf, fvwf, gvwfp, fvwfp, &
    !$omp        gc1, fc1, kap1, coef, triad, w6j1, w6j3, w6j5, w6j7, &
    !$omp        ssum, crsecs, y1, jj1, pun_unit)
    do nei2 = 1, ne2
      e2 = earr2(nei2)
      kout = 0
      if (nei2 == ne) kout = 1
      !$omp critical
      write(pun_unit, '(a,f15.5)') ' e2=', e2
      if (kout == 1) write(pun_unit, '(/a)') ' continuum wavefunctions'
      !$omp end critical

      call wafu(e2, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, &
                gcwf, fcwf, norm, 2, nonrel, pun_unit, lval, coef, kout)

      ! Valence energy loop
      do nei = 1, ne
        neip = nei2 + 1 - nei
        if (irange(nei, nei2) == 0) cycle
        kout1 = 0
        if (nonrel == 1 .and. irange(nei, nei2) == ne * ne) kout1 = 1

        e = earr(nei)
        ep = earr(neip)
        !$omp critical
        if (kout == 1) write(pun_unit, '(/a)') ' valence wavefunctions'
        write(pun_unit, '(t20,a,f15.5)') ' e =', e
        !$omp end critical

        call wafu(e, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, &
                  gvwf, fvwf, norm, 1, nonrel, pun_unit, lval, coef, kout)
        call wafu(ep, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, &
                  gvwfp, fvwfp, norm, 1, nonrel, pun_unit, lval, coef, 0)

        ! Kappa loop (parallelized over kap and kapp)
        do kap = -lval-1, lval
          if (kap == 0) cycle
          y = sqrt(kap * kap - coef)
          if (nonrel == 1) y = real(rl(kap) + 1, 8)
          jj = rj(kap) + 1
          lkap = rl(kap)
          !$omp critical
          if (kout1 == 1) write(pun_unit, '(a,i4)') ' kappa=', kap
          !$omp end critical

          do kapp = -lval-1, lval
            if (kapp == 0) cycle
            yp = sqrt(kapp * kapp - coef)
            if (nonrel == 1) yp = real(rl(kapp) + 1, 8)
            jjp = rj(kapp) + 1
            lkapp = rl(kapp)
            !$omp critical
            if (kout1 == 1) write(pun_unit, '(a,i4)') ' kappa"=', kapp
            !$omp end critical

            ! Initialize integrals
            idgg(:, :) = 0.0d0
            idfg(:, :) = 0.0d0
            idgf(:, :) = 0.0d0
            idff(:, :) = 0.0d0
            iegg(:, :) = 0.0d0
            iefg(:, :) = 0.0d0
            iegf(:, :) = 0.0d0
            ieff(:, :) = 0.0d0
            cdgg(:, :) = 0.0d0
            cdfg(:, :) = 0.0d0
            cdgf(:, :) = 0.0d0
            cdff(:, :) = 0.0d0
            cegg(:, :) = 0.0d0
            cefg(:, :) = 0.0d0
            cegf(:, :) = 0.0d0
            ceff(:, :) = 0.0d0

            ! Compute Coulomb and exchange integrals
            do lam = 0, lmax
              t1v = triad(kap1, kap, lam)
              tm1mv = triad(-kap1, -kap, lam)
              t1vp = triad(kap1, kapp, lam)
              tm1mvp = triad(-kap1, -kapp, lam)

              do kap2 = -lmax-1, lmax
                if (kap2 == 0) cycle
                y2 = sqrt(kap2 * kap2 - coef)
                if (nonrel == 1) y2 = real(rl(kap2) + 1, 8)
                tv2 = triad(kap, kap2, lam)
                tmvm2 = triad(-kap, -kap2, lam)
                tvp2 = triad(kapp, kap2, lam)
                tmvpm2 = triad(-kapp, -kap2, lam)

                if (abs(t1v * tvp2) > tiny) then
                  call cei(gc1, gvwf(:, kap), gvwfp(:, kapp), gcwf(:, kap2), &
                           rs, nws, lam, y1, y, yp, y2, dx, rws, result)
                  idgg(kap2, lam) = result
                  cdgg(kap2, lam) = t1v * tvp2
                end if
                if (abs(t1v * tmvpm2) > tiny) then
                  call cei(gc1, gvwf(:, kap), fvwfp(:, kapp), fcwf(:, kap2), &
                           rs, nws, lam, y1, y, yp, y2, dx, rws, result)
                  idgf(kap2, lam) = result
                  cdgf(kap2, lam) = t1v * tmvpm2
                end if
                if (abs(tm1mv * tvp2) > tiny) then
                  call cei(fc1, fvwf(:, kap), gvwfp(:, kapp), gcwf(:, kap2), &
                           rs, nws, lam, y1, y, yp, y2, dx, rws, result)
                  idfg(kap2, lam) = result
                  cdfg(kap2, lam) = tm1mv * tvp2
                end if
                if (abs(tm1mv * tmvpm2) > tiny) then
                  call cei(fc1, fvwf(:, kap), fvwfp(:, kapp), fcwf(:, kap2), &
                           rs, nws, lam, y1, y, yp, y2, dx, rws, result)
                  idff(kap2, lam) = result
                  cdff(kap2, lam) = tm1mv * tmvpm2
                end if
                if (abs(t1vp * tv2) > tiny) then
                  call cei(gc1, gvwfp(:, kapp), gvwf(:, kap), gcwf(:, kap2), &
                           rs, nws, lam, y1, yp, y, y2, dx, rws, result)
                  iegg(kap2, lam) = result
                  cegg(kap2, lam) = t1vp * tv2
                end if
                if (abs(t1vp * tmvm2) > tiny) then
                  call cei(gc1, gvwfp(:, kapp), fvwf(:, kap), fcwf(:, kap2), &
                           rs, nws, lam, y1, yp, y, y2, dx, rws, result)
                  iegf(kap2, lam) = result
                  cegf(kap2, lam) = t1vp * tmvm2
                end if
                if (abs(tm1mvp * tv2) > tiny) then
                  call cei(fc1, fvwfp(:, kapp), gvwf(:, kap), gcwf(:, kap2), &
                           rs, nws, lam, y1, yp, y, y2, dx, rws, result)
                  iefg(kap2, lam) = result
                  cefg(kap2, lam) = tm1mvp * tv2
                end if
                if (abs(tm1mvp * tmvm2) > tiny) then
                  call cei(fc1, fvwfp(:, kapp), fvwf(:, kap), fcwf(:, kap2), &
                           rs, nws, lam, y1, yp, y, y2, dx, rws, result)
                  ieff(kap2, lam) = result
                  ceff(kap2, lam) = tm1mvp * tmvm2
                end if
              end do
            end do

            ! Sum for direct term
            sigd = 0.0d0
            do lam = 0, lmax
              llam = 2 * lam + 1
              do kap2 = -lmax-1, lmax
                if (kap2 == 0) cycle
                jj2 = rj(kap2) + 1
                sum1 = idgg(kap2, lam) * cdgg(kap2, lam) + &
                       idgf(kap2, lam) * cdgf(kap2, lam) + &
                       idfg(kap2, lam) * cdfg(kap2, lam) + &
                       idff(kap2, lam) * cdff(kap2, lam)
                sigd = sigd + jj1 * jj2 * sum1 * sum1 / llam
              end do
            end do
            sd(kap, kapp) = sigd

            ! Sum for cross term
            sigde = 0.0d0
            do lam = 0, lmax
              llam = 2 * lam + 1
              do lamp = 0, lmax
                llamp = 2 * lam + 1
                do kap2 = -lmax-1, lmax
                  if (kap2 == 0) cycle
                  jj2 = rj(kap2) + 1
                  l2 = int(jj2 * 0.5d0 + small) - 1
                  sum1 = idgg(kap2, lam) * cdgg(kap2, lam) + &
                         idgf(kap2, lam) * cdgf(kap2, lam) + &
                         idfg(kap2, lam) * cdfg(kap2, lam) + &
                         idff(kap2, lam) * cdff(kap2, lam)
                  sum2 = iegg(kap2, lamp) * cegg(kap2, lamp) + &
                         iegf(kap2, lamp) * cegf(kap2, lamp) + &
                         iefg(kap2, lamp) * cefg(kap2, lamp) + &
                         ieff(kap2, lamp) * ceff(kap2, lamp)
                  select case (jj - 1)
                  case (1)
                    contr = jj1 * jj2 * sum1 * sum2 * w6j1(lp, l2, lam, lamp)
                  case (3)
                    contr = jj1 * jj2 * sum1 * sum2 * w6j3(lp, l2, lam, lamp)
                  case (5)
                    contr = jj1 * jj2 * sum1 * sum2 * w6j5(lp, l2, lam, lamp)
                  case (7)
                    contr = jj1 * jj2 * sum1 * sum2 * w6j7(lp, l2, lam, lamp)
                  case default
                    contr = 0.0d0
                  end select
                  sigde = sigde + contr
                end do
              end do
            end do
            sde(kap, kapp) = sigde

            ! Update shared array in critical section
            !$omp critical
            ssum(kap, kapp, nei) = 2.0d0 * (sigd - sigde)
            !$omp end critical
          end do
        end do
      end do

      ! Construct symmetrized cross sections (parallelized)
      !$omp do schedule(dynamic)
      do i = 1, ne
        ip = nei2 + 1 - i
        if (irange(i, nei2) == 0) cycle
        do kap = -lval-1, lval
          if (kap == 0) cycle
          crsecs(kap, kap, i, nei2) = 0.5d0 * (ssum(kap, kap, i) + ssum(kap, kap, ip))
          do kapp = -lval-1, kap-1
            if (kapp == 0) cycle
            crsecs(kap, kapp, i, nei2) = ssum(kap, kapp, i) + ssum(kapp, kap, ip)
          end do
        end do
      end do
      !$omp end do
    end do
    !$omp end parallel do

    ! Write output
    write(mat_unit, '(1x,a/40a1/1x,a/5a1/1x,a/f12.5/1x,a/2f6.2,i5/ &
      &1x,a/i1/1x,a/f8.6/)') &
      'RELATIVISTIC CORE-VALENCE-VALENCE AUGER MATRIXELEMENTS', &
      name, &
      'CORE STATE', &
      norb1, &
      'ENERGY OF THE CORE STATE', ec1, &
      'E0, DE, NE FOR VALENCE BAND', e0, de, ne, &
      'MAXIMAL L QUANTUMNUMBER FOR VALENCE STATES', lval, &
      'LATTICE CONSTANT', a0

    write(mat_unit, '(1x,a2,2x,a1,3x,a2,4x,a2,4x,a1,2x,28(i6,i3,4x))') &
      'i2', 'i', 'e2', 'ep', 'e', &
      (((kap, kapp), kapp=-lval-1,kap), kap=-lval-1,-1), &
      (((kap, kapp), kapp=-lval-1,-1), ((kap, kapp), kapp=1,kap), kap=1,lval)

    do i2 = 1, ne2
      do i = 1, ne
        ip = i2 + 1 - i
        if (irange(i, i2) == 0) cycle
        write(mat_unit, '(2i3,3f6.2,28d13.5)') &
          i2, i, earr2(i2), earr(ip), earr(i), &
          ((crsecs(kap, kapp, i, i2), kapp=-lval-1,kap), kap=-lval-1,-1), &
          ((crsecs(kap, kapp, i, i2), kapp=-lval-1,-1), &
           (crsecs(kap, kapp, i, i2), kapp=1,kap), kap=1,lval)
      end do
    end do

    ! Close files
    close(pun_unit)
    close(po_unit)
    close(wf_unit)
    close(mat_unit)
  end do

  ! Clean up
  close(in_unit)
  close(9)
  stop
contains
  subroutine triadini(triad, rl, rlb, rj, pun)
    use parameters
    implicit none
    real(8), intent(out) :: triad(-lmax-1:lmax, -lmax-1:lmax, 0:lmax2)
    integer, intent(in) :: rl(-lmax-1:lmax), rlb(-lmax-1:lmax), rj(-lmax-1:lmax)
    integer, intent(in) :: pun
    real(8) :: thrcof(ndim), sixcof(ndim)
    real(8) :: nul, j1, j2, xl1, xl2, fact, xlmi1, xlma1, xlmat
    real(8) :: xlmi2, xlma2
    integer :: kap1, kap2, lam, l1, l2, lami1, lama1, lami2, lama2
    integer :: lami, lama, i1, i2, ier
    integer :: iwigtes = 0

    nul = 0.0d0
    triad = nul

    do kap1 = -lmax-1, lmax
      if (kap1 == 0) cycle
      l1 = rl(kap1)
      xl1 = real(l1, 8)
      j1 = rj(kap1) * half
      do kap2 = kap1, lmax
        if (kap2 == 0) cycle
        l2 = rl(kap2)
        xl2 = real(l2, 8)
        j2 = rj(kap2) * half
        fact = sqrt((2.0d0 * xl1 + 1.0d0) * (2.0d0 * xl2 + 1.0d0))

        call rec3jj(thrcof, xl1, xl2, nul, nul, xlmi1, xlma1, xlmat, ndim, ier)
        if (ier < 0) cycle

        call rec6j(sixcof, xl1, xl2, half, j2, j1, xlmi2, xlma2, xlmat, ndim, ier)
        if (ier < 0) cycle

        lami1 = int(xlmi1 + small)
        lami2 = int(xlmi2 + small)
        lama1 = int(xlma1 + small)
        lama2 = int(xlma2 + small)
        lami = max(lami1, lami2)
        lama = min(lama1, lama2)

        do lam = lami, lama
          if (mod(lam - lami1, 2) == 1) cycle
          i1 = lam - lami1 + 1
          i2 = lam - lami2 + 1
          triad(kap1, kap2, lam) = fact * thrcof(i1) * sixcof(i2)
          triad(kap2, kap1, lam) = triad(kap1, kap2, lam)
          if (iwigtes == 1) then
            write(pun, '(2i4,f5.1,5x,2i4,f5.1,5x,i4,10x,d17.10)') &
              kap1, l1, j1, kap2, l2, j2, lam, triad(kap1, kap2, lam)
          end if
        end do
      end do
    end do
  end subroutine triadini

  subroutine wigini(w6j, j1, j3, pun)
    use parameters
    implicit none
    real(8), intent(out) :: w6j(0:lvmax, 0:lmax, 0:lmax, 0:lmax)
    real(8), intent(in) :: j1, j3
    integer, intent(in) :: pun
    real(8) :: sixcof(ndim)
    real(8) :: nul, j, j2, xlam, xlammi, xlamma, am
    integer :: l, l2, lam, lamp, lami1, lama1, lami2, lama2, lami, lama
    integer :: lammi, lamma, ier, n
    integer :: iwigtes = 0

    nul = 0.0d0
    w6j = nul

    lami1 = int(abs(j1 - j3) + small)
    lama1 = int(j1 + j3 + small)

    do l = 0, lvmax
      j = (2.0d0 * l + 1.0d0) * half
      do l2 = 0, lmax
        j2 = (2.0d0 * l2 + 1.0d0) * half
        lami2 = int(abs(j - j2) + small)
        lama2 = int(j + j2 + small)
        lami = max(lami1, lami2)
        lama = min(lama1, lama2)
        do lam = 0, lmax
          if (lam < lami .or. lam > lama) cycle
          xlam = real(lam, 8)
          call rec6j(sixcof, j1, j, xlam, j2, j3, xlammi, xlamma, xlmat, ndim, ier)
          if (ier < 0) cycle
          lammi = int(xlammi + small)
          lamma = int(xlamma + small)
          do lamp = lammi, lamma
            n = lam + lamp + 1
            am = real(mod(n, 2), 8)
            w6j(l, l2, lam, lamp) = (1.0d0 - 2.0d0 * am) * sixcof(lamp - lammi + 1)
            if (iwigtes == 1) then
              write(pun, '(i4,2f5.1,5x,i4,2f5.1,10x,d17.10)') &
                lamp, j1, j, lam, j2, j3, w6j(l, l2, lam, lamp)
            end if
          end do
        end do
      end do
    end do
  end subroutine wigini

  subroutine wafu(en, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, &
                  p, q, norm, nval, nonrel, pun, lval, coef, kout)
    use parameters
    implicit none
    real(8), intent(in) :: en, x0, dx, rmt, rws, coef
    integer, intent(in) :: nmt, nws, norm, nval, nonrel, pun, lval, kout
    integer, intent(in) :: rl(-lmax-1:lmax), rlb(-lmax-1:lmax), rj(-lmax-1:lmax)
    real(8), intent(in) :: v(nrad), rs(nrad)
    real(8), intent(out) :: p(nrad, -lmax-1:lmax), q(nrad, -lmax-1:lmax)
    real(8) :: rr(nrad), vint(-lmax-1:lmax), eta(-lmax-1:lmax)
    real(8) :: tanx(-lmax-1:lmax), ratx(-lmax-1:lmax)
    real(8) :: fb(0:lmax+1), fn(0:lmax+1), fb1(0:lmax+1), fn1(0:lmax+1)
    real(8) :: clight, sign, ekappa, yk, yk2, sk, tane, a1, b1, cose, sine
    real(8) :: r, x1, x2, qint
    integer :: kap, l, lb, lact, i

    if (nonrel == 1) then
      clight = chuge
    else
      clight = csmall
    end if

    if (en > 0.0d0) then
      sign = 1.0d0
      ekappa = sqrt(en)
    else
      sign = -1.0d0
      ekappa = sqrt(-en)
    end if

    call sbf1(en, rmt, fb, fn)

    if (nval == 1) then
      lact = lval
    else
      lact = lmax
    end if

    do kap = -lact-1, lact
      if (kap == 0) cycle
      l = rl(kap)
      lb = rlb(kap)
      if (kap > 0) then
        sk = ekappa
      else
        sk = -sign * ekappa
      end if
      yk = sqrt(kap * kap - coef)
      if (nonrel == 1) then
        yk = real(l + 1, 8)
        lb = l + 1
        sk = -sign * ekappa
      end if
      yk2 = yk + yk

      call comdir(en, kap, v, nmt, nws, dx, x0, q(:, kap), p(:, kap), ratx(kap), nonrel)

      tane = (ratx(kap) * fb(l) - sk * fb(lb)) / (ratx(kap) * fn(l) - sk * fn(lb))
      tanx(kap) = tane

      if (nval == 1) then
        a1 = sk * ekappa * (fn(lb) - fb(lb) / tane) * rmt / q(nmt, kap) / clight
        b1 = ekappa * (fn(l) - fb(l) / tane) * rmt / p(nmt, kap)
      else
        eta(kap) = atan(tane)
        cose = cos(eta(kap))
        sine = sin(eta(kap))
        a1 = sk * (cose * fb(lb) - sine * fn(lb)) * rmt / q(nmt, kap) / clight
        b1 = (cose * fb(l) - sine * fn(l)) * rmt / p(nmt, kap)
      end if

      do i = 1, nmt
        p(i, kap) = p(i, kap) * b1
        q(i, kap) = q(i, kap) * a1
        rr(i) = p(i, kap) * p(i, kap) + q(i, kap) * q(i, kap)
      end do

      do i = nmt + 1, nws + 1
        r = rs(i)
        call sbf1(en, r, fb1, fn1)
        if (nval == 1) then
          q(i, kap) = sk * ekappa * (fn1(lb) - fb1(lb) / tane) * r / clight
          p(i, kap) = ekappa * (fn1(l) - fb1(l) / tane) * r
        else
          q(i, kap) = sk * (cose * fb1(lb) - sine * fn1(lb)) * r / clight
          p(i, kap) = (cose * fb1(l) - sine * fn1(l)) * r
        end if
        rr(i) = p(i, kap) * p(i, kap) + q(i, kap) * q(i, kap)
      end do

      if (nval == 2) cycle

      if (norm == 0) then
        vint(kap) = sintg(yk2, rr, rs, dx, nmt)
      else
        x1 = sintg(yk2, rr, rs, dx, nws - 1)
        x2 = sintg(yk2, rr, rs, dx, nws)
        vint(kap) = x1 + (x2 - x1) * (rws - rs(nws - 1)) / (rs(nws) - rs(nws - 1))
      end if
      qint = sqrt(1.0d0 / vint(kap))
      p(:, kap) = p(:, kap) * qint
      q(:, kap) = q(:, kap) * qint
    end do

    if (kout == 1) then
      if (nval == 1) then
        write(pun, '(/1x,a)') 'CF/G RATIO'
        write(pun, '(2x,12(5x,i3,5x))') (kap, kap=-lact-1,-1), (kap, kap=1,lact)
        write(pun, '(2x,12d13.5)') (ratx(kap), kap=-lact-1,-1), (ratx(kap), kap=1,lact)
        write(pun, '(/1x,a)') 'TANGENT PHASESHIFTS'
        write(pun, '(2x,12(5x,i3,5x))') (kap, kap=-lact-1,-1), (kap, kap=1,lact)
        write(pun, '(2x,12d13.5)') (tanx(kap), kap=-lact-1,-1), (tanx(kap), kap=1,lact)
      else
        write(pun, '(/1x,a)') 'CF/G RATIO'
        write(pun, '(2x,12(5x,i3,5x))') (kap, kap=-1,-lact-1,-1)
        write(pun, '(2x,12d13.5)') (ratx(kap), kap=-1,-lact-1,-1)
        write(pun, '(15x,11(6x,i2,5x))') (kap, kap=1,lact)
        write(pun, '(15x,11d13.5)') (ratx(kap), kap=1,lact)
        write(pun, '(/1x,a)') 'PHASESHIFTS'
        write(pun, '(2x,12(5x,i3,5x))') (kap, kap=-1,-lact-1,-1)
        write(pun, '(2x,12d13.5)') (eta(kap), kap=-1,-lact-1,-1)
        write(pun, '(15x,11(6x,i2,5x))') (kap, kap=1,lact)
        write(pun, '(15x,11d13.5)') (eta(kap), kap=1,lact)
      end if
    end if
  end subroutine wafu

  subroutine comdir(e1, kappa, za, nrc, nnk, dx, x0, q, p, ratfg, nonrel)
    use parameters
    implicit none
    real(8), intent(in) :: e1, dx, x0, za(nrad)
    integer, intent(in) :: kappa, nrc, nnk, nonrel
    real(8), intent(out) :: q(nrad), p(nrad), ratfg
    real(8) :: bgx(nrad), sxk(4), sxm(4), pp(nrad), qp(nrad)
    real(8) :: test, c, cin, hoc, u, tc, xk, x, xc, bgc, wc, uc, t
    integer :: pun, jri, n, ik, lkap, i
    test = 1.0d5
    pun = 99

    bgx(1:nnk) = za(1:nnk)
    xk = real(kappa, 8)
    jri = nrc
    x = x0
    tc = exp(x)
    bgx(1) = -bgx(1) * tc

    if (nonrel == 1) then
      if (kappa < 0) then
        lkap = -kappa - 1
      else
        lkap = kappa
      end if
      xk = real(-lkap - 1, 8)
      u = -0.5d0 * bgx(1) / (lkap + 1.0d0)
      p(1) = 1.0d-20
      q(1) = u * 1.0d-20
    else
      c = csmall
      cin = 1.0d0 / (c * c)
      hoc = bgx(1) / c
      if (abs(hoc / xk) <= 0.05d0) then
        u = (xk + abs(xk)) / hoc - 0.5d0 * hoc / abs(xk)
      else
        u = (xk + sqrt(xk * xk - hoc * hoc)) / hoc
      end if
      p(1) = 1.0d-20
      q(1) = c * u * 1.0d-20
    end if

    if (nonrel /= 1) then
      pp(1) = tc * (cin * (e1 - bgx(1)) + 1.0d0) * q(1) - xk * p(1)
    else
      pp(1) = tc * q(1) - xk * p(1)
    end if
    qp(1) = xk * q(1) - tc * (e1 - bgx(1)) * p(1)

    n = 1
    do while (n < 6)
      ik = 0
      xc = x
      bgc = bgx(n)
      wc = q(n)
      uc = p(n)
      do
        ik = ik + 1
        t = exp(xc)
        if (nonrel /= 1) then
          sxk(ik) = dx * (-xk * uc + t * wc * (cin * (e1 - bgc) + 1.0d0))
        else
          sxk(ik) = dx * (-xk * uc + t * wc)
        end if
        sxm(ik) = dx * (xk * wc - t * (e1 - bgc) * uc)
        select case (ik)
        case (1)
          xc = xc + 0.5d0 * dx
          uc = uc + 0.5d0 * sxk(1)
          wc = wc + 0.5d0 * sxm(1)
          bgc = 0.5d0 * (bgc + bgx(n + 1))
        case (2)
          uc = uc + 0.5d0 * (sxk(2) - sxk(1))
          wc = wc + 0.5d0 * (sxm(2) - sxm(1))
        case (3)
          xc = xc + 0.5d0 * dx
          uc = uc + sxk(3) - 0.5d0 * sxk(2)
          wc = wc + sxm(3) - 0.5d0 * sxm(2)
          bgc = bgx(n + 1)
        case (4)
          q(n + 1) = q(n) + (sxm(1) + 2.0d0 * sxm(2) + 2.0d0 * sxm(3) + sxm(4)) / 6.0d0
          p(n + 1) = p(n) + (sxk(1) + 2.0d0 * sxk(2) + 2.0d0 * sxk(3) + sxk(4)) / 6.0d0
          if (nonrel /= 1) then
            pp(n + 1) = t * q(n + 1) * (cin * (e1 - bgc) + 1.0d0) - xk * p(n + 1)
          else
            pp(n + 1) = t * q(n + 1) - xk * p(n + 1)
          end if
          qp(n + 1) = xk * q(n + 1) - t * (e1 - bgc) * p(n + 1)
          exit
        end select
      end do
      x = x + dx
      n = n + 1
    end do

    do while (n < nnk)
      x = x + dx
      t = exp(x)
      unp = p(n - 5) + 0.3d0 * dx * (11.0d0 * pp(n) - 14.0d0 * pp(n - 1) + &
            26.0d0 * pp(n - 2) - 14.0d0 * pp(n - 3) + 11.0d0 * pp(n - 4))
      wnp = q(n - 5) + 0.3d0 * dx * (11.0d0 * qp(n) - 14.0d0 * qp(n - 1) + &
            26.0d0 * qp(n - 2) - 14.0d0 * qp(n - 3) + 11.0d0 * qp(n - 4))
      nit = 0
      do
        if (nonrel /= 1) then
          pp(n + 1) = t * (cin * (e1 - bgx(n + 1)) + 1.0d0) * wnp - xk * unp
        else
          pp(n + 1) = t * wnp - xk * unp
        end if
        qp(n + 1) = xk * wnp - t * (e1 - bgx(n + 1)) * unp
        unp2 = p(n - 3) + (7.0d0 * pp(n + 1) + 32.0d0 * pp(n) + &
               12.0d0 * pp(n - 1) + 32.0d0 * pp(n - 2) + 7.0d0 * pp(n - 3)) * &
               2.0d0 * dx / 45.0d0
        wnp2 = q(n - 3) + (7.0d0 * qp(n + 1) + 32.0d0 * qp(n) + &
               12.0d0 * qp(n - 1) + 32.0d0 * qp(n - 2) + 7.0d0 * qp(n - 3)) * &
               2.0d0 * dx / 45.0d0
        if (abs(test * (unp2 - unp)) <= abs(unp2) .and. &
            abs(test * (wnp2 - wnp)) <= abs(wnp2)) then
          q(n + 1) = wnp2
          p(n + 1) = unp2
          exit
        end if
        nit = nit + 1
        if (nit > 5) then
          q(n + 1) = wnp2
          p(n + 1) = unp2
          exit
        end if
        wnp = wnp2
        unp = unp2
      end do
      n = n + 1
    end do
    ratfg = q(jri) / p(jri)
  end subroutine comdir

  subroutine cei(f2, f4, f1, f3, r, n, lambda, l2, l4, l1, l3, dx, rws, result)
    use parameters
    implicit none
    real(8), intent(in) :: f1(nrad), f2(nrad), f3(nrad), f4(nrad), r(nrad)
    integer, intent(in) :: n, lambda
    real(8), intent(in) :: l1, l2, l3, l4, dx, rws
    real(8), intent(out) :: result
    real(8) :: rl(ndim), rmlp1(ndim), arg1(ndim), arg2(ndim), arg3(ndim), arg4(ndim)
    real(8) :: ll
    integer :: m, i

    m = n
    do i = 1, m
      rl(i) = r(i) ** lambda
      rmlp1(i) = 1.0d0 / (rl(i) * r(i))
      arg1(i) = f2(i) * f4(i) * rl(i)
    end do
    ll = l2 + l4 + lambda
    call dintg(ll, 1, arg1, arg2, r, dx, m, rws)
    ll = ll + 1
    arg2(1:m) = arg2(1:m) * f1(1:m) * f3(1:m) * rmlp1(1:m)
    ll = ll + l1 + l3 - lambda - 1
    do i = 1, m
      arg3(i) = f2(i) * f4(i) * rmlp1(i)
    end do
    ll = l2 + l4 - lambda - 1
    call dintg(ll, -1, arg3, arg4, r, dx, m, rws)
    ll = ll + 1
    arg4(1:m) = arg4(1:m) * f1(1:m) * f3(1:m) * rl(1:m)
    ll = ll + l1 + l3 + lambda
    arg2(1:m) = arg2(1:m) + arg4(1:m)
    result = sintg(ll, arg2, r, dx, m)
  end subroutine cei

  real(8) function sintg(ll, fct, r, dx, n)
    use parameters
    implicit none
    real(8), intent(in) :: ll, fct(nrad), r(nrad), dx
    integer, intent(in) :: n
    real(8) :: sum, rll, fact, x0, x1, x2
    integer :: i

    if (ll + 1.0d0 <= tiny) then
      write(9, *) ' sintg: ll+1=', ll + 1
      sum = 0.0d0
    else
      rll = r(1) ** ll
      if (rll <= tiny) then
        sum = 0.0d0
      else
        fact = fct(1) / rll
        sum = fact * (r(1) ** (ll + 1)) / (ll + 1)
      end if
    end if
    x0 = fct(1) * r(1)
    do i = 3, n, 2
      x1 = fct(i - 1) * r(i - 1)
      x2 = fct(i) * r(i)
      sum = sum + dx * (x0 + 4.0d0 * x1 + x2)
      x0 = x2
    end do
    sintg = sum / 3.0d0
    if (mod(n, 2) == 0) then
      sintg = sintg + (fct(n) + fct(n - 1)) / 2.0d0 * (r(n) - r(n - 1))
    end if
  end function sintg

  subroutine dintg(ll, idir, fct, yint, r, dx, n, rws)
    use parameters
    implicit none
    real(8), intent(in) :: ll, fct(nrad), r(nrad), dx, rws
    integer, intent(in) :: idir, n
    real(8), intent(out) :: yint(nrad)
    real(8) :: sum, rll, fact, x0, x1, corr
    integer :: i, j

    if (idir > 0) then
      if (ll + 1.0d0 <= tiny) then
        write(9, *) ' dintg: ll+1=', ll + 1
        sum = 0.0d0
      else
        rll = r(1) ** ll
        if (rll <= tiny) then
          sum = 0.0d0
        else
          fact = fct(1) / rll
          sum = fact * (r(1) ** (ll + 1)) / (ll + 1)
        end if
      end if
      yint(1) = sum
      x0 = fct(1) * r(1)
      do i = 2, n
        x1 = fct(i) * r(i)
        sum = sum + dx * (x0 + x1) / 2.0d0
        yint(i) = sum
        x0 = x1
      end do
    else
      yint(n) = 0.0d0
      sum = 0.0d0
      x1 = fct(n) * r(n)
      do i = n - 1, 1, -1
        x0 = fct(i) * r(i)
        sum = sum + dx * (x0 + x1) / 2.0d0
        yint(i) = sum
        if (i == n - 3) then
          call inter1(r(n - 3), yint(n - 3), 4, 1, rws, corr)
          sum = sum - corr
          yint(n - 3:n) = yint(n - 3:n) - corr
        end if
        x1 = x0
      end do
    end if
  end subroutine dintg

  subroutine inter1(r, p, n, id, rs, ps)
    implicit none
    real(8), intent(in) :: r(n), p(n), rs
    integer, intent(in) :: n, id
    real(8), intent(out) :: ps
    real(8) :: term, denom
    integer :: i, j

    ps = 0.0d0
    do j = 1, n, id
      term = 1.0d0
      denom = 1.0d0
      do i = 1, n, id
        if (i == j) cycle
        denom = denom * (r(j) - r(i))
        term = term * (rs - r(i))
      end do
      ps = ps + term * p(j) / denom
    end do
  end subroutine inter1

  subroutine sbf1(e, r, fb, fn)
    use parameters
    implicit none
    real(8), intent(in) :: e, r
    real(8), intent(out) :: fb(0:lmax+1), fn(0:lmax+1)
    real(8) :: eroot, x1, x2, sinx, cosx, tl
    integer :: l

    if (e > 0.0d0) then
      eroot = sqrt(e)
      x1 = eroot * r
      x2 = x1 * x1
      sinx = sin(x1)
      cosx = cos(x1)
      fb(0) = sinx / x1
      fb(1) = sinx / x2 - cosx / x1
      fn(0) = -cosx / x1
      fn(1) = -cosx / x2 - sinx / x1
      do l = 2, lmax + 1
        tl = real(2 * l - 1, 8)
        fb(l) = tl * fb(l - 1) / x1 - fb(l - 2)
        fn(l) = tl * fn(l - 1) / x1 - fn(l - 2)
      end do
    else
      eroot = sqrt(-e)
      x1 = eroot * r
      x2 = x1 * x1
      sinx = sinh(x1)
      cosx = cosh(x1)
      fb(0) = sinx / x1
      fb(1) = -sinx / x2 + cosx / x1
      fn(0) = cosx / x1
      fn(1) = -cosx / x2 + sinx / x1
      do l = 2, lmax + 1
        tl = real(2 * l - 1, 8)
        fb(l) = -tl * fb(l - 1) / x1 + fb(l - 2)
        fn(l) = -tl * fn(l - 1) / x1 + fn(l - 2)
      end do
    end if
  end subroutine sbf1

  subroutine sbf(e, r, fb, fn)
    use parameters
    implicit none
    real(8), intent(in) :: e, r
    real(8), intent(out) :: fb(0:lmax+1), fn(0:lmax+1)
    complex(8) :: cfb(0:lmax+1), cfn(0:lmax+1)
    complex(8) :: ecompl, eroot, x1, sinx, cosx, cimmi
    real(8) :: x2, tl
    integer :: l

    cimmi = cmplx(0.0d0, -1.0d0, 8)
    ecompl = cmplx(e, 0.0d0, 8)
    eroot = sqrt(ecompl)
    x1 = eroot * r
    x2 = real(x1 * x1, 8)
    sinx = sin(x1)
    cosx = cos(x1)
    cfb(0) = sinx / x1
    cfb(1) = sinx / x2 - cosx / x1
    cfn(0) = -cosx / x1
    cfn(1) = -cosx / x2 - sinx / x1
    do l = 2, lmax + 1
      tl = real(2 * l - 1, 8)
      cfb(l) = tl * cfb(l - 1) / x1 - cfb(l - 2)
      cfn(l) = tl * cfn(l - 1) / x1 - cfn(l - 2)
    end do
    if (e < 0.0d0) then
      cfn(0) = cfn(0) * cimmi
      do l = 1, lmax + 1
        cfb(l) = cfb(l) * cimmi ** l
        cfn(l) = cfn(l) * cimmi ** (l + 1)
      end do
    end if
    fb = real(cfb, 8)
    fn = real(cfn, 8)
  end subroutine sbf
end program rcvvmat

! External subroutine interfaces
interface
  subroutine rec3jj(thrcof, l2, l3, m2, m3, l1min, l1max, lmatch, ndim, ier)
    implicit none
    real(8), intent(out) :: thrcof(*)
    real(8), intent(in) :: l2, l3, m2, m3
    real(8), intent(out) :: l1min, l1max, lmatch
    integer, intent(in) :: ndim
    integer, intent(out) :: ier
  end subroutine rec3jj

  subroutine rec6j(sixcof, l2, l3, l4, l5, l6, l1min, l1max, lmatch, ndim, ier)
    implicit none
    real(8), intent(out) :: sixcof(*)
    real(8), intent(in) :: l2, l3, l4, l5, l6
    real(8), intent(out) :: l1min, l1max, lmatch
    integer, intent(in) :: ndim
    integer, intent(out) :: ier
  end subroutine rec6j
end interface
