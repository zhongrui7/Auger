module parameters
  implicit none
  ! Parameters from l.par (assumed values; adjust as needed)
  integer, parameter :: lmax = 6
  integer, parameter :: lvmax = 3
  integer, parameter :: lmax2 = 2 * lmax + 2
  ! Common constants
  integer, parameter :: nrad = 250
  integer, parameter :: nemax = 50
  integer, parameter :: ndim = 100
  real(kind=8), parameter :: tiny = 1.0d-10
  real(kind=8), parameter :: small = 0.01d0
  real(kind=8), parameter :: half = 0.5d0
  real(kind=8), parameter :: pi = 3.1415926535897932d0
  real(kind=8), parameter :: ev = 13.606d0
  real(kind=8), parameter :: csmall = 274.0746d0
  real(kind=8), parameter :: chuge = 1.0d8
  real(kind=8), parameter :: zero = 0.0d0
  real(kind=8), parameter :: one = 1.0d0
  real(kind=8), parameter :: two = 2.0d0
  real(kind=8), parameter :: three = 3.0d0
  real(kind=8), parameter :: srtiny = 1.0d-5
  real(kind=8), parameter :: huge = 1.0d10
  real(kind=8), parameter :: srhuge = 1.0d5
end module parameters

program rccvmat
  use parameters
  implicit none

  ! Arrays
  real(kind=8) :: v(nrad), rs(nrad), xs(nrad), earr(nemax)
  real(kind=8) :: gc1(nrad), fc1(nrad), gc3(nrad), fc3(nrad)
  real(kind=8) :: gvwf(nrad, -lmax-1:lmax), fvwf(nrad, -lmax-1:lmax)
  real(kind=8) :: gcwf(nrad, -lmax-1:lmax), fcwf(nrad, -lmax-1:lmax)
  real(kind=8) :: sd(-lvmax-1:lvmax, nemax), se(-lvmax-1:lvmax, nemax)
  real(kind=8) :: sde(-lvmax-1:lvmax, nemax), ssum(-lvmax-1:lvmax, nemax)
  real(kind=8) :: zd(-lvmax-1:lvmax, nemax), ze(-lvmax-1:lvmax, nemax)
  real(kind=8) :: zde(-lvmax-1:lvmax, nemax), zsum(-lvmax-1:lvmax, nemax)
  real(kind=8) :: idgg(-lmax-1:lmax, 0:lmax), idgf(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: idfg(-lmax-1:lmax, 0:lmax), idff(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: iegg(-lmax-1:lmax, 0:lmax), iegf(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: iefg(-lmax-1:lmax, 0:lmax), ieff(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: cdgg(-lmax-1:lmax, 0:lmax), cdgf(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: cdfg(-lmax-1:lmax, 0:lmax), cdff(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: cegg(-lmax-1:lmax, 0:lmax), cegf(-lmax-1:lmax, 0:lmax)
  real(kind=8) :: cefg(-lmax-1:lmax, 0:lmax), ceff(-lmax-1:lmax, 0:lmax)
  integer :: rl(-lmax-1:lmax), rlb(-lmax-1:lmax), rj(-lmax-1:lmax)
  integer :: triad(-lmax-1:lmax, -lmax-1:lmax, 0:lmax2)
  real(kind=8) :: w6j(0:lvmax, 0:lmax, 0:lmax, 0:lmax)

  ! Scalars
  integer :: in = 1, wf = 2, po = 3, mat = 7, pun = 8
  integer :: zatom, iwigtes, lval, kap1, kap3, nonrel, norm, nmt, nws
  integer :: ne, nei, kap, lam, kap2, jj, l, lkap, llam, jj2, l2, lamp, llamp
  real(kind=8) :: e0, de, clight, x0, dx, rws, a0, zeromt, coef
  real(kind=8) :: x, signv, rmt, ec1, ec3, e, e2, dec, y, y1, y3
  real(kind=8) :: jj1, jj3, xj1, xj3, y2, t13, tm1m3, t1v, tm1mv, tv2
  real(kind=8) :: tmvm2, t32, tm3m2, result, sigd, sige, sigde, sum1, sum2
  real(kind=8) :: contr, result1, result2, dummy
  character(len=1) :: name(5), norb1(5), norb3(5), norb(5)
  character(len=30) :: file

  ! Initialize arrays
  v = 0.0d0; rs = 0.0d0; xs = 0.0d0; earr = 0.0d0
  gc1 = 0.0d0; fc1 = 0.0d0; gc3 = 0.0d0; fc3 = 0.0d0
  gvwf = 0.0d0; fydia = 0.0d0; gcwf = 0.0d0; fcwf = 0.0d0
  sd = 0.0d0; se = 0.0d0; sde = 0.0d0; ssum = 0.0d0
  zd = 0.0d0; ze = 0.0d0; zde = 0.0d0; zsum = 0.0d0
  idgg = 0.0d0; idgf = 0.0d0; idfg = 0.0d0; idff = 0.0d0
  iegg = 0.0d0; iegf = 0.0d0; iefg = 0.0d0; ieff = 0.0d0
  cdgg = 0.0d0; cdgf = 0.0d0; cdfg = 0.0d0; cdff = 0.0d0
  cegg = 0.0d0; cegf = 0.0d0; cefg = 0.0d0; ceff = 0.0d0
  rl = 0; rlb = 0; rj = 0
  triad = 0; w6j = 0.0d0

  ! Open input and error files
  open(unit=in, file='rccvmat.in', status='old')
  open(unit=9, file='rccvmat.err', status='unknown')

  ! Read output file name and open it
  read(in, '(a30)') file
  open(unit=pun, file=trim(file), status='unknown')

  ! Fill relativistic quantum numbers
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

  ! Read control parameter for Wigner coefficients
  read(in, *) iwigtes
  if (iwigtes == 1) then
     read(in, *) lval, kap1, kap3
     goto 499
  end if

  ! Read maximum angular momentum quantum number for valence states
  read(in, *) lval
  write(pun, 200) lval
200 format(' LVAL=', t20, i4, /)

  ! Read energy panel parameters and set up energy array
  read(in, *) e0, de, ne
  do nei = 1, ne
     earr(nei) = e0 + (nei - 1) * de
  end do
  write(pun, 201) (earr(nei), nei=1, ne)
201 format(' ENERGY ARRAY:', /, (10f7.2), /)

  ! Check for non-relativistic limit
  read(in, *) nonrel
  if (nonrel == 1) then
     clight = chuge
  else
     clight = csmall
  end if

  ! Read normalization parameter for valence wavefunctions
  read(in, *) norm
  write(pun, 202) norm
202 format(' NORMALIZATION OF VALENCE WF. (0 FOR MT, 1 FOR WS):', i4, /)

  ! Read core state identifiers
  read(in, '(5a1)') norb1
  read(in, '(5a1)') norb3

  ! Read potential file and parameters
  read(in, '(a30)') file
  open(unit=po, file=trim(file), status='old')
  read(po, '(5a1)') name
  read(po, *) zatom, x0, nmt, dx, rws, a0, zeromt
  read(po, *) v(1:nmt)

  coef = 2.0d0 * zatom / clight
  coef = coef * coef

  ! Set up logarithmic and radial mesh, transform potential
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

  write(pun, 203) zatom, x0, dx, nmt, rmt, nws, rws, zeromt
203 format(' ZATOM =', t20, i4, /, &
       '    X0 =', t20, f5.1, /, &
       '    DX =', t20, f10.7, /, &
       '   NMT =', t20, i4, /, &
       '   RMT =', t20, f10.7, /, &
       '   NWS =', t20, i4, /, &
       '   RWS =', t20, f10.7, /, &
       'ZEROMT =', t20, f5.1, //, &
       'RADIAL MESH AND POTENTIAL*R', /)
  write(pun, 204) (rs(j), v(j) * rs(j), j=1, nws)
204 format(4d20.10)

  ! Read core wavefunctions
  read(in, '(a30)') file
  open(unit=wf, file=trim(file), status='old')

21 continue
  read(wf, '(5a1)') norb
  read(wf, *) ec1
  read(wf, *) kap1
  read(wf, *) nwf
  do i = 1, nwf
     if (nonrel == 1) then
        read(wf, *) dummy, gc1(i)
        fc1(i) = 0.0d0
     else
        read(wf, *) dummy, gc1(i), fc1(i)
     end if
  end do
  do l = 1, 5
     if (norb(l) /= norb1(l)) goto 21
  end do

22 continue
  read(wf, '(5a1)') norb(1:3)
  read(wf, *) ec3
  read(wf, *) kap3
  read(wf, *) nwf
  do i = 1, nwf
     if (nonrel == 1) then
        read(wf, *) dummy, gc3(i)
        fc3(i) = 0.0d0
     else
        read(wf, *) dummy, gc3(i), fc3(i)
     end if
  end do
  do l = 1, 5
     if (norb(l) /= norb3(l)) goto 22
  end do

  write(pun, 205) ec1, kap1
205 format(//, ' EC1 =', t20, f12.6, /, &
          'KAP1 =', t20, i4 4, //, &
          'CORE WAVEFUNCTION NO. 1', /)
  write(pun, 206) (rs(j), gc1(j), fc1(j), j=1, nws)
206 format(3d20.10)
  write(pun, 207) ec3, kap3
207 format(//, ' EC3 =', t20, f12.6, /, &
          'KAP3 =', t20, i4, //, &
          'CORE WAVEFUNCTION NO. 3', /)
  write(pun, 206) (rs(j), gc3(j), fc3(j), j=1, nws)

  read(in, '(a30)') file
  open(unit=mat, file=trim(file), status='unknown')

  write(pun, 208)
208 format(//, ' ****************** END INPUT *******************', //)

499 continue
  l1 = rl(kap1)
  l3 = rl(kap3)
  j1 = rj(kap1)
  j3 = rj(kap3)
  jj1 = j1 + 1
  jj3 = j3 + 1
  xj1 = j1 * half
  xj3 = j3 * half
  y1 = sqrt(kap1 * kap1 - coef)
  y3 = sqrt(kap3 * kap3 - coef)
  if (nonrel == 1) y1 = real(l1 + 1, kind=8)
  if (nonrel == 1) y3 = real(l3 + 1, kind=8)

  if (iwigtes == 1) write(pun, 500)
500 format(/, 'TEST FOR TRIADS', //)

  call triadini(triad, rl, rlb, rj, pun, iwigtes)

  if (iwigtes == 1) write(pun, 501)
501 format(//, 'TEST FOR WIGNER 6J SYMBOLS', //)

  call wigini(w6j, xj1, xj3, pun, iwigtes)

  if (iwigtes == 1) stop ' w3j, w6j test o.k.'

  dec = 2.0d0 * (ec3 - ec1)

  do nei = 1, ne
     write(6, *) nei

     ! Compute valence wavefunctions
     e = earr(nei)
     write(pun, *) ' valence wavefunctions'
     write(pun, *) ' energy=', e
     call wafu(e, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, gvwf, fvwf, &
               norm, 1, nonrel, pun, lval, coef)

     ! Compute continuum wavefunctions
     e2 = e + dec
     write(pun, *) ' continuum wavefunctions'
     write(pun, *) ' energy=', e2
     call wafu(e2, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, gcwf, fcwf, &
               norm, 2, nonrel, pun, lval, coef)

     do kap = -lval-1, lval
        if (kap == 0) cycle
        y = sqrt(kap * kap - coef)
        if (nonrel == 1) y = real(rl(kap) + 1, kind=8)
        jj = rj(kap) + 1
        l = int(jj * 0.5 + small) - 1
        lkap = rl(kap)

        write(pun, *) ' kappa=', kap

        do lam = 0, lmax
           do kap2 = -lmax-1, lmax
              idgg(kap2, lam) = 0.0d0
              idfg(kap2, lam) = 0.0d0
              idgf(kap2, lam) = 0.0d0
              idff(kap2, lam) = 0.0d0
              iegg(kap2, lam) = 0.0d0
              iefg(kap2, lam) = 0.0d0
              iegf(kap2, lam) = 0.0d0
              ieff(kap2, lam) = 0.0d0
           end do
        end do

        do lam = 0, lmax
           t13 = triad(kap1, kap3, lam)
           tm1m3 = triad(-kap1, -kap3, lam)
           t1v = triad(kap1, kap, lam)
           tm1mv = triad(-kap1, -kap, lam)

           do kap2 = -lmax-1, lmax
              if (kap2 == 0) cycle
              y2 = sqrt(kap2 * kap2 - coef)
              if (nonrel == 1) y2 = real(rl(kap2) + 1, kind=8)

              tv2 = triad(kap, kap2, lam)
              tmvm2 = triad(-kap, -kap2, lam)
              t32 = triad(kap3, kap2, lam)
              tm3m2 = triad(-kap3, -kap2, lam)

              if (abs(t13 * tv2) > tiny) then
                 call cei(gc1, gc3, gvwf(:, kap), gcwf(:, kap2), rs, nws, lam, &
                          y1, y3, y, y2, dx, rws, result)
                 idgg(kap2, lam) = result
                 cdgg(kap2, lam) = t13 * tv2
              end if

              if (abs(t13 * tmvm2) > tiny) then
                 call cei(gc1, gc3, fvwf(:, kap), fcwf(:, kap2), rs, nws, lam, &
                          y1, y3, y, y2, dx, rws, result)
                 idgf(kap2, lam) = result
                 cdgf(kap2, lam) = t13 * tmvm2
              end if

              if (abs(tm1m3 * tv2) > tiny) then
                 call cei(fc1, fc3, gvwf(:, kap), gcwf(:, kap2), rs, nws, lam, &
                          y1, y3, y, y2, dx, rws, result)
                 idfg(kap2, lam) = result
                 cdfg(kap2, lam) = tm1m3 * tv2
              end if

              if (abs(tm1m3 * tmvm2) > tiny) then
                 call cei(fc1, fc3, fvwf(:, kap), fcwf(:, kap2), rs, nws, lam, &
                          y1, y3, y, y2, dx, rws, result)
                 idff(kap2, lam) = result
                 cdff(kap2, lam) = tm1m3 * tmvm2
              end if

              if (abs(t1v * t32) > tiny) then
                 call cei(gc1, gvwf(:, kap), gc3, gcwf(:, kap2), rs, nws, lam, &
                          y1, y, y3, y2, dx, rws, result)
                 iegg(kap2, lam) = result
                 cegg(kap2, lam) = t1v * t32
              end if

              if (abs(t1v * tm3m2) > tiny) then
                 call cei(gc1, gvwf(:, kap), fc3, fcwf(:, kap2), rs, nws, lam, &
                          y1, y, y3, y2, dx, rws, result)
                 iegf(kap2, lam) = result
                 cegf(kap2, lam) = t1v * tm3m2
              end if

              if (abs(tm1mv * t32) > tiny) then
                 call cei(fc1, fvwf(:, kap), gc3, gcwf(:, kap2), rs, nws, lam, &
                          y1, y, y3, y2, dx, rws, result)
                 iefg(kap2, lam) = result
                 cefg(kap2, lam) = tm1mv * t32
              end if

              if (abs(tm1mv * tm3m2) > tiny) then
                 call cei(fc1, fvwf(:, kap), fc3, fcwf(:, kap2), rs, nws, lam, &
                          y1, y, y3, y2, dx, rws, result)
                 ieff(kap2, lam) = result
                 ceff(kap2, lam) = tm1mv * tm3m2
              end if
           end do
        end do

        ! Sum up direct and exchange terms
        sigd = 0.0d0
        sige = 0.0d0

        do lam = 0, lmax
           llam = 2 * lam + 1
           do kap2 = -lmax-1, lmax
              if (kap2 == 0) cycle
              jj2 = rj(kap2) + 1

              sum1 = idgg(kap2, lam) * cdgg(kap2, lam) + &
                     idgf(kap2, lam) * cdgf(kap2, lam) + &
                     idfg(kap2, lam) * cdfg(kap2, lam) + &
                     idff(kap2, lam) * cdff(kap2, lam)
              sum2 = iegg(kap2, lam) * cegg(kap2, lam) + &
                     iegf(kap2, lam) * cegf(kap2, lam) + &
                     iefg(kap2, lam) * cefg(kap2, lam) + &
                     ieff(kap2, lam) * ceff(kap2, lam)

              sigd = sigd + jj1 * jj2 * jj3 * sum1 * sum1 / llam
              sige = sige + jj1 * jj2 * jj3 * sum2 * sum2 / llam

              if (nonrel == 1 .and. nei == 1) then
                 if (abs(sum1) > tiny) then
                    write(pun, *) ' nonzero for direct term'
                    write(pun, *) lam, kap2, idgg(kap2, lam)
                 end if
                 if (abs(sum2) > tiny) then
                    write(pun, *) ' nonzero for exchange term'
                    write(pun, *) lam, kap2, iegg(kap2, lam)
                 end if
              end if
           end do
        end do

        sd(kap, nei) = sigd
        se(kap, nei) = sige

        ! Sum up cross term
        sigde = 0.0d0
        do lam = 0, lmax
           llam = 2 * lam + 1
           do lamp = 0, lmax
              llamp = 2 * lam + 1
              do kap2 = -lmax-1, lmax
                 if (kap2 == 0) cycle
                 jj2 = rj(kap2) + 1
                 l2 = int(jj2 * 0.5 + small) - 1

                 sum1 = idgg(kap2, lam) * cdgg(kap2, lam) + &
                        idgf(kap2, lam) * cdgf(kap2, lam) + &
                        idfg(kap2, lam) * cdfg(kap2, lam) + &
                        idff(kap2, lam) * cdff(kap2, lam)
                 sum2 = iegg(kap2, lamp) * cegg(kap2, lamp) + &
                        iegf(kap2, lamp) * cegf(kap2, lamp) + &
                        iefg(kap2, lamp) * cefg(kap2, lamp) + &
                        ieff(kap2, lamp) * ceff(kap2, lamp)

                 contr = jj1 * jj2 * jj3 * sum1 * sum2 * w6j(l, l2, lam, lamp)
                 sigde = sigde + contr

                 if (nonrel == 1 .and. nei == 1) then
                    if (abs(contr) > tiny) then
                       write(pun, *) ' nonzero for cross term'
                       write(pun, *) lam, lamp, kap2
                       write(pun, *) idgg(kap2, lam), iegg(kap2, lamp)
                    end if
                 end if
              end do
           end do
        end do

        sde(kap, nei) = sigde
        ssum(kap, nei) = sigd + sige - 2.0d0 * sigde

        ! Non-relativistic limit check
        if (nonrel == 1 .and. kap1 == -1 .and. kap3 == -1) then
           call cei(gc1, gc3, gvwf(:, kap), gcwf(:, kap), rs, nws, 0, &
                    y1, y3, y, y, dx, rws, result1)
           zd(kap, nei) = 2.0d0 * result1 * result1

           call cei(gc1, gvwf(:, kap), gc3, gcwf(:, kap), rs, nws, lkap, &
                    y1, y, y3, y, dx, rws, result2)
           ze(kap, nei) = 2.0d0 * result2 * result2 / (2.0d0 * lkap + 1)**2

           zde(kap, nei) = result1 * result2 / (2.0d0 * lkap + 1)

           zsum(kap, nei) = zd(kap, nei) + ze(kap, nei) - 2.0d0 * zde(kap, nei)
        end if
     end do
  end do

  if (nonrel == 1 .and. kap1 == -1 .and. kap3 == -1) then
     write(pun, 209)
209  format(//, 'CHECK NON-RELATIVISTIC LIMIT', /)
     write(pun, 210)
     write(pun, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
     do i = 1, ne
        write(pun, 222) earr(i), (zd(kap, i), kap=-lval-1, -1), &
                        (zd(kap, i), kap=1, lval)
     end do

     write(pun, 211)
     write(pun, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
     do i = 1, ne
        write(pun, 222) earr(i), (ze(kap, i), kap=-lval-1, -1), &
                        (ze(kap, i), kap=1, lval)
     end do

     write(pun, 212)
     write(pun, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
     do i = 1, ne
        write(pun, 222) earr(i), (zde(kap, i), kap=-lval-1, -1), &
                        (zde(kap, i), kap=1, lval)
     end do
     write(pun, *)
  end if

  write(pun, 210)
210 format(//, 'DIRECT TERM', /)
  write(pun, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
  do i = 1, ne
     write(pun, 222) earr(i), (sd(kap, i), kap=-lval-1, -1), &
                     (sd(kap, i), kap=1, lval)
  end do

  write(pun, 211)
211 format(//, 'EXCHANGE TERM', /)
  write(pun, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
  do i = 1, ne
     write(pun, 222) earr(i), (se(kap, i), kap=-lval-1, -1), &
                     (se(kap, i), kap=1, lval)
  end do

  write(pun, 212)
212 format(//, 'CROSS TERM', /)
  write(pun, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
  do i = 1, ne
     write(pun, 222) earr(i), (sde(kap, i), kap=-lval-1, -1), &
                     (sde(kap, i), kap=1, lval)
  end do

  write(mat, 220) name, norb1, norb3, dec, lval, a0
220 format(' RELATIVISTIC CORE-CORE-VALENCE AUGER MATRIXELEMENTS', /, &
         5a1, /, &
         ' FINAL CORE STATE', /, &
         5a1, /, &
         ' INITIAL CORE STATE', /, &
         5a1, /, &
         ' TRANSITION ENERGY', /, &
         f12.5, /, &
         ' MAXIMAL L QUANTUMNUMBER FOR VALENCE STATES', /, &
         i1, /, &
         ' LATTICE CONSTANT', /, &
         f10.6, /)
  write(mat, 221) (kap, kap=-lval-1, -1), (kap, kap=1, lval)
221 format(3x, 'e', 2x, 12(6x, i2, 5x))
  do i = 1, ne
     write(mat, 222) earr(i), (ssum(kap, i), kap=-lval-1, -1), &
                     (ssum(kap, i), kap=1, lval)
  end do
222 format(f6.2, 12d13.5)

  stop
contains

  subroutine triadini(triad, rl, rlb, rj, pun, iwigtes)
    implicit none
    integer, intent(inout) :: triad(-lmax-1:lmax, -lmax-1:lmax, 0:lmax2)
    integer, intent(in) :: rl(-lmax-1:lmax), rlb(-lmax-1:lmax), rj(-lmax-1:lmax)
    integer, intent(in) :: pun, iwigtes
    real(kind=8) :: thrcof(ndim), sixcof(ndim)
    real(kind=8) :: j1, j2, fact, xl1, xl2, result
    integer :: kap1, kap2, lam, l1, l2, lami1, lami2, lama1, lama2, lami, lama, i1, i2, ier

    triad = 0
    do kap1 = -lmax-1, lmax
       if (kap1 == 0) cycle
       l1 = rl(kap1)
       xl1 = real(l1, kind=8)
       j1 = rj(kap1) * half
       do kap2 = kap1, lmax
          if (kap2 == 0) cycle
          l2 = rl(kap2)
          xl2 = real(l2, kind=8)
          j2 = rj(kap2) * half

          fact = sqrt((2.0d0 * xl1 + 1.0d0) * (2.0d0 * xl2 + 1.0d0))

          call rec3jj(thrcof, xl1, xl2, zero, zero, xlmi1, xlma1, xlmat, ndim, ier)
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
             result = fact * thrcof(i1) * sixcof(i2)
             triad(kap1, kap2, lam) = result
             triad(kap2, kap1, lam) = result

             if (iwigtes == 1) then
                write(pun, 100) kap1, l1, j1, kap2, l2, j2, lam, result * result
             end if
          end do

          if (iwigtes == 1) write(pun, *) '               *****************'
       end do
       if (iwigtes == 1) write(pun, *) '   ******************************************** '
    end do

100 format(2i4, f5.1, 5x, 2i4, f5.1, 5x, i4, 10x, d17.10)
  end subroutine triadini

  subroutine wigini(w6j, j1, j3, pun, iwigtes)
    implicit none
    real(kind=8), intent(out) :: w6j(0:lvmax, 0:lmax, 0:lmax, 0:lmax)
    real(kind=8), intent(in) :: j1, j3
    integer, intent(in) :: pun, iwigtes
    real(kind=8) :: sixcof(ndim), j, j2, xlam, result, am
    integer :: l, l2, lam, lamp, lami1, lama1, lami2, lama2, lami, lama, lammi, lamma, n, ier

    w6j = zero
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
             xlam = real(lam, kind=8)

             call rec6j(sixcof, j1, j, xlam, j2, j3, xlammi, xlamma, xlmat, ndim, ier)
             if (ier < 0) cycle

             lammi = int(xlammi + small)
             lamma = int(xlamma + small)

             do lamp = lammi, lamma
                n = lam + lamp + 1
                am = real(mod(n, 2), kind=8)
                result = sixcof(lamp - lammi + 1)
                w6j(l, l2, lam, lamp) = (1.0d0 - 2.0d0 * am) * result

                if (iwigtes == 1) then
                   write(pun, 100) lamp, j1, j, lam, j2, j3, result
                end if
             end do
          end do
       end do
    end do

100 format(i4, 2f5.1, 5x, i4, 2f5.1, 10x, d17.10)
  end subroutine wigini

  subroutine wafu(en, rl, rlb, rj, v, rs, x0, dx, rmt, nmt, rws, nws, p, q, norm, nval, nonrel, pun, lval, coef)
    implicit none
    real(kind=8), intent(in) :: en, v(nrad), rs(nrad), x0, dx, rmt, rws, coef
    integer, intent(in) :: nmt, nws, norm, nval, nonrel, pun, lval
    integer, intent(in) :: rl(-lmax-1:lmax), rlb(-lmax-1:lmax), rj(-lmax-1:lmax)
    real(kind=8), intent(out) :: p(nrad, -lmax-1:lmax), q(nrad, -lmax-1:lmax)
    real(kind=8) :: rr(nrad), vint(-lmax-1:lmax), eta(-lmax-1:lmax)
    real(kind=8) :: tanx(-lmax-1:lmax), ratx(-lmax-1:lmax)
    real(kind=8) :: fb(0:lmax+1), fn(0:lmax+1), fb1(0:lmax+1), fn1(0:lmax+1)
    real(kind=8) :: clight, sign, ekappa, sk, yk, yk2, tane, ratfg, a1, b1, qint, x1, x2, r, cose, sine
    integer :: kap, l, lb, lact, i

    if (nonrel == 1) then
       clight = chuge
    else
       clight = csmall
    end if

    if (en > 0) then
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
          yk = real(l + 1, kind=8)
          lb = l + 1
          sk = -sign * ekappa
       end if

       yk2 = yk + yk

       call comdir(en, kap, v, nmt, nws, dx, x0, q(:, kap), p(:, kap), ratfg, nonrel)

       tane = (ratfg * fb(l) - sk * fb(lb)) / (ratfg * fn(l) - sk * fn(lb))
       ratx(kap) = ratfg
       tanx(kap) = tane

       ! Valence or continuum normalization
       if (nval == 1) then
          a1 = sk * ekappa * (fn(lb) - fb(lb) / tane) * rmt / q(nmt, kap) / clight
          b1 = ekappa * (fn(l) - fb(l) / tane) * rmt / p(nmt, kap)
       else
          eta(kap) = atan(tane)
310       cose = cos(eta(kap))
          sine = sin(eta(kap))
          a1 = sk * (cose * fb(lb) - sine * fn(lb)) * rmt / q(nmt, kap) / clight
          b1 = (cose * fb(l) - sine * fn(l)) * rmt / p(nmt, kap)
       end if

       do i = 1, nmt
          p(i, kap) = p(i, kap) * b1
          q(i, kap) = q(i, kap) * a1
          rr(i) = p(i, kap) * p(i, kap) + q(i, kap) * q(i, kap)
       end do

       ! Set up wavefunctions beyond muffin tin radius
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

       ! Normalize valence wavefunctions
       if (norm == 0) then
          vint(kap) = sintg(yk2, rr, rs, dx, nmt)
       else
          x1 = sintg(yk2, rr, rs, dx, nws - 1)
          x2 = sintg(yk2, rr, rs, dx, nws)
          vint(kap) = x1 + (x2 - x1) * (rws - rs(nws - 1)) / (rs(nws) - rs(nws - 1))
       end if
       qint = sqrt(1.0d0 / vint(kap))
       do i = 1, nrad
          p(i, kap) = p(i, kap) * qint
          q(i, kap) = q(i, kap) * qint
       end do
    end do

    if (nval == 1) then
       write(pun, 90)
90     format(/, 'CF/G RATIO', /)
       write(pun, 101) (kap, kap=-lact-1, -1), (kap, kap=1, lact)
       write(pun, 102) (ratx(kap), kap=-lact-1, -1), (ratx(kap), kap=1, lact)
       write(pun, 100)
100    format(/, 'TANGENT PHASESHIFTS', /)
       write(pun, 101) (kap, kap=-lact-1, -1), (kap, kap=1, lact)
101    format(2x, 12(5x, i3, 5x))
       write(pun, 102) (tanx(kap), kap=-lact-1, -1), (tanx(kap), kap=1, lact)
102    format(2x, 12d13.5)
    else
       write(pun, 90)
       write(pun, 101) (kap, kap=-1, -lact-1, -1)
       write(pun, 102) (ratx(kap), kap=-1, -lact-1, -1)
       write(pun, 104) (kap, kap=1, lact)
       write(pun, 105) (ratx(kap), kap=1, lact)
       write(pun, 103)
103    format(/, 'PHASESHIFTS', /)
       write(pun, 101) (kap, kap=-1, -lact-1, -1)
       write(pun, 102) (eta(kap), kap=-1, -lact-1, -1)
       write(pun, 104) (kap, kap=1, lact)
       write(pun, 105) (eta(kap), kap=1, lact)
104    format(15x, 11(6x, i2, 5x))
105    format(15x, 11d13.5)
    end if
  end subroutine wafu

  subroutine comdir(e1, kappa, za, nrc, nnk, dx, x0, q, p, ratfg, nonrel)
    implicit none
    real(kind=8), intent(in) :: e1, za(nrad), dx, x0
    integer, intent(in) :: kappa, nrc, nnk, nonrel
    real(kind=8), intent(out) :: q(nrad), p(nrad), ratfg
    real(kind=8) :: bgx(nrad), sxk(4), sxm(4), pp(nrad), qp(nrad)
    real(kind=8) :: test, c, cin, hoc, u, tc, t, xk, x, xc, bgc, wc, uc
    real(kind=8) :: unp, wnp, unp2, wnp2, z2, stval
    integer :: pun, jri, n, ik, nit, lkap, i

    test = 1.0d5
    pun = 99

    bgx(1:nnk) = za(1:nnk)
    xk = real(kappa, kind=8)
    jri = nrc
    e = e1
    stval = x0
    z2 = -bgx(1) * exp(stval)
    tc = exp(stval)

    if (nonrel == 1) then
       if (kappa < 0) then
          lkap = -kappa - 1
       else
          lkap = kappa
       end if
       kappa = -lkap - 1
       xk = real(kappa, kind=8)
       u = -0.5 * z2 / (lkap + 1.0)
       p(1) = 1.0d-20
       q(1) = u * 1.0d-20
    else
       c = csmall
       cin = 1.0d0 / (c * c)
       hoc = z2 / c
       if (abs(hoc / xk) <= 0.05) then
          u = (xk + abs(xk)) / hoc - 0.5 * hoc / abs(xk)
       else
          u = (xk + sqrt(xk * xk - hoc * hoc)) / hoc
       end if
       p(1) = 1.0d-20
       q(1) = c * u * 1.0d-20
    end if

    if (nonrel /= 1) then
       pp(1) = tc * (cin * (e - bgx(1)) + 1.0d0) * q(1) - xk * p(1)
    else
       pp(1) = tc * q(1) - xk * p(1)
    end if

    qp(1) = xk * q(1) - tc * (e - bgx(1)) * p(1)

    x = stval
    n = 1
25  ik = 0
    xc = x
    bgc = bgx(n)
    wc = q(n)
    uc = p(n)
20  ik = ik + 1
    t = exp(xc)
    if (nonrel /= 1) then
       sxk(ik) = dx * (-xk * uc + t * wc * (cin * (e - bgc) + 1.0d0))
    else
       sxk(ik) = dx * (-xk * uc + t * wc)
    end if

    sxm(ik) = dx * (xk * wc - t * (e - bgc) * uc)
    select case (ik)
    case (1)
       xc = xc + 0.5 * dx
       uc = uc + 0.5 * sxk(1)
       wc = wc + 0.5 * sxm(1)
       bgc = 0.5 * (bgc + bgx(n + 1))
       goto 20
    case (2)
       uc = uc + 0.5 * (sxk(2) - sxk(1))
       wc = wc + 0.5 * (sxm(2) - sxm(1))
       goto 20
    case (3)
       xc = xc + 0.5 * dx
       uc = uc + sxk(3) - 0.5 * sxk(2)
       wc = wc + sxm(3) - 0.5 * sxm(2)
       bgc = bgx(n + 1)
       goto 20
    case (4)
       q(n + 1) = q(n) + (sxm(1) + 2.0d0 * sxm(2) + 2.0d0 * sxm(3) + sxm(4)) / 6.0d0
       p(n + 1) = p(n) + (sxk(1) + 2.0d0 * sxk(2) + 2.0d0 * sxk(3) + sxk(4)) / 6.0d0

       if (nonrel /= 1) then
          pp(n + 1) = t * q(n + 1) * (cin * (e - bgc) + 1.0d0) - xk * p(n + 1)
       else
          pp(n + 1) = t * q(n + 1) - xk * p(n + 1)
       end if

       qp(n + 1) = xk * q(n + 1) - t * (e - bgc) * p(n + 1)
       x = x + dx
       n = n + 1
       if (n <= 6) goto 25
    end select

    x = x + dx
    t = exp(x)
27  unp = p(n - 5) + 0.3 * dx * (11.0d0 * pp(n) - 14.0d0 * pp(n - 1) + &
                                26.0d0 * pp(n - 2) - 14.0d0 * pp(n - 3) + &
                                11.0d0 * pp(n - 4))
    wnp = q(n - 5) + 0.3 * dx * (11.0d0 * qp(n) - 14.0d0 * qp(n - 1) + &
                                26.0d0 * qp(n - 2) - 14.0d0 * qp(n - 3) + &
                                11.0d0 * qp(n - 4))
    nit = 0

33  if (nonrel /= 1) then
       pp(n + 1) = t * (cin * (e - bgx(n + 1)) + 1.0d0) * wnp - xk * unp
    else
       pp(n + 1) = t * wnp - xk * unp
    end if

    qp(n + 1) = xk * wnp - t * (e - bgx(n + 1)) * unp
    unp2 = p(n - 3) + (7.0d0 * pp(n + 1) + 32.0d0 * pp(n) + &
                      12.0d0 * pp(n - 1) + 32.0d0 * pp(n - 2) + &
                      7.0d0 * pp(n - 3)) * 2.0d0 * dx / 45.0d0
    wnp2 = q(n - 3) + (7.0d0 * qp(n + 1) + 32.0d0 * qp(n) + &
                      12.0d0 * qp(n - 1) + 32.0d0 * qp(n - 2) + &
                      7.0d0 * qp(n - 3)) * 2.0d0 * dx / 45.0d0

    if (abs(test * (unp2 - unp)) <= abs(unp2)) then
       if (abs(test * (wnp2 - wnp)) <= abs(wnp2)) then
          q(n + 1) = wnp2
          p(n + 1) = unp2
          n = n + 1
          if (n <= nnk) goto 26
          ratfg = q(jri) / p(jri)
          return
       end if
    end if

    if (nit < 5) then
       nit = nit + 1
       wnp = wnp2
       unp = unp2
       goto 33
    end if

    q(n + 1) = wnp2
    p(n + 1) = unp2
    n = n + 1
    if (n <= nnk) goto 26
    ratfg = q(jri) / p(jri)
  end subroutine comdir

  subroutine cei(f2, f4, f1, f3, r, n, lambda, l2, l4, l1, l3, dx, rws, result)
    implicit none
    real(kind=8), intent(in) :: f1(nrad), f2(nrad), f3(nrad), f4(nrad), r(nrad)
    integer, intent(in) :: n, lambda
    real(kind=8), intent(in) :: l1, l2, l3, l4, dx, rws
    real(kind=8), intent(out) :: result
    real(kind=8) :: rl(ndim), rmlp1(ndim), arg1(ndim), arg2(ndim), arg3(ndim), arg4(ndim)
    real(kind=8) :: ll
    integer :: m, i

    m = n
    do i = 1, m
       rl(i) = r(i)**lambda
       rmlp1(i) = 1.0d0 / (rl(i) * r(i))
       arg1(i) = f2(i) * f4(i) * rl(i)
    end do

    ll = l2 + l4 + lambda
    call dintg(ll, 1, arg1, arg2, r, dx, m, rws)

    ll = ll + 1
    do i = 1, m
       arg2(i) = arg2(i) * f1(i) * f3(i) * rmlp1(i)
    end do

    ll = ll + l1 + l3 - lambda - 1
    do i = 1, m
       arg3(i) = f2(i) * f4(i) * rmlp1(i)
    end do

    ll = l2 + l4 - lambda - 1
    call dintg(ll, -1, arg3, arg4, r, dx, m, rws)

    ll = ll + 1
    do i = 1, m
       arg4(i) = arg4(i) * f1(i) * f3(i) * rl(i)
    end do

    ll = ll + l1 + l3 + lambda
    do i = 1, m
       arg2(i) = arg2(i) + arg4(i)
    end do

    result = sintg(ll, arg2, r, dx, m)
  end subroutine cei

  real(kind=8) function sintg(ll, fct, r, dx, n)
    implicit none
    real(kind=8), intent(in) :: ll, fct(nrad), r(nrad), dx
    integer, intent(in) :: n
    real(kind=8) :: sum, rll, fact, x0, x1, x2
    integer :: i

    if (ll + 1.0d0 <= tiny) then
       write(9, *) ' sintg: ll+1=', ll + 1
       sum = 0.0d0
    else
       rll = r(1)**ll
       if (rll <= tiny) then
          sum = 0.0d0
       else
          fact = fct(1) / rll
          sum = fact * (r(1)**(ll + 1)) / (ll + 1)
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
    implicit none
    real(kind=8), intent(in) :: ll, fct(nrad), r(nrad), dx, rws
    integer, intent(in) :: idir, n
    real(kind=8), intent(out) :: yint(nrad)
    real(kind=8) :: sum, rll, fact, x0, x1, corr
    integer :: i, j

    if (idir > 0) then
       if (ll + 1.0d0 <= tiny) then
          write(9, *) ' dintg: ll+1=', ll + 1
          sum = 0.0d0
       else
          rll = r(1)**ll
          if (rll <= tiny) then
             sum = 0.0d0
          else
             fact = fct(1) / rll
             sum = fact * (r(1)**(ll + 1)) / (ll + 1)
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
             do j = n - 3, n
                yint(j) = yint(j) - corr
             end do
          end if
          x1 = x0
       end do
    end if
  end subroutine dintg

  subroutine inter1(r, p, n, id, rs, ps)
    implicit none
    real(kind=8), intent(in) :: r(n), p(n), rs
    integer, intent(in) :: n, id
    real(kind=8), intent(out) :: ps
    real(kind=8) :: term, denom
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
    implicit none
    real(kind=8), intent(in) :: e, r
    real(kind=8), intent(out) :: fb(0:lmax+1), fn(0:lmax+1)
    real(kind=8) :: eroot, x1, x2, sinx, cosx, tl
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
          tl = real(2 * l - 1, kind=8)
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
          tl = real(2 * l - 1, kind=8)
          fb(l) = -tl * fb(l - 1) / x1 + fb(l - 2)
          fn(l) = -tl * fn(l - 1) / x1 + fn(l - 2)
       end do
    end if
  end subroutine sbf1

end program rccvmat

! External subroutines assumed to be available
interface
!  subroutine rec3jj(thrcof, l2, l3, m2, m3, l1min, l1max, lmatch, ndim, ier)
!    implicit none
!    real(kind=8), intent(out) :: thrcof(*), l1min, l1max, lmatch
!    real(kind=8), intent(in) :: l2, l3, m2, m3
!    integer, intent(in) :: ndim
!    integer, intent(out) :: ier
    SUBROUTINE rec3jj(Thrcof,L2,L3,M2,M3,L1min,L1max,Lmatch,Ndim,Ier)
! =====================================================================
!
!  j1-recursion of 3j-coefficients
! recursive evaluation of 3j- and
! 6j-coefficients.  k. schulten, r.g. gordon.
! ref. in comp. phys. commun. 11 (1976) 269
!
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a1 , a1s , a2 , a2s , c1 , c1old , c2 , cnorm , denom , dv , eps , half , huge , oldfac , one , ratio , sign1 , sign2 ,  &
        & srhuge , srtiny
   REAL*8 sum1 , sum2 , sumbac , sumfor , sumuni , Thrcof , three , thresh , tiny , two , x , x1 , x2 , x3 , y , y1 , y2 , y3 ,    &
        & zero
   INTEGER i , Ier , index , l1cmax , l1cmin , lstep , n , Ndim , nfin , nfinp1 , nfinp2 , nfinp3 , nlim , nstep2
!*** End of declarations inserted by SPAG
   REAL*8 l1 , L2 , L3 , m1 , M2 , M3 , L1min , L1max , newfac , Lmatch
   DIMENSION Thrcof(Ndim)
   INTEGER :: spag_nextblock_1
!
   DATA zero , eps , half , one/0D0 , 0.01D0 , 0.5D0 , 1D0/
   DATA two , three/2D0 , 3D0/
!
!  routine to generate set of 3j-coefficients   l1  l2  l3
!                                               m1  m2  m3
!  by recursion from l1min = max0(/l2-l3/,/m1/)
!                 to l1max =       l2+l3
!
!  the resulting 3j-coefficients are stored as thrcof(l1-l1min+1)
!
!  for a discussion of the recursion equation used see k. schulten
!  and r.g. gordon, j. math. phys. 16, 1961-1970 (1975), ibid. 16,
!  1971-1988 (1975)
!  for the sake of numerical stability the recursion will proceed
!  simultaneously forwards and backwards, starting from l1min and
!  l1max, respectively.
!  lmatch is the l1-value at which forward and backward recursion are
!  matched.
!  ndim is the length of the array thrcof.
!
!  ier is set to -1 if l2-/m2/ or l3-/m3/ less than zero or not
!                   integer, in which case all 3j-coefficients
!                   vanish.
!  ier is set to -2 if number of possible 3j-coefficients exceeds ndim
!  ier is set to non-negative number otherwise (0 + times of rescaling)
!
!  tiny should be set close to the smallest positive floating point
!  number which is representable on the computer.  srtiny is square
!  root of tiny .
!
   DATA tiny , srtiny/1.0D-10 , 1.0D-05/
!
!  huge should be set close to largest positive floating point
!  number which is representable on the computer.  srhuge is
!  square root of huge .
!
   DATA huge , srhuge/1.0D10 , 1.0D05/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         Lmatch = zero
         m1 = -M2 - M3
!
!  check relative magnitude of l- and m-values
         IF ( L2-dabs(M2)+eps>=0 ) THEN
            IF ( L3-dabs(M3)+eps>=0 ) THEN
               IF ( dmod(L2+dabs(M2)+eps,one)<eps+eps ) THEN
                  IF ( dmod(L3+dabs(M3)+eps,one)<eps+eps ) THEN
!
!
!
!  limits for l1
!
                     L1min = dmax1(dabs(L2-L3),dabs(m1))
                     L1max = L2 + L3
                     IF ( L1min<L1max-eps ) THEN
!
!
!
                        Ier = 0
                        nfin = idint(L1max-L1min+one+eps)
                        IF ( Ndim<nfin ) THEN
!
!  dimension of thrcof not large enough to hold all the coefficients
!  required
!
                           Ier = -2
                           WRITE (9,99001) L2 , L3 , m1 , M2 , M3 , nfin , Ndim
99001                      FORMAT (///1x,'3j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,'exceed storage provided  (',i4,',',i4,')')
                           RETURN
                        ELSE
!
!
!  starting forward recursion from l1min taking nstep1 steps
!
                           l1 = L1min
                           Thrcof(1) = srtiny
                           sum1 = (l1+l1+one)*tiny
!
!
                           lstep = 1
                           spag_nextblock_1 = 2
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                     ELSEIF ( L1min<L1max+eps ) THEN
!
!
!  this is reached in case that l1 can take only one value,
!  i.e. l1min = l1max
!
                        Ier = 0
                        Thrcof(1) = (-one)**idint(dabs(L2+M2-L3+M3)+eps)/dsqrt(L1min+L2+L3+one)
                        l1cmin = L1min
                        l1cmax = L1max
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
!
!  this is reached if l2-/m2/ and l3-/m3/  less than zero or not integer
!
         Ier = -1
         WRITE (9,99002) L2 , L3 , m1 , M2 , M3
99002    FORMAT (///1x,'3j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,'do not satisfy the condition l2-/m2/ and l3-/m3/ ge zero ',   &
                &'and integer')
         RETURN
      CASE (2)
         lstep = lstep + 1
         l1 = l1 + one
!
!
         oldfac = newfac
         a1 = (l1+L2+L3+one)*(l1-L2+L3)*(l1+L2-L3)*(-l1+L2+L3+one)
         a2 = (l1+m1)*(l1-m1)
         newfac = dsqrt(a1*a2)
         IF ( l1<one+eps ) THEN
!
!  if l1 = 1  (l1-1) has to be factored out of dv, hence
!
            c1 = -(l1+l1-one)*l1*(M3-M2)/newfac
         ELSE
!
!
            dv = -L2*(L2+one)*m1 + L3*(L3+one)*m1 + l1*(l1-one)*(M3-M2)
            denom = (l1-one)*newfac
!
!
            IF ( lstep>2 ) c1old = dabs(c1)
            c1 = -(l1+l1-one)*dv/denom
         ENDIF
!
         IF ( lstep>2 ) THEN
!
!
            c2 = -l1*oldfac/denom
!
!  recursion to the next 3j-coefficient x
!
            x = c1*Thrcof(lstep-1) + c2*Thrcof(lstep-2)
            Thrcof(lstep) = x
            sumfor = sum1
            sum1 = sum1 + (l1+l1+one)*x*x
            IF ( lstep/=nfin ) THEN
!
!  see if last unnormalized 3j-coefficient exceeds srhuge
!
               IF ( dabs(x)>=srhuge ) THEN
!
!  this is reached if last 3j-coefficient larger than srhuge
!  so that the recursion series thrcof(1), ... , thrcof(lstep)
!  has to be rescaled to prevent overflow
!
                  Ier = Ier + 1
                  DO i = 1 , lstep
                     IF ( dabs(Thrcof(i))<srtiny ) Thrcof(i) = zero
                     Thrcof(i) = Thrcof(i)/srhuge
                  ENDDO
                  sum1 = sum1/huge
                  sumfor = sumfor/huge
                  x = x/srhuge
               ENDIF
!
!  as long as /c1/ is decreasing the recursion proceeds towards
!  increasing 3j-values and, hence, is numerically stable.  once
!  an increase of /c1/ is detected the recursion direction is
!  reversed.
!
               IF ( c1old>dabs(c1) ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDIF
!
!
!  keep three 3j-coefficients around lmatch for comparision with
!  backward recursion.
!
            Lmatch = l1 - 1
            x1 = x
            x2 = Thrcof(lstep-1)
            x3 = Thrcof(lstep-2)
            nstep2 = nfin - lstep + 3
!
!
!
!
!  starting backward recursion from l1max taking nstep2 steps, so
!  that forward and backward recursion overlap at three points
!  l1 = lmatch+1, lmatch, lmatch-1.
!
            nfinp1 = nfin + 1
            nfinp2 = nfin + 2
            nfinp3 = nfin + 3
            l1 = L1max
            Thrcof(nfin) = srtiny
            sum2 = tiny*(l1+l1+one)
!
            l1 = l1 + two
            lstep = 1
            DO
               lstep = lstep + 1
               l1 = l1 - one
!
               oldfac = newfac
               a1s = (l1+L2+L3)*(l1-L2+L3-one)*(l1+L2-L3-one)*(-l1+L2+L3+two)
               a2s = (l1+m1-one)*(l1-m1-one)
               newfac = dsqrt(a1s*a2s)
!
               dv = -L2*(L2+one)*m1 + L3*(L3+one)*m1 + l1*(l1-one)*(M3-M2)
!
               denom = l1*newfac
               c1 = -(l1+l1-one)*dv/denom
               IF ( lstep>2 ) THEN
!
!
                  c2 = -(l1-one)*oldfac/denom
!
!  recursion to the next 3j-coefficient y
!
                  y = c1*Thrcof(nfinp2-lstep) + c2*Thrcof(nfinp3-lstep)
!
                  IF ( lstep==nstep2 ) THEN
!
!
!  the forward recursion 3j-coefficients x1, x2, x3 are to be matched
!  with the corresponding backward recursion values y1, y2, y3.
!
                     y3 = y
                     y2 = Thrcof(nfinp2-lstep)
                     y1 = Thrcof(nfinp3-lstep)
!
!
!  determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
!  with minimal error.
!
                     ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
                     nlim = nfin - nstep2 + 1
!
                     IF ( dabs(ratio)<one ) THEN
!
                        nlim = nlim + 1
                        ratio = one/ratio
                        DO n = nlim , nfin
                           Thrcof(n) = ratio*Thrcof(n)
                        ENDDO
                        sumuni = sumfor + ratio*ratio*sumbac
                     ELSE
!
                        DO n = 1 , nlim
                           Thrcof(n) = ratio*Thrcof(n)
                        ENDDO
                        sumuni = ratio*ratio*sumfor + sumbac
                     ENDIF
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ELSE
!
                     Thrcof(nfinp1-lstep) = y
                     sumbac = sum2
                     sum2 = sum2 + (l1+l1-three)*y*y
!
!  see if last unnormalized 3j-coefficient exceeds srhuge
!
                     IF ( dabs(y)>=srhuge ) THEN
!
!  this is reached if last 3j-coefficient larger than srhuge
!  so that the recursion series thrcof(nfin), ... ,thrcof(nfin-lstep+1)
!  has to be rescaled to prevent overflow
!
                        Ier = Ier + 1
                        DO i = 1 , lstep
                           index = nfin - i + 1
                           IF ( dabs(Thrcof(index))<srtiny ) Thrcof(index) = zero
                           Thrcof(index) = Thrcof(index)/srhuge
                        ENDDO
                        sum2 = sum2/huge
!
!
                        sumbac = sumbac/huge
                     ENDIF
                  ENDIF
               ELSE
!
!  if l1 = l1max + 1  the third term in the recursion formula vanishes
!
                  y = srtiny*c1
                  Thrcof(nfin-1) = y
                  sumbac = sum2
!
                  sum2 = sum2 + tiny*(l1+l1-three)*c1*c1
               ENDIF
            ENDDO
         ELSE
!
!
!  if l1 = l1min + 1  the third term in the recursion equation vanishes
!  , hence
            x = srtiny*c1
            Thrcof(2) = x
            sum1 = sum1 + tiny*(l1+l1+one)*c1*c1
            IF ( lstep/=nfin ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!
            sumuni = sum1
         ENDIF
         spag_nextblock_1 = 3
      CASE (3)
!
!
!  normalize 3j-coefficients
!
         cnorm = one/dsqrt(sumuni)
!
!  sign convention for last 3j-coefficient determines overall phase
!
         sign1 = dsign(one,Thrcof(nfin))
         sign2 = (-one)**idint(dabs(L2+M2-L3+M3)+eps)
         IF ( sign1*sign2<=0 ) cnorm = -cnorm
!
         IF ( dabs(cnorm)<one ) THEN
!
            thresh = tiny/dabs(cnorm)
            DO n = 1 , nfin
               IF ( dabs(Thrcof(n))<thresh ) Thrcof(n) = zero
               Thrcof(n) = cnorm*Thrcof(n)
            ENDDO
            EXIT SPAG_DispatchLoop_1
         ENDIF
!
         DO n = 1 , nfin
            Thrcof(n) = cnorm*Thrcof(n)
         ENDDO
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
  end subroutine rec3jj

!  subroutine rec6j(sixcof, l2, l3, l4, l5, l6, l1min, l1max, lmatch, ndim, ier)
!    implicit none
!    real(kind=8), intent(out) :: sixcof(*), l1min, l1max, lmatch
!    real(kind=8), intent(in) :: l2, l3, l4, l5, l6
!    integer, intent(in) :: ndim
!    integer, intent(out) :: ier
!  end subroutine rec6j
  SUBROUTINE rec6j(Sixcof,L2,L3,L4,L5,L6,L1min,L1max,Lmatch,Ndim,Ier)
! ===================================================================
!
!  j1-recursion of 6j-coefficients
! recursive evaluation of 3j- and
! 6j-coefficients.  k. schulten, r.g. gordon.
! ref. in comp. phys. commun. 11 (1976) 269
!
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a1 , a1s , a2 , a2s , c1 , c1old , c2 , cnorm , denom , dv , eps , half , huge , oldfac , one , ratio , sign1 , sign2 ,  &
        & Sixcof , srhuge
   REAL*8 srtiny , sum1 , sum2 , sumbac , sumfor , sumuni , three , thresh , tiny , two , x , x1 , x2 , x3 , y , y1 , y2 , y3 ,    &
        & zero
   INTEGER i , Ier , index , l1cmax , l1cmin , lstep , n , Ndim , nfin , nfinp1 , nfinp2 , nfinp3 , nlim , nstep2
!*** End of declarations inserted by SPAG
   REAL*8 l1 , L2 , L3 , L4 , L5 , L6 , L1min , L1max , newfac , Lmatch
   DIMENSION Sixcof(Ndim)
   INTEGER :: spag_nextblock_1
!
   DATA zero , eps , half , one/0D0 , 0.01D0 , 0.5D0 , 1D0/
   DATA two , three/2D0 , 3D0/
!
!  routine to generate the set of 6j-coefficients    l1 l2 l3
!                                                    l4 l5 l6
!  by recursion from l1min = max0(/l2-l3/,/l5-l6/)
!                 to l1max = min0( l2+l3 , l5+l6 )
!
!  the resulting 6j-coefficients are stored as sixcof(l1-l1min+1).
!
!  for a discussion of the recursion equation used see k. schulten
!  and r.g. gordon, j. math. phys. 16, 1961-1970 (1975), ibid. 16,
!  1971-1988 (1975)
!  for the sake of numerical stability the recursion will proceed
!  simultaneously forward and backwards, starting at l1min and
!  l1max, respectively.
!  lmatch is the l1-value at which forward and backward recursion
!  are matched.  its value will be returned, though it is not of
!  actual use.
!
!  ndim is the length of the array sixcof.
!
!  ier is set to -1 if either no 6j-coefficient satisfies triangular
!                  condition or if l2+l3+l5+l6 or l2+l4+l6 not
!                  integer, in which case all 6j-coefficients
!                  vanish.
!  ier is set to -2 if number of possible 6j-coefficients exceeds ndim
!
!  tiny should be set close to the smallest positive floating point
!  number which is representable on the computer.  srtiny is square
!  root of tiny .
!
   DATA tiny , srtiny/1.0D-10 , 1.0D-05/
!
!  huge should be set close to largest positive floating point
!  number which is representable on the computer.  srhuge is
!  square root of huge .
   DATA huge , srhuge/1.0D10 , 1.0D05/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         Lmatch = zero
!
!
!
!  check if 6j-coefficients obey selection rules
!
         IF ( dmod(L2+L3+L5+L6+eps,one)<eps+eps ) THEN
            IF ( dmod(L4+L2+L6+eps,one)<eps+eps ) THEN
!
               IF ( L4+L2>=L6 ) THEN
                  IF ( L4-L2+L6>=0 ) THEN
                     IF ( -L4+L2+L6>=0 ) THEN
!
                        IF ( L4+L3>=L5 ) THEN
                           IF ( L4-L3+L5>=0 ) THEN
                              IF ( -L4+L3+L5>=0 ) THEN
!
!  limits for l1
!
                                 L1min = dmax1(dabs(L2-L3),dabs(L5-L6))
                                 L1max = dmin1(L2+L3,L5+L6)
                                 IF ( L1min<L1max-eps ) THEN
!
!
                                    Ier = 0
                                    nfin = idint(L1max-L1min+one+eps)
                                    IF ( Ndim<nfin ) THEN
!
!  this is reached if array sixcof not large enough to hold all
!  6j - coefficients required
!
                                       Ier = -2
                                       WRITE (9,99001) L2 , L3 , L4 , L5 , L6 , nfin , Ndim
99001                                  FORMAT (///1x,'6j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,'exceed storage provided  (',i4, &
                                         &',',i4,')')
                                       RETURN
                                    ELSE
!
!
!
!
!
!  start of forward recursion
!
!
                                       l1 = L1min
                                       Sixcof(1) = srtiny
                                       sum1 = (l1+l1+one)*tiny
!
                                       lstep = 1
                                       spag_nextblock_1 = 2
                                       CYCLE SPAG_DispatchLoop_1
                                    ENDIF
                                 ELSEIF ( L1min<L1max+eps ) THEN
!
!
!  this is reached in case that l1 can take only one value
!
                                    Ier = 0
                                    Sixcof(1) = (-one)**idint(L2+L3+L5+L6+eps)/dsqrt((L1min+L1min+one)*(L4+L4+one))
                                    l1cmin = L1min
                                    l1cmax = L1max
                                    RETURN
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
!
!  this is reached if triangular condition not satisfied
!  or if l2+l3+l5+l6 or l2+l4+l6 not integer
!
         Ier = -1
         WRITE (9,99002) L2 , L3 , L4 , L5 , L6
99002    FORMAT (///1x,'6j-coefficients',9x,'l1',2F7.1,4x,'do not satisfy triangular conditions or'/20x,3F7.1,4x,                  &
                &'l2+l3+l5+l6 or l2+l4+l6 not integer')
         RETURN
      CASE (2)
         lstep = lstep + 1
         l1 = l1 + one
!
         oldfac = newfac
         a1 = (l1+L2+L3+one)*(l1-L2+L3)*(l1+L2-L3)*(-l1+L2+L3+one)
         a2 = (l1+L5+L6+one)*(l1-L5+L6)*(l1+L5-L6)*(-l1+L5+L6+one)
         newfac = dsqrt(a1*a2)
!
         IF ( l1<one+eps ) THEN
!
!  if l1 = 1   (l1 - 1) has to be factored out of dv, hence
!
            c1 = -two*(L2*(L2+one)+L5*(L5+one)-L4*(L4+one))/newfac
         ELSE
!
            dv = two*(L2*(L2+one)*L5*(L5+one)+L3*(L3+one)*L6*(L6+one)-l1*(l1-one)*L4*(L4+one))                                     &
               & - (L2*(L2+one)+L3*(L3+one)-l1*(l1-one))*(L5*(L5+one)+L6*(L6+one)-l1*(l1-one))
!
            denom = (l1-one)*newfac
!
!
            IF ( lstep>2 ) c1old = dabs(c1)
            c1 = -(l1+l1-one)*dv/denom
         ENDIF
!
         IF ( lstep>2 ) THEN
!
!
            c2 = -l1*oldfac/denom
!
!  recursion to the next 6j - coefficient x
!
            x = c1*Sixcof(lstep-1) + c2*Sixcof(lstep-2)
            Sixcof(lstep) = x
!
            sumfor = sum1
            sum1 = sum1 + (l1+l1+one)*x*x
            IF ( lstep/=nfin ) THEN
!
!  see if last unnormalized 6j-coefficient exceeds srhuge
!
               IF ( dabs(x)>=srhuge ) THEN
!
!  this is reached if last 6j-coefficient larger than srhuge
!  so that the recursion series sixcof(1), ... ,sixcof(lstep)
!  has to be rescaled to prevent overflow
!
                  Ier = Ier + 1
                  DO i = 1 , lstep
                     IF ( dabs(Sixcof(i))<srtiny ) Sixcof(i) = zero
                     Sixcof(i) = Sixcof(i)/srhuge
                  ENDDO
                  sum1 = sum1/huge
                  sumfor = sumfor/huge
                  x = x/srhuge
               ENDIF
!
!
!  as long as the coefficient /c1/ is decreasing the recursion proceeds
!  towards increasing 6j-values and, hence, is numerically stable.
!  once an increase of /c1/ is detected, the recursion direction is
!  reversed.
!
               IF ( c1old>dabs(c1) ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDIF
!
!
!  keep three 6j-coefficients around lmatch for comparision later
!  with backward recursion.
!
            Lmatch = l1 - 1
            x1 = x
            x2 = Sixcof(lstep-1)
            x3 = Sixcof(lstep-2)
!
!
!
!  starting backward recursion from l1max taking nstep2 steps, so
!  that forward and backward recursion overlap at the three points
!  l1 = lmatch+1, lmatch, lmatch-1.
!
            nfinp1 = nfin + 1
            nfinp2 = nfin + 2
            nfinp3 = nfin + 3
            nstep2 = nfin - lstep + 3
            l1 = L1max
!
            Sixcof(nfin) = srtiny
            sum2 = (l1+l1+one)*tiny
!
!
            l1 = l1 + two
            lstep = 1
            SPAG_Loop_2_1: DO
               lstep = lstep + 1
               l1 = l1 - one
!
               oldfac = newfac
               a1s = (l1+L2+L3)*(l1-L2+L3-one)*(l1+L2-L3-one)*(-l1+L2+L3+two)
               a2s = (l1+L5+L6)*(l1-L5+L6-one)*(l1+L5-L6-one)*(-l1+L5+L6+two)
               newfac = dsqrt(a1s*a2s)
!
               dv = two*(L2*(L2+one)*L5*(L5+one)+L3*(L3+one)*L6*(L6+one)-l1*(l1-one)*L4*(L4+one))                                  &
                  & - (L2*(L2+one)+L3*(L3+one)-l1*(l1-one))*(L5*(L5+one)+L6*(L6+one)-l1*(l1-one))
!
               denom = l1*newfac
               c1 = -(l1+l1-one)*dv/denom
               IF ( lstep>2 ) THEN
!
!
                  c2 = -(l1-one)*oldfac/denom
!
!  recursion to the next 6j - coefficient y
!
                  y = c1*Sixcof(nfinp2-lstep) + c2*Sixcof(nfinp3-lstep)
                  IF ( lstep==nstep2 ) EXIT SPAG_Loop_2_1
                  Sixcof(nfinp1-lstep) = y
                  sumbac = sum2
                  sum2 = sum2 + (l1+l1-three)*y*y
!
!  see if last unnormalized 6j-coefficient exceeds srhuge
!
                  IF ( dabs(y)>=srhuge ) THEN
!
!  this is reached if last 6j-coefficient larger than srhuge
!  so that the recursion series sixcof(nfin), ... ,sixcof(nfin-lstep+1)
!  has to be rescaled to prevent overflow
!
                     Ier = Ier + 1
                     DO i = 1 , lstep
                        index = nfin - i + 1
                        IF ( dabs(Sixcof(index))<srtiny ) Sixcof(index) = zero
                        Sixcof(index) = Sixcof(index)/srhuge
                     ENDDO
                     sumbac = sumbac/huge
!
                     sum2 = sum2/huge
                  ENDIF
               ELSE
!
!  if l1 = l1max + 1 the third term in the recursion equation vanishes
!
                  y = srtiny*c1
                  Sixcof(nfin-1) = y
                  IF ( lstep==nstep2 ) EXIT SPAG_Loop_2_1
                  sumbac = sum2
                  sum2 = sum2 + (l1+l1-three)*c1*c1*tiny
               ENDIF
            ENDDO SPAG_Loop_2_1
!
!
!  the forward recursion 6j-coefficients x1, x2, x3 are to be matched
!  with the corresponding backward recursion values y1, y2, y3.
!
            y3 = y
            y2 = Sixcof(nfinp2-lstep)
            y1 = Sixcof(nfinp3-lstep)
!
!
!  determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
!  with minimal error.
!
            ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
            nlim = nfin - nstep2 + 1
!
            IF ( dabs(ratio)<one ) THEN
!
               nlim = nlim + 1
               ratio = one/ratio
               DO n = nlim , nfin
                  Sixcof(n) = ratio*Sixcof(n)
               ENDDO
               sumuni = sumfor + ratio*ratio*sumbac
            ELSE
!
               DO n = 1 , nlim
                  Sixcof(n) = ratio*Sixcof(n)
               ENDDO
               sumuni = ratio*ratio*sumfor + sumbac
            ENDIF
         ELSE
!
!  if l1 = l1min + 1 the third term in recursion equation vanishes
!
            x = srtiny*c1
            Sixcof(2) = x
            sum1 = sum1 + tiny*(l1+l1+one)*c1*c1
!
            IF ( lstep/=nfin ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!
            sumuni = sum1
         ENDIF
!
!
!  normalize 6j-coefficients
!
         cnorm = one/dsqrt((L4+L4+one)*sumuni)
!
!  sign convention for last 6j-coeff. determines overall phase
!
         sign1 = dsign(one,Sixcof(nfin))
         sign2 = (-one)**idint(L2+L3+L5+L6+eps)
         IF ( sign1*sign2<=0 ) cnorm = -cnorm
!
         IF ( dabs(cnorm)<one ) THEN
!
            thresh = tiny/dabs(cnorm)
            DO n = 1 , nfin
               IF ( dabs(Sixcof(n))<thresh ) Sixcof(n) = zero
               Sixcof(n) = cnorm*Sixcof(n)
            ENDDO
            EXIT SPAG_DispatchLoop_1
         ENDIF
!
         DO n = 1 , nfin
            Sixcof(n) = cnorm*Sixcof(n)
         ENDDO
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
END SUBROUTINE rec6j

end interface
