module constants
   implicit none
   integer, parameter :: noat = 2
   integer, parameter :: nrad = 400
   integer, parameter :: iorb = 30
   integer, parameter :: irep = 7
   real(8), parameter :: pi = 4.0_8 * atan(1.0_8)
   real(8), parameter :: pi4 = 4.0_8 * pi
   real(8), parameter :: third = 1.0_8 / 3.0_8
   real(8), parameter :: twoth = 2.0_8 * third
   real(8), parameter :: tofp = 3.0_8 / pi4
end module constants

module common_data
   use constants
   implicit none
   real(8) :: den(iorb), dq1(iorb), dfl(iorb), dv(nrad), dr(nrad), &
              dp(nrad), dq(nrad), dpas, z, tets, test, testv, &
              dgc(nrad,iorb), dpc(nrad,iorb), xx(nrad), yy(nrad), &
              ra(noat,nrad), za(noat,nrad), stval, delx, &
              rhotot(nrad,noat), char(nrad,noat), rhoc(nrad,noat), &
              rhoold(nrad,noat), thfpot(nrad), rws, rws1, &
              vcc(noat), vew(noat), tcore(noat), vold(nrad), &
              d(nrad), p(nrad), r(nrad)
   integer :: nqn(iorb), nql(iorb), nk(iorb), nmax(iorb), nel(iorb), &
              norb, icore, nstop, nes, np, nuc, icoul, &
              nt(noat), nz(noat), nnk(noat), nrc(noat), nrws(noat), &
              nrws1(noat), itot, iex, iscfat, nitpot, itatom, &
              iskip, iunit8, iunt14, iturn, itcore, jturn, itval, &
              iemode, icone
   character(4) :: titre(iorb), bar(10)
   real(8) :: rc(noat), vc(noat)
   complex(8) :: sctamp(irep)
   real(8) :: conc, confac, aa, conc1, avol, ec
   integer :: lattyp
   character(80) :: title(20)
   real(8) :: dep(10), deq(10), db, dvc, dsal, dk, dm
   real(8) :: dpno(4,iorb), dqno(4,iorb)
end module common_data

program rcore
   use common_data
   use constants
   implicit none
   real(8) :: cwf(iorb), wa, qval, tail, rlast, anorm, bnorm, eone, ev, evcor
   integer :: i, j, k, natom, ipot, moat, last, last1, imax0, kk, ielec, nitoe, iend
   character(80) :: aa

   write(*, *) ' input?'
   read(*, '(a80)') aa
   open(unit=3, file=trim(aa), status='old', action='read')

   read(3, '(a80)') aa
   open(unit=8, file=trim(aa), status='old', action='read')
   read(3, '(a80)') aa
   open(unit=6, file=trim(aa), status='replace', action='write')
   read(3, '(a80)') aa
   open(unit=7, file=trim(aa), status='replace', action='write')

   read(3, *) ipot, iskip, moat

   if (iunit8 == 1) then
      natom = moat
      call insld(natom)
   else
      select case (ipot)
         case (0)
            call cpain
         case (1)
            call apwin
         case (2)
            call kkrin
      end select
      natom = moat
      call insld(natom)
   end if

   wa = -1.0_8
   z = nz(natom)
   itatom = 0
   if (iscfat == 0) norb = icore

   do i = 1, nrad
      rr = exp(stval + (i-1)*delx)
      rhotot(i,natom) = 0.0_8
      rhoc(i,natom) = 0.0_8
      dr(i) = rr
      thfpot(i) = 0.0_8
      if (ipot < 0) dv(i) = fpot(rr, z, wa)
   end do

   if (ipot >= 0) then
      dp(1:nrad) = rhoold(1:nrad,natom)

      if (iemode == 0) then
         last = nrc(natom)
         rlast = rc(natom)
         qval = rsimp(dp, dr, rc(natom), nrc(natom))
      else if (iemode == 1) then
         last = nrws1(natom)
         rlast = rws1
         qval = rsimp3(dp, dr, rws1, nrws1(natom))
      else
         last = nrws(natom)
         rlast = rws
         qval = rsimp3(dp, dr, rws, nrws(natom))
      end if

      if (icoul == 2) qval = 0.0_8
      if (icoul /= 0) then
         thfpot(1:nrad) = -(nz(natom) - qval) / dr(1:nrad)
      end if

      tail = 0.0_8
      if (icoul /= 0) tail = (nz(natom) - qval) / rlast
      last1 = last + 1

      do i = 1, last
         dv(i) = za(natom,i) / 2.0_8 - vc(natom) / 2.0_8 - vcc(natom) / 2.0_8
      end do
      do i = last1, nrad
         dv(i) = thfpot(i) + tail
      end do

      write(6, '( " Potential for type", i2)') natom
      write(6, '(4d20.10)') (dr(k), dv(k), k=1,nrad)
   end if

   do i = 1, nrad
      xx(i) = 0.0_8
      do j = 1, norb
         dgc(i,j) = 0.0_8
         dpc(i,j) = 0.0_8
      end do
   end do

   imax0 = 0
   np = nrad

   do j = 1, norb
      tets = test
      do
         dp(1:nrad) = 0.0_8
         dq(1:nrad) = 0.0_8
         call resld(nqn(j), nql(j), nk(j), imax, den(j), dfl(j), dq1(j), j, natom)
         if (imax > imax0) imax0 = imax
         if (nstop == 0) exit
         tets = test * 10.0_8
      end do

      dp(1) = 0.0_8
      dq(1) = 0.0_8
      do i = 1, imax
         dgc(i,j) = dp(i)
         dpc(i,j) = dq(i)
         xx(i) = xx(i) + nel(j) * (dp(i)**2 + dq(i)**2)
         if (ipot < 0) rhotot(i,natom) = xx(i)
         if (j <= icore) rhoc(i,natom) = xx(i)
      end do
   end do

   np = imax0
   if (iscfat /= 0) then
      call poisat(natom)
      if (nstop == 0) cycle
   end if

   rinf = dr(np-1)
   if (iskip /= 0) then
      write(6, '(1x, 10x, "one electron energies"//)')
      do i = 1, norb
         write(6, '(10x, i2, a4, 5x, e15.8)') nqn(i), titre(i), den(i)
      end do
      write(6, '(10x//10x, "orthogonality relations"//)')
   end if

   do i = 1, norb
      do j = i, norb
         if (nql(i) /= nql(j) .or. nk(i) /= nk(j)) cycle
         yy(1:np) = dpc(1:np,i) * dpc(1:np,j) + dgc(1:np,i) * dgc(1:np,j)
         if (iemode == 0) then
            anorm = rsimp(yy, dr, rc(natom), nrc(natom))
         else if (iemode == 1) then
            anorm = rsimp3(yy, dr, rws1, nrws1(natom))
         else
            anorm = rsimp3(yy, dr, rws, nrws(natom))
         end if
         if (i == j) cwf(i) = anorm
         if (iskip /= 0) then
            bnorm = rsimp(yy, dr, rinf, np)
            write(6, '(10x, "(", i2, a4, ",", i2, a4, ")", 3x, 2e15.8)') &
               nqn(i), titre(i), nqn(j), titre(j), bnorm, anorm
         end if
      end do
   end do

   xx(1:np) = rhoc(1:np,natom)
   if (iemode == 0) then
      anorm = rsimp(xx, dr, rc(natom), nrc(natom))
   else if (iemode == 1) then
      anorm = rsimp3(xx, dr, rws1, nrws1(natom))
   else
      anorm = rsimp3(xx, dr, rws, nrws(natom))
   end if
   if (iturn >= itcore .and. jturn >= itval) then
      write(6, '(10x//10x, "number of core-electrons in atomic sphere = ", e15.8)') anorm
   end if

   if (iscfat /= 0) then
      call totec(natom)
      cycle
   end if

   eone = 0.0_8
   do i = 1, icore
      anorm = merge(cwf(i), 1.0_8, icone == 1)
      eone = eone + anorm * nel(i) * den(i)
   end do
   eone = 2.0_8 * eone

   dq(1) = 0.0_8
   nkk1 = nnk(natom)
   do k = 1, nkk1
      dp(k) = (za(natom,k) - vc(natom) - vcc(natom)) * rhoc(k,natom)
   end do
   if (iemode == 0) then
      ev = rsimp(dp, dr, rc(natom), nrc(natom))
   else
      ev = rsimp3(dp, dr, rws1, nrws1(natom))
   end if

   evcor = 0.0_8
   if (icoul /= 0) then
      do i = last1, nrad
         if (rhoc(i,natom) == 0.0_8) exit
      end do
      iend = i - 1
      do i = last1, iend
         evcor = evcor + (dv(i-1) + dv(i)) * (dr(i) - dr(i-1))
      end do
   end if

   tcore(natom) = eone - ev - evcor
   if (iskip /= 0) then
      write(6, '(15x, "e-one core ", e20.12/15x, "pot -term ", e20.12/15x, &
                "<V-correction> ", e20.12/15x, "<e-c> ", e20.12/)') &
         eone, ev, evcor, tcore(natom)
   end if

   write(7, '(4d20.12)') (rhoc(j,natom), j=1,nnk(natom))
   do i = 1, norb
      write(7, '(i1, a4)') nqn(i), titre(i)
      write(7, '(e15.8)') den(i)
      write(7, '(5i4)') nk(i)
      write(7, '(5i4)') nnk(natom)
      do j = 1, nnk(natom)
         write(7, '(3e20.10)') dr(j), dgc(j,i), dpc(j,i)
      end do
   end do
end do

close(3)
close(6)
close(7)
close(8)
end program rcore

subroutine apwin
   use common_data
   use constants
   implicit none
   integer :: n, k, kx, k1, kmax, i, j, u, ilat
   real(8) :: vbar

   iunit8 = 1
   read(8, '(20a4)') title
   read(8, *) aa0, conc, ilat
   read(8, *) stval, xdel

   if (stval == 0.0_8) stval = -8.8_8
   if (xdel == 0.0_8) xdel = 0.05_8

   conc1 = conc
   u = 0
   if (iskip /= 0) then
      write(6, '(1x, 20x, 20a4//)') title
      write(6, '(10x, "a0 = ", f10.5, " Charge = ", f10.5/)') aa0, conc
   end if

   do n = 1, noat
      read(8, *) nz(n), nt(n), rc(n), vc(n)
      if (iskip /= 0) write(6, '(10x, "z = ", i4, " nk = ", i4, &
          " rc = ", f10.5, " vc = ", f10.5/)') nz(n), nt(n), rc(n), vc(n)
      kx = nt(n)
      read(8, '(4e20.8)') (ra(n,i), za(n,i), i=1,kx)
   end do

   call fitrad(noat)
   do n = 1, noat
      k1 = nrc(n)
      rc(n) = ra(n,k1)
   end do

   call mike(ilat)
   do n = 1, noat
      kmax = nt(n)
      do k = 1, kmax
         if (ra(n,k) > rws) exit
      end do
      nrws(n) = k - 1
      nt(n) = k
      do k = 1, kmax
         if (ra(n,k) > rws1) exit
      end do
      nrws1(n) = k - 1
   end do

   vbar = conc * vc(1) + (1.0_8 - conc) * vc(2)
   vc(1:2) = vbar

   if (iskip /= 0) then
      write(6, '(1x, 10x/10x///"-----new input"///)')
      do n = 1, noat
         kmax = nt(n)
         write(6, '(72x/12x, 4("r", 14x, "v", 14x))')
         write(6, '(1x, i3, 8e15.7)') &
            (i, (ra(n,i+j), za(n,i+j), j=u,3), i=1,kmax,4)
      end do
   end if
end subroutine apwin

subroutine cpain
   use common_data
   use constants
   implicit none
   integer :: n, i, kmax, k
   real(8) :: vbar

   iunit8 = 1
   read(8, '(20a4)') title
   read(8, *) aa0, conc, ilat
   read(8, *) stval, xdel

   if (stval == 0.0_8) stval = -8.8_8
   if (xdel == 0.0_8) xdel = 0.05_8

   conc1 = conc
   do n = 1, noat
      nt(n) = n
      read(8, *) nz(n), nt(n), nrc(n), rc(n), vc(n)
      kmax = nt(n)
      read(8, '(5d15.8)') (za(n,i), i=1,kmax)
      do i = 1, kmax
         ra(n,i) = exp(stval + xdel * (i-1))
         za(n,i) = za(n,i) / ra(n,i)
      end do
   end do

   call mike(ilat)
   do n = 1, noat
      kmax = nt(n)
      do k = 1, kmax
         if (ra(n,k) > rws) exit
      end do
      nrws(n) = k - 1
      nt(n) = k
      do k = 1, kmax
         if (ra(n,k) > rws1) exit
      end do
      nrws1(n) = k - 1
   end do

   vbar = conc * vc(1) + (1.0_8 - conc) * vc(2)
   vc(1:2) = vbar

   if (iskip >= 2) then
      write(6, '(1x, 10x, 20a4//)') title
      write(6, '(10x, "lattice constant ", f10.5/10x, "concentration ", f10.5/&
                10x, "muffin tin zero ", f10.5//)') aa0, conc, vc(1)
      write(6, '(10x, "conversion factor to rydberg ", f10.5/)') confac
      do n = 1, noat
         write(6, '(1x, 10x, "potential for scatterer ", i4//)') n
         write(6, '(72x/12x, 4("r", 14x, "v", 14x))')
         kmax = nt(n)
         write(6, '(1x, i3, 8e15.7)') &
            (i, (ra(n,i+j), za(n,i+j), j=0,3), i=1,kmax,4)
      end do
   end if
end subroutine cpain

function fpot(r, z, wa) result(val)
   use common_data
   implicit none
   real(8), intent(in) :: r, z, wa
   real(8) :: val, wc, wd, we
   wc = sqrt((r * (z + wa)**(1.0_8/3.0_8)) / 0.8853_8)
   wd = wc * (0.60112_8 * wc + 1.81061_8) + 1.0_8
   we = wc * (wc * (wc * (wc * (0.04793_8 * wc + 0.21465_8) + 0.77112_8) &
        + 1.39515_8) + 1.81061_8) + 1.0_8
   wc = (z + wa) * (wd / we)**2 - wa
   val = -wc / r
end function fpot

subroutine inouh(dp, dq, dr, dq1, dfl, dv, z, test, nuc, nstop, jc)
   use common_data
   use constants
   implicit none
   real(8), intent(inout) :: dp(nrad), dq(nrad)
   real(8), intent(in) :: dr(nrad), dq1, dfl, dv, z, test
   integer, intent(in) :: nuc, jc
   integer, intent(out) :: nstop
   real(8) :: dval, deva1, deva2, deva3, dbe, dsum, dqr, dpr, dm_local, dval_local
   integer :: i, m, j

   dp(1:10) = 0.0_8
   dq(1:10) = 0.0_8

   if (nuc > 0) then
      dval = dv + z * (3.0_8 - dr(1)**2 / dr(nuc)**2) / (dr(nuc) + dr(nuc))
      deva1 = 0.0_8
      deva2 = (dval - 3.0_8 * z / (dr(nuc) + dr(nuc))) / dvc - dd
      deva3 = z / (dr(nuc)**3 * dsal)
      if (dk < 0) then
         dp(10) = dq1
      else
         dq(10) = dq1
      end if
   else
      dval = z / dvc
      deva1 = -dval
      deva2 = dv / dvc + dval / dr(1) - dd
      deva3 = 0.0_8
      dbe = merge((dk - dfl) / dval, dval / (dk + dfl), dk <= 0)
      dq(10) = dq1
      dp(10) = dbe * dq1
   end if

   do i = 1, 5
      dp(i) = dp(10)
      dq(i) = dq(10)
      dep(i) = dp(i) * dfl
      deq(i) = dq(i) * dfl
   end do

   m = 1
   do while (m <= 20)
      dm_local = m + dfl
      dsum = dm_local**2 - dk**2 + deva1**2
      dqr = (dsal - deva2) * dq(m+9) - deva3 * dq(m+7)
      dpr = deva2 * dp(m+9) + deva3 * dp(m+7)
      dval_local = ((dm_local - dk) * dqr - deva1 * dpr) / dsum
      dsum = ((dm_local + dk) * dpr + deva1 * dqr) / dsum

      j = -1
      do i = 1, 5
         dpr = dr(i)**m
         dqr = dsum * dpr
         dpr = dval_local * dpr
         if (m > 1) then
            if (abs(dpr / dp(i)) <= test .and. abs(dqr / dq(i)) <= test) j = 1
         end if
         dp(i) = dp(i) + dpr
         dq(i) = dq(i) + dqr
         dep(i) = dep(i) + dpr * dm_local
         deq(i) = deq(i) + dqr * dm_local
      end do

      if (j == 1) exit
      dp(m+10) = dval_local
      dq(m+10) = dsum
      m = m + 1
   end do

   if (m > 20) nstop = 45
   dpno(1:4,jc) = dp(10:13)
   dqno(1:4,jc) = dq(10:13)
end subroutine inouh

subroutine insld(natom)
   use common_data
   use constants
   implicit none
   integer, intent(in) :: natom
   integer :: i, ll, lp1, j, ielec
   character(4) :: oplus(4) = ['s1/2', 'p3/2', 'd5/2', 'f7/2']
   character(4) :: ominus(3) = ['p1/2', 'd3/2', 'f5/2']

   dpas = delx
   nes = 150
   testv = 1.0e-07_8
   test = 1.0e-10_8
   nuc = 0
   dval = 0.0_8

   read(3, '(10a4)') bar(1:10)
   read(3, *) norb, icore, nitoe, itot, iex, iscfat, nitpot, icoul

   if (nitoe /= 0) nes = nitoe
   if (nitpot == 0) nitpot = 100
   z = nz(natom)
   np = nrad
   itatom = 0

   if (iskip /= 0) then
      select case (icoul)
         case (0)
            write(6, '(10x, "atom ", i4, " no tail"/)') natom
         case (1)
            write(6, '(10x, "atom ", i4, " screened coulomb tail"/)') natom
         case (2)
            write(6, '(10x, "atom ", i4, " unscreened coulomb tail"/)') natom
      end select
      write(6, '(1x, " atomic calculation for scatterer", i4//)') natom
      write(6, '(10x, "number of orbitals:", i4/10x, "number of core orbitals:", i4//)') &
         norb, icore
   end if

   do i = 1, norb
      read(3, *) den(i), nqn(i), nk(i), nel(i)
      nql(i) = abs(nk(i))
      if (nk(i) < 0) nql(i) = nql(i) - 1
      dfl(i) = sqrt(nk(i)**2 - dval)
      if (nk(i) >= 0) then
         do ll = 1, 3
            if (ll == nk(i)) then
               titre(i) = ominus(ll)
               exit
            end if
         end do
      else
         do ll = 1, 4
            lp1 = -ll
            if (lp1 == nk(i)) then
               titre(i) = oplus(ll)
               exit
            end if
         end do
      end if
      if (iskip /= 0) write(6, '(2i6, 5x, a4, i5, e15.8)') &
         nqn(i), nk(i), titre(i), nel(i), den(i)
   end do

   ielec = sum(nel(1:norb))
   nmax(1:norb) = np
   do i = 1, norb
      j = nqn(i) - nql(i)
      l = merge(-1, 1, mod(j, 2) == 0)
      dq1(i) = l * nk(i) / abs(nk(i))
   end do

   if (iscfat /= 0 .and. ielec /= nz(natom)) then
      write(6, '(10x, "get the input right, dummy!"/10x, " z = ", i4, &
                " number of electrons ", i4//)') nz(natom), ielec
      stop 'input screwed up'
   end if
end subroutine insld

subroutine inth(dp, dq, dv, dr)
   use common_data
   implicit none
   real(8), intent(inout) :: dp, dq
   real(8), intent(in) :: dv, dr
   real(8) :: dkoef1 = 475.0_8 / 502.0_8
   real(8) :: dkoef2 = 27.0_8 / 502.0_8
   real(8) :: dpr, dqr, dsum

   dpr = dp + dm * (251.0_8 * dep(1) + 2616.0_8 * dep(3) + 1901.0_8 * dep(5) &
        - (1274.0_8 * dep(2) + 2774.0_8 * dep(4)))
   dqr = dq + dm * (251.0_8 * deq(1) + 2616.0_8 * deq(3) + 1901.0_8 * deq(5) &
        - (1274.0_8 * deq(2) + 2774.0_8 * deq(4)))

   dep(1:4) = dep(2:5)
   deq(1:4) = deq(2:5)

   dsum = (db - dv / dvc) * dr
   dep(5) = -dk * dpr + (dsal * dr + dsum) * dqr
   deq(5) = dk * dqr - dsum * dpr

   dp = dp + dm * (106.0_8 * dep(2) + 646.0_8 * dep(4) + 251.0_8 * dep(5) &
        - (19.0_8 * dep(1) + 264.0_8 * dep(3)))
   dq = dq + dm * (106.0_8 * deq(2) + 646.0_8 * deq(4) + 251.0_8 * deq(5) &
        - (19.0_8 * deq(1) + 264.0_8 * deq(3)))

   dp = dkoef1 * dp + dkoef2 * dpr
   dq = dkoef1 * dq + dkoef2 * dqr
   dep(5) = -dk * dp + (dsal * dr + dsum) * dq
   deq(5) = dk * dq - dsum * dp
end subroutine inth

subroutine kkrin
   use common_data
   use constants
   implicit none
   integer :: n, k, nnk, nnc, iform, iold, kmax
   real(8) :: xrc, x0, vbar, rtemp(3), vtemp(3), dps

   write(6, *) ' enter kkrin!'
   iunit8 = 1
   read(8, '(20a4)') title
   read(8, *) aa0, conc, ilat, iform, iemode, iold
   read(8, *) stval, xdel

   write(6, '(20a4)') title
   write(6, *) aa0, conc, ilat, iform, iemode, iold
   write(6, *) stval, xdel

   if (stval == 0.0_8) stval = -8.8_8
   if (xdel == 0.0_8) xdel = 0.05_8

   conc1 = conc
   call mike(ilat)

   do n = 1, noat
      read(8, *) nt(n), nz(n), nnk(n), nrc(n), rc(n), vc(n), vcc(n)
      write(6, *) nt(n), nz(n), nnk(n), nrc(n), rc(n), vc(n), vcc(n)
      nnk = nnk(n)
      nnc = nrc(n)
      select case (iform)
         case (1)
            read(8, '(4d15.8)') (za(n,k), k=1,nnc)
         case (0)
            read(8, '(4d20.12)') (za(n,k), k=1,nnc)
         case (2)
            read(8, *) (za(n,k), k=1,nnc)
      end select

      xrc = log(rc(n))
      x0 = xrc - (nnc - 1) * xdel
      stval = x0
      if (iold <= 0 .and. iemode == 0) xdel = (-stval + log(rws1)) / (nnk - 6)

      do k = 1, nnk
         ra(n,k) = exp(x0)
         za(n,k) = za(n,k) / ra(n,k)
         x0 = x0 + xdel
      end do
      if (iold < 2) then
         do k = nnc + 1, nnk
            rtemp = ra(n,k-4:k-2)
            vtemp = za(n,k-4:k-2)
            call interp(rtemp, vtemp, 3, ra(n,k), za(n,k), dps, .false.)
         end do
      end if
   end do

   vbar = conc * vc(1) + (1.0_8 - conc) * vc(2)
   vc(1:2) = vbar

   do n = 1, noat
      kmax = nnk(n)
      do k = 1, kmax
         if (ra(n,k) > rws) exit
      end do
      nrws(n) = k - 1
      if (iemode /= 0) then
         nnk(n) = k
         do k = 1, kmax
            if (ra(n,k) > rws1) exit
         end do
         nrws1(n) = k - 1
      else
         nrws1(n) = nrc(n)
         rws1 = rc(n)
      end if
   end do

   write(6, *) ' kkrin left!'
end subroutine kkrin

subroutine poisat(natom)
   use common_data
   use constants
   implicit none
   integer, intent(in) :: natom
   real(8) :: s1(nrad), s2(nrad), s3, rs, rhos, rhox, vexks, vex, del0, del1, alpha
   integer :: i, np1, mft, np11

   real(8) :: beta, xsl
   beta(rs) = 1.0_8 + 0.0545_8 * rs * log(1.0_8 + 11.4_8 / rs)
   xsl = -6.0_8 * (3.0_8 / (8.0_8 * pi)) ** third

   vold(1:np) = dv(1:np)
   s1(1:np) = 0.0_8
   s2(1:np) = 0.0_8

   do i = 3, np
      if (rhotot(i,natom) == 0.0_8) exit
   end do
   np1 = i - 1

   s1(1:np1) = rhotot(1:np1,natom) / dr(1:np1)
   mft = np1 - 1
   do i = 2, mft
      l = mft - i + 1
      s2(l) = s2(l+1) + (s1(l) + s1(l+1)) * (dr(l+1) - dr(l))
   end do

   s3 = rhotot(1,natom) * dr(1)
   za(natom,1) = s2(1) + (s3 - 2.0_8 * z) / dr(1)
   do i = 2, np1
      s3 = s3 + (rhotot(i,natom) + rhotot(i-1,natom)) * (dr(i) - dr(i-1))
      za(natom,i) = s2(i) + (s3 - 2.0_8 * z) / dr(i)
   end do

   np11 = np1 + 1
   za(natom,np11:np) = (s3 - 2.0_8 * z) / dr(np11:np)

   do i = 2, np1
      rhos = rhotot(i,natom) / (dr(i)**2 * pi4)
      rhox = rhos ** third
      rs = (3.0_8 / (pi4 * rhos)) ** third
      vexks = twoth * xsl * rhox
      vex = merge(vexks * beta(rs), vexks, iex == 0)
      za(natom,i) = za(natom,i) + vex
   end do

   do i = 1, np
      dv(i) = za(natom,i) / 2.0_8
      if (dv(i) > -1.0_8 / dr(i)) dv(i) = -1.0_8 / dr(i)
   end do

   nstop = 1
   itatom = itatom + 1
   if (itatom > nitpot) return

   del0 = abs((dv(1) - vold(1)) / dv(1))
   do i = 2, np
      del1 = abs((dv(i) - vold(i)) / dv(i))
      if (del1 > del0) del0 = del1
   end do

   alpha = 0.75_8
   if (iskip /= 0) write(6, '(10x, "iteration ", i4, " delta V-max ", e15.8)') &
      itatom, del0
   if (del0 > testv) nstop = 0
   dv(1:np) = alpha * vold(1:np) + (1.0_8 - alpha) * dv(1:np)
   if (nstop == 0) write(6, '(10x, "iteration ", i4, " delta V-max ", e15.8)') &
      itatom, del0
end subroutine poisat

subroutine resld(nqn, nql, nk, imax, de, dfl, dq1, jc, natom)
   use common_data
   use constants
   implicit none
   integer, intent(in) :: nqn, nql, nk, jc, natom
   integer, intent(out) :: imax
   real(8), intent(inout) :: de, dfl, dq1
   real(8) :: elim, val, dval, dsum, dqm, dpm, dpq, dd, dbe, epriv
   integer :: i, j, k, nd, nodes, ies, imm, imat

   dkoef = 1.0_8 / 720.0_8
   nstop = 0
   dvc = 137.0373_8
   dsal = dvc + dvc

   epriv = de
   imm = 0
   ies = 0
   dk = nk
   lll = (nql * (nql + 1)) / 2
   nd = 0
   nodes = nqn - nql

   if (lll == 0) then
      elim = -z * z / (1.5_8 * nqn * nqn)
   else
      elim = dv(1) + lll / (dr(1) * dr(1))
      do i = 2, np
         val = dv(i) + lll / (dr(i) * dr(i))
         if (val < elim) elim = val
      end do
      if (elim > 0.0_8) then
         nstop = 17
         write(6, '(5x, "nstop = ", i4, "  2*v+l*(l+1)/r**2 is positive"/)') nstop
         stop 'in resld'
      end if
   end if

   if (iskip >= 2) write(6, '(" de=", d20.10, "elim=", d20.10)') de, elim
   if (de <= elim) de = elim * 0.5_8

   do
      do i = 7, np, 2
         imat = np + 1 - i
         if ((dv(imat) + lll / (dr(imat) * dr(imat)) - de) <= 0.0_8) exit
      end do
      if (imat > 5) exit
      de = de * 0.5_8
      if (de >= -test .or. nd > nodes) then
         nstop = 28
         write(6, '(5x, "nstop = ", i4, " 2*v+l*(l+1)/r**2-2*e is positive"/)') nstop
         return
      end if
   end do

   if (iskip >= 2) write(6, '(" de=", d20.10)') de
   db = de / dvc
   call inouh(dp, dq, dr, dq1, dfl, dv(1), z, test, nuc, nstop, jc)
   if (nstop /= 0) then
      write(6, '(5x, "nstop = ", i4, " dexpansion at the origin does not converge"/)') nstop
      return
   end if

   nd = 1
   do i = 1, 5
      dval = dr(i) ** dfl
      if (i > 1 .and. dp(i-1) /= 0.0_8 .and. (dp(i) / dp(i-1)) <= 0.0_8) nd = nd + 1
      dp(i) = dp(i) * dval
      dq(i) = dq(i) * dval
      dep(i) = dep(i) * dval
      deq(i) = deq(i) * dval
   end do

   k = -1 + 2 * mod(nodes, 2)
   if (dp(1) * k <= 0.0_8 .or. k * nk * dq(1) < 0.0_8) then
      nstop = 53
      write(6, '(5x, "nstop ", i4, " dexpansion error at the origin"/)') nstop
      return
   end if

   dm = dpas * dkoef
   do i = 6, imat
      dp(i) = dp(i-1)
      dq(i) = dq(i-1)
      call inth(dp(i), dq(i), dv(i), dr(i))
      if (dp(i-1) /= 0.0_8 .and. (dp(i) / dp(i-1)) <= 0.0_8) then
         nd = nd + 1
         if (nd > nodes) exit
      end if
   end do

   if (nd < nodes) then
      de = 0.8_8 * de
      if (de >= -test) then
         nstop = 206
         write(6, '(5x, "nstop = ", i4, " number of nodes is too small"/)') nstop
         return
      end if
      cycle
   else if (nd > nodes) then
      de = 1.2_8 * de
      if (de <= elim) then
         nstop = 210
         write(6, '(5x, "nstop = ", i4, " number of nodes is too large"/)') nstop
         return
      end if
      cycle
   end if

   dqm = dq(imat)
   dpm = dp(imat)
   if (imm == 0) then
      do i = 1, np, 2
         imax = np + 1 - i
         if (((dv(imax) - de) * dr(imax)**2) <= 300.0_8) exit
      end do
   end if

   dd = sqrt(-de * (2.0_8 + db / dvc))
   dpq = -dd / (dsal + db)
   dm = -dm
   do i = 1, 5
      j = imax + 1 - i
      dp(j) = exp(-dd * dr(j))
      dep(i) = -dd * dp(j) * dr(j)
      dq(j) = dpq * dp(j)
      deq(i) = dpq * dep(i)
   end do

   m = imax - 5
   do i = imat, m
      j = m + imat - i
      dp(j) = dp(j+1)
      dq(j) = dq(j+1)
      call inth(dp(j), dq(j), dv(j), dr(j))
   end do

   dval = dpm / dp(imat)
   if (dval <= 0.0_8) then
      nstop = 312
      write(6, '(5x, "nstop = ", i4, " sign error for the large component"/)') nstop
      return
   end if

   dp(imat:imax) = dp(imat:imax) * dval
   dq(imat:imax) = dq(imat:imax) * dval

   dsum = 0.0_8
   if (dp(1) /= 0.0_8) then
      dsum = 3.0_8 * dr(1) * (dp(1)**2 + dq(1)**2) / (dpas * (dfl + dfl + 1.0_8))
   end if
   do i = 3, imax, 2
      dsum = dsum + dr(i) * (dp(i)**2 + dq(i)**2) + &
             4.0_8 * dr(i-1) * (dp(i-1)**2 + dq(i-1)**2) + &
             dr(i-2) * (dp(i-2)**2 + dq(i-2)**2)
   end do
   dsum = dpas * (dsum + dr(imat) * (dqm**2 - dq(imat)**2)) / 3.0_8

   dbe = dp(imat) * (dqm - dq(imat)) * dvc / dsum
   if (abs(dbe / de) <= test) then
      dsum = sqrt(dsum)
      dq1 = dq1 / dsum
      dp(1:imax) = dp(1:imax) / dsum
      dq(1:imax) = dq(1:imax) / dsum
      dpno(1:4,jc) = dpno(1:4,jc) / dsum
      dqno(1:4,jc) = dqno(1:4,jc) / dsum
      if (imax < np) then
         dp(imax+1:np) = 0.0_8
         dq(imax+1:np) = 0.0_8
      end if
      nstop = 0
      return
   end if

   dval = de + dbe
   if (dval >= 0.0_8) then
      dbe = dbe * 0.5_8
      if (abs(dbe / de) > test) then
         de = de + dbe
         cycle
      end if
      nstop = 345
      write(6, '(5x, "nstop = ", i4, " energy converged to zero"/)') nstop
      return
   end if

   de = dval
   if (abs(de - epriv) < test) then
      dsum = sqrt(dsum)
      dq1 = dq1 / dsum
      dp(1:imax) = dp(1:imax) / dsum
      dq(1:imax) = dq(1:imax) / dsum
      dpno(1:4,jc) = dpno(1:4,jc) / dsum
      dqno(1:4,jc) = dqno(1:4,jc) / dsum
      if (imax < np) then
         dp(imax+1:np) = 0.0_8
         dq(imax+1:np) = 0.0_8
      end if
      nstop = 0
      return
   end if

   epriv = de
   if (iskip == 2) write(6, '(" orbital = ", i3, "iteration ", i3, " test", d13.5, &
      3x, "e0 = ", e15.8, 3x, "e1 = ", e15.8)') jc, ies, test, epriv, de
   if (abs(dbe / de) <= 0.1_8) imm = 1
   ies = ies + 1
   if (ies > nes) then
      nstop = 362
      write(6, '(5x, "nstop = ", i4, "number of iterations is too large"/)') nstop
      write(6, *) ies, nes
      stop
   end if
end subroutine resld

subroutine totec(natom)
   use common_data
   use constants
   implicit none
   integer, intent(in) :: natom
   real(8) :: eone, ecore, ucore, u1, t, u2, envxc, enexc, eatom, rho, rhox, rs, vxc, epsxc
   integer :: i, inf, inf1, iprint

   real(8) :: beta, eps, xks
   beta(rs) = 1.0_8 + 0.0545_8 * rs * log(1.0_8 + 11.4_8 / rs)
   eps(rs) = -0.0666_8 * ((1.0_8 + (rs / 11.4_8)**3) * log(1.0_8 + 11.4_8 / rs) &
             + 0.5_8 * rs / 11.4_8 - (rs / 11.4_8)**2 - 1.0_8 / 3.0_8)
   xks = -4.0_8 * (3.0_8 / (8.0_8 * pi)) ** third

   write(6, '(10x///10x, "separated atom total energy for scatterer [ryd]", i4/)') natom
   write(6, '(10x, a//)') merge('kohn-sham exchange', 'hedin-lundqvist exchange-correlation', iex /= 0)

   do i = 3, np
      if (rhotot(i,natom) == 0.0_8) exit
   end do
   inf1 = i - 1
   rinf1 = dr(inf1)
   inf = np
   rinf = dr(inf)

   eone = 2.0_8 * sum(nel(1:norb) * den(1:norb))
   ecore = 2.0_8 * sum(nel(1:icore) * den(1:icore))

   dq(1:np) = dv(1:np) * rhoc(1:np,natom)
   dp(1:np) = dv(1:np) * rhotot(1:np,natom)
   ucore = 2.0_8 * rsimp(dq, dr, rinf, inf)
   u1 = 2.0_8 * rsimp(dp, dr, rinf, inf)
   t = eone - u1
   tcore(natom) = ecore - ucore

   dp(2:np) = rhotot(2:np,natom) / dr(2:np)
   u2 = -2.0_8 * nz(natom) * rsimp(dp, dr, rinf, inf)

   iprint = 0
   do i = 2, inf
      rho = rhotot(i,natom) / (pi4 * dr(i)**2)
      rhox = rho ** third
      rs = (3.0_8 / (pi4 * rho)) ** third
      if (rs >= 9.0_8) then
         if (iprint == 0) write(6, '(10x, "cut-off at =", f10.5, " rs = ", e10.5//)') dr(i), rs
         iprint = 1
         vxc = xks * rhox
         epsxc = 0.75_8 * vxc
      else if (iex /= 0) then
         vxc = xks * rhox
         epsxc = 0.75_8 * vxc
      else
         vxc = xks * rhox * beta(rs)
         epsxc = 0.75_8 * xks * rhox + eps(rs)
      end if
      dp(i) = vxc * rhotot(i,natom)
      dq(i) = epsxc * rhotot(i,natom)
   end do

   envxc = rsimp(dp, dr, rinf1, inf1)
   enexc = rsimp(dq, dr, rinf1, inf1)
   u = (u1 + u2 - envxc) / 2.0_8
   eatom = t + u + enexc

   write(6, '(15x, "<E>       ", e20.10/15x, "<t>       ", e20.10/15x, &
             "<t-core>  ", e20.10/15x, "<u>       ", e20.10/15x, &
             "<E-xc>    ", e20.10/15x, "<Vxc>     ", e20.10//)') &
      eatom, t, tcore(natom), u, enexc, envxc

   if (iskip /= 0) then
      write(6, '(10x, "one electron energies "///)')
      do i = 1, norb
         write(6, '(15x, i2, a4, e20.8)') nqn(i), titre(i), den(i)
      end do
   end if
end subroutine totec

function rsimp(f, r, rn, irn) result(val)
   use common_data
   use constants
   implicit none
   real(8), intent(in) :: f(nrad), r(nrad), rn
   integer, intent(in) :: irn
   real(8) :: val, s, dx
   integer :: i, ieven, isw, np, nl

   dx = xdel
   val = 0.0_8
   if (irn <= 2) return

   ieven = irn / 2 * 2
   isw = merge(1, 0, ieven == irn)
   np = irn - isw
   s = f(1) * r(1) + f(np) * r(np)
   nl = np - 1
   do i = 2, nl, 2
      s = s + 4.0_8 * f(i) * r(i)
   end do
   nl = nl - 1
   if (nl >= 3) then
      do i = 3, nl, 2
         s = s + 2.0_8 * f(i) * r(i)
      end do
   end if
   s = s * dx / 3.0_8
   val = merge(s + (f(irn) * r(irn) + f(irn-1) * r(irn-1)) * 0.5_8 * dx, s, isw == 1)
end function rsimp

function rsimp3(f, r, rn, jrn) result(val)
   use common_data
   use constants
   implicit none
   real(8), intent(in) :: f(nrad), r(nrad), rn
   integer, intent(in) :: jrn
   real(8) :: val, rad(7), x1(7), s, dx
   integer :: iturn, irn, isw, np, nl, i

   dx = xdel
   do iturn = 1, 7
      irn = jrn + iturn - 6
      x1(iturn) = r(irn)
      s = 0.0_8
      if (irn <= 2) then
         rad(iturn) = 0.0_8
         cycle
      end if
      isw = merge(1, 0, irn / 2 * 2 == irn)
      np = irn - isw
      s = f(1) * r(1) + f(np) * r(np)
      nl = np - 1
      do i = 2, nl, 2
         s = s + 4.0_8 * f(i) * r(i)
      end do
      nl = nl - 1
      if (nl >= 3) then
         do i = 3, nl, 2
            s = s + 2.0_8 * f(i) * r(i)
         end do
      end if
      s = s * dx / 3.0_8
      rad(iturn) = merge(s + (f(irn) * r(irn) + f(irn-1) * r(irn-1)) * 0.5_8 * dx, s, isw == 1)
   end do

   call interp(x1, rad, 7, rn, val, val1, .false.)
end function rsimp3

subroutine fitrad(np)
   use common_data
   use constants
   implicit none
   integer, intent(in) :: np
   integer :: n, l, k, i, kk, kmax, kmax1, kmax3
   real(8) :: x, t, dps

   delx = xdel
   do n = 1, noat
      kmax = nnk(n)
      kmax1 = kmax + 1
      kmax3 = kmax1 - 7
      r(1) = 0.0_8
      p(1) = -2.0_8 * nz(n)
      do l = 1, kmax
         p(l+1) = za(n,l) * ra(n,l)
         r(l+1) = ra(n,l)
      end do
      i = 0
      x = stval - delx
      do while (i <= np)
         x = x + delx
         i = i + 1
         t = exp(x)
         do k = 1, kmax1
            if (r(k) > t) exit
         end do
         if (k <= kmax1) then
            kk = merge(1, merge(k-3, kmax3, k < kmax3), k <= 4)
            call interp(r(kk), p(kk), 7, t, za(n,i), dps, .false.)
            za(n,i) = za(n,i) / t
         end if
      end do
   end do
end subroutine fitrad

subroutine mike(ilat)
   use common_data
   use constants
   implicit none
   integer, intent(in) :: ilat

   lattyp = ilat
   select case (lattyp)
      case (0)
         avol = aa0 ** 3
         ec = 3.1166857_8
         rws1 = aa0 / 2.0_8
         write(6, '(10x, "simple cubic lattice"/)')
      case (1)
         avol = aa0 ** 3 / 4.0_8
         ec = 4.8320664_8
         rws1 = aa0 / sqrt(8.0_8)
         write(6, '(10x, "fcc lattice"/)')
      case (2)
         avol = aa0 ** 3 / 4.0_8
         ec = 4.085521_8
         rws1 = aa0 * sqrt(3.0_8 / 16.0_8)
         write(6, '(10x, "bcc lattice"/)')
   end select

   rws = (tofp * avol) ** third
   confac = (2.0_8 * pi / aa0) ** 2
end subroutine mike

subroutine interp(r, p, n, rs, ps, dps, deriv)
   implicit none
   integer, intent(in) :: n
   real(8), intent(in) :: r(n), p(n), rs
   real(8), intent(out) :: ps, dps
   logical, intent(in) :: deriv
   real(8) :: term, denom, dterm, dterm1
   integer :: i, j, k

   ps = 0.0_8
   dps = 0.0_8
   do j = 1, n
      term = 1.0_8
      denom = 1.0_8
      dterm = 0.0_8
      do i = 1, n
         if (i == j) cycle
         denom = denom * (r(j) - r(i))
         term = term * (rs - r(i))
         if (.not. deriv) cycle
         dterm1 = 1.0_8
         do k = 1, n
            if (k == j .or. k == i) cycle
            dterm1 = dterm1 * (rs - r(k))
         end do
         dterm = dterm + dterm1
      end do
      ps = ps + term * p(j) / denom
      if (deriv) dps = dps + dterm * p(j) / denom
   end do
end subroutine interp
