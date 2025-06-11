module aes_constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: ifilcbwf0 = 60
  character(len=40), parameter :: routine = 'aes'
  character(len=1), parameter :: shell(5) = ['k', 'l', 'm', 'n', 'o']
  character(len=3), parameter :: subsh(0:4) = ['1  ', '2,3', '4,5', '6,7', '8,9']
  character(len=3), parameter :: subshp(0:4) = ['1  ', '23 ', '45 ', '67 ', '89 ']
  character(len=4), parameter :: strme(12) = ['me1 ', 'me2 ', 'me3 ', 'me4 ', 'me5 ', 'me6 ', &
                                             'me7 ', 'me8 ', 'me9 ', 'me10', 'me11', 'me12']
end module aes_constants

subroutine aes(iprint, tsst, msst, mssq, tauq, mezz, mezj, gcor, fcor, &
               ecor, szcor, kapcor, mm05cor, nkpcor, ikmcor, izero, &
               itxray, bcor, bcors, ncstmax)
  use aes_constants
  use mod_energy, only: etab, nemax, efermi
  use mod_rmesh, only: nrmax, jrws, r2drdi
  use mod_calcmode, only: orbpol
  use mod_angmom, only: nkm, nmemax, nkmmax, nlmax, cgc, nlinq, nkmq, nlq
  use mod_sites, only: nq, iqat, nqmax
  use mod_files, only: datset, system, ldatset, ifilbuildbot, wrbuildbot, taudir
  use mod_types, only: soctl, ntmax, lcxray, ncxray, imt, nlt, nt, txt_t, ctl
  use mod_constants, only: ry_ev, c_au
  use omp_lib
  implicit none

  ! Dummy arguments
  integer, intent(in) :: iprint, itxray, ncstmax
  real(dp), intent(in) :: bcor(ntmax), bcors(ntmax), ecor(ncstmax)
  real(dp), intent(in) :: fcor(nrmax, 2, ncstmax), gcor(nrmax, 2, ncstmax)
  real(dp), intent(in) :: szcor(ncstmax)
  integer, intent(in) :: ikmcor(ncstmax, 2), izero(ncstmax), kapcor(ncstmax)
  integer, intent(in) :: mm05cor(ncstmax), nkpcor(ncstmax)
  complex(dp), intent(in) :: mezz(nkmmax, nkmmax, ntmax, nmemax)
  complex(dp), intent(in) :: mezj(nkmmax, nkmmax, ntmax, nmemax)
  complex(dp), intent(in) :: mssq(nkmmax, nkmmax, nqmax)
  complex(dp), intent(in) :: msst(nkmmax, nkmmax, ntmax)
  complex(dp), intent(in) :: tauq(nkmmax, nkmmax, nqmax)
  complex(dp), intent(inout) :: tsst(nkmmax, nkmmax, ntmax)

  ! Local variables
  real(dp), allocatable :: amecif(:, :, :), amecig(:, :, :), rint(:)
  real(dp) :: de, deme, ee, ei, epsd, imefin, jmc, mj, p3, q4
  real(dp) :: rj, sk, wa, wb, wc, wd, we, wf, wgte, xnorm(2)
  logical :: calcint, getirrsol, printame, printmele, readmele
  character(len=2) :: cl
  complex(dp), allocatable :: difme(:, :, :), ea, eb, ebot, ec, ed
  complex(dp), allocatable :: eme(:), etabfin(:, :), me(:, :, :, :, :)
  complex(dp), allocatable :: mmed(:, :, :), mmee(:, :, :), mmetildd(:, :, :)
  complex(dp), allocatable :: mmetilde(:, :, :), p, rme, ssst(:, :, :)
  complex(dp), allocatable :: taub(:, :), tauc(:, :), tautleed(:, :), iaes(:, :, :)
  complex(dp) :: waes(2), wtmp(2)
  character(len=80) :: filnam, spec
  integer :: i, ia_err, ib, ic, icst, icstme, id, ie, ieb, iec, ied
  integer :: ieme, ieme30, ieme40, iepanel, iepath, ifil, ifilb, ifilc
  integer :: ifild, ifilfin, ifilme, ifilval, iflag_tau_dir_inquire
  integer :: igrid_taudir(2), ii, ikmb, ikmc, ikmd, il, ilam, ilamp
  integer :: im, iq, is, it, itx, j, k, kap, lam, lam1, lam2, lam3
  integer :: lammagcpl(nkmmax), lamp, lamp1, lamp2, lamp3, mjm05, ms
  integer :: n, ncst, ne, nebot, nefin, neint, neme, nemeinp, nememax
  integer :: nepanel_taudir, nepath_taudir, netab_taudir(2)
  integer :: nkmb, nkmc, nkmd, nl, nlm, ntxrsgrp
  complex(dp) :: eryd = cmplx(999999.0_dp, 999999.0_dp, dp)
  integer :: iepanel = 1, iepath = 1

  ! Allocate arrays
  allocate(iaes(ncstmax, 2, 2*nemax-1), rint(nrmax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: iaes')
  allocate(mmed(nkmmax, nkmmax, nkmmax), mmee(nkmmax, nkmmax, nkmmax), &
           lammagcpl(nkmmax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: mmed/mmee')
  allocate(taub(nkmmax, nkmmax), tauc(nkmmax, nkmmax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: taub/tauc')
  allocate(mmetildd(nkmmax, nkmmax, nkmmax), mmetilde(nkmmax, nkmmax, nkmmax), &
           stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: mmetildd/mmetilde')
  allocate(tautleed(nkmmax, nkmmax), difme(nkmmax, nkmmax, nkmmax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: tautleed/difme')
  allocate(amecif(nkmmax, nkmmax, 2*nlmax), amecig(nkmmax, nkmmax, 2*nlmax), &
           stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: amecif/amecig')
  allocate(etabfin(2*nemax-1, 2), eme(nememax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: etabfin/eme')
  allocate(me(nkmmax, nkmmax, nkmmax, nememax, nememax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: me')
  allocate(ssst(nkmmax, nkmmax, ntmax), stat=ia_err)
  if (ia_err /= 0) call stop_message(routine, 'alloc: ssst')

  ! Initialize arrays
  iaes = (0.0_dp, 0.0_dp)
  mmed = (0.0_dp, 0.0_dp)
  mmee = (0.0_dp, 0.0_dp)
  mmetildd = (0.0_dp, 0.0_dp)
  mmetilde = (0.0_dp, 0.0_dp)
  difme = (0.0_dp, 0.0_dp)
  tautleed = (0.0_dp, 0.0_dp)
  taub = (0.0_dp, 0.0_dp)
  tauc = (0.0_dp, 0.0_dp)
  me = (0.0_dp, 0.0_dp)
  ssst = (0.0_dp, 0.0_dp)

  ! Read input
  call input_find_section('task', 1)
  call section_set_integer('it', itxray, 1, 0)
  it = itxray
  call section_set_integer('neme', nemeinp, 10, 0)
  nememax = nemeinp
  call section_find_keyword('printme', printmele)
  call section_find_keyword('printame', printame)
  call section_get_core_level_info(cl, ncxray(it), lcxray(it))

  iq = iqat(1, it)
  im = imt(it)

  ! Construct file name
  if (ldatset /= 0) then
    filnam = trim(datset) // trim(txt_t(it)) // '.'
  else
    filnam = trim(txt_t(it)) // '.'
  end if
  filnam = trim(filnam) // shell(ncxray(it))
  spec = '  ' // shell(ncxray(it))
  if (ncxray(it) /= 1) then
    spec = trim(spec) // subsh(lcxray(it))
    filnam = trim(filnam) // subshp(lcxray(it))
  end if
  filnam = trim(filnam) // '.aes'
  spec = trim(spec) // ' - aes spectrum of ' // trim(txt_t(it)) // ' in ' // trim(system)

  ! Open output file
  open(unit=7, file=trim(filnam))
  write(6, '(10x,a,a,/)') 'spec-file :  ( 7) ', trim(filnam)
  write(6, '(a)') trim(spec)
  write(6, '(4x,a,i2,a,i2,a,i2,/)') ' core quantum-numbers  for  it=', it, &
                                    ':   n=', ncxray(it), '  l=', lcxray(it)

  ! Increase angular momentum expansion by 1
  nl = 0
  do iq = 1, nq
    nlq(iq) = nlq(iq) + 1
    nl = max(nl, nlq(iq))
    if (nlq(iq) > nlmax) call stop_message(routine, 'nlq(iq)>nlmax')
    nkmq(iq) = 2 * nlq(iq)**2
    nlinq(iq) = 2 * nlq(iq) * (2 * nlq(iq) - 1)
  end do
  do itx = 1, nt
    nlt(itx) = nlq(iqat(1, itx))
    ctl(it, nlt(itx)) = ctl(it, nlt(itx) - 1)
    soctl(it, nlt(itx)) = soctl(it, nlt(itx) - 1)
  end do

  ! Set up angular momentum coupling
  lam = 0
  do k = 1, 2 * nlmax - 1
    if (mod(k, 2) == 0) then
      l = k / 2
      kap = l
    else
      l = (k - 1) / 2
      kap = -l - 1
    end if
    sk = real(sign(1, kap), dp)
    rj = real(l, dp) - sk / 2.0_dp
    do mjm05 = nint(-rj - 0.5_dp), nint(rj - 0.5_dp)
      mj = real(mjm05, dp) + 0.5_dp
      lam = lam + 1
      jmc = real(l, dp) + sk / 2.0_dp
      if (abs(mj) < l) then
        lammagcpl(lam) = nint(2 * l * (jmc + 0.5_dp) + jmc + mj + 1)
      else
        lammagcpl(lam) = 0
      end if
    end do
  end do
  nlm = nl**2
  nkm = 2 * nlm
  nkmb = 2 * (nl - 1)**2
  nkmc = nkmb
  nkmd = nkm

  ! Read tau to get vb-info
  ifilval = 99
  open(unit=ifilval, status='scratch', form='unformatted', access='direct', &
       recl=(nkmmax * nkmmax * 2 * 8))
  call tau_dir_inquire(taudir, nepanel_taudir, nepath_taudir, igrid_taudir, &
                       netab_taudir, iflag_tau_dir_inquire)
  ne = netab_taudir(1)
  if (nq > 1) call stop_message(routine, 'nq > 1')
  iq = 1
  !$omp parallel do private(ie, eryd)
  do ie = 1, ne
    eryd = etab(ie, 1)
    call tau_dir_read(taudir, eryd, iepanel, iepath, ie, mssq, tauq)
    etab(ie, 1) = eryd
    !$omp critical
    write(ifilval, rec=ie) tauq(:, :, iq)
    !$omp end critical
  end do
  !$omp end parallel do
  ebot = etab(1, 1)

  ! Write header
  call wrhead(7, filnam, 'sp-aes    ', ne)
  ntxrsgrp = 1
  write(7, '(a10,i10)') 'ntxrsgrp  ', ntxrsgrp
  write(7, '(a10,a)') 'spectrum  ', trim(spec)
  write(7, '(a10,i10)') 'it        ', it
  write(7, '(a10,i10)') 'ncxray    ', ncxray(it)
  write(7, '(a10,i10)') 'lcxray    ', lcxray(it)

  ! Prepare core states
  ncst = 4 * lcxray(it) + 2
  if (ncst > ncstmax) then
    write(6, '(//,60("*"),/,10x,a,/,10x,a,i3,a,i3)') ' stop in <aes>', &
      'lcxray=', lcxray(it), ' too large for ncstmax=', ncstmax
    call stop_message(routine, ' ')
  end if
  call core(iprint, gcor, fcor, ecor, szcor, kapcor, mm05cor, nkpcor, &
            ikmcor, izero, itxray, bcor, bcors, ncstmax)

  do ifil = 6, 7
    write(ifil, '(//,a,i4,/,a,20i4)') ' core states :', ncst, &
      ' nkpcor:', (nkpcor(icst), icst=1, ncst)
    write(ifil, '(/,a)') ' icst  n   l  kap  mue  ikm     norm   e(ry)       e(ev)    <sigma_z>  i0'
  end do

  !$omp parallel do private(icst, k, n, rint, xnorm)
  do icst = 1, ncst
    do k = 1, nkpcor(icst)
      do n = 1, jrws(im)
        rint(n) = r2drdi(n, im) * (gcor(n, k, icst)**2 + fcor(n, k, icst)**2)
      end do
      call rradint(im, rint, xnorm(k))
    end do
    !$omp critical
    do ifil = 6, 7
      write(ifil, '(5i4,a,i4,f12.6,f12.4,f12.3,f12.4,i5)') icst, ncxray(it), &
        lcxray(it), kapcor(icst), (2 * mm05cor(icst) + 1), ikmcor(icst, 1), &
        xnorm(1), ecor(icst), ecor(icst) * ry_ev, szcor(icst), izero(icst)
      if (nkpcor(icst) == 2) write(ifil, '(22x,i4,f12.6)') ikmcor(icst, 2), xnorm(2)
    end do
    !$omp end critical
  end do
  !$omp end parallel do

  write(6, '(//,10x,a,i5,/)') 'number of tabulated energies', ne
  write(7, '(/,a,i5,/,a,i5,/)') ' number of energies', ne, ' output format ifmt', 1

  ! Set up energy table for matrix elements
  neme = min(nemeinp, ne)
  deme = (efermi - real(ebot, dp)) / real(neme - 1, dp)
  calcint = .false.
  getirrsol = .false.
  do ie = 1, neme
    eme(ie) = ebot + deme * real(ie - 1, dp)
  end do

  ! Open files
  ifilfin = 98
  ifild = ifilcbwf0
  ifilb = ifilcbwf0 + 1
  ifilc = ifilcbwf0 + 2
  open(unit=ifilfin, status='scratch', form='unformatted', access='direct', &
       recl=(4 + 2 * (2 + 4 * nrmax)) * 2 * 8)
  open(unit=ifild, status='scratch', form='unformatted', access='direct', &
       recl=(4 + 2 * (2 + 4 * nrmax)) * 2 * 8)
  open(unit=ifilb, status='scratch', form='unformatted', access='direct', &
       recl=(4 + 2 * (2 + 4 * nrmax)) * 2 * 8)
  open(unit=ifilc, status='scratch', form='unformatted', access='direct', &
       recl=(4 + 2 * (2 + 4 * nrmax)) * 2 * 8)

  de = real(etab(2, 1) - etab(1, 1), dp)
  write(6, '(a,i3,/,a,f8.3,a,f8.3,a,/,a,f8.3,a,f8.3,a,/,a,f8.3,a,f8.3,a)') &
    ' tau:  number of energies read in  ne=', ne, &
    '       from   e(1)  =', real(etab(1, 1), dp), ' ry    (', real(etab(1, 1), dp) * ry_ev, ' ev)', &
    '       to     e(ne) =', real(etab(ne, 1), dp), ' ry    (', real(etab(ne, 1), dp) * ry_ev, ' ev)', &
    '       with   de    =', de, ' ry    (', de * ry_ev, ' ev)'

  ! Calculate angular matrix elements
  write(6, *) 'angular matrix elements will be calculated'
  call amecoul(amecig, amecif, iprint, nlmax, nkmmax)
  if (printame) then
    do il = 1, 2 * nlmax
      write(6, *) 'angular matrix elements for lr=', il
      call rmatstruct('amecig', amecig(1, 1, il), nkmd, nkmmax, 0, 0, 1, 1e-8_dp, 6)
      call rmatstruct('amecif', amecif(1, 1, il), nkmd, nkmmax, 0, 0, 1, 1e-8_dp, 6)
    end do
  end if

  ! Calculate SPR-AES spectrum
  write(6, *)
  write(6, *) 'calculating spr-aes spectrum'
  write(6, *)

  !$omp parallel do private(icst, ifilme, strfile, readmele, iec, ieb, i, j, icstme, ib, ic, id, rme, &
  !$omp                     ec, eb, ea, ee, ei, ed, p, ieme30, ieme40, p3, q4, wa, wb, wc, wd, we, wf, &
  !$omp                     ikmb, ikmc, ikmd, mmed, mmee, ied, nebot, neint, wgte, tautleed, taub, tauc, &
  !$omp                     mmetildd, mmetilde, difme, lamp, lamp1, lamp2, lamp3, lam1, lam2, lam3, waes, &
  !$omp                     wtmp, ms, epsd) firstprivate(ssst)
  do icst = 1, ncst
    write(6, *) '       core state ', icst
    write(6, *) '       -------------'
    me = (0.0_dp, 0.0_dp)
    mmed = (0.0_dp, 0.0_dp)
    mmee = (0.0_dp, 0.0_dp)

    ! Check for matrix element file
    ifilme = 23
    inquire(file='me', exist=readmele)
    if (readmele) then
      strfile = 'me'
    else
      strfile = strme(icst)
      inquire(file=trim(strfile), exist=readmele)
    end if
    printmele = .not. readmele
    if (readmele) then
      open(unit=ifilme, file=trim(strfile))
      write(6, *)
      write(6, *) '     reading meles from file !!!!!!!!!!! '
      do i = 1, 6
        read(ifilme, *)
      end do
      iec = 0
      do
        iec = iec + 1
        if (iec > neme) exit
        ieb = 0
        do
          ieb = ieb + 1
          do i = 1, 2
            read(ifilme, *)
          end do
          do j = 1, 1000000
            read(ifilme, '(i9,i5,i5,i4,2x,2(g21.14))', iostat=ia_err) icstme, ib, ic, id, rme
            if (ia_err /= 0) exit
            if (icstme == icst) me(ib, ic, id, ieb, iec) = rme
          end do
          if (ieb >= neme) exit
        end do
      end do
      close(ifilme)
      write(6, '(a,i3)') ' me: wave functions read in for neme=', neme
      write(6, '(a,f8.3,a,f8.3,a,/,a,f8.3,a,f8.3,a,/,a,f8.3,a,f8.3,a)') &
        ' me:   from   ebot  =', real(ebot, dp), ' ry    (', real(ebot, dp) * ry_ev, ' ev)', &
        '       to     efermi=', efermi, ' ry    (', efermi * ry_ev, ' ev)', &
        '       with   deme  =', deme, ' ry    (', deme * ry_ev, ' ev)')
    else
      ! Calculate matrix elements
      write(6, *)
      write(6, *) ' calculating me for (iea, ieb, iec, ied)'
      do iec = 1, neme
        ec = ebot + deme * real(iec - 1, dp)
        call ssite(1, 1, ifilc, calcint, getirrsol, ec, p, iprint, nkm, &
                   tsst, msst, ssst, mezz, mezj, orbpol)
        do ieb = 1, neme
          eb = ebot + deme * real(ieb - 1, dp)
          ea = ecor(icst)
          ee = real(eb + ec - ea, dp)
          ei = imag(eb)
          ed = cmplx(ee, ei, dp)
          if (iprint >= 0) write(6, '(28x,a,3(i3,a),i3,a)') '(', icst, ',', ieb, ',', iec, ',', iec + ieb - 1, ')')
          call ssite(1, 1, ifilb, calcint, getirrsol, eb, p, iprint, nkm, &
                     tsst, msst, ssst, mezz, mezj, orbpol)
          call ssite(1, 1, ifild, calcint, getirrsol, ed, p, iprint, nkm, &
                     tsst, msst, ssst, mezz, mezj, orbpol)
          call mecoul(icst, gcor, fcor, mm05cor, nkpcor, ikmcor, lcxray, 1, nkmb, ifilb, &
                      1, nkmc, ifilc, 1, nkmd, ifild, it, me, ieb, iec, amecig, amecif, &
                      nememax, ncstmax)
        end do
      end do
      ! Write matrix elements
      if (printmele) then
        ifilme = 24
        strfile = strme(icst)
        open(unit=ifilme, file=trim(strfile))
        write(6, *)
        write(6, *) ' matrix elements m(c,v,v,l) output in file me'
        write(ifilme, '(a)') '    coulomb matrix elements m(c,v,v,l)'
        write(ifilme, '(10x,a,i2)') 'k_c   = 1,', ncst
        write(ifilme, '(10x,a,i3)') 'k_v   = 1,', nkmb
        write(ifilme, '(10x,a,i3)') 'k_l   = 1,', nkmd
        write(ifilme, *)
        write(ifilme, '(2a)') '*****************************************', &
                              '*****************************'
        do iec = 1, neme
          do ieb = 1, neme
            ec = ebot + deme * real(iec - 1, dp)
            eb = ebot + deme * real(ieb - 1, dp)
            write(ifilme, '(a6,f10.6,a6,f10.6)') ' eb = ', real(eb, dp), ' ec = ', real(ec, dp)
            write(ifilme, '(6x,a17)') 'k_c k_vb k_vc k_l'
            do id = 1, nkmd
              do ic = 1, nkmc
                do ib = 1, nkmb
                  if (abs(me(ib, ic, id, ieb, iec)) > 1.0e-14_dp) then
                    write(ifilme, '(i9,i5,i5,i4,2x,2(g21.14))') icst, ib, ic, id, me(ib, ic, id, ieb, iec)
                  end if
                end do
              end do
            end do
            write(ifilme, '(2a)') '*****************************************', &
                                  '*****************************'
          end do
        end do
        close(ifilme)
      end if
    end if

    ! Set up energy table for final states
    imefin = imag(etab(1, 1))
    nefin = 2 * ne - 1
    do ie = nefin, 1, -1
      etabfin(ie, 1) = cmplx(2 * efermi - ecor(icst), imefin, dp) - &
                       cmplx(real(nefin - ie, dp) * de, 0.0_dp, dp)
    end do

    ! Calculate spectrum
    do ied = 1, 2 * ne - 1
      write(6, *) '  ied=', ied, '  ed=', real(etabfin(ied, 1), dp)
      if (ied <= ne) then
        nebot = 1
        neint = ied
      else
        nebot = ied + 1 - ne
        neint = ne
      end if
      ed = etabfin(ied, 1)
      do ieb = nebot, neint
        iec = ied + 1 - ieb
        wgte = merge(de, 2.0_dp * de, ieb == nebot .or. ieb == neint)
        eb = ebot + real(ieb - 1, dp) * de
        ec = ebot + real(iec - 1, dp) * de

        ! Interpolate matrix elements
        ieme30 = neme
        ieme40 = neme
        do ieme = 1, neme
          if (ieme30 == neme .and. real(eme(ieme), dp) > real(eb, dp)) &
            ieme30 = max(2, ieme - 1)
          if (ieme40 == neme .and. real(eme(ieme), dp) > real(ec, dp)) &
            ieme40 = max(2, ieme - 1)
        end do
        ieme30 = min(ieme30, neme - 1)
        ieme40 = min(ieme40, neme - 1)
        p3 = (real(eb, dp) - real(eme(ieme30), dp)) / deme
        q4 = (real(ec, dp) - real(eme(ieme40), dp)) / deme
        wa = q4 * (q4 - 1.0_dp) / 2.0_dp
        wb = p3 * (p3 - 1.0_dp) / 2.0_dp
        wc = 1.0_dp + p3 * q4 - p3**2 - q4**2
        wd = p3 * (p3 - 2.0_dp * q4 + 1.0_dp) / 2.0_dp
        we = q4 * (q4 - 2.0_dp * p3 + 1.0_dp) / 2.0_dp
        wf = p3 * q4

        do ikmb = 1, nkmb
          do ikmc = 1, nkmc
            do ikmd = 1, nkmd
              mmed(ikmb, ikmc, ikmd) = &
                me(ikmb, ikmc, ikmd, ieme30 + 0, ieme40 - 1) * wa + &
                me(ikmb, ikmc, ikmd, ieme30 - 1, ieme40 + 0) * wb + &
                me(ikmb, ikmc, ikmd, ieme30 + 0, ieme40 + 0) * wc + &
                me(ikmb, ikmc, ikmd, ieme30 + 1, ieme40 + 0) * wd + &
                me(ikmb, ikmc, ikmd, ieme30 + 0, ieme40 + 1) * we + &
                me(ikmb, ikmc, ikmd, ieme30 + 1, ieme40 + 1) * wf
              mmee(ikmb, ikmc, ikmd) = &
                me(ikmb, ikmc, ikmd, ieme40 + 0, ieme30 - 1) * wa + &
                me(ikmb, ikmc, ikmd, ieme40 - 1, ieme30 + 0) * wb + &
                me(ikmb, ikmc, ikmd, ieme40 + 0, ieme30 + 0) * wc + &
                me(ikmb, ikmc, ikmd, ieme40 + 1, ieme30 + 0) * wd + &
                me(ikmb, ikmc, ikmd, ieme40 + 0, ieme30 + 1) * we + &
                me(ikmb, ikmc, ikmd, ieme40 + 1, ieme30 + 1) * wf
            end do
          end do
        end do

        ! Get tau for final states
        eryd = etabfin(ied, 1)
        tsst = (0.0_dp, 0.0_dp)
        call ssite(1, 1, ifilfin, calcint, getirrsol, eryd, p, iprint, nkm, &
                   tsst, msst, ssst, mezz, mezj, orbpol)
        tautleed = tsst(:, :, it)
        read(ifilval, rec=ieb) taub
        read(ifilval, rec=iec) tauc

        ! Calculate tilde matrices
        mmetildd = (0.0_dp, 0.0_dp)
        mmetilde = (0.0_dp, 0.0_dp)
        do ikmb = 1, nkmb
          do ikmc = 1, nkmc
            do ilam = 1, nkmd
              do ilamp = 1, nkmd
                mmetildd(ikmb, ikmc, ilam) = mmetildd(ikmb, ikmc, ilam) + &
                  conjg(tautleed(ilamp, ilam)) * mmed(ikmb, ikmc, ilamp)
                mmetilde(ikmb, ikmc, ilam) = mmetilde(ikmb, ikmc, ilam) + &
                  conjg(tautleed(ilamp, ilam)) * mmee(ikmb, ikmc, ilamp)
              end do
            end do
          end do
        end do

        ! Calculate difference matrix
        difme = (0.0_dp, 0.0_dp)
        do lam1 = 1, nkmb
          do lam2 = 1, nkmc
            do lam3 = 1, nkmd
              difme(lam1, lam2, lam3) = mmetildd(lam1, lam2, lam3) - mmetilde(lam2, lam1, lam3)
            end do
          end do
        end do

        ! Sum over valence states
        wtmp = (0.0_dp, 0.0_dp)
        do lamp2 = 1, nkmb
          do lamp3 = 1, nkmb
            do lamp = 1, nkmc
              do lamp1 = 1, nkmc
                waes = (0.0_dp, 0.0_dp)
                do lam1 = 1, nkmd
                  do lam2 = 1, nkmd
                    if (lam1 == lam2 .or. lam2 == lammagcpl(lam1)) then
                      waes(1) = waes(1) + cgc(lam1, 1) * cgc(lam2, 1) * &
                                difme(lamp2, lamp, lam1) * conjg(difme(lamp3, lamp1, lam2))
                      waes(2) = waes(2) + cgc(lam1, 2) * cgc(lam2, 2) * &
                                difme(lamp2, lamp, lam1) * conjg(difme(lamp3, lamp1, lam2))
                    end if
                  end do
                end do
                do ms = 1, 2
                  wtmp(ms) = wtmp(ms) + imag(tauc(lamp, lamp1)) * &
                             imag(taub(lamp2, lamp3)) * waes(ms)
                end do
              end do
            end do
          end do
        end do
        do ms = 1, 2
          iaes(icst, ms, ied) = iaes(icst, ms, ied) + wgte * wtmp(ms)
        end do
      end do

      ! Apply energy-dependent factor
      epsd = 16.0_dp * (real(etabfin(ied, 1), dp) + c_au**2) / &
             (2.0_dp * real(etabfin(ied, 1), dp) + c_au**2)
      do ms = 1, 2
        iaes(icst, ms, ied) = iaes(icst, ms, ied) * epsd
      end do

      ! Write spectrum
      if (iprint > 0) then
        write(6, *) iaes(icst, 1, ied)
        write(6, *) iaes(icst, 2, ied)
      end if
      write(7, '(1x,e12.5,2x,e12.5,a)') real(etabfin(ied, 1), dp), &
        real(etabfin(ied, 1) - etabfin(nefin, 1), dp)
      write(7, '(20x,6(1x,e18.10))') (iaes(icst, is, ied), is=1, 2), &
        iaes(icst, 1, ied) + iaes(icst, 2, ied)

      ! Buildbot output
      if (wrbuildbot .and. ied <= 3) then
        write(ifilbuildbot, '(a,i5,a,i3,/,a,i5,a,2f10.6,/,1pe22.14)') &
          '# buildbot: ' // trim(routine) // ': aes intensity for it =', it, &
          '  icst =', icst, &
          '#          energy  ied =', ied, ' etabfin = ', etabfin(ied, 1), &
          (iaes(icst, is, ied), is=1, 2)
      end if
    end do
  end do
  !$omp end parallel do

  ! Clean up
  close(ifilval)
  close(ifilfin)
  close(ifild)
  close(ifilb)
  close(ifilc)
  deallocate(iaes, mmed, mmee, rint, taub, tauc, mmetildd, mmetilde)
  deallocate(tautleed, difme, lammagcpl, me, amecif, amecig)
  deallocate(etabfin, eme, ssst)
end subroutine aes
