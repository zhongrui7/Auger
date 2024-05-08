program rcvvaug
!
! **********************************************************
! CALCULATE RELATIVISTIC CORE-VALENCE-VALENCE AUGER SPECTRA
! **********************************************************
!
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NEMAX = 500 , NEMAX2 = 1000 , N1M = 300 , NEMF = 15 , NEMF2 = 30 , KM = 3 , KMTOT = 28
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: A0 , BLA , CONV , DE , DE2 , DUM , DUMMY , E0 , E1RYD , E2MAX , E2MAXEV , E2MIN , ECORE , EE , EEP , EF , EMAX ,&
                 & EMIN , ENERGY , ET0 , FL , FL1 , FL2 , HWC , HWLEV , HWLF , HWLR , HWS , PDIAGG , PMAX , PSUMMA , PSUMTOT ,     &
                 & PTOTT , SMAX , SNORM , U
   real(REAL64) , dimension(nemax,-km-1:km,-km-1:km) :: CRS
   real(REAL64) , dimension(nemf,nemf2,-km-1:km,-km-1:km) :: CRSFILE
   real(REAL64) , dimension(nemax,-km-1:km) :: DOS , UDOS
   real(REAL64) , dimension(nemax) :: E , Y
   real(REAL64) , dimension(n1m) :: E1 , P1DIAG , P1TOT , PHELP
   real(REAL64) , dimension(nemax2) :: E2 , PDIAG , PTOT
   real(REAL64) , save :: EFAC , PI , TINY
   real(REAL64) , dimension(nemf) :: EFILE
   real(REAL64) , dimension(nemf2) :: EFILE2
   character(50) :: FILEN , NAME
   integer :: I , IDOC , IE2 , IFST , II , ILAST , IMINB , IN , IP , J , JDOS , JJ , JMAT , JP , JPRT , JSPECT , K , KAP , KAPP ,  &
            & KP , KREL , KTEST , KTOT , L , LMAX , MAXD , N1 , NE , NE2 , NECOUNT , NEF , NUME , NUME2 , NVAL
   integer , dimension(kmtot) :: IKAP , IKAPP , IORD
   integer , dimension(nemf,nemf2) :: IRANGE
   integer , external :: IRIGHT , LEFTIND
   character(10) :: NAMORB
   real(REAL64) , dimension(nemax2,kmtot) :: P
   real(REAL64) , dimension(n1m,4) :: P1
   real(REAL64) , dimension(kmtot) :: PSUM
   real(REAL64) , external :: SIMPSON
   character(78) :: TITLE
   real(REAL64) , dimension(-km-1:km) :: Z
   external BROAD , DOSBROAD , MATINTER , ORDER
!
! End of declarations rewritten by SPAG
!
   data pi/3.1415926536/ , efac/13.606/ , tiny/1.0D-5/
!
!
!  open i/o files  *****************
!
   in = 1
   jdos = 2
   jmat = 3
   jprt = 7
   jspect = 8
!
   open (unit=in,file='rcvvaug.in')
!
! dos file
!
   read (in,99001) filen
   write (*,99001) filen
   open (unit=jdos,file=filen)
!
!
! matrix element file
!
   read (in,99001) filen
   write (*,99001) filen
   open (unit=jmat,file=filen)
!
!
! print file
!
   read (in,99001) filen
   write (*,99001) filen
   open (unit=jprt,file=filen)
!
!
! spectrum file
!
   read (in,99001) filen
   write (*,99001) filen
   open (unit=jspect,file=filen)
!
! ******************************
!
!
   read (in,99002) title
   write (*,99002) title
!
! spectrometer resolution (fwhm in eV)
!
   read (in,*) hws
   write (*,*) hws
!
!
! valence band lifetime broadening (fwhm in eV)
! (DOS will be broadened before convolution,
!  i.e. valence band holes are considered to be independent)
!
   read (in,*) hwlev
   write (*,*) hwlev
   hwlr = hwlev/efac
!
!
! core hole lifetime broadening (fwhm in eV)
!
   read (in,*) hwc
   write (*,*) hwc
!
!
! two-hole final state lifetime broadening (fwhm in eV)
!
   read (in,*) hwlf
   write (*,*) hwlf
!
!
! relativistic or nonrelativistic densities of states (1/0)
!
   read (in,*) krel
   write (*,*) krel
!
!
! shift the unbroadened spectrum downwards
!
   read (in,*) u
   write (*,*) u
!
!
! ***********  read in matrixelements ***************
!
   read (jmat,*)
   read (jmat,99001) name
   read (jmat,*)
   read (jmat,99003) namorb
   read (jmat,*)
   read (jmat,*) ecore
   read (jmat,*)
   read (jmat,*) e0 , de , nume
   read (jmat,*)
   read (jmat,*) lmax
   read (jmat,*)
   read (jmat,*) a0
   read (jmat,*)
   read (jmat,*)
!
   write (jprt,*) ' Matrixelements'
   write (jprt,99001) name
   write (jprt,99003) namorb
   write (jprt,99004) ((kap,kapp,kapp=-lmax-1,kap),kap=-lmax-1,-1) , ((kap,kapp,kapp=-lmax-1,-1),(kap,kapp,kapp=1,kap),kap=1,lmax)
99004 format (/1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,28(i6,i3,'    ')/)
!
   write (jspect,99002) title
   write (jspect,99001) name
   write (jspect,99010) hws , hwlev , hwc , hwlf
99010 format ('   HWS=',f5.2,' eV','    HWL=',f5.2,' eV','   HWC=',f5.2,' eV','   HWLF=',f5.2,' eV')
!
   conv = 2.0*pi/a0
   conv = conv*conv
   if ( krel==0 ) conv = 1.0
!
   nume2 = 2*nume - 1
   necount = nume*nume
   do i = 1 , nume
      efile(i) = e0 + de*(i-1)
   enddo
!
   et0 = 2.0*e0 - ecore
   do i = 1 , nume2
      efile2(i) = et0 + de*(i-1)
      do j = 1 , nume
         jp = i + 1 - j
         irange(j,i) = 1
         if ( jp<1 .or. jp>nume ) then
            irange(j,i) = 0
            cycle
         endif
!
         read (jmat,*) ii , jj , dum , eep , ee , ((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,kap),kap=-lmax-1,-1) ,                    &
                     & ((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,-1),(crsfile(jj,ii,kap,kapp),kapp=1,kap),kap=1,lmax)
!
         write (jprt,99005) i , j , efile2(i) , efile(jp) , efile(j) , ((crsfile(j,i,kap,kapp),kapp=-lmax-1,kap),kap=-lmax-1,-1) , &
                          & ((crsfile(j,i,kap,kapp),kapp=-lmax-1,-1),(crsfile(j,i,kap,kapp),kapp=1,kap),kap=1,lmax)
99005    format (2I3,3F6.2,28D13.5)
!
         if ( ii/=i .or. jj/=j .or. dabs(eep-efile(jp))>tiny .or. dabs(ee-efile(j))>tiny ) then
            write (6,*) i , j , efile(jp) , efile(j)
            write (6,*) ii , jj , eep , ee
            stop ' shit happened'
         endif
!
      enddo
   enddo
!
!
!
! ***************  read in densities of states  ***************
!
!
   read (jdos,99001) name
   read (jdos,*)
   read (jdos,*) ef
   ef = conv*ef
   read (jdos,*) lmax
   read (jdos,*) maxd
   read (jdos,*)
!
   ii = 0
   do i = 1 , maxd
!
      if ( krel==0 ) then
         read (jdos,*) energy , (z(l),l=0,lmax) , bla
         if ( bla>=tiny .or. ii/=0 ) then
            ii = ii + 1
            e(ii) = energy*conv
            udos(ii,-1) = z(0)
            do l = 1 , lmax
               fl = dfloat(2*l+1)
               fl = fl + fl
               fl1 = dfloat(2*l)
               fl2 = dfloat(2*l+2)
               fl1 = fl1/fl
               fl2 = fl2/fl
               udos(ii,l) = z(l)*fl1
               udos(ii,-l-1) = z(l)*fl2
            enddo
         endif
      else
         read (jdos,*) energy , bla , z(-1) , (z(l),z(-l-1),l=1,lmax)
         if ( bla>=tiny .or. ii/=0 ) then
            ii = ii + 1
            e(ii) = energy*conv
            do k = -lmax - 1 , lmax
               udos(ii,k) = z(k)
            enddo
         endif
      endif
!
   enddo
   ne = ii
!
   write (jprt,*)
   write (jprt,*) ' Densities of states'
   do i = 1 , ne
      write (jprt,99006) e(i) , udos(i,-1) , (udos(i,l),udos(i,-l-1),l=1,lmax)
   enddo
!
!
   nef = iright(e,ne,ef)
!
! check energy ranges:
   if ( e(1)<efile(1) .or. e(nef)>efile(nume) ) then
      write (6,*) ' energy range of matrix elements must exceed' , ' that of DOS!'
      stop
   endif
!
!
!  life time broadening of DOS **********************
!
   call dosbroad(ef,nef,e,udos,hwlr,dos,lmax)
!
   write (jprt,*)
   write (jprt,*) ' Densities of states after broadening'
   do i = 1 , nef
      write (jprt,99006) e(i) , dos(i,-1) , (dos(i,l),dos(i,-l-1),l=1,lmax)
   enddo
!
   ii = 0
   do k = -lmax - 1 , lmax
      if ( k/=0 ) then
         do kp = -lmax - 1 , k
            if ( kp/=0 ) then
               ii = ii + 1
               ikap(ii) = k
               ikapp(ii) = kp
            endif
         enddo
      endif
   enddo
   ktot = ii
   write (6,*) ' ktot=' , ktot
!
! determine boundaries for energy loop ***************
!
   e1ryd = ecore
   e2min = 2.0*e(1) - e1ryd
   e2max = 2.0*e(nef) - e1ryd
!
! set up final energy scale e2 (double grid has been set for array e2)
   de = e(2) - e(1)
   de2 = e(3) - e(1)
!      ne2=iidint((e2max-e2min)/de2)+1
   ne2 = int((e2max-e2min)/de2) + 1
   ee = e2min - de2
   do i = 1 , ne2
      ee = ee + de2
      e2(i) = ee
   enddo
!
! index values for test output of matrix elements:
   ktest = ne2/6
!
!
! start loop over final energies *********************
!
   do ie2 = 1 , ne2
      write (6,*) ie2 , (e2(ie2)+e1ryd)/2.
!
! valence band energy boundaries for this final energy:
      emax = dmin1(e1ryd+e2(ie2)-e(1),e(nef))
      emin = dmax1(e1ryd+e2(ie2)-e(nef),e(1))
      ilast = iright(e,nef,emax)
      ifst = leftind(e,nef,emin)
      nval = ilast - ifst + 1
!
! ********************************************************************
! interpolate matrix elements to DOS energy grid for this final energy
! ********************************************************************
!
! documentation:
      idoc = 0
      if ( mod(ie2,ktest)==0 ) idoc = 1
      if ( ie2==1 .or. ie2==2 ) idoc = 1
      if ( ie2==ne2-1 .or. ie2==ne2 ) idoc = 1
!
      call matinter(crsfile,efile,nume,efile2,nume2,nemf,nemf2,lmax,irange,e2(ie2),emin,emax,e,nemax,crs,ifst,ilast,idoc,jprt)
!
! *******************************************************************
!                   put together spectrum:
! *******************************************************************
!
!
      ii = 0
      do k = -lmax - 1 , lmax
         if ( k/=0 ) then
!
            do kp = -lmax - 1 , k
               if ( kp/=0 ) then
!
                  ii = ii + 1
!
                  do i = ifst , ilast
                     ip = ilast + ifst - i
                     y(i) = crs(i,k,kp)*dos(i,k)*dos(ip,kp)
                  enddo
                  p(ie2,ii) = simpson(y(ifst),nval,de)
               endif
!
            enddo
         endif
      enddo
!
!
! sum up diagonal contributions and total spectrum:
!
      pdiagg = 0.
      ptott = 0.
      do ii = 1 , ktot
         if ( ikap(ii)==ikapp(ii) ) pdiagg = pdiagg + p(ie2,ii)
         ptott = ptott + p(ie2,ii)
      enddo
      pdiag(ie2) = pdiagg
      ptot(ie2) = ptott
!
   enddo
!
!
! end of loop over final energies **********************************
!
!
! determine relative magnitudes of partial spectra:
!
   psumtot = 0.
   do ie2 = 1 , ne2
      psumtot = psumtot + ptot(ie2)
   enddo
   do ii = 1 , ktot
      psumma = 0.
      do ie2 = 1 , ne2
         psumma = psumma + p(ie2,ii)
      enddo
      psum(ii) = psumma/psumtot
   enddo
!
! maximum of total unbroadened spectrum:
!
   pmax = 0.
   do i = 1 , ne2
      if ( ptot(i)>pmax ) pmax = ptot(i)
   enddo
!
! ordering of partial spectra according to magnitude:
!
   call order(psum(1),ktot,iord)
!
!
! ********************************************************************
! BROADENING DUE TO LIFE TIME OF CORE HOLE AND SPECTROMETER RESOLUTION
! ********************************************************************
! (only the most important contributions are broadened)
!
   snorm = 0.
   write (6,*) ' broad'
   call broad(e2(1),ptot,nemax2,ne2,iminb,e1(1),p1tot,phelp,n1m,n1,snorm,smax,hws,hwlf,hwc,u)
   do i = 1 , 4
      write (6,*) ' broad'
      call broad(e2(1),p(1,iord(i)),nemax2,ne2,iminb,e1(1),p1(1,i),phelp,n1m,n1,snorm,dummy,hws,hwlf,hwc,u)
   enddo
   write (6,*) ' broad'
   call broad(e2(1),pdiag,nemax2,ne2,iminb,e1(1),p1diag,phelp,n1m,n1,snorm,dummy,hws,hwlf,hwc,u)
!
! ********************************************************************
!                      OUTPUT OF SPECTRA
! ********************************************************************
!
   e2maxev = e2max*efac
!
!
! unbroadened spectra
!
   write (jprt,*) ' unbroadened spectra'
   write (jprt,99007) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
   write (jprt,99008) (e2(i)-e2max,ptot(i),(p(i,iord(ii)),ii=1,4),pdiag(i),i=1,ne2)
!
!
! broadened spectra
!
   write (jspect,99003) namorb
   write (jspect,99006) ecore*efac
   write (jspect,99007) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
   write (jspect,99008) (e1(i)-e2maxev,p1tot(i)/snorm,(p1(i,ii)/snorm,ii=1,4),p1diag(i)/snorm,i=iminb,n1)
!
   write (jspect,99009) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
99009 format (4x,'e2',5x,4x,'Total',4(4x,i2,1x,i2),4x,'Diag.')
   write (jspect,99006) (e1(i)-e2maxev,p1tot(i),(p1(i,ii),ii=1,4),p1diag(i),i=iminb,n1)
!
!
   stop ' FORTRAN stop'
!
!
99001 format (a50)
99002 format (a78)
99003 format (a10)
99006 format (f11.5,6F9.3)
99007 format (4x,'e2',5x,5x,'Total',3x,4(5x,i2,1x,i2,3x),5x,'Diag.')
99008 format (f11.5,6E13.5)
end program RCVVAUG
!
!
subroutine dosbroad(ef,nef,e,udos,hwlr,dos,lmax)
! ====================================================
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NEMAX = 500 , KM = 3
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(in) :: EF
   integer :: NEF
   real(REAL64) , intent(in) , dimension(nemax) :: E
   real(REAL64) , intent(in) , dimension(nemax,-km-1:km) :: UDOS
   real(REAL64) , intent(in) :: HWLR
   real(REAL64) , intent(out) , dimension(nemax,-km-1:km) :: DOS
   integer , intent(in) :: LMAX
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) , dimension(nemax) :: A
   real(REAL64) , external :: CONVLO
   real(REAL64) :: ESTEP , GAM
   real(REAL64) , save :: HWLR0 , PI
   integer :: I , J , K
   integer , save :: ILOR
!
! End of declarations rewritten by SPAG
!
!
! life-time broadening of valence band states with Lorentzian
! only states below Efermi enter the integral
! DOS is required on an equidistant energy mesh!
!
   data pi/3.141596/ , ilor/2/ , hwlr0/0.1/
!
   write (6,*) ' I am in dosbroad'
   if ( hwlr<hwlr0 ) then
      do k = -lmax - 1 , lmax
         if ( k/=0 ) then
            do j = 1 , nef
               dos(j,k) = udos(j,k)
            enddo
         endif
      enddo
      write (6,*) ' I am out of dosbroad'
      return
   endif
!
   estep = e(2) - e(1)
   do k = -lmax - 1 , lmax
      if ( k/=0 ) then
         do i = 1 , nef
            do j = 1 , nef
               a(j) = udos(j,k)
            enddo
! halfwidth of Lor. increases linearly with binding E.:
            if ( ilor==1 ) gam = hwlr*(ef-e(i))/(ef-e(1))/2.0E0
! halfwidth of Lor. increases quadr.:
            if ( ilor==2 ) gam = hwlr*((ef-e(i))/(ef-e(1)))**2/2.0E0
            dos(i,k) = convlo(gam,a(1),i,estep,1,nef,pi)
         enddo
      endif
   enddo
!
   write (6,*) ' Out of dosbroad'
end subroutine DOSBROAD
!
!
subroutine broad(ein,a0,n0,max0,mine1,e1,a1,ai,n1m,n1,snorm,smax,hws,hwl,hwc,u)
! ==============================================================
!
!     interpolate to plot-grid and convolute with lorentzian
!     for life-time broadening of core hole + valence band holes
!     and gaussian spectrometer resolution
!
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NGC = 100 , NMAX = 500
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: N0
   integer , intent(in) :: N1M
   real(REAL64) , intent(in) , dimension(n0) :: EIN
   real(REAL64) , intent(inout) , dimension(n0) :: A0
   integer , intent(in) :: MAX0
   integer , intent(inout) :: MINE1
   real(REAL64) , intent(inout) , dimension(n1m) :: E1
   real(REAL64) , intent(inout) , dimension(n1m) :: A1
   real(REAL64) , dimension(n1m) :: AI
   integer , intent(inout) :: N1
   real(REAL64) , intent(inout) :: SNORM
   real(REAL64) , intent(inout) :: SMAX
   real(REAL64) , intent(in) :: HWS
   real(REAL64) , intent(in) :: HWL
   real(REAL64) , intent(in) :: HWC
   real(REAL64) , intent(in) :: U
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: AFAC , BFAC , DE0 , DELE , E00 , E0MAX , EMIN0 , EVMAX , EVMIN , GAM
   real(REAL64) , dimension(-ngc:ngc) :: AH , GC
   real(REAL64) , save :: CFAC , ESTEP , HWC0 , PI
   real(REAL64) , external :: CONVGAU , CONVLO
   real(REAL64) , dimension(nmax) :: E0
   integer :: I , I0 , J , MAX , NLST
   integer , save :: IX0
   external INTER
!
! End of declarations rewritten by SPAG
!
   data hwc0/5.E-03/ , cfac/13.606/
   data estep/0.05E0/
   data ix0/201/ , pi/3.141596/
!
!     convert from rydberg to ev:
!
   write (6,*) ' I am in broad'
   do i = 1 , max0
      e0(i) = ein(i)*cfac - u
   enddo
   de0 = e0(2) - e0(1)
   e0max = e0(max0)
   e00 = e0max
   do i = max0 + 1 , nmax
      e00 = e00 + de0
      e0(i) = e00
      a0(i) = 0.0
      if ( e00-e0max>u ) then
         call SPAG_BLOCK_1
         return
      endif
   enddo
   call SPAG_BLOCK_1
contains
   subroutine SPAG_BLOCK_1
      max = i
!udo  this next line is not in the orginal program
      nlst = n1
!
      if ( snorm<=0.E0 ) then
!
         evmin = e0(1) - 1.
         evmax = e0(max) + 1.
         spag_loop_1_2: do
            n1 = idnint((evmax-evmin)/estep) + 1
            if ( n1>n1m ) then
               estep = estep + estep
               cycle
            endif
! create nice numbers for energy:
            i0 = idnint(evmin/estep)
            evmin = i0*estep
!
!     set up plot-grid and determine min-energy
!
            nlst = n1
            write (6,*) ' I know what nlst is:' , nlst
            do i = 1 , n1
               e1(i) = evmin + (i-1)*estep
               if ( e1(i)>e0(max) .and. nlst==n1 ) nlst = i
            enddo
            spag_loop_2_1: do i = 1 , n0
               if ( a0(i)>0.E0 ) exit spag_loop_2_1
            enddo spag_loop_2_1
            exit spag_loop_1_2
         enddo spag_loop_1_2
         emin0 = e0(i)
         spag_loop_1_3: do i = 1 , n1
            if ( e1(i)>emin0 ) exit spag_loop_1_3
         enddo spag_loop_1_3
         mine1 = i
      endif
!
!
!     interpolate
!
      do i = mine1 , nlst
         j = 0
         spag_loop_2_4: do
            j = j + 1
            if ( e0(j)>=e1(i) .or. j>=max ) then
               j = j - 2
               if ( j<1 ) j = 1
               if ( j>max-3 ) j = max - 3
               call inter(e0(j),a0(j),4,e1(i),ai(i))
               exit spag_loop_2_4
            endif
         enddo spag_loop_2_4
      enddo
      do i = 1 , mine1 - 1
         ai(i) = 0.
      enddo
      do i = nlst + 1 , n1
         ai(i) = 0.
      enddo
!
      do i = mine1 , n1
!
!     determine halfwidth for life-time broadening
!     (2 final valence holes + 1 initial core hole)
!
         dele = (e1(nlst)-e1(i))/(e1(nlst)-e1(mine1))
         gam = (hwl*dele*dele+hwc)/2.0E0
!
!     convolute now
!
         a1(i) = convlo(gam,ai,i,estep,mine1,n1,pi)
      enddo
!
!     overwrite unbroadened data in ai with life-time
!     broadened spectrum
!
      do i = 1 , n1
         ai(i) = a1(i)
      enddo
      if ( hws/=0. ) then
!
!        set up gaussian convolution function
!
         afac = -alog(5.E-1)/(hws/2.0E0)**2
         bfac = sqrt(afac/pi)
         do i = -ngc , ngc
            gc(i) = exp(-afac*(i*estep)**2)*bfac
         enddo
!
!        convolute now
!
         do i = mine1 , n1
            a1(i) = convgau(ai,i,gc,ah,estep,ngc,n1)
         enddo
      endif
!
      if ( snorm<=0.E0 ) then
!
!     determine maximum and set equal to 100.
!
         smax = 0.E0
         do i = mine1 , n1
            if ( a1(i)>smax ) smax = a1(i)
         enddo
         snorm = 100.D0/smax
      endif
!
!
      do i = mine1 , n1
         a1(i) = a1(i)*snorm
!
      enddo
!
      write (6,*) ' Out of broad'
   end subroutine SPAG_BLOCK_1
end subroutine BROAD
!
!
function convlo(gam,a1,ind,estep,mine1,n1,pi)
! ==================================================================
   use iso_fortran_env
   implicit none
!
! Function and Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: N1
   real(REAL64) :: CONVLO
   real(REAL64) , intent(in) :: GAM
   real(REAL64) , intent(in) , dimension(n1) :: A1
   integer , intent(in) :: IND
   real(REAL64) , intent(in) :: ESTEP
   integer , intent(in) :: MINE1
   real(REAL64) , intent(in) :: PI
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: ATM , ATP , SUM
   integer :: I
!
! End of declarations rewritten by SPAG
!
!
!     do the convolution integral by trapezoidal rule
!
   write (6,*) ' Out of convlo'
   sum = 0.0
   do i = mine1 , n1
      atp = atan(estep*((i-ind)+0.5)/gam)
      atm = atan(estep*((i-ind)-0.5)/gam)
      sum = sum + a1(i)*(atp-atm)
   enddo
   convlo = sum/pi
   write (6,*) ' Out of convlo'
end function CONVLO
!
!
function convgau(a1,ind,gc,ah,estep,ngc,n1)
! ================================================================
   use iso_fortran_env
   implicit none
!
! Function and Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NGC
   integer , intent(in) :: N1
   real(REAL64) :: CONVGAU
   real(REAL64) , intent(in) , dimension(n1) :: A1
   integer , intent(in) :: IND
   real(REAL64) , intent(in) , dimension(-ngc:ngc) :: GC
   real(REAL64) , intent(inout) , dimension(-ngc:ngc) :: AH
   real(REAL64) , intent(in) :: ESTEP
!
! Local variable declarations rewritten by SPAG
!
   integer :: I , IMX
   real(REAL64) :: SUM
!
! End of declarations rewritten by SPAG
!
!
!     do the convolution integral by trapezoidal rule
!
!
   write (6,*) ' In convgau'
   do i = -ngc , ngc
      ah(i) = 0.D0
   enddo
   do i = 0 , ngc
      imx = ind + i
      if ( imx>n1 ) then
         call SPAG_BLOCK_1
         return
      endif
      ah(i) = a1(imx)*gc(i)
   enddo
   call SPAG_BLOCK_1
contains
   subroutine SPAG_BLOCK_1
      do i = -1 , -ngc , -1
         imx = ind + i
         if ( imx<1 ) then
            call SPAG_BLOCK_2
            return
         endif
         ah(i) = a1(imx)*gc(i)
      enddo
      call SPAG_BLOCK_2
   end subroutine SPAG_BLOCK_1
   subroutine SPAG_BLOCK_2
      sum = 0.D0
      do i = -ngc , ngc
         sum = sum + ah(i)
      enddo
!
      convgau = sum*estep
!
      write (6,*) ' Out of convgau'
   end subroutine SPAG_BLOCK_2
end function CONVGAU
!
!
subroutine matinter(crsfile,efile,nume,efile2,nume2,nemf,nemf2,km,irange,e2,emin,emax,e,nemax,crs,ifst,ilast,idoc,jprt)
! ========================================================
!
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: KML = 3 , NM = 50
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NEMF
   integer , intent(in) :: NEMF2
   integer , intent(in) :: NEMAX
   real(REAL64) , intent(in) , dimension(nemf,nemf2,-kml-1:kml,-kml-1:kml) :: CRSFILE
   real(REAL64) , dimension(nemf) :: EFILE
   integer :: NUME
   real(REAL64) , dimension(nemf2) :: EFILE2
   integer :: NUME2
   integer , intent(in) :: KM
   integer , intent(in) , dimension(nemf,nemf2) :: IRANGE
   real(REAL64) :: E2
   real(REAL64) :: EMIN
   real(REAL64) :: EMAX
   real(REAL64) , dimension(nemax) :: E
   real(REAL64) , dimension(nemax,-kml-1:kml,-kml-1:kml) :: CRS
   integer , intent(in) :: IFST
   integer , intent(in) :: ILAST
   integer , intent(in) :: IDOC
   integer , intent(in) :: JPRT
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) , dimension(nm) :: CHELP , EHELP
   real(REAL64) , dimension(nm,-kml-1:kml,-kml-1:kml) :: CRS1
   real(REAL64) :: EBAR
   integer :: I , I2 , IE , IEXTRA , IF , IF2L , IF2R , IFL , IFR , INSLEFT , INSRIGHT , K , KP , L , LP , NF , NINT
   integer , external :: IRIGHT , LEFTIND
   external LAGRANGE
!
! End of declarations rewritten by SPAG
!
!
! interpolation of CVV matrix elements on DOS energy grid for
! a special final state energy
!
! input: crsfile      matrix elements in the of input file
!        efile        valence band energy of input file
!        nume         number of VB energies actually read in
!        efile2       final state energy of input file
!        nume2        number of FS energies read in
!        nemf,nemf2   array dimensions
!        irange       =1, if crs has been calculated at a special
!                     e and e2; 0 otherwise
!        e2           final state energy
!        e            DOS energy grid
!        nemax        dimension of e
!        ifst, ilast  indices of e between which crs should be
!                     interpolated
!        idoc         if =1: documentation
!
! output: crs
!
!
   write (6,*) ' In matinter'
!
! find efile2 points immediately left and right
! of the actual e2 (if2l,if2r):
   if2l = leftind(efile2,nume2,e2)
   if2r = iright(efile2,nume2,e2)
!
! valence band energy range: find efile points immediately
! left of emin and right of emax (ifl,ifr):
   ifl = leftind(efile,nume,emin)
   ifr = iright(efile,nume,emax)
!
! for interpolation lateron, we may have to calculate the boundary
! values of the crs in a different way (to understand the whole
! thing, draw a diagram e vs. e2)
   insleft = 0
   insright = 0
   if ( irange(ifl,if2r)==0 ) insleft = 1
   if ( irange(ifr,if2l)==0 ) insright = 1
!
! now calculate the crs's for e2 and all values of efile
! between ifl and ifr and store them in crs1:
!
! boundary values (interpolation along boundaries of e/e2 range):
   nint = 5
   if ( nume<nint ) nint = nume
   if ( insleft==1 ) then
      ebar = (efile2(nume2)+efile2(1))/2.
      ehelp(ifl) = e2 - ebar + efile(1)
!        if(ehelp(ifl).ge.emin)then
!          ifl=ifl+1
!          goto 295
!        endif
      do l = -km - 1 , km
         if ( l/=0 ) then
            do lp = -km - 1 , l
               if ( lp/=0 ) then
                  do i = nume , 1 , -1
                     i2 = nume2 - nume + i
                     chelp(i) = crsfile(i,i2,l,lp)
                  enddo
                  call lagrange(efile(1),chelp(1),nume,nint,ehelp(ifl),crs1(ifl,l,lp),iextra)
               endif
            enddo
         endif
      enddo
      if ( idoc==1 ) write (jprt,99001)
99001 format (' non-grid value for lower e boundary')
   endif
   if ( insright==1 ) then
      ehelp(ifr) = e2 - efile2(1) + efile(1)
!        if(ehelp(ifr).le.emax)then
!          ifl=ifl-1
!          goto 297
!        endif
      do l = -km - 1 , km
         if ( l/=0 ) then
            do lp = -km - 1 , l
               if ( lp/=0 ) then
                  do i = 1 , nume
                     chelp(i) = crsfile(i,i,l,lp)
                  enddo
                  call lagrange(efile(1),chelp(1),nume,nint,ehelp(ifr),crs1(ifr,l,lp),iextra)
               endif
            enddo
         endif
      enddo
      if ( idoc==1 ) write (jprt,99002)
99002 format (' non-grid value for upper e boundary')
   endif
   if ( idoc==1 ) write (jprt,99003) ifl , ifr
99003 format (' ifl,ifr: ',2I4)
!
! get the rest of crs1:
   do if = ifl , ifr
      if ( if/=ifl .or. insleft/=1 ) then
         if ( if/=ifr .or. insright/=1 ) then
            ehelp(if) = efile(if)
            do l = -km - 1 , km
               if ( l/=0 ) then
                  do lp = -km - 1 , l
                     if ( lp/=0 ) then
                        do i = if , if + nume - 1
                           chelp(i) = crsfile(if,i,l,lp)
                        enddo
                        call lagrange(efile2(if),chelp(if),nume,nint,e2,crs1(if,l,lp),iextra)
                     endif
                  enddo
               endif
            enddo
         endif
      endif
   enddo
!
   if ( idoc==1 ) then
      write (jprt,99004)
99004 format (' crs1:')
      do i = ifl , ifr
         write (jprt,99005) i , ehelp(i) , ((crs1(i,l,lp),lp=-km-1,l),l=-km-1,-1) ,                                                &
                          & ((crs1(i,l,lp),lp=-km-1,-1),(crs1(i,l,lp),lp=1,l),l=1,km)
99005    format (i3,f10.5,28E13.5)
      enddo
   endif
 
!
! now we have a function crs1 on a mesh efile and we want a function
! crs on a mesh e and that's all
   nf = ifr - ifl + 1
   nint = 6
   if ( nint>nf ) nint = nf
   do ie = ifst , ilast
      do l = -km - 1 , km
         if ( l/=0 ) then
            do lp = -km - 1 , l
               if ( lp/=0 ) then
                  do i = ifl , ifr
                     chelp(i) = crs1(i,l,lp)
                  enddo
                  call lagrange(efile(ifl),chelp(ifl),nf,nint,e(ie),crs(ie,l,lp),iextra)
                  if ( iextra==1 .and. idoc==1 ) write (jprt,99006) ie , e(ie) , ((crs(ie,k,kp),kp=-km-1,k),k=-km-1,-1) ,          &
                     & ((crs(ie,k,kp),kp=-km-1,-1),(crs(ie,k,kp),kp=1,k),k=1,km)
99006             format (' extrapolation: crs, ie=',i3,' e(ie)=',f10.5/28E13.5)
               endif
            enddo
         endif
      enddo
   enddo
!
! test output:
   if ( idoc==1 ) then
      write (jprt,99007) e2
99007 format (' e2: ',f10.5)
      do ie = ifst , ilast
         write (jprt,99008) ie , e(ie) , ((crs(ie,l,lp),lp=-km-1,l),l=-km-1,-1) , ((crs(ie,l,lp),lp=-km-1,-1),(crs(ie,l,lp),lp=1,l)&
                          & ,l=1,km)
99008    format (i3,f10.5,28E13.4)
      enddo
      write (jprt,99009)
99009 format (120('-'))
   endif
   write (6,*) ' out of matinter'
end subroutine MATINTER
!
!
subroutine lagrange(x,y,ndim,npkt,xval,yinter,iextra)
! =========================================================
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NDIM
   real(REAL64) , dimension(ndim) :: X
   real(REAL64) , dimension(ndim) :: Y
   integer :: NPKT
   real(REAL64) :: XVAL
   real(REAL64) :: YINTER
   integer , intent(out) :: IEXTRA
!
! Local variable declarations rewritten by SPAG
!
   integer :: IFST , IGR , ILEFT
   external INTER
!
! End of declarations rewritten by SPAG
!
!
! Interpolation mit Komfort
! iextra = 0, wenn keine Extrapolation
!          1, wenn Extrapolation nach oben
!         -1, wenn Extrapolation nach unten
!
   write (6,*) ' in lagrange'
   iextra = 0
   if ( xval<x(1) ) iextra = -1
   if ( xval>x(ndim) ) iextra = 1
!
! Suche 1. Punkt rechts von xval, falls es ihn gibt:
   igr = 0
   spag_loop_1_1: do
      igr = igr + 1
      if ( x(igr)==xval ) then
         yinter = y(igr)
         return
      endif
      if ( x(igr)>=xval .or. igr>=ndim ) then
!
! Bestimme 1. Punkt f. Interpolation
         ileft = npkt/2
         if ( mod(npkt,2)/=0 ) ileft = ileft + 1
         ifst = igr - ileft
!
         if ( ifst<=0 ) ifst = 1
         if ( ifst>ndim-npkt+1 ) ifst = ndim - npkt + 1
!
! Lasst uns interpolieren:
         call inter(x(ifst),y(ifst),npkt,xval,yinter)
! Und das war's dann auch schon, Leute!
         write (6,*) ' Out of Lagrange'
         exit spag_loop_1_1
      endif
   enddo spag_loop_1_1
end subroutine LAGRANGE
!
!
subroutine inter(r,p,n,rs,ps)
! =================================
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: N
   real(REAL64) , intent(in) , dimension(n) :: R
   real(REAL64) , intent(in) , dimension(n) :: P
   real(REAL64) , intent(in) :: RS
   real(REAL64) , intent(inout) :: PS
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: DENOM , TERM
   integer :: I , J
!
! End of declarations rewritten by SPAG
!
!
!     interpolate via lagrange
!
   write (6,*) ' In inter'
   ps = 0.E0
   do j = 1 , n
      term = 1.E0
      denom = 1.E0
      do i = 1 , n
         if ( i/=j ) then
            denom = denom*(r(j)-r(i))
            term = term*(rs-r(i))
         endif
      enddo
      ps = ps + term*p(j)/denom
   enddo
   write (6,*) ' Out of inter'
end subroutine INTER
!
!
function leftind(field,ndim,value)
   use iso_fortran_env
   implicit none
!
! Function and Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NDIM
   integer :: LEFTIND
   real(REAL64) , intent(in) , dimension(ndim) :: FIELD
   real(REAL64) , intent(in) :: VALUE
!
! Local variable declarations rewritten by SPAG
!
   integer :: I
!
! End of declarations rewritten by SPAG
!
! ==============================================
! index of the field element immediately left of value
! if value is lower than field(1), 1 is returned
!      implict real*8 (a-h,o-z)
   do i = ndim , 1 , -1
      if ( field(i)<value ) then
         call SPAG_BLOCK_1
         return
      endif
   enddo
   call SPAG_BLOCK_1
contains
   subroutine SPAG_BLOCK_1
      leftind = i
      if ( leftind<1 ) leftind = 1
   end subroutine SPAG_BLOCK_1
end function LEFTIND
!
!
function iright(field,ndim,value)
   use iso_fortran_env
   implicit none
!
! Function and Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NDIM
   integer :: IRIGHT
   real(REAL64) , intent(in) , dimension(ndim) :: FIELD
   real(REAL64) , intent(in) :: VALUE
!
! Local variable declarations rewritten by SPAG
!
   integer :: I
!
! End of declarations rewritten by SPAG
!
! =============================================
! index of the field element immediately right of value
! if value is greater than fieldhndim), ndim is returned
!      implict real*8 (a-h,o-z)
   do i = 1 , ndim
      if ( field(i)>value ) then
         call SPAG_BLOCK_1
         return
      endif
   enddo
   call SPAG_BLOCK_1
contains
   subroutine SPAG_BLOCK_1
      iright = i
      if ( iright>ndim ) iright = ndim
   end subroutine SPAG_BLOCK_1
end function IRIGHT
!
!
function simpson(fct,n,dx)
! ===============================================
   use iso_fortran_env
   implicit none
!
! Function and Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: N
   real(REAL64) :: SIMPSON
   real(REAL64) , intent(in) , dimension(n) :: FCT
   real(REAL64) , intent(in) :: DX
!
! Local variable declarations rewritten by SPAG
!
   integer :: I
   real(REAL64) :: SUM , X0 , X1 , X2
!
! End of declarations rewritten by SPAG
!
! integration of fct using simpson integration
! equidistant point mesh
   sum = 0.
   x0 = fct(1)
!
   do i = 3 , n , 2
      x1 = fct(i-1)
      x2 = fct(i)
      sum = sum + dx*(x0+4*x1+x2)
      x0 = x2
   enddo
   simpson = sum/3.
   if ( mod(n,2)==0 ) simpson = simpson + (fct(n)+fct(n-1))*dx/2.
end function SIMPSON
!
!
subroutine order(f,n,ior)
! =============================
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: N
   real(REAL64) , intent(inout) , dimension(n) :: F
   integer , intent(inout) , dimension(n) :: IOR
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: G
   integer :: I , J , L
!
! End of declarations rewritten by SPAG
!
!
! f is to be arranged in decreasing order
   do i = 1 , n
      ior(i) = i
   enddo
!
   do i = 1 , n - 1
      do j = i + 1 , n
         if ( f(j)>f(i) ) then
            g = f(j)
            f(j) = f(i)
            f(i) = g
            l = ior(j)
            ior(j) = ior(i)
            ior(i) = l
         endif
      enddo
   enddo
end subroutine ORDER
