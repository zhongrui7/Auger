PROGRAM rcvvaug
!
! **********************************************************
! CALCULATE RELATIVISTIC CORE-VALENCE-VALENCE AUGER SPECTRA
! **********************************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a0 , bla , conv , crs , crsfile , de , de2 , dos , dum , dummy , e , e0 , e1 , e1ryd , e2 , e2max , e2maxev , e2min ,    &
        & ecore , ee
   REAL*8 eep , ef , efac , efile , efile2 , emax , emin , energy , et0 , fl , fl1 , fl2 , hwc , hwlev , hwlf , hwlr , hws , p ,   &
        & p1 , p1diag
   REAL*8 p1tot , pdiag , pdiagg , phelp , pi , pmax , psum , psumma , psumtot , ptot , ptott , simpson , smax , snorm , tiny , u ,&
        & udos , y , z
   INTEGER i , idoc , ie2 , ifst , ii , ikap , ikapp , ilast , iminb , in , iord , ip , irange , iright , j , jdos , jj , jmat ,   &
         & jp , jprt
   INTEGER jspect , k , kap , kapp , KM , KMTOT , kp , krel , ktest , ktot , l , leftind , lmax , maxd , n1 , N1M , ne , ne2 ,     &
         & necount , nef
   INTEGER NEMAX , NEMAX2 , NEMF , NEMF2 , nume , nume2 , nval
!*** End of declarations inserted by SPAG
!
   PARAMETER (NEMAX=500,NEMAX2=1000,N1M=300,NEMF=15,NEMF2=30)
   PARAMETER (KM=3,KMTOT=28)
!
   CHARACTER name*50 , title*78 , namorb*10 , filen*50
!
   DIMENSION e1(N1M) , e(NEMAX) , efile(NEMF) , efile2(NEMF2) , e2(NEMAX2) , irange(NEMF,NEMF2)
   DIMENSION z(-KM-1:KM) , udos(NEMAX,-KM-1:KM) , dos(NEMAX,-KM-1:KM)
   DIMENSION crsfile(NEMF,NEMF2,-KM-1:KM,-KM-1:KM) , crs(NEMAX,-KM-1:KM,-KM-1:KM)
   DIMENSION p(NEMAX2,KMTOT) , psum(KMTOT)
   DIMENSION pdiag(NEMAX2) , ptot(NEMAX2)
   DIMENSION iord(KMTOT) , ikap(KMTOT) , ikapp(KMTOT)
   DIMENSION y(NEMAX) , phelp(N1M)
   DIMENSION p1(N1M,4) , p1tot(N1M) , p1diag(N1M)
!
   DATA pi/3.1415926536/ , efac/13.606/ , tiny/1.0D-5/
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
   OPEN (UNIT=in,FILE='rcvvaug.in')
!
! dos file
!
   READ (in,99001) filen
   WRITE (*,99001) filen
   OPEN (UNIT=jdos,FILE=filen)
!
!
! matrix element file
!
   READ (in,99001) filen
   WRITE (*,99001) filen
   OPEN (UNIT=jmat,FILE=filen)
!
!
! print file
!
   READ (in,99001) filen
   WRITE (*,99001) filen
   OPEN (UNIT=jprt,FILE=filen)
!
!
! spectrum file
!
   READ (in,99001) filen
   WRITE (*,99001) filen
   OPEN (UNIT=jspect,FILE=filen)
!
! ******************************
!
!
   READ (in,99002) title
   WRITE (*,99002) title
!
! spectrometer resolution (fwhm in eV)
!
   READ (in,*) hws
   WRITE (*,*) hws
!
!
! valence band lifetime broadening (fwhm in eV)
! (DOS will be broadened before convolution,
!  i.e. valence band holes are considered to be independent)
!
   READ (in,*) hwlev
   WRITE (*,*) hwlev
   hwlr = hwlev/efac
!
!
! core hole lifetime broadening (fwhm in eV)
!
   READ (in,*) hwc
   WRITE (*,*) hwc
!
!
! two-hole final state lifetime broadening (fwhm in eV)
!
   READ (in,*) hwlf
   WRITE (*,*) hwlf
!
!
! relativistic or nonrelativistic densities of states (1/0)
!
   READ (in,*) krel
   WRITE (*,*) krel
!
!
! shift the unbroadened spectrum downwards
!
   READ (in,*) u
   WRITE (*,*) u
!
!
! ***********  read in matrixelements ***************
!
   READ (jmat,*)
   READ (jmat,99001) name
   READ (jmat,*)
   READ (jmat,99003) namorb
   READ (jmat,*)
   READ (jmat,*) ecore
   READ (jmat,*)
   READ (jmat,*) e0 , de , nume
   READ (jmat,*)
   READ (jmat,*) lmax
   READ (jmat,*)
   READ (jmat,*) a0
   READ (jmat,*)
   READ (jmat,*)
!
   WRITE (jprt,*) ' Matrixelements'
   WRITE (jprt,99001) name
   WRITE (jprt,99003) namorb
   WRITE (jprt,99004) ((kap,kapp,kapp=-lmax-1,kap),kap=-lmax-1,-1) , ((kap,kapp,kapp=-lmax-1,-1),(kap,kapp,kapp=1,kap),kap=1,lmax)
99004 FORMAT (/1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,28(i6,i3,'    ')/)
!
   WRITE (jspect,99002) title
   WRITE (jspect,99001) name
   WRITE (jspect,99010) hws , hwlev , hwc , hwlf
99010 FORMAT ('   HWS=',f5.2,' eV','    HWL=',f5.2,' eV','   HWC=',f5.2,' eV','   HWLF=',f5.2,' eV')
!
   conv = 2.0*pi/a0
   conv = conv*conv
   IF ( krel==0 ) conv = 1.0
!
   nume2 = 2*nume - 1
   necount = nume*nume
   DO i = 1 , nume
      efile(i) = e0 + de*(i-1)
   ENDDO
!
   et0 = 2.0*e0 - ecore
   DO i = 1 , nume2
      efile2(i) = et0 + de*(i-1)
      DO j = 1 , nume
         jp = i + 1 - j
         irange(j,i) = 1
         IF ( jp<1 .OR. jp>nume ) THEN
            irange(j,i) = 0
            CYCLE
         ENDIF
!
         READ (jmat,*) ii , jj , dum , eep , ee , ((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,kap),kap=-lmax-1,-1) ,                    &
                     & ((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,-1),(crsfile(jj,ii,kap,kapp),kapp=1,kap),kap=1,lmax)
!
         WRITE (jprt,99005) i , j , efile2(i) , efile(jp) , efile(j) , ((crsfile(j,i,kap,kapp),kapp=-lmax-1,kap),kap=-lmax-1,-1) , &
                          & ((crsfile(j,i,kap,kapp),kapp=-lmax-1,-1),(crsfile(j,i,kap,kapp),kapp=1,kap),kap=1,lmax)
99005    FORMAT (2I3,3F6.2,28D13.5)
!
         IF ( ii/=i .OR. jj/=j .OR. dabs(eep-efile(jp))>tiny .OR. dabs(ee-efile(j))>tiny ) THEN
            WRITE (6,*) i , j , efile(jp) , efile(j)
            WRITE (6,*) ii , jj , eep , ee
            STOP ' shit happened'
         ENDIF
!
      ENDDO
   ENDDO
!
!
!
! ***************  read in densities of states  ***************
!
!
   READ (jdos,99001) name
   READ (jdos,*)
   READ (jdos,*) ef
   ef = conv*ef
   READ (jdos,*) lmax
   READ (jdos,*) maxd
   READ (jdos,*)
!
   ii = 0
   DO i = 1 , maxd
!
      IF ( krel==0 ) THEN
         READ (jdos,*) energy , (z(l),l=0,lmax) , bla
         IF ( bla>=tiny .OR. ii/=0 ) THEN
            ii = ii + 1
            e(ii) = energy*conv
            udos(ii,-1) = z(0)
            DO l = 1 , lmax
               fl = dfloat(2*l+1)
               fl = fl + fl
               fl1 = dfloat(2*l)
               fl2 = dfloat(2*l+2)
               fl1 = fl1/fl
               fl2 = fl2/fl
               udos(ii,l) = z(l)*fl1
               udos(ii,-l-1) = z(l)*fl2
            ENDDO
         ENDIF
      ELSE
         READ (jdos,*) energy , bla , z(-1) , (z(l),z(-l-1),l=1,lmax)
         IF ( bla>=tiny .OR. ii/=0 ) THEN
            ii = ii + 1
            e(ii) = energy*conv
            DO k = -lmax - 1 , lmax
               udos(ii,k) = z(k)
            ENDDO
         ENDIF
      ENDIF
!
   ENDDO
   ne = ii
!
   WRITE (jprt,*)
   WRITE (jprt,*) ' Densities of states'
   DO i = 1 , ne
      WRITE (jprt,99006) e(i) , udos(i,-1) , (udos(i,l),udos(i,-l-1),l=1,lmax)
   ENDDO
!
!
   nef = iright(e,ne,ef)
!
! check energy ranges:
   IF ( e(1)<efile(1) .OR. e(nef)>efile(nume) ) THEN
      WRITE (6,*) ' energy range of matrix elements must exceed' , ' that of DOS!'
      STOP
   ENDIF
!
!
!  life time broadening of DOS **********************
!
   CALL dosbroad(ef,nef,e,udos,hwlr,dos,lmax)
!
   WRITE (jprt,*)
   WRITE (jprt,*) ' Densities of states after broadening'
   DO i = 1 , nef
      WRITE (jprt,99006) e(i) , dos(i,-1) , (dos(i,l),dos(i,-l-1),l=1,lmax)
   ENDDO
!
   ii = 0
   DO k = -lmax - 1 , lmax
      IF ( k/=0 ) THEN
         DO kp = -lmax - 1 , k
            IF ( kp/=0 ) THEN
               ii = ii + 1
               ikap(ii) = k
               ikapp(ii) = kp
            ENDIF
         ENDDO
      ENDIF
   ENDDO
   ktot = ii
   WRITE (6,*) ' ktot=' , ktot
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
   ne2 = iidint((e2max-e2min)/de2) + 1
   ee = e2min - de2
   DO i = 1 , ne2
      ee = ee + de2
      e2(i) = ee
   ENDDO
!
! index values for test output of matrix elements:
   ktest = ne2/6
!
!
! start loop over final energies *********************
!
   DO ie2 = 1 , ne2
      WRITE (6,*) ie2 , (e2(ie2)+e1ryd)/2.
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
      IF ( mod(ie2,ktest)==0 ) idoc = 1
      IF ( ie2==1 .OR. ie2==2 ) idoc = 1
      IF ( ie2==ne2-1 .OR. ie2==ne2 ) idoc = 1
!
      CALL matinter(crsfile,efile,nume,efile2,nume2,NEMF,NEMF2,lmax,irange,e2(ie2),emin,emax,e,NEMAX,crs,ifst,ilast,idoc,jprt)
!
! *******************************************************************
!                   put together spectrum:
! *******************************************************************
!
!
      ii = 0
      DO k = -lmax - 1 , lmax
         IF ( k/=0 ) THEN
!
            DO kp = -lmax - 1 , k
               IF ( kp/=0 ) THEN
!
                  ii = ii + 1
!
                  DO i = ifst , ilast
                     ip = ilast + ifst - i
                     y(i) = crs(i,k,kp)*dos(i,k)*dos(ip,kp)
                  ENDDO
                  p(ie2,ii) = simpson(y(ifst),nval,de)
               ENDIF
!
            ENDDO
         ENDIF
      ENDDO
!
!
! sum up diagonal contributions and total spectrum:
!
      pdiagg = 0.
      ptott = 0.
      DO ii = 1 , ktot
         IF ( ikap(ii)==ikapp(ii) ) pdiagg = pdiagg + p(ie2,ii)
         ptott = ptott + p(ie2,ii)
      ENDDO
      pdiag(ie2) = pdiagg
      ptot(ie2) = ptott
!
   ENDDO
!
!
! end of loop over final energies **********************************
!
!
! determine relative magnitudes of partial spectra:
!
   psumtot = 0.
   DO ie2 = 1 , ne2
      psumtot = psumtot + ptot(ie2)
   ENDDO
   DO ii = 1 , ktot
      psumma = 0.
      DO ie2 = 1 , ne2
         psumma = psumma + p(ie2,ii)
      ENDDO
      psum(ii) = psumma/psumtot
   ENDDO
!
! maximum of total unbroadened spectrum:
!
   pmax = 0.
   DO i = 1 , ne2
      IF ( ptot(i)>pmax ) pmax = ptot(i)
   ENDDO
!
! ordering of partial spectra according to magnitude:
!
   CALL order(psum(1),ktot,iord)
!
!
! ********************************************************************
! BROADENING DUE TO LIFE TIME OF CORE HOLE AND SPECTROMETER RESOLUTION
! ********************************************************************
! (only the most important contributions are broadened)
!
   snorm = 0.
   WRITE (6,*) ' broad'
   CALL broad(e2(1),ptot,NEMAX2,ne2,iminb,e1(1),p1tot,phelp,N1M,n1,snorm,smax,hws,hwlf,hwc,u)
   DO i = 1 , 4
      WRITE (6,*) ' broad'
      CALL broad(e2(1),p(1,iord(i)),NEMAX2,ne2,iminb,e1(1),p1(1,i),phelp,N1M,n1,snorm,dummy,hws,hwlf,hwc,u)
   ENDDO
   WRITE (6,*) ' broad'
   CALL broad(e2(1),pdiag,NEMAX2,ne2,iminb,e1(1),p1diag,phelp,N1M,n1,snorm,dummy,hws,hwlf,hwc,u)
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
   WRITE (jprt,*) ' unbroadened spectra'
   WRITE (jprt,99007) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
   WRITE (jprt,99008) (e2(i)-e2max,ptot(i),(p(i,iord(ii)),ii=1,4),pdiag(i),i=1,ne2)
!
!
! broadened spectra
!
   WRITE (jspect,99003) namorb
   WRITE (jspect,99006) ecore*efac
   WRITE (jspect,99007) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
   WRITE (jspect,99008) (e1(i)-e2maxev,p1tot(i)/snorm,(p1(i,ii)/snorm,ii=1,4),p1diag(i)/snorm,i=iminb,n1)
!
   WRITE (jspect,99009) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
99009 FORMAT (4x,'e2',5x,4x,'Total',4(4x,i2,1x,i2),4x,'Diag.')
   WRITE (jspect,99006) (e1(i)-e2maxev,p1tot(i),(p1(i,ii),ii=1,4),p1diag(i),i=iminb,n1)
!
!
   STOP ' FORTRAN stop'
!
!
99001 FORMAT (a50)
99002 FORMAT (a78)
99003 FORMAT (a10)
99006 FORMAT (f11.5,6F9.3)
99007 FORMAT (4x,'e2',5x,5x,'Total',3x,4(5x,i2,1x,i2,3x),5x,'Diag.')
99008 FORMAT (f11.5,6E13.5)
END PROGRAM rcvvaug
!*==DOSBROAD.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dosbroad(Ef,Nef,E,Udos,Hwlr,Dos,Lmax)
! ====================================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a , convlo , Dos , E , Ef , estep , gam , Hwlr , hwlr0 , pi , Udos
   INTEGER i , ilor , j , k , KM , Lmax , Nef , NEMAX
!*** End of declarations inserted by SPAG
!
! life-time broadening of valence band states with Lorentzian
! only states below Efermi enter the integral
! DOS is required on an equidistant energy mesh!
!
   PARAMETER (NEMAX=500)
   PARAMETER (KM=3)
!
   DIMENSION E(NEMAX) , Udos(NEMAX,-KM-1:KM) , Dos(NEMAX,-KM-1:KM)
   DIMENSION a(NEMAX)
!
   DATA pi/3.141596/ , ilor/2/ , hwlr0/0.1/
!
   WRITE (6,*) ' I am in dosbroad'
   IF ( Hwlr<hwlr0 ) THEN
      DO k = -Lmax - 1 , Lmax
         IF ( k/=0 ) THEN
            DO j = 1 , Nef
               Dos(j,k) = Udos(j,k)
            ENDDO
         ENDIF
      ENDDO
      WRITE (6,*) ' I am out of dosbroad'
      RETURN
   ENDIF
!
   estep = E(2) - E(1)
   DO k = -Lmax - 1 , Lmax
      IF ( k/=0 ) THEN
         DO i = 1 , Nef
            DO j = 1 , Nef
               a(j) = Udos(j,k)
            ENDDO
! halfwidth of Lor. increases linearly with binding E.:
            IF ( ilor==1 ) gam = Hwlr*(Ef-E(i))/(Ef-E(1))/2.0E0
! halfwidth of Lor. increases quadr.:
            IF ( ilor==2 ) gam = Hwlr*((Ef-E(i))/(Ef-E(1)))**2/2.0E0
            Dos(i,k) = convlo(gam,a(1),i,estep,1,Nef,pi)
         ENDDO
      ENDIF
   ENDDO
!
   WRITE (6,*) ' Out of dosbroad'
END SUBROUTINE dosbroad
!*==BROAD.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE broad(Ein,A0,N0,Max0,Mine1,E1,A1,Ai,N1m,N1,Snorm,Smax,Hws,Hwl,Hwc,U)
! ==============================================================
!
!     interpolate to plot-grid and convolute with lorentzian
!     for life-time broadening of core hole + valence band holes
!     and gaussian spectrometer resolution
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A0 , A1 , afac , ah , Ai , bfac , cfac , convgau , convlo , de0 , dele , e0 , e00 , e0max , E1 , Ein , emin0 , estep ,   &
        & evmax , evmin
   REAL*8 gam , gc , Hwc , hwc0 , Hwl , Hws , pi , Smax , Snorm , U
   INTEGER i , i0 , ix0 , j , max , Max0 , Mine1 , N0 , N1 , N1m , NGC , nlst , NMAX
!*** End of declarations inserted by SPAG
   PARAMETER (NGC=100,NMAX=500)
   DIMENSION Ein(N0) , A0(N0) , e0(NMAX) , E1(N1m) , A1(N1m) , Ai(N1m)
   DIMENSION gc(-NGC:NGC) , ah(-NGC:NGC)
   INTEGER :: spag_nextblock_1
   DATA hwc0/5.E-03/ , cfac/13.606/
   DATA estep/0.05E0/
   DATA ix0/201/ , pi/3.141596/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!     convert from rydberg to ev:
!
         WRITE (6,*) ' I am in broad'
         DO i = 1 , Max0
            e0(i) = Ein(i)*cfac - U
         ENDDO
         de0 = e0(2) - e0(1)
         e0max = e0(Max0)
         e00 = e0max
         DO i = Max0 + 1 , NMAX
            e00 = e00 + de0
            e0(i) = e00
            A0(i) = 0.0
            IF ( e00-e0max>U ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
         max = i
!udo  this next line is not in the orginal program
         nlst = N1
!
         IF ( Snorm<=0.E0 ) THEN
!
            evmin = e0(1) - 1.
            evmax = e0(max) + 1.
            SPAG_Loop_2_2: DO
               N1 = idnint((evmax-evmin)/estep) + 1
               IF ( N1>N1m ) THEN
                  estep = estep + estep
                  CYCLE
               ENDIF
! create nice numbers for energy:
               i0 = idnint(evmin/estep)
               evmin = i0*estep
!
!     set up plot-grid and determine min-energy
!
               nlst = N1
               WRITE (6,*) ' I know what nlst is:' , nlst
               DO i = 1 , N1
                  E1(i) = evmin + (i-1)*estep
                  IF ( E1(i)>e0(max) .AND. nlst==N1 ) nlst = i
               ENDDO
               SPAG_Loop_3_1: DO i = 1 , N0
                  IF ( A0(i)>0.E0 ) EXIT SPAG_Loop_3_1
               ENDDO SPAG_Loop_3_1
               EXIT SPAG_Loop_2_2
            ENDDO SPAG_Loop_2_2
            emin0 = e0(i)
            SPAG_Loop_2_3: DO i = 1 , N1
               IF ( E1(i)>emin0 ) EXIT SPAG_Loop_2_3
            ENDDO SPAG_Loop_2_3
            Mine1 = i
         ENDIF
!
!
!     interpolate
!
         DO i = Mine1 , nlst
            j = 0
            SPAG_Loop_3_4: DO
               j = j + 1
               IF ( e0(j)>=E1(i) .OR. j>=max ) THEN
                  j = j - 2
                  IF ( j<1 ) j = 1
                  IF ( j>max-3 ) j = max - 3
                  CALL inter(e0(j),A0(j),4,E1(i),Ai(i))
                  EXIT SPAG_Loop_3_4
               ENDIF
            ENDDO SPAG_Loop_3_4
         ENDDO
         DO i = 1 , Mine1 - 1
            Ai(i) = 0.
         ENDDO
         DO i = nlst + 1 , N1
            Ai(i) = 0.
         ENDDO
!
         DO i = Mine1 , N1
!
!     determine halfwidth for life-time broadening
!     (2 final valence holes + 1 initial core hole)
!
            dele = (E1(nlst)-E1(i))/(E1(nlst)-E1(Mine1))
            gam = (Hwl*dele*dele+Hwc)/2.0E0
!
!     convolute now
!
            A1(i) = convlo(gam,Ai,i,estep,Mine1,N1,pi)
         ENDDO
!
!     overwrite unbroadened data in ai with life-time
!     broadened spectrum
!
         DO i = 1 , N1
            Ai(i) = A1(i)
         ENDDO
         IF ( Hws/=0. ) THEN
!
!        set up gaussian convolution function
!
            afac = -alog(5.E-1)/(Hws/2.0E0)**2
            bfac = sqrt(afac/pi)
            DO i = -NGC , NGC
               gc(i) = exp(-afac*(i*estep)**2)*bfac
            ENDDO
!
!        convolute now
!
            DO i = Mine1 , N1
               A1(i) = convgau(Ai,i,gc,ah,estep,NGC,N1)
            ENDDO
         ENDIF
!
         IF ( Snorm<=0.E0 ) THEN
!
!     determine maximum and set equal to 100.
!
            Smax = 0.E0
            DO i = Mine1 , N1
               IF ( A1(i)>Smax ) Smax = A1(i)
            ENDDO
            Snorm = 100.D0/Smax
         ENDIF
!
!
         DO i = Mine1 , N1
            A1(i) = A1(i)*Snorm
!
         ENDDO
!
         WRITE (6,*) ' Out of broad'
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE broad
!*==CONVLO.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION convlo(Gam,A1,Ind,Estep,Mine1,N1,Pi)
! ==================================================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A1 , atm , atp , Estep , Gam , Pi , sum
   INTEGER i , Ind , Mine1 , N1
!*** End of declarations inserted by SPAG
!
!     do the convolution integral by trapezoidal rule
!
   DIMENSION A1(N1)
   WRITE (6,*) ' Out of convlo'
   sum = 0.0
   DO i = Mine1 , N1
      atp = atan(Estep*((i-Ind)+0.5)/Gam)
      atm = atan(Estep*((i-Ind)-0.5)/Gam)
      sum = sum + A1(i)*(atp-atm)
   ENDDO
   convlo = sum/Pi
   WRITE (6,*) ' Out of convlo'
END FUNCTION convlo
!*==CONVGAU.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION convgau(A1,Ind,Gc,Ah,Estep,Ngc,N1)
! ================================================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A1 , Ah , Estep , Gc , sum
   INTEGER i , imx , Ind , N1 , Ngc
!*** End of declarations inserted by SPAG
!
!     do the convolution integral by trapezoidal rule
!
   DIMENSION A1(N1) , Gc(-Ngc:Ngc) , Ah(-Ngc:Ngc)
!
   WRITE (6,*) ' In convgau'
   DO i = -Ngc , Ngc
      Ah(i) = 0.D0
   ENDDO
   DO i = 0 , Ngc
      imx = Ind + i
      IF ( imx>N1 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      Ah(i) = A1(imx)*Gc(i)
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      DO i = -1 , -Ngc , -1
         imx = Ind + i
         IF ( imx<1 ) THEN
            CALL spag_block_2
            RETURN
         ENDIF
         Ah(i) = A1(imx)*Gc(i)
      ENDDO
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
      sum = 0.D0
      DO i = -Ngc , Ngc
         sum = sum + Ah(i)
      ENDDO
!
      convgau = sum*Estep
!
      WRITE (6,*) ' Out of convgau'
   END SUBROUTINE spag_block_2
END FUNCTION convgau
!*==MATINTER.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE matinter(Crsfile,Efile,Nume,Efile2,Nume2,Nemf,Nemf2,Km,Irange,E2,Emin,Emax,E,Nemax,Crs,Ifst,Ilast,Idoc,Jprt)
! ========================================================
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 chelp , Crs , crs1 , Crsfile , E , E2 , ebar , Efile , Efile2 , ehelp , Emax , Emin
   INTEGER i , i2 , Idoc , ie , iextra , if , if2l , if2r , ifl , ifr , Ifst , Ilast , insleft , insright , Irange , iright ,      &
         & Jprt , k , Km , KML
   INTEGER kp , l , leftind , lp , Nemax , Nemf , Nemf2 , nf , nint , NM , Nume , Nume2
!*** End of declarations inserted by SPAG
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
   PARAMETER (KML=3,NM=50)
!
   DIMENSION Crsfile(Nemf,Nemf2,-KML-1:KML,-KML-1:KML) , Efile(Nemf) , Efile2(Nemf2) , Irange(Nemf,Nemf2) , E(Nemax) ,             &
           & Crs(Nemax,-KML-1:KML,-KML-1:KML)
   DIMENSION ehelp(NM) , chelp(NM) , crs1(NM,-KML-1:KML,-KML-1:KML)
   WRITE (6,*) ' In matinter'
!
! find efile2 points immediately left and right
! of the actual e2 (if2l,if2r):
   if2l = leftind(Efile2,Nume2,E2)
   if2r = iright(Efile2,Nume2,E2)
!
! valence band energy range: find efile points immediately
! left of emin and right of emax (ifl,ifr):
   ifl = leftind(Efile,Nume,Emin)
   ifr = iright(Efile,Nume,Emax)
!
! for interpolation lateron, we may have to calculate the boundary
! values of the crs in a different way (to understand the whole
! thing, draw a diagram e vs. e2)
   insleft = 0
   insright = 0
   IF ( Irange(ifl,if2r)==0 ) insleft = 1
   IF ( Irange(ifr,if2l)==0 ) insright = 1
!
! now calculate the crs's for e2 and all values of efile
! between ifl and ifr and store them in crs1:
!
! boundary values (interpolation along boundaries of e/e2 range):
   nint = 5
   IF ( Nume<nint ) nint = Nume
   IF ( insleft==1 ) THEN
      ebar = (Efile2(Nume2)+Efile2(1))/2.
      ehelp(ifl) = E2 - ebar + Efile(1)
!        if(ehelp(ifl).ge.emin)then
!          ifl=ifl+1
!          goto 295
!        endif
      DO l = -Km - 1 , Km
         IF ( l/=0 ) THEN
            DO lp = -Km - 1 , l
               IF ( lp/=0 ) THEN
                  DO i = Nume , 1 , -1
                     i2 = Nume2 - Nume + i
                     chelp(i) = Crsfile(i,i2,l,lp)
                  ENDDO
                  CALL lagrange(Efile(1),chelp(1),Nume,nint,ehelp(ifl),crs1(ifl,l,lp),iextra)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IF ( Idoc==1 ) WRITE (Jprt,99001)
99001 FORMAT (' non-grid value for lower e boundary')
   ENDIF
   IF ( insright==1 ) THEN
      ehelp(ifr) = E2 - Efile2(1) + Efile(1)
!        if(ehelp(ifr).le.emax)then
!          ifl=ifl-1
!          goto 297
!        endif
      DO l = -Km - 1 , Km
         IF ( l/=0 ) THEN
            DO lp = -Km - 1 , l
               IF ( lp/=0 ) THEN
                  DO i = 1 , Nume
                     chelp(i) = Crsfile(i,i,l,lp)
                  ENDDO
                  CALL lagrange(Efile(1),chelp(1),Nume,nint,ehelp(ifr),crs1(ifr,l,lp),iextra)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IF ( Idoc==1 ) WRITE (Jprt,99002)
99002 FORMAT (' non-grid value for upper e boundary')
   ENDIF
   IF ( Idoc==1 ) WRITE (Jprt,99003) ifl , ifr
99003 FORMAT (' ifl,ifr: ',2I4)
!
! get the rest of crs1:
   DO if = ifl , ifr
      IF ( if/=ifl .OR. insleft/=1 ) THEN
         IF ( if/=ifr .OR. insright/=1 ) THEN
            ehelp(if) = Efile(if)
            DO l = -Km - 1 , Km
               IF ( l/=0 ) THEN
                  DO lp = -Km - 1 , l
                     IF ( lp/=0 ) THEN
                        DO i = if , if + Nume - 1
                           chelp(i) = Crsfile(if,i,l,lp)
                        ENDDO
                        CALL lagrange(Efile2(if),chelp(if),Nume,nint,E2,crs1(if,l,lp),iextra)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDIF
   ENDDO
!
   IF ( Idoc==1 ) THEN
      WRITE (Jprt,99004)
99004 FORMAT (' crs1:')
      DO i = ifl , ifr
         WRITE (Jprt,99005) i , ehelp(i) , ((crs1(i,l,lp),lp=-Km-1,l),l=-Km-1,-1) ,                                                &
                          & ((crs1(i,l,lp),lp=-Km-1,-1),(crs1(i,l,lp),lp=1,l),l=1,Km)
99005    FORMAT (i3,f10.5,28E13.5)
      ENDDO
   ENDIF
 
!
! now we have a function crs1 on a mesh efile and we want a function
! crs on a mesh e and that's all
   nf = ifr - ifl + 1
   nint = 6
   IF ( nint>nf ) nint = nf
   DO ie = Ifst , Ilast
      DO l = -Km - 1 , Km
         IF ( l/=0 ) THEN
            DO lp = -Km - 1 , l
               IF ( lp/=0 ) THEN
                  DO i = ifl , ifr
                     chelp(i) = crs1(i,l,lp)
                  ENDDO
                  CALL lagrange(Efile(ifl),chelp(ifl),nf,nint,E(ie),Crs(ie,l,lp),iextra)
                  IF ( iextra==1 .AND. Idoc==1 ) WRITE (Jprt,99006) ie , E(ie) , ((Crs(ie,k,kp),kp=-Km-1,k),k=-Km-1,-1) ,          &
                     & ((Crs(ie,k,kp),kp=-Km-1,-1),(Crs(ie,k,kp),kp=1,k),k=1,Km)
99006             FORMAT (' extrapolation: crs, ie=',i3,' e(ie)=',f10.5/28E13.5)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
!
! test output:
   IF ( Idoc==1 ) THEN
      WRITE (Jprt,99007) E2
99007 FORMAT (' e2: ',f10.5)
      DO ie = Ifst , Ilast
         WRITE (Jprt,99008) ie , E(ie) , ((Crs(ie,l,lp),lp=-Km-1,l),l=-Km-1,-1) , ((Crs(ie,l,lp),lp=-Km-1,-1),(Crs(ie,l,lp),lp=1,l)&
                          & ,l=1,Km)
99008    FORMAT (i3,f10.5,28E13.4)
      ENDDO
      WRITE (Jprt,99009)
99009 FORMAT (120('-'))
   ENDIF
   WRITE (6,*) ' out of matinter'
END SUBROUTINE matinter
!*==LAGRANGE.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE lagrange(X,Y,Ndim,Npkt,Xval,Yinter,Iextra)
! =========================================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   INTEGER Iextra , ifst , igr , ileft , Ndim , Npkt
   REAL*8 X , Xval , Y , Yinter
!*** End of declarations inserted by SPAG
   DIMENSION X(Ndim) , Y(Ndim)
!
! Interpolation mit Komfort
! iextra = 0, wenn keine Extrapolation
!          1, wenn Extrapolation nach oben
!         -1, wenn Extrapolation nach unten
!
   WRITE (6,*) ' in lagrange'
   Iextra = 0
   IF ( Xval<X(1) ) Iextra = -1
   IF ( Xval>X(Ndim) ) Iextra = 1
!
! Suche 1. Punkt rechts von xval, falls es ihn gibt:
   igr = 0
   SPAG_Loop_1_1: DO
      igr = igr + 1
      IF ( X(igr)==Xval ) THEN
         Yinter = Y(igr)
         RETURN
      ENDIF
      IF ( X(igr)>=Xval .OR. igr>=Ndim ) THEN
!
! Bestimme 1. Punkt f. Interpolation
         ileft = Npkt/2
         IF ( mod(Npkt,2)/=0 ) ileft = ileft + 1
         ifst = igr - ileft
!
         IF ( ifst<=0 ) ifst = 1
         IF ( ifst>Ndim-Npkt+1 ) ifst = Ndim - Npkt + 1
!
! Lasst uns interpolieren:
         CALL inter(X(ifst),Y(ifst),Npkt,Xval,Yinter)
! Und das war's dann auch schon, Leute!
         WRITE (6,*) ' Out of Lagrange'
         EXIT SPAG_Loop_1_1
      ENDIF
   ENDDO SPAG_Loop_1_1
END SUBROUTINE lagrange
!*==INTER.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE inter(R,P,N,Rs,Ps)
! =================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 denom , P , Ps , R , Rs , term
   INTEGER i , j , N
!*** End of declarations inserted by SPAG
!
!     interpolate via lagrange
!
   DIMENSION R(N) , P(N)
   WRITE (6,*) ' In inter'
   Ps = 0.E0
   DO j = 1 , N
      term = 1.E0
      denom = 1.E0
      DO i = 1 , N
         IF ( i/=j ) THEN
            denom = denom*(R(j)-R(i))
            term = term*(Rs-R(i))
         ENDIF
      ENDDO
      Ps = Ps + term*P(j)/denom
   ENDDO
   WRITE (6,*) ' Out of inter'
END SUBROUTINE inter
!*==LEFTIND.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
INTEGER FUNCTION leftind(Field,Ndim,Value)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   INTEGER i , Ndim
!*** End of declarations inserted by SPAG
! ==============================================
! index of the field element immediately left of value
! if value is lower than field(1), 1 is returned
!      implict real*8 (a-h,o-z)
   REAL*8 Field , Value
   DIMENSION Field(Ndim)
   DO i = Ndim , 1 , -1
      IF ( Field(i)<Value ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      leftind = i
      IF ( leftind<1 ) leftind = 1
   END SUBROUTINE spag_block_1
END FUNCTION leftind
!*==IRIGHT.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
INTEGER FUNCTION iright(Field,Ndim,Value)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   INTEGER i , Ndim
!*** End of declarations inserted by SPAG
! =============================================
! index of the field element immediately right of value
! if value is greater than fieldhndim), ndim is returned
!      implict real*8 (a-h,o-z)
   REAL*8 Field , Value
   DIMENSION Field(Ndim)
   DO i = 1 , Ndim
      IF ( Field(i)>Value ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      iright = i
      IF ( iright>Ndim ) iright = Ndim
   END SUBROUTINE spag_block_1
END FUNCTION iright
!*==SIMPSON.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION simpson(Fct,N,Dx)
! ===============================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Dx , Fct , sum , x0 , x1 , x2
   INTEGER i , N
!*** End of declarations inserted by SPAG
! integration of fct using simpson integration
! equidistant point mesh
   DIMENSION Fct(N)
   sum = 0.
   x0 = Fct(1)
!
   DO i = 3 , N , 2
      x1 = Fct(i-1)
      x2 = Fct(i)
      sum = sum + Dx*(x0+4*x1+x2)
      x0 = x2
   ENDDO
   simpson = sum/3.
   IF ( mod(N,2)==0 ) simpson = simpson + (Fct(N)+Fct(N-1))*Dx/2.
END FUNCTION simpson
!*==ORDER.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE order(F,N,Ior)
! =============================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 F , g
   INTEGER i , Ior , j , l , N
!*** End of declarations inserted by SPAG
   DIMENSION F(N) , Ior(N)
!
! f is to be arranged in decreasing order
   DO i = 1 , N
      Ior(i) = i
   ENDDO
!
   DO i = 1 , N - 1
      DO j = i + 1 , N
         IF ( F(j)>F(i) ) THEN
            g = F(j)
            F(j) = F(i)
            F(i) = g
            l = Ior(j)
            Ior(j) = Ior(i)
            Ior(i) = l
         ENDIF
      ENDDO
   ENDDO
END SUBROUTINE order
