!*==RCVVAUG.spg  processed by SPAG 6.72Dc https://fortran.uk/plusfortonline.php at 02:03 on 30 Oct 2023
      program RCVVAUG
!
! **********************************************************
! CALCULATE RELATIVISTIC CORE-VALENCE-VALENCE AUGER SPECTRA
! **********************************************************
!
!
      implicit none
!*--RCVVAUG13
!*** Start of declarations inserted by SPAG
      real*8 a0 , bla , conv , crs , crsfile , de , de2 , dos , dum ,   &
           & dummy , e , e0 , e1 , e1ryd , e2 , e2max , e2maxev ,       &
           & e2min , ecore , ee
      real*8 eep , ef , efac , efile , efile2 , emax , emin , energy ,  &
           & et0 , fl , fl1 , fl2 , hwc , hwlev , hwlf , hwlr , hws ,   &
           & p , p1 , p1diag
      real*8 p1tot , pdiag , pdiagg , phelp , pi , pmax , psum ,        &
           & psumma , psumtot , ptot , ptott , SIMPSON , smax , snorm , &
           & tiny , u , udos , y , z
      integer i , idoc , ie2 , ifst , ii , ikap , ikapp , ilast ,       &
            & iminb , in , iord , ip , irange , IRIGHT , j , jdos , jj ,&
            & jmat , jp , jprt
      integer jspect , k , kap , kapp , KM , KMTOT , kp , krel , ktest ,&
            & ktot , l , LEFTIND , lmax , maxd , n1 , N1M , ne , ne2 ,  &
            & necount , nef
      integer NEMAX , NEMAX2 , NEMF , NEMF2 , nume , nume2 , nval
!*** End of declarations inserted by SPAG
!
      parameter (NEMAX=500,NEMAX2=1000,N1M=300,NEMF=15,NEMF2=30)
      parameter (KM=3,KMTOT=28)
!
      character name*50 , title*78 , namorb*10 , filen*50
!
      dimension e1(N1M) , e(NEMAX) , efile(NEMF) , efile2(NEMF2) ,      &
              & e2(NEMAX2) , irange(NEMF,NEMF2)
      dimension z(-KM-1:KM) , udos(NEMAX,-KM-1:KM) , dos(NEMAX,-KM-1:KM)
      dimension crsfile(NEMF,NEMF2,-KM-1:KM,-KM-1:KM) ,                 &
              & crs(NEMAX,-KM-1:KM,-KM-1:KM)
      dimension p(NEMAX2,KMTOT) , psum(KMTOT)
      dimension pdiag(NEMAX2) , ptot(NEMAX2)
      dimension iord(KMTOT) , ikap(KMTOT) , ikapp(KMTOT)
      dimension y(NEMAX) , phelp(N1M)
      dimension p1(N1M,4) , p1tot(N1M) , p1diag(N1M)
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
      read (in,99005) filen
      write (*,99005) filen
      open (unit=jdos,file=filen)
!
!
! matrix element file
!
      read (in,99005) filen
      write (*,99005) filen
      open (unit=jmat,file=filen)
!
!
! print file
!
      read (in,99005) filen
      write (*,99005) filen
      open (unit=jprt,file=filen)
!
!
! spectrum file
!
      read (in,99005) filen
      write (*,99005) filen
      open (unit=jspect,file=filen)
!
! ******************************
!
!
      read (in,99006) title
      write (*,99006) title
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
      read (jmat,99005) name
      read (jmat,*)
      read (jmat,99007) namorb
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
      write (jprt,99005) name
      write (jprt,99007) namorb
      write (jprt,99001) ((kap,kapp,kapp=-lmax-1,kap),kap=-lmax-1,-1) , &
                       & ((kap,kapp,kapp=-lmax-1,-1),                   &
                       & (kap,kapp,kapp=1,kap),kap=1,lmax)
99001 format (/1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,28(i6,i3,'    ')&
            & /)
!
      write (jspect,99006) title
      write (jspect,99005) name
      write (jspect,99002) hws , hwlev , hwc , hwlf
99002 format ('   HWS=',f5.2,' eV','    HWL=',f5.2,' eV','   HWC=',f5.2,&
             &' eV','   HWLF=',f5.2,' eV')
!
      conv = 2.0*pi/a0
      conv = conv*conv
      if ( krel.EQ.0 ) conv = 1.0
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
            if ( jp.LT.1 .OR. jp.GT.nume ) then
               irange(j,i) = 0
               goto 50
            endif
!
            read (jmat,*) ii , jj , dum , eep , ee ,                    &
                        & ((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,kap),  &
                        & kap=-lmax-1,-1) ,                             &
                        & ((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,-1),   &
                        & (crsfile(jj,ii,kap,kapp),kapp=1,kap),kap=1,   &
                        & lmax)
!
            write (jprt,99003) i , j , efile2(i) , efile(jp) , efile(j) &
                             & ,                                        &
                             & ((crsfile(j,i,kap,kapp),kapp=-lmax-1,kap)&
                             & ,kap=-lmax-1,-1) ,                       &
                             & ((crsfile(j,i,kap,kapp),kapp=-lmax-1,-1),&
                             & (crsfile(j,i,kap,kapp),kapp=1,kap),kap=1,&
                             & lmax)
99003       format (2I3,3F6.2,28D13.5)
!
            if ( ii.NE.i .OR. jj.NE.j .OR. dabs(eep-efile(jp))          &
               & .GT.tiny .OR. dabs(ee-efile(j)).GT.tiny ) then
               write (6,*) i , j , efile(jp) , efile(j)
               write (6,*) ii , jj , eep , ee
               stop ' shit happened'
            endif
!
 50      enddo
      enddo
!
!
!
! ***************  read in densities of states  ***************
!
!
      read (jdos,99005) name
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
         if ( krel.EQ.0 ) then
            read (jdos,*) energy , (z(l),l=0,lmax) , bla
            if ( bla.GE.tiny .OR. ii.NE.0 ) then
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
            if ( bla.GE.tiny .OR. ii.NE.0 ) then
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
         write (jprt,99008) e(i) , udos(i,-1) ,                         &
                          & (udos(i,l),udos(i,-l-1),l=1,lmax)
      enddo
!
!
      nef = IRIGHT(e,ne,ef)
!
! check energy ranges:
      if ( e(1).LT.efile(1) .OR. e(nef).GT.efile(nume) ) then
         write (6,*) ' energy range of matrix elements must exceed' ,   &
                    &' that of DOS!'
         stop
      endif
!
!
!  life time broadening of DOS **********************
!
      call DOSBROAD(ef,nef,e,udos,hwlr,dos,lmax)
!
      write (jprt,*)
      write (jprt,*) ' Densities of states after broadening'
      do i = 1 , nef
         write (jprt,99008) e(i) , dos(i,-1) ,                          &
                          & (dos(i,l),dos(i,-l-1),l=1,lmax)
      enddo
!
      ii = 0
      do k = -lmax - 1 , lmax
         if ( k.NE.0 ) then
            do kp = -lmax - 1 , k
               if ( kp.NE.0 ) then
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
      ne2 = iidint((e2max-e2min)/de2) + 1
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
         ilast = IRIGHT(e,nef,emax)
         ifst = LEFTIND(e,nef,emin)
         nval = ilast - ifst + 1
!
! ********************************************************************
! interpolate matrix elements to DOS energy grid for this final energy
! ********************************************************************
!
! documentation:
         idoc = 0
         if ( mod(ie2,ktest).EQ.0 ) idoc = 1
         if ( ie2.EQ.1 .OR. ie2.EQ.2 ) idoc = 1
         if ( ie2.EQ.ne2-1 .OR. ie2.EQ.ne2 ) idoc = 1
!
         call MATINTER(crsfile,efile,nume,efile2,nume2,NEMF,NEMF2,lmax, &
                     & irange,e2(ie2),emin,emax,e,NEMAX,crs,ifst,ilast, &
                     & idoc,jprt)
!
! *******************************************************************
!                   put together spectrum:
! *******************************************************************
!
!
         ii = 0
         do k = -lmax - 1 , lmax
            if ( k.NE.0 ) then
!
               do kp = -lmax - 1 , k
                  if ( kp.NE.0 ) then
!
                     ii = ii + 1
!
                     do i = ifst , ilast
                        ip = ilast + ifst - i
                        y(i) = crs(i,k,kp)*dos(i,k)*dos(ip,kp)
                     enddo
                     p(ie2,ii) = SIMPSON(y(ifst),nval,de)
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
            if ( ikap(ii).EQ.ikapp(ii) ) pdiagg = pdiagg + p(ie2,ii)
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
         if ( ptot(i).GT.pmax ) pmax = ptot(i)
      enddo
!
! ordering of partial spectra according to magnitude:
!
      call ORDER(psum(1),ktot,iord)
!
!
! ********************************************************************
! BROADENING DUE TO LIFE TIME OF CORE HOLE AND SPECTROMETER RESOLUTION
! ********************************************************************
! (only the most important contributions are broadened)
!
      snorm = 0.
      write (6,*) ' broad'
      call BROAD(e2(1),ptot,NEMAX2,ne2,iminb,e1(1),p1tot,phelp,N1M,n1,  &
               & snorm,smax,hws,hwlf,hwc,u)
      do i = 1 , 4
         write (6,*) ' broad'
         call BROAD(e2(1),p(1,iord(i)),NEMAX2,ne2,iminb,e1(1),p1(1,i),  &
                  & phelp,N1M,n1,snorm,dummy,hws,hwlf,hwc,u)
      enddo
      write (6,*) ' broad'
      call BROAD(e2(1),pdiag,NEMAX2,ne2,iminb,e1(1),p1diag,phelp,N1M,n1,&
               & snorm,dummy,hws,hwlf,hwc,u)
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
      write (jprt,99009) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
      write (jprt,99010) (e2(i)-e2max,ptot(i),(p(i,iord(ii)),ii=1,4),   &
                       & pdiag(i),i=1,ne2)
!
!
! broadened spectra
!
      write (jspect,99007) namorb
      write (jspect,99008) ecore*efac
      write (jspect,99009) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
      write (jspect,99010) (e1(i)-e2maxev,p1tot(i)/snorm,(p1(i,ii)/snorm&
                         & ,ii=1,4),p1diag(i)/snorm,i=iminb,n1)
!
      write (jspect,99004) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
99004 format (4x,'e2',5x,4x,'Total',4(4x,i2,1x,i2),4x,'Diag.')
      write (jspect,99008) (e1(i)-e2maxev,p1tot(i),(p1(i,ii),ii=1,4),   &
                         & p1diag(i),i=iminb,n1)
!
!
      stop ' FORTRAN stop'
!
!
99005 format (a50)
99006 format (a78)
99007 format (a10)
99008 format (f11.5,6F9.3)
99009 format (4x,'e2',5x,5x,'Total',3x,4(5x,i2,1x,i2,3x),5x,'Diag.')
99010 format (f11.5,6E13.5)
      end
!*==DOSBROAD.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      subroutine DOSBROAD(Ef,Nef,E,Udos,Hwlr,Dos,Lmax)
! ====================================================
      implicit none
!*--DOSBROAD480
!*** Start of declarations inserted by SPAG
      real*8 a , CONVLO , Dos , E , Ef , estep , gam , Hwlr , hwlr0 ,   &
           & pi , Udos
      integer i , ilor , j , k , KM , Lmax , Nef , NEMAX
!*** End of declarations inserted by SPAG
!
! life-time broadening of valence band states with Lorentzian
! only states below Efermi enter the integral
! DOS is required on an equidistant energy mesh!
!
      parameter (NEMAX=500)
      parameter (KM=3)
!
      dimension E(NEMAX) , Udos(NEMAX,-KM-1:KM) , Dos(NEMAX,-KM-1:KM)
      dimension a(NEMAX)
!
      data pi/3.141596/ , ilor/2/ , hwlr0/0.1/
!
      write (6,*) ' I am in dosbroad'
      if ( Hwlr.LT.hwlr0 ) then
         do k = -Lmax - 1 , Lmax
            if ( k.NE.0 ) then
               do j = 1 , Nef
                  Dos(j,k) = Udos(j,k)
               enddo
            endif
         enddo
         write (6,*) ' I am out of dosbroad'
         return
      endif
!
      estep = E(2) - E(1)
      do k = -Lmax - 1 , Lmax
         if ( k.NE.0 ) then
            do i = 1 , Nef
               do j = 1 , Nef
                  a(j) = Udos(j,k)
               enddo
! halfwidth of Lor. increases linearly with binding E.:
               if ( ilor.EQ.1 ) gam = Hwlr*(Ef-E(i))/(Ef-E(1))/2.0E0
! halfwidth of Lor. increases quadr.:
               if ( ilor.EQ.2 ) gam = Hwlr*((Ef-E(i))/(Ef-E(1)))        &
                                    & **2/2.0E0
               Dos(i,k) = CONVLO(gam,a(1),i,estep,1,Nef,pi)
            enddo
         endif
      enddo
!
      write (6,*) ' Out of dosbroad'
      end
!*==BROAD.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      subroutine BROAD(Ein,A0,N0,max0,Mine1,E1,A1,Ai,N1m,N1,Snorm,Smax, &
                     & Hws,Hwl,Hwc,U)
! ==============================================================
!
!     interpolate to plot-grid and convolute with lorentzian
!     for life-time broadening of core hole + valence band holes
!     and gaussian spectrometer resolution
!
      implicit none
!*--BROAD541
!*** Start of declarations inserted by SPAG
      real*8 A0 , A1 , afac , ah , Ai , bfac , cfac , CONVGAU , CONVLO ,&
           & de0 , dele , e0 , e00 , e0max , E1 , Ein , emin0 , estep , &
           & evmax , evmin
      real*8 gam , gc , Hwc , hwc0 , Hwl , Hws , pi , Smax , Snorm , U
      integer i , i0 , ix0 , j , max , max0 , Mine1 , N0 , N1 , N1m ,   &
            & NGC , nlst , NMAX
!*** End of declarations inserted by SPAG
      parameter (NGC=100,NMAX=500)
      dimension Ein(N0) , A0(N0) , e0(NMAX) , E1(N1m) , A1(N1m) ,       &
              & Ai(N1m)
      dimension gc(-NGC:NGC) , ah(-NGC:NGC)
      data hwc0/5.E-03/ , cfac/13.606/
      data estep/0.05E0/
      data ix0/201/ , pi/3.141596/
!
!     convert from rydberg to ev:
!
      write (6,*) ' I am in broad'
      do i = 1 , max0
         e0(i) = Ein(i)*cfac - U
      enddo
      de0 = e0(2) - e0(1)
      e0max = e0(max0)
      e00 = e0max
      do i = max0 + 1 , NMAX
         e00 = e00 + de0
         e0(i) = e00
         A0(i) = 0.0
         if ( e00-e0max.GT.U ) goto 100
      enddo
 100  max = i
!udo  this next line is not in the orginal program
      nlst = N1
!
      if ( Snorm.LE.0.E0 ) then
!
         evmin = e0(1) - 1.
         evmax = e0(max) + 1.
 150     N1 = idnint((evmax-evmin)/estep) + 1
         if ( N1.GT.N1m ) then
            estep = estep + estep
            goto 150
         endif
! create nice numbers for energy:
         i0 = idnint(evmin/estep)
         evmin = i0*estep
!
!     set up plot-grid and determine min-energy
!
         nlst = N1
         write (6,*) ' I know what nlst is:' , nlst
         do i = 1 , N1
            E1(i) = evmin + (i-1)*estep
            if ( E1(i).GT.e0(max) .AND. nlst.EQ.N1 ) nlst = i
         enddo
         do i = 1 , N0
            if ( A0(i).GT.0.E0 ) goto 200
         enddo
 200     emin0 = e0(i)
         do i = 1 , N1
            if ( E1(i).GT.emin0 ) goto 250
         enddo
 250     Mine1 = i
      endif
!
!
!     interpolate
!
      do i = Mine1 , nlst
         j = 0
 300     j = j + 1
         if ( e0(j).LT.E1(i) .AND. j.LT.max ) goto 300
         j = j - 2
         if ( j.LT.1 ) j = 1
         if ( j.GT.max-3 ) j = max - 3
         call INTER(e0(j),A0(j),4,E1(i),Ai(i))
      enddo
      do i = 1 , Mine1 - 1
         Ai(i) = 0.
      enddo
      do i = nlst + 1 , N1
         Ai(i) = 0.
      enddo
!
      do i = Mine1 , N1
!
!     determine halfwidth for life-time broadening
!     (2 final valence holes + 1 initial core hole)
!
         dele = (E1(nlst)-E1(i))/(E1(nlst)-E1(Mine1))
         gam = (Hwl*dele*dele+Hwc)/2.0E0
!
!     convolute now
!
         A1(i) = CONVLO(gam,Ai,i,estep,Mine1,N1,pi)
      enddo
!
!     overwrite unbroadened data in ai with life-time
!     broadened spectrum
!
      do i = 1 , N1
         Ai(i) = A1(i)
      enddo
      if ( Hws.NE.0. ) then
!
!        set up gaussian convolution function
!
         afac = -alog(5.E-1)/(Hws/2.0E0)**2
         bfac = sqrt(afac/pi)
         do i = -NGC , NGC
            gc(i) = exp(-afac*(i*estep)**2)*bfac
         enddo
!
!        convolute now
!
         do i = Mine1 , N1
            A1(i) = CONVGAU(Ai,i,gc,ah,estep,NGC,N1)
         enddo
      endif
!
      if ( Snorm.LE.0.E0 ) then
!
!     determine maximum and set equal to 100.
!
         Smax = 0.E0
         do i = Mine1 , N1
            if ( A1(i).GT.Smax ) Smax = A1(i)
         enddo
         Snorm = 100.D0/Smax
      endif
!
!
      do i = Mine1 , N1
         A1(i) = A1(i)*Snorm
!
      enddo
!
      write (6,*) ' Out of broad'
      end
!*==CONVLO.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      double precision function CONVLO(Gam,A1,Ind,Estep,Mine1,N1,Pi)
! ==================================================================
      implicit none
!*--CONVLO686
!*** Start of declarations inserted by SPAG
      real*8 A1 , atm , atp , Estep , Gam , Pi , sum
      integer i , Ind , Mine1 , N1
!*** End of declarations inserted by SPAG
!
!     do the convolution integral by trapezoidal rule
!
      dimension A1(N1)
      write (6,*) ' Out of convlo'
      sum = 0.0
      do i = Mine1 , N1
         atp = atan(Estep*((i-Ind)+0.5)/Gam)
         atm = atan(Estep*((i-Ind)-0.5)/Gam)
         sum = sum + A1(i)*(atp-atm)
      enddo
      CONVLO = sum/Pi
      write (6,*) ' Out of convlo'
      end
!*==CONVGAU.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      double precision function CONVGAU(A1,Ind,Gc,Ah,Estep,Ngc,N1)
! ================================================================
      implicit none
!*--CONVGAU709
!*** Start of declarations inserted by SPAG
      real*8 A1 , Ah , Estep , Gc , sum
      integer i , imx , Ind , N1 , Ngc
!*** End of declarations inserted by SPAG
!
!     do the convolution integral by trapezoidal rule
!
      dimension A1(N1) , Gc(-Ngc:Ngc) , Ah(-Ngc:Ngc)
!
      write (6,*) ' In convgau'
      do i = -Ngc , Ngc
         Ah(i) = 0.D0
      enddo
      do i = 0 , Ngc
         imx = Ind + i
         if ( imx.GT.N1 ) goto 100
         Ah(i) = A1(imx)*Gc(i)
      enddo
 100  do i = -1 , -Ngc , -1
         imx = Ind + i
         if ( imx.LT.1 ) goto 200
         Ah(i) = A1(imx)*Gc(i)
      enddo
 200  sum = 0.D0
      do i = -Ngc , Ngc
         sum = sum + Ah(i)
      enddo
!
      CONVGAU = sum*Estep
!
      write (6,*) ' Out of convgau'
      end
!*==MATINTER.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      subroutine MATINTER(Crsfile,Efile,Nume,Efile2,Nume2,Nemf,Nemf2,Km,&
                        & Irange,E2,Emin,Emax,E,Nemax,Crs,Ifst,Ilast,   &
                        & Idoc,Jprt)
! ========================================================
!
      implicit none
!*--MATINTER749
!*** Start of declarations inserted by SPAG
      real*8 chelp , Crs , crs1 , Crsfile , E , E2 , ebar , Efile ,     &
           & Efile2 , ehelp , Emax , Emin
      integer i , i2 , Idoc , ie , iextra , if , if2l , if2r , ifl ,    &
            & ifr , Ifst , Ilast , insleft , insright , Irange ,        &
            & IRIGHT , Jprt , k , Km , KML
      integer kp , l , LEFTIND , lp , Nemax , Nemf , Nemf2 , nf , nint ,&
            & NM , Nume , Nume2
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
      parameter (KML=3,NM=50)
!
      dimension Crsfile(Nemf,Nemf2,-KML-1:KML,-KML-1:KML) , Efile(Nemf) &
              & , Efile2(Nemf2) , Irange(Nemf,Nemf2) , E(Nemax) ,       &
              & Crs(Nemax,-KML-1:KML,-KML-1:KML)
      dimension ehelp(NM) , chelp(NM) , crs1(NM,-KML-1:KML,-KML-1:KML)
      write (6,*) ' In matinter'
!
! find efile2 points immediately left and right
! of the actual e2 (if2l,if2r):
      if2l = LEFTIND(Efile2,Nume2,E2)
      if2r = IRIGHT(Efile2,Nume2,E2)
!
! valence band energy range: find efile points immediately
! left of emin and right of emax (ifl,ifr):
      ifl = LEFTIND(Efile,Nume,Emin)
      ifr = IRIGHT(Efile,Nume,Emax)
!
! for interpolation lateron, we may have to calculate the boundary
! values of the crs in a different way (to understand the whole
! thing, draw a diagram e vs. e2)
      insleft = 0
      insright = 0
      if ( Irange(ifl,if2r).EQ.0 ) insleft = 1
      if ( Irange(ifr,if2l).EQ.0 ) insright = 1
!
! now calculate the crs's for e2 and all values of efile
! between ifl and ifr and store them in crs1:
!
! boundary values (interpolation along boundaries of e/e2 range):
      nint = 5
      if ( Nume.LT.nint ) nint = Nume
      if ( insleft.EQ.1 ) then
         ebar = (Efile2(Nume2)+Efile2(1))/2.
         ehelp(ifl) = E2 - ebar + Efile(1)
!        if(ehelp(ifl).ge.emin)then
!          ifl=ifl+1
!          goto 295
!        endif
         do l = -Km - 1 , Km
            if ( l.NE.0 ) then
               do lp = -Km - 1 , l
                  if ( lp.NE.0 ) then
                     do i = Nume , 1 , -1
                        i2 = Nume2 - Nume + i
                        chelp(i) = Crsfile(i,i2,l,lp)
                     enddo
                     call LAGRANGE(Efile(1),chelp(1),Nume,nint,         &
                                 & ehelp(ifl),crs1(ifl,l,lp),iextra)
                  endif
               enddo
            endif
         enddo
         if ( Idoc.EQ.1 ) write (Jprt,99001)
99001    format (' non-grid value for lower e boundary')
      endif
      if ( insright.EQ.1 ) then
         ehelp(ifr) = E2 - Efile2(1) + Efile(1)
!        if(ehelp(ifr).le.emax)then
!          ifl=ifl-1
!          goto 297
!        endif
         do l = -Km - 1 , Km
            if ( l.NE.0 ) then
               do lp = -Km - 1 , l
                  if ( lp.NE.0 ) then
                     do i = 1 , Nume
                        chelp(i) = Crsfile(i,i,l,lp)
                     enddo
                     call LAGRANGE(Efile(1),chelp(1),Nume,nint,         &
                                 & ehelp(ifr),crs1(ifr,l,lp),iextra)
                  endif
               enddo
            endif
         enddo
         if ( Idoc.EQ.1 ) write (Jprt,99002)
99002    format (' non-grid value for upper e boundary')
      endif
      if ( Idoc.EQ.1 ) write (Jprt,99003) ifl , ifr
99003 format (' ifl,ifr: ',2I4)
!
! get the rest of crs1:
      do if = ifl , ifr
         if ( if.NE.ifl .OR. insleft.NE.1 ) then
            if ( if.NE.ifr .OR. insright.NE.1 ) then
               ehelp(if) = Efile(if)
               do l = -Km - 1 , Km
                  if ( l.NE.0 ) then
                     do lp = -Km - 1 , l
                        if ( lp.NE.0 ) then
                           do i = if , if + Nume - 1
                              chelp(i) = Crsfile(if,i,l,lp)
                           enddo
                           call LAGRANGE(Efile2(if),chelp(if),Nume,nint,&
                            & E2,crs1(if,l,lp),iextra)
                        endif
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
!
      if ( Idoc.EQ.1 ) then
         write (Jprt,99004)
99004    format (' crs1:')
         do i = ifl , ifr
            write (Jprt,99005) i , ehelp(i) ,                           &
                             & ((crs1(i,l,lp),lp=-Km-1,l),l=-Km-1,-1) , &
                             & ((crs1(i,l,lp),lp=-Km-1,-1),             &
                             & (crs1(i,l,lp),lp=1,l),l=1,Km)
99005       format (i3,f10.5,28E13.5)
         enddo
      endif

!
! now we have a function crs1 on a mesh efile and we want a function
! crs on a mesh e and that's all
      nf = ifr - ifl + 1
      nint = 6
      if ( nint.GT.nf ) nint = nf
      do ie = Ifst , Ilast
         do l = -Km - 1 , Km
            if ( l.NE.0 ) then
               do lp = -Km - 1 , l
                  if ( lp.NE.0 ) then
                     do i = ifl , ifr
                        chelp(i) = crs1(i,l,lp)
                     enddo
                     call LAGRANGE(Efile(ifl),chelp(ifl),nf,nint,E(ie), &
                                 & Crs(ie,l,lp),iextra)
                     if ( iextra.EQ.1 .AND. Idoc.EQ.1 )                 &
                        & write (Jprt,99006) ie , E(ie) ,               &
                        & ((Crs(ie,k,kp),kp=-Km-1,k),k=-Km-1,-1) ,      &
                        & ((Crs(ie,k,kp),kp=-Km-1,-1),                  &
                        & (Crs(ie,k,kp),kp=1,k),k=1,Km)
99006                format (' extrapolation: crs, ie=',i3,' e(ie)=',   &
                           & f10.5/28E13.5)
                  endif
               enddo
            endif
         enddo
      enddo
!
! test output:
      if ( Idoc.EQ.1 ) then
         write (Jprt,99007) E2
99007    format (' e2: ',f10.5)
         do ie = Ifst , Ilast
            write (Jprt,99008) ie , E(ie) ,                             &
                             & ((Crs(ie,l,lp),lp=-Km-1,l),l=-Km-1,-1) , &
                             & ((Crs(ie,l,lp),lp=-Km-1,-1),             &
                             & (Crs(ie,l,lp),lp=1,l),l=1,Km)
99008       format (i3,f10.5,28E13.4)
         enddo
         write (Jprt,99009)
99009    format (120('-'))
      endif
      write (6,*) ' out of matinter'
      end
!*==LAGRANGE.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      subroutine LAGRANGE(X,Y,Ndim,Npkt,Xval,Yinter,Iextra)
! =========================================================
      implicit none
!*--LAGRANGE944
!*** Start of declarations inserted by SPAG
      integer Iextra , ifst , igr , ileft , Ndim , Npkt
      real*8 X , Xval , Y , Yinter
!*** End of declarations inserted by SPAG
      dimension X(Ndim) , Y(Ndim)
!
! Interpolation mit Komfort
! iextra = 0, wenn keine Extrapolation
!          1, wenn Extrapolation nach oben
!         -1, wenn Extrapolation nach unten
!
      write (6,*) ' in lagrange'
      Iextra = 0
      if ( Xval.LT.X(1) ) Iextra = -1
      if ( Xval.GT.X(Ndim) ) Iextra = 1
!
! Suche 1. Punkt rechts von xval, falls es ihn gibt:
      igr = 0
 100  igr = igr + 1
      if ( X(igr).EQ.Xval ) then
         Yinter = Y(igr)
         return
      endif
      if ( X(igr).LT.Xval .AND. igr.LT.Ndim ) goto 100
!
! Bestimme 1. Punkt f. Interpolation
      ileft = Npkt/2
      if ( mod(Npkt,2).NE.0 ) ileft = ileft + 1
      ifst = igr - ileft
!
      if ( ifst.LE.0 ) ifst = 1
      if ( ifst.GT.Ndim-Npkt+1 ) ifst = Ndim - Npkt + 1
!
! Lasst uns interpolieren:
      call INTER(X(ifst),Y(ifst),Npkt,Xval,Yinter)
! Und das war's dann auch schon, Leute!
      write (6,*) ' Out of Lagrange'
      end
!*==INTER.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      subroutine INTER(R,P,N,Rs,Ps)
! =================================
      implicit none
!*--INTER987
!*** Start of declarations inserted by SPAG
      real*8 denom , P , Ps , R , Rs , term
      integer i , j , N
!*** End of declarations inserted by SPAG
!
!     interpolate via lagrange
!
      dimension R(N) , P(N)
      write (6,*) ' In inter'
      Ps = 0.E0
      do j = 1 , N
         term = 1.E0
         denom = 1.E0
         do i = 1 , N
            if ( i.NE.j ) then
               denom = denom*(R(j)-R(i))
               term = term*(Rs-R(i))
            endif
         enddo
         Ps = Ps + term*P(j)/denom
      enddo
      write (6,*) ' Out of inter'
      end
!*==LEFTIND.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      integer function LEFTIND(Field,Ndim,value)
      implicit none
!*--LEFTIND1014
!*** Start of declarations inserted by SPAG
      integer i , Ndim
!*** End of declarations inserted by SPAG
! ==============================================
! index of the field element immediately left of value
! if value is lower than field(1), 1 is returned
!      implict real*8 (a-h,o-z)
      real*8 Field , value
      dimension Field(Ndim)
      do i = Ndim , 1 , -1
         if ( Field(i).LT.value ) goto 100
      enddo
 100  LEFTIND = i
      if ( LEFTIND.LT.1 ) LEFTIND = 1
      end
!*==IRIGHT.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      integer function IRIGHT(Field,Ndim,value)
      implicit none
!*--IRIGHT1033
!*** Start of declarations inserted by SPAG
      integer i , Ndim
!*** End of declarations inserted by SPAG
! =============================================
! index of the field element immediately right of value
! if value is greater than fieldhndim), ndim is returned
!      implict real*8 (a-h,o-z)
      real*8 Field , value
      dimension Field(Ndim)
      do i = 1 , Ndim
         if ( Field(i).GT.value ) goto 100
      enddo
 100  IRIGHT = i
      if ( IRIGHT.GT.Ndim ) IRIGHT = Ndim
      end
!*==SIMPSON.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      double precision function SIMPSON(Fct,N,Dx)
! ===============================================
      implicit none
!*--SIMPSON1053
!*** Start of declarations inserted by SPAG
      real*8 Dx , Fct , sum , x0 , x1 , x2
      integer i , N
!*** End of declarations inserted by SPAG
! integration of fct using simpson integration
! equidistant point mesh
      dimension Fct(N)
      sum = 0.
      x0 = Fct(1)
!
      do i = 3 , N , 2
         x1 = Fct(i-1)
         x2 = Fct(i)
         sum = sum + Dx*(x0+4*x1+x2)
         x0 = x2
      enddo
      SIMPSON = sum/3.
      if ( mod(N,2).EQ.0 ) SIMPSON = SIMPSON + (Fct(N)+Fct(N-1))*Dx/2.
      end
!*==ORDER.spg  processed by SPAG 6.72Dc at 02:03 on 30 Oct 2023
      subroutine ORDER(F,N,ior)
! =============================
      implicit none
!*--ORDER1077
!*** Start of declarations inserted by SPAG
      real*8 F , g
      integer i , ior , j , l , N
!*** End of declarations inserted by SPAG
      dimension F(N) , ior(N)
!
! f is to be arranged in decreasing order
      do i = 1 , N
         ior(i) = i
      enddo
!
      do i = 1 , N - 1
         do j = i + 1 , N
            if ( F(j).GT.F(i) ) then
               g = F(j)
               F(j) = F(i)
               F(i) = g
               l = ior(j)
               ior(j) = ior(i)
               ior(i) = l
            endif
         enddo
      enddo
      end
!~~~~~~~~rcvvmat.f~~~begin~~~~~~~~~~~~!
!*==RCVVMAT.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      PROGRAM RCVVMAT
!
!
! ********************************************************************
!  RELATIVISTIC MATRIXELEMENTS FOR CORE-VALENCE-VALENCE AUGER SPECTRA
! ********************************************************************
!
! the highest possible values for momentum, l are set in l.par
!                   (lmax=2*lvmax+lcore+2)
!
      IMPLICIT NONE
!*--RCVVMAT13
!*** Start of declarations inserted by SPAG
      REAL*8 a0 , chuge , clight , coef , contr , crsecs , csmall , de ,&
           & dummy , dx , e , e0 , e2 , earr , earr2 , ec1 , emin2 ,    &
           & ep , ev , fc1
      REAL*8 fcwf , fvwf , fvwfp , gc1 , gcwf , gvwf , gvwfp , half ,   &
           & pi , result , rmt , rs , rws , sd , sde , sigd , sigde ,   &
           & signv , small , ssum
      REAL*8 sum1 , sum2 , tiny , v , w1v , w1vp , w3j , w6j , wv2 ,    &
           & wvp2 , x , x0 , xj1 , xs , y , y1 , y2 , yp , zeromt
      INTEGER i , i2 , ii , ip , irange , j , j1 , jj , jj1 , jj2 ,     &
            & jjp , jp , kap , kap1 , kap2 , kapp , kout , kout1 , l ,  &
            & l1
      INTEGER l2 , lam , lamp , lkap , lkapp , llam , llamp , lmax ,    &
            & lmax2 , lp , lval , lvmax , ne , ne2 , nei , nei2 , neip ,&
            & NEMAX , NEMAX2 , nmt
      INTEGER nonrel , norm , NRAD , nwf , nws
!*** End of declarations inserted by SPAG
!
      PARAMETER (NRAD=250)
      PARAMETER (NEMAX=15,NEMAX2=2*NEMAX)
      INCLUDE 'l.par'
!
      DIMENSION v(NRAD) , rs(NRAD) , xs(NRAD)
      DIMENSION earr(NEMAX) , earr2(NEMAX2) , irange(NEMAX,NEMAX2)
!
      DIMENSION gc1(NRAD) , fc1(NRAD)
      DIMENSION gvwf(NRAD,-lmax-1:lmax) , fvwf(NRAD,-lmax-1:lmax)
      DIMENSION gvwfp(NRAD,-lmax-1:lmax) , fvwfp(NRAD,-lmax-1:lmax)
      DIMENSION gcwf(NRAD,-lmax-1:lmax) , fcwf(NRAD,-lmax-1:lmax)
!
      DIMENSION sd(-lvmax-1:lvmax,-lvmax-1:lvmax)
      DIMENSION sde(-lvmax-1:lvmax,-lvmax-1:lvmax)
      DIMENSION ssum(-lvmax-1:lvmax,-lvmax-1:lvmax,NEMAX)
      DIMENSION crsecs(-lvmax-1:lvmax,-lvmax-1:lvmax,NEMAX,NEMAX2)
!
      REAL*8 idgg(-lmax-1:lmax,0:lmax) , idgf(-lmax-1:lmax,0:lmax) ,    &
           & idfg(-lmax-1:lmax,0:lmax) , idff(-lmax-1:lmax,0:lmax) ,    &
           & iegg(-lmax-1:lmax,0:lmax) , iegf(-lmax-1:lmax,0:lmax) ,    &
           & iefg(-lmax-1:lmax,0:lmax) , ieff(-lmax-1:lmax,0:lmax)
!
      INTEGER rl(-lmax-1:lmax) , rlb(-lmax-1:lmax) , rj(-lmax-1:lmax)
!
      DIMENSION w3j(-lmax-1:lmax,-lmax-1:lmax,0:lmax2)
      DIMENSION w6j(0:lvmax,0:lvmax,0:lmax,0:lmax,0:lmax)
!
      INTEGER in , wf , po , mat , pun , zatom
!
      CHARACTER*1 name(40) , norb1(5) , norb(5)
      CHARACTER*50 filen
!
      DATA tiny , small , half/1.D-10 , 0.01D0 , 0.5D0/
      DATA pi/3.1415926535897932D0/
      DATA ev/13.606/
      DATA csmall/274.0746/ , chuge/1.D+08/
!
! ********************************************************************
!
      in = 1
      wf = 2
      po = 3
      mat = 7
      pun = 8
!
      OPEN (in,FILE='rcvvmat.in',STATUS='old')
      OPEN (9,FILE='rcvvmat.err',STATUS='unknown')
!
!
 100  READ (in,99015,END=99999) filen
      OPEN (pun,FILE=filen,STATUS='unknown')
!
      READ (in,99015) filen
      OPEN (po,FILE=filen,STATUS='unknown')
!
      READ (in,99015) filen
      OPEN (wf,FILE=filen,STATUS='unknown')
!
      READ (in,99015) filen
      OPEN (mat,FILE=filen,STATUS='unknown')
!
!  fill up fields for the relativistic quantumnumbers
!
      DO kap = -lmax - 1 , -1
         l = -kap - 1
         rl(kap) = l
         rlb(kap) = l + 1
         rj(kap) = 2*l + 1
      ENDDO
      DO kap = 1 , lmax
         l = kap
         rl(kap) = l
         rlb(kap) = l - 1
         rj(kap) = 2*l - 1
      ENDDO
!
! read in maximum angular momentum quantumnumber for valence states
!
      READ (in,*) lval
      WRITE (pun,99001) lval
99001 FORMAT (1x,'LVAL=',t20,i4/)
!
! read in starting energy, increment for the energy panel, number of
! energy points for valence band
! (with respect to the muffin tin zero)
! the energy is measured in ryd
!
      READ (in,*) e0 , de , ne
!
! check facility for non-relativistic limit
!
      READ (in,*) nonrel
      IF ( nonrel.EQ.1 ) THEN
         clight = chuge
      ELSE
         clight = csmall
      ENDIF
!
!  normalization of the valence wavefunctions :
!                  norm eq 0  within the muffin tin
!                  norm eq 1  within the wigner-seitz sphere
!
      READ (in,*) norm
      WRITE (pun,99002) norm
99002 FORMAT (1x,'NORMALIZATION OF VALENCE WF. (0 FOR MT, 1 FOR WS):',  &
            & i4/)
!
! read in identifier of core state
!
      READ (in,99016) (norb1(i),i=1,5)
!
      READ (in,*)
!
! read in the name of the species,
!         the atomic number,
!         the first point on the logarithmic mesh,
!         the increment on the log. mesh,
!         the number of point corresponding to the muffin tin zero,
!         the wigner-seitz radius,
!         the lattice constant
!         the value of the muffin tin zero and
!         the potential in form of r*V(r)
!
      READ (po,99016) (name(i),i=1,40)
      READ (po,*) zatom , x0 , nmt , dx , rws , a0 , zeromt
      READ (po,*) (v(j),j=1,nmt)
!
      coef = 2.*zatom/clight
      coef = coef*coef
!
! set up the logarithmic and radial mesh
! transform the potential data to v(R)
!
      x = x0
      xs(1) = x0
      rs(1) = DEXP(x)
      signv = 1.
      IF ( v(1).GT.0. ) signv = -1.
      v(1) = signv*v(1)/rs(1) - zeromt
!
      DO j = 2 , nmt
         x = x + dx
         xs(j) = x
         rs(j) = DEXP(x)
         v(j) = signv*v(j)/rs(j) - zeromt
      ENDDO
      rmt = rs(nmt)
!
      nws = 0
      DO j = nmt + 1 , NRAD
         x = x + dx
         xs(j) = x
         rs(j) = DEXP(x)
         IF ( nws.EQ.0 .AND. rws.LT.rs(j) ) nws = j
         v(j) = 0.D0
      ENDDO
!
      WRITE (*,99016) (name(i),i=1,40)
      WRITE (pun,99016) (name(i),i=1,40)
      WRITE (pun,*)
      WRITE (pun,99003) zatom , x0 , dx , nmt , rmt , nws , rws , zeromt
99003 FORMAT (1x,' ZATOM =',t20,i4/1x,'    X0 =',t20,f5.1/1x,'    DX =',&
            & t20,f10.7/1x,'   NMT =',t20,i4/1x,'   RMT =',t20,f10.7/1x,&
             &'   NWS =',t20,i4/1x,'   RWS =',t20,f10.7/1x,'ZEROMT =',  &
            & t20,f5.1//1x,'RADIAL MESH AND POTENTIAL*R'/)
      WRITE (pun,99004) (rs(j),v(j)*rs(j),j=1,nws)
99004 FORMAT (4D20.10)
!
! read in the core wavefunctions as r*WF(r) and core energies (hartree)
! (the same radial mesh is supposed to be used as for the potential)
!
!
 200  READ (wf,99016) (norb(i),i=1,5)
      READ (wf,*) ec1
      ec1 = 2.*ec1
      READ (wf,*) kap1
      READ (wf,*) nwf
!
      IF ( nonrel.EQ.1 ) THEN
         DO i = 1 , nwf
            READ (wf,*) dummy , gc1(i)
            fc1(i) = 0.D0
         ENDDO
      ELSE
         DO i = 1 , nwf
            READ (wf,*) dummy , gc1(i) , fc1(i)
         ENDDO
      ENDIF
!
      IF ( norb(1).NE.norb1(1) .OR. norb(2).NE.norb1(2) .OR. norb(3)    &
         & .NE.norb1(3) ) GOTO 200
!
      WRITE (pun,99005) (norb(i),i=1,5) , ec1 , kap1
99005 FORMAT (//1x,'CORE WAVEFUNCTION'/1x,'ORB. =',t20,5A1/1x,' EC1 =', &
            & t20,f12.6/1x,'KAP1 =',t20,i4//)
      WRITE (pun,99006) (rs(j),gc1(j),fc1(j),j=1,nws)
99006 FORMAT (3D20.10)
!
      WRITE (pun,99007)
99007 FORMAT (//' ****************** END INPUT *******************'//)
!
!  set up energy panel for valence band :   earr
!            and       for continuum    :   earr2
!
      ii = 0
      DO i = 1 , ne
         earr(i) = e0 + (i-1)*de
      ENDDO
      ne2 = 2*ne - 1
      emin2 = 2.*e0 - ec1
      WRITE (pun,99008) (earr(i),i=1,ne)
99008 FORMAT (1x,'ENERGY PLANE'//(7x,15F7.2)/)
!
      DO i = 1 , ne2
         earr2(i) = emin2 + (i-1)*de
         DO j = 1 , ne
            jp = i + 1 - j
            IF ( jp.GE.1 .AND. jp.LE.ne ) THEN
               ii = ii + 1
               irange(j,i) = ii
            ELSE
               irange(j,i) = 0
            ENDIF
         ENDDO
         WRITE (pun,99009) earr2(i) , (irange(j,i),j=1,ne)
99009    FORMAT (f7.2,15(2x,i3,2x))
      ENDDO
!
      WRITE (pun,*)
!
!  initialize the angular integration coefficients
!
      l1 = rl(kap1)
      j1 = rj(kap1)
      xj1 = j1*half
      jj1 = j1 + 1
      y1 = DSQRT(kap1*kap1-coef)
      IF ( nonrel.EQ.1 ) y1 = DFLOAT(l1+1)
!
      CALL WIG3J(w3j,pun)
      CALL WIG6J(w6j,xj1,pun)
!
!
!         ************************************
!         * loop for the energy in continuum *
!         ************************************
!
      DO nei2 = 1 , ne2
         WRITE (6,*) nei2
!
! compute continuum wavefunctions
!
         e2 = earr2(nei2)
!
         kout = 0
         IF ( nei2.EQ.ne ) kout = 1
!
         WRITE (pun,99010) e2
99010    FORMAT (' e2=',f15.5)
         IF ( kout.EQ.1 ) WRITE (pun,*)
         IF ( kout.EQ.1 ) WRITE (pun,*) ' continuum wavefunctions'
!
         CALL WAFU(e2,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gcwf,fcwf,   &
                 & norm,2,nonrel,pun,lval,coef,kout)
!
!
!         ***************************************
!         * loop for the energy in valence band *
!         ***************************************
!
         DO nei = 1 , ne
            neip = nei2 + 1 - nei
            IF ( irange(nei,nei2).NE.0 ) THEN
!
               kout1 = 0
               IF ( nonrel.EQ.1 .AND. irange(nei,nei2).EQ.ne*ne )       &
                  & kout1 = 1
!
! compute valence wavefunctions
!
               e = earr(nei)
               ep = earr(neip)
!
               IF ( kout.EQ.1 ) WRITE (pun,*)
               IF ( kout.EQ.1 ) WRITE (pun,*) ' valence wavefunctions'
               WRITE (pun,99011) e
99011          FORMAT (t20,' e =',f15.5)
!
               CALL WAFU(e,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gvwf,   &
                       & fvwf,norm,1,nonrel,pun,lval,coef,kout)
!
               CALL WAFU(ep,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gvwfp, &
                       & fvwfp,norm,1,nonrel,pun,lval,coef,0)
!
!         ******************
!         * loop for kappa *
!         ******************
!
               DO kap = -lval - 1 , lval
                  IF ( kap.NE.0 ) THEN
                     y = DSQRT(kap*kap-coef)
                     IF ( nonrel.EQ.1 ) y = DFLOAT(rl(kap)+1)
                     jj = rj(kap) + 1
                     l = IDINT(jj*0.5+small) - 1
                     lkap = rl(kap)
!
                     IF ( kout1.EQ.1 ) WRITE (pun,*) ' kappa=' , kap
!
!         *******************
!         * loop for kappa' *
!         *******************
!
                     DO kapp = -lval - 1 , lval
                        IF ( kapp.NE.0 ) THEN
                           yp = DSQRT(kapp*kapp-coef)
                           IF ( nonrel.EQ.1 ) yp = DFLOAT(rl(kapp)+1)
                           jjp = rj(kapp) + 1
                           lp = IDINT(jjp*0.5+small) - 1
                           lkapp = rl(kapp)
!
                           IF ( kout1.EQ.1 ) WRITE (pun,*) ' kappa"=' , &
                              & kapp
!
!  compute coulomb and exchange integrals for non-vanishing angular
!  integration coefficients
!
                           DO lam = 0 , lmax
                              DO kap2 = -lmax - 1 , lmax
                                 idgg(kap2,lam) = 0.
                                 idfg(kap2,lam) = 0.
                                 idgf(kap2,lam) = 0.
                                 idff(kap2,lam) = 0.
                                 iegg(kap2,lam) = 0.
                                 iefg(kap2,lam) = 0.
                                 iegf(kap2,lam) = 0.
                                 ieff(kap2,lam) = 0.
                              ENDDO
                           ENDDO
!
                           DO lam = 0 , lmax
!
                              w1v = w3j(kap1,kap,lam)
                              w1vp = w3j(kap1,kapp,lam)
!
                              DO kap2 = -lmax - 1 , lmax
                                 IF ( kap2.NE.0 ) THEN
                                    y2 = DSQRT(kap2*kap2-coef)
                                    IF ( nonrel.EQ.1 )                  &
                                     & y2 = DFLOAT(rl(kap2)+1)
!
                                    wv2 = w3j(kap,kap2,lam)
                                    wvp2 = w3j(kapp,kap2,lam)
!
                                    IF ( DABS(w1v*wvp2).GT.tiny ) THEN
                                       CALL CEI(gc1,gvwf(1,kap),        &
                                        & gvwfp(1,kapp),gcwf(1,kap2),rs,&
                                        & nws,lam,y1,y,yp,y2,dx,rws,    &
                                        & result)
                                       idgg(kap2,lam) = result
                                       CALL CEI(gc1,gvwf(1,kap),        &
                                        & fvwfp(1,kapp),fcwf(1,kap2),rs,&
                                        & nws,lam,y1,y,yp,y2,dx,rws,    &
                                        & result)
                                       idgf(kap2,lam) = result
                                       CALL CEI(fc1,fvwf(1,kap),        &
                                        & gvwfp(1,kapp),gcwf(1,kap2),rs,&
                                        & nws,lam,y1,y,yp,y2,dx,rws,    &
                                        & result)
                                       idfg(kap2,lam) = result
                                       CALL CEI(fc1,fvwf(1,kap),        &
                                        & fvwfp(1,kapp),fcwf(1,kap2),rs,&
                                        & nws,lam,y1,y,yp,y2,dx,rws,    &
                                        & result)
                                       idff(kap2,lam) = result
                                    ENDIF
!
                                    IF ( DABS(w1vp*wv2).GT.tiny ) THEN
                                       CALL CEI(gc1,gvwfp(1,kapp),      &
                                        & gvwf(1,kap),gcwf(1,kap2),rs,  &
                                        & nws,lam,y1,yp,y,y2,dx,rws,    &
                                        & result)
                                       iegg(kap2,lam) = result
                                       CALL CEI(gc1,gvwfp(1,kapp),      &
                                        & fvwf(1,kap),fcwf(1,kap2),rs,  &
                                        & nws,lam,y1,yp,y,y2,dx,rws,    &
                                        & result)
                                       iegf(kap2,lam) = result
                                       CALL CEI(fc1,fvwfp(1,kapp),      &
                                        & gvwf(1,kap),gcwf(1,kap2),rs,  &
                                        & nws,lam,y1,yp,y,y2,dx,rws,    &
                                        & result)
                                       iefg(kap2,lam) = result
                                       CALL CEI(fc1,fvwfp(1,kapp),      &
                                        & fvwf(1,kap),fcwf(1,kap2),rs,  &
                                        & nws,lam,y1,yp,y,y2,dx,rws,    &
                                        & result)
                                       ieff(kap2,lam) = result
                                    ENDIF
                                 ENDIF
!
                              ENDDO
                           ENDDO
!
!
! ************* sum up for the direct term
!
!
                           sigd = 0.
!
                           DO lam = 0 , lmax
                              llam = 2*lam + 1
                              w1v = w3j(kap1,kap,lam)
                              w1v = w1v*w1v
!
                              DO kap2 = -lmax - 1 , lmax
                                 IF ( kap2.NE.0 ) THEN
                                    jj2 = rj(kap2) + 1
                                    wvp2 = w3j(kapp,kap2,lam)
                                    wvp2 = wvp2*wvp2
!
                                    sum1 = idgg(kap2,lam)               &
                                     & + idgf(kap2,lam) + idfg(kap2,lam)&
                                     & + idff(kap2,lam)
!
                                    sigd = sigd +                       &
                                     & jj1*jj2*w1v*wvp2*sum1*sum1/llam
!
                                    IF ( kout1.EQ.1 .AND. DABS(sum1)    &
                                     & .GT.tiny ) THEN
                                       WRITE (pun,*)                    &
                                         & ' nonzero for direct term'
                                       WRITE (pun,*) lam , kap2 , sum1
                                    ENDIF
                                 ENDIF
!
                              ENDDO
                           ENDDO
!
                           sd(kap,kapp) = sigd
!
!
! ************* sum up for the cross term
!
!
                           sigde = 0.
!
                           DO lam = 0 , lmax
                              llam = 2*lam + 1
                              w1v = w3j(kap1,kap,lam)
!
                              DO lamp = 0 , lmax
                                 llamp = 2*lam + 1
                                 w1vp = w3j(kap1,kapp,lamp)
!
                                 DO kap2 = -lmax - 1 , lmax
                                    IF ( kap2.NE.0 ) THEN
                                       jj2 = rj(kap2) + 1
                                       l2 = IDINT(jj2*0.5+small) - 1
                                       wv2 = w3j(kap,kap2,lamp)
                                       wvp2 = w3j(kapp,kap2,lam)
!
                                       sum1 = idgg(kap2,lam)            &
                                        & + idgf(kap2,lam)              &
                                        & + idfg(kap2,lam)              &
                                        & + idff(kap2,lam)
                                       sum2 = iegg(kap2,lamp)           &
                                        & + iegf(kap2,lamp)             &
                                        & + iefg(kap2,lamp)             &
                                        & + ieff(kap2,lamp)
!
                                       contr = jj1*jj2*w1v*w1vp*wv2*wvp2&
                                        & *sum1*sum2*w6j(l,lp,l2,lam,   &
                                        & lamp)
                                       sigde = sigde + contr
!
                                       IF ( kout1.EQ.1 .AND. DABS(contr)&
                                        & .GT.tiny ) THEN
                                         WRITE (pun,*)                  &
                                           & ' nonzero for cross term'
                                         WRITE (pun,*) lam , lamp ,     &
                                          & kap2 , sum1*sum2
                                       ENDIF
                                    ENDIF
!
                                 ENDDO
                              ENDDO
                           ENDDO
!
                           sde(kap,kapp) = sigde
!
                           ssum(kap,kapp,nei) = 2.0*(sigd-sigde)
                        ENDIF
!
                     ENDDO
                  ENDIF
               ENDDO
!
!
!         ***************************
!         * end of loop for kappa's *
!         ***************************
!
               IF ( kout1.EQ.1 ) THEN
!
                  WRITE (pun,*)
                  WRITE (pun,*) ' DIRECT TERM'
                  DO kap = -lval - 1 , lval
                     IF ( kap.NE.0 ) THEN
                        WRITE (pun,99017) (kap,kapp,kapp=-lval-1,-1) ,  &
                             & (kap,kapp,kapp=1,lval)
                        WRITE (pun,99018) (sd(kap,kapp),kapp=-lval-1,-1)&
                             & , (sd(kap,kapp),kapp=1,lval)
                     ENDIF
                  ENDDO
                  WRITE (pun,*)
                  WRITE (pun,*) ' CROSS TERM'
                  DO kap = -lval - 1 , lval
                     IF ( kap.NE.0 ) THEN
                        WRITE (pun,99017) (kap,kapp,kapp=-lval-1,-1) ,  &
                             & (kap,kapp,kapp=1,lval)
                        WRITE (pun,99018)                               &
                             & (sde(kap,kapp),kapp=-lval-1,-1) ,        &
                             & (sde(kap,kapp),kapp=1,lval)
                     ENDIF
                  ENDDO
!
               ENDIF
            ENDIF
!
         ENDDO
!
!  ************************************
!  construct symmetrized cross sections
!  ************************************
!
         DO i = 1 , ne
            ip = nei2 + 1 - i
            IF ( irange(i,nei2).NE.0 ) THEN
               DO kap = -lval - 1 , lval
                  IF ( kap.NE.0 ) THEN
                     crsecs(kap,kap,i,nei2)                             &
                      & = (ssum(kap,kap,i)+ssum(kap,kap,ip))*0.5
                     DO kapp = -lval - 1 , kap - 1
                        IF ( kapp.NE.0 ) crsecs(kap,kapp,i,nei2)        &
                           & = ssum(kap,kapp,i) + ssum(kapp,kap,ip)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
!
!
      ENDDO
!
!         ********************************
!         * end of loop for the energies *
!         ********************************
!
      WRITE (mat,99012) (name(i),i=1,40) , (norb1(i),i=1,5) , ec1 , e0 ,&
                      & de , ne , lval , a0
99012 FORMAT (1x,                                                       &
     &'RELATIVISTIC CORE-VALENCE-VALENCE AUGER                MATRIXELEM&
     &ENTS'/40A1/1x,'CORE STATE'/5A1/1x,                                &
    & 'ENERGY OF THE CORE STATE'/f12.5/1x,                              &
    & 'E0, DE, NE FOR VALENCE BAND'/2F6.2,i5/1x,                        &
     &'MAXIMAL L QUANTUMNUMBER FOR VALENCE STATES'/i1/1x,               &
     &'LATTICE CONSTANT'/f10.6/)
!
      WRITE (mat,99013) (((kap,kapp),kapp=-lval-1,kap),kap=-lval-1,-1) ,&
                      & (((kap,kapp),kapp=-lval-1,-1),                  &
                      & ((kap,kapp),kapp=1,kap),kap=1,lval)
99013 FORMAT (1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,28(i6,i3,'    '))
!
      DO i2 = 1 , ne2
         DO i = 1 , ne
            ip = i2 + 1 - i
            IF ( irange(i,i2).NE.0 ) THEN
               WRITE (mat,99014) i2 , i , earr2(i2) , earr(ip) , earr(i)&
                               & , ((crsecs(kap,kapp,i,i2),kapp=-lval-1,&
                               & kap),kap=-lval-1,-1) ,                 &
                               & ((crsecs(kap,kapp,i,i2),kapp=-lval-1,  &
                               & -1),(crsecs(kap,kapp,i,i2),kapp=1,kap),&
                               & kap=1,lval)
99014          FORMAT (2I3,3F6.2,28D13.5)
            ENDIF
         ENDDO
      ENDDO
!
      CLOSE (pun)
      CLOSE (po)
      CLOSE (wf)
      CLOSE (mat)
      GOTO 100
!
!
99015 FORMAT (a50)
99016 FORMAT (40A1)
99017 FORMAT (7(4x,i2,1x,i2,4x))
99018 FORMAT (7D13.5)
!
99999 END
!*==WAFU.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE WAFU(En,Rl,Rlb,Rj,V,Rs,X0,Dx,Rmt,Nmt,Rws,Nws,P,Q,Norm, &
                    & Nval,Nonrel,Pun,Lval,Coef,Kout)
!================================================
!
      IMPLICIT NONE
!*--WAFU639
!*** Start of declarations inserted by SPAG
      REAL*8 a1 , b1 , chuge , clight , Coef , cose , csmall , Dx ,     &
           & ekappa , En , eta , fb , fb1 , fn , fn1 , P , Q , qint ,   &
           & r , ratfg
      REAL*8 ratx , Rmt , rr , Rs , Rws , sign , sine , SINTG , sk ,    &
           & tane , tanx , V , vint , X0 , x1 , x2 , yk , yk2
      INTEGER i , kap , Kout , l , lact , lb , lmax , Lval , Nmt ,      &
            & Nonrel , Norm , NRAD , Nval , Nws
!*** End of declarations inserted by SPAG
!
      PARAMETER (NRAD=250)
      INCLUDE 'l.par'
!
      DIMENSION V(NRAD) , Rs(NRAD)
!
      DIMENSION P(NRAD,-lmax-1:lmax) , Q(NRAD,-lmax-1:lmax)
      DIMENSION rr(NRAD) , vint(-lmax-1:lmax) , eta(-lmax-1:lmax)
      DIMENSION tanx(-lmax-1:lmax) , ratx(-lmax-1:lmax)
      DIMENSION fb(0:lmax+1) , fn(0:lmax+1) , fb1(0:lmax+1) ,           &
              & fn1(0:lmax+1)
!
      INTEGER Rl(-lmax-1:lmax) , Rlb(-lmax-1:lmax) , Rj(-lmax-1:lmax)
      INTEGER Pun
!
      DATA csmall/274.0746/ , chuge/1.D+08/
!
!
      IF ( Nonrel.EQ.1 ) THEN
         clight = chuge
      ELSE
         clight = csmall
      ENDIF
!
      IF ( En.GT.0 ) THEN
         sign = 1.
         ekappa = DSQRT(En)
      ELSE
         sign = -1.
         ekappa = DSQRT(-En)
      ENDIF
!
      CALL SBF1(En,Rmt,fb,fn)
!
      IF ( Nval.EQ.1 ) lact = Lval
      IF ( Nval.EQ.2 ) lact = lmax
!
      DO kap = -lact - 1 , lact
         IF ( kap.NE.0 ) THEN
!
            l = Rl(kap)
            lb = Rlb(kap)
            IF ( kap.GT.0 ) sk = ekappa
            IF ( kap.LT.0 ) sk = -sign*ekappa
            yk = DSQRT(kap*kap-Coef)
!
            IF ( Nonrel.EQ.1 ) THEN
               yk = DFLOAT(l+1)
               lb = l + 1
               sk = -sign*ekappa
            ENDIF
!
            yk2 = yk + yk
!
            CALL COMDIR(En,kap,V,Nmt,Nws,Dx,X0,Q(1,kap),P(1,kap),ratfg, &
                      & Nonrel)
!
            tane = (ratfg*fb(l)-sk*fb(lb))/(ratfg*fn(l)-sk*fn(lb))
            ratx(kap) = ratfg
            tanx(kap) = tane
!
! valence or continuum normalization
!
            IF ( Nval.EQ.1 ) THEN
               a1 = sk*ekappa*(fn(lb)-fb(lb)/tane)*Rmt/Q(Nmt,kap)/clight
               b1 = ekappa*(fn(l)-fb(l)/tane)*Rmt/P(Nmt,kap)
            ELSE
               eta(kap) = DATAN(tane)
               cose = DCOS(eta(kap))
               sine = DSIN(eta(kap))
               a1 = sk*(cose*fb(lb)-sine*fn(lb))*Rmt/Q(Nmt,kap)/clight
               b1 = (cose*fb(l)-sine*fn(l))*Rmt/P(Nmt,kap)
!     if(b1.lt.0.) then
!       if(eta(kap).gt.pi) then
!       eta(kap)=eta(kap)-pi
!       else
!       eta(kap)=eta(kap)+pi
!       endif
!       goto 310
!     endif
            ENDIF
!
!
            DO i = 1 , Nmt
               P(i,kap) = P(i,kap)*b1
               Q(i,kap) = Q(i,kap)*a1
               rr(i) = P(i,kap)*P(i,kap) + Q(i,kap)*Q(i,kap)
            ENDDO
!
!
!  set up the wavefunctions beyond the muffin tin radius
!
            DO i = Nmt + 1 , Nws + 1
               r = Rs(i)
               CALL SBF1(En,r,fb1,fn1)
               IF ( Nval.EQ.1 ) THEN
                  Q(i,kap) = sk*ekappa*(fn1(lb)-fb1(lb)/tane)*r/clight
                  P(i,kap) = ekappa*(fn1(l)-fb1(l)/tane)*r
               ELSE
                  Q(i,kap) = sk*(cose*fb1(lb)-sine*fn1(lb))*r/clight
                  P(i,kap) = (cose*fb1(l)-sine*fn1(l))*r
               ENDIF
               rr(i) = P(i,kap)*P(i,kap) + Q(i,kap)*Q(i,kap)
            ENDDO
!
!
            IF ( Nval.NE.2 ) THEN
!
! normalize valence wavefuntions according to norm
!
               IF ( Norm.EQ.0 ) THEN
                  vint(kap) = SINTG(yk2,rr,Rs,Dx,Nmt)
               ELSE
                  x1 = SINTG(yk2,rr,Rs,Dx,Nws-1)
                  x2 = SINTG(yk2,rr,Rs,Dx,Nws)
                  vint(kap) = x1 + (x2-x1)*(Rws-Rs(Nws-1))              &
                            & /(Rs(Nws)-Rs(Nws-1))
               ENDIF
               qint = DSQRT(1./vint(kap))
               DO i = 1 , NRAD
                  P(i,kap) = P(i,kap)*qint
                  Q(i,kap) = Q(i,kap)*qint
               ENDDO
            ENDIF
         ENDIF
!
!
      ENDDO
!
!
      IF ( Kout.EQ.1 ) THEN
!
         IF ( Nval.EQ.1 ) THEN
!
            WRITE (Pun,99003)
            WRITE (Pun,99004) (kap,kap=-lact-1,-1) , (kap,kap=1,lact)
            WRITE (Pun,99005) (ratx(kap),kap=-lact-1,-1) ,              &
                            & (ratx(kap),kap=1,lact)
            WRITE (Pun,99001)
99001       FORMAT (/1x,'TANGENT PHASESHIFTS'/)
            WRITE (Pun,99004) (kap,kap=-lact-1,-1) , (kap,kap=1,lact)
            WRITE (Pun,99005) (tanx(kap),kap=-lact-1,-1) ,              &
                            & (tanx(kap),kap=1,lact)
!
         ELSE
!
            WRITE (Pun,99003)
            WRITE (Pun,99004) (kap,kap=-1,-lact-1,-1)
            WRITE (Pun,99005) (ratx(kap),kap=-1,-lact-1,-1)
            WRITE (Pun,99006) (kap,kap=1,lact)
            WRITE (Pun,99007) (ratx(kap),kap=1,lact)
            WRITE (Pun,99002)
99002       FORMAT (/1x,'PHASESHIFTS'/)
            WRITE (Pun,99004) (kap,kap=-1,-lact-1,-1)
            WRITE (Pun,99005) (eta(kap),kap=-1,-lact-1,-1)
            WRITE (Pun,99006) (kap,kap=1,lact)
            WRITE (Pun,99007) (eta(kap),kap=1,lact)
!
         ENDIF
!
      ENDIF
99003 FORMAT (/1x,'CF/G RATIO'/)
99004 FORMAT (2x,12(5x,i3,5x))
99005 FORMAT (2x,12D13.5)
99006 FORMAT (15x,11(6x,i2,5x))
99007 FORMAT (15x,11D13.5)
!
      END
!*==COMDIR.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE COMDIR(E1,Kappa,Za,Nrc,Nnk,Dx,X0,Q,P,Ratfg,Nonrel)
!======================================================================
!
!   integration of relativistic radial dirac equations by milne
!   method and calculation of ratio cf/g
!
      IMPLICIT NONE
!*--COMDIR825
!*** Start of declarations inserted by SPAG
      REAL*8 bgc , bgx , c , chuge , cin , csmall , Dx , e , E1 , hoc , &
           & P , pp , Q , qp , Ratfg , stval , sxk , sxm , t , tc
      REAL*8 test , u , uc , unp , unp2 , wc , wnp , wnp2 , x , X0 ,    &
           & xc , xk , z2 , Za
      INTEGER i , ik , jri , kap , Kappa , lkap , n , nit , Nnk ,       &
            & Nonrel , NRAD , Nrc
!*** End of declarations inserted by SPAG
!
      PARAMETER (NRAD=250)
!
      DIMENSION bgx(NRAD) , sxk(4) , sxm(4) , P(NRAD) , Q(NRAD) ,       &
              & pp(NRAD) , qp(NRAD) , Za(NRAD)
      INTEGER pun
      DATA test/1.E+05/ , pun/99/
      DATA csmall/274.0746/ , chuge/1.D+16/
!
      DO i = 1 , Nnk
         bgx(i) = Za(i)
      ENDDO
      kap = Kappa
      xk = DFLOAT(kap)
      jri = Nrc
      e = E1
      stval = X0
      z2 = -bgx(1)*EXP(stval)
      tc = DEXP(stval)
!
!
      IF ( Nonrel.EQ.1 ) THEN
         IF ( kap.LT.0 ) lkap = -kap - 1
         IF ( kap.GT.0 ) lkap = kap
         kap = -lkap - 1
         xk = DFLOAT(kap)
         u = -0.5*z2/(lkap+1.0)
         P(1) = 1.E-20
         Q(1) = u*1.E-20
!
      ELSE
!
         c = csmall
         cin = 1./(c*c)
!
         hoc = z2/c
         IF ( DABS(hoc/xk).LE.0.05 ) THEN
            u = (xk+DABS(xk))/hoc - 0.5*hoc/DABS(xk)
         ELSE
            u = (xk+DSQRT(xk*xk-hoc*hoc))/hoc
         ENDIF
         P(1) = 1.0E-20
         Q(1) = c*u*1.0E-20
!
      ENDIF
!
!
!     eq. 4.92
!
      IF ( Nonrel.NE.1 ) THEN
         pp(1) = tc*(cin*(e-bgx(1))+1.)*Q(1) - xk*P(1)
      ELSE
         pp(1) = tc*Q(1) - xk*P(1)
      ENDIF
!
!     eq. 4.90
!
      qp(1) = xk*Q(1) - tc*(e-bgx(1))*P(1)
!
!     eq. 4.91
!     runge-kutta method
!
      x = stval
      n = 1
 100  ik = 0
      xc = x
      bgc = bgx(n)
      wc = Q(n)
      uc = P(n)
 200  ik = ik + 1
      t = DEXP(xc)
!
      IF ( Nonrel.NE.1 ) THEN
         sxk(ik) = Dx*(-xk*uc+t*wc*(cin*(e-bgc)+1.))
      ELSE
         sxk(ik) = Dx*(-xk*uc+t*wc)
      ENDIF
!
      sxm(ik) = Dx*(xk*wc-t*(e-bgc)*uc)
      IF ( ik.EQ.2 ) THEN
         uc = uc + 0.5*(sxk(2)-sxk(1))
         wc = wc + 0.5*(sxm(2)-sxm(1))
         GOTO 200
      ELSEIF ( ik.EQ.3 ) THEN
         xc = xc + 0.5*Dx
         uc = uc + sxk(3) - 0.5*sxk(2)
         wc = wc + sxm(3) - 0.5*sxm(2)
         bgc = bgx(n+1)
         GOTO 200
      ELSEIF ( ik.EQ.4 ) THEN
         Q(n+1) = Q(n) + (sxm(1)+2.0*sxm(2)+2.0*sxm(3)+sxm(4))/6.0
         P(n+1) = P(n) + (sxk(1)+2.0*sxk(2)+2.0*sxk(3)+sxk(4))/6.0
!
         IF ( Nonrel.NE.1 ) THEN
            pp(n+1) = t*Q(n+1)*(cin*(e-bgc)+1.0) - xk*P(n+1)
         ELSE
            pp(n+1) = t*Q(n+1) - xk*P(n+1)
         ENDIF
!
         qp(n+1) = xk*Q(n+1) - t*(e-bgc)*P(n+1)
         x = x + Dx
         n = n + 1
         IF ( n.LT.6 ) GOTO 100
      ELSE
         xc = xc + 0.5*Dx
         uc = uc + 0.5*sxk(1)
         wc = wc + 0.5*sxm(1)
         bgc = 0.5*(bgc+bgx(n+1))
         GOTO 200
      ENDIF
!
!     milne method
!
 300  x = x + Dx
      t = DEXP(x)
!
      unp = P(n-5) + 0.3*Dx*(11.*pp(n)-14.*pp(n-1)+26.*pp(n-2)          &
          & -14.*pp(n-3)+11.0*pp(n-4))
!
      wnp = Q(n-5) + 0.3*Dx*(11.0*qp(n)-14.0*qp(n-1)+26.0*qp(n-2)       &
          & -14.0*qp(n-3)+11.0*qp(n-4))
      nit = 0
!
!
 400  IF ( Nonrel.NE.1 ) THEN
         pp(n+1) = t*(cin*(e-bgx(n+1))+1.0)*wnp - xk*unp
      ELSE
         pp(n+1) = t*wnp - xk*unp
      ENDIF
!
      qp(n+1) = xk*wnp - t*(e-bgx(n+1))*unp
!
      unp2 = P(n-3) + (7.0*pp(n+1)+32.0*pp(n)+12.0*pp(n-1)+32.0*pp(n-2) &
           & +7.0*pp(n-3))*2.0*Dx/45.0
!
      wnp2 = Q(n-3) + (7.0*qp(n+1)+32.0*qp(n)+12.0*qp(n-1)+32.0*qp(n-2) &
           & +7.0*qp(n-3))*2.0*Dx/45.0
!
      IF ( DABS(test*(unp2-unp)).LE.DABS(unp2) ) THEN
         IF ( DABS(test*(wnp2-wnp)).LE.DABS(wnp2) ) GOTO 500
      ENDIF
!
      IF ( nit.NE.5 ) THEN
         nit = nit + 1
         wnp = wnp2
         unp = unp2
         GOTO 400
      ENDIF
!
 500  Q(n+1) = wnp2
      P(n+1) = unp2
      n = n + 1
      IF ( n.LT.Nnk ) GOTO 300
      Ratfg = Q(jri)/P(jri)
!
      END
!*==CEI.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE CEI(F2,F4,F1,F3,R,N,Lambda,L2,L4,L1,L3,Dx,Rws,Result)
! ====================================================================
      IMPLICIT NONE
!*--CEI994
!*** Start of declarations inserted by SPAG
      REAL*8 arg1 , arg2 , arg3 , arg4 , Dx , F1 , F2 , F3 , F4 , R ,   &
           & Result , rl , rmlp1 , Rws , SINTG
      INTEGER i , Lambda , m , N , NDIM , NRAD
!*** End of declarations inserted by SPAG
!
!     calculate coulomb and exchange integrals
!     (not too) closely following zare, wood and company
!
!        p. marksteiner
!
!  result = int( int( f2(R) *f4(R) * R_<**lambda/R_>**(lambda+1) *
!                     f1(R')*f3(R') ))
!
      PARAMETER (NRAD=250,NDIM=NRAD)
!
      DIMENSION rl(NDIM) , rmlp1(NDIM)
      DIMENSION arg1(NDIM) , arg2(NDIM) , arg3(NDIM) , arg4(NDIM)
!
      DIMENSION F1(1) , F2(1) , F3(1) , F4(1) , R(1)
!
      REAL*8 L1 , L2 , L3 , L4 , ll
!
      m = N
      DO i = 1 , m
         rl(i) = R(i)**Lambda
         rmlp1(i) = 1./(rl(i)*R(i))
         arg1(i) = F2(i)*F4(i)*rl(i)
      ENDDO
!
      ll = L2 + L4 + Lambda
      CALL DINTG(ll,1,arg1,arg2,R,Dx,m,Rws)
! arg2 contains int from r=0 to r(i)
      ll = ll + 1
!
      DO i = 1 , m
         arg2(i) = arg2(i)*F1(i)*F3(i)*rmlp1(i)
      ENDDO
!
      ll = ll + L1 + L3 - Lambda - 1
!
      DO i = 1 , m
         arg3(i) = F2(i)*F4(i)*rmlp1(i)
      ENDDO
!
      ll = L2 + L4 - Lambda - 1
      CALL DINTG(ll,-1,arg3,arg4,R,Dx,m,Rws)
! arg4 contains int. from r(i) to rws
      ll = ll + 1
!
      DO i = 1 , m
         arg4(i) = arg4(i)*F1(i)*F3(i)*rl(i)
      ENDDO
!
      ll = ll + L1 + L3 + Lambda
      DO i = 1 , m
         arg2(i) = arg2(i) + arg4(i)
      ENDDO
! ll = l1+l2+l3+l4 for both arg2 and arg4
      Result = SINTG(ll,arg2,R,Dx,m)
!
      END
!*==SINTG.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      DOUBLE PRECISION FUNCTION SINTG(Ll,Fct,R,Dx,N)
! ========================================
!
! integration of fct over r using simpson integration
! r(i)=r0*exp(dx*(i-1))
!
      IMPLICIT NONE
!*--SINTG1065
!*** Start of declarations inserted by SPAG
      REAL*8 Dx , fact , Fct , R , rll , sum , tiny , x0 , x1 , x2
      INTEGER i , N
!*** End of declarations inserted by SPAG
!
      DIMENSION Fct(1) , R(1)
      REAL*8 Ll
!
      DATA tiny/1.0D-20/
!
      IF ( Ll+1.0.LE.tiny ) THEN
         WRITE (9,*) ' sintg: ll+1=' , Ll + 1
         sum = 0.
      ELSE
         rll = R(1)**Ll
         IF ( rll.LE.tiny ) THEN
            sum = 0.
         ELSE
            fact = Fct(1)/rll
            sum = fact*(R(1)**(Ll+1))/(Ll+1)
         ENDIF
      ENDIF
!
      x0 = Fct(1)*R(1)
      DO i = 3 , N , 2
         x1 = Fct(i-1)*R(i-1)
         x2 = Fct(i)*R(i)
         sum = sum + Dx*(x0+4.D0*x1+x2)
         x0 = x2
      ENDDO
      SINTG = sum/3.D0
      IF ( MOD(N,2).EQ.0 ) SINTG = SINTG + (Fct(N)+Fct(N-1))            &
                                 & /2.D0*(R(N)-R(N-1))
      END
!*==DINTG.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE DINTG(Ll,Idir,Fct,Yint,R,Dx,N,Rws)
! =================================================
!
! integration of fct(i) to yield yint(i)
! r(i)=r0*exp(dx*(i-1))
!
      IMPLICIT NONE
!*--DINTG1108
!*** Start of declarations inserted by SPAG
      REAL*8 corr , Dx , fact , Fct , R , rll , Rws , sum , tiny , x0 , &
           & x1 , Yint
      INTEGER i , Idir , j , N
!*** End of declarations inserted by SPAG
!
      DIMENSION Fct(1) , Yint(1) , R(1)
      REAL*8 Ll
!
      DATA tiny/1.0D-20/
!
      IF ( Idir.GT.0 ) THEN
!
         IF ( Ll+1.0.LE.tiny ) THEN
            WRITE (9,*) ' dintg: ll+1=' , Ll + 1
            sum = 0.
         ELSE
            rll = R(1)**Ll
            IF ( rll.LE.tiny ) THEN
               sum = 0.
            ELSE
               fact = Fct(1)/rll
               sum = fact*(R(1)**(Ll+1))/(Ll+1)
            ENDIF
         ENDIF
         Yint(1) = sum
!
         x0 = Fct(1)*R(1)
         DO i = 2 , N
            x1 = Fct(i)*R(i)
! trapezoidal rule:
            sum = sum + Dx*(x0+x1)/2.D0
            Yint(i) = sum
            x0 = x1
         ENDDO
!
      ELSE
!
         Yint(N) = 0.
         sum = 0.
         x1 = Fct(N)*R(N)
!
         DO i = N - 1 , 1 , -1
            x0 = Fct(i)*R(i)
            sum = sum + Dx*(x0+x1)/2.D0
            Yint(i) = sum
            IF ( i.EQ.N-3 ) THEN
               CALL INTER1(R(N-3),Yint(N-3),4,1,Rws,corr)
               sum = sum - corr
               DO j = N - 3 , N
                  Yint(j) = Yint(j) - corr
               ENDDO
            ENDIF
            x1 = x0
         ENDDO
!
      ENDIF
      END
!*==INTER1.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE INTER1(R,P,N,Id,Rs,Ps)
! =====================================
      IMPLICIT NONE
!*--INTER11171
!*** Start of declarations inserted by SPAG
      REAL*8 denom , P , Ps , R , Rs , term
      INTEGER i , Id , j , N
!*** End of declarations inserted by SPAG
!
!     interpolate via lagrange
!
      DIMENSION R(N) , P(N)
      Ps = 0.D0
      DO j = 1 , N , Id
         term = 1.D0
         denom = 1.D0
         DO i = 1 , N , Id
            IF ( i.NE.j ) THEN
               denom = denom*(R(j)-R(i))
               term = term*(Rs-R(i))
            ENDIF
         ENDDO
         Ps = Ps + term*P(j)/denom
      ENDDO
      END
!*==SBF.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE SBF(E,R,Fb,Fn)
! =============================
!
! spherical bessel and neuman functions for            e > 0
! modified spherical bessel and neuman functions for   e < 0
!
      IMPLICIT NONE
!*--SBF1201
!*** Start of declarations inserted by SPAG
      REAL*8 E , Fb , Fn , R , tl
      INTEGER l , lmax
!*** End of declarations inserted by SPAG
      INCLUDE 'l.par'
!
      COMPLEX*16 cfb(0:lmax+1) , cfn(0:lmax+1)
      COMPLEX*16 ecompl , eroot , x1 , x2 , sinx , cosx , cimmi
      DIMENSION Fb(0:lmax+1) , Fn(0:lmax+1)
!
      DATA cimmi/(0.0D0,-1.0D0)/
!
      ecompl = DCMPLX(E,0.D0)
      eroot = CDSQRT(ecompl)
      x1 = eroot*R
      x2 = x1*x1
      sinx = CDSIN(x1)
      cosx = CDCOS(x1)
!
      cfb(0) = sinx/x1
      cfb(1) = sinx/x2 - cosx/x1
      cfn(0) = -cosx/x1
      cfn(1) = -cosx/x2 - sinx/x1
!
      DO l = 2 , lmax + 1
         tl = DFLOAT(2*l-1)
         cfb(l) = tl*cfb(l-1)/x1 - cfb(l-2)
         cfn(l) = tl*cfn(l-1)/x1 - cfn(l-2)
      ENDDO
!
      IF ( E.LT.0. ) THEN
         cfn(0) = cfn(0)*cimmi
         DO l = 1 , lmax + 1
            cfb(l) = cfb(l)*cimmi**l
            cfn(l) = cfn(l)*cimmi**(l+1)
         ENDDO
      ENDIF
!
      DO l = 0 , lmax + 1
         Fb(l) = DREAL(cfb(l))
         Fn(l) = DREAL(cfn(l))
      ENDDO
!
      END
!*==SBF1.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE SBF1(E,R,Fb,Fn)
! =============================
!
! spherical bessel and neuman functions for            e > 0
! modified spherical bessel and neuman functions for   e < 0
!
      IMPLICIT NONE
!*--SBF11254
!*** Start of declarations inserted by SPAG
      REAL*8 cosx , E , eroot , Fb , Fn , R , sinx , tl , x1 , x2
      INTEGER l , lmax
!*** End of declarations inserted by SPAG
      INCLUDE 'l.par'
!
      DIMENSION Fb(0:lmax+1) , Fn(0:lmax+1)
!
!
!     shyp(x)=0.5d0*(dexp(x)-dexp(-x))
!     chyp(x)=0.5d0*(dexp(x)+dexp(-x))
!
      IF ( E.GT.0. ) THEN
!
         eroot = DSQRT(E)
         x1 = eroot*R
         x2 = x1*x1
         sinx = DSIN(x1)
         cosx = DCOS(x1)
!
         Fb(0) = sinx/x1
         Fb(1) = sinx/x2 - cosx/x1
         Fn(0) = -cosx/x1
         Fn(1) = -cosx/x2 - sinx/x1
!
         DO l = 2 , lmax + 1
            tl = DFLOAT(2*l-1)
            Fb(l) = tl*Fb(l-1)/x1 - Fb(l-2)
            Fn(l) = tl*Fn(l-1)/x1 - Fn(l-2)
         ENDDO
!
      ELSE
!
         eroot = DSQRT(-E)
         x1 = eroot*R
         x2 = x1*x1
         sinx = DSINH(x1)
         cosx = DCOSH(x1)
!
         Fb(0) = sinx/x1
         Fb(1) = -sinx/x2 + cosx/x1
         Fn(0) = cosx/x1
         Fn(1) = -cosx/x2 + sinx/x1
!
         DO l = 2 , lmax + 1
            tl = DFLOAT(2*l-1)
            Fb(l) = -tl*Fb(l-1)/x1 + Fb(l-2)
            Fn(l) = -tl*Fn(l-1)/x1 + Fn(l-2)
         ENDDO
!
      ENDIF
!
      END
!*==WIG3J.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE WIG3J(W3j,Pun)
! ===================================
!
      IMPLICIT NONE
!*--WIG3J1313
!*** Start of declarations inserted by SPAG
      REAL*8 half , halfm , result , small , thrcof , W3j , xlmat ,     &
           & xlmax , xlmin
      INTEGER i1 , ier , iwigtes , kap1 , kap2 , l1 , l2 , lam , lamax ,&
            & lamin , lmax , lmax2 , NDIM
!*** End of declarations inserted by SPAG
!
      PARAMETER (NDIM=100)
      INCLUDE 'l.par'
!
      DIMENSION W3j(-lmax-1:lmax,-lmax-1:lmax,0:lmax2)
      DIMENSION thrcof(NDIM)
      REAL*8 nul , j1 , j2
      INTEGER Pun
!
      DATA nul , small , half , halfm/0.D0 , 0.01D0 , 0.5D0 , -0.5D0/
      DATA iwigtes/0/
!
!
!                                   j1  lambda  j2
!    < kappa1 lambda kappa2 > =  (                  ) =
!                                  1/2    0    -1/2
!
!
!                                  lambda  j2   j1
!                             =  (                  )
!                                    0    -1/2  1/2
!
! for kappa1 = -lmax-1,-lmax,...,lmax-1,lmax
! for kappa2 = -lmax-1,-lmax,...,lmax-1,lmax
!
! convention for storage :  w3j(kappa1,kappa2,lambda)
!
!
      DO kap1 = -lmax - 1 , lmax
         DO kap2 = -lmax - 1 , lmax
            DO lam = 0 , lmax2
               W3j(kap1,kap2,lam) = nul
            ENDDO
         ENDDO
      ENDDO
!
      DO kap1 = -lmax - 1 , lmax
!
         IF ( kap1.NE.0 ) THEN
            IF ( kap1.LT.0 ) THEN
               j1 = -kap1 - half
               l1 = -kap1 - 1
            ELSE
               j1 = kap1 - half
               l1 = kap1
            ENDIF
!
            DO kap2 = -lmax - 1 , lmax
!
               IF ( kap2.NE.0 ) THEN
                  IF ( kap2.LT.0 ) THEN
                     j2 = -kap2 - half
                     l2 = -kap2 - 1
                  ELSE
                     j2 = kap2 - half
                     l2 = kap2
                  ENDIF
!
                  CALL REC3JJ(thrcof,j2,j1,halfm,half,xlmin,xlmax,xlmat,&
                            & NDIM,ier)
                  IF ( ier.GE.0 ) THEN
!
                     lamin = IDINT(xlmin+small)
                     lamax = IDINT(xlmax+small)
!
                     DO lam = lamin , lamax
!
                        i1 = lam - lamin + 1
                        IF ( MOD(lam-lamin,2).NE.1 ) THEN
                           W3j(kap1,kap2,lam) = thrcof(i1)
!
                           IF ( iwigtes.EQ.1 ) THEN
                              WRITE (Pun,99001) kap1 , l1 , j1 , kap2 , &
                                   & l2 , j2 , lam , result
!
99001                         FORMAT (2I4,f5.1,5x,2I4,f5.1,5x,i4,10x,   &
                                    & d17.10)
                           ENDIF
                        ENDIF
!
                     ENDDO
                  ENDIF
               ENDIF
!
            ENDDO
         ENDIF
      ENDDO
!
      END
!*==WIG6J.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE WIG6J(W6j,J1,Pun)
! ================================
!
      IMPLICIT NONE
!*--WIG6J1414
!*** Start of declarations inserted by SPAG
      REAL*8 coeff , half , sixcof , small , W6j , xlam , xlamma ,      &
           & xlammi , xlmat
      INTEGER ier , iwigtes , j3 , l , l2 , lam , lama , lama1 , lama2 ,&
            & lami , lami1 , lami2 , lamma , lammi , lamp , lmax , lp , &
            & lvmax , NDIM
!*** End of declarations inserted by SPAG
!
      PARAMETER (NDIM=100)
      INCLUDE 'l.par'
!
      DIMENSION W6j(0:lvmax,0:lvmax,0:lmax,0:lmax,0:lmax)
      DIMENSION sixcof(NDIM)
      INTEGER Pun
      REAL*8 nul , j , jp , J1 , j2
!
      DATA nul , small , half/0.D0 , 0.01D0 , 0.5D0/
      DATA iwigtes/0/
!
!
!  wigner 6j coefficients
!
!            j'  lambda'  j2         lambda'  j1  j
!          [                 ]  =  [                 ]
!            j   lambda   j1         lambda   j2  j'
!
! multiplied by (-)**(lambda+lambda'+1)
!
!
! for given j1
!
! convention for storage : w6j(l,l',l2,lambda,lambda') ,
!
!  with  j =(2*l +1)/2        l =0,1,...,lvmax
!        j'=(2*l'+1)/2        l'=0,1,...,lvmax
!        j2=(2*l2+1)/2        l2=0,1,...,lmax
!
!
      DO l = 0 , lvmax
         DO lp = 0 , lvmax
            DO l2 = 0 , lmax
               DO lam = 0 , lmax
                  DO lamp = 0 , lmax
                     W6j(l,lp,l2,lam,lamp) = nul
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!
      DO l = 0 , lvmax
         j = (2.*l+1.)*half
         DO lp = 0 , lvmax
            jp = (2.*lp+1.)*half
!
            lami1 = IDINT(DABS(J1-jp)+small)
            lama1 = IDINT(J1+jp+small)
!
            DO l2 = 0 , lmax
               j2 = (2.*l2+1.)*half
!
               lami2 = IDINT(DABS(j-j2)+small)
               lama2 = IDINT(j+j2+small)
               lami = MAX(lami1,lami2)
               lama = MIN(lama1,lama2)
!
               DO lam = 0 , lmax
                  IF ( lam.GE.lami .AND. lam.LE.lama ) THEN
                     xlam = DFLOAT(lam)
!
                     CALL REC6J(sixcof,J1,j,xlam,j2,jp,xlammi,xlamma,   &
                              & xlmat,NDIM,ier)
                     IF ( ier.GE.0 ) THEN
!
                        lammi = IDINT(xlammi+small)
                        lamma = IDINT(xlamma+small)
!
                        DO lamp = lammi , lamma
                           coeff = 1.0 - 2.0*MOD(lam+lamp+1,2)
                           W6j(l,lp,l2,lam,lamp)                        &
                            & = coeff*sixcof(lamp-lammi+1)
!
                           IF ( iwigtes.EQ.1 ) THEN
                              WRITE (Pun,99001) lamp , J1 , j , lam ,   &
                                   & j2 , j3 , sixcof(lamp-lammi+1)
!
99001                         FORMAT (i4,2F5.1,5x,i4,2F5.1,10x,d17.10)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
!
               ENDDO
!
            ENDDO
         ENDDO
      ENDDO
!
      END
!*==REC3JJ.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE REC3JJ(Thrcof,L2,L3,M2,M3,L1min,L1max,Lmatch,Ndim,Ier)
! =====================================================================
!
!  j1-recursion of 3j-coefficients
! recursive evaluation of 3j- and
! 6j-coefficients.  k. schulten, r.g. gordon.
! ref. in comp. phys. commun. 11 (1976) 269
!
!
      IMPLICIT NONE
!*--REC3JJ1526
!*** Start of declarations inserted by SPAG
      REAL*8 a1 , a1s , a2 , a2s , c1 , c1old , c2 , cnorm , denom ,    &
           & dv , eps , half , huge , oldfac , one , ratio , sign1 ,    &
           & sign2 , srhuge , srtiny
      REAL*8 sum1 , sum2 , sumbac , sumfor , sumuni , Thrcof , three ,  &
           & thresh , tiny , two , x , x1 , x2 , x3 , y , y1 , y2 , y3 ,&
           & zero
      INTEGER i , Ier , index , l1cmax , l1cmin , lstep , n , Ndim ,    &
            & nfin , nfinp1 , nfinp2 , nfinp3 , nlim , nstep2
!*** End of declarations inserted by SPAG
      REAL*8 l1 , L2 , L3 , m1 , M2 , M3 , L1min , L1max , newfac ,     &
           & Lmatch
      DIMENSION Thrcof(Ndim)
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
!
      Lmatch = zero
      m1 = -M2 - M3
!
!  check relative magnitude of l- and m-values
      IF ( L2-DABS(M2)+eps.GE.0 ) THEN
         IF ( L3-DABS(M3)+eps.GE.0 ) THEN
            IF ( DMOD(L2+DABS(M2)+eps,one).LT.eps+eps ) THEN
               IF ( DMOD(L3+DABS(M3)+eps,one).LT.eps+eps ) THEN
!
!
!
!  limits for l1
!
                  L1min = DMAX1(DABS(L2-L3),DABS(m1))
                  L1max = L2 + L3
                  IF ( L1min.LT.L1max-eps ) THEN
!
!
!
                     Ier = 0
                     nfin = IDINT(L1max-L1min+one+eps)
                     IF ( Ndim.LT.nfin ) THEN
!
!  dimension of thrcof not large enough to hold all the coefficients
!  required
!
                        Ier = -2
                        WRITE (9,99001) L2 , L3 , m1 , M2 , M3 , nfin , &
                                      & Ndim
99001                   FORMAT (///1x,'3j-coefficients',9x,'l1',        &
                              & 2F7.1/20x,3F7.1,4x,                     &
                               &'exceed storage provided  (',i4,',',i4, &
                               &')')
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
                        GOTO 100
                     ENDIF
                  ELSEIF ( L1min.LT.L1max+eps ) THEN
!
!
!  this is reached in case that l1 can take only one value,
!  i.e. l1min = l1max
!
                     Ier = 0
                     Thrcof(1) = (-one)**IDINT(DABS(L2+M2-L3+M3)+eps)   &
                               & /DSQRT(L1min+L2+L3+one)
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
99002 FORMAT (///1x,'3j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,       &
            &'do not satisfy the condition l2-/m2/ and l3-/m3/ ge zero '&
           & ,'and integer')
      RETURN
 100  lstep = lstep + 1
      l1 = l1 + one
!
!
      oldfac = newfac
      a1 = (l1+L2+L3+one)*(l1-L2+L3)*(l1+L2-L3)*(-l1+L2+L3+one)
      a2 = (l1+m1)*(l1-m1)
      newfac = DSQRT(a1*a2)
      IF ( l1.LT.one+eps ) THEN
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
         IF ( lstep.GT.2 ) c1old = DABS(c1)
         c1 = -(l1+l1-one)*dv/denom
      ENDIF
!
      IF ( lstep.GT.2 ) THEN
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
         IF ( lstep.NE.nfin ) THEN
!
!  see if last unnormalized 3j-coefficient exceeds srhuge
!
            IF ( DABS(x).GE.srhuge ) THEN
!
!  this is reached if last 3j-coefficient larger than srhuge
!  so that the recursion series thrcof(1), ... , thrcof(lstep)
!  has to be rescaled to prevent overflow
!
               Ier = Ier + 1
               DO i = 1 , lstep
                  IF ( DABS(Thrcof(i)).LT.srtiny ) Thrcof(i) = zero
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
            IF ( c1old.GT.DABS(c1) ) GOTO 100
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
 150     lstep = lstep + 1
         l1 = l1 - one
!
         oldfac = newfac
         a1s = (l1+L2+L3)*(l1-L2+L3-one)*(l1+L2-L3-one)*(-l1+L2+L3+two)
         a2s = (l1+m1-one)*(l1-m1-one)
         newfac = DSQRT(a1s*a2s)
!
         dv = -L2*(L2+one)*m1 + L3*(L3+one)*m1 + l1*(l1-one)*(M3-M2)
!
         denom = l1*newfac
         c1 = -(l1+l1-one)*dv/denom
         IF ( lstep.GT.2 ) THEN
!
!
            c2 = -(l1-one)*oldfac/denom
!
!  recursion to the next 3j-coefficient y
!
            y = c1*Thrcof(nfinp2-lstep) + c2*Thrcof(nfinp3-lstep)
!
            IF ( lstep.EQ.nstep2 ) THEN
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
               IF ( DABS(ratio).LT.one ) THEN
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
            ELSE
!
               Thrcof(nfinp1-lstep) = y
               sumbac = sum2
               sum2 = sum2 + (l1+l1-three)*y*y
!
!  see if last unnormalized 3j-coefficient exceeds srhuge
!
               IF ( DABS(y).GE.srhuge ) THEN
!
!  this is reached if last 3j-coefficient larger than srhuge
!  so that the recursion series thrcof(nfin), ... ,thrcof(nfin-lstep+1)
!  has to be rescaled to prevent overflow
!
                  Ier = Ier + 1
                  DO i = 1 , lstep
                     index = nfin - i + 1
                     IF ( DABS(Thrcof(index)).LT.srtiny ) Thrcof(index) &
                        & = zero
                     Thrcof(index) = Thrcof(index)/srhuge
                  ENDDO
                  sum2 = sum2/huge
!
!
                  sumbac = sumbac/huge
               ENDIF
               GOTO 150
            ENDIF
         ELSE
!
!  if l1 = l1max + 1  the third term in the recursion formula vanishes
!
            y = srtiny*c1
            Thrcof(nfin-1) = y
            sumbac = sum2
            sum2 = sum2 + tiny*(l1+l1-three)*c1*c1
!
            GOTO 150
         ENDIF
      ELSE
!
!
!  if l1 = l1min + 1  the third term in the recursion equation vanishes
!  , hence
         x = srtiny*c1
         Thrcof(2) = x
         sum1 = sum1 + tiny*(l1+l1+one)*c1*c1
         IF ( lstep.NE.nfin ) GOTO 100
!
         sumuni = sum1
      ENDIF
!
!
!  normalize 3j-coefficients
!
      cnorm = one/DSQRT(sumuni)
!
!  sign convention for last 3j-coefficient determines overall phase
!
      sign1 = DSIGN(one,Thrcof(nfin))
      sign2 = (-one)**IDINT(DABS(L2+M2-L3+M3)+eps)
      IF ( sign1*sign2.LE.0 ) cnorm = -cnorm
!
      IF ( DABS(cnorm).LT.one ) THEN
!
         thresh = tiny/DABS(cnorm)
         DO n = 1 , nfin
            IF ( DABS(Thrcof(n)).LT.thresh ) Thrcof(n) = zero
            Thrcof(n) = cnorm*Thrcof(n)
         ENDDO
         GOTO 99999
      ENDIF
!
      DO n = 1 , nfin
         Thrcof(n) = cnorm*Thrcof(n)
      ENDDO
      RETURN
!
99999 END
!*==REC6J.spg  processed by SPAG 6.72Dc at 04:58 on  3 Nov 2023
      SUBROUTINE REC6J(Sixcof,L2,L3,L4,L5,L6,L1min,L1max,Lmatch,Ndim,   &
                     & Ier)
! ===================================================================
!
!  j1-recursion of 6j-coefficients
! recursive evaluation of 3j- and
! 6j-coefficients.  k. schulten, r.g. gordon.
! ref. in comp. phys. commun. 11 (1976) 269
!
!
      IMPLICIT NONE
!*--REC6J1886
!*** Start of declarations inserted by SPAG
      REAL*8 a1 , a1s , a2 , a2s , c1 , c1old , c2 , cnorm , denom ,    &
           & dv , eps , half , huge , oldfac , one , ratio , sign1 ,    &
           & sign2 , Sixcof , srhuge
      REAL*8 srtiny , sum1 , sum2 , sumbac , sumfor , sumuni , three ,  &
           & thresh , tiny , two , x , x1 , x2 , x3 , y , y1 , y2 , y3 ,&
           & zero
      INTEGER i , Ier , index , l1cmax , l1cmin , lstep , n , Ndim ,    &
            & nfin , nfinp1 , nfinp2 , nfinp3 , nlim , nstep2
!*** End of declarations inserted by SPAG
      REAL*8 l1 , L2 , L3 , L4 , L5 , L6 , L1min , L1max , newfac ,     &
           & Lmatch
      DIMENSION Sixcof(Ndim)
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
!
      Lmatch = zero
!
!
!
!  check if 6j-coefficients obey selection rules
!
      IF ( DMOD(L2+L3+L5+L6+eps,one).LT.eps+eps ) THEN
         IF ( DMOD(L4+L2+L6+eps,one).LT.eps+eps ) THEN
!
            IF ( L4+L2.GE.L6 ) THEN
               IF ( L4-L2+L6.GE.0 ) THEN
                  IF ( -L4+L2+L6.GE.0 ) THEN
!
                     IF ( L4+L3.GE.L5 ) THEN
                        IF ( L4-L3+L5.GE.0 ) THEN
                           IF ( -L4+L3+L5.GE.0 ) THEN
!
!  limits for l1
!
                              L1min = DMAX1(DABS(L2-L3),DABS(L5-L6))
                              L1max = DMIN1(L2+L3,L5+L6)
                              IF ( L1min.LT.L1max-eps ) THEN
!
!
                                 Ier = 0
                                 nfin = IDINT(L1max-L1min+one+eps)
                                 IF ( Ndim.LT.nfin ) THEN
!
!  this is reached if array sixcof not large enough to hold all
!  6j - coefficients required
!
                                    Ier = -2
                                    WRITE (9,99001) L2 , L3 , L4 , L5 , &
                                     & L6 , nfin , Ndim
99001                               FORMAT (///1x,'6j-coefficients',9x, &
                                      &'l1',2F7.1/20x,3F7.1,4x,         &
                                      &'exceed storage provided  (',i4, &
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
                                    GOTO 100
                                 ENDIF
                              ELSEIF ( L1min.LT.L1max+eps ) THEN
!
!
!  this is reached in case that l1 can take only one value
!
                                 Ier = 0
                                 Sixcof(1) = (-one)                     &
                                  & **IDINT(L2+L3+L5+L6+eps)            &
                                  & /DSQRT((L1min+L1min+one)*(L4+L4+one)&
                                  & )
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
99002 FORMAT (///1x,'6j-coefficients',9x,'l1',2F7.1,4x,                 &
             &'do not satisfy triangular conditions or'/20x,3F7.1,4x,   &
             &'l2+l3+l5+l6 or l2+l4+l6 not integer')
      RETURN
 100  lstep = lstep + 1
      l1 = l1 + one
!
      oldfac = newfac
      a1 = (l1+L2+L3+one)*(l1-L2+L3)*(l1+L2-L3)*(-l1+L2+L3+one)
      a2 = (l1+L5+L6+one)*(l1-L5+L6)*(l1+L5-L6)*(-l1+L5+L6+one)
      newfac = DSQRT(a1*a2)
!
      IF ( l1.LT.one+eps ) THEN
!
!  if l1 = 1   (l1 - 1) has to be factored out of dv, hence
!
         c1 = -two*(L2*(L2+one)+L5*(L5+one)-L4*(L4+one))/newfac
      ELSE
!
         dv = two*(L2*(L2+one)*L5*(L5+one)+L3*(L3+one)*L6*(L6+one)      &
            & -l1*(l1-one)*L4*(L4+one))                                 &
            & - (L2*(L2+one)+L3*(L3+one)-l1*(l1-one))                   &
            & *(L5*(L5+one)+L6*(L6+one)-l1*(l1-one))
!
         denom = (l1-one)*newfac
!
!
         IF ( lstep.GT.2 ) c1old = DABS(c1)
         c1 = -(l1+l1-one)*dv/denom
      ENDIF
!
      IF ( lstep.GT.2 ) THEN
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
         IF ( lstep.NE.nfin ) THEN
!
!  see if last unnormalized 6j-coefficient exceeds srhuge
!
            IF ( DABS(x).GE.srhuge ) THEN
!
!  this is reached if last 6j-coefficient larger than srhuge
!  so that the recursion series sixcof(1), ... ,sixcof(lstep)
!  has to be rescaled to prevent overflow
!
               Ier = Ier + 1
               DO i = 1 , lstep
                  IF ( DABS(Sixcof(i)).LT.srtiny ) Sixcof(i) = zero
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
            IF ( c1old.GT.DABS(c1) ) GOTO 100
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
 150     lstep = lstep + 1
         l1 = l1 - one
!
         oldfac = newfac
         a1s = (l1+L2+L3)*(l1-L2+L3-one)*(l1+L2-L3-one)*(-l1+L2+L3+two)
         a2s = (l1+L5+L6)*(l1-L5+L6-one)*(l1+L5-L6-one)*(-l1+L5+L6+two)
         newfac = DSQRT(a1s*a2s)
!
         dv = two*(L2*(L2+one)*L5*(L5+one)+L3*(L3+one)*L6*(L6+one)      &
            & -l1*(l1-one)*L4*(L4+one))                                 &
            & - (L2*(L2+one)+L3*(L3+one)-l1*(l1-one))                   &
            & *(L5*(L5+one)+L6*(L6+one)-l1*(l1-one))
!
         denom = l1*newfac
         c1 = -(l1+l1-one)*dv/denom
         IF ( lstep.GT.2 ) THEN
!
!
            c2 = -(l1-one)*oldfac/denom
!
!  recursion to the next 6j - coefficient y
!
            y = c1*Sixcof(nfinp2-lstep) + c2*Sixcof(nfinp3-lstep)
            IF ( lstep.NE.nstep2 ) THEN
               Sixcof(nfinp1-lstep) = y
               sumbac = sum2
               sum2 = sum2 + (l1+l1-three)*y*y
!
!  see if last unnormalized 6j-coefficient exceeds srhuge
!
               IF ( DABS(y).GE.srhuge ) THEN
!
!  this is reached if last 6j-coefficient larger than srhuge
!  so that the recursion series sixcof(nfin), ... ,sixcof(nfin-lstep+1)
!  has to be rescaled to prevent overflow
!
                  Ier = Ier + 1
                  DO i = 1 , lstep
                     index = nfin - i + 1
                     IF ( DABS(Sixcof(index)).LT.srtiny ) Sixcof(index) &
                        & = zero
                     Sixcof(index) = Sixcof(index)/srhuge
                  ENDDO
                  sumbac = sumbac/huge
!
                  sum2 = sum2/huge
               ENDIF
               GOTO 150
            ENDIF
         ELSE
!
!  if l1 = l1max + 1 the third term in the recursion equation vanishes
!
            y = srtiny*c1
            Sixcof(nfin-1) = y
            IF ( lstep.NE.nstep2 ) THEN
               sumbac = sum2
               sum2 = sum2 + (l1+l1-three)*c1*c1*tiny
               GOTO 150
            ENDIF
         ENDIF
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
         IF ( DABS(ratio).LT.one ) THEN
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
         IF ( lstep.NE.nfin ) GOTO 100
!
         sumuni = sum1
      ENDIF
!
!
!  normalize 6j-coefficients
!
      cnorm = one/DSQRT((L4+L4+one)*sumuni)
!
!  sign convention for last 6j-coeff. determines overall phase
!
      sign1 = DSIGN(one,Sixcof(nfin))
      sign2 = (-one)**IDINT(L2+L3+L5+L6+eps)
      IF ( sign1*sign2.LE.0 ) cnorm = -cnorm
!
      IF ( DABS(cnorm).LT.one ) THEN
!
         thresh = tiny/DABS(cnorm)
         DO n = 1 , nfin
            IF ( DABS(Sixcof(n)).LT.thresh ) Sixcof(n) = zero
            Sixcof(n) = cnorm*Sixcof(n)
         ENDDO
         GOTO 99999
      ENDIF
!
      DO n = 1 , nfin
         Sixcof(n) = cnorm*Sixcof(n)
      ENDDO
      RETURN
!
99999 END

