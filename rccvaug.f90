PROGRAM rccvaug
!
! calculate relativistic core-core-valence spectrum
!
!     subroutines needed: ccvopen,synspec,broad
!
!
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a0 , bla , cfac , conv , decore , e0 , e1 , ef , em0 , fl , fl1 , fl2 , hwc , hws , hwv , pdos , phelp , pi , pmat ,     &
        & pspect
   REAL*8 ptot , shelp , stot , z
   INTEGER i , j , k , krel , l , lmax , LMAX0 , maxd , maxm , mine1 , n1 , NE0 , NE1 , NEM0
!*** End of declarations inserted by SPAG
   PARAMETER (NE0=500,NE1=301,NEM0=71)
   PARAMETER (LMAX0=3)
!
   CHARACTER*80 head
   CHARACTER*1 name(10) , iddos(10) , idcore1(5) , idcore3(5)
!
   DIMENSION e0(NE0) , e1(NE1) , em0(NEM0)
!
   DIMENSION pdos(NE0,-LMAX0-1:LMAX0)
   DIMENSION pspect(NE0,-LMAX0-1:LMAX0) , stot(NE0) , shelp(NE0)
   DIMENSION ptot(NE1) , phelp(NE1)
   DIMENSION pmat(NEM0,-LMAX0-1:LMAX0) , z(-LMAX0-1:LMAX0)
   REAL*8 norm
!
   INTEGER dos , mat , tty , prt , plot
!
   DATA cfac/13.606E0/ , norm/0.0E0/
   DATA pi/3.141596/
!
!     say hello....
!
   PRINT 99001
   PRINT 99002
99002 FORMAT ('           this is relativistic CCV Auger  ')
   PRINT 99001
!
! open files
!
   CALL ccvopen(dos,mat,tty,prt,plot)
!
! enquiry.....
!
!   identifier
!
   READ (tty,99003) head
99003 FORMAT (a80)
   WRITE (plot,99004) head
99004 FORMAT (1x,a80)
!
!   spectrometer resolution (fwhm in Ev)
!
   READ (tty,*) hws
!
!   halfwidth of lorentzian for valence band lifetime broadening
!
   READ (tty,*) hwv
!
!
!   halfwidth of lorentzian for core-state lifetime broadening
!
   READ (tty,*) hwc
!
!   non-relativistic or relativistic densities of states  (0/1)
!
   READ (tty,*) krel
!
!
! ***************  read in matrix elements  ***************
!
!   (energy is given in ryd)
!
! in relativistic case the ordering of kappa's in the dos file:
!          -1   1  -2   2  -3   3  -4
!                                          but in the mat file:
!          -4  -3  -2  -1   1   2   3
!
!
   READ (mat,*)
   READ (mat,99005) (name(j),j=1,10)
   READ (mat,*)
   READ (mat,99005) (idcore1(j),j=1,5)
   READ (mat,*)
   READ (mat,99005) (idcore3(j),j=1,5)
   READ (mat,*)
   READ (mat,*) decore
   READ (mat,*)
   READ (mat,*) lmax
   READ (mat,*)
   READ (mat,*) a0
   READ (mat,*)
   READ (mat,*)
!
   conv = 2.*pi/a0
   conv = conv*conv
   IF ( krel==0 ) conv = 1.
!
   i = 0
   DO
      i = i + 1
      READ (mat,*,END=100) em0(i) , (pmat(i,k),k=-lmax-1,-1) , (pmat(i,k),k=1,lmax)
   ENDDO
 100  maxm = i - 1
!
!
! ***************  read in densities of states  ***************
!
!
   READ (dos,99005) (iddos(j),j=1,10)
   READ (dos,*)
   READ (dos,*) ef
   ef = conv*ef
   READ (dos,*) lmax
   READ (dos,*) maxd
   READ (dos,*)
!
   DO i = 1 , maxd
!
      IF ( krel==0 ) THEN
         READ (dos,*) e0(i) , (z(l),l=0,lmax) , bla
         e0(i) = e0(i)*conv
         pdos(i,-1) = z(0)
         DO l = 1 , lmax
            fl = dfloat(2*l+1)
            fl = fl + fl
            fl1 = dfloat(2*l)
            fl2 = dfloat(2*l+2)
            fl1 = fl1/fl
            fl2 = fl2/fl
            pdos(i,l) = z(l)*fl1
            pdos(i,-l-1) = z(l)*fl2
         ENDDO
      ELSE
         READ (dos,*) e0(i) , bla , pdos(i,-1) , (pdos(i,l),pdos(i,-l-1),l=1,lmax)
         e0(i) = e0(i)*conv
      ENDIF
!
   ENDDO
!
!
! printout
!
   WRITE (prt,99007) (iddos(j),j=1,10) , (idcore1(j),j=1,5) , (idcore3(j),j=1,5) , decore , ef , hws , hwv , hwc
99007 FORMAT (' Relativistic CCV Auger spectrum of ',10A1/'                   final core state ',                                  &
             &5A1/'                 initial core state ',5A1/'                  transition energy ',f10.3,                         &
             &' ryd'/'                       fermi-energy ',f10.3,' ryd'/'              spectrometer (fwhm) : ',f10.3,             &
             &'  Ev'/'         valence life-time (fwhm) : ',f10.3,'  Ev'/'      core state life-time (fwhm) : ',f10.3,'  Ev'/)
   WRITE (prt,99008)
99008 FORMAT (/' Matrixelements'/)
   DO i = 1 , maxm
      WRITE (prt,99011) em0(i) , (pmat(i,k),k=-lmax-1,-1) , (pmat(i,k),k=1,lmax)
   ENDDO
   WRITE (prt,99009)
99009 FORMAT (/' DOS'/)
   DO i = 1 , maxd
      WRITE (prt,99010) e0(i) , pdos(i,-1) , (pdos(i,l),pdos(i,-l-1),l=1,lmax)
99010 FORMAT (10F11.5)
   ENDDO
!
!
! ******************  put together spectrum
!
   PRINT * , 'starting synspec'
   CALL synspec(e0,pdos,NE0,maxd,pspect,stot,em0,pmat,NEM0,maxm,LMAX0,lmax,prt,1)
!
   WRITE (prt,*) '  total unbroadened spectrum'
   DO i = 1 , maxd
      WRITE (prt,99011) (e0(i)-ef)*cfac , stot(i)
   ENDDO
!
! *****************************
!
!     convolute spectrum with lorentzian in order to
!     account for many-body effects and gaussian
!     spectrometer resolution
!
!
   PRINT * , 'starting broad'
   CALL broad(e0,stot,NE0,maxd,mine1,prt,hwc,hwv,hws,e1,ptot,phelp,NE1,ef,cfac,norm)
!
!
! *********** broadening for partial spectra
!
!
   DO k = -lmax - 1 , lmax
      IF ( k/=0 ) THEN
!
         DO i = 1 , NE0
            shelp(i) = pspect(i,k)
         ENDDO
         PRINT * , 'starting broad'
         CALL broad(e0,shelp,NE0,maxd,mine1,prt,hwc,hwv,hws,e1,pspect(1,k),phelp,NE1,ef,cfac,norm)
      ENDIF
   ENDDO
!
!
! ************** output
!
   n1 = NE1 - mine1 + 1
   WRITE (plot,99006) (idcore1(j),j=1,5) , (idcore3(j),j=1,5) , decore*cfac
99006 FORMAT (2x,5A1,2x,5A1/f15.5)
   WRITE (plot,99013) (l,-l-1,l=1,lmax)
99013 FORMAT (5x,'e',5x,4x,'Total',4x,5x,'-1',6x,6(5x,i2,6x))
   DO i = mine1 , NE1
      WRITE (plot,99011) e1(i) , ptot(i)/norm , pspect(i,-1)/norm , (pspect(i,l)/norm,pspect(i,-l-1)/norm,l=1,lmax)
   ENDDO
!
   WRITE (plot,99014) (l,-l-1,l=1,lmax)
99014 FORMAT (5x,'e',5x,4x,'Total',5x,'-1',2x,6(5x,i2,2x))
   DO i = mine1 , NE1
      WRITE (plot,99012) e1(i) , ptot(i) , pspect(i,-1) , (pspect(i,l),pspect(i,-l-1),l=1,lmax)
99012 FORMAT (f11.5,8F9.3)
   ENDDO
!
   STOP 'end rccvaug'
!
99001 FORMAT (' *****************************************************')
99005 FORMAT (10A1)
99011 FORMAT (f11.5,8E13.5)
END PROGRAM rccvaug
!*==CCVOPEN.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE ccvopen(Dos,Mat,Tty,Prt,Plot)
   IMPLICIT NONE
!
! query the user or the ctl file for a name to be used to
! construct file names for the job, then open all the
! required files.
!
   CHARACTER*40 name
   INTEGER Dos , Mat , Tty , Prt , Plot
!
! initialize i/o units
!
   Dos = 1
   Mat = 2
   Prt = 3
   Tty = 5
   Plot = 7
!
   OPEN (UNIT=Tty,FILE='rccvaug.in')
!
!            dos file
!
   READ (Tty,99001) name
   OPEN (UNIT=Dos,FILE=name)
!
!            matrix element file
!
   READ (Tty,99001) name
   OPEN (UNIT=Mat,FILE=name)
!
!            printer file
!
   READ (Tty,99001) name
   OPEN (UNIT=Prt,FILE=name)
!
!           spectrum file
!
   READ (Tty,99001) name
   OPEN (UNIT=Plot,FILE=name)
!
99001 FORMAT (a40)
END SUBROUTINE ccvopen
!*==SYNSPEC.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE synspec(E0,Pdos,Ne0,Maxd,Pspect,Stot,Em0,Pmat,Nem0,Maxm,Ldim,Lmax,Prt,Ipr)
!
!     put together xps-spectrum on dos mesh
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 E0 , Em0 , Pdos , Pmat , Pspect , Stot , sum , xmat , ylag
   INTEGER i , iex , Ipr , iw , j , k , Ldim , Lmax , Maxd , Maxm , nc , Ne0 , Nem0
!*** End of declarations inserted by SPAG
!
   DIMENSION E0(Ne0) , Em0(Nem0)
!
   DIMENSION Pdos(Ne0,-Ldim-1:Ldim)
   DIMENSION Pspect(Ne0,-Ldim-1:Ldim) , Stot(Ne0)
   DIMENSION Pmat(Nem0,-Ldim-1:Ldim)
!
   INTEGER Prt
!
   DO k = -Lmax - 1 , Lmax
      IF ( k/=0 ) THEN
         DO i = 1 , Maxd
            xmat = ylag(E0(i),Em0,Pmat(1,k),0,3,Maxm,iex)
            Pspect(i,k) = xmat*Pdos(i,k)
         ENDDO
      ENDIF
   ENDDO
!
   IF ( Ipr>0 ) WRITE (Prt,99001)
!
99001 FORMAT (/' synspec: total unbroadened spectrum'//                                                                            &
             &'   e     intensity      e     intensity                         e     intensity      e     intensity'/)
   DO i = 1 , Maxd
      sum = 0.E0
      DO k = -Lmax - 1 , Lmax
         IF ( k==0 ) THEN
         ENDIF
         sum = sum + Pspect(i,k)
      ENDDO
      Stot(i) = sum
   ENDDO
   IF ( Ipr<=0 ) RETURN
   nc = 4
   iw = Maxd/nc + 1
   DO i = 1 , iw
      k = nc*iw + i
      IF ( k>Maxd ) k = Maxd
      WRITE (Prt,99002) (E0(j),Stot(j),j=i,k,iw)
99002 FORMAT (1x,4(0pf6.3,1pe12.4,3x))
   ENDDO
   RETURN
!
END SUBROUTINE synspec
!*==YLAG.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION ylag(Xi,X,Y,Ind1,N1,Imax,Iex)
!
!     program authors a.a.brooks and e.c.long,
!     computing technology center, union carbide corp.,
!     nuclear div., oak ridge, tenn.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     lagrangian interpolation
!     xi is interpolated entry into x-array
!     n is the order of lagrangian interpolation
!     y is array from which ylag is obtained by interpolation
!     ind is the min-i for x(i).gt.xi ==> if ind=0, x-array is searched
!     imax is max index of x-and y-arrays
!     extrapolation can occur ==> iex=-1 or +1
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 d , p , s , X , xd , Xi , Y
   INTEGER i , Iex , Imax , ind , Ind1 , inl , inu , j , n , N1
!*** End of declarations inserted by SPAG
!
   DIMENSION X(Imax) , Y(Imax)
   ind = Ind1
   n = N1
   Iex = 0
   IF ( n>Imax ) THEN
      n = Imax
      Iex = n
   ENDIF
   IF ( ind>0 ) THEN
      CALL spag_block_2
      RETURN
   ENDIF
   DO j = 1 , Imax
      IF ( Xi<X(j) ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      IF ( Xi==X(j) ) THEN
         CALL spag_block_6
         RETURN
      ENDIF
   ENDDO
   Iex = 1
   CALL spag_block_3
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
      ind = j
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
      IF ( ind<=1 ) Iex = -1
      inl = ind - (n+1)/2
      IF ( inl<=0 ) inl = 1
      inu = inl + n - 1
      IF ( inu<=Imax ) THEN
         CALL spag_block_4
         RETURN
      ENDIF
      CALL spag_block_3
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
      inl = Imax - n + 1
      inu = Imax
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
      s = 0.D0
      p = 1.D0
      DO j = inl , inu
         p = p*(Xi-X(j))
         d = 1.D0
         DO i = inl , inu
            IF ( i/=j ) THEN
               xd = X(j)
            ELSE
               xd = Xi
            ENDIF
            d = d*(xd-X(i))
         ENDDO
!      if(i.gt.55) write(6,*) i,j,d
         s = s + Y(j)/d
      ENDDO
      ylag = s*p
      CALL spag_block_5
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
      RETURN
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
      ylag = Y(j)
      CALL spag_block_5
      RETURN
   END SUBROUTINE spag_block_6
END FUNCTION ylag
!*==BROAD.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE broad(E0,A0,N0,Max,Mine1,Prt,Hwc,Hwv,Hws,E1,A1,Ai,N1,Ef,Cfac,Snorm)
!
!     interpolate to plot-grid and convolute with lorentzian
!     (energy dependent halfwidth ) for life-time broadening
!     and gaussian spectrometer resolution.
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A0 , A1 , afac , ah , Ai , bfac , Cfac , convgau , convlo , delta , delta2 , E0 , E1 , Ef , emin0 , estep , gam , gc ,   &
        & Hwc , hwl
   REAL*8 Hws , Hwv , pi , rate1 , smax , Snorm , xmin , ylag
   INTEGER i , iex , ix0 , Max , Mine1 , N0 , N1 , NGC
!*** End of declarations inserted by SPAG
!
   PARAMETER (NGC=100)
!
   DIMENSION E0(N0) , A0(N0) , E1(N1) , A1(N1) , Ai(N1)
   DIMENSION gc(-NGC:NGC) , ah(-NGC:NGC)
   SAVE emin0
!
   INTEGER Prt
!
   DATA xmin/ - 14.E0/ , estep/0.05E0/
   DATA ix0/281/ , pi/3.141596/
!
   IF ( Snorm<=0.E0 ) THEN
!
!     set fermi energy equal to zero and convert from
!     rydberg to ev
!
      DO i = 1 , Max
         E0(i) = (E0(i)-Ef)*Cfac
      ENDDO
!
!     set up plot-grid and determine min-energy
!
      DO i = 1 , N1
         E1(i) = xmin + (i-1)*estep
      ENDDO
      SPAG_Loop_1_1: DO i = 1 , N1
         IF ( A0(i)>0.E0 ) EXIT SPAG_Loop_1_1
      ENDDO SPAG_Loop_1_1
      emin0 = E0(i)
      SPAG_Loop_1_2: DO i = 1 , N1
         IF ( E1(i)>emin0 ) EXIT SPAG_Loop_1_2
      ENDDO SPAG_Loop_1_2
      Mine1 = i
   ENDIF
!
!
!     interpolate
!
   DO i = Mine1 , N1
      Ai(i) = ylag(E1(i),E0,A0,0,3,Max,iex)
   ENDDO
   DO i = 1 , Mine1 - 1
      Ai(i) = 0.0E0
   ENDDO
   DO i = ix0 + 1 , N1
      Ai(i) = 0.E0
   ENDDO
!
   DO i = Mine1 , N1
!
!     determine halfwidth for life-time broadening (core+core+valence)
!
!     valence band broadening decreases qudratically
!     from hwl at the band bottom up to zero at the fermi energy
!
!     core state broadening is constant
!
      rate1 = E1(i)/E1(1)
      gam = (hwl*rate1*rate1+Hwc+Hwc)/2.0
      IF ( E1(i)>0.E0 ) gam = Hwc
!
!     convolute now
!
      A1(i) = convlo(gam,Ai,i,estep,Mine1,N1,pi)
   ENDDO
!
!     overwrite unbroadened data in ai with life-time
!     broadened spectrum
!
   DO i = Mine1 , N1
      Ai(i) = A1(i)
   ENDDO
!
!     set up gaussian convolution function
!
   afac = -dlog(0.5D0)/(Hws/2.0)**2
   bfac = dsqrt(afac/pi)
   DO i = -NGC , NGC
      delta = i*estep
      delta2 = delta*delta
      gc(i) = dexp(-afac*delta2)*bfac
   ENDDO
!
!     convolute now
!
   DO i = Mine1 , N1
      A1(i) = convgau(Ai,i,gc,ah,estep,NGC,N1)
   ENDDO
!
   IF ( Snorm<=0.E0 ) THEN
!
!     determine maximum and set equal to 100.
!
      smax = 0.E0
      DO i = Mine1 , N1
         IF ( A1(i)>smax ) smax = A1(i)
      ENDDO
      Snorm = 100./smax
   ENDIF
!
!
   DO i = Mine1 , N1
      A1(i) = A1(i)*Snorm
   ENDDO
!
END SUBROUTINE broad
!*==CONVLO.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION convlo(Gam,A1,Ind,Estep,Mine1,N1,Pi)
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A1 , atm , atp , Estep , Gam , Pi , sum
   INTEGER i , Ind , Mine1 , N1
!*** End of declarations inserted by SPAG
!
!     do the convolution integral by trapezoidal rule
!
   DIMENSION A1(N1)
!
   sum = 0.0
   DO i = Mine1 , N1
      atp = datan(Estep*((i-Ind)+0.5)/Gam)
      atm = datan(Estep*((i-Ind)-0.5)/Gam)
      sum = sum + A1(i)*(atp-atm)
   ENDDO
   convlo = sum/Pi
END FUNCTION convlo
!*==CONVGAU.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION convgau(A1,Ind,Gc,Ah,Estep,Ngc,N1)
!
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
   END SUBROUTINE spag_block_2
!
END FUNCTION convgau
