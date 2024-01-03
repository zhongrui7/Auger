PROGRAM spag_program_1
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a0 , cfac , conv , conv1 , e0 , ef , pdos , pi
   INTEGER i , kap , l , lmax , LMAXP , ne , NE0
!*** End of declarations inserted by SPAG
   PARAMETER (NE0=500)
   PARAMETER (LMAXP=3)
!
   CHARACTER*80 filen , iddos
   CHARACTER*1 dosunit
!
   DIMENSION e0(NE0)
   DIMENSION pdos(NE0,-LMAXP-1:LMAXP)
!
   DATA cfac/13.606/
   DATA pi/3.14159265359/
!
   OPEN (UNIT=1,FILE='refdos.in')
   DO
!
!            input dos file
!
      READ (1,99001,END=99999) filen
      OPEN (UNIT=2,FILE=filen)
!
!            output dos file
!
      READ (1,99001) filen
      OPEN (UNIT=8,FILE=filen)
!
      READ (1,*) a0
!
      READ (1,*)
!
      conv1 = 2.*pi/a0
      conv1 = conv1*conv1
!
!
! ***************  read in densities of states  ***************
!
! in relativistic case the ordering of kappa's in the dos file:
!          -1   1  -2   2  -3   3  -4
!
!
      READ (2,99001) iddos
      READ (2,99003) dosunit
99003 FORMAT (a1)
      READ (2,*) ef
      READ (2,*) lmax
      READ (2,*) ne
      READ (2,*)
!
      IF ( dosunit=='d' ) conv = conv1*cfac
      IF ( dosunit=='r' ) conv = cfac
      IF ( dosunit=='e' ) conv = 1.0
!
      DO i = 1 , ne
!
         READ (2,*) e0(i) , pdos(i,0) , pdos(i,-1) , (pdos(i,l),pdos(i,-l-1),l=1,lmax)
!
         e0(i) = (e0(i)-ef)*conv
         DO kap = -lmax - 1 , lmax
            pdos(i,kap) = pdos(i,kap)/cfac
         ENDDO
!
      ENDDO
!
      WRITE (8,99002) iddos
99002 FORMAT (1x,a80)
      WRITE (8,*) ne
!
      DO i = 1 , ne
!
         WRITE (8,99004) e0(i) , (pdos(i,kap),kap=-lmax-1,-1) , (pdos(i,kap),kap=1,lmax) , pdos(i,0)
99004    FORMAT (9F10.5)
!
      ENDDO
!
      CLOSE (2)
      CLOSE (8)
   ENDDO
!
99001 FORMAT (a80)
99999 END PROGRAM spag_program_1
