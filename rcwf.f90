PROGRAM rcore
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 anorm , Bar , bnorm , Char , cwf , Delx , Den , Dfl , Dgc , Dp , Dpas , Dpc , Dq , Dq1 , Dr , Dv , eone , ev , evcor ,   &
        & fpot
   INTEGER i , Icone , Icore , Icoul , Iemode , iend , Iex , imax , imax0 , IORB , ipot , Iscfat , Iskip , Itatom , Itcore , Itot ,&
         & Iturn , Itval , Iunit8 , Iunt14
   INTEGER j , Jturn , k , kk , last , last1 , moat , natom , Nel , Nes , Nitpot , Nk , nkk1 , Nmax , Nnk , NOAT , Norb , Np ,     &
         & Nql , Nqn
   INTEGER NRAD , Nrc , Nrws , Nrws1 , Nstop , Nt , Nuc , Nz
   REAL*8 qval , Ra , Rc , Rhoc , Rhoold , Rhotot , rinf , rlast , rr , rsimp , rsimp3 , Rws , Rws1 , Stval , tail , Tcore , Test ,&
        & Testv , Tets , Thfpot
   REAL*8 Titre , Vc , Vcc , Vew , wa , Xx , Yy , Z , Za
!*** End of declarations inserted by SPAG
!
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
   PARAMETER (IORB=30)
!
   COMMON /bla   / Den(IORB) , Dq1(IORB) , Dfl(IORB) , Nqn(IORB) , Nql(IORB) , Nk(IORB) , Nmax(IORB) , Nel(IORB) , Norb , Icore
   COMMON /dira  / Dv(NRAD) , Dr(NRAD) , Dp(NRAD) , Dq(NRAD) , Dpas , Z , Nstop , Nes , Tets , Np , Nuc
   COMMON /ps2   / Titre(IORB) , Bar(10) , Test , Testv
   COMMON /shoot / Dgc(NRAD,IORB) , Dpc(NRAD,IORB)
   COMMON /ajf   / Icoul
   COMMON /oo    / Xx(NRAD) , Yy(NRAD)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /xmesh / Stval , Delx
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nnk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /atom  / Itot , Iex , Iscfat , Nitpot , Itatom
   COMMON /ay    / Rhotot(NRAD,NOAT) , Char(NRAD,NOAT) , Rhoc(NRAD,NOAT) , Rhoold(NRAD,NOAT)
   COMMON /short / Iskip
   COMMON /tom   / Thfpot(NRAD)
   COMMON /files / Iunit8 , Iunt14
   COMMON /kiva  / Iturn , Itcore , Jturn , Itval , Iemode
   COMMON /ws    / Rws , Nrws(NOAT)
   COMMON /ws1   / Rws1 , Nrws1(NOAT)
   COMMON /eshift/ Vcc(NOAT)
   COMMON /ewald / Vew(NOAT)
   COMMON /sumec / Icone
   COMMON /etot  / Tcore(NOAT)
   DIMENSION cwf(IORB)
   CHARACTER aa*80
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!     ***************************************************************
!     * calculate fully relativistically one-electron orbitals and  *
!     * total energy for free atoms or core one-electron orbitals   *
!     * and core charge density                                     *
!     * calculate kinetic energy contribution from the core         *
!     ***************************************************************
!
         WRITE (*,*) ' input?'
         READ (*,99001) aa
         OPEN (UNIT=3,FILE=aa,STATUS='old')
!
         READ (3,99001) aa
         OPEN (UNIT=8,FILE=aa,STATUS='old')
         READ (3,99001) aa
         OPEN (UNIT=6,FILE=aa,STATUS='unknown')
         READ (3,99001) aa
         OPEN (UNIT=7,FILE=aa,STATUS='unknown')
!
         READ (3,99006) ipot , Iskip , moat
!
         IF ( Iunit8/=1 ) THEN
!
            IF ( ipot==0 ) CALL cpain
            IF ( ipot==1 ) CALL apwin
            IF ( ipot==2 ) CALL kkrin
         ENDIF
!
         natom = moat
!
         CALL insld(natom)
!
         wa = -1.0D00
         Z = Nz(natom)
         Itatom = 0
         IF ( Iscfat==0 ) Norb = Icore
!
         DO i = 1 , NRAD
            rr = dexp(Stval+(i-1)*Delx)
            Rhotot(i,natom) = 0.0D00
            Rhoc(i,natom) = 0.0D00
            Dr(i) = rr
            Thfpot(i) = 0.0D00
            IF ( ipot<0 ) Dv(i) = fpot(rr,Z,wa)
         ENDDO
!
!     add tail to the muffin-tin potential
!
         IF ( ipot>=0 ) THEN
!
            DO i = 1 , NRAD
               Dp(i) = Rhoold(i,natom)
            ENDDO
!
            IF ( Iemode==0 ) THEN
               last = Nrc(natom)
               rlast = Rc(natom)
               qval = rsimp(Dp,Dr,Rc(natom),Nrc(natom))
!
            ELSEIF ( Iemode==2 ) THEN
!
               last = Nrws(natom)
               rlast = Rws
               qval = rsimp3(Dp,Dr,Rws,Nrws(natom))
            ELSE
               last = Nrws1(natom)
               rlast = Rws1
               qval = rsimp3(Dp,Dr,Rws1,Nrws1(natom))
            ENDIF
!
            IF ( Icoul==2 ) qval = 0.0D00
            IF ( Icoul/=0 ) THEN
!
!              ..  V(r)/2  -Vc/2               r < R
!     Vat(r) = |
!     [a.u.]   |
!              |   0                                   icoul = 0
!              |  -(Z-Q)[ 1/r + 1/R ]          r > R   icoul = 1
!              .. -Z[ 1/r + 1/R]                       icoul = 0
!
               DO i = 1 , NRAD
                  Thfpot(i) = -(Nz(natom)-qval)/Dr(i)
               ENDDO
            ENDIF
!
            tail = 0.0D00
!
            IF ( Icoul/=0 ) tail = (Nz(natom)-qval)/rlast
!
            last1 = last + 1
!
            DO i = 1 , last
               Dv(i) = Za(natom,i)/2.0D00 - Vc(natom)/2.D00 - Vcc(natom)/2.0D00
            ENDDO
!
!
            DO i = last1 , NRAD
               Dv(i) = Thfpot(i) + tail
            ENDDO
!
            WRITE (6,'( '' Potential for type'',i2)') natom
            WRITE (6,'(4d20.10)') (Dr(k),Dv(k),k=1,NRAD)
         ENDIF
         DO
!
!     now do core part
!
            DO i = 1 , NRAD
               Xx(i) = 0.0D00
!
               DO j = 1 , Norb
                  Dgc(i,j) = 0.0D00
                  Dpc(i,j) = 0.0D00
               ENDDO
            ENDDO
!
            imax0 = 0
            Np = NRAD
!
            DO j = 1 , Norb
               Tets = Test
               SPAG_Loop_4_1: DO
                  DO kk = 1 , NRAD
                     Dp(kk) = 0.0D00
                     Dq(kk) = 0.0D00
                  ENDDO
!
                  CALL resld(Nqn(j),Nql(j),Nk(j),imax,Den(j),Dfl(j),Dq1(j),j,natom)
!
                  IF ( imax>imax0 ) imax0 = imax
!
                  IF ( Nstop==0 ) THEN
!
!     calculate (core) charge density
!
                     Dp(1) = 0.0D00
                     Dq(1) = 0.0D00
!
                     DO i = 1 , imax
                        Dgc(i,j) = Dp(i)
                        Dpc(i,j) = Dq(i)
                        Xx(i) = Xx(i) + Nel(j)*(Dp(i)*Dp(i)+Dq(i)*Dq(i))
                        IF ( ipot<0 ) Rhotot(i,natom) = Xx(i)
                        IF ( j<=Icore ) Rhoc(i,natom) = Xx(i)
                     ENDDO
                     EXIT SPAG_Loop_4_1
                  ELSE
                     Tets = Test*10.0D00
                  ENDIF
               ENDDO SPAG_Loop_4_1
!
            ENDDO
!
!     solve atomic problem selfconsistently, if wanted
!
            Np = imax0
 
            IF ( Iscfat==0 ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!
!
            CALL poisat(natom)
            IF ( Nstop/=0 ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
!
!     write out the precious results
!
         rinf = Dr(Np-1)
!
         IF ( Iskip/=0 ) THEN
            WRITE (6,99002)
!
99002       FORMAT ('1',10x,'one electron energies'//)
            DO i = 1 , Norb
               WRITE (6,99003) Nqn(i) , Titre(i) , Den(i)
99003          FORMAT (10x,i2,a4,5x,e15.8)
            ENDDO
!
!     orthonormality relations
!
            WRITE (6,99004)
99004       FORMAT (10x//10x,'orthogonality relations'//)
         ENDIF
         DO i = 1 , Norb
            DO j = i , Norb
               IF ( Nql(i)==Nql(j) ) THEN
                  IF ( Nk(i)==Nk(j) ) THEN
!
                     DO k = 1 , Np
                        Yy(k) = Dpc(k,i)*Dpc(k,j) + Dgc(k,i)*Dgc(k,j)
                     ENDDO
!
                     IF ( Iemode==0 ) anorm = rsimp(Yy,Dr,Rc(natom),Nrc(natom))
                     IF ( Iemode==1 ) anorm = rsimp3(Yy,Dr,Rws1,Nrws1(natom))
                     IF ( Iemode==2 ) anorm = rsimp3(Yy,Dr,Rws,Nrws(natom))
!
                     IF ( i==j ) cwf(i) = anorm
!
                     IF ( Iskip/=0 ) THEN
!
                        bnorm = rsimp(Yy,Dr,rinf,Np)
                        WRITE (6,99005) Nqn(i) , Titre(i) , Nqn(j) , Titre(j) , bnorm , anorm
99005                   FORMAT (10x,'(',i2,a4,',',i2,a4,')',3x,2E15.8)
                     ENDIF
                  ENDIF
               ENDIF
!
            ENDDO
         ENDDO
!
!     calculate number of core electrons
!
         DO k = 1 , Np
            Xx(k) = Rhoc(k,natom)
         ENDDO
!
         IF ( Iemode==0 ) anorm = rsimp(Xx,Dr,Rc(natom),Nrc(natom))
         IF ( Iemode==1 ) anorm = rsimp3(Xx,Dr,Rws1,Nrws1(natom))
         IF ( Iemode==2 ) anorm = rsimp3(Xx,Dr,Rws,Nrws(natom))
         IF ( Iturn>=Itcore ) THEN
            IF ( Jturn>=Itval ) THEN
               WRITE (6,99007) anorm
99007          FORMAT (10x//10x,'number of core-electrons in atomic sphere = ',e15.8)
            ENDIF
         ENDIF
!
         IF ( Iscfat==0 ) THEN
!
!     calculate core part of total energy
!
            eone = 0.D00
!
            DO i = 1 , Icore
               anorm = 1.0D00
               IF ( Icone==1 ) anorm = cwf(i)
               eone = eone + anorm*Nel(i)*Den(i)
            ENDDO
            eone = 2.0D00*eone
!
            Dq(1) = 0.0D00
            nkk1 = Nnk(natom)
            DO k = 1 , nkk1
               Dp(k) = (Za(natom,k)-Vc(natom)-Vcc(natom))*Rhoc(k,natom)
            ENDDO
!
            IF ( Iemode==0 ) ev = rsimp(Dp,Dr,Rc(natom),Nrc(natom))
            IF ( Iemode==1 ) ev = rsimp3(Dp,Dr,Rws1,Nrws1(natom))
!
            evcor = 0.0D00
            IF ( Icoul/=0 ) THEN
!
               SPAG_Loop_2_2: DO i = last1 , NRAD
                  IF ( Rhoc(i,natom)==0. ) EXIT SPAG_Loop_2_2
               ENDDO SPAG_Loop_2_2
               iend = i - 1
!
               DO i = last1 , iend
                  evcor = evcor + (Dv(i-1)+Dv(i))*(Dr(i)-Dr(i-1))
               ENDDO
            ENDIF
!
            Tcore(natom) = eone - ev - evcor
!
            IF ( Iskip/=0 ) THEN
               WRITE (6,99008) eone , ev , evcor , Tcore(natom)
99008          FORMAT (15x,'e-one core ',e20.12/15x,'pot -term ',e20.12/15x,' <V-correction> ',e20.12/15x,' <e-c> ',e20.12/)
!
!
!  write out core charge densities
!
!      do j=1,nnk(natom)
!        write(7,6010) dr(j),rhoc(j,natom)
!      end do
               WRITE (7,'(4d20.12)') (Rhoc(j,natom),j=1,Nnk(natom))
!
!
!  write out core wavefunctions
!
!
               DO i = 1 , Norb
                  WRITE (7,99009) Nqn(i) , Titre(i)
99009             FORMAT (i1,a4)
                  WRITE (7,99010) Den(i)
99010             FORMAT (e15.8)
                  WRITE (7,99006) Nk(i)
                  WRITE (7,99006) Nnk(natom)
                  DO j = 1 , Nnk(natom)
                     WRITE (7,99011) Dr(j) , Dgc(j,i) , Dpc(j,i)
99011                FORMAT (3E20.10)
                  ENDDO
               ENDDO
            ENDIF
         ELSE
!
!     calculate total energy for free atom
!
            CALL totec(natom)
         ENDIF
!
!
!
         STOP
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99001 FORMAT (a80)
99006 FORMAT (5I4)
END PROGRAM rcore
!
!*==APWIN.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE apwin
!
!     ******************************
!     * potential input apw type   *
!     ******************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Aa , Conc , Conc1 , Confac , P , R , Ra , Rc , Rws , Rws1 , Stval , Title , vbar , Vc , Xdel , Za
   INTEGER i , ilat , imode , IREP , Iskip , Iunit8 , Iunt14 , j , k , k1 , kmax , kx , n , Nk , NOAT , NRAD , Nrc , Nrws , Nrws1 ,&
         & Nt
   INTEGER Nz
!*** End of declarations inserted by SPAG
 
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
   PARAMETER (IREP=7)
!
   INTEGER u
   COMPLEX*16 Sctamp
   COMMON /short / Iskip
   COMMON /dolly / Sctamp(IREP) , Conc
   COMMON /samfac/ Confac
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /oo    / P(NRAD) , R(NRAD)
   COMMON /xmesh / Stval , Xdel
   COMMON /lat   / Aa , Conc1
   COMMON /tit   / Title
   COMMON /files / Iunit8 , Iunt14
   COMMON /ws    / Rws , Nrws(NOAT)
   COMMON /ws1   / Rws1 , Nrws1(NOAT)
!
   DIMENSION Title(20)
   INTEGER :: spag_nextblock_1
!
   Iunit8 = 1
   imode = 0
!
   READ (8,99001) Title
!
!
99001 FORMAT (20A4)
   READ (8,99002) Aa , Conc , ilat
   READ (8,99002) Stval , Xdel
!
   IF ( Stval==0. ) Stval = -8.8D00
   IF ( Xdel==0. ) Xdel = 0.05D00
!
   Conc1 = Conc
   u = 0
   IF ( Iskip/=0 ) THEN
      WRITE (6,99008) Title
99008 FORMAT ('1',20x,20A4//)
      WRITE (6,99009) Aa , Conc
99009 FORMAT (10x,'a0 = ',f10.5,' conc = ',f10.5/)
   ENDIF
!
!
   DO n = 1 , NOAT
      READ (8,99003) Nz(n) , Nk(n) , Rc(n) , Vc(n)
99003 FORMAT (2I4,2F10.5)
      IF ( Iskip/=0 ) WRITE (6,99010) Nz(n) , Nk(n) , Rc(n) , Vc(n)
99010 FORMAT (10x,'z = ',i4,' nk = ',i4,' rc = ',f10.5,' vc = ',f10.5/)
      kx = Nk(n)
      READ (8,99004) (Ra(n,i),Za(n,i),i=1,kx)
99004 FORMAT (4E20.8)
   ENDDO
!
   CALL fitrad(imode)
!
   DO n = 1 , NOAT
      k1 = Nrc(n)
      Rc(n) = Ra(n,k1)
   ENDDO
!
   CALL mike(ilat)
 
!
   DO n = 1 , NOAT
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
            kmax = Nk(n)
            DO k = 1 , kmax
               IF ( Ra(n,k)>Rws ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            spag_nextblock_1 = 2
         CASE (2)
            Nrws(n) = k - 1
            Nk(n) = k
            DO k = 1 , kmax
               IF ( Ra(n,k)>Rws1 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            spag_nextblock_1 = 3
         CASE (3)
            Nrws1(n) = k - 1
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
   ENDDO
!
   vbar = Conc*Vc(1) + (1-Conc)*Vc(2)
   Vc(1) = vbar
   Vc(2) = vbar
!
   IF ( Iskip/=0 ) THEN
!
      WRITE (6,99005)
99005 FORMAT ('1',10x/10x///'-----new input'///)
      DO n = 1 , NOAT
         kmax = Nk(n)
         WRITE (6,99006)
99006    FORMAT (72x/12x,4('r',14x,'v',14x))
         WRITE (6,99007) (i,(Ra(n,i+j),Za(n,i+j),j=u,3),i=1,kmax,4)
99007    FORMAT (1x,i3,8E15.7)
      ENDDO
   ENDIF
   RETURN
99002 FORMAT (2F10.5,i4)
END SUBROUTINE apwin
!
!*==CPAIN.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE cpain
!
!     ****************************************************************
!     *       input routine for potentials in the cpa mode           *
!     ****************************************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Aa , Conc , Conc1 , Confac , Ra , Rc , Rws , Rws1 , Stval , Title , vbar , Vc , Xdel , Za
   INTEGER i , ii , ilat , IREP , Iskip , Iunit8 , Iunt14 , j , k , kmax , n , Nk , NOAT , NRAD , Nrc , Nrws , Nrws1 , Nt , Nz
!*** End of declarations inserted by SPAG
!
   PARAMETER (IREP=7)
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
!
   COMPLEX*16 Sctamp
!
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /samfac/ Confac
   COMMON /dolly / Sctamp(IREP) , Conc
   COMMON /short / Iskip
   COMMON /lat   / Aa , Conc1
   COMMON /tit   / Title
   COMMON /xmesh / Stval , Xdel
   COMMON /files / Iunit8 , Iunt14
   COMMON /ws    / Rws , Nrws(NOAT)
   COMMON /ws1   / Rws1 , Nrws1(NOAT)
!
!      dimension nsymbl(noat)
   DIMENSION Title(20)
   INTEGER :: spag_nextblock_1
!
   Iunit8 = 1
!
   READ (8,99001) Title
!
99001 FORMAT (20A4)
   READ (8,99002) Aa , Conc , ilat
   READ (8,99002) Stval , Xdel
!
   IF ( Stval==0. ) Stval = -8.8D00
   IF ( Xdel==0. ) Xdel = 0.05D00
!
   Conc1 = Conc
   DO n = 1 , NOAT
      Nt(n) = n
      READ (8,99003) Nz(n) , Nk(n) , Nrc(n) , Rc(n) , Vc(n)
99003 FORMAT (4x,3I4,2F11.7)
      kmax = Nk(n)
      READ (8,99004) (Za(n,i),i=1,kmax)
99004 FORMAT (5D15.8)
!
      DO i = 1 , kmax
         Ra(n,i) = dexp(Stval+Xdel*(i-1))
         Za(n,i) = Za(n,i)/Ra(n,i)
      ENDDO
   ENDDO
!
   CALL mike(ilat)
!
   DO n = 1 , NOAT
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
            kmax = Nk(n)
            DO k = 1 , kmax
               IF ( Ra(n,k)>Rws ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            spag_nextblock_1 = 2
         CASE (2)
            Nrws(n) = k - 1
            Nk(n) = k
            DO k = 1 , kmax
               IF ( Ra(n,k)>Rws1 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            spag_nextblock_1 = 3
         CASE (3)
            Nrws1(n) = k - 1
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
   ENDDO
!
   vbar = Conc*Vc(1) + (1-Conc)*Vc(2)
!
   Vc(1) = vbar
   Vc(2) = vbar
!
   IF ( Iskip>=2 ) THEN
!
      WRITE (6,99005) Title
99005 FORMAT ('1',10x,20A4//)
      ii = 0
      WRITE (6,99006) Aa , Conc , Vc(1)
99006 FORMAT (10x,'lattice constant ',f10.5/10x,'concentration ',f10.5/10x,'muffin tin zero ',f10.5//)
      WRITE (6,99010) Confac
99010 FORMAT (10x,'conversion factor to rydberg ',f10.5/)
      ii = 0
      DO n = 1 , NOAT
         WRITE (6,99007) n
99007    FORMAT ('1',10x,'potential for scatterer ',i4//)
         WRITE (6,99008)
99008    FORMAT (72x/12x,4('r',14x,'v',14x))
         kmax = Nk(n)
         WRITE (6,99009) (i,(Ra(n,i+j),Za(n,i+j),j=ii,3),i=1,kmax,4)
99009    FORMAT (1x,i3,8E15.7)
      ENDDO
   ENDIF
   RETURN
99002 FORMAT (2F10.5,i4)
END SUBROUTINE cpain
!
!*==FPOT.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
FUNCTION fpot(R,Z,Wa)
!
!     *************************************************
!     *  thomas fermi potential at the radial point r *
!     *  z          -  atomic number                  *
!     *  wa         -  ionicity - 1                   *
!     *************************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 fpot , R , Wa , wc , wd , we , Z
!*** End of declarations inserted by SPAG
   wc = dsqrt((R*(Z+Wa)**(1./3.))/0.8853)
   wd = wc*(0.60112*wc+1.81061) + 1.
   we = wc*(wc*(wc*(wc*(0.04793*wc+0.21465)+0.77112)+1.39515)+1.81061) + 1.
   wc = (Z+Wa)*(wd/we)**2 - Wa
   fpot = -wc/R
END FUNCTION fpot
!
!*==INOUH.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE inouh(Dp,Dq,Dr,Dq1,Dfl,Dv,Z,Test,Nuc,Nstop,Jc)
!
!     ****************************************
!     * start values for outward integration *
!     ****************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 dbe , Dd , Dep , Deq , deva1 , deva2 , deva3 , Dfl , Dk , Dm , Dp , Dpno , dpr , Dq , Dq1 , Dqno , dqr , Dr , Dsal , dsum
   REAL*8 Dv , dval , Dvc , Test , Z
   INTEGER i , IORB , j , Jc , m , NRAD , Nstop , Nuc
!*** End of declarations inserted by SPAG
   PARAMETER (IORB=30)
   PARAMETER (NRAD=400)
!
   COMMON /ps1   / Dep(10) , Deq(10) , Dd , Dvc , Dsal , Dk , Dm
   COMMON /trois / Dpno(4,IORB) , Dqno(4,IORB)
   DIMENSION Dp(NRAD) , Dq(NRAD) , Dr(NRAD)
!
   DO i = 1 , 10
      Dp(i) = 0.0D00
      Dq(i) = 0.0D00
   ENDDO
!
   IF ( Nuc<=0 ) THEN
!
      dval = Z/Dvc
      deva1 = -dval
      deva2 = Dv/Dvc + dval/Dr(1) - Dd
      deva3 = 0.0D00
!
      IF ( Dk<=0 ) THEN
!
         dbe = (Dk-Dfl)/dval
      ELSE
         dbe = dval/(Dk+Dfl)
      ENDIF
!
      Dq(10) = Dq1
      Dp(10) = dbe*Dq1
   ELSE
!
      dval = Dv + Z*(3.0D00-Dr(1)*Dr(1)/(Dr(Nuc)*Dr(Nuc)))/(Dr(Nuc)+Dr(Nuc))
      deva1 = 0.0D00
      deva2 = (dval-3.0D00*Z/(Dr(Nuc)+Dr(Nuc)))/Dvc - Dd
      deva3 = Z/(Dr(Nuc)*Dr(Nuc)*Dr(Nuc)*Dsal)
!
      IF ( Dk<=0 ) THEN
!
         Dp(10) = Dq1
      ELSE
!
         Dq(10) = Dq1
      ENDIF
   ENDIF
!
   DO i = 1 , 5
      Dp(i) = Dp(10)
      Dq(i) = Dq(10)
      Dep(i) = Dp(i)*Dfl
      Deq(i) = Dq(i)*Dfl
   ENDDO
!
   m = 1
   DO
      Dm = m + Dfl
      dsum = Dm*Dm - Dk*Dk + deva1*deva1
      dqr = (Dsal-deva2)*Dq(m+9) - deva3*Dq(m+7)
      dpr = deva2*Dp(m+9) + deva3*Dp(m+7)
      dval = ((Dm-Dk)*dqr-deva1*dpr)/dsum
      dsum = ((Dm+Dk)*dpr+deva1*dqr)/dsum
!
      j = -1
!
      DO i = 1 , 5
         dpr = Dr(i)**m
         dqr = dsum*dpr
         dpr = dval*dpr
!
         IF ( m/=1 ) THEN
            IF ( dabs(dpr/Dp(i))<=Test .AND. abs(dqr/Dq(i))<=Test ) j = 1
         ENDIF
!
         Dp(i) = Dp(i) + dpr
         Dq(i) = Dq(i) + dqr
         Dep(i) = Dep(i) + dpr*Dm
         Deq(i) = Deq(i) + dqr*Dm
      ENDDO
!
      IF ( j==1 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      Dp(m+10) = dval
      Dq(m+10) = dsum
      m = m + 1
!
      IF ( m>20 ) THEN
!
         Nstop = 45
         CALL spag_block_1
         RETURN
      ENDIF
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      DO i = 1 , 4
         Dpno(i,Jc) = Dp(i+9)
         Dqno(i,Jc) = Dq(i+9)
      ENDDO
   END SUBROUTINE spag_block_1
END SUBROUTINE inouh
!
!*==INSLD.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE insld(Natom)
!
!     ************************************
!     * input subroutine for core-module *
!     ************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Bar , Delx , Den , Dfl , Dp , Dpas , Dq , Dq1 , Dr , Dv , dval , ominus , oplus , Rc , Stval , Test , Testv , Tets ,     &
        & Titre , Vc
   REAL*8 Z , Za
   INTEGER i , Icore , Icoul , ielec , Iex , IORB , Iscfat , Iskip , Itatom , Itot , j , l , ll , lp1 , Natom , Nel , Nes , nitoe ,&
         & Nitpot , Nk
   INTEGER Nkk , Nmax , NOAT , Norb , Np , Nql , Nqn , NRAD , Nrc , Nstop , Nt , Nuc , Nz
!*** End of declarations inserted by SPAG
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
   PARAMETER (IORB=30)
!
!     z        -  atomic number
!     norb     -  number of orbitals
!     icore    -  number of core orbitals
!     nitoe    -  number of iterations to adjust one electron
!                 energies (default=15)
!     itot     -  : 0 no total energy for atom and core
!                 : 1 total energy for atom and core
!     iex      -  : 0 gunnarsson-lunqvist exchange-correlation
!              -  : 1 kohn-sham exchange
!     iscfat   -  : 0 atomic problem is n o t solved selfconsistently
!                     during rkkr-cpa-scf iterations
!                 : 1 atomic problem is solved selfconsistently
!                     to start rkkr-cpa-scf
!     nitpot   -  number of iterations for the atomic potential,
!                 (default=100)
!     dpas     -  logarithmic increment
!     test     -  precission for one electron energies
!     den      -  one electron energy [a.u.]
!     nqn      -  principal quantum number
!     nk       -  kapa
!     nel      -  orbital occupation
!
   COMMON /bla   / Den(IORB) , Dq1(IORB) , Dfl(IORB) , Nqn(IORB) , Nql(IORB) , Nk(IORB) , Nmax(IORB) , Nel(IORB) , Norb , Icore
   COMMON /dira  / Dv(NRAD) , Dr(NRAD) , Dp(NRAD) , Dq(NRAD) , Dpas , Z , Nstop , Nes , Tets , Np , Nuc
   COMMON /ps2   / Titre(IORB) , Bar(10) , Test , Testv
   COMMON /atom  / Itot , Iex , Iscfat , Nitpot , Itatom
   COMMON /xmesh / Stval , Delx
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nkk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /short / Iskip
   COMMON /ajf   / Icoul
   DIMENSION oplus(4) , ominus(3)
   INTEGER :: spag_nextblock_1
!
   DATA oplus/'s1/2' , 'p3/2' , 'd5/2' , 'f7/2'/
   DATA ominus/'p1/2' , 'd3/2' , 'f5/2'/
!
!
   Dpas = Delx
   Nes = 150
   Testv = 1.E-07
   Test = 1.E-10
   Nuc = 0
   dval = 0.0D00
!
   READ (3,99001) (Bar(i),i=1,10)
!
99001 FORMAT (20A4)
   READ (3,99002) Norb , Icore , nitoe , Itot , Iex , Iscfat , Nitpot , Icoul
99002 FORMAT (8I4)
!
   IF ( nitoe==0 ) Nes = nitoe
   IF ( Nitpot==0 ) Nitpot = 100
   Z = Nz(Natom)
   Np = NRAD
   Itatom = 0
!
   IF ( Iskip/=0 ) THEN
      IF ( Icoul==0 ) WRITE (6,99008) Natom
99008 FORMAT (10x,'atom ',i4,' no tail'/)
      IF ( Icoul==1 ) WRITE (6,99009) Natom
99009 FORMAT (10x,'atom ',i4,' screened coulomb tail'/)
      IF ( Icoul==2 ) WRITE (6,99010) Natom
99010 FORMAT (10x,'atom ',i4,' unscreened coulomb tail'/)
      WRITE (6,99003) Natom
99003 FORMAT ('1',' atomic calculation for scatterer',i4//)
      WRITE (6,99006) Norb , Icore
99006 FORMAT (10x,'number of orbitals:',i4/10x,'number of core orbitals:',i4//)
   ENDIF
!
!     read orbital information
!
   DO i = 1 , Norb
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
!
            READ (3,99004) Den(i) , Nqn(i) , Nk(i) , Nel(i)
99004       FORMAT (e15.8,3I4)
!
            Nql(i) = iabs(Nk(i))
            IF ( Nk(i)<0 ) Nql(i) = Nql(i) - 1
            Dfl(i) = Nk(i)*Nk(i)
            Dfl(i) = dsqrt(Dfl(i)-dval)
            IF ( Nk(i)>=0 ) THEN
!
!     j  = l - 1/2
!
               DO ll = 1 , 3
                  IF ( ll==Nk(i) ) THEN
                     Titre(i) = ominus(ll)
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ENDDO
            ENDIF
!
!     j = l + 1/2
!
            DO ll = 1 , 4
               lp1 = -ll
               IF ( lp1==Nk(i) ) THEN
                  Titre(i) = oplus(ll)
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            spag_nextblock_1 = 2
         CASE (2)
            IF ( Iskip/=0 ) WRITE (6,99005) Nqn(i) , Nk(i) , Titre(i) , Nel(i) , Den(i)
99005       FORMAT (2I6,5x,a4,i5,e15.8)
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
   ENDDO
!
   ielec = 0
   DO i = 1 , Norb
      ielec = ielec + Nel(i)
      Nmax(i) = Np
      l = 1
      j = Nqn(i) - Nql(i)
      IF ( (j-2*(j/2))==0 ) l = -l
      Dq1(i) = l*Nk(i)/iabs(Nk(i))
   ENDDO
   IF ( Iscfat/=0 ) THEN
      IF ( ielec/=Nz(Natom) ) THEN
         WRITE (6,99007) Nz(Natom) , ielec
99007    FORMAT (10x,'get the input right, dummy!'/10x,' z = ',i4,' number of electrons ',i4//)
         STOP 'input screwed up'
      ENDIF
   ENDIF
!
   RETURN
END SUBROUTINE insld
!
!*==INTH.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE inth(Dp,Dq,Dv,Dr)
!
!     *****************************************************
!     * integration of the dirac equation using a 5 point *
!     * adams method                                      *
!     *****************************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Db , Dep , Deq , Dk , dkoef1 , dkoef2 , Dm , Dp , dpr , Dq , dqr , Dr , Dsal , dsum , Dv , Dvc
   INTEGER i
!*** End of declarations inserted by SPAG
   COMMON /ps1   / Dep(10) , Deq(10) , Db , Dvc , Dsal , Dk , Dm
!     dpas    -    delx, logarithmic increment for radius
!     dm      -    dpas/720
!     dkoef1  -    475/502
!     dkoef2  -    27/502
!
   dkoef1 = 475.0D00/502.0D00
   dkoef2 = 27.0D00/502.0D00
!
   dpr = Dp + Dm*((251.0D00*Dep(1)+2616.0D00*Dep(3)+1901.0D00*Dep(5))-(1274.0D00*Dep(2)+2774.0D00*Dep(4)))
   dqr = Dq + Dm*((251.0D00*Deq(1)+2616.0D00*Deq(3)+1901.0D00*Deq(5))-(1274.0D00*Deq(2)+2774.0D00*Deq(4)))
!
   DO i = 2 , 5
      Dep(i-1) = Dep(i)
      Deq(i-1) = Deq(i)
   ENDDO
!
   dsum = (Db-Dv/Dvc)*Dr
   Dep(5) = -Dk*dpr + (Dsal*Dr+dsum)*dqr
   Deq(5) = Dk*dqr - dsum*dpr
!
   Dp = Dp + Dm*((106.0D00*Dep(2)+646.0D00*Dep(4)+251.0D00*Dep(5))-(19.0D00*Dep(1)+264.0D00*Dep(3)))
   Dq = Dq + Dm*((106.0D00*Deq(2)+646.0D00*Deq(4)+251.0D00*Deq(5))-(19.0D00*Deq(1)+264.0D00*Deq(3)))
!
   Dp = dkoef1*Dp + dkoef2*dpr
   Dq = dkoef1*Dq + dkoef2*dqr
   Dep(5) = -Dk*Dp + (Dsal*Dr+dsum)*Dq
   Deq(5) = Dk*Dq - dsum*Dp
!
END SUBROUTINE inth
!
!*==KKRIN.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE kkrin
!
!     ******************************************************************
!     *                                                                *
!     * read in kkr-type potential input                               *
!     *                                                                *
!     ******************************************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Aa , Avol , Conc , Conc1 , Confac , dps , Ec , Ra , Rc , rtemp , Rws , Rws1 , Stval , Title , vbar , Vc , Vcc , vtemp ,  &
        & x0 , Xdel
   REAL*8 xrc , Za
   INTEGER iemode , iform , Ilat , iold , IREP , Iunit8 , Iunt14 , k , kk , kmax , n , Nk , nnc , nnk , NOAT , NRAD , Nrc , Nrws , &
         & Nrws1 , Nt
   INTEGER Nz
!*** End of declarations inserted by SPAG
!
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
   PARAMETER (IREP=7)
!
   COMPLEX*16 Sctamp
   COMMON /dolly / Sctamp(IREP) , Conc
   COMMON /samfac/ Confac
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /lat   / Aa , Conc1
   COMMON /xmesh / Stval , Xdel
   COMMON /tit   / Title
   COMMON /files / Iunit8 , Iunt14
   COMMON /ws    / Rws , Nrws(NOAT)
   COMMON /ws1   / Rws1 , Nrws1(NOAT)
   COMMON /eshift/ Vcc(NOAT)
   COMMON /kitty / Avol , Ec , Ilat
!
   DIMENSION Title(20) , rtemp(3) , vtemp(3)
   INTEGER :: spag_nextblock_1
!
   WRITE (6,*) ' enter kkrin!'
   Iunit8 = 1
!
   READ (8,99001) Title
   READ (8,*) Aa , Conc , Ilat , iform , iemode , iold
   READ (8,*) Stval , Xdel
!
   WRITE (6,99001) Title
   WRITE (6,*) Aa , Conc , Ilat , iform , iemode , iold
   WRITE (6,*) Stval , Xdel
!
   IF ( Stval==0. ) Stval = -8.8D00
   IF ( Xdel==0. ) Xdel = 0.05D00
!
   Conc1 = Conc
!
   CALL mike(Ilat)
!
   DO n = 1 , NOAT
      READ (8,*) Nt(n) , Nz(n) , Nk(n) , Nrc(n) , Rc(n) , Vc(n) , Vcc(n)
      WRITE (6,*) Nt(n) , Nz(n) , Nk(n) , Nrc(n) , Rc(n) , Vc(n) , Vcc(n)
      nnk = Nk(n)
      nnc = Nrc(n)
!
!     generate radius and read in potential
!
      IF ( iform==1 ) READ (8,99004) (Za(n,k),k=1,nnc)
99004 FORMAT (4D15.8)
      IF ( iform==0 ) READ (8,99005) (Za(n,k),k=1,nnc)
99005 FORMAT (4D20.12)
      IF ( iform==2 ) READ (8,*) (Za(n,k),k=1,nnc)
!
      xrc = dlog(Rc(n))
      x0 = xrc - (nnc-1)*Xdel
      Stval = x0
      IF ( iold<=0 ) THEN
!
         IF ( iemode==0 ) Xdel = (-Stval+dlog(Rws1))/(nnk-6)
      ENDIF
!
      DO k = 1 , nnk
         Ra(n,k) = dexp(x0)
         Za(n,k) = Za(n,k)/Ra(n,k)
         x0 = x0 + Xdel
      ENDDO
      IF ( iold<2 ) THEN
         DO k = nnc + 1 , nnk
            DO kk = 1 , 3
               rtemp(kk) = Ra(n,k-4+kk)
               vtemp(kk) = Za(n,k-4+kk)
            ENDDO
            CALL interp(rtemp,vtemp,3,Ra(n,k),Za(n,k),dps,.FALSE.)
         ENDDO
      ENDIF
   ENDDO
!
   vbar = Conc*Vc(1) + (1-Conc)*Vc(2)
   Vc(1) = vbar
   Vc(2) = vbar
!
   DO n = 1 , NOAT
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
            kmax = Nk(n)
            DO k = 1 , kmax
               IF ( Ra(n,k)>Rws ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            spag_nextblock_1 = 2
         CASE (2)
            Nrws(n) = k - 1
!
            IF ( iemode==0 ) THEN
!
               Nrws1(n) = Nrc(n)
               Rws1 = Rc(n)
            ELSE
!
               Nk(n) = k
               SPAG_Loop_3_1: DO k = 1 , kmax
                  IF ( Ra(n,k)>Rws1 ) EXIT SPAG_Loop_3_1
               ENDDO SPAG_Loop_3_1
               Nrws1(n) = k - 1
            ENDIF
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
!
   ENDDO
!
   WRITE (6,*) ' kkrin left!'
   RETURN
!
99001 FORMAT (20A4)
99002 FORMAT (f11.6,f9.5,4I4)
99003 FORMAT (4I4,3F11.7)
END SUBROUTINE kkrin
!
!*==POISAT.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE poisat(Natom)
!
!     ***************************************
!     *    calculation of atomic potential  *
!     ***************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 alpha , Bar , beta , Char , D , del0 , del1 , dl0 , Dp , Dpas , Dq , Dr , Dv , pi , pi4 , Ra , Rc , Rhoc , Rhoold , rhos
   REAL*8 Rhotot , rhox , rs , s1 , s2 , s3 , Test , Testv , Tets , third , Titre , twoth , Vc , vex , vexks , Vold , xsl , Z , Za
   INTEGER i , Iex , IORB , Iscfat , Iskip , Itatom , Itot , l , mft , Natom , Nes , Nitpot , Nk , NOAT , Np , np1 , np11 , NRAD , &
         & Nrc , Nstop
   INTEGER Nt , Nuc , Nz
!*** End of declarations inserted by SPAG
!
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
   PARAMETER (IORB=30)
!
   COMMON /atom  / Itot , Iex , Iscfat , Nitpot , Itatom
   COMMON /dira  / Dv(NRAD) , Dr(NRAD) , Dp(NRAD) , Dq(NRAD) , Dpas , Z , Nstop , Nes , Tets , Np , Nuc
   COMMON /ps2   / Titre(IORB) , Bar(10) , Test , Testv
   COMMON /oo    / D(NRAD) , Vold(NRAD)
   COMMON /ay    / Rhotot(NRAD,NOAT) , Char(NRAD,NOAT) , Rhoc(NRAD,NOAT) , Rhoold(NRAD,NOAT)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /short / Iskip
   DIMENSION s1(NRAD) , s2(NRAD)
!
!     set functions for gunnarsson-lundquvist exchange-correlation
!     phys.rev.b13, 4274 (1976)
!
   beta(rs) = 1.0D00 + 0.0545D00*rs*dlog(1.0D00+11.4D00/rs)
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!
         pi = datan(1.0D00)*4.0D00
         pi4 = 4.0D00*pi
         third = 1.0D00/3.0D00
         twoth = 2.0D00*third
!
         xsl = -6.0D00*(3.0D00/(8.0D00*pi))**third
!
!     coulomb potential with poisson equation
!
         DO i = 1 , Np
            Vold(i) = Dv(i)
            s1(i) = 0.0D00
            s2(i) = 0.0D00
         ENDDO
!
         DO i = 3 , Np
            IF ( Rhotot(i,Natom)==0.0 ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
         np1 = i - 1
!
         DO i = 1 , np1
            s1(i) = Rhotot(i,Natom)/Dr(i)
         ENDDO
!
         mft = np1 - 1
         DO i = 2 , mft
            l = mft - i + 1
            s2(l) = s2(l+1) + (s1(l)+s1(l+1))*(Dr(l+1)-Dr(l))
         ENDDO
!
         s3 = Rhotot(1,Natom)*Dr(1)
         Za(Natom,1) = s2(1) + (s3-2.0D00*Z)/Dr(1)
!
         DO i = 2 , np1
            s3 = s3 + (Rhotot(i,Natom)+Rhotot(i-1,Natom))*(Dr(i)-Dr(i-1))
            Za(Natom,i) = s2(i) + (s3-2.0D00*Z)/Dr(i)
         ENDDO
!
         np11 = np1 + 1
         DO i = np11 , Np
            Za(Natom,i) = (s3-2.0D00*Z)/Dr(i)
         ENDDO
!
!     exchange-correlation potential
!
         DO i = 2 , np1
            rhos = Rhotot(i,Natom)/(Dr(i)*Dr(i)*pi4)
            rhox = rhos**third
            rs = (3.0D00/(pi4*rhos))**third
            vexks = twoth*xsl*rhox
            vex = vexks*beta(rs)
            IF ( Iex/=0 ) vex = vexks
            Za(Natom,i) = Za(Natom,i) + vex
         ENDDO
!
         DO i = 1 , Np
            Dv(i) = Za(Natom,i)/2.0D00
            dl0 = -1.0D00/Dr(i)
            IF ( Dv(i)>dl0 ) Dv(i) = dl0
         ENDDO
!
!     check potential
!
         Nstop = 1
         Itatom = Itatom + 1
         IF ( Itatom<=Nitpot ) THEN
!
            del0 = dabs((Dv(1)-Vold(1))/Dv(1))
            DO i = 2 , Np
               del1 = dabs((Dv(i)-Vold(i))/Dv(i))
               IF ( del1>del0 ) del0 = del1
            ENDDO
!
            alpha = 0.75
!
            IF ( Iskip/=0 ) WRITE (6,99001) Itatom , del0
!
            IF ( del0>Testv ) Nstop = 0
            DO i = 1 , Np
               Dv(i) = alpha*Vold(i) + (1.-alpha)*Dv(i)
            ENDDO
            IF ( Nstop/=0 ) WRITE (6,99001) Itatom , del0
         ENDIF
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99001 FORMAT (10x,'iteration ',i4,' delta V-max ',e15.8)
!
END SUBROUTINE poisat
!
!*==RESLD.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE resld(Nqn,Nql,Nk,Imax,De,Dfl,Dq1,Jc,Natom)
!
!     *************************
!     * dirac equation [a.u.] *
!     *************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Db , dbe , dd , De , Dep , Deq , Dfl , Dk , dkoef , Dm , Dp , Dpas , dpm , Dpno , dpq , Dq , Dq1 , dqm , Dqno , Dr
   REAL*8 Dsal , dsum , Dv , dval , Dvc , elim , epriv , Ra , Rc , Test , val , Vc , Z , Za
   INTEGER i , ies , Iex , imat , Imax , imm , IORB , Iscfat , Iskip , Itatom , Itot , j , Jc , k , lll , m , Natom , nd , Nes ,   &
         & Nitpot
   INTEGER Nk , Nnk , NOAT , nodes , Np , Nql , Nqn , NRAD , Nrc , Nstop , Nt , Nuc , Nz
!*** End of declarations inserted by SPAG
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
   PARAMETER (IORB=30)
   COMMON /dira  / Dv(NRAD) , Dr(NRAD) , Dp(NRAD) , Dq(NRAD) , Dpas , Z , Nstop , Nes , Test , Np , Nuc
   COMMON /ps1   / Dep(10) , Deq(10) , Db , Dvc , Dsal , Dk , Dm
   COMMON /trois / Dpno(4,IORB) , Dqno(4,IORB)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nnk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /atom  / Itot , Iex , Iscfat , Nitpot , Itatom
   COMMON /short / Iskip
   INTEGER :: spag_nextblock_1
   INTEGER :: spag_nextblock_2
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!
         dkoef = 1.0D00/720.0D00
!
         Nstop = 0
         Dvc = 137.0373
         Dsal = Dvc + Dvc
!
         epriv = De
         imm = 0
         ies = 0
         Dk = Nk
         lll = (Nql*(Nql+1))/2
         nd = 0
         nodes = Nqn - Nql
         IF ( lll/=0 ) THEN
            elim = Dv(1) + lll/(Dr(1)*Dr(1))
!
            DO i = 2 , Np
               val = Dv(i) + lll/(Dr(i)*Dr(i))
               IF ( val<=elim ) elim = val
            ENDDO
!
            IF ( elim>=0 ) THEN
!
               Nstop = 17
               WRITE (6,99001) Nstop
99001          FORMAT (5x,'nstop = ',i4,'  2*v+l*(l+1)/r**2 is positive'/)
               STOP 'in resld'
            ENDIF
         ELSE
!
!
            elim = -Z*Z/(1.5D00*Nqn*Nqn)
         ENDIF
!
         IF ( Iskip>=2 ) WRITE (6,'('' de='',d20.10,''elim='',d20.10)') De , elim
         IF ( De<=elim ) De = elim*0.5D00
         spag_nextblock_1 = 2
      CASE (2)
         DO WHILE ( imm/=1 )
            spag_nextblock_2 = 1
            SPAG_DispatchLoop_2: DO
               SELECT CASE (spag_nextblock_2)
               CASE (1)
!
                  DO i = 7 , Np , 2
                     imat = Np + 1 - i
                     IF ( (Dv(imat)+lll/(Dr(imat)*Dr(imat))-De)<=0. ) THEN
                        spag_nextblock_2 = 2
                        CYCLE SPAG_DispatchLoop_2
                     ENDIF
                  ENDDO
                  spag_nextblock_2 = 2
               CASE (2)
                  IF ( imat>5 ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
!
                  De = De*0.5D00
                  IF ( De>=-Test .OR. nd>nodes ) THEN
                     Nstop = 28
                     WRITE (6,99002) Nstop
99002                FORMAT (5x,'nstop = ',i4,' 2*v+l*(l+1)/r**2-2*e is positive'/)
                     RETURN
                  ENDIF
                  EXIT SPAG_DispatchLoop_2
               END SELECT
            ENDDO SPAG_DispatchLoop_2
         ENDDO
         spag_nextblock_1 = 3
      CASE (3)
!
!     get start values for outward integration
!
         IF ( Iskip>=2 ) WRITE (6,'('' de='',d20.10)') De
         Db = De/Dvc
         CALL inouh(Dp,Dq,Dr,Dq1,Dfl,Dv(1),Z,Test,Nuc,Nstop,Jc)
         IF ( Nstop/=0 ) THEN
            Nstop = 36
            WRITE (6,99003) Nstop
99003       FORMAT (5x,'nstop = ',i4,' dexpansion at the origin does not converge'/)
            RETURN
         ELSE
!
!     calculate number of nodes for the large component
!
            nd = 1
!
            DO i = 1 , 5
               dval = Dr(i)**Dfl
               IF ( i/=1 ) THEN
                  IF ( Dp(i-1)/=0.0D00 ) THEN
                     IF ( (Dp(i)/Dp(i-1))<=0.0D00 ) nd = nd + 1
                  ENDIF
               ENDIF
!
               Dp(i) = Dp(i)*dval
               Dq(i) = Dq(i)*dval
               Dep(i) = Dep(i)*dval
               Deq(i) = Deq(i)*dval
            ENDDO
            k = -1 + 2*(nodes-2*(nodes/2))
            IF ( (Dp(1)*k)>0. ) THEN
!
               IF ( (k*Nk*Dq(1))>=0. ) THEN
                  Dm = Dpas*dkoef
!
!     do now outward integration
!
                  DO i = 6 , imat
                     Dp(i) = Dp(i-1)
                     Dq(i) = Dq(i-1)
!
                     CALL inth(Dp(i),Dq(i),Dv(i),Dr(i))
!
                     IF ( Dp(i-1)/=0.0D00 ) THEN
                        IF ( (Dp(i)/Dp(i-1))<=0.0D00 ) THEN
                           nd = nd + 1
                           IF ( nd>nodes ) THEN
                              spag_nextblock_1 = 4
                              CYCLE SPAG_DispatchLoop_1
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
!
                  IF ( nd==nodes ) THEN
!
!     start values for inward integration, use descleaux's
!     magic boundary: 300
!
                     dqm = Dq(imat)
                     dpm = Dp(imat)
                     IF ( imm/=1 ) THEN
!
                        SPAG_Loop_2_1: DO i = 1 , Np , 2
                           Imax = Np + 1 - i
                           IF ( ((Dv(Imax)-De)*Dr(Imax)*Dr(Imax))<=300. ) EXIT SPAG_Loop_2_1
                        ENDDO SPAG_Loop_2_1
                     ENDIF
!
                     dd = dsqrt(-De*(2.0D00+Db/Dvc))
!
                     dpq = -dd/(Dsal+Db)
                     Dm = -Dm
                     DO i = 1 , 5
                        j = Imax + 1 - i
                        Dp(j) = dexp(-dd*Dr(j))
                        Dep(i) = -dd*Dp(j)*Dr(j)
                        Dq(j) = dpq*Dp(j)
                        Deq(i) = dpq*Dep(i)
                     ENDDO
                     m = Imax - 5
!
!     do now inward integration
!
                     DO i = imat , m
                        j = m + imat - i
                        Dp(j) = Dp(j+1)
                        Dq(j) = Dq(j+1)
                        CALL inth(Dp(j),Dq(j),Dv(j),Dr(j))
                     ENDDO
!
!     check left and right large components
!
                     dval = dpm/Dp(imat)
                     IF ( dval>0. ) THEN
!
                        DO i = imat , Imax
                           Dp(i) = Dp(i)*dval
                           Dq(i) = Dq(i)*dval
                        ENDDO
!
!     calculate the norm
!
                        dsum = 0.0D00
                        IF ( Dp(1)/=0.0 ) dsum = 3.0D00*Dr(1)*(Dp(1)**2+Dq(1)**2)/(Dpas*(Dfl+Dfl+1.))
!
                        DO i = 3 , Imax , 2
                           dsum = dsum + Dr(i)*(Dp(i)**2+Dq(i)**2) + 4.*Dr(i-1)*(Dp(i-1)**2+Dq(i-1)**2) + Dr(i-2)                  &
                                & *(Dp(i-2)**2+Dq(i-2)**2)
                        ENDDO
                        dsum = Dpas*(dsum+Dr(imat)*(dqm*dqm-Dq(imat)*Dq(imat)))/3.0D00
!
!     modify one-electron energy
!
                        dbe = Dp(imat)*(dqm-Dq(imat))*Dvc/dsum
                        imm = 0
                        val = dabs(dbe/De)
                        IF ( val<=Test ) THEN
                           spag_nextblock_1 = 5
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        DO
                           dval = De + dbe
!
                           IF ( dval<0. ) THEN
!
                              De = dval
                              IF ( dabs(De-epriv)<Test ) THEN
                                 spag_nextblock_1 = 5
                                 CYCLE SPAG_DispatchLoop_1
                              ENDIF
                              epriv = De
                              IF ( Iskip==2 ) WRITE (6,99004) Jc , ies , Test , epriv , De
99004                         FORMAT (' orbital = ',i3,'iteration ',i3,' test',d13.5,3x,'e0 = ',e15.8,3x,'e1 = ',e15.8)
                              IF ( val<=0.1 ) imm = 1
                              ies = ies + 1
                              IF ( ies<=Nes ) THEN
                                 spag_nextblock_1 = 2
                                 CYCLE SPAG_DispatchLoop_1
                              ENDIF
                              Nstop = 362
                              WRITE (6,99005) Nstop
99005                         FORMAT (5x,'nstop = ',i4,'number of iterations is too large'/)
                              WRITE (6,*) ies , Nes
                              STOP
                           ELSE
                              dbe = dbe*0.5D00
                              val = val*0.5D00
                              IF ( val<=Test ) THEN
                                 Nstop = 345
                                 WRITE (6,99006) Nstop
99006                            FORMAT (5x,'nstop = ',i4,' energy converged to zero'/)
                                 RETURN
                              ENDIF
                           ENDIF
                        ENDDO
                     ELSE
                        Nstop = 312
                        WRITE (6,99007) Nstop
99007                   FORMAT (5x,'nstop = ',i4,' sign error for the large component'/)
                        RETURN
                     ENDIF
                  ELSE
                     De = 0.8*De
                     IF ( De<-Test ) THEN
                        spag_nextblock_1 = 2
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     Nstop = 206
                     WRITE (6,99008) Nstop
99008                FORMAT (5x,'nstop = ',i4,' number of nodes is too small'/)
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
            Nstop = 53
            WRITE (6,99009) Nstop
99009       FORMAT (5x,'nstop ',i4,' dexpansion error at the origin'/)
            RETURN
         ENDIF
      CASE (4)
!
         De = 1.2*De
         IF ( De>elim ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         Nstop = 210
         WRITE (6,99010) Nstop
99010    FORMAT (5x,'nstop = ',i4,' number of nodes is too large'/)
         RETURN
      CASE (5)
!      return
!
!     renormalize wavefunction
!
         dsum = dsqrt(dsum)
         Dq1 = Dq1/dsum
!
         DO i = 1 , Imax
            Dp(i) = Dp(i)/dsum
            Dq(i) = Dq(i)/dsum
         ENDDO
!
         DO i = 1 , 4
            Dpno(i,Jc) = Dpno(i,Jc)/dsum
            Dqno(i,Jc) = Dqno(i,Jc)/dsum
         ENDDO
!
         IF ( Imax/=Np ) THEN
            j = Imax + 1
!
            DO i = j , Np
               Dp(i) = 0.0D00
               Dq(i) = 0.0D00
            ENDDO
         ENDIF
!
         Nstop = 0
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE resld
!
!*==TOTEC.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE totec(Natom)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 Bar , beta , Char , Den , Dfl , Dp , Dpas , Dq , Dq1 , Dr , Dv , eatom , ecore , enexc , envxc , eone , eps , epsxc ,    &
        & pi , pi4
   REAL*8 Ra , Rc , rho , Rhoc , Rhoold , Rhotot , rhox , rinf , rinf1 , rs , rsimp , Rws , Rws1 , t , Tcore , Test , Testv ,      &
        & Tets , third , Title
   REAL*8 Titre , u , u1 , u2 , ucore , Vc , Vcc , Vew , vxc , xks , Xx , Yy , Z , Za
   INTEGER i , Icone , Icore , Iemode , Iex , inf , inf1 , IORB , iprint , Iscfat , Iskip , Itatom , Itcore , Itot , Iturn ,       &
         & Itval , Jturn , Natom , Nel , Nes
   INTEGER Nitpot , Nk , Nkk , Nmax , NOAT , Norb , Np , Nql , Nqn , NRAD , Nrc , Nrws , Nrws1 , Nstop , Nt , Nuc , Nz
!*** End of declarations inserted by SPAG
!
   PARAMETER (NOAT=2)
   PARAMETER (IORB=30)
   PARAMETER (NRAD=400)
!
!     ***********************************************
!     *  total energy for free atoms                *
!     *        in  [ryd]   !!!                      *
!     ***********************************************
!
   COMMON /oo    / Xx(NRAD) , Yy(NRAD)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /bla   / Den(IORB) , Dq1(IORB) , Dfl(IORB) , Nqn(IORB) , Nql(IORB) , Nk(IORB) , Nmax(IORB) , Nel(IORB) , Norb , Icore
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /tit   / Title(20)
   COMMON /ps2   / Titre(IORB) , Bar(10) , Test , Testv
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nkk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /ay    / Rhotot(NRAD,NOAT) , Char(NRAD,NOAT) , Rhoc(NRAD,NOAT) , Rhoold(NRAD,NOAT)
   COMMON /etot  / Tcore(NOAT)
   COMMON /atom  / Itot , Iex , Iscfat , Nitpot , Itatom
   COMMON /short / Iskip
   COMMON /ws    / Rws , Nrws(NOAT)
   COMMON /ws1   / Rws1 , Nrws1(NOAT)
   COMMON /kiva  / Iturn , Itcore , Jturn , Itval , Iemode
   COMMON /sumec / Icone
   COMMON /eshift/ Vcc(NOAT)
   COMMON /ewald / Vew(NOAT)
!
   COMMON /dira  / Dv(NRAD) , Dr(NRAD) , Dp(NRAD) , Dq(NRAD) , Dpas , Z , Nstop , Nes , Tets , Np , Nuc
!
!
!     set functions for local density functional:
!     gunnarson-lundqvist, phys.rev.b13, 4274 (1976)
!
   beta(rs) = 1.0D00 + 0.0545D00*rs*dlog(1.0D00+11.4D00/rs)
   eps(rs) = -0.0666D00*((1.0D00+(rs/11.4D00)**3)*dlog(1.0D00+11.4D00*rs)+0.5D00*rs/11.4D00-(rs/11.4D00)**2-1.0D00/3.0D00)
   INTEGER :: spag_nextblock_1
   INTEGER :: spag_nextblock_2
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         pi = datan(1.0D00)*4.0D00
         pi4 = 4.0D00*pi
         third = 1.0D00/3.0D00
         xks = -4.0D00*(3.0D00/(8.0D00*pi))**third
!
!     ****************************************************
!     *      atomic reference total energy               *
!     ****************************************************
!
         WRITE (6,99002) Natom
!
99002    FORMAT (10x///10x,'separated atom total energy for scatterer        [ryd]',i4/)
         IF ( Iex/=0 ) WRITE (6,99003)
99003    FORMAT (10x,'kohn-sham exchange'//)
         IF ( Iex==0 ) WRITE (6,99004)
99004    FORMAT (10x,'hedin-lundqvist exchange-correlation'/)
!
         DO i = 3 , Np
            IF ( Rhotot(i,Natom)==0. ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
         inf1 = i - 1
         rinf1 = Dr(inf)
         inf = Np
         rinf = Dr(inf)
!
         eone = 0.0D00
         DO i = 1 , Norb
            eone = eone + Nel(i)*Den(i)
         ENDDO
         eone = eone*2.0D00
!
         ecore = 0.0D00
         DO i = 1 , Icore
            ecore = ecore + Nel(i)*Den(i)
         ENDDO
         ecore = 2.0D00*ecore
!
         DO i = 1 , Np
            Dq(i) = Dv(i)*Rhoc(i,Natom)
            Dp(i) = Dv(i)*Rhotot(i,Natom)
         ENDDO
!
         ucore = 2.0D00*rsimp(Dq,Dr,rinf,inf)
         u1 = 2.0D00*rsimp(Dp,Dr,rinf,inf)
         t = eone - u1
         Tcore(Natom) = ecore - ucore
!
 
         DO i = 2 , Np
            Dp(i) = Rhotot(i,Natom)/Dr(i)
         ENDDO
!
         u2 = -2.0D00*Nz(Natom)*rsimp(Dp,Dr,rinf,inf)
!
!     exchange-correlation terms
!
         iprint = 0
!
         DO i = 2 , inf
            spag_nextblock_2 = 1
            SPAG_DispatchLoop_2: DO
               SELECT CASE (spag_nextblock_2)
               CASE (1)
                  rho = Rhotot(i,Natom)/(pi4*Dr(i)*Dr(i))
                  rhox = rho**third
                  rs = (3.0D00/(pi4*rho))**third
!     if(dr(i).gt.rc(natom)) go to 53
                  IF ( rs>=9.0 ) THEN
                     IF ( iprint==0 ) WRITE (6,99001) Dr(i) , rs
99001                FORMAT (10x,'cut-off at =',f10.5,' rs = ',e10.5//)
                     iprint = 1
!
                  ELSEIF ( Iex==0 ) THEN
!
!     gunnarsson-lundquvist exchange-correlation
!
                     vxc = xks*rhox*beta(rs)
                     epsxc = 0.75D00*xks*rhox + eps(rs)
                     spag_nextblock_2 = 2
                     CYCLE SPAG_DispatchLoop_2
                  ENDIF
!
!     kohn-sham exchange
!
                  vxc = xks*rhox
                  epsxc = 0.75D00*vxc
                  spag_nextblock_2 = 2
               CASE (2)
!
                  Dp(i) = vxc*Rhotot(i,Natom)
                  Dq(i) = epsxc*Rhotot(i,Natom)
                  EXIT SPAG_DispatchLoop_2
               END SELECT
            ENDDO SPAG_DispatchLoop_2
         ENDDO
!
!
         envxc = rsimp(Dp,Dr,rinf1,inf1)
         enexc = rsimp(Dq,Dr,rinf1,inf1)
!
!
         u = (u1+u2-envxc)/2.0D00
!
         eatom = t + u + enexc
!
         WRITE (6,99005) eatom , t , Tcore(Natom) , u , enexc , envxc
99005    FORMAT (15x,'<E>       ',e20.10/15x,'<t>       ',e20.10/15x,'<t-core>  ',e20.10/15x,'<u>       ',e20.10/15x,'<E-xc>    ', &
               & e20.10/15x,'<Vxc>     ',e20.10//)
!
         IF ( Iskip/=0 ) THEN
            WRITE (6,99006)
99006       FORMAT (10x,'one electron energies '///)
            DO i = 1 , Norb
               WRITE (6,99007) Nqn(i) , Titre(i) , Den(i)
99007          FORMAT (15x,i2,a4,e20.8)
            ENDDO
         ENDIF
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE totec
!
!*==RSIMP.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
FUNCTION rsimp(F,R,Rn,Irn)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 dx , F , R , Rn , rsimp , s , Stval , Xdel
   INTEGER i , ieven , Irn , isw , nl , np , NRAD
!*** End of declarations inserted by SPAG
   PARAMETER (NRAD=400)
!
!     **************************************
!     *   radial integration via simpson   *
!     **************************************
!
   COMMON /xmesh / Stval , Xdel
!
!     deltax=0.05 for r in logar. scale
!
   DIMENSION F(NRAD) , R(NRAD)
!
   dx = Xdel
!
   isw = 0
   rsimp = 0.0D00
   IF ( Irn<=2 ) RETURN
   ieven = (Irn/2)*2
   IF ( ieven==Irn ) isw = 1
   np = Irn - isw
   s = F(1)*R(1) + F(np)*R(np)
   nl = np - 1
   DO i = 2 , nl , 2
      s = s + 4.0D00*F(i)*R(i)
   ENDDO
   nl = nl - 1
   IF ( nl>=3 ) THEN
      DO i = 3 , nl , 2
         s = s + 2.0D00*F(i)*R(i)
      ENDDO
   ENDIF
   s = s*dx/3.0D00
   IF ( isw==1 ) THEN
      rsimp = s + (F(Irn)*R(Irn)+F(Irn-1)*R(Irn-1))*0.5D00*dx
      RETURN
   ENDIF
   rsimp = s
   RETURN
END FUNCTION rsimp
!
!*==RSIMP3.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
FUNCTION rsimp3(F,R,Rn,Jrn)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 dx , F , R , rad , Rn , rsimp3 , s , Stval , val , val1 , x1 , Xdel
   INTEGER i , irn , isw , iturn , Jrn , nl , np , NRAD
!*** End of declarations inserted by SPAG
   PARAMETER (NRAD=400)
!
!     **************************************
!     *   radial integration via simpson   *
!     *   with interpolation               *
!     **************************************
!
   COMMON /xmesh / Stval , Xdel
   DIMENSION rad(7) , x1(7)
!
!     deltax for r in logar. scale
!
   DIMENSION F(NRAD) , R(NRAD)
!
   dx = Xdel
!
   DO iturn = 1 , 7
      irn = Jrn + iturn - 6
      x1(iturn) = R(irn)
!
      isw = 0
      rsimp3 = 0.0D00
      IF ( irn<=2 ) RETURN
      IF ( irn/2*2==irn ) isw = 1
      np = irn - isw
      s = F(1)*R(1) + F(np)*R(np)
      nl = np - 1
      DO i = 2 , nl , 2
         s = s + 4.0D00*F(i)*R(i)
      ENDDO
      nl = nl - 1
      IF ( nl>=3 ) THEN
         DO i = 3 , nl , 2
            s = s + 2.0D00*F(i)*R(i)
         ENDDO
      ENDIF
      s = s*dx/3.0D00
      IF ( isw==1 ) THEN
         rad(iturn) = s + (F(irn)*R(irn)+F(irn-1)*R(irn-1))*0.5D00*dx
      ELSE
         rad(iturn) = s
      ENDIF
   ENDDO
!
   CALL interp(x1,rad,7,Rn,val,val1,.FALSE.)
!
   rsimp3 = val
END FUNCTION rsimp3
!
!*==FITRAD.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE fitrad(Np)
!
!     ******************************************************************
!     *                                                                *
!     * fits apw. radius to logarithmic scale                          *
!     *                                                                *
!     ******************************************************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 delx , dps , P , R , Ra , Rc , Rws , Stval , t , Vc , x , Xdel , Za
   INTEGER i , k , kk , kmax , kmax1 , kmax3 , l , n , Nk , NOAT , Np , NRAD , Nrc , Nrws , Nt , Nz
!*** End of declarations inserted by SPAG
!
   PARAMETER (NOAT=2)
   PARAMETER (NRAD=400)
!
   COMMON /dd    / Nt(NOAT) , Nz(NOAT) , Nk(NOAT) , Nrc(NOAT) , Rc(NOAT) , Vc(NOAT)
   COMMON /ee    / Za(NOAT,NRAD)
   COMMON /eee   / Ra(NOAT,NRAD)
   COMMON /oo    / P(NRAD) , R(NRAD)
   COMMON /xmesh / Stval , Xdel
   COMMON /ws    / Rws , Nrws(NOAT)
   INTEGER :: spag_nextblock_1
!
   delx = Xdel
!
   DO n = 1 , NOAT
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
            kmax = Nk(n)
            kmax1 = kmax + 1
            kmax3 = kmax1 - 7
!
            R(1) = 0.
            P(1) = -2.0D00*Nz(n)
!
            DO l = 1 , kmax
               P(l+1) = Za(n,l)*Ra(n,l)
               R(l+1) = Ra(n,l)
            ENDDO
!
            i = 0
!
            x = Stval - delx
            spag_nextblock_1 = 2
         CASE (2)
            x = x + delx
            i = i + 1
            t = dexp(x)
            IF ( i>Np ) EXIT SPAG_DispatchLoop_1
!
            DO k = 1 , kmax1
               IF ( R(k)>t ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
            kk = kmax3
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         CASE (3)
            IF ( k<=4 ) THEN
               kk = 1
            ELSEIF ( k>=kmax3 ) THEN
               kk = kmax3
            ELSE
               kk = k - 3
            ENDIF
            spag_nextblock_1 = 4
         CASE (4)
!
            CALL interp(R(kk),P(kk),7,t,Za(n,i),dps,.FALSE.)
!
!     if(t.gt.rws) za(n,i)=0.0d 00
            Za(n,i) = Za(n,i)/t
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
   ENDDO
!
END SUBROUTINE fitrad
!
!*==MIKE.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE mike(Ilat)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A0 , Avol , Conc , Confac , Ec , pi , Rws , Rws1 , third , tofp
   INTEGER Ilat , Lattyp , NOAT , Nrws , Nrws1
!*** End of declarations inserted by SPAG
   PARAMETER (NOAT=2)
   COMMON /ws    / Rws , Nrws(NOAT)
   COMMON /ws1   / Rws1 , Nrws1(NOAT)
   COMMON /kitty / Avol , Ec , Lattyp
   COMMON /lat   / A0 , Conc
   COMMON /samfac/ Confac
!
!     ec is the ewald constant as listed in:
!         j.f.janak, phys.rev.b9, 3985 (1974)
!     avol volume of unit cell
!
   third = 1.0D00/3.0D00
   pi = datan(1.0D00)*4.0D00
   tofp = 3.0D00/(4.0D00*pi)
!
   Lattyp = Ilat
   IF ( Lattyp==1 ) THEN
!
!     lattyp=1, fcc
!
      Avol = (A0**3.0D00)/4.0D00
      Ec = 4.8320664D00
      Rws1 = A0/dsqrt(8.0D00)
      WRITE (6,99002)
99002 FORMAT (10x,'fcc lattice'/)
   ELSEIF ( Lattyp==2 ) THEN
!
!     lattyp=2, bcc
!
      Avol = (A0**3.0D00)/4.0D00
      Ec = 4.085521D00
      Rws1 = A0*dsqrt(3.0D00/16.0D00)
      WRITE (6,99003)
99003 FORMAT (10x,'bcc lattice'/)
   ELSE
!
!     lattyp=0, simple cubic
!
      Avol = A0**3.0D00
      Ec = 3.1166857D00
      Rws1 = A0/2.0D00
      WRITE (6,99001)
99001 FORMAT (10x,'simple cubic lattice'/)
   ENDIF
!
!     set wigner-seitz-radius and conversion factor
!     for dimensionless units [d.u.]
!
   Rws = (tofp*Avol)**third
   Confac = (2.0D00*pi/A0)**2.0D00
!
   RETURN
END SUBROUTINE mike
!
!*==INTERP.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
SUBROUTINE interp(R,P,N,Rs,Ps,Dps,Deriv)
!
!     ****************************
!     *                          *
!     * interpolate via lagrange *
!     *                          *
!     ****************************
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 denom , Dps , dterm , dterm1 , P , Ps , R , Rs , term
   INTEGER i , j , k , N
!*** End of declarations inserted by SPAG
   LOGICAL Deriv , nodriv
   DIMENSION R(N) , P(N)
   nodriv = .NOT.Deriv
   Ps = 0.D00
   Dps = 0.D00
   DO j = 1 , N
      term = 1.D00
      denom = 1.D00
      dterm = 0.D00
      DO i = 1 , N
         IF ( i/=j ) THEN
            denom = denom*(R(j)-R(i))
            term = term*(Rs-R(i))
            IF ( .NOT.(nodriv) ) THEN
               dterm1 = 1.D00
               DO k = 1 , N
                  IF ( k/=j .AND. k/=i ) dterm1 = dterm1*(Rs-R(k))
               ENDDO
               dterm = dterm + dterm1
            ENDIF
         ENDIF
      ENDDO
      IF ( .NOT.(nodriv) ) Dps = Dps + dterm*P(j)/denom
      Ps = Ps + term*P(j)/denom
   ENDDO
END SUBROUTINE interp
