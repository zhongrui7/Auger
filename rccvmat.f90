!*==RCCVMAT.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM rccvmat
!
!
! ********************************************************************
!   RELATIVISTIC MATRIXELEMENTS FOR CORE-CORE-VALENCE AUGER SPECTRA
! ********************************************************************
!
! the highest possible values for momentum indices, l are set in l.par
! for core-core-valence lmax=2*lcore+lvmax+2
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a0 , cdff , cdfg , cdgf , cdgg , ceff , cefg , cegf , cegg , chuge , clight , coef , contr , csmall , de , dec , dummy , &
        & dx , e , e0
   REAL*8 e2 , earr , ec1 , ec3 , ev , fc1 , fc3 , fcwf , fvwf , gc1 , gc3 , gcwf , gvwf , half , pi , result , result1 , result2 ,&
        & rmt , rs
   REAL*8 rws , sd , sde , se , sigd , sigde , sige , signv , small , ssum , sum1 , sum2 , t13 , t1v , t32 , tiny , tm1m3 , tm1mv ,&
        & tm3m2 , tmvm2
   REAL*8 triad , tv2 , v , w6j , x , x0 , xj1 , xj3 , xs , y , y1 , y2 , y3 , zd , zde , ze , zeromt , zsum
   INTEGER i , iwigtes , j , j1 , j3 , jj , jj1 , jj2 , jj3 , kap , kap1 , kap2 , kap3 , l , l1 , l2 , l3 , lam , lamp , lkap
   INTEGER llam , llamp , lval , ne , nei , NEMAX , nmt , nonrel , norm , NRAD , nwf , nws
!*** End of declarations inserted by SPAG
!
   PARAMETER (NRAD=250)
   PARAMETER (NEMAX=50)
   INCLUDE 'l.par'
!
   DIMENSION v(NRAD) , rs(NRAD) , xs(NRAD) , earr(NEMAX)
!
   DIMENSION gc1(NRAD) , fc1(NRAD) , gc3(NRAD) , fc3(NRAD)
   DIMENSION gvwf(NRAD,-LMAX-1:LMAX) , fvwf(NRAD,-LMAX-1:LMAX)
   DIMENSION gcwf(NRAD,-LMAX-1:LMAX) , fcwf(NRAD,-LMAX-1:LMAX)
!
   DIMENSION sd(-LVMAX-1:LVMAX,NEMAX) , se(-LVMAX-1:LVMAX,NEMAX)
   DIMENSION sde(-LVMAX-1:LVMAX,NEMAX) , ssum(-LVMAX-1:LVMAX,NEMAX)
!
   DIMENSION zd(-LVMAX-1:LVMAX,NEMAX) , ze(-LVMAX-1:LVMAX,NEMAX)
   DIMENSION zde(-LVMAX-1:LVMAX,NEMAX) , zsum(-LVMAX-1:LVMAX,NEMAX)
!
   REAL*8 idgg(-LMAX-1:LMAX,0:LMAX)
   REAL*8 idgf(-LMAX-1:LMAX,0:LMAX)
   REAL*8 idfg(-LMAX-1:LMAX,0:LMAX)
   REAL*8 idff(-LMAX-1:LMAX,0:LMAX)
   REAL*8 iegg(-LMAX-1:LMAX,0:LMAX)
   REAL*8 iegf(-LMAX-1:LMAX,0:LMAX)
   REAL*8 iefg(-LMAX-1:LMAX,0:LMAX)
   REAL*8 ieff(-LMAX-1:LMAX,0:LMAX)
!
   DIMENSION cdgg(-LMAX-1:LMAX,0:LMAX) , cdgf(-LMAX-1:LMAX,0:LMAX) , cdfg(-LMAX-1:LMAX,0:LMAX) , cdff(-LMAX-1:LMAX,0:LMAX) ,       &
           & cegg(-LMAX-1:LMAX,0:LMAX) , cegf(-LMAX-1:LMAX,0:LMAX) , cefg(-LMAX-1:LMAX,0:LMAX) , ceff(-LMAX-1:LMAX,0:LMAX)
!
   INTEGER rl(-LMAX-1:LMAX) , rlb(-LMAX-1:LMAX) , rj(-LMAX-1:LMAX)
   DIMENSION triad(-LMAX-1:LMAX,-LMAX-1:LMAX,0:LMAX2)
   DIMENSION w6j(0:LVMAX,0:LMAX,0:LMAX,0:LMAX)
!
   INTEGER in , wf , po , mat , pun , zatom
!
   CHARACTER*1 name(5) , norb1(5) , norb3(5) , norb(5)
   CHARACTER*30 file
   INTEGER :: spag_nextblock_1
!
   DATA tiny , small , half/1.D-10 , 0.01D0 , 0.5D0/
   DATA pi/3.1415926535897932D0/
   DATA ev/13.606/
   DATA csmall/274.0746/ , chuge/1.D+08/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
! ********************************************************************
!
         in = 1
         wf = 2
         po = 3
         mat = 7
         pun = 8
!
         OPEN (in,FILE='rccvmat.in')
         OPEN (9,FILE='rccvmat.err')
!
         READ (in,99004) file
         OPEN (pun,FILE=file)
!
!
!  fill up fields for the relativistic quantumnumbers
!
!
         DO kap = -LMAX - 1 , -1
            l = -kap - 1
            rl(kap) = l
            rlb(kap) = l + 1
            rj(kap) = 2*l + 1
         ENDDO
         DO kap = 1 , LMAX
            l = kap
            rl(kap) = l
            rlb(kap) = l - 1
            rj(kap) = 2*l - 1
         ENDDO
!
!  control parameter for wigner coefficients
!
         READ (in,*) iwigtes
!
         IF ( iwigtes==1 ) THEN
            READ (in,*) lval , kap1 , kap3
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!
! read in maximum angular momentum quantumnumber for valence states
!
         READ (in,*) lval
         WRITE (pun,99001) lval
99001    FORMAT (1x,'LVAL=',t20,i4/)
!
!
! read in starting energy, increment for the energy panel, number of
! energy points and set up the energy panel
! (with respect to the muffin tin zero)
! the energy is measured in ryd
!
         READ (in,*) e0 , de , ne
         DO nei = 1 , ne
            earr(nei) = e0 + (nei-1)*de
         ENDDO
         WRITE (pun,99002) (earr(nei),nei=1,ne)
99002    FORMAT (1x,'ENERGY ARRAY:'/(10F7.2)/)
!
!
! check facility for non-relativistic limit
!
         READ (in,*) nonrel
         IF ( nonrel==1 ) THEN
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
         WRITE (pun,99003) norm
99003    FORMAT (1x,'NORMALIZATION OF VALENCE WF. (0 FOR MT, 1 FOR WS):',i4/)
!
!
! read in identifiers of core states
!
         READ (in,99005) (norb1(i),i=1,5)
         READ (in,99005) (norb3(i),i=1,5)
!
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
         READ (in,99004) file
         OPEN (po,FILE=file)
!
         READ (po,99005) (name(i),i=1,5)
         READ (po,*) zatom , x0 , nmt , dx , rws , a0 , zeromt
         READ (po,*) (v(j),j=1,nmt)
!
         coef = 2.*zatom/clight
         coef = coef*coef
!
!
! set up the logarithmic and radial mesh
! transform the potential data to v(R)
!
         x = x0
         xs(1) = x0
         rs(1) = dexp(x)
         signv = 1.
         IF ( v(1)>0. ) signv = -1.
         v(1) = signv*v(1)/rs(1) - zeromt
!
         DO j = 2 , nmt
            x = x + dx
            xs(j) = x
            rs(j) = dexp(x)
            v(j) = signv*v(j)/rs(j) - zeromt
         ENDDO
         rmt = rs(nmt)
!
         nws = 0
         DO j = nmt + 1 , NRAD
            x = x + dx
            xs(j) = x
            rs(j) = dexp(x)
            IF ( nws==0 .AND. rws<rs(j) ) nws = j
            v(j) = 0.D0
         ENDDO
!
         WRITE (pun,99006) zatom , x0 , dx , nmt , rmt , nws , rws , zeromt
99006    FORMAT (1x,' ZATOM =',t20,i4/1x,'    X0 =',t20,f5.1/1x,'    DX =',t20,f10.7/1x,'   NMT =',t20,i4/1x,'   RMT =',t20,       &
               & f10.7/1x,'   NWS =',t20,i4/1x,'   RWS =',t20,f10.7/1x,'ZEROMT =',t20,f5.1//1x,'RADIAL MESH AND POTENTIAL*R'/)
         WRITE (pun,99007) (rs(j),v(j)*rs(j),j=1,nws)
99007    FORMAT (4D20.10)
!
! read in the core wavefunctions as r*WF(r) and core energies (hartree)
! (the same radial mesh is supposed to be used as for the potential)
!
         READ (in,99004) file
         OPEN (wf,FILE=file)
         SPAG_Loop_2_1: DO
!
            READ (wf,99005) (norb(i),i=1,5)
            READ (wf,*) ec1
            READ (wf,*) kap1
            READ (wf,*) nwf
            DO i = 1 , nwf
               IF ( nonrel==1 ) THEN
                  READ (wf,*) dummy , gc1(i)
                  fc1(i) = 0.D0
               ELSE
                  READ (wf,*) dummy , gc1(i) , fc1(i)
               ENDIF
            ENDDO
            DO l = 1 , 5
               IF ( norb(l)/=norb1(l) ) CYCLE SPAG_Loop_2_1
            ENDDO
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDDO SPAG_Loop_2_1
         spag_nextblock_1 = 2
      CASE (2)
         SPAG_Loop_2_2: DO
            READ (wf,99005) (norb(i),i=1,3)
            READ (wf,*) ec3
            READ (wf,*) kap3
            READ (wf,*) nwf
            DO i = 1 , nwf
               IF ( nonrel==1 ) THEN
                  READ (wf,*) dummy , gc3(i)
                  fc3(i) = 0.D0
               ELSE
                  READ (wf,*) dummy , gc3(i) , fc3(i)
               ENDIF
            ENDDO
            DO l = 1 , 5
               IF ( norb(l)/=norb3(l) ) CYCLE SPAG_Loop_2_2
            ENDDO
!
            WRITE (pun,99008) ec1 , kap1
99008       FORMAT (//1x,' EC1 =',t20,f12.6/1x,'KAP1 =',t20,i4//1x,'CORE WAVEFUNCTION NO. 1'/)
            WRITE (pun,99009) (rs(j),gc1(j),fc1(j),j=1,nws)
            WRITE (pun,99010) ec3 , kap3
99010       FORMAT (//1x,' EC3 =',t20,f12.6/1x,'KAP3 =',t20,i4//1x,'CORE WAVEFUNCTION NO. 3'/)
            WRITE (pun,99009) (rs(j),gc3(j),fc3(j),j=1,nws)
!
            READ (in,99004) file
            OPEN (mat,FILE=file)
!
!
            WRITE (pun,99011)
99011       FORMAT (//' ****************** END INPUT *******************'//)
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ENDDO SPAG_Loop_2_2
         spag_nextblock_1 = 3
      CASE (3)
!
!
!  initialize the angular integration coefficients
!
!
         l1 = rl(kap1)
         l3 = rl(kap3)
         j1 = rj(kap1)
         j3 = rj(kap3)
         jj1 = j1 + 1
         jj3 = j3 + 1
         xj1 = j1*half
         xj3 = j3*half
         y1 = dsqrt(kap1*kap1-coef)
         y3 = dsqrt(kap3*kap3-coef)
         IF ( nonrel==1 ) y1 = dfloat(l1+1)
         IF ( nonrel==1 ) y3 = dfloat(l3+1)
!
         IF ( iwigtes==1 ) WRITE (pun,99012)
99012    FORMAT (/'TEST FOR TRIADS'//)
!
         CALL triadini(triad,rl,rlb,rj,pun,iwigtes)
!
         IF ( iwigtes==1 ) WRITE (pun,99013)
99013    FORMAT (//'TEST FOR WIGNER 6J SYMBOLS'//)
!
         CALL wigini(w6j,xj1,xj3,pun,iwigtes)
!
         IF ( iwigtes==1 ) STOP ' w3j, w6j test o.k.'
!
!
!
!         ***********************
!         * loop for the energy *
!         ***********************
!
!
         dec = 2.0*(ec3-ec1)
!
         DO nei = 1 , ne
            WRITE (6,*) nei
!
! compute valence wavefunctions
!
            e = earr(nei)
            WRITE (pun,*) ' valence wavefunctions'
            WRITE (pun,*) ' energy=' , e
            CALL wafu(e,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gvwf,fvwf,norm,1,nonrel,pun,lval,coef)
!
! compute continuum wavefunctions
!
            e2 = e + dec
            WRITE (pun,*) ' continuum wavefunctions'
            WRITE (pun,*) ' energy=' , e2
            CALL wafu(e2,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gcwf,fcwf,norm,2,nonrel,pun,lval,coef)
!
!
!         ******************
!         * loop for kappa *
!         ******************
!
            DO kap = -lval - 1 , lval
               IF ( kap/=0 ) THEN
                  y = dsqrt(kap*kap-coef)
                  IF ( nonrel==1 ) y = dfloat(rl(kap)+1)
                  jj = rj(kap) + 1
                  l = idint(jj*0.5+small) - 1
                  lkap = rl(kap)
!
                  WRITE (pun,*) ' kappa=' , kap
!
!
!  compute coulomb and exchange integrals for non-vanishing angular
!  integration coefficients
!
!
                  DO lam = 0 , LMAX
                     DO kap2 = -LMAX - 1 , LMAX
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
!
                  DO lam = 0 , LMAX
!
                     t13 = triad(kap1,kap3,lam)
                     tm1m3 = triad(-kap1,-kap3,lam)
                     t1v = triad(kap1,kap,lam)
                     tm1mv = triad(-kap1,-kap,lam)
!
                     DO kap2 = -LMAX - 1 , LMAX
                        IF ( kap2/=0 ) THEN
                           y2 = dsqrt(kap2*kap2-coef)
                           IF ( nonrel==1 ) y2 = dfloat(rl(kap2)+1)
!
                           tv2 = triad(kap,kap2,lam)
                           tmvm2 = triad(-kap,-kap2,lam)
                           t32 = triad(kap3,kap2,lam)
                           tm3m2 = triad(-kap3,-kap2,lam)
!
                           IF ( dabs(t13*tv2)>tiny ) THEN
                              CALL cei(gc1,gc3,gvwf(1,kap),gcwf(1,kap2),rs,nws,lam,y1,y3,y,y2,dx,rws,result)
                              idgg(kap2,lam) = result
                              cdgg(kap2,lam) = t13*tv2
                           ENDIF
!
                           IF ( dabs(t13*tmvm2)>tiny ) THEN
                              CALL cei(gc1,gc3,fvwf(1,kap),fcwf(1,kap2),rs,nws,lam,y1,y3,y,y2,dx,rws,result)
                              idgf(kap2,lam) = result
                              cdgf(kap2,lam) = t13*tmvm2
                           ENDIF
!
                           IF ( dabs(tm1m3*tv2)>tiny ) THEN
                              CALL cei(fc1,fc3,gvwf(1,kap),gcwf(1,kap2),rs,nws,lam,y1,y3,y,y2,dx,rws,result)
                              idfg(kap2,lam) = result
                              cdfg(kap2,lam) = tm1m3*tv2
                           ENDIF
!
                           IF ( dabs(tm1m3*tmvm2)>tiny ) THEN
                              CALL cei(fc1,fc3,fvwf(1,kap),fcwf(1,kap2),rs,nws,lam,y1,y3,y,y2,dx,rws,result)
                              idff(kap2,lam) = result
                              cdff(kap2,lam) = tm1m3*tmvm2
                           ENDIF
!
                           IF ( dabs(t1v*t32)>tiny ) THEN
                              CALL cei(gc1,gvwf(1,kap),gc3,gcwf(1,kap2),rs,nws,lam,y1,y,y3,y2,dx,rws,result)
                              iegg(kap2,lam) = result
                              cegg(kap2,lam) = t1v*t32
                           ENDIF
!
                           IF ( dabs(t1v*tm3m2)>tiny ) THEN
                              CALL cei(gc1,gvwf(1,kap),fc3,fcwf(1,kap2),rs,nws,lam,y1,y,y3,y2,dx,rws,result)
                              iegf(kap2,lam) = result
                              cegf(kap2,lam) = t1v*tm3m2
                           ENDIF
!
                           IF ( dabs(tm1mv*t32)>tiny ) THEN
                              CALL cei(fc1,fvwf(1,kap),gc3,gcwf(1,kap2),rs,nws,lam,y1,y,y3,y2,dx,rws,result)
                              iefg(kap2,lam) = result
                              cefg(kap2,lam) = tm1mv*t32
                           ENDIF
!
                           IF ( dabs(tm1mv*tm3m2)>tiny ) THEN
                              CALL cei(fc1,fvwf(1,kap),fc3,fcwf(1,kap2),rs,nws,lam,y1,y,y3,y2,dx,rws,result)
                              ieff(kap2,lam) = result
                              ceff(kap2,lam) = tm1mv*tm3m2
                           ENDIF
                        ENDIF
!
                     ENDDO
                  ENDDO
!
!
! ************* sum up for the direct and exchange terms
!
!
                  sigd = 0.
                  sige = 0.
!
                  DO lam = 0 , LMAX
                     llam = 2*lam + 1
!
                     DO kap2 = -LMAX - 1 , LMAX
                        IF ( kap2/=0 ) THEN
                           jj2 = rj(kap2) + 1
!
                           sum1 = idgg(kap2,lam)*cdgg(kap2,lam) + idgf(kap2,lam)*cdgf(kap2,lam) + idfg(kap2,lam)*cdfg(kap2,lam)    &
                                & + idff(kap2,lam)*cdff(kap2,lam)
                           sum2 = iegg(kap2,lam)*cegg(kap2,lam) + iegf(kap2,lam)*cegf(kap2,lam) + iefg(kap2,lam)*cefg(kap2,lam)    &
                                & + ieff(kap2,lam)*ceff(kap2,lam)
!
                           sigd = sigd + jj1*jj2*jj3*sum1*sum1/llam
                           sige = sige + jj1*jj2*jj3*sum2*sum2/llam
!
                           IF ( nonrel==1 .AND. nei==1 ) THEN
                              IF ( dabs(sum1)>tiny ) THEN
                                 WRITE (pun,*) ' nonzero for direct term'
                                 WRITE (pun,*) lam , kap2 , idgg(kap2,lam)
                              ENDIF
                              IF ( dabs(sum2)>tiny ) THEN
                                 WRITE (pun,*) ' nonzero for exchange term'
                                 WRITE (pun,*) lam , kap2 , iegg(kap2,lam)
                              ENDIF
                           ENDIF
                        ENDIF
!
                     ENDDO
                  ENDDO
!
                  sd(kap,nei) = sigd
                  se(kap,nei) = sige
!
!
! ************* sum up for the cross term
!
!
                  sigde = 0.
!
                  DO lam = 0 , LMAX
                     llam = 2*lam + 1
!
                     DO lamp = 0 , LMAX
                        llamp = 2*lam + 1
!
                        DO kap2 = -LMAX - 1 , LMAX
                           IF ( kap2/=0 ) THEN
                              jj2 = rj(kap2) + 1
                              l2 = idint(jj2*0.5+small) - 1
!
                              sum1 = idgg(kap2,lam)*cdgg(kap2,lam) + idgf(kap2,lam)*cdgf(kap2,lam) + idfg(kap2,lam)*cdfg(kap2,lam) &
                                   & + idff(kap2,lam)*cdff(kap2,lam)
                              sum2 = iegg(kap2,lamp)*cegg(kap2,lamp) + iegf(kap2,lamp)*cegf(kap2,lamp) + iefg(kap2,lamp)           &
                                   & *cefg(kap2,lamp) + ieff(kap2,lamp)*ceff(kap2,lamp)
!
                              contr = jj1*jj2*jj3*sum1*sum2*w6j(l,l2,lam,lamp)
                              sigde = sigde + contr
!
                              IF ( nonrel==1 .AND. nei==1 ) THEN
                                 IF ( dabs(contr)>tiny ) THEN
                                    WRITE (pun,*) ' nonzero for cross term'
                                    WRITE (pun,*) lam , lamp , kap2
                                    WRITE (pun,*) idgg(kap2,lam) , iegg(kap2,lamp)
                                 ENDIF
                              ENDIF
                           ENDIF
!
                        ENDDO
                     ENDDO
                  ENDDO
!
                  sde(kap,nei) = sigde
!
                  ssum(kap,nei) = sigd + sige - 2.0*sigde
!
!
!  easy check for the non-relativistic limit if nonrel eq 1 and
!  if kap1=kap3=-1
!
                  IF ( nonrel==1 .AND. kap1==-1 .AND. kap3==-1 ) THEN
                     CALL cei(gc1,gc3,gvwf(1,kap),gcwf(1,kap),rs,nws,0,y1,y3,y,y,dx,rws,result1)
                     zd(kap,nei) = 2.D0*result1*result1
!
                     CALL cei(gc1,gvwf(1,kap),gc3,gcwf(1,kap),rs,nws,lkap,y1,y,y3,y,dx,rws,result2)
                     ze(kap,nei) = 2.D0*result2*result2/(2.*lkap+1.)**2
!
                     zde(kap,nei) = result1*result2/(2.*lkap+1.)
!
                     zsum(kap,nei) = zd(kap,nei) + ze(kap,nei) - 2.D0*zde(kap,nei)
                  ENDIF
               ENDIF
!
            ENDDO
!
!
!         *************************
!         * end of loop for kappa *
!         *************************
!
         ENDDO
!
!
!         ******************************
!         * end of loop for the energy *
!         ******************************
!
!
!
         IF ( nonrel==1 .AND. kap1==-1 .AND. kap3==-1 ) THEN
!
            WRITE (pun,99014)
99014       FORMAT (//1x,'CHECK NON-RELATIVISTIC LIMIT'/)
            WRITE (pun,99015)
            WRITE (pun,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
            DO i = 1 , ne
               WRITE (pun,99020) earr(i) , (zd(kap,i),kap=-lval-1,-1) , (zd(kap,i),kap=1,lval)
            ENDDO
!
!
            WRITE (pun,99016)
            WRITE (pun,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
            DO i = 1 , ne
               WRITE (pun,99020) earr(i) , (ze(kap,i),kap=-lval-1,-1) , (ze(kap,i),kap=1,lval)
            ENDDO
!
!
            WRITE (pun,99017)
            WRITE (pun,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
            DO i = 1 , ne
               WRITE (pun,99020) earr(i) , (zde(kap,i),kap=-lval-1,-1) , (zde(kap,i),kap=1,lval)
            ENDDO
            WRITE (pun,*)
!
         ENDIF
!
!
         WRITE (pun,99015)
         WRITE (pun,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
         DO i = 1 , ne
            WRITE (pun,99020) earr(i) , (sd(kap,i),kap=-lval-1,-1) , (sd(kap,i),kap=1,lval)
         ENDDO
!
!
         WRITE (pun,99016)
         WRITE (pun,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
         DO i = 1 , ne
            WRITE (pun,99020) earr(i) , (se(kap,i),kap=-lval-1,-1) , (se(kap,i),kap=1,lval)
         ENDDO
!
!
         WRITE (pun,99017)
         WRITE (pun,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
         DO i = 1 , ne
            WRITE (pun,99020) earr(i) , (sde(kap,i),kap=-lval-1,-1) , (sde(kap,i),kap=1,lval)
         ENDDO
!
!
         WRITE (mat,99018) (name(i),i=1,5) , (norb1(i),i=1,5) , (norb3(i),i=1,5) , dec , lval , a0
99018    FORMAT (1x,'RELATIVISTIC CORE-CORE-VALENCE AUGER MATRIXELEMENTS'/5A1/1x,'FINAL CORE STATE'/5A1/1x,                        &
               & 'INITIAL CORE STATE'/5A1/1x,'TRANSITION ENERGY'/f12.5/1x,'MAXIMAL L QUANTUMNUMBER FOR VALENCE STATES'/i1/1x,      &
                &'LATTICE CONSTANT'/f10.6/)
         WRITE (mat,99019) (kap,kap=-lval-1,-1) , (kap,kap=1,lval)
         DO i = 1 , ne
            WRITE (mat,99020) earr(i) , (ssum(kap,i),kap=-lval-1,-1) , (ssum(kap,i),kap=1,lval)
         ENDDO
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99004 FORMAT (a30)
99005 FORMAT (5A1)
99009 FORMAT (3D20.10)
99015 FORMAT (//1x,'DIRECT TERM'/)
99016 FORMAT (//1x,'EXCHANGE TERM'/)
99017 FORMAT (//1x,'CROSS TERM'/)
99019 FORMAT (3x,'e',2x,12(6x,i2,5x))
99020 FORMAT (f6.2,12D13.5)
!
!
END PROGRAM rccvmat
!*==TRIADINI.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE triadini(Triad,Rl,Rlb,Rj,Pun,Iwigtes)
! ====================================================
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 fact , half , result , sixcof , small , thrcof , tiny , Triad , xl1 , xl2 , xlma1 , xlma2 , xlmat , xlmi1 , xlmi2
   INTEGER i1 , i2 , ier , Iwigtes , kap1 , kap2 , l1 , l2 , lam , lama , lama1 , lama2 , lami , lami1 , lami2 , NDIM
!*** End of declarations inserted by SPAG
!
   PARAMETER (NDIM=100)
   INCLUDE 'l.par'
!
   INTEGER Rl(-LMAX-1:LMAX) , Rlb(-LMAX-1:LMAX) , Rj(-LMAX-1:LMAX)
   DIMENSION Triad(-LMAX-1:LMAX,-LMAX-1:LMAX,0:LMAX2)
   DIMENSION thrcof(NDIM) , sixcof(NDIM)
   REAL*8 nul , j1 , j2
   INTEGER Pun
!
   DATA nul , tiny , small , half/0.D0 , 1.D-10 , 0.01D0 , 0.5D0/
!
!
!
!                       < kappa1 lambda kappa2 >
!
! for kappa1 = -lmax-1,-lmax,...,lmax-1,lmax
! for kappa2 = -lmax-1,-lmax,...,lmax-1,lmax
!
! convention for storage :  triad(kappa1,kappa2,lambda)
!
!
   DO kap1 = -LMAX - 1 , LMAX
      DO kap2 = -LMAX - 1 , LMAX
         DO lam = 0 , LMAX2
            Triad(kap1,kap2,lam) = nul
         ENDDO
      ENDDO
   ENDDO
!
   DO kap1 = -LMAX - 1 , LMAX
      IF ( kap1/=0 ) THEN
         l1 = Rl(kap1)
         xl1 = dfloat(l1)
         j1 = Rj(kap1)*half
         DO kap2 = kap1 , LMAX
            IF ( kap2/=0 ) THEN
               l2 = Rl(kap2)
               xl2 = dfloat(l2)
               j2 = Rj(kap2)*half
!
               fact = dsqrt((2.*xl1+1.)*(2.*xl2+1.))
!
               CALL rec3jj(thrcof,xl1,xl2,nul,nul,xlmi1,xlma1,xlmat,NDIM,ier)
               IF ( ier>=0 ) THEN
!
                  CALL rec6j(sixcof,xl1,xl2,half,j2,j1,xlmi2,xlma2,xlmat,NDIM,ier)
                  IF ( ier>=0 ) THEN
!
                     lami1 = idint(xlmi1+small)
                     lami2 = idint(xlmi2+small)
                     lama1 = idint(xlma1+small)
                     lama2 = idint(xlma2+small)
                     lami = max(lami1,lami2)
                     lama = min(lama1,lama2)
!
                     DO lam = lami , lama
                        IF ( mod(lam-lami1,2)/=1 ) THEN
                           i1 = lam - lami1 + 1
                           i2 = lam - lami2 + 1
                           result = fact*thrcof(i1)*sixcof(i2)
                           Triad(kap1,kap2,lam) = result
                           Triad(kap2,kap1,lam) = result
!
                           IF ( Iwigtes==1 ) THEN
                              WRITE (Pun,99001) kap1 , l1 , j1 , kap2 , l2 , j2 , lam , result*result
!
99001                         FORMAT (2I4,f5.1,5x,2I4,f5.1,5x,i4,10x,d17.10)
                           ENDIF
                        ENDIF
!
                     ENDDO
!
                     IF ( Iwigtes==1 ) WRITE (Pun,*) '               *****************'
                  ENDIF
               ENDIF
            ENDIF
!
         ENDDO
!
         IF ( Iwigtes==1 ) WRITE (Pun,*) '   ******************************************** '
      ENDIF
!
   ENDDO
!
END SUBROUTINE triadini
!*==WIGINI.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE wigini(W6j,J1,J3,Pun,Iwigtes)
! ================================================
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 am , half , result , sixcof , small , tiny , W6j , xlam , xlamma , xlammi , xlmat
   INTEGER ier , Iwigtes , l , l2 , lam , lama , lama1 , lama2 , lami , lami1 , lami2 , lamma , lammi , lamp , n , NDIM
!*** End of declarations inserted by SPAG
!
   PARAMETER (NDIM=100)
   INCLUDE 'l.par'
!
   DIMENSION W6j(0:LVMAX,0:LMAX,0:LMAX,0:LMAX)
   DIMENSION sixcof(NDIM)
   INTEGER Pun
   REAL*8 nul , j , J1 , j2 , J3
!
   DATA nul , tiny , small , half/0.D0 , 1.D-10 , 0.01D0 , 0.5D0/
!
!
!  wigner 6j coefficients
!
!            j3  lambda'  j2         lambda'  j1  j
!                               =
!            j   lambda   j1         lambda   j2  j3
!
! multiplied by (-)**(lambda+lambda'+1)
!
!
! for given j1,j3
!
! convention for storage : w6j(l,l2,lambda,lambda') ,
!
!  with  j=(2*l+1)/2        l=0,1,...,lvmax
!        j2=(2*l2+1)/2      l2=0,1,...,lmax
!
!
   DO l = 0 , LVMAX
      DO l2 = 0 , LMAX
         DO lam = 0 , LMAX
            DO lamp = 0 , LMAX
               W6j(l,l2,lam,lamp) = nul
            ENDDO
         ENDDO
      ENDDO
   ENDDO
!
   lami1 = idint(dabs(J1-J3)+small)
   lama1 = idint(J1+J3+small)
!
   DO l = 0 , LVMAX
      j = (2.*l+1.)*half
      DO l2 = 0 , LMAX
         j2 = (2.*l2+1.)*half
         lami2 = idint(dabs(j-j2)+small)
         lama2 = idint(j+j2+small)
         lami = max(lami1,lami2)
         lama = min(lama1,lama2)
         DO lam = 0 , LMAX
            IF ( lam>=lami .AND. lam<=lama ) THEN
               xlam = dfloat(lam)
!
               CALL rec6j(sixcof,J1,j,xlam,j2,J3,xlammi,xlamma,xlmat,NDIM,ier)
               IF ( ier>=0 ) THEN
!
                  lammi = idint(xlammi+small)
                  lamma = idint(xlamma+small)
!
                  DO lamp = lammi , lamma
                     n = lam + lamp + 1
                     am = dfloat(mod(n,2))
                     result = sixcof(lamp-lammi+1)
                     W6j(l,l2,lam,lamp) = (1.0-2.0*am)*result
!
                     IF ( Iwigtes==1 ) THEN
                        WRITE (Pun,99001) lamp , J1 , j , lam , j2 , J3 , result
 
!
99001                   FORMAT (i4,2F5.1,5x,i4,2F5.1,10x,d17.10)
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
!
         ENDDO
      ENDDO
   ENDDO
!
END SUBROUTINE wigini
!*==WAFU.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE wafu(En,Rl,Rlb,Rj,V,Rs,X0,Dx,Rmt,Nmt,Rws,Nws,P,Q,Norm,Nval,Nonrel,Pun,Lval,Coef)
!================================================
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a1 , b1 , chuge , clight , Coef , cose , csmall , Dx , ekappa , En , eta , fb , fb1 , fn , fn1 , P , Q , qint , r , ratfg
   REAL*8 ratx , Rmt , rr , Rs , Rws , sign , sine , sintg , sk , tane , tanx , V , vint , X0 , x1 , x2 , yk , yk2
   INTEGER i , kap , l , lact , lb , Lval , Nmt , Nonrel , Norm , NRAD , Nval , Nws
!*** End of declarations inserted by SPAG
!
   PARAMETER (NRAD=250)
   INCLUDE 'l.par'
!
   DIMENSION V(NRAD) , Rs(NRAD)
!
   DIMENSION P(NRAD,-LMAX-1:LMAX) , Q(NRAD,-LMAX-1:LMAX)
   DIMENSION rr(NRAD) , vint(-LMAX-1:LMAX) , eta(-LMAX-1:LMAX)
   DIMENSION tanx(-LMAX-1:LMAX) , ratx(-LMAX-1:LMAX)
   DIMENSION fb(0:LMAX+1) , fn(0:LMAX+1) , fb1(0:LMAX+1) , fn1(0:LMAX+1)
!
   INTEGER Rl(-LMAX-1:LMAX) , Rlb(-LMAX-1:LMAX) , Rj(-LMAX-1:LMAX)
   INTEGER Pun
!
   DATA csmall/274.0746/ , chuge/1.D+08/
!
!
   IF ( Nonrel==1 ) THEN
      clight = chuge
   ELSE
      clight = csmall
   ENDIF
!
   IF ( En>0 ) THEN
      sign = 1.
      ekappa = dsqrt(En)
   ELSE
      sign = -1.
      ekappa = dsqrt(-En)
   ENDIF
!
   CALL sbf1(En,Rmt,fb,fn)
!
   IF ( Nval==1 ) lact = Lval
   IF ( Nval==2 ) lact = LMAX
!
   DO kap = -lact - 1 , lact
      IF ( kap/=0 ) THEN
!
         l = Rl(kap)
         lb = Rlb(kap)
         IF ( kap>0 ) sk = ekappa
         IF ( kap<0 ) sk = -sign*ekappa
         yk = dsqrt(kap*kap-Coef)
!
         IF ( Nonrel==1 ) THEN
            yk = dfloat(l+1)
            lb = l + 1
            sk = -sign*ekappa
         ENDIF
!
         yk2 = yk + yk
!
         CALL comdir(En,kap,V,Nmt,Nws,Dx,X0,Q(1,kap),P(1,kap),ratfg,Nonrel)
!
         tane = (ratfg*fb(l)-sk*fb(lb))/(ratfg*fn(l)-sk*fn(lb))
         ratx(kap) = ratfg
         tanx(kap) = tane
!
! valence or continuum normalization
!
         IF ( Nval==1 ) THEN
            a1 = sk*ekappa*(fn(lb)-fb(lb)/tane)*Rmt/Q(Nmt,kap)/clight
            b1 = ekappa*(fn(l)-fb(l)/tane)*Rmt/P(Nmt,kap)
         ELSE
            eta(kap) = datan(tane)
            cose = dcos(eta(kap))
            sine = dsin(eta(kap))
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
            CALL sbf1(En,r,fb1,fn1)
            IF ( Nval==1 ) THEN
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
         IF ( Nval/=2 ) THEN
!
! normalize valence wavefuntions according to norm
!
            IF ( Norm==0 ) THEN
               vint(kap) = sintg(yk2,rr,Rs,Dx,Nmt)
            ELSE
               x1 = sintg(yk2,rr,Rs,Dx,Nws-1)
               x2 = sintg(yk2,rr,Rs,Dx,Nws)
               vint(kap) = x1 + (x2-x1)*(Rws-Rs(Nws-1))/(Rs(Nws)-Rs(Nws-1))
            ENDIF
            qint = dsqrt(1./vint(kap))
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
   IF ( Nval==1 ) THEN
!
      WRITE (Pun,99001)
      WRITE (Pun,99003) (kap,kap=-lact-1,-1) , (kap,kap=1,lact)
      WRITE (Pun,99004) (ratx(kap),kap=-lact-1,-1) , (ratx(kap),kap=1,lact)
      WRITE (Pun,99002)
99002 FORMAT (/1x,'TANGENT PHASESHIFTS'/)
      WRITE (Pun,99003) (kap,kap=-lact-1,-1) , (kap,kap=1,lact)
      WRITE (Pun,99004) (tanx(kap),kap=-lact-1,-1) , (tanx(kap),kap=1,lact)
!
   ELSE
!
      WRITE (Pun,99001)
      WRITE (Pun,99003) (kap,kap=-1,-lact-1,-1)
      WRITE (Pun,99004) (ratx(kap),kap=-1,-lact-1,-1)
      WRITE (Pun,99006) (kap,kap=1,lact)
      WRITE (Pun,99007) (ratx(kap),kap=1,lact)
      WRITE (Pun,99005)
99005 FORMAT (/1x,'PHASESHIFTS'/)
      WRITE (Pun,99003) (kap,kap=-1,-lact-1,-1)
      WRITE (Pun,99004) (eta(kap),kap=-1,-lact-1,-1)
      WRITE (Pun,99006) (kap,kap=1,lact)
      WRITE (Pun,99007) (eta(kap),kap=1,lact)
!
   ENDIF
99001 FORMAT (/1x,'CF/G RATIO'/)
99003 FORMAT (2x,12(5x,i3,5x))
99004 FORMAT (2x,12D13.5)
99006 FORMAT (15x,11(6x,i2,5x))
99007 FORMAT (15x,11D13.5)
!
END SUBROUTINE wafu
!*==COMDIR.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE comdir(E1,Kappa,Za,Nrc,Nnk,Dx,X0,Q,P,Ratfg,Nonrel)
!======================================================================
!
!   integration of relativistic radial dirac equations by milne
!   method and calculation of ratio cf/g
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 bgc , bgx , c , chuge , cin , csmall , Dx , e , E1 , hoc , P , pp , Q , qp , Ratfg , stval , sxk , sxm , t , tc
   REAL*8 test , u , uc , unp , unp2 , wc , wnp , wnp2 , x , X0 , xc , xk , z2 , Za
   INTEGER i , ik , jri , kap , Kappa , lkap , n , nit , Nnk , Nonrel , NRAD , Nrc
!*** End of declarations inserted by SPAG
!
   PARAMETER (NRAD=250)
!
   DIMENSION bgx(NRAD) , sxk(4) , sxm(4) , P(NRAD) , Q(NRAD) , pp(NRAD) , qp(NRAD) , Za(NRAD)
   INTEGER pun
   INTEGER :: spag_nextblock_1
   DATA test/1.E+05/ , pun/99/
   DATA csmall/274.0746/ , chuge/1.D+16/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         DO i = 1 , Nnk
            bgx(i) = Za(i)
         ENDDO
         kap = Kappa
         xk = dfloat(kap)
         jri = Nrc
         e = E1
         stval = X0
         z2 = -bgx(1)*exp(stval)
         tc = dexp(stval)
!
!
         IF ( Nonrel==1 ) THEN
            IF ( kap<0 ) lkap = -kap - 1
            IF ( kap>0 ) lkap = kap
            kap = -lkap - 1
            xk = dfloat(kap)
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
            IF ( dabs(hoc/xk)<=0.05 ) THEN
               u = (xk+dabs(xk))/hoc - 0.5*hoc/dabs(xk)
            ELSE
               u = (xk+dsqrt(xk*xk-hoc*hoc))/hoc
            ENDIF
            P(1) = 1.0E-20
            Q(1) = c*u*1.0E-20
!
         ENDIF
!
!
!     eq. 4.92
!
         IF ( Nonrel/=1 ) THEN
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
         spag_nextblock_1 = 2
      CASE (2)
         ik = 0
         xc = x
         bgc = bgx(n)
         wc = Q(n)
         uc = P(n)
         DO
            ik = ik + 1
            t = dexp(xc)
!
            IF ( Nonrel/=1 ) THEN
               sxk(ik) = Dx*(-xk*uc+t*wc*(cin*(e-bgc)+1.))
            ELSE
               sxk(ik) = Dx*(-xk*uc+t*wc)
            ENDIF
!
            sxm(ik) = Dx*(xk*wc-t*(e-bgc)*uc)
            IF ( ik==2 ) THEN
               uc = uc + 0.5*(sxk(2)-sxk(1))
               wc = wc + 0.5*(sxm(2)-sxm(1))
            ELSEIF ( ik==3 ) THEN
               xc = xc + 0.5*Dx
               uc = uc + sxk(3) - 0.5*sxk(2)
               wc = wc + sxm(3) - 0.5*sxm(2)
               bgc = bgx(n+1)
            ELSEIF ( ik==4 ) THEN
               Q(n+1) = Q(n) + (sxm(1)+2.0*sxm(2)+2.0*sxm(3)+sxm(4))/6.0
               P(n+1) = P(n) + (sxk(1)+2.0*sxk(2)+2.0*sxk(3)+sxk(4))/6.0
!
               IF ( Nonrel/=1 ) THEN
                  pp(n+1) = t*Q(n+1)*(cin*(e-bgc)+1.0) - xk*P(n+1)
               ELSE
                  pp(n+1) = t*Q(n+1) - xk*P(n+1)
               ENDIF
!
               qp(n+1) = xk*Q(n+1) - t*(e-bgc)*P(n+1)
               x = x + Dx
               n = n + 1
               IF ( n>=6 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ELSE
               xc = xc + 0.5*Dx
               uc = uc + 0.5*sxk(1)
               wc = wc + 0.5*sxm(1)
               bgc = 0.5*(bgc+bgx(n+1))
            ENDIF
         ENDDO
         spag_nextblock_1 = 3
      CASE (3)
!
!     milne method
!
         x = x + Dx
         t = dexp(x)
!
         unp = P(n-5) + 0.3*Dx*(11.*pp(n)-14.*pp(n-1)+26.*pp(n-2)-14.*pp(n-3)+11.0*pp(n-4))
!
         wnp = Q(n-5) + 0.3*Dx*(11.0*qp(n)-14.0*qp(n-1)+26.0*qp(n-2)-14.0*qp(n-3)+11.0*qp(n-4))
         nit = 0
         spag_nextblock_1 = 4
      CASE (4)
!
!
         IF ( Nonrel/=1 ) THEN
            pp(n+1) = t*(cin*(e-bgx(n+1))+1.0)*wnp - xk*unp
         ELSE
            pp(n+1) = t*wnp - xk*unp
         ENDIF
!
         qp(n+1) = xk*wnp - t*(e-bgx(n+1))*unp
!
         unp2 = P(n-3) + (7.0*pp(n+1)+32.0*pp(n)+12.0*pp(n-1)+32.0*pp(n-2)+7.0*pp(n-3))*2.0*Dx/45.0
!
         wnp2 = Q(n-3) + (7.0*qp(n+1)+32.0*qp(n)+12.0*qp(n-1)+32.0*qp(n-2)+7.0*qp(n-3))*2.0*Dx/45.0
!
         IF ( dabs(test*(unp2-unp))<=dabs(unp2) ) THEN
            IF ( dabs(test*(wnp2-wnp))<=dabs(wnp2) ) THEN
               spag_nextblock_1 = 5
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
!
         IF ( nit/=5 ) THEN
            nit = nit + 1
            wnp = wnp2
            unp = unp2
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
!
         Q(n+1) = wnp2
         P(n+1) = unp2
         n = n + 1
         IF ( n<Nnk ) THEN
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         Ratfg = Q(jri)/P(jri)
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
END SUBROUTINE comdir
!*==CEI.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cei(F2,F4,F1,F3,R,N,Lambda,L2,L4,L1,L3,Dx,Rws,Result)
! ====================================================================
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 arg1 , arg2 , arg3 , arg4 , Dx , F1 , F2 , F3 , F4 , R , Result , rl , rmlp1 , Rws , sintg
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
   CALL dintg(ll,1,arg1,arg2,R,Dx,m,Rws)
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
   CALL dintg(ll,-1,arg3,arg4,R,Dx,m,Rws)
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
   Result = sintg(ll,arg2,R,Dx,m)
!
END SUBROUTINE cei
!*==SINTG.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
DOUBLE PRECISION FUNCTION sintg(Ll,Fct,R,Dx,N)
! ========================================
!
! integration of fct over r using simpson integration
! r(i)=r0*exp(dx*(i-1))
!
   IMPLICIT NONE
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
   IF ( Ll+1.0<=tiny ) THEN
      WRITE (9,*) ' sintg: ll+1=' , Ll + 1
      sum = 0.
   ELSE
      rll = R(1)**Ll
      IF ( rll<=tiny ) THEN
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
   sintg = sum/3.D0
   IF ( mod(N,2)==0 ) sintg = sintg + (Fct(N)+Fct(N-1))/2.D0*(R(N)-R(N-1))
END FUNCTION sintg
!*==DINTG.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dintg(Ll,Idir,Fct,Yint,R,Dx,N,Rws)
! =================================================
!
! integration of fct(i) to yield yint(i)
! r(i)=r0*exp(dx*(i-1))
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 corr , Dx , fact , Fct , R , rll , Rws , sum , tiny , x0 , x1 , Yint
   INTEGER i , Idir , j , N
!*** End of declarations inserted by SPAG
!
   DIMENSION Fct(1) , Yint(1) , R(1)
   REAL*8 Ll
!
   DATA tiny/1.0D-20/
!
   IF ( Idir>0 ) THEN
!
      IF ( Ll+1.0<=tiny ) THEN
         WRITE (9,*) ' dintg: ll+1=' , Ll + 1
         sum = 0.
      ELSE
         rll = R(1)**Ll
         IF ( rll<=tiny ) THEN
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
         IF ( i==N-3 ) THEN
            CALL inter1(R(N-3),Yint(N-3),4,1,Rws,corr)
            sum = sum - corr
            DO j = N - 3 , N
               Yint(j) = Yint(j) - corr
            ENDDO
         ENDIF
         x1 = x0
      ENDDO
!
   ENDIF
END SUBROUTINE dintg
!*==INTER1.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE inter1(R,P,N,Id,Rs,Ps)
! =====================================
   IMPLICIT NONE
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
         IF ( i/=j ) THEN
            denom = denom*(R(j)-R(i))
            term = term*(Rs-R(i))
         ENDIF
      ENDDO
      Ps = Ps + term*P(j)/denom
   ENDDO
END SUBROUTINE inter1
!*==SBF.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE sbf(E,R,Fb,Fn)
! =============================
!
! spherical bessel and neuman functions for            e > 0
! modified spherical bessel and neuman functions for   e < 0
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 E , Fb , Fn , R , tl
   INTEGER l
!*** End of declarations inserted by SPAG
   INCLUDE 'l.par'
!
   COMPLEX*16 cfb(0:LMAX+1) , cfn(0:LMAX+1)
   COMPLEX*16 ecompl , eroot , x1 , x2 , sinx , cosx , cimmi
   DIMENSION Fb(0:LMAX+1) , Fn(0:LMAX+1)
!
   DATA cimmi/(0.0D0,-1.0D0)/
!
   ecompl = dcmplx(E,0.D0)
   eroot = cdsqrt(ecompl)
   x1 = eroot*R
   x2 = x1*x1
   sinx = cdsin(x1)
   cosx = cdcos(x1)
!
   cfb(0) = sinx/x1
   cfb(1) = sinx/x2 - cosx/x1
   cfn(0) = -cosx/x1
   cfn(1) = -cosx/x2 - sinx/x1
!
   DO l = 2 , LMAX + 1
      tl = dfloat(2*l-1)
      cfb(l) = tl*cfb(l-1)/x1 - cfb(l-2)
      cfn(l) = tl*cfn(l-1)/x1 - cfn(l-2)
   ENDDO
!
   IF ( E<0. ) THEN
      cfn(0) = cfn(0)*cimmi
      DO l = 1 , LMAX + 1
         cfb(l) = cfb(l)*cimmi**l
         cfn(l) = cfn(l)*cimmi**(l+1)
      ENDDO
   ENDIF
!
   DO l = 0 , LMAX + 1
      Fb(l) = dreal(cfb(l))
      Fn(l) = dreal(cfn(l))
   ENDDO
!
END SUBROUTINE sbf
!*==SBF1.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE sbf1(E,R,Fb,Fn)
! =============================
!
! spherical bessel and neuman functions for            e > 0
! modified spherical bessel and neuman functions for   e < 0
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 cosx , E , eroot , Fb , Fn , R , sinx , tl , x1 , x2
   INTEGER l
!*** End of declarations inserted by SPAG
   INCLUDE 'l.par'
!
   DIMENSION Fb(0:LMAX+1) , Fn(0:LMAX+1)
!
!
   IF ( E>0. ) THEN
!
      eroot = dsqrt(E)
      x1 = eroot*R
      x2 = x1*x1
      sinx = dsin(x1)
      cosx = dcos(x1)
!
      Fb(0) = sinx/x1
      Fb(1) = sinx/x2 - cosx/x1
      Fn(0) = -cosx/x1
      Fn(1) = -cosx/x2 - sinx/x1
!
      DO l = 2 , LMAX + 1
         tl = dfloat(2*l-1)
         Fb(l) = tl*Fb(l-1)/x1 - Fb(l-2)
         Fn(l) = tl*Fn(l-1)/x1 - Fn(l-2)
      ENDDO
!
   ELSE
!
      eroot = dsqrt(-E)
      x1 = eroot*R
      x2 = x1*x1
      sinx = dsinh(x1)
      cosx = dcosh(x1)
!
      Fb(0) = sinx/x1
      Fb(1) = -sinx/x2 + cosx/x1
      Fn(0) = cosx/x1
      Fn(1) = -cosx/x2 + sinx/x1
!
      DO l = 2 , LMAX + 1
         tl = dfloat(2*l-1)
         Fb(l) = -tl*Fb(l-1)/x1 + Fb(l-2)
         Fn(l) = -tl*Fn(l-1)/x1 + Fn(l-2)
      ENDDO
!
   ENDIF
!
END SUBROUTINE sbf1
!*==REC3JJ.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
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
END SUBROUTINE rec3jj
!*==REC6J.f90 processed by SPAG 8.02DA 11:10  3 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
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
