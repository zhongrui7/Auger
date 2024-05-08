program rcvvmat
!
! ********************************************************************
!  RELATIVISTIC MATRIXELEMENTS FOR CORE-VALENCE-VALENCE AUGER SPECTRA
! ********************************************************************
!
! the highest possible values for momentum, l are set in l.par
!                   (lmax=2*lvmax+lcore+2)
!
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NRAD = 250 , NEMAX = 15 , NEMAX2 = 2*nemax , LVMAX = 2 , LCORE = 2 , LMAX = lcore + 2*lvmax + 2 ,        &
                        & LMAX2 = 2*lmax
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: A0 , CLIGHT , COEF , CONTR , DE , DUMMY , DX , E , E0 , E2 , EC1 , EMIN2 , EP , RESULT , RMT , RWS , SIGD ,     &
                 & SIGDE , SIGNV , SUM1 , SUM2 , W1V , W1VP , WV2 , WVP2 , X , X0 , XJ1 , Y , Y1 , Y2 , YP , ZEROMT
   real(REAL64) , save :: CHUGE , CSMALL , EV , HALF , PI , SMALL , TINY
   real(REAL64) , dimension(-lvmax-1:lvmax,-lvmax-1:lvmax,nemax,nemax2) :: CRSECS
   real(REAL64) , dimension(nemax) :: EARR
   real(REAL64) , dimension(nemax2) :: EARR2
   real(REAL64) , dimension(nrad) :: FC1 , GC1 , RS , V , XS
   real(REAL64) , dimension(nrad,-lmax-1:lmax) :: FCWF , FVWF , FVWFP , GCWF , GVWF , GVWFP
   character(50) :: FILEN
   integer :: I , I2 , II , IN , IP , J , J1 , JJ , JJ1 , JJ2 , JJP , JP , KAP , KAP1 , KAP2 , KAPP , KOUT , KOUT1 , L , L1 , L2 , &
            & LAM , LAMP , LKAP , LKAPP , LLAM , LLAMP , LP , LVAL , MAT , NE , NE2 , NEI , NEI2 , NEIP , NMT , NONREL , NORM ,    &
            & NWF , NWS , PO , PUN , WF , ZATOM
   real(REAL64) , dimension(-lmax-1:lmax,0:lmax) :: IDFF , IDFG , IDGF , IDGG , IEFF , IEFG , IEGF , IEGG
   integer , dimension(nemax,nemax2) :: IRANGE
   character(1) , dimension(40) :: NAME
   character(1) , dimension(5) :: NORB , NORB1
   integer , dimension(-lmax-1:lmax) :: RJ , RL , RLB
   real(REAL64) , dimension(-lvmax-1:lvmax,-lvmax-1:lvmax) :: SD , SDE
   real(REAL64) , dimension(-lvmax-1:lvmax,-lvmax-1:lvmax,nemax) :: SSUM
   real(REAL64) , dimension(-lmax-1:lmax,-lmax-1:lmax,0:lmax2) :: W3J
   real(REAL64) , dimension(0:lvmax,0:lvmax,0:lmax,0:lmax,0:lmax) :: W6J
   external CEI , WAFU , WIG3J , WIG6J
!
! End of declarations rewritten by SPAG
!
!      include 'l.par'
 
   integer :: SPAG_NextBlock_1
!
   data tiny , small , half/1.D-10 , 0.01D0 , 0.5D0/
   data pi/3.1415926535897932D0/
   data ev/13.606/
   data csmall/274.0746/ , chuge/1.D+08/
   SPAG_NextBlock_1 = 1
   spag_dispatchloop_1: do
      select case (SPAG_NextBlock_1)
      case (1)
!
! ********************************************************************
!
         in = 1
         wf = 2
         po = 3
         mat = 7
         pun = 8
!
         open (in,file='rcvvmat.in',status='old')
         open (9,file='rcvvmat.err',status='unknown')
         SPAG_NextBlock_1 = 2
      case (2)
!
!
         read (in,99001,end=99999) filen
         open (pun,file=filen,status='unknown')
!
         read (in,99001) filen
         open (po,file=filen,status='unknown')
!
         read (in,99001) filen
         open (wf,file=filen,status='unknown')
!
         read (in,99001) filen
         open (mat,file=filen,status='unknown')
!
!  fill up fields for the relativistic quantumnumbers
!
         do kap = -lmax - 1 , -1
            l = -kap - 1
            rl(kap) = l
            rlb(kap) = l + 1
            rj(kap) = 2*l + 1
         enddo
         do kap = 1 , lmax
            l = kap
            rl(kap) = l
            rlb(kap) = l - 1
            rj(kap) = 2*l - 1
         enddo
!
! read in maximum angular momentum quantumnumber for valence states
!
         read (in,*) lval
         write (pun,99002) lval
99002    format (1x,'LVAL=',t20,i4/)
!
! read in starting energy, increment for the energy panel, number of
! energy points for valence band
! (with respect to the muffin tin zero)
! the energy is measured in ryd
!
         read (in,*) e0 , de , ne
!
! check facility for non-relativistic limit
!
         read (in,*) nonrel
         if ( nonrel==1 ) then
            clight = chuge
         else
            clight = csmall
         endif
!
!  normalization of the valence wavefunctions :
!                  norm eq 0  within the muffin tin
!                  norm eq 1  within the wigner-seitz sphere
!
         read (in,*) norm
         write (pun,99003) norm
99003    format (1x,'NORMALIZATION OF VALENCE WF. (0 FOR MT, 1 FOR WS):',i4/)
!
! read in identifier of core state
!
         read (in,99004) (norb1(i),i=1,5)
!
         read (in,*)
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
         read (po,99004) (name(i),i=1,40)
         read (po,*) zatom , x0 , nmt , dx , rws , a0 , zeromt
         read (po,*) (v(j),j=1,nmt)
!
         coef = 2.*zatom/clight
         coef = coef*coef
!
! set up the logarithmic and radial mesh
! transform the potential data to v(R)
!
         x = x0
         xs(1) = x0
         rs(1) = dexp(x)
         signv = 1.
         if ( v(1)>0. ) signv = -1.
         v(1) = signv*v(1)/rs(1) - zeromt
!
         do j = 2 , nmt
            x = x + dx
            xs(j) = x
            rs(j) = dexp(x)
            v(j) = signv*v(j)/rs(j) - zeromt
         enddo
         rmt = rs(nmt)
!
         nws = 0
         do j = nmt + 1 , nrad
            x = x + dx
            xs(j) = x
            rs(j) = dexp(x)
            if ( nws==0 .and. rws<rs(j) ) nws = j
            v(j) = 0.D0
         enddo
!
         write (*,99004) (name(i),i=1,40)
         write (pun,99004) (name(i),i=1,40)
         write (pun,*)
         write (pun,99005) zatom , x0 , dx , nmt , rmt , nws , rws , zeromt
99005    format (1x,' ZATOM =',t20,i4/1x,'    X0 =',t20,f5.1/1x,'    DX =',t20,f10.7/1x,'   NMT =',t20,i4/1x,'   RMT =',t20,       &
               & f10.7/1x,'   NWS =',t20,i4/1x,'   RWS =',t20,f10.7/1x,'ZEROMT =',t20,f5.1//1x,'RADIAL MESH AND POTENTIAL*R'/)
         write (pun,99006) (rs(j),v(j)*rs(j),j=1,nws)
99006    format (4D20.10)
         do
!
! read in the core wavefunctions as r*WF(r) and core energies (hartree)
! (the same radial mesh is supposed to be used as for the potential)
!
!
            read (wf,99004) (norb(i),i=1,5)
            read (wf,*) ec1
            ec1 = 2.*ec1
            read (wf,*) kap1
            read (wf,*) nwf
!
            if ( nonrel==1 ) then
               do i = 1 , nwf
                  read (wf,*) dummy , gc1(i)
                  fc1(i) = 0.D0
               enddo
            else
               do i = 1 , nwf
                  read (wf,*) dummy , gc1(i) , fc1(i)
               enddo
            endif
!
            if ( norb(1)==norb1(1) .and. norb(2)==norb1(2) .and. norb(3)==norb1(3) ) then
!
               write (pun,99007) (norb(i),i=1,5) , ec1 , kap1
99007          format (//1x,'CORE WAVEFUNCTION'/1x,'ORB. =',t20,5A1/1x,' EC1 =',t20,f12.6/1x,'KAP1 =',t20,i4//)
               write (pun,99008) (rs(j),gc1(j),fc1(j),j=1,nws)
99008          format (3D20.10)
!
               write (pun,99009)
99009          format (//' ****************** END INPUT *******************'//)
!
!  set up energy panel for valence band :   earr
!            and       for continuum    :   earr2
!
               ii = 0
               do i = 1 , ne
                  earr(i) = e0 + (i-1)*de
               enddo
               ne2 = 2*ne - 1
               emin2 = 2.*e0 - ec1
               write (pun,99010) (earr(i),i=1,ne)
99010          format (1x,'ENERGY PLANE'//(7x,15F7.2)/)
!
               do i = 1 , ne2
                  earr2(i) = emin2 + (i-1)*de
                  do j = 1 , ne
                     jp = i + 1 - j
                     if ( jp>=1 .and. jp<=ne ) then
                        ii = ii + 1
                        irange(j,i) = ii
                     else
                        irange(j,i) = 0
                     endif
                  enddo
                  write (pun,99011) earr2(i) , (irange(j,i),j=1,ne)
99011             format (f7.2,15(2x,i3,2x))
               enddo
!
               write (pun,*)
!
!  initialize the angular integration coefficients
!
               l1 = rl(kap1)
               j1 = rj(kap1)
               xj1 = j1*half
               jj1 = j1 + 1
               y1 = dsqrt(kap1*kap1-coef)
               if ( nonrel==1 ) y1 = dfloat(l1+1)
!
               call wig3j(w3j,pun)
               call wig6j(w6j,xj1,pun)
!
!
!         ************************************
!         * loop for the energy in continuum *
!         ************************************
!
               do nei2 = 1 , ne2
                  write (6,*) nei2
!
! compute continuum wavefunctions
!
                  e2 = earr2(nei2)
!
                  kout = 0
                  if ( nei2==ne ) kout = 1
!
                  write (pun,99012) e2
99012             format (' e2=',f15.5)
                  if ( kout==1 ) write (pun,*)
                  if ( kout==1 ) write (pun,*) ' continuum wavefunctions'
!
                  call wafu(e2,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gcwf,fcwf,norm,2,nonrel,pun,lval,coef,kout)
!
!
!         ***************************************
!         * loop for the energy in valence band *
!         ***************************************
!
                  do nei = 1 , ne
                     neip = nei2 + 1 - nei
                     if ( irange(nei,nei2)/=0 ) then
!
                        kout1 = 0
                        if ( nonrel==1 .and. irange(nei,nei2)==ne*ne ) kout1 = 1
!
! compute valence wavefunctions
!
                        e = earr(nei)
                        ep = earr(neip)
!
                        if ( kout==1 ) write (pun,*)
                        if ( kout==1 ) write (pun,*) ' valence wavefunctions'
                        write (pun,99013) e
99013                   format (t20,' e =',f15.5)
!
                        call wafu(e,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gvwf,fvwf,norm,1,nonrel,pun,lval,coef,kout)
!
                        call wafu(ep,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,gvwfp,fvwfp,norm,1,nonrel,pun,lval,coef,0)
!
!         ******************
!         * loop for kappa *
!         ******************
!
                        do kap = -lval - 1 , lval
                           if ( kap/=0 ) then
                              y = dsqrt(kap*kap-coef)
                              if ( nonrel==1 ) y = dfloat(rl(kap)+1)
                              jj = rj(kap) + 1
                              l = idint(jj*0.5+small) - 1
                              lkap = rl(kap)
!
                              if ( kout1==1 ) write (pun,*) ' kappa=' , kap
!
!         *******************
!         * loop for kappa' *
!         *******************
!
                              do kapp = -lval - 1 , lval
                                 if ( kapp/=0 ) then
                                    yp = dsqrt(kapp*kapp-coef)
                                    if ( nonrel==1 ) yp = dfloat(rl(kapp)+1)
                                    jjp = rj(kapp) + 1
                                    lp = idint(jjp*0.5+small) - 1
                                    lkapp = rl(kapp)
!
                                    if ( kout1==1 ) write (pun,*) ' kappa"=' , kapp
!
!  compute coulomb and exchange integrals for non-vanishing angular
!  integration coefficients
!
                                    do lam = 0 , lmax
                                       do kap2 = -lmax - 1 , lmax
                                         idgg(kap2,lam) = 0.
                                         idfg(kap2,lam) = 0.
                                         idgf(kap2,lam) = 0.
                                         idff(kap2,lam) = 0.
                                         iegg(kap2,lam) = 0.
                                         iefg(kap2,lam) = 0.
                                         iegf(kap2,lam) = 0.
                                         ieff(kap2,lam) = 0.
                                       enddo
                                    enddo
!
                                    do lam = 0 , lmax
!
                                       w1v = w3j(kap1,kap,lam)
                                       w1vp = w3j(kap1,kapp,lam)
!
                                       do kap2 = -lmax - 1 , lmax
                                         if ( kap2/=0 ) then
                                         y2 = dsqrt(kap2*kap2-coef)
                                         if ( nonrel==1 ) y2 = dfloat(rl(kap2)+1)
!
                                         wv2 = w3j(kap,kap2,lam)
                                         wvp2 = w3j(kapp,kap2,lam)
!
                                         if ( dabs(w1v*wvp2)>tiny ) then
                                         call cei(gc1,gvwf(1,kap),gvwfp(1,kapp),gcwf(1,kap2),rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                                         idgg(kap2,lam) = result
                                         call cei(gc1,gvwf(1,kap),fvwfp(1,kapp),fcwf(1,kap2),rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                                         idgf(kap2,lam) = result
                                         call cei(fc1,fvwf(1,kap),gvwfp(1,kapp),gcwf(1,kap2),rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                                         idfg(kap2,lam) = result
                                         call cei(fc1,fvwf(1,kap),fvwfp(1,kapp),fcwf(1,kap2),rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                                         idff(kap2,lam) = result
                                         endif
!
                                         if ( dabs(w1vp*wv2)>tiny ) then
                                         call cei(gc1,gvwfp(1,kapp),gvwf(1,kap),gcwf(1,kap2),rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                                         iegg(kap2,lam) = result
                                         call cei(gc1,gvwfp(1,kapp),fvwf(1,kap),fcwf(1,kap2),rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                                         iegf(kap2,lam) = result
                                         call cei(fc1,fvwfp(1,kapp),gvwf(1,kap),gcwf(1,kap2),rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                                         iefg(kap2,lam) = result
                                         call cei(fc1,fvwfp(1,kapp),fvwf(1,kap),fcwf(1,kap2),rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                                         ieff(kap2,lam) = result
                                         endif
                                         endif
!
                                       enddo
                                    enddo
!
!
! ************* sum up for the direct term
!
!
                                    sigd = 0.
!
                                    do lam = 0 , lmax
                                       llam = 2*lam + 1
                                       w1v = w3j(kap1,kap,lam)
                                       w1v = w1v*w1v
!
                                       do kap2 = -lmax - 1 , lmax
                                         if ( kap2/=0 ) then
                                         jj2 = rj(kap2) + 1
                                         wvp2 = w3j(kapp,kap2,lam)
                                         wvp2 = wvp2*wvp2
!
                                         sum1 = idgg(kap2,lam) + idgf(kap2,lam) + idfg(kap2,lam) + idff(kap2,lam)
!
                                         sigd = sigd + jj1*jj2*w1v*wvp2*sum1*sum1/llam
!
                                         if ( kout1==1 .and. dabs(sum1)>tiny ) then
                                         write (pun,*) ' nonzero for direct term'
                                         write (pun,*) lam , kap2 , sum1
                                         endif
                                         endif
!
                                       enddo
                                    enddo
!
                                    sd(kap,kapp) = sigd
!
!
! ************* sum up for the cross term
!
!
                                    sigde = 0.
!
                                    do lam = 0 , lmax
                                       llam = 2*lam + 1
                                       w1v = w3j(kap1,kap,lam)
!
                                       do lamp = 0 , lmax
                                         llamp = 2*lam + 1
                                         w1vp = w3j(kap1,kapp,lamp)
!
                                         do kap2 = -lmax - 1 , lmax
                                         if ( kap2/=0 ) then
                                         jj2 = rj(kap2) + 1
                                         l2 = idint(jj2*0.5+small) - 1
                                         wv2 = w3j(kap,kap2,lamp)
                                         wvp2 = w3j(kapp,kap2,lam)
!
                                         sum1 = idgg(kap2,lam) + idgf(kap2,lam) + idfg(kap2,lam) + idff(kap2,lam)
                                         sum2 = iegg(kap2,lamp) + iegf(kap2,lamp) + iefg(kap2,lamp) + ieff(kap2,lamp)
!
                                         contr = jj1*jj2*w1v*w1vp*wv2*wvp2*sum1*sum2*w6j(l,lp,l2,lam,lamp)
                                         sigde = sigde + contr
!
                                         if ( kout1==1 .and. dabs(contr)>tiny ) then
                                         write (pun,*) ' nonzero for cross term'
                                         write (pun,*) lam , lamp , kap2 , sum1*sum2
                                         endif
                                         endif
!
                                         enddo
                                       enddo
                                    enddo
!
                                    sde(kap,kapp) = sigde
!
                                    ssum(kap,kapp,nei) = 2.0*(sigd-sigde)
                                 endif
!
                              enddo
                           endif
                        enddo
!
!
!         ***************************
!         * end of loop for kappa's *
!         ***************************
!
                        if ( kout1==1 ) then
!
                           write (pun,*)
                           write (pun,*) ' DIRECT TERM'
                           do kap = -lval - 1 , lval
                              if ( kap/=0 ) then
                                 write (pun,99014) (kap,kapp,kapp=-lval-1,-1) , (kap,kapp,kapp=1,lval)
                                 write (pun,99015) (sd(kap,kapp),kapp=-lval-1,-1) , (sd(kap,kapp),kapp=1,lval)
                              endif
                           enddo
                           write (pun,*)
                           write (pun,*) ' CROSS TERM'
                           do kap = -lval - 1 , lval
                              if ( kap/=0 ) then
                                 write (pun,99014) (kap,kapp,kapp=-lval-1,-1) , (kap,kapp,kapp=1,lval)
                                 write (pun,99015) (sde(kap,kapp),kapp=-lval-1,-1) , (sde(kap,kapp),kapp=1,lval)
                              endif
                           enddo
!
                        endif
                     endif
!
                  enddo
!
!  ************************************
!  construct symmetrized cross sections
!  ************************************
!
                  do i = 1 , ne
                     ip = nei2 + 1 - i
                     if ( irange(i,nei2)/=0 ) then
                        do kap = -lval - 1 , lval
                           if ( kap/=0 ) then
                              crsecs(kap,kap,i,nei2) = (ssum(kap,kap,i)+ssum(kap,kap,ip))*0.5
                              do kapp = -lval - 1 , kap - 1
                                 if ( kapp/=0 ) crsecs(kap,kapp,i,nei2) = ssum(kap,kapp,i) + ssum(kapp,kap,ip)
                              enddo
                           endif
                        enddo
                     endif
                  enddo
!
!
               enddo
!
!         ********************************
!         * end of loop for the energies *
!         ********************************
!
               write (mat,99016) (name(i),i=1,40) , (norb1(i),i=1,5) , ec1 , e0 , de , ne , lval , a0
99016          format (1x,'RELATIVISTIC CORE-VALENCE-VALENCE AUGER                MATRIXELEMENTS'/40A1/1x,'CORE STATE'/5A1/1x,     &
                      &'ENERGY OF THE CORE STATE'/f12.5/1x,'E0, DE, NE FOR VALENCE BAND'/2F6.2,i5/1x,                              &
                      &'MAXIMAL L QUANTUMNUMBER FOR VALENCE STATES'/i1/1x,'LATTICE CONSTANT'/f10.6/)
!
               write (mat,99017) (((kap,kapp),kapp=-lval-1,kap),kap=-lval-1,-1) ,                                                  &
                               & (((kap,kapp),kapp=-lval-1,-1),((kap,kapp),kapp=1,kap),kap=1,lval)
99017          format (1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,28(i6,i3,'    '))
!
               do i2 = 1 , ne2
                  do i = 1 , ne
                     ip = i2 + 1 - i
                     if ( irange(i,i2)/=0 ) then
                        write (mat,99018) i2 , i , earr2(i2) , earr(ip) , earr(i) , ((crsecs(kap,kapp,i,i2),kapp=-lval-1,kap),     &
                             & kap=-lval-1,-1) , ((crsecs(kap,kapp,i,i2),kapp=-lval-1,-1),(crsecs(kap,kapp,i,i2),kapp=1,kap),kap=1,&
                             & lval)
99018                   format (2I3,3F6.2,28D13.5)
                     endif
                  enddo
               enddo
!
               close (pun)
               close (po)
               close (wf)
               close (mat)
               SPAG_NextBlock_1 = 2
               cycle spag_dispatchloop_1
            endif
         enddo
         exit spag_dispatchloop_1
      end select
   enddo spag_dispatchloop_1
!
!
99001 format (a50)
99004 format (40A1)
99014 format (7(4x,i2,1x,i2,4x))
99015 format (7D13.5)
!
99999 end program RCVVMAT
!
!
subroutine wafu(en,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,p,q,norm,nval,nonrel,pun,lval,coef,kout)
!================================================
!
   use iso_fortran_env
   implicit none
   include 'l.par'
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NRAD = 250
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) :: EN
   integer , intent(in) , dimension(-lmax-1:lmax) :: RL
   integer , intent(in) , dimension(-lmax-1:lmax) :: RLB
   integer , dimension(-lmax-1:lmax) :: RJ
   real(REAL64) , dimension(nrad) :: V
   real(REAL64) , dimension(nrad) :: RS
   real(REAL64) :: X0
   real(REAL64) :: DX
   real(REAL64) :: RMT
   integer :: NMT
   real(REAL64) , intent(in) :: RWS
   integer :: NWS
   real(REAL64) , intent(inout) , dimension(nrad,-lmax-1:lmax) :: P
   real(REAL64) , intent(inout) , dimension(nrad,-lmax-1:lmax) :: Q
   integer , intent(in) :: NORM
   integer , intent(in) :: NVAL
   integer :: NONREL
   integer , intent(in) :: PUN
   integer , intent(in) :: LVAL
   real(REAL64) , intent(in) :: COEF
   integer , intent(in) :: KOUT
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: A1 , B1 , CLIGHT , COSE , EKAPPA , QINT , R , RATFG , SIGN , SINE , SK , TANE , X1 , X2 , YK , YK2
   real(REAL64) , save :: CHUGE , CSMALL
   real(REAL64) , dimension(-lmax-1:lmax) :: ETA , RATX , TANX , VINT
   real(REAL64) , dimension(0:lmax+1) :: FB , FB1 , FN , FN1
   integer :: I , KAP , L , LACT , LB
   real(REAL64) , dimension(nrad) :: RR
   real(REAL64) , external :: SINTG
   external COMDIR , SBF1
!
! End of declarations rewritten by SPAG
!
   data csmall/274.0746/ , chuge/1.D+08/
!
!
   if ( nonrel==1 ) then
      clight = chuge
   else
      clight = csmall
   endif
!
   if ( en>0 ) then
      sign = 1.
      ekappa = dsqrt(en)
   else
      sign = -1.
      ekappa = dsqrt(-en)
   endif
!
   call sbf1(en,rmt,fb,fn)
!
   if ( nval==1 ) lact = lval
   if ( nval==2 ) lact = lmax
!
   do kap = -lact - 1 , lact
      if ( kap/=0 ) then
!
         l = rl(kap)
         lb = rlb(kap)
         if ( kap>0 ) sk = ekappa
         if ( kap<0 ) sk = -sign*ekappa
         yk = dsqrt(kap*kap-coef)
!
         if ( nonrel==1 ) then
            yk = dfloat(l+1)
            lb = l + 1
            sk = -sign*ekappa
         endif
!
         yk2 = yk + yk
!
         call comdir(en,kap,v,nmt,nws,dx,x0,q(1,kap),p(1,kap),ratfg,nonrel)
!
         tane = (ratfg*fb(l)-sk*fb(lb))/(ratfg*fn(l)-sk*fn(lb))
         ratx(kap) = ratfg
         tanx(kap) = tane
!
! valence or continuum normalization
!
         if ( nval==1 ) then
            a1 = sk*ekappa*(fn(lb)-fb(lb)/tane)*rmt/q(nmt,kap)/clight
            b1 = ekappa*(fn(l)-fb(l)/tane)*rmt/p(nmt,kap)
         else
            eta(kap) = datan(tane)
            cose = dcos(eta(kap))
            sine = dsin(eta(kap))
            a1 = sk*(cose*fb(lb)-sine*fn(lb))*rmt/q(nmt,kap)/clight
            b1 = (cose*fb(l)-sine*fn(l))*rmt/p(nmt,kap)
!     if(b1.lt.0.) then
!       if(eta(kap).gt.pi) then
!       eta(kap)=eta(kap)-pi
!       else
!       eta(kap)=eta(kap)+pi
!       endif
!       goto 310
!     endif
         endif
!
!
         do i = 1 , nmt
            p(i,kap) = p(i,kap)*b1
            q(i,kap) = q(i,kap)*a1
            rr(i) = p(i,kap)*p(i,kap) + q(i,kap)*q(i,kap)
         enddo
!
!
!  set up the wavefunctions beyond the muffin tin radius
!
         do i = nmt + 1 , nws + 1
            r = rs(i)
            call sbf1(en,r,fb1,fn1)
            if ( nval==1 ) then
               q(i,kap) = sk*ekappa*(fn1(lb)-fb1(lb)/tane)*r/clight
               p(i,kap) = ekappa*(fn1(l)-fb1(l)/tane)*r
            else
               q(i,kap) = sk*(cose*fb1(lb)-sine*fn1(lb))*r/clight
               p(i,kap) = (cose*fb1(l)-sine*fn1(l))*r
            endif
            rr(i) = p(i,kap)*p(i,kap) + q(i,kap)*q(i,kap)
         enddo
!
!
         if ( nval/=2 ) then
!
! normalize valence wavefuntions according to norm
!
            if ( norm==0 ) then
               vint(kap) = sintg(yk2,rr,rs,dx,nmt)
            else
               x1 = sintg(yk2,rr,rs,dx,nws-1)
               x2 = sintg(yk2,rr,rs,dx,nws)
               vint(kap) = x1 + (x2-x1)*(rws-rs(nws-1))/(rs(nws)-rs(nws-1))
            endif
            qint = dsqrt(1./vint(kap))
            do i = 1 , nrad
               p(i,kap) = p(i,kap)*qint
               q(i,kap) = q(i,kap)*qint
            enddo
         endif
      endif
!
!
   enddo
!
!
   if ( kout==1 ) then
!
      if ( nval==1 ) then
!
         write (pun,99001)
         write (pun,99003) (kap,kap=-lact-1,-1) , (kap,kap=1,lact)
         write (pun,99004) (ratx(kap),kap=-lact-1,-1) , (ratx(kap),kap=1,lact)
         write (pun,99002)
99002    format (/1x,'TANGENT PHASESHIFTS'/)
         write (pun,99003) (kap,kap=-lact-1,-1) , (kap,kap=1,lact)
         write (pun,99004) (tanx(kap),kap=-lact-1,-1) , (tanx(kap),kap=1,lact)
!
      else
!
         write (pun,99001)
         write (pun,99003) (kap,kap=-1,-lact-1,-1)
         write (pun,99004) (ratx(kap),kap=-1,-lact-1,-1)
         write (pun,99006) (kap,kap=1,lact)
         write (pun,99007) (ratx(kap),kap=1,lact)
         write (pun,99005)
99005    format (/1x,'PHASESHIFTS'/)
         write (pun,99003) (kap,kap=-1,-lact-1,-1)
         write (pun,99004) (eta(kap),kap=-1,-lact-1,-1)
         write (pun,99006) (kap,kap=1,lact)
         write (pun,99007) (eta(kap),kap=1,lact)
!
      endif
!
   endif
99001 format (/1x,'CF/G RATIO'/)
99003 format (2x,12(5x,i3,5x))
99004 format (2x,12D13.5)
99006 format (15x,11(6x,i2,5x))
99007 format (15x,11D13.5)
!
end subroutine WAFU
!
!
subroutine comdir(e1,kappa,za,nrc,nnk,dx,x0,q,p,ratfg,nonrel)
!======================================================================
!
!   integration of relativistic radial dirac equations by milne
!   method and calculation of ratio cf/g
!
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NRAD = 250
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(in) :: E1
   integer , intent(in) :: KAPPA
   real(REAL64) , intent(in) , dimension(nrad) :: ZA
   integer , intent(in) :: NRC
   integer , intent(in) :: NNK
   real(REAL64) , intent(in) :: DX
   real(REAL64) , intent(in) :: X0
   real(REAL64) , intent(inout) , dimension(nrad) :: Q
   real(REAL64) , intent(inout) , dimension(nrad) :: P
   real(REAL64) , intent(out) :: RATFG
   integer , intent(in) :: NONREL
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: BGC , C , CIN , E , HOC , STVAL , T , TC , U , UC , UNP , UNP2 , WC , WNP , WNP2 , X , XC , XK , Z2
   real(REAL64) , dimension(nrad) :: BGX , PP , QP
   real(REAL64) , save :: CHUGE , CSMALL , TEST
   integer :: I , IK , JRI , KAP , LKAP , N , NIT
   integer , save :: PUN
   real(REAL64) , dimension(4) :: SXK , SXM
!
! End of declarations rewritten by SPAG
!
   integer :: SPAG_NextBlock_1
   data test/1.E+05/ , pun/99/
   data csmall/274.0746/ , chuge/1.D+16/
   SPAG_NextBlock_1 = 1
   spag_dispatchloop_1: do
      select case (SPAG_NextBlock_1)
      case (1)
!
         do i = 1 , nnk
            bgx(i) = za(i)
         enddo
         kap = kappa
         xk = dfloat(kap)
         jri = nrc
         e = e1
         stval = x0
         z2 = -bgx(1)*exp(stval)
         tc = dexp(stval)
!
!
         if ( nonrel==1 ) then
            if ( kap<0 ) lkap = -kap - 1
            if ( kap>0 ) lkap = kap
            kap = -lkap - 1
            xk = dfloat(kap)
            u = -0.5*z2/(lkap+1.0)
            p(1) = 1.E-20
            q(1) = u*1.E-20
!
         else
!
            c = csmall
            cin = 1./(c*c)
!
            hoc = z2/c
            if ( dabs(hoc/xk)<=0.05 ) then
               u = (xk+dabs(xk))/hoc - 0.5*hoc/dabs(xk)
            else
               u = (xk+dsqrt(xk*xk-hoc*hoc))/hoc
            endif
            p(1) = 1.0E-20
            q(1) = c*u*1.0E-20
!
         endif
!
!     eq. 4.92
         if ( nonrel/=1 ) then
            pp(1) = tc*(cin*(e-bgx(1))+1.)*q(1) - xk*p(1)
         else
            pp(1) = tc*q(1) - xk*p(1)
         endif
!     eq. 4.90
!
         qp(1) = xk*q(1) - tc*(e-bgx(1))*p(1)
!
!     eq. 4.91
!     runge-kutta method
!
         x = stval
         n = 1
         SPAG_NextBlock_1 = 2
      case (2)
         ik = 0
         xc = x
         bgc = bgx(n)
         wc = q(n)
         uc = p(n)
         do
            ik = ik + 1
            t = dexp(xc)
!
            if ( nonrel/=1 ) then
               sxk(ik) = dx*(-xk*uc+t*wc*(cin*(e-bgc)+1.))
            else
               sxk(ik) = dx*(-xk*uc+t*wc)
            endif
!
            sxm(ik) = dx*(xk*wc-t*(e-bgc)*uc)
            if ( ik==2 ) then
               uc = uc + 0.5*(sxk(2)-sxk(1))
               wc = wc + 0.5*(sxm(2)-sxm(1))
            elseif ( ik==3 ) then
               xc = xc + 0.5*dx
               uc = uc + sxk(3) - 0.5*sxk(2)
               wc = wc + sxm(3) - 0.5*sxm(2)
               bgc = bgx(n+1)
            elseif ( ik==4 ) then
               q(n+1) = q(n) + (sxm(1)+2.0*sxm(2)+2.0*sxm(3)+sxm(4))/6.0
               p(n+1) = p(n) + (sxk(1)+2.0*sxk(2)+2.0*sxk(3)+sxk(4))/6.0
!
               if ( nonrel/=1 ) then
                  pp(n+1) = t*q(n+1)*(cin*(e-bgc)+1.0) - xk*p(n+1)
               else
                  pp(n+1) = t*q(n+1) - xk*p(n+1)
               endif
!
               qp(n+1) = xk*q(n+1) - t*(e-bgc)*p(n+1)
               x = x + dx
               n = n + 1
               if ( n>=6 ) then
                  SPAG_NextBlock_1 = 3
                  cycle spag_dispatchloop_1
               endif
               SPAG_NextBlock_1 = 2
               cycle spag_dispatchloop_1
            else
               xc = xc + 0.5*dx
               uc = uc + 0.5*sxk(1)
               wc = wc + 0.5*sxm(1)
               bgc = 0.5*(bgc+bgx(n+1))
            endif
         enddo
         SPAG_NextBlock_1 = 3
      case (3)
!
!     milne method
!
         x = x + dx
         t = dexp(x)
!
         unp = p(n-5) + 0.3*dx*(11.*pp(n)-14.*pp(n-1)+26.*pp(n-2)-14.*pp(n-3)+11.0*pp(n-4))
!
         wnp = q(n-5) + 0.3*dx*(11.0*qp(n)-14.0*qp(n-1)+26.0*qp(n-2)-14.0*qp(n-3)+11.0*qp(n-4))
         nit = 0
         SPAG_NextBlock_1 = 4
      case (4)
!
!
         if ( nonrel/=1 ) then
            pp(n+1) = t*(cin*(e-bgx(n+1))+1.0)*wnp - xk*unp
         else
            pp(n+1) = t*wnp - xk*unp
         endif
!
         qp(n+1) = xk*wnp - t*(e-bgx(n+1))*unp
!
         unp2 = p(n-3) + (7.0*pp(n+1)+32.0*pp(n)+12.0*pp(n-1)+32.0*pp(n-2)+7.0*pp(n-3))*2.0*dx/45.0
!
         wnp2 = q(n-3) + (7.0*qp(n+1)+32.0*qp(n)+12.0*qp(n-1)+32.0*qp(n-2)+7.0*qp(n-3))*2.0*dx/45.0
!
         if ( dabs(test*(unp2-unp))<=dabs(unp2) ) then
            if ( dabs(test*(wnp2-wnp))<=dabs(wnp2) ) then
               SPAG_NextBlock_1 = 5
               cycle spag_dispatchloop_1
            endif
         endif
!
         if ( nit/=5 ) then
            nit = nit + 1
            wnp = wnp2
            unp = unp2
            SPAG_NextBlock_1 = 4
            cycle spag_dispatchloop_1
         endif
         SPAG_NextBlock_1 = 5
      case (5)
!
         q(n+1) = wnp2
         p(n+1) = unp2
         n = n + 1
         if ( n<nnk ) then
            SPAG_NextBlock_1 = 3
            cycle spag_dispatchloop_1
         endif
         ratfg = q(jri)/p(jri)
         exit spag_dispatchloop_1
      end select
   enddo spag_dispatchloop_1
!
end subroutine COMDIR
!
!
subroutine cei(f2,f4,f1,f3,r,n,lambda,l2,l4,l1,l3,dx,rws,result)
! ====================================================================
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NRAD = 250 , NDIM = nrad
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(in) , dimension(1) :: F2
   real(REAL64) , intent(in) , dimension(1) :: F4
   real(REAL64) , intent(in) , dimension(1) :: F1
   real(REAL64) , intent(in) , dimension(1) :: F3
   real(REAL64) , dimension(1) :: R
   integer , intent(in) :: N
   integer , intent(in) :: LAMBDA
   real(REAL64) , intent(in) :: L2
   real(REAL64) , intent(in) :: L4
   real(REAL64) , intent(in) :: L1
   real(REAL64) , intent(in) :: L3
   real(REAL64) :: DX
   real(REAL64) :: RWS
   real(REAL64) , intent(out) :: RESULT
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) , dimension(ndim) :: ARG1 , ARG2 , ARG3 , ARG4 , RL , RMLP1
   integer :: I , M
   real(REAL64) :: LL
   real(REAL64) , external :: SINTG
   external DINTG
!
! End of declarations rewritten by SPAG
!
!
!     calculate coulomb and exchange integrals
!     (not too) closely following zare, wood and company
!
!        p. marksteiner
!
!  result = int( int( f2(R) *f4(R) * R_<**lambda/R_>**(lambda+1) *
!                     f1(R')*f3(R') ))
   m = n
   do i = 1 , m
      rl(i) = r(i)**lambda
      rmlp1(i) = 1./(rl(i)*r(i))
      arg1(i) = f2(i)*f4(i)*rl(i)
   enddo
!
   ll = l2 + l4 + lambda
   call dintg(ll,1,arg1,arg2,r,dx,m,rws)
! arg2 contains int from r=0 to r(i)
   ll = ll + 1
!
   do i = 1 , m
      arg2(i) = arg2(i)*f1(i)*f3(i)*rmlp1(i)
   enddo
!
   ll = ll + l1 + l3 - lambda - 1
!
   do i = 1 , m
      arg3(i) = f2(i)*f4(i)*rmlp1(i)
   enddo
!
   ll = l2 + l4 - lambda - 1
   call dintg(ll,-1,arg3,arg4,r,dx,m,rws)
! arg4 contains int. from r(i) to rws
   ll = ll + 1
!
   do i = 1 , m
      arg4(i) = arg4(i)*f1(i)*f3(i)*rl(i)
   enddo
!
   ll = ll + l1 + l3 + lambda
   do i = 1 , m
      arg2(i) = arg2(i) + arg4(i)
   enddo
! ll = l1+l2+l3+l4 for both arg2 and arg4
   result = sintg(ll,arg2,r,dx,m)
!
end subroutine CEI
!
!
function sintg(ll,fct,r,dx,n)
! ========================================
!
! integration of fct over r using simpson integration
! r(i)=r0*exp(dx*(i-1))
!
   use iso_fortran_env
   implicit none
!
! Function and Dummy argument declarations rewritten by SPAG
!
   real(REAL64) :: SINTG
   real(REAL64) , intent(in) :: LL
   real(REAL64) , intent(in) , dimension(1) :: FCT
   real(REAL64) , intent(in) , dimension(1) :: R
   real(REAL64) , intent(in) :: DX
   integer , intent(in) :: N
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: FACT , RLL , SUM , X0 , X1 , X2
   integer :: I
   real(REAL64) , save :: TINY
!
! End of declarations rewritten by SPAG
!
!
!
   data tiny/1.0D-20/
!
   if ( ll+1.0<=tiny ) then
      write (9,*) ' sintg: ll+1=' , ll + 1
      sum = 0.
   else
      rll = r(1)**ll
      if ( rll<=tiny ) then
         sum = 0.
      else
         fact = fct(1)/rll
         sum = fact*(r(1)**(ll+1))/(ll+1)
      endif
   endif
!
   x0 = fct(1)*r(1)
   do i = 3 , n , 2
      x1 = fct(i-1)*r(i-1)
      x2 = fct(i)*r(i)
      sum = sum + dx*(x0+4.D0*x1+x2)
      x0 = x2
   enddo
   sintg = sum/3.D0
   if ( mod(n,2)==0 ) sintg = sintg + (fct(n)+fct(n-1))/2.D0*(r(n)-r(n-1))
end function SINTG
!
!
subroutine dintg(ll,idir,fct,yint,r,dx,n,rws)
! =================================================
!
! integration of fct(i) to yield yint(i)
! r(i)=r0*exp(dx*(i-1))
!
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(in) :: LL
   integer , intent(in) :: IDIR
   real(REAL64) , intent(in) , dimension(1) :: FCT
   real(REAL64) , intent(inout) , dimension(1) :: YINT
   real(REAL64) , dimension(1) :: R
   real(REAL64) , intent(in) :: DX
   integer , intent(in) :: N
   real(REAL64) :: RWS
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: CORR , FACT , RLL , SUM , X0 , X1
   integer :: I , J
   real(REAL64) , save :: TINY
   external INTER1
!
! End of declarations rewritten by SPAG
!
   data tiny/1.0D-20/
!
   if ( idir>0 ) then
!
      if ( ll+1.0<=tiny ) then
         write (9,*) ' dintg: ll+1=' , ll + 1
         sum = 0.
      else
         rll = r(1)**ll
         if ( rll<=tiny ) then
            sum = 0.
         else
            fact = fct(1)/rll
            sum = fact*(r(1)**(ll+1))/(ll+1)
         endif
      endif
      yint(1) = sum
!
      x0 = fct(1)*r(1)
      do i = 2 , n
         x1 = fct(i)*r(i)
! trapezoidal rule:
         sum = sum + dx*(x0+x1)/2.D0
         yint(i) = sum
         x0 = x1
      enddo
!
   else
!
      yint(n) = 0.
      sum = 0.
      x1 = fct(n)*r(n)
!
      do i = n - 1 , 1 , -1
         x0 = fct(i)*r(i)
         sum = sum + dx*(x0+x1)/2.D0
         yint(i) = sum
         if ( i==n-3 ) then
            call inter1(r(n-3),yint(n-3),4,1,rws,corr)
            sum = sum - corr
            do j = n - 3 , n
               yint(j) = yint(j) - corr
            enddo
         endif
         x1 = x0
      enddo
!
   endif
end subroutine DINTG
!
!
subroutine inter1(r,p,n,id,rs,ps)
! =====================================
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: N
   real(REAL64) , intent(in) , dimension(n) :: R
   real(REAL64) , intent(in) , dimension(n) :: P
   integer , intent(in) :: ID
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
   ps = 0.D0
   do j = 1 , n , id
      term = 1.D0
      denom = 1.D0
      do i = 1 , n , id
         if ( i/=j ) then
            denom = denom*(r(j)-r(i))
            term = term*(rs-r(i))
         endif
      enddo
      ps = ps + term*p(j)/denom
   enddo
end subroutine INTER1
!
!
subroutine sbf(e,r,fb,fn)
! =============================
!
! spherical bessel and neuman functions for            e > 0
! modified spherical bessel and neuman functions for   e < 0
!
   use iso_fortran_env
   implicit none
   include 'l.par'
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(in) :: E
   real(REAL64) , intent(in) :: R
   real(REAL64) , intent(out) , dimension(0:lmax+1) :: FB
   real(REAL64) , intent(out) , dimension(0:lmax+1) :: FN
!
! Local variable declarations rewritten by SPAG
!
   complex(REAL64) , dimension(0:lmax+1) :: CFB , CFN
   complex(REAL64) , save :: CIMMI
   complex(REAL64) :: COSX , ECOMPL , EROOT , SINX , X1 , X2
   integer :: L
   real(REAL64) :: TL
!
! End of declarations rewritten by SPAG
!
   data cimmi/(0.0D0,-1.0D0)/
!
   ecompl = dcmplx(e,0.D0)
   eroot = cdsqrt(ecompl)
   x1 = eroot*r
   x2 = x1*x1
   sinx = cdsin(x1)
   cosx = cdcos(x1)
!
   cfb(0) = sinx/x1
   cfb(1) = sinx/x2 - cosx/x1
   cfn(0) = -cosx/x1
   cfn(1) = -cosx/x2 - sinx/x1
!
   do l = 2 , lmax + 1
      tl = dfloat(2*l-1)
      cfb(l) = tl*cfb(l-1)/x1 - cfb(l-2)
      cfn(l) = tl*cfn(l-1)/x1 - cfn(l-2)
   enddo
!
   if ( e<0. ) then
      cfn(0) = cfn(0)*cimmi
      do l = 1 , lmax + 1
         cfb(l) = cfb(l)*cimmi**l
         cfn(l) = cfn(l)*cimmi**(l+1)
      enddo
   endif
!
   do l = 0 , lmax + 1
      fb(l) = dreal(cfb(l))
      fn(l) = dreal(cfn(l))
   enddo
!
end subroutine SBF
!
!
subroutine sbf1(e,r,fb,fn)
! =============================
!
! spherical bessel and neuman functions for            e > 0
! modified spherical bessel and neuman functions for   e < 0
!
   use iso_fortran_env
   implicit none
   include 'l.par'
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(in) :: E
   real(REAL64) , intent(in) :: R
   real(REAL64) , intent(inout) , dimension(0:lmax+1) :: FB
   real(REAL64) , intent(inout) , dimension(0:lmax+1) :: FN
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: COSX , EROOT , SINX , TL , X1 , X2
   integer :: L
!
! End of declarations rewritten by SPAG
!     shyp(x)=0.5d0*(dexp(x)-dexp(-x))
!     chyp(x)=0.5d0*(dexp(x)+dexp(-x))
!
   if ( e>0. ) then
!
      eroot = dsqrt(e)
      x1 = eroot*r
      x2 = x1*x1
      sinx = dsin(x1)
      cosx = dcos(x1)
!
      fb(0) = sinx/x1
      fb(1) = sinx/x2 - cosx/x1
      fn(0) = -cosx/x1
      fn(1) = -cosx/x2 - sinx/x1
!
      do l = 2 , lmax + 1
         tl = dfloat(2*l-1)
         fb(l) = tl*fb(l-1)/x1 - fb(l-2)
         fn(l) = tl*fn(l-1)/x1 - fn(l-2)
      enddo
!
   else
!
      eroot = dsqrt(-e)
      x1 = eroot*r
      x2 = x1*x1
      sinx = dsinh(x1)
      cosx = dcosh(x1)
!
      fb(0) = sinx/x1
      fb(1) = -sinx/x2 + cosx/x1
      fn(0) = cosx/x1
      fn(1) = -cosx/x2 + sinx/x1
!
      do l = 2 , lmax + 1
         tl = dfloat(2*l-1)
         fb(l) = -tl*fb(l-1)/x1 + fb(l-2)
         fn(l) = -tl*fn(l-1)/x1 + fn(l-2)
      enddo
!
   endif
!
end subroutine SBF1
!
!
subroutine wig3j(w3j,pun)
! ===================================
!
   use iso_fortran_env
   implicit none
   include 'l.par'
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NDIM = 100
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(out) , dimension(-lmax-1:lmax,-lmax-1:lmax,0:lmax2) :: W3J
   integer , intent(in) :: PUN
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) , save :: HALF , HALFM , SMALL
   integer :: I1 , IER , KAP1 , KAP2 , L1 , L2 , LAM , LAMAX , LAMIN
   integer , save :: IWIGTES
   real(REAL64) :: J1 , J2
   real(REAL64) , save :: NUL
   real(REAL64) :: RESULT , XLMAT , XLMAX , XLMIN
   real(REAL64) , dimension(ndim) :: THRCOF
   external REC3JJ
!
! End of declarations rewritten by SPAG
!
   data nul , small , half , halfm/0.D0 , 0.01D0 , 0.5D0 , -0.5D0/
   data iwigtes/0/
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
   do kap1 = -lmax - 1 , lmax
      do kap2 = -lmax - 1 , lmax
         do lam = 0 , lmax2
            w3j(kap1,kap2,lam) = nul
         enddo
      enddo
   enddo
!
   do kap1 = -lmax - 1 , lmax
!
      if ( kap1/=0 ) then
         if ( kap1<0 ) then
            j1 = -kap1 - half
            l1 = -kap1 - 1
         else
            j1 = kap1 - half
            l1 = kap1
         endif
!
         do kap2 = -lmax - 1 , lmax
!
            if ( kap2/=0 ) then
               if ( kap2<0 ) then
                  j2 = -kap2 - half
                  l2 = -kap2 - 1
               else
                  j2 = kap2 - half
                  l2 = kap2
               endif
!
               call rec3jj(thrcof,j2,j1,halfm,half,xlmin,xlmax,xlmat,ndim,ier)
               if ( ier>=0 ) then
!
                  lamin = idint(xlmin+small)
                  lamax = idint(xlmax+small)
!
                  do lam = lamin , lamax
!
                     i1 = lam - lamin + 1
                     if ( mod(lam-lamin,2)/=1 ) then
                        w3j(kap1,kap2,lam) = thrcof(i1)
!
                        if ( iwigtes==1 ) then
                           write (pun,99001) kap1 , l1 , j1 , kap2 , l2 , j2 , lam , result
!
99001                      format (2I4,f5.1,5x,2I4,f5.1,5x,i4,10x,d17.10)
                        endif
                     endif
!
                  enddo
               endif
            endif
!
         enddo
      endif
   enddo
!
end subroutine WIG3J
!
!
subroutine wig6j(w6j,j1,pun)
! ================================
!
   use iso_fortran_env
   implicit none
   include 'l.par'
!
! PARAMETER definitions rewritten by SPAG
!
   integer , parameter :: NDIM = 100
!
! Dummy argument declarations rewritten by SPAG
!
   real(REAL64) , intent(out) , dimension(0:lvmax,0:lvmax,0:lmax,0:lmax,0:lmax) :: W6J
   real(REAL64) :: J1
   integer , intent(in) :: PUN
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: COEFF , XLAM , XLAMMA , XLAMMI , XLMAT
   real(REAL64) , save :: HALF , SMALL
   integer :: IER , J3 , L , L2 , LAM , LAMA , LAMA1 , LAMA2 , LAMI , LAMI1 , LAMI2 , LAMMA , LAMMI , LAMP , LP
   integer , save :: IWIGTES
   real(REAL64) :: J , J2 , JP
   real(REAL64) , save :: NUL
   real(REAL64) , dimension(ndim) :: SIXCOF
   external REC6J
!
! End of declarations rewritten by SPAG
!
   data nul , small , half/0.D0 , 0.01D0 , 0.5D0/
   data iwigtes/0/
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
   do l = 0 , lvmax
      do lp = 0 , lvmax
         do l2 = 0 , lmax
            do lam = 0 , lmax
               do lamp = 0 , lmax
                  w6j(l,lp,l2,lam,lamp) = nul
               enddo
            enddo
         enddo
      enddo
   enddo
!
!
   do l = 0 , lvmax
      j = (2.*l+1.)*half
      do lp = 0 , lvmax
         jp = (2.*lp+1.)*half
!
         lami1 = idint(dabs(j1-jp)+small)
         lama1 = idint(j1+jp+small)
!
         do l2 = 0 , lmax
            j2 = (2.*l2+1.)*half
!
            lami2 = idint(dabs(j-j2)+small)
            lama2 = idint(j+j2+small)
            lami = max(lami1,lami2)
            lama = min(lama1,lama2)
!
            do lam = 0 , lmax
               if ( lam>=lami .and. lam<=lama ) then
                  xlam = dfloat(lam)
!
                  call rec6j(sixcof,j1,j,xlam,j2,jp,xlammi,xlamma,xlmat,ndim,ier)
                  if ( ier>=0 ) then
!
                     lammi = idint(xlammi+small)
                     lamma = idint(xlamma+small)
!
                     do lamp = lammi , lamma
                        coeff = 1.0 - 2.0*mod(lam+lamp+1,2)
                        w6j(l,lp,l2,lam,lamp) = coeff*sixcof(lamp-lammi+1)
!
                        if ( iwigtes==1 ) then
                           write (pun,99001) lamp , j1 , j , lam , j2 , j3 , sixcof(lamp-lammi+1)
!
99001                      format (i4,2F5.1,5x,i4,2F5.1,10x,d17.10)
                        endif
                     enddo
                  endif
               endif
!
            enddo
!
         enddo
      enddo
   enddo
!
end subroutine WIG6J
!
!
subroutine rec3jj(thrcof,l2,l3,m2,m3,l1min,l1max,lmatch,ndim,ier)
! =====================================================================
!
!  j1-recursion of 3j-coefficients recursive evaluation of 3j- and
! 6j-coefficients.  k. schulten, r.g. gordon.
! ref. in comp. phys. commun. 11 (1976) 269
!
!
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NDIM
   real(REAL64) , intent(inout) , dimension(ndim) :: THRCOF
   real(REAL64) , intent(in) :: L2
   real(REAL64) , intent(in) :: L3
   real(REAL64) , intent(in) :: M2
   real(REAL64) , intent(in) :: M3
   real(REAL64) , intent(inout) :: L1MIN
   real(REAL64) , intent(inout) :: L1MAX
   real(REAL64) , intent(out) :: LMATCH
   integer , intent(inout) :: IER
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: A1 , A1S , A2 , A2S , C1 , C1OLD , C2 , CNORM , DENOM , DV , OLDFAC , RATIO , SIGN1 , SIGN2 , SUM1 , SUM2 ,     &
                 & SUMBAC , SUMFOR , SUMUNI , THRESH , X , X1 , X2 , X3 , Y , Y1 , Y2 , Y3
   real(REAL64) , save :: EPS , HALF , HUGE , ONE , SRHUGE , SRTINY , THREE , TINY , TWO , ZERO
   integer :: I , INDEX , L1CMAX , L1CMIN , LSTEP , N , NFIN , NFINP1 , NFINP2 , NFINP3 , NLIM , NSTEP2
   real(REAL64) :: L1 , M1 , NEWFAC
!
! End of declarations rewritten by SPAG
!
   integer :: SPAG_NextBlock_1
!
   data zero , eps , half , one/0D0 , 0.01D0 , 0.5D0 , 1D0/
   data two , three/2D0 , 3D0/
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
   data tiny , srtiny/1.0D-10 , 1.0D-05/
!
!  huge should be set close to largest positive floating point
!  number which is representable on the computer.  srhuge is
!  square root of huge .
!
   data huge , srhuge/1.0D10 , 1.0D05/
   SPAG_NextBlock_1 = 1
   spag_dispatchloop_1: do
      select case (SPAG_NextBlock_1)
      case (1)
!
         lmatch = zero
         m1 = -m2 - m3
!
!  check relative magnitude of l- and m-values
         if ( l2-dabs(m2)+eps>=0 ) then
            if ( l3-dabs(m3)+eps>=0 ) then
               if ( dmod(l2+dabs(m2)+eps,one)<eps+eps ) then
                  if ( dmod(l3+dabs(m3)+eps,one)<eps+eps ) then
!
!  limits for l1
!
                     l1min = dmax1(dabs(l2-l3),dabs(m1))
                     l1max = l2 + l3
                     if ( l1min<l1max-eps ) then
!
!
!
                        ier = 0
                        nfin = idint(l1max-l1min+one+eps)
                        if ( ndim<nfin ) then
!
!  dimension of thrcof not large enough to hold all the coefficients
!  required
!
                           ier = -2
                           write (9,99001) l2 , l3 , m1 , m2 , m3 , nfin , ndim
99001                      format (///1x,'3j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,'exceed storage provided  (',i4,',',i4,')')
                           return
                        else
!
!
!  starting forward recursion from l1min taking nstep1 steps
!
                           l1 = l1min
                           thrcof(1) = srtiny
                           sum1 = (l1+l1+one)*tiny
!
!
                           lstep = 1
                           SPAG_NextBlock_1 = 2
                           cycle spag_dispatchloop_1
                        endif
                     elseif ( l1min<l1max+eps ) then
!
!
!  this is reached in case that l1 can take only one value,
!  i.e. l1min = l1max
!
                        ier = 0
                        thrcof(1) = (-one)**idint(dabs(l2+m2-l3+m3)+eps)/dsqrt(l1min+l2+l3+one)
                        l1cmin = l1min
                        l1cmax = l1max
                        return
                     endif
                  endif
               endif
            endif
         endif
!
!
!  this is reached if l2-/m2/ and l3-/m3/  less than zero or not integer
!
         ier = -1
         write (9,99002) l2 , l3 , m1 , m2 , m3
99002    format (///1x,'3j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,'do not satisfy the condition l2-/m2/ and l3-/m3/ ge zero ',   &
                &'and integer')
         return
      case (2)
         lstep = lstep + 1
         l1 = l1 + one
!
!
         oldfac = newfac
         a1 = (l1+l2+l3+one)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+one)
         a2 = (l1+m1)*(l1-m1)
         newfac = dsqrt(a1*a2)
         if ( l1<one+eps ) then
!
!  if l1 = 1  (l1-1) has to be factored out of dv, hence
!
            c1 = -(l1+l1-one)*l1*(m3-m2)/newfac
         else
!
!
            dv = -l2*(l2+one)*m1 + l3*(l3+one)*m1 + l1*(l1-one)*(m3-m2)
            denom = (l1-one)*newfac
!
!
            if ( lstep>2 ) c1old = dabs(c1)
            c1 = -(l1+l1-one)*dv/denom
         endif
!
         if ( lstep>2 ) then
!
!
            c2 = -l1*oldfac/denom
!
!  recursion to the next 3j-coefficient x
!
            x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
            thrcof(lstep) = x
            sumfor = sum1
            sum1 = sum1 + (l1+l1+one)*x*x
            if ( lstep/=nfin ) then
!
!  see if last unnormalized 3j-coefficient exceeds srhuge
!
               if ( dabs(x)>=srhuge ) then
!
!  this is reached if last 3j-coefficient larger than srhuge
!  so that the recursion series thrcof(1), ... , thrcof(lstep)
!  has to be rescaled to prevent overflow
!
                  ier = ier + 1
                  do i = 1 , lstep
                     if ( dabs(thrcof(i))<srtiny ) thrcof(i) = zero
                     thrcof(i) = thrcof(i)/srhuge
                  enddo
                  sum1 = sum1/huge
                  sumfor = sumfor/huge
                  x = x/srhuge
               endif
!
!  as long as /c1/ is decreasing the recursion proceeds towards
!  increasing 3j-values and, hence, is numerically stable.  once
!  an increase of /c1/ is detected the recursion direction is
!  reversed.
!
               if ( c1old>dabs(c1) ) then
                  SPAG_NextBlock_1 = 2
                  cycle spag_dispatchloop_1
               endif
            endif
!
!
!  keep three 3j-coefficients around lmatch for comparision with
!  backward recursion.
!
            lmatch = l1 - 1
            x1 = x
            x2 = thrcof(lstep-1)
            x3 = thrcof(lstep-2)
            nstep2 = nfin - lstep + 3
!
!  starting backward recursion from l1max taking nstep2 steps, so
!  that forward and backward recursion overlap at three points
!  l1 = lmatch+1, lmatch, lmatch-1.
!
            nfinp1 = nfin + 1
            nfinp2 = nfin + 2
            nfinp3 = nfin + 3
            l1 = l1max
            thrcof(nfin) = srtiny
            sum2 = tiny*(l1+l1+one)
!
            l1 = l1 + two
            lstep = 1
            do
               lstep = lstep + 1
               l1 = l1 - one
!
               oldfac = newfac
               a1s = (l1+l2+l3)*(l1-l2+l3-one)*(l1+l2-l3-one)*(-l1+l2+l3+two)
               a2s = (l1+m1-one)*(l1-m1-one)
               newfac = dsqrt(a1s*a2s)
!
               dv = -l2*(l2+one)*m1 + l3*(l3+one)*m1 + l1*(l1-one)*(m3-m2)
!
               denom = l1*newfac
               c1 = -(l1+l1-one)*dv/denom
               if ( lstep>2 ) then
!
!
                  c2 = -(l1-one)*oldfac/denom
!
!  recursion to the next 3j-coefficient y
!
                  y = c1*thrcof(nfinp2-lstep) + c2*thrcof(nfinp3-lstep)
!
                  if ( lstep==nstep2 ) then
!
!
!  the forward recursion 3j-coefficients x1, x2, x3 are to be matched
!  with the corresponding backward recursion values y1, y2, y3.
!
                     y3 = y
                     y2 = thrcof(nfinp2-lstep)
                     y1 = thrcof(nfinp3-lstep)
!
!
!  determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
!  with minimal error.
!
                     ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
                     nlim = nfin - nstep2 + 1
!
                     if ( dabs(ratio)<one ) then
!
                        nlim = nlim + 1
                        ratio = one/ratio
                        do n = nlim , nfin
                           thrcof(n) = ratio*thrcof(n)
                        enddo
                        sumuni = sumfor + ratio*ratio*sumbac
                     else
!
                        do n = 1 , nlim
                           thrcof(n) = ratio*thrcof(n)
                        enddo
                        sumuni = ratio*ratio*sumfor + sumbac
                     endif
                     SPAG_NextBlock_1 = 3
                     cycle spag_dispatchloop_1
                  else
!
                     thrcof(nfinp1-lstep) = y
                     sumbac = sum2
                     sum2 = sum2 + (l1+l1-three)*y*y
!
!  see if last unnormalized 3j-coefficient exceeds srhuge
!
                     if ( dabs(y)>=srhuge ) then
!
!  this is reached if last 3j-coefficient larger than srhuge
!  so that the recursion series thrcof(nfin), ... ,thrcof(nfin-lstep+1)
!  has to be rescaled to prevent overflow
!
                        ier = ier + 1
                        do i = 1 , lstep
                           index = nfin - i + 1
                           if ( dabs(thrcof(index))<srtiny ) thrcof(index) = zero
                           thrcof(index) = thrcof(index)/srhuge
                        enddo
                        sum2 = sum2/huge
!
!
                        sumbac = sumbac/huge
                     endif
                  endif
               else
!
!  if l1 = l1max + 1  the third term in the recursion formula vanishes
!
                  y = srtiny*c1
                  thrcof(nfin-1) = y
                  sumbac = sum2
!
                  sum2 = sum2 + tiny*(l1+l1-three)*c1*c1
               endif
            enddo
         else
!
!
!  if l1 = l1min + 1  the third term in the recursion equation vanishes
!  , hence
            x = srtiny*c1
            thrcof(2) = x
            sum1 = sum1 + tiny*(l1+l1+one)*c1*c1
            if ( lstep/=nfin ) then
               SPAG_NextBlock_1 = 2
               cycle spag_dispatchloop_1
            endif
!
            sumuni = sum1
         endif
         SPAG_NextBlock_1 = 3
      case (3)
!
!
!  normalize 3j-coefficients
!
         cnorm = one/dsqrt(sumuni)
!
!  sign convention for last 3j-coefficient determines overall phase
!
         sign1 = dsign(one,thrcof(nfin))
         sign2 = (-one)**idint(dabs(l2+m2-l3+m3)+eps)
         if ( sign1*sign2<=0 ) cnorm = -cnorm
!
         if ( dabs(cnorm)<one ) then
!
            thresh = tiny/dabs(cnorm)
            do n = 1 , nfin
               if ( dabs(thrcof(n))<thresh ) thrcof(n) = zero
               thrcof(n) = cnorm*thrcof(n)
            enddo
            exit spag_dispatchloop_1
         endif
!
         do n = 1 , nfin
            thrcof(n) = cnorm*thrcof(n)
         enddo
         return
      end select
   enddo spag_dispatchloop_1
!
end subroutine REC3JJ
!
!
subroutine rec6j(sixcof,l2,l3,l4,l5,l6,l1min,l1max,lmatch,ndim,ier)
! ===================================================================
!
!  j1-recursion of 6j-coefficients
! recursive evaluation of 3j- and
! 6j-coefficients.  k. schulten, r.g. gordon.
! ref. in comp. phys. commun. 11 (1976) 269
!
!
   use iso_fortran_env
   implicit none
!
! Dummy argument declarations rewritten by SPAG
!
   integer , intent(in) :: NDIM
   real(REAL64) , intent(inout) , dimension(ndim) :: SIXCOF
   real(REAL64) , intent(in) :: L2
   real(REAL64) , intent(in) :: L3
   real(REAL64) , intent(in) :: L4
   real(REAL64) , intent(in) :: L5
   real(REAL64) , intent(in) :: L6
   real(REAL64) , intent(inout) :: L1MIN
   real(REAL64) , intent(inout) :: L1MAX
   real(REAL64) , intent(out) :: LMATCH
   integer , intent(inout) :: IER
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) :: A1 , A1S , A2 , A2S , C1 , C1OLD , C2 , CNORM , DENOM , DV , OLDFAC , RATIO , SIGN1 , SIGN2 , SUM1 , SUM2 ,     &
                 & SUMBAC , SUMFOR , SUMUNI , THRESH , X , X1 , X2 , X3 , Y , Y1 , Y2 , Y3
   real(REAL64) , save :: EPS , HALF , HUGE , ONE , SRHUGE , SRTINY , THREE , TINY , TWO , ZERO
   integer :: I , INDEX , L1CMAX , L1CMIN , LSTEP , N , NFIN , NFINP1 , NFINP2 , NFINP3 , NLIM , NSTEP2
   real(REAL64) :: L1 , NEWFAC
!
! End of declarations rewritten by SPAG
!
   integer :: SPAG_NextBlock_1
!
   data zero , eps , half , one/0D0 , 0.01D0 , 0.5D0 , 1D0/
   data two , three/2D0 , 3D0/
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
   data tiny , srtiny/1.0D-10 , 1.0D-05/
!
!  huge should be set close to largest positive floating point
!  number which is representable on the computer.  srhuge is
!  square root of huge .
   data huge , srhuge/1.0D10 , 1.0D05/
   SPAG_NextBlock_1 = 1
   spag_dispatchloop_1: do
      select case (SPAG_NextBlock_1)
      case (1)
!
         lmatch = zero
!
!
!
!  check if 6j-coefficients obey selection rules
!
         if ( dmod(l2+l3+l5+l6+eps,one)<eps+eps ) then
            if ( dmod(l4+l2+l6+eps,one)<eps+eps ) then
!
               if ( l4+l2>=l6 ) then
                  if ( l4-l2+l6>=0 ) then
                     if ( -l4+l2+l6>=0 ) then
!
                        if ( l4+l3>=l5 ) then
                           if ( l4-l3+l5>=0 ) then
                              if ( -l4+l3+l5>=0 ) then
!
!  limits for l1
!
                                 l1min = dmax1(dabs(l2-l3),dabs(l5-l6))
                                 l1max = dmin1(l2+l3,l5+l6)
                                 if ( l1min<l1max-eps ) then
!
!
                                    ier = 0
                                    nfin = idint(l1max-l1min+one+eps)
                                    if ( ndim<nfin ) then
!
!  this is reached if array sixcof not large enough to hold all
!  6j - coefficients required
!
                                       ier = -2
                                       write (9,99001) l2 , l3 , l4 , l5 , l6 , nfin , ndim
99001                                  format (///1x,'6j-coefficients',9x,'l1',2F7.1/20x,3F7.1,4x,'exceed storage provided  (',i4, &
                                         &',',i4,')')
                                       return
                                    else
!
!
!
!
!
!  start of forward recursion
!
!
                                       l1 = l1min
                                       sixcof(1) = srtiny
                                       sum1 = (l1+l1+one)*tiny
!
                                       lstep = 1
                                       SPAG_NextBlock_1 = 2
                                       cycle spag_dispatchloop_1
                                    endif
                                 elseif ( l1min<l1max+eps ) then
!
!
!  this is reached in case that l1 can take only one value
!
                                    ier = 0
                                    sixcof(1) = (-one)**idint(l2+l3+l5+l6+eps)/dsqrt((l1min+l1min+one)*(l4+l4+one))
                                    l1cmin = l1min
                                    l1cmax = l1max
                                    return
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
!
!
!  this is reached if triangular condition not satisfied
!  or if l2+l3+l5+l6 or l2+l4+l6 not integer
!
         ier = -1
         write (9,99002) l2 , l3 , l4 , l5 , l6
99002    format (///1x,'6j-coefficients',9x,'l1',2F7.1,4x,'do not satisfy triangular conditions or'/20x,3F7.1,4x,                  &
                &'l2+l3+l5+l6 or l2+l4+l6 not integer')
         return
      case (2)
         lstep = lstep + 1
         l1 = l1 + one
!
         oldfac = newfac
         a1 = (l1+l2+l3+one)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+one)
         a2 = (l1+l5+l6+one)*(l1-l5+l6)*(l1+l5-l6)*(-l1+l5+l6+one)
         newfac = dsqrt(a1*a2)
!
         if ( l1<one+eps ) then
!
!  if l1 = 1   (l1 - 1) has to be factored out of dv, hence
!
            c1 = -two*(l2*(l2+one)+l5*(l5+one)-l4*(l4+one))/newfac
         else
!
            dv = two*(l2*(l2+one)*l5*(l5+one)+l3*(l3+one)*l6*(l6+one)-l1*(l1-one)*l4*(l4+one))                                     &
               & - (l2*(l2+one)+l3*(l3+one)-l1*(l1-one))*(l5*(l5+one)+l6*(l6+one)-l1*(l1-one))
!
            denom = (l1-one)*newfac
!
!
            if ( lstep>2 ) c1old = dabs(c1)
            c1 = -(l1+l1-one)*dv/denom
         endif
!
         if ( lstep>2 ) then
!
!
            c2 = -l1*oldfac/denom
!
!  recursion to the next 6j - coefficient x
!
            x = c1*sixcof(lstep-1) + c2*sixcof(lstep-2)
            sixcof(lstep) = x
!
            sumfor = sum1
            sum1 = sum1 + (l1+l1+one)*x*x
            if ( lstep/=nfin ) then
!
!  see if last unnormalized 6j-coefficient exceeds srhuge
!
               if ( dabs(x)>=srhuge ) then
!
!  this is reached if last 6j-coefficient larger than srhuge
!  so that the recursion series sixcof(1), ... ,sixcof(lstep)
!  has to be rescaled to prevent overflow
!
                  ier = ier + 1
                  do i = 1 , lstep
                     if ( dabs(sixcof(i))<srtiny ) sixcof(i) = zero
                     sixcof(i) = sixcof(i)/srhuge
                  enddo
                  sum1 = sum1/huge
                  sumfor = sumfor/huge
                  x = x/srhuge
               endif
!
!
!  as long as the coefficient /c1/ is decreasing the recursion proceeds
!  towards increasing 6j-values and, hence, is numerically stable.
!  once an increase of /c1/ is detected, the recursion direction is
!  reversed.
!
               if ( c1old>dabs(c1) ) then
                  SPAG_NextBlock_1 = 2
                  cycle spag_dispatchloop_1
               endif
            endif
!
!
!  keep three 6j-coefficients around lmatch for comparision later
!  with backward recursion.
!
            lmatch = l1 - 1
            x1 = x
            x2 = sixcof(lstep-1)
            x3 = sixcof(lstep-2)
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
            l1 = l1max
!
            sixcof(nfin) = srtiny
            sum2 = (l1+l1+one)*tiny
!
!
            l1 = l1 + two
            lstep = 1
            spag_loop_2_1: do
               lstep = lstep + 1
               l1 = l1 - one
!
               oldfac = newfac
               a1s = (l1+l2+l3)*(l1-l2+l3-one)*(l1+l2-l3-one)*(-l1+l2+l3+two)
               a2s = (l1+l5+l6)*(l1-l5+l6-one)*(l1+l5-l6-one)*(-l1+l5+l6+two)
               newfac = dsqrt(a1s*a2s)
!
               dv = two*(l2*(l2+one)*l5*(l5+one)+l3*(l3+one)*l6*(l6+one)-l1*(l1-one)*l4*(l4+one))                                  &
                  & - (l2*(l2+one)+l3*(l3+one)-l1*(l1-one))*(l5*(l5+one)+l6*(l6+one)-l1*(l1-one))
!
               denom = l1*newfac
               c1 = -(l1+l1-one)*dv/denom
               if ( lstep>2 ) then
!
!
                  c2 = -(l1-one)*oldfac/denom
!
!  recursion to the next 6j - coefficient y
!
                  y = c1*sixcof(nfinp2-lstep) + c2*sixcof(nfinp3-lstep)
                  if ( lstep==nstep2 ) exit spag_loop_2_1
                  sixcof(nfinp1-lstep) = y
                  sumbac = sum2
                  sum2 = sum2 + (l1+l1-three)*y*y
!
!  see if last unnormalized 6j-coefficient exceeds srhuge
!
                  if ( dabs(y)>=srhuge ) then
!
!  this is reached if last 6j-coefficient larger than srhuge
!  so that the recursion series sixcof(nfin), ... ,sixcof(nfin-lstep+1)
!  has to be rescaled to prevent overflow
!
                     ier = ier + 1
                     do i = 1 , lstep
                        index = nfin - i + 1
                        if ( dabs(sixcof(index))<srtiny ) sixcof(index) = zero
                        sixcof(index) = sixcof(index)/srhuge
                     enddo
                     sumbac = sumbac/huge
!
                     sum2 = sum2/huge
                  endif
               else
!
!  if l1 = l1max + 1 the third term in the recursion equation vanishes
!
                  y = srtiny*c1
                  sixcof(nfin-1) = y
                  if ( lstep==nstep2 ) exit spag_loop_2_1
                  sumbac = sum2
                  sum2 = sum2 + (l1+l1-three)*c1*c1*tiny
               endif
            enddo spag_loop_2_1
!
!
!  the forward recursion 6j-coefficients x1, x2, x3 are to be matched
!  with the corresponding backward recursion values y1, y2, y3.
!
            y3 = y
            y2 = sixcof(nfinp2-lstep)
            y1 = sixcof(nfinp3-lstep)
!
!
!  determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
!  with minimal error.
!
            ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
            nlim = nfin - nstep2 + 1
!
            if ( dabs(ratio)<one ) then
!
               nlim = nlim + 1
               ratio = one/ratio
               do n = nlim , nfin
                  sixcof(n) = ratio*sixcof(n)
               enddo
               sumuni = sumfor + ratio*ratio*sumbac
            else
!
               do n = 1 , nlim
                  sixcof(n) = ratio*sixcof(n)
               enddo
               sumuni = ratio*ratio*sumfor + sumbac
            endif
         else
!
!  if l1 = l1min + 1 the third term in recursion equation vanishes
!
            x = srtiny*c1
            sixcof(2) = x
            sum1 = sum1 + tiny*(l1+l1+one)*c1*c1
!
            if ( lstep/=nfin ) then
               SPAG_NextBlock_1 = 2
               cycle spag_dispatchloop_1
            endif
!
            sumuni = sum1
         endif
!
!
!  normalize 6j-coefficients
!
         cnorm = one/dsqrt((l4+l4+one)*sumuni)
!
!  sign convention for last 6j-coeff. determines overall phase
!
         sign1 = dsign(one,sixcof(nfin))
         sign2 = (-one)**idint(l2+l3+l5+l6+eps)
         if ( sign1*sign2<=0 ) cnorm = -cnorm
!
         if ( dabs(cnorm)<one ) then
!
            thresh = tiny/dabs(cnorm)
            do n = 1 , nfin
               if ( dabs(sixcof(n))<thresh ) sixcof(n) = zero
               sixcof(n) = cnorm*sixcof(n)
            enddo
            exit spag_dispatchloop_1
         endif
!
         do n = 1 , nfin
            sixcof(n) = cnorm*sixcof(n)
         enddo
         return
      end select
   enddo spag_dispatchloop_1
!
end subroutine REC6J
