      program rcvvmat
c
c
c ********************************************************************
c  RELATIVISTIC MATRIXELEMENTS FOR CORE-VALENCE-VALENCE AUGER SPECTRA
c ********************************************************************
c
c the highest possible values for momentum, l are set in l.par
c                   (lmax=2*lvmax+lcore+2)
c
      implicit real*8 (a-h,o-z)
c
      parameter (nrad=250)
      parameter (nemax=15, nemax2=2*nemax)
      include 'l.par'
c
      dimension v(nrad),rs(nrad),xs(nrad)
      dimension earr(nemax),earr2(nemax2),irange(nemax,nemax2)
c
      dimension gc1(nrad),fc1(nrad)
      dimension gvwf(nrad,-lmax-1:lmax),fvwf(nrad,-lmax-1:lmax)
      dimension gvwfp(nrad,-lmax-1:lmax),fvwfp(nrad,-lmax-1:lmax)
      dimension gcwf(nrad,-lmax-1:lmax),fcwf(nrad,-lmax-1:lmax)
c
      dimension sd(-lvmax-1:lvmax,-lvmax-1:lvmax)
      dimension sde(-lvmax-1:lvmax,-lvmax-1:lvmax)
      dimension ssum(-lvmax-1:lvmax,-lvmax-1:lvmax,nemax)
      dimension crsecs(-lvmax-1:lvmax,-lvmax-1:lvmax,nemax,nemax2)
c
      real*8 idgg(-lmax-1:lmax,0:lmax), idgf(-lmax-1:lmax,0:lmax),
     $       idfg(-lmax-1:lmax,0:lmax), idff(-lmax-1:lmax,0:lmax),
     $       iegg(-lmax-1:lmax,0:lmax), iegf(-lmax-1:lmax,0:lmax),
     $       iefg(-lmax-1:lmax,0:lmax), ieff(-lmax-1:lmax,0:lmax)
c
      integer rl(-lmax-1:lmax),rlb(-lmax-1:lmax),rj(-lmax-1:lmax)
c
      dimension w3j(-lmax-1:lmax,-lmax-1:lmax,0:lmax2)
      dimension w6j(0:lvmax,0:lvmax,0:lmax,0:lmax,0:lmax)
c
      integer in,wf,po,mat,pun,zatom
c
      character*1 name(40),norb1(5),norb(5)
      character*50 filen
c
      data tiny,small,half/1.d-10,0.01d0,0.5d0/
      data pi/3.1415926535897932d0/
      data ev/13.606/
      data csmall/274.0746/,chuge/1.d+08/
c
c ********************************************************************
c
      in=1
      wf=2
      po=3
      mat=7
      pun=8
c
      open(in,file='rcvvmat.in',status='old')
      open(9,file='rcvvmat.err',status='unknown')
c
 1000 continue
c
      read(in,250,end=1001) filen
      open(pun,file=filen,status='unknown')
c
      read(in,250) filen
      open(po,file=filen,status='unknown')
c
      read(in,250) filen
      open(wf,file=filen,status='unknown')
c
      read(in,250) filen
      open(mat,file=filen,status='unknown')
c
c  fill up fields for the relativistic quantumnumbers
c
      do kap=-lmax-1,-1
        l=-kap-1
        rl(kap) =l
        rlb(kap)=l+1
        rj(kap) =2*l+1
      end do
      do kap=1,lmax
        l=kap
        rl(kap) =l
        rlb(kap)=l-1
        rj(kap) =2*l-1
      end do
c
c read in maximum angular momentum quantumnumber for valence states
c
      read(in,*) lval
      write(pun,200) lval
c
c read in starting energy, increment for the energy panel, number of
c energy points for valence band 
c (with respect to the muffin tin zero)
c the energy is measured in ryd
c
      read(in,*) e0,de,ne
c
c check facility for non-relativistic limit
c
      read(in,*) nonrel
      if(nonrel.eq.1) then
        clight=chuge
      else
        clight=csmall
      endif
c
c  normalization of the valence wavefunctions :
c                  norm eq 0  within the muffin tin
c                  norm eq 1  within the wigner-seitz sphere
c
      read(in,*) norm
      write(pun,202) norm
c       
c read in identifier of core state
c
      read(in,251) (norb1(i),i=1,5)
c
      read(in,*)
c
c read in the name of the species,
c         the atomic number,
c         the first point on the logarithmic mesh,
c         the increment on the log. mesh,
c         the number of point corresponding to the muffin tin zero,
c         the wigner-seitz radius,
c         the lattice constant
c         the value of the muffin tin zero and
c         the potential in form of r*V(r)
c
      read(po,251) (name(i),i=1,40)
      read(po,*) zatom,x0,nmt,dx,rws,a0,zeromt
      read(po,*) (v(j),j=1,nmt)
c
      coef=2.*zatom/clight
      coef=coef*coef
c
c set up the logarithmic and radial mesh
c transform the potential data to v(R)
c
      x=x0
      xs(1)=x0
      rs(1)=dexp(x)
      signv=1.
      if(v(1).gt.0.) signv=-1.
      v(1)=signv*v(1)/rs(1)-zeromt
c
      do j=2,nmt
        x=x+dx
        xs(j)=x
        rs(j)=dexp(x)
        v(j)=signv*v(j)/rs(j)-zeromt
      end do
      rmt=rs(nmt)
c
      nws=0
      do j=nmt+1,nrad
        x=x+dx
        xs(j)=x
        rs(j)=dexp(x)
        if(nws.eq.0.and.rws.lt.rs(j)) nws=j
        v(j)=0.d0
      end do
c
      write(*,251) (name(i),i=1,40)
      write(pun,251) (name(i),i=1,40)
      write(pun,*)
      write(pun,203) zatom,x0,dx,nmt,rmt,nws,rws,zeromt
      write(pun,204) (rs(j),v(j)*rs(j),j=1,nws)
c
c read in the core wavefunctions as r*WF(r) and core energies (hartree)
c (the same radial mesh is supposed to be used as for the potential)
c
   21 continue
c
      read(wf,251) (norb(i),i=1,5)
      read(wf,*) ec1
      ec1=2.*ec1
      read(wf,*) kap1
      read(wf,*) nwf
c
      if(nonrel.eq.1) then
        do i=1,nwf
          read(wf,*) dummy,gc1(i)
          fc1(i)=0.d0
        end do
      else
        do i=1,nwf
          read(wf,*) dummy,gc1(i),fc1(i)
        end do
      endif
c
      if(norb(1).ne.norb1(1).or.norb(2).ne.norb1(2)
     $   .or.norb(3).ne.norb1(3)) goto 21
c
      write(pun,205) (norb(i),i=1,5),ec1,kap1
      write(pun,206) (rs(j),gc1(j),fc1(j),j=1,nws)
c
      write(pun,208)
c
c  set up energy panel for valence band :   earr
c            and       for continuum    :   earr2
c
      ii=0
      do i=1,ne
        earr(i)=e0+(i-1)*de
      end do
      ne2=2*ne-1
      emin2=2.*e0-ec1
      write(pun,201) (earr(i),i=1,ne)
c
      do i=1,ne2
        earr2(i)=emin2+(i-1)*de
        do j=1,ne
          jp=i+1-j
          if(jp.ge.1.and.jp.le.ne) then
            ii=ii+1
            irange(j,i)=ii
          else
            irange(j,i)=0
          endif
        end do
        write(pun,220) earr2(i),(irange(j,i),j=1,ne)
      end do
c
      write(pun,*)
c
c  initialize the angular integration coefficients
c
      l1=rl(kap1)
      j1=rj(kap1)
      xj1=j1*half
      jj1=j1+1
      y1=dsqrt(kap1*kap1-coef)
      if(nonrel.eq.1) y1=dfloat(l1+1)
c
      call wig3j(w3j,pun)
      call wig6j(w6j,xj1,pun)
c
c
c         ************************************
c         * loop for the energy in continuum *
c         ************************************
c
      do nei2=1,ne2
        write(6,*) nei2
c
c compute continuum wavefunctions
c
        e2=earr2(nei2)
c
        kout=0
        if(nei2.eq.ne) kout=1
c
        write(pun,270) e2
        if(kout.eq.1) write(pun,*) 
        if(kout.eq.1) write(pun,*) ' continuum wavefunctions'
c
        call wafu(e2,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,
     >            gcwf,fcwf,norm,2,nonrel,pun,lval,coef,kout)
c
c
c         ***************************************
c         * loop for the energy in valence band *
c         ***************************************
c
        do 31 nei=1,ne
          neip=nei2+1-nei
          if(irange(nei,nei2).eq.0) goto 31
c
          kout1=0
          if(nonrel.eq.1.and.irange(nei,nei2).eq.ne*ne) kout1=1
c
c compute valence wavefunctions
c
          e=earr(nei)
          ep=earr(neip)
c
          if(kout.eq.1) write(pun,*) 
          if(kout.eq.1) write(pun,*) ' valence wavefunctions'
          write(pun,271) e
c
          call wafu(e,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,
     >              gvwf,fvwf,norm,1,nonrel,pun,lval,coef,kout)
c
          call wafu(ep,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,
     >              gvwfp,fvwfp,norm,1,nonrel,pun,lval,coef,0)
c
c         ******************
c         * loop for kappa *
c         ******************
c
          do 35 kap=-lval-1,lval
            if(kap.eq.0) goto 35
            y=dsqrt(kap*kap-coef)
            if(nonrel.eq.1) y=dfloat(rl(kap)+1)
            jj=rj(kap)+1
            l=idint(jj*0.5+small)-1
            lkap=rl(kap)
c
            if(kout1.eq.1) write(pun,*) ' kappa=',kap
c
c         *******************
c         * loop for kappa' *
c         *******************
c
            do 36 kapp=-lval-1,lval
              if(kapp.eq.0) goto 36
              yp=dsqrt(kapp*kapp-coef)
              if(nonrel.eq.1) yp=dfloat(rl(kapp)+1)
              jjp=rj(kapp)+1
              lp=idint(jjp*0.5+small)-1
              lkapp=rl(kapp)
c
              if(kout1.eq.1) write(pun,*) ' kappa"=',kapp
c
c  compute coulomb and exchange integrals for non-vanishing angular
c  integration coefficients
c
              do lam=0,lmax
                do kap2=-lmax-1,lmax
                  idgg(kap2,lam)=0.
                  idfg(kap2,lam)=0.
                  idgf(kap2,lam)=0.
                  idff(kap2,lam)=0.
                  iegg(kap2,lam)=0.
                  iefg(kap2,lam)=0.
                  iegf(kap2,lam)=0.
                  ieff(kap2,lam)=0.
                end do
              end do
c                                                                      
              do lam=0,lmax
c
                w1v=w3j(kap1,kap,lam)
                w1vp=w3j(kap1,kapp,lam)
c
                do 45 kap2=-lmax-1,lmax
                  if(kap2.eq.0) goto 45
                  y2=dsqrt(kap2*kap2-coef)
                  if(nonrel.eq.1) y2=dfloat(rl(kap2)+1)
c
                  wv2=w3j(kap,kap2,lam)
                  wvp2=w3j(kapp,kap2,lam)
c
                  if(dabs(w1v*wvp2).gt.tiny) then
                    call cei(gc1,gvwf(1,kap),gvwfp(1,kapp),gcwf(1,kap2),
     >                      rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                    idgg(kap2,lam)=result
                    call cei(gc1,gvwf(1,kap),fvwfp(1,kapp),fcwf(1,kap2),
     >                       rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                    idgf(kap2,lam)=result
                    call cei(fc1,fvwf(1,kap),gvwfp(1,kapp),gcwf(1,kap2),
     >                       rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                    idfg(kap2,lam)=result
                    call cei(fc1,fvwf(1,kap),fvwfp(1,kapp),fcwf(1,kap2),
     >                       rs,nws,lam,y1,y,yp,y2,dx,rws,result)
                    idff(kap2,lam)=result
                  endif
c
                  if(dabs(w1vp*wv2).gt.tiny) then
                    call cei(gc1,gvwfp(1,kapp),gvwf(1,kap),gcwf(1,kap2),
     >                       rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                    iegg(kap2,lam)=result
                    call cei(gc1,gvwfp(1,kapp),fvwf(1,kap),fcwf(1,kap2),
     >                       rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                    iegf(kap2,lam)=result
                    call cei(fc1,fvwfp(1,kapp),gvwf(1,kap),gcwf(1,kap2),
     >                       rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                    iefg(kap2,lam)=result
                    call cei(fc1,fvwfp(1,kapp),fvwf(1,kap),fcwf(1,kap2),
     >                       rs,nws,lam,y1,yp,y,y2,dx,rws,result)
                    ieff(kap2,lam)=result
                  endif
c
  45            continue
              end do
c
c
c ************* sum up for the direct term
c
c
              sigd=0.
c
              do lam=0,lmax
                llam=2*lam+1
                w1v=w3j(kap1,kap,lam)
                w1v=w1v*w1v
c
                do 55 kap2=-lmax-1,lmax
                  if(kap2.eq.0) goto 55
                  jj2=rj(kap2)+1
                  wvp2=w3j(kapp,kap2,lam)
                  wvp2=wvp2*wvp2
c
                  sum1=idgg(kap2,lam)+idgf(kap2,lam)+
     >                 idfg(kap2,lam)+idff(kap2,lam)
c
                  sigd=sigd+jj1*jj2*w1v*wvp2*sum1*sum1/llam
c
                  if(kout1.eq.1.and.dabs(sum1).gt.tiny) then
                    write(pun,*)  ' nonzero for direct term'
                    write(pun,*)  lam, kap2, sum1
                  endif
c
   55           continue
              end do
c
              sd(kap,kapp)=sigd
c
c
c ************* sum up for the cross term
c
c
              sigde=0.
c
              do lam=0,lmax
                llam=2*lam+1
                w1v=w3j(kap1,kap,lam)
c
                do lamp=0,lmax
                  llamp=2*lam+1
                  w1vp=w3j(kap1,kapp,lamp)
c
                  do 65 kap2=-lmax-1,lmax
                    if(kap2.eq.0) goto 65
                    jj2=rj(kap2)+1
                    l2=idint(jj2*0.5+small)-1
                    wv2=w3j(kap,kap2,lamp)
                    wvp2=w3j(kapp,kap2,lam)
c
                    sum1=idgg(kap2,lam)+idgf(kap2,lam)+
     >                   idfg(kap2,lam)+idff(kap2,lam)
                    sum2=iegg(kap2,lamp)+iegf(kap2,lamp)+
     >                   iefg(kap2,lamp)+ieff(kap2,lamp)
c
                    contr=jj1*jj2*w1v*w1vp*wv2*wvp2*
     >                    sum1*sum2*w6j(l,lp,l2,lam,lamp)
                    sigde=sigde+contr
c
                    if(kout1.eq.1.and.dabs(contr).gt.tiny) then
                      write(pun,*)  ' nonzero for cross term'
                      write(pun,*)  lam, lamp, kap2, sum1*sum2
                    endif
c
   65             continue
                end do
              end do
c
              sde(kap,kapp)=sigde
c
              ssum(kap,kapp,nei)=2.0*(sigd-sigde)
c
   36       continue
   35     continue
c
c
c         ***************************
c         * end of loop for kappa's *
c         ***************************
c
          if(kout1.eq.1) then
c
            write(pun,*) 
            write(pun,*) ' DIRECT TERM'
            do 70 kap=-lval-1,lval
              if(kap.eq.0) goto 70
              write(pun,230) (kap,kapp,kapp=-lval-1,-1),
     >                       (kap,kapp,kapp=1,lval)
              write(pun,231) (sd(kap,kapp),kapp=-lval-1,-1),
     >                       (sd(kap,kapp),kapp=1,lval)
   70       continue
            write(pun,*) 
            write(pun,*) ' CROSS TERM'
            do 71 kap=-lval-1,lval
              if(kap.eq.0) goto 71
              write(pun,230) (kap,kapp,kapp=-lval-1,-1),
     >                       (kap,kapp,kapp=1,lval)
              write(pun,231) (sde(kap,kapp),kapp=-lval-1,-1),
     >                       (sde(kap,kapp),kapp=1,lval)
   71       continue
c
          endif
c
   31   continue
c
c  ************************************ 
c  construct symmetrized cross sections
c  ************************************
c
        do 66 i=1,ne
          ip=nei2+1-i
          if(irange(i,nei2).eq.0) goto 66
          do 67 kap=-lval-1,lval
            if(kap.eq.0) goto 67
            crsecs(kap,kap,i,nei2)=(ssum(kap,kap,i)+ssum(kap,kap,ip))
     >                              *0.5
            do 68 kapp=-lval-1,kap-1
              if(kapp.eq.0) goto 68
              crsecs(kap,kapp,i,nei2)=ssum(kap,kapp,i)+ssum(kapp,kap,ip)
   68       continue
   67     continue
   66   continue
c
c
      end do
c
c         ********************************
c         * end of loop for the energies *
c         ********************************
c
      write(mat,225) (name(i),i=1,40),(norb1(i),i=1,5),
     >               ec1,e0,de,ne,lval,a0
c
      write(mat,221) 
     $ ( ( (kap,kapp), kapp=-lval-1,kap ), kap=-lval-1,-1 ),
     $ ( ( (kap,kapp), kapp=-lval-1,-1  ),
     $   ( (kap,kapp), kapp=1,kap       ), kap=1,lval)
c
      do 99 i2=1,ne2
        do 98 i=1,ne
          ip=i2+1-i
          if(irange(i,i2).eq.0) goto 98
          write(mat,222) i2,i,earr2(i2),earr(ip),earr(i),
     > ( (crsecs(kap,kapp,i,i2),kapp=-lval-1,kap), kap=-lval-1,-1),
     > ( (crsecs(kap,kapp,i,i2),kapp=-lval-1,-1 ),
     >   (crsecs(kap,kapp,i,i2),kapp=1,kap      ), kap=1,lval)
   98   continue
   99 continue
c
c
  250 format(a50)
  200 format(1x,'LVAL=',t20,i4/)
  202 format(1x,'NORMALIZATION OF VALENCE WF. (0 FOR MT, 1 FOR WS):',
     $i4/)
  251 format(40a1)
  203 format(1x,' ZATOM =',t20,i4/
     $       1x,'    X0 =',t20,f5.1/
     $       1x,'    DX =',t20,f10.7/
     $       1x,'   NMT =',t20,i4/
     $       1x,'   RMT =',t20,f10.7/
     $       1x,'   NWS =',t20,i4/
     $       1x,'   RWS =',t20,f10.7/
     $       1x,'ZEROMT =',t20,f5.1//
     $       1x,'RADIAL MESH AND POTENTIAL*R'/)
  204 format(4d20.10)
  205 format(//1x,'CORE WAVEFUNCTION'/
     $         1x,'ORB. =',t20,5a1/
     $         1x,' EC1 =',t20,f12.6/
     $         1x,'KAP1 =',t20,i4//)
  206 format(3d20.10)
  208 format(//' ****************** END INPUT *******************'//)
  201 format(1x,'ENERGY PLANE'//(7x,15f7.2)/)
  220 format(f7.2,15(2x,i3,2x))
  270   format(' e2=',f15.5)
  271     format(t20,' e =',f15.5)
  230 format(7(4x,i2,1x,i2,4x))
  231 format(7d13.5)
  225 format(1x,'RELATIVISTIC CORE-VALENCE-VALENCE AUGER 
     $MATRIXELEMENTS'/
     $       40a1/
     $       1x,'CORE STATE'/
     $       5a1/
     $       1x,'ENERGY OF THE CORE STATE'/
     $       f12.5/
     $       1x,'E0, DE, NE FOR VALENCE BAND'/
     $       2f6.2,i5/
     $       1x,'MAXIMAL L QUANTUMNUMBER FOR VALENCE STATES'/
     $       i1/
     $       1x,'LATTICE CONSTANT'/
     $       f10.6/)
  221 format(1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,
     $28(i6,i3,'    '))
  222 format(2i3,3f6.2,28d13.5)
c
      close(pun)
      close(po)
      close(wf)
      close(mat)
      goto 1000
c
 1001 stop
      end
      subroutine wafu(en,rl,rlb,rj,v,rs,x0,dx,rmt,nmt,rws,nws,
     $p,q,norm,nval,nonrel,pun,lval,coef,kout)
c================================================
c
      implicit real*8 (a-h,o-z)
c
      parameter (nrad=250)
      include 'l.par'
c
      dimension v(nrad),rs(nrad)
c
      dimension p(nrad,-lmax-1:lmax),q(nrad,-lmax-1:lmax)
      dimension rr(nrad),vint(-lmax-1:lmax),eta(-lmax-1:lmax)
      dimension tanx(-lmax-1:lmax),ratx(-lmax-1:lmax)
      dimension fb(0:lmax+1),fn(0:lmax+1),fb1(0:lmax+1),fn1(0:lmax+1)
c
      integer rl(-lmax-1:lmax),rlb(-lmax-1:lmax),rj(-lmax-1:lmax)
      integer pun
c
      data csmall/274.0746/,chuge/1.d+08/
c
c
      if(nonrel.eq.1) then
      clight=chuge
      else
      clight=csmall
      endif
c
      if(en.gt.0) then
      sign=1.
      ekappa = dsqrt ( en )
      else
      sign=-1.
      ekappa = dsqrt ( -en )
      endif
c
      call sbf1(en,rmt,fb,fn)
c
      if(nval.eq.1) lact=lval
      if(nval.eq.2) lact=lmax
c
      do 10 kap=-lact-1,lact
      if(kap.eq.0) goto 10
c
      l=rl(kap)
      lb=rlb(kap)
      if(kap.gt.0)  sk= ekappa
      if(kap.lt.0)  sk= - sign * ekappa
      yk=dsqrt(kap*kap-coef)
c
      if(nonrel.eq.1) then
      yk=dfloat(l+1)
      lb=l+1
      sk=-sign*ekappa
      endif
c
      yk2=yk+yk
c
      call comdir(en,kap,v,nmt,nws,dx,x0,q(1,kap),p(1,kap),ratfg,
     $            nonrel)
c
      tane=(ratfg*fb(l)-sk*fb(lb)) / (ratfg*fn(l)-sk*fn(lb))
      ratx(kap)=ratfg
      tanx(kap)=tane
c
c valence or continuum normalization
c
      if(nval.eq.1) then
      a1=sk*ekappa*(fn(lb) - fb(lb)/tane)*rmt/q(nmt,kap)/clight
      b1=   ekappa*(fn(l ) - fb(l )/tane)*rmt/p(nmt,kap)
      else
      eta(kap)=datan(tane)
  310 cose=dcos(eta(kap))
      sine=dsin(eta(kap))
      a1=sk*(cose*fb(lb)-sine*fn(lb))*rmt/q(nmt,kap)/clight
      b1=   (cose*fb(l )-sine*fn(l ))*rmt/p(nmt,kap)
c     if(b1.lt.0.) then
c       if(eta(kap).gt.pi) then
c       eta(kap)=eta(kap)-pi
c       else
c       eta(kap)=eta(kap)+pi
c       endif
c       goto 310
c     endif
      endif
c
c
      do 15 i=1,nmt
      p(i,kap)=p(i,kap)*b1
      q(i,kap)=q(i,kap)*a1
      rr(i)=p(i,kap)*p(i,kap)+q(i,kap)*q(i,kap)
   15 continue
c
c
c  set up the wavefunctions beyond the muffin tin radius
c
      do 20 i=nmt+1,nws+1
      r=rs(i)
      call sbf1(en,r,fb1,fn1)
      if(nval.eq.1) then
      q(i,kap)=sk*ekappa*(fn1(lb) - fb1(lb)/tane)*r/clight
      p(i,kap)=   ekappa*(fn1(l ) - fb1(l )/tane)*r
      else
      q(i,kap)=sk*(cose*fb1(lb) - sine*fn1(lb))*r/clight
      p(i,kap)=   (cose*fb1(l ) - sine*fn1(l ))*r
      endif
      rr(i)=p(i,kap)*p(i,kap)+q(i,kap)*q(i,kap)
   20 continue
c
c
      if(nval.eq.2) goto 10
c
c normalize valence wavefuntions according to norm
c
      if(norm.eq.0) then
      vint(kap)=sintg(yk2,rr,rs,dx,nmt)
      else
      x1=sintg(yk2,rr,rs,dx,nws-1)
      x2=sintg(yk2,rr,rs,dx,nws)
      vint(kap)=x1+(x2-x1)*(rws-rs(nws-1))/(rs(nws)-rs(nws-1))
      endif
      qint=dsqrt(1./vint(kap))
      do 25 i=1,nrad
      p(i,kap)=p(i,kap)*qint
      q(i,kap)=q(i,kap)*qint
   25 continue
c
c
   10 continue
c
c
      if(kout.eq.1) then
c
      if(nval.eq.1) then
c
      write(pun,90)
   90 format(/1x,'CF/G RATIO'/)
      write(pun,101) (kap,kap=-lact-1,-1),(kap,kap=1,lact)
      write(pun,102) (ratx(kap),kap=-lact-1,-1),
     $               (ratx(kap),kap=1,lact)
      write(pun,100)
  100 format(/1x,'TANGENT PHASESHIFTS'/)
      write(pun,101) (kap,kap=-lact-1,-1),(kap,kap=1,lact)
  101 format(2x,12(5x,i3,5x))
      write(pun,102) (tanx(kap),kap=-lact-1,-1),
     $               (tanx(kap),kap=1,lact)
  102 format(2x,12d13.5)
c
      else
c
      write(pun,90)
      write(pun,101) (kap,kap=-1,-lact-1,-1)
      write(pun,102) (ratx(kap),kap=-1,-lact-1,-1)
      write(pun,104) (kap,kap=1,lact)
      write(pun,105) (ratx(kap),kap=1,lact)
      write(pun,103)
  103 format(/1x,'PHASESHIFTS'/)
      write(pun,101) (kap,kap=-1,-lact-1,-1)
      write(pun,102) (eta(kap),kap=-1,-lact-1,-1)
      write(pun,104) (kap,kap=1,lact)
      write(pun,105) (eta(kap),kap=1,lact)
  104 format(15x,11(6x,i2,5x))
  105 format(15x,11d13.5)
c
      endif
c
      endif
c
      return
      end
      subroutine comdir(e1,kappa,za,nrc,nnk,dx,x0,q,p,ratfg,nonrel)
c======================================================================
c
c   integration of relativistic radial dirac equations by milne
c   method and calculation of ratio cf/g
c
      implicit real*8 (a-h,o-z)
c
      parameter (nrad=250)
c
      dimension bgx(nrad),sxk(4),sxm(4),p(nrad),q(nrad),
     *          pp(nrad),qp(nrad),za(nrad)
      integer pun
      data test/1.e+05/,pun/99/
      data csmall/274.0746/,chuge/1.d+16/
c
      do 200 i=1,nnk
      bgx(i)=za(i)
  200 continue
      kap=kappa
      xk=dfloat(kap)
      jri=nrc
      e=e1
      stval = x0
      z2=-bgx(1)*exp(stval)
      tc=dexp(stval)
c
c
      if(nonrel.eq.1) then
      if(kap.lt.0) lkap=-kap-1
      if(kap.gt.0) lkap=kap
      kap=-lkap-1
      xk=dfloat(kap)
      u=-0.5*z2/(lkap+1.0)
      p(1)=1.e-20
      q(1)=u*1.e-20
c
      else
c
      c=csmall
      cin=1./(c*c)
c
      hoc=z2/c
      if (dabs(hoc/xk)-0.05) 5,5,6
5     u=(xk+dabs(xk))/hoc-0.5*hoc/dabs(xk)
      go to 7
6     u=(xk+dsqrt(xk*xk-hoc*hoc))/hoc
7     p(1)=1.0e-20
      q(1)=c*u*1.0e-20
c
      endif
c      
c
c     eq. 4.92
c
      if(nonrel.ne.1) then
      pp(1)=tc*(cin*(e-bgx(1))+1.)*q(1)-xk*p(1)
      else
      pp(1)=tc*q(1)-xk*p(1)
      endif
c
c     eq. 4.90
c
      qp(1)=xk*q(1)-tc*(e-bgx(1))*p(1)
c
c     eq. 4.91
c     runge-kutta method
c
11    x=stval
      n=1
25    ik=0
      xc=x
      bgc=bgx(n)
      wc=q(n)
      uc=p(n)
20    ik=ik+1
      t=dexp(xc)
12    continue
c
      if(nonrel.ne.1) then 
      sxk(ik)=dx*(-xk*uc+t*wc*(cin*(e-bgc)+1.))
      else
      sxk(ik)=dx*(-xk*uc+t*wc)
      endif
c
      sxm(ik)=dx*(xk*wc-t*(e-bgc)*uc)
15    go to (16,17,18,19),ik
16    xc=xc+0.5*dx
      uc=uc+0.5*sxk(1)
      wc=wc+0.5*sxm(1)
      bgc=0.5*(bgc+bgx(n+1))
      go to 20
17    uc=uc+0.5*(sxk(2)-sxk(1))
      wc=wc+0.5*(sxm(2)-sxm(1))
      go to 20
   18 xc=xc+0.5*dx
      uc=uc+sxk(3)-0.5*sxk(2)
      wc=wc+sxm(3)-0.5*sxm(2)
      bgc=bgx(n+1)
      go to 20
19    q(n+1) = q(n)+(sxm(1)+2.0*sxm(2)+2.0*sxm(3)
     *        +sxm(4))/6.0
      p(n+1) = p(n)+(sxk(1)+2.0*sxk(2)+2.0*sxk(3)
     *        +sxk(4))/6.0
c
      if(nonrel.ne.1) then
      pp(n+1)= t*q(n+1)*(cin*(e-bgc)+1.0)-xk*p(n+1)
      else
      pp(n+1)= t*q(n+1)-xk*p(n+1)
      endif
c
      qp(n+1)= xk*q(n+1)-t*(e-bgc)*p(n+1)
24    x=x+dx
      n=n+1
      if (n-6) 25,26,26
c
c     milne method
c
26    x=x+dx
      t=dexp(x)
c
27    unp = p(n-5)+0.3*dx*(11.*pp(n)-14.*pp(n-1)
     *      +26.*pp(n-2)-14.*pp(n-3)+11.0*pp(n-4))
c
      wnp = q(n-5)+0.3*dx*(11.0*qp(n)-14.0*qp(n-1)
     *      +26.0*qp(n-2)-14.0*qp(n-3)+11.0*qp(n-4))
      nit=0
c
33    continue
c
      if(nonrel.ne.1) then
      pp(n+1)=t*(cin*(e-bgx(n+1))+1.0)*wnp-xk*unp
      else
      pp(n+1)=t*wnp-xk*unp
      endif
c      
      qp(n+1)=xk*wnp-t*(e-bgx(n+1))*unp
c
      unp2 = p(n-3)+(7.0*pp(n+1)+32.0*pp(n)
     *       +12.0*pp(n-1)+32.0*pp(n-2)
     *       +7.0*pp(n-3))*2.0*dx/45.0
c
      wnp2 = q(n-3)+(7.0*qp(n+1)+32.0*qp(n)
     *       +12.0*qp(n-1)+32.0*qp(n-2)
     *       +7.0*qp(n-3))*2.0*dx/45.0
c
      if(dabs(test*(unp2-unp))-dabs(unp2)) 30,30,31
   30 if(dabs(test*(wnp2-wnp))-dabs(wnp2)) 32,32,31
c
31    if (nit-5) 81,32,81
81    nit=nit+1
      wnp=wnp2
      unp=unp2
      go to 33
c
32    q(n+1)=wnp2
      p(n+1)=unp2
      n=n+1
      if (n-nnk) 26,34,34
   34 ratfg=q(jri)/p(jri)
c
      return
      end
      subroutine cei(f2,f4,f1,f3,r,n,lambda,l2,l4,l1,l3,dx,rws,result)
c ====================================================================
      implicit real*8 (a-h,o-z)
c
c     calculate coulomb and exchange integrals
c     (not too) closely following zare, wood and company
c
c        p. marksteiner
c
c  result = int( int( f2(R) *f4(R) * R_<**lambda/R_>**(lambda+1) *
c                     f1(R')*f3(R') ))
c
      parameter (nrad=250,ndim=nrad)
c
      dimension rl(ndim),rmlp1(ndim)
      dimension arg1(ndim),arg2(ndim),arg3(ndim),arg4(ndim)
c
      dimension f1(1),f2(1),f3(1),f4(1),r(1)
c
      real*8 l1,l2,l3,l4,ll
c
      m=n
      do 1 i=1,m
      rl(i)=r(i)**lambda
      rmlp1(i)=1./(rl(i)*r(i))
      arg1(i)=f2(i)*f4(i)*rl(i)
    1 continue
c
      ll=l2+l4+lambda
      call dintg(ll,1,arg1,arg2,r,dx,m,rws)
c arg2 contains int from r=0 to r(i)
      ll=ll+1
c
      do 2 i=1,m
      arg2(i)=arg2(i)*f1(i)*f3(i)*rmlp1(i)
    2 continue
c
      ll=ll+l1+l3-lambda-1
c
      do 3 i=1,m
      arg3(i)=f2(i)*f4(i)*rmlp1(i)
    3 continue
c
      ll=l2+l4-lambda-1
      call dintg(ll,-1,arg3,arg4,r,dx,m,rws)
c arg4 contains int. from r(i) to rws
      ll=ll+1
c
      do 4 i=1,m
      arg4(i)=arg4(i)*f1(i)*f3(i)*rl(i)
    4 continue
c
      ll=ll+l1+l3+lambda
      do 5 i=1,m
      arg2(i)=arg2(i)+arg4(i)
    5 continue
c ll = l1+l2+l3+l4 for both arg2 and arg4
      result=sintg(ll,arg2,r,dx,m)
c
      return
      end
      double precision function sintg(ll,fct,r,dx,n)
c ========================================
c
c integration of fct over r using simpson integration
c r(i)=r0*exp(dx*(i-1))
c
      implicit real*8 (a-h,o-z)
c
      dimension fct(1),r(1)
      real*8 ll
c
      data tiny/1.0d-20/
c
      if(ll+1.0.le.tiny)then
        write(9,*) ' sintg: ll+1=',ll+1
        sum=0.
      else
        rll=r(1)**ll
        if(rll.le.tiny)then
          sum=0.
        else
          fact=fct(1)/rll
          sum=fact*(r(1)**(ll+1))/(ll+1)
        endif
      endif
c
      x0=fct(1)*r(1)
      do 100 i=3,n,2
      x1=fct(i-1)*r(i-1)
      x2=fct(i)*r(i)
      sum=sum+dx*(x0+4.d0*x1+x2)
      x0=x2
  100 continue
      sintg=sum/3.d0
      if(mod(n,2).eq.0)then
        sintg=sintg+(fct(n)+fct(n-1))/2.d0*(r(n)-r(n-1))
      endif
      return
      end
      subroutine dintg(ll,idir,fct,yint,r,dx,n,rws)
c =================================================
c
c integration of fct(i) to yield yint(i)
c r(i)=r0*exp(dx*(i-1))
c
      implicit real*8 (a-h,o-z)
c
      dimension fct(1),yint(1),r(1)
      real*8 ll
c
      data tiny/1.0d-20/
c
      if(idir.gt.0)then
c
      if(ll+1.0.le.tiny)then
        write(9,*) ' dintg: ll+1=',ll+1
        sum=0.
      else
        rll=r(1)**ll
        if(rll.le.tiny)then
          sum=0.
        else
          fact=fct(1)/rll
          sum=fact*(r(1)**(ll+1))/(ll+1)
        endif
      endif
      yint(1)=sum
c
      x0=fct(1)*r(1)
      do 100 i=2,n
      x1=fct(i)*r(i)
c trapezoidal rule:
      sum=sum+dx*(x0+x1)/2.d0
      yint(i)=sum
      x0=x1
  100 continue
c
      else
c
      yint(n)=0.
      sum=0.
      x1=fct(n)*r(n)
c
      do 500 i=n-1,1,-1
      x0=fct(i)*r(i)
      sum=sum+dx*(x0+x1)/2.d0
      yint(i)=sum
      if(i.eq.n-3)then
        call inter1(r(n-3),yint(n-3),4,1,rws,corr)
        sum=sum-corr
        do 490 j=n-3,n
  490   yint(j)=yint(j)-corr
      endif
      x1=x0
  500 continue
c
      endif
      return
      end
      subroutine inter1(r,p,n,id,rs,ps)
c =====================================
      implicit real*8 (a-h,o-z)
c
c     interpolate via lagrange
c
      dimension r(n),p(n)
      ps=0.d0
      do 1 j=1,n,id
      term=1.d0
      denom=1.d0
      do 2 i=1,n,id
      if(i.eq.j) go to 2
      denom=denom*(r(j)-r(i))
      term=term*(rs-r(i))
    2 continue
    1 ps=ps+term *p(j)/denom
      return
      end
      subroutine sbf(e,r,fb,fn)
c =============================
c
c spherical bessel and neuman functions for            e > 0
c modified spherical bessel and neuman functions for   e < 0
c
      implicit real*8 (a-h,o-z)
      include 'l.par'
c
      complex*16 cfb(0:lmax+1),cfn(0:lmax+1)
      complex*16 ecompl,eroot,x1,x2,sinx,cosx,cimmi
      dimension  fb(0:lmax+1),fn(0:lmax+1)
c
      data cimmi/(0.0d0,-1.0d0)/
c
      ecompl=dcmplx(e,0.d0)
      eroot=cdsqrt(ecompl)
      x1=eroot*r
      x2=x1*x1
      sinx=cdsin(x1)
      cosx=cdcos(x1)
c
      cfb(0)= sinx/x1
      cfb(1)= sinx/x2 - cosx/x1
      cfn(0)=-cosx/x1
      cfn(1)=-cosx/x2 - sinx/x1
c
      do 5 l=2,lmax+1
      tl=dfloat(2*l-1)
      cfb(l)=tl*cfb(l-1)/x1-cfb(l-2)
    5 cfn(l)=tl*cfn(l-1)/x1-cfn(l-2)
c
      if(e.lt.0.) then
      cfn(0)=cfn(0)*cimmi
      do 10 l=1,lmax+1
      cfb(l)=cfb(l)*cimmi**l
   10 cfn(l)=cfn(l)*cimmi**(l+1)
      endif
c
      do 15 l=0,lmax+1
      fb(l)=dreal(cfb(l))
   15 fn(l)=dreal(cfn(l))
c
      return
      end
      subroutine sbf1(e,r,fb,fn)
c =============================
c
c spherical bessel and neuman functions for            e > 0
c modified spherical bessel and neuman functions for   e < 0
c
      implicit real*8 (a-h,o-z)
      include 'l.par'
c
      dimension  fb(0:lmax+1),fn(0:lmax+1)
c
c
c     shyp(x)=0.5d0*(dexp(x)-dexp(-x))
c     chyp(x)=0.5d0*(dexp(x)+dexp(-x))
c
      if(e.gt.0.) then
c
      eroot=dsqrt(e)
      x1=eroot*r
      x2=x1*x1
      sinx=dsin(x1)
      cosx=dcos(x1)
c
      fb(0)= sinx/x1
      fb(1)= sinx/x2 - cosx/x1
      fn(0)=-cosx/x1
      fn(1)=-cosx/x2 - sinx/x1
c
      do 5 l=2,lmax+1
      tl=dfloat(2*l-1)
      fb(l)=tl*fb(l-1)/x1-fb(l-2)
    5 fn(l)=tl*fn(l-1)/x1-fn(l-2)
c
      else
c
      eroot=dsqrt(-e)
      x1=eroot*r
      x2=x1*x1
      sinx=dsinh(x1)
      cosx=dcosh(x1)
c
      fb(0)= sinx/x1
      fb(1)=-sinx/x2 + cosx/x1
      fn(0)= cosx/x1
      fn(1)=-cosx/x2 + sinx/x1
c
      do 10 l=2,lmax+1
      tl=dfloat(2*l-1)
      fb(l)=-tl*fb(l-1)/x1+fb(l-2)
   10 fn(l)=-tl*fn(l-1)/x1+fn(l-2)
c
      endif
c
      return
      end
      subroutine wig3j(w3j,pun)
c ===================================
c
      implicit real*8 (a-h,o-z)
c
      parameter (ndim=100)
      include 'l.par'
c
      dimension w3j(-lmax-1:lmax,-lmax-1:lmax,0:lmax2)
      dimension thrcof(ndim)
      real*8 nul,j1,j2
      integer pun
c
      data nul,small,half,halfm/0.d0,0.01d0,0.5d0,-0.5d0/
      data iwigtes/0/
c
c                                
c                                   j1  lambda  j2 
c    < kappa1 lambda kappa2 > =  (                  ) =
c                                  1/2    0    -1/2
c
c
c                                  lambda  j2   j1
c                             =  (                  )  
c                                    0    -1/2  1/2
c
c for kappa1 = -lmax-1,-lmax,...,lmax-1,lmax
c for kappa2 = -lmax-1,-lmax,...,lmax-1,lmax
c
c convention for storage :  w3j(kappa1,kappa2,lambda)
c
c
      do kap1=-lmax-1,lmax
        do kap2=-lmax-1,lmax
          do lam =0,lmax2
            w3j(kap1,kap2,lam)=nul
          end do
        end do
      end do
c
      do 10 kap1=-lmax-1,lmax
c
        if(kap1.eq.0) goto 10
        if(kap1.lt.0) then
           j1=-kap1-half
           l1=-kap1-1
        else
           j1= kap1-half
           l1=kap1
        endif
c
        do 15 kap2=-lmax-1,lmax
c
          if(kap2.eq.0) goto 15
          if(kap2.lt.0) then
             j2=-kap2-half
             l2=-kap2-1
          else
             j2= kap2-half
             l2=kap2
          endif
c
          call rec3jj(thrcof,j2,j1,halfm,half,xlmin,xlmax,xlmat,ndim,
     >                ier)
          if(ier.lt.0) goto 15
c
          lamin=idint(xlmin+small)
          lamax=idint(xlmax+small)
c
          do 20 lam=lamin,lamax
c
            i1=lam-lamin+1
            if(mod(lam-lamin,2).eq.1) goto 20
            w3j(kap1,kap2,lam)=thrcof(i1)
c
            if(iwigtes.eq.1) then
              write(pun,100) kap1,l1,j1,kap2,l2,j2,lam,result
            endif
c
   20     continue
c
   15   continue
   10 continue
c
  100 format(2i4,f5.1,5x,2i4,f5.1,5x,i4,10x,d17.10)
c
      return
      end
      subroutine wig6j(w6j,j1,pun)
c ================================
c
      implicit real*8 (a-h,o-z)
c
      parameter (ndim=100)
      include 'l.par'
c
      dimension w6j(0:lvmax,0:lvmax,0:lmax,0:lmax,0:lmax)
      dimension sixcof(ndim)
      integer pun
      real*8 nul,j,jp,j1,j2
c
      data nul,small,half/0.d0,0.01d0,0.5d0/
      data iwigtes/0/
c
c
c  wigner 6j coefficients
c
c            j'  lambda'  j2         lambda'  j1  j
c          [                 ]  =  [                 ]
c            j   lambda   j1         lambda   j2  j'
c
c multiplied by (-)**(lambda+lambda'+1)
c
c
c for given j1
c
c convention for storage : w6j(l,l',l2,lambda,lambda') ,
c
c  with  j =(2*l +1)/2        l =0,1,...,lvmax
c        j'=(2*l'+1)/2        l'=0,1,...,lvmax
c        j2=(2*l2+1)/2        l2=0,1,...,lmax
c
c
      do l=0,lvmax
        do lp=0,lvmax
          do l2=0,lmax
            do lam=0,lmax
              do lamp=0,lmax
                w6j(l,lp,l2,lam,lamp)=nul
              end do
            end do
          end do
        end do
      end do
c
c
      do l=0,lvmax
        j=(2.*l+1.)*half
        do lp=0,lvmax
          jp=(2.*lp+1.)*half
c
          lami1=idint(dabs(j1-jp)+small)
          lama1=idint(j1+jp+small)
c
          do l2=0,lmax
            j2=(2.*l2+1.)*half
c
            lami2=idint(dabs(j-j2)+small)
            lama2=idint(j+j2+small)
            lami=max(lami1,lami2)
            lama=min(lama1,lama2)
c
            do 20 lam=0,lmax
              if(lam.lt.lami.or.lam.gt.lama) goto 20
              xlam=dfloat(lam)
c
              call rec6j(sixcof,j1,j,xlam,j2,jp,
     >                   xlammi,xlamma,xlmat,ndim,ier)
              if(ier.lt.0) goto 20
c
              lammi=idint(xlammi+small)
              lamma=idint(xlamma+small)
c
              do lamp=lammi,lamma
                coeff=1.0-2.0*mod(lam+lamp+1,2)
                w6j(l,lp,l2,lam,lamp)=coeff*sixcof(lamp-lammi+1)
c
                if(iwigtes.eq.1) then
                  write(pun,100) lamp,j1,j,lam,j2,j3,
     >                           sixcof(lamp-lammi+1)
                endif
              end do
c
   20       continue
c
          end do
        end do
      end do
c
  100 format(i4,2f5.1,5x,i4,2f5.1,10x,d17.10)
c
      return
      end
      subroutine rec3jj(thrcof,l2,l3,m2,m3,l1min,l1max,lmatch,ndim,ier)
c =====================================================================
c
c  j1-recursion of 3j-coefficients
c recursive evaluation of 3j- and
c 6j-coefficients.  k. schulten, r.g. gordon.
c ref. in comp. phys. commun. 11 (1976) 269
c
c
      implicit real * 8  (a-h)
      implicit real * 8  (o-z)
      real * 8   l1, l2, l3, m1, m2, m3, l1min, l1max, newfac, lmatch
      dimension thrcof(ndim)
c
      data  zero,eps,half,one /0d0,0.01d0,0.5d0,1d0/
      data two, three /2d0,3d0/
c
c  routine to generate set of 3j-coefficients   l1  l2  l3
c                                               m1  m2  m3
c  by recursion from l1min = max0(/l2-l3/,/m1/)
c                 to l1max =       l2+l3
c
c  the resulting 3j-coefficients are stored as thrcof(l1-l1min+1)
c
c  for a discussion of the recursion equation used see k. schulten
c  and r.g. gordon, j. math. phys. 16, 1961-1970 (1975), ibid. 16,
c  1971-1988 (1975)
c  for the sake of numerical stability the recursion will proceed
c  simultaneously forwards and backwards, starting from l1min and
c  l1max, respectively.
c  lmatch is the l1-value at which forward and backward recursion are
c  matched.
c  ndim is the length of the array thrcof.
c
c  ier is set to -1 if l2-/m2/ or l3-/m3/ less than zero or not
c                   integer, in which case all 3j-coefficients
c                   vanish.
c  ier is set to -2 if number of possible 3j-coefficients exceeds ndim
c  ier is set to non-negative number otherwise (0 + times of rescaling)
c
c  tiny should be set close to the smallest positive floating point
c  number which is representable on the computer.  srtiny is square
c  root of tiny .
c
      data tiny, srtiny  /1.0 d-10,1.0 d-05/
c
c  huge should be set close to largest positive floating point
c  number which is representable on the computer.  srhuge is
c  square root of huge .
c
      data  huge, srhuge /1.0 d 10,1.0 d 05/
c
      lmatch = zero
      m1 = - m2 - m3
c
c  check relative magnitude of l- and m-values
      if(l2-dabs(m2)+eps) 5,1,1
    1 if(l3-dabs(m3)+eps) 5,2,2
    2 if(dmod(l2+dabs(m2)+eps,one).ge.eps+eps)   go to 5
      if(dmod(l3+dabs(m3)+eps,one).ge.eps+eps)   go to 5
c
c
c
c  limits for l1
c
      l1min = dmax1(dabs(l2-l3),dabs(m1))
      l1max = l2 + l3
      if(l1min.lt.l1max-eps)   go to 20
      if(l1min.lt.l1max+eps)   go to 10
c
c
c  this is reached if l2-/m2/ and l3-/m3/  less than zero or not integer
c
    5 ier = - 1
      write(9,51) l2, l3, m1, m2, m3
   51 format(///1x,15h3j-coefficients,9x,2hl1,2f7.1/20x,3f7.1,4x,
     1 57hdo not satisfy the condition l2-/m2/ and l3-/m3/ ge zero ,
     2 11hand integer)
      return
c
c
c  this is reached in case that l1 can take only one value,
c  i.e. l1min = l1max
c
   10 ier = 0
      thrcof(1) = (-one) ** idint(dabs(l2+m2-l3+m3)+eps) /
     1 dsqrt(l1min + l2 + l3 + one)
      l1cmin = l1min
      l1cmax = l1max
      return
c
c
c
   20 ier = 0
      nfin = idint(l1max-l1min+one+eps)
      if(ndim-nfin)  21, 23, 23
c
c  dimension of thrcof not large enough to hold all the coefficients
c  required
c
   21 ier = - 2
      write(9,22)  l2, l3, m1, m2, m3, nfin, ndim
   22 format(///1x,15h3j-coefficients,9x,2hl1,2f7.1/20x,3f7.1,4x,
     1 26hexceed storage provided  (,i4,1h,,i4,1h))
      return
c
c
c  starting forward recursion from l1min taking nstep1 steps
c
   23 l1 = l1min
      thrcof(1) = srtiny
      sum1 = (l1+l1+one) * tiny
c
c
      lstep = 1
   30 lstep = lstep + 1
      l1 = l1 + one
c
c
      oldfac = newfac
      a1 = (l1+l2+l3+one) * (l1-l2+l3) * (l1+l2-l3) * (-l1+l2+l3+one)
      a2 = (l1+m1) * (l1-m1)
      newfac = dsqrt(a1*a2)
      if(l1.lt.one+eps)   go to 40
c
c
      dv = - l2*(l2+one) * m1 + l3*(l3+one) * m1 + l1*(l1-one) * (m3-m2)
      denom = (l1-one) * newfac
c
      if(lstep-2)  32, 32, 31
c
   31 c1old = dabs(c1)
   32 c1 = - (l1+l1-one) * dv / denom
      go to 50
c
c  if l1 = 1  (l1-1) has to be factored out of dv, hence
c
   40 c1 = - (l1+l1-one) * l1 * (m3-m2) / newfac
c
   50 if(lstep.gt.2)   go to 60
c
c
c  if l1 = l1min + 1  the third term in the recursion equation vanishes
c  , hence
      x = srtiny * c1
      thrcof(2) = x
      sum1 = sum1 + tiny * (l1+l1+one) * c1*c1
      if(lstep.eq.nfin)   go to 220
      go to 30
c
c
   60 c2 = - l1 * oldfac / denom
c
c  recursion to the next 3j-coefficient x
c
      x = c1 * thrcof(lstep-1) + c2 * thrcof(lstep-2)
      thrcof(lstep) = x
      sumfor = sum1
      sum1 = sum1 + (l1+l1+one) * x*x
      if(lstep.eq.nfin)   go to 100
c
c  see if last unnormalized 3j-coefficient exceeds srhuge
c
      if(dabs(x).lt.srhuge)   go to 80
c
c  this is reached if last 3j-coefficient larger than srhuge
c  so that the recursion series thrcof(1), ... , thrcof(lstep)
c  has to be rescaled to prevent overflow
c
      ier = ier + 1
      do 70 i=1,lstep
      if(dabs(thrcof(i)).lt.srtiny)   thrcof(i) = zero
   70 thrcof(i) = thrcof(i) / srhuge
      sum1 = sum1 / huge
      sumfor = sumfor / huge
      x = x / srhuge
c
c  as long as /c1/ is decreasing the recursion proceeds towards
c  increasing 3j-values and, hence, is numerically stable.  once
c  an increase of /c1/ is detected the recursion direction is
c  reversed.
c
   80 if(c1old-dabs(c1))   100, 100, 30
c
c
c  keep three 3j-coefficients around lmatch for comparision with
c  backward recursion.
c
  100 lmatch = l1 - 1
      x1 = x
      x2 = thrcof(lstep-1)
      x3 = thrcof(lstep-2)
      nstep2 = nfin - lstep + 3
c
c
c
c
c  starting backward recursion from l1max taking nstep2 steps, so
c  that forward and backward recursion overlap at three points
c  l1 = lmatch+1, lmatch, lmatch-1.
c
      nfinp1 = nfin + 1
      nfinp2 = nfin + 2
      nfinp3 = nfin + 3
      l1 = l1max
      thrcof(nfin) = srtiny
      sum2 = tiny * (l1+l1+one)
c
      l1 = l1 + two
      lstep = 1
  110 lstep = lstep + 1
      l1 = l1 - one
c
      oldfac = newfac
      a1s = (l1+l2+l3)*(l1-l2+l3-one)*(l1+l2-l3-one)*(-l1+l2+l3+two)
      a2s = (l1+m1-one) * (l1-m1-one)
      newfac = dsqrt(a1s*a2s)
c
      dv = - l2*(l2+one) * m1 + l3*(l3+one) * m1 + l1*(l1-one) * (m3-m2)
c
      denom = l1 * newfac
      c1 = - (l1+l1-one) * dv / denom
      if(lstep.gt.2)   go to 120
c
c  if l1 = l1max + 1  the third term in the recursion formula vanishes
c
      y = srtiny * c1
      thrcof(nfin-1) = y
      sumbac = sum2
      sum2 = sum2 + tiny * (l1+l1-three) * c1*c1
c
      go to 110
c
c
  120 c2 = - (l1 - one) * oldfac / denom
c
c  recursion to the next 3j-coefficient y
c
      y = c1 * thrcof(nfinp2-lstep) + c2 * thrcof(nfinp3-lstep)
c
      if(lstep.eq.nstep2)   go to 200
c
      thrcof(nfinp1-lstep) = y
      sumbac = sum2
      sum2 = sum2 + (l1+l1-three) * y*y
c
c  see if last unnormalized 3j-coefficient exceeds srhuge
c
      if(dabs(y).lt.srhuge)   go to 110
c
c  this is reached if last 3j-coefficient larger than srhuge
c  so that the recursion series thrcof(nfin), ... ,thrcof(nfin-lstep+1)
c  has to be rescaled to prevent overflow
c
      ier = ier + 1
      do 130 i=1,lstep
      index = nfin - i + 1
      if(dabs(thrcof(index)).lt.srtiny)   thrcof(index) = zero
  130 thrcof(index) = thrcof(index) / srhuge
      sum2 = sum2 / huge
      sumbac = sumbac / huge
c
c
      go to 110
c
c
c  the forward recursion 3j-coefficients x1, x2, x3 are to be matched
c  with the corresponding backward recursion values y1, y2, y3.
c
  200 y3 = y
      y2 = thrcof(nfinp2-lstep)
      y1 = thrcof(nfinp3-lstep)
c
c
c  determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
c  with minimal error.
c
      ratio = ( x1*y1 + x2*y2 + x3*y3 ) / ( x1*x1 + x2*x2 + x3*x3 )
      nlim = nfin - nstep2 + 1
c
      if(dabs(ratio).lt.one)   go to 211
c
      do 210 n=1,nlim
  210 thrcof(n) = ratio * thrcof(n)
      sumuni = ratio * ratio * sumfor + sumbac
      go to 230
c
  211 nlim = nlim + 1
      ratio = one / ratio
      do 212 n=nlim,nfin
  212 thrcof(n) = ratio * thrcof(n)
      sumuni = sumfor + ratio*ratio*sumbac
      go to 230
c
  220 sumuni = sum1
c
c
c  normalize 3j-coefficients
c
  230 cnorm = one / dsqrt(sumuni)
c
c  sign convention for last 3j-coefficient determines overall phase
c
      sign1 = dsign(one,thrcof(nfin))
      sign2 = (-one) ** idint(dabs(l2+m2-l3+m3)+eps)
      if(sign1*sign2) 235,235,236
  235 cnorm = - cnorm
c
  236 if(dabs(cnorm).lt.one)   go to 250
c
      do 240 n=1,nfin
  240 thrcof(n) = cnorm * thrcof(n)
      return
c
  250 thresh = tiny / dabs(cnorm)
      do 251 n=1,nfin
      if(dabs(thrcof(n)).lt.thresh)   thrcof(n) = zero
  251 thrcof(n) = cnorm * thrcof(n)
c
      return
      end
      subroutine rec6j(sixcof,l2,l3,l4,l5,l6,l1min,l1max,lmatch,ndim,
     1                 ier)
c ===================================================================
c
c  j1-recursion of 6j-coefficients
c recursive evaluation of 3j- and
c 6j-coefficients.  k. schulten, r.g. gordon.
c ref. in comp. phys. commun. 11 (1976) 269
c
c
      implicit real * 8  (a-h)
      implicit real * 8  (o-z)
      real * 8 l1,l2,l3,l4,l5,l6,l1min,l1max,newfac,lmatch
      dimension sixcof(ndim)
c
      data  zero,eps,half,one /0d0,0.01d0,0.5d0,1d0/
      data  two,three /2d0,3d0/
c
c  routine to generate the set of 6j-coefficients    l1 l2 l3
c                                                    l4 l5 l6
c  by recursion from l1min = max0(/l2-l3/,/l5-l6/)
c                 to l1max = min0( l2+l3 , l5+l6 )
c
c  the resulting 6j-coefficients are stored as sixcof(l1-l1min+1).
c
c  for a discussion of the recursion equation used see k. schulten
c  and r.g. gordon, j. math. phys. 16, 1961-1970 (1975), ibid. 16,
c  1971-1988 (1975)
c  for the sake of numerical stability the recursion will proceed
c  simultaneously forward and backwards, starting at l1min and
c  l1max, respectively.
c  lmatch is the l1-value at which forward and backward recursion
c  are matched.  its value will be returned, though it is not of
c  actual use.
c
c  ndim is the length of the array sixcof.
c
c  ier is set to -1 if either no 6j-coefficient satisfies triangular
c                  condition or if l2+l3+l5+l6 or l2+l4+l6 not
c                  integer, in which case all 6j-coefficients
c                  vanish.
c  ier is set to -2 if number of possible 6j-coefficients exceeds ndim
c
c  tiny should be set close to the smallest positive floating point
c  number which is representable on the computer.  srtiny is square
c  root of tiny .
c
      data  tiny, srtiny  /1.0 d-10, 1.0 d-05/
c
c  huge should be set close to largest positive floating point
c  number which is representable on the computer.  srhuge is
c  square root of huge .
      data  huge, srhuge  /1.0 d 10, 1.0 d 05/
c
      lmatch = zero
c
c
c
c  check if 6j-coefficients obey selection rules
c
      if(dmod(l2+l3+l5+l6+eps,one).ge.eps+eps)   go to 7
      if(dmod(l4+l2+l6+eps,one).ge.eps+eps)   go to 7
c
      if( l4+l2-l6)  7,1,1
    1 if( l4-l2+l6)  7,2,2
    2 if(-l4+l2+l6)  7,3,3
c
    3 if( l4+l3-l5)  7,4,4
    4 if( l4-l3+l5)  7,5,5
    5 if(-l4+l3+l5)  7,6,6
c
c  limits for l1
c
    6 l1min = dmax1(dabs(l2-l3),dabs(l5-l6))
      l1max = dmin1(l2+l3,l5+l6)
      if(l1min.lt.l1max-eps)   go to 20
      if(l1min.lt.l1max+eps)   go to 10
c
c
c  this is reached if triangular condition not satisfied
c  or if l2+l3+l5+l6 or l2+l4+l6 not integer
c
    7 ier = - 1
      write(9,71) l2, l3, l4, l5, l6
   71 format(///1x,15h6j-coefficients,9x,2hl1,2f7.1,4x,
     1 39hdo not satisfy triangular conditions or/20x,3f7.1,4x,
     2 35hl2+l3+l5+l6 or l2+l4+l6 not integer)
      return
c
c
c  this is reached in case that l1 can take only one value
c
   10 ier = 0
      sixcof(1) = (-one) ** idint(l2+l3+l5+l6+eps) /
     1            dsqrt((l1min+l1min+one)*(l4+l4+one))
      l1cmin = l1min
      l1cmax = l1max
      return
c
c
   20 ier = 0
      nfin = idint(l1max-l1min+one+eps)
      if(ndim-nfin)   21, 23, 23
c
c  this is reached if array sixcof not large enough to hold all
c  6j - coefficients required
c
   21 ier = - 2
      write(9,22)  l2, l3, l4, l5, l6, nfin, ndim
   22 format(///1x,'6j-coefficients',9x,'l1',2f7.1/20x,3f7.1,4x,
     1 'exceed storage provided  (',i4,',',i4,')')
      return
c
c
c
c
c
c  start of forward recursion
c
c
   23 l1 = l1min
      sixcof(1) = srtiny
      sum1 = (l1+l1+one) * tiny
c
      lstep = 1
   30 lstep = lstep + 1
      l1 = l1 + one
c
      oldfac = newfac
      a1 = (l1+l2+l3+one) * (l1-l2+l3) * (l1+l2-l3) * (-l1+l2+l3+one)
      a2 = (l1+l5+l6+one) * (l1-l5+l6) * (l1+l5-l6) * (-l1+l5+l6+one)
      newfac = dsqrt(a1*a2)
c
      if(l1.lt.one+eps)   go to 40
c
      dv = two * ( l2*(l2+one)*l5*(l5+one) + l3*(l3+one)*l6*(l6+one)
     1           - l1*(l1-one)*l4*(l4+one) )
     2                   - (l2*(l2+one) + l3*(l3+one) - l1*(l1-one))
     3                   * (l5*(l5+one) + l6*(l6+one) - l1*(l1-one))
c
      denom = (l1-one) * newfac
c
      if(lstep-2)  32, 32, 31
c
   31 c1old = dabs(c1)
   32 c1 = - (l1+l1-one) * dv / denom
      go to 50
c
c  if l1 = 1   (l1 - 1) has to be factored out of dv, hence
c
   40 c1 = - two * ( l2*(l2+one) + l5*(l5+one) - l4*(l4+one) )
     1 / newfac
c
   50 if(lstep.gt.2)   go to 60
c
c  if l1 = l1min + 1 the third term in recursion equation vanishes
c
      x = srtiny * c1
      sixcof(2) = x
      sum1 = sum1 + tiny * (l1+l1+one) * c1 * c1
c
      if(lstep.eq.nfin)   go to 220
      go to 30
c
c
   60 c2 = - l1 * oldfac / denom
c
c  recursion to the next 6j - coefficient x
c
      x = c1 * sixcof(lstep-1) + c2 * sixcof(lstep-2)
      sixcof(lstep) = x
c
      sumfor = sum1
      sum1 = sum1 + (l1+l1+one) * x * x
      if(lstep.eq.nfin)   go to 100
c
c  see if last unnormalized 6j-coefficient exceeds srhuge
c
      if(dabs(x).lt.srhuge)   go to 80
c
c  this is reached if last 6j-coefficient larger than srhuge
c  so that the recursion series sixcof(1), ... ,sixcof(lstep)
c  has to be rescaled to prevent overflow
c
      ier = ier + 1
      do 70 i=1,lstep
      if(dabs(sixcof(i)).lt.srtiny)   sixcof(i) = zero
   70 sixcof(i) = sixcof(i) / srhuge
      sum1 = sum1 / huge
      sumfor = sumfor / huge
      x = x / srhuge
c
c
c  as long as the coefficient /c1/ is decreasing the recursion proceeds
c  towards increasing 6j-values and, hence, is numerically stable.
c  once an increase of /c1/ is detected, the recursion direction is
c  reversed.
c
   80 if(c1old-dabs(c1))   100, 100, 30
c
c
c  keep three 6j-coefficients around lmatch for comparision later
c  with backward recursion.
c
  100 lmatch = l1 - 1
      x1 = x
      x2 = sixcof(lstep-1)
      x3 = sixcof(lstep-2)
c
c
c
c  starting backward recursion from l1max taking nstep2 steps, so
c  that forward and backward recursion overlap at the three points
c  l1 = lmatch+1, lmatch, lmatch-1.
c
      nfinp1 = nfin + 1
      nfinp2 = nfin + 2
      nfinp3 = nfin + 3
      nstep2 = nfin - lstep + 3
      l1 = l1max
c
      sixcof(nfin) = srtiny
      sum2 = (l1+l1+one) * tiny
c
c
      l1 = l1 + two
      lstep = 1
  110 lstep = lstep + 1
      l1 = l1 - one
c
      oldfac = newfac
      a1s = (l1+l2+l3)*(l1-l2+l3-one)*(l1+l2-l3-one)*(-l1+l2+l3+two)
      a2s = (l1+l5+l6)*(l1-l5+l6-one)*(l1+l5-l6-one)*(-l1+l5+l6+two)
      newfac = dsqrt(a1s*a2s)
c
      dv = two * ( l2*(l2+one)*l5*(l5+one) + l3*(l3+one)*l6*(l6+one)
     1           - l1*(l1-one)*l4*(l4+one) )
     2                   - (l2*(l2+one) + l3*(l3+one) - l1*(l1-one))
     3                   * (l5*(l5+one) + l6*(l6+one) - l1*(l1-one))
c
      denom = l1 * newfac
      c1 = - (l1+l1-one) * dv / denom
      if(lstep.gt.2)   go to 120
c
c  if l1 = l1max + 1 the third term in the recursion equation vanishes
c
      y = srtiny * c1
      sixcof(nfin-1) = y
      if(lstep.eq.nstep2)   go to 200
      sumbac = sum2
      sum2 = sum2 + (l1+l1-three) * c1 * c1 * tiny
      go to 110
c
c
  120 c2 = - (l1-one) * oldfac / denom
c
c  recursion to the next 6j - coefficient y
c
      y = c1 * sixcof(nfinp2-lstep) + c2 * sixcof(nfinp3-lstep)
      if(lstep.eq.nstep2)   go to 200
      sixcof(nfinp1-lstep) = y
      sumbac = sum2
      sum2 = sum2 + (l1+l1-three) * y * y
c
c  see if last unnormalized 6j-coefficient exceeds srhuge
c
      if(dabs(y).lt.srhuge)   go to 110
c
c  this is reached if last 6j-coefficient larger than srhuge
c  so that the recursion series sixcof(nfin), ... ,sixcof(nfin-lstep+1)
c  has to be rescaled to prevent overflow
c
      ier = ier + 1
      do 130 i=1,lstep
      index = nfin-i+1
      if(dabs(sixcof(index)).lt.srtiny)   sixcof(index) = zero
  130 sixcof(index) = sixcof(index) / srhuge
      sumbac = sumbac / huge
      sum2 = sum2 / huge
c
      go to 110
c
c
c  the forward recursion 6j-coefficients x1, x2, x3 are to be matched
c  with the corresponding backward recursion values y1, y2, y3.
c
  200 y3 = y
      y2 = sixcof(nfinp2-lstep)
      y1 = sixcof(nfinp3-lstep)
c
c
c  determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
c  with minimal error.
c
      ratio = ( x1*y1 + x2*y2 + x3*y3 ) / ( x1*x1 + x2*x2 + x3*x3 )
      nlim = nfin - nstep2 + 1
c
      if(dabs(ratio).lt.one)   go to 211
c
      do 210 n=1,nlim
  210 sixcof(n) = ratio * sixcof(n)
      sumuni = ratio * ratio * sumfor + sumbac
      go to 230
c
  211 nlim = nlim + 1
      ratio = one / ratio
      do 212 n=nlim,nfin
  212 sixcof(n) = ratio * sixcof(n)
      sumuni = sumfor + ratio*ratio*sumbac
      go to 230
c
  220 sumuni = sum1
c
c
c  normalize 6j-coefficients
c
  230 cnorm = one / dsqrt((l4+l4+one)*sumuni)
c
c  sign convention for last 6j-coeff. determines overall phase
c
      sign1 = dsign(one,sixcof(nfin))
      sign2 = (-one) ** idint(l2+l3+l5+l6+eps)
      if(sign1*sign2) 235,235,236
  235 cnorm = - cnorm
c
  236 if(dabs(cnorm).lt.one)   go to 250
c
      do 240 n=1,nfin
  240 sixcof(n) = cnorm * sixcof(n)
      return
c
  250 thresh = tiny / dabs(cnorm)
      do 251 n=1,nfin
      if(dabs(sixcof(n)).lt.thresh)   sixcof(n) = zero
  251 sixcof(n) = cnorm * sixcof(n)
c
      return
      end
