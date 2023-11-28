      program rcore
c
      implicit real*8(a-h,o-z)
c
      parameter(noat=2)
      parameter(nrad=400)
      parameter(iorb=30)
c
      common/bla/den(iorb),dq1(iorb),dfl(iorb),nqn(iorb),nql(iorb),
     *           nk(iorb),nmax(iorb),nel(iorb),norb,icore
      common/dira/dv(nrad),dr(nrad),dp(nrad),dq(nrad),dpas,z,nstop,
     *            nes,tets,np,nuc
      common/ps2/titre(iorb),bar(10),test,testv
      common/shoot/dgc(nrad,iorb),dpc(nrad,iorb)
      common/ajf/icoul
      common/oo/ xx(nrad),yy(nrad)
      common/eee/ ra(noat,nrad)
      common/ee/ za(noat,nrad)
      common/xmesh/stval,delx
      common/dd/ nt(noat),nz(noat),nnk(noat),nrc(noat),rc(noat),
     *           vc(noat)
      common/atom/itot,iex,iscfat,nitpot,itatom
      common/ay/rhotot(nrad,noat),char(nrad,noat),rhoc(nrad,noat),
     *          rhoold(nrad,noat)
      common/short/iskip
      common/tom/thfpot(nrad)
      common/files/iunit8,iunt14
      common/kiva/iturn,itcore,jturn,itval,iemode
      common/ws/rws,nrws(noat)
      common/ws1/rws1,nrws1(noat)
      common/eshift/vcc(noat)
      common/ewald/vew(noat)
      common/sumec/icone
      common/etot/tcore(noat)
      dimension cwf(iorb)
      character aa*20, potin*20, potout*20,rcwf_out*20
c
c     ***************************************************************
c     * calculate fully relativistically one-electron orbitals and  *
c     * total energy for free atoms or core one-electron orbitals   *
c     * and core charge density                                     *
c     * calculate kinetic energy contribution from the core         *
c     ***************************************************************
c
      write(*,*) ' input?'
      read(*,80) aa
      open(unit=3,file=aa,status='old')
c
      read(3,80) potin
   80 format(a80)
      open(unit=8,file=potin,status='old')
      read(3,80) aa
      open(unit=6,file=potout,status='unknown')
      read(3,80) rcwf_out
      open(unit=7,file=rcwf_out,status='unknown')
c
      read(3,5) ipot,iskip,moat
c
      if(iunit8.eq.1) go to 9
c
      if(ipot.eq.0) call cpain
      if(ipot.eq.1) call apwin
      if(ipot.eq.2) call kkrin
c
    9 natom=moat
c
      call insld(natom)
c
      wa=-1.0d 00
      z=nz(natom)
      itatom=0
      if(iscfat.eq.0) norb=icore
c
    8 do 10 i=1,nrad
      rr=dexp(stval+(i-1)*delx)
      rhotot(i,natom)=0.0d 00
      rhoc(i,natom)=0.0d 00
      dr(i)=rr
      thfpot(i)=0.0d 00
      if(ipot.lt.0)   dv(i)=fpot(rr,z,wa)
   10 continue
c
c     add tail to the muffin-tin potential
c
      if(ipot.lt.0)   go to 25
c
      do 11 i=1,nrad
   11 dp(i)=rhoold(i,natom)
c
   12 if(iemode.ne.0) go to 13
      last=nrc(natom)
      rlast=rc(natom)
      qval=rsimp(dp,dr,rc(natom),nrc(natom))
      go to 15
c
   13 if(iemode.eq.2) go to 14
      last=nrws1(natom)
      rlast=rws1
      qval=rsimp3(dp,dr,rws1,nrws1(natom))
      go to 15
c
   14 last=nrws(natom)
      rlast=rws
      qval=rsimp3(dp,dr,rws,nrws(natom))
c
   15 if(icoul.eq.2) qval=0.0d 00
      if(icoul.eq.0) go to 18
c
c              ..  V(r)/2  -Vc/2               r < R
c     Vat(r) = |
c     [a.u.]   |
c              |   0                                   icoul = 0
c              |  -(Z-Q)[ 1/r + 1/R ]          r > R   icoul = 1
c              .. -Z[ 1/r + 1/R]                       icoul = 0
c
      do 17 i=1,nrad
   17 thfpot(i)= -(nz(natom)-qval)/dr(i)
c
   18 tail= 0.0d 00
c
      if(icoul.ne.0) tail=(nz(natom)-qval)/rlast
c
      last1=last+1
c
      do 19 i=1,last
      dv(i) = za(natom,i)/2.0d 00 -vc(natom)/2.d 00 -vcc(natom)/2.0d 00
   19 continue
c
c
      do 20 i=last1,nrad
      dv(i) = thfpot(i) + tail
   20 continue
c
      write(6,'( '' Potential for type'',i2)') natom
      write(6,'(4d20.10)') (dr(k),dv(k),k=1,nrad)
c
c     now do core part
c
   25 do 27 i=1,nrad
      xx(i)=0.0d 00
c
      do 27 j=1,norb
      dgc(i,j)=0.0d 00
      dpc(i,j)=0.0d 00
   27 continue
c
      imax0=0
      np=nrad
c
      do 100  j=1,norb
      tets=test
   30 do 31 kk=1,nrad
      dp(kk)=0.0d 00
   31 dq(kk)=0.0d 00
c
      call resld (nqn(j),nql(j),nk(j),imax,den(j),dfl(j),dq1(j),j,
     *            natom )
c
      if(imax.gt.imax0) imax0=imax
c
      if (nstop.eq.0) go to 40
      tets=test*10.0d 00
      go to 30
c
c     calculate (core) charge density
c
   40 dp(1)=0.0d 00
      dq(1)=0.0d 00
c
      do 50  i=1,imax
      dgc(i,j)=dp(i)
      dpc(i,j)=dq(i)
      xx(i)=xx(i)+nel(j)*(dp(i)*dp(i)+dq(i)*dq(i))
      if(ipot.lt.0) rhotot(i,natom)=xx(i)
      if(j.gt.icore) go to 50
      rhoc(i,natom)=xx(i)
   50 continue
c
  100 continue
c
c     solve atomic problem selfconsistently, if wanted
c
      np=imax0

      if(iscfat.eq.0) go to 300
c
c
      call poisat(natom)
      if(nstop.eq.0) go to 25
c
c     write out the precious results
c
  300 rinf=dr(np-1)
c
      if(iskip.eq.0) go to 425
      write(6,1)
      do 420 i=1,norb
  420 write(6,2) nqn(i),titre(i),den(i)
c
c     orthonormality relations
c
      write(6,3)
  425 do 450 i=1,norb
      do 450 j=i,norb
      if(nql(i).ne.nql(j)) go to 450
      if(nk(i).ne.nk(j)) go to 450
c
      do 440 k=1,np
      yy(k)=dpc(k,i)*dpc(k,j)+dgc(k,i)*dgc(k,j)
  440 continue
c
      if(iemode.eq.0) anorm=rsimp(yy,dr,rc(natom),nrc(natom))
      if(iemode.eq.1) anorm=rsimp3(yy,dr,rws1,nrws1(natom))
      if(iemode.eq.2) anorm=rsimp3(yy,dr,rws,nrws(natom))
c
      if(i.eq.j) cwf(i)=anorm
c
      if(iskip.eq.0) go to 450
c
      bnorm=rsimp(yy,dr,rinf,np)
      write(6,4) nqn(i),titre(i),nqn(j),titre(j),bnorm,anorm
c
  450 continue
c
c     calculate number of core electrons
c
      do 460 k=1,np
  460 xx(k)=rhoc(k,natom)
c
      if(iemode.eq.0) anorm=rsimp(xx,dr,rc(natom),nrc(natom))
      if(iemode.eq.1) anorm=rsimp3(xx,dr,rws1,nrws1(natom))
      if(iemode.eq.2) anorm=rsimp3(xx,dr,rws,nrws(natom))
      if(iturn.lt.itcore) go to 465
      if(jturn.lt.itval) go to 465
      write(6,6) anorm
c
  465 if(iscfat.eq.0) go to 500
c
c     calculate total energy for free atom
c
      call totec(natom)
      go to 1000
c
c     calculate core part of total energy
c
  500 eone=0.d 00
c
      do 700 i=1,icore
      anorm=1.0d 00
      if(icone.eq.1) anorm=cwf(i)
  700 eone = eone + anorm*nel(i)*den(i)
      eone = 2.0d 00*eone
c
      dq(1)=0.0d 00
      nkk1=nnk(natom)
      do 750 k=1,nkk1
      dp(k) = (za(natom,k)-vc(natom)-vcc(natom))*rhoc(k,natom)
  750 continue
c
      if(iemode.eq.0) ev = rsimp(dp,dr,rc(natom),nrc(natom))
      if(iemode.eq.1) ev = rsimp3(dp,dr,rws1,nrws1(natom))
c
      evcor=0.0d 00
      if(icoul.eq.0) go to 800
c
      do 760 i=last1,nrad
      if(rhoc(i,natom).eq.0.) go to 770
  760 continue
  770 iend=i-1
c
      do 780 i=last1,iend
  780 evcor = evcor + (dv(i-1)+dv(i))*(dr(i)-dr(i-1))
c
  800 tcore(natom) = eone - ev - evcor
c
      if(iskip.eq.0) go to 1000
      write(6,7) eone,ev,evcor,tcore(natom)
c
c
c  write out core charge densities
c
c      do j=1,nnk(natom)
c        write(7,6010) dr(j),rhoc(j,natom)
c      end do
        write(7,'(4d20.12)') (rhoc(j,natom),j=1,nnk(natom))
c
c
c  write out core wavefunctions
c
c
      do 6000 i=1,norb
      write(7,6008) nqn(i),titre(i)
      write(7,6009) den(i)
      write(7,5) nk(i)
      write(7,5) nnk(natom)
      do 6001 j=1,nnk(natom)
      write(7,6010) dr(j),dgc(j,i),dpc(j,i)
 6001 continue
 6000 continue
c
c
 1000 continue
c
      stop
c
    1 format(1h1,10x,'one electron energies'//)
    2 format(10x,i2,a4,5x,e15.8)
    3 format(10x//10x,'orthogonality relations'//)
    4 format(10x,'(',i2,a4,',',i2,a4,')',3x,2e15.8)
    5 format(5i4)
    6 format(10x//10x,'number of core-electrons in atomic sphere = ',
     *e15.8)
    7 format(15x,'e-one core ',e20.12/15x,'pot -term ',e20.12/
     *      15x,' <V-correction> ',e20.12/
     *      15x,' <e-c> ',e20.12/)
 6008 format(i1,a4)
 6009 format(e15.8)
 6010 format(3e20.10)
      end
      subroutine apwin
c
c     ******************************
c     * potential input apw type   *
c     ******************************
c
      implicit real*8(a-h,o-z)

      parameter(noat=2)
      parameter(nrad=400)
      parameter(irep=7)
c
      integer u
      complex*16 sctamp
      common/short/iskip
      common/dolly/sctamp(irep),conc
      common/samfac/confac
      common/dd/nt(noat),nz(noat),nk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/ee/za(noat,nrad)
      common/eee/ra(noat,nrad)
      common/oo/p(nrad),r(nrad)
      common/xmesh/stval,xdel
      common/lat/aa,conc1
      common/tit/title
      common/files/iunit8,iunt14
      common/ws/rws,nrws(noat)
      common/ws1/rws1,nrws1(noat)
c
      dimension title(20)
c
      iunit8=1
      imode=0
c
      read(8,1) title
      read(8,2) aa,conc,ilat
      read(8,2) stval,xdel
c
      if(stval.eq.0.) stval=-8.8d 00
      if(xdel.eq.0.) xdel=0.05d 00
c
      conc1=conc
      u=0
      if(iskip.eq.0) go to 50
      write(6,8)title
      write(6,9) aa,conc
c
c
   50 do 120 n=1,noat
      read(8,3) nz(n),nk(n),rc(n),vc(n)
      if(iskip.ne.0) write(6,10) nz(n),nk(n),rc(n),vc(n)
      kx=nk(n)
      read(8,4) (ra(n,i),za(n,i),i=1,kx)
  120 continue
c
      call fitrad(imode)
c
      do 190 n=1,noat
      k1=nrc(n)
  190 rc(n)=ra(n,k1)
c
      call mike(ilat)

c
      do 150 n=1,noat
      kmax=nk(n)
      do 122 k=1,kmax
      if(ra(n,k).gt.rws) go to 130
  122 continue
  130 nrws(n)=k-1
      nk(n)=k
      do 140 k=1,kmax
      if(ra(n,k).gt.rws1) go to 145
  140 continue
  145 nrws1(n)=k-1
  150 continue
c
      vbar=conc*vc(1)+(1-conc)*vc(2)
      vc(1)=vbar
      vc(2)=vbar
c
      if(iskip.eq.0) go to 350
c
      write(6,5)
      do 300 n=1,noat
      kmax=nk(n)
      write(6,6)
      write(6,7) (i,(ra(n,i+j),za(n,i+j),j=u,3),i=1,kmax,4)
  300 continue
  350 return
c
c
    1 format(20a4)
    2 format(2f10.5,i4)
    3 format(2i4,2f10.5)
    4 format(4e20.8)
    5 format(1h1,10x/10x///'-----new input'///)
    6 format(72x/12x,4('r',14x,'v',14x))
    7 format(1x,i3,8e15.7)
    8 format(1h1,20x,20a4//)
    9 format(10x,'a0 = ',f10.5,' conc = ',f10.5/)
   10 format(10x,'z = ',i4,' nk = ',i4,' rc = ',f10.5,' vc = ',f10.5/)
      end
      subroutine cpain
c
c     ****************************************************************
c     *       input routine for potentials in the cpa mode           *
c     ****************************************************************
c
      implicit real*8(a-h,o-z)
c
      parameter(irep=7)
      parameter(noat=2)
      parameter(nrad=400)
c
      complex*16 sctamp
c
      common/dd/nt(noat),nz(noat),nk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/ee/za(noat,nrad)
      common/eee/ra(noat,nrad)
      common/samfac/confac
      common/dolly/sctamp(irep),conc
      common/short/iskip
      common/lat/aa,conc1
      common/tit/title
      common/xmesh/stval,xdel
      common/files/iunit8,iunt14
      common/ws/rws,nrws(noat)
      common/ws1/rws1,nrws1(noat)
c
c      dimension nsymbl(noat)
      dimension title(20)
c
      iunit8=1
c
      read(8,1) title
      read(8,2) aa,conc,ilat
      read(8,2) stval,xdel
c
      if(stval.eq.0.) stval=-8.8d 00
      if(xdel.eq.0.) xdel=0.05d 00
c
      conc1=conc
      do 100 n=1,noat
      nt(n)=n
      read(8,3) nz(n),nk(n),nrc(n),rc(n),vc(n)
      kmax=nk(n)
      read(8,4) (za(n,i),i=1,kmax)
c
      do 50 i=1,kmax
      ra(n,i)=dexp(stval+ xdel*(i-1))
      za(n,i)=za(n,i)/ra(n,i)
   50 continue
  100 continue
c
      call mike(ilat)
c
      do 150 n=1,noat
      kmax=nk(n)
      do 120 k=1,kmax
      if(ra(n,k).gt.rws) go to 130
  120 continue
  130 nrws(n)=k-1
      nk(n)=k
      do 140 k=1,kmax
      if(ra(n,k).gt.rws1) go to 145
  140 continue
  145 nrws1(n)=k-1
  150 continue
c
      vbar=conc*vc(1)+(1-conc)*vc(2)
c
      vc(1)=vbar
      vc(2)=vbar
c
      if(iskip.lt.2) go to 300
c
      write(6,5) title
      ii=0
      write(6,6) aa,conc,vc(1)
      write(6,10) confac
      ii=0
      do 200 n=1,noat
      write(6,7) n
      write(6,8)
      kmax=nk(n)
      write(6,9) (i,(ra(n,i+j),za(n,i+j),j=ii,3),i=1,kmax,4)
  200 continue
  300 return
c
    1 format(20a4)
    2 format(2f10.5,i4)
    3 format(4x,3i4,2f11.7)
    4 format(5d15.8)
    5 format(1h1,10x, 20a4//)
    6 format(10x,'lattice constant ',f10.5/
     '10x,'concentration ',f10.5/
     *10x,'muffin tin zero ',f10.5//)
    7 format(1h1,10x,'potential for scatterer ',i4//)
    8 format(72x/12x,4('r',14x,'v',14x))
    9 format(1x,i3,8e15.7)
   10 format(10x,'conversion factor to rydberg ',f10.5/)
      end

      function fpot(r,z,wa)
c
c     *************************************************
c     *  thomas fermi potential at the radial point r *
c     *  z          -  atomic number                  *
c     *  wa         -  ionicity - 1                   *
c     *************************************************
c
      implicit real*8(a-h,o-z)
      wc=dsqrt((r*(z+wa)**(1./3.))/0.8853)
      wd=wc*(0.60112*wc+1.81061)+1.
      we=wc*(wc*(wc*(wc*(0.04793*wc+0.21465)+0.77112)+1.39515)+1.81061)+
     11.
      wc=(z+wa)*(wd/we)**2-wa
      fpot=-wc/r
      return
      end
      subroutine inouh (dp,dq,dr,dq1,dfl,dv,z,test,nuc,nstop,jc)
c
c     ****************************************
c     * start values for outward integration *
c     ****************************************
c
      implicit real*8(a-h,o-z)
      parameter(iorb=30)
      parameter(nrad=400)
c
      common/ps1/dep(10),deq(10),dd,dvc,dsal,dk,dm
      common/trois/ dpno(4,iorb),dqno(4,iorb)
      dimension dp(nrad),dq(nrad),dr(nrad)
c
      do 3 i=1,10
      dp(i)=0.0d 00
 3    dq(i)=0.0d 00
c
      if (nuc) 5,5,21
c
 5    dval=z/dvc
      deva1=-dval
      deva2=dv/dvc+dval/dr(1)-dd
      deva3=0.0d 00
c
      if (dk) 9,9,11
c
 9    dbe=(dk-dfl)/dval
      go to 15
 11   dbe=dval/(dk+dfl)
c
 15   dq(10)=dq1
      dp(10)=dbe*dq1
      go to 39
c
 21   dval=dv+z*(3.0d 00-dr(1)*dr(1)/(dr(nuc)*dr(nuc)))
     *           /(dr(nuc)+dr(nuc))
      deva1=0.0d 00
      deva2=(dval-3.0d 00*z/(dr(nuc)+dr(nuc)))/dvc-dd
      deva3=z/(dr(nuc)*dr(nuc)*dr(nuc)*dsal)
c
      if (dk) 33,33,35
c
 33   dp(10)=dq1
      go to 39
c
 35   dq(10)=dq1
c
 39   do 40 i=1,5
      dp(i)=dp(10)
      dq(i)=dq(10)
      dep(i)=dp(i)*dfl
 40   deq(i)=dq(i)*dfl
c
      m=1
 41   dm=m+dfl
      dsum=dm*dm-dk*dk+deva1*deva1
      dqr=(dsal-deva2)*dq(m+9)-deva3*dq(m+7)
      dpr=deva2*dp(m+9)+deva3*dp(m+7)
      dval=((dm-dk)*dqr-deva1*dpr)/dsum
      dsum=((dm+dk)*dpr+deva1*dqr)/dsum
c
      j=-1
c
      do 44 i=1,5
      dpr=dr(i)**m
      dqr=dsum*dpr
      dpr=dval*dpr
c
      if (m.eq.1) go to 43
      if (dabs(dpr/dp(i)).le.test.and.abs(dqr/dq(i)).le.test) j=1
c
 43   dp(i)=dp(i)+dpr
      dq(i)=dq(i)+dqr
      dep(i)=dep(i)+dpr*dm
 44   deq(i)=deq(i)+dqr*dm
c
      if (j.eq.1) go to 99
      dp(m+10)=dval
      dq(m+10)=dsum
      m=m+1
c
      if (m.le.20) go to 41
c
      nstop=45
 99   do 101 i=1,4
      dpno(i,jc)=dp(i+9)
 101  dqno(i,jc)=dq(i+9)
 999  return
      end
      subroutine insld(natom)
c
c     ************************************
c     * input subroutine for core-module *
c     ************************************
c
      implicit real*8(a-h,o-z)
      parameter(noat=2)
      parameter(nrad=400)
      parameter(iorb=30)
c
c     z        -  atomic number
c     norb     -  number of orbitals
c     icore    -  number of core orbitals
c     nitoe    -  number of iterations to adjust one electron
c                 energies (default=15)
c     itot     -  : 0 no total energy for atom and core
c                 : 1 total energy for atom and core
c     iex      -  : 0 gunnarsson-lunqvist exchange-correlation
c              -  : 1 kohn-sham exchange
c     iscfat   -  : 0 atomic problem is n o t solved selfconsistently
c                     during rkkr-cpa-scf iterations
c                 : 1 atomic problem is solved selfconsistently
c                     to start rkkr-cpa-scf
c     nitpot   -  number of iterations for the atomic potential,
c                 (default=100)
c     dpas     -  logarithmic increment
c     test     -  precission for one electron energies
c     den      -  one electron energy [a.u.]
c     nqn      -  principal quantum number
c     nk       -  kapa
c     nel      -  orbital occupation
c
      common/bla/den(iorb),dq1(iorb),dfl(iorb),nqn(iorb),nql(iorb),
     *           nk(iorb),nmax(iorb),nel(iorb),norb,icore
      common/dira/dv(nrad),dr(nrad),dp(nrad),dq(nrad),dpas,z,
     *            nstop,nes,tets,np,nuc
      common/ps2/ titre(iorb),bar(10),test,testv
      common/atom/itot,iex,iscfat,nitpot,itatom
      common/xmesh/stval,delx
      common/ee/za(noat,nrad)
      common/dd/nt(noat),nz(noat),nkk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/short/iskip
      common/ajf/icoul
      dimension oplus(4),ominus(3)
c
      data oplus/4hs1/2,4hp3/2,4hd5/2,4hf7/2/
      data ominus/4hp1/2,4hd3/2,4hf5/2/
c
c
      dpas = delx
      nes=150
      testv=1.e-07
      test=1.e-10
      nuc=0
      dval=0.0d 00
c
      read(3,1000) (bar(i),i=1,10)
      read(3,1001) norb,icore,nitoe,itot,iex,iscfat,nitpot,icoul
c
      if(nitoe.eq.0) nes=nitoe
      if(nitpot.eq.0) nitpot=100
      z=nz(natom)
      np= nrad
      itatom=0
c
      if(iskip.eq.0) go to 10
      if(icoul.eq.0) write(6,1007) natom
      if(icoul.eq.1) write(6,1008) natom
      if(icoul.eq.2) write(6,1009) natom
      write(6,1002) natom
      write(6,1005) norb,icore
c
c     read orbital information
c
   10 do 50 i=1,norb
c
      read(3,1003) den(i),nqn(i),nk(i),nel(i)
c
      nql(i)=iabs(nk(i))
      if (nk(i).lt.0) nql(i)=nql(i)-1
      dfl(i)=nk(i)*nk(i)
      dfl(i)=dsqrt(dfl(i)-dval)
      if(nk(i).lt.0) go to 25
c
c     j  = l - 1/2
c
      do 20 ll=1,3
      if(ll.ne.nk(i)) go to 20
      titre(i)=ominus(ll)
      go to 40
   20 continue
c
c     j = l + 1/2
c
   25 do 30 ll=1,4
      lp1=-ll
      if(lp1.ne.nk(i)) go to 30
      titre(i)=oplus(ll)
      go to 40
   30 continue
   40 if(iskip.ne.0) write (6,1004) nqn(i),nk(i),titre(i),nel(i),den(i)
   50 continue
c
      ielec=0
      do 100 i=1,norb
      ielec=ielec+nel(i)
      nmax(i)=np
      l=1
      j=nqn(i)-nql(i)
      if((j-2*(j/2)).eq.0) l=-l
      dq1(i)=l*nk(i)/iabs(nk(i))
  100 continue
      if(iscfat.eq.0) go to 200
      if(ielec.ne.nz(natom)) go to 300
c
  200 return
  300 write(6,1006) nz(natom),ielec
      stop 'input screwed up'
c
 1000 format(20a4)
 1001 format(8i4)
 1002 format(1h1,' atomic calculation for scatterer',i4//)
 1003 format(e15.8,3i4)
 1004 format(2i6,5x,a4,i5,e15.8)
 1005 format(10x,'number of orbitals:',i4/
     *       10x,'number of core orbitals:',i4//)
 1006 format(10x,'get the input right, dummy!'/
     *       10x,' z = ',i4,' number of electrons ',i4//)
 1007 format(10x,'atom ',i4,' no tail'/)
 1008 format(10x,'atom ',i4,' screened coulomb tail'/)
 1009 format(10x,'atom ',i4,' unscreened coulomb tail'/)
      end
      subroutine inth (dp,dq,dv,dr)
c
c     *****************************************************
c     * integration of the dirac equation using a 5 point *
c     * adams method                                      *
c     *****************************************************
c
      implicit real*8(a-h,o-z)
      common/ps1/dep(10),deq(10),db,dvc,dsal,dk,dm
c     dpas    -    delx, logarithmic increment for radius
c     dm      -    dpas/720
c     dkoef1  -    475/502
c     dkoef2  -    27/502
c
      dkoef1 = 475.0d 00/502.0d 00
      dkoef2 = 27.0d 00/502.0d 00
c
      dpr=dp+dm*((251.0d 00*dep(1)+2616.0d 00*dep(3)+1901.0d 00*dep(5))
     *    -(1274.0d 00*dep(2)+2774.0d 00*dep(4)))
      dqr=dq+dm*((251.0d 00*deq(1)+2616.0d 00*deq(3)+1901.0d 00*deq(5))
     *    -(1274.0d 00*deq(2)+2774.0d 00*deq(4)))
c
      do 13 i=2,5
      dep(i-1)=dep(i)
   13 deq(i-1)=deq(i)
c
      dsum=(db-dv/dvc)*dr
      dep(5)=-dk*dpr+(dsal*dr+dsum)*dqr
      deq(5)=dk*dqr-dsum*dpr
c
      dp=dp+dm*((106.0d 00*dep(2)+646.0d 00*dep(4)+251.0d 00*dep(5))
     *   -(19.0d 00*dep(1)+264.0d 00*dep(3)))
      dq=dq+dm*((106.0d 00*deq(2)+646.0d 00*deq(4)+251.0d 00*deq(5))
     *   -(19.0d 00*deq(1)+264.0d 00*deq(3)))
c
      dp=dkoef1*dp+dkoef2*dpr
      dq=dkoef1*dq+dkoef2*dqr
      dep(5)=-dk*dp+(dsal*dr+dsum)*dq
      deq(5)=dk*dq-dsum*dp
c
      return
      end
      subroutine kkrin
c
c     ******************************************************************
c     *                                                                *
c     * read in kkr-type potential input                               *
c     *                                                                *
c     ******************************************************************
c
      implicit real*8(a-h,o-z)
c
      parameter(noat=2)
      parameter(nrad=400)
      parameter(irep=7)
c
      complex*16 sctamp
      common/dolly/sctamp(irep),conc
      common/samfac/confac
      common/dd/nt(noat),nz(noat),nk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/ee/za(noat,nrad)
      common/eee/ra(noat,nrad)
      common/lat/aa,conc1
      common/xmesh/stval,xdel
      common/tit/title
      common/files/iunit8,iunt14
      common/ws/rws,nrws(noat)
      common/ws1/rws1,nrws1(noat)
      common/eshift/vcc(noat)
      common/kitty/avol,ec,ilat
c
      dimension title(20),rtemp(3),vtemp(3)
c
      write(6,*) ' enter kkrin!'
      iunit8=1
c
      read(8,10005) title
      read(8,*) aa,conc,ilat,iform,iemode,iold
      read(8,*) stval,xdel
c
      write(6,10005) title
      write(6,*) aa,conc,ilat,iform,iemode,iold
      write(6,*) stval,xdel
c
      if(stval.eq.0.) stval=-8.8d 00
      if(xdel.eq.0.) xdel=0.05d 00
c
      conc1=conc
c
      call mike(ilat)
c
      do 105 n=1,noat
      read(8,*) nt(n),nz(n),nk(n),nrc(n),rc(n),vc(n),vcc(n)
      write(6,*) nt(n),nz(n),nk(n),nrc(n),rc(n),vc(n),vcc(n)
      nnk=nk(n)
      nnc=nrc(n)
c
c     generate radius and read in potential
c
      if(iform.eq.1) read(8,10015) (za(n,k),k=1,nnc)
      if(iform.eq.0) read(8,10016) (za(n,k),k=1,nnc)
      if(iform.eq.2) read(8,*) (za(n,k),k=1,nnc)
c
      xrc=dlog(rc(n))
      x0=xrc-(nnc-1)*xdel
      stval=x0
      if(iold.gt.0) go to 102
c
      if(iemode.eq.0) xdel=(-stval+dlog(rws1))/(nnk-6)
c
  102 do 103 k=1,nnk
      ra(n,k)=dexp(x0)
      za(n,k)=za(n,k)/ra(n,k)
  103 x0=x0+xdel
      if(iold.ge.2) goto 105
      do k=nnc+1,nnk
        do kk=1,3
          rtemp(kk)=ra(n,k-4+kk)
          vtemp(kk)=za(n,k-4+kk)
        end do
        call interp(rtemp,vtemp,3,ra(n,k),za(n,k),dps,.false.)
      end do
  105 continue
c
      vbar=conc*vc(1)+(1-conc)*vc(2)
      vc(1)=vbar
      vc(2)=vbar
c
      do 150 n=1,noat
      kmax=nk(n)
      do 120 k=1,kmax
      if(ra(n,k).gt.rws) go to 130
  120 continue
  130 nrws(n)=k-1
c
      if(iemode.eq.0) go to 147
c
      nk(n)=k
      do 140 k=1,kmax
      if(ra(n,k).gt.rws1) go to 145
  140 continue
  145 nrws1(n)=k    - 1
      go to 150
c
  147 nrws1(n)=nrc(n)
      rws1=rc(n)
c
  150 continue
c
      write(6,*) ' kkrin left!'
      return
c
10005 format(20a4)
10007 format(f11.6,f9.5,4i4)
10014 format(4i4,3f11.7)
10015 format(4d15.8)
10016 format(4d20.12)
      end
      subroutine poisat(natom)
c
c     ***************************************
c     *    calculation of atomic potential  *
c     ***************************************
c
      implicit real*8(a-h,o-z)
c
      parameter(noat=2)
      parameter(nrad=400)
      parameter(iorb=30)
c
      common/atom/itot,iex,iscfat,nitpot,itatom
      common/dira/dv(nrad),dr(nrad),dp(nrad),dq(nrad),dpas,z,nstop,
     *            nes,tets,np,nuc
      common/ps2/titre(iorb),bar(10),test,testv
      common/oo/d(nrad),vold(nrad)
      common/ay/rhotot(nrad,noat),char(nrad,noat),rhoc(nrad,noat),
     *          rhoold(nrad,noat)
      common/eee/ra(noat,nrad)
      common/ee/za(noat,nrad)
      common/dd/nt(noat),nz(noat),nk(noat),nrc(noat),rc(noat),vc(noat)
      common/short/iskip
      dimension s1(nrad),s2(nrad)
c
c     set functions for gunnarsson-lundquvist exchange-correlation
c     phys.rev.b13, 4274 (1976)
c
      beta(rs)=1.0d 00 + 0.0545d 00*rs
     $        *dlog(1.0d 00 + 11.4d 00/rs)
c
c
      pi=datan(1.0d 00)*4.0d 00
      pi4=4.0d 00*pi
      third=1.0d 00/3.0d 00
      twoth=2.0d 00*third
c
      xsl= -6.0d 00*( 3.0d 00/(8.0d 00*pi) )**third
c
c     coulomb potential with poisson equation
c
      do 4 i=1,np
      vold(i)=dv(i)
      s1(i)=0.0d 00
      s2(i)=0.0d 00
    4 continue
c
      do 3 i=3,np
      if(rhotot(i,natom).eq.0.0) go to 5
    3 continue
    5 np1=i-1
c
      do 16 i=1,np1
   16 s1(i)=rhotot(i,natom)/dr(i)
c
      mft=np1-1
      do 6 i=2,mft
      l=mft-i+1
      s2(l)=s2(l+1)+(s1(l)+s1(l+1))*(dr(l+1)-dr(l))
    6 continue
c
    9 s3=rhotot(1,natom)*dr(1)
      za(natom,1)=s2(1)+(s3-2.0d 00*z)/dr(1)
c
      do 10 i=2,np1
      s3=s3+(rhotot(i,natom)+rhotot(i-1,natom))*(dr(i)-dr(i-1))
      za(natom,i)=s2(i)+(s3-2.0d 00*z)/dr(i)
   10 continue
c
      np11=np1+1
      do 15 i=np11,np
   15 za(natom,i)=(s3-2.0d 00*z)/dr(i)
c
c     exchange-correlation potential
c
      do 40 i=2,np1
      rhos=rhotot(i,natom)/(dr(i)*dr(i)*pi4)
      rhox=rhos**third
      rs=(3.0d 00/(pi4*rhos))**third
      vexks=twoth*xsl*rhox
      vex= vexks*beta(rs)
      if(iex.ne.0) vex=vexks
      za(natom,i)= za(natom,i) + vex
   40 continue
c
      do 42 i=1,np
      dv(i)=za(natom,i)/2.0d 00
      dl0=-1.0d 00/dr(i)
      if(dv(i).gt.dl0) dv(i)=dl0
   42 continue
c
c     check potential
c
      nstop=1
      itatom=itatom+1
      if(itatom.gt.nitpot) go to 60
c
      del0=dabs( (dv(1)-vold(1))/dv(1) )
      do 45 i=2,np
      del1= dabs( (dv(i)-vold(i))/dv(i) )
      if(del1.le.del0) go to 45
      del0=del1
   45 continue
c
      alpha=0.75
c
      if(iskip.ne.0) write(6,1) itatom,del0
    1 format(10x,'iteration ',i4,' delta V-max ',e15.8)
c
      if(del0.gt.testv) nstop=0
      do 50 i=1,np
      dv(i)=alpha*vold(i)+(1.-alpha)*dv(i)
   50 continue
      if(nstop.eq.0) go to 60
      write(6,1) itatom,del0
   60 return
c
      end
      subroutine resld (nqn,nql,nk,imax,de,dfl,dq1,jc,natom)
c
c     *************************
c     * dirac equation [a.u.] *
c     *************************
c
      implicit real*8(a-h,o-z)
      parameter(noat=2)
      parameter(nrad=400)
      parameter(iorb=30)
      common/dira/dv(nrad),dr(nrad),dp(nrad),dq(nrad),dpas,z,
     *            nstop,nes,test,np,nuc
      common/ps1/dep(10),deq(10),db,dvc,dsal,dk,dm
      common/trois/ dpno(4,iorb),dqno(4,iorb)
      common/ee/za(noat,nrad)
      common/eee/ra(noat,nrad)
      common/dd/nt(noat),nz(noat),nnk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/atom/itot,iex,iscfat,nitpot,itatom
      common/short/iskip
c
c
      dkoef=1.0d 00/720.0d 00
c
      nstop=0
      dvc=137.0373
      dsal=dvc+dvc
c
      epriv=de
      imm=0
      ies=0
      dk=nk
      lll=(nql*(nql+1))/2
      nd=0
      nodes=nqn-nql
      if (lll.ne.0) go to 11
c
      elim=-z*z/(1.5d 00*nqn*nqn)
c
      go to 19
   11 elim=dv(1)+lll/(dr(1)*dr(1))
c
      do 15 i=2,np
      val=dv(i)+lll/(dr(i)*dr(i))
      if (val.le.elim) elim=val
   15 continue
c
      if(elim) 19,17,17
c
   17 nstop=17
      write(6,1000) nstop
 1000 format(5x,'nstop = ',i4,'  2*v+l*(l+1)/r**2 is positive'/)
      stop 'in resld'
c
   19 if(iskip.ge.2)
     >write(6,'('' de='',d20.10,''elim='',d20.10)') de,elim
      if(de.le.elim) de=elim*0.5d 00
   21 if (imm.eq.1) go to 35
c
      do 25 i=7,np,2
      imat=np+1-i
      if ((dv(imat)+lll/(dr(imat)*dr(imat))-de).le.0.) go to 26
   25 continue
   26 if (imat.gt.5) go to 35
c
      de=de*0.5d 00
      if(de.lt.-test.and.nd.le.nodes) go to 21
   28 nstop=28
      write(6,1001) nstop
 1001 format(5x,'nstop = ',i4,' 2*v+l*(l+1)/r**2-2*e is positive'/)
      return
c
c     get start values for outward integration
c
   35 if(iskip.ge.2) write(6,'('' de='',d20.10)') de
      db=de/dvc
      call inouh (dp,dq,dr,dq1,dfl,dv(1),z,test,nuc,nstop,jc)
      if (nstop) 36,47,36
   36 nstop=36
      write(6,1002) nstop
 1002 format(5x,'nstop = ',i4,
     *' dexpansion at the origin does not converge'/)
      return
c
c     calculate number of nodes for the large component
c
   47 nd=1
c
      do 51 i=1,5
      dval=dr(i)**dfl
      if (i.eq.1) go to 50
      if (dp(i-1).eq.0.0d 00) go to 50
      if ((dp(i)/dp(i-1)).gt.0.0d 00) go to 50
      nd=nd+1
   50 continue
c
      dp(i)=dp(i)*dval
      dq(i)=dq(i)*dval
      dep(i)=dep(i)*dval
   51 deq(i)=deq(i)*dval
      k=-1+2*(nodes-2*(nodes/2))
      if ((dp(1)*k).gt.0.) go to 54
   53 nstop=53
      write(6,1003) nstop
 1003 format(5x,'nstop ',i4,' dexpansion error at the origin'/)
      return
c
   54 if ((k*nk*dq(1)).lt.0.) go to 53
      dm=dpas*dkoef
c
c     do now outward integration
c
      do 195 i=6,imat
      dp(i)=dp(i-1)
      dq(i)=dq(i-1)
c
      call inth (dp(i),dq(i),dv(i),dr(i))
c
      if (dp(i-1).eq.0.0d 00) go to 195
      if ((dp(i)/dp(i-1)).gt.0.0d 00) go to 195
      nd=nd+1
      if(nd.gt.nodes) go to 209
  195 continue
c
      if (nd.eq.nodes) go to 240
      de=0.8*de
      if(de.lt.-test) go to 21
  206 nstop=206
      write(6,1004) nstop
 1004 format(5x,'nstop = ',i4,' number of nodes is too small'/)
      return
c
  209 de=1.2*de
      if(de.gt.elim) go to 21
  210 nstop=210
      write(6,1005) nstop
 1005 format(5x,'nstop = ',i4,' number of nodes is too large'/)
      return
c
c     start values for inward integration, use descleaux's
c     magic boundary: 300
c
  240 dqm=dq(imat)
      dpm=dp(imat)
      if (imm.eq.1) go to 258
c
      do 255 i=1,np,2
      imax=np+1-i
      if (((dv(imax)-de)*dr(imax)*dr(imax)).le.300.) go to 258
  255 continue
c
  258 dd=dsqrt(-de*(2.0d 00+db/dvc))
c
      dpq=-dd/(dsal+db)
      dm=-dm
      do 277 i=1,5
      j=imax+1-i
      dp(j)=dexp(-dd*dr(j))
      dep(i)=-dd*dp(j)*dr(j)
      dq(j)=dpq*dp(j)
  277 deq(i)=dpq*dep(i)
      m=imax-5
c
c     do now inward integration
c
      do 301 i=imat,m
      j=m+imat-i
      dp(j)=dp(j+1)
      dq(j)=dq(j+1)
  301 call inth (dp(j),dq(j),dv(j),dr(j))
c
c     check left and right large components
c
      dval=dpm/dp(imat)
      if (dval.gt.0.) go to 313
      nstop=312
      write(6,1006) nstop
 1006 format(5x,'nstop = ',i4,' sign error for the large component'/)
      return
c
  313 do 315 i=imat,imax
      dp(i)=dp(i)*dval
  315 dq(i)=dq(i)*dval
c
c     calculate the norm
c
      dsum=0.0d 00
      if(dp(1).eq.0.0) go to 332
      dsum=3.0d 00*dr(1)*(dp(1)**2+dq(1)**2)/(dpas*(dfl+dfl+1.))
c
  332 do 333 i=3,imax,2
  333 dsum=dsum+dr(i)*(dp(i)**2+dq(i)**2)+4.*dr(i-1)*(dp(i-1)**2+dq(i-1)
     1**2)+dr(i-2)*(dp(i-2)**2+dq(i-2)**2)
      dsum=dpas*(dsum+dr(imat)*(dqm*dqm-dq(imat)*dq(imat)))/3.0d 00
c
c     modify one-electron energy
c
      dbe=dp(imat)*(dqm-dq(imat))*dvc/dsum
      imm=0
      val=dabs(dbe/de)
      if (val.le.test) go to 365
  340 dval=de+dbe
c
      if(dval.lt.0.) go to 360
      dbe=dbe*0.5d 00
      val=val*0.5d 00
      if (val.gt.test) go to 340
  345 nstop=345
      write(6,1007) nstop
 1007 format(5x,'nstop = ',i4,' energy converged to zero'/)
      return
c
  360 de=dval
      if(dabs(de-epriv).lt.test) go to 365
      epriv=de
      if(iskip.eq.2) write(6,1009) jc,ies,test,epriv,de
 1009 format(' orbital = ',i3,'iteration ',i3,' test',d13.5,
     *       3x,'e0 = ',e15.8,3x,'e1 = ',e15.8)
      if (val.le.0.1) imm=1
      ies=ies+1
      if(ies.le.nes) go to 21
  362 nstop=362
      write(6,1008) nstop
 1008 format(5x,'nstop = ',i4,'number of iterations is too large'/)
      write(6,*) ies,nes
      stop
c      return
c
c     renormalize wavefunction
c
  365 dsum=dsqrt(dsum)
      dq1=dq1/dsum
c
      do 367 i=1,imax
      dp(i)=dp(i)/dsum
  367 dq(i)=dq(i)/dsum
c
      do 368 i=1,4
      dpno(i,jc)=dpno(i,jc)/dsum
  368 dqno(i,jc)=dqno(i,jc)/dsum
c
      if(imax.eq.np) go to 398
      j=imax+1
c
      do 371 i=j,np
      dp(i)=0.0d 00
  371 dq(i)=0.0d 00
c
  398 nstop=0
  399 return
      end
      subroutine totec(natom)
      implicit real*8(a-h,o-z)
c
      parameter(noat=2)
      parameter(iorb=30)
      parameter(nrad=400)
c
c     ***********************************************
c     *  total energy for free atoms                *
c     *        in  [ryd]   !!!                      *
c     ***********************************************
c
      common/oo/xx(nrad),yy(nrad)
      common/ee/za(noat,nrad)
      common/bla/den(iorb),dq1(iorb),dfl(iorb),nqn(iorb),nql(iorb),
     *           nk(iorb),nmax(iorb),nel(iorb),norb,icore
      common/eee/ra(noat,nrad)
      common/tit/title(20)
      common/ps2/titre(iorb),bar(10),test,testv
      common/dd/nt(noat),nz(noat),nkk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/ay/rhotot(nrad,noat),char(nrad,noat),rhoc(nrad,noat),
     *          rhoold(nrad,noat)
      common/etot/tcore(noat)
      common/atom/itot,iex,iscfat,nitpot,itatom
      common/short/iskip
      common/ws/rws,nrws(noat)
      common/ws1/rws1,nrws1(noat)
      common/kiva/iturn,itcore,jturn,itval,iemode
      common/sumec/icone
      common/eshift/vcc(noat)
      common/ewald/vew(noat)
c
      common/dira/dv(nrad),dr(nrad),dp(nrad),dq(nrad),dpas,z,nstop,
     *            nes,tets,np,nuc
c
c
c     set functions for local density functional:
c     gunnarson-lundqvist, phys.rev.b13, 4274 (1976)
c
      beta(rs) = 1.0d 00 + 0.0545d 00*rs*dlog(1.0d 00 + 11.4d 00/rs)
      eps(rs) = - 0.0666d 00*( (1.0d 00 + (rs/11.4d 00)**3)
     &        *dlog(1.0d 00 + 11.4d 00*rs) + 0.5d 00*rs/11.4d 00
     &        - (rs/11.4d 00)**2 - 1.0d 00/3.0d 00)
c
      pi=datan(1.0d 00)*4.0d 00
      pi4=4.0d 00*pi
      third=1.0d 00/3.0d 00
      xks=-4.0d 00*( 3.0d 00/(8.0d 00*pi) )**third
c
c     ****************************************************
c     *      atomic reference total energy               *
c     ****************************************************
c
      write(6,1) natom
      if(iex.ne.0) write(6,2)
      if(iex.eq.0) write(6,3)
c
      do 20 i=3,np
      if(rhotot(i,natom).eq.0.) go to 25
   20 continue
   25 inf1=i-1
      rinf1=dr(inf)
      inf=np
      rinf=dr(inf)
c
      eone=0.0d 00
      do 40 i=1,norb
   40 eone = eone+nel(i)*den(i)
      eone = eone*2.0d 00
c
      ecore=0.0d 00
      do 45 i=1,icore
   45 ecore = ecore+nel(i)*den(i)
      ecore = 2.0d 00*ecore
c
      do 50 i=1,np
      dq(i) = dv(i)*rhoc(i,natom)
   50 dp(i) = dv(i)*rhotot(i,natom)
c
      ucore = 2.0d 00*rsimp(dq,dr,rinf,inf)
      u1 = 2.0d 00*rsimp(dp,dr,rinf,inf)
      t = eone - u1
      tcore(natom) = ecore - ucore
c

      do 52 i=2,np
   52 dp(i) = rhotot(i,natom)/dr(i)
c
      u2 = -2.0d 00*nz(natom)*rsimp(dp,dr,rinf,inf)
c
c     exchange-correlation terms
c
      iprint=0
c
      do 60 i=2,inf
      rho=rhotot(i,natom)/(pi4*dr(i)*dr(i))
      rhox=rho**third
      rs = (3.0d 00/(pi4*rho))**third
c     if(dr(i).gt.rc(natom)) go to 53
      if(rs.lt.9.0) go to 51
      if(iprint.eq.0) write(6,7) dr(i),rs
    7 format(10x,'cut-off at =',f10.5,' rs = ',e13.5//)
      iprint=1
      go to 53
c
   51 if(iex.eq.0) go to 55
c
c     kohn-sham exchange
c
   53 vxc=xks*rhox
      epsxc=0.75d 00*vxc
      go to 57
c
c     gunnarsson-lundquvist exchange-correlation
c
   55 vxc=xks*rhox*beta(rs)
      epsxc=0.75d 00*xks*rhox + eps(rs)
c
   57 dp(i)=vxc*rhotot(i,natom)
      dq(i)=epsxc*rhotot(i,natom)
   60 continue
c
c
      envxc=rsimp(dp,dr,rinf1,inf1)
      enexc=rsimp(dq,dr,rinf1,inf1)
c
c
      u = (u1 + u2 - envxc)/2.0d 00
c
      eatom = t + u + enexc
c
      write(6,4) eatom,t,tcore(natom),u,enexc,envxc
c
      if(iskip.eq.0) go to 200
      write(6,5)
      do 80 i=1,norb
   80 write(6,6) nqn(i),titre(i),den(i)
  200 return
c
    1 format(10x///10x,'separated atom total energy for scatterer
     * [ryd]',i4/)
    2 format(10x,'kohn-sham exchange'//)
    3 format(10x,'hedin-lundqvist exchange-correlation'/)
    4 format(15x,'<E>       ',e20.10/
     *       15x,'<t>       ',e20.10/
     *       15x,'<t-core>  ',e20.10/
     *       15x,'<u>       ',e20.10/
     *       15x,'<E-xc>    ',e20.10/
     *       15x,'<Vxc>     ',e20.10//)
    5 format(10x,'one electron energies '///)
    6 format(15x,i2,a4,e20.8)
      end
      function rsimp(f,r,rn,irn)
      implicit real*8(a-h,o-z)
      parameter(nrad=400)
c
c     **************************************
c     *   radial integration via simpson   *
c     **************************************
c
      common/xmesh/stval,xdel
c
c     deltax=0.05 for r in logar. scale
c
      dimension f(nrad),r(nrad)
c
      dx=xdel
c
      isw=0
      rsimp=0.0d 00
      if(irn.le.2) return
      ieven=(irn/2)*2
      if(ieven.eq.irn) isw=1
      np=irn-isw
      s=f(1)*r(1)+f(np)*r(np)
      nl=np-1
      do 5 i=2,nl,2
    5 s=s+4.0d 00*f(i)*r(i)
      nl=nl-1
      if(nl.lt.3) goto 15
      do 10 i=3,nl,2
   10 s=s+2.0d 00*f(i)*r(i)
   15 s=s*dx/3.0d 00
      if(isw.eq.1) goto 30
      rsimp=s
      return
   30 rsimp=s+(f(irn)*r(irn)+f(irn-1)*r(irn-1))*0.5d 00*dx
      return
      end
      function rsimp3(f,r,rn,jrn)
      implicit real*8(a-h,o-z)
      parameter(nrad=400)
c
c     **************************************
c     *   radial integration via simpson   *
c     *   with interpolation               *
c     **************************************
c
      common/xmesh/stval,xdel
      dimension rad(7),x1(7)
c
c     deltax for r in logar. scale
c
      dimension f(nrad),r(nrad)
c
      dx=xdel
c
      do 50 iturn=1,7
      irn= jrn+iturn-6
      x1(iturn)=r(irn)
c
      isw=0
      rsimp3=0.0d 00
      if(irn.le.2) return
      if(irn/2*2.eq.irn) isw=1
      np=irn-isw
      s=f(1)*r(1)+f(np)*r(np)
      nl=np-1
      do 5 i=2,nl,2
    5 s=s+4.0d 00*f(i)*r(i)
      nl=nl-1
      if(nl.lt.3) go to 15
      do 10 i=3,nl,2
   10 s=s+2.0d 00*f(i)*r(i)
   15 s=s*dx/3.0d 00
      if(isw.eq.1) go to 17
      rad(iturn)=s
      go to 50
   17 rad(iturn)=s+(f(irn)*r(irn)+f(irn-1)*r(irn-1))*0.5d 00*dx
   50 continue
c
      call interp(x1,rad,7,rn,val,val1,.false.)
c
      rsimp3=val
      return
      end
      subroutine fitrad(np)
c
c     ******************************************************************
c     *                                                                *
c     * fits apw. radius to logarithmic scale                          *
c     *                                                                *
c     ******************************************************************
c
      implicit real*8(a-h,o-z)
c
      parameter(noat=2)
      parameter(nrad=400)
c
      common/dd/nt(noat),nz(noat),nk(noat),nrc(noat),rc(noat),
     *          vc(noat)
      common/ee/za(noat,nrad)
      common/eee/ra(noat,nrad)
      common/oo/p(nrad),r(nrad)
      common/xmesh/stval,xdel
      common/ws/rws,nrws(noat)
c
      delx=xdel
c
      do 250 n=1,noat
      kmax=nk(n)
      kmax1=kmax+1
      kmax3=kmax1-7
c
      r(1)=0.
      p(1)=-2.0d 00*nz(n)
c
      do 50 l=1,kmax
      p(l+1)=za(n,l)*ra(n,l)
   50 r(l+1)=ra(n,l)
c
      i=0
c
      x=stval-delx
   60 x=x+delx
      i=i+1
      t=dexp(x)
      if(i.gt.np) go to 250
c
      do 70 k=1,kmax1
      if(r(k).gt.t) go to 80
   70 continue
c      go to 150
      go to 105
   80 if(k-4) 90,90,100
   90 kk=1
      go to 110
  100 if(k.ge.kmax3) go to 105
      kk=k-3
      go to 110
  105 kk=kmax3
c
  110 call interp(r(kk),p(kk),7,t,za(n,i),dps,.false.)
c
c     if(t.gt.rws) za(n,i)=0.0d 00
      za(n,i)= za(n,i)/t
      go to 60
  250 continue
c
  500 return
      end
      subroutine mike(ilat)
      implicit real*8(a-h,o-z)
      parameter(noat=2)
      common/ws/rws,nrws(noat)
      common/ws1/rws1,nrws1(noat)
      common/kitty/avol,ec,lattyp
      common/lat/a0,conc
      common/samfac/confac
c
c     ec is the ewald constant as listed in:
c         j.f.janak, phys.rev.b9, 3985 (1974)
c     avol volume of unit cell
c
      third=1.0d 00/3.0d 00
      pi=datan(1.0d 00)*4.0 d 00
      tofp=3.0d 00/(4.0d 00*pi)
c
      lattyp=ilat
      if(lattyp.eq.1) go to 10
      if(lattyp.eq.2) go to 20
c
c     lattyp=0, simple cubic
c
      avol=a0**3.0d 00
      ec=3.1166857d 00
      rws1=a0/2.0d 00
      write(6,1)
      go to 30
c
c     lattyp=1, fcc
c
   10 avol=(a0**3.0d 00)/4.0d 00
      ec=4.8320664d 00
      rws1=a0/dsqrt(8.0d 00)
      write(6,2)
      go to 30
c
c     lattyp=2, bcc
c
   20 avol=(a0**3.0d 00)/4.0d 00
      ec=4.085521d 00
      rws1=a0*dsqrt(3.0d 00/16.0d 00)
      write(6,3)
c
c     set wigner-seitz-radius and conversion factor
c     for dimensionless units [d.u.]
c
   30 rws=(tofp*avol)**third
      confac=(2.0d 00*pi/a0)**2.0d 00
c
      return
    1 format(10x,'simple cubic lattice'/)
    2 format(10x,'fcc lattice'/)
    3 format(10x,'bcc lattice'/)
      end
      subroutine interp(r,p,n,rs,ps,dps,deriv)
c
c     ****************************
c     *                          *
c     * interpolate via lagrange *
c     *                          *
c     ****************************
c
      implicit real*8(a-h,o-z)
      logical deriv,nodriv
      dimension r(n),p(n)
      nodriv=.not.deriv
      ps=0.d 00
      dps=0.d 00
      do 1 j=1,n
      term=1.d 00
      denom=1.d 00
      dterm=0.d 00
      do 2 i=1,n
      if(i.eq.j) go to 2
      denom=denom*(r(j)-r(i))
      term=term*(rs-r(i))
      if(nodriv) go to 2
      dterm1=1.d 00
      do 3 k=1,n
      if(k.eq.j.or.k.eq.i) go to 3
      dterm1=dterm1*(rs-r(k))
    3 continue
      dterm=dterm+dterm1
    2 continue
      if(nodriv) go to 1
      dps=dps+dterm*p(j)/denom
    1 ps=ps+term *p(j)/denom
      return
      end
