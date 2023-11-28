      program rcvvaug
c
c **********************************************************
C CALCULATE RELATIVISTIC CORE-VALENCE-VALENCE AUGER SPECTRA
c **********************************************************
c
      implicit real*8 (a-h,o-z)
c
      parameter (nemax=500,nemax2=1000,n1m=300,nemf=15,nemf2=30)
      parameter (km=3, kmtot=28 )
c
      character name*50,title*78,namorb*10,filen*50
c
      dimension e1(n1m),e(nemax),efile(nemf),efile2(nemf2),
     $e2(nemax2),irange(nemf,nemf2)
      dimension z(-km-1:km),udos(nemax,-km-1:km),dos(nemax,-km-1:km)
      dimension crsfile(nemf,nemf2,-km-1:km,-km-1:km),
     $crs(nemax,-km-1:km,-km-1:km)
      dimension p(nemax2,kmtot),psum(kmtot)
      dimension pdiag(nemax2),ptot(nemax2)
      dimension iord(kmtot),ikap(kmtot),ikapp(kmtot)
      dimension y(nemax),phelp(n1m)
      dimension p1(n1m,4),p1tot(n1m),p1diag(n1m)
c
      data pi/3.1415926536/,efac/13.606/,tiny/1.0d-5/
c
c
c  open i/o files  *****************
c
      in=1
      jdos=2
      jmat=3
      jprt=7
      jspect=8
c
      open(unit=in,file='rcvvaug.in')
c
c dos file
c
      read(in,1) filen
      write(*,1) filen
      open(unit=jdos,file=filen)
c
c 
c matrix element file
c
      read(in,1) filen
      write(*,1) filen
      open(unit=jmat,file=filen)
c
c 
c print file
c
      read(in,1) filen
      write(*,1) filen
      open(unit=jprt,file=filen)
c
c 
c spectrum file
c
      read(in,1) filen
      write(*,1) filen
      open(unit=jspect,file=filen)
c
c ******************************
c
c
      read(in,2) title
      write(*,2) title
c
c spectrometer resolution (fwhm in eV)
c
      read(in,*) hws
      write(*,*) hws
c
c 
c valence band lifetime broadening (fwhm in eV)
c (DOS will be broadened before convolution, 
c  i.e. valence band holes are considered to be independent) 
c
      read(in,*) hwlev
      write(*,*) hwlev
      hwlr=hwlev/efac
c
c
c core hole lifetime broadening (fwhm in eV)
c
      read(in,*) hwc 
      write(*,*) hwc
c
c
c two-hole final state lifetime broadening (fwhm in eV)
c
      read(in,*) hwlf
      write(*,*) hwlf
c
c
c relativistic or nonrelativistic densities of states (1/0)
c
      read(in,*) krel
      write(*,*) krel
c      
c
c shift the unbroadened spectrum downwards
c
      read(in,*) u
      write(*,*) u
c      
c 
c ***********  read in matrixelements ***************
c
      read(jmat,*)
      read(jmat,1) name
      read(jmat,*)
      read(jmat,3) namorb
      read(jmat,*)
      read(jmat,*) ecore
      read(jmat,*)
      read(jmat,*) e0,de,nume
      read(jmat,*)
      read(jmat,*) lmax
      read(jmat,*)
      read(jmat,*) a0
      read(jmat,*)
      read(jmat,*)
c
      write(jprt,*) ' Matrixelements'
      write(jprt,1) name
      write(jprt,3) namorb 
      write(jprt,4) ((kap,kapp,kapp=-lmax-1,kap),kap=-lmax-1,-1),
     $((kap,kapp,kapp=-lmax-1,-1),(kap,kapp,kapp=1,kap),kap=1,lmax)
c
      write(jspect,2) title
      write(jspect,1) name
      write(jspect,10) hws,hwlev,hwc,hwlf
c
      conv=2.0*pi/a0
      conv=conv*conv
      if(krel.eq.0) conv=1.0
c
      nume2=2*nume-1
      necount=nume*nume
      do 100 i=1,nume
      efile(i)=e0+de*(i-1)
  100 continue
c
      et0=2.0*e0-ecore
      do 105 i=1,nume2
      efile2(i)=et0+de*(i-1)
      do 106 j=1,nume
      jp=i+1-j
      irange(j,i)=1
      if(jp.lt.1.or.jp.gt.nume) then
      irange(j,i)=0
      goto 106
      endif
c
      read(jmat,*) ii,jj,dum,eep,ee,
     $((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,kap),kap=-lmax-1,-1),
     $((crsfile(jj,ii,kap,kapp),kapp=-lmax-1,-1),
     $ (crsfile(jj,ii,kap,kapp),kapp=1,kap),kap=1,lmax)
c
      write(jprt,5) i,j,efile2(i),efile(jp),efile(j),
     $((crsfile(j,i,kap,kapp),kapp=-lmax-1,kap),kap=-lmax-1,-1),
     $((crsfile(j,i,kap,kapp),kapp=-lmax-1,-1),
     $ (crsfile(j,i,kap,kapp),kapp=1,kap),kap=1,lmax)
c
      if(ii.ne.i.or.jj.ne.j.or.dabs(eep-efile(jp)).gt.tiny.or.
     $   dabs(ee-efile(j)).gt.tiny) then
      write(6,*) i,j,efile(jp),efile(j)
      write(6,*) ii,jj,eep,ee
      stop ' shit happened'
      endif
c
  106 continue
  105 continue   
c
c
c
c ***************  read in densities of states  ***************
c
c
      read(jdos,1) name
      read(jdos,*) 
      read(jdos,*) ef
      ef=conv*ef
      read(jdos,*) lmax
      read(jdos,*) maxd
      read(jdos,*)
c
      ii=0
      do 110 i=1,maxd
c
      if(krel.eq.0) then
      read(jdos,*) energy,(z(l),l=0,lmax),bla
      if(bla.lt.tiny.and.ii.eq.0) goto 110
      ii=ii+1
      e(ii)=energy*conv
      udos(ii,-1)=z(0)
      do 115 l=1,lmax
      fl=dfloat(2*l+1)
      fl=fl+fl
      fl1=dfloat(2*l)
      fl2=dfloat(2*l+2)
      fl1=fl1/fl
      fl2=fl2/fl
      udos(ii,l)   =z(l)*fl1
      udos(ii,-l-1)=z(l)*fl2
  115 continue
      else
      read(jdos,*) energy,bla,z(-1),
     $            (z(l),z(-l-1),l=1,lmax)
      if(bla.lt.tiny.and.ii.eq.0) goto 110
      ii=ii+1
      e(ii)=energy*conv
      do 116 k=-lmax-1,lmax
      udos(ii,k)=z(k)
  116 continue
      endif
c
  110 continue
      ne=ii
c      
      write(jprt,*)
      write(jprt,*) ' Densities of states'
      do 120 i=1,ne
      write(jprt,6) e(i),udos(i,-1),
     $            (udos(i,l),udos(i,-l-1),l=1,lmax)
  120 continue    
c  
c
      nef=iright(e,ne,ef)
c  
c check energy ranges:
      if(e(1).lt.efile(1).or.e(nef).gt.efile(nume))then
        write(6,*) ' energy range of matrix elements must exceed',
     *  ' that of DOS!'
        stop
      endif
c
c
c  life time broadening of DOS **********************
c
      call dosbroad(ef,nef,e,udos,hwlr,dos,lmax)
c
      write(jprt,*)
      write(jprt,*) ' Densities of states after broadening'
      do 125 i=1,nef
      write(jprt,6) e(i),dos(i,-1),
     $            (dos(i,l),dos(i,-l-1),l=1,lmax)
  125 continue    
c
      ii=0
      do 127 k=-lmax-1,lmax
      if(k.eq.0) goto 127 
      do 128 kp=-lmax-1,k
      if(kp.eq.0) goto 128
      ii=ii+1 
      ikap(ii)=k
      ikapp(ii)=kp
  128 continue
  127 continue
      ktot=ii
      write(6,*) ' ktot=',ktot
c
c determine boundaries for energy loop ***************
c
      e1ryd=ecore
      e2min=2.0*e(1)-e1ryd
      e2max=2.0*e(nef)-e1ryd
c
c set up final energy scale e2 (double grid has been set for array e2)
      de=e(2)-e(1)
      de2=e(3)-e(1)
      ne2=iidint((e2max-e2min)/de2)+1
      ee=e2min-de2
      do 130 i=1,ne2
      ee=ee+de2
      e2(i)=ee
  130 continue
c
c index values for test output of matrix elements:
      ktest=ne2/6
c
c
c start loop over final energies *********************
c
      do 150 ie2=1,ne2
      write(6,*) ie2,(e2(ie2)+e1ryd)/2.
c
c valence band energy boundaries for this final energy:
      emax=dmin1(e1ryd+e2(ie2)-e(1),e(nef))
      emin=dmax1(e1ryd+e2(ie2)-e(nef),e(1))
      ilast=iright(e,nef,emax)
      ifst=leftind(e,nef,emin)
      nval=ilast-ifst+1
c
c ********************************************************************
c interpolate matrix elements to DOS energy grid for this final energy
c ********************************************************************
c
c documentation:
      idoc=0
      if(mod(ie2,ktest).eq.0)idoc=1
      if(ie2.eq.1.or.ie2.eq.2)idoc=1
      if(ie2.eq.ne2-1.or.ie2.eq.ne2)idoc=1
c
      call matinter(crsfile,efile,nume,efile2,nume2,
     *nemf,nemf2,lmax,irange,e2(ie2),emin,emax,e,nemax,crs,
     *ifst,ilast,idoc,jprt)
c
c *******************************************************************
c                   put together spectrum:
c *******************************************************************
c
c
      ii=0
      do 160 k=-lmax-1,lmax
      if(k.eq.0) goto 160
c
      do 165 kp=-lmax-1,k
      if(kp.eq.0) goto 165
c
      ii=ii+1
c
      do 170 i=ifst,ilast
      ip=ilast+ifst-i
  170 y(i)=crs(i,k,kp)*dos(i,k)*dos(ip,kp)
      p(ie2,ii)=simpson(y(ifst),nval,de)
c
  165 continue
  160 continue
c
c
c sum up diagonal contributions and total spectrum:
c
      pdiagg=0.
      ptott=0.
      do 180 ii=1,ktot
      if(ikap(ii).eq.ikapp(ii)) pdiagg=pdiagg+p(ie2,ii)
      ptott=ptott+p(ie2,ii)
  180 continue
      pdiag(ie2)=pdiagg
      ptot(ie2)=ptott
c
  150 continue
c
c
c end of loop over final energies **********************************
c
c
c determine relative magnitudes of partial spectra:
c
      psumtot=0.
      do 190 ie2=1,ne2
  190 psumtot=psumtot+ptot(ie2)
      do 200 ii=1,ktot
      psumma=0.
      do 201 ie2=1,ne2
  201 psumma=psumma+p(ie2,ii)
  200 psum(ii)=psumma/psumtot
c
c maximum of total unbroadened spectrum:
c
      pmax=0.
      do 205 i=1,ne2
  205 if(ptot(i).gt.pmax) pmax=ptot(i)
c
c ordering of partial spectra according to magnitude:
c
      call order(psum(1),ktot,iord)
c
c
c ********************************************************************
c BROADENING DUE TO LIFE TIME OF CORE HOLE AND SPECTROMETER RESOLUTION
c ********************************************************************
c (only the most important contributions are broadened)
c
      snorm=0.
      write(6,*) ' broad'
      call broad(e2(1),ptot,nemax2,ne2,iminb,
     *e1(1),p1tot,phelp,n1m,n1,snorm,smax,hws,hwlf,hwc,u)
      do 210 i=1,4
      write(6,*) ' broad'
  210 call broad(e2(1),p(1,iord(i)),nemax2,ne2,iminb,
     *e1(1),p1(1,i),phelp,n1m,n1,snorm,dummy,hws,hwlf,hwc,u)
      write(6,*) ' broad'
      call broad(e2(1),pdiag,nemax2,ne2,iminb,
     *e1(1),p1diag,phelp,n1m,n1,snorm,dummy,hws,hwlf,hwc,u)
c
c ********************************************************************
C                      OUTPUT OF SPECTRA
c ********************************************************************
c
      e2maxev=e2max*efac
c
c
c unbroadened spectra 
c
      write(jprt,*) ' unbroadened spectra' 
      write(jprt,7) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
      write(jprt,8) (e2(i)-e2max,ptot(i),
     $(p(i,iord(ii)),ii=1,4),pdiag(i),i=1,ne2)
c
c
c broadened spectra
c
      write(jspect,3) namorb
      write(jspect,6) ecore*efac
      write(jspect,7) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
      write(jspect,8) (e1(i)-e2maxev,p1tot(i)/snorm,
     $(p1(i,ii)/snorm,ii=1,4),p1diag(i)/snorm,i=iminb,n1)
c
      write(jspect,9) (ikap(iord(i)),ikapp(iord(i)),i=1,4)
      write(jspect,6) (e1(i)-e2maxev,p1tot(i),
     $(p1(i,ii),ii=1,4),p1diag(i),i=iminb,n1)
c
c
    1 format(a50)
    2 format(a78)
    3 format(a10)
    4 format(/1x,'i2',2x,'i',3x,'e2',4x,'ep',4x,'e',2x,
     $28(i6,i3,'    ')/)
    5 format(2i3,3f6.2,28d13.5)    
    6 format(f11.5,6f9.3)
    7 format(4x,'e2',5x,5x,'Total',3x,4(5x,i2,1x,i2,3x),5x,'Diag.')
    8 format(f11.5,6e13.5)
    9 format(4x,'e2',5x,4x,'Total',4(4x,i2,1x,i2),4x,'Diag.')
   10 format('   HWS=',f5.2,' eV','    HWL=',f5.2,' eV',
     <       '   HWC=',f5.2,' eV','   HWLF=',f5.2,' eV')
c
c
      stop' FORTRAN stop'
      end
      subroutine dosbroad(ef,nef,e,udos,hwlr,dos,lmax)
c ====================================================
      implicit real*8 (a-h,o-z)
c
c life-time broadening of valence band states with Lorentzian
c only states below Efermi enter the integral
c DOS is required on an equidistant energy mesh!
c
      parameter (nemax=500)
      parameter (km=3)
c
      dimension e(nemax),udos(nemax,-km-1:km),dos(nemax,-km-1:km)
      dimension a(nemax)
c
      data pi/3.141596/,ilor/2/,hwlr0/0.1/
c
      write(6,*)' I am in dosbroad'
      if(hwlr.lt.hwlr0) then
      do 200 k=-lmax-1,lmax
      if(k.eq.0) goto 200
      do 201 j=1,nef 
      dos(j,k)=udos(j,k)
  201 continue
  200 continue
      write(6,*)' I am out of dosbroad'
      return
      endif
c
      estep=e(2)-e(1)
      do 100 k=-lmax-1,lmax
      if(k.eq.0) goto 100
      do 90 i=1,nef
      do 80 j=1,nef
      a(j)=udos(j,k)
   80 continue
      if(ilor.eq.1)then
c halfwidth of Lor. increases linearly with binding E.:
      gam=hwlr*(ef-e(i))/(ef-e(1))/2.0E0
      endif
      if(ilor.eq.2)then
c halfwidth of Lor. increases quadr.:
      gam=hwlr*((ef-e(i))/(ef-e(1)))**2/2.0E0
      endif
      dos(i,k)=convlo(gam,a(1),i,estep,1,nef,pi)
   90 continue
  100 continue
c 
      write(6,*)' Out of dosbroad'
      return
      end
      subroutine broad(ein,a0,n0,max0,mine1,
     ,                 e1,a1,ai,n1m,n1,snorm,smax,hws,hwl,hwc,u)
c ==============================================================
c
c     interpolate to plot-grid and convolute with lorentzian
c     for life-time broadening of core hole + valence band holes
c     and gaussian spectrometer resolution
c
      implicit real*8 (a-h,o-z)
      parameter (ngc=100,nmax=500)
      dimension ein(n0),a0(n0),e0(nmax),e1(n1m),a1(n1m),ai(n1m)
      dimension gc(-ngc:ngc),ah(-ngc:ngc)
      data hwc0/5.e-03/,cfac/13.606/
      data estep /0.05e0/
      data ix0/201/,pi/3.141596/
c
c     convert from rydberg to ev:
c
      write(6,*)' I am in broad'
      do i=1,max0
        e0(i)=ein(i)*cfac-u
      end do
      de0=e0(2)-e0(1)
      e0max=e0(max0)
      e00=e0max 
      do i=max0+1,nmax
        e00=e00+de0
        e0(i)=e00
        a0(i)=0.0
        if(e00-e0max.gt.u) goto 10
      end do
   10 max=i
cudo  this next line is not in the orginal program
      nlst=n1
c
      if(snorm.gt.0.e0) goto 200
c
      evmin=e0(1)-1.
      evmax=e0(max)+1.
    5 n1=idnint((evmax-evmin)/estep)+1
      if(n1.gt.n1m)then
        estep=estep+estep
        goto 5
      endif
c create nice numbers for energy:
      i0=idnint(evmin/estep)
      evmin=i0*estep
c
c     set up plot-grid and determine min-energy
c
      nlst=n1
      write(6,*)' I know what nlst is:',nlst
      do 20 i=1,n1
      e1(i)=evmin+(i-1)*estep
   20 if(e1(i).gt.e0(max).and.nlst.eq.n1)nlst=i
   22 continue
      do 25 i=1,n0
      if(a0(i).gt.0.e0) goto 27
   25 continue
   27 emin0=e0(i)
      do 30 i=1,n1
      if(e1(i).gt.emin0) goto 40
   30 continue
   40 mine1=i
c
  200 continue
c
c     interpolate
c
      do 50 i=mine1,nlst
      j=0
   45 j=j+1
      if(e0(j).lt.e1(i).and.j.lt.max)goto45
      j=j-2
      if(j.lt.1)j=1
      if(j.gt.max-3)j=max-3
      call inter(e0(j),a0(j),4,e1(i),ai(i))
   50 continue
      do i=1,mine1-1
        ai(i)=0.
      end do
      do i=nlst+1,n1
        ai(i)=0.
      end do
c
      do 90 i=mine1,n1
c
c     determine halfwidth for life-time broadening 
c     (2 final valence holes + 1 initial core hole)
c
      dele=(e1(nlst)-e1(i))/(e1(nlst)-e1(mine1))
      gam=(hwl*dele*dele+hwc)/2.0e0
c
c     convolute now
c
   90 a1(i)=convlo(gam,ai,i,estep,mine1,n1,pi)
c
c     overwrite unbroadened data in ai with life-time
c     broadened spectrum
c
      do 95 i=1,n1
   95 ai(i)=a1(i)
      if(hws.ne.0.)then
c
c        set up gaussian convolution function
c
         afac=-alog(5.e-1)/(hws/2.0e0)**2
         bfac=sqrt(afac/pi)
         do 80 i=-ngc,ngc
            gc(i)=exp(-afac*(i*estep)**2)*bfac
   80    continue
c
c        convolute now
c
         do 85 i=mine1,n1
            a1(i)=convgau(ai,i,gc,ah,estep,ngc,n1)
   85    continue
      endif
C
      if(snorm.gt.0.e0) goto 300
c
c     determine maximum and set equal to 100.
c
      smax=0.e0
      do 100 i=mine1,n1
      if(a1(i).gt.smax) smax=a1(i)
  100 continue
      snorm=100.d0/smax
c
  300 continue
c
      do 110 i=mine1,n1
      a1(i)=a1(i)*snorm
c
  110 continue
c
      write(6,*)' Out of broad'
      return
      end
      double precision function convlo(gam,a1,ind,estep,mine1,n1,pi)
c ==================================================================
      implicit real*8 (a-h,o-z)
c
c     do the convolution integral by trapezoidal rule
C
      dimension a1(n1)
      write(6,*)' Out of convlo'
      sum=0.0
      do 10 i=mine1,n1
      atp=atan(estep*((i-ind)+0.5)/gam)
      atm=atan(estep*((i-ind)-0.5)/gam)
   10 sum=sum+a1(i)*(atp-atm)
      convlo=sum/pi
      write(6,*)' Out of convlo'
      return
      end
      double precision function convgau(a1,ind,gc,ah,estep,ngc,n1)
c ================================================================
      implicit real*8 (a-h,o-z)
c
c     do the convolution integral by trapezoidal rule
c
      dimension a1(n1), gc(-ngc:ngc), ah(-ngc:ngc)
c
      write(6,*)' In convgau'
      do 10 i=-ngc,ngc
   10 ah(i)=0.d0
      do 20 i=0,ngc
      imx=ind+i
      if(imx.gt.n1) goto 30
   20 ah(i)=a1(imx)*gc(i)
   30 do 40 i=-1,-ngc,-1
      imx=ind+i
      if(imx.lt.1) goto 50
   40 ah(i)=a1(imx)*gc(i)
   50 sum=0.d0
      do 60 i=-ngc,ngc
   60 sum=sum+ah(i)
c
      convgau=sum*estep
c
      write(6,*)' Out of convgau'
      return
      end
      subroutine matinter(crsfile,efile,nume,efile2,nume2,
     *nemf,nemf2,km,irange,e2,emin,emax,e,nemax,crs,
     *ifst,ilast,idoc,jprt)
c ========================================================
c
      implicit real*8(a-h,o-z)
c
c interpolation of CVV matrix elements on DOS energy grid for
c a special final state energy
c
c input: crsfile      matrix elements in the of input file
c        efile        valence band energy of input file
c        nume         number of VB energies actually read in
c        efile2       final state energy of input file
c        nume2        number of FS energies read in
c        nemf,nemf2   array dimensions
c        irange       =1, if crs has been calculated at a special
c                     e and e2; 0 otherwise
c        e2           final state energy
c        e            DOS energy grid
c        nemax        dimension of e
c        ifst, ilast  indices of e between which crs should be
c                     interpolated
c        idoc         if =1: documentation
c
c output: crs
c
      parameter(kml=3,nm=50)
c
      dimension crsfile(nemf,nemf2,-kml-1:kml,-kml-1:kml),
     *efile(nemf),efile2(nemf2),
     *irange(nemf,nemf2),e(nemax),crs(nemax,-kml-1:kml,-kml-1:kml)
      dimension ehelp(nm),chelp(nm),crs1(nm,-kml-1:kml,-kml-1:kml)
      write(6,*)' In matinter'
c
c find efile2 points immediately left and right
c of the actual e2 (if2l,if2r):
      if2l=leftind(efile2,nume2,e2)
      if2r=iright(efile2,nume2,e2)
c
c valence band energy range: find efile points immediately
c left of emin and right of emax (ifl,ifr):
      ifl=leftind(efile,nume,emin)
      ifr=iright(efile,nume,emax)
c
c for interpolation lateron, we may have to calculate the boundary
c values of the crs in a different way (to understand the whole
c thing, draw a diagram e vs. e2)
      insleft=0
      insright=0
      if(irange(ifl,if2r).eq.0)insleft=1
      if(irange(ifr,if2l).eq.0)insright=1
c
c now calculate the crs's for e2 and all values of efile
c between ifl and ifr and store them in crs1:
c
c boundary values (interpolation along boundaries of e/e2 range):
      nint=5
      if(nume.lt.nint)nint=nume
      if(insleft.eq.1)then
        ebar=(efile2(nume2)+efile2(1))/2.
  295   ehelp(ifl)=e2-ebar+efile(1)
c        if(ehelp(ifl).ge.emin)then
c          ifl=ifl+1
c          goto 295
c        endif
        do 305 l=-km-1,km
         if(l.eq.0) goto 305
           do 304 lp=-km-1,l
            if(lp.eq.0) goto 304
              do 300 i=nume,1,-1
                 i2=nume2-nume+i
                 chelp(i)=crsfile(i,i2,l,lp)
  300         continue
              call lagrange(efile(1),chelp(1),nume,nint,ehelp(ifl),
     *        crs1(ifl,l,lp),iextra)
  304      continue
  305   continue
        if(idoc.eq.1)write(jprt,307)
  307   format(' non-grid value for lower e boundary')
      endif
      if(insright.eq.1)then
  297   ehelp(ifr)=e2-efile2(1)+efile(1)
c        if(ehelp(ifr).le.emax)then
c          ifl=ifl-1
c          goto 297
c        endif
        do 315 l=-km-1,km
         if(l.eq.0) goto 315
         do 314 lp=-km-1,l
          if(lp.eq.0) goto 314
          do 310 i=1,nume
  310      chelp(i)=crsfile(i,i,l,lp)
           call lagrange(efile(1),chelp(1),nume,nint,ehelp(ifr),
     *     crs1(ifr,l,lp),iextra)
  314    continue
  315   continue
        if(idoc.eq.1)write(jprt,317)
  317   format(' non-grid value for upper e boundary')
      endif
      if(idoc.eq.1)write(jprt,320)ifl,ifr
  320 format(' ifl,ifr: ',2i4)
c
c get the rest of crs1:
      do 202 if=ifl,ifr
       if(if.eq.ifl.and.insleft.eq.1)goto 202
       if(if.eq.ifr.and.insright.eq.1)goto 202
       ehelp(if)=efile(if)
        do 201 l=-km-1,km
         if(l.eq.0) goto 201
          do 200 lp=-km-1,l
           if(lp.eq.0) goto 200
           do 190 i=if,if+nume-1
  190      chelp(i)=crsfile(if,i,l,lp)
           call lagrange(efile2(if),chelp(if),nume,nint,
     *     e2,crs1(if,l,lp),iextra)
  200     continue
  201   continue
  202 continue
c
      if(idoc.eq.1)then
        write(jprt,204)
  204   format(' crs1:')
        do 207 i=ifl,ifr
        write(jprt,206)i,ehelp(i),
     $  ((crs1(i,l,lp),lp=-km-1,l),l=-km-1,-1),
     $  ((crs1(i,l,lp),lp=-km-1,-1),(crs1(i,l,lp),lp=1,l),l=1,km)
  207   continue  
  206   format(i3,f10.5,28e13.5)
      endif

c
c now we have a function crs1 on a mesh efile and we want a function
c crs on a mesh e and that's all
      nf=ifr-ifl+1
      nint=6
      if(nint.gt.nf)nint=nf
      do 210 ie=ifst,ilast
       do 209 l=-km-1,km
        if(l.eq.0) goto 209
         do 208 lp=-km-1,l
          if(lp.eq.0) goto 208
           do 205 i=ifl,ifr
  205      chelp(i)=crs1(i,l,lp)
          call lagrange(efile(ifl),chelp(ifl),nf,nint,e(ie),
     *    crs(ie,l,lp),iextra)
          if(iextra.eq.1.and.idoc.eq.1)write(jprt,220)ie,e(ie),
     $    ((crs(ie,k,kp),kp=-km-1,k),k=-km-1,-1),
     $    ((crs(ie,k,kp),kp=-km-1,-1),(crs(ie,k,kp),kp=1,k),k=1,km)
  208    continue
  209  continue
  210 continue
  220 format(' extrapolation: crs, ie=',i3,' e(ie)=',f10.5/
     *28e13.5)
c
c test output:
      if(idoc.eq.1)then
        write(jprt,230)e2
        do 225 ie=ifst,ilast
        write(jprt,250)ie,e(ie),
     $  ((crs(ie,l,lp),lp=-km-1,l),l=-km-1,-1),
     $  ((crs(ie,l,lp),lp=-km-1,-1),(crs(ie,l,lp),lp=1,l),l=1,km)
  225   continue
        write(jprt,260)
      endif
  230 format(' e2: ',f10.5)
  250 format(i3,f10.5,28e13.4)
  260 format(120('-'))
      write(6,*)' out of matinter'
      return
      end
      subroutine lagrange(x,y,ndim,npkt,xval,yinter,iextra)
c =========================================================
      implicit real*8 (a-h,o-z)
      dimension x(ndim),y(ndim)
c
c Interpolation mit Komfort
c iextra = 0, wenn keine Extrapolation
c          1, wenn Extrapolation nach oben
c         -1, wenn Extrapolation nach unten
c
      write(6,*)' in lagrange'
      iextra=0
      if(xval.lt.x(1))iextra=-1
      if(xval.gt.x(ndim))iextra=1
c
c Suche 1. Punkt rechts von xval, falls es ihn gibt:
      igr=0
   10 igr=igr+1
      if(x(igr).eq.xval)then
        yinter=y(igr)
        return
      endif
      if(x(igr).lt.xval.and.igr.lt.ndim)goto 10
c
c Bestimme 1. Punkt f. Interpolation
      ileft=npkt/2
      if(mod(npkt,2).ne.0)ileft=ileft+1
      ifst=igr-ileft
c
      if(ifst.le.0)ifst=1
      if(ifst.gt.ndim-npkt+1)ifst=ndim-npkt+1
c
c Lasst uns interpolieren:
      call inter(x(ifst),y(ifst),npkt,xval,yinter)
c Und das war's dann auch schon, Leute!
      write(6,*)' Out of Lagrange'
      return
      end
      subroutine inter(r,p,n,rs,ps)
c =================================
      implicit real*8 (a-h,o-z)
c
c     interpolate via lagrange
c
      dimension r(n),p(n)
      write(6,*)' In inter'
      ps=0.e0
      do 1 j=1,n
      term=1.e0
      denom=1.e0
      do 2 i=1,n
      if(i.eq.j) go to 2
      denom=denom*(r(j)-r(i))
      term=term*(rs-r(i))
    2 continue
    1 ps=ps+term *p(j)/denom
      write(6,*)' Out of inter'
      return
      end
      integer function leftind(field,ndim,value)
c ==============================================
c index of the field element immediately left of value
c if value is lower than field(1), 1 is returned
c      implict real*8 (a-h,o-z)
      real*8 field,value
      dimension field(ndim)
      do 10 i=ndim,1,-1
   10 if(field(i).lt.value)goto 20
   20 leftind=i
      if(leftind.lt.1)leftind=1
      return
      end
      integer function iright(field,ndim,value)
c =============================================
c index of the field element immediately right of value
c if value is greater than fieldhndim), ndim is returned
c      implict real*8 (a-h,o-z)
      real*8 field,value
      dimension field(ndim)
      do 10 i=1,ndim
   10 if(field(i).gt.value)goto 20
   20 iright=i
      if(iright.gt.ndim)iright=ndim
      return
      end
      double precision function simpson(fct,n,dx)
c ===============================================
      implicit real*8 (a-h,o-z)
c integration of fct using simpson integration
C equidistant point mesh
      dimension fct(n)
      sum=0.
      x0=fct(1)
c
      do 100 i=3,n,2
      x1=fct(i-1)
      x2=fct(i)
      sum=sum+dx*(x0+4*x1+x2)
      x0=x2
  100 continue
      simpson=sum/3.
      if(mod(n,2).eq.0)simpson=simpson+(fct(n)+fct(n-1))*dx/2.
      return
      end
      subroutine order(f,n,ior)
c =============================
      implicit real*8(a-h,o-z)
      dimension f(n),ior(n)
c
c f is to be arranged in decreasing order
      do 10 i=1,n
   10 ior(i)=i
c
      do 20 i=1,n-1
      do 20 j=i+1,n
      if(f(j).gt.f(i))then
        g=f(j)
        f(j)=f(i)
        f(i)=g
        l=ior(j)
        ior(j)=ior(i)
        ior(i)=l
      endif
   20 continue
      return
      end
