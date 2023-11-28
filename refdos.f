      implicit real*8 (a-h,o-z)
      parameter (ne0=500)
      parameter (lmaxp=3)
c
      character*80 filen,iddos
      character*1  dosunit
c
      dimension e0(ne0)
      dimension pdos(ne0,-lmaxp-1:lmaxp)
c
      data cfac/13.606/
      data pi/3.14159265359/
c
      open (unit=1,file='refdos.in')      
 1000 continue
c
c            input dos file 
c
      read (1,10,end=1001) filen
      open (unit=2,file=filen)
c
c            output dos file 
c
      read (1,10) filen
      open (unit=8,file=filen)
c
      read(1,*) a0
c
      read(1,*) 
c
      conv1=2.*pi/a0
      conv1=conv1*conv1
c
c
c ***************  read in densities of states  ***************
c
c in relativistic case the ordering of kappa's in the dos file:
c          -1   1  -2   2  -3   3  -4
c
c
      read(2,10) iddos
      read(2,12) dosunit
      read(2,*)  ef
      read(2,*)  lmax
      read(2,*)  ne
      read(2,*)
c
      if(dosunit.eq.'d') conv=conv1*cfac
      if(dosunit.eq.'r') conv=cfac
      if(dosunit.eq.'e') conv=1.0
c
      do i=1,ne
c
         read(2,*) e0(i),pdos(i,0),pdos(i,-1),
     >               (pdos(i,l),pdos(i,-l-1),l=1,lmax)
c
         e0(i)=(e0(i)-ef)*conv
         do kap=-lmax-1,lmax
           pdos(i,kap)=pdos(i,kap)/cfac
         end do
c
      end do
c      
      write(8,11) iddos
      write(8,*) ne
c
      do i=1,ne
c
        write(8,13) e0(i),(pdos(i,kap),kap=-lmax-1,-1),
     >              (pdos(i,kap),kap=1,lmax),pdos(i,0)
c
      end do      
c
   10 format(a80)
   11 format(1x,a80)
   12 format(a1)
   13 format(9f10.5)
c
      close(2)
      close(8)
      goto 1000
 1001 stop 
      end
