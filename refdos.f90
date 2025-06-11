program dos_conversion
   implicit none

   ! Parameters
   integer, parameter :: ne0 = 500
   integer, parameter :: lmaxp = 3
   real(8), parameter :: cfac = 13.606
   real(8), parameter :: pi = 3.14159265359

   ! Variable declarations
   real(8) :: a0, conv, conv1, e0(ne0), ef, pdos(ne0, -lmaxp-1:lmaxp)
   integer :: i, kap, l, lmax, ne
   character(80) :: filen, iddos
   character(1) :: dosunit

   ! Main loop
   do
      ! Input DOS file
      open(unit=1, file='refdos.in', status='old', action='read')
      read(1, '(a80)', iostat=i) filen
      if (i /= 0) exit  ! Exit on end-of-file
      open(unit=2, file=trim(filen), status='old', action='read')

      ! Output DOS file
      read(1, '(a80)') filen
      open(unit=8, file=trim(filen), status='replace', action='write')

      read(1, *) a0
      read(1, *)  ! Skip a blank line

      conv1 = 2.0_8 * pi / a0
      conv1 = conv1 * conv1

      ! Read the density of states
      read(2, '(a80)') iddos
      read(2, '(a1)') dosunit
      read(2, *) ef
      read(2, *) lmax
      read(2, *) ne
      read(2, *)  ! Skip a blank line

      ! Set conversion factor
      select case (dosunit)
         case ('d')
            conv = conv1 * cfac
         case ('r')
            conv = cfac
         case ('e')
            conv = 1.0_8
      end select

      ! Read and process DOS data
      do i = 1, ne
         read(2, *) e0(i), pdos(i,0), pdos(i,-1), &
                    (pdos(i,l), pdos(i,-l-1), l=1,lmax)

         e0(i) = (e0(i) - ef) * conv
         do kap = -lmax - 1, lmax
            pdos(i,kap) = pdos(i,kap) / cfac
         end do
      end do

      ! Write output
      write(8, '(1x,a80)') iddos
      write(8, *) ne

      do i = 1, ne
         write(8, '(9f10.5)') e0(i), (pdos(i,kap), kap=-lmax-1,-1), &
                              (pdos(i,kap), kap=1,lmax), pdos(i,0)
      end do

      close(2)
      close(8)
   end do

   close(1)
end program dos_conversion
