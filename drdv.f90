!*==DVDR.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dvdr(R,Rv,Z,Nmt,R2dvdr)
!     ===============
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 a , R , R2dvdr , r2v , Rv , Z
   INTEGER i , NA , Nmt , NRAD
!*** End of declarations inserted by SPAG
!
!     calculate derivative of the potential times r**2
!     r**2*dv/dr=d(r**2*v)/dr-2.r*v
!     finally override array rv with potential
!
!     j.redinger april 1985
!
   PARAMETER (NRAD=400)
   PARAMETER (NA=NRAD*3)
!
   DIMENSION R(NRAD) , Rv(NRAD) , r2v(NRAD) , R2dvdr(NRAD) , a(NA)
!
   DO i = 1 , Nmt
      r2v(i) = R(i)*Rv(i)
   ENDDO
!
   CALL derspl(Nmt,R,r2v,R2dvdr,a)
!
   DO i = 1 , Nmt
      R2dvdr(i) = R2dvdr(i) - 2.D0*Rv(i)
      Rv(i) = Rv(i)/R(i)
   ENDDO
!
   DO i = Nmt + 1 , NRAD
      R2dvdr(i) = 0.D0
   ENDDO
END SUBROUTINE dvdr
!*==DERSPL.f90 processed by SPAG 8.02DA 00:54  4 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE derspl(N,X,F,D,A)
!     =================
!
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL*8 A , D , F , h1 , h2 , p , X
   INTEGER i , j , k , N
!*** End of declarations inserted by SPAG
!
   DIMENSION X(N) , F(N) , D(N) , A(N)
!
!     f(i) are the function values at the points x(i) for i=1,n
!     and the spline derivatives d(i) are found.
!     the dimension of a must not be less than 3*n.
!
   DO i = 2 , N
      IF ( X(i)<=X(i-1) ) THEN
         WRITE (6,99001) i
99001    FORMAT (' return from derspl  ',i3,' out of order')
         A(1) = 1.D0
         RETURN
      ENDIF
   ENDDO
   DO i = 1 , N
      j = 2
      IF ( i/=1 ) THEN
         j = N - 1
         IF ( i/=N ) THEN
            h1 = 1.D0/(X(i)-X(i-1))
            h2 = 1.D0/(X(i+1)-X(i))
            A(3*i-2) = h1
            A(3*i-1) = 2.D0*(h1+h2)
            A(3*i) = h2
            D(i) = 3*(F(i+1)*h2*h2+F(i)*(h1*h1-h2*h2)-F(i-1)*h1*h1)
            CYCLE
         ENDIF
      ENDIF
      h1 = 1.D0/(X(j)-X(j-1))
      h2 = 1.D0/(X(j+1)-X(j))
      A(3*i-2) = h1*h1
      A(3*i-1) = h1*h1 - h2*h2
      A(3*i) = -h2*h2
      D(i) = 2.D0*(F(j)*(h2*h2*h2+h1*h1*h1)-F(j+1)*h2*h2*h2-F(j-1)*h1*h1*h1)
   ENDDO
   p = A(4)/A(1)
   A(5) = A(5) - p*A(2)
   A(6) = A(6) - p*A(3)
   D(2) = D(2) - p*D(1)
   DO i = 3 , N
      k = 3*i - 4
      p = A(k+2)/A(k)
      A(k+3) = A(k+3) - p*A(k+1)
      D(i) = D(i) - p*D(i-1)
      IF ( i==N-1 ) THEN
         p = A(k+5)/A(k)
         A(k+5) = A(k+6) - p*A(k+1)
         A(k+6) = A(k+7)
         D(N) = D(N) - p*D(N-2)
      ENDIF
   ENDDO
   D(N) = D(N)/A(3*N-1)
   DO i = 3 , N
      j = N + 2 - i
      D(j) = (D(j)-A(3*j)*D(j+1))/A(3*j-1)
   ENDDO
   D(1) = (D(1)-D(2)*A(2)-D(3)*A(3))/A(1)
   A(1) = 0.D0
END SUBROUTINE derspl
