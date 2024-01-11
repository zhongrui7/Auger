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
!
!*==SINTG.f90 
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
!
!*==DINTG.f90 
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
!
!*==INTER1.f90 
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
