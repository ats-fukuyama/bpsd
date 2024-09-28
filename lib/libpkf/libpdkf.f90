MODULE libpdkf

  USE task_kinds,ONLY: rkind
  USE task_constants,ONLY: PI

  INTEGER:: N1
  REAL(rkind):: G1,G2,G3,G4,G5

  PRIVATE
  PUBLIC pdkf
  PUBLIC pdkf_A

CONTAINS

  !     ****** FUNCTION FOR INTEGRO-DIFFERENTIAL ANALYSIS ******
  !     --- plasma dispersion kernel function ---

  FUNCTION pdkf(X,ALPHA,BETA,RKY,M,N)
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,ALPHA,BETA,RKY
    INTEGER,INTENT(IN):: M,N
    COMPLEX(rkind):: pdkf
    INTEGER:: L

    N1=N
    CALL EUL(X,ALPHA,BETA,RKY,pdkf,M,L)

    RETURN
  END FUNCTION pdkf

!     *****  EULER TRANSFOMATION  *****

  SUBROUTINE EUL(X,ALPHA,BETA,RKY,CS,M,L)

    USE libgrf
      IMPLICIT NONE
      REAL(rkind),PARAMETER:: HP=0.5D0*PI
      INTEGER,PARAMETER:: LMAX=1000
      REAL(rkind),INTENT(IN):: X,RKY,ALPHA,BETA
      INTEGER,INTENT(IN):: M
      INTEGER,INTENT(OUT):: L
      COMPLEX(rkind),INTENT(OUT):: CS
      REAL(rkind),DIMENSION(LMAX)::  A,B
      INTEGER:: ILST,K
      REAL(rkind):: H0,EPS,SR1,SI1,SR,SI,ESR,ESI,SR2,SI2,PARITY,SKR,SKI
      REAL(rkind):: TMIN,YMIN,EYMAX
      INTEGER:: nx,nxmax
      REAL(rkind),ALLOCATABLE:: XX(:),YY(:,:)
      REAL(rkind):: xstep

      G2=HP
      G3=X/BETA
      G4=0.5D0*ALPHA*BETA
      G5=RKY*BETA

      TMIN=(G3*G3/(G4*G4+G5*G5))**0.25D0
      YMIN=0.5D0*(G3*G3/TMIN**2+(G4*G4+G5*G5)*TMIN**2)
      EYMAX=EXP(-YMIN)
      H0=0.5D0
      EPS=1.D-20
      ILST=0

      nxmax=999
      xstep=2.0/(nxmax+1)
      ALLOCATE(XX(nxmax),YY(nxmax,2))
      G1=DBLE(0)
      DO nx=1,nxmax
         xx(nx)=-1.D0+xstep*nx
         yy(nx,1)=FUNR(xx(nx),1.0-xx(nx),1.D0+xx(nx))
         yy(nx,2)=FUNI(xx(nx),1.0-xx(nx),1.D0+xx(nx))
      END DO
      CALL PAGES
      CALL grd1d(0,xx,yy,nxmax,nxmax,2,'@funr/i vs x@')
      CALL PAGEE
      DEALLOCATE(XX,YY)
      
      SR1=0.D0
      SI1=0.D0
      IF(M.NE.0) THEN 
         DO K=M-1,0,-1
            G1=DBLE(K)
            CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS,ILST)
            SR1=SR1+SR
            SI1=SI1+SI
         END DO
      ENDIF

      G1=DBLE(M)
      CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS,ILST)
      A(1)=SR
      SR2=0.5D0*SR
      B(1)=SI
      SI2=0.5D0*SI
      PARITY=-1.D0
      L=0

   30 CONTINUE
      L=L+1
      IF(L.GE.LMAX) GOTO 9000
      G1=DBLE(M+L)
      CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS,ILST)
      A(L+1)=SR*PARITY
      B(L+1)=SI*PARITY
      DO K=L,1,-1    
         A(K)=A(K+1)-A(K)
         B(K)=B(K+1)-B(K)
      END DO
      SKR=A(1)*PARITY*0.5D0**(L+1)
      SR2=SR2+SKR
      SKI=B(1)*PARITY*0.5D0**(L+1)
      SI2=SI2+SKI
      PARITY=-PARITY
      IF(DABS(SKR).GT.EPS.OR.DABS(SKI).GT.EPS) GOTO 30

      SR=(SR1+SR2)/SQRT(4.D0*G2)
      SI=(SI1+SI2)/SQRT(4.D0*G2)
      CS=DCMPLX(SR,SI)

      RETURN

 9000 WRITE(6,*) ' ## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
      RETURN
    END SUBROUTINE EUL


!     *****  REAL PART  *****

    FUNCTION FUNR(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNR
      REAL(rkind):: Y1,T,T2,YY,AN2

      Y1=XM
      IF(INT(G1).EQ.0) THEN 
         T=0.5D0*G2*XP
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
!         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
         IF(ABS(YY).LT.352.D0) THEN
            SELECT CASE(N1)
            CASE(1)
               AN2=1.D0
            CASE(2)
               AN2=T
            CASE DEFAULT
               WRITE(6,'(A,I6)') 'XX pdkf-funr: undefined N1: N1=',N1
               STOP
            END SELECT
            FUNR=AN2*0.5*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR=0.D0
         ENDIF
      ELSE
         T=G2*(X+2.D0*G1)
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            SELECT CASE(N1)
            CASE(1)
               AN2=1.D0
            CASE(2)
               AN2=T
            CASE DEFAULT
               WRITE(6,'(A,I6)') 'XX pdkf-funr: undefined N1: N1=',N1
               STOP
            END SELECT
            FUNR=AN2*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR=0.D0
         ENDIF
      ENDIF
      RETURN
    END FUNCTION FUNR

!     *****  IMAG PART  *****

    FUNCTION FUNI(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNI
      REAL(rkind):: Y1,Y2,T,T2,YY,AN2

      Y1=X
      Y2=XM
      T=G2*(XP+2.D0*G1)
      T2=T*T
      YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
      IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
         SELECT CASE(N1)
         CASE(1)
            AN2=1.D0
         CASE(2)
            AN2=T
         CASE DEFAULT
            WRITE(6,'(A,I6)') 'XX pdkf-funi: undefined N1: N1=',N1
            STOP
         END SELECT
         FUNI=AN2*G2*DEXP(YY)*DSIN(T)
      ELSE
         FUNI=0.D0
      ENDIF
      RETURN
    END FUNCTION FUNI

  !     ****** FUNCTION FOR INTEGRO-DIFFERENTIAL ANALYSIS ******
  !     --- plasma dispersion kernel function ---

  FUNCTION pdkf_A(X,ALPHA,BETA,RKY,M,N)
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,ALPHA,BETA,RKY
    INTEGER,INTENT(IN):: M,N
    COMPLEX(rkind):: pdkf_A
    INTEGER:: L

    N1=N
    CALL EUL_A(X,ALPHA,BETA,RKY,pdkf_A,M,L)

    RETURN
  END FUNCTION pdkf_A

!     *****  EULER TRANSFOMATION  *****

  SUBROUTINE EUL_A(X,ALPHA,BETA,RKY,CS,M,L)

      IMPLICIT NONE
      REAL(rkind),PARAMETER:: HP=0.5D0*PI
      INTEGER,PARAMETER:: LMAX=1000
      REAL(rkind),INTENT(IN):: X,RKY,ALPHA,BETA
      INTEGER,INTENT(IN):: M
      INTEGER,INTENT(OUT):: L
      COMPLEX(rkind),INTENT(OUT):: CS
      REAL(rkind),DIMENSION(LMAX)::  A,B
      INTEGER:: ILST,K
      REAL(rkind):: H0,EPS,SR1,SI1,SR,SI,ESR,ESI,SR2,SI2,PARITY,SKR,SKI
      REAL(rkind):: TMIN,YMIN,EYMAX

      G2=HP
      G3=X/BETA
      G4=0.5D0*ALPHA*BETA
      G5=RKY*BETA

      TMIN=(G3*G3/(G4*G4+G5*G5))**0.25D0
      YMIN=0.5D0*(G3*G3/TMIN**2+(G4*G4+G5*G5)*TMIN**2)
      EYMAX=EXP(-YMIN)
      WRITE(21,'(A,3ES12.4)') 'G3,G4,G5       =',G3,G4,G5
      WRITE(21,'(A,3ES12.4)') 'TMIN,YMIN,EYMAX=',TMIN,YMIN,EYMAX
      H0=0.5D0
      EPS=1.D-10
      ILST=0

      SR1=0.D0
      SI1=0.D0
      IF(M.NE.0) THEN 
         DO K=M-1,0,-1
            G1=DBLE(K)
            CALL DEFTC2_A(SR,SI,ESR,ESI,H0,EPS,ILST)
            SR1=SR1+SR
            SI1=SI1+SI
         END DO
      ENDIF

      G1=DBLE(M)
      CALL DEFTC2_A(SR,SI,ESR,ESI,H0,EPS,ILST)
      A(1)=SR
      SR2=0.5D0*SR
      B(1)=SI
      SI2=0.5D0*SI
      PARITY=-1.D0
      L=0

   30 CONTINUE
      L=L+1
      IF(L.GE.LMAX) GOTO 9000
      G1=DBLE(M+L)
      CALL DEFTC2_A(SR,SI,ESR,ESI,H0,EPS,ILST)
      WRITE(21,'(A,I4,6es12.4)') 'ML:',M+L,X,G3,SR,SI
      A(L+1)=SR*PARITY
      B(L+1)=SI*PARITY
      DO K=L,1,-1    
         A(K)=A(K+1)-A(K)
         B(K)=B(K+1)-B(K)
      END DO
      SKR=A(1)*PARITY*0.5D0**(L+1)
      SR2=SR2+SKR
      SKI=B(1)*PARITY*0.5D0**(L+1)
      SI2=SI2+SKI
      PARITY=-PARITY
      IF(DABS(SKR).GT.EPS.OR.DABS(SKI).GT.EPS) GOTO 30

      SR=(SR1+SR2)/SQRT(4.D0*G2)
      SI=(SI1+SI2)/SQRT(4.D0*G2)
      CS=DCMPLX(SR,SI)

      RETURN

 9000 WRITE(6,*) ' ## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
      RETURN
    END SUBROUTINE EUL_A


!     *****  REAL PART  *****

    FUNCTION FUNR_A(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNR_A
      REAL(rkind):: Y1,T,T2,YY,AN2

      Y1=XM
      IF(INT(G1).EQ.0) THEN 
         T=0.5D0*G2*XP
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
!         IF(ABS(YY).LT.352.D0) THEN
            SELECT CASE(N1)
            CASE(1)
               AN2=1.D0
            CASE(2)
               AN2=T
            CASE DEFAULT
               WRITE(6,'(A,I6)') 'XX pdkf-funr_a: undefined N1: N1=',N1
               STOP
            END SELECT
            FUNR_A=AN2*0.5*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR_A=0.D0
         ENDIF
         WRITE(21,'(A,4ES12.4)') 'XP,T,YY,FUNR_A=',XP,T,YY,FUNR_A
      ELSE
         T=G2*(X+2.D0*G1)
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            SELECT CASE(N1)
            CASE(1)
               AN2=1.D0
            CASE(2)
               AN2=T
            CASE DEFAULT
               WRITE(6,'(A,I6)') 'XX pdkf-funr_a: undefined N1: N1=',N1
               STOP
            END SELECT
            FUNR_A=AN2*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR_A=0.D0
         ENDIF
         WRITE(21,'(A,4ES12.4)') 'X, T,YY,FUNR_A=',X,T,YY,FUNR_A
      ENDIF
      RETURN
    END FUNCTION FUNR_A

!     *****  IMAG PART  *****

    FUNCTION FUNI_A(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNI_A
      REAL(rkind):: Y1,Y2,T,T2,YY,AN2

      Y1=X
      Y2=XM
      T=G2*(XP+2.D0*G1)
      T2=T*T
      YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
      IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
         SELECT CASE(N1)
         CASE(1)
            AN2=1.D0
         CASE(2)
            AN2=T
         CASE DEFAULT
            WRITE(6,'(A,I6)') 'XX pdkf-funi_i: undefined N1: N1=',N1
            STOP
         END SELECT
         FUNI_A=AN2*G2*DEXP(YY)*DSIN(T)
      ELSE
         FUNI_A=0.D0
      ENDIF
      RETURN
    END FUNCTION FUNI_A

!     *****  DOUBLE EXPONENTIAL FORMULA  *****

    SUBROUTINE DEFTC2(CSR,CSI,ESR,ESI,H0,EPS,ILST)

!        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (-1.D0, +1.D0)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)

      IMPLICIT NONE
      REAL(rkind),INTENT(OUT):: CSR,CSI ! Integral
      REAL(rkind),INTENT(OUT):: ESR,ESI ! Estimated error
      REAL(rkind),INTENT(IN)::  H0      ! Initial step size
      REAL(rkind),INTENT(IN)::  EPS     ! Convergence thrshold
      INTEGER,INTENT(IN)::  ILST    ! print out control: 0 for no print out
      REAL(rkind),PARAMETER:: HP=1.5707963267948966192D0

      REAL(rkind):: EPS1,H,X,CSRP,CSIP,ATPR,ATPI,ATMR,ATMI,EPSI
      REAL(rkind):: HN,HC,HS,CC,XM,XP,CTR,CTI,ATR,ATI
      INTEGER:: N,NP,NM,NMIN,IND,ND

      EPS1=EPS**0.75
      H=H0
      X=0.D0
      CSR=HP*H*FUNR_A(X,1.D0-X,1.D0+X)
      CSRP=0.D0
      CSI=HP*H*FUNI_A(X,1.D0-X,1.D0+X)
      CSIP=0.D0
      N=0
      NP=0
      NM=0
      NMIN=1

    5 IND=0
      ATPR=1.D0
      ATPI=1.D0
      ATMR=1.D0
      ATMI=1.D0
      ND=2
      EPSI=MAX(EPS1*H,2.D-17)
      IF(N.EQ.0) ND=1

   10 N=N+ND
      HN=DBLE(N)*H
      HC=HP*H*COSH(HN)
      IF(IND.NE.1) THEN
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR_A(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI_A(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         NP=NP+1
         ATR=ATPR
         ATPR=ABS(CTR)
         ATI=ATPI
         ATPI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATPR)**2+(ATI+ATPI)**2) &
                 .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
               IF(IND.EQ.-1) GO TO 100
               IND=1
            ENDIF
         ENDIF
      ENDIF

      IF(IND.NE.-1) THEN
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR_A(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI_A(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         NM=NM+1
         ATR=ATMR
         ATMR=ABS(CTR)
         ATI=ATMI
         ATMI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATMR)**2+(ATI+ATMI)**2) &
                 .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
               IF(IND.EQ.1) GO TO 100
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10

  100 CONTINUE
      ESR=ABS(CSR-CSRP)
      CSRP=CSR
      ESI=ABS(CSI-CSIP)
      CSIP=CSI
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) THEN
            WRITE(6,601) H,NP,NM,CSR
            WRITE(6,604) CSI
         ENDIF
         IF(H.LT.H0) THEN
            WRITE(6,602) H,NP,NM,CSR,ESR
            WRITE(6,605) CSI,ESI
         ENDIF
      ENDIF
      IF(SQRT(ESR*ESR+ESI*ESI) &
        .LT.EPS1*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) GOTO 200

      IF(N.GT.1000) THEN
         WRITE(6,*) 'XX DEFTC2: Loop count exceeds 1000'
         STOP
      ENDIF

      H=0.5D0*H
      CSR=0.5D0*CSR
      CSI=0.5D0*CSI
      NMIN=N/2
      N=-1
      GOTO 5

  200 RETURN

  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  604 FORMAT(1H ,13X,16X,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  605 FORMAT(1H ,13X,16X,1PD24.15,1PD14.5)
    END SUBROUTINE DEFTC2
    
!     *****  DOUBLE EXPONENTIAL FORMULA  *****

    SUBROUTINE DEFTC2_A(CSR,CSI,ESR,ESI,H0,EPS,ILST)

!        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (-1.D0, +1.D0)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)

      USE libgrf
      IMPLICIT NONE
      REAL(rkind),INTENT(OUT):: CSR,CSI ! Integral
      REAL(rkind),INTENT(OUT):: ESR,ESI ! Estimated error
      REAL(rkind),INTENT(IN)::  H0      ! Initial step size
      REAL(rkind),INTENT(IN)::  EPS     ! Convergence thrshold
      INTEGER,INTENT(IN)::  ILST    ! print out control: 0 for no print out
      REAL(rkind),PARAMETER:: HP=1.5707963267948966192D0

      REAL(rkind):: EPS1,H,X,CSRP,CSIP,ATPR,ATPI,ATMR,ATMI,EPSI
      REAL(rkind):: HN,HC,HS,CC,XM,XP,CTR,CTI,ATR,ATI
      INTEGER:: N,NP,NM,NMIN,IND,ND

      INTEGER:: IG,I
      REAL(rkind):: XGR(10000),YGR(10000,2)

      EPS1=EPS**0.75
      H=H0
      X=0.D0
      CSR=HP*H*FUNR_A(X,1.D0-X,1.D0+X)
      CSRP=0.D0
      CSI=HP*H*FUNI_A(X,1.D0-X,1.D0+X)
      CSIP=0.D0
      N=0
      NP=0
      NM=0
      NMIN=1
      ig=0

    5 IND=0
      ATPR=1.D0
      ATPI=1.D0
      ATMR=1.D0
      ATMI=1.D0
      ND=2
      EPSI=MAX(EPS1*H,2.D-17)
      IF(N.EQ.0) ND=1

   10 N=N+ND
      HN=DBLE(N)*H
      HC=HP*H*COSH(HN)
      IF(IND.NE.1) THEN
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR_A(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI_A(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         IG=IG+1
         XGR(IG)=G2*(X+2.D0*G1)
         YGR(IG,1)=CTR
         YGR(IG,2)=CTI
         NP=NP+1
         ATR=ATPR
         ATPR=ABS(CTR)
         ATI=ATPI
         ATPI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATPR)**2+(ATI+ATPI)**2) &
                 .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
               IF(IND.EQ.-1) GO TO 100
               IND=1
            ENDIF
         ENDIF
      ENDIF

      IF(IND.NE.-1) THEN
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR_A(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI_A(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         IG=IG+1
         XGR(IG)=X
         YGR(IG,1)=CTR
         YGR(IG,2)=CTI
         NM=NM+1
         ATR=ATMR
         ATMR=ABS(CTR)
         ATI=ATMI
         ATMI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATMR)**2+(ATI+ATMI)**2) &
                 .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
               IF(IND.EQ.1) GO TO 100
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10

  100 CONTINUE
      ESR=ABS(CSR-CSRP)
      CSRP=CSR
      ESI=ABS(CSI-CSIP)
      CSIP=CSI
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) THEN
            WRITE(6,601) H,NP,NM,CSR
            WRITE(6,604) CSI
         ENDIF
         IF(H.LT.H0) THEN
            WRITE(6,602) H,NP,NM,CSR,ESR
            WRITE(6,605) CSI,ESI
         ENDIF
      ENDIF
      IF(SQRT(ESR*ESR+ESI*ESI) &
        .LT.EPS1*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) GOTO 200

      IF(N.GT.1000) THEN
         WRITE(6,*) 'XX DEFTC2: Loop count exceeds 1000'
         STOP
      ENDIF

      H=0.5D0*H
      CSR=0.5D0*CSR
      CSI=0.5D0*CSI
      NMIN=N/2
      N=-1
      GOTO 5

200   CONTINUE
      WRITE(6,'(A,I8,3ES12.4)') &
           'DEFTC2_A: ',IG,CSR,CSI
      IF(ABS(CSR).GT.1.D-8) THEN
         DO I=1,IG
            WRITE(6,'(I6,3ES12.4)') I,XGR(I),YGR(I,1),YGR(I,2)
         END DO
         CALL PAGES
         CALL grd1d(0,XGR,YGR,10000,IG,2,'@f vs x@')
         CALL PAGEE
      END IF
      RETURN

  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  604 FORMAT(1H ,13X,16X,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  605 FORMAT(1H ,13X,16X,1PD24.15,1PD14.5)
    END SUBROUTINE DEFTC2_A
END MODULe libpdkf
