MODULE wicomm_g
  USE task_kinds,ONLY: rkind
  INTEGER:: N1
  REAL(rkind):: G1,G2,G3,G4,G5
END MODULE wicomm_g

MODULE libpdkf
  USE task_kinds,ONLY: rkind
  USE task_constants,ONLY: PI

  INTEGER:: np_tau
  REAL(rkind):: xi,eta,rnu

  PRIVATE
  PUBLIC func_pdkf_r
  PUBLIC func_pdkf_i
  PUBLIC func_pdkf_rab
  PUBLIC func_pdkf_iab
  PUBLIC pdkf_hift
  PUBLIC pdkf_ft
  PUBLIC pdkf_eul

CONTAINS

  !   *****  REAL PART  *****

  FUNCTION func_pdkf_r(tau)
    
    !   real part fucntion of tau

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau
    REAL(rkind):: func_pdkf_r
    REAL(rkind):: arg

    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_r=tau**np_tau*DEXP(arg)*DCOS(rnu*tau)
    ELSE
       func_pdkf_r=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_r

  !   *****  IMAGINARY PART  *****

  FUNCTION func_pdkf_i(tau)
    
    !   imaginary part fucntion of tau

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau
    REAL(rkind):: func_pdkf_i
    REAL(rkind):: arg

    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_i=tau**np_tau*DEXP(arg)*DSIN(rnu*tau)
    ELSE
       func_pdkf_i=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_i

  !   *****  REAL PART  *****

  FUNCTION func_pdkf_rab(tau,Btau,tauA)
    
    ! real part fucntion of tau, tau_max-tau, tau-tau_min

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau,Btau,tauA
    REAL(rkind):: func_pdkf_rab
    REAL(rkind):: arg,dummy

    dummy=Btau
    dummy=tauA
    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_rab=tau**np_tau*DEXP(arg)*DCOS(rnu*tau)
    ELSE
       func_pdkf_rab=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_rab

  !     *****  IMAGINARY PART  *****

  FUNCTION func_pdkf_iab(tau,Btau,tauA)
    
    ! imaginary part fucntion of tau, tau_max-tau, tau-tau_min

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau,Btau,tauA
    REAL(rkind):: func_pdkf_iab
    REAL(rkind):: arg,dummy

    dummy=Btau
    dummy=tauA
    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_iab=tau**np_tau*DEXP(arg)*DSIN(rnu*tau)
    ELSE
       func_pdkf_iab=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_iab

  !   *** plasma dispersion kernel function (zero-infinity DE integral) ***

  FUNCTION pdkf_hift(xi_,eta_,rnu_,np_tau_)
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi_,eta_,rnu_
    INTEGER,INTENT(IN):: np_tau_
    COMPLEX(rkind):: pdkf_hift
    REAL(rkind):: H0,EPS,SR,SI,ESR,ESI
    INTEGER:: ILST,M

    xi=xi_
    eta=eta_
    rnu=rnu_
    np_tau=np_tau_
    
    H0=0.5D0
    EPS=1.D-8
    ILST=0

    CALL DEHIFT(SR,ESR,H0,EPS,ILST,func_pdkf_r,'func_pdkf_r')
    CALL DEHIFT(SI,ESI,H0,EPS,ILST,func_pdkf_i,'func_pdkf_i')
    pdkf_hift=DCMPLX(SR,SI)
    RETURN
  END FUNCTION pdkf_hift

  !   *** plasma dispersion kernel function (finite DE integral) ***

  FUNCTION pdkf_ft(xi_,eta_,rnu_,np_tau_)
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi_,eta_,rnu_
    INTEGER,INTENT(IN):: np_tau_
    COMPLEX(rkind):: pdkf_ft
    REAL(rkind):: H0,EPS,tau_min,tau_max,SR,SI,ESR,ESI
    INTEGER:: ILST,M

    xi=xi_
    eta=eta_
    rnu=rnu_
    np_tau=np_tau_
    
    H0=0.5D0
    EPS=1.D-8
    ILST=0

    tau_min=xi/18.D0       ! xi/tau<18
    IF(eta.EQ.0.D0) THEN
       tau_max=0.D0
    ELSE
       tau_max=18.D0/eta      ! eta*tau<18
    END IF

    IF(tau_max.EQ.0.D0) THEN
       CALL DEHIFT(SR,ESR,H0,EPS,ILST,func_pdkf_r,'func_pdkf_r')
       CALL DEHIFT(SI,ESI,H0,EPS,ILST,func_pdkf_i,'func_pdkf_i')
    ELSE
       CALL DEFTAB(tau_min,tau_max,SR,ESR,H0,EPS,ILST,func_pdkf_rab, &
            'func_pdkf_rab')
       CALL DEFTAB(tau_min,tau_max,SI,ESI,H0,EPS,ILST,func_pdkf_iab, &
            'func_pdkf_iab')
    END IF
    pdkf_ft=DCMPLX(SR,SI)
    RETURN
  END FUNCTION pdkf_ft

  !   ***** plasma dispersion kernel function (Euler transformation) *****
  !             real part: 0-pai/2, pai/2-3*pai/2, 3*pai/2-5*pai/2,...
  !             imag part: 0-pai, pai-2*pai, 2*pai-3*pai,...

  FUNCTION pdkf_eul(xi,eta,rnu,np_tau)

    USE task_constants,ONLY: PI
    USE libde
      USE wicomm_G
      IMPLICIT NONE
      REAL(rkind),PARAMETER:: HP=0.5D0*PI
      INTEGER,PARAMETER:: LMAX=1000
      REAL(rkind),INTENT(IN):: xi,eta,rnu
      INTEGER,INTENT(IN):: np_tau
      COMPLEX(rkind):: pdkf_eul
      REAL(rkind),DIMENSION(LMAX)::  A,B
      INTEGER:: ILST,K,L,M
      REAL(rkind):: H0,EPS,SR1,SI1,SR,SI,ESR,ESI,SR2,SI2,PARITY,SKR,SKI

      n1=np_tau
      G2=HP
      G3=xi
      G4=eta
      G5=0.D0

      H0=0.5D0
      EPS=1.D-8
      ILST=0
      M=5

      SR1=0.D0
      SI1=0.D0
      IF(M.NE.0) THEN 
         DO K=M-1,0,-1
            IF(K.EQ.0) H0=0.05D0
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
      pdkf_eul=CMPLX(SR,SI,rkind)

      RETURN

 9000 WRITE(6,*) ' ## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
      RETURN
    END FUNCTION pdkf_eul


!     *****  REAL PART  *****

    FUNCTION FUNR(X,XM,XP)

      USE wicomm_g
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNR
      REAL(rkind):: Y1,T,T2,YY,AN2

      Y1=XM
      IF(INT(G1).EQ.0) THEN 
         T=0.5D0*G2*XP
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            IF(N1.EQ.1) THEN
               AN2=1.D0
            ELSEIF(N1.EQ.2) THEN
               AN2=T
            ENDIF
               FUNR=AN2*0.5*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR=0.D0
         ENDIF 
      ELSE
         T=G2*(X+2.D0*G1)
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            IF(N1.EQ.1) THEN
               AN2=1.D0
            ELSEIF(N1.EQ.2) THEN
               AN2=T
            ENDIF
               FUNR=AN2*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR=0.D0
         ENDIF
      ENDIF
      RETURN
    END FUNCTION FUNR

!     *****  IMAG PART  *****

    FUNCTION FUNI(X,XM,XP)

      USE wicomm_g
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
         IF(N1.EQ.1) THEN
            AN2=1.D0
         ELSEIF(N1.EQ.2) THEN
            AN2=T
         ENDIF
         FUNI=AN2*G2*DEXP(YY)*DSIN(T)
      ELSE
         FUNI=0.D0
      ENDIF
      RETURN
    ENDFUNCTION FUNI

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
!      INTERFACE
!         FUNCTION FUNR(X,XM,XP)
!           USE wicomm,ONLY: rkind,ikind
!           REAL(rkind):: FUNR
!           REAL(rkind),INTENT(IN):: X,XM,XP
!         END FUNCTION FUNR
!         FUNCTION FUNI(X,XM,XP)
!           USE wicomm,ONLY: rkind,ikind
!           REAL(rkind):: FUNI
!           REAL(rkind),INTENT(IN):: X,XM,XP
!         END FUNCTION FUNI
!      END INTERFACE
      REAL(rkind),PARAMETER:: HP=1.5707963267948966192D0

      REAL(rkind):: EPS1,H,X,CSRP,CSIP,ATPR,ATPI,ATMR,ATMI,EPSI
      REAL(rkind):: HN,HC,HS,CC,XM,XP,CTR,CTI,ATR,ATI
      INTEGER:: N,NP,NM,NMIN,IND,ND

      EPS1=EPS**0.75
      H=H0
      X=0.D0
      CSR=HP*H*FUNR(X,1.D0-X,1.D0+X)
      CSRP=0.D0
      CSI=HP*H*FUNI(X,1.D0-X,1.D0+X)
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
         CTR=HC*FUNR(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI(X,XM,XP)*CC*CC
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
         CTR=HC*FUNR(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI(X,XM,XP)*CC*CC
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
END MODULe libpdkf
