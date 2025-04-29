!     $Id$

!     ***** CALCULATE EM FIELD IN A WAVEGUIDE *****
!
!     KA 11 : COAXIAL        TEM
!        12 :                TE11
!        13 :                TM01
!        14 :                multi
!        15 :                TE11 17 degree rotation
!        21 : CIRCULAR       TE11
!        22 :                TM01
!        23 :                TE01
!        24 :                TE11 17 degree rotation
!        25 :                TE11 RHS rotation
!        26 :                TE11 LHS rotation
!        27 :                multi
!        31 : RECTANGULAR    TE01
!        32 :                TM11
!        33 :                TE11
!        34 :                multi
!        4m : PARALLEL PLATE CIRCULAR
!
!     ***** CALCULATE E FIELD IN A WAVEGUIDE *****

subroutine WFEFWG(NSD,CEWG)

  use wfcomm
  implicit none
  integer    :: NB,NSD,ND1,ND2,I
  real(8)    :: ASD(3),X,Y,Z
  complex(8) :: CEWG(0:NMDM),CEW(3,0:NMDM,3),CBW(3,0:NMDM)
  
  NB=-KASID(NSD)
  
  ND1=NDSID(1,NSD)
  ND2=NDSID(2,NSD)
  ASD(1)=XND(ND2)-XND(ND1)
  ASD(2)=YND(ND2)-YND(ND1)
  ASD(3)=ZND(ND2)-ZND(ND1)
  
  X=XND(ND1)
  Y=YND(ND1)
  Z=ZND(ND1)
  call WFEMWG(NB,X,Y,Z,CEW(1,0,1),CBW)
  
  X=0.5D0*(XND(ND1)+XND(ND2))
  Y=0.5D0*(YND(ND1)+YND(ND2))
  Z=0.5D0*(ZND(ND1)+ZND(ND2))
  call WFEMWG(NB,X,Y,Z,CEW(1,0,2),CBW)
  
  X=XND(ND2)
  Y=YND(ND2)
  Z=ZND(ND2)
  call WFEMWG(NB,X,Y,Z,CEW(1,0,3),CBW)
  
  do I=0,NMBDY(NB)
     CEWG(I)=ASD(1)*(CEW(1,I,1)+2*CEW(1,I,2)+CEW(1,I,3))/4.D0 &
          &          +ASD(2)*(CEW(2,I,1)+2*CEW(2,I,2)+CEW(2,I,3))/4.D0 &
          &          +ASD(3)*(CEW(3,I,1)+2*CEW(3,I,2)+CEW(3,I,3))/4.D0
  end do
  return
end subroutine WFEFWG

!     ***** CALCULATE B FIELD IN A WAVEGUIDE *****

subroutine WFBFWG(ND,CBWG)

  use wfcomm
  implicit none
  integer    :: NB,ND
  real(8)    :: X,Y,Z
  complex(8) :: CEW(3,0:NMDM),CBWG(3,0:NMDM)

  NB=-KANOD(ND)
  
  X=XND(ND)
  Y=YND(ND)
  Z=ZND(ND)
  call WFEMWG(NB,X,Y,Z,CEW,CBWG)
  return
end subroutine WFBFWG

!     ***** CALCULATE EM FIELD IN A WAVEGUIDE *****

SUBROUTINE WFEMWG(NB,X,Y,Z,CEW,CBW)

  use wfcomm
  implicit none
  integer    :: NB,KA,KA10,I,N1,N2
  real(8)    :: ROMEG,SPWRX,DX,X,DY,Y,DZ,Z,XN1,YN1,ZN1
  real(8)    :: XN2,YN2,ZN2,XN3,YN3,ZN3,RWG1,RWG2,R1,R2
  real(8)    :: RWG,WG1,WG2,D1,D2,ZWG,RL
  complex(8) :: CE,CBX,CBY,CBZ,CB,CEW,CBW
  
  DIMENSION CEW(3,0:NMDM),CBW(3,0:NMDM)
  DIMENSION CE(3,0:NMDM),CB(3,0:NMDM)
  
  ROMEG=2.D0*PI*RF*1.D6
  SPWRX=SQRT(PWRBDY(NB))
  
  DX=X-XPBDY(NB)
  DY=Y-YPBDY(NB)
  DZ=Z-ZPBDY(NB)
  XN1=XNBDY(1,NB)
  YN1=YNBDY(1,NB)
  ZN1=ZNBDY(1,NB)
  XN2=XNBDY(2,NB)
  YN2=YNBDY(2,NB)
  ZN2=ZNBDY(2,NB)
  XN3=XNBDY(3,NB)
  YN3=YNBDY(3,NB)
  ZN3=ZNBDY(3,NB)
  
  KA=KABDY(NB)
  KA10=KA/10
  
  IF(KA10.LE.3) THEN
     
!     ----- COAXIAL WG ------
     
     IF(KA10.EQ.1) THEN
        RWG1=SZBDY(1,NB)
        RWG2=SZBDY(2,NB)
        R1=XN1*DX+YN1*DY+ZN1*DZ
        R2=XN2*DX+YN2*DY+ZN2*DZ
        
!     ----- COAXIAL WG TEM------
        
        IF(KA.EQ.11) THEN
           CALL WFEWGA(0,0,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,0),CB(1,0))
           
!     ----- COAXIAL WG TE11 ------

        ELSEIF(KA.EQ.12) THEN
           CALL WFEWGA(0,1,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,0),CB(1,0))
           
!     ----- COAXIAL WG TM01 ------
           
        ELSEIF(KA.EQ.13) THEN
           CALL WFEWGA(1,0,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,0),CB(1,0))
           
!     ----- COAXIAL WG multi ------

        ELSEIF(KA.EQ.14) THEN
           CALL WFEWGA(0,0,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,0),CB(1,0))
           CALL WFEWGA(0,1,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,1),CB(1,1))
           CALL WFEWGA(1,0,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,2),CB(1,2))
           DO I=1,3
              CE(I,0)=CE(I,0)+CE(I,1)+CE(I,2)
              CB(I,0)=CB(I,0)+CB(I,1)+CB(I,2)
           ENDDO

!     ----- COAXIAL WG TE11 17 degree------

        ELSEIF(KA.EQ.15) THEN
           CALL WFEWGA(0,1,1,17.D0,RWG1,RWG2,R1,R2,ROMEG,&
                &                     CE(1,0),CB(1,0))
        ENDIF

        CALL WFEWGA(0,0,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
             &                  CE(1,1),CB(1,1))
        CALL WFEWGA(0,1,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
             &                  CE(1,2),CB(1,2))
        CALL WFEWGA(0,1,1,90.D0,RWG1,RWG2,R1,R2,ROMEG,&
             &                  CE(1,3),CB(1,3))
        CALL WFEWGA(1,0,1,0.D0,RWG1,RWG2,R1,R2,ROMEG,&
             &                  CE(1,4),CB(1,4))

!     ----- CIRCULAR WG ------

     ELSEIF(KA10.EQ.2) THEN
        RWG=SZBDY(1,NB)
        R1=XN1*DX+YN1*DY+ZN1*DZ
        R2=XN2*DX+YN2*DY+ZN2*DZ
        
!     ----- CIRCULAR WG TE11 ------

        IF(KA.EQ.21) THEN
           CALL WFEWGC(0,1,1,0.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))

!     ----- CIRCULAR WG TM01 ------

        ELSEIF(KA.EQ.22) THEN
           CALL WFEWGC(1,0,1,0.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))

!     ----- CIRCULAR WG TE01 ------

        ELSEIF(KA.EQ.23) THEN
           CALL WFEWGC(0,0,1,0.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))

!     ----- CIRCULAR WG TE11 17 degree ------

        ELSEIF(KA.EQ.24) THEN
           CALL WFEWGC(0,1,1,17.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))

!     ----- CIRCULAR WG TE11 RHS rotation ------

        ELSEIF(KA.EQ.25) THEN
           CALL WFEWGC(0,1,1, 0.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))
           CALL WFEWGC(0,1,1,90.D0,RWG,R1,R2,ROMEG,CE(1,1),CB(1,1))
           DO I=1,3
              CE(I,0)=(CE(I,0)-CII*CE(I,1))/SQRT(2.D0)
              CB(I,0)=(CB(I,0)-CII*CB(I,1))/SQRT(2.D0)
           ENDDO
           
!     ----- CIRCULAR WG TE11 LHS rotation ------

        ELSEIF(KA.EQ.26) THEN
           CALL WFEWGC(0,1,1, 0.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))
           CALL WFEWGC(0,1,1,90.D0,RWG,R1,R2,ROMEG,CE(1,1),CB(1,1))
           DO I=1,3
              CE(I,0)=(CE(I,0)+CII*CE(I,1))/SQRT(2.D0)
              CB(I,0)=(CB(I,0)+CII*CB(I,1))/SQRT(2.D0)
           ENDDO

!     ----- CIRCULAR WG multi

        ELSEIF(KA.EQ.27) THEN
           CALL WFEWGC(0,1,1,0.D0,RWG,R1,R2,ROMEG,CE(1,0),CB(1,0))
           CALL WFEWGC(1,0,1,0.D0,RWG,R1,R2,ROMEG,CE(1,1),CB(1,1))
           CALL WFEWGC(0,0,1,0.D0,RWG,R1,R2,ROMEG,CE(1,2),CB(1,2))
           DO I=1,3
              CE(I,0)=CE(I,0)+CE(I,1)+CE(I,2)
              CB(I,0)=CB(I,0)+CB(I,1)+CB(I,2)
           ENDDO
        ENDIF
        
        CALL WFEWGC(0,1,1, 0.D0,RWG,R1,R2,ROMEG,CE(1,1),CB(1,1))
        CALL WFEWGC(0,1,1,90.D0,RWG,R1,R2,ROMEG,CE(1,2),CB(1,2))
        CALL WFEWGC(1,0,1, 0.D0,RWG,R1,R2,ROMEG,CE(1,3),CB(1,3))
        CALL WFEWGC(0,0,1, 0.D0,RWG,R1,R2,ROMEG,CE(1,4),CB(1,4))

!     ----- RECTANGULAR WG ------
        
     ELSEIF(KA10.EQ.3) THEN
        WG1=SZBDY(1,NB)
        WG2=SZBDY(2,NB)
        D1=XN1*DX+YN1*DY+ZN1*DZ+0.5D0*WG1
        D2=XN2*DX+YN2*DY+ZN2*DZ+0.5D0*WG2
        
!     ----- RECTANGULAR WG TE01 ------
        
        IF(KA.EQ.31) THEN
           CALL WFEWGR(0,0,1,WG1,WG2,D1,D2,ROMEG,CE(1,0),CB(1,0))
           
!     ----- RECTANGULAR WG TM11 ------

        ELSEIF(KA.EQ.32) THEN
           CALL WFEWGR(1,1,1,WG1,WG2,D1,D2,ROMEG,CE(1,0),CB(1,0))

!     ----- RECTANGULAR WG TE11 ------

        ELSEIF(KA.EQ.33) THEN
           CALL WFEWGR(0,1,1,WG1,WG2,D1,D2,ROMEG,CE(1,0),CB(1,0))
           
!     ----- RECTANGULAR WG multi ------

        ELSEIF(KA.EQ.34) THEN
           CALL WFEWGR(0,0,1,WG1,WG2,D1,D2,ROMEG,CE(1,0),CB(1,0))
           CALL WFEWGR(1,1,1,WG1,WG2,D1,D2,ROMEG,CE(1,1),CB(1,1))
           CALL WFEWGR(0,1,1,WG1,WG2,D1,D2,ROMEG,CE(1,2),CB(1,2))
           
           DO I=1,3
              CE(I,0)=CE(I,0)+CE(I,1)+CE(I,2)
              CB(I,0)=CB(I,0)+CB(I,1)+CB(I,2)
           ENDDO
        ENDIF
        
        CALL WFEWGR(0,0,1,WG1,WG2,D1,D2,ROMEG,CE(1,1),CB(1,1))
        CALL WFEWGR(1,1,1,WG1,WG2,D1,D2,ROMEG,CE(1,2),CB(1,2))
        CALL WFEWGR(0,1,1,WG1,WG2,D1,D2,ROMEG,CE(1,3),CB(1,3))
        
     ENDIF
     
     CBX=-CB(2,0)
     CBY=+CB(1,0)
     CBZ= 0.D0
     
     CEW(1,0)=( CE(1,0)*XN1+CE(2,0)*XN2+CE(3,0)*XN3)*SPWRX
     CEW(2,0)=( CE(1,0)*YN1+CE(2,0)*YN2+CE(3,0)*YN3)*SPWRX
     CEW(3,0)=( CE(1,0)*ZN1+CE(2,0)*ZN2+CE(3,0)*ZN3)*SPWRX
     CBW(1,0)=( CBX*XN1+CBY*XN2+CBZ*XN3)*SPWRX
     CBW(2,0)=( CBX*YN1+CBY*YN2+CBZ*YN3)*SPWRX
     CBW(3,0)=( CBX*ZN1+CBY*ZN2+CBZ*ZN3)*SPWRX
     
     DO I=1,NMBDY(NB)
        CEW(1,I)=( CE(1,I)*XN1+CE(2,I)*XN2+CE(3,I)*XN3)
        CEW(2,I)=( CE(1,I)*YN1+CE(2,I)*YN2+CE(3,I)*YN3)
        CEW(3,I)=( CE(1,I)*ZN1+CE(2,I)*ZN2+CE(3,I)*ZN3)
        CBX= CB(2,I)
        CBY=-CB(1,I)
        CBZ= 0.D0
        CBW(1,I)=( CBX*XN1+CBY*XN2+CBZ*XN3)
        CBW(2,I)=( CBX*YN1+CBY*YN2+CBZ*YN3)
        CBW(3,I)=( CBX*ZN1+CBY*ZN2+CBZ*ZN3)
     ENDDO

!     ----- PARALLEL PLATE WG ------

  ELSEIF(KA10.EQ.4) THEN
     N1=KA-KA10*10
     N2=0
     ZWG=SZBDY(1,NB)
     DX=X
     DY=Y
     DZ=Z-ZPBDY(NB)+0.5D0*ZWG
     RL=SQRT(DX*DX+DY*DY)
     CALL WFEWGP(N1,N2,ZWG,DX,DY,DZ,ROMEG,CE(1,0),CB(1,0))
     
     CBX= DY/RL*CB(3,0)
     CBY=              -DX/RL*CB(3,0)
     CBZ= DX/RL*CB(2,0)-DY/RL*CB(1,0)
     
     IF(XNBDY(1,NB).GT.0.D0) THEN
        CEW(1,0)= CE(1,0)*SPWRX
        CEW(2,0)= CE(2,0)*SPWRX
        CEW(3,0)= CE(3,0)*SPWRX
        CBW(1,0)= CBX*SPWRX
        CBW(2,0)= CBY*SPWRX
        CBW(3,0)= CBZ*SPWRX
        CEW(1,1)= CONJG(CE(1,0))
        CEW(2,1)= CONJG(CE(2,0))
        CEW(3,1)= CONJG(CE(3,0))
        CBW(1,1)=-CONJG(CBX)
        CBW(2,1)=-CONJG(CBY)
        CBW(3,1)=-CONJG(CBZ)
     ELSE
        CEW(1,0)= CONJG(CE(1,0))*SPWRX
        CEW(2,0)= CONJG(CE(2,0))*SPWRX
        CEW(3,0)= CONJG(CE(3,0))*SPWRX
        CBW(1,0)= CONJG(CBX)*SPWRX
        CBW(2,0)= CONJG(CBY)*SPWRX
        CBW(3,0)= CONJG(CBZ)*SPWRX
        CEW(1,1)= CE(1,0)
        CEW(2,1)= CE(2,0)
        CEW(3,1)= CE(3,0)
        CBW(1,1)=-CBX
        CBW(2,1)=-CBY
        CBW(3,1)=-CBZ
     ENDIF
     
  ELSE
     if(nrank.eq.0) WRITE(6,*) 'XX WFEFWG: UNDEFINED NB,KA=',NB,KA
     STOP
  ENDIF

!      WRITE(6,'(A,I6,1P3E12.4)') & 
!     &     'NB,X,Y,Z=',NB,X,Y,Z
!      WRITE(6,'(A,1P6E12.4)') &
!     &     '    CEW=',CEW(1,0),CEW(2,0),CEW(3,0)
!      WRITE(6,'(A,1P6E12.4)') &
!     &     '    CBW=',CBW(1,0),CBW(2,0),CBW(3,0)

  RETURN
END SUBROUTINE WFEMWG

!     ***** FUNCTION BESSJY *****

FUNCTION BESSJY(RK)

  use wfcomm
  implicit none
  integer :: N1,NBES
  real(8) :: BESSJY,ARGA,RK,RBESA,ARGB,RBESB,BYB,BESSY,BJB,BESSJ,BYA,BJA
  COMMON /WFEFN1/ RBESB,RBESA,NBES
  
  N1=NBES
  ARGA=RK*RBESA
  ARGB=RK*RBESB
  BYB=BESSY(N1,ARGB)
  BJB=BESSJ(N1,ARGB)
  BYA=BESSY(N1,ARGA)
  BJA=BESSJ(N1,ARGA)
  BESSJY=BYB*BJA-BJB*BYA

  RETURN
END FUNCTION BESSJY

!     ***** FUNCTION BESSJYD *****

FUNCTION BESSJYD(RK)

  use wfcomm
  implicit none
  integer :: N1,NBES
  real(8) :: BESSJYD,ARGA,RK,RBESA,ARGB,RBESB,BDYB,BESSY,BDJB,BESSJ,BDYA,BDJA
  COMMON /WFEFN1/ RBESB,RBESA,NBES

  N1=NBES
  ARGA=RK*RBESA
  ARGB=RK*RBESB
  BDYB=N1*BESSY(N1,ARGB)/ARGB-BESSY(N1+1,ARGB)
  BDJB=N1*BESSJ(N1,ARGB)/ARGB-BESSJ(N1+1,ARGB)
  BDYA=N1*BESSY(N1,ARGA)/ARGA-BESSY(N1+1,ARGA)
  BDJA=N1*BESSJ(N1,ARGA)/ARGA-BESSJ(N1+1,ARGA)
  BESSJYD=BDYB*BDJA-BDJB*BDYA
! WRITE(6,'(A,I5,1P4E12.4)') 'N1,RK,ARGA,ARGB,BESSJYD=', &
!   &                         N1,RK,ARGA,ARGB,BESSJYD
  RETURN
END FUNCTION BESSJYD

!************************************************************************

SUBROUTINE NEWTON(XS,DXS,XR,FUNC,EPS,ILL)

  implicit none
  integer :: ILMAX,ILL
  real(8) :: DXS,DX,X,XS,Y,FUNC,EPS,XR,XN,YN
  EXTERNAL FUNC

  ILMAX=ILL
  IF(ILMAX.EQ.0) ILMAX=30
  ILL=0
  IF(ABS(DXS).LE.1.E-15) THEN
     ILL=900
     RETURN
  ENDIF
  
  DX=DXS
  X=XS
  Y=FUNC(X)
  IF(ABS(Y).LE.EPS) THEN
     XR=X
     RETURN
  ENDIF
  XN=X+DX
  
10 CONTINUE
  ILL=ILL+1
  IF(ILL.GT.ILMAX) THEN
     ILL=901
     XR=XN
     RETURN
  ENDIF
  YN=FUNC(XN)
  IF(ABS(YN).LE.EPS) THEN
     XR=XN
     RETURN
  ENDIF
  IF(ABS(YN-Y).LE.0.D0) THEN
     ILL=902
     RETURN
  ENDIF
  X=XN
  XN=X-DX*Y/(YN-Y)
  Y=YN
  IF(ABS(XN-X).LE.EPS) THEN
     XR=XN
     RETURN
  ENDIF
  GOTO 10
END SUBROUTINE NEWTON

!     ***** CO-AXIAL WAVEGUIDE *****

SUBROUTINE WFEWGA(MODE,N1,N2,TH,RWG1,RWG2,R1,R2,ROMEG,CE,CB)

  use wfcomm
  implicit none
  integer    :: NBES,N1,MODE,N2,ILL
  real(8)    :: EPS,RBESB,RWG1,RBESA,RWG2,RVIMP,FACTOR,AMPE,RL,R1,R2,THETA
  real(8)    :: XS,DX,X0,RKC,RK2,ROMEG,ARGB,ARGA,BDYB,BESSY,BDJB,BESSJ,BESSA,BESSB
  real(8)    :: ARG,BESS,FBESS,DBESS,TH,BYB,BJB,BDYA,BDJA,BESSDA,BESSDB
  real(8)    :: BESSJY,BESSJYD
  complex(8) :: CERHO,CEPHI,CEZ,CBRHO,CBPHI,CE(3),CBZ,CB(3),CK,CFACTOR,CAMPH,CAMPE

  COMMON /WFEFN1/ RBESB,RBESA,NBES
  EXTERNAL BESSJY,BESSJYD
  DATA EPS/1.D-6/
  
  RBESB=RWG1
  RBESA=RWG2
  NBES=N1
  
  IF(MODE.EQ.0) THEN
     IF(N1.EQ.0) THEN
        RVIMP=SQRT(RMU0/EPS0)
        FACTOR=PI*RVIMP*LOG(RWG1/RWG2)
        AMPE=RVIMP/SQRT(FACTOR)
        
        RL=SQRT(R1*R1+R2*R2)
        THETA=ATAN2(R2,R1)
        CERHO= AMPE/RL
        CEPHI= 0.D0
        CEZ  = 0.D0
        CBRHO= 0.D0
        CBPHI= AMPE/RL
        CBZ  = 0.D0
        
        CE(1)=CERHO*COS(THETA)-CEPHI*SIN(THETA)
        CE(2)=CERHO*SIN(THETA)+CEPHI*COS(THETA)
        CE(3)=CEZ
        CB(1)=(CBRHO*COS(THETA)-CBPHI*SIN(THETA))
        CB(2)=(CBRHO*SIN(THETA)+CBPHI*COS(THETA))
        CB(3)=CBZ
        RETURN
     ELSE
        IF(N2.LE.0) THEN
           if (nrank.eq.0) WRITE(6,*) 'XX WFEWGA: INVALID N2.LE.0'
           STOP
        ENDIF
        IF(N2.EQ.1) THEN
           XS=2.D0*N1/(RBESB+RBESA)
        ELSE
           XS=PI*(N2-1)/(RBESB-RBESA)
        ENDIF
        DX=1.D-6
        ILL=0
        CALL NEWTON(XS,DX,X0,BESSJYD,EPS,ILL)
        IF(ILL.GE.30) THEN
           if (nrank.eq.0) WRITE(6,*) 'XX NEWTON: ILL=',ILL
           STOP
        ENDIF
        
        RKC=X0
        RK2=(ROMEG/VC)**2-RKC**2
        IF(RK2.LT.0.D0) THEN
           CK=CII*SQRT(-RK2)
        ELSE
           CK=SQRT(RK2)
        ENDIF
        ARGB=RKC*RBESB
        ARGA=RKC*RBESA
        
        BDYB=N1*BESSY(N1,ARGB)/ARGB-BESSY(N1+1,ARGB)
        BDJB=N1*BESSJ(N1,ARGB)/ARGB-BESSJ(N1+1,ARGB)
        BESSA=BDYB*BESSJ(N1,ARGA)-BDJB*BESSY(N1,ARGA)
        BESSB=BDYB*BESSJ(N1,ARGB)-BDJB*BESSY(N1,ARGB)
         
        CFACTOR=0.5D0*PI*ROMEG*RMU0*CK/RKC**4 &
             &              *(ARGB**2*(1.D0-N1**2/ARGB**2)*BESSB**2 &
             &               -ARGA**2*(1.D0-N1**2/ARGA**2)*BESSA**2)
        IF(N1.NE.0) CFACTOR=0.5D0*CFACTOR
        CAMPH=1.D0/SQRT(ABS(CFACTOR))
        CAMPE=ROMEG*RMU0*CAMPH/CK
        
        RL=SQRT(R1*R1+R2*R2)
        ARG=RKC*RL
        BESS=BDYB*BESSJ(N1,ARG)-BDJB*BESSY(N1,ARG)
        FBESS=N1*BESS/ARG
        DBESS=FBESS-(BDYB*BESSJ(N1+1,ARG)-BDJB*BESSY(N1+1,ARG))
        
        THETA=ATAN2(R2,R1)
        CERHO=CAMPE*CK/RKC*FBESS*SIN(N1*THETA+TH*PI/180.D0)
        CEPHI=CAMPE*CK/RKC*DBESS*COS(N1*THETA+TH*PI/180.D0)
        CEZ  =0.D0
        CBRHO=  -CK*VC*CEPHI/ROMEG
        CBPHI=   CK*VC*CERHO/ROMEG
        CBZ  =CII*CK*VC*CAMPE/ROMEG*BESS*COS(N1*THETA+TH*PI/180.D0)
     ENDIF
  ELSE
     IF(N2.LE.0) THEN
        if (nrank.eq.0) WRITE(6,*) 'XX WFEWGA: INVALID N2.LE.0'
        STOP
     ENDIF
     XS=PI*N2/(RBESB-RBESA)
     
     DX=1.D-6
     ILL=0
     CALL NEWTON(XS,DX,X0,BESSJY,EPS,ILL)
     IF(ILL.GE.30) THEN
        if (nrank.eq.0) WRITE(6,*) 'XX NEWTON: ILL=',ILL
        STOP
     ENDIF
     
     RKC=X0
     RK2=(ROMEG/VC)**2-RKC**2
     IF(RK2.LT.0.D0) THEN
        CK=CII*SQRT(-RK2)
     ELSE
        CK=SQRT(RK2)
     ENDIF
     ARGB=RKC*RBESB
     ARGA=RKC*RBESA
     BYB=BESSY(N1,ARGB)
     BJB=BESSJ(N1,ARGB)
     BDYA=N1*BESSY(N1,ARGA)/ARGA-BESSY(N1+1,ARGA)
     BDJA=N1*BESSJ(N1,ARGA)/ARGA-BESSJ(N1+1,ARGA)
     BDYB=N1*BESSY(N1,ARGB)/ARGB-BESSY(N1+1,ARGB)
     BDJB=N1*BESSJ(N1,ARGB)/ARGB-BESSJ(N1+1,ARGB)
     BESSDA=BYB*BDJA-BJB*BDYA
     BESSDB=BYB*BDJB-BJB*BDYB
     
     CFACTOR=0.5D0*PI*ROMEG*EPS0*CK/RKC**4 &
          &           *(ARGB**2*BESSDB**2 &
          &            -ARGA**2*BESSDA**2)
     IF(N1.NE.0) CFACTOR=0.5D0*CFACTOR
     CAMPE=1.D0/SQRT(ABS(CFACTOR))
     
     RL=SQRT(R1*R1+R2*R2)
     ARG=RKC*RL
     BESS=BYB*BESSJ(N1,ARG)-BJB*BESSY(N1,ARG)
     FBESS=N1*BESS/ARG
     DBESS=FBESS-(BYB*BESSJ(N1+1,ARG)-BJB*BESSY(N1+1,ARG))
     
     THETA=ATAN2(R2,R1)
     CERHO=   -CAMPE*CK*DBESS/RKC*COS(N1*THETA+TH*PI/180.D0)
     CEPHI=    CAMPE*CK*FBESS/RKC*SIN(N1*THETA+TH*PI/180.D0)
     CEZ  = CII*CAMPE*    BESS    *COS(N1*THETA+TH*PI/180.D0)
     CBRHO= -ROMEG*CEPHI/(VC*CK)
     CBPHI=  ROMEG*CERHO/(VC*CK)
     CBRHO=  0.D0
  ENDIF
  CE(1)=CERHO*COS(THETA)-CEPHI*SIN(THETA)
  CE(2)=CERHO*SIN(THETA)+CEPHI*COS(THETA)
  CE(3)=CEZ
  CB(1)=CBRHO*COS(THETA)-CBPHI*SIN(THETA)
  CB(2)=CBRHO*SIN(THETA)+CBPHI*COS(THETA)
  CB(3)=CBZ
  RETURN
END SUBROUTINE WFEWGA

!     ***** CIRCULAR WAVEGUIDE *****

SUBROUTINE WFEWGC(MODE,N1,N2,TH,RWG,R1,R2,ROMEG,CE,CB)

  use wfcomm
  implicit none
  integer :: MODE,N1,N2
  real(8) :: ZFJ(0:5,6),ZDJ(0:5,6),RKC,RWG,ARG,ROMEG,BESSJ
  real(8) :: RL,R1,R2,BESS,DBESS,FBESS,THETA,TH
  complex(8) :: CK,CFACTOR,CAMPH,CAMPE,CERHO,CEPHI,CEZ,CBRHO,CBPHI,CBZ,CE(3),CB(3)

  DATA ZFJ/&
       & 2.40483D0, 3.83171D0, 5.13562D0, 6.38016D0, 7.58834D0, 8.77148D0,&
       & 5.52008D0, 7.01559D0, 8.41724D0, 9.76102D0,11.06471D0,12.33860D0,&
       & 8.65373D0,10.17347D0,11.61984D0,13.01520D0,14.37254D0,15.70017D0,&
       &11.79153D0,13.32369D0,14.79595D0,16.22347D0,17.61597D0,18.98013D0,&
       &14.93092D0,16.47063D0,17.95982D0,19.40942D0,20.82693D0,22.21780D0,&
       &18.07106D0,19.61586D0,21.11700D0,22.58273D0,24.01902D0,25.43034D0/
  DATA ZDJ/&
       & 3.83171D0, 1.84118D0, 3.05424D0, 4.20119D0, 5.31755D0, 6.41562D0,&
       & 7.01559D0, 5.33144D0, 6.70613D0, 8.01524D0, 9.28238D0,10.51986D0,&
       &10.17347D0, 8.53632D0, 9.96947D0,11.34592D0,12.68191D0,13.98719D0,&
       &13.32369D0,11.70600D0,13.17037D0,14.58585D0,15.96411D0,17.31284D0,&
       &16.47063D0,14.86359D0,16.34752D0,17.78875D0,19.19603D0,20.57551D0,&
       &19.61586D0,18.01555D0,19.51291D0,20.97248D0,22.40103D0,23.80358D0/
  
  IF(MODE.EQ.0) THEN
     RKC=ZDJ(N1,N2)/RWG
     ARG=(ROMEG/VC)**2-RKC**2
     IF(ARG.LT.0.D0) THEN
        CK=CII*SQRT(-ARG)
     ELSE
        CK=SQRT(ARG)
     ENDIF
     
     ARG=ZDJ(N1,N2)
     CFACTOR=0.5D0*PI*ROMEG*RMU0*CK/RKC**4 &
          &           *ARG**2*(1.D0-N1**2/ARG**2)*BESSJ(N1,ARG)**2
     IF(N1.NE.0) CFACTOR=0.5D0*CFACTOR
     CAMPH=1.D0/SQRT(ABS(CFACTOR))
     CAMPE=ROMEG*RMU0*CAMPH/CK
     
     RL=SQRT(R1*R1+R2*R2)
     ARG=RL*RKC
     BESS=BESSJ(N1,ARG)
     IF(ARG.LE.1.D-6) THEN
        IF(N1.EQ.1) THEN
           DBESS=0.5D0
           FBESS=0.5D0
        ELSE
           DBESS=0.D0
           FBESS=0.D0
        ENDIF
     ELSE
        DBESS=N1*BESS/ARG-BESSJ(N1+1,ARG)
        FBESS=N1*BESS/ARG
     ENDIF
     THETA=ATAN2(R2,R1)
     CERHO=CAMPE*CK*FBESS/RKC*SIN(N1*THETA+TH*PI/180.D0)
     CEPHI=CAMPE*CK*DBESS/RKC*COS(N1*THETA+TH*PI/180.D0)
     CEZ  =0.D0
     CBRHO=  -CK*VC*CEPHI/ROMEG
     CBPHI=   CK*VC*CERHO/ROMEG
     CBZ  =CII*CK*VC*CAMPE*BESS*COS(N1*THETA+TH*PI/180.D0)/ROMEG
  ELSE
     RKC=ZFJ(N1,N2)/RWG
     ARG=(ROMEG/VC)**2-RKC**2
     IF(ARG.LT.0.D0) THEN
        CK=CII*SQRT(-ARG)
     ELSE
        CK=SQRT(ARG)
     ENDIF
     
     ARG=ZFJ(N1,N2)
     DBESS=-BESSJ(N1+1,ARG)
     CFACTOR=0.5D0*PI*ROMEG*EPS0*CK/RKC**4 &
          &          *ARG**2*DBESS**2
     IF(N1.NE.0) CFACTOR=0.5D0*CFACTOR
     CAMPE=1.D0/SQRT(ABS(CFACTOR))
     
     RL=SQRT(R1*R1+R2*R2)
     ARG=ZFJ(N1,N2)*RL/RWG
     BESS=BESSJ(N1,ARG)
     IF(ARG.LE.1.D-6) THEN
        IF(N1.EQ.1) THEN
           DBESS=0.5D0
           FBESS=0.5D0
        ELSE
           DBESS=0.D0
           FBESS=0.D0
        ENDIF
     ELSE
        DBESS=N1*BESS/ARG-BESSJ(N1+1,ARG)
        FBESS=N1*BESS/ARG
     ENDIF
     THETA=ATAN2(R2,R1)
     CERHO=   -CAMPE*CK*DBESS/RKC*COS(N1*THETA+TH*PI/180.D0)
     CEPHI=    CAMPE*CK*FBESS/RKC*SIN(N1*THETA+TH*PI/180.D0)
     CEZ  = CII*CAMPE*    BESS   *COS(N1*THETA+TH*PI/180.D0)
     CBRHO= -ROMEG*CEPHI/(VC*CK)
     CBPHI=  ROMEG*CERHO/(VC*CK)
     CBRHO=  0.D0
  ENDIF
  CE(1)=CERHO*COS(THETA)-CEPHI*SIN(THETA)
  CE(2)=CERHO*SIN(THETA)+CEPHI*COS(THETA)
  CE(3)=CEZ
  CB(1)=CBRHO*COS(THETA)-CBPHI*SIN(THETA)
  CB(2)=CBRHO*SIN(THETA)+CBPHI*COS(THETA)
  CB(3)=CBZ
  RETURN
END SUBROUTINE WFEWGC


!     ***** RECTANGULAR WAVEGUIDE *****

SUBROUTINE WFEWGR(MODE,N1,N2,WG1,WG2,D1,D2,ROMEG,CE,CB)

  use wfcomm
  implicit none
  integer :: N1,N2,MODE
  real(8) :: RKC,WG1,WG2,ARG,ROMEG,ARG1,D1,ARG2,AMPH,D2
  complex(8) :: CK,CFACTOR,CAMPE,CE(3),CB(3)

  RKC=SQRT((N1*PI/WG1)**2+(N2*PI/WG2)**2)
  ARG=(ROMEG/VC)**2-RKC**2
  IF(ARG.LT.0.D0) THEN
     CK=CII*SQRT(-ARG)
  ELSE
     CK=SQRT(ARG)
  ENDIF
  
  ARG1=N1*PI*D1/WG1
  ARG2=N2*PI*D2/WG2
  
  IF(MODE.EQ.0) THEN
     CFACTOR=0.5D0*ROMEG*RMU0*CK/RKC**2 *WG1*WG2
     IF(N1.NE.0) CFACTOR=0.5D0*CFACTOR
     IF(N2.NE.0) CFACTOR=0.5D0*CFACTOR
     AMPH=1.D0/SQRT(ABS(CFACTOR))
     CAMPE=ROMEG*RMU0*AMPH/CK
     CE(1)= CAMPE*N2*PI/WG2*CK/RKC**2*COS(ARG1)*SIN(ARG2)
     CE(2)=-CAMPE*N1*PI/WG1*CK/RKC**2*SIN(ARG1)*COS(ARG2)
     CE(3)= 0.D0
     CB(1)=  -CK*VC*CE(2)/ROMEG
     CB(2)=   CK*VC*CE(1)/ROMEG
     CB(3)=CII*CK*VC*CAMPE/ROMEG          *COS(ARG1)*COS(ARG2)
  ELSE
     CFACTOR=0.5D0*ROMEG*EPS0*CK/RKC**2 *WG1*WG2
     IF(N1.NE.0) CFACTOR=0.5D0*CFACTOR
     IF(N2.NE.0) CFACTOR=0.5D0*CFACTOR
     CAMPE=1.D0/SQRT(ABS(CFACTOR))
     CE(1)=   -CAMPE*N1*PI/WG1*CK/RKC**2*COS(ARG1)*SIN(ARG2)
     CE(2)=   -CAMPE*N2*PI/WG2*CK/RKC**2*SIN(ARG1)*COS(ARG2)
     CE(3)= CII*CAMPE                   *SIN(ARG1)*SIN(ARG2)
     CB(1)=-ROMEG*CE(2)/(VC*CK)
     CB(2)= ROMEG*CE(1)/(VC*CK)
     CB(3)= 0.D0
  ENDIF
  RETURN
END SUBROUTINE WFEWGR

!     ***** CIRCULAR PARALLEL PLATE *****

SUBROUTINE WFEWGP(N1,N2,ZWG,DX,DY,DZ,ROMEG,CE,CB)

!     TEM ONLY, N1: azimuthal mode number, N2: axial mode number

  use wfcomm
  implicit none
  integer :: N2,N1
  real(8) :: RK0,ROMEG,RKZ,ZWG,RKR2,RKR,FACTOR,AMPE,RL,DX,DY,PHI
  real(8) :: ARG,FBESJ,BESSJ,DBESJ,FBESN,BESSY,DBESN,DZ 
  complex(8) :: CE(3),CB(3),CFH,CDH,CEPHI,CBZ,CEZ,CBPHI

  RK0=ROMEG/VC
  RKZ=N2*PI/ZWG
  RKR2=RK0**2-RKZ**2
  IF(RKR2.LT.0.D0) THEN
     if (nrank.eq.0) WRITE(6,*) 'XX WFEWGP: IMAGINARY WAVE NUMBER'
     if (nrank.eq.0) WRITE(6,*) 'XX WFEWGP: RF_CUTOFF=',RKZ*VC/(2.D0*PI)
     STOP
  ENDIF
  RKR=SQRT(RKR2)
  FACTOR=2.D0*ZWG/(ROMEG*RMU0)
  IF(N1.NE.0) FACTOR=0.5D0*FACTOR
  IF(N2.NE.0) FACTOR=0.5D0*FACTOR
  AMPE=1.D0/SQRT(FACTOR)
!      WRITE(6,'(A,1PE12.4)') 'AMPE=',AMPE
  
  RL=SQRT(DX*DX+DY*DY)
  PHI=ATAN2(DY,DX)
  
  ARG=RKR*RL
  FBESJ=BESSJ(N1,ARG)
  DBESJ=N1*FBESJ/ARG-BESSJ(N1+1,ARG)
  FBESN=BESSY(N1,ARG)
  DBESN=N1*FBESN/ARG-BESSY(N1+1,ARG)
  CFH=FBESJ+CII*FBESN
  CDH=DBESJ+CII*DBESN
!  WRITE(6,'(A,1P3E12.4)') 'BES=',ARG,FBESJ,FBESN
  IF(N1.EQ.0) THEN
     CEPHI=0.D0
     CBZ=0.D0
  ELSE
     CEPHI=-RKZ*RL/N1 *CFH*SIN(N1*PHI)*SIN(RKZ*DZ)
     CBZ=CII*VC/ROMEG *RKR*RKZ/N1 *(2.D0*CFH+ARG*CDH) &
          &                         *SIN(N1*PHI)*SIN(RKZ*DZ)
  ENDIF
  CEZ=                    CFH*COS(N1*PHI)*COS(RKZ*DZ)
  CBPHI= CII*VC/ROMEG *RKR*CDH*COS(N1*PHI)*COS(RKZ*DZ)
  
  CE(1)=-AMPE*CEPHI*SIN(PHI)
  CE(2)= AMPE*CEPHI*COS(PHI)
  CE(3)= AMPE*CEZ
  CB(1)=-AMPE*CBPHI*SIN(PHI)
  CB(2)= AMPE*CBPHI*COS(PHI)
  CB(3)= AMPE*CBZ
  RETURN
END SUBROUTINE WFEWGP


!     ******* SET WG BOUNDARY *******

SUBROUTINE SETWGB(IERR)

  use wfcomm
  implicit none
  integer :: IERR,M100,M10,M1,NB,NBP,NSRF,NE,IN0,IN1,IN2,IN3,NN0,NN1,NN2,NN3,NN
  real(8) :: RL,EPS

  DATA EPS/1.D-6/

  IERR=0
  M100=model_prof/100
  M10=model_prof/10
  M1=MOD(model_prof,10)
  
  IF(M100.EQ.1) THEN
     IF(M1.EQ.0) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        ZPBDY(1)=ZNDMAX
        YNBDY(2,1)=-1.D0
        ZNBDY(3,1)=-1.D0
     ELSEIF(M1.EQ.1) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        ZPBDY(1)=0.D0
        YNBDY(2,1)= 1.D0
        ZNBDY(3,1)= 1.D0
     ELSEIF(M1.EQ.2.OR.M1.EQ.4) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        ZPBDY(1)=ZNDMAX
        ZPBDY(2)=0.D0
        YNBDY(2,1)=-1.D0
        ZNBDY(3,1)=-1.D0
        YNBDY(2,2)= 1.D0
        ZNBDY(3,2)= 1.D0
     ELSEIF(M1.EQ.3.OR.M1.EQ.5) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        ZPBDY(1)=0.D0
        ZPBDY(2)=ZNDMAX
        YNBDY(2,1)= 1.D0
        ZNBDY(3,1)= 1.D0
        YNBDY(2,2)=-1.D0
        ZNBDY(3,2)=-1.D0
     ENDIF
     DO NB=1,NBMAX
        KDBDY(NB)='coaxial WG'
        KABDY(NB)=M10
        NMBDY(NB)=4
!        PWRBDY(NB)=1.D0
        PHABDY(NB)=0.D0
        XPBDY(NB)=0.D0
        YPBDY(NB)=0.D0
!        ZPBDY(NB)=0.D0
        XNBDY(1,NB)= 1.D0
        YNBDY(1,NB)= 0.D0
        ZNBDY(1,NB)= 0.D0
        XNBDY(2,NB)= 0.D0
!        YNBDY(2,NB)= 1.D0
        ZNBDY(2,NB)= 0.D0
        XNBDY(3,NB)= 0.D0
        YNBDY(3,NB)= 0.D0
!        ZNBDY(3,NB)= 1.D0
        SZBDY(1,NB)= RNDMAX
        SZBDY(2,NB)= RNDMIN
     ENDDO

  ELSEIF(M100.EQ.2) THEN
     IF(M1.EQ.0) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        ZPBDY(1)=ZNDMAX
        YNBDY(2,1)=-1.D0
        ZNBDY(3,1)=-1.D0
     ELSEIF(M1.EQ.1) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        ZPBDY(1)=0.D0
        YNBDY(2,1)= 1.D0
        ZNBDY(3,1)= 1.D0
     ELSEIF(M1.EQ.2.OR.M1.EQ.4) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        ZPBDY(1)=ZNDMAX
        ZPBDY(2)=0.D0
        YNBDY(2,1)=-1.D0
        ZNBDY(3,1)=-1.D0
        YNBDY(2,2)= 1.D0
        ZNBDY(3,2)= 1.D0
     ELSEIF(M1.EQ.3.OR.M1.EQ.5) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        ZPBDY(1)=0.D0
        ZPBDY(2)=ZNDMAX
        YNBDY(2,1)= 1.D0
        ZNBDY(3,1)= 1.D0
        YNBDY(2,2)=-1.D0
        ZNBDY(3,2)=-1.D0
     ENDIF
     DO NB=1,NBMAX
        KDBDY(NB)='circular WG'
        KABDY(NB)=M10
        NMBDY(NB)=4
!        PWRBDY(NB)=1.D0
        PHABDY(NB)=0.D0
        XPBDY(NB)=0.D0
        YPBDY(NB)=0.D0
!        ZPBDY(NB)=0.D0
        XNBDY(1,NB)= 1.D0
        YNBDY(1,NB)= 0.D0
        ZNBDY(1,NB)= 0.D0
        XNBDY(2,NB)= 0.D0
!        YNBDY(2,NB)= 1.D0
        ZNBDY(2,NB)= 0.D0
        XNBDY(3,NB)= 0.D0
        YNBDY(3,NB)= 0.D0
!        ZNBDY(3,NB)= 1.D0
        SZBDY(1,NB)= RNDMAX
        SZBDY(2,NB)= 0.0D0
     ENDDO
     
  ELSEIF(M100.EQ.3) THEN
     IF(M1.EQ.0) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        ZPBDY(1)=ZNDMAX
        XNBDY(2,1)= 1.D0
        ZNBDY(3,1)=-1.D0
     ELSEIF(M1.EQ.1) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        ZPBDY(1)=0.D0
        XNBDY(2,1)=-1.D0
        ZNBDY(3,1)= 1.D0
     ELSEIF(M1.EQ.2.OR.M1.EQ.4) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        ZPBDY(1)=ZNDMAX
        ZPBDY(2)=0.D0
        XNBDY(2,1)= 1.D0
        ZNBDY(3,1)=-1.D0
        XNBDY(2,2)=-1.D0
        ZNBDY(3,2)= 1.D0
     ELSEIF(M1.EQ.3.OR.M1.EQ.5) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        ZPBDY(1)=0.D0
        ZPBDY(2)=ZNDMAX
        XNBDY(2,1)=-1.D0
        ZNBDY(3,1)= 1.D0
        XNBDY(2,2)= 1.D0
        ZNBDY(3,2)=-1.D0
     ENDIF
     DO NB=1,NBMAX
        KDBDY(NB)='rectangular WG'
        KABDY(NB)=M10
        NMBDY(NB)=3
!        PWRBDY(NB)=1.D0
        PHABDY(NB)=0.D0
        XPBDY(NB)=0.D0
        YPBDY(NB)=0.D0
!        ZPBDY(NB)=0.D0
        XNBDY(1,NB)= 0.D0
        YNBDY(1,NB)= 1.D0
        ZNBDY(1,NB)= 0.D0
!        XNBDY(2,NB)= 1.D0
        YNBDY(2,NB)= 0.D0
        ZNBDY(2,NB)= 0.D0
        XNBDY(3,NB)= 0.D0
        YNBDY(3,NB)= 0.D0
!        ZNBDY(3,NB)= 1.D0
        SZBDY(1,NB)= YNDMAX-YNDMIN
        SZBDY(2,NB)= XNDMAX-XNDMIN
     ENDDO
     
  ELSEIF(M100.EQ.4) THEN
     IF(M1.EQ.0) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        XPBDY(1)=XNDMAX
        ZPBDY(1)=0.5D0*(ZNDMIN+ZNDMAX)
        XNBDY(1,1)=-1.D0
        SZBDY(1,1)=ZNDMAX-ZNDMIN
     ELSEIF(M1.EQ.1) THEN
        NBMAX=1
        PWRBDY(1)=1.D0
        XPBDY(1)=0.05D0
        ZPBDY(1)=0.5D0*(ZNDMIN+ZNDMAX)
        XNBDY(1,1)= 1.D0
        SZBDY(1,1)=ZNDMAX-ZNDMIN
     ELSEIF(M1.EQ.2.OR.M1.EQ.4) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        XPBDY(1)=XNDMAX
        XPBDY(2)=0.05D0
        ZPBDY(1)=0.5D0*(ZNDMIN+ZNDMAX)
        ZPBDY(2)=0.5D0*(ZNDMIN+ZNDMAX)
        XNBDY(1,1)=-1.D0
        XNBDY(1,2)= 1.D0
        SZBDY(1,1)=ZNDMAX-ZNDMIN
        SZBDY(1,2)=ZNDMAX-ZNDMIN
     ELSEIF(M1.EQ.3.OR.M1.EQ.5) THEN
        NBMAX=2
        PWRBDY(1)=1.D0
        PWRBDY(2)=0.D0
        XPBDY(1)=0.05D0
        XPBDY(2)=XNDMAX
        ZPBDY(1)=0.5D0*(ZNDMIN+ZNDMAX)
        ZPBDY(2)=0.5D0*(ZNDMIN+ZNDMAX)
        XNBDY(1,1)= 1.D0
        XNBDY(1,2)=-1.D0
        SZBDY(1,1)=ZNDMAX-ZNDMIN
        SZBDY(1,2)=ZNDMAX-ZNDMIN
     ENDIF
     DO NB=1,NBMAX
        KDBDY(NB)='parallel plate WG'
        KABDY(NB)=M10
        NMBDY(NB)=1
!        PWRBDY(NB)=1.D0
        PHABDY(NB)=0.D0
!        XPBDY(NB)=0.D0
        YPBDY(NB)=0.D0
!        ZPBDY(NB)=0.D0
!        XNBDY(1,NB)= 0.D0
        YNBDY(1,NB)= 1.D0
        ZNBDY(1,NB)= 0.D0
        XNBDY(2,NB)= 1.D0
        YNBDY(2,NB)= 0.D0
        ZNBDY(2,NB)= 0.D0
        XNBDY(3,NB)= 0.D0
        YNBDY(3,NB)= 0.D0
        ZNBDY(3,NB)= 1.D0
     ENDDO
  ENDIF

  DO NB=1,NBMAX
     NBP=0
     IF(NB.EQ.NBMAX.AND.(M1.EQ.4.OR.M1.EQ.5)) THEN
        KABDY(NB)=4
        DO NSRF=1,NSRFMAX
           NE=NESRF(NSRF)
           IN0=INSRF(NSRF)
           IN1=MOD(IN0,4)+1
           IN2=MOD(IN0+1,4)+1
           IN3=MOD(IN0+2,4)+1
           NN0=NDELM(IN0,NE)
           NN1=NDELM(IN1,NE)
           NN2=NDELM(IN2,NE)
           NN3=NDELM(IN3,NE)
           IF(ABS(ZND(NN1)-ZPBDY(NB)).LE.EPS.AND.&
                &            ABS(ZND(NN2)-ZPBDY(NB)).LE.EPS.AND.&
                &            ABS(ZND(NN3)-ZPBDY(NB)).LE.EPS) THEN
              NBP=NBP+1
              IF(NBP.LE.NBPM) THEN
!                 WRITE(6,'(A,4I8)') 'NB,NBP,NE,NN0=',NB,NBP,NE,NN0
                 NENBP(NBP,NB)=NE
                 NDNBP(NBP,NB)=NN0
              ENDIF
           ENDIF
        ENDDO
     ELSE
        IF(M100.LE.3) THEN
           DO NN=1,NNMAX
              IF(ABS(ZND(NN)-ZPBDY(NB)).LE.EPS) THEN
                 NBP=NBP+1
                 IF(NBP.LE.NBPM) THEN
                    NDNBP(NBP,NB)=NN
                 ENDIF
              ENDIF
           ENDDO
        ELSE
           DO NN=1,NNMAX
              RL=SQRT(XND(NN)**2+YND(NN)**2)
              IF((XNBDY(1,NB).GT.0.D0.AND.&
                   &                 RL.LE.XPBDY(NB)+0.001D0).OR.&
                   &               (XNBDY(1,NB).LT.0.D0.AND.&
                   &                 RL.GE.XPBDY(NB)-0.001D0)) THEN
                 NBP=NBP+1
                 IF(NBP.LE.NBPM) THEN
                    NDNBP(NBP,NB)=NN
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDIF
     NBPMAX(NB)=NBP
     IF(NBP.GT.NBPM) THEN
        if(nrank.eq.0) WRITE(6,*) 'XX SETWGB: NBP.GT.NBPM: NB,NBP,NBPM=',&
             &                                             NB,NBP,NBPM
        IERR=1
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE SETWGB
