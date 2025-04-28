! trprof.f90

MODULE trprof

  USE trcomm,ONLY: rkind
  PRIVATE

  !  *** Define fixed profile of density and temperature ***

  INTEGER:: ntime_ngiven_max                ! number of time points
  INTEGER:: ndata_ngiven_max                ! number of coef data
  REAL(rkind):: rho_min_ngiven,rho_max_ngiven  ! range of fixed profile
  REAL(rkind),ALLOCATABLE:: time_ngiven(:)     ! time points t_i
  REAL(rkind),ALLOCATABLE:: coef_ngiven(:,:)   ! coef data for t_i<= t <t_{i+1}
                                            ! temperature profile
  INTEGER:: ntime_tgiven_max                ! number of time points
  INTEGER:: ndata_tgiven_max                ! number of coef data
  REAL(rkind) :: rho_min_tgiven,rho_max_tgiven  ! range of fixed profile
  REAL(rkind),ALLOCATABLE:: time_tgiven(:)     ! time points t_i
  REAL(rkind),ALLOCATABLE:: coef_tgiven(:,:)   ! coef data for t_i<= t <t_{i+1}

  PUBLIC tr_prof
  PUBLIC tr_prof_impurity
  PUBLIC tr_prof_current

  PUBLIC tr_set_ngiven  ! set coef matrix for n
  PUBLIC tr_set_tgiven  ! set coef matrix for nT
  PUBLIC tr_prof_ngiven ! set fixed density profile
  PUBLIC tr_prof_tgiven ! set fixed temperature profile
  PUBLIC tr_prep_ngiven ! read fixed density pfofile parameters
  PUBLIC tr_prep_tgiven ! read fixed temperature profile parameters

CONTAINS

!     ***********************************************************

!           SET INITIAL PROFILE

!     ***********************************************************

  SUBROUTINE tr_prof

      USE trcomm
      USE libfio
      USE libspl1d
      IMPLICIT NONE
      INTEGER:: NR, NS, NF
      REAL(rkind) :: PROF,qsurf,qaxis
      REAL(rkind),ALLOCATABLE:: &
           rs_prof(:),rn_prof(:),rdn_prof(:),uprof(:,:)
      INTEGER:: nrmax_prof,ierr,i
      REAL(rkind):: R1,RN1

      ! *** number of radial mesh ***
      
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      
      ! *** set radial mesh: RG: grid, RM: center  *** 
      
      DO NR=1,NRMAX
         RG(NR) = DBLE(NR)*DR
         RM(NR) =(DBLE(NR)-0.5D0)*DR
      END DO

      ! *** initialize radial variables ***

      DO NR=1,NRMAX
         VTOR(NR)=0.D0
         VPAR(NR)=0.D0
         VPRP(NR)=0.D0
         VPOL(NR)=0.D0
         WROT(NR)=0.D0
         DO NS=1,NSTM
            RN(NR,NS)=0.D0
            RT(NR,NS)=0.D0
            RU(NR,NS)=0.D0
         END DO
         DO NF=1,NFMAX
            RW(NR,NF)=0.D0
            RNF(NR,NF)=0.D0
            RTF(NR,NF)=0.D0
         END DO
      END DO

      ! *** read intial profile data from knam_prof ***
      
      SELECT CASE(model_prof)
      CASE(11)
         CALL FROPEN(21,knam_prof,1,0,'PN',ierr)
         IF(ierr.NE.0) THEN
            WRITE(6,'(A,I8)') &
                 'XX file open error: knam_prof: ierr=',ierr
            STOP
         END IF

         ! --- count number of data ---
         
         I=0
         READ(21,'(A)')
100      CONTINUE
         I=I+1
         READ(21,*,ERR=190,END=200) R1,RN1
         GO TO 100
190      WRITE(6,*) 'XX prof error'
         STOP
200      WRITE(6,'(A,I6)') '## prof_data size=',I-1
         nrmax_prof=I-1
         ALLOCATE(rs_prof(nrmax_prof),rn_prof(nrmax_prof))
         REWIND(21)

         ! --- read electron density profile data ---
         
         READ(21,'(A)')
300      CONTINUE
         DO I=1,nrmax_prof
            READ(21,*,ERR=390,END=400) rs_prof(I),rn_prof(I)
            rn_prof(I)=rn_PROF(I)*1.D-20
            WRITE(6,'(A,I8,2ES12.4)') 'read:',I,rs_prof(I),rn_prof(I)
         END DO
         GOTO 400
390      WRITE(6,*) 'XX prof data error'
         STOP
               
         ! --- set spline data for electron density profile ---
         
400      CONTINUE
         WRITE(6,*) nrmax_prof,rs_prof(1),rs_prof(nrmax_prof)
         ALLOCATE(rdn_prof(nrmax_prof))
         ALLOCATE(uprof(4,nrmax_prof))
         rdn_prof(1)=0.D0
         CALL SPL1D(rs_prof,rn_prof,rdn_prof,uprof,nrmax_prof,1,ierr)
         WRITE(6,*) nrmax_prof,rs_prof(1),rs_prof(nrmax_prof)
         IF(ierr.NE.0) THEN
            WRITE(6,*) 'XX prof spline error !'
            STOP
         END IF

         DO nr=1,nrmax
            ! --- set interpolated electron density profile ---
            CALL SPL1DF(RM(nr),RN(nr,1),rs_prof,uprof,nrmax_prof,ierr)
            IF(ierr.NE.0) THEN
               WRITE(6,*) nr,RM(nr),rs_prof(1),rs_prof(nrmax_prof)
               WRITE(6,*) 'XX prof splinef error !',ierr
               STOP
            END IF
            ! --- set non-electron density profile ---
            DO ns=2,nsmax
               RN(nr,ns)=PN(ns)/PN(1)*RN(nr,1)
            END DO
            WRITE(6,'(A,I6,5ES12.4)') &
                 'prof: ',nr,RM(nr),RN(nr,1),RN(nr,2),RN(nr,3),RN(nr,4)
            ! --- set temperature and rotation profile ---
            DO ns=1,nsmax
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1(NS))**PROFT2(NS)
               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1(NS))**PROFU2(NS)
               RU(NR,NS) = (PU(NS)-PUS(NS))*PROF+PUS(NS)
            END DO
         END DO

      CASE default
            ! --- set default profiles ---

            DO nr=1,nrmax
            DO ns=1,nsmax
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1(NS))**PROFN2(NS)
               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1(NS))**PROFT2(NS)
               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1(NS))**PROFU2(NS)
               RU(NR,NS) = (PU(NS)-PUS(NS))*PROF+PUS(NS)
            END DO
            END DO
      END SELECT

      DO nr=1,nrmax
         PEX(NR,1:NSM) = 0.D0
         SEX(NR,1:NSM) = 0.D0
         PRF(NR,1:NSM) = 0.D0
         RNF(NR,1:NFMAX) = 0.D0
         PBM(NR)=0.D0
         WROT(NR)=0.D0
         VTOR(NR)=0.D0

         IF(MDLEQ0.EQ.1) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1(7))**PROFU2(7)
            RN(NR,7) = (PN(7)-PNS(7))*PROF+PNS(7)
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1(8))**PROFU2(8)
            RN(NR,8) = (PN(8)-PNS(8))*PROF+PNS(8)
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         ELSE
            ANNU(NR)=0.D0
         ENDIF

         RW(NR,1:NFMAX) = 0.D0

         SUMPBM=SUMPBM+PBM(NR)
      ENDDO

      SELECT CASE(model_ngiven)
      CASE(1)
         CALL tr_prep_ngiven
         IF(time_ngiven(1).LE.0.D0) THEN
            DO nr=1,nrmax
               CALL tr_prof_ngiven(rm(nr),t,rn(nr,1))
               DO ns=2,nsmax
                  rn(nr,ns)=pn(ns)/(pz(ns)*pn(1))*rn(nr,1)
                  IF(ns.EQ.2) WRITE(6,'(I6,3ES12.4)') nr,rm(nr),rn(nr,1),rn(nr,2)
               END DO
            END DO
            CALL tr_prof_ngiven(1.D0,t,pns(1))
            pnss(1)=pns(1)
            DO ns=2,nsmax
               pns(ns)=pn(ns)/(pz(ns)*pn(1))*pns(1)
               pnss(ns)=pns(ns)
            END DO
         END IF
      CASE(2)
         CALL tr_prep_ngiven
         IF(time_ngiven(1).LE.0.D0) THEN
            DO nr=1,nrmax
               IF((rm(nr).GE.rho_min_ngiven).AND. &
                  (rm(nr).LE.rho_max_ngiven)) THEN
                  CALL tr_prof_ngiven(rm(nr),t,rn(nr,1))
                  DO ns=2,nsmax
                     rn(nr,ns)=pn(ns)/(pz(ns)*pn(1))*rn(nr,1)
                  END DO
               END IF
            END DO
            CALL tr_prof_ngiven(1.D0,t,pns(1))
            pnss(1)=pns(1)
            DO ns=2,nsmax
               pns(ns)=pn(ns)/(pz(ns)*pn(1))*pns(1)
               pnss(ns)=pns(1)
            END DO
         END IF
      END SELECT
      SELECT CASE(model_tgiven)
      CASE(1)
         CALL tr_prep_tgiven
         IF(time_tgiven(1).LE.0.D0) THEN
            DO nr=1,nrmax
               CALL tr_prof_tgiven(rm(nr),t,rt(nr,1))
               DO ns=2,nsmax
                  rt(nr,ns)=rt(nr,1)
               END DO
            END DO
            CALL tr_prof_tgiven(1.D0,t,pts(1))
            DO ns=2,nsmax
               pts(ns)=pts(1)
            END DO
         END IF
      CASE(2)
         CALL tr_prep_tgiven
         IF(time_tgiven(1).LE.0.D0) THEN
            DO nr=1,nrmax
               IF((rm(nr).GE.rho_min_tgiven).AND. &
                  (rm(nr).LE.rho_max_tgiven)) THEN
                  CALL tr_prof_tgiven(rm(nr),t,rt(nr,1))
                  DO ns=2,nsmax
                     rt(nr,ns)=rt(nr,1)
                  END DO
               END IF
            END DO
            CALL tr_prof_tgiven(1.D0,t,pts(1))
            DO ns=2,nsmax
               pts(ns)=pts(1)
            END DO
         END IF
      END SELECT

      qsurf=5.D0*(RA/RR)*(BB/RIPS)
      qaxis=0.5D0*qsurf
      DO nr=1,nrmax
         qpinv(nr)=1.D0/(qaxis+(qsurf-qaxis)*RG(nr)**2)
      END DO

    END SUBROUTINE tr_prof


    SUBROUTINE tr_prof_impurity

      USE trcomm
      IMPLICIT NONE
      REAL(rkind)   :: ANEAVE, ANI, ANZ, DILUTE, TE, TRZEC,TRZEFE
      INTEGER NR
      EXTERNAL TRZEC,TRZEFE
      
!     *** CALCULATE ANEAVE and ANC, ANFE ***

      ANEAVE=SUM(RN(1:NRMAX,1)*RM(1:NRMAX))*2.D0*DR
      SELECT CASE(MDLIMP)
      CASE(0)
         ANC (1:NRMAX)=0.D0
         ANFE(1:NRMAX)=0.D0
      CASE(1,3)
         DO NR=1,NRMAX
            ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC *1.D-2*RN(NR,1)
            ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE*1.D-2*RN(NR,1)
         END DO
      CASE(2,4)
         ANC (1:NRMAX)=PNC *RN(1:NRMAX,1)
         ANFE(1:NRMAX)=PNFE*RN(1:NRMAX,1)
      END SELECT

!     *** CALCULATE PZC,PZFE ***

      DO NR=1,NRMAX
         TE=RT(NR,1)
         PZC(NR)=TRZEC(TE)
         PZFE(NR)=TRZEFE(TE)
      ENDDO

!     *** Dilution of ION due to IMPURITY DENSITY ***

         DO NR=1,NRMAX
            ANI = SUM(PZ(2:NSMAX)*RN(NR,2:NSMAX))     ! main ion charge density
            ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)   ! imputity ion charge den
            DILUTE = 1.D0-ANZ/ANI                     ! dilution factor
            IF(DILUTE.LT.0.D0) THEN
               WRITE(6,*) 'XX trprof: negative DILUTE: reduce PNC/PNFE'
               STOP
            END IF
            RN(NR,2:NSMAX) = RN(NR,2:NSMAX)*DILUTE    ! main ions diluted
         ENDDO
         PNSS(1)=PNS(1)
         PNSS(2:NSMAX)=PNS(2:NSMAX)*DILUTE
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1)=PNSA(1)
            PNSSA(2:NSMAX)=PNSA(2:NSMAX)*DILUTE
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      CALL TRZEFF

    END SUBROUTINE tr_prof_impurity

!     *** CALCULATE PROFILE OF AJ(R) ***

    SUBROUTINE tr_prof_current

      USE trcomm
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind), DIMENSION(NRMAX) :: DSRHO
      REAL(rkind) :: FACT,FACTOR0,FACTORM,FACTORP,PROF
      REAL(rkind) :: SUML

!     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***

         DO NR=1,NRMAX
            IF((1.D0-RM(NR)**ABS(PROFJ1)).LE.0.D0) THEN
               PROF=0.D0
            ELSE
               PROF= (1.D0-RM(NR)**ABS(PROFJ1))**ABS(PROFJ2)
            ENDIF
            AJOH(NR)= PROF
            AJ(NR)  = PROF
         ENDDO

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO

         RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)=FACT*RDP(1:NRMAX)
         RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
         AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)  =AJOH(1:NRMAX)
         BP(1:NRMAX)  =FACT*BP(1:NRMAX)
         QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX) &
              /(4.D0*PI**2*RDPVRHOG(1:NRMAX))
         
         Q0=(20.D0*QP(1)-23.D0*QP(2)+8.D0*QP(3))/5.D0

!     *** THIS MODEL ASSUMES CONSTANT EZ ***

      IF(PROFJ1.LE.0.D0.OR.MDLNCL.EQ.1) THEN
         CALL TRZEFF
         CALL TRCFET
         IF(PROFJ1.GT.0.D0.AND.MDLNCL.EQ.1) GOTO 2000

         AJOH(1:NRMAX)=1.D0/ETA(1:NRMAX)
         AJ(1:NRMAX)  =1.D0/ETA(1:NRMAX)

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO
         BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)  =FACT*RDP(1:NRMAX)
         RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
         AJOH(1:NRMAX) =FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)   =AJOH(1:NRMAX)
         AJTOR(1:NRMAX)=FACT*AJTOR(1:NRMAX)
         BP(1:NRMAX)   =FACT*BP(1:NRMAX)
         QP(1:NRMAX)   =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      ENDIF
 2000 CONTINUE
      SUML=0.D0
      DO NR=1,NRMAX
         SUML=SUML+RDP(NR)*DR
         RPSI(NR)=SUML
         BPRHO(NR)=BP(NR)
      ENDDO

    END SUBROUTINE tr_prof_current

  !     ***** Routine for fixed density profile *****
      
  SUBROUTINE tr_set_ngiven(nr,time)

    USE trcomm
    USE trcomx
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr
    REAL(rkind):: rn_local
    REAL(rkind),INTENT(IN):: time
    INTEGER:: NS,NEQ,NW

    IF(model_ngiven.EQ.0) RETURN
    IF(time.LE.time_ngiven(1)) return
    IF(model_ngiven.EQ.2) THEN
       IF((rm(nr).LT.rho_min_ngiven).OR. &
            (rm(nr).GT.rho_max_ngiven)) RETURN
    END IF
    CALL tr_prof_ngiven(rm(nr),time,rn_local)
    NEQ=NEA(1,1) ! NEQ of electron density equation
    DO NW=1,NEQMAX
       A(NEQ,NW,NR) = 0.D0
       B(NEQ,NW,NR) = 0.D0
       C(NEQ,NW,NR) = 0.D0
    END DO
    D(NEQ,NR)=0.D0
!    B(NEQ,NEQ,NR)=-1.D0/tau_ngiven
!    D(NEQ,NR)=rn_local/tau_ngiven
    RD(NEQ,NR)=1.D0
    DO NS=2,NSMAX
       NEQ=NEA(NS,1) ! NEQ of density equation
       DO NW=1,NEQMAX
          A(NEQ,NW,NR) = 0.D0
          B(NEQ,NW,NR) = 0.D0
          C(NEQ,NW,NR) = 0.D0
       END DO
       D(NEQ,NR)=0.D0
!       B(NEQ,NEQ,NR)=-1.D0/tau_ngiven
!       D(NEQ,NR)=pn(ns)/(pz(ns)*pn(1))*rn_local/tau_ngiven
       RD(NEQ,NR)=1.D0
    END DO
    RETURN
  END SUBROUTINE tr_set_ngiven

  !     ***** Routine for fixed temperature profile *****
      
  SUBROUTINE tr_set_tgiven(nr,time)

    USE trcomm
    USE trcomx
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr
    REAL(rkind),INTENT(IN):: time
    REAL(rkind):: rt_local
    INTEGER:: NS,NEQ,NW

    IF(model_tgiven.EQ.0) RETURN
    IF(time.LE.time_tgiven(1)) return
    IF(model_tgiven.EQ.2) THEN
       IF((rm(nr).LT.rho_min_tgiven).OR. &
            (rm(nr).GT.rho_max_tgiven)) RETURN
    END IF
    CALL tr_prof_tgiven(rm(nr),time,rt_local)
    NEQ=NEA(1,1) ! NEQ of electron density equation
    DO NW=1,NEQMAX
       A(NEQ,NW,NR) = 0.D0
       B(NEQ,NW,NR) = 0.D0
       C(NEQ,NW,NR) = 0.D0
    END DO
    D(NEQ,NR)=0.D0
!    B(NEQ,NEQ,NR)=-1.D0/tau_tgiven
!    D(NEQ,NR)=rt_local/tau_tgiven
    RD(NEQ,NR)=1.D0
    DO NS=2,NSMAX
       NEQ=NEA(NS,2) ! NEQ of temperature equation
       DO NW=1,NEQMAX
          A(NEQ,NW,NR) = 0.D0
          B(NEQ,NW,NR) = 0.D0
          C(NEQ,NW,NR) = 0.D0
       END DO
       D(NEQ,NR)=0.D0
!       B(NEQ,NEQ,NR)=-1.D0/tau_tgiven
!       D(NEQ,NR)=rt_local/tau_tgiven
       RD(NEQ,NR)=1.D0
    END DO
    RETURN
  END SUBROUTINE tr_set_tgiven
      
  ! *** set fixed density profile ***
  
  SUBROUTINE tr_prof_ngiven(rho,time,rn_local)
  
    USE trcomm
    IMPLICIT NONE    
    REAL(rkind),INTENT(IN):: rho,time
    REAL(rkind),INTENT(OUT):: rn_local
    REAL(rkind):: tr_func_ngiven
    REAL(rkind),ALLOCATABLE:: coef(:)
    REAL(rkind):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LT.time_ngiven(1)) THEN
       RETURN
    ELSE IF (time.GE.time_ngiven(ntime_ngiven_max)) THEN
       id=ntime_ngiven_max
    ELSE
       DO ntime=1,ntime_ngiven_max-1
          IF(time.GE.time_ngiven(ntime).AND. &
               time.LT.time_ngiven(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_ngiven_max))
    IF(id.EQ.ntime_ngiven_max) THEN ! after time_ngiven(ntime_ngiven_max)
       DO i=0,ndata_ngiven_max
          coef(i)=coef_ngiven(i,ntime_ngiven_max)
       END DO
    ELSE ! between time_ngiven(id) and time_ngiven(id+1)
       factor=(time-time_ngiven(id)) &
             /(time_ngiven(id+1)-time_ngiven(id))
       DO i=0,ndata_ngiven_max
          coef(i)=(1.D0-factor)*coef_ngiven(i,id) &
                        +factor*coef_ngiven(i,id+1)
       END DO
    END IF

    ! --- set local density profile ---
    
    rn_local=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    rn_local=rn_local*1.D-20
    IF(rn_local.LE.0.D0) rn_local=1.D-8
    RETURN
  END SUBROUTINE tr_prof_ngiven

  ! *** set fixed temperature profile ***
  
  SUBROUTINE tr_prof_tgiven(rho,time,rt_local)
  
    USE trcomm
    IMPLICIT NONE    
    REAL(rkind),INTENT(IN):: rho,time
    REAL(rkind),INTENT(OUT):: rt_local
    REAL(rkind):: tr_func_tgiven
    REAL(rkind),ALLOCATABLE:: coef(:)
    REAL(rkind):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LE.time_tgiven(1)) THEN
       RETURN
    ELSE IF (time.GE.time_tgiven(ntime_tgiven_max)) THEN
       id=ntime_tgiven_max
    ELSE
       DO ntime=1,ntime_tgiven_max-1
          IF(time.GE.time_tgiven(ntime).AND. &
               time.LE.time_tgiven(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_tgiven_max))
    IF(id.EQ.0) THEN ! before time_ngiven(1)
       DO i=0,ndata_tgiven_max
          coef(i)=coef_tgiven(i,1)
       END DO
    ELSE IF(id.EQ.ntime_tgiven_max) THEN ! after time_ngiven(ntime_ngiven_max)
       DO i=0,ndata_tgiven_max
          coef(i)=coef_tgiven(i,ntime_tgiven_max)
       END DO
    ELSE ! between time_ngiven(id) and time_ngiven(id+1)
       factor=(time-time_tgiven(id)) &
             /(time_tgiven(id+1)-time_tgiven(id))
       DO i=0,ndata_tgiven_max
          coef(i)=(1.D0-factor)*coef_tgiven(i,id) &
                        +factor*coef_tgiven(i,id+1)
       END DO
    END IF

    ! --- set temperature profile ---
    
    rt_local=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    rt_local=rt_local*1.D-3
    IF(rt_local.LE.0.D0) rt_local=3.D-5
    RETURN
  END SUBROUTINE tr_prof_tgiven

  ! *** read density profile data from file ***

  SUBROUTINE tr_prep_ngiven
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_ngiven,1,0,'fn',ierr)
    READ(NFL,*) ntime_ngiven_max,ndata_ngiven_max,rho_min_ngiven,rho_max_ngiven
    IF(ALLOCATED(time_ngiven)) DEALLOCATE(time_ngiven)
    IF(ALLOCATED(coef_ngiven)) DEALLOCATE(coef_ngiven)
    ALLOCATE(time_ngiven(ntime_ngiven_max))
    ALLOCATE(coef_ngiven(0:ndata_ngiven_max,ntime_ngiven_max))
    DO ntime=1,ntime_ngiven_max
       READ(NFL,*) time_ngiven(ntime)
       READ(NFL,*) (coef_ngiven(ndata,ntime),ndata=0,ndata_ngiven_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_ngiven
    
  ! *** read temperature profile data from file ***

  SUBROUTINE tr_prep_tgiven
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_tgiven,1,0,'ft',ierr)
    READ(NFL,*) ntime_tgiven_max,ndata_tgiven_max,rho_min_tgiven,rho_max_tgiven
    IF(ALLOCATED(time_tgiven)) DEALLOCATE(time_tgiven)
    IF(ALLOCATED(coef_tgiven)) DEALLOCATE(coef_tgiven)
    ALLOCATE(time_tgiven(ntime_tgiven_max))
    ALLOCATE(coef_tgiven(0:ndata_tgiven_max,ntime_tgiven_max))
    DO ntime=1,ntime_tgiven_max
       READ(NFL,*) time_tgiven(ntime)
       READ(NFL,*) (coef_tgiven(ndata,ntime),ndata=0,ndata_tgiven_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_tgiven
  
END MODULE trprof
