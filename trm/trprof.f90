! trprof.f90

MODULE trprof

  PRIVATE
  PUBLIC tr_prof
  PUBLIC tr_prof_impurity
  PUBLIC tr_prof_current

CONTAINS

!     ***********************************************************

!           SET INITIAL PROFILE

!     ***********************************************************

  SUBROUTINE tr_prof

      USE trcomm
      USE trfixed
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
      END SELECT

      DO nr=1,nrmax
         SELECT CASE(model_prof)
         CASE(11)
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
            
         CASE default
            ! --- set default profiles ---

            DO ns=1,nsmax
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1(NS))**PROFN2(NS)
               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1(NS))**PROFT2(NS)
               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1(NS))**PROFU2(NS)
               RU(NR,NS) = (PU(NS)-PUS(NS))*PROF+PUS(NS)
            END DO
         END SELECT
         
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

      SELECT CASE(model_nfixed)
      CASE(1)
         CALL tr_prep_nfixed
         IF(time_nfixed(1).LE.0.D0) THEN
            DO nr=1,nrmax
               CALL tr_prof_nfixed(rm(nr),t,rn(nr,1))
               DO ns=2,nsmax
                  rn(nr,ns)=pn(ns)/(pz(ns)*pn(1))*rn(nr,1)
                  IF(ns.EQ.2) WRITE(6,'(I6,3ES12.4)') nr,rm(nr),rn(nr,1),rn(nr,2)
               END DO
            END DO
            CALL tr_prof_nfixed(1.D0,t,pns(1))
            pnss(1)=pns(1)
            DO ns=2,nsmax
               pns(ns)=pn(ns)/(pz(ns)*pn(1))*pns(1)
               pnss(ns)=pns(ns)
            END DO
         END IF
      CASE(2)
         CALL tr_prep_nfixed
         IF(time_nfixed(1).LE.0.D0) THEN
            DO nr=1,nrmax
               IF((rm(nr).GE.rho_min_nfixed).AND. &
                  (rm(nr).LE.rho_max_nfixed)) THEN
                  CALL tr_prof_nfixed(rm(nr),t,rn(nr,1))
                  DO ns=2,nsmax
                     rn(nr,ns)=pn(ns)/(pz(ns)*pn(1))*rn(nr,1)
                  END DO
               END IF
            END DO
            CALL tr_prof_nfixed(1.D0,t,pns(1))
            pnss(1)=pns(1)
            DO ns=2,nsmax
               pns(ns)=pn(ns)/(pz(ns)*pn(1))*pns(1)
               pnss(ns)=pns(1)
            END DO
         END IF
      END SELECT
      SELECT CASE(model_tfixed)
      CASE(1)
         CALL tr_prep_tfixed
         IF(time_tfixed(1).LE.0.D0) THEN
            DO nr=1,nrmax
               CALL tr_prof_tfixed(rm(nr),t,rt(nr,1))
               DO ns=2,nsmax
                  rt(nr,ns)=rt(nr,1)
               END DO
            END DO
            CALL tr_prof_tfixed(1.D0,t,pts(1))
            DO ns=2,nsmax
               pts(ns)=pts(1)
            END DO
         END IF
      CASE(2)
         CALL tr_prep_tfixed
         IF(time_tfixed(1).LE.0.D0) THEN
            DO nr=1,nrmax
               IF((rm(nr).GE.rho_min_tfixed).AND. &
                  (rm(nr).LE.rho_max_tfixed)) THEN
                  CALL tr_prof_tfixed(rm(nr),t,rt(nr,1))
                  DO ns=2,nsmax
                     rt(nr,ns)=rt(nr,1)
                  END DO
               END IF
            END DO
            CALL tr_prof_tfixed(1.D0,t,pts(1))
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

  END MODULE trprof
