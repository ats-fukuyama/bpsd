! fpcaldrrem.f90

MODULE fpcaldrrem

  ! *******************************************
  !     Radial transport due to EM turbulence
  ! *******************************************

  PRIVATE
  PUBLIC fp_caldrr_em   ! modeld=2
  PUBLIC fp_caldrr_em2  ! modeld=3

CONTAINS

  SUBROUTINE fp_caldrr_em

    USE fpcomm
    USE fpdebug
    USE plprof
    USE libell
    IMPLICIT NONE
    INTEGER:: NSA,NSSA,NR,NP,NTH,ierr
    
    REAL(rkind):: factor_c,rhon,ql,rs,factor_r,pfpl,vfpl,factor_t
    REAL(rkind):: sinml,cosml,raml,rk,factor_th

    !---- Calculation of EM diffusion coefficient ----

    DO NSA=NSASTART,NSAEND
       NSSA=NS_NSA(NSA)

       factor_c=Pi*rd_em_k_ratio**2*rd_em_amp**2

       DO NR=NRSTART,NRENDWG
          RHON=RG(NR)
          CALL pl_qprf(rhon,ql)
          RS=RA*RG(NR)
          factor_t=SQRT(rs/RR**3)/ql

          DO NP=NPSTART,NPEND
             PFPL=PM(NP,NSSA)*PTFP0(NSA)
             VFPL=PFPL/AMFP(NSA)
             factor_t=VFPL*SQRT(rs/RR**3)/ql
             factor_r=VFPL*(AMP*PA(NSSA)/(AEE*PZ(NSSA)*BB))

             DO NTH=1,NTHMAX
                SINML=SINM(NTH)
                COSML=COSM(NTH)
                RAML=SINML**2*(1+EPSRM(NR))
                RK=SQRT((1.D0-RAML*(1.D0-EPSRM(NR)))/(2.D0*EPSRM(NR)*RAML))
                IF(RK.LT.1.D0) THEN
                   factor_th=0.5D0*Pi*SINML**4/(COSML**2*ELLFC(RK,IERR))
                   IF(ierr.NE.0) THEN
                      WRITE(6,*) 'XX fp_cal_rd_em: ELLFC error, ierr=',ierr
                      STOP
                   END IF
                ELSE
                   factor_th=0.D0
                END IF
!                  WRITE(61,'(4I4,5ES12.4)') NTH,NP,NR,NSA, &
!                       SINML,RAML,RK,ELLFC(RK,IERR),factor_th
                DRR(NTH,NP,NR,NSA)=factor_c*factor_r**2*factor_t*factor_th
!                  WRITE(61,'(16X,5ES12.4)') &
!                       factor_c,factor_r,factor_t,factor_th,DRR(NTH,NP,NR,NSA)
                FRR(NTH,NP,NR,NSA)=0.D0
                IF(NTH.EQ.3.AND.NP.EQ.3) &
                     WRITE(6,*) 'DRR: NR,NSA:',NR,NSA, &
                                 factor_c,factor_r,factor_t,factor_th
             END DO ! NTH
          END DO ! NP
       END DO ! NR
    END DO ! NSA

    CALL fp_debug_drr
    RETURN
  END SUBROUTINE fp_caldrr_em

  SUBROUTINE fp_caldrr_em2

    USE fpcomm
    USE fpdebug
      USE plprof
      USE libell
      IMPLICIT NONE
      integer:: NSA,NSSA,NR,NP,NTH,ierr
      REAL(rkind):: factor_c,rhon,ql,rs,factor_r,pfpl,vfpl,factor_t
      REAL(rkind):: sinml,cosml,raml,rk,factor_th

      CALL fp_debug_drr

!---- Calculation of EM diffusion coefficient ----

      DO NSA=NSASTART,NSAEND
         NSSA=NS_NSA(NSA)

         factor_c=Pi*rd_em_k_ratio**2*rd_em_amp**2

         DO NR=NRSTART,NRENDWG
            RHON=RG(NR)
            CALL pl_qprf(rhon,ql)
            RS=RA*RG(NR)
            factor_t=SQRT(rs/RR**3)/ql

            DO NP=NPSTART,NPEND
               PFPL=PM(NP,NSSA)*PTFP0(NSA)
               VFPL=PFPL/AMFP(NSA)
               factor_t=VFPL*SQRT(rs/RR**3)/ql
               factor_r=VFPL*(AMP*PA(NSSA)/(AEE*PZ(NSSA)*BB))

               DO NTH=1,NTHMAX
                  SINML=SINM(NTH)
                  COSML=COSM(NTH)
                  RAML=SINML**2*(1+EPSRM(NR))
                  RK=SQRT((1.D0-RAML*(1.D0-EPSRM(NR)))/(2.D0*EPSRM(NR)*RAML))
                  IF(RK.LT.1.D0) THEN
                     factor_th=0.5D0*Pi*SINML**4/(COSML**2*ELLFC(RK,IERR))
                     IF(ierr.NE.0) THEN
                        WRITE(6,*) 'XX fp_cal_rd_em: ELLFC error, ierr=',ierr
                        STOp
                     END IF
                  ELSE
                     factor_th=0.D0
                  END IF
!                  WRITE(61,'(4I4,5ES12.4)') NTH,NP,NR,NSA, &
!                       SINML,RAML,RK,ELLFC(RK,IERR),factor_th
                  DRR(NTH,NP,NR,NSA)=factor_c*factor_r**2*factor_t*factor_th
!                  WRITE(61,'(16X,5ES12.4)') &
!                       factor_c,factor_r,factor_t,factor_th,DRR(NTH,NP,NR,NSA)
                  FRR(NTH,NP,NR,NSA)=0.D0
                  IF(NTH.EQ.3.AND.NP.EQ.3) &
                     WRITE(6,*) 'DRR: NR,NSA:',NR,NSA,factor_c,factor_r,factor_t,factor_th
               END DO ! NTH
            END DO ! NP
         END DO ! NR
      END DO ! NSA

      CALL fp_debug_drr
      RETURN
    END SUBROUTINE fp_caldrr_em2

  END MODULE fpcaldrrem
