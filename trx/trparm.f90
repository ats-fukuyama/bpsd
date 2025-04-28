! trparm.f90

MODULE trparm

  USE trcomm_parm
  NAMELIST /TR/ &
       RR,RA,RB,RKAP,RDLT,BB,RIPS,RIPE,RHOA, &
       NSMAX,NSZMAX,NSNMAX, &
       PM,PZ,PN,PNS,PT,PTS,PU,PUS, &
       MDLIMP,PNC,PNFE,PNNU,PNNUS, &
       PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
       PROFNU1,PROFNU2,PROFJ1,PROFJ2,ALP, &
       model_prof,knam_prof,knam_ngiven,knam_tgiven, &
       AD0,AV0,CNP,CNH,CDP,CDH,CNN,CDW, &
       MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF, &
       MDLWLD,MDLDW,MDLCD05, &
       CHP,CK0,CK1,CWEB,CALF,CKALFA,CKBETA,CKGUMA, &
       DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST, &
       EPSLTR,LMAXTR,TPRST,MDLST,IZERO, &
       nnfmax,model_nnf, &
       NNBMAX,NECMAX,NLHMAX,NICMAX,NPELMAX,NPSCMAX, &
       model_nnb,MDLEC,MDLLH,MDLIC,MDLPEL,MDLPSC, &
       PNBIN,PNBR0,PNBRW,PNBCD,PNBVY,PNBVW,PNBENG,PNBRTG, &
       ns_nnb,nrmax_nnb, &
       PECIN,PECR0,PECRW,PECCD,PECTOE,PECNPR, &
       PLHIN,PLHR0,PLHRW,PLHCD,PLHTOE,PLHNPR, &
       PICIN,PICR0,PICRW,PICCD,PICTOE,PICNPR, &
       PELIN,PELR0,PELRW,PELRAD,PELVEL,PELTIM, &
       pellet_time_start,pellet_time_interval,PELPAT, &
       number_of_pellet_repeat, &
       PSCIN,PSCR0,PSCRW,NSPSC,PSCIN_MAX, &
       MDLCD,PBSCD, &
       MDLPR,SYNCABS,SYNCSELF, &
       MDLEDGE,CSPRS, &
       MODELG,NTEQIT,MODEP,MDLJQ,MDLFLX, &
       MDLNI, &
       MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI, &
       MDLER,MDLNCL,NSLMAX, &
       MDLELM,ELMWID,ELMDUR,ELMNRD,ELMTRD,ELMENH, &
       MDLTC,MDLPCK,model_ngiven,model_tgiven, &
       KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS
  PRIVATE
  PUBLIC tr_parm
  PUBLIC tr_nlin
  PUBLIC tr_save_parm
  PUBLIC tr_load_parm

CONTAINS

!     ***********************************************************

!           PARAMETER INPUT

!     ***********************************************************

  SUBROUTINE tr_parm(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

      USE trcomm
      USE libkio
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER,INTENT(OUT):: IERR

    1 CALL TASK_PARM(MODE,'TR',KIN,tr_nlin,trplst,IERR)
      IF(IERR.NE.0) RETURN

      CALL TRCHEK(IERR)
      NTMAX_SAVE=NTMAX
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
    END SUBROUTINE tr_parm

!     ****** INPUT NAMELIST ******

    SUBROUTINE tr_nlin(NID,IST,IERR)

      USE trcomm
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NID
      INTEGER,INTENT(OUT):: IST, IERR

      IF(NID.GE.0) THEN
         READ(NID,TR,IOSTAT=IST,ERR=9800,END=9900)
         NTMAX_SAVE=NTMAX
      ENDIF
      IST=0
      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
    END SUBROUTINE tr_nlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE trplst

      WRITE(6,601)
      RETURN

  601 FORMAT(' ','# &TR : RR,RA,RB,RKAP,RDLT,BB,RIPS,RIPE,RHOA'/ &
             ' ',8X,'(PM,PZ,PN,PNS,PT,PTS:NSM)'/ &
             ' ',8X,'PNC,PNFE,PNNU,PNNUS'/ &
             ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/ &
             ' ',8X,'PROFJ1,PROFJ2,ALP'/ &
             ' ',8X,'CK0,CK1,CNP,CNH,CDP,CDH,CNN,CDW,CSPRS'/ &
             ' ',8X,'CWEB,CALF,CKALFA,CKBETA,MDLKNC,MDLTPF'/ &
             ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/ &
             ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/ &
             ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/ &
             ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG'/ &
             ' ',8X,'PNBCD,MDLNB,NSNNB,NRNBMAX'/ &
             ' ',8X,'PECIN,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/ &
             ' ',8X,'PLHIN,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/ &
             ' ',8X,'PICIN,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/ &
             ' ',8X,'PELIN,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/ &
             ' ',8X,'PELTIM,PELPAT'/ &
             ' ',8X,'pellet_time_start,pellet_time_interval'/ &
             ' ',8X,'number_of_pellet_repeat'/ &
             ' ',8X,'MDLPR,SYNCABS,SYNCSELF,MODELG,NTEQIT'/&
             ' ',8X,'MDEDGE,MDLIMP,model_prof,knam_prof,'/ &
             ' ',8X,'MDLXP,MDLUF,MDLNCL,MDLWLD,MDLFLX,MDLER,MDCD05'/ &
             ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/ &
             ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDIR,KUFDEV,KUFDCG'/ &
             ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDLTC,MDLPCK'/ &
             ' ',8X,'KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS,'/ &
             ' ',8X,'MDLPSC,NPSCMAX,PSCIN,PSCR0,PSCRW,NSPSC,PSCIN_MAX,'/ &
             ' ',8X,'knam_ngiven,knam_tgiven'/ &
             ' ',8X,'model_ngiven,model_tgiven')
    END SUBROUTINE trplst

!     ***** CHECK INPUT PARAMETERS *****

    SUBROUTINE trchek(IERR)

      USE trcomm
      IMPLICIT NONE
      INTEGER, INTENT(OUT):: IERR


      IERR=0

      IF(NRMAX.LT.1) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX =',NRMAX
         IERR=1
      ENDIF

      IF(NTMAX.LT.0.OR.NTMAX/NGTSTP.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM*NGTSTP
         IERR=1
      ENDIF

      RETURN
    END SUBROUTINE trchek


!     ****** save NAMELIST tr ******

  SUBROUTINE tr_save_parm(nfc)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nfc

    WRITE(nfc,TR)
    RETURN
  END SUBROUTINE tr_save_parm

!     ****** load NAMELIST tr ******

  SUBROUTINE tr_load_parm(nfc,ierr)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nfc
    INTEGER,INTENT(OUT) :: ierr
    INTEGER:: ist

    CALL tr_nlin(nfc,ist,ierr)
    IF(ierr.EQ.8) GOTO 9800
    IF(ierr.EQ.9) GOTO 9900
    ierr=0
    WRITE(6,'(A,I5)') '## tr_load_parm: namelist file loaded.'
    RETURN

9800 CONTINUE
    ierr=1
    WRITE(6,'(A,I5)') 'XX tr_load_parm: read error: ist=',ist
    RETURN

9900 CONTINUE
    ierr=2
    WRITE(6,'(A)') 'XX tr_load_parm: file open error:'
    RETURN
  END SUBROUTINE tr_load_parm
    
END MODULE trparm
