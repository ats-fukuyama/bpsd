! bpsd_libfio.f90

MODULE bpsd_libfio

  PRIVATE
  PUBLIC fropen,fwopen

CONTAINS

!     ***** OPEN FILE FOR READ *****

      SUBROUTINE FROPEN(NFL,KNAMFL,MODEF,MODEP,KPR,IERR,CONVERT)

!     INPUT:
!        NFL    : FILE DEVICE NUMBER
!        KNAMFL : FILE NAME
!        MODEF  : 0 : UNFORMATTED
!                 1 : FORMATTED
!        MODEP  : 0 : WITHOUT PROMPT
!                 1 : WITH FILE NAME INPUT
!        KPR    : PROMPT
!        CONVERT: OPTIONAL: 'BIG_ENDIAN','LITTLE_ENDIAN','NATIVE'

!     OUTPUT:
!        IERR   : ERROR CODE
!                 0 : NO ERROR
!                 1 : EOF IN READ FILE NAME
!                 2 : CANCEL WITH BLANK FILE NAME
!                 3 : UNDEFINED MODEP
!                 4 : FILE NAME ERROR
!                 5 : UNDEFINED MODEF
!                 6 : OLD FILE OPEN ERROR
!                 7 : OLD FILE NOT FOUND
!                 8 : EMPTY FILE NAME

      IMPLICIT NONE

      INTEGER(4),       INTENT(IN)   :: NFL, MODEF, MODEP
      INTEGER(4),       INTENT(OUT)  :: IERR
      CHARACTER(LEN=*), INTENT(IN)   :: KNAMFL
      CHARACTER(LEN=*), INTENT(IN)   :: KPR
      CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: CONVERT
      CHARACTER(LEN=256):: CONVERT_
      INTEGER(4)        :: IST
      CHARACTER(LEN=80) :: KNAM
      LOGICAL           :: LEX

      IF(PRESENT(CONVERT)) THEN
         CONVERT_=CONVERT
      ELSE
         CONVERT_="native"
      END IF

      KNAM=KNAMFL
      IF(MODEP.EQ.0) THEN
         IF(LEN_TRIM(KNAM).EQ.0) GOTO  9008
      ELSEIF(MODEP.EQ.1) THEN
    1    WRITE(6,*) '#',KPR,'> INPUT : LOAD FILE NAME : ',TRIM(KNAM)
         READ(5,'(A80)',ERR=1,END=9001) KNAM
         IF(KNAM(1:2).EQ.'  ') GOTO 9002
      ELSE
         WRITE(6,*) 'XX FROPEN: UNKNOWN MODEP : MODEP=',MODEP
         GOTO 9003
      ENDIF

      INQUIRE(FILE=KNAM,EXIST=LEX,ERR=9004)
      IF(LEX) THEN
         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAM,IOSTAT=IST,STATUS='OLD',ERR=20, &
     &           FORM='UNFORMATTED',CONVERT=CONVERT_)
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAM,IOSTAT=IST,STATUS='OLD',ERR=20, &
     &           FORM='FORMATTED',CONVERT=CONVERT_)
         ELSE
            WRITE(6,*) 'XX FROPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# OLD FILE (',TRIM(KNAM),') IS ASSIGNED FOR INPUT.'
         GOTO 9000

   20    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ELSE
         WRITE(6,*) 'XX FILE (',TRIM(KNAM),') NOT FOUND'
         GOTO 9007
      ENDIF

 9000 IERR=0
      RETURN

 9001 IERR=1
      RETURN

 9002 IERR=2
      RETURN

 9003 IERR=3
      RETURN

 9004 IERR=4
      RETURN

 9005 IERR=5
      RETURN

 9006 IERR=6
      RETURN

 9007 IERR=7
      RETURN

 9008 IERR=8
      RETURN
      END SUBROUTINE FROPEN

!     ***** OPEN FILE FOR WRITE *****

      SUBROUTINE FWOPEN(NFL,KNAMFL,MODEF,MODEP,KPR,IERR)

!     INPUT:
!        NFL    : FILE DEVICE NUMBER
!        KNAMFL : FILE NAME
!        MODEF  : 0 : UNFORMATTED
!                 1 : FORMATTED
!        MODEP  : 0 : WITHOUT PROMPT, ALWAYS OVERWRITE
!                 1 : WITHOUT PROMPT, CONFIRM, IF FILE EXISTS
!                 2 : WITHOUT PROMPT, ASK NEW NAME, IF FILE EXISTS
!                 3 : WITHOUT PROMPT, ERROR, IF FILE EXISTS
!                 4 : WITH FILE NAME INPUT, ALWAYS OVERWRITE
!                 5 : WITH FILE NAME INPUT, CONFIRM, IF FILE EXISTS
!                 6 : WITH FILE NAME INPUT, ASK NEW NAME, IF FILE EXISTS
!                 7 : WITH FILE NAME INPUT, ERROR, IF FILE EXISTS

!     OUTPUT:
!        IERR   : ERROR CODE
!                 0 : NO ERROR
!                 1 : EOF IN READ FILE NAME
!                 2 : CANCEL WITH BLANK FILE NAME
!                 3 : UNDEFINED MODEP
!                 4 : FILE NAME ERROR
!                 5 : UNDEFINED MODEF
!                 6 : OLD FILE OPEN ERROR
!                 7 : OLD FILE EXISTS
!                 8 : EMPTY FILE NAME

      IMPLICIT NONE

      INTEGER(4),       INTENT(IN)   :: NFL, MODEF, MODEP
      INTEGER(4),       INTENT(OUT)  :: IERR
      CHARACTER(LEN=*), INTENT(IN)   :: KNAMFL
      CHARACTER(LEN=*), INTENT(IN)   :: KPR
      INTEGER(4)        :: MODEPI, MODEPII, IST
      CHARACTER(LEN=80) :: KNAM
      CHARACTER(LEN=1)  :: KID
      LOGICAL           :: LEX

      MODEPI=MODEP

      KNAM=KNAMFL

 1000 IF(MODEPI.LE.3) THEN
         IF(LEN_TRIM(KNAM).EQ.0) GOTO 9008
      ELSE
    1    WRITE(6,*) '#',KPR,'> INPUT : SAVE FILE NAME : ',TRIM(KNAM)
         READ(5,'(A80)',ERR=1,END=9001) KNAM
         IF(KNAM(1:2).EQ.'  ') GOTO 9002
      ENDIF

      INQUIRE(FILE=KNAM,EXIST=LEX,ERR=9004)

      IF(LEX) THEN
         MODEPII=MOD(MODEPI,4)
         IF(MODEPII.EQ.0) THEN
            WRITE(6,*) '# OLD FILE (',TRIM(KNAM), &
     &                 ') WILL BE OVERWRITTEN'
         ELSEIF(MODEPII.EQ.1) THEN
    3       WRITE(6,*) '# OLD FILE (',TRIM(KNAM), &
     &                 ') IS GOING TO BE OVERWRITTEN'
            WRITE(6,*) '  ARE YOU SURE ? (Y/N)'
            READ(5,'(A1)',ERR=3,END=9001) KID
            CALL GUCPTL(KID)
            IF(KID.EQ.'N') GOTO 9007
         ELSEIF(MODEPII.EQ.2) THEN
            MODEPI=1
            GOTO 1000
         ELSEIF(MODEPII.EQ.3) THEN
            WRITE(6,*) 'XX FWOPEN: FILE ALREADY EXISTS.'
            GOTO 9007
         ELSE
            WRITE(6,*) 'XX FWOPEN: UNKNOWN MODEP : MODEP=',MODEP
            GOTO 9003
         ENDIF

         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAM,IOSTAT=IST,STATUS='OLD',ERR=10, &
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAM,IOSTAT=IST,STATUS='OLD',ERR=10, &
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FEOPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# OLD FILE (',TRIM(KNAM), &
     &                 ') IS ASSIGNED FOR OUTPUT.'
         GOTO 9000

   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ELSE
         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAM,IOSTAT=IST,STATUS='NEW',ERR=20, &
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAM,IOSTAT=IST,STATUS='NEW',ERR=20, &
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FEOPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# NEW FILE (',TRIM(KNAM),') IS CREATED FOR OUTPUT.'
         GOTO 9000

   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ENDIF

 9000 IERR=0
      RETURN

 9001 IERR=1
      RETURN

 9002 IERR=2
      RETURN

 9003 IERR=3
      RETURN

 9004 IERR=4
      RETURN

 9005 IERR=5
      RETURN

 9006 IERR=6
      RETURN

 9007 IERR=7
      RETURN

 9008 IERR=8
      RETURN
    END SUBROUTINE FWOPEN
  END MODULE bpsd_libfio
