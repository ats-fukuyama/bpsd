! bpsd_libchar.f90

MODULE bpsd_libchar

  PRIVATE
  PUBLIC toupper
  PUBLIC tolower
  
CONTAINS

!***************************************************************
!
!   Convert Strings to Upper Case
!
!***************************************************************

SUBROUTINE TOUPPER(KTEXT)

  implicit none
  character(len=*), INTENT(INOUT) ::  KTEXT

  INTEGER :: NCHAR, I, ID

  NCHAR = LEN(KTEXT)
  DO I = 1, NCHAR
     ID=IACHAR(KTEXT(I:I))
     IF(ID >= 97 .AND. ID <= 122) ID = ID - 32
     KTEXT(I:I)=ACHAR(ID)
  END DO

  RETURN
END SUBROUTINE TOUPPER

!***************************************************************
!
!   Convert Strings to Lower Case
!
!***************************************************************

SUBROUTINE TOLOWER(KTEXT)

  implicit none
  character(len=*), INTENT(INOUT) ::  KTEXT

  INTEGER :: NCHAR, I, ID

  NCHAR = LEN(KTEXT)
  DO I = 1, NCHAR
     ID=IACHAR(KTEXT(I:I))
     IF(ID >= 65 .AND. ID <= 90) ID = ID + 32
     KTEXT(I:I)=ACHAR(ID)
  END DO

  RETURN
END SUBROUTINE TOLOWER
END MODULE bpsd_libchar
