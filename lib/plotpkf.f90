! plotpkf.f90

PROGRAM plotpkf
  USE task_kinds
  USE task_constants
  IMPLICIT NONE
  INTEGER:: mode=0

  CALL GSOPEN

1 CONTINUE
  
  WRITE(6,'(A)') '## mode_pkf: 1:pdkf 2:pdkf2 3:pgkf 9:QUIT'
  READ(5,*,ERR=1,END=9000) mode

  SELECT CASE(mode)
  CASE(1)
     CALL plot_pdkf
  CASE(2)
!     CALL plot_pdkf2
  CASE(3)
!     CALL plot_pgkf
  CASE(9)
     GO TO 9000
  CASE DEFAULT
     WRITE(6,*) 'XX unknown mode[1,2,3,9]: mode=',mode
  END SELECT
  GO TO 1

9000 CONTINUE
  CALL GSCLOS
  STOP

CONTAINS

  ! *** plot function pdkf ***
  
  SUBROUTINE plot_pdkf

    USE libpdkf
    USE libpdkf_eul
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nxmax,nymax
    REAL(rkind):: xmin,xmax
    COMPLEX(rkind):: cf_
    REAL(rkind),ALLOCATABLE:: xa(:),ya(:),fa(:,:),f2d(:,:,:)
    INTEGER:: mode
  
1   CONTINUE
  
    WRITE(6,'(A)') '## mode_pdkf: 1:pdkf_eul 2:pdkf_ft 3:pdkf_hift  9:QUIT'
    READ(5,*,ERR=1,END=9000) mode

    SELECT CASE(mode)
    CASE(1)
       CALL plot_pdkf_eul
    CASE(2)
       CALL plot_pdkf_ft
    CASE(3)
       CALL plot_pdkf_hift
    CASE(9)
       GO TO 9000
    CASE DEFAULT
       WRITE(6,*) 'XX unknown mode[1,2,3,9]: mode=',mode
    END SELECT
    GO TO 1
    
9000 CONTINUE
    RETURN
  END SUBROUTINE plot_pdkf

  SUBROUTINE plot_pdkf_eul

    USE libpdkf_eul
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nxmax,nymax,nzmax    ! x=xi, y=eta, z=nu
    REAL(rkind):: xmin,xmax,ymin,ymax,zmin,zmax
    COMPLEX(rkind):: cf_
    REAL(rkind),ALLOCATABLE:: xa(:),ya(:),za(:),fa(:,:),f2d(:,:,:)
    REAL(rkind):: dx,dy,dz
    INTEGER:: mode,n1,nn,nf,nm,nx,ntau,ny,i
    CHARACTER(LEN=80):: title
  
    nxmax=101
    nymax=1
    
    xmin=0.01D0
    xmax=1.D0
    rky=0.D0
    n1=1
    ntau=0
    nf=0
    nn=0
    nm=0
    
    WRITE(6,'(A)') '## n1=1:2'
10  CONTINUE
    WRITE(6,'(A/A,I4,3ES12.4,I4,2ES12.4,I6)') &
          '## INPUT : alpha,beta,rky,n1,xmin,xmax,nxmax:', &
          '##   ', mode,alpha,beta,rky,n1,xmin,xmax,nxmax
     READ(5,*,END=9000,ERR=10) mode,alpha,beta,rky,n1,xmin,xmax,nxmax

     IF(nxmax.LT.1) GO TO 10
     
     IF(n1.LT.1.OR.N1.GT.2) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pkrf: n1 must be in [1:2]: n1=',n1
        GO TO 10
     END IF

     ALLOCATE(xa(nxmax),fa(nxmax,4))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
        cf_=pdkf_eul(alpha*xa(nx),beta,1.D0,n1)
        fa(nx,1)=REAL(cf_)
        fa(nx,2)=AIMAG(cf_)
        cf_=pdkf_ft(xa(nx),beta,1.D0,n1)
        fa(nx,3)=REAL(cf_)
        fa(nx,4)=AIMAG(cf_)
        WRITE(6,'(A,5ES12.4)') 'x,f=',xa(nx),(fa(nx,i),i=1,4)
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF vs x: n1=',n1,'@'
     CALL grd1d(0,xa,fa,nxmax,nxmax,4,title)
     CALL PAGEE
     DEALLOCATE(xa,fa)
     GO TO 10


9000 CONTINUE
     RETURN
   END SUBROUTINE plot_pdkf
END PROGRAM plotpkf
