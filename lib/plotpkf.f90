! plotpkf.f90

PROGRAM plotpkf
  USE task_kinds
  USE task_constants
  IMPLICIT NONE
  INTEGER:: mode=0

  CALL GSOPEN

1 CONTINUE
  
  WRITE(6,'(A)') '## mode: 1:pdkf 2:pdkf2 3:pgkf 9:QUIT'
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
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nxmax,nymax
    REAL(rkind):: xmin,xmax
    COMPLEX(rkind):: cf_
    REAL(rkind),ALLOCATABLE:: xa(:),ya(:),fa(:,:),f2d(:,:,:)
    REAL(rkind):: alpha,beta,rky,dx,dy,alpha_min,alpha_max
    REAL(rkind):: beta_min,beta_max
    INTEGER:: mode,n1,nn,nf,nm,nx,ntau,ny
    CHARACTER(LEN=80):: title
  
    nxmax=101
    nymax=101
    xmin=0.01D0
    xmax=1.D0
    alpha=0.1D0
    alpha_min=0.1D0
    alpha_max=1000.D0
    !  beta=1.0D0
    beta=15.0D0
    beta_min=0.001D0
    beta_max=1000.D0
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

     ALLOCATE(xa(nxmax),fa(nxmax,2))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
!        cf_=pdkf_eul(alpha*xa(nx),beta,1.D0,n1)
        cf_=pdkf_ft(xa(nx),beta,1.D0,n1)
        WRITE(21,'(A,3ES12.4)') 'x,cf        =',xa(nx),cf_
        fa(nx,1)=REAL(cf_)
        fa(nx,2)=AIMAG(cf_)
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF vs x: n1=',n1,'@'
     CALL grd1d(0,xa,fa,nxmax,nxmax,2,title)
     CALL PAGEE
     DEALLOCATE(xa,fa)
     GO TO 10


9000 CONTINUE
     RETURN
   END SUBROUTINE plot_pdkf
END PROGRAM plotpkf
