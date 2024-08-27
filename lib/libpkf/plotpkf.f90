! plotpkf.f90

PROGRAM plotpkf
  USE task_kinds
  USE task_constants
  USE libpdkf
  USE libpdkf2
  USE libpgkf
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
  
  CALL GSOPEN

  nxmax=101
  nymax=101
  xmin=0.001D0
  xmax=1000.D0
  nxmax=2
  xmax=0.006D0
  xmax=0.007D0
  alpha=0.1D0
  alpha_min=0.1D0
  alpha_max=1000.D0
!  beta=1.0D0
  beta=300.0D0
  beta_min=0.001D0
  beta_max=1000.D0
  rky=0.D0
  n1=1
  ntau=0
  nf=0
  nn=0
  nm=0
  mode=0

1 CONTINUE
  
  WRITE(6,'(A)') '## mode: 0,1,2:pdkf 10:pdkf2 20:pgkf 9:QUIT'
  READ(5,*,ERR=1,END=9000) mode

  SELECT CASE(mode)
  CASE(0) ! pdkf
     WRITE(6,'(A)') '## n1=1:2'
10   CONTINUE
     WRITE(6,'(A/A,3ES12.4,I4,2ES12.4,I6)') &
          '## INPUT : alpha,beta,rky,n1,xmin,xmax,nxmax:', &
          '##   ', alpha,beta,rky,n1,xmin,xmax,nxmax
     READ(5,*,END=1) alpha,beta,rky,n1,xmin,xmax,nxmax

     IF(nxmax.LT.1) GO TO 1
     
     IF(n1.LT.1.OR.N1.GT.2) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pkrf: n1 must be in [1:2]: n1=',n1
        GO TO 10
     END IF

     ALLOCATE(xa(nxmax),fa(nxmax,2))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
        cf_=pdkf(xa(nx),alpha,beta,rky,5,n1)
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

  CASE(1) ! pdkf 2D
     WRITE(6,'(A)') '## n1=1:2, alpha>0, x>0'
11   CONTINUE
     WRITE(6,'(A/A,2ES12.4,I4,2ES12.4,I6,2ES12.4,I6)') &
          '## INPUT : beta,rky,n1,xmin,xmax,nxmax,alpha_min,alpha_max,nymax:',&
          '##   ', beta,rky,n1,xmin,xmax,nxmax,alpha_min,alpha_max,nymax
     READ(5,*,ERR=11,END=1) &
          beta,rky,n1,xmin,xmax,nxmax,alpha_min,alpha_max,nymax

     IF(nxmax.LT.1.OR.nymax.LT.1) GO TO 11
     
     IF(n1.LT.1.OR.N1.GT.2) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pdkf: n1 must be in [1:2]: n1=',n1
        GO TO 11
     END IF

     ALLOCATE(xa(nxmax),ya(nymax),f2d(nxmax,nymax,2))

     dx=(log10(xmax)-log10(xmin))/(nxmax-1)
     dy=(log10(alpha_max)-log10(alpha_min))/(nymax-1)
     DO nx=1,nxmax
        xa(nx)=log10(xmin)+dx*(nx-1)
        DO ny=1,nymax
           ya(ny)=log10(alpha_min)+dy*(ny-1)
           cf_=pdkf(10.D0**xa(nx),10.D0**ya(ny),beta,rky,5,n1)
           f2d(nx,ny,1)=REAL(cf_)
           f2d(nx,ny,2)=AIMAG(cf_)
        END DO
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF vs x, alpha: n1=',n1,'@'
     CALL grd2d(0,xa,ya,f2d,nxmax,nxmax,nymax,title,3,ASPECT=1.D0)
     CALL PAGEE
     DEALLOCATE(xa,ya,f2d)
     GO TO 11

  CASE(2) ! pdkf 2D
     WRITE(6,'(A)') '## n1=1:2, beta>0, x>0'
12   CONTINUE
     WRITE(6,'(A/A,2ES12.4,I4,2ES12.4,I6,2ES12.4,I6)') &
          '## INPUT : alpha,rky,n1,xmin,xmax,nxmax,beta_min,beta_max,nymax:',&
          '##   ', alpha,rky,n1,xmin,xmax,nxmax,beta_min,beta_max,nymax
     READ(5,*,ERR=11,END=1) &
          alpha,rky,n1,xmin,xmax,nxmax,beta_min,beta_max,nymax

     IF(nxmax.LT.1.OR.nymax.LT.1) GO TO 11
     
     IF(n1.LT.1.OR.N1.GT.2) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pdkf: n1 must be in [1:2]: n1=',n1
        GO TO 11
     END IF

     ALLOCATE(xa(nxmax),ya(nymax),f2d(nxmax,nymax,2))

     dx=(log10(xmax)-log10(xmin))/(nxmax-1)
     dy=(log10(beta_max)-log10(beta_min))/(nymax-1)
     DO nx=1,nxmax
        xa(nx)=log10(xmin)+dx*(nx-1)
        DO ny=1,nymax
           ya(ny)=log10(beta_min)+dy*(ny-1)
           cf_=pdkf(10.D0**xa(nx),alpha,10.D0**ya(ny),rky,5,n1)
           f2d(nx,ny,1)=REAL(cf_)
           f2d(nx,ny,2)=AIMAG(cf_)
        END DO
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF real vs x, beta: n1=',n1,'@'
     CALL grd2d(0,xa,ya,f2d(:,:,1),nxmax,nxmax,nymax,title,3, &
          ASPECT=1.D0,NLMAX=12)
     CALL PAGEE
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF imag vs x, beta: n1=',n1,'@'
     CALL grd2d(0,xa,ya,f2d(:,:,2),nxmax,nxmax,nymax,title,3, &
          ASPECT=1.D0,NLMAX=12)
     CALL PAGEE
     DEALLOCATE(xa,ya,f2d)
     GO TO 12

  CASE(10) ! pgkf
     WRITE(6,'(A)') '## nf=0:4, xmax<=10.D0'
20   CONTINUE
     WRITE(6,'(A/A,3I4,2ES12.4,I6)') &
          '## INPUT : nf,nm,nn,xmin,xmax,nxmax:', &
          '##   ', nf,nm,nn,xmin,xmax,nxmax
     READ(5,*,END=1) nf,nm,nn,xmin,xmax,nxmax
     IF(xmax.gt.10.D0) xmax=10.D0

     IF(nxmax.LT.1) GO TO 1
     
     IF(nf.LT.0.OR.nf.GT.4) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pkrf: nf must be in [0:4]: n4=',nf
        GO TO 20
     END IF

     ALLOCATE(xa(nxmax),fa(nxmax,2))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
        cf_=pgkf(xa(nx),nf,nn,nm)
        fa(nx,1)=REAL(cf_)
        fa(nx,2)=AIMAG(cf_)
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF vs x: n1=',n1,'@'
     CALL grd1d(0,xa,fa,nxmax,nxmax,2,title)
     CALL PAGEE
     DEALLOCATE(xa,fa)
     GO TO 20
     
  CASE(30) ! pdkf2
     WRITE(6,'(A)') '## n1=1:2'
30   CONTINUE
     WRITE(6,'(A/A,2ES12.4,I4,2ES12.4,I6)') &
          '## INPUT : alpha,beta,ntau,xmin,xmax,nxmax:', &
          '##   ', alpha,beta,ntau,xmin,xmax,nxmax
     READ(5,*,END=1) alpha,beta,ntau,xmin,xmax,nxmax

     IF(nxmax.LT.1) GO TO 1
     
     ALLOCATE(xa(nxmax),fa(nxmax,2))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
        cf_=pdkf2(xa(nx),alpha,beta,ntau)
        fa(nx,1)=REAL(cf_)
        fa(nx,2)=AIMAG(cf_)
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i2,A)') '@PDKFs vs x: ntau=',ntau,'@'
     CALL grd1d(0,xa,fa,nxmax,nxmax,2,title)
     CALL PAGEE
     DEALLOCATE(xa,fa)
     GO TO 30

  CASE(9)
     GO TO 9000
  END SELECT
  GO TO 1

9000 CALL GSCLOS
  STOP
END PROGRAM plotpkf
