! testpkf.f90

PROGRAM testpkf
  USE task_kinds
  USE task_constants
  USE libpdkf
  USE libpdkf2
  USE libpgkf
  USE libgrf
  IMPLICIT NONE

  COMPLEX(rkind):: cf_
  REAL(rkind):: x,f
  REAL(rkind):: alpha,beta,rky
  INTEGER:: mode,ntype,n1,m,nx,ny
  CHARACTER(LEN=80):: title
  INTEGER:: nxmax,nymax
  REAL(rkind):: xmin,xmax,dx,dy
  REAL(rkind),ALLOCATABLE:: xa(:),ya(:),fa(:,:),f2d(:,:,:)
  REAL(rkind):: alpha_min,alpha_max
  REAL(rkind):: beta_min,beta_max
  
  CALL GSOPEN

  mode=0
  n1=1
  m=5
  ntype=1
  alpha=300.D0
  beta=0.1D0
  rky=0.D0
  x=0.005
  xmin=0.001D0
  xmax=0.01D0
  nxmax=101
  alpha_min=100.D0
  alpha_max=1000.D0
  nymax=101
  
1 CONTINUE

  WRITE(6,'(A)') '## ntype: 1:PDKF 2:PDKF2 3:PGKF 9:end'
  READ(5,*,ERR=1,END=9000) ntype

  SELECT CASE(ntype)
  CASE(1) !  PDKF
     WRITE(6,'(A)') '## PDKF:'
11   CONTINUE

     WRITE(6,'(A)') '## mode: 0:point 1:1D 2:1D+ 3:2D 9:end'
     READ(5,*,ERR=11,END=1) mode

     SELECT CASE(mode)
     CASE(0)
        WRITE(6,'(A)') '## point:'
21      CONTINUE
        WRITE(6,'(A)') '## n1,m,x,alpha,beta,rky'
        WRITE(6,'(A,2I4,4ES12.4)') '##', n1,m,x,alpha,beta,rky
        READ(5,*,ERR=21,END=11) n1,m,x,alpha,beta,rky
        IF(n1.EQ.0) GO TO 11
        IF(n1.EQ.1.OR.n1.EQ.2) THEN
           cf_=pdkf(x,alpha,beta,rky,m,n1)
           WRITE(6,'(A,3ES12.4)') 'xa,cf       =',x,cf_
        ELSE
           WRITE(6,*) 'XX n1 must be 1 or 2.'
        END IF
        GO TO 21
        
     CASE(1)
        WRITE(6,'(A)') '## 1D single alpha scan:'
31      CONTINUE

        WRITE(6,'(A/A,2I4,5ES12.4,I6)') &
          '## INPUT : n1,m,alpha,beta,rky,xmin,xmax,nxmax:', &
          '## ', n1,m,alpha,beta,rky,xmin,xmax,nxmax
        READ(5,*,ERR=31,END=11) n1,m,alpha,beta,rky,xmin,xmax,nxmax

        IF(nxmax.LT.1) GO TO 31
     
        IF(n1.EQ.1.OR.n1.EQ.2) THEN
           ALLOCATE(xa(nxmax),fa(nxmax,2))

           dx=(xmax-xmin)/(nxmax-1)
           DO nx=1,nxmax
              xa(nx)=xmin+dx*(nx-1)
              cf_=pdkf(xa(nx),alpha,beta,rky,m,n1)
              fa(nx,1)=REAL(cf_)
              fa(nx,2)=AIMAG(cf_)
           END DO
  
           CALL PAGES
           WRITE(title,'(A,i1,A)') '@PDKF vs x: n1=',n1,'@'
           CALL grd1d(0,xa,fa,nxmax,nxmax,2,title)
           CALL PAGEE
           DEALLOCATE(xa,fa)
        ELSE
           WRITE(6,*) 'XX n1 must be 1 or 2.'
        END IF
        GO TO 31
        
     CASE(2)
        WRITE(6,'(A)') '## 1D multi alpha scan:'
41      CONTINUE
        
        WRITE(6,'(A)') &
          '## INPUT : n1,m,nymax,nxmax/alpha_min,alpha_max,beta,rky,xmin,xmax:'
        WRITE(6,'(A,4I6/A,6ES12.4)') &
             '## ',n1,m,nymax,nxmax, &
             '   ',alpha_min,alpha_max,beta,rky,xmin,xmax
        READ(5,*,ERR=41,END=11) &
             n1,m,nymax,nxmax,alpha_min,alpha_max,beta,rky,xmin,xmax

        IF(nxmax.LT.1.OR.nymax.LT.1) GO TO 41
     
        IF(n1.EQ.1.OR.n1.EQ.2) THEN
           ALLOCATE(xa(nxmax),ya(nymax),f2d(nxmax,nymax,2))

           dx=(log10(xmax)-log10(xmin))/(nxmax-1)
           dy=(log10(alpha_max)-log10(alpha_min))/(nymax-1)
           DO nx=1,nxmax
              xa(nx)=log10(xmin)+dx*(nx-1)
              DO ny=1,nymax
                 ya(ny)=log10(alpha_min)+dy*(ny-1)
                 cf_=pdkf(10.D0**xa(nx),10.D0**ya(ny),beta,rky,m,n1)
                 f2d(nx,ny,1)=REAL(cf_)
                 f2d(nx,ny,2)=AIMAG(cf_)
              END DO
           END DO
  
           CALL PAGES
           WRITE(title,'(A,i1,A)') '@PDKF real vs x, alpha: n1=',n1,'@'
           CALL grd2d(0,xa,ya,f2d(:,:,1),nxmax,nxmax,nymax,title,3, &
                ASPECT=1.D0,NLMAX=12)
           CALL PAGEE
           CALL PAGES
           WRITE(title,'(A,i1,A)') '@PDKF imag vs x, beta: n1=',n1,'@'
           CALL grd2d(0,xa,ya,f2d(:,:,2),nxmax,nxmax,nymax,title,3, &
                ASPECT=1.D0,NLMAX=12)
           CALL PAGEE
           DEALLOCATE(xa,ya,f2d)
        ELSE
           WRITE(6,*) 'XX n1 must be 1 or 2.'
        END IF
     CASE(9)
        GOTO 21
     CASE DEFAULT
        WRITE(6,'(A,I4)') 'XX undefined mode=',mode
     END SELECT
     GO TO 41
           
  CASE(9)
     GO TO 9000
  CASE DEFAULT
     WRITE(6,'(A,I4)') 'XX undefined mode=',mode
  END SELECT
  GO TO 1

9000 CONTINUE
  CLOSE(21)
  CALL GSCLOS
END PROGRAM
     
