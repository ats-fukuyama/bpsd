! testsgima.f90

MODULE sigmav_nrl_dt_comm
  USE bpsd_kinds
  REAL(rkind):: T_bulk_dt
END MODULE sigmav_nrl_dt_comm

PROGRAM testsigma
  USE bpsd_kinds
  USE bpsd_constants
  USE libgrf
  IMPLICIT NONE

  REAL(rkind),DIMENSION(10):: &
       tempa=(/ 1.D0, 2.D0, 5.D0, 10.D0, 20.D0, &
               50.D0, 100.D0, 200.D0, 500.D0, 1000.D0 /)
  REAL(rkind),DIMENSION(10):: &
       sva=(/ 5.5D-21, 2.6D-19, 1.3D-17, 1.1D-16, 4.2D-16, &
             8.7D-16, 8.5D-16, 6.3D-16, 3.7D-16, 2.7D-16 /)
  REAL(rkind),DIMENSION(10):: templa,svla
  REAL(rkind),DIMENSION(4,10):: usvla

  REAL(rkind),DIMENSION(10):: x
  REAL(rkind),DIMENSION(10,3):: y
  INTEGER:: i

  CALL GSOPEN
  
  CALL sigmavm_nrl_table_dt
  
  DO i=1,10
     x(i)=LOG10(tempa(i))
     y(i,1)=LOG10(sva(i))
     y(i,2)=LOG10(sigmavm_int_dt(tempa(i)))
     y(i,3)=LOG10(sigmavm_low_dt(tempa(i)))
     WRITE(6,'(I4,4ES12.4)') i,x(i),y(i,1),y(i,2),y(i,3)
  END DO

  CALL PAGES
  CALL grd1d(0,x,y,10,10,3,'@temp vs sigmav@',3)
  CALL PAGEE

  CALL GSCLOS
  STOP

CONTAINS

  ! Crosssection for energy E
  
  FUNCTION sigma_nrl_dt(E)
    IMPLICIT NONE
    REAL(rkind),PARAMETER:: A1=45.95D0
    REAL(rkind),PARAMETER:: A2= 5.02D4
    REAL(rkind),PARAMETER:: A3=1.368D-2
    REAL(rkind),PARAMETER:: A4=1.076D0
    REAL(rkind),PARAMETER:: A5= 409.D0
    ! E: energy [keV]
    REAL(rkind),INTENT(IN):: E ! energy [keV]
    REAL(rkind):: sigma_nrl_dt  ! [barn=1.D-28 m^2]

    sigma_nrl_dt=(A5+A2/((A4-A3*E)**2+1))/(E*EXP(A1/SQRT(E))-1.D0)

    RETURN
  END FUNCTION sigma_nrl_dt

  ! Reaction rate for energy E with Maxwellian of temperature T
    
  FUNCTION sigmav_nrl_dt(E)
    USE sigmav_nrl_dt_comm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: E
    REAL(rkind):: sigmav_nrl_dt

    sigmav_nrl_dt=sigma_nrl_dt(E)*E*EXP(-E/T_bulk_dt)

    RETURN
  END FUNCTION sigmav_nrl_dt

  ! Reaction rate for Maxwellian with temperature T
    
  FUNCTION sigmavm_nrl_dt(T)
    USE sigmav_nrl_dt_comm
    USE sigmav_nrl_dt_comm
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: T
    REAL(rkind):: sigmavm_nrl_dt
    REAL(rkind):: rmu,factor,cs,es,H0,EPS
    INTEGER:: ILST

    ! T: temperature [keV]

    ! <sigma v>=(8/\pi \mu (T)^3)^{1/2}
    !           \int_0^\infty \sigma(e) E \exp(-E/T)

    ! \mu = m1*m2/(m1+m2) 
      
    RMU=AMD*AMT/(AMD+AMT)
    FACTOR=SQRT(8.D0/(PI*RMU*(T)**3))
    H0=0.1D0
    EPS=1.D-6
    ILST=0
    CALL DEHIFT(CS,ES,H0,EPS,ILST,sigmav_nrl_dt,'sigmav_nrl_dt')
    sigmavm_nrl_dt=FACTOR*CS

    RETURN
  END FUNCTION sigmavm_nrl_dt

  ! set spline coefficients usvla for reaction rate for discrete energies

  SUBROUTINE sigmavm_nrl_table_dt
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),DIMENSION(10):: dsvla
    INTEGER:: ierr,n

    DO n=1,10
       templa(n)=LOG10(tempa(n))
       svla(n)=LOG10(sva(n))
       WRITE(6,'(I4,4ES12.4)') n,tempa(n),templa(n),sva(n),svla(n)
    END DO

    CALL SPL1D(templa,svla,dsvla,usvla,10,0,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX spl1d svl: ierr=',ierr
       STOP
    END IF

    RETURN
  END SUBROUTINE sigmavm_nrl_table_dt

  ! interpolate sigma_v_maxwellian

  FUNCTION sigmavm_int_dt(temp)

    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: temp
    REAL(rkind):: sigmavm_int_dt
    REAL(rkind):: svl
    INTEGER:: ierr

    CALL SPL1DF(LOG10(temp),svl,templa,usvla,10,ierr)
!    WRITE(6,'(5ES12.4)') temp,LOG10(temp),svl,templa(1),templa(10)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX sigmavm_int_dt: SPL1DF error: ierr=',ierr
       STOP
    END IF

    sigmavm_int_dt=10.D0**svl
    RETURN
  END FUNCTION sigmavm_int_dt

  ! analytic sigmavm

  FUNCTION sigmavm_low_dt(temp)

    IMPLICIT NONE
    REAL(rkind), INTENT(IN):: temp
    REAL(rkind):: sigmavm_low_dt
    REAL(rkind):: temp3

    IF(temp.LE.0.D0) THEN
       sigmavm_low_dt=0.D0
       RETURN
    END IF
    temp3=1.D0/temp**(1.D0/3.D0)
    sigmavm_low_dt=2.33D-14*temp3**2*EXP(-18.76D0*temp3)
    RETURN
  END FUNCTION sigmavm_low_dt

END PROGRAM testsigma
      
