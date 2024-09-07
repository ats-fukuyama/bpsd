! libsigma.f90

MODULE libsigma
  USE bpsd_kinds
  USE bpsd_constants

  REAL(rkind),DIMENSION(10):: &
       tempa=(/ 1.D0, 2.D0, 5.D0, 10.D0, 20.D0, &
               50.D0, 100.D0, 200.D0, 500.D0, 1000.D0 /)
  REAL(rkind),DIMENSION(10):: &
       sigmavma_dt=(/ 5.5D-21, 2.6D-19, 1.3D-17, 1.1D-16, 4.2D-16, &
             8.7D-16, 8.5D-16, 6.3D-16, 3.7D-16, 2.7D-16 /)
  REAL(rkind),DIMENSION(10):: sigmavmal_dt

  REAL(rkind),DIMENSION(10):: tempal
  REAL(rkind),DIMENSION(4,10):: usigmavmal_dt

  PRIVATE

  PUBLIC set_usigmavmal_dt
  PUBLIC sigmavm_dt

CONTAINS

  ! set spline coefficients usvla for reaction rate for discrete energies

  SUBROUTINE set_usigmavmal_dt
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),DIMENSION(10):: dsigmavmal_dt
    INTEGER:: ierr,n

    DO n=1,10
       tempal(n)=LOG10(tempa(n))
       sigmavmal_dt(n)=LOG10(sigmavma_dt(n))
    END DO

    CALL SPL1D(tempal,sigmavmal_dt,dsigmavmal_dt,usigmavmal_dt,10,0,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX spl1d sigmavmal_dt: ierr=',ierr
       STOP
    END IF

    RETURN
  END SUBROUTINE set_usigmavmal_dt

  ! interpolate sigma_v_maxwellian

  FUNCTION sigmavm_dt(temp)

    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: temp
    REAL(rkind):: templ
    REAL(rkind):: sigmavm_dt
    REAL(rkind):: sigmavml_dt
    INTEGER:: ierr

    templ=LOG10(temp)
    IF(templ.LE.tempal(1)) THEN
       sigmavm_dt=0.D0
    ELSE IF(templ.GE.tempal(10)) THEN
       sigmavm_dt=sigmavma_dt(10)
    ELSE
       CALL SPL1DF(templ,sigmavml_dt,tempal,usigmavmal_dt,10,ierr)
!       WRITE(6,'(5ES12.4)') temp,templ,sigmavml_dt,tempal(1),tempal(10)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,I4)') 'XX sigmavm_dt: SPL1DF error: ierr=',ierr
          STOP
       END IF
       sigmavm_dt=10.D0**sigmavml_dt
    END IF
    RETURN
  END FUNCTION sigmavm_dt
  
END MODULE libsigma
      
