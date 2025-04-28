! trgiven.f90

MODULE trgiven

  !  *** Calculate given profile of density and temperature ***

  USE task_kinds,ONLY: dp

  PRIVATE
                                            ! density profile
  INTEGER:: ntime_ngiven_max                ! number of time points
  INTEGER:: ndata_ngiven_max                ! number of coef data
  REAL(dp),PUBLIC:: rho_min_ngiven,rho_max_ngiven  ! range of given profile
  REAL(dp),PUBLIC,ALLOCATABLE:: time_ngiven(:)     ! time points t_i
  REAL(dp),ALLOCATABLE:: coef_ngiven(:,:)   ! coef data for t_i<= t <t_{i+1}
                                            ! temperature profile
  INTEGER:: ntime_tgiven_max                ! number of time points
  INTEGER:: ndata_tgiven_max                ! number of coef data
  REAL(dp),PUBLIC :: rho_min_tgiven,rho_max_tgiven  ! range of given profile
  REAL(dp),PUBLIC,ALLOCATABLE:: time_tgiven(:)     ! time points t_i
  REAL(dp),ALLOCATABLE:: coef_tgiven(:,:)   ! coef data for t_i<= t <t_{i+1}

  PUBLIC tr_set_ngiven  ! set coef matrix for n
  PUBLIC tr_set_tgiven  ! set coef matrix for nT
  PUBLIC tr_prof_ngiven ! set given density profile
  PUBLIC tr_prof_tgiven ! set given temperature profile
  PUBLIC tr_prep_ngiven ! read given density pfofile parameters
  PUBLIC tr_prep_tgiven ! read given temperature profile parameters

CONTAINS

  !     ***** Routine for given density profile *****
      
  SUBROUTINE tr_set_ngiven(nr,time)

    USE trcomm
    USE trcomx
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr
    REAL(rkind),INTENT(IN):: time
    REAL(rkind):: rn_local
    INTEGER:: NS,NEQ,NW

    IF(model_ngiven.EQ.0) RETURN
    IF(time.LE.time_ngiven(1)) return
    IF(model_ngiven.EQ.2) THEN
       IF((rm(nr).LT.rho_min_ngiven).OR. &
            (rm(nr).GT.rho_max_ngiven)) RETURN
    END IF
    CALL tr_prof_ngiven(rm(nr),time,rn_local)
    NEQ=NEA(1,1) ! NEQ of electron density equation
    DO NW=1,NEQMAX
       A(NEQ,NW,NR) = 0.D0
       B(NEQ,NW,NR) = 0.D0
       C(NEQ,NW,NR) = 0.D0
    END DO
    D(NEQ,NR)=0.D0
!    B(NEQ,NEQ,NR)=-1.D0/tau_ngiven
!    D(NEQ,NR)=rn_local/tau_ngiven
    RD(NEQ,NR)=1.D0
    DO NS=2,NSMAX
       NEQ=NEA(NS,1) ! NEQ of density equation
       DO NW=1,NEQMAX
          A(NEQ,NW,NR) = 0.D0
          B(NEQ,NW,NR) = 0.D0
          C(NEQ,NW,NR) = 0.D0
       END DO
       D(NEQ,NR)=0.D0
!       B(NEQ,NEQ,NR)=-1.D0/tau_ngiven
!       D(NEQ,NR)=pn(ns)/(pz(ns)*pn(1))*rn_local/tau_ngiven
       RD(NEQ,NR)=1.D0
    END DO
    RETURN
  END SUBROUTINE tr_set_ngiven
      
  !     ***** Routine for given temperature profile *****
      
  SUBROUTINE tr_set_tgiven(nr,time)

    USE trcomm
    USE trcomx
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr
    REAL(rkind),INTENT(IN):: time
    REAL(rkind):: rt_local
    INTEGER:: NS,NEQ,NW

    IF(model_tgiven.EQ.0) RETURN
    IF(time.LE.time_tgiven(1)) return
    IF(model_tgiven.EQ.2) THEN
       IF((rm(nr).LT.rho_min_tgiven).OR. &
            (rm(nr).GT.rho_max_tgiven)) RETURN
    END IF
    CALL tr_prof_tgiven(rm(nr),time,rt_local)
    NEQ=NEA(1,1) ! NEQ of electron density equation
    DO NW=1,NEQMAX
       A(NEQ,NW,NR) = 0.D0
       B(NEQ,NW,NR) = 0.D0
       C(NEQ,NW,NR) = 0.D0
    END DO
    D(NEQ,NR)=0.D0
!    B(NEQ,NEQ,NR)=-1.D0/tau_tgiven
!    D(NEQ,NR)=rt_local/tau_tgiven
    RD(NEQ,NR)=1.D0
    DO NS=2,NSMAX
       NEQ=NEA(NS,2) ! NEQ of temperature equation
       DO NW=1,NEQMAX
          A(NEQ,NW,NR) = 0.D0
          B(NEQ,NW,NR) = 0.D0
          C(NEQ,NW,NR) = 0.D0
       END DO
       D(NEQ,NR)=0.D0
!       B(NEQ,NEQ,NR)=-1.D0/tau_tgiven
!       D(NEQ,NR)=rt_local/tau_tgiven
       RD(NEQ,NR)=1.D0
    END DO
    RETURN
  END SUBROUTINE tr_set_tgiven
      
  ! *** set given density profile ***
  
  SUBROUTINE tr_prof_ngiven(rho,time,rn)
  
    USE task_kinds,ONLY: dp
    IMPLICIT NONE    
    REAL(dp),INTENT(IN):: rho,time
    REAL(dp),INTENT(OUT):: rn
    REAL(dp):: tr_func_ngiven
    REAL(dp),ALLOCATABLE:: coef(:)
    REAL(dp):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LT.time_ngiven(1)) THEN
       RETURN
    ELSE IF (time.GE.time_ngiven(ntime_ngiven_max)) THEN
       id=ntime_ngiven_max
    ELSE
       DO ntime=1,ntime_ngiven_max-1
          IF(time.GE.time_ngiven(ntime).AND. &
               time.LT.time_ngiven(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_ngiven_max))
    IF(id.EQ.ntime_ngiven_max) THEN ! after time_ngiven(ntime_ngiven_max)
       DO i=0,ndata_ngiven_max
          coef(i)=coef_ngiven(i,ntime_ngiven_max)
       END DO
    ELSE ! between time_ngiven(id) and time_ngiven(id+1)
       factor=(time-time_ngiven(id)) &
             /(time_ngiven(id+1)-time_ngiven(id))
       DO i=0,ndata_ngiven_max
          coef(i)=(1.D0-factor)*coef_ngiven(i,id) &
                        +factor*coef_ngiven(i,id+1)
       END DO
    END IF

    ! --- set local density profile ---
    
    rn=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    rn=rn*1.D-20
    IF(rn.LE.0.D0) rn=1.D-8
    RETURN
  END SUBROUTINE tr_prof_ngiven

  ! *** set given temperature profile ***
  
  SUBROUTINE tr_prof_tgiven(rho,time,rt)
  
    USE task_kinds,ONLY: dp
    USE trcomm,ONLY: model_tgiven
    IMPLICIT NONE    
    REAL(dp),INTENT(IN):: rho,time
    REAL(dp),INTENT(OUT):: rt
    REAL(dp):: tr_func_tgiven
    REAL(dp),ALLOCATABLE:: coef(:)
    REAL(dp):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LE.time_tgiven(1)) THEN
       RETURN
    ELSE IF (time.GE.time_tgiven(ntime_tgiven_max)) THEN
       id=ntime_tgiven_max
    ELSE
       DO ntime=1,ntime_tgiven_max-1
          IF(time.GE.time_tgiven(ntime).AND. &
               time.LE.time_tgiven(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_tgiven_max))
    IF(id.EQ.0) THEN ! before time_ngiven(1)
       DO i=0,ndata_tgiven_max
          coef(i)=coef_tgiven(i,1)
       END DO
    ELSE IF(id.EQ.ntime_tgiven_max) THEN ! after time_ngiven(ntime_ngiven_max)
       DO i=0,ndata_tgiven_max
          coef(i)=coef_tgiven(i,ntime_tgiven_max)
       END DO
    ELSE ! between time_ngiven(id) and time_ngiven(id+1)
       factor=(time-time_tgiven(id)) &
             /(time_tgiven(id+1)-time_tgiven(id))
       DO i=0,ndata_tgiven_max
          coef(i)=(1.D0-factor)*coef_tgiven(i,id) &
                        +factor*coef_tgiven(i,id+1)
       END DO
    END IF

    ! --- set temperature profile ---
    
    rt=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    rt=rt*1.D-3
    RETURN
  END SUBROUTINE tr_prof_tgiven

  ! *** read density profile data from file ***

  SUBROUTINE tr_prep_ngiven
    USE trcomm,ONLY: knam_ngiven
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_ngiven,1,0,'fn',ierr)
    READ(NFL,*) ntime_ngiven_max,ndata_ngiven_max,rho_min_ngiven,rho_max_ngiven
    IF(ALLOCATED(time_ngiven)) DEALLOCATE(time_ngiven)
    IF(ALLOCATED(coef_ngiven)) DEALLOCATE(coef_ngiven)
    ALLOCATE(time_ngiven(ntime_ngiven_max))
    ALLOCATE(coef_ngiven(0:ndata_ngiven_max,ntime_ngiven_max))
    DO ntime=1,ntime_ngiven_max
       READ(NFL,*) time_ngiven(ntime)
       READ(NFL,*) (coef_ngiven(ndata,ntime),ndata=0,ndata_ngiven_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_ngiven
    
  ! *** read temperature profile data from file ***

  SUBROUTINE tr_prep_tgiven
    USE trcomm,ONLY: knam_tgiven
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_tgiven,1,0,'ft',ierr)
    READ(NFL,*) ntime_tgiven_max,ndata_tgiven_max,rho_min_tgiven,rho_max_tgiven
    IF(ALLOCATED(time_tgiven)) DEALLOCATE(time_tgiven)
    IF(ALLOCATED(coef_tgiven)) DEALLOCATE(coef_tgiven)
    ALLOCATE(time_tgiven(ntime_tgiven_max))
    ALLOCATE(coef_tgiven(0:ndata_tgiven_max,ntime_tgiven_max))
    DO ntime=1,ntime_tgiven_max
       READ(NFL,*) time_tgiven(ntime)
       READ(NFL,*) (coef_tgiven(ndata,ntime),ndata=0,ndata_tgiven_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_tgiven
  
END MODULE trgiven
