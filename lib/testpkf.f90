! testpkf.f90

PROGRAM testpkf
  USE task_kinds
  USE task_constants
  USE libpdkf
  USE libpgkf
  USE libgrf
  IMPLICIT NONE
  REAL(rkind):: xi,eta,rnu,xi_min,xi_max
  INTEGER:: np_tau,ntau_max
  COMPLEX(rkind):: cf1,cf2,cf3

  xi=1.D0
  eta=0.D0
  rnu=1.D0
  np_tau=0

  xi_min=1.D-2
  xi_max=1.D+2
  ntau_max=101
  
1 CONTINUE
  WRITE(6,'(A,3ES12.4,I4)') '## xi,eta,rnu,np=',xi,eta,rnu,np_tau
  READ(5,*,ERR=1,END=9000) xi,eta,rnu,np_tau

!  cf1=pdkf_hift(xi,eta,rnu,np_tau)
!  cf2=pdkf_ft(xi,eta,rnu,np_tau)
  cf3=pdkf_eul(xi,eta,rnu,np_tau)
  WRITE(6,'(A,4ES12.4)') 'Re: ',REAL(cf1),REAL(cf2),REAL(cf3)
  WRITE(6,'(A,4ES12.4)') 'Im: ',AIMAG(cf1),AIMAg(cf2),AIMAG(cf3)
  GO to 1

9000 CONTINUE
  STOP
  
END PROGRAM testpkf
