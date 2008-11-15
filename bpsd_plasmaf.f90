!     $Id$
!=======================================================================
module bpsd_plasmaf

  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_set_plasmaf,bpsd_get_plasmaf, &
       & bpsd_save_plasmaf,bpsd_load_plasmaf
  private

  logical, save :: bpsd_plasmafx_init_flag = .TRUE.
  type(bpsd_data1Dx_type), save :: plasmafx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_plasmafx
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    plasmafx%status=0
    plasmafx%dataName='plasmaf'
    plasmafx%ndmax=0

    bpsd_plasmafx_init_flag = .FALSE.

    return
  end subroutine bpsd_init_plasmafx

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_plasmaf_data(data,n1,n2)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_plasmaf_data),DIMENSION(:,:),POINTER,INTENT(INOUT):: data
    INTEGER(ikind),INTENT(IN):: n1,n2

    IF(ASSOCIATED(data)) THEN
       IF(n1.LE.0.or.n2.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1).OR.n2.NE.SIZE(data,2)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1,n2))
       END IF
    ELSE
       IF(n1.GT.0.AND.n2.GT.0)  THEN
          ALLOCATE(data(n1,n2))
       ENDIF
    ENDIF
  END SUBROUTINE bpsd_adjust_plasmaf_data

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_setup_plasmaf_kdata
!-----------------------------------------------------------------------
    IMPLICIT NONE
    integer(ikind):: nd

    do nd=0,plasmafx%ndmax-1,5
       plasmafx%kid(nd+1)='plasmaf%pn'
       plasmafx%kid(nd+2)='plasmaf%pt'
       plasmafx%kid(nd+3)='plasmaf%ptpr'
       plasmafx%kid(nd+4)='plasmaf%ptpp'
       plasmafx%kid(nd+5)='plasmaf%pu'
       plasmafx%kunit(nd+1)='10^20/m^3'
       plasmafx%kunit(nd+2)='eV'
       plasmafx%kunit(nd+3)='eV'
       plasmafx%kunit(nd+4)='eV'
       plasmafx%kunit(nd+5)='m/s'
    enddo
    plasmafx%kid(plasmafx%ndmax)='plasmaf%qinv'
    plasmafx%kunit(plasmafx%ndmax)=' '
    RETURN
  END SUBROUTINE bpsd_setup_plasmaf_kdata

!-----------------------------------------------------------------------
  subroutine bpsd_set_plasmaf(plasmaf_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_plasmaf_type),intent(in):: plasmaf_in
    integer,intent(out) :: ierr
    integer :: ns,nr,nd

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    plasmafx%nrmax=plasmaf_in%nrmax
    plasmafx%ndmax=plasmaf_in%nsmax*5+1
    CALL bpsd_adjust_karray(plasmafx%kid,plasmafx%ndmax)
    CALL bpsd_adjust_karray(plasmafx%kunit,plasmafx%ndmax)
    CALL bpsd_adjust_array1D(plasmafx%rho,plasmafx%nrmax)
    CALL bpsd_adjust_array2D(plasmafx%data,plasmafx%nrmax,plasmafx%ndmax)

    CALL bpsd_setup_plasmaf_kdata

    plasmafx%time = plasmaf_in%time
    do nr=1,plasmaf_in%nrmax
       plasmafx%rho(nr) = plasmaf_in%rho(nr)
       do ns=1,plasmaf_in%nsmax
          nd=5*(ns-1)
          plasmafx%data(nr,nd+1) = plasmaf_in%data(nr,ns)%pn
          plasmafx%data(nr,nd+2) = plasmaf_in%data(nr,ns)%pt
          plasmafx%data(nr,nd+3) = plasmaf_in%data(nr,ns)%ptpr
          plasmafx%data(nr,nd+4) = plasmaf_in%data(nr,ns)%ptpp
          plasmafx%data(nr,nd+5) = plasmaf_in%data(nr,ns)%pu
       enddo
       plasmafx%data(nr,plasmafx%ndmax) = plasmaf_in%qinv(nr)
    enddo
    CALL DATE_AND_TIME(plasmafx%created_date, &
         &             plasmafx%created_time, &
         &             plasmafx%created_timezone)

    if(plasmafx%status.ge.3) then
       plasmafx%status=3
    else
       plasmafx%status=2
    endif

    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_set_plasmaf'
       write(6,*) '---- plasmafx%rho'
       write(6,'(1P5E12.4)') &
            &        (plasmafx%rho(nr),nr=1,plasmafx%nrmax)
       do nd=1,plasmafx%ndmax
          write(6,*) '---- ',plasmafx%kid(nd)
          write(6,'(1P5E12.4)') &
               &           (plasmafx%data(nr,nd),nr=1,plasmafx%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_set_plasmaf

!-----------------------------------------------------------------------
  subroutine bpsd_get_plasmaf(plasmaf_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_plasmaf_type),intent(out) :: plasmaf_out
    integer,intent(out) :: ierr
    integer :: nr, nd, ns, mode
    real(8) :: s
    real(8), dimension(:), pointer :: v

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    if(plasmafx%status.eq.0) then
       write(6,*) 'XX bpsd_get_plasmaf: no space allocated to plasmafx%data'
       ierr=1
       return
    endif

    if(plasmafx%status.eq.1) then
       write(6,*) 'XX bpsd_get_plasmaf: no data in plasmafx%data'
       ierr=2
       return
    endif

    if(plasmaf_out%nrmax.eq.0) then
       mode=0
       plasmaf_out%nrmax = plasmafx%nrmax
    else
       mode=1
    endif
    plasmaf_out%nsmax = (plasmafx%ndmax-1)/5

    CALL bpsd_adjust_array1D(plasmaf_out%rho,plasmaf_out%nrmax)
    CALL bpsd_adjust_array1D(plasmaf_out%qinv,plasmaf_out%nrmax)
    CALL bpsd_adjust_plasmaf_data(plasmaf_out%data,plasmaf_out%nrmax, &
         &                                         plasmaf_out%nsmax)

    if(mode.eq.0) then
       plasmaf_out%time  = plasmafx%time
       do nr=1,plasmafx%nrmax
          plasmaf_out%rho(nr)=plasmafx%rho(nr)
          do ns=1,plasmaf_out%nsmax
             nd=5*(ns-1)
             plasmaf_out%data(nr,ns)%pn  =plasmafx%data(nr,nd+1)
             plasmaf_out%data(nr,ns)%pt  =plasmafx%data(nr,nd+2)
             plasmaf_out%data(nr,ns)%ptpr=plasmafx%data(nr,nd+3)
             plasmaf_out%data(nr,ns)%ptpp=plasmafx%data(nr,nd+4)
             plasmaf_out%data(nr,ns)%pu  =plasmafx%data(nr,nd+5)
          enddo
          plasmaf_out%qinv(nr)=plasmafx%data(nr,plasmafx%ndmax)
       enddo
       ierr=0
       return
    endif

    if(plasmafx%status.eq.2) then
       CALL bpsd_adjust_array3D(plasmafx%spline,4,plasmafx%nrmax, &
            &                                     plasmafx%ndmax)
       plasmafx%status=3
    endif

    if(plasmafx%status.eq.3) then
       CALL bpsd_adjust_array1D(plasmafx%s,plasmafx%nrmax)
       do nr=1,plasmafx%nrmax
          plasmafx%s(nr)=plasmafx%rho(nr)**2
       enddo
       do nd=1,plasmafx%ndmax
          call bpsd_spl1D(plasmafx,nd,ierr)
       enddo
       plasmafx%status=4
    endif
    !
    allocate(v(plasmafx%ndmax))
    do nr=1,plasmaf_out%nrmax
       s = plasmaf_out%rho(nr)**2
       do nd=1,plasmafx%ndmax
          call bpsd_spl1DF(s,v(nd),plasmafx,nd,ierr)
       enddo
       do ns=1,plasmaf_out%nsmax
          nd=5*(ns-1)
          plasmaf_out%data(nr,ns)%pn   = v(nd+1)
          plasmaf_out%data(nr,ns)%pt   = v(nd+2)
          plasmaf_out%data(nr,ns)%ptpr = v(nd+3)
          plasmaf_out%data(nr,ns)%ptpp = v(nd+4)
          plasmaf_out%data(nr,ns)%pu   = v(nd+5)
       enddo
       plasmaf_out%qinv(nr)  = v(plasmafx%ndmax)
    enddo
    deallocate(v)
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_plasmaf'
       write(6,*) '---- plasmafx%rho'
       write(6,'(1P5E12.4)') &
            &        (plasmaf_out%rho(nr),nr=1,plasmaf_out%nrmax)
       do ns=1,plasmaf_out%nsmax
          write(6,*) '---- plasmafx%pn(',ns,')'
          write(6,'(1P5E12.4)') &
               &     (plasmaf_out%data(nr,ns)%pn,nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%pt(',ns,')'
          write(6,'(1P5E12.4)') &
               &     (plasmaf_out%data(nr,ns)%pt,nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%ptpr(',ns,')'
          write(6,'(1P5E12.4)') &
               &     (plasmaf_out%data(nr,ns)%ptpr,nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%ptpp(',ns,')'
          write(6,'(1P5E12.4)') &
               &     (plasmaf_out%data(nr,ns)%ptpp,nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%pu(',ns,')'
          write(6,'(1P5E12.4)') &
               &     (plasmaf_out%data(nr,ns)%pu,nr=1,plasmaf_out%nrmax)
       enddo
       write(6,*) '---- plasmafx%qinv'
       write(6,'(1P5E12.4)') &
            &        (plasmaf_out%qinv(nr),nr=1,plasmaf_out%nrmax)
    endif
    return
  end subroutine bpsd_get_plasmaf

!-----------------------------------------------------------------------
  subroutine bpsd_save_plasmaf(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    if(plasmafx%status.gt.1) call bpsd_save_data1Dx(fid,plasmafx,ierr)
    return

  end subroutine bpsd_save_plasmaf

!-----------------------------------------------------------------------
  subroutine bpsd_load_plasmaf(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data1Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nr,nd

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    plasmafx%dataName=datax%dataName
    plasmafx%time = datax%time
    plasmafx%nrmax=datax%nrmax
    plasmafx%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(plasmafx%kid,plasmafx%ndmax)
    CALL bpsd_adjust_karray(plasmafx%kunit,plasmafx%ndmax)
    CALL bpsd_adjust_array1D(plasmafx%rho,plasmafx%nrmax)
    CALL bpsd_adjust_array2D(plasmafx%data,plasmafx%nrmax, &
         &                                 plasmafx%ndmax)

    do nr=1,plasmafx%nrmax
       plasmafx%rho(nr) = datax%rho(nr)
       do nd=1,plasmafx%ndmax
          plasmafx%data(nr,nd) = datax%data(nr,nd)
       enddo
    enddo
    do nd=1,plasmafx%ndmax
       plasmafx%kid(nd)=datax%kid(nd)
       plasmafx%kunit(nd)=datax%kunit(nd)
    enddo
    plasmafx%created_date = datax%created_date
    plasmafx%created_time = datax%created_time
    plasmafx%created_timezone = datax%created_timezone

    if(plasmafx%status.ge.3) then
       plasmafx%status=3
    else
       plasmafx%status=2
    endif

    ierr=0
    return

  end subroutine bpsd_load_plasmaf

end module bpsd_plasmaf
