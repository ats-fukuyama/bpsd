!     $Id$
!=======================================================================
module bpsd_shot

  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_set_shot,bpsd_get_shot, &
       & bpsd_save_shot,bpsd_load_shot
  private

  logical, save :: bpsd_shotx_init_flag = .TRUE.
  type(bpsd_shotx_type), save :: shotx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_shotx_init
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    shotx%status = 2
    shotx%dataName = 'shot'
    bpsd_shotx_init_flag = .FALSE.

    return
  end subroutine bpsd_shotx_init

!-----------------------------------------------------------------------
  subroutine bpsd_set_shot(shot_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_shot_type) :: shot_in
    integer :: ierr

    if(bpsd_shotx_init_flag) call bpsd_shotx_init

    shotx%deviceID = shot_in%deviceID
    shotx%shotID = shot_in%shotID
    shotx%modelID = shot_in%modelID
    CALL DATE_AND_TIME(shotx%created_date, &
         &             shotx%created_time, &
         &             shotx%created_timezone)
    shotx%status = 2
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,'(A/A32,A32/A32,I12/A32,I12)') &
            &        '-- bpsd_set_shot', &
            &        'shot%deviceID: ',shotx%deviceID, &
            &        'shot%shotID  : ',shotx%shotID, &
            &        'shot%modelID : ',shotx%modelID
    endif
    return
  end subroutine bpsd_set_shot

!-----------------------------------------------------------------------
  subroutine bpsd_get_shot(shot_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_shot_type) :: shot_out
    integer :: ierr

    if(bpsd_shotx_init_flag) call bpsd_shotx_init

    if(shotx%status.eq.1) then
       write(6,*) 'XX bpsd_get_shot: no data in shot'
       ierr=2
       return
    endif

    shot_out%deviceID = shotx%deviceID
    shot_out%shotID = shotx%shotID
    shot_out%modelID = shotx%modelID
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,'(A/A32,A32/A32,I12/A32,I12)') &
            &        '-- bpsd_get_shot', &
            &        'shot%deviceID: ',shotx%deviceID, &
            &        'shot%shotID  : ',shotx%shotID, &
            &        'shot%modelID : ',shotx%modelID
    endif
    return
  end subroutine bpsd_get_shot

!-----------------------------------------------------------------------
  subroutine bpsd_save_shot(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_shotx_init_flag) call bpsd_shotx_init

    if(shotx%status.gt.1) call bpsd_save_shotx(fid,shotx,ierr)
    return

  end subroutine bpsd_save_shot

!-----------------------------------------------------------------------
  subroutine bpsd_load_shot(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_shotx_type),intent(in) :: datax
    integer,intent(out) :: ierr

    if(bpsd_shotx_init_flag) call bpsd_shotx_init

    shotx%deviceID = datax%deviceID
    shotx%shotID = datax%shotID
    shotx%modelID = datax%modelID
    shotx%created_date = datax%created_date
    shotx%created_time = datax%created_time
    shotx%created_timezone = datax%created_timezone
    shotx%status=2
    ierr=0
    return
  end subroutine bpsd_load_shot

end module bpsd_shot
