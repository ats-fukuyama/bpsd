c     $Id$
c=======================================================================
      module bpsd_shot
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_shot,bpsd_get_shot,
     &       bpsd_save_shot,bpsd_load_shot
      private
c
      logical, save :: bpsd_shotx_init_flag = .TRUE.
      type(bpsd_shotx_type), save :: shotx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_shotx_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      shotx%dataName = 'shot'
      bpsd_shotx_init_flag = .FALSE.
c
      return
      end subroutine bpsd_shotx_init
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_shot(shot_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_shot_type) :: shot_in
      integer :: ierr
c
      if(bpsd_shotx_init_flag) call bpsd_shotx_init
c
      shotx%deviceID = shot_in%deviceID
      shotx%shotID = shot_in%shotID
      shotx%modelID = shot_in%modelID
      shotx%status = 2
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,'(A/A32,A32/A32,I12/A32,I12)')
     &        '-- bpsd_set_shot',
     &        'shot%deviceID: ',shotx%deviceID,
     &        'shot%shotID  : ',shotx%shotID,
     &        'shot%modelID : ',shotx%modelID
      endif
      return
      end subroutine bpsd_set_shot
c
c-----------------------------------------------------------------------
      subroutine bpsd_get_shot(shot_out,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_shot_type) :: shot_out
      integer :: ierr
c
      if(bpsd_shotx_init_flag) call bpsd_shotx_init
c
      if(shotx%status.eq.1) then
         write(6,*) 'XX bpsd_get_shot: no data in shot'
         ierr=2
         return
      endif
c
      shot_out%deviceID = shotx%deviceID
      shot_out%shotID = shotx%shotID
      shot_out%modelID = shotx%modelID
      ierr = 0
c
      if(bpsd_debug_flag) then
         write(6,'(A/A32,A32/A32,I12/A32,I12)')
     &        '-- bpsd_get_shot',
     &        'shot%deviceID: ',shotx%deviceID,
     &        'shot%shotID  : ',shotx%shotID,
     &        'shot%modelID : ',shotx%modelID
      endif
      return
      end subroutine bpsd_get_shot
c
c-----------------------------------------------------------------------
      subroutine bpsd_save_shot(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      integer,intent(in) :: fid
      integer,intent(out) :: ierr
c
      if(bpsd_shotx_init_flag) call bpsd_shotx_init
c
      if(shotx%status.gt.1) 
     &     call bpsd_save_shotx(fid,shotx,ierr)
      return
c
      end subroutine bpsd_save_shot
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_shot(datax,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_subs
      implicit none
      type(bpsd_shotx_type),intent(in) :: datax
      integer,intent(out) :: ierr
c
      if(bpsd_shotx_init_flag) call bpsd_shotx_init
c
      shotx%deviceID = datax%deviceID
      shotx%shotID = datax%shotID
      shotx%modelID = datax%modelID
      shotx%status=2
      ierr=0
      return
c
      end subroutine bpsd_load_shot
c
      end module bpsd_shot
