c     $Id$
c=======================================================================
      module bpsd_shot
c
      use bpsd_flags
      use bpsd_types
      use bpsd_types_internal
      public bpsd_set_shot,bpsd_get_shot
      private
c
      logical, save :: bpsd_shot_init_flag = .TRUE.
      type(bpsd_shotx_type), save :: shotx
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_shot_init
c-----------------------------------------------------------------------
      use bpsd_subs
      implicit none
c
      bpsd_shot_init_flag = .FALSE.
c
      return
      end subroutine bpsd_shot_init
c
c-----------------------------------------------------------------------
      subroutine bpsd_set_shot(shot_in,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_shot_type) :: shot_in
      integer :: ierr
c
      if(bpsd_shot_init_flag) call bpsd_shot_init
c
      shotx%deviceID = shot_in%deviceID
      shotx%dataName = 'shot'
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
      use bpsd_types_internal
      use bpsd_subs
      implicit none
      type(bpsd_shot_type) :: shot_out
      integer :: ierr
c
      if(bpsd_shot_init_flag) call bpsd_shot_init
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
      end module bpsd_shot
