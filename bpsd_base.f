c     $Id$
c=======================================================================
      module bpsd
c
      use bpsd_kinds
      use bpsd_constants
      use bpsd_flags
      use bpsd_types
      use bpsd_shot
      use bpsd_device
      use bpsd_species
      use bpsd_equ1D
      use bpsd_metric1D
      use bpsd_plasmaf
c
      interface bpsd_set_data
         module procedure bpsd_set_shot, 
     &                    bpsd_set_device, 
     &                    bpsd_set_species, 
     &                    bpsd_set_equ1D, 
     &                    bpsd_set_metric1D, 
     &                    bpsd_set_plasmaf
      end interface bpsd_set_data
c
      interface bpsd_get_data
         module procedure bpsd_get_shot, 
     &                    bpsd_get_device, 
     &                    bpsd_get_species, 
     &                    bpsd_get_equ1D, 
     &                    bpsd_get_metric1D, 
     &                    bpsd_get_plasmaf
      end interface bpsd_get_data
c
      contains
c
c-----------------------------------------------------------------------
      subroutine bpsd_save(ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(out):: ierr
      integer:: fid
c
      fid=31
c
      call fwopen(fid,'bpsd.data',0,0,'bpsd',ierr)
      if(ierr.ne.0) return
c
      write(fid) '##bspd.data version 0.01'
c
      call bpsd_save_shot(fid,ierr)
      if(ierr.ne.0) return
      call bpsd_save_device(fid,ierr)
      if(ierr.ne.0) return
      call bpsd_save_equ1D(fid,ierr)
      if(ierr.ne.0) return
      call bpsd_save_metric1D(fid,ierr)
      if(ierr.ne.0) return
      call bpsd_save_species(fid,ierr)
      if(ierr.ne.0) return
      call bpsd_save_plasmaf(fid,ierr)
      if(ierr.ne.0) return
c
      write(fid) 'data:end'
c
      close(fid)
c      
      return
      end subroutine bpsd_save
c
c-----------------------------------------------------------------------
      subroutine bpsd_load(ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(out):: ierr
      integer:: fid
      character(len=24):: line
      character(len=8):: datadim
c
      fid=32
c
      call fropen(fid,'bpsd.data',0,0,'bpsd',ierr)
      if(ierr.ne.0) return
c
      read(fid,IOSTAT=ierr,ERR=8,END=9) line
      if(line(1:12).ne.'##bspd.data ') then
         ierr=1001
         goto 9000
      endif
      if(line(21:24).ne.'0.01') then
         ierr=1002
         goto 9000
      endif
c
    1 continue
         read(fid,IOSTAT=ierr,ERR=8,END=9) datadim

         if(datadim.eq.'data:str') then
            call bpsd_load_shotx(fid,ierr)
            if(ierr.ne.0) then
               ierr=2000+ierr
               goto 9000
            endif
         endif

         if(datadim.eq.'data::0D') then
            call bpsd_load_data0Dx(fid,ierr)
            if(ierr.ne.0) then
               ierr=3000+ierr
               goto 9000
            endif
         endif

         if(datadim.eq.'data::1D') then
            call bpsd_load_data1Dx(fid,ierr)
            if(ierr.ne.0) then
               ierr=4000+ierr
               goto 9000
            endif
         endif

         if(datadim.eq.'data::2D') then
            call bpsd_load_data2Dx(fid,ierr)
            if(ierr.ne.0) then
               ierr=5000+ierr
               goto 9000
            endif
         endif

         if(datadim.eq.'data::3D') then
            call bpsd_load_data3Dx(fid,ierr)
            if(ierr.ne.0) then
               ierr=6000+ierr
               goto 9000
            endif
         endif

         if(datadim.eq.'data:end') then
            ierr=0
            goto 9000
         endif
      go to 1
!
    8 ierr=8000
      goto 9000
    9 ierr=9000
!
 9000 close(fid)
      return
      end subroutine bpsd_load
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_shotx(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_shotx_type):: datax
      integer,intent(out):: ierr
c
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%deviceID
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%shotID,datax%modelID

      if(datax%dataName(1:4).eq.'shot') 
     &     call bpsd_load_shot(datax,ierr)

      ierr=0
      return

    8 ierr=8000+ierr
      return
    9 ierr=9000+ierr
      return
      end subroutine bpsd_load_shotx
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_data0Dx(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data0Dx_type):: datax
      integer,intent(out):: ierr
c
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%ndmax
      allocate(datax%kid(datax%ndmax))
      allocate(datax%data(datax%ndmax))
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

      if(datax%dataName(1:6).eq.'device') 
     &     call bpsd_load_device(datax,ierr)
      if(datax%dataName(1:7).eq.'species') 
     &     call bpsd_load_species(datax,ierr)

      deallocate(datax%data)
      deallocate(datax%kid)

      ierr=0
      return

    8 ierr=8000+ierr
      return
    9 ierr=9000+ierr
      return
      end subroutine bpsd_load_data0Dx
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_data1Dx(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data1Dx_type):: datax
      integer,intent(out):: ierr
c
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%nrmax,datax%ndmax
      allocate(datax%kid(datax%ndmax))
      allocate(datax%s(datax%nrmax))
      allocate(datax%data(datax%nrmax,datax%ndmax))
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%s
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

      if(datax%dataName(1:5).eq.'equ1D') 
     &     call bpsd_load_equ1D(datax,ierr)
      if(datax%dataName(1:8).eq.'metric1D') 
     &     call bpsd_load_metric1D(datax,ierr)
      if(datax%dataName(1:7).eq.'plasmaf') 
     &     call bpsd_load_plasmaf(datax,ierr)

      deallocate(datax%data)
      deallocate(datax%s)
      deallocate(datax%kid)

      ierr=0
      return

    8 ierr=8000+ierr
      return
    9 ierr=9000+ierr
      return
      end subroutine bpsd_load_data1Dx
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_data2Dx(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data2Dx_type):: datax
      integer,intent(out):: ierr
c
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
      read(fid,IOSTAT=ierr,ERR=8,END=9) 
     &     datax%nrmax,datax%nthmax,datax%ndmax
      allocate(datax%kid(datax%ndmax))
      allocate(datax%s(datax%nrmax))
      allocate(datax%th(datax%nthmax))
      allocate(datax%data(datax%nrmax,datax%nthmax,datax%ndmax))
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%s
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%th
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

!      if(datax%dataName(1:7).eq.'plasmaf') 
!     &     call bpsd_load_plasmaf(datax,ierr)

      deallocate(datax%data)
      deallocate(datax%th)
      deallocate(datax%s)
      deallocate(datax%kid)

      ierr=0
      return

    8 ierr=8000+ierr
      return
    9 ierr=9000+ierr
      return
      end subroutine bpsd_load_data2Dx
c
c-----------------------------------------------------------------------
      subroutine bpsd_load_data3Dx(fid,ierr)
c-----------------------------------------------------------------------
c
      use bpsd_types_internal
      implicit none
      integer,intent(in):: fid
      type(bpsd_data3Dx_type):: datax
      integer,intent(out):: ierr
c
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
      read(fid,IOSTAT=ierr,ERR=8,END=9) 
     &     datax%nrmax,datax%nthmax,datax%nphmax,datax%ndmax
      allocate(datax%kid(datax%ndmax))
      allocate(datax%s(datax%nrmax))
      allocate(datax%th(datax%nthmax))
      allocate(datax%ph(datax%nphmax))
      allocate(datax%data(datax%nrmax,datax%nthmax,datax%nphmax,
     &         datax%ndmax))
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%s
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%th
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%ph
      read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

!      if(datax%dataName(1:7).eq.'plasmaf') 
!     &     call bpsd_load_plasmaf(datax,ierr)

      deallocate(datax%data)
      deallocate(datax%ph)
      deallocate(datax%th)
      deallocate(datax%s)
      deallocate(datax%kid)

      ierr=0
      return

    8 ierr=8000+ierr
      return
    9 ierr=9000+ierr
      return
      end subroutine bpsd_load_data3Dx
c
      end module bpsd
