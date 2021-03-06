! bpsd_trmatrix.f90

module bpsd_trmatrix

  use bpsd_kinds
  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_put_trmatrix,bpsd_get_trmatrix, &
         bpsd_save_trmatrix,bpsd_load_trmatrix
  private

  logical, save :: bpsd_trmatrixx_init_flag = .TRUE.
  type(bpsd_data1Dx_type), save :: trmatrixx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_trmatrixx
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    trmatrixx%status=0
    trmatrixx%dataName='trmatrix'
    trmatrixx%ndmax=0

    bpsd_trmatrixx_init_flag = .FALSE.

    return
  end subroutine bpsd_init_trmatrixx

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_trmatrix_data(data,n1,n2)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_trmatrix_data),DIMENSION(:,:),POINTER:: data
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
  END SUBROUTINE bpsd_adjust_trmatrix_data

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_setup_trmatrix_kdata
!-----------------------------------------------------------------------
    IMPLICIT NONE
    integer(ikind):: nd

    do nd=0,trmatrixx%ndmax-2,5
       trmatrixx%kid(nd+1)='trmatrix%nip'
       trmatrixx%kid(nd+2)='trmatrix%nim'
       trmatrixx%kid(nd+3)='trmatrix%ncx'
       trmatrixx%kid(nd+4)='trmatrix%Pec'
       trmatrixx%kid(nd+5)='trmatrix%Plh'
       trmatrixx%kid(nd+6)='trmatrix%Pic'
       trmatrixx%kid(nd+7)='trmatrix%Pbr'
       trmatrixx%kid(nd+8)='trmatrix%Pcy'
       trmatrixx%kid(nd+9)='trmatrix%Plr'
       trmatrixx%kunit(nd+1)='1/(m^3 s)'
       trmatrixx%kunit(nd+2)='1/(m^3 s)'
       trmatrixx%kunit(nd+3)='1/(m^3 s)'
       trmatrixx%kunit(nd+4)='W1/m^3'
       trmatrixx%kunit(nd+5)='W1/m^3'
       trmatrixx%kunit(nd+6)='W1/m^3'
       trmatrixx%kunit(nd+7)='W1/m^3'
       trmatrixx%kunit(nd+8)='W1/m^3'
       trmatrixx%kunit(nd+9)='W1/m^3'
    enddo
    RETURN
  END SUBROUTINE bpsd_setup_trmatrix_kdata

!-----------------------------------------------------------------------
  subroutine bpsd_put_trmatrix(trmatrix_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_trmatrix_type),intent(in):: trmatrix_in
    integer,intent(out) :: ierr
    integer :: ns,nr,nd

    if(bpsd_trmatrixx_init_flag) call bpsd_init_trmatrixx

    trmatrixx%nrmax=trmatrix_in%nrmax
    trmatrixx%ndmax=trmatrix_in%nsmax*9
    CALL bpsd_adjust_karray(trmatrixx%kid,trmatrixx%ndmax)
    CALL bpsd_adjust_karray(trmatrixx%kunit,trmatrixx%ndmax)
    CALL bpsd_adjust_array1D(trmatrixx%rho,trmatrixx%nrmax)
    CALL bpsd_adjust_array2D(trmatrixx%data,trmatrixx%nrmax,trmatrixx%ndmax)

    CALL bpsd_setup_trmatrix_kdata

    trmatrixx%time = trmatrix_in%time
    do nr=1,trmatrix_in%nrmax
       trmatrixx%rho(nr) = trmatrix_in%rho(nr)
       do ns=1,trmatrix_in%nsmax
          nd=6*(ns-1)
          trmatrixx%data(nr,nd+1) = trmatrix_in%data(nr,ns)%Dn
          trmatrixx%data(nr,nd+2) = trmatrix_in%data(nr,ns)%Dp
          trmatrixx%data(nr,nd+3) = trmatrix_in%data(nr,ns)%DT
          trmatrixx%data(nr,nd+4) = trmatrix_in%data(nr,ns)%un
          trmatrixx%data(nr,nd+5) = trmatrix_in%data(nr,ns)%up
          trmatrixx%data(nr,nd+6) = trmatrix_in%data(nr,ns)%uT
       enddo
    enddo
    CALL DATE_AND_TIME(trmatrixx%created_date, &
                       trmatrixx%created_time, &
                       trmatrixx%created_timezone)

    if(trmatrixx%status.ge.3) then
       trmatrixx%status=3
    else
       trmatrixx%status=2
    endif

    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_trmatrix'
       write(6,*) '---- trmatrixx%rho'
       write(6,'(1P5E12.4)') &
                     (trmatrixx%rho(nr),nr=1,trmatrixx%nrmax)
       do nd=1,trmatrixx%ndmax
          write(6,*) '---- ',trmatrixx%kid(nd)
          write(6,'(1P5E12.4)') &
                           (trmatrixx%data(nr,nd),nr=1,trmatrixx%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_put_trmatrix

!-----------------------------------------------------------------------
  subroutine bpsd_get_trmatrix(trmatrix_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_trmatrix_type),intent(inout) :: trmatrix_out
    integer,intent(out) :: ierr
    integer :: nr, nd, ns, mode
    real(dp), dimension(:), pointer :: v

    if(bpsd_trmatrixx_init_flag) call bpsd_init_trmatrixx

    if(trmatrixx%status.eq.0) then
       write(6,*) 'XX bpsd_get_trmatrix: no space allocated to trmatrixx%data'
       ierr=1
       return
    endif

    if(trmatrixx%status.eq.1) then
       write(6,*) 'XX bpsd_get_trmatrix: no data in trmatrixx%data'
       ierr=2
       return
    endif

    if(trmatrix_out%nrmax.eq.0) then
       mode=0
       trmatrix_out%nrmax = trmatrixx%nrmax
    else
       mode=1
    endif
    trmatrix_out%nsmax = (trmatrixx%ndmax-1)/6

    CALL bpsd_adjust_array1D(trmatrix_out%rho,trmatrix_out%nrmax)
    CALL bpsd_adjust_trmatrix_data(trmatrix_out%data,trmatrix_out%nrmax, &
                                                     trmatrix_out%nsmax)

    if(mode.eq.0) then
       trmatrix_out%time  = trmatrixx%time
       do nr=1,trmatrixx%nrmax
          trmatrix_out%rho(nr)=trmatrixx%rho(nr)
          do ns=1,trmatrix_out%nsmax
             nd=6*(ns-1)
             trmatrix_out%data(nr,ns)%Dn =trmatrixx%data(nr,nd+1)
             trmatrix_out%data(nr,ns)%Dp =trmatrixx%data(nr,nd+2)
             trmatrix_out%data(nr,ns)%DT =trmatrixx%data(nr,nd+3)
             trmatrix_out%data(nr,ns)%un =trmatrixx%data(nr,nd+4)
             trmatrix_out%data(nr,ns)%up =trmatrixx%data(nr,nd+5)
             trmatrix_out%data(nr,ns)%uT =trmatrixx%data(nr,nd+6)
          enddo
       enddo
       ierr=0
       return
    endif

    if(trmatrixx%status.eq.2) then
       CALL bpsd_adjust_array3D(trmatrixx%spline,6,trmatrixx%nrmax, &
                                                   trmatrixx%ndmax)
       trmatrixx%status=3
    endif

    if(trmatrixx%status.eq.3) then
       CALL bpsd_adjust_array1D(trmatrixx%s,trmatrixx%nrmax)
       do nr=1,trmatrixx%nrmax
          trmatrixx%s(nr)=trmatrixx%rho(nr)**2
       enddo
       do nd=1,trmatrixx%ndmax
          call bpsd_spl1D(trmatrixx,nd,ierr)
       enddo
       trmatrixx%status=4
    endif

    allocate(v(trmatrixx%ndmax))

    do nr=1,trmatrix_out%nrmax
       do nd=1,trmatrixx%ndmax
          call bpsd_spl1DF(trmatrix_out%rho(nr),v(nd),trmatrixx,nd,ierr)
       enddo
       do ns=1,trmatrix_out%nsmax
          nd=6*(ns-1)
          trmatrix_out%data(nr,ns)%Dn  = v(nd+1)
          trmatrix_out%data(nr,ns)%Dp  = v(nd+2)
          trmatrix_out%data(nr,ns)%DT  = v(nd+3)
          trmatrix_out%data(nr,ns)%un  = v(nd+4)
          trmatrix_out%data(nr,ns)%up  = v(nd+5)
          trmatrix_out%data(nr,ns)%uT  = v(nd+6)
       enddo
    enddo
    deallocate(v)
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_trmatrix'
       write(6,*) '---- trmatrixx%rho'
       write(6,'(1P5E12.4)') &
                     (trmatrix_out%rho(nr),nr=1,trmatrix_out%nrmax)
       do ns=1,trmatrix_out%nsmax
          write(6,*) '---- trmatrixx%Dn(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trmatrix_out%data(nr,ns)%Dn,nr=1,trmatrix_out%nrmax)
          write(6,*) '---- trmatrixx%Dp(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trmatrix_out%data(nr,ns)%Dp,nr=1,trmatrix_out%nrmax)
          write(6,*) '---- trmatrixx%DT(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trmatrix_out%data(nr,ns)%DT,nr=1,trmatrix_out%nrmax)
          write(6,*) '---- trmatrixx%un(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trmatrix_out%data(nr,ns)%un,nr=1,trmatrix_out%nrmax)
          write(6,*) '---- trmatrixx%up(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trmatrix_out%data(nr,ns)%up,nr=1,trmatrix_out%nrmax)
          write(6,*) '---- trmatrixx%uT(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trmatrix_out%data(nr,ns)%uT,nr=1,trmatrix_out%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_get_trmatrix

!-----------------------------------------------------------------------
  subroutine bpsd_save_trmatrix(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_trmatrixx_init_flag) call bpsd_init_trmatrixx

    if(trmatrixx%status.gt.1) call bpsd_save_data1Dx(fid,trmatrixx,ierr)
    return

  end subroutine bpsd_save_trmatrix

!-----------------------------------------------------------------------
  subroutine bpsd_load_trmatrix(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data1Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nr,nd

    if(bpsd_trmatrixx_init_flag) call bpsd_init_trmatrixx

    trmatrixx%dataName=datax%dataName
    trmatrixx%time = datax%time
    trmatrixx%nrmax=datax%nrmax
    trmatrixx%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(trmatrixx%kid,trmatrixx%ndmax)
    CALL bpsd_adjust_karray(trmatrixx%kunit,trmatrixx%ndmax)
    CALL bpsd_adjust_array1D(trmatrixx%rho,trmatrixx%nrmax)
    CALL bpsd_adjust_array2D(trmatrixx%data,trmatrixx%nrmax, &
                                            trmatrixx%ndmax)

    do nr=1,trmatrixx%nrmax
       trmatrixx%rho(nr) = datax%rho(nr)
       do nd=1,trmatrixx%ndmax
          trmatrixx%data(nr,nd) = datax%data(nr,nd)
       enddo
    enddo
    do nd=1,trmatrixx%ndmax
       trmatrixx%kid(nd)=datax%kid(nd)
       trmatrixx%kunit(nd)=datax%kunit(nd)
    enddo
    trmatrixx%created_date = datax%created_date
    trmatrixx%created_time = datax%created_time
    trmatrixx%created_timezone = datax%created_timezone

    if(trmatrixx%status.ge.3) then
       trmatrixx%status=3
    else
       trmatrixx%status=2
    endif

    ierr=0
    return

  end subroutine bpsd_load_trmatrix

end module bpsd_trmatrix
