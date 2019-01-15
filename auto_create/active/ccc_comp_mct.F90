!xmlinsert(list[0])
module {ccc}_comp_mct

!xmlinsert(list[1])

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: {ccc}_init_mct
  public :: {ccc}_run_mct
  public :: {ccc}_final_mct

!xmlinsert(list[2])

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: {ccc}_SetgsMap_mct
  private :: {ccc}_import_mct
  private :: {ccc}_export_mct
  private :: {ccc}_domain_mct
!xmlinsert(list[3])
!
!================================================================================
CONTAINS
!================================================================================
!xmlinsert(list[4])
!xmlinsert(list[5])
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)                 :: EClock
    type(seq_cdata), intent(inout)              :: cdata_{c}
!xmlinsert(list[6])
    type(mct_aVect), intent(inout)              :: x2{c}_{c}
    type(mct_aVect), intent(inout)              :: {c}2x_{c}   
    character(len=*), optional,   intent(IN)    :: NLFilename ! Namelist filename
    !
    ! Locals
    !
    type(mct_gsMap), pointer   :: gsMap_{ccc}
    type(mct_gGrid), pointer   :: dom_{c}
    integer :: lsize
!xmlinsert(list[7])
    !--------------------------------------------------------------------------
    ! Determine attribute vector indices
    !--------------------------------------------------------------------------

    call {ccc}_cpl_indices_set()
    
!xmlinsert(list[8])

       if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
!{list} key="runtypemap" value="0"
       else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
!{list} key="runtypemap" value="1"
       else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
!{list} key="runtypemap" value="2"
       else
          write(iulog,*) '{ccc}_comp_mct: ERROR: unknown starttype'
          call shr_sys_abort()
       end if
!xmlinsert(list[9])
       !
       ! Initialize MCT gsMap, domain and attribute vectors
       !
       call {ccc}_SetgsMap_mct( mpicom_{ccc}, {CCC}ID, gsMap_{ccc} )
       lsize = mct_gsMap_lsize(gsMap_{ccc}, mpicom_{ccc})
       !
       ! Initialize MCT domain 
       !
       call {ccc}_domain_mct( lsize, gsMap_{ccc}, dom_{c} )
       !
       ! Initialize MCT attribute vectors
       !
       call mct_aVect_init({c}2x_{c}, rList=seq_flds_{c}2x_fields, lsize=lsize)
       call mct_aVect_zero({c}2x_{c})
       
       call mct_aVect_init(x2{c}_{c}, rList=seq_flds_x2{c}_fields, lsize=lsize) 
       call mct_aVect_zero(x2{c}_{c})
!xmlinsert(list[10])
       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
!xmlinsert(list[11])

 end subroutine {ccc}_init_mct

!================================================================================
!xmlinsert(list[12])
!xmlinsert(list[13])
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_{c}
    type(mct_aVect)             , intent(inout) :: x2{c}_{c}
    type(mct_aVect)             , intent(inout) :: {c}2x_{c}

!xmlinsert(list[14])
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

!xmlinsert(list[15])

!xmlinsert(list[16])

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)

!xmlinsert(list[17])

#if (defined _MEMTRACE)
    if(my_task == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':end::',lbnum) 
       call memmon_reset_addr()
    endif
#endif

  end subroutine {ccc}_run_mct

  subroutine {ccc}_final_mct( )
!xmlinsert(list[18])    
!xmlinsert(list[19])
  end subroutine {ccc}_final_mct

  subroutine {ccc}_SetGSMap_mct( mpicom_{ccc}, {CCC}ID, gsMap_{ccc} )

!xmlinsert(list[20])
    integer        , intent(in)    :: mpicom_{ccc}
    integer        , intent(in)    :: {CCC}ID
    type(mct_gsMap), intent(inout) :: gsMap_{ccc}
    integer, allocatable :: gindex(:)
!xmlinsert(list[21])
    call mct_gsMap_init( gsMap_{ccc}, gindex, mpicom_{ccc}, {CCC}ID, nlcols, ngcols)
    deallocate(gindex)
  end subroutine atm_SetgsMap_mct

!xmlinsert(list[22])  
    type(mct_aVect)   , intent(inout) :: x2{c}_{c}
!xmlinsert(list[23])
 end subroutine {ccc}_import_mct

!xmlinsert(list[export_mct])
    type(mct_aVect)   , intent(inout) :: {c}2x_{c}


!<list> key="export"
    {c}2x_{c}%rAttr(index_{c}2x_{x1}_{x2}   ,{x}) = {ccc_out}%{x3}(i){x4}
!</list>

!xmlinsert(list[24])
 end subroutine {ccc}_export_mct

 subroutine {ccc}_domain_mct( lsize, gsMap_{c}, dom_{c})
!xmlinsert(list[domain_use])
    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_{c}
    type(mct_ggrid), intent(inout) :: dom_{c}
!xmlinsert(list[domain_local_var])
    call mct_gGrid_init( GGrid=dom_{c}, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )
!xmlinsert(list[19_mct_aVect_zero])
    allocate(data(lsize))
    call mct_gsMap_orderedPoints(gsMap_{c}, {my_task}, idata)
    call mct_gGrid_importIAttr(dom_{c},'GlobGridNum',idata,lsize)

    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_{c},"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_{c},"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_{c},"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_{c},"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_{c},"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_{c},"frac",data,lsize)
!xmlinsert(list[20_domain_computing])
    deallocate(data)
    deallocate(idata)
!xmlinsert(list[deallocate_local])
  end subroutine {ccc}_domain_mct
!xmlinsert(list[21_other_function])
end module {ccc}_comp_mct
