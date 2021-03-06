module map_iceocn_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of ocn-ice.
!       
!
! Author: R. Jacob, M. Vertenstein
! Revision History:
! 30Mar06 - P. Worley - added optional arguments to MCT_Rearrange call
! 13Apr06 - M. Vertenstein - cleaned up interfaces 
!
!---------------------------------------------------------------------



  use shr_sys_mod
  use mct_mod
  use seq_cdata_mod
  use seq_comm_mct

  implicit none
  save
  private  ! except


!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_ocn2ice_init_mct
  public :: map_ice2ocn_init_mct
  public :: map_ocn2ice_mct
  public :: map_ice2ocn_mct



!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_ice2ocn
  type(mct_rearr), private :: Re_ocn2ice

#ifdef CPP_VECTOR
    logical :: usevector = .true.
#else
    logical :: usevector = .false.
#endif

#ifdef SYSUNICOS
    logical :: usealltoall = .true.
#else
    logical :: usealltoall = .false.
#endif

  character(*),parameter :: subName = '(map_iceocn_mct)'






!=======================================================================
   contains
!=======================================================================
  subroutine map_ice2ocn_init_mct( cdata_i, cdata_o)

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_i
    type(seq_cdata),intent(in) :: cdata_o

    ! 
    ! Local variables
    !

    integer                  :: icesize, ocnsize  ! global grid sizes

    type(mct_gsMap), pointer :: gsMap_i           ! ice gsMap
    type(mct_gsMap), pointer :: gsMap_o           ! ocn gsMap

    integer                  :: mpicom            ! communicator spanning atm and lnd



    character(*),parameter :: subName = '(map_ice2ocn_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_i, gsMap=gsMap_i)
    call seq_cdata_setptrs(cdata_o, gsMap=gsMap_o)

    call seq_cdata_setptrs(cdata_o, mpicom=mpicom)

    icesize = mct_gsMap_gsize(gsMap_i) 
    ocnsize = mct_gsMap_gsize(gsMap_o) 
    if (icesize /= ocnsize) then 
      write(logunit,*) "(map_ice2ocn_init_mct) ice and ocn are different." 
      write(logunit,*) "(map_ice2ocn_init_mct) Must be same size. Exiting." 
      call shr_sys_abort(subName // "different size") 
    endif 







    call mct_rearr_init(gsMap_i, gsMap_o, mpicom, Re_ice2ocn)



  end subroutine map_ice2ocn_init_mct

!=======================================================================
  subroutine map_ocn2ice_init_mct( cdata_o, cdata_i )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_o
    type(seq_cdata),intent(in) :: cdata_i
    ! 
    ! Local variables
    !

    integer                  :: icesize, ocnsize  ! global grid sizes

    type(mct_gsMap), pointer :: gsMap_i           ! ice gsMap
    type(mct_gsMap), pointer :: gsMap_o           ! ocn gsMap
    type(mct_gGrid), pointer :: dom_o             ! ocn domain
    type(mct_gGrid), pointer :: dom_i             ! ice domain
    integer                  :: mpicom            ! communicator spanning ice and ocn
    character(*),parameter :: subName = '(map_ocn2ice_init_mct) '

    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! ocn areas from mapping file
    type(mct_aVect)          :: areadst           ! atm areas from mapping file


    !--------------------------------------------------
    call seq_cdata_setptrs(cdata_o, gsMap=gsMap_o, dom=dom_o)
    call seq_cdata_setptrs(cdata_i, gsMap=gsMap_i, dom=dom_i)





    call seq_cdata_setptrs(cdata_o, mpicom=mpicom)




    ! Initialize ocn -> ice mapping or rearranger
    icesize = mct_gsMap_gsize(gsMap_i)
    ocnsize = mct_gsMap_gsize(gsMap_o)
    if (icesize /= ocnsize) then
      write(logunit,*) "(map_ocn2ice_init_mct) ocn and ice are different."
      write(logunit,*) "(map_ocn2ice_init_mct) Must be same size. Exiting."
      call shr_sys_abort(subName // "different size")
    endif

    ! Initialize rearranger

    call mct_rearr_init(gsMap_o, gsMap_i, mpicom, Re_ocn2ice)




    ! Set ice aream to ocn aream
    ! Note that "aream" attribute of dom_o%data is set in routine
    ! map_ocn2ice_init_mct

    lsize = mct_gsMap_lsize(gsMap_o, mpicom)
    call mct_aVect_init(areasrc, rList="aream", lsize=lsize )

    lsize = mct_gsMap_lsize(gsMap_i, mpicom)
    call mct_aVect_init(areadst, rList="aream", lsize=lsize )
    call mct_aVect_zero(areadst)

    km = mct_aVect_indexRA(dom_o%data,"aream")
    ka = mct_aVect_indexRA(areasrc   ,"aream")
    areasrc%rAttr(ka,:) = dom_o%data%rAttr(km,:)

    call mct_rearr_rearrange(areasrc, areadst, Re_ocn2ice, VECTOR=usevector, ALLTOALL=usealltoall)

    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_i%data,"aream")
    dom_i%data%rAttr(km,:) = areadst%rAttr(ka,:)

    call mct_aVect_clean(areasrc)
    call mct_aVect_clean(areadst)

  end subroutine map_ocn2ice_init_mct

!=======================================================================
  subroutine map_ice2ocn_mct( cdata_i, av_i, cdata_o, av_o, fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata_i
    type(mct_aVect), intent(in)  :: av_i
    type(seq_cdata), intent(in)  :: cdata_o
    type(mct_aVect), intent(out) :: av_o
    character(len=*),intent(in), optional :: statelist
    character(len=*),intent(in), optional :: fluxlist








    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_i, av_o, Re_ice2ocn, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_i, av_o, Re_ice2ocn, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_i, av_o, Re_ice2ocn, VECTOR=usevector, ALLTOALL=usealltoall)
    end if




  end subroutine map_ice2ocn_mct

!=======================================================================
  subroutine map_ocn2ice_mct( cdata_o, av_o, cdata_i, av_i, &
                              fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata) ,intent(in)  :: cdata_o
    type(mct_aVect) ,intent(in)  :: av_o
    type(seq_cdata) ,intent(in)  :: cdata_i
    type(mct_aVect) ,intent(out) :: av_i

    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist



    character(*),parameter :: subName = '(map_ocn2ice_mct) '



    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_o, av_i, Re_ocn2ice, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_o, av_i, Re_ocn2ice, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_o, av_i, Re_ocn2ice, VECTOR=usevector, ALLTOALL=usealltoall)
    end if


  end subroutine map_ocn2ice_mct


end module map_iceocn_mct
