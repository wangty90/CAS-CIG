module map_rofocn_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of ocn-rof.
!       
!
! Author: R. Jacob, M. Vertenstein
! Revision History:
! 30Mar06 - P. Worley - added optional arguments to MCT_Rearrange call
! 13Apr06 - M. Vertenstein - cleaned up interfaces 
!
!---------------------------------------------------------------------




  use shr_sys_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod
  use seq_cdata_mod
  use seq_infodata_mod

  implicit none

  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------


  public :: map_rof2ocn_init_mct
  public :: map_rof2ocn_mct


!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------


  type(mct_rearr), private :: Re_rof2ocn
  type(mct_sMatp), private :: sMatp_Fr2o

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

  logical, private :: samegrid_mapr2o

  character(*),parameter :: subName = '(map_rofocn_mct)'





!=======================================================================
   contains
!=======================================================================
  subroutine map_rof2ocn_init_mct( cdata_r, cdata_o)

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_r
    type(seq_cdata),intent(in) :: cdata_o

    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata


    type(mct_gsMap), pointer :: gsMap_r           ! rof gsMap
    type(mct_gsMap), pointer :: gsMap_o           ! ocn gsMap
    type(mct_gGrid), pointer :: dom_r             ! rof domain
    type(mct_gGrid), pointer :: dom_o             ! ocn domain

    integer                  :: mpicom            ! communicator spanning atm and lnd

    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! atm areas from mapping file
    type(mct_aVect)          :: areadst           ! lnd areas set to atm areas
    type(mct_rearr)          :: Re_ocn2rof


    character(*),parameter :: subName = '(map_rof2ocn_init_mct) '
    !--------------------------------------------------
    call seq_cdata_setptrs(cdata_r, gsMap=gsMap_r, dom=dom_r)
    call seq_cdata_setptrs(cdata_o, gsMap=gsMap_o, dom=dom_o)


    call seq_cdata_setptrs(cdata_r, mpicom=mpicom, infodata=infodata)
    call seq_infodata_GetData( infodata, samegrid_ro=samegrid_mapr2o)






    if (samegrid_mapr2o) then

       call mct_rearr_init(gsMap_r, gsMap_o, mpicom, Re_rof2ocn)
       call mct_rearr_init(gsMap_o, gsMap_r, mpicom, Re_ocn2rof)

       call mct_rearr_rearrange_fldlist(dom_o%data, dom_r%data, Re_ocn2rof, VECTOR=usevector, &
          ALLTOALL=usealltoall, fldlist='aream')

       call mct_rearr_clean(Re_ocn2rof)


    else

       ! Initialize rof->ocn mapping or rearranging

       lsize = mct_gsMap_lsize(gsMap_r, mpicom)
       call mct_aVect_init( areasrc, rList="aream", lsize=lsize )

       call shr_mct_sMatPInitnc(sMatp_Fr2o, gsMap_r, gsMap_o, "seq_maps.rc", &
              "rof2ocnFmapname:","rof2ocnFmaptype:",mpicom, &
               areasrc=areasrc)

      ! Determine rof grid areas from mapping files

       km = mct_aVect_indexRA(dom_r%data,"aream", perrWith=subName)
       ka = mct_aVect_indexRA(areasrc   ,"aream", perrWith=subName)
       dom_r%data%rAttr(km,:) = areasrc%rAttr(ka,:)

       call mct_aVect_clean(areasrc)

    endif






  end subroutine map_rof2ocn_init_mct

!=======================================================================

!=======================================================================

  subroutine map_rof2ocn_mct( cdata_r, r2x_r, cdata_o, r2x_o)

    type(seq_cdata),intent(in) :: cdata_r
    type(mct_aVect),intent(in) :: r2x_r
    type(seq_cdata),intent(in) :: cdata_o
    type(mct_aVect),intent(out):: r2x_o









    if (samegrid_mapr2o) then

          call mct_rearr_rearrange(r2x_r, r2x_o, Re_rof2ocn, VECTOR=usevector, &
             ALLTOALL=usealltoall)

    else

!tcx    call mct_aVect_zero(r2x_o)
       call mct_sMat_avMult(r2x_r, sMatp_Fr2o, r2x_o, VECTOR=usevector)

    endif




  end subroutine map_rof2ocn_mct

!=======================================================================


end module map_rofocn_mct
