module map_{ccc1}{ccc2}_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of {ccc2}-{ccc1}.
!       
!
! Author: R. Jacob, M. Vertenstein
! Revision History:
! 30Mar06 - P. Worley - added optional arguments to MCT_Rearrange call
! 13Apr06 - M. Vertenstein - cleaned up interfaces 
!
!---------------------------------------------------------------------
!<list> key="use" number=1 
  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod

  use seq_comm_mct, only : logunit, loglevel
  use seq_cdata_mod
  use seq_flds_indices
  use seq_infodata_mod
  use seq_map_mod
  use m_die

  implicit none
  save
  private  ! except
!</list>

!<list> key="use" number=2
  use shr_kind_mod     , only: R8 => SHR_KIND_R8
  use shr_sys_mod

  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_cdata_mod
  use seq_comm_mct, only: logunit, loglevel
  use seq_infodata_mod
  use seq_map_mod
  use m_die

  implicit none
  save
  private  ! except
!</list>

!<list> key="use" number=3
  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod      ,only: CL => SHR_KIND_CL
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod

  use seq_comm_mct, only : logunit, loglevel
  use seq_cdata_mod
  use seq_flds_indices
  use seq_infodata_mod
  use seq_map_mod
  use m_die

  implicit none
  save
  private  ! except
!</list>

!<list> key="use" number=4
  use shr_sys_mod
  use mct_mod
  use seq_cdata_mod
  use seq_comm_mct

  implicit none
  save
  private  ! except
!</list>

!<list> key="use" number=5
  use shr_sys_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod
  use seq_cdata_mod
  use seq_infodata_mod

  implicit none

  private  ! except
!</list>

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------
!<list> key="public_interface" number=1

  public :: map_{ccc2}2{ccc1}_init_mct
  public :: map_{ccc1}2{ccc2}_init_mct
  public :: map_{ccc2}2{ccc1}_mct
  public :: map_{ccc1}2{ccc2}_mct

!</list>

!<list> key="public_interface" number=2

  public :: map_{ccc1}2{ccc2}_init_mct
  public :: map_{ccc1}2{ccc2}_mct

!</list>

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------
!<list> key="private_data" number=1
  type(mct_rearr), private :: Re_{ccc2}2{ccc1}
  type(mct_rearr), private :: Re_{ccc1}2{ccc2}
  type(mct_sMatp), private :: sMatp_F{c1}2{c2}
  type(mct_sMatp), private :: sMatp_S{c1}2{c2}
  type(mct_sMatp), private :: sMatp_F{c2}2{c1}
  type(mct_sMatp), private :: sMatp_S{c2}2{c1}

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
  logical, private :: samegrid_map{c1}2{c2}
!</list>

!<list> key="private_data" number=2
  type(mct_rearr), private :: Re_{ccc1}2{ccc2}
  type(mct_rearr), private :: Re_{ccc2}2{ccc1}

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

  character(*),parameter :: subName = '(map_{ccc1}{ccc2}_mct)'

!</list>

!<list> key="private_data" number=3
  type(mct_rearr), private :: Re_{ccc1}2{ccc2}
  type(mct_sMatp), private :: sMatp_F{c1}2{c2}

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

  logical, private :: samegrid_map{c1}2{c2}

  character(*),parameter :: subName = '(map_{ccc1}{ccc2}_mct)'

!</list>

!<list> key="vector_mapping" number=1
  !--- for vector mapping ---
  real(R8), private, allocatable :: slon_{c1}(:),clon_{c1}(:),slat_{c1}(:),clat_{c1}(:)
  real(R8), private, allocatable :: slon_{c2}(:),clon_{c2}(:),slat_{c2}(:),clat_{c2}(:)
  type(mct_aVect), private :: av3d_{c1}, av3d_{c2}
!</list>

!<list> key="deg2rad" number=1
  real(R8),parameter,private :: shr_const_deg2rad = shr_const_pi/180.0_R8  ! deg to rads
!</list>

!<list> key="iogrid" number=1
  integer         :: ni_a   ! number of longitudes on input grid
  integer         :: nj_a   ! number of latitudes  on input grid
  integer         :: ni_o   ! number of longitudes on output grid
  integer         :: nj_o   ! number of latitudes  on output grid
!</list>

!=======================================================================
   contains
!=======================================================================
!<list> key="init1" number=1
  subroutine map_{ccc1}2{ccc2}_init_mct( cdata_{c1}, cdata_{c2})

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_{c1}
    type(seq_cdata),intent(in) :: cdata_{c2}
!</list>

    ! 
    ! Local variables
    !
!{list} key="infodata" number=1
!    type(seq_infodata_type), pointer :: infodata
!

!{list} key="global_grid_sizes" number=1
!    integer                  :: {ccc1}size, {ccc2}size  ! global grid sizes
!

    type(mct_gsMap), pointer :: gsMap_{c1}           ! {ccc1} gsMap
    type(mct_gsMap), pointer :: gsMap_{c2}           ! {ccc2} gsMap
!<list> key="dom" number=1
    type(mct_gGrid), pointer :: dom_{c1}             ! {ccc1} domain
    type(mct_gGrid), pointer :: dom_{c2}             ! {ccc2} domain
!</list>

    integer                  :: mpicom            ! communicator spanning atm and lnd
!<list> key="areasrc_dst" number=1
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! atm areas from mapping file
    type(mct_aVect)          :: areadst           ! lnd areas set to atm areas
!</list>

!<list> key="areasrc_dst" number=2
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! atm areas from mapping file
    type(mct_aVect)          :: areadst           ! lnd areas set to atm areas
    type(mct_rearr)          :: Re_{ccc2}2{ccc1}
!</list>

!<list> key="grids"  number=1
    integer                  :: klon,klat,n
!</list>

    character(*),parameter :: subName = '(map_{ccc1}2{ccc2}_init_mct) '
    !--------------------------------------------------
!<list> key="dom" number=1
    call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1}, dom=dom_{c1})
    call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2}, dom=dom_{c2})
!</list>

!<list> key="dom" number=2
    call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1})
    call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2})
!</list>

!{list} key="infodata" number=1
!    call seq_cdata_setptrs(cdata_{c1}, mpicom=mpicom, infodata=infodata)

!{list} key="infodata" number=2
!    call seq_cdata_setptrs(cdata_{c1}, mpicom=mpicom)

!{list} key="infodata" number=3
!    call seq_cdata_setptrs(cdata_{c2}, mpicom=mpicom)

!{list} key="infodata" number=1
!    call seq_infodata_GetData( infodata, samegrid_{c1}{c2}=samegrid_map{c1}2{c2})
!

!<list> key="without_samegrid" number=1
    {ccc1}size = mct_gsMap_gsize(gsMap_{c1}) 
    {ccc2}size = mct_gsMap_gsize(gsMap_{c2}) 
    if ({ccc1}size /= {ccc2}size) then 
      write(logunit,*) "(map_{ccc1}2{ccc2}_init_mct) {ccc1} and {ccc2} are different." 
      write(logunit,*) "(map_{ccc1}2{ccc2}_init_mct) Must be same size. Exiting." 
      call shr_sys_abort(subName // "different size") 
    endif 
!</list>

!<list> key="before_samegrid" number=1
    lsize = mct_gsMap_lsize(gsMap_{c1}, mpicom)
    call mct_aVect_init( areasrc, rList="aream", lsize=lsize )      
    lsize = mct_gsMap_lsize(gsMap_{c2}, mpicom)
    call mct_aVect_init( areadst, rList="aream", lsize=lsize )
!</list>

!<list> key="samegrid" number=1
    if (samegrid_map{c1}2{c2}) then
       
       call mct_rearr_init(gsMap_{c1}, gsMap_{c2}, mpicom, Re_{ccc1}2{ccc2})
       
       ! --- copy {ccc1} area into {ccc2} aream
 
       ka = mct_aVect_indexRA(areasrc   , "aream")
       km = mct_aVect_indexRa(dom_{c1}%data, "aream" )
       areasrc%rAttr(ka,:) = dom_{c1}%data%rAttr(km,:)
       call mct_rearr_rearrange(areasrc, areadst, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
          ALLTOALL=usealltoall)
    else

       call shr_mct_sMatPInitnc(sMatp_F{c1}2{c2},gsMap_{c1},gsMap_{c2},"seq_maps.rc", &
          "{ccc1}2{ccc2}Fmapname:","{ccc1}2{ccc2}Fmaptype:",mpicom, &
          areasrc=areasrc, areadst=areadst)

       call shr_mct_sMatPInitnc(sMatp_S{c1}2{c2},gsMap_{c1},gsMap_{c2},"seq_maps.rc", &
          "{ccc1}2{ccc2}Smapname:","{ccc1}2{ccc2}Smaptype:",mpicom)

    endif
!</list>

!<list> key="samegrid" number=2
    if (samegrid_map{c1}2{c2}) then

       call mct_rearr_init(gsMap_{c1}, gsMap_{c2}, mpicom, Re_{ccc1}2{ccc2})

    else

       call shr_mct_sMatPInitnc(sMatp_F{c1}2{c2},gsMap_{c1},gsMap_{c2},"seq_maps.rc", &
            "{ccc1}2{ccc2}Fmapname:","{ccc1}2{ccc2}Fmaptype:",mpicom)
       call shr_mct_sMatPInitnc(sMatp_S{c1}2{c2},gsMap_{c1},gsMap_{c2},"seq_maps.rc", &
            "{ccc1}2{ccc2}Smapname:","{ccc1}2{ccc2}Smaptype:",mpicom)

    endif
!</list>

!<list> key="samegrid" number=3
    if (samegrid_map{c1}2{c2}) then

       call mct_rearr_init(gsMap_{c1}, gsMap_{c2}, mpicom, Re_{ccc1}2{ccc2})

       !--- want to copy atm area into atm and ocn aream ---

       lsize = mct_gsMap_lsize(gsMap_{c1}, mpicom)
       call mct_aVect_init( areasrc, rList="aream", lsize=lsize )

       lsize = mct_gsMap_lsize(gsMap_{c2}, mpicom)
       call mct_aVect_init( areadst, rList="aream", lsize=lsize )

       ka = mct_aVect_indexRa(dom_{c1}%data, "area" )
       km = mct_aVect_indexRA(areasrc   , "aream")
       areasrc%rAttr(km,:) = dom_{c1}%data%rAttr(ka,:)

       call mct_rearr_rearrange(areasrc, areadst, Re_{ccc1}2{ccc2}, VECTOR=usevector,&
          ALLTOALL=usealltoall)

       ka = mct_aVect_indexRA(areasrc   ,"aream")
       km = mct_aVect_indexRA(dom_{c1}%data,"aream")
       dom_{c1}%data%rAttr(km,:) = areasrc%rAttr(ka,:)

       ka = mct_aVect_indexRA(areadst   ,"aream")
       km = mct_aVect_indexRA(dom_{c2}%data,"aream")
       dom_{c2}%data%rAttr(km,:) = areadst%rAttr(ka,:)

       call mct_aVect_clean(areasrc)
       call mct_aVect_clean(areadst)

    else

       call shr_mct_sMatPInitnc(sMatp_F{c1}2{c2},gsMap_{c1},gsMap_{c2},"seq_maps.rc", &
          "{ccc1}2{ccc2}Fmapname:","{ccc1}2{ccc2}Fmaptype:",mpicom,&
          ni_i=ni_{c1},nj_i=nj_{c1},ni_o=ni_{c2},nj_o=nj_{c2})

       call shr_mct_sMatPInitnc(sMatp_S{c1}2{c2},gsMap_{c1},gsMap_{c2},"seq_maps.rc", &
          "{ccc1}2{ccc2}Smapname:","{ccc1}2{ccc2}Smaptype:",mpicom)

       !--- compute these up front for vector mapping ---
       lsize = mct_gsMap_lsize(gsMap_{c1}, mpicom)
       call mct_avect_init(av3d_{c1},rList='ux:uy:uz',lsize=lsize)
       allocate(slon_{c1}(lsize),clon_{c1}(lsize),slat_{c1}(lsize),clat_{c1}(lsize))
       klon = mct_aVect_indexRa(dom_{c1}%data, "lon" )
       klat = mct_aVect_indexRa(dom_{c1}%data, "lat" )
       do n = 1,lsize
          slon_{c1}(n) = sin(dom_{c1}%data%rAttr(klon,n)*shr_const_deg2rad)
          clon_{c1}(n) = cos(dom_{c1}%data%rAttr(klon,n)*shr_const_deg2rad)
          slat_{c1}(n) = sin(dom_{c1}%data%rAttr(klat,n)*shr_const_deg2rad)
          clat_{c1}(n) = cos(dom_{c1}%data%rAttr(klat,n)*shr_const_deg2rad)
       enddo

       lsize = mct_gsMap_lsize(gsMap_{c2}, mpicom)
       call mct_avect_init(av3d_{c2},rList='ux:uy:uz',lsize=lsize)
       allocate(slon_{c2}(lsize),clon_{c2}(lsize),slat_{c2}(lsize),clat_{c2}(lsize))
       klon = mct_aVect_indexRa(dom_{c2}%data, "lon" )
       klat = mct_aVect_indexRa(dom_{c2}%data, "lat" )
       do n = 1,lsize
          slon_o(n) = sin(dom_{c2}%data%rAttr(klon,n)*shr_const_deg2rad)
          clon_o(n) = cos(dom_{c2}%data%rAttr(klon,n)*shr_const_deg2rad)
          slat_o(n) = sin(dom_{c2}%data%rAttr(klat,n)*shr_const_deg2rad)
          clat_o(n) = cos(dom_{c2}%data%rAttr(klat,n)*shr_const_deg2rad)
       enddo

    endif

!</list>

!<list> key="samegrid" number=4
    if (samegrid_map{c1}2{c2}) then

       call mct_rearr_init(gsMap_{c1}, gsMap_{c2}, mpicom, Re_{ccc1}2{ccc2})
       call mct_rearr_init(gsMap_{c2}, gsMap_{c1}, mpicom, Re_{ccc2}2{ccc1})

       call mct_rearr_rearrange_fldlist(dom_{c2}%data, dom_{c1}%data, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
          ALLTOALL=usealltoall, fldlist='aream')

       call mct_rearr_clean(Re_{ccc2}2{ccc1})


    else

       ! Initialize rof->ocn mapping or rearranging

       lsize = mct_gsMap_lsize(gsMap_{c1}, mpicom)
       call mct_aVect_init( areasrc, rList="aream", lsize=lsize )

       call shr_mct_sMatPInitnc(sMatp_F{c1}2{c2}, gsMap_{c1}, gsMap_{c2}, "seq_maps.rc", &
              "{ccc1}2{ccc2}Fmapname:","{ccc1}2{ccc2}Fmaptype:",mpicom, &
               areasrc=areasrc)

      ! Determine rof grid areas from mapping files

       km = mct_aVect_indexRA(dom_{c1}%data,"aream", perrWith=subName)
       ka = mct_aVect_indexRA(areasrc   ,"aream", perrWith=subName)
       dom_{c1}%data%rAttr(km,:) = areasrc%rAttr(ka,:)

       call mct_aVect_clean(areasrc)

    endif

!</list>

!<list> key="samegrid" number=5
    km = mct_aVect_indexRA(dom_{c1}%data,"aream")
    ka = mct_aVect_indexRA(areasrc   ,"aream")
    areasrc%rAttr(ka,:) = dom_{c1}%data%rAttr(km,:)

    ! tcraig, 6/2010, initialize areadst with dom_g km because areasrc is only
    ! over land
    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_{c2}%data,"aream")
    areadst%rAttr(ka,:) = dom_{c2}%data%rAttr(km,:)

    call mct_rearr_rearrange(areasrc, areadst, Re_{ccc1}2{ccc2}, VECTOR=usevector, ALLTOALL=usealltoall)
!</list>

!<list> key="init_rearranger" number=1
    call mct_rearr_init(gsMap_{c1}, gsMap_{c2}, mpicom, Re_{ccc1}2{ccc2})
!</list>


!<list> key="before_samegrid" number=1
    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_{c2}%data,"aream")
    dom_{c2}%data%rAttr(km,:) = areadst%rAttr(ka,:)
    call mct_aVect_clean(areasrc)      
    call mct_aVect_clean(areadst)      
!</list>

  end subroutine map_{ccc1}2{ccc2}_init_mct

!=======================================================================
!<list> key="public_interface" number=1
  subroutine map_{ccc2}2{ccc1}_init_mct( cdata_{c2}, cdata_{c1} )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_{c2}
    type(seq_cdata),intent(in) :: cdata_{c1}
    ! 
    ! Local variables
    !
!<list> key="infodata" number=1
    type(seq_infodata_type), pointer :: infodata
!</list>

!{list} key="global_grid_sizes" number=1
!    integer                  :: {ccc1}size, {ccc2}size  ! global grid sizes
!

    type(mct_gsMap), pointer :: gsMap_{c1}           ! {ccc1} gsMap
    type(mct_gsMap), pointer :: gsMap_{c2}           ! {ccc2} gsMap
    type(mct_gGrid), pointer :: dom_{c2}             ! {ccc2} domain
    type(mct_gGrid), pointer :: dom_{c1}             ! {ccc1} domain
    integer                  :: mpicom            ! communicator spanning {ccc1} and {ccc2}
    character(*),parameter :: subName = '(map_{ccc2}2{ccc1}_init_mct) '
!<list> key="init2" number=1
    integer                  :: kf,iam,ierr,lsize
!</list>

!<list> key="init2" number=2
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! ocn areas from mapping file
    type(mct_aVect)          :: areadst           ! atm areas from mapping file
!</list>

!<list> key="mapr2o" number=1
  logical, private :: samegrid_mapr2o
!</list>

    !--------------------------------------------------
!<list> key="dom2" number=1
    call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2}, dom=dom_{c2})
    call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1}, dom=dom_{c1})
!</list>

!<list> key="dom2" number=2
    call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2})
    call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1})
!</list>

!<list> key="dom2" number=3
    call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2}, dom=dom_{c2})
    call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1})
!</list>

!<list> key="infodata" number=1
    call seq_cdata_setptrs(cdata_{c1}, mpicom=mpicom,infodata=infodata)

!</list>

!<list> key="infodata" number=2
    call seq_cdata_setptrs(cdata_{c1}, mpicom=mpicom)

!</list>

!<list> key="infodata" number=3
    call seq_cdata_setptrs(cdata_{c2}, mpicom=mpicom)

!</list>

!<list> key="init2" number=1
    call mpi_comm_rank(mpicom,iam,ierr)
!</list>

!<list> key="infodata" number=1    
    call seq_infodata_GetData( infodata, samegrid_{c1}{c2}=samegrid_map{c1}2{c2})
!</list>

    ! Initialize {ccc2} -> {ccc1} mapping or rearranger
!<list> key="without_samegrid" number=1
    {ccc1}size = mct_gsMap_gsize(gsMap_{c1})
    {ccc2}size = mct_gsMap_gsize(gsMap_{c2})
    if ({ccc1}size /= {ccc2}size) then
      write(logunit,*) "(map_{ccc2}2{ccc1}_init_mct) {ccc2} and {ccc1} are different."
      write(logunit,*) "(map_{ccc2}2{ccc1}_init_mct) Must be same size. Exiting."
      call shr_sys_abort(subName // "different size")
    endif
!</list>

!<list> key="init_rearranger" number=1
    ! Initialize rearranger

    call mct_rearr_init(gsMap_{c2}, gsMap_{c1}, mpicom, Re_{ccc2}2{ccc1})
!</list>

!<list> key="samegrid" number=1
    if (samegrid_map{c1}2{c2}) then

       call mct_rearr_init(gsMap_{c2}, gsMap_{c1}, mpicom, Re_{ccc2}2{ccc1})
    else
       call shr_mct_sMatPInitnc(sMatp_F{c2}2{c1}, gsMap_{c2}, gsMap_{c1}, "seq_maps.rc", &
            "{ccc2}2{ccc1}Fmapname:", "{ccc2}2{ccc1}Fmaptype:", mpicom)

       call shr_mct_sMatPInitnc(sMatp_S{c2}2{c1}, gsMap_{c2}, gsMap_{c1}, "seq_maps.rc", &
            "{ccc2}2{ccc1}Smapname:", "{ccc2}2{ccc1}Smaptype:", mpicom)

    endif
!</list>

!<list> key="samegrid" number=2
    if (samegrid_map{c1}2{c2}) then

       call mct_rearr_init(gsMap_{c2}, gsMap_{c1}, mpicom, Re_{ccc2}2{ccc1})
    else
       call shr_mct_sMatPInitnc(sMatp_F{c2}2{c1}, gsMap_{c2}, gsMap_{c1}, "seq_maps.rc", &
            "{ccc2}2{ccc1}Fmapname:", "{ccc2}2{ccc1}Fmaptype:", mpicom)

       call shr_mct_sMatPInitnc(sMatp_S{c2}2{c1}, gsMap_{c2}, gsMap_{c1}, "seq_maps.rc", &
            "{ccc2}2{ccc1}Smapname:", "{ccc2}2{ccc1}Smaptype:", mpicom)

    endif
!</list>

!<list> key="samegrid" number=3
    ! Initialize {ccc2}->{ccc1} mapping or rearranging

    if (samegrid_map{c1}2{c2}) then

       call mct_rearr_init(gsMap_{c2}, gsMap_{c1}, mpicom, Re_{ccc2}2{ccc1})

    else

       lsize = mct_gsMap_lsize(gsMap_{c2}, mpicom)
       call mct_aVect_init( areasrc, rList="aream", lsize=lsize )

       lsize = mct_gsMap_lsize(gsMap_{c1}, mpicom)
       call mct_aVect_init( areadst, rList="aream", lsize=lsize )

       call shr_mct_sMatPInitnc(sMatp_F{c2}2{c1}, gsMap_{c2}, gsMap_{c1}, "seq_maps.rc", &
            "{ccc2}2{ccc1}Fmapname:", "{ccc2}2{ccc1}Fmaptype:", mpicom, &
            areasrc=areasrc, areadst=areadst)

       call shr_mct_sMatPInitnc(sMatp_S{c2}2{c1}, gsMap_{c2}, gsMap_{c1}, "seq_maps.rc", &
            "{ccc2}2{ccc1}Smapname:", "{ccc2}2{ccc1}Smaptype:", mpicom)

       !--- copy aream from mapping files

       km = mct_aVect_indexRA(dom_{c2}%data,"aream")
       ka = mct_aVect_indexRA(areasrc   ,"aream")
       dom_{c2}%data%rAttr(km,:) = areasrc%rAttr(ka,:)

       km = mct_aVect_indexRA(dom_{c1}%data,"aream")
       ka = mct_aVect_indexRA(areadst   ,"aream")
       dom_{c1}%data%rAttr(km,:) = areadst%rAttr(ka,:)

       call mct_aVect_clean(areasrc)
       call mct_aVect_clean(areadst)

    endif
!</list>

!<list> key="samegrid" number=6
    ! Set {ccc1} aream to {ccc2} aream
    ! Note that "aream" attribute of dom_{c2}%data is set in routine
    ! map_{ccc2}2{ccc1}_init_mct

    lsize = mct_gsMap_lsize(gsMap_{c2}, mpicom)
    call mct_aVect_init(areasrc, rList="aream", lsize=lsize )

    lsize = mct_gsMap_lsize(gsMap_{c1}, mpicom)
    call mct_aVect_init(areadst, rList="aream", lsize=lsize )
    call mct_aVect_zero(areadst)

    km = mct_aVect_indexRA(dom_{c2}%data,"aream")
    ka = mct_aVect_indexRA(areasrc   ,"aream")
    areasrc%rAttr(ka,:) = dom_{c2}%data%rAttr(km,:)

    call mct_rearr_rearrange(areasrc, areadst, Re_{ccc2}2{ccc1}, VECTOR=usevector, ALLTOALL=usealltoall)

    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_{c1}%data,"aream")
    dom_{c1}%data%rAttr(km,:) = areadst%rAttr(ka,:)

    call mct_aVect_clean(areasrc)
    call mct_aVect_clean(areadst)
!</list>

  end subroutine map_{ccc2}2{ccc1}_init_mct
!</list>

!=======================================================================
!<list> key="map1" number=1
  subroutine map_{ccc1}2{ccc2}_mct( cdata_{c1}, av_{c1}, cdata_{c2}, av_{c2}, fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata_{c1}
    type(mct_aVect), intent(in)  :: av_{c1}
    type(seq_cdata), intent(in)  :: cdata_{c2}
    type(mct_aVect), intent(out) :: av_{c2}
    character(len=*),intent(in), optional :: statelist
    character(len=*),intent(in), optional :: fluxlist
!</list>

!<list> key="map1" number=2
  subroutine map_rof2ocn_mct( cdata_{c1}, r2x_{c1}, cdata_{c2}, r2x_{c2})

    type(seq_cdata),intent(in) :: cdata_{c1}
    type(mct_aVect),intent(in) :: r2x_{c1}
    type(seq_cdata),intent(in) :: cdata_{c2}
    type(mct_aVect),intent(out):: r2x_{c2}

!</list>

!<list> key="local" number=1
    !
    ! Local
    !
    integer                :: lsize
    type(mct_aVect)        :: av_l_f     ! temporary flux attribute vector
    type(mct_aVect)        :: av_l_s     ! temporary state attribute vector
    character(*),parameter :: subName = '(map_{ccc1}2{ccc2}_mct) '
!</list>

!<list> key="local" number=2
    !
    ! Local variables
    !
    integer         :: lsize      ! temporary
    type(mct_aVect) :: av_i_f     ! temporary flux attribute vector
    type(mct_aVect) :: av_i_s     ! temporary state attribute vector
    type(mct_aVect) :: av_a_fl    ! temporary av for rearranging
    type(mct_aVect) :: av_i_fl    ! temporary av for rearranging
    character(*),parameter :: subName = '(map_{ccc1}2{ccc2}_mct) '
!</list>

!<list> key="local" number=3

    !
    ! Local Variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_ggrid), pointer :: dom_o
    type(mct_ggrid), pointer :: dom_a
    type(mct_gsMap), pointer :: gsmap_a
    integer                  :: mpicom
    integer                  :: lsize
    integer                  :: ku,kv,kux,kuy,kuz,n
    real(r8)                 :: ue,un,ur,ux,uy,uz,speed
    real(r8)                 :: urmaxl,urmax,uravgl,uravg,spavgl,spavg
    integer(in)              :: my_task,ierr,urcnt,urcntl
    type(mct_aVect)          :: av_o_f     ! temporary flux attribute vector
    type(mct_aVect)          :: av_o_s     ! temporary state attribute vector
    logical,save             :: write_info = .true.
    character(CL)            :: vect_map   ! vector mapping option
    character(*),parameter :: subName = '(map_atm2ocn_mct) '
!</list>

!<list> key="samegrid" number=1
    if (samegrid_map{c1}2{c2}) then
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
              call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
              call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, ALLTOALL=usealltoall)
       end if
    else
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
            call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_F{c1}2{c2}, rList=fluxlist, donorm=.false.)
          end if
          if (present(statelist)) then
            call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_S{c1}2{c2}, rList=statelist, donorm=.false.)
          end if
       else
          call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_F{c1}2{c2}, donorm=.false.)
       endif
    endif
!</list>

!<list> key="samegrid" number=2
    if (samegrid_map{c1}2{c2}) then
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
              call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
              call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, ALLTOALL=usealltoall)
       end if
    else
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
            call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_F{c1}2{c2}, rList=fluxlist, donorm=.false.)
          end if
          if (present(statelist)) then
            call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_S{c1}2{c2}, rList=statelist, donorm=.false.)
          end if
       else
          call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_F{c1}2{c2}, donorm=.false.)
       endif
    endif
!</list>

!<list> key="samegrid" number=3
    if (samegrid_map{c1}2{c2}) then
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
              call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
              call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, ALLTOALL=usealltoall)
       end if
    else
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
            call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_F{c1}2{c2}, rList=fluxlist, donorm=.false.)
          end if
          if (present(statelist)) then
            call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_S{c1}2{c2}, rList=statelist, donorm=.false.)
          end if
       else
          call seq_map_avNorm(av_{c1}, av_{c2}, sMatp_F{c1}2{c2}, donorm=.false.)
       endif
!</list>

!<list> key="samegrid" number=6
    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, ALLTOALL=usealltoall)
    end if
!</list>

!<list> key="samegrid" number=4
    if (samegrid_map{c1}2{c2}) then

          call mct_rearr_rearrange(r2x_{c1}, r2x_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
             ALLTOALL=usealltoall)

    else

!tcx    call mct_aVect_zero(r2x_o)
       call mct_sMat_avMult(r2x_{c1}, sMatp_F{c1}2{c2}, r2x_{c2}, VECTOR=usevector)

    endif

!</list>

!<list> key="samegrid" number=5
    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_{c1}, av_{c2}, Re_{ccc1}2{ccc2}, VECTOR=usevector, ALLTOALL=usealltoall)
    end if
!</list>

!<list> key="samegrid" number=3
      ! Correct {c1}->{c2} vector mapping near NP if appropriate

      if (trim(vect_map) == 'npfix') then
         if (write_info) write(logunit,*) subname,' :calling1 ',trim(vect_map)
         ku = mct_aVect_indexRA(av_{c1}, 'Sa_u', perrwith='quiet')
         kv = mct_aVect_indexRA(av_{c1}, 'Sa_v', perrwith='quiet')
         if (ku /= 0 .and. kv /= 0) then
            call seq_cdata_setptrs(cdata_{c2}, dom=dom_{c2})
            call seq_cdata_setptrs(cdata_{c1}, dom=dom_{c1}, gsmap=gsmap_{c1})
            call seq_cdata_setptrs(cdata_{c1}, mpicom=mpicom)
            call map_npfixNew4R(av_{c1}, av_{c2}, 'Sa_u', 'Sa_v', gsmap_{c1}, dom_{c1}, dom_{c2}, ni_{c1}, nj_{c1}, mpicom)
         end if
         write_info = .false.
      end if

      if (trim(vect_map(1:6)) == 'cart3d') then
         ku = mct_aVect_indexRA(av_{c1}, 'Sa_u', perrwith='quiet')
         kv = mct_aVect_indexRA(av_{c1}, 'Sa_v', perrwith='quiet')
         if (ku /= 0 .and. kv /= 0) then
            kux = mct_aVect_indexRA(av3d_{c1},'ux')
            kuy = mct_aVect_indexRA(av3d_{c1},'uy')
            kuz = mct_aVect_indexRA(av3d_{c1},'uz')
            lsize = mct_aVect_lsize(av_{c1})
            do n = 1,lsize
               ur = 0.0_r8
               ue = av_{c1}%rAttr(ku,n)
               un = av_{c1}%rAttr(kv,n)
               ux = clon_{c1}(n)*clat_{c1}(n)*ur - clon_{c1}(n)*slat_{c1}(n)*un - slon_{c1}(n)*ue
               uy = slon_{c1}(n)*clon_{c1}(n)*ur - slon_{c1}(n)*slat_{c1}(n)*un + clon_{c1}(n)*ue
               uz = slat_{c1}(n)          *ur + clat_{c1}(n)          *un
               av3d_{c1}%rAttr(kux,n) = ux
               av3d_{c1}%rAttr(kuy,n) = uy
               av3d_{c1}%rAttr(kuz,n) = uz
            enddo
            call mct_sMat_avMult(av3d_{c1}, sMatp_S{c1}2{c2}, av3d_{c2}, VECTOR=usevector)

            kux = mct_aVect_indexRA(av3d_{c2},'ux')
            kuy = mct_aVect_indexRA(av3d_{c2},'uy')
            kuz = mct_aVect_indexRA(av3d_{c2},'uz')
            lsize = mct_aVect_lsize(av_{c2})
            urmaxl = -1.0_r8
            uravgl = 0.0_r8
            urcntl = 0
            spavgl = 0.0_r8
            do n = 1,lsize
               ux = av3d_{c2}%rAttr(kux,n)
               uy = av3d_{c2}%rAttr(kuy,n)
               uz = av3d_{c2}%rAttr(kuz,n)
               ue = -slon_{c2}(n)          *ux + clon_{c2}(n)          *uy
               un = -clon_{c2}(n)*slat_{c2}(n)*ux - slon_{c2}(n)*slat_{c2}(n)*uy + clat_{c2}(n)*uz
               ur =  clon_{c2}(n)*clat_{c2}(n)*ux + slon_{c2}(n)*clat_{c2}(n)*uy - slat_{c2}(n)*uz
               speed = sqrt(ur*ur + ue*ue + un*un)
               if (trim(vect_map) == 'cart3d_diag' .or. trim(vect_map) == 'cart3d_uvw_diag') then
                  if (speed /= 0.0_r8) then
                     urmaxl = max(urmaxl,abs(ur))
                     uravgl = uravgl + abs(ur)
                     spavgl = spavgl + speed
                     urcntl = urcntl + 1
                  endif
               endif
               if (vect_map(1:10) == 'cart3d_uvw') then
                  if (write_info) write(logunit,*) subname,' :calling3 ',trim(vect_map)
                  !--- this adds ur to ue and un, while preserving u/v angle and
                  !total speed ---
                  if (un == 0.0_R8) then
                     !--- if ue is also 0.0 then just give speed to ue, this is
                     !arbitrary ---
                     av_{c2}%rAttr(ku,n) = sign(speed,ue)
                     av_{c2}%rAttr(kv,n) = 0.0_r8
                  else if (ue == 0.0_R8) then
                     av_{c2}%rAttr(ku,n) = 0.0_r8
                     av_{c2}%rAttr(kv,n) = sign(speed,un)
                  else
                     av_{c2}%rAttr(ku,n) = sign(speed/sqrt(1.0_r8 + ((un*un)/(ue*ue))),ue)
                     av_{c2}%rAttr(kv,n) = sign(speed/sqrt(1.0_r8 + ((ue*ue)/(un*un))),un)
                  endif
               else
                  if (write_info) write(logunit,*) subname,' :calling4 ',trim(vect_map)
                  !--- this ignores ur ---
                  av_{c2}%rAttr(ku,n) = ue
                  av_{c2}%rAttr(kv,n) = un
               endif
               write_info = .false.
            enddo
            if (trim(vect_map) == 'cart3d_diag' .or. trim(vect_map) == 'cart3d_uvw_diag') then
               call seq_cdata_setptrs(cdata_{c1}, mpicom=mpicom)
               call mpi_comm_rank(mpicom,my_task,ierr)
               call shr_mpi_max(urmaxl,urmax,mpicom,'urmax')
               call shr_mpi_sum(uravgl,uravg,mpicom,'uravg')
               call shr_mpi_sum(spavgl,spavg,mpicom,'spavg')
               call shr_mpi_sum(urcntl,urcnt,mpicom,'urcnt')
               if (my_task == 0 .and. urcnt > 0) then
                  uravg = uravg / urcnt
                  spavg = spavg / urcnt
                  write(logunit,*) trim(subname),' cart3d uravg,urmax,spavg = ',uravg,urmax,spavg
               endif
            endif
         endif  ! ku,kv
      endif  ! cart3d

    endif

!</list>

  end subroutine map_{ccc1}2{ccc2}_mct

!=======================================================================
!<list> key="public_interface" number=1
  subroutine map_{ccc2}2{ccc1}_mct( cdata_{c2}, av_{c2}, cdata_{c1}, av_{c1}, &
!{list} key="fraction" number=1
!                              fractions_{c2}, fractions_{c1}, &
!
                              fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata) ,intent(in)  :: cdata_{c2}
    type(mct_aVect) ,intent(in)  :: av_{c2}
    type(seq_cdata) ,intent(in)  :: cdata_{c1}
    type(mct_aVect) ,intent(out) :: av_{c1}
!<list> key="fraction" number=1
    type(mct_aVect) ,intent(in), optional :: fractions_{c2}
    type(mct_aVect) ,intent(in), optional :: fractions_{c1}
!</list>    

    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
!<list> key="local" number=1
    !
    ! Local
    !
    integer  :: i,j,kl,lsize,numats,ier
    type(mct_aVect)          :: av_a_f     ! temporary flux attribute vector
    type(mct_aVect)          :: av_a_s     ! temporary state attribute vector
!</list>

!<list> key="local" number=2
    !
    ! Local Variables
    !
    integer              :: n,ki,i,j       ! indices
    integer              :: numats,ier     ! number of attributes
    integer              :: lsize          ! size of attribute vector
    type(mct_aVect)      :: av_a_f, av_a_s ! temporary
!</list>

!<list> key="local" number=3
    !
    ! Local variables
    !
    type(mct_aVect)          :: av_a_f, av_a_s
    integer                  :: k{c1}, k{c2}, kSo_t
    integer                  :: numats,i,j,ier
    integer                  :: lsize
    real(R8),allocatable     :: recip(:)
    integer                  :: rcnt
    real(R8)                 :: rmax,rsum,rval
!</list>

    character(*),parameter :: subName = '(map_{ccc2}2{ccc1}_mct) '
!<list> key="samegrid" number=1
    if (samegrid_map{c1}2{c2}) then
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, ALLTOALL=usealltoall)
       end if
    else
        if (present(fractions_{c2})) then
          if (present(fluxlist) .or. present(statelist)) then
             if (present(fluxlist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, fractions_{c2}, '{c2}frin', fractions_{c1}, '{c2}frin', rList=fluxlist)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_S{c2}2{c1}, fractions_{c2}, '{c2}frin', fractions_{c1}, '{c2}frin', rList=statelist)
             end if
          else
             call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, fractions_{c2}, '{c2}frin', fractions_{c1}, '{c2}frin')
          endif
       else
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, rList=fluxlist, donorm=.false.)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_S{c2}2{c1}, rList=statelist, donorm=.false.)
             end if
          else
             ! --- default is flux mapping
             call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, donorm=.false.)
          endif
       endif
    endif
!</list>

!<list> key="samegrid" number=2
    if (samegrid_map{c1}2{c2}) then
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, ALLTOALL=usealltoall)
       end if
    else
        if (present(fractions_{c2})) then
          if (present(fluxlist) .or. present(statelist)) then
             if (present(fluxlist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, fractions_{c2}, '{c2}fr{c1}c', fractions_{c1}, '{c2}fr{c1}c', rList=fluxlist)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_S{c2}2{c1}, fractions_{c2}, '{c2}fr{c1}c', fractions_{c1}, '{c2}fr{c1}c', rList=statelist)
             end if
          else
             call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, fractions_{c2}, '{c2}fr{c1}c', fractions_{c1}, '{c2}fr{c1}c')
          endif
       else
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, rList=fluxlist, donorm=.false.)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_S{c2}2{c1}, rList=statelist, donorm=.false.)
             end if
          else
             ! --- default is flux mapping
             call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, donorm=.false.)
          endif
       endif
    endif
!</list>

!<list> key="samegrid" number=3
    if (samegrid_map{c1}2{c2}) then
       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, ALLTOALL=usealltoall)
       end if
    else
        if (present(fractions_{c2})) then
          if (present(fluxlist) .or. present(statelist)) then
             if (present(fluxlist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, fractions_{c2}, '{c2}fr{c1}c', fractions_{c1}, '{c2}fr{c1}c', rList=fluxlist)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_S{c2}2{c1}, fractions_{c2}, '{c2}fr{c1}c', fractions_{c1}, '{c2}fr{c1}c', rList=statelist)
             end if
          else
             call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, fractions_{c2}, '{c2}fr{c1}c', fractions_{c1}, '{c2}fr{c1}c')
          endif
       else
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, rList=fluxlist, donorm=.false.)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_S{c2}2{c1}, rList=statelist, donorm=.false.)
             end if
          else
             ! --- default is flux mapping
             call seq_map_avNorm(av_{c2}, av_{c1}, sMatp_F{c2}2{c1}, donorm=.false.)
          endif
       endif
    endif
!</list>

!<list> key="samegrid" number=6
    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, ALLTOALL=usealltoall)
    end if
!</list>

!<list> key="samegrid" number=5
    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_{c2}, av_{c1}, Re_{ccc2}2{ccc1}, VECTOR=usevector, ALLTOALL=usealltoall)
    end if

!</list>

  end subroutine map_{ccc2}2{ccc1}_mct
!</list>

!<list> key="npfixnew4r" number=1
 subroutine map_npFixNew4R(buni,buno,fld1,fld2,gsmapi,domi,domo,ni_i,nj_i,mpicom)

   !===============================================================================
   !    Correct the north pole mapping of velocity fields from the atm to ocn
   !    grids.  This assumes the input grid is a regular lat/lon with the north
   !    pole surrounded by the last latitude line in the input array.  The
   !    longitudes in the last latitude must be ordered and equally spaced.
   !
   !    4R is a low memory version of 3R.
   !    This version (New4R) is the same as 3R except it uses a lot less memory
   !    and is a bit faster.  Like 3R, it saves data between calls and so
   !    assumes the input grid remains constant for all calls.  This is bfb
   !    with 3R on bluevista as of 2/15/07.
   !
   !    !REVISION HISTORY:
   !    2007-Feb-12 - T. Craig -- modified New3R to reduce memory
   !    2007-Apr-27 - M. Vertenstein - implemented in sequential system
   !===============================================================================

#include <mpif.h>
   !
   ! Arguments
   !
   type(mct_Avect),intent(in)   :: buni    ! input  attribute vec
   type(mct_Avect),intent(out)  :: buno    ! output attribute vec
   character(*)   ,intent(in)   :: fld1    ! name of first input field
   character(*)   ,intent(in)   :: fld2    ! name of second input field
   integer        ,intent(in)   :: ni_i    ! number of longitudes in input grid
   integer        ,intent(in)   :: nj_i    ! number of latitudes in input grid
   type(mct_gsMap),pointer      :: gsmapi  ! input gsmap
   type(mct_gGrid),pointer      :: domi    ! input domain
   type(mct_gGrid),pointer      :: domo    ! output domain
   integer        ,intent(in)   :: mpicom  ! mpi communicator group
   !
   ! Local Variables
   !
   integer(IN)  :: n,m                           ! generic indices
   integer(IN)  :: n1,n2,n3                      ! generic indices
   integer(IN)  :: kui,kvi                       ! field indices
   integer(IN)  :: kuo,kvo                       ! field indices
   integer(IN)  :: kin                           ! index index
   integer(IN)  :: nmin,nmax                     ! indices of highest latitude in input
   integer(IN)  :: npts                          ! local number of points in an aV
   integer(IN)  :: num                           ! number of points at highest latitude
   integer(IN)  :: kloni                         ! longitude index on input domain
   integer(IN)  :: klati                         ! latitude index on input domain
   integer(IN)  :: klono                         ! longitude index on output domain
   integer(IN)  :: klato                         ! latitude index on output domain
   integer(IN)  :: index                         ! index value
   real(R8)     :: rindex                        ! index value
   real(R8)     :: latmax                        ! value of highest latitude
   real(R8)     :: olon,olat                     ! output bundle lon/lat
   real(R8)     :: ilon,ilat                     ! input bundle lon/lat
   real(R8)     :: npu,npv                       ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2                 ! angles for trig functions 
   real(R8),allocatable,save :: ilon1(:)         ! lon of input grid at highest latitude
   real(R8),allocatable,save :: ilat1(:)         ! lat of input grid at highest latitude
   real(R8),allocatable,save :: ilon2(:)         ! lon of input grid at highest latitude
   real(R8),allocatable,save :: ilat2(:)         ! lat of input grid at highest latitude
   real(R8),allocatable      :: rarray(:)        ! temporary array
   real(R8),allocatable      :: rarray2(:,:)     ! temporary array
   real(R8)     :: w1,w2,w3,w4                   ! weights
   real(R8)     :: f1,f2,f3,f4                   ! function values
   real(R8)     :: alpha,beta                    ! used to generate weights
   real(R8)     :: rtmp                          ! real temporary
   real(R8),allocatable :: lData(:,:)            ! last lat local input bundle data
                                                 ! also compressed global data
   real(R8)   ,allocatable,save :: alphafound(:) ! list of found alphas
   real(R8)   ,allocatable,save :: betafound(:)  ! list of found betas
   integer(IN),allocatable,save :: mfound(:)     ! list of found ms
   integer(IN),allocatable,save :: nfound(:)     ! list of found ns
   real(R8)   ,allocatable      :: rfound(:)     ! temporary for copy
   integer(IN),allocatable      :: ifound(:)     ! temporary for copy
   integer(IN),save             :: cntfound      ! number found
   integer(IN),save             :: cntf_tot      ! cntfound total
   integer(IN)  :: cnt                           ! loop counter
   logical      :: found                         ! search for new interpolation
   integer(IN)  :: rcode                         ! error code
   integer(IN)  :: np1                           ! n+1 or tmp
   real(R8)     :: ilon1x                        ! tmp
   integer(IN),pointer,save  :: gindex(:)        ! global index
   logical,save :: first_call = .true.           ! flags 1st invocation of routine
   integer(IN),save   :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9
   integer(IN) :: mype

   !--- formats ---
   character(*),parameter :: subName = '(map_npFixNew4R) '
   character(*),parameter :: F00 = "('(map_npFixNew4R) ',8a)"
   character(*),parameter :: F01 = "('(map_npFixNew4R) ',a,i12)"
   !-------------------------------------------------------------------------------

   call MPI_COMM_RANK(mpicom,mype,rcode)

   kui   = mct_aVect_indexRA(buni,fld1,perrWith=subName)
   kvi   = mct_aVect_indexRA(buni,fld2,perrWith=subName)
   kuo   = mct_aVect_indexRA(buno,fld1,perrWith=subName)
   kvo   = mct_aVect_indexRA(buno,fld2,perrWith=subName)

!  tcx 3/19/08, don't use GlobGridNum, it's not set properly in models
!   kin   = mct_aVect_indexIA(domi%data,"GlobGridNum",perrWith=subName)
   klati = mct_aVect_indexRA(domi%data,"lat"        ,perrWith=subName)
   kloni = mct_aVect_indexRA(domi%data,"lon"        ,perrWith=subName)
   klato = mct_aVect_indexRA(domo%data,"lat"        ,perrWith=subName)
   klono = mct_aVect_indexRA(domo%data,"lon"        ,perrWith=subName)

   ! ni_i, nj_i should be read in from mapping file

   nmin = (ni_i)*(nj_i-1) + 1
   nmax = ni_i*nj_i
   num  = ni_i

   !---------------------------------------------------------------------------
   ! Initialization only
   !---------------------------------------------------------------------------

   if (first_call) then

     if (loglevel > 0) write(logunit,F00) " compute bilinear weights & indicies for NP region."

     allocate(ilon1(num))
     allocate(ilon2(num))
     allocate(ilat1(num))
     allocate(ilat2(num))

     ilon1 = 0._r8
     ilon2 = 0._r8
     ilat1 = 0._r8
     ilat2 = 0._r8

     npts = mct_aVect_lSize(domi%data)
     call mct_gsMap_orderedPoints(gsMapi, mype, gindex)
     do m=1,npts
       if (gindex(m).ge.nmin) then               ! are on highest latitude
         n = gindex(m) - nmin + 1                ! n is 1->ni_i lon index on highest latitude
         rtmp = domi%data%rAttr(kloni,m)      ! rtmp is longitude value on highest latitude
         ilon1(n) = mod(rtmp+360._R8,360._R8) ! ilon1(n) is longitude val mapped from 0->360
         ilat1(n) = domi%data%rAttr(klati,m)  ! ilat1(n) values should all be the same (i.e. highest lat)
       endif
     enddo

     !--- all gather local data, MPI_SUM is low memory and simple
     !--- but is a performance penalty compared to gatherv and copy
     !--- or a fancy send/recv

     allocate(rarray(num))
     rarray = ilat1
     call MPI_ALLREDUCE(rarray,ilat1,num,MPI_REAL8,MPI_SUM,mpicom,rcode)
     if (rcode.ne.0) then
       write(logunit,*) trim(subName),' ilat1 rcode error ',rcode
       call shr_sys_abort()
     endif
     rarray = ilon1
     call MPI_ALLREDUCE(rarray,ilon1,num,MPI_REAL8,MPI_SUM,mpicom,rcode)
     if (rcode.ne.0) then
       write(logunit,*) trim(subName),' ilon1 rcode error ',rcode
       call shr_sys_abort()
     endif

     do n = 1,num
       np1 = mod(n,num)+1
       ilat2(n) = ilat1(np1)
       ilon2(n) = ilon1(np1)
       if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360._R8
     enddo

     do n = 1,num
       if (ilat1(n) /= ilat2(n)) then
          write(logunit,*) trim(subname),' ERROR: ilat1 ne ilat2 ',n,ilat1(n),ilat2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilat1 ne ilat2')
       endif
       if (ilon2(n) < ilon1(n)) then
          write(logunit,*) trim(subname),' ERROR: ilon2 lt ilon1 ',n,ilon1(n),ilon2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilon2 ilon1 error')
       endif
       ! tcraig check that lon diffs are reasonable 4x average dlon seems like
       ! reasonable limit
       if (ilon2(n) - ilon1(n) > (360.0_R8/(num*1.0_R8))*4.0) then
          write(logunit,*) trim(subname),' ERROR: ilon1,ilon2 ',n,ilon1(n),ilon2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilon2 ilon1 size diff')
       endif
     enddo

     latmax = maxval(ilat1)

     !--- compute weights and save them ---

     npts = mct_aVect_lSize(buno)
     allocate(mfound(npts),nfound(npts),alphafound(npts),betafound(npts))
     cntfound = 0
     do m = 1,npts
       olat = domo%data%rAttr(klato,m)
       if (olat >= latmax) then
         rtmp = domo%data%rAttr(klono,m)
         olon = mod(rtmp,360._R8)
         n = 1
         found = .false.
         do while (n <= num .and. .not.found )
           if (    olon         >= ilon1(n) .and. olon < ilon2(n) .or.   &
                   olon+360._R8 >= ilon1(n) .and. olon < ilon2(n)) then


!tcx ilat2==ilat1 so don't average
!--->        ilat = (ilat1(n) + ilat2(n)) * 0.5_R8
             ilat = ilat1(n)
             if (ilon2(n) == ilon1(n)) then
               alpha = 0.5_R8
             else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
               alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
             else if (olon+360._R8>= ilon1(n) .and. olon < ilon2(n)) then
               alpha = (olon+360._R8 - ilon1(n)) / (ilon2(n) - ilon1(n))
             else
               write(logunit,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
             endif
             if (ilat >= 90._R8) then
               beta  = 1.0_R8
             else
               beta  = (olat - ilat) / (90._R8 - ilat)
             endif

             cntfound = cntfound + 1
             mfound(cntfound) = m
             nfound(cntfound) = n
             alphafound(cntfound) = alpha
             betafound(cntfound) = beta
             found = .true.

           endif
           n = n + 1     ! normal increment
         enddo
         if ( .not.found ) then
           write(logunit,*) subName,' ERROR: found = false ',found,m,olon,olat
         endif
       endif
     end do

     allocate(ifound(npts))
     ifound(1:cntfound) = mfound(1:cntfound)
     deallocate(mfound)
     if (cntfound > 0) then
        allocate(mfound(cntfound))
        mfound(1:cntfound) = ifound(1:cntfound)
     endif

     ifound(1:cntfound) = nfound(1:cntfound)
     deallocate(nfound)
     if (cntfound > 0) then
        allocate(nfound(cntfound))
        nfound(1:cntfound) = ifound(1:cntfound)
     endif
     deallocate(ifound)

     allocate(rfound(npts))
     rfound(1:cntfound) = alphafound(1:cntfound)
     deallocate(alphafound)
     if (cntfound > 0) then
        allocate(alphafound(cntfound))
        alphafound(1:cntfound) = rfound(1:cntfound)
     endif

     rfound(1:cntfound) = betafound(1:cntfound)
     deallocate(betafound)
     if (cntfound > 0) then
        allocate(betafound(cntfound))
        betafound(1:cntfound) = rfound(1:cntfound)
     endif
     deallocate(rfound)

     call MPI_ALLREDUCE(cntfound,cntf_tot,1,MPI_INTEGER,MPI_SUM,mpicom,rcode)
     if (mype == 0) then
        write(logunit,F01) ' total npfix points found = ',cntf_tot
     endif

     first_call = .false.

   endif

   !---------------------------------------------------------------------------
   ! Return if there is nothing to do; must be no points on any pes
   ! If there are any npfix points, all pes must continue to the np u,v calc
   !---------------------------------------------------------------------------

   if (cntf_tot < 1) then
      return
   endif

   !---------------------------------------------------------------------------
   ! Non-initialization, run-time fix
   !---------------------------------------------------------------------------

   !--- barrier not required but interesting for timing. ---
   !  call shr_mpi_barrier(mpicom,subName//" barrier")

   !--- extract index, u, v from buni ---

   allocate(lData(3,num))
   lData = 0._R8
   npts = mct_aVect_lSize(buni)
   do n=1,npts
     if (gindex(n).ge.nmin) then
       m = gindex(n) - nmin + 1
       lData(1,m) = gindex(n)
       lData(2,m) = buni%rAttr(kui,n)
       lData(3,m) = buni%rAttr(kvi,n)
     endif
   enddo

   !--- all gather local data, MPI_SUM is low memory and simple
   !--- but is a performance penalty compared to gatherv and copy
   !--- KLUDGE - this should be looked at when it becomes a performance/memory
   !--- penalty

   allocate(rarray2(3,num))
   rarray2=lData
   call MPI_ALLREDUCE(rarray2,lData,3*num,MPI_REAL8,MPI_SUM,mpicom,rcode)
   deallocate(rarray2)

   if (rcode.ne.0) then
     write(logunit,*) trim(subName),' rcode error ',rcode
     call shr_sys_abort()
   endif

   do n2=1,num
     if (lData(1,n2).lt.0.1_R8) then
       write(logunit,*) trim(subName),' error allreduce ',n2
     endif
   enddo

   !--- compute npu, npv (pole data) and initialize ilon,ilat arrays ---

   npu = 0._R8
   npv = 0._R8
   do n = 1,num
     theta1 = ilon1(n)*shr_const_deg2rad
     npu = npu + cos(theta1)*lData(2,n) &
               - sin(theta1)*lData(3,n)
     npv = npv + sin(theta1)*lData(2,n) &
               + cos(theta1)*lData(3,n)
   enddo
   npu = npu / real(num,R8)
   npv = npv / real(num,R8)

   !--- compute updated pole vectors ---

!DIR$ CONCURRENT
   do cnt = 1,cntfound
      m     = mfound(cnt)
      n     = nfound(cnt)
      np1   = mod(n,num)+1
      alpha = alphafound(cnt)
      beta  = betafound(cnt)

      w1 = (1.0_R8-alpha)*(1.0_R8-beta)
      w2 = (    alpha)*(1.0_R8-beta)
      w3 = (    alpha)*(    beta)
      w4 = (1.0_R8-alpha)*(    beta)

      theta1 = ilon1(n)*shr_const_deg2rad
      theta2 = ilon2(n)*shr_const_deg2rad

      f1 = lData(2,n)
      f2 = lData(2,np1)
      f3 =  cos(theta1)*npu + sin(theta1)*npv
      f4 =  cos(theta2)*npu + sin(theta2)*npv
      rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
      buno%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

      f1 = lData(3,n)
      f2 = lData(3,np1)
      f3 = -sin(theta1)*npu + cos(theta1)*npv
      f4 = -sin(theta2)*npu + cos(theta2)*npv
      rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
      buno%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
   enddo

   deallocate(lData)
   first_call = .false.

end subroutine map_npFixNew4R
!</list>

end module map_{ccc1}{ccc2}_mct
