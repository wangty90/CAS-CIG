!===============================================================================
! SVN $Id: seq_avdata_mod.F90 18516 2009-09-25 22:54:10Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/branch_tags/t3148b_tags/t3148b02_drvseq3_1_48/driver/seq_avdata_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_avdata_mod -- provides use access to public cpl7 aVect, domain,
!     and fraction data.
!
! !DESCRIPTION:
!
!    provides use access to public cpl7 aVect, domain, and fraction data.
!
! !REMARKS:
!
!    use access to public cpl7 aVect, domain, and fraction info is to avoid 
!    excessively long routine arg lists, eg. for history & restart modules.
!    Note: while cpl7's non-main program ("driver") routines CAN access this
!    data by use'ing this module, they SHOULD access it via agrument lists 
!    if it is reasonable to do so.  Do the right thing.
!
! !REVISION HISTORY:
!     2009-Sep-25 - B. Kauffman - initial version
!
! !INTERFACE: ------------------------------------------------------------------

module seq_avdata_mod

! !USES:

   use shr_kind_mod  ,only: IN => SHR_KIND_IN
   use mct_mod           ! mct_ wrappers for mct lib

   use seq_cdata_mod     ! "cdata" type & methods (domain + decomp + infodata in one datatype)
   use seq_infodata_mod  ! "infodata" gathers various control flags into one datatype

   implicit none

   public  ! default is public

! !PUBLIC DATA MEMBERS:

   !----------------------------------------------------------------------------
   ! Infodata: inter-model control flags, domain info
   !----------------------------------------------------------------------------

   type (seq_infodata_type) :: infodata ! single instance for cpl and all comps

   !----------------------------------------------------------------------------
   ! cdata types: contains pointers to domain info + component ID + infobuffer
   !----------------------------------------------------------------------------

!{list} key="listc"
!  type (seq_cdata) :: cdata_{c}{c}
!wangty modify
#ifdef wrf
   type (seq_cdata) :: cdata_ww ! juanxiong he for wrf/ccsm4 coupling
   type (seq_cdata) :: cdata_cc ! juanxiong he for wrf/cam coupling, cam end
   type (seq_cdata) :: cdata_mm ! juanxiong he for wrf/cam coupling, wrf end
   type (seq_cdata) :: cdata_gege ! juanxiong he for geatm/cam coupling, geatm end
   type (seq_cdata) :: cdata_caca ! juanxiong he for geatm/cam coupling, cam end
#endif

!{list} key="listc"
!  type (seq_cdata) :: cdata_{c}x
!wangty modify
#ifdef wrf
   type (seq_cdata) :: cdata_wx ! juanxiong he for wrf/ccsm4 coupling
   type (seq_cdata) :: cdata_cx ! juanxiong he for wrf/cam coupling, cam end
   type (seq_cdata) :: cdata_mx ! juanxiong he for wrf/cam coupling, wrf end
   type (seq_cdata) :: cdata_cax ! juanxiong he for geatm/cam coupling, cam end
   type (seq_cdata) :: cdata_gex ! juanxiong he for geatm/cam coupling, geatm end
#endif
   !----------------------------------------------------------------------------
   ! domain info: coords, fractions, decomps, area correction factors
   !----------------------------------------------------------------------------

   !--- domain coords, area, mask  (MCT General Grids) --

!{list} key="listc"
!  type(mct_gGrid) :: dom_{c}{c}
!wangty modify
#ifdef wrf
   type(mct_gGrid) :: dom_ww      ! wrf domain,juanxiong he for wrf/ccsm4 coupling
   type(mct_gGrid) :: dom_cc      ! cam domain,juanxiong he for wrf/cam coupling, cam end
   type(mct_gGrid) :: dom_mm      ! wrf domain,juanxiong he for wrf/cam coupling, wrf grid wrf end
   type(mct_gGrid) :: dom_caca      ! cam domain,juanxiong he for geatm/cam coupling, cam end
   type(mct_gGrid) :: dom_gege      ! wrf domain,juanxiong he for geatm/cam coupling, geatm grid geatm end
#endif

!{list} key="listc"
!  type(mct_gGrid) :: dom_{c}x
!wangty modify
#ifdef wrf
   type(mct_gGrid) :: dom_wx      ! wrf domain for wrf/ccsm4 coupling, juanxiong he
   type(mct_gGrid) :: dom_cx      ! cam domain for wrf/cam coupling, cam end, juanxiong he
   type(mct_gGrid) :: dom_mx      ! wrf domain for wrf/cam coupling, wrf end, juanxiong he
   type(mct_gGrid) :: dom_cax      ! cam domain for geatm/cam coupling, cam end, juanxiong he
   type(mct_gGrid) :: dom_gex      ! wrf domain for geatm/cam coupling, geatm end, juanxiong he
#endif

   !--- domain fractions (only defined on cpl pes) ---

!{list} key="fractionc"
!  type(mct_aVect) :: fractions_{c}x   
!wangty modify
#ifdef wrf
   type(mct_aVect) :: fractions_wx   ! Fractions on wrf grid, juanxiong he
   type(mct_aVect) :: fractions_gex   ! Fractions on geatm grid, juanxiong he
#endif

   !----------------------------------------------------------------------------
   ! State/flux field bundles (MCT attribute vectors)
   !----------------------------------------------------------------------------

!{list} key="cimport"
!  type(mct_aVect) :: x2{c}_{c}{c}

!{list} key="listc"
!  type(mct_aVect) :: {c}2x_{c}{c}

!{list} key="cimport"
!  type(mct_aVect) :: x2{c}_{c}x

!{list} key="listc"
!  type(mct_aVect) :: {c}2x_{c}x

!{list} key="mergec"
!  type(mct_aVect) :: {c}2x_{co}x

   type(mct_aVect) :: xao_ox    ! Atm-ocn fluxes, ocn grid, cpl pes - defined in flux_ao gc 
   type(mct_aVect) :: xao_ax    ! Atm-ocn fluxes, atm grid, cpl pes - defined in flux_ao gc 
   type(mct_accum) :: r2xacc_rx ! Rof export, rof grid, cpl pes - defined in driver
   type(mct_accum) :: x2oacc_ox ! Ocn import, ocn grid, cpl pes - defined in driver

   integer(IN)     :: r2xacc_rx_cnt ! r2xacc_rx: number of time samples accumulated
   integer(IN)     :: x2oacc_ox_cnt ! x2oacc_ox: number of time samples accumulated
!wangty modify
#ifdef wrf
  !--------------------------------------------------------------------------------
  ! juanxiong he for geatm/ccsm4 coupling
  !--------------------------------------------------------------------------------
   type(mct_aVect) :: chem2x_lx    ! geatm export, lnd grid, cpl pes - defined in mrg_x2l
   type(mct_aVect) :: chem2x_ix    ! geatm export, ice grid, cpl pes - defined in mrg_x2i
   type(mct_aVect) :: chem2x_ox    ! geatm export, ocn grid, cpl pes - defined in mrg_x2o
   type(mct_aVect) :: xchemo_ox    ! geatm-ocn fluxes, ocn grid, cpl pes - defined in flux_wo gc
   type(mct_aVect) :: xchemo_gex    ! geatm-ocn fluxes, wrf grid, cpl pes - defined in flux_wo gc

   type(mct_aVect) :: l2x_chemx    ! lnd export, geatm grid, cpl pes - defined in mrg_x2ge
   type(mct_aVect) :: l2x_chemchem    ! lnd export, geatm grid, cpl pes - defined in mrg_x2ge
   type(mct_aVect) :: x2l_chemx    ! lnd import, geatm grid, defunct 
   type(mct_aVect) :: x2l_chemchem   ! lnd import, geatm grid, defunct 
   type(mct_aVect) :: i2x_chemx    ! Ice export, geatm grid, cpl pes - defined in mrg_x2ge
   type(mct_aVect) :: i2x_chemchem    ! Ice export, geatm grid, cpl pes - defined in mrg_x2ge
   type(mct_aVect) :: x2i_chemx    ! Ice import, geatm grid, defunct
   type(mct_aVect) :: x2i_chemchem   ! Ice import, geatm grid, defunct

   type(mct_aVect) :: o2x_chemx    ! ocn export, geatm grid, cpl pes - defined in mrg_x2ge
   type(mct_aVect) :: o2x_chemchem    ! ocn export, geatm grid, cpl pes - defined in mrg_x2ge
   type(mct_accum) :: x2oacc_chemx ! Ocn import, ocn grid, cpl pes - defined in driver
   type(mct_aVect) :: x2o_chemx    ! onc import, geatm grid, defunct
   type(mct_aVect) :: x2o_chemchem  ! onc import, geatm grid, defunct

   ! cam/geatm coupling
   ! cam to geatm
   type(mct_aVect) :: ca2x_caca1    ! Cam export, cam grid, cam end
   type(mct_aVect) :: ca2x_caca2    ! Cam export, cam grid, cam end
   type(mct_aVect) :: ca2x_cax1    ! Cam export, cam grid, rearrange, cpl end
   type(mct_aVect) :: ca2x_cax2    ! Cam export, cam grid, rearrange, cpl end
   type(mct_aVect) :: ca2x_chemx1    ! Cam export, wrf grid, defined in map_wrfcam, cpl end
   type(mct_aVect) :: ca2x_chemx2    ! Cam export, wrf grid, defined in map_wrfcam, cpl end
   type(mct_aVect) :: x2chem_chemx1    ! wrf import, wrf grid, merge, cpl end
   type(mct_aVect) :: x2chem_chemx2    ! wrf import, wrf grid, merge, cpl end
   type(mct_aVect) :: x2chem_chemchem1    ! wrf import, wrf grid, rearrange, wrf end
   type(mct_aVect) :: x2chem_chemchem2    ! wrf import, wrf grid, rearrange, wrf end
   ! geatm to cam
   type(mct_aVect) :: x2ca_caca1    ! Cam import, cam grid, rearrange, cam end
   type(mct_aVect) :: x2ca_caca2    ! Cam import, cam grid, rearrange, cam end
   type(mct_aVect) :: x2ca_cax1    ! Cam import, cam grid, merge, cpl end
   type(mct_aVect) :: x2ca_cax2    ! Cam import, cam grid, merge, cpl end
   type(mct_aVect) :: chem2x_cax    ! Cam import, cam grid, defined in map_camwrf, cpl end
   type(mct_aVect) :: chem2x_chemx    ! wrf export, wrf grid, rearrange, cpl end
   type(mct_aVect) :: chem2x_chemchem    ! wrf export, wrf grid, wrf end

  !--------------------------------------------------------------------------------
  ! juanxiong he for geatm/ccsm4 coupling
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  ! juanxiong he for wrf/cam/ccsm4 coupling
  !--------------------------------------------------------------------------------
   type(mct_aVect) :: x2w_ww    ! Wrf import, wrf grid, atm pes - defined in wrf gc
   type(mct_aVect) :: w2x_ww    ! Wrf export, wrf grid, atm pes - defined in wrf gc
   type(mct_aVect) :: x2w_wx    ! Wrf import, wrf grid, cpl pes - defined in map_wrfwrf
   type(mct_aVect) :: w2x_wx    ! Wrf export, wrf grid, cpl pes - defined in map_wrfwrf
   type(mct_aVect) :: w2x_lx    ! Wrf export, lnd grid, cpl pes - defined in mrg_x2l
   type(mct_aVect) :: w2x_ix    ! Wrf export, ice grid, cpl pes - defined in mrg_x2i
   type(mct_aVect) :: w2x_ox    ! Wrf export, ocn grid, cpl pes - defined in mrg_x2o
   type(mct_aVect) :: xwo_ox    ! Wrf-ocn fluxes, ocn grid, cpl pes - defined in flux_wo gc
   type(mct_aVect) :: xwo_wx    ! Wrf-ocn fluxes, wrf grid, cpl pes - defined in flux_wo gc

   type(mct_aVect) :: l2x_wx    ! lnd export, wrf grid, cpl pes - defined in mrg_x2w
   type(mct_aVect) :: l2x_ww    ! lnd export, wrf grid, cpl pes - defined in mrg_x2w
   type(mct_aVect) :: x2l_wx    ! lnd import, wrf grid, defunct
   type(mct_aVect) :: x2l_ww    ! lnd import, wrf grid, defunct

   type(mct_aVect) :: i2x_wx    ! Ice export, wrf grid, cpl pes - defined in mrg_x2w
   type(mct_aVect) :: i2x_ww    ! Ice export, wrf grid, cpl pes - defined in mrg_x2w
   type(mct_aVect) :: x2i_wx    ! Ice import, wrf grid, defunct
   type(mct_aVect) :: x2i_ww    ! Ice import, wrf grid, defunct

   type(mct_aVect) :: o2x_wx    ! ocn export, wrf grid, cpl pes - defined in mrg_x2w
   type(mct_aVect) :: o2x_ww    ! ocn export, wrf grid, cpl pes - defined in mrg_x2w
   type(mct_accum) :: x2oacc_wx ! Ocn import, ocn grid, cpl pes - defined in driver
   type(mct_aVect) :: x2o_wx    ! onc import, wrf grid, defunct
   type(mct_aVect) :: x2o_ww    ! onc import, wrf grid, defunct

   ! cam/wrf coupling
   ! cam to wrf
   type(mct_aVect) :: c2x_cc1    ! Cam export, cam grid, cam end
   type(mct_aVect) :: c2x_cc2    ! Cam export, cam grid, cam end
   type(mct_aVect) :: c2x_cx1    ! Cam export, cam grid, rearrange, cpl end
   type(mct_aVect) :: c2x_cx2    ! Cam export, cam grid, rearrange, cpl end
   type(mct_aVect) :: c2x_mx1    ! Cam export, wrf grid, defined in map_wrfcam, cpl end
   type(mct_aVect) :: c2x_mx2    ! Cam export, wrf grid, defined in map_wrfcam, cpl end
   type(mct_aVect) :: x2m_mx1    ! wrf import, wrf grid, merge, cpl end
   type(mct_aVect) :: x2m_mx2    ! wrf import, wrf grid, merge, cpl end
   type(mct_aVect) :: x2m_mm1    ! wrf import, wrf grid, rearrange, wrf end
   type(mct_aVect) :: x2m_mm2    ! wrf import, wrf grid, rearrange, wrf end
   ! wrf to cam
   type(mct_aVect) :: x2c_cc1    ! Cam import, cam grid, rearrange, cam end
   type(mct_aVect) :: x2c_cc2    ! Cam import, cam grid, rearrange, cam end
   type(mct_aVect) :: x2c_cx1    ! Cam import, cam grid, merge, cpl end
   type(mct_aVect) :: x2c_cx2    ! Cam import, cam grid, merge, cpl end
   type(mct_aVect) :: m2x_cx    ! Cam import, cam grid, defined in map_camwrf, cpl end
   type(mct_aVect) :: m2x_mx    ! wrf export, wrf grid, rearrange, cpl end
   type(mct_aVect) :: m2x_mm    ! wrf export, wrf grid, wrf end

  !--------------------------------------------------------------------------------
  ! juanxiong he for wrf/cam/ccsm4 coupling
  !--------------------------------------------------------------------------------
#endif

! !PUBLIC MEMBER FUNCTIONS

   ! no public routines

! !PUBLIC TYPES:

   ! no public types

end module seq_avdata_mod

!===============================================================================
