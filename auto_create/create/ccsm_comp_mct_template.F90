
module ccsm_comp_mod

#define NEW_BUDGET
#define DEBUG

!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CCSM4/cpl7. Can have different
!          land, sea-ice, and ocean models plugged in at compile-time.
!          These models can be either: stub, dead, data, or active
!          components or some combination of the above.
!
!               stub -------- Do nothing.
!               dead -------- Send analytic data back.
!               data -------- Send data back interpolated from input
!               files.
!               prognostic -- Prognostically simulate the given
!               component.
!
! Method: Call appropriate initialization, run (time-stepping), and
!         finalization routines.
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! share code & libs
   !----------------------------------------------------------------------------
   use shr_kind_mod,      only: r8 => SHR_KIND_R8
   use shr_kind_mod,      only: cs => SHR_KIND_CS
   use shr_kind_mod,      only: cl => SHR_KIND_CL
   use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
   use shr_const_mod,     only: shr_const_cday
   use shr_file_mod,      only: shr_file_setLogLevel, shr_file_setLogUnit
   use shr_file_mod,      only: shr_file_setIO, shr_file_getUnit
   use shr_scam_mod,      only: shr_scam_checkSurface
   use shr_map_mod,       only: shr_map_setDopole
   use shr_mpi_mod,       only: shr_mpi_min, shr_mpi_max
   use shr_mem_mod,       only: shr_mem_init, shr_mem_getusage
   use shr_cal_mod,       only: shr_cal_date2ymd
   use shr_orb_mod,       only: shr_orb_params
   use mct_mod            ! mct_ wrappers for mct lib
   use perf_mod
   use ESMF_Mod

   !----------------------------------------------------------------------------
   ! component model interfaces (init, run, final methods)
   !----------------------------------------------------------------------------
!wangty modify
#ifdef wrf
   use wrf_comp_mct, only: wrf_init_mct, wrf_run_mct, wrf_final_mct  ! juanxiong he
   use geatm_comp_mct, only: geatm_init_mct, geatm_run_mct, geatm_final_mct  ! juanxiong he
#endif

!{list} key="components" number=1
!   use {ccc}_comp_mct, only: {ccc}_init_mct, {ccc}_run_mct, {ccc}_final_mct

#ifdef ESMF_INTERFACE
   use esmfshr_attribute_mod
!{list} key="components" number=1
!   use {ccc}_comp_mct, only: {ccc}_register

#endif

   !----------------------------------------------------------------------------
   ! cpl7 modules
   !----------------------------------------------------------------------------

   !--- modules with public read/write data ---
   use seq_avdata_mod    ! drv aVects & associated domain, fraction, cdata
   use seq_diag_mct      ! diagnostic routines

   !--- other ---
   use seq_flds_indices  ! drv aVect field indicies, routines
   use seq_comm_mct      ! mpi comm data & routines, plus logunit and loglevel
   use seq_timemgr_mod   ! clock & alarm routines
   use seq_infodata_mod  ! "infodata" gathers various control flags into one datatype
   use seq_cdata_mod     ! "cdata" type & methods (domain + decomp + infodata in one datatype)
   use seq_domain_mct    ! domain related routines
   use seq_flux_mct      ! flux calc routines
   use seq_frac_mct      ! domain fraction routines
   use seq_rest_mod      ! restart file routines
   use seq_hist_mod      ! history file routines
   use seq_io_mod        ! i/o subroutins

   !--- merging routines ---
!wangty modify
#ifdef wrf
   use mrg_x2w_mct      ! wrf merge gridded component (merging, history) by juanxiong he
   use mrg_x2c_mct      ! cam merge gridded component (merging, history) by juanxiong he
   use mrg_x2ca_mct      ! cam merge gridded component (merging, history) by juanxiong he
   use mrg_x2ge_mct      ! wrf merge gridded component (merging, history) by juanxiong he
#endif
!{list} key="merge" number=1
!   use mrg_x2{c}_mct      ! {ccc} merge gridded component

   !--- mapping routines ---
!wangty modify
#ifdef wrf
   use map_camcam_mct   ! juanxiong he
   use map_wrfwrf_mct   ! juanxiong he
   use map_wrfcam_mct   ! juanxiong he
   use map_gcamcam_mct   ! juanxiong he
   use map_geagea_mct   ! juanxiong he
   use map_geacam_mct   ! juanxiong he
#endif
!{list} key="mapping" number=1
!   use map_{ccc1}{ccc2}_mct   !{ccc1} to {ccc2} mapping routines

   implicit none

   private

   public ccsm_pre_init, ccsm_init, ccsm_run, ccsm_final
#ifdef ESMF_INTERFACE
   public ccsm_comp_register
#endif
   public timing_dir, mpicom_GLOID

#include <mpif.h>

   !----------------------------------------------------------------------------
   ! domains & related
   !----------------------------------------------------------------------------
 
   !--- domain decomps (MCT Global Seg Maps) ---
!{list} key="gsMap_components" number=1
!   type(mct_gsMap)  :: gsMap_{c}{c}   !on component {c} pes

!wangty modify
#ifdef wrf
   type(mct_gsMap)  :: gsMap_ww   !juanxiong he
   type(mct_gsMap)  :: gsMap_cc   !juanxiong he
   type(mct_gsMap)  :: gsMap_mm   !juanxiong he
   type(mct_gsMap)  :: gsMap_gege   !juanxiong he
   type(mct_gsMap)  :: gsMap_caca   !juanxiong he
#endif
!{list} key="gsMap_cpl" number=1
!   type(mct_gsMap)  :: gsMap_{c}x   !on cpl pes

!wangty modify
#ifdef wrf
   type(mct_gsMap)  :: gsMap_cx  !juanxiong he
   type(mct_gsMap)  :: gsMap_wx  !juanxiong he
   type(mct_gsMap)  :: gsMap_mx   !juanxiong he
   type(mct_gsMap)  :: gsMap_cax  !juanxiong he
   type(mct_gsMap)  :: gsMap_gex  !juanxiong he
#endif
   !--- domain area correction factors (only defined on cpl pes) ---
!{list} key="correction_factors" number=1
!   real(r8),pointer :: drv2mdl_{c}{c}(:), mdl2drv_{c}{c}(:)

!wangty modify
#ifdef wrf
   real(r8),pointer :: drv2mdl_ww(:), mdl2drv_ww(:)   ! juanxiong he
   real(r8),pointer :: drv2mdl_ge(:), mdl2drv_ge(:)   ! juanxiong he
#endif

   !--- domain equivalent 2d grid size ---
!{list} key="2dgrid_size" number=1
!   integer          :: {ccc}_nx, {ccc}_ny

!wangty modify
#ifdef wrf
   integer          :: wrf_nx, wrf_ny  ! juanxiong he
   integer          :: geatm_nx, geatm_ny  ! juanxiong he
#endif
   !----------------------------------------------------------------------------
   ! time management
   !----------------------------------------------------------------------------

   type (seq_timemgr_type), SAVE :: seq_SyncClock ! array of all clocks & alarm
   type (ESMF_Clock), SAVE       :: EClock_d      ! driver clock
!{list} key="Eclock" number=1
!  type (ESMF_Clock), SAVE       :: EClock_{c}

!wangty modify
#ifdef wrf
   type (ESMF_Clock)       :: EClock_w  ! juanxiong he for wrf/ccsm4 coupling
   type (ESMF_Clock)       :: EClock_ge  ! juanxiong he for wrf/ccsm4 coupling
#endif
   logical  :: restart_alarm          ! restart alarm
   logical  :: history_alarm          ! history alarm
   logical  :: histavg_alarm          ! history alarm
   logical  :: stop_alarm             ! stop alarm
!{list} key="compoents_alarm" number=1
!  logical  :: {ccc}run_alarm          ! {ccc} run alarm

!{list} key="compoents_nextalarm" number=1
!  logical  :: {ccc}next_alarm         ! {ccc} run alarm on next timestep

!wangty modify
#ifdef wrf
   logical  :: wrfrun_alarm           ! wrf run alarm, juanxiong he
   logical  :: gearun_alarm           ! gea run alarm, juanxiong he
#endif
   logical  :: tprof_alarm            ! timing profile alarm
   logical  :: t1hr_alarm             ! alarm every hour
   logical  :: t2hr_alarm             ! alarm every two hours
   logical  :: t3hr_alarm             ! alarm every three hours
   logical  :: t6hr_alarm             ! alarm every six hours
   logical  :: t12hr_alarm            ! alarm every twelve hours
   logical  :: t24hr_alarm            ! alarm every twentyfour hours

   real(r8) :: days_per_year = 365.0  ! days per year

   integer  :: dtime                  ! dt of one coupling interval
   integer  :: ncpl                   ! number of coupling intervals per day
   integer  :: ymd                    ! Current date (YYYYMMDD)
   integer  :: year                   ! Current date (YYYY)
   integer  :: month                  ! Current date (MM)
   integer  :: day                    ! Current date (DD)
   integer  :: tod                    ! Current time of day (seconds)
   character(CL) :: orb_mode          ! orbital mode
   integer  :: orb_iyear              ! orbital year
   integer  :: orb_iyear_align        ! associated with model year
   integer  :: orb_cyear              ! orbital year for current orbital computation
   integer  :: orb_nyear              ! orbital year associated with currrent model year
   real(r8) :: orb_eccen              ! orbital eccentricity
   real(r8) :: orb_obliq              ! obliquity in degrees
   real(r8) :: orb_mvelp              ! moving vernal equinox long
   real(r8) :: orb_obliqr             ! Earths obliquity in rad
   real(r8) :: orb_lambm0             ! Mean long of perihelion at vernal equinox (radians)
   real(r8) :: orb_mvelpp             ! moving vernal equinox long

   !--- for documenting speed of the model ---
   character( 8) :: dstr              ! date string
   character(10) :: tstr              ! time string
   integer       :: begStep, endStep  ! Begining and ending step number
   real(r8)      :: simDays           ! Number of simulated days
   real(r8)      :: SYPD              ! Simulated years per day
   real(r8)      :: Time_begin        ! Start time
   real(r8)      :: Time_end          ! Ending time
   real(r8)      :: Time_bstep        ! Start time
   real(r8)      :: Time_estep        ! Ending time
   real(r8)      :: dtstep            ! delta time
   real(r8)      :: dtstep_acc        ! dtstep accumulator
   integer       :: dtstep_cnt        ! dtstep counter
   character(CL) :: timing_file       ! Local path to tprof filename
   character(CL) :: timing_dir        ! timing directory
   character(CL) :: tchkpt_dir        ! timing checkpoint directory

   !----------------------------------------------------------------------------
   ! control flags
   !----------------------------------------------------------------------------

!{list} key="components_present" number=1
!  logical  :: {ccc}_present            ! .true.  => {ccc} is present

!wangty modify
#ifdef wrf
   logical  :: wrf_present            ! .true.  => wrf is present, juanxiong he
   logical  :: geatm_present            ! .true.  => gea is present, juanxiong he
#endif

!{list} key="components_prognostic" number=1
!  logical  :: {ccc}_prognostic         ! .true.  => {ccc} comp expects input

!{list} key="prognostic_expect_others" number=1
!  logical  :: {ccc1}{ccc2}_prognostic         ! .true.  => {ccc1} comp expects {ccc2} input

!wangty modify
#ifdef wrf
   logical  :: wrf_prognostic         ! .true.  => wrf comp expects input, juanxiong he
   logical  :: geatm_prognostic         ! .true.  => gea comp expects input, juanxiong he
#endif
   logical  :: dead_comps             ! .true.  => dead components

   logical  :: single_column          ! scm mode logical
   real(r8) :: scmlon                 ! single column lon
   real(r8) :: scmlat                 ! single column lat
   logical  :: aqua_planet            ! aqua planet mode
   real(r8) :: nextsw_cday            ! radiation control
!{list} key="aerosol_data" number=1
!   logical  :: {ccc}_aero               ! {ccc} provides aerosol data

   real(r8) :: flux_epbalfact         ! precip factor

!<list> key="ocn_def" number=1
   logical  :: {ccccc}_tight_coupling   ! couple {ccc} as frequently as lnd & ice
   logical  :: skip_{ccccc}_run         ! skip the {ccccc} model first pass
   logical  :: cpl2{ccc}_first          ! use to call initial cpl2{ccc} timer
   character(CS) :: a{c}flux_grid       ! grid for a/{c} flux calc: atm xor {ccc}
!</list>

   logical  :: run_barriers           ! barrier the component run calls


   logical       :: read_restart      ! local read restart flag
   character(CL) :: rest_file         ! restart file path + filename

   logical  :: domain_check           ! .true.  => check consistency of domains
   logical  :: shr_map_dopole         ! logical for dopole in shr_map_mod

   !--- history & budgets ---
   logical       :: do_budgets        ! heat/water budgets on
   logical       :: do_histinit       ! initial hist file
   logical       :: do_histavg        ! histavg on or off
!<list> key="aux_files" number=1
   logical       :: do_hist_{c}2x       ! create aux files: {c}2x
   character(CL) :: hist_{c}2x_flds     = '{value}'
!</list>

!<list> key="aux_files_24hr" number=1
   logical       :: do_hist_{c}2x24hr   ! create aux files: {c}2x 24hr
   character(CL) :: hist_{c}2x24hr_flds     = '{value}'
!</list>

!<list> key="aux_files_3hr" number=1
   logical       :: do_hist_{c}2x3hrp   ! create aux files: {c}2x 3hr precip
   logical       :: do_hist_{c}2x3hr    ! create aux files: {c}2x 3hr states
   character(CL) :: hist_{c}2x3hrp_flds     = '{value}'
   character(CL) :: hist_{c}2x3hr_flds     = '{value}'
!</list>

   integer  :: budget_inst            ! instantaneous budget flag
   integer  :: budget_daily           ! daily budget flag
   integer  :: budget_month           ! monthly budget flag
   integer  :: budget_ann             ! annual budget flag
   integer  :: budget_ltann           ! long term budget flag for end of year writing
   integer  :: budget_ltend           ! long term budget flag for end of run writing

   ! --- other ---
   integer  :: k1,k2,k3               ! aVect field indices
!{list} key="alarm_count" number=1
!   integer  :: {ccc}run_count           ! number of times {ccc} run alarm went on

   logical  :: exists                 ! true if file exists
   integer  :: ierr                   ! MPI error return
   integer  :: rc                     ! return code
   logical  :: cdf64                  ! true => use 64 bit addressing in netCDF files

   character(*), parameter :: NLFileName = "drv_in"  ! input namelist filename

   integer  :: info_debug = 0         ! local info_debug level

   !----------------------------------------------------------------------------
   ! memory monitoring
   !----------------------------------------------------------------------------
   real(r8) :: msize,msize0,msize1     ! memory size (high water)
   real(r8) :: mrss ,mrss0 ,mrss1      ! resident size (current memory use)

   !----------------------------------------------------------------------------
   ! threading control
   !----------------------------------------------------------------------------
   integer  :: nthreads_GLOID         ! OMP global number of threads
   integer  :: nthreads_CPLID         ! OMP cpl number of threads
!{list} key="nthreads_ID" number=1
!   integer  :: nthreads_{ccc}ID         ! OMP {ccc} number of threads

!wangty modify
#ifdef wrf
   integer  :: nthreads_WRFID         ! OMP wrf number of threads, juanxiong he
   integer  :: nthreads_GEAID         ! OMP wrf number of threads, juanxiong he
#endif
   integer  :: pethreads_GLOID        ! OMP number of threads per task

!{list} key="nthreads_CPLID" number=1
!   integer  :: nthreads_CPL{ccc}ID      ! OMP cpl-{ccc} number of threads

   logical  :: drv_threading          ! driver threading control

   !----------------------------------------------------------------------------
   ! communicator groups and related
   !----------------------------------------------------------------------------
   integer  :: Global_Comm
   integer  :: mpicom_GLOID           ! MPI global communicator
   integer  :: mpicom_CPLID           ! MPI cpl communicator
!{list} key="mpicom" number=1
!   integer  :: mpicom_{ccc}ID           ! MPI {ccc} communicator

!wangty modify
#ifdef wrf
   integer  :: mpicom_WRFID           ! MPI wrf communicator, juanxiong he
   integer  :: mpicom_GEAID           ! MPI wrf communicator, juanxiong he
#endif
!{list} key="mpicom_cpl" number=1
!   integer  :: mpicom_CPL{ccc}ID        ! MPI cpl-{ccc} communicator

!wangty modify
#ifdef wrf
   integer  :: mpicom_CPLWRFID        ! MPI cpl-wrf communicator, juanxiong he
   integer  :: mpicom_CPLGEAID        ! MPI cpl-wrf communicator, juanxiong he
#endif
   logical  :: iamroot_GLOID          ! GLOID masterproc
   logical  :: iamroot_CPLID          ! CPLID masterproc
!{list} key="masterproc" number=1
!   logical  :: iamroot_{ccc}ID          ! {ccc}ID masterproc

!wangty modify
#ifdef wrf
   logical  :: iamroot_WRFID          ! WRFID masterproc, juanxiong he
   logical  :: iamroot_GEAID          ! GEAID masterproc, juanxiong he
#endif
   logical  :: iamin_CPLID            ! pe associated with CPLID
!{list} key="iamin" number=1
!   logical  :: iamin_{ccc}ID            ! pe associated with {ccc}ID

!wangty modify
#ifdef wrf
   logical  :: iamin_WRFID            ! pe associated with WRFID, juanxiong he
   logical  :: iamin_GEAID            ! pe associated with WRFID, juanxiong he
#endif
!{list} key="iamin_cpl" number=1
!   logical  :: iamin_CPL{ccc}ID         ! pe associated with CPL{ccc}ID

!wangty modify
#ifdef wrf
   logical  :: iamin_CPLWRFID         ! pe associated with CPLWRFID, juanxiong he
   logical  :: iamin_CPLGEAID         ! pe associated with CPLWRFID, juanxiong he
#endif
   character(CL) :: complist          ! list of comps on this pe
   integer  :: iam_GLOID              ! pe number in global id
!wangty modify
#ifdef wrf
   integer, pointer :: atm_petlist(:), lnd_petlist(:), ice_petlist(:), ocn_petlist(:), glc_petlist(:), &
                       wrf_petlist(:), gea_petlist(:) ! juanxiong he
#else
!{list} key="petlist" number=1
!  integer, pointer :: {ccc}_petlist(:)

#endif
   !----------------------------------------------------------------------------
   ! formats
   !----------------------------------------------------------------------------
   character(*), parameter :: subname = '(seq_mct_drv)'
   character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
   character(*), parameter :: F0L = "('"//subname//" : ', A, L6 )"
   character(*), parameter :: F0I = "('"//subname//" : ', A, 2i8 )"
   character(*), parameter :: F0R = "('"//subname//" : ', A, 2g23.15 )"
   character(*), parameter :: FormatA = '(A,": =============== ", A41, " ===============")'
   character(*), parameter :: FormatD = '(A,": =============== ", A20,2I8,5x,   " ===============")'
   character(*), parameter :: FormatR = '(A,": =============== ", A31,F9.3,1x,  " ===============")'
   character(*), parameter :: FormatQ = '(A,": =============== ", A20,2F10.2,1x," ===============")'
!wangty modify
#ifdef wrf
   logical :: twoway_coupling ! juanxiong he
   integer :: twoway_nudging ! juanxiong he
   integer :: integration_phase ! juanxiong he
#endif
!===============================================================================
contains
!===============================================================================

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

subroutine ccsm_pre_init()

   !--------------------------------------------------------------------------
   ! Initialize MCT and MPI communicators and IO
   !--------------------------------------------------------------------------

   call mpi_init(ierr)
   call shr_mpi_chkerr(ierr,subname//' mpi_init')

   Global_Comm=MPI_COMM_WORLD
   call seq_io_init1(NLFileName, Global_Comm)
   if (Global_Comm /= MPI_COMM_NULL) then
!wangty modify
#ifdef wrf
        call seq_comm_init(Global_Comm, NLFileName, atm_petlist=atm_petlist, lnd_petlist=lnd_petlist, &
           ice_petlist=ice_petlist, ocn_petlist=ocn_petlist, glc_petlist=glc_petlist, &
           wrf_petlist=wrf_petlist, gea_petlist=gea_petlist) ! juanxiong he
#else
!<list> key="petlist" number=1
     call seq_comm_init(Global_Comm, NLFileName, {ccc}_petlist={ccc}_petlist, &
     )
!</list>

#endif
   end if

   call seq_io_init2()

   !--- set task based threading counts ---
   call seq_comm_setptrs(GLOID,pethreads=pethreads_GLOID,iam=iam_GLOID)
   call seq_comm_setnthreads(pethreads_GLOID)

   !--- get some general data ---
   call seq_comm_setptrs(GLOID,mpicom=mpicom_GLOID,iamroot=iamroot_GLOID,nthreads=nthreads_GLOID)
   call seq_comm_setptrs(CPLID,mpicom=mpicom_CPLID,iamroot=iamroot_CPLID,nthreads=nthreads_CPLID)
!{list} key="mpicom" number=1
!call seq_comm_setptrs({CCC}ID,mpicom=mpicom_{CCC}ID,iamroot=iamroot_{CCC}ID,nthreads=nthreads_{CCC}ID)

!wangty modify
#ifdef wrf
   call seq_comm_setptrs(WRFID,mpicom=mpicom_WRFID,iamroot=iamroot_WRFID,nthreads=nthreads_WRFID) ! juanxiong he
   call seq_comm_setptrs(GEAID,mpicom=mpicom_GEAID,iamroot=iamroot_GEAID,nthreads=nthreads_GEAID) ! juanxiong he
#endif
!{list} key="mpicom_cpl" number=1
!   call seq_comm_setptrs(CPL{ccc}ID,mpicom=mpicom_CPL{ccc}ID,nthreads=nthreads_CPL{ccc}ID)

!wangty modify
#ifdef wrf
   call seq_comm_setptrs(CPLWRFID,mpicom=mpicom_CPLWRFID,nthreads=nthreads_CPLWRFID) !juanxiong he
   call seq_comm_setptrs(CPLGEAID,mpicom=mpicom_CPLGEAID,nthreads=nthreads_CPLGEAID) !juanxiong he
#endif
   iamin_CPLID    = seq_comm_iamin(CPLID)
!{list} key="iamin" number=1
!   iamin_{ccc}ID    = seq_comm_iamin({ccc}ID)

!wangty modify
#ifdef wrf
   iamin_CPLWRFID = seq_comm_iamin(CPLWRFID) !juanxiong he
   iamin_CPLGEAID = seq_comm_iamin(CPLGEAID) !juanxiong he
#endif
   complist = " "
   if (iamin_CPLID) complist = trim(complist)//' cpl'
!{list} key="iamin" number=1
!   if (iamin_{ccc}ID) complist = trim(complist)//' {ccc}'

!wangty modify
#ifdef wrf
   if (iamin_WRFID) complist = trim(complist)//' wrf' !juanxiong he
   if (iamin_GEAID) complist = trim(complist)//' gea' !juanxiong he
#endif
   !--------------------------------------------------------------------------
   ! Set logging parameters both for shr code and locally
   !--------------------------------------------------------------------------

   if (iamroot_CPLID) then
      inquire(file='cpl_modelio.nml',exist=exists)
      if (exists) then
         logunit = shr_file_getUnit()
         call shr_file_setIO('cpl_modelio.nml',logunit)
         call shr_file_setLogUnit(logunit)
         loglevel = 1
         call shr_file_setLogLevel(loglevel)
      endif
   else
      loglevel = 0
      call shr_file_setLogLevel(loglevel)
   endif

   !--------------------------------------------------------------------------
   ! Log info about the environment settings
   !--------------------------------------------------------------------------

   if (iamroot_CPLID) then
#ifdef USE_ESMF_LIB
      write(logunit,'(2A)') subname,' USE_ESMF_LIB is set'
#else
      write(logunit,'(2A)') subname,' USE_ESMF_LIB is NOT set, using esmf_wrf_timemgr'
#endif
#ifdef MCT_INTERFACE
      write(logunit,'(2A)') subname,' MCT_INTERFACE is set' 
#endif
#ifdef ESMF_INTERFACE
      write(logunit,'(2A)') subname,' ESMF_INTERFACE is set'
#endif
   endif


end subroutine ccsm_pre_init

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

subroutine ccsm_init()

 101  format( A, 2i8, 12A, A, F8.2, A, F8.2 )
 102  format( A, 2i8, A, 5L3 )
 103  format( 5A )
 104  format( A, 2i8)
 105  format( A, 2i8, A, f10.2, A, f10.2, A, A, i5, A, A)
 106  format( A, f23.12)

   !--------------------------------------------------------------------------
   ! Print Model heading and copyright message
   !--------------------------------------------------------------------------

   if (iamroot_CPLID) call seq_ccsm_printlogheader()

   !-----------------------------------------------------------------------------
   ! Timer initialization (has to be after mpi init)
   !-----------------------------------------------------------------------------

   call t_initf(NLFileName, LogPrint=.false., mpicom=mpicom_GLOID, &
                MasterTask=iamroot_GLOID)

   if (iamin_CPLID) then
      call seq_io_cpl_init()
   endif

   call t_startf('DRIVER_INIT')

   !-----------------------------------------------------------------------------
   ! Memory test
   !-----------------------------------------------------------------------------
   call shr_mem_init(prt=.true.)

   !-----------------------------------------------------------------------------
   ! Initialize coupled field indices
   !-----------------------------------------------------------------------------

   call seq_flds_set()
   call seq_flds_indices_set( )

   !-----------------------------------------------------------------------------
   ! Initialize infodata
   !-----------------------------------------------------------------------------

   call seq_infodata_init(infodata,nlfilename,GLOID)
   if (iamroot_CPLID) then
      write(logunit,*) ' '
      write(logunit,'(2A)') 'Status of infodata after seq_infodata_init'
      call seq_infodata_print( infodata )
      write(logunit,*) ' '
   endif

   call seq_infodata_GetData(infodata,read_restart=read_restart, restart_file=rest_file, &
        timing_dir=timing_dir, tchkpt_dir=tchkpt_dir)
!wangty modify
#ifdef wrf
   call seq_infodata_GetData(infodata, info_debug=info_debug, atm_present=atm_present, &
        lnd_present=lnd_present, ice_present=ice_present, ocn_present=ocn_present, &
        glc_present=glc_present, sno_present=sno_present, wrf_present=wrf_present, & !juanxiong he
        geatm_present=geatm_present, & !juanxiong he
        single_column=single_column, aqua_planet=aqua_planet, &
        ocean_tight_coupling=ocean_tight_coupling, drv_threading=drv_threading)
#else
!{list} key="components_present" number=1
!   call seq_infodata_GetData(infodata, info_debug=info_debug, {ccc}_present={ccc}_present, &

#endif
   call seq_infodata_GetData(infodata, do_histinit=do_histinit)
   call seq_infodata_GetData(infodata, do_budgets=do_budgets, budget_inst=budget_inst, &
        budget_daily=budget_daily, budget_month=budget_month, budget_ann=budget_ann, &
        budget_ltann=budget_ltann, budget_ltend=budget_ltend)

!<list> key="getdata" number=1
   call seq_infodata_GetData(infodata, &
        histaux_{c1}2x    =do_hist_{c1}2x    , histaux_{c3}2x3hr =do_hist_{c3}2x3hr , &
        histaux_{c3}2x3hrp=do_hist_{c3}2x3hrp, histaux_{c2}2x24hr=do_hist_{c2}2x24hr, &
       )
!</list>

   call seq_infodata_GetData(infodata, run_barriers = run_barriers)

   call seq_infodata_GetData(infodata, aoflux_grid=aoflux_grid)

   call seq_infodata_GetData(infodata, shr_map_dopole=shr_map_dopole)
   call shr_map_setDopole(shr_map_dopole)

   !-----------------------------------------------------------------------------
   ! Test Threading Setup in driver, happens to be valid on all pes for all IDs
   !-----------------------------------------------------------------------------

   if (drv_threading) then
      if (iamroot_GLOID) write(logunit,*) ' '
      if (iamroot_GLOID) write(logunit,'(2A)    ') subname,' Test Threading in driver'
      call seq_comm_setnthreads(nthreads_GLOID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLOID = ',nthreads_GLOID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_CPLID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_CPLID = ',nthreads_CPLID,seq_comm_getnthreads()
!{list} key="nthreads_ID"
!      call seq_comm_setnthreads(nthreads_{ccc}ID)
!      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_{ccc}ID = ',nthreads_{ccc}ID,seq_comm_getnthreads()
      if (iamroot_GLOID) write(logunit,*) ' '
!wangty modify
#ifdef wrf
      call seq_comm_setnthreads(nthreads_WRFID) ! juanxiong he
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_WRFID = ',nthreads_WRFID,seq_comm_getnthreads()
      if (iamroot_GLOID) write(logunit,*)
      call seq_comm_setnthreads(nthreads_GEAID) ! juanxiong he
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GEAID =',nthreads_GEAID,seq_comm_getnthreads()
      if (iamroot_GLOID) write(logunit,*)
#endif
      call seq_comm_setnthreads(nthreads_GLOID)
   endif

   !-----------------------------------------------------------------------------
   ! Setup cdata types, call on all pes so the ID is set on all pes even
   ! though other data may be invalid
   !-----------------------------------------------------------------------------
!{list} key="gsMap_cpl"
!   call seq_cdata_init(cdata_{c}x, CPLID, dom_{c}x, gsMap_{c}x, infodata, 'cdata_{c}x' )
!wangty modify
#ifdef wrf
   call seq_cdata_init(cdata_wx, CPLID, dom_wx, gsMap_wx, infodata, 'cdata_wx' ) !juanxiong he, wrf/ccsm4 coupling
   call seq_cdata_init(cdata_mx, CPLID, dom_mx, gsMap_mx, infodata, 'cdata_mx' ) !juanxiong he, wrf/cam coupling, wrf end
   call seq_cdata_init(cdata_cx, CPLID, dom_cx, gsMap_cx, infodata, 'cdata_cx' ) !juanxiong he, wrf/cam coupling, cam end
   call seq_cdata_init(cdata_gex, CPLID, dom_gex, gsMap_gex, infodata, 'cdata_gex' )  !juanxiong he, geatm/ccsm4 coupling
   call seq_cdata_init(cdata_cax, CPLID, dom_cax, gsMap_cax, infodata, 'cdata_cax' )  !juanxiong he, geatm/cam coupling, cam end
#endif

!{list} key="gsMap_components"
!   call seq_cdata_init(cdata_{cc}, ATMID, dom_{cc}, gsMap_{cc}, infodata, 'cdata_{cc}')
!wangty modify
#ifdef wrf
   call seq_cdata_init(cdata_ww, WRFID, dom_ww, gsMap_ww, infodata, 'cdata_ww' ) !juanxiong he, wrf/ccsm4 coupling
   call seq_cdata_init(cdata_mm, WRFID, dom_mm, gsMap_mm, infodata, 'cdata_mm' ) !juanxiong he, wrf/cam coupling, wrf end
   call seq_cdata_init(cdata_cc, ATMID, dom_cc, gsMap_cc, infodata, 'cdata_cc' ) !juanxiong he, wrf/cam coupling, cam end
   call seq_cdata_init(cdata_gege, GEAID, dom_gege, gsMap_gege, infodata, 'cdata_gege' )   !juanxiong he, geatm/cam coupling, geatm end
   call seq_cdata_init(cdata_caca, ATMID, dom_caca, gsMap_caca, infodata, 'cdata_caca' )   !juanxiong he, geatm/cam coupling, cam end
#endif
   !-----------------------------------------------------------------------------
   ! Initialize time manager
   !-----------------------------------------------------------------------------
!wangty modify
#ifdef wrf
   call seq_timemgr_clockInit(seq_SyncClock,nlfilename,read_restart,rest_file,mpicom_gloid, &
      EClock_d, EClock_a, EClock_l, EClock_o, EClock_i, Eclock_g, Eclock_w, Eclock_ge)  ! juanxiong he for wrf/ccsm4 coupling
#else
!{list} key="Eclock"
!   call seq_timemgr_clockInit(seq_SyncClock,nlfilename,read_restart,rest_file,mpicom_gloid, EClock_d, &
!      EClock_{c}, &
!   )
#endif
   if (iamroot_CPLID) then
       call seq_timemgr_clockPrint(seq_SyncClock)
   endif

   call seq_infodata_getData(infodata,orb_iyear=orb_iyear,orb_iyear_align=orb_iyear_align, &
      orb_mode=orb_mode)
   if (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd)
      call shr_cal_date2ymd(ymd,year,month,day)
      orb_cyear = orb_iyear + (year - orb_iyear_align)
      call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                          orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
      call seq_infodata_putData(infodata,orb_eccen=orb_eccen,orb_obliqr=orb_obliqr, &
           orb_lambm0=orb_lambm0,orb_mvelpp=orb_mvelpp)
   endif
!wangty modify
#ifdef wrf
   call seq_infodata_putData(infodata,atm_phase=1,lnd_phase=1,ocn_phase=1,ice_phase=1,glc_phase=1,&
                             wrf_phase=1, geatm_phase=1) !juanxiong he
#else
!<list> key="components" number=1
   call seq_infodata_putData(infodata, &
      {ccc}_phase=1, &
   )
!</list>

#endif
   !-----------------------------------------------------------------------------
   ! If in single column mode overwrite lnd,ocn,ice_present flags according to
   ! focndomain file in ocn_in namelist
   ! SCAM can reset the lnd_present, ice_present and ocn_present flags
   !-----------------------------------------------------------------------------

   if (.not.aqua_planet .and. single_column) then
      call seq_infodata_getData( infodata, scmlon=scmlon, scmlat=scmlat)
      call shr_scam_checkSurface(scmlon, scmlat, OCNID, mpicom_OCNID, &
!{list} key="overwrite_present"
!           {ccc}_present={ccc}_present, &
!      )
!      call seq_infodata_putData( infodata, & 
!           {ccc}_present={ccc}_present, &
!      )
   endif

   !-----------------------------------------------------------------------------
   ! Component Initialization
   ! Note that within each component initialization, the relevant x_pxresent
   ! flag
   ! part of CCSMInit (contained as a pointer in cdata_xc) can be modified
   ! By default, all these flags are set to true
   ! The atm can reset the lnd_present, ice_present and ocn_present flags based
   ! on aqua_planet, ideal_phys and adiabatic modes
   ! The stub components will reset the present flags to false, all other
   ! components will set them to true for the purposes of symmetry
   !-----------------------------------------------------------------------------

   call t_startf('driver_init_comps')
   if ( iamroot_CPLID )then
      write(logunit,*) ' '
!wangty modify
#ifdef wrf
      write(logunit,F00) 'Initialize each component: atm, wrf, gea, lnd, ocn, and ice' ! juanxiong he
#else
!{list} key="components"
!      write(logunit,F00) 'Initialize each component: {ccc}'
#endif
      call shr_sys_flush(logunit)
   endif
   call t_adj_detailf(+2)

!{list} key='component_init' number=1
!   !-----------------------------------------------------------------------------
!   ! Initialization {ccc} component
!   !-----------------------------------------------------------------------------


!   if (iamin_CPL{ccc}ID) then
!      call seq_infodata_exchange(infodata,CPL{ccc}ID,'cpl2{ccc}_init')
!   endif
!   if (iamin_{ccc}ID .and. {ccc}_present) then
!      if (seq_comm_iamroot({ccc}ID)) write(logunit,F00) 'Initialize {ccc} component'
!      call shr_sys_flush(logunit)
!      if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!!wangty modify
!#ifdef wrf
!      call atm_init_mct( EClock_a, cdata_aa, cdata_cc, cdata_caca, x2a_aa, a2x_aa, &
!                         x2c_cc1, x2c_cc2, c2x_cc1, c2x_cc2, &  ! juanxiong he for wrf
!                         x2ca_caca1, x2ca_caca2, ca2x_caca1, ca2x_caca2, & ! juanxiong he for geatm
!                         twoway_coupling, twoway_nudging, NLFilename=NLFilename ) !juanxiong he
!#else
!      call {ccc}_init_mct( EClock_{c}, cdata_{cc}, x2{c}_{cc}, {c}2x_{cc}, NLFilename=NLFilename )
!#endif
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!   endif
!   if (iamin_CPL{ccc}ID) then
!      call seq_infodata_exchange(infodata,CPL{ccc}ID,'{ccc}2cpl_init')
!   endif
!wangty modify
#ifdef wrf
   !-----------------------------------------------------------------------------
   ! Initialization WRF component, juanxiong he
   !-----------------------------------------------------------------------------
   if (iamin_CPLWRFID) then
      call seq_infodata_exchange(infodata,CPLWRFID,'cpl2wrf_init')
   endif
   if (iamin_WRFID .and. wrf_present) then
      if (seq_comm_iamroot(WRFID)) write(logunit,F00) 'Initialize wrf component'
      call shr_sys_flush(logunit)
      if (drv_threading) call seq_comm_setnthreads(nthreads_WRFID)
      call wrf_init_mct( EClock_w, cdata_ww, cdata_mm, x2w_ww, w2x_ww, x2m_mm1, x2m_mm2, &
                         m2x_mm, twoway_coupling, twoway_nudging, NLFilename=NLFilename )
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif
   if (iamin_CPLWRFID) then
      call seq_infodata_exchange(infodata,CPLWRFID,'wrf2cpl_init')
   endif

   !-----------------------------------------------------------------------------
   ! Initialization GEATM component, juanxiong he
   !-----------------------------------------------------------------------------

   if (iamin_CPLGEAID) then
      call seq_infodata_exchange(infodata,CPLGEAID,'cpl2gea_init')
   endif
   if (iamin_GEAID .and. geatm_present) then
      if (seq_comm_iamroot(GEAID)) write(logunit,F00) 'Initialize gea component'
      call shr_sys_flush(logunit)
      if (drv_threading) call seq_comm_setnthreads(nthreads_GEAID)
      call geatm_init_mct( EClock_ge, cdata_gege, &
                           x2chem_chemchem1,x2chem_chemchem2, chem2x_chemchem )
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif
   if (iamin_CPLGEAID) then
      call seq_infodata_exchange(infodata,CPLGEAID,'gea2cpl_init')
   endif
#endif
!{list} key="component_init" number=2
!   !-----------------------------------------------------------------------------
!   ! Initialization {ccc} component
!   !-----------------------------------------------------------------------------

!   if (iamin_CPL{ccc}ID) then
!      call seq_infodata_exchange(infodata,CPL{ccc}ID,'cpl2{ccc}_init')
!   endif
!   if (iamin_{ccc}ID .and. {ccc}_present) then
!      if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!      if (seq_comm_iamroot({ccc}ID)) write(logunit,F00) 'Initialize {ccc} component'
!      call shr_sys_flush(logunit)
!      call lnd_init_mct({allccc})
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!    endif
!   if (iamin_CPL{ccc}ID) then
!      call seq_infodata_exchange(infodata,CPL{ccc}ID,'{ccc}2cpl_init')
!   endif

   call t_adj_detailf(-2)

   call t_stopf  ('driver_init_comps')

   !-----------------------------------------------------------------------------
   ! Determine final settings for presence of land, ice and ocean and the
   ! prognostic flags
   !-----------------------------------------------------------------------------

!{list} key="iamin_cpl"
!   if (iamin_CPL{ccc}ID) then
!      call seq_infodata_exchange(infodata,CPL{ccc}ID,'cpl2{ccc}_init')
!   endif

!wangty modify
#ifdef wrf
   if (iamin_CPLWRFID) then  !juanxiong he
      call seq_infodata_exchange(infodata,CPLWRFID,'cpl2wrf_init')
   endif
   if (iamin_CPLGEAID) then  !juanxiong he
      call seq_infodata_exchange(infodata,CPLGEAID,'cpl2gea_init')
   endif
#endif
   if ( iamroot_CPLID) then
      write(logunit,F00) 'Determine final settings for presence of surface components'
      call shr_sys_flush(logunit)
   endif
  
   call seq_infodata_getData(infodata, &
!{list} key="components_present"
!        {ccc}_present={ccc}_present, &
!wangty modify
#ifdef wrf
        wrf_present=wrf_present, & !juaxniong he
        geatm_present=geatm_present, & !juaxniong he
#endif
!{list} key="components_prognostic"
!        {ccc}_prognostic={ccc}_prognostic, &
!wangty modify
#ifdef wrf
        wrf_prognostic=wrf_prognostic, & !juanxiong he
        geatm_prognostic=geatm_prognostic, & !juanxiong he
#endif
!{list} key="prognostic_expect_others"
!        {ccc1}{ccc2}_prognostic={ccc1}{ccc2}_prognostic, &
        dead_comps=dead_comps, &
!{list} key="2dgrid_size"
!        {ccc}_nx={ccc}_nx, {ccc}_ny={ccc}_ny, &
!wangty modify
#ifdef wrf
        wrf_nx=wrf_nx, wrf_ny=wrf_ny, & ! juanxiong he
        geatm_nx=geatm_nx, geatm_ny=geatm_ny, & ! juanxiong he
#endif
        cpl_cdf64=cdf64, &
!{list} key="aerosol_data"
!        {ccc}_aero={ccc}_aero 
!        )

!{list} key="present_test"
!   if (.not. {ccc}_present) then
!      call shr_sys_abort('{ccc} must be present')
!   endif
!{list} key="prognostic_expect_others"
!   if ({ccc1}{ccc2}_prognostic .and. .not.{ccc2}_present) then
!      if (iamroot_CPLID) then
!         write(logunit,F00) 'WARNING: {ccc1}{ccc2}_prognostic is TRUE but {ccc2}_present is FALSE'
!         call shr_sys_flush(logunit)
!      endif
!   endif
!{list} key="present_prognostic_test"
!   if ({ccc}_prognostic .and. .not.{ccc}_present) then
!      call shr_sys_abort('if prognostic {ccc} must also have ocn present')
!   endif
!wangty modify
#ifdef wrf
   if (wrf_prognostic .and. .not.wrf_present) then  !juanxiong he
      call shr_sys_abort('if prognostic wrf must also have wrf present')
   endif
   if (geatm_prognostic .and. .not.geatm_present) then  !juanxiong he
      call shr_sys_abort('if prognostic geatm must also have geatm present')
   endif
#endif
! tcx remove temporarily for development
!   if (glc_prognostic .and. .not.sno_present) then
!      call shr_sys_abort('if prognostic glc must also have sno present')
!   endif
!   if (sno_prognostic .and. .not.glc_present) then
!      call shr_sys_abort('if prognostic sno must also have glc present')
!   endif

   !-----------------------------------------------------------------------------
   ! Set domain check and other flag
   !-----------------------------------------------------------------------------

   domain_check = .true.
   if (single_column         ) domain_check = .false.
   if (dead_comps            ) domain_check = .false.

   ! set skip_ocean_run flag, used primarily for ocn run on first timestep
   ! use reading a restart as a surrogate from whether this is a startup run

!{list} key="ocn_def"
!   skip_{ccccc}_run = .true.
!   if ( read_restart) skip_{ccccc}_run = .false.
!   {ccc}run_count = 0
!   cpl2{ccc}_first = .true.

   do_histavg = .true.
   if (seq_timemgr_histavg_type == seq_timemgr_type_never) then
      do_histavg = .false.
   endif

   !-----------------------------------------------------------------------------
   ! Write output
   ! NOTE- assume that runoff will only be mapped from land to ocean if
   !       prognostic ocean is true
   !-----------------------------------------------------------------------------

   if (iamroot_CPLID) then
      write(logunit,*  )' '
      write(logunit,F00)'After component initialization:'
!{list} key="components_present"
!      write(logunit,F0L)'{ccc} model present     = ',{ccc}_present
!wangty modify
#ifdef wrf
      write(logunit,F0L)'wrf model present     = ',wrf_present ! juanxiong he
      write(logunit,F0L)'geatm  model present  = ',geatm_present ! juanxiong he
#endif
!{list} key="components_prognostic"
!      write(logunit,F0L)'{ccc} model prognostic  = ',{ccc}_prognostic
!wangty modify
#ifdef wrf
      write(logunit,F0L)'wrf model prognostic  = ',wrf_prognostic ! juanxiong he
      write(logunit,F0L)'geatm model prognostic  = ',geatm_prognostic ! juanxiong he
#endif
!{list} key="prognostic_expect_others"
!      write(logunit,F0L)'{ccc1} {ccc2}   prognostic  = ',{ccc1}{ccc2}_prognostic
      write(logunit,F0L)'dead components       = ',dead_comps
      write(logunit,F0L)'domain_check          = ',domain_check
!{list} key="2dgrid_size"
!      write(logunit,F0I)'{ccc}_nx,{ccc}_ny         = ',{ccc}_nx,{ccc}_ny
!wangty modify
#ifdef wrf
      write(logunit,F0I)'wrf_nx,wrf_ny         = ',wrf_nx,wrf_ny !juanxiong he
      write(logunit,F0I)'geatm_nx,geatm_ny         = ',geatm_nx,geatm_ny !juanxiong he
#endif
!{list} key="ocn_def"
!      write(logunit,F0L)'skip init {ccccc} run   = ',skip_{ccccc}_run
!      write(logunit,F0L)'{ccccc} tight coupling  = ',{cccc}_tight_coupling
      write(logunit,F0L)'cpl_cdf64             = ',cdf64
      write(logunit,F0L)'do_histavg            = ',do_histavg
!{list} key="aerosol_data"
!      write(logunit,F0L)'{ccc}_aero              = ',{ccc}_aero
      write(logunit,*  )' '
      call shr_sys_flush(logunit)
   endif

   !-----------------------------------------------------------------------------
   ! Need to initialize aream, set it to area for now until maps are read
   !   in some cases, maps are not read at all !!
   ! Need to initialize ascale
   ! Entire domain must have reasonable values before calling xxx2xxx init
   ! NOTE (tcx) : use cdata%dom instead of dom% due to seg fault on bluevista I,
   ! why?
   !-----------------------------------------------------------------------------
!{list} key="init_aream_ascale" number=1
!   if (iamin_{ccc}ID .and. {ccc}_present) then
!      if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!      k1 = mct_aVect_indexRa(cdata_{c}{c}%dom%data,"area"  ,perrWith='{c}{c} area ')
!      k2 = mct_aVect_indexRa(cdata_{c}{c}%dom%data,"aream" ,perrWith='{c}{c} aream')
!      k3 = mct_aVect_indexRa(cdata_{c}{c}%dom%data,"ascale",perrWith='{c}{c} ascale')
!      cdata_{c}{c}%dom%data%rAttr(k2,:) = cdata_{c}{c}%dom%data%rAttr(k1,:)
!      cdata_{c}{c}%dom%data%rAttr(k3,:) = 1.0_r8
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!   endif
!{list} key="init_aream_ascale" number=2
!    if (iamin_{comp}ID .and. {ccc}_present) then
!      if (drv_threading) call seq_comm_setnthreads(nthreads_{comp}ID)
!      k1 = mct_aVect_indexRa(cdata_{c}{c}%dom%data,"area"  ,perrWith='{c}{c} area ')
!      k2 = mct_aVect_indexRa(cdata_{c}{c}%dom%data,"aream" ,perrWith='{c}{c} aream')
!      k3 = mct_aVect_indexRa(cdata_{c}{c}%dom%data,"ascale",perrWith='{c}{c} ascale')
!      cdata_{c}{c}%dom%data%rAttr(k2,:) = cdata_{c}{c}%dom%data%rAttr(k1,:)
!      cdata_{c}{c}%dom%data%rAttr(k3,:) = 1.0_r8
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!   endif
   !-----------------------------------------------------------------------------
   ! Initialize driver rearrangers and AVs on driver
   ! Initialize cdata_*x data
   ! Zero out x2*_** in case it never gets used then it'll produce zeros in
   ! diags
   !-----------------------------------------------------------------------------
   call t_startf('driver_init_xxx2xxx')
!{list} key="init_driver" number=1
!   if (iamin_CPL{ccc}ID .and. {ccc}_present) then
!      call map_{ccc}2{ccc}_init_mct({allccc1})
!      call mct_avect_zero(x2{c}_{c}{c})
!      call map_{ccc}{c}2{ccc}x_mct({allccc2})
!   endif
!{list} key="init_driver" number=2
!   if (iamin_CPL{comp}ID .and. {ccc}_present) then
!      call map_{ccc}2{ccc}_init_mct({allccc1})
!      call mct_avect_init({c}2xacc_{c}x%data, [c]2x_{c}x, mct_aVect_lsize({c}2x_{c}x))
!      call mct_accum_zero({c}2xacc_{c}x)
!      {c}2xacc_{c}x_cnt = 0
!   endif
!{list} key="init_driver" number=3
!   if (iamin_CPL{ccc}ID .and. {ccc}_present) then
!      call map_{ccc}2{ccc}_init_mct({allccc1})
!      call mct_avect_zero(x2{c}_{cc})
!      call map_{ccc}{c}2{ccc}x_mct(all{ccc2})
!      call mct_avect_init(x2{c}acc_{c}x%data, x2{c}_{c}x, mct_aVect_lsize(x2{c}_{c}x))
!      call mct_accum_zero(x2{c}acc_{c}x)
!      x2{c}acc_{c}x_cnt = 0
!   endif

!wangty modify
#ifdef wrf
   !-----------------------------------------------------------------------------
   ! GEATM and WRF, Juanxiong He
   !-----------------------------------------------------------------------------

   if (iamin_CPLATMID .and. atm_present) then !juanxiong he, for wrf/cam coupling, cam import
      call map_cam2cam_init_mct(cdata_cc, x2c_cc1, c2x_cc1, ATMID, &
                                cdata_cx, x2c_cx1, c2x_cx1, CPLID, CPLATMID)
      call map_cam2cam_init_mct(cdata_cc, x2c_cc2, c2x_cc2, ATMID, &
                                cdata_cx, x2c_cx2, c2x_cx2, CPLID, CPLATMID)
      call mct_avect_zero(x2c_cc1)
      call mct_avect_zero(x2c_cc2)
      call map_cama2camx_mct(cdata_cc, x2c_cc1, cdata_cx, x2c_cx1)  ! obviously, cam doesn't need data from wrf to drive its run
      call map_cama2camx_mct(cdata_cc, x2c_cc2, cdata_cx, x2c_cx2)  ! obviously, cam doesn't need data from wrf to drive its run
   endif

   if (iamin_CPLWRFID .and. wrf_present) then !juanxiong he, for wrf/cam coupling, wrf import
      call map_wrf2wrf_init_mct(cdata_mm, x2m_mm1, m2x_mm, WRFID, &
                                cdata_mx, x2m_mx1, m2x_mx, CPLID, CPLWRFID)
      call map_wrf2wrf_init_mct(cdata_mm, x2m_mm2, m2x_mm, WRFID, &
                                cdata_mx, x2m_mx2, m2x_mx, CPLID, CPLWRFID)
      call mct_avect_zero(x2m_mm1)
      call mct_avect_zero(x2m_mm2)
      call map_wrfw2wrfx_mct(cdata_mm, x2m_mm1, cdata_mx, x2m_mx1)   ! prepare import
      call map_wrfw2wrfx_mct(cdata_mm, x2m_mm2, cdata_mx, x2m_mx2)   ! prepare import
   endif

   if (iamin_CPLATMID .and. atm_present) then !juanxiong he, for gea/cam coupling, cam import
      call map_gcam2cam_init_mct(cdata_caca, x2ca_caca1, ca2x_caca1, ATMID, &
                                 cdata_cax, x2ca_cax1, ca2x_cax1, CPLID, CPLATMID)
      call map_gcam2cam_init_mct(cdata_caca, x2ca_caca2, ca2x_caca2, ATMID, &
                                 cdata_cax, x2ca_cax2, ca2x_cax2, CPLID, CPLATMID)
      call mct_avect_zero(x2ca_caca1)
      call mct_avect_zero(x2ca_caca2)
      call map_gcama2camx_mct(cdata_caca, x2ca_caca1, cdata_cax, x2ca_cax1)
      call map_gcama2camx_mct(cdata_caca, x2ca_caca2, cdata_cax, x2ca_cax2)
   endif

   if (iamin_CPLGEAID .and. geatm_present) then !juanxiong he, for geatm/cam coupling, geatm import
      call map_gea2gea_init_mct(cdata_gege, x2chem_chemchem1, chem2x_chemchem, GEAID, &
                                cdata_gex, x2chem_chemx1, chem2x_chemx, CPLID, CPLGEAID)
      call map_gea2gea_init_mct(cdata_gege, x2chem_chemchem2, chem2x_chemchem, GEAID, &
                                cdata_gex, x2chem_chemx2, chem2x_chemx, CPLID, CPLGEAID)
      call mct_avect_zero(x2chem_chemchem1)
      call mct_avect_zero(x2chem_chemchem2)
      call map_geaw2geax_mct(cdata_gege, x2chem_chemchem1, cdata_gex, x2chem_chemx1)   ! prepare import
      call map_geaw2geax_mct(cdata_gege, x2chem_chemchem2, cdata_gex, x2chem_chemx2)   ! prepare import
   endif
   !-----------------------------------------------------------------------------
   ! GEATM and WRF, Juanxiong He
   !-----------------------------------------------------------------------------
#endif
   call t_stopf  ('driver_init_xxx2xxx')

   !-----------------------------------------------------------------------------
   ! Remainder of initialization
   !-----------------------------------------------------------------------------

   if (iamin_CPLID) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

      !-----------------------------------------------------------------------------
      ! Allocate attribute vectors for merge components
      !-----------------------------------------------------------------------------

      if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing merge components'
!{list} key="init_remainder" number=1
!                       call mrg_x2{c}_init_mct( {allccc} )
!{list} key="init_remainder" number=2
!      if ({ccc}_present) call mrg_x2{c}_init_mct( {allccc} )
!wangty modify
#ifdef wrf
                       call mrg_x2w_init_mct( cdata_mx, c2x_mx1 )  !juanxiong he
                       call mrg_x2w_init_mct( cdata_mx, c2x_mx2 )  !juanxiong he
                       call mrg_x2c_init_mct( cdata_cx, m2x_cx )  !juanxiong he
                       call mrg_x2ge_init_mct( cdata_gex, ca2x_chemx1 )
!juanxiong he
                       call mrg_x2ge_init_mct( cdata_gex, ca2x_chemx2 )
!juanxiong he
                       call mrg_x2ca_init_mct( cdata_cax, chem2x_cax )
!juanxiong he
#endif

      !-----------------------------------------------------------------------------
      ! Initialize mapping
      ! Read aream into domains!
      !-----------------------------------------------------------------------------


      call t_startf('driver_init_maps')

!{list} key="init_map" number=1
!      if ({ccc2}_present) then
!         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing {ccc1}/{ccc2} mapping'
!         call map_{ccc1}2{ccc2}_init_mct(cdata_{c1}x, cdata_{c2}x)
!         call map_{ccc2}2{ccc1}_init_mct(cdata_{c2}x, cdata_{c1}x)
!      endif
!{list} key="init_map" number=2
!      if ({ccc1}_present .and. {ccc2}_present) then
!         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing {ccc1}/{ccc2} mapping'
!         call map_{ccc1}2{ccc2}_init_mct(cdata_{c1}x, cdata_{c2}x)
!         call map_{ccc2}2{ccc1}_init_mct(cdata_{c2}x, cdata_{c1}x)
!      endif
!{list} key="init_map" number=3
!      if ({ccc2}_present) then
!         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing {ccc1}/{ccc2} mapping'
!         call map_{ccc2}2{ccc1}_init_mct(cdata_{c2}x, cdata_{c1}x)
!      endif
!{list} key="init_map" number=4
!      if (rof_present .and. ocnrof_prognostic) then
!         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing rof/ocn mapping'
!         call map_rof2ocn_init_mct(cdata_rx, cdata_ox)
!      endif
!wangty modify
#ifdef wrf
      if (wrf_present) then     ! juanxiong he
         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing wrf/cam mapping'
         call map_cam2wrf_init_mct(cdata_cx, cdata_mx)
         call map_wrf2cam_init_mct(cdata_mx, cdata_cx)
      endif
      if (geatm_present) then     ! juanxiong he
         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing geatm/cam mapping'
         call map_cam2gea_init_mct(cdata_cax, cdata_gex)
         call map_gea2cam_init_mct(cdata_gex, cdata_cax)
      endif
#endif
      call t_stopf  ('driver_init_maps')

      !-----------------------------------------------------------------------------
      ! Check domains if appropriate
      ! This must be done after the mappers are initialized since
      ! checking is done on each processor and not with a global gather
      !-----------------------------------------------------------------------------

      if (domain_check) then
         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Performing domain checking'
!{list} key="domain_check"
!         call seq_domain_check_mct( cdata_{c}x, &                                    
!         )
      endif

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif ! iamin_CPLID

   !-----------------------------------------------------------------------------
   ! Map  dom_*x to dom_** in case any domain fields have been updated on cpl
   ! pes
   ! Initialize area corrections based on aream (read in map_init) and area
   ! Area correct component initialization output fields
   ! Map initial component AVs from component to coupler pes
   !-----------------------------------------------------------------------------
!{list} key="mapping"
!   if (iamin_CPL{ccc}ID .and. {ccc}_present) then
!      call map_{ccc}x2{ccc}{c1}_mct({allccc1})
!      if (iamin_{ccc}ID) then
!         call domain_areafactinit_mct(cdata_{cc},mdl2drv_{cc},drv2mdl_{cc},'areafact_{c}')
!         call mct_avect_vecmult({c}2x_{cc},mdl2drv_{cc},seq_flds_{c}2x_fluxes)
!      endif
!      call map_{ccc}{c1}2{ccc}x_mct({allccc2})
!   endif

   !-----------------------------------------------------------------------------
   ! global sum diagnostics for IC data
   !-----------------------------------------------------------------------------
   if (iamin_CPLID .and. info_debug > 1) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!{list} key="diagnostics_sum"
!      if ({ccc}_present) call seq_diag_avect_mct(cdata_{c}x,{c}2x_{c}x,'recv {ccc} IC')
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   end if

   !-----------------------------------------------------------------------------
   ! Initialize fractions
   !-----------------------------------------------------------------------------

   if (iamin_CPLID) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing fractions'
!{list} key="init_fractions"
!      call seq_frac_init(cdata_{c1}x, 
!                         {ccc1}_present, 
!                         dead_comps, &
!                         fractions_{c2}x
!                        )
      if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Setting fractions '
!{list} key="set_fractions"
!      call seq_frac_set({c1}2x_{c1}x, 
!                        cdata_{c2}x, 
!                        {ccc1}_present, 
!                        fractions_{c3}x
!                        )

      !-----------------------------------------------------------------------------
      ! Initialize atm/ocn flux component and compute ocean albedos
      !-----------------------------------------------------------------------------
!{list} key="flux_albedos"
!      if ({ccc2}_present) then
!         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Initializing {ccc1}/{ccc2} flux component'
!         call mct_aVect_init(x{c1}{c2}_{c2}x, rList=seq_flds_x{c1}{c2}_fields, lsize=mct_aVect_lsize({c2}2x_{c2}x))
!         call mct_aVect_zero(x{c1}{c2}_{c2}x)
!         if (trim({c1}{c2}flux_grid) == '{ccc2}') then
!            call seq_flux_init_mct(cdata_[c2]x,fractions_{c2}x)
!         elseif (trim({c1}{c2}flux_grid) == '{ccc1}') then
!            call seq_flux_init_mct(cdata_{c1}x,fractions_{c1}x)
!         elseif (trim({c1}{c2}flux_grid) == 'exch') then
!            call seq_flux_initexch_mct(cdata_{c1}x,cdata_{c2}x)
!         endif
!         call seq_flux_{ccc2}alb_mct(cdata_{c2}x,x{c1}{c2}_{c2}x,fractions_{c2}x)
!      endif

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif

   !-----------------------------------------------------------------------------
   ! Recalculate initial solar. Merge atmosphere input state and run atmospheric radiation
   ! tcx - for initialization only?
   !-----------------------------------------------------------------------------
!{list} key="compdriver_init"
!   call t_startf('driver_init_{ccc}init')

!   if ({ccc}_prognostic) then
      if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!{list} key="call_map" number="1"
!         if ({ccc2}_present) then
!            if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling map_{ccc2}2{ccc1}_mct'
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x, &
!                                  fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                  fluxlist=seq_flds_{c2}2x_fluxes, statelist=seq_flds_{c2}2x_states )
!         endif
!{list} key="call_map" number="2"
!         if ({ccc2}_present) then
!            if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling map_{ccc2}2{ccc1}_mct for mapping {c2}2x_{c2}x to {c2}2x_{c1}x'
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x, &
!                                  fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                  statelist=seq_flds_{c2}2x_states )
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x, &
!                                  fluxlist=seq_flds_{c2}2x_fluxes )
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, x{c1}{c2}_{c2}x, cdata_{c1}x, x{c1}{c2}_{c1}x, &
!                                  fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                  statelist=seq_flds_x{c1}{c2}_albedo )
!            if (trim(aoflux_grid) == '{ccc2}') then
!               if ( seq_comm_iamroot(CPLID)) &
!                  write(logunit,F00) 'Calling map_{ccc2}2{ccc1}_mct for mapping x{c1}{c2}_{c2}x to x{c1}{c2}_{c1}x'
!               call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, x{c1}{c2}_{c2}x, cdata_{c1}x, x{c1}{c2}_{c1}x, &
!                                     fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                     fluxlist=seq_flds_x{c1}{c2}_fluxes, statelist=seq_flds_x{c1}{c2}_states )
!            endif
!            if (trim(aoflux_grid) == '{ccc1}') then
!               if ( seq_comm_iamroot(CPLID)) &
!                  write(logunit,F00) 'Calling map_{ccc1}2{ccc2}_mct for mapping x{c1}{c2}_{c1}x to x{c1}{c2}_{c2}x'
!! tcraig: this mapping has to be done with area overlap mapping for all fields
!! due to the masking of the x{c1}{c2}_{c1}x data and the fact that states are mapped with
!! bilinear mapping currently
!!               call map_{ccc1}2{ccc2}_mct( cdata_{c1}x, x{c1}{c2}_{c1}x, cdata_{c2}x, x{c1}{c2}_{c2}x, &
!!                                     fluxlist=seq_flds_x{c1}{c2}_fluxes, statelist=seq_flds_x{c1}{c2}_states )
!               call map_{ccc1}2{ccc2}_mct( cdata_{c1}x, x{c1}{c2}_{c1}x, cdata_{c2}x, x{c1}{c2}_{c2}x, &
!                                     fluxlist=seq_flds_x{c1}{c2}_states//":"//seq_flds_x{c1}{c2}_fluxes)
!            endif
!         endif

!{list} key="call_mrg" 
!         if ({ccc2}_present .or. {ccc3}_present) then
!            if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling mrg_x2{c1}_run_mct'
!            call mrg_x2{c1}_run_mct( cdata_{c1}x, {c2}2x_{c1}x, {c3}2x_{c1}x, x{c1}{c3}_{c1}x, {c4}2x_{c1}x, fractions_{c1}x, x2{c1}_{c1}x )
!         endif

!{list} key="compdriver_init"
!         if (info_debug > 1) call seq_diag_avect_mct(cdata_{c}x,x2{c}_{c}x,'send {ccc} IC2')
!         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling {ccc}_init_mct'
!
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!      endif
!
!      if (iamin_CPL{ccc}ID) then
!         call map_{ccc}x2{ccc}{c}_mct( cdata_{c}x, x2{c}_{c}x, cdata_{c}{c}, x2{c}_{c}{c})
!         call seq_infodata_exchange(infodata,CPL{ccc}ID,'cpl2{ccc}_init')
!      endif
!   endif  

!{list} key="compdriver_init"
!   if (iamin_{ccc}ID) then
!      call t_adj_detailf(+2)
!      if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!      call seq_infodata_putData(infodata,{ccc}_phase=2)
!      call mct_avect_vecmult(x2{c}_{c}{c},drv2mdl_{c}{c},seq_flds_x2{c}_fluxes)
!!wangty modify
!#ifdef wrf
!      call atm_init_mct( EClock_a, cdata_aa, cdata_cc, cdata_caca, x2a_aa, a2x_aa, &
!                         x2c_cc1, x2c_cc2,  c2x_cc1, c2x_cc2, &
!                         x2ca_caca1, x2ca_caca2,  ca2x_caca1, ca2x_caca2, &
!                         twoway_coupling, twoway_nudging) !juanxiong he
!#else
!      call {ccc}_init_mct( EClock_{c}, cdata_{c}{c}, x2{c}_{c}{c}, {c}2x_{c}{c})
!#endif
!      call mct_avect_vecmult({c}2x_{c}{c},mdl2drv_{c}{c},seq_flds_{c}2x_fluxes)
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!      call t_adj_detailf(-2)
!   endif
!   if (iamin_CPL{ccc}ID) then
!      call map_{ccc}{c}2{ccc}x_mct( cdata_{c}{c}, {c}2x_{c}{c}, cdata_{c}x, {c}2x_{c}x)
!!wangty modify
!#ifdef wrf
!      call map_cama2camx_mct(cdata_cc, c2x_cc1, cdata_cx, c2x_cx1)  !juanxiong he
!      call map_cama2camx_mct(cdata_cc, c2x_cc2, cdata_cx, c2x_cx2)  !juanxiong he
!      call map_gcama2camx_mct(cdata_caca, ca2x_caca1, cdata_cax, ca2x_cax1) !juanxiong he
!      call map_gcama2camx_mct(cdata_caca, ca2x_caca2, cdata_cax, ca2x_cax2) !juanxiong he
!#endif
!      call seq_infodata_exchange(infodata,CPL{ccc}ID,'{ccc}2cpl_init')
!   endif
!   if (iamin_CPLID) then
!      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!      if (info_debug > 1) call seq_diag_avect_mct(cdata_{c}x,{c}2x_{c}x,'recv {ccc} IC2')
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!   endif

!wangty modify
#ifdef wrf
   !------------------------------------------------------------------
   ! wrf  and geatm initial, second phase, juanxiong he
   !------------------------------------------------------------------

   if (iamin_CPLID) then
!            call seq_hist_write(EClock_d)
            call map_cam2wrf_mct( cdata_cx, c2x_cx1, cdata_mx, c2x_mx1 , &
                                    statelist=seq_flds_c2x_states ) ! juanxiong he
            call map_cam2wrf_mct( cdata_cx, c2x_cx2, cdata_mx, c2x_mx2 , &
                                    statelist=seq_flds_c2x_states ) ! juanxiong he
            call mrg_x2w_run_mct( cdata_mx, c2x_mx1, x2m_mx1)  ! juanxiong he
            call mrg_x2w_run_mct( cdata_mx, c2x_mx2, x2m_mx2)  ! juanxiong he
            call map_cam2gea_mct( cdata_cax, ca2x_cax1, cdata_gex, ca2x_chemx1 , &
                                  statelist=seq_flds_ca2x_states ) ! juanxiong he
            call map_cam2gea_mct( cdata_cax, ca2x_cax2, cdata_gex, ca2x_chemx2 , &
                                  statelist=seq_flds_ca2x_states ) ! juanxiong he
            call mrg_x2ge_run_mct( cdata_gex, ca2x_chemx1, x2chem_chemx1)  ! juanxiong he
            call mrg_x2ge_run_mct( cdata_gex, ca2x_chemx2, x2chem_chemx2)  ! juanxiong he
   endif

   if (iamin_CPLWRFID .and. wrf_prognostic) then
            call t_drvstartf ('driver_c2w_wrfx2wrfa',barrier=mpicom_CPLWRFID)
            call map_wrfx2wrfw_mct( cdata_mx, x2m_mx1, cdata_mm, x2m_mm1)   ! juaxniong he
            call map_wrfx2wrfw_mct( cdata_mx, x2m_mx2, cdata_mm, x2m_mm2)   ! juaxniong he
            call t_drvstopf  ('driver_c2w_wrfx2wrfa')
   endif

   if (iamin_WRFID) then
      call t_adj_detailf(+2)
      if (drv_threading) call seq_comm_setnthreads(nthreads_WRFID)
      call seq_infodata_putData(infodata,wrf_phase=2)
      call wrf_init_mct( EClock_w, cdata_ww, cdata_mm, x2w_ww, w2x_ww, x2m_mm1, x2m_mm2, m2x_mm, &
                         twoway_coupling, twoway_nudging )    ! juanxiong he
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      twoway_coupling = .False.
      call t_adj_detailf(-2)
   endif

   if (iamin_CPLGEAID .and. geatm_prognostic) then
            call t_drvstartf ('driver_ca2ge_geax2geaa',barrier=mpicom_CPLGEAID)
            call map_geax2geaw_mct( cdata_gex, x2chem_chemx1, cdata_gege, x2chem_chemchem1)   !juaxniong he
            call map_geax2geaw_mct( cdata_gex, x2chem_chemx2, cdata_gege, x2chem_chemchem2)   !juaxniong he
            call t_drvstopf  ('driver_ca2ge_geax2geaa')
   endif

   if (iamin_GEAID) then
      call t_adj_detailf(+2)
      if (drv_threading) call seq_comm_setnthreads(nthreads_GEAID)
      call seq_infodata_putData(infodata,geatm_phase=2)
      call geatm_init_mct( EClock_ge, cdata_gege, &
                           x2chem_chemchem1, x2chem_chemchem2, chem2x_chemchem ) ! juanxiong he
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      call t_adj_detailf(-2)
   endif
#endif
!{list} key="compdriver_init"
!   call t_stopf  ('driver_init_{ccc}init')

   !-----------------------------------------------------------------------------
   ! Read driver restart file, overwrite anything previously sent or computed
   !-----------------------------------------------------------------------------

   call t_startf('driver_init_readrestart')
   call seq_diag_zero_mct(mode='all')
   if (read_restart) call seq_rest_read(rest_file)
   call t_stopf  ('driver_init_readrestart')

   if (do_histinit) then
      if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod )
            write(logunit,104) ' Write history file at ',ymd,tod
            call shr_sys_flush(logunit)
         endif
         call seq_hist_write(EClock_d)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   endif

   if ( iamroot_CPLID )then
      write(logunit,*) ' '
      write(logunit,F00) 'Model initialization complete '
      write(logunit,*) ' '
      call shr_sys_flush(logunit)
   endif

   call t_stopf  ('DRIVER_INIT')

end subroutine ccsm_init

!===============================================================================

subroutine ccsm_run()

 101  format( A, 2i8, 12A, A, F8.2, A, F8.2 )
 102  format( A, 2i8, A, 5L3 )
 103  format( 5A )
 104  format( A, 2i8)
 105  format( A, 2i8, A, f10.2, A, f10.2, A, A, i5, A, A)
 106  format( A, f23.12)

!{list} key="putdata"
!   call seq_infodata_putData(infodata,{ccc}_phase=1)
   call seq_timemgr_EClockGetData( EClock_d, stepno=begstep)
   call seq_timemgr_EClockGetData( EClock_d, dtime=dtime)
   ncpl = 86400/dtime
   dtstep_acc = 0._r8
   dtstep_cnt = 0
   stop_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
   if (seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_datestop)) then
      if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,103) subname,' NOTE: Stopping from alarm STOP DATE'
          write(logunit,*) ' '
      endif
      stop_alarm = .true.
   endif

   !----------------------------------------------------------
   ! Beginning of basic time step loop
   !----------------------------------------------------------

   call t_startf ('DRIVER_RUN_LOOP_BSTART')
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf ('DRIVER_RUN_LOOP_BSTART')
   Time_begin = mpi_wtime()
   Time_bstep = mpi_wtime()
   do while ( .not. stop_alarm)

      call t_startf('DRIVER_RUN_LOOP')
      call t_drvstartf ('DRIVER_CLOCK_ADVANCE',cplrun=.true.)

      !----------------------------------------------------------
      ! Advance sync clock time (this is time that models should have before
      ! they return to the driver).  Write timestamp and run alarm status
      !----------------------------------------------------------

      call seq_timemgr_clockAdvance( seq_SyncClock)
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod )
      call shr_cal_date2ymd(ymd,year,month,day)
      stop_alarm    = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
!wangty modify
#ifdef wrf
      wrfrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_wrfrun) !juanxiong he
      gearun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_gearun) !juanxiong he
#endif
!{list} key="compoents_alarm"
!      {ccc}run_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_{ccc}run)
!{list} key="compoents_nextalarm"
!      {ccc}next_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_{ccc}next)
      restart_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_restart)
      history_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_history)
      histavg_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_histavg)
      tprof_alarm   = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_tprof)

      ! this probably belongs in seq_timemgr somewhere using proper clocks
      t1hr_alarm = .false.
      t2hr_alarm = .false.
      t3hr_alarm = .false.
      t6hr_alarm = .false.
      t12hr_alarm = .false.
      t24hr_alarm = .false.
      if (mod(tod, 3600) == 0) t1hr_alarm = .true.
      if (mod(tod, 7200) == 0) t2hr_alarm = .true.
      if (mod(tod,10800) == 0) t3hr_alarm = .true.
      if (mod(tod,21600) == 0) t6hr_alarm = .true.
      if (mod(tod,43200) == 0) t12hr_alarm = .true.
      if (tod == 0) t24hr_alarm = .true.

!{list} key="alarm_putdata"
!      call seq_infodata_putData(infodata, {ccc}run_alarm={ccc}run_alarm)

      if (seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_datestop)) then
         if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,103) subname,' NOTE: Stopping from alarm STOP DATE'
             write(logunit,*) ' '
         endif
         stop_alarm = .true.
      endif

      ! update the orbital data as needed
      if (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
         orb_nyear =  orb_iyear + (year - orb_iyear_align)
         if (orb_nyear /= orb_cyear) then
            orb_cyear = orb_nyear
            call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                                orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
            call seq_infodata_putData(infodata,orb_eccen=orb_eccen,orb_obliqr=orb_obliqr, &
                 orb_lambm0=orb_lambm0,orb_mvelpp=orb_mvelpp)
         endif
      endif

!{list} key="compoents_nextalarm"
!      ! override {ccc}run_alarm and {ccc}next_alarm for first {ccc} run
!      ! skip_{ccccc}_run is initialized above to true if it's a startup
!      ! if it's not a startup, ignore all of this
!      ! stop the overide on the second {ccc}run_alarm
!
!      if ({ccc}run_alarm) {ccc}run_count = {ccc}run_count + 1
!      if ({ccc}run_count > 1) skip_{ccccc}_run = .false.
!      if (skip_{ccccc}_run) then
!         {ccc}run_alarm = .false.
!         {ccc}next_alarm = .false.
!      endif

!{list} key="alarm_state"
!      if (iamroot_CPLID) then
!         if (loglevel > 1) then
!!wangty modify
!#ifdef wrf
!            write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
!               ' aliog run alarms = ',  atmrun_alarm, lndrun_alarm, &
!                         icerun_alarm, ocnrun_alarm, glcrun_alarm, wrfrun_alarm, gearun_alarm !juanxiong he
!#else
!            write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
!               ' aliog run alarms = ',  {ccc}run_alarm
!#endif
!            write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
!               ' 1.2.3.6.12 run alarms = ',  t1hr_alarm, t2hr_alarm, &
!                         t3hr_alarm, t6hr_alarm, t12hr_alarm, t24hr_alarm
!            call shr_sys_flush(logunit)
!         endif
!      endif

      call t_drvstopf  ('DRIVER_CLOCK_ADVANCE',cplrun=.true.)

!{list} key="model_prep"
!      !----------------------------------------------------------
!      ! {ccc1}/{ccc2} PREP
!      ! Map for {ccc2} prep and {ccc3}{ccc1} flux
!      !----------------------------------------------------------
!
!      if (iamin_CPLID .and. ({ccc2}_present.or.{ccc1}_present)) then
!         call t_drvstartf ('DRIVER_{ccc1}PREP',cplrun=.true.,barrier=mpicom_CPLID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!         call t_drvstartf ('driver_{ccc1}prep_{ccc3}2{ccc1}',barrier=mpicom_CPLID)
!         call map_{ccc3}2{ccc1}_mct( cdata_{c3}x, {c3}2x_{c3}x, cdata_{c1}x, {c3}2x_{c1}x, &
!                               fluxlist=seq_flds_{c3}2x_fluxes, statelist=seq_flds_{c3}2x_states )
!         call t_drvstopf  ('driver_{ccc1}prep_{ccc3}2{ccc1}')
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_{ccc1}PREP',cplrun=.true.)
!      endif

!{list} key="model_setup" number="1"
!      !----------------------------------------------------------
!      ! {ccc1} SETUP
!      !----------------------------------------------------------
!
!      if ({ccc1}_present .and. {ccc1}run_alarm) then
!
!         !----------------------------------------------------
!         ! "startup" wait
!         !----------------------------------------------------
!
!         if (iamin_CPL{ccc1}ID .and. cpl2{ccc1}_first) then
!            ! want to know the time the ocean pes waited for the cpl pes
!            !   at the first ocnrun_alarm, min ocean wait is wait time
!            ! do not use t_barrierf here since it can be "off", use mpi_barrier
!            if (iamin_{ccc1}ID) call t_drvstartf ('DRIVER_C2O_INITWAIT')
!            call mpi_barrier(mpicom_CPL{ccc1}ID,ierr)
!            if (iamin_{ccc1}ID) call t_drvstopf  ('DRIVER_C2O_INITWAIT')
!            cpl2{ccc1}_first = .false.
!         endif
!
!         !----------------------------------------------------
!         ! {ccc1} prep
!         !----------------------------------------------------
!
!         if (iamin_CPLID .and. {ccc1}_prognostic) then
!            call t_drvstartf ('DRIVER_{ccc1}PREP',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            call t_drvstartf ('driver_{ccc1}prep_avg',barrier=mpicom_CPLID)
!            ! temporary formation of average
!!           call mct_accum_average(x2{c1}acc_{c1}x)
!            if (x2{c1}acc_{c1}x_cnt > 0) then
!               x2{c1}acc_{c1}x%data%rAttr = x2{c1}acc_{c1}x%data%rAttr / (x2{c1}acc_{c1}x_cnt*1.0_r8)
!            endif
!            x2{c1}acc_{c1}x_cnt = 0
!            call t_drvstopf  ('driver_{ccc1}prep_avg')
!
!            if ({ccc2}_present .and. {ccc1}rof_prognostic) then
!               if ({c2}2xacc_{c2}x_cnt > 0) then
!                  call t_drvstartf ('driver_{ccc1}prep_ravg',barrier=mpicom_CPLID)
!                  {c2}2xacc_{c2}x%data%rAttr = {c2}2xacc_{c2}x%data%rAttr / ({c2}2xacc_{c2}x_cnt*1.0_r8)
!                  {c2}2xacc_{c2}x_cnt = 0
!                  call t_drvstopf ('driver_{ccc1}prep_ravg')
!                  call t_drvstartf ('driver_{ccc1}prep_{ccc2}2{ccc1}',barrier=mpicom_CPLID)
!                  call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2xacc_{c2}x%data, cdata_{c1}x, {c2}2x_{c1}x )
!                  if (do_hist_{c2}2x) then
!                     call seq_hist_writeaux(EClock_d,'{c2}2xacc','dom{c2}',cdata_{c2}x,{c2}2xacc_{c2}x%data, &
!                          {ccc2}_nx,{ccc2}_ny,1)
!                  endif
!                  call t_drvstopf  ('driver_{ccc1}prep_{ccc2}2{ccc1}')
!                  call t_drvstartf ('driver_{ccc1}prep_{ccc2}copy',barrier=mpicom_CPLID)
!                  call mct_aVect_copy(aVin={c2}2x_{c1}x, aVout=x2{c1}acc_{c1}x%data)
!                  call t_drvstopf  ('driver_{ccc1}prep_{ccc2}copy')
!               endif
!            endif
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc1}prep_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c1}x,x2{c1}acc_{c1}x%data,'send {ccc1}')
!               call t_drvstopf  ('driver_{ccc1}prep_diagav')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc1}PREP',cplrun=.true.)
!         endif
!
!
!         !----------------------------------------------------
!         ! cpl -> {ccc1}
!         !----------------------------------------------------
!
!         if (iamin_CPL{ccc1}ID .and. {ccc1}_prognostic) then
!            call t_drvstartf ('DRIVER_C2{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call t_drvstartf ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}',barrier=mpicom_CPLOCNID)
!            call map_{ccc1}x2{ccc1}{c1}_mct( cdata_{c1}x, x2{c1}acc_{c1}x%data, cdata_{c1}{c1}, x2{c1}_{c1}{c1})
!            call t_drvstopf  ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}')
!            call t_drvstartf ('driver_c2{c1}_infoexch',barrier=mpicom_CPL{ccc1}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc1}ID,'cpl2{ccc1}_run')
!            call t_drvstopf  ('driver_c2{c1}_infoexch')
!            call t_drvstopf  ('DRIVER_C2{c1}')
!         endif
!
!      endif

!{list} key="model_setup" number="2"
!      !----------------------------------------------------------
!      ! {ccc1} SETUP
!      !----------------------------------------------------------
!
!      if ({ccc1}_present .and. {ccc1}run_alarm) then
!
!         !----------------------------------------------------
!         ! {ccc1} prep
!         !----------------------------------------------------
!
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER_{ccc1}PREP',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!
!            if ({ccc1}_prognostic) then
!               call t_drvstartf ('driver_{ccc1}prep_{ccc2}2{ccc1}',barrier=mpicom_CPLID)
!               call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x )
!               call t_drvstopf  ('driver_{ccc1}prep_{ccc2}2{ccc1}')
!               call t_drvstartf ('driver_{ccc1}prep_mrgx2{c1}',barrier=mpicom_CPLID)
!               call mrg_x2{c1}_run_mct( cdata_{c1}x, {c2}2x_{c1}x, x2{c1}_{c1}x )
!               call t_drvstopf  ('driver_{ccc1}prep_mrgx2{c1}')
!               if (info_debug > 1) then
!                  call t_drvstartf ('driver_{ccc1}prep_diagav',barrier=mpicom_CPLID)
!                  call seq_diag_avect_mct(cdata_{c1}x,x2{c1}_{c1}x,'send {ccc1}')
!                  call t_drvstopf  ('driver_{ccc1}prep_diagav')
!               endif
!            endif
!
!            if ({ccc3}_present .and. {ccc4}_prognostic) then
!               if (info_debug > 1) then
!                  call t_drvstartf ('driver_{ccc1}prep_diagav',barrier=mpicom_CPLID)
!                  call seq_diag_avect_mct(cdata_{c4}x,x2{c3}_{c4}x,'send {ccc4}')
!                  call t_drvstopf  ('driver_{ccc1}prep_diagav')
!               endif
!            endif
!
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc1}PREP',cplrun=.true.)
!         endif
!
!         !----------------------------------------------------
!         ! cpl -> {ccc1}
!         !----------------------------------------------------
!
!         if (iamin_CPL{ccc1}ID) then
!            call t_drvstartf ('DRIVER_C2{c1}',barrier=mpicom_CPL{ccc1}ID)
!            if ({ccc1}_prognostic) then
!               call t_drvstartf ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}',barrier=mpicom_CPL{ccc1}ID)
!               call map_{ccc1}x2{ccc1}{c1}_mct( cdata_{c1}x, x2{c1}_{c1}x, cdata_{c1}{c1}, x2{c1}_{c1}{c1})
!               call t_drvstopf  ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}')
!            endif
!            if ({ccc3}_present .and. {ccc4}_prognostic) then
!               call t_drvstartf ('driver_c2{c1}_{ccc4}x2{ccc4}{c4}',barrier=mpicom_CPL{ccc1}ID)
!               call map_{ccc4}x2{ccc4}{c4}_mct( cdata_{c4}x, x2{c4}_{c4}x, cdata_{c4}{c4}, x2{c4}_{c4}{c4})
!               call t_drvstopf  ('driver_c2{c1}_{ccc4}x2{ccc4}{c4}')
!            endif
!            if ({ccc1}_prognostic .or. {ccc4}_prognostic) then
!               call t_drvstartf ('driver_c2{c1}_infoexch',barrier=mpicom_CPL{ccc1}ID)
!               call seq_infodata_exchange(infodata,CPL{ccc1}ID,'cpl2{ccc1}_run')
!               call t_drvstopf  ('driver_c2{c1}_infoexch')
!            endif
!            call t_drvstopf  ('DRIVER_C2{c1}')
!         endif
!
!      endif

!{list} key="model_setup" number="3"
!      !----------------------------------------------------------
!      ! {ccc1} SETUP
!      !----------------------------------------------------------
!
!      if ({ccc1}_present .and. {ccc1}run_alarm) then
!
!         !----------------------------------------------------
!         ! {ccc1} prep
!         !----------------------------------------------------
!
!         if (iamin_CPLID .and. {ccc1}_prognostic) then
!            call t_drvstartf ('DRIVER_{ccc1}PREP',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            call t_drvstartf ('driver_{ccc1}prep_{ccc2}2{ccc1}',barrier=mpicom_CPLID)
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x )
!            call t_drvstopf  ('driver_{ccc1}prep_{ccc2}2{ccc1}')
!
!            call t_drvstartf ('driver_{ccc1}prep_{ccc3}2{ccc1}',barrier=mpicom_CPLID)
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c3}2x_{c2}x, cdata_{c1}x, {c3}2x_{c1}x )
!!tcx fails            call map_{ccc3}2{ccc1}_mct( cdata_{c3}x, {c3}2x_{c3}x, cdata_{c1}x, {c3}2x_{c1}x )
!            call t_drvstopf  ('driver_{ccc1}prep_{ccc3}2{ccc1}')
!
!            call t_drvstartf ('driver_{ccc1}prep_mrgx2{c1}',barrier=mpicom_CPLID)
!            call mrg_x2{c1}_run_mct( cdata_{c1}x, {c3}2x_{c1}x, {c2}2x_{c1}x, x2{c1}_{c1}x )
!            call t_drvstopf  ('driver_{ccc1}prep_mrgx2{c1}')
!
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc1}prep_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c1}x,x2{c1}_{c1}x,'send {ccc1}')
!               call t_drvstopf  ('driver_{ccc1}prep_diagav')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc1}PREP',cplrun=.true.)
!         endif
!
!         !----------------------------------------------------
!         ! cpl -> {ccc1}
!         !----------------------------------------------------
!
!         if (iamin_CPLICEID .and. {ccc1}_prognostic) then
!            call t_drvstartf ('DRIVER_C2{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call t_drvstartf ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call map_{ccc1}x2{ccc1}{c1}_mct( cdata_{c1}x, x2{c1}_{c1}x, cdata_{c1}{c1}, x2{c1}_{c1}{c1})
!            call t_drvstopf  ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}')
!            call t_drvstartf ('driver_c2{c1}_infoexch',barrier=mpicom_CPL{ccc1}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc1}ID,'cpl2{ccc1}_run')
!            call t_drvstopf  ('driver_c2{c1}_infoexch')
!            call t_drvstopf  ('DRIVER_C2{c1}')
!         endif
!
!      endif

!{list} key="model_run" number="1"
!      !----------------------------------------------------------
!      ! Run {ccc} Model
!      !----------------------------------------------------------
!
!      if ({ccc}_present .and. {ccc}run_alarm .and. iamin_{ccc}ID) then
!         if (run_barriers) then
!            call t_drvstartf ('DRIVER_{ccc}_RUN_BARRIER')
!            call mpi_barrier(mpicom_{ccc}ID,ierr)
!            call t_drvstopf ('DRIVER_{ccc}_RUN_BARRIER')
!         endif
!         call t_drvstartf ('DRIVER_{ccc}_RUN',barrier=mpicom_{ccc}ID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!         if ({ccc}_prognostic) call mct_avect_vecmult(x2{c}_{c}{c},drv2mdl_{c}{c},seq_flds_x2{c}_fluxes)
!         call ocn_run_mct( EClock_{c}, cdata_{c}{c}, x2{c}_{c}{c}, {c}2x_{c}{c})
!         call mct_avect_vecmult({c}2x_{c}{c},mdl2drv_{c}{c},seq_flds_{c}2x_fluxes)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_{ccc}_RUN')
!      endif

!{list} key="model_run" number="2"
!      !----------------------------------------------------------
!      ! Run {ccc1} Model
!      !----------------------------------------------------------
!
!      if (({ccc1}_present.or.{ccc2}_present.or.{ccc3}_present) .and. &
!           {ccc1}run_alarm .and. iamin_{ccc1}ID) then
!         if (run_barriers) then
!            call t_drvstartf ('DRIVER_{ccc1}_RUN_BARRIER')
!            call mpi_barrier(mpicom_{ccc1}ID,ierr)
!            call t_drvstopf ('DRIVER_{ccc1}_RUN_BARRIER')
!         endif
!         call t_drvstartf ('DRIVER_{ccc1}_RUN',barrier=mpicom_{ccc1}ID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc1}ID)
!         if ({ccc1}_prognostic) then
!            call mct_avect_vecmult(x2{c1}_{c1}{c1},drv2mdl_{c1}{c1},seq_flds_x2{c1}_fluxes)
!         endif
!         if ({ccc3}_prognostic) then
!            call mct_avect_vecmult(x2{c3}_{c3}{c3},drv2mdl_{c3}{c3},seq_flds_x2{c3}_fluxes)
!         endif
!         call {ccc1}_run_mct( EClock_{c1}, cdata_{c1}{c1}, x2{c1}_{c1}{c1}, {c1}2x_{c1}{c1}, &
!                                     cdata_{c2}{c2},         {c2}2x_{c2}{c2}, &
!                                     cdata_{c3}{c3}, x2{c3}_{c3}{c3}, {c3}2x_{c3}{c3})
!         call mct_avect_vecmult({c1}2x_{c1}{c1},mdl2drv_{c1}{c1},seq_flds_{c1}2x_fluxes)
!         if ({ccc2}_present) then
!            call mct_avect_vecmult({c2}2x_{c2}{c2},mdl2drv_{c2}{c2},seq_flds_{c2}2x_fluxes)
!         endif
!         if ({ccc3}_present) then
!            call mct_avect_vecmult({c3}2x_{c3}{c3},mdl2drv_{c3}{c3},seq_flds_{c3}2x_fluxes)
!         endif
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_{ccc1}_RUN')
!      endif

!{list} key="ocn_def"
!      !----------------------------------------------------------
!      ! {ccc} -> cpl, tight coupling (sequential type mode)
!      !----------------------------------------------------------
!
!      if ({ccccc}_tight_coupling) then
!      if ({ccc}_present .and. {ccc}next_alarm) then
!         if (iamin_CPL{ccc}ID) then
!            call t_drvstartf ('DRIVER_{c}2C',barrier=mpicom_CPL{ccc}ID)
!            call t_drvstartf ('driver_{c}2c_{ccc}{c}2{ccc}x',barrier=mpicom_CPL{ccc}ID)
!            call map_{ccc}{c}2{ccc}x_mct( cdata_{c}{c}, {c}2x_{c}{c}, cdata_{c}x, {c}2x_{c}x)
!            call t_drvstopf  ('driver_{c}2c_{ccc}{c}2{ccc}x')
!            call t_drvstartf ('driver_{c}2c_infoexch',barrier=mpicom_CPL{ccc}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc}ID,'{ccc}2cpl_run')
!            call t_drvstopf  ('driver_{c}2c_infoexch')
!            call t_drvstopf  ('DRIVER_{c}2C')
!         endif
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER_{ccc}POST',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc}post_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c}x,{c}2x_{c}x,'recv {ccc}')
!               call t_drvstopf  ('driver_{ccc}post_diagav')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc}POST',cplrun=.true.)
!         endif
!      endif
!      endif

!{list} key="model_cpl" number="1"
!      if ({ccc1}_present .and. iamin_CPLID) then
!         call t_drvstartf ('DRIVER_{ccc2}{ccc1}P',cplrun=.true.,barrier=mpicom_CPLID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!         if ({ccc1}_prognostic) then
!
!            ! Map {ccc3} to {ccc1}
!            if ({ccc3}_present) then
!               call t_drvstartf ('driver_{ccc2}{ccc1}p_{ccc3}2{ccc1}',barrier=mpicom_CPLID)
!               call map_{ccc3}2{ccc1}_mct( cdata_{c3}x, {c3}2x_{c3}x, cdata_{c1}x, {c3}2x_{c1}x)
!               call t_drvstopf  ('driver_{ccc2}{ccc1}p_{ccc3}2{ccc1}')
!            endif
!
!            ! Merge {ccc1} inputs
!            call t_drvstartf ('driver_{ccc2}{ccc1}p_mrgx2{c1}',barrier=mpicom_CPLID)
!            call mrg_x2{c1}_run_mct( cdata_{c1}x, {c2}2x_{c1}x, {c3}2x_{c1}x, x{c2}{c1}_{c1}x, fractions_{c1}x, x2{c1}_{c1}x )
!            call t_drvstopf  ('driver_{ccc2}{ccc1}p_mrgx2{c1}')
!
!            ! Accumulate {ccc1} inputs
!            ! Form partial sum of tavg ocn inputs (virtual "send" to {ccc1})
!            call t_drvstartf ('driver_{ccc2}{ccc1}p_accum',barrier=mpicom_CPLID)
!!     !         call mct_accum_accumulate(x2o_ox, x2oacc_ox)
!            if (x2{c1}acc_{c1}x_cnt == 0) then
!               x2{c1}acc_{c1}x%data%rAttr = x2{c1}_{c1}x%rAttr
!            else
!               x2{c1}acc_{c1}x%data%rAttr = x2{c1}acc_{c1}x%data%rAttr + x2{c1}_{c1}x%rAttr
!            endif
!            x2{c1}acc_{c1}x_cnt = x2{c1}acc_{c1}x_cnt + 1
!            call t_drvstopf  ('driver_{ccc2}{ccc1}p_accum')
!         endif
!
!         ! Compute {ccc2}/{ccc1} fluxes (virtual "recv" from {ccc1})
!         call t_drvstartf ('driver_{ccc2}{ccc1}p_flux',barrier=mpicom_CPLID)
!         if (trim({c2}{c1}flux_grid) == '{ccc1}') then
!            call seq_flux_{ccc2}{ccc1}_mct( cdata_{c1}x, {c2}2x_{c1}x, {c1}2x_{c1}x, x{c2{c1}_{c1}x)
!            call seq_flux_{ccc1}alb_mct(cdata_{c1}x,x{c2}{c1}_{c1}x,fractions_{c1}x)
!         else if (trim({c2}{c1}flux_grid) == '{ccc2}') then
!            call seq_flux_{ccc2}{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, {c1}2x_{c2}x, x{c2}[c1]_{c2}x)
!            call seq_flux_{ccc1}alb_mct(cdata_{c1}x,x{c2}{c1}_{c1}x,fractions_{c1}x)
!         else if (trim({c2}{c1}flux_grid) == 'exch') then
!            call seq_flux_{ccc2}{ccc1}exch_mct( cdata_{c2}x, cdata_{c1}x, {c2}2x_{c2}x, {c1}2x_{c1}x, x{c2}{c1}_{c2}x, x{c2}{c1}_{c1}x, &
!                                       fractions_{c2}x, fractions_{c1}x)
!            call seq_flux_{ccc1}alb_mct(cdata_{c1}x,x{c2}{c1}_{c1}x,fractions_{c1}x)
!         endif  ! {c2}{c1}flux_grid
!         call t_drvstopf  ('driver_{ccc2}{ccc1}p_flux')
!
!         if (trim({c2}{c1}flux_grid) == '{ccc2}') then
!            call t_drvstartf ('driver_{ccc2}{ccc1}p_{ccc2}2{ccc1}',barrier=mpicom_CPLID)
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, x{c2}{c1}_{c2}x, cdata_{c1}x, x{c2}{c1}_{c1}x, &
!                                   fluxlist=seq_flds_x{c2}{c1}_states//":"//seq_flds_x{c2}{c1}_fluxes)
!            call t_drvstopf  ('driver_{ccc2}{ccc1}p_{ccc2}2{ccc1}')
!         endif
!
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_{ccc2}{ccc1}P',cplrun=.true.)
!      endif

!{list} key="model_cpl" number="2"
!      !----------------------------------------------------------
!      ! {ccc1} -> cpl
!      !----------------------------------------------------------
!
!      if (({ccc1}_present.or.{ccc2}_present.or.{ccc3}_present) .and. {ccc1}run_alarm) then
!
!         if (iamin_CPL{ccc1}ID) then
!            call t_drvstartf ('DRIVER_{c1}2C',barrier=mpicom_CPL{ccc1}ID)
!            if ({ccc1}_present) then
!               call t_drvstartf ('driver_{c1}2c_{ccc1}{c1}2{ccc1}x',barrier=mpicom_CPL{ccc1}ID)
!               call map_{ccc1}{c1}2{ccc1}x_mct( cdata_{c1}{c1}, {c1}2x_{c1}{c1}, cdata_{c1}x, {c1}2x_{c1}x)
!               call t_drvstopf  ('driver_{c1}2c_{ccc1}{c1}2{ccc1}x')
!            endif
!            if ({ccc2}_present) then
!               call t_drvstartf ('driver_{c1}2c_{ccc2}{c2}2{ccc2}x',barrier=mpicom_CPL{ccc1}ID)
!               call map_{ccc2}{c2}2{ccc2}x_mct( cdata_{c2}{c2}, {c2}2x_{c2}{c2}, cdata_{c2}x, {c2}2x_{c2}x)
!               call t_drvstopf  ('driver_{c1}2c_{ccc2}{c2}2{ccc2}x')
!            endif
!            if ({ccc3}_present .and. {ccc4}_prognostic .and. {ccc4}run_alarm) then
!               call t_drvstartf ('driver_{c1}2c_{ccc3}{c3}2{ccc3}x',barrier=mpicom_CPL{ccc1}ID)
!               call map_{ccc3}{c3}2{ccc3}x_mct( cdata_{c3}{c3}, {c3}2x_{c3}{c3}, cdata_{c3}x, {c3}2x_{c3}x)
!               call t_drvstopf  ('driver_{c1}2c_{ccc3}{c3}2{ccc3}x')
!            endif
!            call t_drvstartf ('driver_{c1}2c_infoexch',barrier=mpicom_CPL{ccc1}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc1}ID,'{ccc1}2cpl_run')
!            call t_drvstopf  ('driver_{c1}2c_infoexch')
!            call t_drvstopf  ('DRIVER_{c1}2C')
!         endif
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER_{ccc1}POST',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc1}post_diagav',barrier=mpicom_CPLID)
!               if ({ccc1}_present) then
!                  call seq_diag_avect_mct(cdata_{c1}x,{c1}2x_{c1}x,'recv {ccc1}')
!               endif
!               if ({ccc2}_present) then
!                  call seq_diag_avect_mct(cdata_{c2}x,{c2}2x_{c2}x,'recv {ccc2}')
!               endif
!               if ({ccc3}_present .and. {ccc4}_prognostic .and. {ccc4}run_alarm) then
!                  call seq_diag_avect_mct(cdata_{c4}x,{c4}2x_{c4}x,'recv {ccc4}')
!               endif
!               call t_drvstopf  ('driver_{ccc1}post_diagav')
!            endif
!            if ({ccc2}_present .and. {ccc5}{ccc2}_prognostic) then
!               call t_drvstartf ('driver_{ccc1}post_raccum',barrier=mpicom_CPLID)
!               ! better to flux correct here if flux_epbalfact varies
!               ! over the accumulation period
!               call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)
!               if ({c2}2xacc_{c2}x_cnt == 0) then
!                  {c2}2xacc_{c2}x%data%rAttr = {c2}2x_{c2}x%rAttr * flux_epbalfact
!               else
!                  {c2}2xacc_{c2}x%data%rAttr = {c2}2xacc_{c2}x%data%rAttr + {c2}2x_{c2}x%rAttr * flux_epbalfact
!               endif
!               {c2}2xacc_{c2}x_cnt = {c2}2xacc_{c2}x_cnt + 1
!               call t_drvstopf ('driver_{ccc1}post_raccum')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc1}POST',cplrun=.true.)
!         endif
!
!      endif   ! run alarm, lnd_present

!{list} key="model_setup" number="4"
!      !----------------------------------------------------------
!      ! {ccc1} SETUP
!      !----------------------------------------------------------
!
!      if ({ccc2}_present .and. {ccc1}run_alarm) then
!         if (iamin_CPLID .and. {ccc1}_prognostic) then
!            call t_drvstartf ('DRIVER_{ccc1}PREP',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!
!            ! Map {ccc2} to {ccc1}
!            call t_drvstartf ('driver_{ccc1}prep_{ccc2}2{ccc1}',barrier=mpicom_CPLID)
!            call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x )
!            call t_drvstopf  ('driver_{ccc1}prep_{ccc2}2{ccc1}')
!
!            ! Merge {ccc1} inputs
!            call t_drvstartf ('driver_{ccc1}prep_mrgx2{c1}',barrier=mpicom_CPLID)
!            call mrg_x2{c1}_run_mct( cdata_{c1}x, {c2}2x_{c1}x, x2{c1}_{c1}x)
!            call t_drvstopf  ('driver_{ccc1}prep_mrgx2{c1}')
!
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc1}prep_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c1}x,x2{c1}_{c1}x,'send {ccc1}')
!               call t_drvstopf  ('driver_{ccc1}prep_diagav')
!            endif
!
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc}PREP',cplrun=.true.)
!         endif
!
!         !----------------------------------------------------
!         ! cpl -> {ccc}
!         !----------------------------------------------------
!
!         if (iamin_CPL{ccc}ID .and. {ccc}_prognostic) then
!            call t_drvstartf ('DRIVER_C2{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call t_drvstartf ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call map_{ccc1}x2{ccc1}{c1}_mct( cdata_{c1}x, x2{c1}_{c1}x, cdata_{c1}{c1}, x2{c1}_{c1}{c1})
!            call t_drvstopf  ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}')
!            call t_drvstartf ('driver_c2{c1}_infoexch',barrier=mpicom_CPL{ccc1}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc1}ID,'cpl2{ccc1}_run')
!            call t_drvstopf  ('driver_c2{c1}_infoexch')
!            call t_drvstopf  ('DRIVER_C2{c1}')
!         endif
!      endif

#if (defined NEW_BUDGET)
      if (iamin_CPLID .and. do_budgets) then
         call t_drvstartf ('DRIVER_BUDGET1',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
!{list} key="budget_fractions" number="1"
!         if ({ccc}_present) then
!            call seq_diag_{ccc}_mct(dom_{c}x, fractions_{c}x, {c}2x_{c}={c}2x_{c}x, x2{c}_{c}=x2{c}_{c}x)
!         endif

!{list} key="budget_fractions" number="2"
!         if ({ccc}_present) then
!            call seq_diag_rtm_mct(dom_{c}x, {c}2x_{c}={c}2x_{c}x)
!         endif

!{list} key="budget_fractions" number="3"
!         if ({ccc1}_present) then
!            call seq_diag_{ccc1}_mct(dom_{c1}x, fractions_{c1}x, {c1}2x_{c1}={c1}2x_{c1}x, x2{c1}_{c1}=x2{c1}_{c1}x, x{c2}{c1}_{c1}=x{c2}{c1}_{c1}x, {c3}2x_{c1}={c3}2x_{c1}x)
!         endif

!{list} key="budget_fractions" number="4"
!         if ({ccc}_present) then
!            call seq_diag_ice_mct(dom_{c}x, fractions_{c}x, x2{c}_{c}=x2{c}_{c}x)
!         endif

         call t_drvstopf  ('DRIVER_BUDGET1',cplrun=.true.,budget=.true.)
      endif
#endif

!{list} key="model_cpl" number="3"
!      !----------------------------------------------------------
!      ! {ccc} -> cpl
!      !----------------------------------------------------------
!
!      if ({ccc}_present .and. {ccc}run_alarm) then
!         if (iamin_CPL{ccc}ID) then
!            call t_drvstartf ('DRIVER_{c}2C',barrier=mpicom_CPL{ccc}ID)
!            call t_drvstartf ('driver_{c}2c_{ccc}{c}2{ccc}x',barrier=mpicom_CPL{ccc}ID)
!            call map_{ccc}{c}2{ccc}x_mct( cdata_{c}{c}, {c}2x_{c}{c}, cdata_{c}x, {c}2x_{c}x)
!            call t_drvstopf  ('driver_{c}2c_{ccc}{c}2{ccc}x')
!            call t_drvstartf ('driver_{c}2c_infoexch',barrier=mpicom_CPL{ccc}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc}ID,'{ccc}2cpl_run')
!            call t_drvstopf  ('driver_{c}2c_infoexch')
!            call t_drvstopf  ('DRIVER_{c}2C')
!         endif
!
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER{ccc}POST',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc}post_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c}x,{c}2x_{c}x,'recv {ccc}')
!               call t_drvstopf  ('driver_{ccc}post_diagav')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc}POST',cplrun=.true.)
!         endif
!
!      endif   ! run alarm, {ccc}_present

!{list} key="set_fractions"
!      !----------------------------------------------------------
!      ! Update fractions based on new {ccc1} fractions
!      !----------------------------------------------------------
!
!      if (iamin_CPLID) then
!         call t_drvstartf ('DRIVER_FRACSET',cplrun=.true.,barrier=mpicom_CPLID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!         call t_drvstartf ('driver_fracset_fracset',barrier=mpicom_CPLID)
!         call seq_frac_set({c1}2x_{c1}x, 
!                           cdata_{c2}x,
!                           {ccc3}_present,
!                           fractions_{c4}x,
!                           )
!         call t_drvstopf  ('driver_fracset_fracset')
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_FRACSET',cplrun=.true.)
!      endif

!{list} key="model_setup" number="5"
!      !----------------------------------------------------------
!      ! {ccc1} SETUP
!      !----------------------------------------------------------
!
!      if ({ccc1}_present .and. {ccc1}run_alarm) then
!
!         !----------------------------------------------------------
!         ! {ccc1} prep
!         !----------------------------------------------------------
!
!         if (iamin_CPLID) then
!         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!         call t_drvstartf ('DRIVER_{ccc1}PREP',cplrun=.true.,barrier=mpicom_CPLID)
!         if ({ccc1}_prognostic) then
!            if ({ccc2}_present) then
!               if (trim({c1}{c2}flux_grid) == '{ccc2}') then
!                  call t_drvstartf ('driver_{ccc1}prep_{ccc2}2{ccc2}2',barrier=mpicom_CPLID)
!                  call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, x{c1}{c2}_{c2}x, cdata_{c1}x, x{c1}{c2}_{c1}x, &
!                                        fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                        fluxlist=seq_flds_x{c1}{c2}_fluxes, statelist=seq_flds_x{c1}{c2}_states )
!                  call t_drvstopf  ('driver_{ccc1}prep_{ccc2}2{ccc1}2')
!               endif
!            endif
!            if ({ccc2}_present) then
!               call t_drvstartf ('driver_{ccc1}prep_{ccc2}2{ccc1}1',barrier=mpicom_CPLID)
!               call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x, {c2}2x_{c1}x, &
!                                     fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                     statelist=seq_flds_{c2}2x_states )
!               call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, {c2}2x_{c2}x, cdata_{c1}x,{c2}2x_{c1}x, &
!                                     fluxlist=seq_flds_{c2}2x_fluxes )
!               call map_{ccc2}2{ccc1}_mct( cdata_{c2}x, x{c1}{c2}_{c2}x, cdata_{c1}x, x{c1}{c2}_{c1}x, &
!                                     fractions_{c2}=fractions_{c2}x, fractions_{c1}=fractions_{c1}x, &
!                                     statelist=seq_flds_x{c1}{c2}_albedo )
!               call t_drvstopf  ('driver_{ccc1}prep_{ccc2}2{ccc1}1')
!            endif
!            if ({ccc3}_present) then
!               call t_drvstartf ('driver_{ccc1}prep_{ccc3}2{ccc1}',barrier=mpicom_CPLID)
!               call map_{ccc3}2{ccc1}_mct( cdata_{c3}x, {c3}2x_{c3}x, cdata_{c1}x, {c3}2x_{c1}x, &
!                                     fractions_{c3}=fractions_{c3}x, fractions_{c1}=fractions_{c1}x, &
!                                     fluxlist=seq_flds_{c3}2x_fluxes, statelist=seq_flds_{c3}2x_states )
!               call t_drvstopf  ('driver_{ccc1}prep_{ccc3}2{ccc1}')
!            endif
!            if ({ccc4}_present) then
!               call t_drvstartf ('driver_{ccc1}prep_{ccc4}2{ccc1}',barrier=mpicom_CPLID)
!               call map_{ccc4}2{ccc1}_mct( cdata_{c4}x, {c4}2x_{c4}x, cdata_{c1}x, {c4}2x_{c1}x , &
!                                     fractions_{c4}=fractions_{c4}x, fractions_{c1}=fractions_{c1}x, &
!                                     fluxlist=seq_flds_{c4}2x_fluxes, statelist=seq_flds_{c4}2x_states )
!               call t_drvstopf  ('driver_{ccc1}prep_{ccc4}2{ccc1}')
!            endif
!            call t_drvstartf ('driver_{ccc1}prep_mrgx2{c1}',barrier=mpicom_CPLID)
!            call mrg_x2{c1}_run_mct( cdata_{c1}x, {c4}2x_{c1}x, {c2}2x_{c1}x, x{c1}{c2}_{c1}x, {c3}2x_{c1}x, fractions_{c1}x, x2{c1}_{c1}x )
!            call t_drvstopf  ('driver_{ccc1}prep_mrgx2{c1}')
!
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc1}prep_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c1}x,x2{c1}_{c1}x,'send {ccc1}')
!               call t_drvstopf  ('driver_{ccc1}prep_diagav')
!            endif
!         endif  ! {ccc1}_prognostic
!         call t_drvstopf  ('DRIVER_{ccc1}PREP',cplrun=.true.)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         endif  ! CPLID
!
!         !----------------------------------------------------------
!         ! cpl -> {ccc1}
!         !----------------------------------------------------------
!
!         if (iamin_CPL{ccc1}ID .and. {ccc1}_prognostic) then
!            call t_drvstartf ('DRIVER_C2{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call t_drvstartf ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}',barrier=mpicom_CPL{ccc1}ID)
!            call map_{ccc1}x2{ccc1}{c1}_mct( cdata_{c1}x, x2{c1}_{c1}x, cdata_{c1}{c1}, x2{c1}_{c1}{c1})
!            call t_drvstopf  ('driver_c2{c1}_{ccc1}x2{ccc1}{c1}')
!            call t_drvstartf ('driver_c2{c1}_infoexch',barrier=mpicom_CPL{ccc1}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc1}ID,'cpl2{cc1}_run')
!            call t_drvstopf  ('driver_c2{c1}_infoexch')
!            call t_drvstopf  ('DRIVER_C2{c1}')
!         endif
!
!      endif

!{list} key="model_run" number=3
!      if ({ccc}_present .and. {ccc}run_alarm .and. iamin_{ccc}ID) then
!         if (run_barriers) then
!            call t_drvstartf ('DRIVER_{ccc}_RUN_BARRIER')
!            call mpi_barrier(mpicom_{ccc}ID,ierr)
!            call t_drvstopf ('DRIVER_{ccc}_RUN_BARRIER')
!         endif
!         call t_drvstartf ('DRIVER_{ccc}_RUN',barrier=mpicom_{ccc}ID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!         if ({ccc}_prognostic) call mct_avect_vecmult(x2{c}_{c}{c},drv2mdl_{c}{c},seq_flds_x2{c}_fluxes)
!!wangty modify
!#ifdef wrf
!         integration_phase = 1 !juanxiong he
!         call atm_run_mct( EClock_a, EClock_w, cdata_aa, cdata_cc, cdata_caca, x2a_aa, a2x_aa, &
!                           x2c_cc1, x2c_cc2, c2x_cc1, c2x_cc2, &
!                           x2ca_caca1, x2ca_caca2, ca2x_caca1, ca2x_caca2, &
!                           twoway_coupling, twoway_nudging, integration_phase )
!!juanxiong he
!#else
!         call {ccc}_run_mct( EClock_{c}, cdata_{c}{c}, x2{c}_{c}{c}, {c}2x_{c}{c})
!#endif
!         call mct_avect_vecmult({c}2x_{c}{c},mdl2drv_{c}{c},seq_flds_{c}2x_fluxes)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_{ccc}_RUN')
!      endif

!{list} key="model_run" number=4
!      !----------------------------------------------------------
!      ! Run {ccc1} Model
!      !----------------------------------------------------------
!
!      if ({ccc1}_present .and. {ccc1}run_alarm .and. iamin_{ccc1}ID) then
!         if (run_barriers) then
!            call t_drvstartf ('DRIVER_{ccc1}_RUN_BARRIER')
!            call mpi_barrier(mpicom_{ccc1}ID,ierr)
!            call t_drvstopf ('DRIVER_{ccc1}_RUN_BARRIER')
!         endif
!         call t_drvstartf ('DRIVER_{ccc1}_RUN',barrier=mpicom_{ccc1}ID)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc1}ID)
!         if ({ccc1}_prognostic) call mct_avect_vecmult(x2{c1}_{c1}{c1},drv2mdl_{c1}{c1},seq_flds_x2{c1}_fluxes)
!         call {ccc1}_run_mct( EClock_{c1}, cdata_{c1}{c1}, x2{c1}_{c1}{c1}, {c1}2x_{c1}{c1})
!         call mct_avect_vecmult({c1}2x_{c1}{c1},mdl2drv_{c1}{c1},seq_flds_{c1}2x_fluxes)
!         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!         call t_drvstopf  ('DRIVER_{ccc1}_RUN')
!      endif
!
!      !----------------------------------------------------------
!      ! {ccc1} -> cpl
!      !----------------------------------------------------------
!
!      if ({ccc1}_present .and. {ccc2}_prognostic .and. {ccc1}run_alarm) then
!
!         if (iamin_CPL{ccc1}ID) then
!            call t_drvstartf ('DRIVER_{c1}2C',barrier=mpicom_CPL{ccc1}ID)
!            call t_drvstartf ('driver_{c1}2c_{ccc1}{c1}2{ccc1}x',barrier=mpicom_CPL{ccc1}ID)
!            call map_{ccc1}{c1}2{ccc1}x_mct( cdata_{c1}{c1}, {c1}2x_{c1}{c1}, cdata_{c1}x, {c1}2x_{c1}x)
!            call t_drvstopf  ('driver_{c1}2c_{ccc1}{c1}2{ccc1}x')
!            call t_drvstartf ('driver_{c1}2c_infoexch',barrier=mpicom_CPL{ccc1}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc1}ID,'{ccc1}2cpl_run')
!            call t_drvstopf  ('driver_{c1}2c_infoexch')
!            call t_drvstopf  ('DRIVER_{c1}2C')
!         endif
!
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER_{ccc1}POST',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc1}post_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c1}x,{c1}2x_{c1}x,'recv {ccc1}')
!               call t_drvstopf  ('driver_{ccc1}post_diagav')
!            endif
!            call t_drvstartf ('driver_{ccc1}post_{ccc1}2{ccc2}',barrier=mpicom_CPLID)
!            call map_{ccc1}2{ccc2}_mct( cdata_{c1}x, {c1}2x_{c1}x, cdata_{c2}x, {c1}2x_{c2}x )
!            call t_drvstopf  ('driver_{ccc1}post_{ccc1}2{ccc2}')
!            call t_drvstartf ('driver_{ccc1}post_mrgx2{c2}',barrier=mpicom_CPLID)
!            call mrg_x2{c2}_run_mct( cdata_{c2}x, {c1}2x_{c2}x, x2{c3}_{c2}x )
!            call t_drvstopf  ('driver_{ccc1}post_mrgx2{c2}')
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc1}POST',cplrun=.true.)
!         endif
!
!      endif   ! run alarm, {ccc1}_present

!wangty modify
#ifdef wrf
      !----------------------------------------------------------
      ! atm -> cpl for wrf, juangxiong he
      !----------------------------------------------------------
      if (atm_present .and. atmrun_alarm) then
         if (iamin_CPLATMID) then
            call t_drvstartf ('DRIVER_A2C',barrier=mpicom_CPLATMID)
            call t_drvstartf ('driver_a2c_atma2atmx',barrier=mpicom_CPLATMID)
            call map_cama2camx_mct(cdata_cc, c2x_cc1, cdata_cx, c2x_cx1) !juanxiong he
            call map_cama2camx_mct(cdata_cc, c2x_cc2, cdata_cx, c2x_cx2) !juanxiong he
            call t_drvstopf  ('driver_a2c_atma2atmx')
            call t_drvstartf ('driver_a2c_infoexch',barrier=mpicom_CPLATMID)
            call seq_infodata_exchange(infodata,CPLATMID,'atm2cpl_run')
            call t_drvstopf  ('driver_a2c_infoexch')
            call t_drvstopf  ('DRIVER_A2C')
         endif
      endif ! run alarm

     !juanxiong he
     if (wrf_present .and. wrfrun_alarm) then

         !----------------------------------------------------------
         ! wrf prep
         !----------------------------------------------------------

         if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call t_drvstartf ('DRIVER_WRFPREP',cplrun=.true.,barrier=mpicom_CPLID)

         call t_drvstartf ('driver_wrfprep_cam2wrf',barrier=mpicom_CPLID)
         call map_cam2wrf_mct( cdata_cx, c2x_cx1, cdata_mx, c2x_mx1 , &
                                    statelist=seq_flds_c2x_states )
         call map_cam2wrf_mct( cdata_cx, c2x_cx2, cdata_mx, c2x_mx2 , &
                                    statelist=seq_flds_c2x_states )
         call t_drvstopf  ('driver_wrfprep_cam2wrf')

         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling mrg_x2w_run_mct'
         call mrg_x2w_run_mct( cdata_mx, c2x_mx1, x2m_mx1)
         call mrg_x2w_run_mct( cdata_mx, c2x_mx2, x2m_mx2)

         call t_drvstopf  ('DRIVER_WRFPREP',cplrun=.true.)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif  ! CPLID

         !----------------------------------------------------------
         ! cpl -> wrf, 3d
         !----------------------------------------------------------

         if (iamin_CPLWRFID .and. wrf_prognostic) then
            call t_drvstartf ('DRIVER_C2W',barrier=mpicom_CPLWRFID)
            call t_drvstartf ('driver_c2w_wrfx2wrfa',barrier=mpicom_CPLWRFID)
            call map_wrfx2wrfw_mct( cdata_mx, x2m_mx1, cdata_mm, x2m_mm1)   ! juaxniong he
            call map_wrfx2wrfw_mct( cdata_mx, x2m_mx2, cdata_mm, x2m_mm2)   ! juaxniong he
            call t_drvstopf  ('driver_c2w_wrfx2wrfa')
            call t_drvstartf ('driver_c2w_infoexch',barrier=mpicom_CPLWRFID)
            call seq_infodata_exchange(infodata,CPLWRFID,'cpl2wrf_run')
            call t_drvstopf  ('driver_c2w_infoexch')
            call t_drvstopf  ('DRIVER_C2W')
         endif

      endif

      !----------------------------------------------------------
      ! RUN wrf model
      !----------------------------------------------------------

      if (wrf_present .and. wrfrun_alarm .and. iamin_WRFID) then
         call t_drvstartf ('DRIVER_WRF_RUN',barrier=mpicom_WRFID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_WRFID)
         call wrf_run_mct( EClock_a, EClock_w, cdata_ww, cdata_mm, x2w_ww, w2x_ww, x2m_mm1, &
                           x2m_mm2, m2x_mm, twoway_coupling, twoway_nudging ) !juanxiong he
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_WRF_RUN')
      endif

      !----------------------------------------------------------
      ! wrf -> cpl
      !----------------------------------------------------------
      if (wrf_present) then
         if (iamin_CPLWRFID.and.twoway_coupling) then
            call t_drvstartf ('DRIVER_M2X',barrier=mpicom_CPLWRFID)
            call t_drvstartf ('driver_m2x_wrfw2wrfx',barrier=mpicom_CPLWRFID)
            call map_wrfw2wrfx_mct(cdata_mm, m2x_mm, cdata_mx, m2x_mx) !juanxiong he
            call t_drvstopf  ('driver_m2x_wrfw2wrfx')
            call t_drvstartf ('driver_m2x_infoexch',barrier=mpicom_CPLWRFID)
            call seq_infodata_exchange(infodata,CPLWRFID,'wrf2cpl_run')
            call t_drvstopf  ('driver_m2x_infoexch')
            call t_drvstopf  ('DRIVER_M2X')
         endif

         if (iamin_CPLID) then
            call t_drvstartf ('DRIVER_WRFPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_wrfpost_diagav',barrier=mpicom_CPLID)
               call seq_diag_avect_mct(cdata_mx,m2x_mx,'recv wrf')
               call t_drvstopf  ('driver_wrfpost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_WRFPOST',cplrun=.true.)
         endif
      endif
      !----------------------------------------------------------
      ! budget with new fractions
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! atm -> cpl for geatm, juanxiong he
      !----------------------------------------------------------
      if (atm_present .and. atmrun_alarm) then
         if (iamin_CPLATMID) then
            call t_drvstartf ('DRIVER_A2C',barrier=mpicom_CPLATMID)
            call t_drvstartf ('driver_a2ca_gatma2atmx',barrier=mpicom_CPLATMID)
            call map_gcama2camx_mct(cdata_caca, ca2x_caca1, cdata_cax, ca2x_cax1)  !juanxiong he
            call map_gcama2camx_mct(cdata_caca, ca2x_caca2, cdata_cax, ca2x_cax2)  !juanxiong he
            call t_drvstopf  ('driver_a2ca_gatma2atmx')
            call t_drvstartf ('driver_a2ca_infoexch',barrier=mpicom_CPLATMID)
            call seq_infodata_exchange(infodata,CPLATMID,'atm2cpl_run')
            call t_drvstopf  ('driver_a2ca_infoexch')
            call t_drvstopf  ('DRIVER_A2C')
         endif
      endif ! run alarm

     !juanxiong he
     if (geatm_present .and. gearun_alarm) then

         !----------------------------------------------------------
         ! gea prep
         !----------------------------------------------------------

         if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call t_drvstartf ('DRIVER_GEAPREP',cplrun=.true.,barrier=mpicom_CPLID)

         call t_drvstartf ('driver_geaprep_cam2gea',barrier=mpicom_CPLID)
         call map_cam2gea_mct( cdata_cax, ca2x_cax1, cdata_gex, ca2x_chemx1 , &
                               statelist=seq_flds_ca2x_states )
         call map_cam2gea_mct( cdata_cax, ca2x_cax2, cdata_gex, ca2x_chemx2 , &
                               statelist=seq_flds_ca2x_states )
         call t_drvstopf  ('driver_geaprep_cam2gea')

         if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling mrg_x2w_run_mct'
         call mrg_x2ge_run_mct( cdata_gex, ca2x_chemx1, x2chem_chemx1)
         call mrg_x2ge_run_mct( cdata_gex, ca2x_chemx2, x2chem_chemx2)

         call t_drvstopf  ('DRIVER_GEAPREP',cplrun=.true.)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif  ! CPLID

         !----------------------------------------------------------
         ! cpl -> gea, 3d
         !----------------------------------------------------------

         if (iamin_CPLGEAID .and. geatm_prognostic) then
            call t_drvstartf ('DRIVER_CA2GE',barrier=mpicom_CPLGEAID)
            call t_drvstartf ('driver_ca2ge_geax2geaa',barrier=mpicom_CPLGEAID)
            call map_geax2geaw_mct( cdata_gex, x2chem_chemx1, cdata_gege, x2chem_chemchem1)   ! juaxniong he
            call map_geax2geaw_mct( cdata_gex, x2chem_chemx2, cdata_gege, x2chem_chemchem2)   ! juaxniong he
            call t_drvstopf  ('driver_ca2ge_geax2geaa')
            call t_drvstartf ('driver_ca2ge_infoexch',barrier=mpicom_CPLGEAID)
            call seq_infodata_exchange(infodata,CPLGEAID,'cpl2gea_run')
            call t_drvstopf  ('driver_ca2ge_infoexch')
            call t_drvstopf  ('DRIVER_CA2GE')
         endif
      endif

      !----------------------------------------------------------
      ! RUN geatm model
      !----------------------------------------------------------

      if (geatm_present .and. gearun_alarm .and. iamin_GEAID) then
         call t_drvstartf ('DRIVER_GEA_RUN',barrier=mpicom_GEAID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GEAID)
         call geatm_run_mct( EClock_a, EClock_ge, cdata_gege, &
                             x2chem_chemchem1, x2chem_chemchem2, chem2x_chemchem ) !juanxiong he
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_GEA_RUN')
      endif

      !----------------------------------------------------------
      ! geatm -> cpl
      !----------------------------------------------------------
      if (geatm_present) then
         if (iamin_CPLGEAID.and.twoway_coupling) then
            call t_drvstartf ('DRIVER_GEA2X',barrier=mpicom_CPLGEAID)
            call t_drvstartf ('driver_chem2x_geaw2geax',barrier=mpicom_CPLGEAID)
            call map_geaw2geax_mct(cdata_gege, chem2x_chemchem, cdata_gex, chem2x_chemx)  !juanxiong he
            call t_drvstopf  ('driver_chem2x_geaw2geax')
            call t_drvstartf ('driver_chem2x_infoexch',barrier=mpicom_CPLGEAID)
            call seq_infodata_exchange(infodata,CPLGEAID,'gea2cpl_run')
            call t_drvstopf  ('driver_chem2x_infoexch')
            call t_drvstopf  ('DRIVER_GEA2X')
         endif

         if (iamin_CPLID) then
            call t_drvstartf ('DRIVER_GEAPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_geapost_diagav',barrier=mpicom_CPLID)
               call seq_diag_avect_mct(cdata_gex,chem2x_chemx,'recv geatm')
               call t_drvstopf  ('driver_geapost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_GEAPOST',cplrun=.true.)
         endif
      endif
      !----------------------------------------------------------
      ! RUN atm model, phase = 2, juanxiong he
      !----------------------------------------------------------

         !----------------------------------------------------------
         ! wrf -> cpl -> atm
         !----------------------------------------------------------

        if (iamin_CPLID.and.twoway_coupling) then

             ! wrf->cpl, juanxiong he
             if (wrf_present) then
               call t_drvstartf ('driver_atmprep_wrf2atm',barrier=mpicom_CPLID)
               call map_wrf2cam_mct( cdata_mx, m2x_mx, cdata_cx, m2x_cx , &
                                     fluxlist=seq_flds_w2x_fluxes, statelist=seq_flds_w2x_states )
               call t_drvstopf  ('driver_atmprep_wrf2atm')
            endif
            call t_drvstartf ('driver_atmprep_mrgx2c',barrier=mpicom_CPLID)
            if ( seq_comm_iamroot(CPLID)) write(logunit,F00) 'Calling mrg_x2c_run_mct'
            call mrg_x2c_run_mct( cdata_cx, m2x_cx, x2c_cx1 )
            call t_drvstopf  ('driver_atmprep_mrgx2c')
         endif

         if (iamin_CPLATMID .and. atm_prognostic.and.twoway_coupling) then
            !cpl->cam, juanxiong he
            call t_drvstartf ('DRIVER_C2A',barrier=mpicom_CPLATMID)
            call t_drvstartf ('driver_c2a_atmx2atma',barrier=mpicom_CPLATMID)
            call map_camx2cama_mct( cdata_cx, x2c_cx1, cdata_cc, x2c_cc1)
            call t_drvstopf  ('driver_c2a_atmx2atma')

            call t_drvstartf ('driver_c2a_infoexch',barrier=mpicom_CPLATMID)
            call seq_infodata_exchange(infodata,CPLATMID,'cpl2atm_run')
            call t_drvstopf  ('driver_c2a_infoexch')
            call t_drvstopf  ('DRIVER_C2A')
         endif
      if (atm_present .and. atmrun_alarm .and. iamin_ATMID) then
         call t_drvstartf ('DRIVER_ATM_RUN',barrier=mpicom_ATMID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_ATMID)
         integration_phase = 2 !juanxiong he
         call atm_run_mct( EClock_a, EClock_w, cdata_aa, cdata_cc, cdata_caca, x2a_aa, a2x_aa, &
                           x2c_cc1, x2c_cc2, c2x_cc1, c2x_cc2, &
                           x2ca_caca1, x2ca_caca2, ca2x_caca1, ca2x_caca2, &
                           twoway_coupling, twoway_nudging, integration_phase )
         call mct_avect_vecmult(a2x_aa,mdl2drv_aa,seq_flds_a2x_fluxes)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_ATM_RUN')
      endif
#endif

!{list} key="model_cpl" number=4
!      !----------------------------------------------------------
!      ! {ccc} -> cpl for ccsm, juanxiong he
!      !----------------------------------------------------------
!
!      if ({ccc}_present .and. {ccc1}run_alarm) then
!         if (iamin_CPL{ccc1}ID) then
!            call t_drvstartf ('DRIVER_{c}2C',barrier=mpicom_CPL{ccc}ID)
!            call t_drvstartf ('driver_{c}2c_{ccc}{c}2{ccc}x',barrier=mpicom_CPL{ccc}ID)
!            call map_{ccc}{c}2{ccc}x_mct( cdata_{c}{c}, {c}2x_{c}{c}, cdata_{c}x, {c}2x_{c}x)
!            call t_drvstopf  ('driver_{c}2c_{ccc}{c}2{ccc}x')
!            call t_drvstartf ('driver_{c}2c_infoexch',barrier=mpicom_CPL{ccc}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc}ID,'{ccc}2cpl_run')
!            call t_drvstopf  ('driver_{c}2c_infoexch')
!            call t_drvstopf  ('DRIVER_{c}2C')
!         endif
!
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER_{ccc}POST',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc}post_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c}x,{c}2x_{c}x,'recv {ccc}')
!               call t_drvstopf  ('driver_{ccc1}post_diagav')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc1}POST',cplrun=.true.)
!         endif
!
!      endif ! run alarm

      !----------------------------------------------------------
      ! budget with new fractions
      !----------------------------------------------------------

#if (defined NEW_BUDGET)
      if (iamin_CPLID .and. do_budgets) then
         call t_drvstartf ('DRIVER_BUDGET2',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
!{list} key="budget_newfractions" number=1
!         if ({ccc}_present) then
!            call seq_diag_{ccc}_mct(dom_{c}x, fractions_{c}x, {c}2x_{c}={c}2x_{c}x, x2{c}_{c}=x2{c}_{c}x)
!         endif

!{list} key="budget_newfractions" number=2
!         if ({ccc}_present) then
!            call seq_diag_{ccc}_mct(dom_{c}x, fractions_{c}x, {c}2x_{c}={c}2x_{c}x)
!         endif

         call t_drvstopf  ('DRIVER_BUDGET2',cplrun=.true.,budget=.true.)
      endif
      if (iamin_CPLID .and. do_budgets) then
         call t_drvstartf ('DRIVER_BUDGET3',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         call seq_diag_accum_mct()
         call t_drvstopf  ('DRIVER_BUDGET3',cplrun=.true.,budget=.true.)
      endif
      if (iamin_CPLID .and. do_budgets) then
         call t_drvstartf ('DRIVER_BUDGETF',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (.not. dead_comps) call seq_diag_print_mct(EClock_d,stop_alarm,budget_inst, &
            budget_daily, budget_month, budget_ann, budget_ltann, budget_ltend)
         call seq_diag_zero_mct(EClock=EClock_d)
         call t_drvstopf  ('DRIVER_BUDGETF',cplrun=.true.,budget=.true.)
      endif
#endif

!{list} key="ocn_def"
!      !----------------------------------------------------------
!      ! {ccc} -> cpl, loose coupling (concurrent type mode)
!      !----------------------------------------------------------
!
!      if (.not.{ccccc}_tight_coupling) then
!      if ({ccc}_present .and. {ccc}next_alarm) then
!         if (iamin_CPL{ccc}ID) then
!            call t_drvstartf ('DRIVER_{c}2C',barrier=mpicom_CPL{ccc}ID)
!            call t_drvstartf ('driver_{c}2c_{ccc}{c}2{ccc}x',barrier=mpicom_CPL{ccc}ID)
!            call map_{ccc}{c}2{ccc}x_mct( cdata_{c}{c}, {c}2x_{c}{c}, cdata_{c}x, {c}2x_{c}x)
!            call t_drvstopf  ('driver_{c}2c_{ccc}{c}2{ccc}x')
!            call t_drvstartf ('driver_{c}2c_infoexch',barrier=mpicom_CPL{ccc}ID)
!            call seq_infodata_exchange(infodata,CPL{ccc}ID,'{cc}2cpl_run')
!            call t_drvstopf  ('driver_{c}2c_infoexch')
!            call t_drvstopf  ('DRIVER_{c}2C')
!         endif
!         if (iamin_CPLID) then
!            call t_drvstartf ('DRIVER_{ccc}POST',cplrun=.true.,barrier=mpicom_CPLID)
!            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
!            if (info_debug > 1) then
!               call t_drvstartf ('driver_{ccc}post_diagav',barrier=mpicom_CPLID)
!               call seq_diag_avect_mct(cdata_{c}x,{c}2x_{c}x,'recv {ccc}')
!               call t_drvstopf  ('driver_{ccc}post_diagav')
!            endif
!            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!            call t_drvstopf  ('DRIVER_{ccc}POST',cplrun=.true.)
!         endif
!      endif
!      endif


      !----------------------------------------------------------
      ! Save driver level restart information
      !----------------------------------------------------------

      if ( restart_alarm .and. iamin_CPLID) then
         call t_drvstartf ('DRIVER_RESTART',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            write(logunit,104) ' Write restart file at ',ymd,tod
            call shr_sys_flush(logunit)
         endif
         call seq_rest_write(EClock_d,seq_SyncClock)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_RESTART',cplrun=.true.)
      endif


      !----------------------------------------------------------
      ! Write history file, only AVs on CPLID
      !----------------------------------------------------------

      if (iamin_CPLID) then

         call t_drvstartf ('DRIVER_HISTORY',cplrun=.true.,barrier=mpicom_CPLID)
         if ( history_alarm) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (iamroot_CPLID) then
               write(logunit,104) ' Write history file at ',ymd,tod
               call shr_sys_flush(logunit)
            endif
            call seq_hist_write(EClock_d)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif

         if (do_histavg) then
            call seq_hist_writeavg(EClock_d,histavg_alarm)
         endif

!{list} key="do_hist" number=1
!         if (do_hist_{c}2x) then
!            if (trim(hist_{c}2x_flds) == 'all') then
!               call seq_hist_writeaux(EClock_d,'{c}2x','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,ncpl)
!            else
!               call seq_hist_writeaux(EClock_d,'{c}2x','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,ncpl,&
!                                                flds=hist_{c}2x_flds)
!            endif
!         endif
!         if (do_hist_{c}2x3hr) then
!            if (trim(hist_{c}2x3hr_flds) == 'all') then
!               call seq_hist_writeaux(EClock_d,'{c}2x3h','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,8,t3hr_alarm)
!            else
!               call seq_hist_writeaux(EClock_d,'{c}2x3h','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,8,t3hr_alarm,&
!                                                 flds=hist_{c}2x3hr_flds)
!            end if
!         endif
!         if (do_hist_{c}2x3hrp) then
!            if (trim(hist_{c}2x3hrp_flds) == 'all') then
!               call seq_hist_writeaux(EClock_d,'{c}2x3h_prec','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,8,t3hr_alarm)
!            else
!               call seq_hist_writeaux(EClock_d,'{c}2x3h_prec','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,8,t3hr_alarm,&
!                                              flds=hist_{c}2x3hrp_flds)
!            end if
!         endif
!         if (do_hist_{c}2x24hr) then
!            call seq_hist_writeaux(EClock_d,'{c}2x1d','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,1,t24hr_alarm)
!         endif

!{list} key="do_hist" number=2
!         if (do_hist_{c}2x) then
!            call seq_hist_writeaux(EClock_d,'{c}2x','dom{c}',cdata_{c}x,{c}2x_{c}x,{ccc}_nx,{ccc}_ny,ncpl)
!         endif

         call t_drvstopf  ('DRIVER_HISTORY',cplrun=.true.)

      end if

      ! --- End timestep clock/timing diagnostics
      call t_drvstartf ('DRIVER_TSTAMP_WRITE',cplrun=.true.)
      if (tod == 0 .or. info_debug > 1) then
         if (iamroot_CPLID) then
            call date_and_time(dstr,tstr)
            Time_estep = mpi_wtime()
            dtstep = time_estep-time_bstep
            dtstep_acc = dtstep_acc + dtstep
            dtstep_cnt = dtstep_cnt + 1
            write(logunit,101) ' tStamp_write: model date = ',ymd,tod, &
               ' wall clock = ',dstr(1:4),'-',dstr(5:6),'-',dstr(7:8),' ',&
                                tstr(1:2),':',tstr(3:4),':',tstr(5:6), &
               ' avg dt = ',dtstep_acc/dtstep_cnt,' dt = ',dtstep
            Time_bstep = mpi_wtime()
            call shr_sys_flush(logunit)
         endif
      endif

!{list} key="memory_write"
!      if (tod == 0) then
!         if (iamroot_CPLID .or. iamroot_{ccc}ID
!            ) then
!            call shr_mem_getusage(msize,mrss)
!            write(logunit,105) ' memory_write: model date = ',ymd,tod, &
!               ' memory = ',mrss,' MB (highwater)    ',msize,' MB (usage)', &
!               '  (pe=',iam_GLOID,' comps=',trim(complist)//')'
!         endif
!      endif

      if (info_debug > 1) then
         if (iamroot_CPLID) then
            call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday)
!            write(logunit,106) ' nextsw_cday = ',nextsw_cday
            write(logunit,*) '  nextsw_cday = ',nextsw_cday
         endif
      endif
      call t_drvstopf  ('DRIVER_TSTAMP_WRITE',cplrun=.true.)

      call t_stopf  ('DRIVER_RUN_LOOP')
      ! --- Write out performance data
      call t_drvstartf  ('DRIVER_TPROF_WRITE',cplrun=.true.)
      if (tprof_alarm) then
         call t_startf("sync1_tprof")
         call mpi_barrier(mpicom_GLOID,ierr)
         call t_stopf("sync1_tprof")

         write(timing_file,'(a,i8.8,a1,i5.5)') trim(tchkpt_dir)//"/ccsm_timing_",ymd,"_",tod
         call t_prf(filename=trim(timing_file), mpicom=mpicom_GLOID, &
                    num_outpe=1)

         call t_startf("sync2_tprof")
         call mpi_barrier(mpicom_GLOID,ierr)
         call t_stopf("sync2_tprof")
      endif
      call t_drvstopf  ('DRIVER_TPROF_WRITE',cplrun=.true.)

   end do   ! driver run loop
   call t_startf ('DRIVER_RUN_LOOP_BSTOP')
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf ('DRIVER_RUN_LOOP_BSTOP')

   Time_end = mpi_wtime()

   !----------------------------------------------------------
   ! Ending of basic time step loop
   !----------------------------------------------------------

end subroutine ccsm_run

!===============================================================================

subroutine ccsm_final()

   !------------------------------------------------------------------------
   ! Finalization of all models
   !------------------------------------------------------------------------

 103  format( 5A )
   ! TODO finalize routines need to be cleaned up

   call t_barrierf ('DRIVER_FINAL_BARRIER', mpicom_GLOID)
   call t_startf ('DRIVER_FINAL')

   call seq_timemgr_EClockGetData( EClock_d, stepno=endstep)
   call shr_mem_getusage(msize,mrss)

!{list} key="model_final" 
!   if (iamin_{ccc}ID) then
!      if (drv_threading) call seq_comm_setnthreads(nthreads_{ccc}ID)
!      call {ccc}_final_mct( )
!      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
!   endif


!wangty modify
#ifdef wrf
   if (iamin_WRFID) then  !juanxiong he
      if (drv_threading) call seq_comm_setnthreads(nthreads_WRFID)
      call wrf_final_mct( )
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif
   if (iamin_GEAID) then  !juanxiong he
      if (drv_threading) call seq_comm_setnthreads(nthreads_GEAID)
      call geatm_final_mct( )
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif
#endif


   !------------------------------------------------------------------------
   ! End the run cleanly
   !------------------------------------------------------------------------

   call seq_io_finalize()

   call shr_mpi_min(msize,msize0,mpicom_GLOID,'driver msize0',all=.true.)
   call shr_mpi_max(msize,msize1,mpicom_GLOID,'driver msize1',all=.true.)
   call shr_mpi_min(mrss,mrss0,mpicom_GLOID,'driver mrss0',all=.true.)
   call shr_mpi_max(mrss,mrss1,mpicom_GLOID,'driver mrss1',all=.true.)
   if ( iamroot_CPLID )then
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod, dtime=dtime)
      write(logunit,'(//)')
      write(logunit,FormatA) subname, 'SUCCESSFUL TERMINATION OF CPL7-CCSM'
      write(logunit,FormatD) subname, '  at YMD,TOD = ',ymd,tod
      simDays = (endStep-begStep)*dtime/(24._r8*3600._r8)
      write(logunit,FormatR) subname, '# simulated days (this run) = ', simDays 
      write(logunit,FormatR) subname, 'compute time (hrs)          = ', (Time_end-Time_begin)/3600._r8
      if ( (Time_end /= Time_begin) .and. (simDays /= 0.0_r8) )then
         SYPD = shr_const_cday*simDays/(days_per_year*(Time_end-Time_begin))
         write(logunit,FormatR) subname, '# simulated years / cmp-day = ', SYPD
      endif
      write(logunit,FormatR) subname,' pes min memory highwater  (MB)  = ',mrss0
      write(logunit,FormatR) subname,' pes max memory highwater  (MB)  = ',mrss1
      write(logunit,FormatR) subname,' pes min memory last usage (MB)  = ',msize0
      write(logunit,FormatR) subname,' pes max memory last usage (MB)  = ',msize1
      write(logunit,'(//)')
      close(logunit)
   endif

   call t_stopf  ('DRIVER_FINAL')
   call t_prf(trim(timing_dir)//'/ccsm_timing', mpicom_GLOID)
   call t_finalizef()

end subroutine ccsm_final


!===============================================================================

subroutine t_drvstartf(string,cplrun,budget,barrier)

   implicit none

   character(len=*),intent(in) :: string
   logical,intent(in),optional :: cplrun
   logical,intent(in),optional :: budget
   integer,intent(in),optional :: barrier

   character(len=128) :: strbar
   character(len=*),parameter :: strcpl = 'DRIVER_CPL_RUN'
   character(len=*),parameter :: strbud = 'DRIVER_BUDGET'
   logical :: lcplrun,lbudget

!-------------------------------------------------------------------------------

   lcplrun  = .false.
   lbudget  = .false.
   if (present(cplrun)) then
      lcplrun = cplrun
   endif
   if (present(budget)) then
      lbudget = budget
   endif

   if (present(barrier)) then
      strbar = trim(string)//'_BARRIER'
      call t_barrierf (trim(strbar), barrier)
   endif

   if (lcplrun) then
      call t_startf   (trim(strcpl))
   endif

   if (lbudget) then
      call t_startf   (trim(strbud))
   endif

   call t_startf   (trim(string))

end subroutine t_drvstartf

!===============================================================================

subroutine t_drvstopf(string,cplrun,budget)

   implicit none

   character(len=*),intent(in) :: string
   logical,intent(in),optional :: cplrun
   logical,intent(in),optional :: budget

   character(len=128) :: strbar
   character(len=*),parameter :: strcpl = 'DRIVER_CPL_RUN'
   character(len=*),parameter :: strbud = 'DRIVER_BUDGET'
   logical :: lcplrun,lbudget

!-------------------------------------------------------------------------------

   lcplrun = .false.
   lbudget = .false.
   if (present(cplrun)) then
      lcplrun = cplrun
   endif
   if (present(budget)) then
      lbudget = budget
   endif

!  strbar = trim(string)//'_BARRIER'

   call t_stopf   (trim(string))

   if (lbudget) then
      call t_stopf   (trim(strbud))
   endif

   if (lcplrun) then
      call t_stopf   (trim(strcpl))
   endif

end subroutine t_drvstopf

!===============================================================================

subroutine seq_ccsm_printlogheader()

  !-----------------------------------------------------------------------
  !
  ! Purpose: Print basic information on what this driver program is
  ! to the logfile.
  !
  !-----------------------------------------------------------------------
  !
  ! Local variables
  !
  implicit none

  character(len=8) :: cdate          ! System date
  character(len=8) :: ctime          ! System time
  integer          :: values(8)
  character        :: date*8, time*10, zone*5

!-------------------------------------------------------------------------------

  call date_and_time (date, time, zone, values)
  cdate(1:2) = date(5:6)
  cdate(3:3) = '/'
  cdate(4:5) = date(7:8)
  cdate(6:6) = '/'
  cdate(7:8) = date(3:4)
  ctime(1:2) = time(1:2)
  ctime(3:3) = ':'
  ctime(4:5) = time(3:4)
  ctime(6:6) = ':'
  ctime(7:8) = time(5:6)
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '     NCAR CPL7 Community Climate System Model (CCSM)  '
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '     (Online documentation is available on the CCSM '
  write(logunit,F00) '      Models page: http://www.ccsm.ucar.edu/models/ '
  write(logunit,F00) '      License information is available as a link from above '
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '                DATE ',cdate, ' TIME ', ctime 
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,*)' '
  write(logunit,*)' '

end subroutine seq_ccsm_printlogheader


#ifdef ESMF_INTERFACE

!===============================================================================

subroutine ccsm_comp_init(comp, importState, exportState, clock, rc)
   type(ESMF_CplComp)   :: comp
   type(ESMF_State)     :: importState, exportState
   type(ESMF_Clock)     :: clock
   integer, intent(out) :: rc

   ! Local variables
   type(ESMF_State)    :: attState
   integer             :: localrc
!{list} key="ESMF_comp"
!   type(ESMF_GridComp) :: {ccc}Comp
   type(ESMF_VM)       :: vm
   character(len=80)   :: str
   integer             :: rc2
   integer, dimension(1) :: rootList

   rc = ESMF_SUCCESS

   !------
   ! Create a state object to which the field level attributes will be
   ! attached, and link the state to the specified component
   !------
   attState = ESMF_StateCreate(name="ccsm_atts", rc=localrc)
   if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to create state for attributes')

   call ESMF_AttributeLink(comp, attState, rc=localrc)
   if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attributes')

!{list} key="model_register"
!   !------
!   ! Create and setup the model components
!   !------
!   call {ccc}_register({ccc}_petlist, comp, {ccc}Comp)

   !------
   ! Process the CESM initialization
   !------
   call ccsm_init()

   !------
   ! Set the application and field level attributes
   !------
   call esmfshr_attribute_appl_init(comp, rc=localrc)
   !call esmfshr_attribute_fields_init(attState, rc=localrc)

   !------
   ! Get the VM and root pet list to be used for the AttributeUpdate call
   !------
! get current
   !call ESMF_VMGetGlobal(vm, rc=localrc)
   call ESMF_VMGetCurrent(vm, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to get VM')

!{list} key="rootlist" number=1
!   rootList(1) = {ccc}_petlist(1)
!
!! use this one
!   call ESMF_AttributeUpdate({ccc}Comp, vm, rootList=atm_petlist, rc=localrc)
!!   call ESMF_AttributeUpdate({ccc}Comp, vm, rootList=rootList, rc=localrc)
!! try using just first pet in petlist
!! /= ESMF_SUCCESS
!! Add ESMF_Finalize w/ abort
!! Macro in ESMF_Attribute - #define DEBUG
!   if (localrc /= ESMF_SUCCESS) then
!      write(logunit,*) ' '
!      write(logunit,*) 'KDS: Error updating {ccc} attributes: ', localrc
!      call ESMF_Finalize()
!      call shr_sys_abort('failed to update {ccc} attributes')
!   endif

!{list} key="rootlist" number=2
!   call ESMF_AttributeUpdate({ccc}Comp, vm, rootList={ccc}_petlist, rc=localrc)

   !------
   ! Write out all of the attributes to the CIM compliant XML file
   !------
   if (iamroot_GLOID) then

      call ESMF_AttributeWrite( &
              comp, &
              convention='CIM 1.0', &
              purpose='Model Component Simulation Description', &
              attwriteflag=ESMF_ATTWRITE_XML, rc=localrc)

   endif

   rc = localrc

end subroutine ccsm_comp_init
!===============================================================================

subroutine ccsm_comp_run(comp, importState, exportState, clock, rc)
   type(ESMF_CplComp)   :: comp
   type(ESMF_State)     :: importState, exportState
   type(ESMF_Clock)     :: clock
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   call ccsm_run()

end subroutine ccsm_comp_run

!===============================================================================

subroutine ccsm_comp_final(comp, importState, exportState, clock, rc)
   type(ESMF_CplComp)   :: comp
   type(ESMF_State)     :: importState, exportState
   type(ESMF_Clock)     :: clock
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   call ccsm_final()

end subroutine ccsm_comp_final


!===============================================================================
!
! This subroutine registers the initialization, run and finalization routines
! for the specified coupler component.
!
subroutine ccsm_comp_register(comp, rc)
     type(ESMF_CplComp) :: comp
     integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   call ESMF_CplCompSetEntryPoint(comp, ESMF_SETINIT, &
                                  userRoutine=ccsm_comp_init, rc=rc)
   call ESMF_CplCompSetEntryPoint(comp, ESMF_SETRUN, &
                                  userRoutine=ccsm_comp_run, rc=rc)
   call ESMF_CplCompSetEntryPoint(comp, ESMF_SETFINAL, &
                                  userRoutine=ccsm_comp_final, rc=rc)

end subroutine ccsm_comp_register

!===============================================================================

#endif

end module ccsm_comp_mod


