module cam_CplIndices
  
  use seq_flds_mod
  use mct_mod
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep

  implicit none

  SAVE
  public                               ! By default make data private
!wangty modify
#ifdef wrf
    !-----------------------------------------------------------------
    ! juanxiong he
    !-----------------------------------------------------------------

  ! drv -> cam, three dimension (flux), juanxiong he

  integer, dimension(1:num_cam_levs) :: index_x2c_Fcxx_dudt      ! heat from radiation tendency (K/s)
  integer, dimension(1:num_cam_levs) :: index_x2c_Fcxx_dvdt      ! heat from pbl tendency (K/s)
  integer, dimension(1:num_cam_levs) :: index_x2c_Fcxx_dtdt      ! heat from cumulus tendency (K/s)
  integer, dimension(1:num_cam_levs) :: index_x2c_Fcxx_dqdt     ! heat from microphysic tendency (K/s)

  integer, dimension(1:num_cam_levs) :: index_x2c_Sx_u3d      ! u 
  integer, dimension(1:num_cam_levs) :: index_x2c_Sx_v3d      ! v
  integer, dimension(1:num_cam_levs) :: index_x2c_Sx_t3d      ! t
  integer, dimension(1:num_cam_levs) :: index_x2c_Sx_q3d      ! q
  integer :: index_x2c_Sx_ps      ! ps

  ! cam -> drv, three dimension (scalar), juanxiong he

  integer, dimension(1:num_cam_levs)   :: index_c2x_Sc_z3d            ! atm level height
  integer, dimension(1:num_cam_levs)    :: index_c2x_Sc_u3d            ! atm level zon wind
  integer, dimension(1:num_cam_levs)    :: index_c2x_Sc_v3d            ! atm level mer wind
  integer, dimension(1:num_cam_levs)    :: index_c2x_Sc_t3d            ! atm level temp
  integer, dimension(1:num_cam_levs)    :: index_c2x_Sc_w3d         ! atm level vert wind
  integer, dimension(1:num_cam_levs)    :: index_c2x_Sc_q3d         ! atm level spec hum
  integer, dimension(1:num_cam_levs)    :: index_c2x_Sc_p3d           ! atm level pressure
  integer, dimension(1:num_cam_levs)  :: index_c2x_Fcxc_utend         ! atm level u wind tendency 
  integer, dimension(1:num_cam_levs)  :: index_c2x_Fcxc_vtend         ! atm level v wind tendency
  integer, dimension(1:num_cam_levs)  :: index_c2x_Fcxc_ttend           ! atm level t tendency
  integer, dimension(1:num_cam_levs)  :: index_c2x_Fcxc_qtend           ! atm level q tendency
  integer :: index_c2x_Sc_ps         ! atm surface pressure
  integer :: index_c2x_Sc_phis         ! atm surface geopotential  height
  integer :: index_c2x_Sc_lat         ! atm latitude
  integer :: index_c2x_Sc_lon         ! atm longitude
  integer :: index_c2x_Sc_ts         !surface temperature
  integer :: index_c2x_Sc_sst ! sst
  integer :: index_c2x_Sc_snowhland  ! snow height over land
  integer :: index_c2x_Sc_snowhice !snow height over ice
  integer :: index_c2x_Sc_seaice  ! seaice
  integer :: index_c2x_Sc_ocnfrac !ocn fraction
  integer, dimension(1:num_soil_layers)  :: index_c2x_Sc_soildepth ! soil layer 1 depth
  integer, dimension(1:num_soil_layers)  :: index_c2x_Sc_soilthick ! soil layer 1 thickness
  integer, dimension(1:num_soil_layers)  :: index_c2x_Sc_soilt ! soil layer 1 temperature
  integer, dimension(1:num_soil_layers)  :: index_c2x_Sc_soilm ! soil layer 1 moisture

! drv -> geatm, three dimension
  ! drv -> cam, three dimension (flux)
  integer, dimension(1:num_cam_levs,1:num_tracers) :: index_x2ca_Fcaxx_tracer ! for radiation

  ! cam -> drv, three dimension (scalar)
  integer, dimension(1:num_cam_levs) :: index_ca2x_Sca_z3d            ! atm level height
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_u3d            ! atm level zon wind
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_v3d            ! atm level mer wind
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_t3d            ! atm level temp
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_qv3d         ! atm level qv
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_qc3d         ! atm level qc
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_qi3d         ! atm level qi
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_p3d          ! atm level pressure
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_rh3d        ! atm level rh
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_taucldi3d     ! atm level taucldi
  integer, dimension(1:num_cam_levs)  :: index_ca2x_Sca_taucldv3d     ! atm level taucldc
  integer :: index_ca2x_Sca_ps         ! atm surface pressure
  integer :: index_ca2x_Sca_phis         ! atm surface geopotential height
  integer :: index_ca2x_Sca_lat         ! atm latitude
  integer :: index_ca2x_Sca_lon         ! atm longittude
  integer :: index_ca2x_Sca_ts         !surface temperature
  integer :: index_ca2x_Sca_sst        ! sst
  integer :: index_ca2x_Sca_snowhland  ! snow height over land
  integer :: index_ca2x_Sca_snowhice !snow height over ice
  integer :: index_ca2x_Sca_seaice  ! seaice
  integer :: index_ca2x_Sca_ocnfrac !ocn fraction
  integer :: index_ca2x_Sca_t2 ! t2
  integer :: index_ca2x_Sca_q2 ! q2
  integer :: index_ca2x_Sca_rh2 ! q2
  integer :: index_ca2x_Sca_u10 ! u10
  integer :: index_ca2x_Sca_v10 ! v10

  integer :: index_ca2x_Sca_ust ! ust
  integer :: index_ca2x_Sca_rmol ! rmol
  integer :: index_ca2x_Sca_pblh ! pblh
  integer :: index_ca2x_Sca_raincv ! raincv
  integer :: index_ca2x_Sca_rainncv ! rainncv
  integer :: index_ca2x_Sca_swdown !swdown
  integer :: index_ca2x_Sca_clflo !clflo
  integer :: index_ca2x_Sca_clfmi !clfmi
  integer :: index_ca2x_Sca_clfhi !clfhi

  integer, dimension(1:num_soil_layers) :: index_ca2x_Sca_soildepth ! soil layer depth
  integer, dimension(1:num_soil_layers)  :: index_ca2x_Sca_soilthick ! soil layer thickness
  integer, dimension(1:num_soil_layers)  :: index_ca2x_Sca_soilt ! soil layer temperature
  integer, dimension(1:num_soil_layers)  :: index_ca2x_Sca_soilm ! soil layer moisture

    !-----------------------------------------------------------------
    ! juanxiong he
    !-----------------------------------------------------------------
#endif
!wangty modify
#ifdef wrf
!--------------------------------------------------------------------------------------
! soildepth/height and soil temperature/moisture for wrf/cam coupling, added
! juanxiong he
!--------------------------------------------------------------------------------------
  integer, dimension(1:num_soil_layers) :: index_x2a_Sx_soildepth ! soil layer depth
  integer, dimension(1:num_soil_layers)  :: index_x2a_Sx_soilthick ! soil layer thickness
  integer, dimension(1:num_soil_layers)  :: index_x2a_Sx_soilt ! soil layer temperature
  integer, dimension(1:num_soil_layers)  :: index_x2a_Sx_soilm ! soil layer moisture
!--------------------------------------------------------------------------------------
! soildepth/height and soil temperature/moisture for wrf/cam coupling, added
! juanxiong he
!--------------------------------------------------------------------------------------
#endif

integer :: index_a2x_Sa_z
integer :: index_a2x_Sa_u
integer :: index_a2x_Sa_v
integer :: index_a2x_Sa_tbot
integer :: index_a2x_Sa_ptem
integer :: index_a2x_Sa_shum
integer :: index_a2x_Sa_dens
integer :: index_a2x_Sa_pbot
integer :: index_a2x_Sa_pslv
integer :: index_a2x_Faxa_lwdn
integer :: index_a2x_Faxa_rainc
integer :: index_a2x_Faxa_rainl
integer :: index_a2x_Faxa_snowc
integer :: index_a2x_Faxa_snowl
integer :: index_a2x_Faxa_swndr
integer :: index_a2x_Faxa_swvdr
integer :: index_a2x_Faxa_swndf
integer :: index_a2x_Faxa_swvdf
integer :: index_a2x_Faxa_swnet
integer :: index_a2x_Faxa_bcphidry
integer :: index_a2x_Faxa_bcphodry
integer :: index_a2x_Faxa_bcphiwet
integer :: index_a2x_Faxa_ocphidry
integer :: index_a2x_Faxa_ocphodry
integer :: index_a2x_Faxa_ocphiwet
integer :: index_a2x_Faxa_dstwet1
integer :: index_a2x_Faxa_dstwet2
integer :: index_a2x_Faxa_dstwet3
integer :: index_a2x_Faxa_dstwet4
integer :: index_a2x_Faxa_dstdry1
integer :: index_a2x_Faxa_dstdry2
integer :: index_a2x_Faxa_dstdry3
integer :: index_a2x_Faxa_dstdry4
integer :: index_a2x_Sa_co2prog
integer :: index_a2x_Sa_co2diag

integer :: index_x2a_Sx_t
integer :: index_x2a_So_t
integer :: index_x2a_Sa_lfrac
integer :: index_x2a_Sa_ifrac
integer :: index_x2a_Sa_ofrac
integer :: index_x2a_Sx_tref
integer :: index_x2a_Sx_qref
integer :: index_x2a_Sx_avsdr
integer :: index_x2a_Sx_anidr
integer :: index_x2a_Sx_avsdf
integer :: index_x2a_Sx_anidf
integer :: index_x2a_Sl_snowh
integer :: index_x2a_Si_snowh
integer :: index_x2a_Sl_fv
integer :: index_x2a_Sl_ram1
integer :: index_x2a_Faxx_taux
integer :: index_x2a_Faxx_tauy
integer :: index_x2a_Faxx_lat
integer :: index_x2a_Faxx_sen
integer :: index_x2a_Faxx_lwup
integer :: index_x2a_Faxx_evap
integer :: index_x2a_Fall_flxdst1
integer :: index_x2a_Fall_flxdst2
integer :: index_x2a_Fall_flxdst3
integer :: index_x2a_Fall_flxdst4
integer :: index_x2a_Faxx_flxvoc1
integer :: index_x2a_Faxx_flxvoc2
integer :: index_x2a_Faxx_flxvoc3
integer :: index_x2a_Faxx_flxvoc4
integer :: index_x2a_Faxx_flxvoc5
integer :: index_x2a_So_ustar
integer :: index_x2a_So_re
integer :: index_x2a_So_ssq
integer :: index_x2a_Sx_ddvel
integer :: index_x2a_Sx_u10
integer :: index_x2a_Faxx_fco2_lnd
integer :: index_x2a_Faxx_fco2_ocn
integer :: index_x2a_Faxx_fdms_ocn

contains


  subroutine {ccc}_cpl_indices_set( )

    type(mct_aVect) :: a2x      ! temporary
    type(mct_aVect) :: x2a      ! temporary
!wangty modify
#ifdef wrf
    type(mct_aVect) :: x2c      ! temporary, juanxiong he
    type(mct_aVect) :: c2x      ! temporary, juanxiong he
    type(mct_aVect) :: ca2x      ! temporary, juanxiong he
    type(mct_aVect) :: x2ca      ! temporary, juanxiong he

    integer :: i,j,k  ! juanxiong he

!-----------------------------------------------------------------------------------------------------------
!    for wrf-cam coupling
!    Juanxiong He, 05/10/2010
!-----------------------------------------------------------------------------------------------------------
#endif

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=1)
    call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=1)


index_a2x_Sa_z = mct_avect_indexra(a2x,'Sa_z')
index_a2x_Sa_u = mct_avect_indexra(a2x,'Sa_u')
index_a2x_Sa_v = mct_avect_indexra(a2x,'Sa_v')
index_a2x_Sa_tbot = mct_avect_indexra(a2x,'Sa_tbot')
index_a2x_Sa_ptem = mct_avect_indexra(a2x,'Sa_ptem')
index_a2x_Sa_shum = mct_avect_indexra(a2x,'Sa_shum')
index_a2x_Sa_dens = mct_avect_indexra(a2x,'Sa_dens')
index_a2x_Sa_pbot = mct_avect_indexra(a2x,'Sa_pbot')
index_a2x_Sa_pslv = mct_avect_indexra(a2x,'Sa_pslv')
index_a2x_Faxa_lwdn = mct_avect_indexra(a2x,'Faxa_lwdn')
index_a2x_Faxa_rainc = mct_avect_indexra(a2x,'Faxa_rainc')
index_a2x_Faxa_rainl = mct_avect_indexra(a2x,'Faxa_rainl')
index_a2x_Faxa_snowc = mct_avect_indexra(a2x,'Faxa_snowc')
index_a2x_Faxa_snowl = mct_avect_indexra(a2x,'Faxa_snowl')
index_a2x_Faxa_swndr = mct_avect_indexra(a2x,'Faxa_swndr')
index_a2x_Faxa_swvdr = mct_avect_indexra(a2x,'Faxa_swvdr')
index_a2x_Faxa_swndf = mct_avect_indexra(a2x,'Faxa_swndf')
index_a2x_Faxa_swvdf = mct_avect_indexra(a2x,'Faxa_swvdf')
index_a2x_Faxa_swnet = mct_avect_indexra(a2x,'Faxa_swnet')
index_a2x_Faxa_bcphidry = mct_avect_indexra(a2x,'Faxa_bcphidry')
index_a2x_Faxa_bcphodry = mct_avect_indexra(a2x,'Faxa_bcphodry')
index_a2x_Faxa_bcphiwet = mct_avect_indexra(a2x,'Faxa_bcphiwet')
index_a2x_Faxa_ocphidry = mct_avect_indexra(a2x,'Faxa_ocphidry')
index_a2x_Faxa_ocphodry = mct_avect_indexra(a2x,'Faxa_ocphodry')
index_a2x_Faxa_ocphiwet = mct_avect_indexra(a2x,'Faxa_ocphiwet')
index_a2x_Faxa_dstwet1 = mct_avect_indexra(a2x,'Faxa_dstwet1')
index_a2x_Faxa_dstwet2 = mct_avect_indexra(a2x,'Faxa_dstwet2')
index_a2x_Faxa_dstwet3 = mct_avect_indexra(a2x,'Faxa_dstwet3')
index_a2x_Faxa_dstwet4 = mct_avect_indexra(a2x,'Faxa_dstwet4')
index_a2x_Faxa_dstdry1 = mct_avect_indexra(a2x,'Faxa_dstdry1')
index_a2x_Faxa_dstdry2 = mct_avect_indexra(a2x,'Faxa_dstdry2')
index_a2x_Faxa_dstdry3 = mct_avect_indexra(a2x,'Faxa_dstdry3')
index_a2x_Faxa_dstdry4 = mct_avect_indexra(a2x,'Faxa_dstdry4')
index_a2x_Sa_co2prog = mct_avect_indexra(a2x,'Sa_co2prog')
index_a2x_Sa_co2diag = mct_avect_indexra(a2x,'Sa_co2diag')
!wangty modify
#ifdef wrf
    call mct_aVect_init(c2x, rList=seq_flds_c2x_fields, lsize=1)  ! for wrf/cam, juanxiong he
    call mct_aVect_init(x2c, rList=seq_flds_x2c_fields, lsize=1)  ! for wrf/cam, juanxiong he

    call mct_aVect_init(ca2x, rList=seq_flds_ca2x_fields, lsize=1)  ! for geatm/cam, juanxiong he
    call mct_aVect_init(x2ca, rList=seq_flds_x2ca_fields, lsize=1)  ! for geatm/cam, juanxiong he
#endif
    index_x2a_Faxx_fco2_lnd = mct_avect_indexra(x2a,'Faxx_fco2_lnd',perrWith='quiet')
    index_x2a_Faxx_fco2_ocn = mct_avect_indexra(x2a,'Faxx_fco2_ocn',perrWith='quiet')

    !index_x2a_Faxx_fdms_ocn =
    !mct_avect_indexra(x2a,'Faxx_fdms_ocn',perrWith='quiet')
    index_x2a_Faxx_fdms_ocn  = mct_avect_indexra(x2a,'Faxx_fdms' ,perrWith='quiet')

    index_x2a_Faxx_flxvoc1  = mct_avect_indexra(x2a,'Fall_flxvoc1' ,perrWith='quiet')
    index_x2a_Faxx_flxvoc2  = mct_avect_indexra(x2a,'Fall_flxvoc2' ,perrWith='quiet')
    index_x2a_Faxx_flxvoc3  = mct_avect_indexra(x2a,'Fall_flxvoc3' ,perrWith='quiet')
    index_x2a_Faxx_flxvoc4  = mct_avect_indexra(x2a,'Fall_flxvoc4' ,perrWith='quiet')
    index_x2a_Faxx_flxvoc5  = mct_avect_indexra(x2a,'Fall_flxvoc5' ,perrWith='quiet')
    if ( lnd_drydep )then
       index_x2a_Sx_ddvel   = mct_avect_indexra(x2a, trim(drydep_fields_token))
    else
       index_x2a_Sx_ddvel   = 0
    end if
!wangty modify
#ifdef wrf
    !-----------------------------------------------------------------
    ! juanxiong he
    !-----------------------------------------------------------------
    do k=1,num_soil_layers
      index_x2a_Sx_soildepth(k) = mct_avect_indexra(x2a,'Sx_soildepth'//slayer(k))  ! soil layer  depth
      index_x2a_Sx_soilthick(k) = mct_avect_indexra(x2a,'Sx_soilthick'//slayer(k))! soil layer  thickness
      index_x2a_Sx_soilt(k) = mct_avect_indexra(x2a,'Sx_soilt'//slayer(k))! soil layer  temperature
      index_x2a_Sx_soilm(k) = mct_avect_indexra(x2a,'Sx_soilm'//slayer(k))! soil layer  moisture
    end do

    do k = 1, num_cam_levs
    index_x2c_Fcxx_dudt(k)     = mct_avect_indexra(x2c,'Fcxx_dudt'//clev(k))
    index_x2c_Fcxx_dvdt(k)      = mct_avect_indexra(x2c,'Fcxx_dvdt'//clev(k))
    index_x2c_Fcxx_dtdt(k)      = mct_avect_indexra(x2c,'Fcxx_dtdt'//clev(k))
    index_x2c_Fcxx_dqdt(k)     = mct_avect_indexra(x2c,'Fcxx_dqdt'//clev(k))
    index_x2c_Sx_u3d(k)      = mct_avect_indexra(x2c,'Sx_u3d'//clev(k))
    index_x2c_Sx_v3d(k)      = mct_avect_indexra(x2c,'Sx_v3d'//clev(k))
    index_x2c_Sx_t3d(k)      = mct_avect_indexra(x2c,'Sx_t3d'//clev(k))
    index_x2c_Sx_q3d(k)      = mct_avect_indexra(x2c,'Sx_q3d'//clev(k))
    enddo
    index_x2c_Sx_ps      = mct_avect_indexra(x2c,'Sx_ps')

    do k = 1, num_cam_levs
    index_c2x_Sc_z3d(k)          = mct_avect_indexra(c2x,'Sc_z3d'//clev(k) )
    index_c2x_Sc_u3d(k)          = mct_avect_indexra(c2x,'Sc_u3d'//clev(k) )
    index_c2x_Sc_v3d(k)       = mct_avect_indexra(c2x,'Sc_v3d'//clev(k) )
    index_c2x_Sc_t3d(k)         = mct_avect_indexra(c2x,'Sc_t3d'//clev(k) )
    index_c2x_Sc_w3d(k)        = mct_avect_indexra(c2x,'Sc_w3d'//clev(k) )
    index_c2x_Sc_q3d(k)          = mct_avect_indexra(c2x,'Sc_q3d'//clev(k) )
    index_c2x_Sc_p3d(k)          = mct_avect_indexra(c2x,'Sc_p3d'//clev(k) )
    index_c2x_Fcxc_utend(k)          = mct_avect_indexra(c2x,'Fcxc_utend'//clev(k) )
    index_c2x_Fcxc_vtend(k)          = mct_avect_indexra(c2x,'Fcxc_vtend'//clev(k) )
    index_c2x_Fcxc_ttend(k)          = mct_avect_indexra(c2x,'Fcxc_ttend'//clev(k) )
    index_c2x_Fcxc_qtend(k)          = mct_avect_indexra(c2x,'Fcxc_qtend'//clev(k) )
    end do
    index_c2x_Sc_ps        = mct_avect_indexra(c2x,'Sc_ps')
    index_c2x_Sc_phis        = mct_avect_indexra(c2x,'Sc_phis')
    index_c2x_Sc_lat        = mct_avect_indexra(c2x,'Sc_lat')
    index_c2x_Sc_lon        = mct_avect_indexra(c2x,'Sc_lon')
    index_c2x_Sc_ts        = mct_avect_indexra(c2x,'Sc_ts')
    index_c2x_Sc_sst        = mct_avect_indexra(c2x,'Sc_sst')
    index_c2x_Sc_snowhland  = mct_avect_indexra(c2x,'Sc_snowhland')
    index_c2x_Sc_snowhice  = mct_avect_indexra(c2x,'Sc_snowhice')
    index_c2x_Sc_seaice  = mct_avect_indexra(c2x,'Sc_seaice')
    index_c2x_Sc_ocnfrac  = mct_avect_indexra(c2x,'Sc_ocnfrac')
    do k = 1, num_soil_layers
    index_c2x_Sc_soildepth(k)   = mct_avect_indexra(c2x,'Sc_soildepth'//slayer(k) )
    index_c2x_Sc_soilthick(k)   = mct_avect_indexra(c2x,'Sc_soilthick'//slayer(k) )
    index_c2x_Sc_soilt(k)   = mct_avect_indexra(c2x,'Sc_soilt'//slayer(k) )
    index_c2x_Sc_soilm(k)   = mct_avect_indexra(c2x,'Sc_soilm'//slayer(k) )
    end do
  ! cam -> drv, three dimension (flux)
    do k = 1, num_cam_levs
    index_ca2x_Sca_z3d(k)          = mct_avect_indexra(ca2x,'Sca_z3d'//clev(k) )
    index_ca2x_Sca_u3d(k)          = mct_avect_indexra(ca2x,'Sca_u3d'//clev(k) )
    index_ca2x_Sca_v3d(k)          = mct_avect_indexra(ca2x,'Sca_v3d'//clev(k) )
    index_ca2x_Sca_t3d(k)          = mct_avect_indexra(ca2x,'Sca_t3d'//clev(k) )
    index_ca2x_Sca_rh3d(k)          = mct_avect_indexra(ca2x,'Sca_rh3d'//clev(k) )
    index_ca2x_Sca_qv3d(k)          = mct_avect_indexra(ca2x,'Sca_qv3d'//clev(k) )
    index_ca2x_Sca_qi3d(k)          = mct_avect_indexra(ca2x,'Sca_qi3d'//clev(k) )
    index_ca2x_Sca_qc3d(k)          = mct_avect_indexra(ca2x,'Sca_qc3d'//clev(k) )
    index_ca2x_Sca_p3d(k)          = mct_avect_indexra(ca2x,'Sca_p3d'//clev(k) ) 
    index_ca2x_Sca_taucldi3d(k)          = mct_avect_indexra(ca2x,'Sca_taucldi3d'//clev(k) )
    index_ca2x_Sca_taucldv3d(k)          = mct_avect_indexra(ca2x,'Sca_taucldv3d'//clev(k) )
    enddo
    index_ca2x_Sca_ps        = mct_avect_indexra(ca2x,'Sca_ps')
    index_ca2x_Sca_phis      = mct_avect_indexra(ca2x,'Sca_phis')
    index_ca2x_Sca_lat       = mct_avect_indexra(ca2x,'Sca_lat')
    index_ca2x_Sca_lon       = mct_avect_indexra(ca2x,'Sca_lon')
    index_ca2x_Sca_ts        = mct_avect_indexra(ca2x,'Sca_ts')
    index_ca2x_Sca_sst        = mct_avect_indexra(ca2x,'Sca_sst')
    index_ca2x_Sca_snowhland  = mct_avect_indexra(ca2x,'Sca_snowhland')
    index_ca2x_Sca_snowhice  = mct_avect_indexra(ca2x,'Sca_snowhice')
    index_ca2x_Sca_seaice  = mct_avect_indexra(ca2x,'Sca_seaice')
    index_ca2x_Sca_ocnfrac  = mct_avect_indexra(ca2x,'Sca_ocnfrac')
    index_ca2x_Sca_t2  = mct_avect_indexra(ca2x,'Sca_t2')
    index_ca2x_Sca_q2  = mct_avect_indexra(ca2x,'Sca_q2')
    index_ca2x_Sca_rh2  = mct_avect_indexra(ca2x,'Sca_rh2')
    index_ca2x_Sca_u10  = mct_avect_indexra(ca2x,'Sca_u10')
    index_ca2x_Sca_v10  = mct_avect_indexra(ca2x,'Sca_v10')
    index_ca2x_Sca_ust  = mct_avect_indexra(ca2x,'Sca_ust')
    index_ca2x_Sca_rmol  = mct_avect_indexra(ca2x,'Sca_rmol')
    index_ca2x_Sca_pblh  = mct_avect_indexra(ca2x,'Sca_pblh')
    index_ca2x_Sca_rainncv  = mct_avect_indexra(ca2x,'Sca_rainncv')
    index_ca2x_Sca_swdown  = mct_avect_indexra(ca2x,'Sca_swdown')
    index_ca2x_Sca_clflo  = mct_avect_indexra(ca2x,'Sca_clflo')
    index_ca2x_Sca_clfmi  = mct_avect_indexra(ca2x,'Sca_clfmi')
    index_ca2x_Sca_clfhi  = mct_avect_indexra(ca2x,'Sca_clfhi')
     do k = 1, num_soil_layers
    index_ca2x_Sca_soildepth(k)   = mct_avect_indexra(ca2x,'Sca_soildepth'//slayer(k) )
    index_ca2x_Sca_soilthick(k)   = mct_avect_indexra(ca2x,'Sca_soilthick'//slayer(k) )
    index_ca2x_Sca_soilt(k)   = mct_avect_indexra(ca2x,'Sca_soilt'//slayer(k) )
    index_ca2x_Sca_soilm(k)   = mct_avect_indexra(ca2x,'Sca_soilm'//slayer(k) )
    end do

  ! drv -> cam, three dimension (scalar)
    do k = 1, num_cam_levs
    do i=1, num_tracers
    index_x2ca_Fcaxx_tracer(k,i)          = mct_avect_indexra(x2ca,'Fcaxx_tracer'//clev(k)//ctracer(i))
    end do
    end do

    !-----------------------------------------------------------------
    ! juanxiong he
    !-----------------------------------------------------------------
#endif

index_x2a_Sx_t = mct_avect_indexra(x2a,'Sx_t')
index_x2a_So_t = mct_avect_indexra(x2a,'So_t')
index_x2a_Sa_lfrac = mct_avect_indexra(x2a,'Sa_lfrac')
index_x2a_Sa_ifrac = mct_avect_indexra(x2a,'Sa_ifrac')
index_x2a_Sa_ofrac = mct_avect_indexra(x2a,'Sa_ofrac')
index_x2a_Sx_tref = mct_avect_indexra(x2a,'Sx_tref')
index_x2a_Sx_qref = mct_avect_indexra(x2a,'Sx_qref')
index_x2a_Sx_avsdr = mct_avect_indexra(x2a,'Sx_avsdr')
index_x2a_Sx_anidr = mct_avect_indexra(x2a,'Sx_anidr')
index_x2a_Sx_avsdf = mct_avect_indexra(x2a,'Sx_avsdf')
index_x2a_Sx_anidf = mct_avect_indexra(x2a,'Sx_anidf')
index_x2a_Sl_snowh = mct_avect_indexra(x2a,'Sl_snowh')
index_x2a_Si_snowh = mct_avect_indexra(x2a,'Si_snowh')
index_x2a_Sl_fv = mct_avect_indexra(x2a,'Sl_fv')
index_x2a_Sl_ram1 = mct_avect_indexra(x2a,'Sl_ram1')
index_x2a_Faxx_taux = mct_avect_indexra(x2a,'Faxx_taux')
index_x2a_Faxx_tauy = mct_avect_indexra(x2a,'Faxx_tauy')
index_x2a_Faxx_lat = mct_avect_indexra(x2a,'Faxx_lat')
index_x2a_Faxx_sen = mct_avect_indexra(x2a,'Faxx_sen')
index_x2a_Faxx_lwup = mct_avect_indexra(x2a,'Faxx_lwup')
index_x2a_Faxx_evap = mct_avect_indexra(x2a,'Faxx_evap')
index_x2a_Fall_flxdst1 = mct_avect_indexra(x2a,'Fall_flxdst1')
index_x2a_Fall_flxdst2 = mct_avect_indexra(x2a,'Fall_flxdst2')
index_x2a_Fall_flxdst3 = mct_avect_indexra(x2a,'Fall_flxdst3')
index_x2a_Fall_flxdst4 = mct_avect_indexra(x2a,'Fall_flxdst4')
index_x2a_Faxx_flxvoc1 = mct_avect_indexra(x2a,'Faxx_flxvoc1')
index_x2a_Faxx_flxvoc2 = mct_avect_indexra(x2a,'Faxx_flxvoc2')
index_x2a_Faxx_flxvoc3 = mct_avect_indexra(x2a,'Faxx_flxvoc3')
index_x2a_Faxx_flxvoc4 = mct_avect_indexra(x2a,'Faxx_flxvoc4')
index_x2a_Faxx_flxvoc5 = mct_avect_indexra(x2a,'Faxx_flxvoc5')
index_x2a_So_ustar = mct_avect_indexra(x2a,'So_ustar')
index_x2a_So_re = mct_avect_indexra(x2a,'So_re')
index_x2a_So_ssq = mct_avect_indexra(x2a,'So_ssq')
index_x2a_Sx_ddvel = mct_avect_indexra(x2a,'Sx_ddvel')
index_x2a_Sx_u10 = mct_avect_indexra(x2a,'Sx_u10')
index_x2a_Faxx_fco2_lnd = mct_avect_indexra(x2a,'Faxx_fco2_lnd')
index_x2a_Faxx_fco2_ocn = mct_avect_indexra(x2a,'Faxx_fco2_ocn')
index_x2a_Faxx_fdms_ocn = mct_avect_indexra(x2a,'Faxx_fdms_ocn')
    index_a2x_Sa_co2prog    = mct_avect_indexra(a2x,'Sa_co2prog',perrWith='quiet')
    index_a2x_Sa_co2diag    = mct_avect_indexra(a2x,'Sa_co2diag',perrWith='quiet')

    call mct_aVect_clean(x2a)
    call mct_aVect_clean(a2x)
    !call mct_aVect_clean(r2x_o) !LPF 20121219

!wangty modify
#ifdef wrf
    call mct_aVect_clean(c2x)
    call mct_aVect_clean(x2c)
    call mct_aVect_clean(ca2x)
    call mct_aVect_clean(x2ca)
#endif

  end subroutine cam_CplIndicesSet

end module cam_CplIndices

