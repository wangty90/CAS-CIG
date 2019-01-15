[start module_use]
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep
[end module_use]
[start var_model_to_cpl]
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
[end var_model_to_cpl]
[start var_cpl_to_model]
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
[end var_cpl_to_model]
[start indices_set_var]
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
[end indices_set_var]
[start init_var_x2c]
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
[end init_var_x2c]
[start init_var_c2x]
    index_a2x_Sa_co2prog    = mct_avect_indexra(a2x,'Sa_co2prog',perrWith='quiet')
    index_a2x_Sa_co2diag    = mct_avect_indexra(a2x,'Sa_co2diag',perrWith='quiet')
[end init_var_c2x]
[start other_clean]
!wangty modify
#ifdef wrf
    call mct_aVect_clean(c2x)
    call mct_aVect_clean(x2c)
    call mct_aVect_clean(ca2x)
    call mct_aVect_clean(x2ca)
#endif
[end other_clean]
