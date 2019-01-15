[start var_model_to_cpl]
  public :: clm_cpl_indices_set        ! Set the coupler indices
[end var_model_to_cpl]
[start indices_set_use]
  use seq_flds_mod  , only: seq_flds_x2l_fields, seq_flds_l2x_fields,     &
                            seq_flds_x2s_fields, seq_flds_s2x_fields,     &
                                                 seq_flds_r2x_fields
  use mct_mod       , only: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                            mct_aVect_clean, mct_avect_nRattr
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep
    implicit none
[end indices_set_use]
[start indices_set_var]
    type(mct_aVect) :: r2x      ! temporary, runoff to coupler
    type(mct_aVect) :: s2x      ! temporary, glacier to coupler
    type(mct_aVect) :: x2s      ! temporary, coupler to glacier
    character(len=32) :: subname = 'clm_cpl_indices_set'  ! subroutine name
[end indices_set_var]
[start other_init]
    call mct_aVect_init(r2x, rList=seq_flds_r2x_fields, lsize=1)
    call mct_aVect_init(x2s, rList=seq_flds_x2s_fields, lsize=1)
    call mct_aVect_init(s2x, rList=seq_flds_s2x_fields, lsize=1)
[end other_init]
[start init_var_x2c]
    index_l2x_Fall_flxvoc1  = mct_avect_indexra(l2x,'Fall_flxvoc1' ,perrwith='quiet')
    index_l2x_Fall_flxvoc2  = mct_avect_indexra(l2x,'Fall_flxvoc2' ,perrwith='quiet')
    index_l2x_Fall_flxvoc3  = mct_avect_indexra(l2x,'Fall_flxvoc3' ,perrwith='quiet')
    index_l2x_Fall_flxvoc4  = mct_avect_indexra(l2x,'Fall_flxvoc4' ,perrwith='quiet')
    index_l2x_Fall_flxvoc5  = mct_avect_indexra(l2x,'Fall_flxvoc5' ,perrwith='quiet')
    index_l2x_Fall_fco2_lnd = mct_avect_indexra(l2x,'Fall_fco2_lnd',perrwith='quiet')
    if ( lnd_drydep )then
       index_l2x_Sl_ddvel = mct_avect_indexra(l2x, trim(drydep_fields_token))
    else
       index_l2x_Sl_ddvel = 0
    end if

    nflds_l2x = mct_avect_nRattr(l2x)

    nflds_r2x = mct_avect_nRattr(r2x)
    ! sno -> drv

#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
    index_s2x_Ss_tsrf01   = mct_avect_indexra(s2x,'Ss_tsrf01')
    index_s2x_Ss_topo01   = mct_avect_indexra(s2x,'Ss_topo01')
    index_s2x_Fgss_qice01 = mct_avect_indexra(s2x,'Fgss_qice01')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
    index_s2x_Ss_tsrf02   = mct_avect_indexra(s2x,'Ss_tsrf02')
    index_s2x_Ss_topo02   = mct_avect_indexra(s2x,'Ss_topo02')
    index_s2x_Ss_tsrf03   = mct_avect_indexra(s2x,'Ss_tsrf03')
    index_s2x_Ss_topo03   = mct_avect_indexra(s2x,'Ss_topo03')
    index_s2x_Fgss_qice02 = mct_avect_indexra(s2x,'Fgss_qice02')
    index_s2x_Fgss_qice03 = mct_avect_indexra(s2x,'Fgss_qice03')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
    index_s2x_Ss_tsrf04   = mct_avect_indexra(s2x,'Ss_tsrf04')
    index_s2x_Ss_topo04   = mct_avect_indexra(s2x,'Ss_topo04')
    index_s2x_Ss_tsrf05   = mct_avect_indexra(s2x,'Ss_tsrf05')
    index_s2x_Ss_topo05   = mct_avect_indexra(s2x,'Ss_topo05')
    index_s2x_Fgss_qice04 = mct_avect_indexra(s2x,'Fgss_qice04')
    index_s2x_Fgss_qice05 = mct_avect_indexra(s2x,'Fgss_qice05')
#endif
#if (defined GLC_NEC_10 )
    index_s2x_Ss_tsrf06   = mct_avect_indexra(s2x,'Ss_tsrf06')
    index_s2x_Ss_topo06   = mct_avect_indexra(s2x,'Ss_topo06')
    index_s2x_Ss_tsrf07   = mct_avect_indexra(s2x,'Ss_tsrf07')
    index_s2x_Ss_topo07   = mct_avect_indexra(s2x,'Ss_topo07')
    index_s2x_Ss_tsrf08   = mct_avect_indexra(s2x,'Ss_tsrf08')
    index_s2x_Ss_topo08   = mct_avect_indexra(s2x,'Ss_topo08')
    index_s2x_Ss_tsrf09   = mct_avect_indexra(s2x,'Ss_tsrf09')
    index_s2x_Ss_topo09   = mct_avect_indexra(s2x,'Ss_topo09')
    index_s2x_Ss_tsrf10   = mct_avect_indexra(s2x,'Ss_tsrf10')
    index_s2x_Ss_topo10   = mct_avect_indexra(s2x,'Ss_topo10')
    index_s2x_Fgss_qice06 = mct_avect_indexra(s2x,'Fgss_qice06')
    index_s2x_Fgss_qice07 = mct_avect_indexra(s2x,'Fgss_qice07')
    index_s2x_Fgss_qice08 = mct_avect_indexra(s2x,'Fgss_qice08')
    index_s2x_Fgss_qice09 = mct_avect_indexra(s2x,'Fgss_qice09')
    index_s2x_Fgss_qice10 = mct_avect_indexra(s2x,'Fgss_qice10')
#endif

    nflds_s2x = mct_avect_nRattr(s2x)
[end init_var_x2c]
[start init_var_c2x]
    index_x2l_Sa_co2prog    = mct_avect_indexra(x2l,'Sa_co2prog',perrwith='quiet')
    index_x2l_Sa_co2diag    = mct_avect_indexra(x2l,'Sa_co2diag',perrwith='quiet')

    nflds_x2l = mct_avect_nRattr(x2l)
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
    index_x2s_Sg_frac01   = mct_avect_indexra(x2s,'Sg_frac01')
    index_x2s_Sg_topo01   = mct_avect_indexra(x2s,'Sg_topo01')
    index_x2s_Fsgg_rofi01 = mct_avect_indexra(x2s,'Fsgg_rofi01')
    index_x2s_Fsgg_rofl01 = mct_avect_indexra(x2s,'Fsgg_rofl01')
    index_x2s_Fsgg_hflx01 = mct_avect_indexra(x2s,'Fsgg_hflx01')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
    index_x2s_Sg_frac02   = mct_avect_indexra(x2s,'Sg_frac02')
    index_x2s_Sg_topo02   = mct_avect_indexra(x2s,'Sg_topo02')
    index_x2s_Sg_frac03   = mct_avect_indexra(x2s,'Sg_frac03')
    index_x2s_Sg_topo03   = mct_avect_indexra(x2s,'Sg_topo03')
    index_x2s_Fsgg_rofi02 = mct_avect_indexra(x2s,'Fsgg_rofi02')
    index_x2s_Fsgg_rofl02 = mct_avect_indexra(x2s,'Fsgg_rofl02')
    index_x2s_Fsgg_hflx02 = mct_avect_indexra(x2s,'Fsgg_hflx02')
    index_x2s_Fsgg_rofi03 = mct_avect_indexra(x2s,'Fsgg_rofi03')
    index_x2s_Fsgg_rofl03 = mct_avect_indexra(x2s,'Fsgg_rofl03')
    index_x2s_Fsgg_hflx03 = mct_avect_indexra(x2s,'Fsgg_hflx03')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
    index_x2s_Sg_frac04   = mct_avect_indexra(x2s,'Sg_frac04')
    index_x2s_Sg_topo04   = mct_avect_indexra(x2s,'Sg_topo04')
    index_x2s_Sg_frac05   = mct_avect_indexra(x2s,'Sg_frac05')
    index_x2s_Sg_topo05   = mct_avect_indexra(x2s,'Sg_topo05')
    index_x2s_Fsgg_rofi04 = mct_avect_indexra(x2s,'Fsgg_rofi04')
    index_x2s_Fsgg_rofl04 = mct_avect_indexra(x2s,'Fsgg_rofl04')
    index_x2s_Fsgg_hflx04 = mct_avect_indexra(x2s,'Fsgg_hflx04')
    index_x2s_Fsgg_rofi05 = mct_avect_indexra(x2s,'Fsgg_rofi05')
    index_x2s_Fsgg_rofl05 = mct_avect_indexra(x2s,'Fsgg_rofl05')
    index_x2s_Fsgg_hflx05 = mct_avect_indexra(x2s,'Fsgg_hflx05')
#endif
#if (defined GLC_NEC_10 )
    index_x2s_Sg_frac06   = mct_avect_indexra(x2s,'Sg_frac06')
    index_x2s_Sg_topo06   = mct_avect_indexra(x2s,'Sg_topo06')
    index_x2s_Sg_frac07   = mct_avect_indexra(x2s,'Sg_frac07')
    index_x2s_Sg_topo07   = mct_avect_indexra(x2s,'Sg_topo07')
    index_x2s_Sg_frac08   = mct_avect_indexra(x2s,'Sg_frac08')
    index_x2s_Sg_topo08   = mct_avect_indexra(x2s,'Sg_topo08')
    index_x2s_Sg_frac09   = mct_avect_indexra(x2s,'Sg_frac09')
    index_x2s_Sg_topo09   = mct_avect_indexra(x2s,'Sg_topo09')
    index_x2s_Sg_frac10   = mct_avect_indexra(x2s,'Sg_frac10')
    index_x2s_Sg_topo10   = mct_avect_indexra(x2s,'Sg_topo10')
    index_x2s_Fsgg_rofi06 = mct_avect_indexra(x2s,'Fsgg_rofi06')
    index_x2s_Fsgg_rofl06 = mct_avect_indexra(x2s,'Fsgg_rofl06')
    index_x2s_Fsgg_hflx06 = mct_avect_indexra(x2s,'Fsgg_hflx06')
    index_x2s_Fsgg_rofi07 = mct_avect_indexra(x2s,'Fsgg_rofi07')
    index_x2s_Fsgg_rofl07 = mct_avect_indexra(x2s,'Fsgg_rofl07')
    index_x2s_Fsgg_hflx07 = mct_avect_indexra(x2s,'Fsgg_hflx07')
    index_x2s_Fsgg_rofi08 = mct_avect_indexra(x2s,'Fsgg_rofi08')
    index_x2s_Fsgg_rofl08 = mct_avect_indexra(x2s,'Fsgg_rofl08')
    index_x2s_Fsgg_hflx08 = mct_avect_indexra(x2s,'Fsgg_hflx08')
    index_x2s_Fsgg_rofi09 = mct_avect_indexra(x2s,'Fsgg_rofi09')
    index_x2s_Fsgg_rofl09 = mct_avect_indexra(x2s,'Fsgg_rofl09')
    index_x2s_Fsgg_hflx09 = mct_avect_indexra(x2s,'Fsgg_hflx09')
    index_x2s_Fsgg_rofi10 = mct_avect_indexra(x2s,'Fsgg_rofi10')
    index_x2s_Fsgg_rofl10 = mct_avect_indexra(x2s,'Fsgg_rofl10')
    index_x2s_Fsgg_hflx10 = mct_avect_indexra(x2s,'Fsgg_hflx10')
#endif
    nflds_x2s = mct_avect_nRattr(x2s)
[end init_var_c2x]
[start other_clean]
    call mct_aVect_clean(r2x)
    call mct_aVect_clean(x2s)
    call mct_aVect_clean(s2x)
[end other_clean]
