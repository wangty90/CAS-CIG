      !----------------------------------------------------------
      ! ocn -> cpl, tight coupling (sequential type mode)
      !----------------------------------------------------------

      if (ocn_present .and. ocnnext_alarm) then
         if (iamin_CPLOCNID) then
            call t_drvstartf ('DRIVER_O2C',barrier=mpicom_CPLOCNID)
            call t_drvstartf ('driver_o2c_ocno2ocnx',barrier=mpicom_CPLOCNID)
            call map_ocno2ocnx_mct( cdata_oo, o2x_oo, cdata_ox, o2x_ox)
            call t_drvstopf  ('driver_o2c_ocno2ocnx')
            call t_drvstartf ('driver_o2c_infoexch',barrier=mpicom_CPLOCNID)
            call seq_infodata_exchange(infodata,CPLOCNID,'ocn2cpl_run')
            call t_drvstopf  ('driver_o2c_infoexch')
            call t_drvstopf  ('DRIVER_O2C')
         endif
         if (iamin_CPLID) then
            call t_drvstartf  ('DRIVER_OCNPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_ocnpost_diagav',barrier=mpicom_CPLID)
               call seq_diag_avect_mct(cdata_ox,o2x_ox,'recv ocn')
               call t_drvstopf  ('driver_ocnpost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_OCNPOST',cplrun=.true.)
         endif
      endif
