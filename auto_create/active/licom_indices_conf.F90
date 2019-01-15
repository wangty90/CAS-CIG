[start init_var_x2c]
    index_o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2' ,perrWith='quiet')

    index_o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')
[end init_var_x2c]
[start init_var_c2x]
    index_x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
    index_x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')
[end init_var_c2x]
