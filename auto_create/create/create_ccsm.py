#!/usr/bin/python
import re
import ConfigParser
import types

srcfile = "ccsm_comp_mct_template.F90"
src = open(srcfile,'r')
lines = src.readlines()
src.close()

config_file_path="conf_ccsm_comp.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)

mutual1 = cf.sections()
count = 1
mutual = {}
for k in mutual1:
        mutual[count] = k
        count = count+1

mutual_value = {}
for k in mutual:
        mutual_value[k] = cf.items(mutual[k])

for k in mutual:
    listwrite = []
    i = 0
    pair = mutual[k][0:3]
    while i < len(lines):
        flag = 1
        line = lines[i]
	keytemp={}
        if '{list}' in line:
            keys = re.findall(r'key\=\"(.*?)\"',line)
            keyslist = cf.options(mutual[k])
            number = re.findall(r'number\=(\d)',line)
	    print k
	    print i
	    print keys
            if keys[0] in keyslist and cf.get(mutual[k],keys[0])[1] == number[0]:
                line = lines[i+1][1:]
	    elif keys[0] in keyslist and (keys[0] == "mapping" or keys[0] == "aux_files" or keys[0] == "aux_files_24hr" or keys[0] == "aux_files_3hr" or keys[0] == "init_map" or keys[0] == "flux_albedos" or keys[0] == "call_map" or keys[0] == "call_mrg" or keys[0] == "model_setup" or keys[0] == "init_fractions" or keys[0] == "set_fractions" or keys[0] == "model_run" or keys[0] == "model_cpl" or keys[0] == "prognostic_expect_others" or keys[0] == "model_prep" or keys[0] == "budget_fractions"):
		keytemp = keys
		line = lines[i+1][1:]
            else:
                line = lines[i+2][1:]
                if not line.strip():
                    flag = 0
            i = i + 2
        elif '<list>' in line:
	    flag = 0
            tmplist = []
            keys = re.findall(r'key\=\"(.*?)\"',line)
            keyslist = cf.options(mutual[k])
            number = re.findall(r'number\=(\d)',line)
            if keys[0] in keyslist and cf.get(mutual[k],keys[0])[1] == number[0]:
	            i = i + 1
                    line = lines[i]
                    while '</list>' not in line:
                            tmplist.append(line)
                            i = i+1
                            line = lines[i]
                    for tl in tmplist:
                            tmpline = tl
                            tmpline = tmpline.replace('{ccc}',pair)
			    tmpline = tmpline.replace('{ccc1}',pair)
                            tmpline = tmpline.replace('{c}',pair[0])
		            tmpline = tmpline.replace('{c1}',pair[0])
                            listwrite.append(tmpline)
                    tmplist = []
	    elif keys[0] in keyslist and (keys[0] == "mapping" or keys[0] == "aux_files" or keys[0] == "aux_files_24hr" or keys[0] == "aux_files_3hr" or keys[0] == "init_map" or keys[0] == "flux_albedos" or keys[0] == "call_map" or keys[0] == "call_mrg" or keys[0] == "model_setup" or keys[0] == "init_fractions" or keys[0] == "set_fractions" or keys[0] == "model_run" or keys[0] == "model_cpl" or keys[0] == "prognostic_expect_others" or keys[0] == "model_prep" or keys[0] == "budget_fractions"):
		    if keys[0] == "mapping":
			    print type(cf.get(mutual[k],keys[0]))
	    else:
		    i = i + 1
		    line = lines[i]
		    while '</list>' not in line:
			    i = i + 1
			    line = lines[i]
        tmpline = line
        tmpline = tmpline.replace('{ccc}',pair)
        tmpline = tmpline.replace('{ccc1}',pair)
        tmpline = tmpline.replace('{c}',pair[0])
        tmpline = tmpline.replace('{c1}',pair[0])
        if flag:
	    if keytemp and keytemp[0] == "mapping":
		print type(cf.get(mutual[k],keytemp[0]))
            listwrite.append(tmpline)
        i = i + 1
    desfile = "mrg_x2"+pair[0]+"_mct.F90"
    des = open(desfile,'w')
    des.writelines(listwrite)
    des.close()
