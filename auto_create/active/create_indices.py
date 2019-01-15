#!/usr/bin/python
import re
import ConfigParser
import types
import os

srcfile = "ccc_CplIndices.F90"
#modelccc = ['ocn']
modconf = "../common.xml"
modcf = ConfigParser.ConfigParser()
modcf.read(modconf)
common = modcf.items("common")
other = modcf.items("other")
modccc=[]
typccc=dict()
for types in common:
        tt = types[0]
        ty = eval(types[1])
        for name in ty:
                modccc.append(name)
                typccc[name]=tt

src = open(srcfile,'r')
lines = src.readlines()
src.close()

for ccc in modccc:
    xmlconf = ccc+"_indices_conf.xml"
    xmlcf = ConfigParser.ConfigParser()
    xmlcf.read(xmlconf)
  
    xmlfile = ccc+"_indices_conf.F90"
    xmlsrc = open(xmlfile,'r')
    xmllines = xmlsrc.readlines()
    xmlsrc.close()

    listwrite = []
    i = 0
    while i < len(lines):
        if '{list}' in lines[i]:
            key = re.findall(r'key=\"(.*?)\"',lines[i])
            if key[0] in xmlcf.options(xmlcf.sections()[0]):
#		print key[0]
		if key[0] == "indices_set":
			number = re.findall(r'number\=(\d)',lines[i])
			if xmlcf.get(xmlcf.sections()[0],key[0])[1] == number[0]:
				line = lines[i+1][1:]
		                listwrite.append(line)
			i = i+2
		elif key[0] == "var_model_to_cpl":
			values = xmlcf.get(xmlcf.sections()[0],key[0]).strip('[').strip(']').split(',')
			x2c_values = []
			for value in values:
				define = value.strip('"').split(':')
				listwrite.append(define[1]+" :: "+define[0]+"\n")
				x2c_values.append(define[0])
			i = i+1
		elif key[0] == "var_cpl_to_model":
                        values = xmlcf.get(xmlcf.sections()[0],key[0]).strip('[').strip(']').split(',')
			c2x_values = []
                        for value in values:
                                define = value.strip('"').split(':')
                                listwrite.append(define[1]+" :: "+define[0]+"\n")
                                c2x_values.append(define[0])
                        i = i+1
		elif key[0] == "init_var_x2c":
			no_x2c_values = xmlcf.get(xmlcf.sections()[0],key[0]).strip('[').strip(']').split(',') 
			for value in x2c_values:
				if value not in no_x2c_values:
					listwrite.append(value+" = mct_avect_indexra("+typccc[ccc][0]+"2x,'"+value.split(typccc[ccc][0]+"2x")[1][1:]+"')\n")
			i = i+1
		elif key[0] == "init_var_c2x":
                        no_c2x_values = xmlcf.get(xmlcf.sections()[0],key[0]).strip('[').strip(']').split(',')
                        for value in c2x_values:
                                if value not in no_c2x_values:
                                        listwrite.append(value+" = mct_avect_indexra(x2"+typccc[ccc][0]+",'"+value.split("x2"+typccc[ccc][0])[1][1:]+"')\n")
                        i = i+1
	    else:
		i = i+1 	
        elif 'xmlinsert' in lines[i]:
            id = re.findall(r'xmlinsert\(list\[(.*?)\]\)',lines[i])[0]
            j = 0
            while j < len(xmllines):
                while j<len(xmllines) and "[start "+id+"]" not in xmllines[j]:
                    j = j + 1
                j = j + 1
                if j>=len(xmllines):
                    break
                while "[end "+id+"]" not in xmllines[j]:
                    listwrite.append(xmllines[j])
                    j = j + 1
                break
            i = i+1
	else :
        	flag = 1
        	line = lines[i]
        	tmpline = line
        	tmpline = tmpline.replace('{ccc}',ccc).replace('{c}',typccc[ccc][0])
        	if flag:
            		listwrite.append(tmpline)
        	i = i + 1
    isExists=os.path.exists("./"+ccc)
    if not isExists:
        os.makedirs("./"+ccc)
    desfile = "./"+ccc+"/"+ccc+"_CplIndices.F90"
    des = open(desfile,'w')
    des.writelines(listwrite)
    des.close()

