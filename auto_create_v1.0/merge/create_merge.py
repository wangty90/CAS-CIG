#!/usr/bin/python
import re
import ConfigParser
import types
import commands

srcfile = "merge/templates/mrg_x2{c}_mct.F90"
src = open(srcfile,'r')
lines = src.readlines()
src.close()

config_file_path="merge/conf_merge.xml"
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
        if '{list}' in line:
            keys = re.findall(r'key\=\"(.*?)\"',line)
            keyslist = cf.options(mutual[k])
            number = re.findall(r'number\=(\d)',line)
            if keys[0] in keyslist and cf.get(mutual[k],keys[0])[1] == number[0]:
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
            listwrite.append(tmpline)
        i = i + 1
    desfile = "mrg_x2"+pair[0]+"_mct.F90"
    des = open(desfile,'w')
    des.writelines(listwrite)
    des.close()

commands.getoutput('mv mrg*.F90 ./generation/drv/driver/.')
