#!/usr/bin/python
import re
import ConfigParser
import types
import commands

config_file_path="../common.xml"
model = ConfigParser.ConfigParser()
model.read(config_file_path)

#common = eval(model.options("common"))
common = model.items("common")
other = model.items("other")

config_file_path="conf_mapping.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)

#print 'options:',common
#print 'options1:',common[0][0]
#mutual = {}
#count = 0
#for k,value1 in common:
#	for j,value2 in common:
#		if k <> j:
#			mutual[count] = k+j
#			count = count+1
#print mutual
#for k ,value in other:
#	print value
#	value1 = eval(value[1:len(value)-1])
#	print value1
#	if type(value1) is types.StringType:
#		mutual[count] = value1
#		count = count+1
#	else:
#		for j in value1:
#			print j
#			mutual[count] = j
#			count = count+1
#print mutual
mutual1 = cf.sections()
count = 1
mutual = {}
for k in mutual1:
	mutual[count] = k
	count = count+1

mutual_value = {}
for k in mutual:
	mutual_value[k] = cf.items(mutual[k])

srcfile = "templates/map_ccc1ccc2_mct.F90"
src = open(srcfile,'r')
lines = src.readlines()
src.close()

for k in mutual:
    listwrite = []
    i = 0
    pair = [mutual[k][0:3],mutual[k][3:6]]
    stack = 0
#    print mutual[k]
    while i < len(lines):
#	print i
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
#	    print i
            tmplist = []
            keys = re.findall(r'key\=\"(.*?)\"',line)
#	    if keys[0] == "public_interface":
#		print i
            keyslist = cf.options(mutual[k])
	    number = re.findall(r'number\=(\d)',line)
            if keys[0] in keyslist and cf.get(mutual[k],keys[0])[1] == number[0]:
#		print i
#            if keys in keyslist and number in cf.get(mutual[k],keys):
#		print "jinlaile",i
		stack = stack+1
	    	while stack > 0:
                	i = i + 1
                	line = lines[i]
			while '</list>' not in line and '<list>' not in line and '{list}' not in line:
	                    	tmplist.append(line)     
				i = i+1
				line = lines[i]
	                for tl in tmplist:
        	               	tmpline = tl
                	       	tmpline = tmpline.replace('{ccc1}',pair[0]).replace('{ccc2}',pair[1])
                        	tmpline = tmpline.replace('{c1}',pair[0][0]).replace('{c2}',pair[1][0])
                        	listwrite.append(tmpline)
			tmplist = []	
			if '</list>' in line:
				stack = stack-1
			elif '<list>' in line:
#				print i
   		                keys = re.findall(r'key\=\"(.*?)\"',line)
           			keyslist = cf.options(mutual[k])
		                number = re.findall(r'number\=(\d)',line)
            			if keys[0] in keyslist and cf.get(mutual[k],keys[0])[1] == number[0]:
					stack = stack+1
				else:
					stack1 = stack+1
					while stack1>stack:
		                                i = i+1
        		                        line = lines[i]
                        		        if '<list>' in line:
                                        		stack1 = stack1+1
                                		elif '</list>' in line:
                                        		stack1 = stack1-1
                                		elif '{list}' in line:
                                        		i = i+2			
			elif '{list}' in line:
		                keys = re.findall(r'key\=\"(.*?)\"',line)
            			keyslist = cf.options(mutual[k])
                                number = re.findall(r'number\=(\d)',line)
            			if keys[0] in keyslist and cf.get(mutual[k],keys[0])[1] == number[0]:
                			line = lines[i+1][1:]
				        tmpline = line
				        tmpline = tmpline.replace('{ccc1}',pair[0]).replace('{ccc2}',pair[1])
				        tmpline = tmpline.replace('{c1}',pair[0][0]).replace('{c2}',pair[1][0])
					listwrite.append(tmpline)
            			else:
                			line = lines[i+2][1:]
            			i = i + 2
				line = lines[i]
	    else:
#		print "no"
#		print keys
#		print stack
		stack1 = stack+1
		while stack1 > stack:
			i = i+1
			line = lines[i]
			if '<list>' in line:
				stack1 = stack1+1
			elif '</list>' in line:
				stack1 = stack1-1
			elif '{list}' in line:
				i = i+2
				line = lines[i][1:]		
#		print i
            i = i + 1
            line = lines[i]
        tmpline = line
        tmpline = tmpline.replace('{ccc1}',pair[0]).replace('{ccc2}',pair[1])
        tmpline = tmpline.replace('{c1}',pair[0][0]).replace('{c2}',pair[1][0])
        if flag:
            listwrite.append(tmpline)
        i = i + 1 
    desfile = "map_"+pair[0]+pair[1]+"_mct.F90"
    des = open(desfile,'w')
    des.writelines(listwrite)
    des.close()
    
commands.getoutput('mv map*.F90 ../generation/drv/driver/.')

