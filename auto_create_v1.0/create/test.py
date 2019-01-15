#!/usr/bin/python
filename="ccsm_pre_init_example.F90"
CCC=['CPL','OCN','ATM','LND','ICE','GLC','WRF','GEA']

src = open(filename, "r+")
listwrite = []
clist = []
while True:
    line = src.readline()
    if '{lists}' in line:  
	#print line
	#listwrite.append(line)
	line = src.readline()
        #print line
	#listwrite.append(line)
	
	if 'CPL' in line:
	    clist = CCC[1:]
	else:
	    clist = CCC 
	
	for c in clist:
	    listwrite.append("   "+line[1:].replace('{CCC}',c).replace('{ccc}',c.lower()))
    elif line:
	listwrite.append(line)
    else:
	break
src.close()

fw = open("abc.F90",'w')
for listw in listwrite:
    fw.writelines(listw)
fw.close()


