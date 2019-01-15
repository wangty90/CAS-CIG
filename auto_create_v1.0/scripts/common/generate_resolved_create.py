#!/usr/bin/python
import re
import ConfigParser
config_file_path="conf.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)
srcfile = eval(cf.get("generate_resolved","srcfile"))
desfile = eval(cf.get("generate_resolved","desfile"))

alcomp = eval(cf.get("generate_resolved","alcomp"))

src = open(srcfile,'r+')
lines = src.readlines()
src.close()

i = 0
listwrite=[]
while i < len(lines):
    line = lines[i]
    if '{list}' in line:
        key = re.findall(r'key\=\"(.*?)\"',line)
#	print key
        mylist = vars()[key[0]]
#	print mylist
        i = i+1
        line = lines[i]
        for c in mylist:
	    if c[2] == 0: 
		tmplist = "if ($COMP_{comp1} != 'none') then"
		tmplist = tmplist.replace('{comp1}',c[0])
	    	listwrite.append(tmplist+'\n')
            tmplist = line[1:].replace('{comp1}',c[0])
            tmplist = tmplist.replace('{comp2}',c[1])
            tmplist = tmplist.replace('{comp4}',c[3])
            listwrite.append(tmplist)
	    if c[2] == 0: 
		tmplist = "endif\n"
	        listwrite.append(tmplist)
    else:
        listwrite.append(line)
    i = i+1
des = open(desfile,'w')
des.writelines(listwrite)
des.close()

