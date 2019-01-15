#!/usr/bin/python
import re
import ConfigParser
config_file_path="conf.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)
srcfile = eval(cf.get("create_newcase","srcfile"))
desfile = eval(cf.get("create_newcase","desfile"))

allcomp = eval(cf.get("create_newcase","allcomp"))
comp = eval(cf.get("create_newcase","comp"))
everycomp = eval(cf.get("create_newcase","everycomp"))
onecomp = eval(cf.get("create_newcase","onecomp"))
alcomp = eval(cf.get("create_newcase","alcomp"))
everycompdir = eval(cf.get("create_newcase","everycompdir"))
dcomp = eval(cf.get("create_newcase","dcomp"))
compmodel = eval(cf.get("create_newcase","compmodel"))
gridequal = eval(cf.get("create_newcase","gridequal"))

src = open(srcfile,'r+')
lines = src.readlines()
src.close()

i = 0
listwrite=[]
while i < len(lines):
    line = lines[i]
    if '{list}' in line:
        key = re.findall(r'key\=\"(.*?)\"',line)
        mylist = vars()[key[0]]
        i = i+1
        line = lines[i]
        for c in mylist:
            if isinstance(c,list):
                tmplist = line[1:].replace('{comp1}',c[0])
                tmplist = tmplist.replace('{comp2}',c[1])
            else:
		tmplist = line[1:].replace('{'+key[0]+'}',c)
	    if c == mylist[len(mylist)-1]:
		if tmplist[len(tmplist)-1] == ',': 
			tmplist = tmplist.replace(',','')
		tmplist = tmplist.replace('||','')
            listwrite.append(" "+tmplist)
    else:
        listwrite.append(line)
    i = i+1
des = open(desfile,'w')
des.writelines(listwrite)
des.close()

