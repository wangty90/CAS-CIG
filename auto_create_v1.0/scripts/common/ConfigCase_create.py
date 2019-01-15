#!/usr/bin/python
import re
import ConfigParser
config_file_path="conf.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)
srcfile = eval(cf.get("ConfigCase","srcfile"))
desfile = eval(cf.get("ConfigCase","desfile"))

optcomp = eval(cf.get("ConfigCase","optcomp"))
alcomp = eval(cf.get("ConfigCase","alcomp"))
allcomp = eval(cf.get("ConfigCase","allcomp"))
ocncomp = eval(cf.get("ConfigCase","ocncomp"))
ocnicecomp = eval(cf.get("ConfigCase","ocnicecomp"))

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
            if isinstance(c,list):
#		print c
                tmplist = line[1:].replace('{comp1}',c[0])
                tmplist = tmplist.replace('{comp2}',c[1])
            else:
		tmplist = line[1:].replace('{'+key[0]+'}',c)
            listwrite.append(" "+tmplist)
    else:
        listwrite.append(line)
    i = i+1
des = open(desfile,'w')
des.writelines(listwrite)
des.close()

