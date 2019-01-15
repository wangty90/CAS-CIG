#!/usr/bin/python
import re
import ConfigParser

# Read all common component models in common.xml
config_file_path="common.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)

atm = eval(cf.get("common","atm"))
lnd = eval(cf.get("common","lnd"))
ice = eval(cf.get("common","ice"))
ocn = eval(cf.get("common","ocn"))
glc = eval(cf.get("common","glc"))

# Read all other component models in common.xml
compmodel=cf.options("other")
src = open("common.xml",'r+')
lines = src.readlines()
src.close()

i = 7
j = 0
listwrite=[]
model=[]
while i < len(lines):
        line = lines[i]
        pos = line.index('=')
        tmpline = line[0:pos-1]
        model.append(tmpline)
#       print model[j]
        pos = line.index('[')
        tmpline = line[pos+2:len(line)-3]
        compmodel[j].append(tmpline.split("','"))
#       print compmodel[0]
        i = i+1
        j = j+1

# Print necessary variables on conf.xml
listwrite.append("[create_newcase]")
listwrite.append("srcfile=\"create_newcase_template\"")
listwrite.append("desfile=\"create_newcase\"")
line = "allcomp = ['ATM','LND','ICE','OCN','CPL','GLC'"
for c in model:
        line = line+",'"+c.upper()+"'"
line = line+']'
listwrite.append(line)
line = "comp = ['ATM','LND','ICE','OCN','CPL','GLC']"
listwrite.append(line)
line = "everycomp = ["
for c in atm:
        line = line+",'"+c+"'"
for c in lnd:
        line = line+",'"+c+"'"
for c in ice:
        line = line+",'"+c+"'"
for c in ocn:
        line = line+",'"+c+"'"
for c in glc:
        line = line+",'"+c+"'"
for c in compmodel:
        for cc in c:
                line = line+",'"+cc[0]+"'"
line = line+']'
listwrite.append(line)
line = "onecomp = ['atm','lnd','ocn','ice','glc']"
listwrite.append(line)
line = "alcomp = ['ATM','LND','ICE','OCN','GLC'"
for c in model:
        line = line+",'"+c.upper()+"'"
line = line+']'
listwrite.append(line)
line = "everycompdir = ["
for c in atm:
        line = line+"['"+c+"','atm/"+c+'/bld/'+c+".cpl7.template'],"
for c in lnd:
        line = line+"['"+c+"','lnd/"+c+'/bld/'+c+".cpl7.template'],"
for c in ice:
        line = line+"['"+c+"','ice/"+c+'/bld/'+c+".cpl7.template'],"
for c in ocn:
        line = line+"['"+c+"','ocn/"+c+'/bld/'+c+".cpl7.template'],"
for c in glc:
        line = line+"['"+c+"','glc/"+c+'/bld/'+c+".cpl7.template'],"
i =0
for c in model:
        for cc in compmodel[i]:
                line = line+"['"+cc[0]+"','"+c+"/"+cc[0]+'/bld/'+cc[0]+".cpl7.template'],"
        i = i+1
line = line.rstrip(',')
line = line+']'
print line
listwrite.append(line)

