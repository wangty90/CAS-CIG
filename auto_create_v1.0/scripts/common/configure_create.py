#!/usr/bin/python
import re
import ConfigParser
config_file_path="conf.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)
srcfile = eval(cf.get("configure","srcfile"))
desfile = eval(cf.get("configure","desfile"))

decomp = eval(cf.get("configure","decomp"))

src = open(srcfile,'r+')
lines = src.readlines()
src.close()

i = 0
listwrite=[]
while i < len(lines):
    line = lines[i]
    if '{list}' in line:
        for c in decomp:
	    tmplist = "         if ($COMP_{comp1} == '{comp2}') then"
	    tmplist = tmplist.replace('{comp1}',c[0])
	    tmplist = tmplist.replace('{comp2}',c[1])
            listwrite.append(tmplist+'\n')
	    tmplist = "            if (${comp3}_AUTO_DECOMP == 'true') then"
	    tmplist = tmplist.replace('{comp3}',c[2])
            listwrite.append(tmplist+'\n')
	    tmplist = "            set config = `./Tools/Templates/generate_{comp4}_decomp.pl -res ${comp1}_GRID -platform $OS -model {comp4} -nproc $NTASKS_{comp1} -thrds $NTHRDS_{comp1} -output all`"
            tmplist = tmplist.replace('{comp1}',c[0])
            tmplist = tmplist.replace('{comp2}',c[1])
            tmplist = tmplist.replace('{comp4}',c[3])
            listwrite.append(tmplist+'\n')
	    tmplist = "            if ($config[1] >= 0) then"
            listwrite.append(tmplist+'\n')
	    tmplist = "                ./xmlchange -file env_build.xml -id {comp3}_BLCKX      -val $config[3]"
	    tmplist = tmplist.replace('{comp3}',c[2])
            listwrite.append(tmplist+'\n')
            tmplist = "                ./xmlchange -file env_build.xml -id {comp3}_BLCKY      -val $config[4]"
            tmplist = tmplist.replace('{comp3}',c[2])
            listwrite.append(tmplist+'\n')
            tmplist = "                ./xmlchange -file env_build.xml -id {comp3}_MXBLCKS    -val $config[5]"
            tmplist = tmplist.replace('{comp3}',c[2])
            listwrite.append(tmplist+'\n')
            tmplist = "                ./xmlchange -file env_build.xml -id {comp3}_DECOMPTYPE -val $config[6]"
            tmplist = tmplist.replace('{comp3}',c[2])
            listwrite.append(tmplist+'\n')
	    tmplist = "            else"
            listwrite.append(tmplist+'\n')
	    tmplist = "               echo \"{comp4} decomp not set for ${comp1}_GRID on $NTASKS_{comp1} x $NTHRDS_{comp1} procs\""
            tmplist = tmplist.replace('{comp4}',c[3])
            tmplist = tmplist.replace('{comp1}',c[0])
            listwrite.append(tmplist+'\n')
            tmplist = "               exit -1"
            listwrite.append(tmplist+'\n')
            tmplist = "            endif"
            listwrite.append(tmplist+'\n')
            tmplist = "            endif"
            listwrite.append(tmplist+'\n')
            tmplist = "         endif"
            listwrite.append(tmplist+'\n')
    else:
        listwrite.append(line)
    i = i+1
des = open(desfile,'w')
des.writelines(listwrite)
des.close()

