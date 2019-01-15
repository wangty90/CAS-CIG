#!/usr/bin/python
import re
import ConfigParser
import types
import commands

config_file_path="./common.xml"
model = ConfigParser.ConfigParser()
model.read(config_file_path)

#common = eval(model.options("common"))
common = model.items("common")
other = model.items("other")

commands.getoutput('cp ./mapping/defaults/conf_mapping.xml.defaults ./mapping/conf_mapping.xml')
commands.getoutput('cp ./map/defaults/conf_map.xml.defaults ./map/conf_map.xml')
commands.getoutput('cp ./merge/defaults/conf_merge.xml.defaults ./merge/conf_merge.xml')
