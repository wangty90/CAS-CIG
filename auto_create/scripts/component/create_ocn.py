#!/usr/bin/python
import re
import ConfigParser
config_file_path="conf.xml"
cf = ConfigParser.ConfigParser()
cf.read(config_file_path)
srcfile = eval(cf.get("create_newcase","srcfile"))
desfile = eval(cf.get("create_newcase","desfile"))

