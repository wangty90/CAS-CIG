#!/usr/bin/python
import re
import ConfigParser
import sys
try: 
  import xml.etree.cElementTree as ET 
except ImportError: 
  import xml.etree.ElementTree as ET 

#Read existant compsets & create new compset
tree = ET.parse("config_compsets.xml")
root = tree.getroot()

tag = 0
for compset in root.findall('compset'):
  name = compset.get('NAME')
  shortname = compset.get('SHORTNAME')
  if sys.argv[1]==name or sys.argv[1]==shortname:
    tag =1

if tag==0:
  print "new compset creating"
  newcompset = open('set.xml').read()
  f=open('config_compsets.xml','r+')
  while newcompset.find('<compset')!= -1:
    front = newcompset.find('<compset')
    rear = newcompset.find('/>')
    tmp = newcompset[front:rear+2]
#2016.8.30
    newcompset = newcompset[rear+1:]
    f.seek(-19,2)
    f.write(tmp)
  f.seek(0,2)
  f.write('\n\n</config_compset> \n')
  f.close()
  
#Read existant grids & create new grid
tree = ET.parse("config_grid.xml")
root = tree.getroot()

tag = 0
for compset in root.findall('horiz_grid'):
  grid = compset.get('GRID')
  shortname = compset.get('SHORTNAME')
  if sys.argv[2]==grid or sys.argv[2]==shortname:
    tag =1

if tag==0:
  print "new grid creating"

#Read existant machines & create new machine
tree = ET.parse("config_machines.xml")
root = tree.getroot()

tag = 0
for compset in root.findall('machine'):
  mach = compset.get('MACH')
  if sys.argv[3]==mach:
    tag =1

if tag==0:
  print "new machine creating"

#Read existant definitions & create new values
tree = ET.parse("config_definition.xml")
root = tree.getroot()

tag = 0
for compset in root.findall('definition')
  mach = compset.get('definition')
  if sys.argv[4]==definition:
    tag =1

if tag==0:
  print "new values creating"

