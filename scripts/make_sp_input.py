#! /usr/bin/python 

import sys, re, os
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import operations as op
import singlePoint as sp
from optparse import OptionParser


def extract_xyz_irc(startstring,endstring,lines,lineNum,mult):
  xyz_info = []
  atoms_num = []
  xyz = []
  x_coords = []
  y_coords = []
  z_coords = []
  if startstring in lines[lineNum]:  
    for lineNum1 in range(lineNum,len(lines)):
      if not endstring in lines[lineNum1]: 
        p = re.search('\d+\s*-?\d+.\d+\s*-?\d+.\d+', str(lines[lineNum1]))
        if p: 
          if len(lines[lineNum1].split()) == 5:
            atoms_num.append(lines[lineNum1].split()[1]) # Atom number
            x_coords.append(lines[lineNum1].split()[2])
            y_coords.append(lines[lineNum1].split()[3])
            z_coords.append(lines[lineNum1].split()[4])
          elif len(lines[lineNum1].split()) == 6:
            atoms_num.append(lines[lineNum1].split()[1]) # Atom number
            x_coords.append(lines[lineNum1].split()[3])
            y_coords.append(lines[lineNum1].split()[4])
            z_coords.append(lines[lineNum1].split()[5])
      else:
        if len(lines[lineNum1].split()) > 4:
          rx_coord = float(lines[lineNum1].split()[8])*mult
        else: 
          rx_coord = 0.00000
        xyz.append(atoms_num)
        xyz.append(x_coords)
        xyz.append(y_coords)
        xyz.append(z_coords)
        xyz_info.append(rx_coord)
        xyz_info.append(xyz)
        break
    return xyz_info

def get_xyz(lines):
  structures = []
  for lineNum in range(0,len(lines)):
    if not  "Distance matrix (angstroms):" in lines[lineNum]:
      if extract_xyz_irc("Input orientation:", "Distance matrix (angstroms):", lines,lineNum,-1):
        structures.append(extract_xyz_irc("Input orientation:", "Distance matrix (angstroms):", lines, lineNum,-1))
    else:
      break
  for lineNum in range(0,len(lines)):
    if not "Forward" in lines[lineNum]:
      if extract_xyz_irc("CURRENT STRUCTURE", "NET REACTION COORDINATE UP TO THIS POINT", lines,lineNum,-1):
        structures.append(extract_xyz_irc("CURRENT STRUCTURE", "NET REACTION COORDINATE UP TO THIS POINT", lines,lineNum,-1))
    else:
      for lineNum0 in range(lineNum,len(lines)):
        if extract_xyz_irc("CURRENT STRUCTURE", "NET REACTION COORDINATE UP TO THIS POINT", lines,lineNum0,1):
          structures.append(extract_xyz_irc("CURRENT STRUCTURE", "NET REACTION COORDINATE UP TO THIS POINT", lines,lineNum0,1))
      break
  structures.sort()
  return structures

parser = OptionParser() 
parser.add_option("-o", "--output_irc", dest="filenameIRC", help="File with IRC output (default: output.dat) Warning: Use join_irc.py to create the correct outputfile ", default = "output.dat") #, metavar="FILE")
parser.add_option("-i", "--input_sp", dest="input_sp", help="File with single points output (default: input_sp.dat)", default = "input_sp.dat") #, metavar="FILE")
parser.add_option("-p", "--numproc", dest="proc_num", help="Number of procesors (default: 4)", default="4") 
parser.add_option("-m", "--memory", dest="mem_num", help="Memory alloted for the computation (default: 2GB)", default="2GB") 
parser.add_option("-a", "--folder", dest="folder", help="Name of the Single Point folder (default: sp)", default = "sp") 
parser.add_option("-s", "--scf", dest="scf", help="The cut off criterium for scf convergence (default: tight)", default = "tight") 
parser.add_option("-c", "--charge", dest="charge", help="charge of the compuation (default: 0)", default = "0") 
parser.add_option("-n", "--multiplicity", dest="multi", help="multiplicity of the computation (default: 1)", default = "1") 
parser.add_option("-z", "--otheroption", dest="other", help="Any other option that should be added", default = " ") 

(options, args) = parser.parse_args()

input_irc = options.filenameIRC
input_sp = options.input_sp
proc_num = options.proc_num
mem_num = options.mem_num
sp_folder = options.folder
scf = options.scf
multi = options.multi
charge = options.charge
other = options.other

filelines = op.read_lines(input_irc)

theory = sp.get_level_of_theory(filelines)
basis = sp.get_basis(filelines)

if not os.path.isdir(sp_folder):
  os.makedirs(sp_folder)
os.chdir(sp_folder)

f = open(input_sp,'w') 

if not os.path.isdir("CHK"):
  os.makedirs("CHK")

xyz_list = get_xyz(filelines)

for strucNum in range(0,len(xyz_list)):
  f.write("%NProcShared="+proc_num+"\n")
  f.write("%Mem="+mem_num+"\n")
  f.write("%chk=CHK/sp_"+str(strucNum).zfill(3)+".chk\n")
  f.write("#P  "+theory+"/"+basis+" pop=NBORead  scf=("+scf+") int=Ultrafine\n\n")
  f.write("Single Point computation for reaction coordinate: "+str(xyz_list[strucNum][0])+"\n\n"+str(charge)+" "+str(multi)+"\n")
  for atmNum in range(0, len(xyz_list[strucNum][1][0])):
    f.write(xyz_list[strucNum][1][0][atmNum]+"    "+xyz_list[strucNum][1][1][atmNum]+"    "+xyz_list[strucNum][1][2][atmNum]+"    "+xyz_list[strucNum][1][3][atmNum]+"\n")
  f.write("\n$nbo bndidx $end\n")
  if strucNum != len(xyz_list)-1:
    #print str(strucNum) + "\t" + str(len(xyz_list))
    f.write("\n--Link1--\n")
  else: 
    print "Input file contaiing "+str(strucNum)+" inputs was written to --->  sp/"+input_sp
    f.write("\n\n")
f.close()


