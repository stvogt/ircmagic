#! /usr/bin/python 

import sys, re, os, commands, shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.integrate import cumtrapz
sys.path.append("/home/stvogt/bin")
import cubman

header =  '''
-----------------------------------------------------------------------
          IrcMagic: An Open-Source IRC proccesing suit
                            IrcMagic 0.1

                            S. Vogt-Geisse
-----------------------------------------------------------------------

'''
def read_lines(outfile):
  if not os.path.exists(outfile):
    print outfile + " does not exist there -> " + os.getcwd()
    sys.exit(2)
  
  outputfile = open(outfile,"r")
  outputLines = outputfile.readlines()
  outputfile.close()
  return outputLines

def simple_plot(x, y):
  plt.plot(x, y, 'ro')
  #plt.axhline(y=0, color='red')
  plt.show()

def fit(x_coord, y_coord,rank):
  fx_coord = map(float, x_coord)
  fy_coord = map(float, y_coord)
  x = np.array(fx_coord)
  y = np.array(fy_coord)
  z = np.polyfit(x,y,rank)

  # creat polynomial function
  p = np.poly1d(z)

  #plot
  xp = np.linspace(frx_coord[0], frx_coord[-1], 100)
  plt.plot(x, y, '.', xp, p(xp), '-')
  plt.show()


def neg_derivative(x_coord, y_coord):
  fx_coord = map(float, x_coord)
  fy_coord = map(float, y_coord)
  x = np.array(fx_coord)
  y = np.array(fy_coord)
  f = -np.diff(y)/np.diff(x)
  return (fx_coord[1:], f)

def integrate(x_coord, y_coord):
  fx_coord = map(float, x_coord)
  fy_coord = map(float, y_coord)
  x = np.array(fx_coord)
  y = np.array(fy_coord)
  work = simps(y, x)
  #work = np.trapz(y, x)
  return work

def unique(List):
    output = []
    for x in List:
      if not x in output:
        output.append(x)
    return output

def print_lists(list1,list2,list1_name,list2_name,filename):
  f = open(filename+'.dat', 'w')
  f.write(header)
  f.write("Generated output from print_lists routine\n\n")
  f.write('%2s %20s %20s' % (" ",  list1_name, list2_name)+'\n' )
  for j in range(0,len(list1)):
      f.write('%2d %20f %20f' % (j+1, round(float(list1[j]),5), round(float(list2[j]),5))+'\n')
  

def general_print(filename, **kwargs):
  f = open(filename, 'w')
  f.write(header)
  f.write("Generated output from print routine\n\n")
  f.write('%5s' % (" "))
  for key in kwargs:
    if key == "Reaction Coordinate":
      iter_range = len(kwargs[key])
      f.write('%20s' % (key))
  for key in kwargs:
    if key != "Reaction Coordinate":
      f.write('%15s' % (key))
  f.write('\n')
  #f.write('%2d' % (j+1)
  #print kwargs
  for j in range(0,iter_range):
    f.write('%2d' % (j+1)) 
    for key in kwargs:
      
      f.write('%20f' % (round(float(kwargs[key][j]),5)))
    f.write('\n')

def num_valence_orbs(atom_list):
  val_dict= {
      "1"  : 1,
      "2"  : 2,
      "3"  : 1,
      "4"  : 2,
      "5"  : 3,
      "6"  : 4,
      "7"  : 5,
      "8"  : 6,
      "9"  : 7,
      "10" : 8,
      "11" : 1, 
      "12" : 2,
      "13" : 3,
      "14" : 4,
      "15" : 5,
      "16" : 6,
      "17" : 7,
      "18" : 8,
      "19" : 1,
      "20" : 2,
      }
  val_elec = 0
  for atom in atom_list:
    val_elec = val_elec + val_dict[atom]
  val_orbs = val_elec/2
  return val_orbs

def general_plot(plotname, ylabel, works, Zeroline, Show, **kwargs):
  bullets = ['bo',  'rs', 'k^', 'g-', 'm-.','bo',  'rs', 'k^', 'g-', 'm-.']
  count = 0
  print "I am using this operation.py"
  for key in kwargs:
    if key == "Reaction Coordinate":
      x = kwargs[key]
      plt.xlabel(key)
  for key in kwargs:
    if key != "Reaction Coordinate":
      plt.plot(x, kwargs[key], bullets[count])
      count = count + 1
  if Zeroline:
    plt.axhline(y=0, color = 'red')
  plt.axvline(x=works[4], color = 'black', linestyle = "dashed" )
  plt.axvline(x=works[5], color = 'black', linestyle = "dashed" )
  plt.ylabel(ylabel)
  plt.savefig(plotname)
  if Show:
    plt.show()


def general_plot_orbs(plotname, ylabel, orblist, zeroline=False, NoSave = False, **kwargs):
  bullets = ['bo',  'rs', 'k^', 'g-', 'm-.','bo',  'rs', 'k^', 'g-', 'm-.']
  keys = []
  count = 0
  occupied = False
  virtual = False
  for key in kwargs:
    if key == "Reaction Coordinate":
      x = kwargs[key]
      plt.xlabel(key)
  for key in kwargs:
    if key != "Reaction Coordinate":
      keys.append(key)

  if all(float(kwargs[key][0]) < 0 for key in kwargs):
    keys.sort(key=int, reverse=True)
    occupied = True
  else:
    keys.sort(key=int)
    virtual = True
    
  for orbNum in orblist:
    if occupied:
      key = str(int(1 + float(keys[0]) - float(orbNum)))
    elif virtual:
      key = str(int(float(orbNum)))
    else:
      print "No orbitals found, check your orbital list"
      sys.exit(1)
    plt.plot(x, kwargs[key], bullets[count])
    if zeroline:
      plt.axhline(y=0, color='black')
    plt.ylabel(ylabel)
    plt.savefig(plotname)
    #plt.savefig(plotname.split('.')[0]+'_'+str(count1)+plotname.split('.')[1])
    count = count + 1
  if not NoSave:
    plt.show()

def general_plot_prop(plotname, ylabel, proplist, works, Zeroline, Show, **kwargs):
  doc = """ 
        Takes as input a diccionary, where the keys are the property to be ploted,
        and the value is a list with the property  value for each point of the reaction
        coordinate.
  """
  bullets = ['bo',  'rs', 'k^', 'g-', 'm-.','bo',  'rs', 'k^', 'g-', 'm-.']
  keys = []
  count = 0
  occupied = False
  virtual = False
  for key in kwargs:
    if key == "Reaction Coordinate":
      x = kwargs[key]
      plt.xlabel(key)

  for key in proplist:
    if key != "Reaction Coordinate":
      plt.plot(x, kwargs[key], bullets[count])
      count = count + 1
  if Zeroline:
    plt.axhline(y=0, color='black')
  plt.ylabel(ylabel)
  plt.axvline(x=works[4], color = 'black', linestyle = "dashed" )
  plt.axvline(x=works[5], color = 'black', linestyle = "dashed" )
  plt.savefig(plotname)
  if Show:
    plt.show()
  
def atom_label(atom_num_list):
  atom_dict= {
      "1"  : "H",
      "2"  : "He",
      "3"  : "Li",
      "4"  : "Be",
      "5"  : "B",
      "6"  : "C",
      "7"  : "N",
      "8"  : "O",
      "9"  : "F",
      "10" : "Ne",
      "11" : "Na", 
      "12" : "Mg",
      "13" : "Al",
      "14" : "Si",
      "15" : "P",
      "16" : "S",
      "17" : "Cl",
      "18" : "Ar",
      }
  atom_label = []
  for atom_num in atom_num_list:
    atom_label.append(atom_dict[atom_num])
  return atom_label

def render_image(cubefile,jmolfile):
  print "Generating image corresponding to file:  "+cubefile
  lines=read_lines(jmolfile)
  f = open("jmol.spt", 'w')
  for line in lines:
    if "%full_file_dir" in line:
      f.write(line.replace('%full_file_dir', os.getcwd()+"/"+cubefile))
    elif "%image_file" in line:
      f.write(line.replace('%image_file', cubefile.split(".")[0]+".png"))
    else:
      f.write(line)
  f.close()
  if not os.path.isdir("orb_images"):
    os.makedirs("orb_images")
  os.chdir("orb_images")
  input0 = "/home/stvogt/.jmol/jmol.sh -n -s ../jmol.spt"
  print os.getcwd()
  print input0
  job0 = commands.getstatusoutput(input0)
  os.chdir('../')

#def generate_cube(orbital,fchk,cube_output,sq=False):
#  cub = "orbital_"+str(orbital)+"_"+cube_output
#  input1="cubegen 0 MO="+str(orbital)+" "+fchk+" "+cub+" 0 h" 
#  print input1
#  job1 = commands.getstatusoutput(input1)
#  if sq:
#    sq_cub = cubman.sq_cube(cub)
#    os.system("rm "+cub)
#    return sq_cub
#  else:
#    return cub
#

def fchk_gen(chk_folder):
    for file_ in os.listdir(chk_folder):
      print "Generating formated checkpoint file for:  "+file_
      input0 = "formchk "+chk_folder+'/'+file_ 
      print input0
      job0 = commands.getstatusoutput(input0)

def cube_files(orb_range, atom_list):
  if os.path.isdir("CHK"):
    f = open("cube.log", "w")
    f.write("Log file for cube generation\n\n")
    if not os.path.isdir("cubes"):
      os.makedirs("cubes")
    #fchk_gen("CHK/")
    for file_ in os.listdir("CHK"):
      print "generating checkpoint file for:  "+file_
      input0 = "formchk "+"CHK/"+file_ 
      print input0
      job0 = commands.getstatusoutput(input0)
      os.chdir("cubes")
      print "Entering this directory --> "+os.getcwd()
      fchk = "../CHK/"+file_.split('.')[0]+".fchk"
      cube_output = file_.split(".")[0]+".cub"
      generate_cubes(orb_range,fchk,cube_output)
      print "ordering orbital cubes into directories"
      for i in orb_range: 
        if not os.path.isdir("orbitals_"+str(i)):
          os.makedirs("orbitals_"+str(i))
        os.chdir("orbitals_"+str(i))
        original = "../orbitals_"+str(i)+"_"+file_.split(".")[0]+".cub_sq" 
        destiny  = "orbitals_"+str(i)+"_"+file_.split(".")[0]+".cub"
        if os.path.isfile(original):
          shutil.move(original, destiny)
        else:
          print "Cube file "+original+" not found, please check your checkpoint file"
          f.write("Cube file "+original+" not found, please check your checkpoint file\n")
        #render_image(destiny,"/home/stvogt/bin/jmol.scripts")
        os.chdir("../")
      #if sum_orbs:
      #  if not os.path.isdir("orb_sum"):
      #    os.makedirs("orb_sum")
      #  os.chdir("orb_sum")
      #  original = "../sum_of_all_cubes_"+cube_output+"_sq"
      #  destiny  = "sum_of_all_cubes_"+cube_output
      #  if os.path.isfile(original):
      #    shutil.move(original, destiny)
      #  else:
      #    print "Cube file "+original+" not found, please check!"
      #    f.write("Cube file "+original+" not found, please check!\n")
      #  render_image(destiny,"/home/stvogt/bin/jmol.scripts")
      #  os.chdir("../")
      os.system("rm *.cub")
      os.chdir("../")
    f.close()  

  else:
    print "No checkpoint folder CHK!"
    sys.exit(1)
#
#def generate_cubes(orbitals,fchk,cube_output,sq=False):
#  count = 0
#  orb_names = []
#  for i in orbitals: 
#		print "Genrating cube for orbital number: "+str(i)
#		cub_i = "orbitals_"+str(i)+"_"+cube_output
#		orb_name.append(cub_i)
#		input1="cubegen 0 MO="+str(i)+" "+fchk+" "+cub_i+" 0 h" 
#		print input1
#		job1 = commands.getstatusoutput(input1)
#	#	if sq:
#	#		sq_cub = cubman.sq_cube(cub_i)
#	#		os.system("rm "+cub_i)
#    #os.system("ulimit -s unlimited")
    #print job1

    #if sum_orbs:
    #  cub0 = "sum_of_all_cubes_"+cube_output
    #  if count == 0:
    #    print "Generating the cube of the sum accumulation"
    #    input2="cubegen 0 MO="+str(i)+" "+fchk+" "+cub0+" 0 h" 
    #    job2 = commands.getstatusoutput(input2)
    #    cubman.sq_cube(cub0)
    #    print input2
    #  else:
    #    print "Summing cubes "+cub0+" with "+cub_i
    #    cubman.sum_total(cub0+"_sq", cub_i+"_sq")
    #  count = count +1


#def valence_density():
#  pass

#def general_plot(plotname, ylabel, works, zeroline=False, blocks = 3, **kwargs):
#  count = 0
#  count1 = 0
#  bullets = ['bo',  'rs', 'k^', 'g-', 'm-.','bo',  'rs', 'k^', 'g-', 'm-.']
#  for key in kwargs:
#    if key == "Reaction Coordinate":
#      x = kwargs[key]
#      plt.xlabel(key)
#  for key in kwargs:
#    if key != "Reaction Coordinate":
#      plt.plot(x, kwargs[key], bullets[count])
#      count = count + 1
#      if count % blocks == 0:
#        count1 = count1+1
#        if zeroline:
#          plt.axhline(y=0, color='black')
#        plt.ylabel(ylabel)
#        plt.savefig(plotname.split('.')[0]+'_'+str(count1))
#        #plt.savefig(plotname.split('.')[0]+'_'+str(count1)+plotname.split('.')[1])
#        plt.show()
#  if len(kwargs) < 3:
#    if zeroline:
#      plt.axhline(y=0, color='black')
#    plt.ylabel(ylabel)
#    plt.savefig(plotname)
#    plt.show()
