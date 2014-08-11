#! /usr/bin/python

import sys, os, glob
import irclib
import singlePoint as sp
import operations as op
from pyDensity import Density

class Irc:
  def __init__(self, outfile):
    self.outfile = outfile
    self.outpath = os.path.abspath(outfile)
    self.lines = op.read_lines(outfile)

  def rxCoord(self):
    rx_coord = irclib.get_rxCoord(self.lines)  
    return rx_coord

  def atoms(self):
    return sp.get_xyz(self.lines)[0]

  def energy(self):
    energy = irclib.all_energies(self.lines)
    return {"Reaction Coordinate": self.rxCoord(), "Energy":energy}

  def force(self):
    coord = op.neg_derivative(self.rxCoord(),self.energy()['Energy'])[0]
    force = op.neg_derivative(self.rxCoord(),self.energy()['Energy'])[1]
    return {"Reaction Coordinate":coord , "Reaction Force":force}

  def ReactionWorks(self):
    y_in = self.force()["Reaction Force"]
    x_in = self.force()["Reaction Coordinate"]
    if float(x_in[0]) > 0:
      x = x_in[::-1] 
      y = y_in[::-1] 
    else: 
      x = x_in 
      y = y_in 
    max_force = max(y)
    min_force = min(y)

    for i in range(0,len(y)):
      if y[i] == min_force:
        w1_cut = i
        w1_rxcoord = x[w1_cut]
      elif y[i] == max_force:
        w3_cut = i
        w3_rxcoord = x[w3_cut]

    for i in range(5,len(y)):
      if y[i] > 0.0:
        w2_cut = i-1
        w2_rxcoord = x[w2_cut]
        break

    w1_y = y[:w1_cut+1]
    w1_x = x[:w1_cut+1]
    w2_y = y[w1_cut:w2_cut+1] ## Due to python list handling
    w2_x = x[w1_cut:w2_cut+1]
    w3_y = y[w2_cut:w3_cut+1]
    w3_x = x[w2_cut:w3_cut+1]
    w4_y = y[w3_cut:]
    w4_x = x[w3_cut:]

    w1 = round(-1.0*(op.integrate(w1_x,w1_y)),2)
    w2 = round(-1.0*(op.integrate(w2_x,w2_y)),2)
    w3 = round(-1.0*(op.integrate(w3_x,w3_y)),2)
    w4 = round(-1.0*(op.integrate(w4_x,w4_y)),2)
    return (w1,w2,w3,w4,w1_rxcoord,w3_rxcoord)

  def distances(self):
    print "----Distances---"
    sigma = []
    Sigma = {}
    all_dis = irclib.all_distances(self.lines)
    num_dist = len(all_dis[0][0])
    Sigma["Reaction Coordinate"] = self.rxCoord()
    for j in range(0, num_dist):
      for coord in range(0,len(all_dis)):
        sigma.append(all_dis[coord][1][j])
      Sigma[all_dis[0][0][j]] = sigma
      print all_dis[0][0][j]
      sigma = []
    return Sigma

  def angles(self):
    print "----Angles---"
    sigma = []
    Sigma = {}
    all_angles = irclib.all_angles(self.lines)
    num_ang = len(all_angles[0][0])
    Sigma["Reaction Coordinate"] = self.rxCoord()
    for j in range(0, num_ang):
      for coord in range(0,len(all_angles)):
        sigma.append(all_angles[coord][1][j])
      Sigma[all_angles[0][0][j]] = sigma
      print all_angles[0][0][j]
      sigma = []
    return Sigma

  def bondOrders(self):
    Bndo = {}
    bndo_list = []
    bondOrder =  irclib.all_bondOrders(self.lines)
    rx_coord = self.rxCoord()
    Bndo["Reaction Coordinate"] = rx_coord
    print "----------- Bonds ----------"
    for key in bondOrder[0]:
      print key
      for coord_num in range(0,len(rx_coord)):
        dicts = bondOrder[coord_num]
        bndo_list.append(dicts[key])
      Bndo[key] = bndo_list  
      bndo_list = []
    return Bndo

  def natCharges(self):
    Charges = {}
    charges_list = []
    charges =  irclib.all_NatCharges(self.lines)
    rx_coord = self.rxCoord()
    Charges["Reaction Coordinate"] = rx_coord
    print "------------Charges-------------"
    for key in charges[0]:
      print key
      for coord_num in range(0,len(rx_coord)):
        dicts = charges[coord_num]
        charges_list.append(dicts[key])
      Charges[key] = charges_list  
      charges_list = []
    return Charges

  def bondOrbital(self):
    irclib.all_bondOrbitals(self.lines)

  def occ_orbtitals(self):
    all_orbs = irclib.all_orbs(self.lines)
    val_orbs = op.num_valence_orbs(self.atoms())
    epsilon = []
    Epsilon = {}
    num_orbs = len(all_orbs[0][0])
    Epsilon["Reaction Coordinate"] = self.rxCoord()
    end = num_orbs
    start = num_orbs - val_orbs
    for j in range(start,end):
      for coord in range(0,len(all_orbs)):
        epsilon.append(all_orbs[coord][0][j])
      Epsilon[str(j)] = epsilon
      epsilon = []
    return Epsilon

  def virt_orbtitals(self):
    all_orbs = irclib.all_orbs(self.lines)
    val_orbs = op.num_valence_orbs(self.atoms())
    epsilon = []
    Epsilon = {}
    num_orbs = len(all_orbs[0][1])
    Epsilon["Reaction Coordinate"] = self.rxCoord()
    end = val_orbs
    start = 1
    for j in range(start,end):
      for coord in range(0,len(all_orbs)):
        epsilon.append(all_orbs[coord][1][j])
      Epsilon[str(j)] = epsilon
      epsilon = []
    return Epsilon

  def all_orbtitals(self):
    all_orbs = irclib.all_orbs(self.lines)
    epsilon = []
    Epsilon = {}
    num_orbs = len(all_orbs[0][2])
    Epsilon["Reaction Coordinate"] = self.rxCoord()
    start = 1
    for j in range(start, num_orbs):
      for coord in range(0,len(all_orbs)):
        epsilon.append(all_orbs[coord][2][j])
      Epsilon[str(j)] = epsilon
      epsilon = []
    return Epsilon

  def symm_orbitals(self):
      symm_orbs_all = irclib.all_symm_orbs_energ(self.lines)[0]
      symm_orbs_occ = sp.get_symm_orbs(self.lines)[0]
      symm_orbs_virt = sp.get_symm_orbs(self.lines)[1]
      print "\nSymmetries of occupied orbitals:"
      for key in sorted(symm_orbs_occ.iterkeys()):
          print str(key) +":  "+ symm_orbs_occ[key]
      count = 0
      print "\nSymmetries of occupied orbitals:"
      for key in sorted(symm_orbs_virt.iterkeys()):
          print str(key) +":  "+ symm_orbs_virt[key]
          count = count +1
          if count == 5:
              break
      print "###########################################"
      symm_orbs_all["Reaction Coordinate"] = self.rxCoord()
      return symm_orbs_all

  def chemPotKoopman(self):
    chemPot = irclib.koopman(self.lines)
    return {"Reaction Coordinate": self.rxCoord(), "Chemical Potential":chemPot}

  def chemPotGeneral(self,orbOcc,orbVirt)
    pass

  def flux(self, chemPot):
    coord = op.neg_derivative(self.rxCoord(),chemPot['Chemical Potential'])[0]
    flux  = op.neg_derivative(self.rxCoord(),chemPot['Chemical Potential'])[1]
    return {'Reaction Coordinate':coord, 'REF' : flux}

  def bondOrderDeriv(self):
    bondOrder = self.bondOrders()
    der_dict = {}
    #coord = op.neg_derivative(self.rxCoord,self.rxCoord)[0]
    for order in bondOrder:
      if order != "Reaction Coordinate":
        neg_der = op.neg_derivative(self.rxCoord(),bondOrder[order])
        der_dict[order] = neg_der[1]
        der_dict["Reaction Coordinate"] = neg_der[0]
    return der_dict
      
  def IP(self, outfile_cat):
    cat_lines = op.read_lines(outfile_cat)
    ip = irclib.get_IP(self.lines,cat_lines)
    return {"Reaction Coordinate": self.rxCoord(), "IP":ip}
    
  def EA(self, outfile_an):
    an_lines = op.read_lines(outfile_an)
    ea = irclib.get_EA(self.lines,an_lines)
    return {"Reaction Coordinate": self.rxCoord(), "EA":ea}

  def chemPotFinitDiff(self,ip, ea):
    chemPotFD = irclib.finite_diff(ip['IP'], ea['EA'])
    return {"Reaction Coordinate": self.rxCoord(), "Chemical Potential":chemPotFD}

  def simplePlot(self,x,y):
    return op.simple_plot(x,y)

  def savePlot(self,plotname, ylabel, zeroline=False, Show=True, **kwargs):
    print "I am using this IrcMagic class!!!!!"
    if not os.path.isdir("figures"):
      os.makedirs("figures")
    os.chdir("figures")
    works = self.ReactionWorks()
    op.general_plot(plotname, ylabel, works, zeroline, Show, **kwargs)
    os.chdir("../")

  def savePlotOrbs(self,plotname, ylabel, orblist, zeroline=False, **kwargs):
    if not os.path.isdir("orbitals"):
      os.makedirs("orbitals")
    os.chdir("orbitals")
    op.general_plot_orbs(plotname, ylabel, orblist, zeroline=False, **kwargs)
    os.chdir("../")

  def savePlotProp(self,plotname, ylabel, proplist, zeroline=False, Show=True,  **kwargs):
    if not os.path.isdir("figures"):
      os.makedirs("figures")
    os.chdir("figures")
    works = self.ReactionWorks()
    op.general_plot_prop(plotname, ylabel, proplist, works,  zeroline, Show, **kwargs)
    os.chdir("../")

  def save(self, filename, **kwargs):
    if not os.path.isdir("data"):
      os.makedirs("data")
    os.chdir("data")
    op.general_print(filename, **kwargs)
    os.chdir("../")

  def generate_cube_file(self,orb_range):
    atom_list = self.atoms()
    op.cube_files(orb_range,atom_list)

  def dual(self):
    if not os.path.isdir("density"):
        os.makedirs("density")
    op.fchk_gen("CHK")
    for file_ in glob.glob("CHK/*.fchk"):
        density =  Density(self.outfile, file_)
        os.chdir("density")
        density.dualFMOA()
        destiny  = "dualFMOA_"+file_.split("_")[1].split(".")[0]+".cub"
        print destiny
        os.system("mv dualFMOA/dualFMOA.cub dualFMOA/"+destiny)
        os.chdir("../")


