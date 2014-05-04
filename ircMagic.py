#! /usr/bin/python

import sys, os
sys.path.append('/home/stvogt/qtc/scripts')
import operations as op
import irc, chem_pot
import singlePoint as sp

class Irc:
  def __init__(self, outfile):
    self.outfile = outfile
    self.lines = op.read_lines(outfile)

  def rxCoord(self):
    rx_coord = chem_pot.get_rxCoord(self.lines)  
    return rx_coord

  def atoms(self):
    return sp.get_xyz(self.lines)[0]

  def energy(self):
    energy = chem_pot.all_energies(self.lines)
    return {"Reaction Coordinate": self.rxCoord(), "Energy":energy}

  def force(self):
    coord = op.neg_derivative(self.rxCoord(),self.energy()['Energy'])[0]
    force = op.neg_derivative(self.rxCoord(),self.energy()['Energy'])[1]
    return {"Reaction Coordinate":coord , "Reaction Force":force}

  def ReactionWorks(self):
    y = self.force()["Reaction Force"]
    x = self.force()["Reaction Coordinate"]
    max_force = max(y)
    min_force = min(y)

    for i in range(0,len(y)):
      if y[i] == min_force:
        w1_cut = i
        w1_rxcoord = x[i]
      elif y[i] == max_force:
        w3_cut = i+1
        w3_rxcoord = x[i+1]

    for i in range(0,len(y)):
      if y[i] > 0.0:
        w2_cut = i
        break

    w1_y = y[:w1_cut]
    w1_x = x[:w1_cut]
    w2_y = y[w1_cut:w2_cut]
    w2_x = x[w1_cut:w2_cut]
    w3_y = y[w2_cut:w3_cut]
    w3_x = x[w2_cut:w3_cut]
    w4_y = y[w3_cut:]
    w4_x = x[w3_cut:]

    w1 = round(op.integrate(w1_x,w1_y),2)
    w2 = round(op.integrate(w2_x,w2_y),2)
    w3 = round(op.integrate(w3_x,w3_y),2)
    w4 = round(op.integrate(w4_x,w4_y),2)
    return (w1,w2,w3,w4,w1_rxcoord,w3_rxcoord)

  def distances(self):
    print "----Distances---"
    sigma = []
    Sigma = {}
    all_dis = chem_pot.all_distances(self.lines)
    num_dist = len(all_dis[0][0])
    Sigma["Reaction Coordinate"] = self.rxCoord()
    for j in range(0, num_dist):
      for coord in range(0,len(all_dis)):
        sigma.append(all_dis[coord][1][j])
      Sigma[all_dis[0][0][j]] = sigma
      #print all_dis[0][0][j]
      sigma = []
    return Sigma

  def bondOrders(self):
    Bndo = {}
    bndo_list = []
    bondOrder =  chem_pot.all_bondOrders(self.lines)
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
    charges =  chem_pot.all_NatCharges(self.lines)
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
    chem_pot.all_bondOrbitals(self.lines)

  def occ_orbtitals(self):
    all_orbs = chem_pot.all_orbs(self.lines)
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
    all_orbs = chem_pot.all_orbs(self.lines)
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
    all_orbs = chem_pot.all_orbs(self.lines)
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

  def chemPotKoopman(self):
    chemPot = chem_pot.koopman(self.lines)
    return {"Reaction Coordinate": self.rxCoord(), "Chemical Potential":chemPot}

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
    ip = chem_pot.get_IP(self.lines,cat_lines)
    return {"Reaction Coordinate": self.rxCoord(), "IP":ip}
    
  def EA(self, outfile_an):
    an_lines = op.read_lines(outfile_an)
    ea = chem_pot.get_EA(self.lines,an_lines)
    return {"Reaction Coordinate": self.rxCoord(), "EA":ea}

  def chemPotFinitDiff(self,ip, ea):
    chemPotFD = chem_pot.finite_diff(ip['IP'], ea['EA'])
    return {"Reaction Coordinate": self.rxCoord(), "Chemical Potential":chemPotFD}

  def simplePlot(self,x,y):
    return op.simple_plot(x,y)

  def savePlot(self,plotname, ylabel, zeroline=False, Show=True, **kwargs):
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

  def generate_cube_file(self,orb_range,sum_orbs=False):
    atom_list = self.atoms()
    op.cube_files(orb_range,atom_list,sum_orbs)
