#! /usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import singlePoint as sp
import nbo

def process_occ_energies(occ_energies):
  processed_energies = []
  for item in occ_energies:
    if len(item) > 10:
      gaussian_fuckup = item.split("-")
      gaussian_fuckup = filter(None, gaussian_fuckup) 
      for item_stick in gaussian_fuckup:
        processed_energies.append("-"+item_stick)
    else: 
      processed_energies.append(item)
  return processed_energies

def get_rxCoord(lines):
  rx_coord = []
  for Nline in range(0,len(lines)):
    line = lines[Nline]
    if "Single Point computation for reaction coordinate: " in line:
      if not "\\" in line:
        if not "Job title:" in line:
          rx_coord.append(line.split()[6])
  return rx_coord

def get_blocks(lines):
  block_lines = []
  blocks = []
  rx_coord = get_rxCoord(lines)
  for coord in rx_coord:
    startline = False
    for Nline1 in range(0,len(lines)):
      line1 = lines[Nline1]
      if "Single Point computation for reaction coordinate: "+coord in line1: 
        if not "\\" in line1:
          startline = True
          for Nline2 in range(Nline1,len(lines)):
            line2 = lines[Nline2]
            block_lines.append(line2)
            if "Normal termination of Gaussian 09" in line2:
              break
      elif startline:
        blocks.append(block_lines)
        block_lines = []
        break

  return blocks

def all_energies(lines):
  energy = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    energy.append(sp.get_energy(blocks[blockNum]))
  energy[:] = [(float(x) - float(energy[0]))*627.509469 for x in energy] 
  return energy

def all_orbs(lines):
  all_orbs = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    all_orbs.append(sp.get_orbitals(blocks[blockNum]))
  return all_orbs

def all_symm_orbs_energ(lines):
    occ_energ_symm = {}
    virt_energ_symm = {}
    blocks = get_blocks(lines)
    for blockNum in range(0,len(blocks)):
        # Get the dictcionary with the number of orbitals and respective symmetries
        Occ_dict_num = sp.get_symm_orbs(blocks[blockNum])[0]
        Virt_dict_num = sp.get_symm_orbs(blocks[blockNum])[1]
        #Get the orbital energies:
        occ_energ = sp.get_orbitals(blocks[blockNum])[0]
        virt_energ = sp.get_orbitals(blocks[blockNum])[1]
        all_orbs = sp.get_orbitals(blocks[blockNum])[2]
        #In the first block initialize the symm_orbital dictcionary
        if blockNum == 0:
            for itemNum in range(1,len(occ_energ)+1):
                symm = Occ_dict_num[itemNum]
                occ_energ_symm[symm] = []
            for itemNum in range(len(occ_energ)+1, len(all_orbs)+1):
                symm = Virt_dict_num[itemNum]
                virt_energ_symm[symm] = []
        # Appending the orbital energies of that point in the IRC with the orbtial enengies 
        for itemNum in range(1,len(occ_energ)+1):
            symm = str(Occ_dict_num[itemNum])
            occ_energ_symm[symm].append(occ_energ[itemNum-1])
        for itemNum in range(len(occ_energ)+1, len(all_orbs)+1):
            #print itemNum
            symm = str(Virt_dict_num[itemNum])
            virt_energ_symm[symm].append(virt_energ[itemNum-len(occ_energ)-1])
    all_energ_symm = dict(occ_energ_symm.items() + virt_energ_symm.items())        
    return(all_energ_symm, occ_energ_symm, virt_energ_symm)

def all_distances(lines):
  all_dist = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    xyz = sp.get_xyz(blocks[blockNum])
    dis = sp.bonddistance(xyz)
    all_dist.append(dis)
  return all_dist

def all_angles(lines):
  all_angles = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    xyz = sp.get_xyz(blocks[blockNum])
    ang = sp.angles(xyz)
    all_angles.append(ang)
  return all_angles

def all_bondOrders(lines):
  all_bndo = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    bndo = nbo.bondOrder(blocks[blockNum])
    all_bndo.append(bndo)
  return all_bndo

def all_NatCharges(lines):
  all_charges = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    charges = nbo.natCharges(blocks[blockNum])
    all_charges.append(charges)
  return all_charges

def all_bondOrbitals(lines):
  orbs = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    natorb = nbo.bondOrbitals(blocks[blockNum])
    for key in natorb:
     print key
     print natorb[key]
     #pass
    break
    orbs.append(charges)
  return orbs


def koopman(lines):
  orbs = all_orbs(lines)
  Mu = []
  
  #The first index are the points along the IRC, in the second index, 0 is occupied, 1 is virtual, and the thrid index are the individual orbitals
  for i in range(0,len(orbs)):
    homo = float(orbs[i][0][-1])
    lumo = float(orbs[i][1][0])
    mu = -0.5*(homo + lumo)*627.509469
    Mu.append(mu)
  return Mu

def get_IP(lines_neut,lines_cat):
  cat_enrg = all_energies(lines_cat)
  neut_enrg = all_energies(lines_neut)
  I = []
  for enrg_item in range(0,len(cat_enrg)):
    ip = float(cat_enrg[enrg_item])-float(neut_enrg[enrg_item])
    I.append(ip)
  return I

def get_EA(lines_neut,lines_an):
  an_enrg = all_energies(lines_an)
  neut_enrg = all_energies(lines_neut)
  A = []
  for enrg_item in range(0,len(an_enrg)):
    ea = float(neut_enrg[enrg_item])-float(an_enrg[enrg_item])
    A.append(ea)
  return A

def finite_diff(ip,ea):
  MU = [] 
  for i in range(0,len(ip)):
    mu = 0.5*(float(ip[i])+float(ea[i]))*627.509469
    MU.append(mu)
  return MU

def get_all_xyz(lines):
  xyz = []
  blocks = get_blocks(lines)
  for blockNum in range(0,len(blocks)):
    xyz.append(sp.get_xyz(blocks[blockNum]))
  return xyz

#def get_orbs(lines):
#  # Extract the orbital energies within these lines
#  All_occ = []
#  All_virt = []
#  rx_coord = []
#  for Nline in range(0,len(lines)):
#    if "Single Point computation for reaction coordinate:" in lines[Nline]:
#      if not "\\" in lines[Nline]:
#        rx_coord.append(lines[Nline].split()[6])
#    if "The electronic state is" in lines[Nline]:
#      Orbitalenergies_temp=[]
#      Virtenergies_temp=[]
#      Occenergies_temp=[]
#      Orbitalenergies=[]
#      Virtenergies=[]
#      Occenergies=[]
#      for Orbline in range(Nline,len(lines)):
#        if not "Condensed to atoms" in lines[Orbline]:
#          linelist=lines[Orbline].split()
#          for lineitem in range(4,len(linelist)):
#            orbenergy=linelist[lineitem]
#            if "occ." in lines[Orbline]: 
#              Occenergies_temp.append(orbenergy)
#            if "virt." in lines[Orbline]: 
#              Virtenergies_temp.append(orbenergy)
#            Orbitalenergies_temp.append(orbenergy)	
#        else:
#          break
#
#      Orbitalenergies=process_occ_energies(Orbitalenergies_temp)
#      Virtenergies=process_occ_energies(Virtenergies_temp)
#      Occenergies=process_occ_energies(Occenergies_temp)
#      
#      All_occ.append(Occenergies)
#      All_virt.append(Virtenergies)
#
#  #return (Occenergies,Virtenergies,Orbitalenergies,All_energies)
#  return (All_occ, All_virt, rx_coord)
#
#def koopmans_old(lines):
#  occ = get_orbs(lines)[0]
#  virt = get_orbs(lines)[1]
#  
#  HOMO = []
#  LUMO = []
#  Mu = []
#  
#  for i in range(0,len(occ)):
#    HOMO.append(occ[i][-1])
#    LUMO.append(virt[i][0])
#    mu = -0.5*(float(occ[i][-1]) + float(virt[i][0]))
#    Mu.append(mu)
#
#  return(Mu, HOMO, LUMO)
#def finite_diff(lines_neut, lines_an, lines_cat):
#  cat_enrg = get_energy(lines_cat)[1]
#  an_enrg = get_energy(lines_an)[1]
#  neut_enrg = get_energy(lines_neut)[1]
#  rx_coord = get_energy(lines_neut)[0]
#  I = []
#  A = [] 
#  MU = [] 
#  for enrg_item in range(0,len(cat_enrg)):
#    ip = float(cat_enrg[enrg_item])-float(neut_enrg[enrg_item])
#    af = -float(an_enrg[enrg_item])+float(neut_enrg[enrg_item])
#    #ip = float(neut_enrg[enrg_item]) - float(cat_enrg[enrg_item])
#    #af = float(neut_enrg[enrg_item]) - float(an_enrg[enrg_item])
#    mu = 0.5*(float(ip)+float(af))
#    I.append(ip)
#    A.append(af)
#    MU.append(mu)
#  return (MU,rx_coord, I, A, cat_enrg, an_enrg, neut_enrg)


