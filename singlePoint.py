import sys,re,math
import operations as op

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

def get_energy(lines):
  for line in lines:
    if "SCF Done:" in line:
      energy = line.split()[4]
      return energy

def get_orbitals(lines):
  # Obtain iteration range within the output file
  All_energies = []
  for Nline in range(0,len(lines)):
    if "The electronic state is" in lines[Nline]:
      startline=Nline+1
    if "Molecular Orbital Coefficients:" in lines[Nline]:
      endline=Nline
      break
    elif "Condensed to atoms"  in lines[Nline]:
      endline=Nline
      break

  # Extract the orbital energies within these lines
  Orbitalenergies_temp=[]
  Virtenergies_temp=[]
  Occenergies_temp=[]

  Orbitalenergies=[]
  Virtenergies=[]
  Occenergies=[]

  for Orbline in range(startline,endline):
    linelist=lines[Orbline].split()
    for lineitem in range(4,len(linelist)):
      orbenergy=linelist[lineitem]
      if "occ." in lines[Orbline]: 
        Occenergies_temp.append(orbenergy)
      if "virt." in lines[Orbline]: 
        Virtenergies_temp.append(orbenergy)
      Orbitalenergies_temp.append(orbenergy)	

  Orbitalenergies=process_occ_energies(Orbitalenergies_temp)
  Virtenergies=process_occ_energies(Virtenergies_temp)
  Occenergies=process_occ_energies(Occenergies_temp)
  
  return (Occenergies,Virtenergies,Orbitalenergies)

def bonddistance(xyz):
  bondLabel = []
  bondDistance = []
  atomLabel = op.atom_label(xyz[0])
  xCoord = xyz[1]
  yCoord = xyz[2]
  zCoord = xyz[3]

  for i in range(0,len(atomLabel)-1):
    for j in range(i,len(atomLabel)):
      if i is not j:
        bondLabel.append(atomLabel[i] + str(i+1) + "-" + atomLabel[j] + str(j+1))
        bondDistance.append(math.sqrt((float(xCoord[i])-float(xCoord[j]))**2 + (float(yCoord[i])-float(yCoord[j]))**2 + (float(zCoord[i])-float(zCoord[j]))**2))
  return(bondLabel,bondDistance)

def angles(xyz):
  angleLabel = []
  angle = []
  atomLabel = op.atom_label(xyz[0])
  xCoord = xyz[1]
  yCoord = xyz[2]
  zCoord = xyz[3]

  for i in range(0,len(atomLabel)-1):
    for j in range(i, len(atomLabel)):
      for k in range(0, len(atomLabel)):
        if i is j:
          continue
        elif i is k:
          continue
        elif j is k:
          continue
        else:
          ex_ki = float(xCoord[k]) - float(xCoord[i])
          ex_kj = float(xCoord[k]) - float(xCoord[j])
          ey_ki = float(yCoord[k]) - float(yCoord[i])
          ey_kj = float(yCoord[k]) - float(yCoord[j])
          ez_ki = float(zCoord[k]) - float(zCoord[i])
          ez_kj = float(zCoord[k]) - float(zCoord[j])
          rki = math.sqrt(ex_ki**2.0 + ey_ki**2.0 + ez_ki **2.0)
          rkj = math.sqrt(ex_kj**2.0 + ey_kj**2.0 + ez_kj **2.0)
          ex_ki = ex_ki/rki
          ey_ki = ey_ki/rki
          ez_ki = ez_ki/rki
          ex_kj = ex_kj/rkj
          ey_kj = ey_kj/rkj
          ez_kj = ez_kj/rkj
          angle.append(180.0/math.pi*math.acos(ex_ki*ex_kj + ey_ki*ey_kj + ez_ki*ez_kj))
          angleLabel.append(atomLabel[i] + str(i+1) + "-" + atomLabel[k] + str(k+1) + "-" + atomLabel[j] + str(j+1))
  return(angleLabel,angle)


def get_xyz(lines):
  xyz_info = []
  atoms_num = []
  xyz = []
  x_coords = []
  y_coords = []
  z_coords = []
  found = False
  for lineNum in range(0,len(lines)):
    if "Standard orientation:" in lines[lineNum]:  
      found = True
      for lineNum1 in range(lineNum+5,len(lines)):
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
          xyz.append(atoms_num)
          xyz.append(x_coords)
          xyz.append(y_coords)
          xyz.append(z_coords)
          break
    elif found:
      break
  return xyz

def get_last_xyz(lines):
  xyz_info = []
  atoms_num = []
  x_coords = []
  y_coords = []
  z_coords = []
  for lineNum in range(0,len(lines)):
    #if "Input orientation:" in lines[lineNum]:  
    if "Standard orientation:" in lines[lineNum]:  
      xyz = []
      atoms_num = []
      x_coords = []
      y_coords = []
      z_coords = []
      for lineNum1 in range(lineNum+5,len(lines)):
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
          xyz.append(atoms_num)
          xyz.append(x_coords)
          xyz.append(y_coords)
          xyz.append(z_coords)
          break
  return xyz

def get_orb_symm(lines):
  symmetries = []
  found = False
  for lineNum in range(0,len(lines)):
    if "Population analysis using the SCF density." in lines[Nline]:
      found = True
      for lineNum1 in range(lineNum+4,len(lines)):
        linelist=lines[lineNum1].split()
        if not Virtual in lines[lineNum1]: 
          if "Occupied" in lines[lineNum1]:
            for lineitem in range(1,len(linelist)):
              symmetries.append(lineitem.strip(')').strip('('))
          else:
            for lineitem in range(0,len(linelist)):
              symmetries.append(lineitem.strip(')').strip('('))

    elif found:
      break
  return symmetries

# make it work for level of theory on any place of that line
def get_level_of_theory(lines):
  for lineNum in range(0,len(lines)):
    if "#P" in lines[lineNum]:
      theory = lines[lineNum].split()[1].split("/")[0]    
      return theory
  print "No level of theory found"

def get_basis(lines):
  for lineNum in range(0,len(lines)):
      if "Standard basis:" in lines[lineNum]:
        basis = lines[lineNum].split()[2]    
  return basis

def get_charge(lines):
  for lineNum in range(0,len(lines)):
      if "Charge =" in lines[lineNum]:
        charge = lines[lineNum].split()[2]    
  return charge

def get_multiplicity(lines):
  for lineNum in range(0,len(lines)):
      if "Multiplicity =" in lines[lineNum]:
        multi = lines[lineNum].split()[5]    
  return multi
