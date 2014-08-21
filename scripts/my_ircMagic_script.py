#! /usr/bin/python

import sys
from ircMagic import Irc

outfile = sys.argv[1]

##Initiate molecule object
HSCHO = Irc(outfile)

#calculates and prints the reaction works
works = HSCHO.ReactionWorks()
print works

# Plot energy profile
energy = HSCHO.energy()
HSCHO.savePlot('energy.png', "Energy", **energy)
HSCHO.save("energy.dat",**energy)


#plot amd save reaction force profile
force = HSCHO.force()
print force
HSCHO.savePlot('force.png',"Force", **force)
HSCHO.save("force.dat",**force)


#plot amd save chemical potential profile
mu = HSCHO.chemPotKoopman()
HSCHO.savePlot('chem_pot.png',"mu", **mu)
HSCHO.save("chem_pot.dat",**mu)

#plot and save REF
flux_k  = HSCHO.flux(mu)
HSCHO.savePlot('flux_k.png',"REF", **flux_k)
HSCHO.save("flux_k.dat",**flux_k)

#Obtain all HSCHO moecular orbital energies
all_orbitals = HSCHO.all_orbtitals()

#Obtain all symm_orbitals
symm_orbitals = HSCHO.symm_orbitals()

# Save orbitlas
HSCHO.save("orbitals.dat", **all_orbitals)

#plot amd save chemical potential profile with symmetry orbitals
mu = HSCHO.chemPotGeneral("13A'",'4A"')
HSCHO.savePlot('chem_pot_13A_4A.png',"mu", **mu1)

#Obtain flux with symmetry chem pot
flux1  = HSCHO.flux(mu1)
HSCHO.savePlot('flux_13A_4A.png',"REF", **flux1)

#Save some orbitals
HSCHO.savePlotProp("orbitals_3A''_13A'.png","Occ. Orbs.",['3A"',"13A'"], **symm_orbitals)

#Obtain Bond orders
bondorders = HSCHO.bondOrders()
HSCHO.save("bndord.dat", **bondorders)
HSCHO.savePlotProp("BondOrder1.png","Bond Order",["S2-C3","C3-O4"] , **bondorders)

#Obtain Bond order derivatives
bnd_der = HSCHO.bondOrderDeriv()
HSCHO.save("bndordder.dat", **bnd_der)
HSCHO.savePlotProp("BondOrderDer.png","Bond Order",["S2-C3","C3-O4"],True, **bnd_der)
HSCHO.savePlotProp("BondOrderDer_H.png","Bond Order",["H1-S2","H1-O4"],True, **bnd_der)

##Obtain Bond distances
distances = HSCHO.distances()
HSCHO.save("distances.dat", **distances)
HSCHO.savePlotProp("Distance_SNO.svg","Distances",["S2-N3","N3-O4"] ,True, **distances)

#Obtain Angles
angles = HSCHO.angles()
HSCHO.save("angels.dat", **angles)
HSCHO.savePlotProp("angles.png","Distances",["C1-F6-H4", "C1-F5-H4"] ,True, **angles)

#Obtain atomic natural charges
charges = HSCHO.natCharges()
HSCHO.save("charges.dat", **charges)
HSCHO.savePlotProp("Charges_SNO.svg","Natural Charge",["S2","N3","O4"] ,True, **charges)


#Obtain cubes for Dual descriptor
HSCHO.dual()

##################################################################################


#Open ionized outputs
outfile_cat = sys.argv[2]
outfile_an = sys.argv[3]

#Compute ionization potential and chemical potential using finite differnece approximation
ip = HSCHO.IP(outfile_cat)
ea = HSCHO.EA(outfile_an)

# Obtain chem pot from finite differences
mu_fd = HSCHO.chemPotFinitDiff(ip, ea)
flux_fd = HSCHO.flux(mu_fd)

# Extract Plots
HSCHO.savePlot('flux_fd.png',"REF",**flux_fd)
HSCHO.savePlot('ip.png',"ip",**ip)
HSCHO.savePlot('ea.png',"ea",**ea)

#separate occupied and virtual orbitals, a bit outdateed better use all the orbitals
occ_orbitals = HSCHO.occ_orbtitals()
virt_orbitals = HSCHO.virt_orbtitals()



