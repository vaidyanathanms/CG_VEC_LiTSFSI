#####################################################
# Extra code to compute Rg                          #
# Version: May-20-2024                              #
#####################################################


# Initialization
units		lj
boundary	p p p
atom_style	full
log 		log.txt
read_restart	archival_*.restart

# Neighbor information
neighbor        0.3 bin
neigh_modify	every 10 one 10000

#####################################################
# Interaction parameters
#####################################################

# Pair Information
include         pair_info_file.dat
pair_modify	shift yes mix arithmetic
kspace_style	pppm 0.00001
dielectric	1

# Pair Coeffs - read from restart file

# Bond Information
bond_style	harmonic
special_bonds   lj/coul 0 0 1

# Bond Coeffs
include         bond_info_file.dat

# Angle Information
angle_style  cosine/squared

# Angle Coeffs
angle_coeff * 3 110   # Set as default to account for the absence of 3-3-3
angle_coeff 1 3 110   # 1 -> 1-1-1
angle_coeff 2 5 70    # 2 -> 1-1-2 -> Stronger angle
angle_coeff 3*5 3 110 # 3-5 -> 1-1-3;1-3-1;1-3-3
angle_coeff 6 3 70    # 6 -> 2-1-3


# Dihedral Information
dihedral_style	none
improper_style 	none


# Initiate/Write Atoms and Other Details
variable        ts equal step
thermo_style	custom step temp press pe evdwl ecoul epair ebond eangle
thermo          10000
dump            main all custom 10000 config_${ts}.lammpstrj id type xu yu zu

# Code to compute Rg
compute cc1 all chunk/atom molecule
compute myChunk all gyration/chunk cc1
fix myrg all ave/time 100 1 100 c_myChunk file tmp.out mode vector

#####################################################
# Equilibration (NVT dynamics at kT=1)
#####################################################

# Main Fixes - NVT
fix 1 all nvt temp 1.0 1.0 0.5

# Run Styles
run_style verlet
restart   5000 restart1 restart2
restart   500000 archival_*.restart
timestep  0.01
run	  30000 
unfix 1
write_restart 	restart_NVT.dat


