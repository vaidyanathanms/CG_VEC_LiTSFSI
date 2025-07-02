# CG_VEC_LiSTSFSI

Files to do CG simulations of VEC-LiSTSFI.

## Required depedencies/software

- Python3.x

- VMD-1.9.4

- Fortran

- LAMMPS

## Systems possible

- System-1: Cations with SIC copolymer comprising VEC and anions.

- System-2: Cations with blend of copolymer comprising VEC/anions and VEC monomers

- System-3: Cations with tertiary blend of VEC homopolymer, anion homopolymer and VEC monomers 

## Notes

Users can choose the number of beads constituting the VEC monomer

## Generating trajectories

Step-1: Load VMD, python and intel modules (if HPC is used)

Step-2: Set input parameters in `genconf.py` and run using `python3 genconf.py`

Step-3: Set `restart` option to 0 for the first run and 1 to subsequent runs.


## Analyzing trajectories

Step-1: To analyze the trajectories, set the parameters in `anainp.txt`

Step-2: Compile using
```
ifx -qopenmp -check all -traceback ana_params.f90 analyze.f90 -o ana.o
```

Step-3: Run using
```
./ana.o anainp.txt
```

