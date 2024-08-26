package require topotools 1.6
topo readlammpsdata VECdata20.dat full
topo numbonds
topo guessangles
topo writelammpsdata VECmain_20.dat full
exit
