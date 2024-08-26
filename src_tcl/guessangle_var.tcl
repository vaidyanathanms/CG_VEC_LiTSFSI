package require topotools 1.6
topo readlammpsdata py_inpname full
topo numbonds
topo guessangles
topo writelammpsdata py_outname full
exit
