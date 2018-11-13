%mem=2000MB
%chk=model_amber2.chk
%nprocshared=1
#p oniom(amber) geom=connectivity

Gaussian input prepared by ASE

0 1
 H-HW-1.0000(PDBName=HW,ResName=WAT,ResNum=1)   0   0.00000    1.00000   0.000000 H
 O-OW--1.000(PDBName=OW,ResName=WAT,ResNum=1)    0   0.00000    -1.00000   0.000000 H

 1 
 2
 

NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000
