%mem=2000MB
%chk=model_amber2.chk
%nprocshared=1
#p oniom(amber) geom=connectivity

Gaussian input prepared by ASE

0 1
H-HO-0.6000(PDBName=HW,ResName=WAT,ResNum=1)   0  -1.00000     1.00000   0.000000 H
O-OH--0.600(PDBName=OW,ResName=WAT,ResNum=1)   0  -1.00000     0.00000   0.000000 H
O-OH--0.600(PDBName=OW,ResName=WAT,ResNum=1)   0   1.00000     0.00000   0.000000 H
H-HO-0.6000(PDBName=HW,ResName=WAT,ResNum=1)   0   1.00000     0.00000   1.000000 H

 1 2 1.0
 2 3 1.0
 3 4 1.0
 4
 

NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000
HrmStr1 HO OH  553.0  0.950
HrmStr1 OH OH  353.0  1.475
HrmBnd1 HO OH OH  60.0  94.8
AmbTrs HO OH OH HO    0   0   0   0  0.000  0.000  1.400  0.000 1.0