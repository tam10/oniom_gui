FN="AmberFort"
gfortran -O2 -g -fPIC -c -o "$FN".o "$FN".f90
gfortran -O2 -g -fPIC -dynamiclib "$FN".o LBFGSB_Opt/Build/*.o -o "$FN".dll
gfortran -O2 -g -fPIC -shared "$FN".o LBFGSB_Opt/Build/*.o -o "$FN".bundle
