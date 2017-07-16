# Compile code:
Go to user directory:

cd /src/Users/dargent/

Toy fit with full PDF:

cd /BsDsKpipi_TD

Toy fit with TD PDF (DsK like) and test with B2DX-Fitter classes:

cd /B2DX-Tests

make -j 2

# Run it:
All options like number of generated events, amplitude model etc. are set in steering files: 

TD_Bs.txt 

(timeFit.txt)

see README in the corresponding directories for more information

./ampFit 0 < TD_Bs.txt 

(./timeFit < timeFit.txt)

# How to use classes from B2DXFitter?

Some classes are already included, for example: RooCubicSplineFun

If you need additional classes:

1) Copy headers to TD-MINT2/Mint/

2) Add correct directory to included files: 

e.g: #include "RooCubicSplineKnot.h"  --> #include "Mint/RooCubicSplineKnot.h"

3) Remove ClassDef(RooCubicSplineFun,1)

4) Copy src files to TD-MINT2/src/Mojito/TimeDependent/src/

5) Rename file from *.cxx to *.cpp

6) Again add correct directory to included files

7) Classes can now be used in your ampFit.cpp code as usual



