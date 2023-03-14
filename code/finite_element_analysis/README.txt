This README.txt file was generated on 2023/03/13 by Michael S. Dimitriyev (michael.dimitriyev@gmail.com).

*******************
* INSTRUCTIONS *
*******************

This directory contains the script uniaxial_FEA_loop.py, which is run via the terminal command
> python uniaxial_FEA_loop.py
and exports the following files:
- my_mesh.xml: contains mesh data
- out_(DISP)/disp.out: displacement field
- out_(DISP)/ed.out: energy density field
- out_(DISP)/energy.out: total energy

This script sequentially solves for the displacement field of a 2D continuum fabric, where the top boundary is displaced by a distance DISP, in units of the fabric's original height H. The current displacement DISP is incremented from 0.0 by a value specified by the variable dy, until a maximum displacement is reached, which is set by the variable DEF_MAX. The results of each sequential stretch are stored in a folder of the form out_(DISP). 
	
VERSION INFORMATION:
    - Python3 (v3.7.3), as packaged in Anaconda3 (release 2019.03).
    - FEniCS (2019.1.0) downloaded via conda-forge (https://fenicsproject.org/download/archive/)


	




