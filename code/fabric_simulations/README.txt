These simulations were developed by M. Dimitriyev and S. Gonzalez.

For questions, email sgonzalez49@gatech.edu or visit matsumoto.gatech.edu for the contact information of current graduate students.

-------------
VERSION INFORMATION
Python code developed for v3.6.13 of the Python Programming Language, using v1.18.1 of NumPy and v1.4.1 of SciPy libraries as packaged in Anaconda3 (https://www.anaconda.com/)

-------------
RUNNING THE SIMULATIONS

To run, execute the following script in terminal: python3 stockinette.py config_file_name.dat
This script works for all four fabric-type simulations, where "stockinette" is subbed out for the desired fabric type.

To run, you must have an initialization file, named "Init_XXXX.dat", and a configuration file that is called in the terminal. The XXXX in the initialization file name can be subbed out for any number or string. That XXXX suffix is included in the configuration file to indicate the initialization file that is being called into the simulation. The initialization file must be in the same folder as the simulation code.

INITIALIZATION FILE
The initialization file has a list of inputs used to construct an initial Bezier curve for the stitch. Example initialization files for each fabric type have been provided. To create new initialization files, you can either add a small amount of random noise to the provided initialization file or use the output of one of the simulations of the same fabric type. To make a new init file from a simulation, take the ResultsOut.dat output file and add in the missing cell dimension. For example, if you ran your simulation to stretch in the x-direction, you will have to add in the x-dimension box width to the ResultsOut.dat file to turn it into an initialization file. The x-dimension is always the second-to-last entry and the the y-dimension is the last entry of the initialization files.

CONFIGURATION FILE
The configuration files include a list of fabric and yarn parameters. It also specifies the prefix of the output file name, the stretching direction, and the initialization file needed to start the simulation. Example configuration files have been included for each fabric type, where the fabrics are made of acrylic yarn and stretched in the x-direction.

The configuration file also tells the simulation the cell dimension to start and end stretching. To derive a stress-v-strain curve from the simulation outputs, this range must include the zero-force simulation box dimension. The zero-force box dimension will be the lowest energy configuration. If stretching the fabric, you should look for simulation outputs where the total energy, listed in EnergyOut.dat, decreases and then increases as the stitch cell box increases in size. This will ensure you have captured the zero-force configuration in your sweep. 

SIMULATION STEP SIZE
The step size the simulation takes to increase the cell dimension between your lower and upper limit is hard-coded in for the different stretching directions. Search for "step_size" in each simulation and you will be able to change them. Generally, they are currently set to generate 20-30 simulation steps between the zero-force and ~0.4 N/mm stress mark for each fabric stretched in each direction with the provided configuration file parameters. A smaller step size gives more accuracy in determining the zero-force configuration.

For compression data, you can set the step_size to a negative number. Please note that the energy quickly diverges upon compression so you should start your simulations near the zero-force state and decrease the stitch cell dimension with a negative step size. If you start in a large (order ~100 J) energy state, the simulation is unlikely  to converge well to a lower energy state.

-------------
ASSOCIATED SCRIPTS

Bezier.py, ConnectingCurve.py, Crossover.py, and mathHelper.py are associated scripts that construct and update the Bezier curve splines used in the simulations. These codes must be in the same folder as the fabric simulations for the simulations to run. Initialized Bezier curves are provided via the initialization files.

-------------
SIMULATION MECHANICS
For an in-depth explanation of the physics utilized and code mechanics, see the SI of "Programming Mechanics in Knitted Materials, Stitch by Stitch" by Singal et al.

The simulations take a stitch constructed of Bezier curve splines, calculate the bending and compression energy along the curves, and iteratively change the shape of the curve to minimize the total energy for the specified stitch cell dimension. The unspecified stitch cell dimension can vary freely as an unclamped edge to reflect the Poisson ratio of the fabric. These simulations are static and would reflect stretching experiments done in the quasi-static regime.

BEZIER SPLINES
The Bezier curve splines constructed in the simulation are labelled according to subfigure (a) in "splinenotation.png". Subfigure (b) shows the total shape of the stitch from a top view (left) and angled perspective (right). 

MINIMIZATION SCHEME
Standard scipy.minimize is used to minimize the total energy, using the sequential least-squares ("SLSQP") method. Due to the large dimension of the energy landscape, the simulation can fall into local minima. 

****To increase the likelihood of finding a global minima, run the simulation multiple times for the same configuration parameters but change the initialization used to start the simulation.****

Typically, four initialization files for a single simulation provides good stability. Garter can have convergence issues for very small step sizes, which can be solved by running more initializations.

--------------
SIMULATION OUTPUTS
When a simulation code is started, it writes the input configuration file into an output configuration file. This allows you to delete the input configuration file and still keep track of the parameters for each simulation set. 

The simulations output ten files for each set of cell dimensions. 

"BendingDensityOut.dat", "BezierOut.dat", "BezierOutCourse.dat", and "BezierOutWale.dat" provide information about the final Bezier curve spline. "BezierOut.dat", "BezierOutCourse.dat", and "BezierOutWale.dat" list control points for the Bezier curves. See the print function in Bezier.py for how the control points are printed. See the simulation code for which splines are printed to which file.

"ResultsOut.dat" provides information that can be used to generate a new initialization. See the INITIALIZATION FILE section of this text for instructions. 

"ContactDensityOut.dat" and "ContactMapOut.dat" provide information about the contact energy density.

"EnergyOut.dat", "StretchOut.dat", and "CellPropsOut.dat" can be used to generate stress v strain data. "CellPropsOut.dat" has information on shearing via the join and cell angles. Currently, the simulation is hard-coded for uniaxial stretching or compression and these angles are set to zero. "EnergyOut.dat" has three entries: the contact energy, then the bending energy, then the total energy. "StretchOut.dat" has three entries: the length of yarn per stitch, the stitch cell x-dimension, and the stitch cell y-dimension.

-------------
ANALYZING RESULTS
To produce stress and strain data, first interpolate the energy v stitch cell dimension in the direction of stretching. To reduce kinks, only interpolate within the low-energy region. The minimum value of the interpolation tells you where the zero-force state is. Record the cell dimension of the zero-force state and use interpolation or linear fitting to get the other cell dimension at zero force.

To get force data from the energy v cell dimension, take a discrete, midpoint derivative. To transform the force into stress, divide by the cell dimension you weren't sweeping over at zero-force. To get strain from the cell dimension data, calculate (a_x-a_x0)/a_x0 where a_x is the x-dimension of the cell and a_x0 is the x-dimension of the cell at zero force. The same calculation can be made for the y-dimension.

Remember to run multiple initializations and choose the minimum energy simulation for a given stitch cell dimension. Before doing your stress strain analysis. This will give the best results.

A Mathematica code to analyze results is included in the sample_analysis folder, as well as sample data for stockinette. The final output of the code is a list of stress and strain data for plotting. 



