1) Running code used for all the simulation of "Paper" where melts of semiflexible polymers on a face centred cubic lattice has been studied in bulk (with periodic boundary conditions). 

Initial conditions: 

1) The polymer length and number of polymers to be simulated have to be set inside the main_ring.cpp

2) The code takes as input file the FCC-lattice of the polymer melts. The coordinates have to be presented as a single row of successive coordinates (3 for each monomer) to be put inside the ICS_RINGS folder with the name IC_RING_$Density_$Bending_rigidity_$Polymer_Length. This code works with the assumption of a mono-disperse solution of polymers.

e.g. Take a system of 2 polymer of 3 monomers each than the input file is as follows: 

0 1 0 0 2 0 0 1 0 3 1 2 ... 

3) The density of the system can be tuned by varying the Lattice Length of the box via the file Include/lattice.hpp

4) Custom initial conditions can be created by using the code inside Initial_Conditions folder. In this case a melt of double folded rings are constructed.  

5) By setting the bool variable to false then the simulation set as initial condition the last saved frame inside the directory DIR if present any. In the case DIR is empty an error will occour. 

Simulation: 

1) The length of the simulation and the time separation to save the frame are contained respectively in the variable MC_steps and stride.
2) The root folder to store the frames and all the data is DIR_MASTER. Specific folder in which each trajectory will be saved is DIR. 



