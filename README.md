# cell-gradients
 
The two subfolders (numerics and simulations) contain two different sets of code. Numerics has the codes used to numerically compute the concentration gradients based on the analytic solutions of the differential equations. Simulations has codes for particle based simulations based on the same equations.

We list out what the codes do:

## NUMERICS
sphere.py is the code that computes the solution for a spherical cell. Corresponds to figure 2.   
cylinder.py does the same for a cylindrical cell. Corresponds to figure 3.   
ellipse.cpp does the same for a spheroidal cell. Corresponds to figure 4.  
sphere_reflect.cpp does the same for a spherical cell with partially absorbing boundaries. Corresponds to figure 5.  

## SIMULATIONS
sphere_pz.cpp is the code that simulates our model in a spherical cell. Corresponds to figure 2.  
cylinder.cpp does the same for a cylindrical cell. Corresponds to figure 3.  
ellipse.cpp does the same for a spheroidal cell. Corresponds to figure 4.  
sphere_reflect.cpp does the same for a spherical cell with partially absorbing boundaries. Corresponds to figure 5.  
sphere_pz_fall.cpp does the same for a spherical cell where the particles being transported on the surface can fall off before they reach the pole. Corresponds to figure 6.  
