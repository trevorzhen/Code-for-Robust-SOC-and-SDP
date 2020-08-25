This folder includes all code used in the numerical experiments in Section 
6 for the paper "Robust optimization for models with uncertain SOC and SDP 
constraints" by Jianzhe Zhen, Frans de Ruiter, Ernst Roos & Dick den Hertog.

Most of the code is authored by Ernst Roos, except for the the files:
core_fme.m      
which is adopted from https://www.mathworks.com/matlabcentral/fileexchange/7957-fourier-motzkin-elimination

polyhedron_copositive.m
which is authored by Areesh Mittal (used for comparison)

VChooseK.c, VChooseK.m, and VChooseK.mex
which are authored by Jan Simon (obtained from Mathworks forum)

--------------------------------------------------------------------------------

The folder contains the following script files with described function:

core_fme.m
This function performs Fourier-Motzkin Elimination on a given linear system of 
inequalities.

find_vertices.m
This function finds all vertices of a given polyhedron.

fme_mve.m
This function extracts the proper constraint matrices to be used for 
Fourier-Motzkin elimination from the adjustable robust formulation of the 
minimum volume circumscribing ellipsoid problem.

polyhedron_copositive.m
This function finds an approximate solution to the minimum volume circumscribing
ellipsoid problem through a copositive programming approach introduced by
Mittal & Hanasusanto (2018).

RandomPolyhedronMittal.m
This function generates a random polyhedron of prescribed size according to the
procedure introduced by Mittal & Hanasusanto (2018).

SOCP_MVE.m
This script applies all provided methods to the selected instances of the
minimum volume circumscribing ellipsoid problem and prints the published
results.

SOCP_MVE_exact.m
This function solves the minimum volume circumscribing ellipsoid problem exactly
by enumerating all the polyhedron's vertices.

SOCP_MVE_full_quadratic.m
This function approximates the minimum volume circumscribing ellipsoid problem
through an adjustable robust linear optimization reformulation that is 
approximated with quadratic decision rules.

SOCP_MVE_linear.m
This function approximates the minimum volume circumscribing ellipsoid problem
through an adjustable robust linear optimization reformulation that is 
approximated with linear decision rules.

VChooseK.c, VChooseK.m, VChooseK.mex
This function creates a matrix, the rows of which are all combinations of
choosing K elements of a vector V without order and withour repetition.
