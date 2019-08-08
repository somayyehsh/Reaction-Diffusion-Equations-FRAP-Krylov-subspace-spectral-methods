# SciCode

This is the solution of PDEs for first-order photobleaching kinetics (Reaction-Diffusion Equations) using Krylov subspace spectral methods written in Matlab. Please, cite both below papers. Feel free to contact me for any additional information at sheikho.physics@gmail.com.


	S. Sheikholeslami, James V. Lambers, “Modeling of first-order 
	photobleaching kinetics using Krylov subspace spectral methods”, 
	Computers and Mathematics with Applications, 2017, 
	https://doi.org/10.1016/j.camwa.2017.10.019.

	S. Sheikholeslami, James V. Lambers, Carley Walker, 
	“Convergence Analysis of Krylov Subspace Spectral Methods for 
	Reaction-Diffusion Equations”, Journal of Scientific Computing, 2018, 
	http://doi.org/10.1007/s10915-018-0824-5.


Files for performance comparisons:

kssmain - helper function called by testopt2
testopt1 - optimized version of test, that performs only one time step
testopt1N - same as testopt1, for Neumann boundary conditions
testopt2 - optimized version of test, that performs time-stepping
crank - Crank-Nicholson
eulerforward - forward Euler
rk4 - 4th-order
methodwrap - used to run test cases using either crank, eulerforward or 
rk4. Sample usage:

[solntime,abserr,relerr]=methodwrap(@crank,N,nsteps,tf,coefs)

methodwrapN - same as methodwrap, for Neumann boundary conditions
refsolnN - same as refsoln, for Neumann boundary conditions
mydct, mydct2, myidct, myidct2 - cosine transform functions for Neumann 
boundary conditions

