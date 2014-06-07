1D-CESE
=======

A 1D CESE C++ code

A very straightforward code, in the code, you can modify the value N and length to refine the "mesh "
#define N 10000
#define length 1

The values of x1 and x2 are the initial and final distance of the mesh. x2-x1 should be equal to length

#define x1 0
#define x2 1

To set the simulation time we should change T
#define T 0.1
To define the heat capacity ratio we use gama. 
#define gama 1.4
Cf stands for the Courant number, used to define the Courant-Friedrichs-Lewy condition
#define Cf 0.75
To calculate du/dx we use central diference method for the right (r) and left (l) side points of the point we are trying to calculate. We use a weight function of the form  W=((r*l^alpha)(l*r^alpha))/(l^alpha+r^alpha)
Alpha =0 for continous flows, alpha = 1,2,3.... for a flow with higher discontinuity, as the discontinuity increases alpha can be increased, but the diffusion also increases. 
#define alpha 6

We use rho1, u1 and p1 to define the density, velocity and pressure of the left side of the sod's tube
#define rho1 1
#define u1 0
#define p1 1
We use rho2, u2 and p2 to define the density, velocity and pressure of the right side of the sod's tube
#define rho2 0.125
#define u2 0
#define p2 0.1

