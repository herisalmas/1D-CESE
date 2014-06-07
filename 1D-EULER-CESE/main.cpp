//
//  main.cpp
//  1D-EULER-CESE
//
//  Created by Heri Salmas on 11/26/13.
//  Copyright (c) 2013 Patita Co. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

#define N 10000
#define length 1
#define x1 0
#define x2 1
#define T 0.1
#define gama 1.4
#define Cf 0.75
#define alpha 6
//initial conditions
#define rho1 1
#define u1 0
#define p1 1

#define rho2 0.125
#define u2 0
#define p2 0.1


double U[N+1][3],U_p[N+1][3],dU_dT[N+1][3],F[N+1][3],dU_dX[N+1][3],dF_dT[N+1][3];

// Define the Courant Condition
double CFL (double U[N+1][3])
{
    int n;
    double max_vel, vel, p, u,dx;
    max_vel=1.0e-4;
    dx= double(x2-x1)/double(N);
    for (n=1;n<=N-1;n++)
    {
        u=U[n][1]/U[n][0];
        p=(gama-1)*(U[n][2]-0.5*U[n][0]*u*u);
        vel = sqrt(gama*p/U[n][0])+fabs(u);
        if (vel>max_vel)max_vel=vel;
    }
    return Cf*dx/max_vel;
 }

// Define W_function in function of alpha
double Wfunct(double x_plus, double x_min)
{
    double W;
    W=(pow(fabs(x_plus),alpha)*x_min+pow(fabs(x_min),alpha)*x_plus)/(pow(fabs(x_plus),alpha)+pow(fabs(x_min),alpha)+1e-20);
    return W;
}

// Define the initial Values
void Initial_Values(double U[N+1][3], double U_p[N+1][3], double dU_dX[N+1][3])
{
    int n, k;
    double dU_dX_plus, dU_dX_min, dx;
    dx= double(x2-x1)/double(N);
    for (n=0; n<N+1; n++) {
        // get the first section of the shock tube
        if (n<=N/2) {
            U[n][0]=rho1;
            U[n][1]=rho1*u1;
            U[n][2]=p1/(gama-1)+0.5*(rho1*u1*u1);
        }
        // get the second section of the shock tube
        else
        {
            U[n][0]=rho2;
            U[n][1]=rho2*u2;
            U[n][2]=p2/(gama-1)+0.5*(rho2*u2*u2);
        }
    }
    for (n=0; n<=N; n++) {
        for (k=0; k<3; k++) {
            U_p[n][k]=U[n][k];
        }
    }
    // get intial dudx
    for (n=0; n<N; n++) {
        for (k=0; k<3; k++) {
            dU_dX_plus=2*(U_p[n+1][k]-U[n][k])/dx;
            dU_dX_min=2*(U[n][k]-U_p[n][k])/dx;
            dU_dX[n][k]=Wfunct(dU_dX_plus,dU_dX_min);
        }
    }
    for(k=0;k<3;k++)
		dU_dX[N][k]=2*(U[N][k]-U_p[N][k])/dx;
}

// Define the boundary conditions
void Boundaries (double U[N+1][3])
{
    int k;
	for(k=0;k<3;k++)
	{
		U[0][k]=U[1][k];
		U[N][k]=U[N-1][k];
	}
}
// Asign F from U // PROBLEMS

void get_F(double U[N+1][3], double F[N+1][3])
{
    int n;
    double u, p;
    for (n=0; n<N+1; n++)
    {
        u=U[n][1]/U[n][0];
        p=(gama-1)*(U[n][2]-0.5*U[n][1]*U[n][1]/U[n][0]);
        F[n][0]=U[n][1];
        F[n][1]=U[n][0]*u*u+p;
        F[n][2]=(U[n][2]+p)*u;
    }
}
// Get Dudt, dfdT
void get_dUdT_dFdT (double U[N+1][3],double dU_dX[N+1][3],double dU_dT[N+1][3],double dF_dT[N+1][3])
{
    int n;
	double u;
	for(n=0;n<=N;n++)
	{
		u=U[n][1]/U[n][0];
		
		dU_dT[n][0]=-dU_dX[n][1];
		dU_dT[n][1]=-(gama-3)*u*u*dU_dX[n][0]/2-(3-gama)*u*dU_dX[n][1]-(gama-1)*dU_dX[n][2];
		dU_dT[n][2]=-((gama-1)*u*u*u-gama*u*U[n][2]/U[n][0])*dU_dX[n][0]-
        (gama*U[n][2]/U[n][0]-3*(gama-1)/2*u*u)*dU_dX[n][1]-gama*u*dU_dX[n][2];
        
		dF_dT[n][0]=dU_dT[n][1];
		dF_dT[n][1]=(gama-3)*u*u*dU_dT[n][0]/2+(3-gama)*u*dU_dT[n][1]+(gama-1)*dU_dT[n][2];
		dF_dT[n][2]=((gama-1)*u*u*u-gama*u*U[n][2]/U[n][0])*dU_dT[n][0]+
        (gama*U[n][2]/U[n][0]-3*(gama-1)/2*u*u)*dU_dT[n][1]+gama*u*dU_dT[n][2];
         }
}

// Get dudx for complete time step
void dU_dX_n(double U[N+1][3],double U_p[N+1][3],double dU_dX[N+1][3],double dU_dT[N+1][3],double dt)
{
	int n,k;
	double dU_dX_plus,dU_dX_min,dx;
	dx=(double(x2-x1))/double(N);
	for(n=1;n<N;n++)
	{
		for(k=0;k<3;k++)
		{
			dU_dX_plus=2*(U_p[n][k]+0.5*dt*dU_dT[n][k]-U[n][k])/dx;
			dU_dX_min=2*(U[n][k]-U_p[n-1][k]-0.5*dt*dU_dT[n-1][k])/dx;
			dU_dX[n][k]=Wfunct(dU_dX_plus, dU_dX_min);
       
		}
	}
}

// Get dudx for half time step

void dU_dX_p(double U[N+1][3],double U_p[N+1][3],double dU_dX[N+1][3],double dU_dT[N+1][3],double dt)
{
	int n,k;
	double dU_dX_plus,dU_dX_min,dx;
	dx=(double(x2-x1))/double(N);
	for(n=0;n<N;n++)
	{
		for(k=0;k<3;k++)
		{
			dU_dX_plus=2*(U[n+1][k]+0.5*dt*dU_dT[n+1][k]-U_p[n][k])/dx;
			dU_dX_min=2*(U_p[n][k]-U[n][k]-0.5*dt*dU_dT[n][k])/dx;
			dU_dX[n][k]=Wfunct(dU_dX_plus, dU_dX_min);
		}
	}
}
void CE_SE(double U[N+1][3],double U_p[N+1][3])
{
	int n,k;
	double dt,t,dx,V1,V2;
	dx=(double(x2-x1))/double(N);
	t=0;
	Initial_Values(U, U_p, dU_dX);
	while(t<T)
	{
		dt=CFL(U);
		t=t+dt;
		get_F(U,F);
		get_dUdT_dFdT(U,dU_dX,dU_dT,dF_dT);
		for(n=0;n<N;n++)
		{
			for(k=0;k<3;k++)
			{
                V1=dx*dU_dX[n][k]/4+dt*F[n][k]/dx+dt*dt*dF_dT[n][k]/(4*dx);
                V2=dx*dU_dX[n+1][k]/4+dt*F[n+1][k]/dx+dt*dt*dF_dT[n+1][k]/(4*dx);
                U_p[n][k]=0.5*(U[n][k]+U[n+1][k]+V1-V2);
			}
		}
		get_F(U_p,F);
		dU_dX_p(U,U_p,dU_dX,dU_dT,dt);
		get_dUdT_dFdT(U_p,dU_dX,dU_dT,dF_dT);
		for(n=1;n<N;n++)
		{
			for(k=0;k<3;k++)
			{
                V1=dx*dU_dX[n-1][k]/4+dt*F[n-1][k]/dx+dt*dt*dF_dT[n-1][k]/(4*dx);
                V2=dx*dU_dX[n][k]/4+dt*F[n][k]/dx+dt*dt*dF_dT[n][k]/(4*dx);
                U[n][k]=0.5*(U_p[n-1][k]+U_p[n][k]+V1-V2);
			}
		}
		Boundaries(U);
		dU_dX_n(U,U_p,dU_dX,dU_dT,dt);
	}
}
// Write the output file
void Output(double U[N+1][3])
{
	int n;
	double x,dx;
    dx=(double(x2-x1))/double(N);
    FILE *fp;
	double rho,u,p;
	fp=fopen("result.dat","w+");
    fprintf(fp,"%2.30s\n %20.60s\n %20.18s\t %2.3d\t %2.18s\t ","TITLE = \"1D-EULER.dat \"","variables = \"x\", \"rho\", \"u\", \"p\"","zone i=",N+1,"f=point\n");
    
    
	for(n=0;n<N+1;n++)
	{
		x=x1+n*dx;
		rho=U[n][0];
		u=U[n][1]/U[n][0];
		p=(gama-1)*(U[n][2]-0.5*U[n][0]*u*u);
		fprintf(fp,"%20.5f\t%20.5f\t%20.5f\t%20.5f\n",x,rho,u,p);
	}
	fclose(fp);
}

// IMPLEMENT THE CODE

int main()
{
    clock_t start,finish;
	double duration;
	start=clock();
    CE_SE(U, U_p);
    Output(U);
    finish=clock();
	duration=(double)(finish-start)/CLOCKS_PER_SEC;
    cout << "Complete in \n" << duration << "s \n";
    return 0;
}