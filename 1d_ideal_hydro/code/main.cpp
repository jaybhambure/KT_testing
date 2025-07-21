#include<iostream>
#include<fstream>
#include<cmath>
#include "macros.h"
#include "pressure_entropy.h"
#include "RK3_evol.h"

using namespace std;

int main()
{

global_N = 0;

double E[Vol];
double M_x[Vol];

int coeff = N/10; // N is total time steps, coeff is for picking out times when to write to file


/********************* INITIAL CONDITIONS ************************/

/*************** From D.Teaney et al's draft, test 1 *****************/

/*
for (int ix = 0; ix < L; ix++)
{
	double x = ix * dx;
	E[ix] = 0.48 * exp(-(x - L*dx/2.0)*(x - L*dx/2.0)/5.0/5.0) + 0.12;
	M_x[ix] = 0.0;
}	
*/

/******************************************************************/



/*************** From D.Teaney et al's draft, test 2 *****************/

/*
for (int ix = 0; ix < L; ix++)
{
	double x = ix * dx;
	E[ix] = 9.6 * exp(-(x - L*dx/2.0)*(x - L*dx/2.0)/5.0/5.0) + 0.06;
	M_x[ix] = 0.0;
}
*/

/******************************************************************/



/**********************************************************/


/*************** From D.Teaney et al's draft, test 3 *****************/

/*
for (int ix = 0; ix < L; ix++)
{
	
double x = (ix - L/2) *dx; 

if (x <= 0.0) 
	E[ix] = 0.48;  // units GeV^-4 
else
	E[ix] = 0.12; 

M_x[ix] = 0.0;

}
*/

/**********************************************************/



/*************** From D. Teaney et al's draft, test 4 *****************/

for (int ix = 0; ix < L; ix++)
{
	
double x = (ix - L/2) *dx; 

if (x <= 0.0) 
	E[ix] = 4.8;  // units GeV^-4 
else
	E[ix] = 0.06; 

M_x[ix] = 0.0;

}


/**********************************************************/


ofstream file;

for (int it = 0; it <= N; it++)
{
double sum_E = 0.0;
double sum_M_x = 0.0;

global_N++; // this variable tracks global time; used in RK3_evol.h to output time at which pressure becomes negative (if it happens)

for (int ix = 0; ix < L; ix++)
{
sum_E = sum_E + E[ix]*dx; // checks energy momentum conservation
sum_M_x = sum_M_x + M_x[ix]*dx;
}

cout<<"it: "<<it<<" sum_E: "<<sum_E<<" sum_M_x: "<<sum_M_x<<endl;

	double v;

	if (it%coeff == 0)
	{
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test1/E_M_"+std::to_string(it)+".txt"); 
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test2/E_M_"+std::to_string(it)+".txt"); 
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test3/E_M_"+std::to_string(it)+".txt");
	file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test4/E_M_"+std::to_string(it)+".txt");
	
	for (int ix = 0; ix<L; ix++)
	{
	v = vel(E[ix], M_x[ix]);
	double gamma = 1.0/sqrt(1.0 - v*v);	
	file<<it*dt<<" "<<(ix - L/2.0)*dx<<" "<<E[ix]<<" "<<M_x[ix]<<" "<<v<<" "<<gamma * v<<endl;
	}
	
	}

	file.close();
	
	/*
	if (it*dt == 47.5) // this particular time was used by D.T. et al
	{
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test1/E_M_"+std::to_string(it)+".txt");
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test2/E_M_"+std::to_string(it)+".txt");
	
	for (int ix = 0; ix<L; ix++)
	{
	v = vel(E[ix], M_x[ix]);
	double gamma = 1.0/sqrt(1.0 - v*v);	
	file<<it*dt<<" "<<(ix - L/2.0)*dx<<" "<<E[ix]<<" "<<M_x[ix]<<" "<<v<<" "<<gamma * v<<endl;
	}

	}

	file.close();
	*/	

	
	if (it*dt == 61.85) // this particular time was used by D.T. et al
	{
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test3/E_M_"+std::to_string(it)+".txt");
	file.open("/home/chandrodoy/Dropbox/Density_frame_codes/1d_ideal_hydro/test4/E_M_"+std::to_string(it)+".txt"); 
	
	for (int ix = 0; ix<L; ix++)
	{
	v = vel(E[ix], M_x[ix]);
	double gamma = 1.0/sqrt(1.0 - v*v);	
	file<<it*dt<<" "<<(ix - L/2.0)*dx<<" "<<E[ix]<<" "<<M_x[ix]<<" "<<v<<" "<<gamma * v<<endl;
	}

	}

	file.close();
	

	/*****************************  RK3 Shu and Osher ***************************/
	
	RK3(E, M_x); // Runge_Kutta evolution
		
	/*****************************  RK3 Shu and Osher ***************************/

cout<<"Time: "<<it*dt<<endl;

}



}


