#include<iostream>
#include<math.h>
#include<fstream>
#include<cmath>
#include "macros.h"
#include "neighbour_cells_2d.h"
#include "pressure_entropy.h"
#include "deriv_minmod.h"
#include "regulate.h"
#include "RK3_evol_x.h"
#include "RK3_evol_y.h"
#include "density_frame_noise.h"
#include "Metropolis.h"

using namespace std;

int main()
{
double E[Vol];
double M_x[Vol];
double M_y[Vol];

double M_in = - 10.0, M_fin = 10.0;

double h = 0.01;
int size = (M_fin - M_in)/h;

long int count_x[Vol][size];
long int count_y[Vol][size];

long int total_x[Vol] = {};
long int total_y[Vol] = {};

double e_avg[Vol] = {};

	for (int L_coord = 0; L_coord < Vol; L_coord++)
	{
		for (int i = 0; i < size; i++)
		{
		count_x[L_coord][i] = 0;
		count_y[L_coord][i] = 0;
		}
	}

int coeff = N/10; // N is total time steps, coeff is for printing to file 

/**** initial conditions: 2d Riemann slab *******/
/*
for (int L_coord = 0; L_coord < Vol; L_coord++)
{
	int j = L_coord/L;
	int i = L_coord - j*L;
	
	double x = (i - L/2)*dx;
	
	double y = (j - L/2)*dy;

	if (abs(x) < 1.0 && abs(y)< 2.0 )
	{
	E[L_coord] = A;
	}
	else
	E[L_coord] = delta;	
	

	M_x[L_coord] = 0.0;
	M_y[L_coord] = 0.0;
}
*/
/**** initial conditions: 2d Riemann slab *******/


/************* initial conditions: fixed energy density, coldstart momentum densities ****/

for (int L_coord = 0; L_coord < Vol; L_coord++)
{
	E[L_coord] = 2.0;
	//E[L_coord] = 5.0;

	M_x[L_coord] = 0.0;
	M_y[L_coord] = 0.0;
}

/************* initial conditions: fixed energy density, coldstart momentum densities ****/


ofstream file;


for (int it = 0; it < N; it++)
{

	/***** check for energy, momentum conservation ******/	

	double sum_E = 0.0; 
	double sum_M_x = 0.0;
	double sum_M_y = 0.0;

	for (int ix = 0; ix < Vol; ix++)
	{
	sum_E = sum_E + E[ix];
	sum_M_x = sum_M_x + M_x[ix];
	sum_M_y = sum_M_y + M_y[ix];
	}

	cout<<"it: "<<it<<" sum_E: "<<sum_E<<" sum_M_x: "<<sum_M_x<<" sum_M_y: "<<sum_M_y<<endl;

	/***** check for energy, momentum conservation ******/	
	
	/*	
	if (it%coeff == 0)
	{
	file.open("/home/chandrodoy/Dropbox/Density_frame_codes/2d/results/it="+std::to_string(it)+".txt");	
	
		for (int ix = 0; ix<L; ix++)
		{
		file<<it*dt<<" "<<(ix - L/2)*dx<<" "<<E[ix]<<" "<<M_x[ix]<<" "<<M_y[ix]<<endl;
		}
	file.close();
	}

	if (it%coeff == 0)
	{
	file.open("/home/chandrodoy/Dropbox/Density_frame_codes/2d/results_contour/it="+std::to_string(it)+".txt");	
	
		for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
		int j = L_coord/L;
		int i = L_coord - j*L;	

		file<<it*dt<<" "<<(i - L/2)*dx<<" "<<(j - L/2)*dy<<" "<<E[L_coord]<<" "<<M_x[L_coord]<<" "<<M_y[L_coord]<<endl;
		}
	file.close();
	}
	*/

	/***************  Runge-Kutta update (RK3 Shu Osher) ***************************/
	
	if (it%2 == 0) // update along x, then y
	{
	RK3_x(E, M_x, M_y);
	RK3_y(E, M_x, M_y);
	}

	else if (it%2==1) // update along y, then x
	{
	RK3_y(E, M_x, M_y);
	RK3_x(E, M_x, M_y);
	}
	
	/***************  Runge-Kutta update (RK3 Shu Osher) ***************************/


	/************ Metropolis update ***********************/

	Von_Neumann(E, M_x, M_y);	
	
	/************ Metropolis update ***********************/
	
		/***** histogram M_x, M_y **********/	
		
		for (int L_coord = 0; L_coord < Vol; L_coord++)
		{

			int pos_x = (M_x[L_coord] - M_in)/h;
			int pos_y = (M_y[L_coord] - M_in)/h;

			count_x[L_coord][pos_x]++;	
			count_y[L_coord][pos_y]++;	
				
		}

		/***** histogram M_x, M_y **********/	
		
		/***** compute average LRF energy density **********/	
		
		for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
			double e = - (d-1.0)/2.0*E[L_coord] + 1.0/2.0 * sqrt( (d+1.0)*(d+1.0) * E[L_coord] * E[L_coord] - 4.0*d* M_x[L_coord]*M_x[L_coord] - 4.0*d* M_y[L_coord]*M_y[L_coord] );
			
			e_avg[L_coord]+= e/(double)N;	
				
		}

		/***** compute average LRF energy density **********/	

cout<<"Time: "<<it*dt<<endl;

}


	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/2d/histogram/L="+std::to_string(L)+"/E_2.txt");	
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/2d/histogram/L="+std::to_string(L)+"/E_5.txt");	
	
	file.open("/home/chandrodoy/Dropbox/Density_frame_codes/2d/histogram/L="+std::to_string(L)+"/E_2_dV_0.25.txt");	
	//file.open("/home/chandrodoy/Dropbox/Density_frame_codes/2d/histogram/L="+std::to_string(L)+"/E_5_dV_0.25.txt");	
	
	double e_avg_Vol = 0.0;

	for (int L_coord = 0; L_coord < Vol; L_coord++)
	{
	e_avg_Vol += e_avg[L_coord]/Vol;	
	cout<<L_coord<<" e: "<<e_avg[L_coord]<<endl;
	}

	cout<<"e_avg: "<<e_avg_Vol<<endl;

	for (int i = 0; i<size; i++)
	{
	file<<M_in + i * h;	

		for (int L_coord = 0; L_coord < L; L_coord++) // Vol or L depending on size
		{
		file<<" "<<count_x[L_coord][i]/(double)N/h;	
		}
		
		for (int L_coord = 0; L_coord < L; L_coord++)
		{
		file<<" "<<count_y[L_coord][i]/(double)N/h;	
		}

		file<<endl;	
	}
	
	file.close();


}


