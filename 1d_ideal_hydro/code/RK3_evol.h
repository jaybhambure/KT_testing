#include "evol_functions.h"

using namespace std;

	void RK3(double *E, double *M_x)
	{
		/****** Method of Shu & Osher *********/	

		double	E_1[Vol] = {}, E_2[Vol] = {}, E_3[Vol] = {};

		double M_x_1[Vol] = {}, M_x_2[Vol] = {}, M_x_3[Vol] = {};	
		
		double res_E, res_M_x;

	int ix, ixp, ixpp, ixm, ixmm;

	int flag = 0;	

	//Step 1
	for (int L_coord = 0; L_coord < Vol; L_coord++)
		{

		ix  = L_coord;	
			
		get_cells_p(ix, &ixp); // gets neighbouring cells (imposes periodic boundary conditions)
		get_cells_p(ixp, &ixpp);
		get_cells_m(ix, &ixm);
		get_cells_m(ixm, &ixmm);

	res_E = F_E(E[ix], E[ixp], E[ixpp], E[ixm], E[ixmm], M_x[ix], M_x[ixp], M_x[ixpp], M_x[ixm], M_x[ixmm]);
	res_M_x = F_M_x(E[ix], E[ixp], E[ixpp], E[ixm], E[ixmm], M_x[ix], M_x[ixp], M_x[ixpp], M_x[ixm], M_x[ixmm]);

	E_1[L_coord] = E[L_coord] + dt/dx * res_E;
	M_x_1[L_coord] = M_x[L_coord] + dt/dx * res_M_x;
	
	double p = pressure(E_1[L_coord], M_x_1[L_coord]);
		
	// checks if pressure is non-negative	
	if (p <= 0.0 || E_1[L_coord] < 0.0)
	{
	flag = 1;
	cout<<"Step 1, pressure negative, pressure: "<<p<<" E_1: "<<E_1[L_coord]<<" "<<L_coord<<" res_E: "<<dt/dx * res_E<<" res_M: "<<dt/dx*res_M_x<<" Mx: "<<M_x_1[L_coord]<<" global_time: "<<global_N<<endl;
	}

		}


	//Step 2
	for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
		ix = L_coord;	

		get_cells_p(ix, &ixp);
		get_cells_p(ixp, &ixpp);
		get_cells_m(ix, &ixm);
		get_cells_m(ixm, &ixmm);

	res_E = F_E(E_1[ix], E_1[ixp], E_1[ixpp], E_1[ixm], E_1[ixmm], M_x_1[ix], M_x_1[ixp], M_x_1[ixpp], M_x_1[ixm], M_x_1[ixmm]);
	res_M_x = F_M_x(E_1[ix], E_1[ixp], E_1[ixpp], E_1[ixm], E_1[ixmm], M_x_1[ix], M_x_1[ixp], M_x_1[ixpp], M_x_1[ixm], M_x_1[ixmm]);

	E_2[L_coord] = 3.0/4.0*E[L_coord] + 1.0/4.0*( E_1[L_coord] + dt/dx * res_E );
	M_x_2[L_coord] = 3.0/4.0*M_x[L_coord] + 1.0/4.0*( M_x_1[L_coord] + dt/dx * res_M_x );
	
	double p = pressure(E_2[L_coord], M_x_2[L_coord]);
	
	if (p <= 0.0 || E_2[L_coord] < 0.0)	
	{
	flag = 1;
	cout<<"Step 2, presure negative"<<" pressure: "<<p<<" E_2: "<<E_2[L_coord]<<" "<<L_coord<<" res_E: "<<dt/dx * res_E<<" res_M: "<<dt/dx*res_M_x<<" Mx: "<<M_x_2[L_coord]<<" global_time: "<<global_N<<endl;
	}

		}
		
	//Step 3
	
	for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
	
		ix = L_coord;

		get_cells_p(ix, &ixp);
		get_cells_p(ixp, &ixpp);
		get_cells_m(ix, &ixm);
		get_cells_m(ixm, &ixmm);

	res_E = F_E(E_2[ix], E_2[ixp], E_2[ixpp], E_2[ixm], E_2[ixmm], M_x_2[ix], M_x_2[ixp], M_x_2[ixpp], M_x_2[ixm], M_x_2[ixmm]);
	res_M_x = F_M_x(E_2[ix], E_2[ixp], E_2[ixpp], E_2[ixm], E_2[ixmm], M_x_2[ix], M_x_2[ixp], M_x_2[ixpp], M_x_2[ixm], M_x_2[ixmm]);

	E_3[L_coord] = 1.0/3.0*E[L_coord] + 2.0/3.0*(E_2[L_coord] + dt/dx * res_E );
	M_x_3[L_coord] = 1.0/3.0*M_x[L_coord] + 2.0/3.0*(M_x_2[L_coord] + dt/dx * res_M_x );

	double p = pressure(E_3[L_coord], M_x_3[L_coord]);
	

		if (p <= 0.0 || E_3[L_coord] < 0.0)
		{
		flag = 1;	
		cout<<"Step 3, presure negative"<<" pressure: "<<p<<" E_3: "<<E_3[L_coord]<<" "<<L_coord<<" res_E: "<<dt/dx * res_E<<" res_M: "<<dt/dx*res_M_x<<" Mx: "<<M_x_3[L_coord]<<" global_time: "<<global_N<<endl;
		}

		}
		

	for (int L_coord = 0; L_coord < Vol; L_coord++)
	{
	E[L_coord] = E_3[L_coord];
	M_x[L_coord] = M_x_3[L_coord];
	}

	
	
	}


