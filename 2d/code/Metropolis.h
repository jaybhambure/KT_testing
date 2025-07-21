#include<random>
using namespace std;

double generateRandom_uniform() 
{
    static thread_local std::mt19937 generator1{ std::random_device{}() };
    std::uniform_real_distribution<double> distribution(0, 1);
    return distribution(generator1);
}



		void update_field(double *E, double *M_x, double *M_y, int L_coord, int update_M)
		{
		
		int neighbour_x, neighbour_y, neighbour_xy;
		
			get_cells_p(L_coord, &neighbour_x, 1, 0);	
			get_cells_p(L_coord, &neighbour_y, 0, 1);	
			get_cells_p(neighbour_x, &neighbour_xy, 0, 1);	
	
	double update_Mx_L_coord = M_x[L_coord];
	double update_Mx_neighbour_x = M_x[neighbour_x];
       	double update_Mx_neighbour_y = M_x[neighbour_y];
       	double update_Mx_neighbour_xy = M_x[neighbour_xy];
	
	double update_My_L_coord = M_y[L_coord];
	double update_My_neighbour_x = M_y[neighbour_x];
       	double update_My_neighbour_y = M_y[neighbour_y];
       	double update_My_neighbour_xy = M_y[neighbour_xy];

	double xi_xx, xi_xy, xi_yx, xi_yy; // noise in density frame, has units of [T * \eta/(dt * dV)]^0.5 = GeV^(d+1) where d = no. of spatial dimensions 

	double T_L_coord, T_neighbour_x, T_neighbour_y, T_neighbour_xy;
	double vx_L_coord, vx_neighbour_x, vx_neighbour_y, vx_neighbour_xy;
	double vy_L_coord, vy_neighbour_x, vy_neighbour_y, vy_neighbour_xy;
	
	get_T_v(E[L_coord], M_x[L_coord], M_y[L_coord], &T_L_coord, &vx_L_coord, &vy_L_coord);
	get_T_v(E[neighbour_x], M_x[neighbour_x], M_y[neighbour_x], &T_neighbour_x, &vx_neighbour_x, &vy_neighbour_x);
	get_T_v(E[neighbour_y], M_x[neighbour_y], M_y[neighbour_y], &T_neighbour_y, &vx_neighbour_y, &vy_neighbour_y);
	get_T_v(E[neighbour_xy], M_x[neighbour_xy], M_y[neighbour_xy], &T_neighbour_xy, &vx_neighbour_xy, &vy_neighbour_xy);

	double T, vx, vy;
	
	T = 1.0/4.0 * (T_L_coord + T_neighbour_x + T_neighbour_y + T_neighbour_xy);
	vx = 1.0/4.0 * (vx_L_coord + vx_neighbour_x + vx_neighbour_y + vx_neighbour_xy);
	vy = 1.0/4.0 * (vy_L_coord + vy_neighbour_x + vy_neighbour_y + vy_neighbour_xy);

	get_noise(T, vx, vy, &xi_xx, &xi_xy, &xi_yx, &xi_yy);
	//get_noise(T_L_coord, vx_L_coord, vy_L_coord, &xi_xx, &xi_xy, &xi_yx, &xi_yy);

	double delta_p_1_x = 0.0, delta_p_2_x = 0.0, delta_p_1_y = 0.0, delta_p_2_y = 0.0; // 1 denotes x-surface, 2 denotes y-surface, x, y denotes which component of momenta
	
	if (update_M == 1) // x-component of momenta
	{
	delta_p_1_x = dt * dy * xi_xx; // units of GeV
	delta_p_2_x = dt * dx * xi_yx; // flux of x component of momenta along direction perpendicular to y i.e., dx
	
	update_Mx_L_coord = M_x[L_coord] + 1.0/dV * ( - 0.5 * delta_p_1_x - 0.5 * delta_p_2_x ); // trial fluxes are divided by dV as Mx is momentum density
	update_Mx_neighbour_x = M_x[neighbour_x] + 1.0/dV * ( 0.5 * delta_p_1_x - 0.5 * delta_p_2_x );
	update_Mx_neighbour_y = M_x[neighbour_y] + 1.0/dV * ( - 0.5 * delta_p_1_x + 0.5 * delta_p_2_x );
	update_Mx_neighbour_xy = M_x[neighbour_xy] + 1.0/dV * ( 0.5 * delta_p_1_x + 0.5 * delta_p_2_x );
	}

	else if (update_M == 2) // y-component of momenta
	{
	delta_p_1_y = dt * dy * xi_xy;
	delta_p_2_y = dt * dx * xi_yy;
	
	update_My_L_coord = M_y[L_coord] + 1.0/dV * ( - 0.5 * delta_p_1_y - 0.5 * delta_p_2_y );
	update_My_neighbour_x = M_y[neighbour_x] + 1.0/dV * ( 0.5 * delta_p_1_y - 0.5 * delta_p_2_y );
	update_My_neighbour_y = M_y[neighbour_y] + 1.0/dV * ( - 0.5 * delta_p_1_y + 0.5 * delta_p_2_y );
	update_My_neighbour_xy = M_y[neighbour_xy] + 1.0/dV * ( 0.5 * delta_p_1_y + 0.5 * delta_p_2_y );
	}

	/****************** check whether this update gives unphysical pressure *********************/
		
	   /*** conformal fluid: p = 1.0/(2*d) * ( - (d-1) E + sqrt( (d+1)^2 E^2 - 4*d \vec{M}^2 ) ); **********/

	int flag = 0;
	
	double magnitude_M_L_coord = sqrt(update_Mx_L_coord * update_Mx_L_coord + update_My_L_coord * update_My_L_coord);
	double magnitude_M_neighbour_x = sqrt(update_Mx_neighbour_x * update_Mx_neighbour_x + update_My_neighbour_x * update_My_neighbour_x);
	double magnitude_M_neighbour_y = sqrt(update_Mx_neighbour_y * update_Mx_neighbour_y + update_My_neighbour_y * update_My_neighbour_y);
	double magnitude_M_neighbour_xy = sqrt(update_Mx_neighbour_xy * update_Mx_neighbour_xy + update_My_neighbour_xy * update_My_neighbour_xy);

	if (E[L_coord] <= magnitude_M_L_coord || E[neighbour_x] <= magnitude_M_neighbour_x || E[neighbour_y] <= magnitude_M_neighbour_y || E[neighbour_xy] <= magnitude_M_neighbour_xy)
		flag = 1;
	
	/********************************************************************************************/

	/******************* if the pressure is negative/zero do not update *********************************/

	if (flag == 1)
	{
	//counter_no_update++;

	update_Mx_L_coord = M_x[L_coord];
	update_Mx_neighbour_x = M_x[neighbour_x];
	update_Mx_neighbour_y = M_x[neighbour_y];
	update_Mx_neighbour_xy = M_x[neighbour_xy];

	update_My_L_coord = M_y[L_coord];
	update_My_neighbour_x = M_y[neighbour_x];
	update_My_neighbour_y = M_y[neighbour_y];
	update_My_neighbour_xy = M_y[neighbour_xy];
	}

	/*********************************************************************************************/



		double entropy_density_old = entropy_density(E[L_coord], M_x[L_coord], M_y[L_coord]) + entropy_density(E[neighbour_x], M_x[neighbour_x], M_y[neighbour_x]) 
					     + entropy_density(E[neighbour_y], M_x[neighbour_y], M_y[neighbour_y]) + entropy_density(E[neighbour_xy], M_x[neighbour_xy], M_y[neighbour_xy]);


		double entropy_density_new = entropy_density(E[L_coord], update_Mx_L_coord, update_My_L_coord) + entropy_density(E[neighbour_x], update_Mx_neighbour_x, update_My_neighbour_x) 
					     + entropy_density(E[neighbour_y], update_Mx_neighbour_y, update_My_neighbour_y) + entropy_density(E[neighbour_xy], update_Mx_neighbour_xy, update_My_neighbour_xy);

		double del_S = (entropy_density_new - entropy_density_old) * dx * dy;

			double probability = min(1.0, exp(del_S));
	
			double random = generateRandom_uniform();

		if (random <= probability)
		{
		M_x[L_coord] = update_Mx_L_coord;	
		M_x[neighbour_x] = update_Mx_neighbour_x;	
		M_x[neighbour_y] = update_Mx_neighbour_y;	
		M_x[neighbour_xy] = update_Mx_neighbour_xy;	

		M_y[L_coord] = update_My_L_coord;	
		M_y[neighbour_x] = update_My_neighbour_x;	
		M_y[neighbour_y] = update_My_neighbour_y;	
		M_y[neighbour_xy] = update_My_neighbour_xy;	

		}

		}



	void Von_Neumann(double *E, double *M_x, double *M_y)
	{

	
			for (int L_coord = 0; L_coord < Vol; L_coord++)
			{
				
				update_field(E, M_x, M_y, L_coord, 1); // 1 means update M_x
				update_field(E, M_x, M_y, L_coord, 2); // 2 means update M_y
				
			}


	}


