using namespace std;

	void regulate_half_cell_values(double *E_temp, double *M_x_temp, double *M_y_temp, int e_x, int e_y)
	{
		
		int neighbour_p, neighbour_pp, neighbour_m, neighbour_mm;	

		int neighbour_plus, neighbour_minus;

		double avg_E, avg_Mx, avg_My;

		double da;

		if (e_x == 1)
		da = dx;

		else 
		da = dy;	

		for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
		get_cells_p(L_coord, &neighbour_p, e_x, e_y);
		get_cells_p(neighbour_p, &neighbour_pp, e_x, e_y);
		
		get_cells_m(L_coord, &neighbour_m, e_x, e_y);
		get_cells_m(neighbour_m, &neighbour_mm, e_x, e_y);
		
		double E_plus_cell_plus_half = E_temp[neighbour_p] - da/2.0 * deriv(E_temp[neighbour_p], E_temp[neighbour_pp], E_temp[L_coord], da);
		double E_minus_cell_plus_half = E_temp[L_coord] + da/2.0 * deriv(E_temp[L_coord], E_temp[neighbour_p], E_temp[neighbour_m], da);

		double E_plus_cell_minus_half = E_temp[L_coord] - da/2.0 * deriv(E_temp[L_coord], E_temp[neighbour_p], E_temp[neighbour_m], da);
		double E_minus_cell_minus_half = E_temp[neighbour_m] + da/2.0 * deriv(E_temp[neighbour_m], E_temp[L_coord], E_temp[neighbour_mm], da);
		
		double Mx_plus_cell_plus_half = M_x_temp[neighbour_p] - da/2.0 * deriv(M_x_temp[neighbour_p], M_x_temp[neighbour_pp], M_x_temp[L_coord], da);
		double Mx_minus_cell_plus_half = M_x_temp[L_coord] + da/2.0 * deriv(M_x_temp[L_coord], M_x_temp[neighbour_p], M_x_temp[neighbour_m], da);
	
		double Mx_plus_cell_minus_half = M_x_temp[L_coord] - da/2.0 * deriv(M_x_temp[L_coord], M_x_temp[neighbour_p], M_x_temp[neighbour_m], da);
		double Mx_minus_cell_minus_half = M_x_temp[neighbour_m] + da/2.0 * deriv(M_x_temp[neighbour_m], M_x_temp[L_coord], M_x_temp[neighbour_mm], da);
		
		double My_plus_cell_plus_half = M_y_temp[neighbour_p] - da/2.0 * deriv(M_y_temp[neighbour_p], M_y_temp[neighbour_pp], M_y_temp[L_coord], da);
		double My_minus_cell_plus_half = M_y_temp[L_coord] + da/2.0 * deriv(M_y_temp[L_coord], M_y_temp[neighbour_p], M_y_temp[neighbour_m], da);
	
		double My_plus_cell_minus_half = M_y_temp[L_coord] - da/2.0 * deriv(M_y_temp[L_coord], M_y_temp[neighbour_p], M_y_temp[neighbour_m], da);
		double My_minus_cell_minus_half = M_y_temp[neighbour_m] + da/2.0 * deriv(M_y_temp[neighbour_m], M_y_temp[L_coord], M_y_temp[neighbour_mm], da);


		double M_plus_cell_plus_half = sqrt(pow(Mx_plus_cell_plus_half, 2.0) + pow(My_plus_cell_plus_half, 2.0));
		double M_minus_cell_plus_half = sqrt(pow(Mx_minus_cell_plus_half, 2.0) + pow(My_minus_cell_plus_half, 2.0));
		double M_plus_cell_minus_half = sqrt(pow(Mx_plus_cell_minus_half, 2.0) + pow(My_plus_cell_minus_half, 2.0));
		double M_minus_cell_minus_half = sqrt(pow(Mx_minus_cell_minus_half, 2.0) + pow(My_minus_cell_minus_half, 2.0));


		if (E_plus_cell_plus_half < M_plus_cell_plus_half) 
		{
		//counter_regulation++;

		get_cells_p(neighbour_p, &neighbour_plus, e_x, e_y);
		get_cells_m(neighbour_p, &neighbour_minus, e_x, e_y);

		avg_E = 1.0/3.0 * (E_temp[neighbour_minus] + E_temp[neighbour_p] + E_temp[neighbour_plus]);	
		avg_Mx = 1.0/3.0 * (M_x_temp[neighbour_minus] + M_x_temp[neighbour_p] + M_x_temp[neighbour_plus]);	
		avg_My = 1.0/3.0 * (M_y_temp[neighbour_minus] + M_y_temp[neighbour_p] + M_y_temp[neighbour_plus]);	

		cout<<"check half: "<<E_temp[neighbour_minus]<<" "<<E_temp[neighbour_p]<<" "<<E_temp[neighbour_plus]<<" "<<avg_E<<" L_coord: "<<L_coord<<endl;

		E_temp[neighbour_minus] = avg_E;	
		E_temp[neighbour_p] = avg_E;	
		E_temp[neighbour_plus] = avg_E;	
		
		M_x_temp[neighbour_minus] = avg_Mx;	
		M_x_temp[neighbour_p] = avg_Mx;	
		M_x_temp[neighbour_plus] = avg_Mx;	

		M_y_temp[neighbour_minus] = avg_My;	
		M_y_temp[neighbour_p] = avg_My;	
		M_y_temp[neighbour_plus] = avg_My;	
		}		
		
		if (E_minus_cell_plus_half < M_minus_cell_plus_half) 
		{
		//counter_regulation++;

		get_cells_p(L_coord, &neighbour_plus, e_x, e_y);
		get_cells_m(L_coord, &neighbour_minus, e_x, e_y);

		avg_E = 1.0/3.0 * (E_temp[neighbour_minus] + E_temp[L_coord] + E_temp[neighbour_plus]);	
		avg_Mx = 1.0/3.0 * (M_x_temp[neighbour_minus] + M_x_temp[L_coord] + M_x_temp[neighbour_plus]);	
		avg_My = 1.0/3.0 * (M_y_temp[neighbour_minus] + M_y_temp[L_coord] + M_y_temp[neighbour_plus]);	

		cout<<"check half: "<<E_temp[neighbour_minus]<<" "<<E_temp[L_coord]<<" "<<E_temp[neighbour_plus]<<" "<<avg_E<<" L_coord: "<<L_coord<<endl;
		
		E_temp[neighbour_minus] = avg_E;	
		E_temp[L_coord] = avg_E;	
		E_temp[neighbour_plus] = avg_E;	
	
		M_x_temp[neighbour_minus] = avg_Mx;	
		M_x_temp[L_coord] = avg_Mx;	
		M_x_temp[neighbour_plus] = avg_Mx;	
		
		M_y_temp[neighbour_minus] = avg_My;	
		M_y_temp[L_coord] = avg_My;	
		M_y_temp[neighbour_plus] = avg_My;	

		}	


		if (E_plus_cell_minus_half < M_plus_cell_minus_half) 
		{
		//counter_regulation++;
		
		get_cells_p(L_coord, &neighbour_plus, e_x, e_y);
		get_cells_m(L_coord, &neighbour_minus, e_x, e_y);

		avg_E = 1.0/3.0 * (E_temp[neighbour_minus] + E_temp[L_coord] + E_temp[neighbour_plus]);	
		avg_Mx = 1.0/3.0 * (M_x_temp[neighbour_minus] + M_x_temp[L_coord] + M_x_temp[neighbour_plus]);	
		avg_My = 1.0/3.0 * (M_y_temp[neighbour_minus] + M_y_temp[L_coord] + M_y_temp[neighbour_plus]);	

		E_temp[neighbour_minus] = avg_E;	
		E_temp[L_coord] = avg_E;	
		E_temp[neighbour_plus] = avg_E;	
		
		M_x_temp[neighbour_minus] = avg_Mx;	
		M_x_temp[L_coord] = avg_Mx;	
		M_x_temp[neighbour_plus] = avg_Mx;	

		M_y_temp[neighbour_minus] = avg_My;	
		M_y_temp[L_coord] = avg_My;	
		M_y_temp[neighbour_plus] = avg_My;	

		}	


		if (E_minus_cell_minus_half < M_minus_cell_minus_half) 
		{
		//counter_regulation++;
		
		get_cells_p(neighbour_m, &neighbour_plus, e_x, e_y);
		get_cells_m(neighbour_m, &neighbour_minus, e_x, e_y);

		avg_E = 1.0/3.0 * (E_temp[neighbour_minus] + E_temp[neighbour_m] + E_temp[neighbour_plus]);	
		avg_Mx = 1.0/3.0 * (M_x_temp[neighbour_minus] + M_x_temp[neighbour_m] + M_x_temp[neighbour_plus]);	
		avg_My = 1.0/3.0 * (M_y_temp[neighbour_minus] + M_y_temp[neighbour_m] + M_y_temp[neighbour_plus]);	

		E_temp[neighbour_minus] = avg_E;	
		E_temp[neighbour_m] = avg_E;	
		E_temp[neighbour_plus] = avg_E;	
		
		M_x_temp[neighbour_minus] = avg_Mx;	
		M_x_temp[neighbour_m] = avg_Mx;	
		M_x_temp[neighbour_plus] = avg_Mx;	

		M_y_temp[neighbour_minus] = avg_My;	
		M_y_temp[neighbour_m] = avg_My;	
		M_y_temp[neighbour_plus] = avg_My;	

		}	


		}

	}




	void regulate_central_cell_values(double *E_temp, double *M_x_temp, double *M_y_temp, int e_x, int e_y)
	{
		int neighbour_p, neighbour_m;

		double M_magnitude;

		for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
			M_magnitude = sqrt(pow(M_x_temp[L_coord], 2.0) + pow(M_y_temp[L_coord], 2.0)); 

			if (E_temp[L_coord] < M_magnitude)
			{

				//counter_regulation++;

				get_cells_p(L_coord, &neighbour_p, e_x, e_y);	
				get_cells_m(L_coord, &neighbour_m, e_x, e_y);	
		
			double avg_E = 1.0/3.0 * (E_temp[neighbour_m] + E_temp[L_coord] + E_temp[neighbour_p]);
			double avg_Mx = 1.0/3.0 * (M_x_temp[neighbour_m] + M_x_temp[L_coord] + M_x_temp[neighbour_p]);
			double avg_My = 1.0/3.0 * (M_y_temp[neighbour_m] + M_y_temp[L_coord] + M_y_temp[neighbour_p]);

			cout<<"E: "<<E_temp[neighbour_m]<<" "<<E_temp[L_coord]<<" "<<E_temp[neighbour_p]<<" "<<avg_E<<endl;

			E_temp[neighbour_m] = avg_E;		
			E_temp[L_coord] = avg_E;		
			E_temp[neighbour_p] = avg_E;		
			
			M_x_temp[neighbour_m] = avg_Mx;		
			M_x_temp[L_coord] = avg_Mx;		
			M_x_temp[neighbour_p] = avg_Mx;
			
			M_y_temp[neighbour_m] = avg_My;		
			M_y_temp[L_coord] = avg_My;		
			M_y_temp[neighbour_p] = avg_My;

			}

		//cout<<"L_coord: "<<L_coord<<" E1: "<<E_1[L_coord]<<" Mx1: "<<M_x_1[L_coord]<<endl;
		}

	}


