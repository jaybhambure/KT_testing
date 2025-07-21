using namespace std;

	void Heun_x(double *E, double *M_x, double *M_y)
	{
	
		/*** Heun's method ****/

		double	E_1[Vol] = {}, E_2[Vol] = {};

		double M_x_1[Vol] = {}, M_x_2[Vol] = {};	
		double M_y_1[Vol] = {}, M_y_2[Vol] = {};	
		
		double res_E, res_M_x, res_M_y;

	int cell_p, cell_pp, cell_m, cell_mm; // cell_p: cell_plus (or to the right of current cell), cell_m: cell_minus (or to the left of current cell)

	int flag = 0;	

	if (regulator == 1)
	regulate_half_cell_values(E, M_x, M_y, 1, 0);	

	//Step 1
	for (int L_coord = 0; L_coord < Vol; L_coord++)
		{

		get_cells_p(L_coord, &cell_p, 1, 0);
		get_cells_p(cell_p, &cell_pp, 1, 0);
		get_cells_m(L_coord, &cell_m, 1, 0);
		get_cells_m(cell_m, &cell_mm, 1, 0);

	res_E = F_E_x(E[L_coord], E[cell_p], E[cell_pp], E[cell_m], E[cell_mm], M_x[L_coord], M_x[cell_p], M_x[cell_pp], M_x[cell_m], M_x[cell_mm], M_y[L_coord], M_y[cell_p], M_y[cell_pp], M_y[cell_m], M_y[cell_mm]);
	res_M_x = F_M_xx(E[L_coord], E[cell_p], E[cell_pp], E[cell_m], E[cell_mm], M_x[L_coord], M_x[cell_p], M_x[cell_pp], M_x[cell_m], M_x[cell_mm], M_y[L_coord], M_y[cell_p], M_y[cell_pp], M_y[cell_m], M_y[cell_mm]);
	res_M_y = F_M_xy(E[L_coord], E[cell_p], E[cell_pp], E[cell_m], E[cell_mm], M_x[L_coord], M_x[cell_p], M_x[cell_pp], M_x[cell_m], M_x[cell_mm], M_y[L_coord], M_y[cell_p], M_y[cell_pp], M_y[cell_m], M_y[cell_mm]);

	E_1[L_coord] = E[L_coord] + dt/dx * res_E;
	M_x_1[L_coord] = M_x[L_coord] + dt/dx * res_M_x;
	M_y_1[L_coord] = M_y[L_coord] + dt/dx * res_M_y;
	

		}

	if (regulator == 1)
	{
	regulate_central_cell_values(E_1, M_x_1, M_y_1, 1, 0);
	regulate_half_cell_values(E_1, M_x_1, M_y_1, 1, 0);
	}

	//Step 2
	for (int L_coord = 0; L_coord < Vol; L_coord++)
		{
		
		get_cells_p(L_coord, &cell_p, 1, 0);
		get_cells_p(cell_p, &cell_pp, 1, 0);
		get_cells_m(L_coord, &cell_m, 1, 0);
		get_cells_m(cell_m, &cell_mm, 1, 0);


res_E = F_E_x(E_1[L_coord], E_1[cell_p], E_1[cell_pp], E_1[cell_m], E_1[cell_mm], M_x_1[L_coord], M_x_1[cell_p], M_x_1[cell_pp], M_x_1[cell_m], M_x_1[cell_mm], M_y_1[L_coord], M_y_1[cell_p], M_y_1[cell_pp], M_y_1[cell_m], M_y_1[cell_mm]);
res_M_x = F_M_xx(E_1[L_coord], E_1[cell_p], E_1[cell_pp], E_1[cell_m], E_1[cell_mm], M_x_1[L_coord], M_x_1[cell_p], M_x_1[cell_pp], M_x_1[cell_m], M_x_1[cell_mm], M_y_1[L_coord], M_y_1[cell_p], M_y_1[cell_pp], M_y_1[cell_m], M_y_1[cell_mm]);
res_M_y = F_M_xy(E_1[L_coord], E_1[cell_p], E_1[cell_pp], E_1[cell_m], E_1[cell_mm], M_x_1[L_coord], M_x_1[cell_p], M_x_1[cell_pp], M_x_1[cell_m], M_x_1[cell_mm], M_y_1[L_coord], M_y_1[cell_p], M_y_1[cell_pp], M_y_1[cell_m], M_y_1[cell_mm]);

	E_2[L_coord] = 1.0/2.0*E[L_coord] + 1.0/2.0*( E_1[L_coord] + dt/dx * res_E );
	M_x_2[L_coord] = 1.0/2.0*M_x[L_coord] + 1.0/2.0*( M_x_1[L_coord] + dt/dx * res_M_x );
	M_y_2[L_coord] = 1.0/2.0*M_y[L_coord] + 1.0/2.0*( M_y_1[L_coord] + dt/dx * res_M_x );
	
		}
		
		
	if (regulator == 1)
	{
	regulate_central_cell_values(E_2, M_x_2, M_y_2, 1, 0);
	regulate_half_cell_values(E_2, M_x_2, M_y_2, 1, 0);
	}

	for (int L_coord = 0; L_coord < Vol; L_coord++)
	{
	E[L_coord] = E_2[L_coord];
	M_x[L_coord] = M_x_2[L_coord];
	M_y[L_coord] = M_y_2[L_coord];
	}

	
	
	}


