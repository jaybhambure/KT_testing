using namespace std;
	
		/***** Kurganov-Tadmor scheme *****/

		/****************** F_E_x *******************/

		double F_E_x(double E_cell, double E_cell_p, double E_cell_pp, double E_cell_m, double E_cell_mm, double M_x_cell, double M_x_cell_p, double M_x_cell_pp, double M_x_cell_m, double M_x_cell_mm, double M_y_cell, double M_y_cell_p, double M_y_cell_pp, double M_y_cell_m, double M_y_cell_mm)
		{

		double H_cell_plus_half;
		double H_cell_minus_half;
		
		double J_cell = M_x_cell;
		double J_cell_p = M_x_cell_p;
		double J_cell_pp = M_x_cell_pp;
		double J_cell_m = M_x_cell_m;
		double J_cell_mm = M_x_cell_mm;

		double a_cell_plus_half;
		double a_cell_minus_half;
		
		double E_plus_cell_plus_half = E_cell_p - dx/2.0 * deriv(E_cell_p, E_cell_pp, E_cell, dx);
		double E_minus_cell_plus_half = E_cell + dx/2.0 * deriv(E_cell, E_cell_p, E_cell_m, dx);
		
		double E_plus_cell_minus_half = E_cell - dx/2.0 * deriv(E_cell, E_cell_p, E_cell_m, dx);
		double E_minus_cell_minus_half = E_cell_m + dx/2.0 * deriv(E_cell_m, E_cell, E_cell_mm, dx);
		
		double Mx_plus_cell_plus_half = M_x_cell_p - dx/2.0 * deriv(M_x_cell_p, M_x_cell_pp, M_x_cell, dx);
		double Mx_minus_cell_plus_half = M_x_cell + dx/2.0 * deriv(M_x_cell, M_x_cell_p, M_x_cell_m, dx);
		
		double Mx_plus_cell_minus_half = M_x_cell - dx/2.0 * deriv(M_x_cell, M_x_cell_p, M_x_cell_m, dx);
		double Mx_minus_cell_minus_half = M_x_cell_m + dx/2.0 * deriv(M_x_cell_m, M_x_cell, M_x_cell_mm, dx);
		
		double My_plus_cell_plus_half = M_y_cell_p - dx/2.0 * deriv(M_y_cell_p, M_y_cell_pp, M_y_cell, dx);
		double My_minus_cell_plus_half = M_y_cell + dx/2.0 * deriv(M_y_cell, M_y_cell_p, M_y_cell_m, dx);
		
		double My_plus_cell_minus_half = M_y_cell - dx/2.0 * deriv(M_y_cell, M_y_cell_p, M_y_cell_m, dx);
		double My_minus_cell_minus_half = M_y_cell_m + dx/2.0 * deriv(M_y_cell_m, M_y_cell, M_y_cell_mm, dx);

		/*************** to compute maximum propagation speed ******************************/

	       	double p;
		double v_plus_cell_plus_half, v_minus_cell_plus_half;
		double v_plus_cell_minus_half, v_minus_cell_minus_half;

		double a_plus_cell_plus_half, a_minus_cell_plus_half;
		double a_plus_cell_minus_half, a_minus_cell_minus_half;

		p = pressure(E_plus_cell_plus_half, Mx_plus_cell_plus_half, My_plus_cell_plus_half);
		v_plus_cell_plus_half = Mx_plus_cell_plus_half/(E_plus_cell_plus_half + p);

		p = pressure(E_minus_cell_plus_half, Mx_minus_cell_plus_half, My_minus_cell_plus_half);
		v_minus_cell_plus_half = Mx_minus_cell_plus_half/(E_minus_cell_plus_half + p);
	
		p = pressure(E_plus_cell_minus_half, Mx_plus_cell_minus_half, My_plus_cell_minus_half);
		v_plus_cell_minus_half = Mx_plus_cell_minus_half/(E_plus_cell_minus_half + p);
			
		p = pressure(E_minus_cell_minus_half, Mx_minus_cell_minus_half, My_minus_cell_minus_half);
		v_minus_cell_minus_half = Mx_minus_cell_minus_half/(E_minus_cell_minus_half + p);

		double a_plus_1 = (v_plus_cell_plus_half + cs)/(1.0 + v_plus_cell_plus_half * cs);
		double a_plus_2 = (v_plus_cell_plus_half - cs)/(1.0 - v_plus_cell_plus_half * cs);
		
		double a_minus_1 = (v_minus_cell_plus_half + cs)/(1.0 + v_minus_cell_plus_half * cs);
		double a_minus_2 = (v_minus_cell_plus_half - cs)/(1.0 - v_minus_cell_plus_half * cs);

		a_cell_plus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		a_plus_1 = (v_plus_cell_minus_half + cs)/(1.0 + v_plus_cell_minus_half * cs);
		a_plus_2 = (v_plus_cell_minus_half - cs)/(1.0 - v_plus_cell_minus_half * cs);
		
		a_minus_1 = (v_minus_cell_minus_half + cs)/(1.0 + v_minus_cell_minus_half * cs);
		a_minus_2 = (v_minus_cell_minus_half - cs)/(1.0 - v_minus_cell_minus_half * cs);

		a_cell_minus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		//a_cell_plus_half = 1.0;
		//a_cell_minus_half = 1.0;

		double J_plus_cell_plus_half = J_cell_p - dx/2.0 * deriv(J_cell_p, J_cell_pp, J_cell, dx);
		double J_minus_cell_plus_half = J_cell + dx/2.0 * deriv(J_cell, J_cell_p, J_cell_m, dx);
		
		double J_plus_cell_minus_half = J_cell - dx/2.0 * deriv(J_cell, J_cell_p, J_cell_m, dx);
		double J_minus_cell_minus_half = J_cell_m + dx/2.0 * deriv(J_cell_m, J_cell, J_cell_mm, dx);

		H_cell_plus_half = 0.5 * (J_plus_cell_plus_half + J_minus_cell_plus_half) - a_cell_plus_half/2.0 * (E_plus_cell_plus_half - E_minus_cell_plus_half);
		H_cell_minus_half = 0.5 * (J_plus_cell_minus_half + J_minus_cell_minus_half) - a_cell_minus_half/2.0 * (E_plus_cell_minus_half - E_minus_cell_minus_half);

		double result = - (H_cell_plus_half - H_cell_minus_half);

		return result;

		}		
	
		/********************************************/

		
		/****************** F_M_xx *******************/

		double F_M_xx(double E_cell, double E_cell_p, double E_cell_pp, double E_cell_m, double E_cell_mm, double M_x_cell, double M_x_cell_p, double M_x_cell_pp, double M_x_cell_m, double M_x_cell_mm, double M_y_cell, double M_y_cell_p, double M_y_cell_pp, double M_y_cell_m, double M_y_cell_mm)
		{
		double H_cell_plus_half;
		double H_cell_minus_half;
		
		double J_cell = pressure(E_cell, M_x_cell, M_y_cell) + M_x_cell * M_x_cell/(E_cell + pressure(E_cell, M_x_cell, M_y_cell));
		double J_cell_p = pressure(E_cell_p, M_x_cell_p, M_y_cell_p) + M_x_cell_p * M_x_cell_p/(E_cell_p + pressure(E_cell_p, M_x_cell_p, M_y_cell_p));
		double J_cell_pp = pressure(E_cell_pp, M_x_cell_pp, M_y_cell_pp) + M_x_cell_pp * M_x_cell_pp/(E_cell_pp + pressure(E_cell_pp, M_x_cell_pp, M_y_cell_pp));
		double J_cell_m = pressure(E_cell_m, M_x_cell_m, M_y_cell_m) + M_x_cell_m * M_x_cell_m/(E_cell_m + pressure(E_cell_m, M_x_cell_m, M_y_cell_m));
		double J_cell_mm = pressure(E_cell_mm, M_x_cell_mm, M_y_cell_mm) + M_x_cell_mm * M_x_cell_mm/(E_cell_mm + pressure(E_cell_mm, M_x_cell_mm, M_y_cell_mm));

		double a_cell_plus_half;
		double a_cell_minus_half;

		/************** compute E and M_x at edge of cells **************/	

		double E_plus_cell_plus_half = E_cell_p - dx/2.0 * deriv(E_cell_p, E_cell_pp, E_cell, dx);
		double E_minus_cell_plus_half = E_cell + dx/2.0 * deriv(E_cell, E_cell_p, E_cell_m, dx);
		
		double E_plus_cell_minus_half = E_cell - dx/2.0 * deriv(E_cell, E_cell_p, E_cell_m, dx);
		double E_minus_cell_minus_half = E_cell_m + dx/2.0 * deriv(E_cell_m, E_cell, E_cell_mm, dx);
	
		double Mx_plus_cell_plus_half = M_x_cell_p - dx/2.0 * deriv(M_x_cell_p, M_x_cell_pp, M_x_cell, dx);
		double Mx_minus_cell_plus_half = M_x_cell + dx/2.0 * deriv(M_x_cell, M_x_cell_p, M_x_cell_m, dx);
		
		double Mx_plus_cell_minus_half = M_x_cell - dx/2.0 * deriv(M_x_cell, M_x_cell_p, M_x_cell_m, dx);
		double Mx_minus_cell_minus_half = M_x_cell_m + dx/2.0 * deriv(M_x_cell_m, M_x_cell, M_x_cell_mm, dx);
		
		double My_plus_cell_plus_half = M_y_cell_p - dx/2.0 * deriv(M_y_cell_p, M_y_cell_pp, M_y_cell, dx);
		double My_minus_cell_plus_half = M_y_cell + dx/2.0 * deriv(M_y_cell, M_y_cell_p, M_y_cell_m, dx);
		
		double My_plus_cell_minus_half = M_y_cell - dx/2.0 * deriv(M_y_cell, M_y_cell_p, M_y_cell_m, dx);
		double My_minus_cell_minus_half = M_y_cell_m + dx/2.0 * deriv(M_y_cell_m, M_y_cell, M_y_cell_mm, dx);
		
		double J_plus_cell_plus_half = J_cell_p - dx/2.0 * deriv(J_cell_p, J_cell_pp, J_cell, dx);
		double J_minus_cell_plus_half = J_cell + dx/2.0 * deriv(J_cell, J_cell_p, J_cell_m, dx);
		
		double J_plus_cell_minus_half = J_cell - dx/2.0 * deriv(J_cell, J_cell_p, J_cell_m, dx);
		double J_minus_cell_minus_half = J_cell_m + dx/2.0 * deriv(J_cell_m, J_cell, J_cell_mm, dx);


		/*************** to compute maximum propagation speed ******************************/

	       	double p;
		double v_plus_cell_plus_half, v_minus_cell_plus_half;
		double v_plus_cell_minus_half, v_minus_cell_minus_half;

		double a_plus_cell_plus_half, a_minus_cell_plus_half;
		double a_plus_cell_minus_half, a_minus_cell_minus_half;

		p = pressure(E_plus_cell_plus_half, Mx_plus_cell_plus_half, My_plus_cell_plus_half);
		v_plus_cell_plus_half = Mx_plus_cell_plus_half/(E_plus_cell_plus_half + p);

		p = pressure(E_minus_cell_plus_half, Mx_minus_cell_plus_half, My_minus_cell_plus_half);
		v_minus_cell_plus_half = Mx_minus_cell_plus_half/(E_minus_cell_plus_half + p);
	
		p = pressure(E_plus_cell_minus_half, Mx_plus_cell_minus_half, My_plus_cell_minus_half);
		v_plus_cell_minus_half = Mx_plus_cell_minus_half/(E_plus_cell_minus_half + p);
			
		p = pressure(E_minus_cell_minus_half, Mx_minus_cell_minus_half, My_minus_cell_minus_half);
		v_minus_cell_minus_half = Mx_minus_cell_minus_half/(E_minus_cell_minus_half + p);

		double a_plus_1 = (v_plus_cell_plus_half + cs)/(1.0 + v_plus_cell_plus_half * cs);
		double a_plus_2 = (v_plus_cell_plus_half - cs)/(1.0 - v_plus_cell_plus_half * cs);
		
		double a_minus_1 = (v_minus_cell_plus_half + cs)/(1.0 + v_minus_cell_plus_half * cs);
		double a_minus_2 = (v_minus_cell_plus_half - cs)/(1.0 - v_minus_cell_plus_half * cs);

		a_cell_plus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		a_plus_1 = (v_plus_cell_minus_half + cs)/(1.0 + v_plus_cell_minus_half * cs);
		a_plus_2 = (v_plus_cell_minus_half - cs)/(1.0 - v_plus_cell_minus_half * cs);
		
		a_minus_1 = (v_minus_cell_minus_half + cs)/(1.0 + v_minus_cell_minus_half * cs);
		a_minus_2 = (v_minus_cell_minus_half - cs)/(1.0 - v_minus_cell_minus_half * cs);

		a_cell_minus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		//a_cell_plus_half = 1.0;
		//a_cell_minus_half = 1.0;

		H_cell_plus_half = 0.5 * (J_plus_cell_plus_half + J_minus_cell_plus_half) - a_cell_plus_half/2.0 * (Mx_plus_cell_plus_half - Mx_minus_cell_plus_half);
		H_cell_minus_half = 0.5 * (J_plus_cell_minus_half + J_minus_cell_minus_half) - a_cell_minus_half/2.0 * (Mx_plus_cell_minus_half - Mx_minus_cell_minus_half);

		double result = - (H_cell_plus_half - H_cell_minus_half);

		return result;

		}		
	
		/********************************************/



		double F_M_xy(double E_cell, double E_cell_p, double E_cell_pp, double E_cell_m, double E_cell_mm, double M_x_cell, double M_x_cell_p, double M_x_cell_pp, double M_x_cell_m, double M_x_cell_mm, double M_y_cell, double M_y_cell_p, double M_y_cell_pp, double M_y_cell_m, double M_y_cell_mm)
		{
		double H_cell_plus_half;
		double H_cell_minus_half;
		
		double J_cell = M_x_cell * M_y_cell/(E_cell + pressure(E_cell, M_x_cell, M_y_cell));
		double J_cell_p = M_x_cell_p * M_y_cell_p/(E_cell_p + pressure(E_cell_p, M_x_cell_p, M_y_cell_p));
		double J_cell_pp = M_x_cell_pp * M_y_cell_pp/(E_cell_pp + pressure(E_cell_pp, M_x_cell_pp, M_y_cell_pp));
		double J_cell_m = M_x_cell_m * M_y_cell_m/(E_cell_m + pressure(E_cell_m, M_x_cell_m, M_y_cell_m));
		double J_cell_mm = M_x_cell_mm * M_y_cell_mm/(E_cell_mm + pressure(E_cell_mm, M_x_cell_mm, M_y_cell_mm));

		double a_cell_plus_half;
		double a_cell_minus_half;

		/************** compute E, M_x, and M_y at edge of cells **************/	

		double E_plus_cell_plus_half = E_cell_p - dx/2.0 * deriv(E_cell_p, E_cell_pp, E_cell, dx);
		double E_minus_cell_plus_half = E_cell + dx/2.0 * deriv(E_cell, E_cell_p, E_cell_m, dx);
		
		double E_plus_cell_minus_half = E_cell - dx/2.0 * deriv(E_cell, E_cell_p, E_cell_m, dx);
		double E_minus_cell_minus_half = E_cell_m + dx/2.0 * deriv(E_cell_m, E_cell, E_cell_mm, dx);
	
		double Mx_plus_cell_plus_half = M_x_cell_p - dx/2.0 * deriv(M_x_cell_p, M_x_cell_pp, M_x_cell, dx);
		double Mx_minus_cell_plus_half = M_x_cell + dx/2.0 * deriv(M_x_cell, M_x_cell_p, M_x_cell_m, dx);
		
		double Mx_plus_cell_minus_half = M_x_cell - dx/2.0 * deriv(M_x_cell, M_x_cell_p, M_x_cell_m, dx);
		double Mx_minus_cell_minus_half = M_x_cell_m + dx/2.0 * deriv(M_x_cell_m, M_x_cell, M_x_cell_mm, dx);
		
		double My_plus_cell_plus_half = M_y_cell_p - dx/2.0 * deriv(M_y_cell_p, M_y_cell_pp, M_y_cell, dx);
		double My_minus_cell_plus_half = M_y_cell + dx/2.0 * deriv(M_y_cell, M_y_cell_p, M_y_cell_m, dx);
		
		double My_plus_cell_minus_half = M_y_cell - dx/2.0 * deriv(M_y_cell, M_y_cell_p, M_y_cell_m, dx);
		double My_minus_cell_minus_half = M_y_cell_m + dx/2.0 * deriv(M_y_cell_m, M_y_cell, M_y_cell_mm, dx);

		double J_plus_cell_plus_half = J_cell_p - dx/2.0 * deriv(J_cell_p, J_cell_pp, J_cell, dx);
		double J_minus_cell_plus_half = J_cell + dx/2.0 * deriv(J_cell, J_cell_p, J_cell_m, dx);
		
		double J_plus_cell_minus_half = J_cell - dx/2.0 * deriv(J_cell, J_cell_p, J_cell_m, dx);
		double J_minus_cell_minus_half = J_cell_m + dx/2.0 * deriv(J_cell_m, J_cell, J_cell_mm, dx);


		/*************** to compute maximum propagation speed ******************************/

	       	double p;
		double v_plus_cell_plus_half, v_minus_cell_plus_half;
		double v_plus_cell_minus_half, v_minus_cell_minus_half;

		double a_plus_cell_plus_half, a_minus_cell_plus_half;
		double a_plus_cell_minus_half, a_minus_cell_minus_half;

		p = pressure(E_plus_cell_plus_half, Mx_plus_cell_plus_half, My_plus_cell_plus_half);
		v_plus_cell_plus_half = Mx_plus_cell_plus_half/(E_plus_cell_plus_half + p);

		p = pressure(E_minus_cell_plus_half, Mx_minus_cell_plus_half, My_plus_cell_plus_half);
		v_minus_cell_plus_half = Mx_minus_cell_plus_half/(E_minus_cell_plus_half + p);
	
		p = pressure(E_plus_cell_minus_half, Mx_plus_cell_minus_half, My_plus_cell_minus_half);
		v_plus_cell_minus_half = Mx_plus_cell_minus_half/(E_plus_cell_minus_half + p);
			
		p = pressure(E_minus_cell_minus_half, Mx_minus_cell_minus_half, My_minus_cell_minus_half);
		v_minus_cell_minus_half = Mx_minus_cell_minus_half/(E_minus_cell_minus_half + p);

		double a_plus_1 = (v_plus_cell_plus_half + cs)/(1.0 + v_plus_cell_plus_half * cs);
		double a_plus_2 = (v_plus_cell_plus_half - cs)/(1.0 - v_plus_cell_plus_half * cs);
		
		double a_minus_1 = (v_minus_cell_plus_half + cs)/(1.0 + v_minus_cell_plus_half * cs);
		double a_minus_2 = (v_minus_cell_plus_half - cs)/(1.0 - v_minus_cell_plus_half * cs);

		a_cell_plus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		a_plus_1 = (v_plus_cell_minus_half + cs)/(1.0 + v_plus_cell_minus_half * cs);
		a_plus_2 = (v_plus_cell_minus_half - cs)/(1.0 - v_plus_cell_minus_half * cs);
		
		a_minus_1 = (v_minus_cell_minus_half + cs)/(1.0 + v_minus_cell_minus_half * cs);
		a_minus_2 = (v_minus_cell_minus_half - cs)/(1.0 - v_minus_cell_minus_half * cs);

		a_cell_minus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		//a_cell_plus_half = 1.0;
		//a_cell_minus_half = 1.0;

		H_cell_plus_half = 0.5 * (J_plus_cell_plus_half + J_minus_cell_plus_half) - a_cell_plus_half/2.0 * (My_plus_cell_plus_half - My_minus_cell_plus_half);
		H_cell_minus_half = 0.5 * (J_plus_cell_minus_half + J_minus_cell_minus_half) - a_cell_minus_half/2.0 * (My_plus_cell_minus_half - My_minus_cell_minus_half);

		double result = - (H_cell_plus_half - H_cell_minus_half);

		return result;

		}		
	
		/********************************************/


