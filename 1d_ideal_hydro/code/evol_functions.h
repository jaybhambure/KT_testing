#include "neighbour_cells_1d.h"

using namespace std;
	
double minmod(double a, double b, double c) 
{
    if ((a > 0.0 && b > 0.0 && c > 0.0)) 
    {
	return std::min(std::min(a, b), c);
    } 

    else if ((a < 0.0 && b < 0.0 && c < 0.0)) 
    {
	  return  std::max(std::max(a, b), c);
    }

    else 
    {
        return 0.0;
    }
}


	double maximum(double a, double b, double c, double d)
	{
	
	return std::max(std::max(a,b), std::max(c,d));
	
	}




		double deriv(double rho_i, double rho_ip, double rho_im)
		{
		double forward_deriv = (rho_ip - rho_i)/dx;
		double backward_deriv = (rho_i - rho_im)/dx;
		double centred_deriv = (rho_ip - rho_im)/2.0/dx;

		return ( minmod(theta * forward_deriv, centred_deriv, theta * backward_deriv) );
		}
	
		

		/********************************* Kurganov Tadmor scheme **********************************/

		/****************** F_E: computes flux of energy density *******************/

		double F_E(double E_ix, double E_ixp, double E_ixpp, double E_ixm, double E_ixmm, double M_x_ix, double M_x_ixp, double M_x_ixpp, double M_x_ixm, double M_x_ixmm)
		{

		double H_i_plus_half;
		double H_i_minus_half;
		double H;
		
		/********** Flux of energy density: J = M_x *************/

		double J_ix = M_x_ix;
		double J_ixp = M_x_ixp;
		double J_ixpp = M_x_ixpp;
		double J_ixm = M_x_ixm;
		double J_ixmm = M_x_ixmm;


		double a_i_plus_half;
		double a_i_minus_half;
		
		/************** compute E, M_x at edge of cells **************/	
		
		double E_plus_i_plus_half = E_ixp - dx/2.0 * deriv(E_ixp, E_ixpp, E_ix);
		double E_minus_i_plus_half = E_ix + dx/2.0 * deriv(E_ix, E_ixp, E_ixm);
		
		double E_plus_i_minus_half = E_ix - dx/2.0 * deriv(E_ix, E_ixp, E_ixm);
		double E_minus_i_minus_half = E_ixm + dx/2.0 * deriv(E_ixm, E_ix, E_ixmm);
		
		double Mx_plus_i_plus_half = M_x_ixp - dx/2.0 * deriv(M_x_ixp, M_x_ixpp, M_x_ix);
		double Mx_minus_i_plus_half = M_x_ix + dx/2.0 * deriv(M_x_ix, M_x_ixp, M_x_ixm);
		
		double Mx_plus_i_minus_half = M_x_ix - dx/2.0 * deriv(M_x_ix, M_x_ixp, M_x_ixm);
		double Mx_minus_i_minus_half = M_x_ixm + dx/2.0 * deriv(M_x_ixm, M_x_ix, M_x_ixmm);


		/*************** to compute maximum propagation speed ******************************/

	       	double p;
		double v_plus_i_plus_half, v_minus_i_plus_half;
		double v_plus_i_minus_half, v_minus_i_minus_half;

		double a_plus_i_plus_half, a_minus_i_plus_half;
		double a_plus_i_minus_half, a_minus_i_minus_half;

		p = pressure(E_plus_i_plus_half, Mx_plus_i_plus_half);
		v_plus_i_plus_half = Mx_plus_i_plus_half/(E_plus_i_plus_half + p);

		p = pressure(E_minus_i_plus_half, Mx_minus_i_plus_half);
		v_minus_i_plus_half = Mx_minus_i_plus_half/(E_minus_i_plus_half + p);
	
		p = pressure(E_plus_i_minus_half, Mx_plus_i_minus_half);
		v_plus_i_minus_half = Mx_plus_i_minus_half/(E_plus_i_minus_half + p);
			
		p = pressure(E_minus_i_minus_half, Mx_minus_i_minus_half);
		v_minus_i_minus_half = Mx_minus_i_minus_half/(E_minus_i_minus_half + p);

		double a_plus_1 = (v_plus_i_plus_half + cs)/(1.0 + v_plus_i_plus_half * cs);
		double a_plus_2 = (v_plus_i_plus_half - cs)/(1.0 - v_plus_i_plus_half * cs);
		
		double a_minus_1 = (v_minus_i_plus_half + cs)/(1.0 + v_minus_i_plus_half * cs);
		double a_minus_2 = (v_minus_i_plus_half - cs)/(1.0 - v_minus_i_plus_half * cs);

		a_i_plus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		a_plus_1 = (v_plus_i_minus_half + cs)/(1.0 + v_plus_i_minus_half * cs);
		a_plus_2 = (v_plus_i_minus_half - cs)/(1.0 - v_plus_i_minus_half * cs);
		
		a_minus_1 = (v_minus_i_minus_half + cs)/(1.0 + v_minus_i_minus_half * cs);
		a_minus_2 = (v_minus_i_minus_half - cs)/(1.0 - v_minus_i_minus_half * cs);

		a_i_minus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		//a_i_plus_half = 1.0;
		//a_i_minus_half = 1.0;

		double J_plus_i_plus_half = J_ixp - dx/2.0 * deriv(J_ixp, J_ixpp, J_ix);
		double J_minus_i_plus_half = J_ix + dx/2.0 * deriv(J_ix, J_ixp, J_ixm);
		
		double J_plus_i_minus_half = J_ix - dx/2.0 * deriv(J_ix, J_ixp, J_ixm);
		double J_minus_i_minus_half = J_ixm + dx/2.0 * deriv(J_ixm, J_ix, J_ixmm);

		H_i_plus_half = 0.5 * (J_plus_i_plus_half + J_minus_i_plus_half) - a_i_plus_half/2.0 * (E_plus_i_plus_half - E_minus_i_plus_half);
		H_i_minus_half = 0.5 * (J_plus_i_minus_half + J_minus_i_minus_half) - a_i_minus_half/2.0 * (E_plus_i_minus_half - E_minus_i_minus_half);

		double result = - (H_i_plus_half - H_i_minus_half);

		return result;

		}		
	
		/********************************************/





		
		/****************** F_M_x: computes flux of momentum density *******************/

		double F_M_x(double E_ix, double E_ixp, double E_ixpp, double E_ixm, double E_ixmm, double M_x_ix, double M_x_ixp, double M_x_ixpp, double M_x_ixm, double M_x_ixmm)
		{
		double H_i_plus_half;
		double H_i_minus_half;
		double H;
		
		/******* Flux of momentum density J = M_x^2/(E + p) + p ***************/

		double J_ix = pressure(E_ix, M_x_ix) + M_x_ix * M_x_ix/(E_ix + pressure(E_ix, M_x_ix)); 
		double J_ixp = pressure(E_ixp, M_x_ixp) + M_x_ixp * M_x_ixp/(E_ixp + pressure(E_ixp, M_x_ixp));
		double J_ixpp = pressure(E_ixpp, M_x_ixpp) + M_x_ixpp * M_x_ixpp/(E_ixpp + pressure(E_ixpp, M_x_ixpp));
		double J_ixm = pressure(E_ixm, M_x_ixm) + M_x_ixm * M_x_ixm/(E_ixm + pressure(E_ixm, M_x_ixm));
		double J_ixmm = pressure(E_ixmm, M_x_ixmm) + M_x_ixmm * M_x_ixmm/(E_ixmm + pressure(E_ixmm, M_x_ixmm));

		double a_i_plus_half;
		double a_i_minus_half;

		/************** compute E, M_x, and J at edge of cells **************/	

		double E_plus_i_plus_half = E_ixp - dx/2.0 * deriv(E_ixp, E_ixpp, E_ix);
		double E_minus_i_plus_half = E_ix + dx/2.0 * deriv(E_ix, E_ixp, E_ixm);
		
		double E_plus_i_minus_half = E_ix - dx/2.0 * deriv(E_ix, E_ixp, E_ixm);
		double E_minus_i_minus_half = E_ixm + dx/2.0 * deriv(E_ixm, E_ix, E_ixmm);
	
		double Mx_plus_i_plus_half = M_x_ixp - dx/2.0 * deriv(M_x_ixp, M_x_ixpp, M_x_ix);
		double Mx_minus_i_plus_half = M_x_ix + dx/2.0 * deriv(M_x_ix, M_x_ixp, M_x_ixm);
		
		double Mx_plus_i_minus_half = M_x_ix - dx/2.0 * deriv(M_x_ix, M_x_ixp, M_x_ixm);
		double Mx_minus_i_minus_half = M_x_ixm + dx/2.0 * deriv(M_x_ixm, M_x_ix, M_x_ixmm);
		
		double J_plus_i_plus_half = J_ixp - dx/2.0 * deriv(J_ixp, J_ixpp, J_ix);
		double J_minus_i_plus_half = J_ix + dx/2.0 * deriv(J_ix, J_ixp, J_ixm);
		
		double J_plus_i_minus_half = J_ix - dx/2.0 * deriv(J_ix, J_ixp, J_ixm);
		double J_minus_i_minus_half = J_ixm + dx/2.0 * deriv(J_ixm, J_ix, J_ixmm);


		/*************** to compute maximum propagation speed ******************************/

	       	double p;
		double v_plus_i_plus_half, v_minus_i_plus_half;
		double v_plus_i_minus_half, v_minus_i_minus_half;

		double a_plus_i_plus_half, a_minus_i_plus_half;
		double a_plus_i_minus_half, a_minus_i_minus_half;

		p = pressure(E_plus_i_plus_half, Mx_plus_i_plus_half);
		v_plus_i_plus_half = Mx_plus_i_plus_half/(E_plus_i_plus_half + p);

		p = pressure(E_minus_i_plus_half, Mx_minus_i_plus_half);
		v_minus_i_plus_half = Mx_minus_i_plus_half/(E_minus_i_plus_half + p);
	
		p = pressure(E_plus_i_minus_half, Mx_plus_i_minus_half);
		v_plus_i_minus_half = Mx_plus_i_minus_half/(E_plus_i_minus_half + p);
			
		p = pressure(E_minus_i_minus_half, Mx_minus_i_minus_half);
		v_minus_i_minus_half = Mx_minus_i_minus_half/(E_minus_i_minus_half + p);

		double a_plus_1 = (v_plus_i_plus_half + cs)/(1.0 + v_plus_i_plus_half * cs);
		double a_plus_2 = (v_plus_i_plus_half - cs)/(1.0 - v_plus_i_plus_half * cs);
		
		double a_minus_1 = (v_minus_i_plus_half + cs)/(1.0 + v_minus_i_plus_half * cs);
		double a_minus_2 = (v_minus_i_plus_half - cs)/(1.0 - v_minus_i_plus_half * cs);

		a_i_plus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		a_plus_1 = (v_plus_i_minus_half + cs)/(1.0 + v_plus_i_minus_half * cs);
		a_plus_2 = (v_plus_i_minus_half - cs)/(1.0 - v_plus_i_minus_half * cs);
		
		a_minus_1 = (v_minus_i_minus_half + cs)/(1.0 + v_minus_i_minus_half * cs);
		a_minus_2 = (v_minus_i_minus_half - cs)/(1.0 - v_minus_i_minus_half * cs);

		a_i_minus_half = maximum(abs(a_plus_1), abs(a_plus_2), abs(a_minus_1), abs(a_minus_2));

		//a_i_plus_half = 1.0;
		//a_i_minus_half = 1.0;

		H_i_plus_half = 0.5 * (J_plus_i_plus_half + J_minus_i_plus_half) - a_i_plus_half/2.0 * (Mx_plus_i_plus_half - Mx_minus_i_plus_half);
		H_i_minus_half = 0.5 * (J_plus_i_minus_half + J_minus_i_minus_half) - a_i_minus_half/2.0 * (Mx_plus_i_minus_half - Mx_minus_i_minus_half);

		double result = - (H_i_plus_half - H_i_minus_half);

		return result;

		}		
	
		/********************************************/


