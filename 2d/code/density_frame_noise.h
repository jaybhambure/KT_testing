#include "random_LRF.h"
#include "matrix_multiply.h"
using namespace std;

double eta_lower_lower[dim][dim]; // the Minkowski metric

double vx, vy;
double vel[d];
double gamma_v;
double Lambda[dim][dim];

double u_upper[dim]; // u^mu
double u_lower[dim]; // u_mu

double l_x[dim], l_y[dim];
double LL[d][d];
double PL;

double Delta_lower_lower[dim][dim]; // Delta with 2 lower indices
double Delta_upper_lower[dim][dim]; // Delta with 1 upper and 1 lower indices

double P_lower_lower[dim][dim];
double lPl_lower_lower[dim][dim];

double kappa_xx[dim][dim];
double kappa_xy[dim][dim];
double kappa_yx[dim][dim];
double kappa_yy[dim][dim];

double cs2 = 1.0/(double)d;

double Kronecker_delta(int a, int b)
{
if (a == b)
	return (1.0);
else
	return (0.0);
}
	

	void metric()	
	{
		for (int mu = 0; mu < dim; mu++)
		{
			for (int nu = 0; nu < dim; nu++)
			{
			    	eta_lower_lower[mu][nu] = pow(-1, Kronecker_delta(0, mu)) * Kronecker_delta(mu, nu);
			}
		}
	}



	void initialize_matrices(double vx, double vy)
	{
	double v = sqrt(vx*vx + vy*vy);
	
	gamma_v = 1.0/sqrt(1.0 - v*v);
	vel[0] = vx; vel[1] = vy;
	
	u_upper[0] = gamma_v; u_upper[1] = gamma_v * vx; u_upper[2] = gamma_v * vy;
	u_lower[0] = - gamma_v; u_lower[1] = gamma_v * vx; u_lower[2] = gamma_v * vy;
	
	PL = 1.0 - cs2 * v*v;

/******** write (2+1)-dimensional Lorentz matrix ************/

Lambda[0][0] = gamma_v; 
Lambda[0][1] = gamma_v * vx; 
Lambda[0][2] = gamma_v * vy;

Lambda[1][0] = gamma_v * vx; 
Lambda[1][1] = 1.0 + (gamma_v * gamma_v)/(1.0 + gamma_v) * vx * vx; 
Lambda[1][2] = (gamma_v * gamma_v)/(1.0 + gamma_v) * vx * vy;

Lambda[2][0] = gamma_v * vy; 
Lambda[2][1] = (gamma_v*gamma_v)/(1.0 + gamma_v) * vx * vy; 
Lambda[2][2] = 1.0 + gamma_v*gamma_v/(1.0 + gamma_v) * vy * vy;

/**********************************************************************/



	/***** initialize l^i_mu ***************/

	for (int mu = 0; mu < dim; mu++)
	{
	l_x[mu] = Kronecker_delta(1, mu) - vx * Kronecker_delta(0, mu); // 0 means t, 1 means x
	l_y[mu] = Kronecker_delta(2, mu) - vy * Kronecker_delta(0, mu); // 0 means t, 2 means y
	}
	
	/***** initialize l^i_mu ***************/




	/*********** initialize Lij *****************/

	for (int i = 0; i < d; i++) // here 0 means x, 1 means 1
	{
	for (int j = 0; j < d; j++)
	{
	 	LL[i][j] = Kronecker_delta(i, j) - vel[i] * vel[j];
	}
	}

	/*********** initialize Lij *****************/


				
	/************ compute Delta_mu_nu ******************/

	for (int mu = 0; mu < dim; mu++)
	{
	for (int nu = 0; nu < dim; nu++)
	{
	Delta_lower_lower[mu][nu] = eta_lower_lower[mu][nu] + u_lower[mu] * u_lower[nu];
	}
	}


	/************ compute Delta^mu_nu ******************/

	for (int mu = 0; mu < dim; mu++)
	{
	for (int nu = 0; nu < dim; nu++)
	{
	Delta_upper_lower[mu][nu] = Kronecker_delta(mu, nu) + u_upper[mu] * u_lower[nu];
	}
	}


	
	/*********** Initialize projector P_mu_nu ******************/

	for (int mu = 0; mu < dim; mu++)
	{
	for (int nu = 0; nu < dim; nu++)
	{
	P_lower_lower[mu][nu] = - cs2 * u_lower[mu] * u_lower[nu] + 1.0/(double)d * Delta_lower_lower[mu][nu];
	}
	}
	
	/*********** Initialize projector P_mu_nu ******************/



	/************ initialize lPl_mu_nu *************/


	for (int mu = 0; mu < dim; mu++)
	{
	for (int nu = 0; nu < dim; nu++)
	{
	lPl_lower_lower[mu][nu] = - cs2/gamma_v/gamma_v * Delta_upper_lower[0][mu] * Delta_upper_lower[0][nu] + 1.0/(double)d * Delta_lower_lower[mu][nu];
	}
	}

	/************ initialize lPl_mu_nu *************/


 /*************************** all ingredients present ****************************************/


	/*************** initialize kappa matrices *****************/


	for (int mu = 0; mu < dim; mu++)
	{
	for (int nu = 0; nu < dim; nu++)
	{
	kappa_xx[mu][nu] = l_x[mu] * l_x[nu] - LL[0][0]/PL * lPl_lower_lower[mu][nu] + LL[0][0]/PL * P_lower_lower[mu][nu];
	kappa_xy[mu][nu] = l_x[mu] * l_y[nu] - LL[0][1]/PL * lPl_lower_lower[mu][nu] + LL[0][1]/PL * P_lower_lower[mu][nu];
	kappa_yx[mu][nu] = l_y[mu] * l_x[nu] - LL[1][0]/PL * lPl_lower_lower[mu][nu] + LL[1][0]/PL * P_lower_lower[mu][nu];
	kappa_yy[mu][nu] = l_y[mu] * l_y[nu] - LL[1][1]/PL * lPl_lower_lower[mu][nu] + LL[1][1]/PL * P_lower_lower[mu][nu];
	}
	}
	
	/*************** initialize kappa matrices *****************/

}


	double lambda_LRF_xx, lambda_LRF_xy, lambda_LRF_yy;

	double lambda_matrix[dim][dim];

	double lambda_mu_nu[dim][dim];

	

	void boost_noise_matrix()
	{
	
	/************* Noise matrix adding a 0 row and column ****************/

	lambda_matrix[0][0] = 0.0;
	lambda_matrix[0][1] = 0.0;
	lambda_matrix[0][2] = 0.0;

	lambda_matrix[1][0] = 0.0;
	lambda_matrix[1][1] = lambda_LRF_xx;
	lambda_matrix[1][2] = lambda_LRF_xy;

	lambda_matrix[2][0] = 0.0;
	lambda_matrix[2][1] = lambda_LRF_xy;
	lambda_matrix[2][2] = lambda_LRF_yy;

	/****************** Lorentz boost matrix ******************/

	double product[dim][dim] = {};

	double get_result[dim*dim];

	multiply(Lambda, lambda_matrix, get_result); // computes Lambda . noise. Lambda^T

	for (int L_coord = 0; L_coord < dim*dim; L_coord++)
	{
	int j = L_coord/dim;
	int i = L_coord - j*dim;
	lambda_mu_nu[i][j] = get_result[L_coord];
	}

	}
		


void get_noise(double T, double vx, double vy, double *xi_xx, double *xi_xy, double *xi_yx, double *xi_yy)
{

metric();

double A = 4.0 *eta*T/dt/dx/dy; // strength of noise

//double A = 4.0 *eta*T; // strength of noise

double v = sqrt(vx*vx + vy*vy);
double phi = atan(vy/vx);

/*************************************************************************************************/

	initialize_matrices(vx, vy);

	random_numbers_LRF(&lambda_LRF_xx, &lambda_LRF_xy, &lambda_LRF_yy, A); // construct noise in LRF, calls function in header file random_LRF.h

	boost_noise_matrix(); // boost noise
	
	double xi_xx_temp = 0.0, xi_xy_temp = 0.0, xi_yx_temp = 0.0, xi_yy_temp = 0.0; // noise in density frame (temporary variables)

		for (int mu = 0; mu < dim; mu++)
		{
		for (int nu = 0; nu < dim; nu++)
		{
		xi_xx_temp = xi_xx_temp + kappa_xx[mu][nu] * lambda_mu_nu[mu][nu];
		xi_xy_temp = xi_xy_temp + kappa_xy[mu][nu] * lambda_mu_nu[mu][nu];
		xi_yx_temp = xi_yx_temp + kappa_yx[mu][nu] * lambda_mu_nu[mu][nu];
		xi_yy_temp = xi_yy_temp + kappa_yy[mu][nu] * lambda_mu_nu[mu][nu];
		}
		}	

		*xi_xx = xi_xx_temp;
		*xi_xy = xi_xy_temp;
		*xi_yx = xi_yx_temp;
		*xi_yy = xi_yy_temp;

}


