#include<iostream>
#include<math.h>
#include<fstream>
#include<cmath>
#include "macros.h"
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
double L[d][d];
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

double delta(int a, int b)
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
			    	eta_lower_lower[mu][nu] = pow(-1, delta(0, mu)) * delta(mu, nu);
			}
		}
	}



	void initialize_matrices(double v, double phi)
	{
	vx = v * cos(phi), vy = v * sin(phi);
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
	l_x[mu] = delta(1, mu) - vx * delta(0, mu); // 0 means t, 1 means x
	l_y[mu] = delta(2, mu) - vy * delta(0, mu); // 0 means t, 2 means y
	}
	
	/***** initialize l^i_mu ***************/


	/*********** initialize Lij *****************/

	for (int i = 0; i < d; i++) // here 0 means x, 1 means 1
	{
	for (int j = 0; j < d; j++)
	{
	 	L[i][j] = delta(i, j) - vel[i]	* vel[j];
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
	Delta_upper_lower[mu][nu] = delta(mu, nu) + u_upper[mu] * u_lower[nu];
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
	kappa_xx[mu][nu] = l_x[mu] * l_x[nu] - L[0][0]/PL * lPl_lower_lower[mu][nu] + L[0][0]/PL * P_lower_lower[mu][nu];
	kappa_xy[mu][nu] = l_x[mu] * l_y[nu] - L[0][1]/PL * lPl_lower_lower[mu][nu] + L[0][1]/PL * P_lower_lower[mu][nu];
	kappa_yx[mu][nu] = l_y[mu] * l_x[nu] - L[1][0]/PL * lPl_lower_lower[mu][nu] + L[1][0]/PL * P_lower_lower[mu][nu];
	kappa_yy[mu][nu] = l_y[mu] * l_y[nu] - L[1][1]/PL * lPl_lower_lower[mu][nu] + L[1][1]/PL * P_lower_lower[mu][nu];
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
		


int main()
{

metric();

double A = 4.0 *eta*T/dt;

double phi_divide = 8.0/3.0;

double v = 0.0, phi = M_PI/(double)phi_divide;

ofstream file;
ofstream file2;
ofstream file3;
ofstream file4;


	file.open("./results/numerics/xi_compare_phi_"+std::to_string(phi_divide)+".txt");
	file2.open("./results/analytics/analytic_"+std::to_string(phi_divide)+".txt");
	file3.open("./results/Landau/analytic_"+std::to_string(phi_divide)+".txt");
	file4.open("./results/Landau/numerics_"+std::to_string(phi_divide)+".txt");

/*************************************************************************************************/

	int N = 1000000;

	int count = 0;

do
{

initialize_matrices(v, phi);

double corr_xx_xx = 0.0, corr_xy_xy = 0.0, corr_yy_yy = 0.0;  // self-correlations
double corr_xx_xy = 0.0, corr_xx_yy = 0.0, corr_xy_yy = 0.0; // cross-correlations

double corr_xx_xx_Landau = 0.0, corr_xy_xy_Landau = 0.0, corr_yy_yy_Landau = 0.0;  // self-correlations
double corr_xx_xy_Landau = 0.0, corr_xx_yy_Landau = 0.0, corr_xy_yy_Landau = 0.0; // cross-correlations

double noise_xx = 0.0;

	for (int i = 1; i <= N; i++)
	{
	random_numbers_LRF(&lambda_LRF_xx, &lambda_LRF_xy, &lambda_LRF_yy, A); // construct noise in LRF

	boost_noise_matrix(); // boost noise
	
	corr_xx_xx_Landau = corr_xx_xx_Landau + lambda_mu_nu[1][1] * lambda_mu_nu[1][1]/N; // noise in Landau frame
	corr_xy_xy_Landau = corr_xy_xy_Landau + lambda_mu_nu[1][2] * lambda_mu_nu[1][2]/N;
	corr_yy_yy_Landau = corr_yy_yy_Landau + lambda_mu_nu[2][2] * lambda_mu_nu[2][2]/N;
	corr_xx_xy_Landau = corr_xx_xy_Landau + lambda_mu_nu[1][1] * lambda_mu_nu[1][2]/N;
	corr_xx_yy_Landau = corr_xx_yy_Landau + lambda_mu_nu[1][1] * lambda_mu_nu[2][2]/N;
	corr_xy_yy_Landau = corr_xy_yy_Landau + lambda_mu_nu[1][2] * lambda_mu_nu[2][2]/N;

	double xi_xx = 0.0, xi_xy = 0.0, xi_yx = 0.0, xi_yy = 0.0; // noise in density frame

		for (int mu = 0; mu < dim; mu++)
		{
		for (int nu = 0; nu < dim; nu++)
		{
		xi_xx = xi_xx + kappa_xx[mu][nu] * lambda_mu_nu[mu][nu];
		xi_xy = xi_xy + kappa_xy[mu][nu] * lambda_mu_nu[mu][nu];
		xi_yx = xi_yx + kappa_yx[mu][nu] * lambda_mu_nu[mu][nu];
		xi_yy = xi_yy + kappa_yy[mu][nu] * lambda_mu_nu[mu][nu];
		}
		}	


	corr_xx_xx = corr_xx_xx + xi_xx * xi_xx/N;
	corr_xy_xy = corr_xy_xy + xi_xy * xi_xy/N;
	corr_yy_yy = corr_yy_yy + xi_yy * xi_yy/N;
	corr_xx_xy = corr_xx_xy + xi_xx * xi_xy/N;
	corr_xx_yy = corr_xx_yy + xi_xx * xi_yy/N;
	corr_xy_yy = corr_xy_yy + xi_xy * xi_yy/N;

	noise_xx += xi_xx/N;

	}

	/**************************** for analytical result ************************************/

	double kappa_xxxx, kappa_xyxy, kappa_yyyy; // self correlations
       	double kappa_xxxy, kappa_xxyy, kappa_xyyy; // cross correlations	

	double LPL[d][d];

		for (int i = 0; i < d; i++)
		{
			for (int j = 0; j < d; j++)
			{
				LPL[i][j] = -1.0/gamma_v/gamma_v * cs2 * vel[i]*vel[j] + L[i][j]/(double)d;
			}
		}	

	double PLPL = (d - 1.0)/(double)d * cs2 * cs2 * v*v * v*v + 1.0/d * (1.0 - cs2 * v*v) * (1.0 - cs2 * v*v);

	kappa_xxxx = 2.0 * eta * 2.0*T/dt *( (L[0][0]*L[0][0] + L[0][0]*L[0][0])/2.0 - 1.0/PL *L[0][0] * LPL[0][0] - 1.0/PL * L[0][0] * LPL[0][0] +  PLPL/PL/PL * L[0][0] * L[0][0] );
	kappa_xyxy = 2.0 * eta * 2.0*T/dt *( (L[0][0]*L[1][1] + L[0][1]*L[1][0])/2.0 - 1.0/PL *L[0][1] * LPL[0][1] - 1.0/PL * L[0][1] * LPL[0][1] +  PLPL/PL/PL * L[0][1] * L[0][1] );
	kappa_yyyy = 2.0 * eta * 2.0*T/dt *( (L[1][1]*L[1][1] + L[1][1]*L[1][1])/2.0 - 1.0/PL *L[1][1] * LPL[1][1] - 1.0/PL * L[1][1] * LPL[1][1] +  PLPL/PL/PL * L[1][1] * L[1][1] );

	kappa_xxxy = 2.0 * eta * 2.0*T/dt *( (L[0][0]*L[0][1] + L[0][1]*L[0][0])/2.0 - 1.0/PL *L[0][0] * LPL[0][1] - 1.0/PL * L[0][1] * LPL[0][0] +  PLPL/PL/PL * L[0][0] * L[0][1] );
	kappa_xxyy = 2.0 * eta * 2.0*T/dt *( (L[0][1]*L[0][1] + L[0][1]*L[0][1])/2.0 - 1.0/PL *L[0][0] * LPL[1][1] - 1.0/PL * L[1][1] * LPL[0][0] +  PLPL/PL/PL * L[0][0] * L[1][1] );
	kappa_xyyy = 2.0 * eta * 2.0*T/dt *( (L[0][1]*L[1][1] + L[0][1]*L[1][1])/2.0 - 1.0/PL *L[0][1] * LPL[1][1] - 1.0/PL * L[1][1] * LPL[0][1] +  PLPL/PL/PL * L[0][1] * L[1][1] );

	double kappa = 1.0/pow(gamma_v, 4.0) * 2.0 * (d - 1.0)/d * T * eta * 1.0/pow(1.0 - cs2 * v*v, 2.0);

	//cout<<endl;

	/*
	cout<<"v: "<<v<<" kappa_xxxx: "<<kappa_xxxx<<" "<<corr_xx_xx<<" "<<noise_xx<<endl;
	cout<<"v: "<<v<<" kappa_xyxy: "<<kappa_xyxy<<" "<<corr_xy_xy<<endl;
	cout<<"v: "<<v<<" kappa_yyyy: "<<kappa_yyyy<<" "<<corr_yy_yy<<endl;
	
	cout<<"v: "<<v<<" kappa_xxxy: "<<kappa_xxxy<<" "<<corr_xx_xy<<endl;
	cout<<"v: "<<v<<" kappa_xxyy: "<<kappa_xxyy<<" "<<corr_xx_yy<<endl;
	cout<<"v: "<<v<<" kappa_xyyy: "<<kappa_xyyy<<" "<<corr_xy_yy<<endl;
	*/

	double kappa_xxxx_Landau, kappa_xyxy_Landau, kappa_yyyy_Landau; // self correlations Landau
       	double kappa_xxxy_Landau, kappa_xxyy_Landau, kappa_xyyy_Landau; // cross correlations Landau	
	
	kappa_xxxx_Landau = 2.0 * eta * 2.0*T/dt * ( (Delta_lower_lower[1][1]*Delta_lower_lower[1][1] + Delta_lower_lower[1][1]*Delta_lower_lower[1][1])/2.0 - 1.0/d * Delta_lower_lower[1][1] * Delta_lower_lower[1][1] );
	kappa_xyxy_Landau = 2.0 * eta * 2.0*T/dt * ( (Delta_lower_lower[1][1]*Delta_lower_lower[2][2] + Delta_lower_lower[1][2]*Delta_lower_lower[2][1])/2.0 - 1.0/d * Delta_lower_lower[1][2] * Delta_lower_lower[1][2] );
	kappa_yyyy_Landau = 2.0 * eta * 2.0*T/dt * ( (Delta_lower_lower[2][2]*Delta_lower_lower[2][2] + Delta_lower_lower[2][2]*Delta_lower_lower[2][2])/2.0 - 1.0/d * Delta_lower_lower[2][2] * Delta_lower_lower[2][2] );

	kappa_xxxy_Landau = 2.0 * eta * 2.0*T/dt * ( (Delta_lower_lower[1][1]*Delta_lower_lower[1][2] + Delta_lower_lower[1][2]*Delta_lower_lower[1][1])/2.0 - 1.0/d * Delta_lower_lower[1][1] * Delta_lower_lower[1][2] );
	kappa_xxyy_Landau = 2.0 * eta * 2.0*T/dt * ( (Delta_lower_lower[1][2]*Delta_lower_lower[1][2] + Delta_lower_lower[1][2]*Delta_lower_lower[1][2])/2.0 - 1.0/d * Delta_lower_lower[1][1] * Delta_lower_lower[2][2] );
	kappa_xyyy_Landau = 2.0 * eta * 2.0*T/dt * ( (Delta_lower_lower[1][2]*Delta_lower_lower[2][2] + Delta_lower_lower[1][2]*Delta_lower_lower[2][2])/2.0 - 1.0/d * Delta_lower_lower[1][2] * Delta_lower_lower[2][2] );

	cout<<"v: "<<v<<" phi: "<<phi<<endl;

	if (count % 2 == 0)
	{
	file<<v<<" "<<phi<<" "<<corr_xx_xx<<" "<<corr_xy_xy<<" "<<corr_yy_yy<<" "<<corr_xx_xy<<" "<<corr_xx_yy<<" "<<corr_xy_yy<<" "<<endl;
       	}

       	file2<<v<<" "<<phi<<" "<<kappa_xxxx<<" "<<kappa_xyxy<<" "<<kappa_yyyy<<" "<<kappa_xxxy<<" "<<kappa_xxyy<<" "<<kappa_xyyy<<" "<<endl;

	if (count % 1 == 0)
	{
	file3<<v<<" "<<phi<<" "<<corr_xx_xx_Landau<<" "<<corr_xy_xy_Landau<<" "<<corr_yy_yy_Landau<<" "<<corr_xx_xy_Landau<<" "<<corr_xx_yy_Landau<<" "<<corr_xy_yy_Landau<<" "<<endl;
       	}

	if (count % 1 == 0)
	{
	file4<<v<<" "<<phi<<" "<<kappa_xxxx_Landau<<" "<<kappa_xyxy_Landau<<" "<<kappa_yyyy_Landau<<" "<<kappa_xxxy_Landau<<" "<<kappa_xxyy_Landau<<" "<<kappa_xyyy_Landau<<" "<<endl;
       	}

	v = v + 0.01;
	
	count++;

	}while (v<=0.9999);

}


