#include <math.h>

double pressure(double E, double M_x, double M_y)
{

	//double p = 1.0/3.0 * (-E + sqrt(4.0*E*E - 3.0*M_x*M_x - 3.0*M_y*M_y));
	double p = 1.0/2.0/(double)d * (- (d-1.0)*E + sqrt( (d+1.0)*(d+1.0) *E*E - 4.0*d * M_x*M_x - 4.0*d * M_y*M_y) );

	return p;
}

double entropy_density(double E, double M_x, double M_y)
{

	double p = pressure(E, M_x, M_y);

	//double gamma_vel = sqrt( (E + p)/4.0/p );
	//double T = pow(M_PI*M_PI * p, 0.25); // p = T^4/pi^2, single component Boltzmann gas in 3d
	//double s = 4.0 * p/T *gamma_vel;
	
	double gamma_vel = sqrt( (E + p)/(d+1.0)/p );
	double T = pow(2 * M_PI * p, 1.0/3.0); // p = T^3/(2 pi), single component Boltzmann gas in 2d
	double s = (d+1.0) * p/T *gamma_vel;

	return s;

}
	

void get_T_v(double E, double M_x, double M_y, double *T, double *vx, double *vy)
{

	double p = pressure(E, M_x, M_y);

	//*T = pow(M_PI*M_PI*p, 0.25);
	
	*T = pow(2 * M_PI * p, 1.0/3.0);
	*vx = M_x/(E + p);
	*vy = M_y/(E + p);

}

