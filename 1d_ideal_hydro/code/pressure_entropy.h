#include <math.h>
using namespace std;

double pressure(double E, double M_x)
{
	double p = 1.0/3.0 * (-E + sqrt(4.0*E*E - 3.0*M_x*M_x));

	return p;
}


	double vel(double E, double M_x)
	{
	double epsilon = pow(10.0, -20.0);
	double result;	

	if (abs(E) < epsilon)
	result = 0.0;

	else
	{	
	result = M_x/(E + pressure(E, M_x));	
	}
	
	return result;
	}


