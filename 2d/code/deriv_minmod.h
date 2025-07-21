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
	


		double deriv(double rho_i, double rho_ip, double rho_im, double delta_a)
		{
		double forward_deriv = (rho_ip - rho_i)/delta_a; // delta_a : dx, dy
		double backward_deriv = (rho_i - rho_im)/delta_a;
		double centred_deriv = (rho_ip - rho_im)/2.0/delta_a;

		return ( minmod(theta * forward_deriv, centred_deriv, theta * backward_deriv) );
		}

