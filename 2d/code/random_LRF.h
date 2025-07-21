#include<random>
using namespace std;


double generateRandom_normal() 
{
    static thread_local std::mt19937 generator2{ std::random_device{}() };
    std::normal_distribution<double> normal(0, 1);
    return normal(generator2);
}


/*
double generateRandom_normal() 
{
    std::random_device rd;
    std::mt19937 generator2(rd());

    std::normal_distribution<double> normal(0, 1);
    return normal(generator2);
}
*/

void random_numbers_uncorrelated(double *g_xx, double *g_xy, double *g_yy, double A)
{

/** generates Gaussian random numbers ***/

*g_xx = sqrt(A) * generateRandom_normal();
*g_xy = sqrt(A/2.0)  * generateRandom_normal();
*g_yy = sqrt(A) * generateRandom_normal();

}

void random_numbers_LRF(double *lambda_xx, double *lambda_xy, double *lambda_yy, double A)
{

double gamma_xx, gamma_xy, gamma_yy;

random_numbers_uncorrelated(&gamma_xx, &gamma_xy, &gamma_yy, A);

double trace = gamma_xx + gamma_yy;

*lambda_xx = gamma_xx - trace/(double)d;
*lambda_xy = gamma_xy;
*lambda_yy = gamma_yy - trace/(double)d;

}
