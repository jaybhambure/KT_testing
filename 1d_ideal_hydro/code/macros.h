int N = 100000;
double dt = 0.001; // perhaps much smaller than required

const int L = 25000;
double dx = 0.01;

const int Vol = L;

double theta = 1.0; // flux limiting parameter

int global_flag;
int global_N;

double cs = 1.0/sqrt(3.0);
