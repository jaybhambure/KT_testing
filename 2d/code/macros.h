int N = 20000000;
const int L = 4;
const int Vol = L*L;

float dt = 0.01;
float dx = 1;
float dy = 1;

//float dx = 0.5;
//float dy = 0.5;

float dV = dx * dy;

float A = 10.0; // in GeV^4 (initial conditions for test problem)
float delta = 0.12;

float theta = 1.0; // flux limiter 

const int d = 2;
const int dim = d+1;
double eta = 0.1;

float cs = 1.0/sqrt(d);

int regulator = 1;

long int counter_no_update = 0;
long int counter_regulation = 0;
