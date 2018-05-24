#include "Initialization.h"

//function to initialize the coordinate array
double* initialize_x(){
    int i;
    static double x[SIZE];
    for(i=0; i<SIZE; i++){
        x[i]=(i*M_PI)/(2.0*(SIZE-1));
    }
    return x;
}

//function to initialize the fields arrays with a gaussian in Pi
void initialization_1(double *x, double *A, double *delta, double *Phi, double *Pi){
    int i;
    for(i=0;i<SIZE;i++){
        Phi[i]=0.0;
        Pi[i]=2.0*eps*exp(-(4.0*pow(tan(x[i]),2.0))/pow(M_PI*sigma,2))/M_PI;
    }
    A_Solver(x, Phi, Pi, A);
    Delta_Solver(x, Phi, Pi, delta);
}

//function to input a custom starting setup
void initialization_custom(double *x, double *A, double *delta, double *Phi, double *Pi){
}
