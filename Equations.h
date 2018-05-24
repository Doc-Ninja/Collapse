#ifndef EQUATIONS_H_INCLUDED
#define EQUATIONS_H_INCLUDED

#include <math.h>


#include "Parameters.h"
//Phi dot computation block
double Phi_dot_PRE (double A,double delta,double Pi);
double Phi_dot (double h, double LL, double L, double R, double RR);
double Phi_dot_LB (double h, double L, double C, double R, double RR);
double Phi_dot_RB (double h, double LL, double L, double C, double R);
double Phi_dot_LLB (double h, double C, double R, double RR, double RRR);
double Phi_dot_RRB (double h, double LLL, double LL, double L, double C);

//Pi dot computation block
double Pi_dot_PRE (double x, double A, double delta, double Phi);
double Pi_dot (double h, double LL, double L, double R, double RR, double x);
double Pi_dot_LB (double h, double L, double C, double R, double RR, double x);
double Pi_dot_LLB (double h, double C, double R, double RR, double RRR, double x);
double Pi_dot_RB (double h, double LL, double L, double C, double R, double x);
double Pi_dot_RRB (double h, double LLL, double LL, double L, double C, double x);

//Delta computation function
double A_prime (double x, double A, double Phi, double Pi);

//A computation Function
double delta_prime (double x, double Phi, double Pi);

#endif // EQUATIONS_H_INCLUDED
