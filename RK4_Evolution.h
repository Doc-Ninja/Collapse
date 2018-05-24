#ifndef RK4_EVOLUTION_H_INCLUDED
#define RK4_EVOLUTION_H_INCLUDED

#include "Parameters.h"
#include "Equations.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void Delta_Solver(double* x, double* Phi, double* Pi, double* delta);
void A_Solver(double* x, double* Phi, double* Pi, double* A);
double* Phi_dot_FULL(double* x, double* A, double* delta, double* Pi);
double* Pi_dot_FULL(double* x, double* A, double* delta, double* Phi);
void evolve(double *x, double *A, double *delta, double *Phi, double *Pi, double dt);
int minpos(double *a);
double dt_cal(double h, double* A, double* delta);
bool Horizon_con(double* A);



#endif // RK4_EVOLUTION_H_INCLUDED
