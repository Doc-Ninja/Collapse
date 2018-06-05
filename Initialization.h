#ifndef INITIALIZATION_H_INCLUDED
#define INITIALIZATION_H_INCLUDED

#define _USE_MATH_DEFINES
#include <math.h>
#include "Parameters.h"
#include "RK4_Evolution.h"



void initialize_x(double*x);
void initialize_vars(double *x, double *A, double *delta, double *Phi, double *Pi, double* t);
void initialization_1(double *x, double *A, double *delta, double *Phi, double *Pi, double* t);
void initialization_custom(double *x, double *A, double *delta, double *Phi, double *Pi, double* t);


#define ERR1(e) {printf("Error: %s\n", nc_strerror(e)); getchar(); exit(EXIT_FAILURE);}


#endif // INITIALIZATION_H_INCLUDED
