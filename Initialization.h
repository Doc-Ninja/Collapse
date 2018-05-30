#ifndef INITIALIZATION_H_INCLUDED
#define INITIALIZATION_H_INCLUDED

#define _USE_MATH_DEFINES
#include <math.h>
#include "Parameters.h"
#include "RK4_Evolution.h"

double* initialize_x();
void initialization_1(double *x, double *A, double *delta, double *Phi, double *Pi);
void initialization_custom(double *x, double *A, double *delta, double *Phi, double *Pi);
void initialize_fields(double *x, double *A, double *delta, double *Phi, double *Pi);



#endif // INITIALIZATION_H_INCLUDED
