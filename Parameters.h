#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include <math.h>
#include <stdio.h>

// definition of boolean type
typedef int bool;
#define true 1
#define false 0

// dimension
#define d 3.0

#define TIME_STRIDE 100


int type;
double eps;
double sigma;
const char* file;

#define SIZE 2048
//int SIZE;
double A_HORIZON;
int STEP_LIMIT;

void load_param();
void f_init();
void p_init();


/****************************
netCDF parameters and defines
****************************/



//dimension defines
#define NDIMS 2





#endif // PARAMETERS_H_INCLUDED
