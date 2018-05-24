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

// initialization parameters
#define eps 45.0
#define sigma 1/16.0



//Simulation parameters
#define SIZE 128
#define A_HORIZON 0.01
#define STEP_LIMIT 5
#define TIME_STRIDE 1
/****************************
netCDF parameters and defines
****************************/

//error handling
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}

//dimension defines
#define NDIMS 2


#endif // PARAMETERS_H_INCLUDED
