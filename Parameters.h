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
#define eps 55.0
#define sigma 1/16.0

//Simulation parameters
#define SIZE 1024
#define A_HORIZON 0.1
#define STEP_LIMIT 1000

/*********************************************
netCDF/data collection  parameters and defines
*********************************************/



//dimension defines
#define NDIMS 2

//Enabling/disabling files
#define EXTRA_PROBE true
#define PROBE true
#define BIGFILE true
#define CHECKPOINT true
#define CONSTRAINT true

//ALL file Time and space strides
#define TIME_STRIDE 100
#define SPACE_STRIDE 100
#define B_SIZE ((SIZE-1)/SPACE_STRIDE)+1

//Time Stride for probe files
#define PROBE_STRIDE 100

//Number of probe files
#define N_PROBE 5

//position of extra_probe
#define X_PROBE 5

//space and time stride for the constraint file
#define C_TIME_STRIDE 100
#define C_SPACE_STRIDE 100
#define C_SIZE  ((SIZE-1)/C_SPACE_STRIDE)+1

//Time steps beetween checkpoints
#define CHECKSTEP 100


//error handling
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); getchar(); close_all(); exit(EXIT_FAILURE);}


#endif // PARAMETERS_H_INCLUDED
