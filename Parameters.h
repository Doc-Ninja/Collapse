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
#define eps 40.0
#define sigma 1/16.0
//0 for custom, 1 for gaussian on Pi
#define initialization 1

//Simulation parameters
#define N 12
#define SIZE ((1<<N)+1)
#define A_HORIZON 0.001
#define STEP_LIMIT 1
#define C_COUR 0.25

/*********************************************
netCDF/data collection  parameters and defines
*********************************************/



//dimension defines
#define NDIMS 2

//Enabling/disabling files
#define EXTRA_PROBE false
#define PROBE false
#define BIGFILE false
#define CHECKPOINT true
#define CONSTRAINT false

//ALL file Time and space strides
#define TIME_STRIDE 100
#define SPACE_STRIDE 1
#define B_SIZE ((SIZE-1)/SPACE_STRIDE)+1

//Time Stride for probe files
#define PROBE_STRIDE 1

//Number of probe files
#define N_PROBE 4

//position of extra_probe
#define X_PROBE 5

//space and time stride for the constraint file
#define C_TIME_STRIDE 100
#define C_SPACE_STRIDE 1
#define C_SIZE  ((SIZE-1)/C_SPACE_STRIDE)+1

//Time steps beetween checkpoints
#define CHECKSTEP 100000


//error handling
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); getchar(); close_all(); exit(EXIT_FAILURE);}


#endif // PARAMETERS_H_INCLUDED
