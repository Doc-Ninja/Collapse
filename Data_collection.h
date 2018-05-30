#ifndef DATA_COLLECTION_H_INCLUDED
#define DATA_COLLECTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#include "Parameters.h"

//net CDF error handling
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}

void initialize_files(double *x);
void data_stamp(int t_step, double t, double* A, double* delta, double* Phi, double* Pi);
void close_all();

#endif //DATA_COLLECTION_H_INCLUDED