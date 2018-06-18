#ifndef DATA_COLLECTION_H_INCLUDED
#define DATA_COLLECTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#include "Parameters.h"

void initialize_files(double *x);
void initialize_bigfile(double *x);
void initialize_probe(double* x);
void initialize_con(double *x);
void initialize_extra(double* x);
void data_stamp(int t_step, double t, double* x, double* A, double* delta, double* Phi, double* Pi, double **Con, double t_past);
void data_stamp_bigfile(int t_step, double t, double* A, double* delta, double* Phi, double* Pi);
void data_stamp_probe(int t_step, double t, double* A, double* delta, double* Phi, double* Pi);
void data_stamp_con(int t_step, double **Con, double t_past);
void data_stamp_extra(int t_step, double t, double* A, double* delta, double* Phi, double* Pi);
void data_stamp_check(int t_step, double t, double* x, double* A, double* delta, double* Phi, double* Pi);
void data_stamp_final(double t, double* x, double* A, double* delta, double* Phi, double* Pi);
void close_all();
void close_bigfile();
void close_probe();
void close_con();
void close_extra();

#endif //DATA_COLLECTION_H_INCLUDED