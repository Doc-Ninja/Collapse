#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Parameters.h"
#include "Equations.h"
#include "RK4_Evolution.h"
#include "Initialization.h"
#include "Data_collection.h"

int main(){
	//initializaion of parameters, variables, fields and files
	//extern int SIZE;
	extern double A_HORIZON;
	extern int STEP_LIMIT;
	int t_step = 0;
	double t = 0.0;
	double dt;
	double* x = initialize_x();
    double h=x[1];
	double* A = calloc(SIZE, sizeof(double));
    double* delta = calloc(SIZE, sizeof(double));
    double* Phi = calloc(SIZE, sizeof(double));
    double* Pi = calloc(SIZE, sizeof(double));
	load_param();
    initialize_fields(x, A, delta, Phi, Pi);
	initialize_files(x);

	//first data collection
	data_stamp(t_step, t, A, delta, Phi, Pi);

	//setup for the first cycle
	dt = dt_cal(h, A, delta);
	t += dt;
	
	//main cycle
	for (t_step = 1; (t_step <= STEP_LIMIT)&&(Horizon_con(A)); t_step++) {
		printf("ITERATION: %d\n", t_step);
		evolve(x, A, delta, Phi, Pi, dt);
		data_stamp(t_step, t, A, delta, Phi, Pi);
		dt = dt_cal(h, A, delta);
		t += dt;
	}

	//closing the open files
	close_all();
	//freeing the fields array positions
	free(x);
	free(A);
	free(delta);
	free(Phi);
	free(Pi);
	if (t_step > STEP_LIMIT)
	  printf("CYCLE LIMIT!");
	getchar();
    return 0;
}
