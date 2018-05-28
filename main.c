#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Parameters.h"
#include "Equations.h"
#include "RK4_Evolution.h"
#include "Initialization.h"
#include "Data_collection.h"

int main(){
	//initializaion of variables, fields and files
	int t_step = 0;
	double t = 0.0;
	double dt;
    double* x=initialize_x();
    double h=x[1];
    double A[SIZE];
    double delta[SIZE];
    double Phi[SIZE];
    double Pi[SIZE];
    initialization_1(x, A, delta, Phi, Pi);
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
	if (t_step > STEP_LIMIT)
	  printf("CYCLE LIMIT!");
	getchar();
    return 0;
}
