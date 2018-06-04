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
	double t_past = 0.0;
	double dt = 0.0;
    double* x=initialize_x();
    double h=x[1];
    double A[SIZE];
    double delta[SIZE];
    double Phi[SIZE];
    double Pi[SIZE];
	double Constr[3][SIZE] = { 0 };
    initialization_1(x, A, delta, Phi, Pi);
	initialize_files(x);

	//first data collection
	data_stamp(t_step, t, x, A, delta, Phi, Pi, Constr, t_past);
	
	//main cycle
	for (t_step = 1; (t_step <= STEP_LIMIT)&&(Horizon_con(A)); t_step++) {
		Con_past(A, dt, t, &t_past);
		Con2_FULL(x, A, delta, Phi, Pi, Constr);
		dt = dt_cal(h, A, delta);
		t += dt;
		printf("ITERATION: %d\n", t_step);
		evolve(x, A, delta, Phi, Pi, dt);
		Con1_FULL(A, dt, Constr);
		Con3_FULL(Constr);
		data_stamp(t_step, t, x, A, delta, Phi, Pi, Constr, t_past);		
	}
	data_stamp_final(t, x, A, delta, Phi, Pi);
	//closing the open files
	close_all();
	if (t_step > STEP_LIMIT)
	  printf("CYCLE LIMIT!");
	if (!(Horizon_con(A))){
		FILE* report = fopen("Output/Report.txt", "w");
		int hpos_n;
		double hpos_x;
		hpos_n = minpos(A);
		hpos_x = x[hpos_n];
		fprintf(report, "Horizon found at x[%d] = %f \n Cycles elapsed: %d", hpos_n, hpos_x, t_step);
		fclose(report);
	}
	getchar();
    return 0;
}
