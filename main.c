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
	double t;
	double dt = 0.0;
	double *x = calloc(SIZE, sizeof(double));
    double *A = calloc(SIZE, sizeof(double));
    double *delta = calloc(SIZE, sizeof(double));
    double *Phi = calloc(SIZE, sizeof(double));
    double *Pi = calloc(SIZE, sizeof(double));
	//double Constr[3][SIZE] = { 0 };
	double **Constr = (double **)calloc(3, sizeof(double*));
	for (int i_c = 0; i_c < 3; i_c++) Constr[i_c] = (double *)calloc(SIZE, sizeof(double));
	
	initialize_vars(x, A, delta, Phi, Pi, &t);
	double h = x[1];
	double t_past = t;
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
		printf("Horizon Found!");
		FILE* report = fopen("Output/Report.txt", "w");
		int hpos_n;
		double hpos_x;
		hpos_n = minpos(A);
		hpos_x = x[hpos_n];
		fprintf(report, "Horizon found at x[%d] = %f \nCycles elapsed: %d\nStarting eps: %f", hpos_n, hpos_x, t_step, eps);
		fclose(report);
	}
	//getchar();
	free(x);
	free(A);
	free(delta);
	free(Pi);
	free(Phi);
	for (int i_c = 0; i_c < 3; i_c++) free(Constr[i_c]);
	free(Constr);
    return 0;
}
