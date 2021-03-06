#include "RK4_Evolution.h"


//Solver for delta, since the derivative does not depend from the function, it becomes a Simpson integral
void Delta_Solver(double* x, double* Phi, double* Pi, double* delta) {
	int i;
	double h = x[1];
	delta[0] = 0.0;
	double k1, k2, k3;
	for (i = 0; i < SIZE - 1; i++) {
		k1 = h * delta_prime(x[i], Phi[i], Pi[i]);
		k2 = h * delta_prime((x[i] + x[i + 1]) / 2.0, (Phi[i] + Phi[i + 1]) / 2.0, (Pi[i] + Pi[i + 1]) / 2.0);
		k3 = h * delta_prime(x[i + 1], Phi[i + 1], Pi[i + 1]);
		delta[i + 1] = delta[i] + (k1 + 4.0*k2 + k3) / 6.0;
	}
}

//Solver for A
void A_Solver(double* x, double* Phi, double* Pi, double* A) {
	int i;
	A[0] = 1.0;
	double k1, k2, k3, k4;
	double h = x[1];
	for (i = 0; i < SIZE - 2; i++) {
		k1 = h * A_prime(x[i], A[i], Phi[i], Pi[i]);
		k2 = h * A_prime((x[i] + x[i + 1]) / 2.0, A[i] + k1 / 2.0, (Phi[i] + Phi[i + 1]) / 2.0, (Pi[i] + Pi[i + 1]) / 2.0);
		k3 = h * A_prime((x[i] + x[i + 1]) / 2.0, A[i] + k2 / 2.0, (Phi[i] + Phi[i + 1]) / 2.0, (Pi[i] + Pi[i + 1]) / 2.0);
		k4 = h * A_prime(x[i + 1], A[i] + k3, Phi[i + 1], Pi[i + 1]);
		A[i + 1] = A[i] + (k1 + 2.0f * k2 + 2.0f * k3 + k4) / 6.0f;
	}
	A[SIZE - 1] = 1.0; //check if it is possible to do without this
}

//Functions to compute Phi_dot and Pi_dot over all the array
double* Phi_dot_FULL(double* x, double* A, double* delta, double* Pi) {
	static double Phi_DOT[SIZE];
	double PRE[SIZE];
	double h = x[1];
	int i;
	for (i = 0; i < SIZE; i++) {
		PRE[i] = Phi_dot_PRE(A[i], delta[i], Pi[i]);
	}
	Phi_DOT[0] = Phi_dot_LLB(h, PRE[0], PRE[1], PRE[2], PRE[3]);
	Phi_DOT[1] = Phi_dot_LB(h, PRE[0], PRE[1], PRE[2], PRE[3]);
	for (i = 2; i < SIZE - 2; i++) {
		Phi_DOT[i] = Phi_dot(h, PRE[i - 2], PRE[i - 1], PRE[i + 1], PRE[i + 2]);
	}
	Phi_DOT[SIZE - 2] = Phi_dot_RB(h, PRE[SIZE - 4], PRE[SIZE - 3], PRE[SIZE - 2], PRE[SIZE - 1]);
	Phi_DOT[SIZE - 1] = Phi_dot_RRB(h, PRE[SIZE - 4], PRE[SIZE - 3], PRE[SIZE - 2], PRE[SIZE - 1]);
	return Phi_DOT;
}
double* Pi_dot_FULL(double* x, double* A, double* delta, double* Phi) {
	static double Pi_DOT[SIZE];
	double PRE[SIZE];
	double h = x[1];
	int i;
	for (i = 0; i < SIZE; i++) {
		PRE[i] = Pi_dot_PRE(x[i], A[i], delta[i], Phi[i]);
	}
	Pi_DOT[1] = Pi_dot_LB(h, PRE[0], PRE[1], PRE[2], PRE[3], x[1]);
	for (i = 2; i < SIZE - 2; i++) {
		Pi_DOT[i] = Pi_dot(h, PRE[i - 2], PRE[i - 1], PRE[i + 1], PRE[i + 2], x[i]);
	}
	Pi_DOT[SIZE - 2] = Pi_dot_RB(h, PRE[SIZE - 4], PRE[SIZE - 3], PRE[SIZE - 2], PRE[SIZE - 1], x[SIZE - 2]);
	Pi_DOT[SIZE - 1] = Pi_dot_RRB(h, PRE[SIZE - 4], PRE[SIZE - 3], PRE[SIZE - 2], PRE[SIZE - 1], x[SIZE - 1]);
	Pi_DOT[0] = Pi_dot_LLB(h, PRE[0], Pi_DOT[1], Pi_DOT[2], Pi_DOT[3], x[0]);

	return Pi_DOT;
}

//variables to store A and dt in the past
double A_past[SIZE] = { 0 };
double A_copy[SIZE] = { 0 };
double dt_past;

//Function to update the stored values for the Constrain computing
void Con_past(double * A, double dt, double t, double *t_past) {
	int i;
	for (i = 0; i < SIZE; i++) {
		A_past[i] = A_copy[i];
		A_copy[i] = A[i];
	}
	dt_past = dt;
	(*t_past) = t;
}

//functions to compute the 2 parts of the constraint over the entire array and their sum
void Con1_FULL(double* A, double dt, double **Con) {
	double h = dt + dt_past;
	int i;
	for (i = 0; i < SIZE; i++)
		Con[0][i] = Con1(A[i], A_past[i], h);
}
void Con2_FULL(double*x, double* A, double* delta, double* Phi, double* Pi, double **Con) {
	int i;
	for (i = 0; i < SIZE; i++)
		Con[1][i] = Con2(x[i], A[i], delta[i], Phi[i], Pi[i]);
}
void Con3_FULL(double **Con) {
	int i;
	for (i = 0; i < SIZE; i++)
		Con[2][i] = Con[0][i] + Con[1][i];
}


//function that evolves all the fields by 1 time step
void evolve(double *x, double *A, double *delta, double *Phi, double *Pi, double dt) {
	//declarations
	int i;
	double* tempPhi = calloc(SIZE, sizeof(double));
	double* tempPi = calloc(SIZE, sizeof(double));
	double* tempkPhi;
	double* tempkPi;
	double* kPhi1 = calloc(SIZE, sizeof(double));
	double* kPhi2 = calloc(SIZE, sizeof(double));
	double* kPhi3 = calloc(SIZE, sizeof(double));
	double* kPhi4 = calloc(SIZE, sizeof(double));
	double* kPi1 = calloc(SIZE, sizeof(double));
	double* kPi2 = calloc(SIZE, sizeof(double));
	double* kPi3 = calloc(SIZE, sizeof(double));
	double* kPi4 = calloc(SIZE, sizeof(double));

	//computation of the first RK4 step coefficients
	tempkPhi = Phi_dot_FULL(x, A, delta, Pi);
	tempkPi = Pi_dot_FULL(x, A, delta, Phi);
	for (i = 0; i < SIZE; i++) {
		kPhi1[i] = dt * tempkPhi[i];
		kPi1[i] = dt * tempkPi[i];
		tempPhi[i] = Phi[i] + kPhi1[i] * 0.5;
		tempPi[i] = Pi[i] + kPi1[i] * 0.5;
	}

	//Updating A and delta 
	A_Solver(x, tempPhi, tempPi, A);
	Delta_Solver(x, tempPhi, tempPi, delta);

	//computation of the second RK4 step coefficients
	tempkPhi = Phi_dot_FULL(x, A, delta, tempPi);
	tempkPi = Pi_dot_FULL(x, A, delta, tempPhi);
	for (i = 0; i < SIZE; i++) {
		kPhi2[i] = dt * tempkPhi[i];
		kPi2[i] = dt * tempkPi[i];
		tempPhi[i] = Phi[i] + kPhi2[i] * 0.5;
		tempPi[i] = Pi[i] + kPi2[i] * 0.5;
	}

	//Updating A and delta 
	A_Solver(x, tempPhi, tempPi, A);
	Delta_Solver(x, tempPhi, tempPi, delta);


	//computation of the third RK4 step coefficients
	tempkPhi = Phi_dot_FULL(x, A, delta, tempPi);
	tempkPi = Pi_dot_FULL(x, A, delta, tempPhi);

	for (i = 0; i < SIZE; i++) {
		kPhi3[i] = dt * tempkPhi[i];
		kPi3[i] = dt * tempkPi[i];
		tempPhi[i] = Phi[i] + kPhi3[i];
		tempPi[i] = Pi[i] + kPi3[i];
	}

	//Updating A and delta
	A_Solver(x, tempPhi, tempPi, A);
	Delta_Solver(x, tempPhi, tempPi, delta);

	//computation of the fourth RK4 step coefficients
	tempkPhi = Phi_dot_FULL(x, A, delta, tempPi);
	tempkPi = Pi_dot_FULL(x, A, delta, tempPhi);
	//updating Phi and Pi
	for (i = 0; i < SIZE; i++) {
		kPhi4[i] = dt * tempkPhi[i];
		kPi4[i] = dt * tempkPi[i];
		Phi[i] = Phi[i] + (kPhi1[i] + 2.0*kPhi2[i] + 2.0*kPhi3[i] + kPhi4[i]) / 6.0;
		Pi[i] = Pi[i] + (kPi1[i] + 2.0*kPi2[i] + 2.0*kPi3[i] + kPi4[i]) / 6.0;
	}

	//Updating A and delta
	A_Solver(x, Phi, Pi, A);
	Delta_Solver(x, Phi, Pi, delta);

	free(kPhi1);
	free(kPhi2);
	free(kPhi3);
	free(kPhi4);
	free(kPi1);
	free(kPi2);
	free(kPi3);
	free(kPi4);
	free(tempPhi);
	free(tempPi);


}

//function to find the minimum entry in an array
int minpos(double *a) {
	int pos = 0;
	int i;
	for (i = 1; i < SIZE; i++) {
		if (a[i] < a[pos])
			pos = i;
	}
	return pos;
}

//function to compute the time step
double dt_cal(double h, double* A, double* delta) {
	double temp_courant[SIZE];
	double dtp;
	int i;
	for (i = 0; i < SIZE; i++) {
		temp_courant[i] = fabs(exp(delta[i])*h / A[i]);
	}
	double courant_con = temp_courant[minpos(temp_courant)];
	dtp = C_COUR* courant_con;
	return dtp;
}

//function to check if the collapse has happened
bool Horizon_con(double* A) {
	bool b;
	double min_A = A[minpos(A)];
	if (min_A < A_HORIZON) 
		b = false;
	else
		b = true;
	return b;
}
