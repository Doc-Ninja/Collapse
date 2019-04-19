/**
This file contains the equations to evolve the system
*/

#include "Equations.h"

//function to compute the coefficients for the derivative in Phi dot (A exp(-delta) Pi)
double Phi_dot_PRE(double A, double delta, double Pi) {
	double PRE = A * exp(-delta)*Pi;
	return PRE;
}

//Functions to compute Phi dot, 4-point precision, with special stencils near the boundaries
double Phi_dot(double h, double LL, double L, double R, double RR) {
	double PHI_dot = (LL - 8.0*L + 8.0*R - RR) / (12.0*h);
	return PHI_dot;
}

double Phi_dot_LB(double h, double L, double C, double R, double RR) {
	double PHI_dot = (-2.0*L - 3.0*C + 6.0*R - 1.0*RR) / (6.0*h);
	return PHI_dot;
}

double Phi_dot_LLB(double h, double C, double R, double RR, double RRR) {
	//double PHI_dot = (-11.0*C + 18.0*R - 9.0*RR + 2.0*RRR) / (6.0*h);
  double PHI_dot=0; 
	return PHI_dot;
}

double Phi_dot_RB(double h, double LL, double L, double C, double R) {
	double PHI_dot = (1.0*LL - 6.0*L + 3.0*C + 2.0*R) / (6.0*h);
	return PHI_dot;
}

double Phi_dot_RRB(double h, double LLL, double LL, double L, double C) {
	//double PHI_dot= (-2.0*LLL +9.0*LL -18.0*L +11.0*C)/(6.0*h);
	double PHI_dot = 0.0;
	return PHI_dot;
}

//function to compute the coefficients for the derivative in Pi dot (tgx^d-1 A exp-delta Phi)
double Pi_dot_PRE(double x, double A, double delta, double Phi) {
	double PRE = pow(tan(x), d - 1.0)*A*exp(-delta)*Phi;
	return PRE;
}

//Functions to compute Pi dot, 4-point precision, with special stencils near the boundaries
double Pi_dot(double h, double LL, double L, double R, double RR, double x) {
	double PI_dot = pow(tan(x), 1.0 - d)*(LL - 8.0*L + 8.0*R - RR) / (12.0*h);
	return PI_dot;
}

double Pi_dot_LB(double h, double L, double C, double R, double RR, double x) {
	double PI_dot = pow(tan(x), 1.0 - d)*(-2.0*L - 3.0*C + 6.0*R - 1.0*RR) / (6.0*h);
	return PI_dot;
}
// this function computes Pi dot in the origin based on the fields in the following positions
double Pi_dot_LLB(double h, double C, double R, double RR, double RRR, double x) {
	double PI_dot = (18.0*R - 9.0*RR + 2.0*RRR) / 11.0;
	return PI_dot;
}

double Pi_dot_RB(double h, double LL, double L, double C, double R, double x) {
	double PI_dot = pow(tan(x), 1.0 - d)*(1.0*LL - 6.0*L + 3.0*C + 2.0*R) / (6.0*h);
	return PI_dot;
}

double Pi_dot_RRB(double h, double LLL, double LL, double L, double C, double x) {
	double PI_dot = 0.0;
	return PI_dot;
}

//function to compute A prime, no derivatives involved
double A_prime(double x, double A, double Phi, double Pi) {
	/*double A_PRIME = -sin(2.0*x)*A*(pow(Phi, 2.0) + pow(Pi, 2.0)) / 2.0;
	if (A != 1.0)
		A_PRIME += 2.0*(d - 2.0 + 2.0*pow(sin(x), 2.0)) * (1.0 - A) / (sin(2.0*x));*/
	double A_PRIME = - sin(x)*cos(x)*A*(pow(Phi, 2.0) + pow(Pi, 2.0));
	if (A != 1.0)
		A_PRIME += (d - 2.0 + 2 * pow(sin(x), 2.0))*(1 - A) / (sin(x)*cos(x));
	return A_PRIME;
}

//function to compute delta prime, no derivatives involved
double delta_prime(double x, double Phi, double Pi) {
	//double delta_PRIME = -sin(2.0*x)*(pow(Phi, 2.0) + pow(Pi, 2.0)) / 2.0;
	double delta_PRIME = -sin(x)*cos(x)*(pow(Phi, 2.0) + pow(Pi, 2.0));
	return delta_PRIME;
}

//functions for the 2 parts of the constraint
double Con1(double A_future, double A_past, double dt) {
	double A_dot = (A_future - A_past) / dt;
	return A_dot;
}
double Con2(double x, double A, double delta, double Phi, double Pi) {
	double A_sq_pt = sin(2.0*x)*pow(A, 2.0)*exp(-delta)*Phi*Pi;
	return A_sq_pt;
}
