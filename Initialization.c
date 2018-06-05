#include "Initialization.h"
#include "netcdf.h"

//function to initialize all the variables
void initialize_vars(double *x, double *A, double *delta, double *Phi, double *Pi, double* t) {
	initialize_x(x);
	if (initialization == 0)
		initialization_custom(x, A, delta, Phi, Pi, t);
	else if (initialization == 1)
		initialization_1(x, A, delta, Phi, Pi, t);
}

//function to initialize the coordinate array
void initialize_x(double*x) {
	int i;
	for (i = 0; i < SIZE; i++) {
		x[i] = (i*M_PI) / (2.0*(SIZE - 1));
	}
}

//function to initialize the fields arrays with a gaussian in Pi
void initialization_1(double *x, double *A, double *delta, double *Phi, double *Pi, double* t) {
	int i;
	for (i = 0; i < SIZE; i++) {
		Phi[i] = 0.0;
		Pi[i] = 2.0*eps*exp(-(4.0*pow(tan(x[i]), 2.0)) / pow(M_PI*sigma, 2)) / M_PI;
	}
	A_Solver(x, Phi, Pi, A);
	Delta_Solver(x, Phi, Pi, delta);
	(*t) = 0.0;
}

//function to input a custom starting setup, in case of grid size not matching between file and simulation settings, it will interpolate (linearly) or subsample accordingly
void initialization_custom(double *x, double *A, double *delta, double *Phi, double *Pi, double *t) {
	//file and data pointer
	int ncid, size_varid, A_varid, delta_varid, Phi_varid, Pi_varid, t_varid;

	//size storage
	int SIZE_f;

	//error handling
	int retval;

	//file opening
	if ((retval = nc_open("Input/start.nc", NC_NOWRITE, &ncid)))
		ERR1(retval);

	//data pointers assignement
	if ((retval = nc_inq_varid(ncid, "time", &t_varid)))
		ERR1(retval);
	if ((retval = nc_inq_varid(ncid, "A", &A_varid)))
		ERR1(retval);
	if ((retval = nc_inq_varid(ncid, "delta", &delta_varid)))
		ERR1(retval);
	if ((retval = nc_inq_varid(ncid, "PHI", &Phi_varid)))
		ERR1(retval);
	if ((retval = nc_inq_varid(ncid, "PI", &Pi_varid)))
		ERR1(retval);
	if ((retval = nc_inq_varid(ncid, "SIZE", &size_varid)))
		ERR1(retval);

	//Size acquisition
	if ((retval = nc_get_var_int(ncid, size_varid, &SIZE_f)))
		ERR1(retval);

	//Comparison of sizes and fields arrays population
	if (SIZE_f == SIZE) {
		//filling of the fields arrays
		if ((retval = nc_get_var_double(ncid, A_varid, A)))
			ERR1(retval);
		if ((retval = nc_get_var_double(ncid, delta_varid, delta)))
			ERR1(retval);
		if ((retval = nc_get_var_double(ncid, Phi_varid, Phi)))
			ERR1(retval);
		if ((retval = nc_get_var_double(ncid, Pi_varid, Pi)))
			ERR1(retval);
	}
	else {
		//allocation of temporary arrays for the variables
		double* A_temp = (double*)calloc(SIZE_f, sizeof(double));
		double* delta_temp = (double*)calloc(SIZE_f, sizeof(double));
		double* Phi_temp = (double*)calloc(SIZE_f, sizeof(double));
		double* Pi_temp = (double*)calloc(SIZE_f, sizeof(double));

		//Acquisition of data
		if ((retval = nc_get_var_double(ncid, A_varid, A_temp)))
			ERR1(retval);
		if ((retval = nc_get_var_double(ncid, delta_varid, delta_temp)))
			ERR1(retval);
		if ((retval = nc_get_var_double(ncid, Phi_varid, Phi_temp)))
			ERR1(retval);
		if ((retval = nc_get_var_double(ncid, Pi_varid, Pi_temp)))
			ERR1(retval);

		if (SIZE < SIZE_f) {
			//subsampling of temp arrays
			int i;
			int p;
			for (i = 0; i < SIZE; i++) {
				p = (i*(SIZE_f - 1)) / (SIZE - 1);
				A[i] = A_temp[p];
				Phi[i] = Phi_temp[p];
				Pi[i] = Pi_temp[p];
				delta[i] = delta_temp[p];
			}
		}
		else {
			//indexes for the 2 cycles and for the internal cycle stop
			int i, p, q, r;
			//doubles for the interpolaton coefficients
			double a, b, c;
			//interpolation
			for (p = 0; p < SIZE_f - 1; p++) {
				q = (p*(SIZE - 1)) / (SIZE_f - 1);
				r = ((p + 1)*(SIZE - 1)) / (SIZE_f - 1);
				c = r - q;
				for (i = q; i < r; i++) {
					a = r - i;
					b = i - q;
					A[i] = (A_temp[p] * a + A_temp[p + 1] * b) / c;
					delta[i] = (delta_temp[p] * a + delta_temp[p + 1] * b) / c;
					Phi[i] = (Phi_temp[p] * a + Phi_temp[p + 1] * b) / c;
					Pi[i] = (Pi_temp[p] * a + Pi_temp[p + 1] * b) / c;
				}
			}
			i = SIZE - 1;
			//last spot special treatment
			A[i] = A_temp[p];
			delta[i] = delta_temp[p];
			Phi[i] = Phi_temp[p];
			Pi[i] = Pi_temp[p];

		}
		//deallocation fo the temporary arrays
		free(A_temp);
		free(delta_temp);
		free(Phi_temp);
		free(Pi_temp);
	}
	//time acquisition
	if ((retval = nc_get_var_double(ncid, t_varid, t)))
		ERR1(retval);

	//file closure
	if ((retval = nc_close(ncid)))
		ERR1(retval);


}
